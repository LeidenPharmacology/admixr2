# =============================================================================
# Inference under misspecification: TIC, and the corrected likelihood-ratio test
# =============================================================================
#
# Eq. (1) is a QUASI-likelihood -- the exact density of a model known to be
# false, since the subjects behind a summary are a mixture over (a, b) and are
# normal only where f is linear in both. Estimates stay consistent (the score has
# zero expectation at Psi for any weight), but three things that all descend from
# the information equality H = J do not:
#
#   standard errors      -> covMethod = "r,s"      (R/adfweight.R)
#   the LRT reference    -> anova.admFit()         (here)
#   the AIC penalty 2p   -> TIC, tr(H^-1 J)        (here)
#
# One matrix product H^-1 J repairs all three. Under correct specification
# J = 2H, every eigenvalue below is 1, tr(H^-1 J) = 2p, and each collapses to
# the textbook form exactly rather than approximately.
#
# See algorithm/adf/HANDOFF-INFERENCE.md and DERIVATION-DOFV.md.

# Ruben's series for P(sum_i lambda_i chi2_1 < x).
#
# Implemented here rather than taking a dependency on CompQuadForm: it is 25
# lines and that package would be a new Imports for them. Imhof's inversion was
# tried and rejected -- its integrand is 0/0 at u = 0 and decays like
# u^(-1-q/2) with an oscillating numerator, so integrate() returns NA at q = 1.
#
# Ruben converges geometrically at small q and is exactly self-checking: with
# all lambda equal every g_r is zero, every correction term vanishes, and it
# returns pchisq(x/lambda, q) identically.
.admRubenP <- function(x, lambda, maxit = 2000L, tol = 1e-10) {
  q  <- length(lambda)
  b  <- min(lambda)
  rl <- 1 - b / lambda                       # zero when all lambda are equal
  a  <- numeric(maxit + 1L)
  a[1L] <- prod(sqrt(b / lambda))
  gr <- vapply(seq_len(maxit), function(r) sum(rl^r), numeric(1))
  p  <- a[1L] * stats::pchisq(x / b, q)
  for (k in seq_len(maxit)) {
    a[k + 1L] <- sum(gr[seq_len(k)] * a[k:1L]) / (2 * k)
    trm <- a[k + 1L] * stats::pchisq(x / b, q + 2 * k)
    p   <- p + trm
    if (is.finite(trm) && abs(trm) < tol * max(p, 1e-300)) break
  }
  min(max(p, 0), 1)
}

# The (H, J) pair a "r,s" fit stored, or NULL. One accessor, so anova.admFit()
# and the TIC block cannot disagree about where it lives -- adirmc uses
# `adirmcExtra` where the other three use `admExtra` (the same split
# plot.admFit resolves).
.admSandwichParts <- function(fit) {
  e <- tryCatch(fit$env, error = function(e) NULL)
  if (is.null(e)) return(NULL)
  x <- (e$admExtra %||% e$adirmcExtra)$sandwich
  if (is.null(x) || is.null(x$H) || is.null(x$J)) return(NULL)
  x
}

# The optimizer-scale parameter names a fit was built on, sandwich or not.
.admFitParNames <- function(fit) {
  e <- tryCatch(fit$env, error = function(e) NULL)
  if (is.null(e)) return(NULL)
  (e$admExtra %||% e$adirmcExtra)$par_names
}

# TIC and the effective parameter count.
#
# CONVENTION CHECK, and it is the one that catches a factor of two: under
# correct specification J = 2H, so tr(H^-1 J) = 2p, p_eff = p, and
# TIC = objective + 2p = AIC exactly. Any implementation failing that is wrong.
#
# The penalty is LARGER than AIC's wherever the model is misspecified (measured
# p_eff = 10.73 against p = 7, a 7.5-unit shift), so AIC systematically
# OVER-SELECTS here -- squarely inside the range where covariate decisions turn.
.admTICStats <- function(sw, objective) {
  if (is.null(sw)) return(NULL)
  Hi <- tryCatch(solve(sw$H), error = function(e) NULL)
  if (is.null(Hi)) return(NULL)
  tr <- sum(diag(Hi %*% sw$J))
  if (!is.finite(tr)) return(NULL)
  list(p_eff = tr / 2, TIC = objective + tr)
}

# The corrected test for ONE nested pair. anova.admFit() orders the fits and
# formats the table; everything statistical is here.
.admLRT <- function(full, reduced) {
  sw_f <- .admSandwichParts(full)
  if (is.null(sw_f))
    stop("anova(): the larger model must be fitted with covMethod = \"r,s\" -- ",
         "the corrected reference distribution is built from its H and J, and ",
         "only that covMethod stores them.", call. = FALSE)
  # NESTING is decided on the FITS' parameter sets, not on the stored
  # covariance's. .admReduceNpdOmega() drops the omega block when the full
  # Hessian is not positive definite, so sw$par_names can be a strict subset of
  # the fit's -- and comparing that against the other fit's FULL set would report
  # a spurious "not nested" on a pair that plainly is.
  nm_f <- .admFitParNames(full)
  nm_r <- .admFitParNames(reduced)
  if (is.null(nm_f) || is.null(nm_r))
    stop("anova(): could not read the fits' parameter names.", call. = FALSE)
  # Nesting is a precondition, not something to approximate around: a non-nested
  # pair is a different problem (Vuong) and must not come back with a p-value.
  if (!all(nm_r %in% nm_f))
    stop("anova(): the models are not nested -- the smaller one has parameters ",
         "the larger does not (", paste(setdiff(nm_r, nm_f), collapse = ", "),
         "). A non-nested comparison needs a Vuong test, not this.",
         call. = FALSE)
  gamma <- setdiff(nm_f, nm_r)
  if (!length(gamma))
    stop("anova(): the two fits have the same parameters; nothing to test.",
         call. = FALSE)
  # A BOUNDARY RESTRICTION IS REFUSED, NOT RESCALED. Fixing a variance to zero
  # puts the null on the edge of the parameter space, where the limit is a
  # chi-bar-squared MIXTURE -- the expansion behind the eigenvalue weights
  # assumes an interior null and does not apply. Rescaling dOFV there produces a
  # finite, plausible p-value from the wrong reference distribution, which is
  # worse than refusing.
  #
  # Only the DIAGONAL is a boundary: Omega_ii >= 0, and .admBuildOptVec stores
  # it as log(Omega_ii) under `logchol_<eta>`. An OFF-diagonal covariance
  # (`chol_<eta_i>_<eta_j>`) may take either sign, so fixing one to zero is an
  # interior restriction and is tested normally.
  .bnd <- grep("^logchol_", gamma, value = TRUE)
  if (length(.bnd))
    stop("anova(): ", paste(.bnd, collapse = ", "), " is a variance, so ",
         "restricting it to zero puts the null on the BOUNDARY of the ",
         "parameter space. The limit is then a chi-bar-squared mixture rather ",
         "than the weighted chi-squares this rescales to, and no rescaling of ",
         "dOFV repairs that. Compare these models on TIC, or simulate under ",
         "the null.", call. = FALSE)
  q <- length(gamma)

  dOFV <- as.numeric(reduced$objective) - as.numeric(full$objective)
  # Nesting guarantees dOFV >= 0, so this costs nothing and catches a failed
  # full-model fit that the optimiser convergence flag does not -- three in
  # 150,000, at -8e33. Never pmax(dOFV, 0) as a matter of course: that turns a
  # failed fit into a silent non-rejection.
  if (!is.finite(dOFV))
    stop("anova(): dOFV is not finite; one of the fits failed.", call. = FALSE)
  if (dOFV < -1e-6)
    stop(sprintf(paste0("anova(): dOFV = %.4g < 0, which nesting forbids. The ",
                        "larger model failed to reach the smaller one's ",
                        "optimum -- refit it; this must not be clamped to 0."),
                 dOFV), call. = FALSE)
  dOFV <- max(dOFV, 0)

  # Both covariances from the LARGER fit, and both SUBSET AFTER INVERTING. That
  # order is the whole point: the inverse's gamma block is the Schur complement,
  # which PROFILES OUT the nuisance parameters, where inverting H[gamma, gamma]
  # would condition on them instead. Measured, the wrong order inflated
  # sum(lambda) 12x (4.66 -> 58.4) -- in the conservative direction, so the test
  # rejected at 0.0002 instead of 0.05 and no covariate would ever have been
  # selected. Pinned in test-inference.R.
  Hi <- tryCatch(solve(sw_f$H), error = function(e) NULL)
  if (is.null(Hi))
    stop("anova(): the larger fit's Hessian is not invertible.", call. = FALSE)
  # By NAME, never by position: .admCovThetaOrder() exists because a positional
  # match silently rotated SEs once, every number finite and plausible. The
  # lookup is into the COVARIANCE's names, which is where a reduced omega block
  # shows up -- testing an omega the sandwich had to drop must say so, not index
  # something else.
  k <- match(gamma, sw_f$par_names)
  if (anyNA(k))
    stop("anova(): ", paste(gamma[is.na(k)], collapse = ", "),
         " is not in the larger fit's stored covariance. The omega block is ",
         "dropped when the full Hessian is not positive definite, so a ",
         "restriction on omega cannot be tested from that fit.", call. = FALSE)
  A <- (2 * Hi)[k, k, drop = FALSE]
  B <- (Hi %*% sw_f$J %*% Hi)[k, k, drop = FALSE]

  lambda <- tryCatch(
    sort(Re(eigen(solve(A) %*% B, only.values = TRUE)$values), decreasing = TRUE),
    error = function(e) NULL)
  if (is.null(lambda) || !all(is.finite(lambda)) || any(lambda <= 0))
    stop("anova(): the eigenvalues of A^-1 B are not usable.", call. = FALSE)

  kap  <- tryCatch(kappa(A, exact = TRUE), error = function(e) Inf)
  warn <- NA_character_
  if (q == 1L) {
    method <- "scaled chi-squared"
    p <- stats::pchisq(dOFV / lambda, 1, lower.tail = FALSE)
  } else if (is.finite(kap) && kap < 1000) {
    method <- "weighted chi-squares"
    p <- 1 - .admRubenP(dOFV, lambda)
  } else {
    # Estimated eigenvalues repel when the block is ill-conditioned. The trace
    # is a smooth function of the same matrices and degrades gracefully -- at
    # kappa = 5387 it was the only form still near nominal. The limit here is
    # identifiability of the tested block, not its dimension.
    method <- "trace approximation"
    warn <- sprintf(paste0("the tested block is ill-conditioned (kappa = %.3g); ",
                           "using the trace approximation instead of the ",
                           "eigenvalue weights."), kap)
    warning("anova(): ", warn, call. = FALSE)
    p <- stats::pchisq(dOFV / (sum(lambda) / q), q, lower.tail = FALSE)
  }
  list(dOFV = dOFV, df = q, lambda = lambda, p = p, method = method,
       gamma = gamma, naive_p = stats::pchisq(dOFV, q, lower.tail = FALSE),
       warning = warn)
}

#' Compare nested aggregate-data fits
#'
#' `anova()` on two or more `admFit` objects gives the likelihood-ratio test with
#' the correction misspecification requires -- no separate function to learn, the
#' ordinary idiom carries the corrected reference distribution.
#'
#' The chi-squared reference distribution is a consequence of the information
#' equality `H = J`, which holds only under correct specification. The aggregate
#' objective is a quasi-likelihood, so `dOFV` is not chi-squared and testing it
#' against one rejects far too often -- measured at **20.2% for a nominal 5%**.
#'
#' The limit is instead a weighted sum of chi-squares, `sum_i lambda_i chi2_1`,
#' with `lambda` the eigenvalues of `A^-1 B`, `A` the naive covariance of the
#' tested block and `B` its sandwich covariance. Under correct specification
#' every `lambda_i` is 1 and this reduces to the ordinary test exactly, so the
#' `Weight` column doubles as a diagnostic: values far from 1 say the uncorrected
#' test would have been badly wrong.
#'
#' The correction is a monotone rescaling of `dOFV`, so at matched actual size it
#' is the SAME test and gives up no power. What it supplies is the threshold: the
#' critical value that actually delivers 5% ranged from 4.1 to 14.0 across
#' configurations, and it cannot be tabulated, because it depends on the
#' between-subject spread and the residual model of the problem at hand. Obtained
#' any other way it needs a simulation under the null for every comparison.
#'
#' The larger model must be fitted with `covMethod = "r,s"`, which is what stores
#' the `H` and `J` the weights are built from.
#'
#' Refused rather than approximated: non-nested comparisons (that is a Vuong
#' test) and boundary restrictions such as `omega^2 = 0`, whose limit is a
#' chi-bar-squared mixture the expansion behind this does not cover.
#'
#' One thing no rescaling repairs: a misspecified model can converge to a
#' non-zero pseudo-true value even where a covariate has no effect (measured
#' `+0.0117`, p < 1e-4, with a proportional residual). That cannot be detected
#' from a single fit.
#'
#' @param object An `admFit`.
#' @param ... One or more further `admFit` objects to compare it with.
#' @return A `data.frame` of class `anova.admFit`, one row per model ordered
#'   smallest to largest, carrying `dOFV`, `Weight`, the corrected `p` and the
#'   uncorrected `p.naive` beside it for contrast.
#' @examples
#' \dontrun{
#' f1 <- nlmixr2(mod_cov, admData(), est = "adgh",
#'               control = adghControl(studies = st, covMethod = "r,s"))
#' f0 <- nlmixr2(mod_nocov, admData(), est = "adgh",
#'               control = adghControl(studies = st, covMethod = "r,s"))
#' anova(f1, f0)
#' }
#' @method anova admFit
#' @export
anova.admFit <- function(object, ...) {
  fits <- c(list(object), list(...))
  ok   <- vapply(fits, function(z) inherits(z, "admFit"), logical(1))
  if (!all(ok))
    stop("anova(): every argument must be an admFit; got ",
         paste(vapply(fits[!ok], function(z) class(z)[1L], character(1)),
               collapse = ", "), ".", call. = FALSE)
  if (length(fits) < 2L)
    stop("anova() on a single admFit has nothing to compare it with. Pass the ",
         "nested pair, e.g. anova(full, reduced).", call. = FALSE)
  npar <- vapply(fits, function(z) length(.admFitParNames(z)), integer(1))
  # Smallest first, the way anova.lme and anova.lm order their tables, so each
  # test row reads "the model above, plus these parameters".
  ord  <- order(npar)
  fits <- fits[ord]; npar <- npar[ord]
  nms  <- make.unique(vapply(fits, function(z) z$env$method %||% "adm",
                             character(1)), sep = "_")

  n   <- length(fits)
  num <- function(nm) vapply(fits, function(z)
    as.numeric(z$env[[nm]] %||% NA_real_), numeric(1))
  out <- data.frame(
    Npar    = npar,
    OBJF    = vapply(fits, function(z) as.numeric(z$objective), numeric(1)),
    AIC     = num("AIC"),
    BIC     = num("BIC"),
    TIC     = num("TIC"),
    Test    = c(NA_character_, sprintf("%d vs %d", seq_len(n - 1L), 2:n)),
    dOFV    = NA_real_,
    Df      = NA_integer_,
    Weight  = NA_character_,
    p       = NA_real_,
    p.naive = NA_real_,
    row.names        = nms,
    stringsAsFactors = FALSE,
    check.names      = FALSE)
  lams <- vector("list", n)
  for (i in seq_len(n - 1L)) {
    r <- .admLRT(fits[[i + 1L]], fits[[i]])
    out$dOFV[i + 1L]    <- r$dOFV
    out$Df[i + 1L]      <- r$df
    out$Weight[i + 1L]  <- paste(sprintf("%.3f", r$lambda), collapse = ", ")
    out$p[i + 1L]       <- r$p
    out$p.naive[i + 1L] <- r$naive_p
    lams[[i + 1L]]      <- r$lambda
  }
  # TIC is dropped rather than printed as a column of NA when no fit carries it.
  if (all(is.na(out$TIC))) out$TIC <- NULL
  attr(out, "lambda")  <- lams
  attr(out, "heading") <- c(
    "Likelihood-ratio test, corrected for misspecification",
    paste("dOFV rescaled by the eigenvalues of A^-1 B;",
          "Weight = 1 is the ordinary chi-squared test"))
  class(out) <- c("anova.admFit", "anova", "data.frame")
  out
}
