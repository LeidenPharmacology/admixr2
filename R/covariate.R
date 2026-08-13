# =============================================================================
# Covariate marginalisation for aggregate data modelling
# =============================================================================
#
# A study's subjects span a distribution of covariate values (`cov_dist`), and
# the aggregate (E, V) it reports is MARGINAL over that distribution. So the
# prediction has to be marginal too:
#
#   mu = E_{a,eta}[f]        V = Cov_{a,eta}(f) + residual
#
# Covariate-induced between-subject variability belongs INSIDE V_pred, because
# that is where the observed V carries it. Three paths compute those moments --
# collapse, u-quantile and the general per-row path -- chosen per study by
# .admCheckCovariates(); see each section below.
#
# NODE METHODS (gl / gh / taylor) WERE REMOVED, deliberately. They score the
# study at fixed covariate values and combine the per-node -2LL linearly. On the
# one pooled (E, V) a publication reports that is sum_k w_k * NLL(same obs),
# i.e. E_a[-2LL] -- an average of log-densities, so -2 log of an unnormalised
# geometric mean of densities, which is not a likelihood and is biased (~58% low
# on the covariate coefficient; validation/covariate-threeway.R). On per-node
# data it IS a likelihood, but only because sum_k c_k n NLL_k is exactly
# sum_k n_k NLL_k with n_k = c_k n -- ordinary multi-study fitting with the
# quadrature weights standing in for stratum sizes. When a publication reports
# strata it also reports their real n_k, which is strictly better, so the
# machinery was redundant where valid and invalid where it was not. Fit stratum
# summaries as ordinary studies, each with its own `cov` and `n`.
#
# (`%||%` is defined in utils.R.)

# Append this study's covariate columns to an rxSolve params frame (matrix or
# data.frame). Returns `mat` unchanged when there is nothing to add.
#
# ONLY names present in `cov_s` are added -- deliberately NOT
# setdiff(mod$params, colnames(mat)) with a 0 default. That blanket zero-fill is
# the pattern this package has been bitten by twice: it clobbered genuine
# hard-coded model constants (qout/vb -> /0 -> an NA objective), and it handed
# the solve lambda = 0 for an ESTIMATED boxCox/yeoJohnson, making the sens
# gradient ~60x wrong for boxCox and NaN for yeoJohnson. Everything that is not
# a covariate is rxSolve's to default from the model itself.
#
# Named .adm* on purpose: the dev-mode daemon payload is collected by a
# /^\.(adm|adfo|adgh|adirmc)/i regex, so a helper named .cov_fills would be
# missing inside every parallel-restart worker.
.admCovCols <- function(mat, mod_params, cov_s, cov_rows = NULL) {
  # PER-ROW values (general path): each simulated subject carries its own
  # covariate, so rxode2 evaluates whatever functional form the model contains.
  if (!is.null(cov_rows)) {
    nms <- setdiff(intersect(colnames(cov_rows), mod_params), colnames(mat))
    if (length(nms) == 0L) return(mat)
    if (nrow(cov_rows) != nrow(mat))
      stop(".admCovCols: cov_rows has ", nrow(cov_rows), " rows but the params ",
           "frame has ", nrow(mat), ". Recycling here would hand subjects the ",
           "wrong covariate values silently.", call. = FALSE)
    add <- as.matrix(cov_rows[, nms, drop = FALSE])
    return(if (is.data.frame(mat)) cbind(mat, as.data.frame(add, check.names = FALSE))
           else                    cbind(mat, add))
  }
  if (is.null(cov_s) || length(cov_s) == 0L) return(mat)
  nms <- setdiff(intersect(names(cov_s), mod_params), colnames(mat))
  if (length(nms) == 0L) return(mat)
  add <- matrix(rep(as.numeric(unlist(cov_s[nms], use.names = FALSE)),
                    each = nrow(mat)),
                nrow(mat), length(nms), dimnames = list(NULL, nms))
  if (is.data.frame(mat)) cbind(mat, as.data.frame(add, check.names = FALSE))
  else                    cbind(mat, add)
}

# =============================================================================
# Covariate collapse (exact marginal moments, no nodes)
# =============================================================================

# Map covariate -> (coefficient theta, eta) from rxode2's OWN mu-reference
# metadata, so nothing has to be declared by hand. muRefCovariateDataFrame gives
# (theta, covariate, covariateParameter); muRefDataFrame gives (theta, eta).
# Joining on `theta` says which eta's mu-referenced argument each covariate
# enters, and with which coefficient.
.admCovMap <- function(ui) {
  cd <- tryCatch(ui$muRefCovariateDataFrame, error = function(e) NULL)
  md <- tryCatch(ui$muRefDataFrame,          error = function(e) NULL)
  if (is.null(cd) || is.null(md)) return(NULL)
  if (!NROW(cd) || !NROW(md)) return(NULL)
  m <- merge(cd, md, by = "theta")
  if (!NROW(m)) return(NULL)
  data.frame(covariate = as.character(m$covariate),
             coef      = as.character(m$covariateParameter),
             eta       = as.character(m$eta),
             stringsAsFactors = FALSE)
}

# Sigma_a for the named covariates, in the order given.
.admCovDistSigma <- function(cov_dist, nms) {
  d   <- length(nms)
  sds <- vapply(nms, function(n) as.numeric(cov_dist[[n]]$sd), numeric(1))
  if (!is.null(cov_dist$Sigma)) return(cov_dist$Sigma[nms, nms, drop = FALSE])
  rho <- cov_dist$rho %||% 0
  S <- diag(sds^2, d)
  if (d >= 2L && abs(rho) > 1e-12)
    for (i in seq_len(d - 1L)) for (j in seq(i + 1L, d))
      S[i, j] <- S[j, i] <- rho * sds[i] * sds[j]
  S
}

# Effective Cholesky for a study whose subjects SPAN a covariate distribution.
#
# With mu-referencing the covariate and the random effect enter the SAME
# argument:  param <- g(theta + theta_cov*a + eta).  So for a ~ N(mu_a, Sigma_a),
#
#     theta + theta_cov*a + eta  ==  theta + theta_cov*mu_a + s,
#     s ~ N(0, Omega + J Sigma_a J'),      J[eta_k, k] = theta_cov_k
#
# EXACTLY. Solving at the covariate MEAN and inflating Omega therefore reproduces
# the marginal moments with one solve set -- no quadrature nodes, no importance
# weights, and the same cost as a fit with no covariate at all.
#
# Exact only for a LINEAR covariate effect on the mu-referenced scale and a
# normal covariate. Returns NULL whenever that cannot be established, so every
# caller falls back to the plain Cholesky rather than silently approximating.
.admCovInflateL <- function(pars, pinfo, s) {
  cd <- s$cov_dist
  if (is.null(cd) || pinfo$n_eta == 0L || is.null(pars$L)) return(NULL)
  cmap <- pinfo$cov_map
  if (is.null(cmap) || !NROW(cmap)) return(NULL)
  eta_nms <- pinfo$eta_col_names
  has_sd  <- vapply(names(cd), function(n)
    is.list(cd[[n]]) && !is.null(cd[[n]]$sd), logical(1))
  keep <- cmap$covariate %in% names(cd)[has_sd] & cmap$eta %in% eta_nms
  use  <- cmap[keep, , drop = FALSE]
  if (!NROW(use)) return(NULL)
  J <- matrix(0, pinfo$n_eta, NROW(use))
  for (k in seq_len(NROW(use))) {
    cf <- pars$struct[[use$coef[k]]]
    if (is.null(cf) || !is.finite(cf)) return(NULL)
    J[match(use$eta[k], eta_nms), k] <- cf
  }
  Sig <- .admCovDistSigma(cd, use$covariate)
  Om  <- tcrossprod(pars$L) + J %*% Sig %*% t(J)
  tryCatch(t(chol(Om)), error = function(e) NULL)
}

# The Cholesky a given study should be simulated with: inflated when the study
# declares a covariate DISTRIBUTION, the plain one otherwise.
#
# Returns NULL -- never a silent fall back to pars$L -- when a study DOES declare
# cov_dist and the inflation cannot be built. Falling back there would solve at
# the covariate mean and understate the variance, which is a plausible-looking
# fit of the wrong model. Callers treat NULL as "cannot evaluate" (Inf / NA).
# .admCheckCovariates() rejects the reachable causes up front, so in practice
# this only fires on a pathological Omega.
.admStudyL <- function(pars, pinfo, s) {
  if (!identical(s$.adm_cov_path, "collapse")) return(pars$L)
  .admCovInflateL(pars, pinfo, s)
}

# Covariate columns for a params frame that stacks `n_blk` blocks of `n_sim`
# rows, each block a PERTURBATION OF THE SAME SUBJECTS (the finite-difference
# frames in .admGrad). The covariate rows are tiled per block so every block sees
# the same subjects' covariates -- which is exactly what makes the difference a
# common-random-numbers one. Tiling with the wrong stride would give each
# perturbation different subjects and quietly turn the gradient into noise.
.admCovColsTiled <- function(mat, mod_params, s, n_sim, n_blk) {
  cr <- s[["cov_rows"]]
  if (!is.null(cr)) cr <- cr[rep(seq_len(n_sim), times = n_blk), , drop = FALSE]
  .admCovCols(mat, mod_params, s[["cov"]], cr)
}

# Attach per-row covariate values to a study for the GENERAL path. Returns the
# study unchanged on the collapse path (where the covariate is held at its mean
# and its variance lives in Omega) and when no distribution is declared.
.admStudyCovRows <- function(s, pinfo, n_row) {
  if (!identical(s$.adm_cov_path, "rows")) return(s)
  s$cov_rows <- .admCovRowsFor(s$cov_dist, n_row, pinfo$n_eta)
  s
}

# Up-front validation for every study that declares a covariate distribution.
#
# Scope, deliberately narrow and enforced rather than assumed: the collapse is
# EXACT only for a covariate that enters as a bare `theta * COV` product inside a
# mu-referenced expression, appears nowhere else in the model, and is normally
# distributed. rxode2 records the first of those in muRefCovariateDataFrame --
# and records it ONLY for that bare product. Measured on 6 covariate models:
#
#   cl <- exp(tcl + tcov*WT + eta.cl)              -> 1 row   (collapsible)
#   cl <- exp(tcl + tcov*log(WT/70) + eta.cl)      -> 0 rows
#   cl <- exp(tcl + eta.cl)*(WT/70)^tcov           -> 0 rows
#   cl <- exp(tcl + eta.cl)*exp(tcov*WT)           -> 0 rows
#   cl <- exp(tcl + eta.cl)*(1 + emax*WT/(ec50+WT))-> 0 rows
#   cl <- exp(tcl + tcov*SEX + eta.cl)             -> 1 row   (but SEX is not normal)
#
# The allometric form is mathematically collapsible (it is linear in log(WT)) and
# is still refused here, because rxode2's metadata cannot distinguish it from the
# Emax form. Refusing is the safe direction: the alternative is a fit that runs,
# converges, and reports an omega that has quietly eaten the covariate spread.
# Refuse `cov_dist` for an estimator that has no covariate path.
#
# Silence here is the dangerous outcome, not an error: every study carries a
# covariate VALUE as well (derived from the distribution when not given), so an
# unwired estimator does not fail -- it solves at the covariate mean and reports
# a fit whose omega has quietly absorbed the between-subject covariate spread.
# Measured on the general path's own test model, that is omega 0.30 -> 0.44.
.admRefuseCovariates <- function(studies, est) {
  has <- vapply(studies, function(s) !is.null(s[["cov_dist"]]), logical(1))
  if (!any(has)) return(invisible(NULL))
  stop("admixr2: `", est, "` does not support covariate marginalisation. ",
       "Stud", if (sum(has) > 1L) "ies " else "y ",
       paste(sQuote(names(studies)[has]), collapse = ", "),
       " declare(s) `cov_dist`. Use `admc`, which marginalises over the ",
       "covariate distribution; running here would silently solve at the ",
       "covariate mean and inflate omega instead.", call. = FALSE)
}

.admCheckCovariates <- function(.ui, pinfo, studies, grad, est = NULL) {
  has <- vapply(studies, function(s) !is.null(s$cov_dist), logical(1))
  if (!any(has)) return(studies)
  bad <- function(...) stop("admixr2: ", ..., call. = FALSE)

  cmap <- pinfo$cov_map
  covs <- tryCatch(.ui$allCovs, error = function(e) character(0))
  for (nm in names(studies)[has]) {
    cd <- studies[[nm]]$cov_dist
    ok_collapse <- pinfo$n_eta > 0L
    ok_uq       <- pinfo$n_eta > 0L
    uq          <- list()

    for (cv in names(cd)) {
      if (!cv %in% covs)
        bad("study '", nm, "' declares `cov_dist` for '", cv,
            "', which the model never reads. Model covariates: ",
            if (length(covs)) paste(covs, collapse = ", ") else "(none)", ".")
      sp <- cd[[cv]]
      if (!is.list(sp)) bad("`cov_dist` for '", cv, "' must be a list.")
      normal <- !is.null(sp$mu) && !is.null(sp$sd) &&
                is.finite(sp$mu) && is.finite(sp$sd) && sp$sd > 0
      lnorm  <- !is.null(sp$meanlog) && !is.null(sp$sdlog) &&
                is.finite(sp$meanlog) && is.finite(sp$sdlog) && sp$sdlog > 0
      disc   <- !is.null(sp$values) && length(sp$values) > 0L &&
                all(is.finite(sp$values))
      userq  <- is.function(sp$quantile)
      if (!(normal || lnorm || disc || userq))
        bad("`cov_dist` for '", cv, "' in study '", nm, "' is not a supported ",
            "distribution. Give one of: `mu`+`sd` (normal), `meanlog`+`sdlog` ",
            "(lognormal), `values`(+`probs`) (discrete/categorical), or ",
            "`quantile` (a function of a uniform).")

      # Fill the solve value from the distribution when the study omits it.
      if (is.null(studies[[nm]][["cov"]][[cv]])) {
        m <- .admCovMeanOf(sp)
        if (!is.null(m) && is.finite(m)) studies[[nm]][["cov"]][[cv]] <- m
      }

      once <- .admNameOccurrence(.ui, cv)[[cv]] == 1L
      pe   <- if (once) .admCovParamEta(.ui, cv, pinfo$eta_col_names) else NULL

      # COLLAPSE: exact closed form. Needs a normal covariate entering as a
      # single bare `theta * COV` term.
      if (!normal || is.null(cmap) || !cv %in% cmap$covariate || !once)
        ok_collapse <- FALSE
      # u-QUANTILE: any functional form, but the covariate's WHOLE effect has to
      # fit in one eta column -- so it must appear exactly once, in a parameter
      # assignment that also carries an eta, with a reference value to measure
      # Delta against.
      if (!once || is.null(pe) || is.null(studies[[nm]][["cov"]][[cv]])) ok_uq <- FALSE
      else uq[[length(uq) + 1L]] <- list(cov = cv, param = pe$param, eta = pe$eta,
                                         idx = match(pe$eta, pinfo$eta_col_names))
    }

    # A GRADIENT is only carried through the general path. On "rows" the
    # covariate is data -- a per-row params column -- so the existing sensitivity
    # directions differentiate exactly the function the NLL evaluates, with no
    # new chain rule. "collapse" moves Omega itself (Omega + J Sigma_a J', with J
    # carrying the covariate coefficient) and "uq" replaces an eta column with
    # u = F_u^-1(Phi(z)); neither derivative exists yet.
    #
    # So the gradient mode is part of what makes a path VALID, not a reason to
    # refuse the fit. Erroring here would make a covariate model fail out of the
    # box, since every estimator defaults to a gradient.
    if (!identical(grad, "none")) { ok_collapse <- FALSE; ok_uq <- FALSE }


    # Most efficient VALID path wins. "rows" assumes nothing at all: every
    # simulated subject carries its own covariate value, so rxode2 evaluates the
    # whole model -- covariate on several parameters, on a parameter with no eta,
    # or interacting with one. It is the only path with no structural
    # precondition, hence the fallback.
    studies[[nm]]$.adm_cov_path <- if (ok_collapse) "collapse"
                                   else if (ok_uq)  "uq" else "rows"
    studies[[nm]]$.adm_cov_collapse <- ok_collapse
    studies[[nm]]$.adm_cov_uq <- if (ok_uq) uq else NULL
  }

  studies
}

# =============================================================================
# General covariate marginalisation (ANY functional form)
# =============================================================================
#
# The collapse above is exact but narrow: it needs the covariate to enter as a
# bare `theta * COV` term so its variance can be folded into Omega. Everything
# else -- (WT/70)^theta, exp(theta*WT), allometric-in-log, Emax, if/else,
# categorical -- goes through here instead.
#
# The trick is that we never have to KNOW the functional form. A covariate
# reaches rxode2 only as a column of the params frame, so if each simulated
# subject carries its OWN covariate value, rxode2 evaluates whatever the model
# contains and the aggregate moments come out of the pooled ensemble. That is
# also why it costs nothing extra for admc: n_sim subjects still means n_sim
# rows, each now carrying its own (covariate, eta) pair rather than a shared
# covariate.

# Quantile function for one covariate spec. Supported:
#   list(mu, sd)              normal
#   list(meanlog, sdlog)      lognormal
#   list(values, probs)       discrete / categorical (probs default to uniform)
#   list(quantile = f(u))     anything else, supplied by the caller
.admCovQuantile <- function(spec, u) {
  if (is.function(spec$quantile)) return(as.numeric(spec$quantile(u)))
  if (!is.null(spec$values)) {
    pr <- spec$probs %||% rep(1 / length(spec$values), length(spec$values))
    pr <- pr / sum(pr)
    return(as.numeric(spec$values)[findInterval(u, cumsum(pr), rightmost.closed = TRUE) + 1L])
  }
  if (!is.null(spec$meanlog)) return(stats::qlnorm(u, spec$meanlog, spec$sdlog))
  stats::qnorm(u, spec$mu, spec$sd)
}

# Deterministic per-row covariate values for `n` simulated subjects.
#
# Deterministic on purpose: the covariate distribution is DATA, not a parameter,
# so the same rows must come back on every objective evaluation or the optimizer
# sees noise. Common random numbers therefore hold with no seed plumbing.
#
# The uniforms are taken from Sobol dimensions AFTER the eta dimensions.
# sobol(n, dim = k)[, 1:j] is exactly sobol(n, dim = j) (verified), so this
# yields dimensions genuinely different from the ones .admMakeZ used for eta.
# Drawing a separate halton/sobol sequence instead would NOT: every low-
# discrepancy family starts from the same base-2 van der Corput sequence, so
# covariate column 1 would have been a copy of eta column 1.
.admCovRowsFor <- function(cov_dist, n, n_eta) {
  nms <- names(cov_dist)
  d   <- length(nms)
  u   <- randtoolbox::sobol(n, dim = n_eta + d)
  if (!is.matrix(u)) u <- matrix(u, nrow = n)
  u   <- u[, n_eta + seq_len(d), drop = FALSE]
  # sobol emits an exact 0 in its first row; qnorm(0) is -Inf.
  u   <- pmin(pmax(u, .Machine$double.eps), 1 - .Machine$double.eps)
  out <- vapply(seq_len(d), function(k) .admCovQuantile(cov_dist[[nms[k]]], u[, k]),
                numeric(n))
  if (!is.matrix(out)) out <- matrix(out, nrow = n)
  colnames(out) <- nms
  out
}

# Deterministic quadrature nodes + weights for one covariate (adgh path).
# Discrete specs enumerate exactly; normal/lognormal use Gauss-Hermite.
.admCovNodesFor <- function(spec, n_nodes) {
  if (!is.null(spec$values)) {
    pr <- spec$probs %||% rep(1 / length(spec$values), length(spec$values))
    return(list(x = as.numeric(spec$values), w = pr / sum(pr)))
  }
  g <- .adghNodes1(n_nodes)                       # standard-normal nodes/weights
  if (!is.null(spec$meanlog)) list(x = exp(spec$meanlog + spec$sdlog * g$x), w = g$w)
  else                        list(x = spec$mu + spec$sd * g$x,              w = g$w)
}

# Product grid over several covariates: every combination, weights multiplied.
.admCovGrid <- function(cov_dist, n_nodes) {
  nms <- names(cov_dist)
  one <- lapply(nms, function(n) .admCovNodesFor(cov_dist[[n]], n_nodes))
  X   <- as.matrix(expand.grid(lapply(one, `[[`, "x"), KEEP.OUT.ATTRS = FALSE))
  W   <- Reduce(`*`, lapply(seq_along(one), function(k)
           rep(rep(one[[k]]$w, each = prod(vapply(one[seq_len(k - 1L)],
                 function(o) length(o$x), integer(1)))),
               length.out = nrow(X))))
  colnames(X) <- nms
  list(X = X, W = W / sum(W))
}

# =============================================================================
# u-quantile marginalisation -- ANY functional form
# =============================================================================
#
# The model sees the covariate and the random effect only through their sum in
# the mu-referenced argument, u = Delta(a) + eta. Delta(a) is arithmetic, so u's
# distribution can be pinned down exactly BEFORE any solve, and the whole solve
# budget spent representing it:
#
#     F_u(t) = E_a[ Phi( (t - Delta(a)) / sd_eta ) ]
#
# Nothing here inspects the functional form. Delta is measured from the model
# itself (see .admCovDelta), so power, exponential, allometric, Emax, if/else
# and categorical covariate effects are all handled identically.

# Which parameter assignment does the covariate enter, and which eta shares it?
# The shift structure the whole method rests on is exactly "covariate and eta in
# the same mu-referenced argument", so this also decides whether it applies.
.admCovParamEta <- function(ui, cov, eta_names) {
  lst <- tryCatch(ui$lstExpr, error = function(e) NULL)
  if (is.null(lst)) return(NULL)
  for (e in lst) {
    if (!is.call(e) || length(e) < 3L) next
    if (!as.character(e[[1L]])[1L] %in% c("<-", "=")) next
    v <- all.vars(e[[3L]])
    if (!cov %in% v) next
    et <- intersect(eta_names, v)
    if (length(et) != 1L) next
    return(list(param = as.character(e[[2L]])[1L], eta = et))
  }
  NULL
}

# Delta(a) MEASURED from the model, not parsed out of it.
#
# One rxSolve at eta = 0 over a covariate grid; the mu-referenced parameter comes
# back as an lhs column, and Delta(a) = log(param(a) / param(a_ref)). Exact to
# machine precision (measured 8e-17 on an allometric model) and indifferent to
# how the covariate effect is written.
#
# Costs one extra rxSolve call (~11 ms) per objective evaluation. Delta depends
# on the covariate coefficients, which move, so it cannot be cached across
# iterations -- but 11 ms against a solve budget in seconds is noise.
# ngrid = 32: Gauss-Hermite over a smooth Delta converges long before this.
# Measured against a 200-node reference, F_u agreed to 3.8e-15 at 16 nodes and
# 2.2e-15 at 24, over covariate spreads from sdlog 0.15 to 0.45. The node count
# costs twice -- once in the probe's rxSolve rows, once in every Newton
# iteration's n_sim x n_node matrices -- so 64 was paying for nothing.
.admCovDelta <- function(rxMod, pars, pinfo, s, cov, param_nm, cores,
                         ngrid = 32L) {
  spec  <- s$cov_dist[[cov]]
  a_ref <- s[["cov"]][[cov]]
  if (is.null(a_ref)) return(NULL)
  an <- .admCovANodes(spec, ngrid)
  a  <- an$x
  aa <- c(a, a_ref)
  pm <- .admMakeParamsList(length(aa), pinfo, 1L)[[1L]]
  for (nm in pinfo$struct_names) pm[, nm] <- pars$struct[[nm]]
  pm <- .admCovCols(pm, rxMod$params, NULL,
                    matrix(aa, ncol = 1L, dimnames = list(NULL, cov)))
  out <- tryCatch(rxode2::rxSolve(rxMod, params = as.data.frame(pm),
                                  events = s$ev_full, cores = cores,
                                  nDisplayProgress = .Machine$integer.max),
                  error = function(e) NULL)
  if (is.null(out) || is.null(out[[param_nm]])) return(NULL)
  v <- out[[param_nm]][out[["time"]] == s$times[[1L]]]
  if (length(v) != length(aa) || any(!is.finite(v)) || any(v <= 0)) return(NULL)
  list(a = a, delta = log(v[seq_along(a)]) - log(v[length(aa)]), w = an$w)
}

# The covariate value the model is SOLVED at when the study does not name one.
#
# Every path needs it: the collapse solves at the covariate mean so the model
# itself produces the mean shift, and u-quantile measures Delta relative to it.
# Deriving it from `cov_dist` rather than making the user restate it removes a
# way for the two to disagree silently.
.admCovMeanOf <- function(spec) {
  if (!is.null(spec$values)) {
    pr <- spec$probs %||% rep(1, length(spec$values))
    return(sum(as.numeric(spec$values) * pr) / sum(pr))
  }
  if (is.function(spec$quantile))
    return(mean(spec$quantile((seq_len(1024L) - 0.5) / 1024L)))
  if (!is.null(spec$meanlog)) return(exp(spec$meanlog + spec$sdlog^2 / 2))
  spec$mu
}

# Spread of a covariate spec, on the scale the identifiability ridge lives on.
.admCovSdOf <- function(spec) {
  if (!is.null(spec$values)) {
    pr <- spec$probs %||% rep(1, length(spec$values)); pr <- pr / sum(pr)
    m <- sum(as.numeric(spec$values) * pr)
    return(sqrt(sum(pr * (as.numeric(spec$values) - m)^2)))
  }
  if (is.function(spec$quantile))
    return(stats::sd(spec$quantile((seq_len(1024L) - 0.5) / 1024L)))
  if (!is.null(spec$meanlog)) return(spec$sdlog)
  spec$sd
}

# Warn when a covariate coefficient cannot be identified from the data supplied.
#
# When a covariate shares a mu-referenced argument with a random effect, the
# model sees only u = Delta(a) + eta, so ONE population determines just two
# quantities -- u's mean and variance -- against three parameters. The likelihood
# is then EXACTLY flat along
#
#     theta'   = theta + (b - b')*mu_a
#     omega'^2 = omega^2 + (b^2 - b'^2)*sd_a^2
#
# (verified: the objective is bit-identical across b from 0.40 to 1.10). Only
# BETWEEN-STUDY variation in the covariate distribution breaks it -- differing
# means break the first equation, differing spreads the second. Note this is
# variation in the DATA, not a between-study random effect: admixr2 has none, and
# a tau^2 would in fact compete with the covariate for the same signal.
#
# A covariate on a parameter with NO random effect is not affected: there is no
# omega for its variance to be absorbed into, so it is identified by shape.
.admWarnCovIdentifiability <- function(.ui, pinfo, studies) {
  has <- vapply(studies, function(s) !is.null(s[["cov_dist"]]), logical(1))
  if (!any(has)) return(invisible(NULL))
  cds <- lapply(studies[has], `[[`, "cov_dist")
  for (cv in unique(unlist(lapply(cds, names)))) {
    # only covariates that share an argument with an eta sit on the ridge
    if (is.null(.admCovParamEta(.ui, cv, pinfo$eta_col_names))) next
    sp <- Filter(Negate(is.null), lapply(cds, `[[`, cv))
    if (length(sp) == 0L) next
    mu <- vapply(sp, function(x) .admCovMeanOf(x) %||% NA_real_, numeric(1))
    sd <- vapply(sp, function(x) .admCovSdOf(x)  %||% NA_real_, numeric(1))
    varies <- (length(unique(signif(mu, 10))) > 1L) ||
              (length(unique(signif(sd, 10))) > 1L)
    if (!varies)
      warning("admixr2: the coefficient on covariate '", cv, "' is not ",
              "identifiable from these data. It enters the same argument as a ",
              "random effect, and every study declaring it has the SAME ",
              "covariate distribution, so the likelihood is exactly flat along ",
              "a trade-off between that coefficient, the corresponding fixed ",
              "effect and omega. Identification needs studies whose covariate ",
              "MEANS or SPREADS differ.", call. = FALSE)
  }
  invisible(NULL)
}

# Nodes + weights for the a-integral inside F_u.
#
# GAUSS-HERMITE, not equal-weight quantiles. Here we are INTEGRATING over the
# covariate, and GH carries the tails through its weights; N equal-weight
# quantiles truncate the covariate at the (0.5/N) and (1-0.5/N) points and drop
# the tail mass entirely. That understates the spread of Delta, and the fit then
# inflates the covariate coefficient to make up the missing variance -- measured
# as tcov 0.80 against a truth of 0.75 on an allometric model, while omega stayed
# correct. The mirror-image mistake (using GH to REPRESENT u's distribution) is
# equally wrong; see .admCovUQuantile, which uses the closed-form CDF.
.admCovANodes <- function(spec, ngrid) {
  if (!is.null(spec$values)) {
    pr <- spec$probs %||% rep(1, length(spec$values))
    return(list(x = as.numeric(spec$values), w = pr / sum(pr)))
  }
  if (is.function(spec$quantile)) {           # arbitrary user CDF: no GH rule
    p <- (seq_len(ngrid) - 0.5) / ngrid
    return(list(x = spec$quantile(p), w = rep(1 / ngrid, ngrid)))
  }
  g <- .adghNodes1(ngrid)                     # standard-normal nodes/weights
  if (!is.null(spec$meanlog))
    list(x = exp(spec$meanlog + spec$sdlog * g$x), w = g$w)
  else
    list(x = spec$mu + spec$sd * g$x, w = g$w)
}

# Inverse CDF of u = Delta(a) + eta, evaluated at the uniforms `p`.
# Phi carries the tails analytically, which a grid of summed draws cannot.
# Inverse CDF of u = Delta(a) + eta at the uniforms `p`, by NEWTON on the exact
# mixture CDF rather than interpolation on a grid.
#
# F_u(t) = sum_j w_j Phi((t - Delta_j)/s) is smooth and strictly increasing, and
# its density is available in closed form, so Newton converges in a handful of
# steps from a moment-matched normal start.
#
# Why not the grid: interpolating F_u on a grid makes u PIECEWISE LINEAR in the
# parameters, because Delta and s move the grid's knots. The objective then has
# small kinks, which (a) is the exact mismatch that stops an analytic gradient
# being consistent with the objective it differentiates, and (b) defeats a
# trust-region derivative-free search, which builds quadratic models. Newton
# gives the true implicit function, so d(u)/d(theta) = -(dF_u/d(theta)) / f_u(u)
# is the derivative of what is actually computed.
#
# Bracketed and clamped: f_u underflows in the far tails, so an undamped Newton
# step there can throw the iterate out of the support entirely.
.admCovUQuantile <- function(dl, sd_eta, p, tol = 1e-11, maxit = 50L) {
  d <- dl$delta; w <- dl$w
  lo <- min(d) - 12 * sd_eta
  hi <- max(d) + 12 * sd_eta
  # moment-matched normal start: exact when Delta is degenerate, close otherwise
  m  <- sum(w * d)
  sdu <- sqrt(sum(w * (d - m)^2) + sd_eta^2)
  u  <- pmin(pmax(m + sdu * stats::qnorm(p), lo), hi)
  for (it in seq_len(maxit)) {
    z  <- outer(u, d, "-") / sd_eta
    r  <- as.numeric(stats::pnorm(z) %*% w) - p
    if (max(abs(r)) < tol) break
    fu <- as.numeric(stats::dnorm(z) %*% w) / sd_eta
    u  <- pmin(pmax(u - r / pmax(fu, 1e-300), lo), hi)
  }
  u
}
# Replace the covariate-affected eta column with draws of u = Delta(a) + eta.
#
# The substitution is by INVERSE TRANSFORM of the very dimension that would have
# produced that eta: u = F_u^-1(Phi(z_j)). So the low-discrepancy quality of the
# existing draws carries over, u gets exactly the right marginal, and its
# independence from the other eta columns is preserved by construction -- which
# a sorted set of quantiles pasted in would have destroyed.
#
# Requires Omega DIAGONAL in the affected row (or a single eta): with a
# correlated Omega, column j is a mixture of several z columns and replacing it
# would break the joint distribution. Enforced by .admCheckCovariates().
.admCovUQEta <- function(eta_mat, z, pars, pinfo, s, rxMod, cores) {
  info <- s$.adm_cov_uq
  if (is.null(info)) return(eta_mat)
  for (k in seq_along(info)) {
    it <- info[[k]]
    dl <- .admCovDelta(rxMod, pars, pinfo, s, it$cov, it$param, cores)
    if (is.null(dl)) return(NULL)                 # never silently fall back
    j  <- match(it$eta, pinfo$eta_col_names)
    sd_eta <- sqrt(tcrossprod(pars$L)[j, j])
    u <- .admCovUQuantile(dl, sd_eta, stats::pnorm(z[, j]))
    if (any(!is.finite(u))) return(NULL)
    eta_mat[, j] <- u
  }
  eta_mat
}

# Refuse the removed node-quadrature inputs, loudly.
#
# gl/gh/taylor covariate marginalisation was removed (see the file header). Old
# code passing `weight`, `cov_method` or a `quadrature` attribute would
# otherwise be accepted and silently fitted as an ordinary unweighted study
# list, which is a DIFFERENT objective -- a plausible number from a model the
# user did not ask for, which is the failure mode this package keeps meeting.
.admRefuseNodeStudies <- function(studies) {
  bad_w <- vapply(studies, function(s) {
    x <- s[["weight"]]
    !is.null(x) && is.numeric(x) && length(x) && !isTRUE(all.equal(as.numeric(x)[[1L]], 1))
  }, logical(1))
  bad_m <- vapply(studies, function(s) {
    x <- s[["cov_method"]]
    !is.null(x) && !identical(as.character(x)[[1L]], "marginal")
  }, logical(1))
  if (!any(bad_w) && !any(bad_m)) return(invisible(NULL))
  which_f <- if (any(bad_m)) "cov_method" else "weight"
  stop("admixr2: node-quadrature covariate marginalisation (gl / gh / taylor) ",
       "was removed.\n",
       "  Stud", if (sum(bad_w | bad_m) > 1L) "ies " else "y ",
       paste(sQuote(names(studies)[bad_w | bad_m]), collapse = ", "),
       " set `", which_f, "`.\n",
       "  For ONE pooled (E, V) per study -- what a publication reports -- give ",
       "the study its\n  `cov_dist` and let admixr2 marginalise over it.\n",
       "  For summaries BY COVARIATE STRATUM, pass the strata as ordinary ",
       "studies, each with\n  its own `cov` and its own `n`: that is the same ",
       "likelihood the node route computed,\n  with the real stratum sizes ",
       "instead of quadrature weights standing in for them.",
       call. = FALSE)
}
