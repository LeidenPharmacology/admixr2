# EXPERIMENTAL: joint individual + aggregate estimation (issue #60).
#
# The joint likelihood across independent studies is a product, so on the -2LL
# scale it is a sum:
#
#   -2LL(theta) = sum over AGGREGATE studies of the marginal MVN term
#               + sum over INDIVIDUAL subjects of the usual NLME term
#
# The asymmetry between the two is correct rather than a compromise: you
# CONDITION on a subject's covariates when you observe them, and you MARGINALISE
# over the study's covariate distribution when you do not. This is the framework
# Phillippo et al. formalise as multilevel network meta-regression.
#
# Why it earns its place: every covariate coefficient an aggregate-only fit can
# estimate comes from a BETWEEN-study contrast, and is therefore confounded with
# whatever else differs between those studies -- a true 0.75 comes back at 1.14
# on a three-trial design whose baselines happen to track mean weight. Individual
# records carry a WITHIN-study contrast, which is the only thing that breaks
# that confounding.
#
# The individual term is delegated to nlmixr2est's FOCEI rather than reimplemented.
# It is evaluated, not optimised: `maxOuterIterations = 0` runs the inner EBE
# problem at the parameters admixr2's optimizer is currently proposing and
# returns the objective there. So the outer optimisation stays admixr2's, and the
# individual contribution is a black-box function of the same parameter vector.
#
# THE INDIVIDUAL RECORDS MUST NOT BE THE SAME SUBJECTS AS AN AGGREGATE STUDY.
# Summing the two terms is a PRODUCT of likelihoods, which is correct only when
# the two sources are independent -- i.e. the individual records come from a
# study that is NOT also present in `studies`. Supply both representations of the
# same trial and those subjects are counted twice: the objective is not a
# likelihood of anything, and every standard error shrinks by about sqrt(2) with
# no data added. Haneuse & Wakefield (2008, JRSS-B 70(1):73-93) treat the
# same-population case by CONDITIONING instead -- the individual sample shrinks
# the support of the aggregate likelihood rather than multiplying an independent
# term -- which is a different construction from this one. Until that is
# implemented the restriction stands: `ipd` is for studies not in `studies`. If
# you have both for one trial, either drop the aggregate summary or reduce its
# `n` by the number of individual subjects.
#
# TWO THINGS INDIVIDUAL DATA DOES NOT FIX, both measured in the ecological
# literature and worth knowing before reaching for it:
#   * It does not rescue a MISSPECIFIED model -- it makes it worse. With an
#     ignored covariate interaction, adding individual records moved the bias
#     from 26.6% to 39.1% and coverage to 70% (Jackson, Best & Richardson 2006,
#     Stat Med 25:2136-2159). Individual data sharpens whatever model you have.
#   * For a CONTINUOUS covariate it may add little over the published SD. In the
#     same study, supplying the within-study variance alone moved the bias from
#     -7% to -3%, with "only a little extra benefit" from coupled individual
#     records. The `cov_dist` route is doing most of the work already.
#
# If individual records ARE available for only a handful of subjects, take the
# covariate EXTREMES rather than a random sample: five well-chosen subjects can
# beat fifty random ones (Glynn, Wakefield, Handcock & Richardson 2008, JRSS-A
# 171(1):179-202). Such a sample cannot be analysed on its own -- it is
# identifiable only jointly, which is what this path is for.
#
# Two things this costs, both accepted deliberately for an experimental path:
#
#   * SPEED. Every objective evaluation runs a full inner optimisation over every
#     individual subject. There is no gradient through it, so the joint fit is
#     derivative-free and needs more evaluations than it would otherwise.
#   * The likelihood CONSTANTS differ between the two terms. FOCEI's objective
#     carries the log(2*pi) terms that admData()'s all-NA DV column deliberately
#     keeps out of the aggregate one. Those constants do not depend on the
#     parameters, so they shift the objective by a fixed amount and leave the
#     argmin -- the estimates -- untouched. They would matter for comparing the
#     joint objective against anything else, so the value reported by a joint fit
#     is NOT comparable with an aggregate-only one. See issue #60.

# A private, reparsed copy of the model for the individual term.
#
# The inits are overwritten on every evaluation, so this must not be the caller's
# ui: an rxUi is an environment and mutating it in place would rewrite the model
# the user handed in. Reparsing from `$fun` is the one construction that is
# unambiguously a fresh object.
.admIpdSetup <- function(ui, ipd) {
  if (!is.data.frame(ipd) || !nrow(ipd))
    stop("admixr2: `ipd` must be a non-empty data frame of individual records.",
         call. = FALSE)
  need <- c("ID", "TIME", "DV")
  if (!all(need %in% names(ipd)))
    stop("admixr2: `ipd` needs column(s) ",
         paste(setdiff(need, names(ipd)), collapse = ", "), ".", call. = FALSE)
  fun <- tryCatch(ui$fun, error = function(e) NULL)
  if (is.null(fun))
    stop("admixr2: cannot take a private copy of the model for the individual ",
         "term (no `$fun`).", call. = FALSE)
  suppressMessages(rxode2::rxode2(fun))
}

# The individual -2LL at the parameters the outer optimizer is proposing.
#
# Returns Inf rather than propagating an error: FOCEI's inner problem can fail at
# a parameter vector the outer search is merely probing, and an Inf is a value the
# optimizer can step away from, where an error would abort the whole fit.
.admIpdFit <- function(pars, pinfo, ui_ipd, ipd, inner_iter = 100L,
                       fast = FALSE) {
  df <- pinfo$iniDf
  est <- tryCatch(.admFullTheta(pars, pinfo), error = function(e) NULL)
  if (is.null(est) || !all(is.finite(est))) return(NULL)
  df$est <- unname(est)
  ui_ipd$iniDf <- df
  tryCatch(suppressWarnings(suppressMessages(
    nlmixr2est::nlmixr2(
      ui_ipd, ipd, est = "focei",
      control = nlmixr2est::foceiControl(
        maxOuterIterations = 0L,          # EVALUATE at these parameters
        maxInnerIterations = as.integer(inner_iter),
        covMethod = "", calcTables = FALSE, print = 0L,
        compress = FALSE, addProp = "combined2", fast = fast))
  )), error = function(e) NULL)
}

.admIpdNLL <- function(pars, pinfo, ui_ipd, ipd, inner_iter = 100L) {
  f <- .admIpdFit(pars, pinfo, ui_ipd, ipd, inner_iter)
  if (is.null(f)) return(Inf)
  v <- tryCatch(as.numeric(f$objective), error = function(e) NA_real_)
  if (!is.finite(v)) Inf else v
}

# Gradient of the individual term.
#
# The joint objective is a SUM, so its gradient is the sum of the two and the
# aggregate half stays exactly as analytic as it ever was.
#
# The individual half is analytic for the structural thetas and the residual
# parameters, and finite-differenced for omega only. That split is not a
# compromise on the theta side -- nlmixr2est's own analytic FOCEI gradient
# supplies those, and they agree with central differences to ~1e-5 relative.
#
# HOW IT IS REACHED. nlmixr2est 6.2.0 still ships the R implementation of the
# analytic outer gradient (`.foceiAnalyticGradSetup` + `.foceiAnalyticGradCore`);
# 7.x deleted the R route in favour of a C++ one driven from inside the fit
# (commit "delete the R analytic gradient"), so there is no equivalent entry
# point there. Hence the version gate: on anything that does not expose the R
# core, this falls back to differencing the whole vector.
#
# Three things this depends on, each of which silently returns NULL if missed:
#   * the ui must carry `fast = TRUE` -- .foceiAnalyticGradSetup gates on it
#     first, so a fit run without it looks exactly like "no analytic gradient
#     available";
#   * the augmented model must be built with the SAME direction set the setup
#     returns, or the core refuses on an ndir mismatch;
#   * a proportional/combined residual declines when a prediction vanishes
#     (min|f| < 1e-6 max|f|), which a long terminal sampling tail will trigger.
#     That is a real refusal rather than a bug, and it is why the FD path stays.
#
# WHY OMEGA IS DIFFERENCED. The analytic vector reports omega as `om.chol.*`, on
# nlmixr2est's Cholesky estimation scale; admixr2's optimizer holds
# log(Omega_ii) on the diagonal and raw L_ij off it. Rotating between the two is
# the step most likely to be silently wrong -- a wrong factor there still
# converges, to the wrong place -- so until it is derived and FD-verified per
# block, those few entries are differenced directly. Omega is typically 1-3
# parameters, so this is 2-6 inner solves rather than 2*npar.
# Memo for the augmented sensitivity model. Named .adm* so the daemon payload
# filter would catch it, and kept OUT of anything serialised: a deserialised
# rxode2 model carries a dead pointer.
.adm_ipd_env <- new.env(parent = emptyenv())

.admIpdAnalyticGrad <- function(fit) {
  ns <- asNamespace("nlmixr2est")
  gt <- function(nm) if (exists(nm, envir = ns, inherits = FALSE))
    get(nm, envir = ns) else NULL
  setup <- gt(".foceiAnalyticGradSetup"); core <- gt(".foceiAnalyticGradCore")
  aug   <- gt(".foceiAnalyticAugModelDirs")
  if (is.null(setup) || is.null(core) || is.null(aug)) return(NULL)
  # The augmented model depends on the MODEL and the direction set, not on the
  # parameters, so it is built once and reused. Rebuilding it per gradient means
  # a symengine expansion and a gcc compile on every optimizer step, which made
  # the analytic route slower than differencing the whole vector.
  memo <- .adm_ipd_env
  tryCatch({
    ui <- fit$finalUi; Om <- fit$omega; ini <- ui$iniDf
    th <- ini[!is.na(ini$ntheta), , drop = FALSE]
    th <- th[order(th$ntheta), , drop = FALSE]
    v  <- fit$theta[th$name]; if (anyNA(v)) v <- th$est
    thVals <- stats::setNames(as.numeric(v), th$name)
    st <- setup(ui, thVals, Om)
    if (is.null(st)) return(NULL)
    key <- paste0(.admPkgKey(), "|", digest::digest(list(ui$lstExpr, st$dir$dirs)))
    am  <- memo[[key]]
    if (is.null(am)) {
      am <- aug(ui, st$dir$dirs)
      if (is.null(am)) return(NULL)
      assign(key, am, envir = memo)
    }
    tol <- tryCatch(gt(".foceiGradSolveTolOr")(ui), error = function(e) 1e-10)
    r <- core(ui, stats::setNames(as.numeric(thVals),
                                  paste0("THETA_", seq_along(thVals), "_")),
              as.matrix(fit$eta[, st$etaNames, drop = FALSE]), fit$eta$ID,
              fit$dataSav, Om, st$ef, st$dir, st$dOiEst, st$tr28, st$omNames,
              tol, interaction = st$interaction, foceType = st$foceType,
              startedEnv = NULL, am = am)
    if (is.null(r) || !all(is.finite(r$g))) return(NULL)
    r$g
  }, error = function(e) NULL)
}

.admIpdGrad <- function(p, pinfo, ui_ipd, ipd, h = NULL, inner_iter = 100L) {
  f <- function(pp) .admIpdNLL(.admUnpack(pp, pinfo), pinfo, ui_ipd, ipd,
                               inner_iter)
  fd_idx <- function(idx, base = NULL) {
    steps <- tryCatch(.admShi21Steps(f, p, idx = idx), error = function(e) NULL)
    if (is.null(steps) || !all(is.finite(steps[idx])))
      steps <- pmax(abs(p), 1) * (h %||% 1e-4)
    out <- if (is.null(base)) numeric(length(p)) else base
    for (k in idx) {
      hk <- steps[[k]]
      a <- b <- p; a[k] <- a[k] + hk; b[k] <- b[k] - hk
      fa <- f(a); fb <- f(b)
      out[k] <- if (is.finite(fa) && is.finite(fb)) (fa - fb) / (2 * hk) else 0
    }
    out
  }
  # the analytic route needs a fit at THESE parameters, with fast = TRUE
  fit <- .admIpdFit(.admUnpack(p, pinfo), pinfo, ui_ipd, ipd, inner_iter,
                    fast = TRUE)
  g <- if (is.null(fit)) NULL else .admIpdAnalyticGrad(fit)
  if (is.null(g)) return(fd_idx(seq_along(p)))       # refused: difference it all

  out <- numeric(length(p))
  nmT <- pinfo$struct_names
  nmS <- pinfo$sigma_names
  # struct thetas: iniDf holds the optimizer value itself (log scale), factor 1
  for (j in seq_along(nmT)) {
    k <- match(nmT[j], names(g))
    if (!is.na(k)) out[j] <- g[[k]]
  }
  # sigma: iniDf reports an SD, the optimizer holds log(sigma^2), so
  # d/d(log sigma^2) = d/d(sigma) * sigma/2 for a variance-role sigma
  role <- .admSigmaRole(pinfo)
  sv   <- .admUnpack(p, pinfo)$sigma_var
  off  <- length(nmT)
  for (j in seq_along(nmS)) {
    k <- match(nmS[j], names(g))
    if (is.na(k)) next
    out[off + j] <- if (identical(role[[j]], "var"))
      g[[k]] * sqrt(unname(sv[nmS[j]])) / 2 else g[[k]]
  }
  # omega: analytic reports om.chol.* on nlmixr2est's Cholesky scale; difference
  # those entries on admixr2's own scale instead of rotating between the two
  om_idx <- setdiff(seq_along(p), seq_len(off + length(nmS)))
  if (length(om_idx)) out <- fd_idx(om_idx, base = out)
  out
}

# Observations an individual dataset contributes to the joint objective.
#
# Counted as records the likelihood actually scores -- EVID == 0 with a non-NA
# DV -- rather than nrow(), which would also count dose rows.
.admIpdNobs <- function(ipd) {
  if (!is.data.frame(ipd) || !nrow(ipd)) return(0L)
  keep <- rep(TRUE, nrow(ipd))
  if (!is.null(ipd$EVID)) keep <- keep & ipd$EVID == 0
  if (!is.null(ipd$DV))   keep <- keep & !is.na(ipd$DV)
  sum(keep)
}
