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
.admIpdNLL <- function(pars, pinfo, ui_ipd, ipd, inner_iter = 100L) {
  df <- pinfo$iniDf
  est <- tryCatch(.admFullTheta(pars, pinfo), error = function(e) NULL)
  if (is.null(est) || !all(is.finite(est))) return(Inf)
  df$est <- unname(est)
  ui_ipd$iniDf <- df
  f <- tryCatch(suppressWarnings(suppressMessages(
    nlmixr2est::nlmixr2(
      ui_ipd, ipd, est = "focei",
      control = nlmixr2est::foceiControl(
        maxOuterIterations = 0L,          # EVALUATE at these parameters
        maxInnerIterations = as.integer(inner_iter),
        covMethod = "", calcTables = FALSE, print = 0L,
        compress = FALSE, addProp = "combined2"))
  )), error = function(e) NULL)
  if (is.null(f)) return(Inf)
  v <- tryCatch(as.numeric(f$objective), error = function(e) NA_real_)
  if (!is.finite(v)) Inf else v
}
