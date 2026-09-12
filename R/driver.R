# Shared machinery for the four nlmixr2Est.* drivers.
#
# admc, adfo, adgh and adirmc differ in how they compute an objective and a
# gradient, and in almost nothing else. Everything that happens once the optimum
# is in hand -- assembling the fit environment, handing it to nlmixr2est, and
# repairing what nlmixr2est does to it on the way back -- was written out four
# times, token for token, differing only in an estimator name and a field name.

# Apply the per-fit settings that every driver copies from its control onto
# `pinfo`, and return the endpoint name alongside.
#
# These four fields travel on `pinfo` rather than as arguments for one reason,
# recorded here once instead of four times: a mirai daemon resolves the restart
# worker from the STALE INSTALLED namespace, so a new formal on
# `.admRestartWorker` (or its three siblings) throws `unused argument` before the
# patched dev body can run -- but `pinfo` is sent by value and carries anything.
#
# Single-sourced because the set GROWS. Each of `sigdig`, `sim_cache_file` and
# `resid_nodes` was added to four drivers by hand, and a field set in three of
# them is not an error anywhere -- it is a silent divergence in whichever
# estimator was missed (`sim_cache_file` decides which compiled model a parallel
# worker reads; `resid_nodes` decides the accuracy a transform-both-sides
# residual is integrated at).
#
# adirmc adds `grad_h` of its own after this returns; it is the only
# estimator whose inner optimiser finite-differences on a separate path.
.admDriverPinfo <- function(.ui, .ctl) {
  pinfo <- .admParseIniDf(.ui$iniDf, .ui)
  pinfo$nDisplayProgress <- .ctl$nDisplayProgress %||% pinfo$nDisplayProgress
  pinfo$sigdig           <- .ctl$sigdig
  # The compiled-model cache path, for the parallel workers. They have no `ui`,
  # so they cannot derive the .admIniKey() component that keys a fix()ed
  # parameter's VALUE -- reading the wrong entry means solving at another model's
  # fixed value, silently. See .admModelCacheFile().
  pinfo$sim_cache_file   <- tryCatch(.admModelCacheFile(.ui), error = function(e) NULL)
  # Residual-quadrature nodes travel pinfo -> arr -> .admResidApply/.admResidDeriv.
  pinfo$resid_nodes      <- .ctl$resid_nodes %||% .ADM_TBS_NODES
  # Nodes for the COVARIATE dimension of the adgh product grid. Read in exactly
  # one place (.admCovGrid via .adghGrid) and, until this line existed, set
  # nowhere -- so the `%||% 7L` fallback there always fired and the covariate
  # integration was fixed at 7 nodes with no way to change it. That is what made
  # the accuracy sweep plateau: raising n_nodes refines only the random-effect
  # dimensions, so the covariate error never moved.
  pinfo$cov_nodes        <- .ctl$cov_nodes %||% 7L
  # NULL for the three estimators that have no eta grid; the joint collapse
  # needs it to price itself against the design it would replace
  pinfo$n_nodes          <- .ctl$n_nodes
  # How the covariate distribution is integrated: "on" (the product grid above,
  # reduced wherever a reduction verifies), "off" (that grid, unreduced) or
  # "sparse" (a Smolyak grid over the same margins). Only adghControl() exposes
  # it, so every other estimator lands on "off" here -- the reductions are
  # quadrature constructions and the other estimators reach their etas through
  # .admMakeZ, so "off" is not a default, it is the only correct value for them.
  pinfo$cov_integration  <- .ctl$cov_integration %||% "off"
  # 3, matching adghControl's default and .adghGrid's own fallback. Level 2 is
  # the axial rule the retired Taylor design was; adghControl always supplies
  # the value so this fallback is unreachable today, but a control that did not
  # would otherwise have built a rule an order of magnitude coarser in silence.
  pinfo$cov_sparse_level <- .ctl$cov_sparse_level %||% 3L
  pinfo
}


# The preamble every driver runs before it touches the studies: name them, refuse
# the retired node grammar, build pinfo, and warn on an unidentifiable covariate
# coefficient. Copy-pasted into all four before this, comments included.
#
# ORDER MATTERS TWICE, which is why it is one function rather than a convention.
.admDriverStudies <- function(.ui, .ctl, est) {
  studies <- .ctl$studies
  if (length(studies) == 0L)
    stop(est, "Control(studies=...) required", call. = FALSE)
  # PARTIAL NAMES ARE THE DANGEROUS CASE, and `is.null(names(studies))` is
  # FALSE for them. `list(pub = admStudy(...), list(E = ..., n = 200))` has
  # names c("pub", ""), so it passed straight through -- and .admMaterialise()
  # rebuilds by name, where `studies[[""]]` is NULL rather than an error and
  # `out[[""]] <- NULL` adds nothing. A two-study fit quietly became a
  # one-study fit, printed "Obs units: 1" and reported estimates from half the
  # data. Duplicates collapse to the first the same way. admStudies() already
  # guards blanks, NAs and duplicates; the plain list() route did not.
  .nm <- names(studies)
  if (is.null(.nm)) .nm <- rep("", length(studies))
  .nm[is.na(.nm)] <- ""
  .blank <- !nzchar(.nm)
  if (any(.blank)) .nm[.blank] <- paste0("study", which(.blank))
  if (anyDuplicated(.nm))
    stop(est, "Control(studies=): duplicate study name(s) ",
         paste(sQuote(unique(.nm[duplicated(.nm)])), collapse = ", "),
         ". Studies are matched by name throughout, so a duplicate silently ",
         "drops all but the first. Give each study its own name.",
         call. = FALSE)
  names(studies) <- .nm
  # Materialise before shared validation reads study fields.
  studies <- .admMaterialise(studies)
  pinfo <- .admDriverPinfo(.ui, .ctl)
  .admWarnCovIdentifiability(.ui, pinfo, studies)
  list(studies = studies, pinfo = pinfo)
}

# Normalise the studies, flatten them to observation units, and report the two
# model-level flags every driver derives from the result.
#
# `multi_out` is MODEL-level (`length(.admOutputVars(ui)) > 1`), not study-level,
# so a multi-endpoint model always tags observations by `cmt` -- the solve cannot
# disambiguate endpoints otherwise. `any_joint` is study-level and gates the
# paths that have no joint implementation.
#
# `tag_cmt` defaults to `multi_out`, which is what three drivers pass. adirmc
# refuses both multi-output and joint studies, and does so BETWEEN the flatten and
# the ev_full build -- so it calls this with `ev_full = FALSE`, runs its refusal,
# and builds ev_full itself. Keeping that as a flag rather than moving the refusal
# is deliberate: an error must fire where it fires now.
.admDriverUnits <- function(studies, .ui, output_var, ev_full = TRUE) {
  for (nm in names(studies))
    studies[[nm]] <- .admNormaliseStudy(studies[[nm]], nm, output_var)
  studies   <- .admFlattenStudies(studies)
  multi_out <- length(.admOutputVars(.ui)) > 1L
  any_joint <- any(vapply(studies, function(u) isTRUE(u$is_joint), logical(1)))
  if (ev_full) studies <- .admBuildEvFull(studies, tag_cmt = multi_out)
  list(studies = studies, multi_out = multi_out, any_joint = any_joint)
}

# Turn a fully-populated `.ret` environment into the returned nlmixr2 fit.
#
# The variation points are ARGUMENTS, never inferred:
#
#   est          "admc" / "adfo" / "adgh" / "adirmc" -- the est string, the
#                $method field, and the objDf row name, which must agree
#   objective    opt$objective for the nloptr drivers, best_nll for adirmc
#   extra_field  "admExtra", or "adirmcExtra" for adirmc. NOT unified: it is a
#                user-visible field on fit$env (plot.admFit resolves both), so
#                renaming it here would be a silent interface change for anyone
#                reading fit$env$adirmcExtra
#   handle_ctl   the nmObjHandleControlObject.* METHOD, passed as a function.
#                Deliberately not generic dispatch: the generic is registered
#                into nlmixr2est's namespace by zzz.R and is not imported here,
#                so dispatching would make this depend on load order for no gain
#   multi_out    adirmc rejects multi-endpoint models, so it passes FALSE and
#                gets admData() exactly as before
#
# ORDER IS PART OF THE CONTRACT, not a stylistic choice:
#
#   * `cov_nms` must be SNAPSHOT BY THE CALLER, before nlmixr2CreateOutputFromUi
#     runs. foceiFitCpp_ re-dimnames the covariance IN PLACE on the same SEXP
#     (`as<NumericMatrix>` wraps, it does not copy), so reading the names back
#     off .cov afterwards restores nothing. This silently no-opped once already.
#   * `.admCovSkip(cov, .ui)` must be computed from the covariance THIS fit
#     built, because parFixedDf$SE is filled positionally against skipCov.
#   * the class rewrite must carry `.foceiEnv` across, or the fit loses the
#     environment nlmixr2est's own methods reach through.
# The dummy frame the post-fit output stage solves over.
#
# A model that READS a covariate needs that covariate present here, or rxode2
# stops with "the following parameter(s) are required for solving: <name>" --
# AFTER the whole optimisation has run, so the fit is lost at the last step. The
# frame is built by the driver, not taken from the caller, so adding the column
# to the data passed to nlmixr2() does not help; it has to happen here.
#
# The VALUE is immaterial to the result -- this row carries the non-NA
# placeholder DV and never enters the reported objective, which every estimator
# overwrites with its own aggregate -2LL. It must merely be finite and not
# something the model divides by, so use the studies' own covariate values when
# they carry any and fall back to 1 rather than 0.
.admDummyData <- function(.ui, multi_out, studies) {
  d  <- if (multi_out) admData(.admEndpointNames(.ui)) else admData()
  cv <- tryCatch(.ui$allCovs, error = function(e) NULL)
  for (nm in cv) {
    # [[ ]] on an ATOMIC s$cov is "subscript out of bounds", raised after the
    # whole fit has run. A study may carry its covariates as a named numeric
    # vector as easily as a list, so read both.
    vals <- unlist(lapply(studies, function(s) {
      cs <- s[["cov"]]
      if (is.null(cs) || !nm %in% names(cs)) NULL else cs[[nm]]
    }), use.names = FALSE)
    vals <- vals[is.finite(vals)]
    d[[nm]] <- if (length(vals)) mean(vals) else 1
  }
  d
}


.admFinaliseFit <- function(.ret, .ui, .ctl, est, objective, ov, studies,
                            cov, cov_nms, multi_out, extra_field, handle_ctl,
                            t_opt, t_cov, t_elapsed) {
  # The stratum resolution the studies were generated at, if any. Recorded
  # HERE because all four drivers pass through this function -- a field set in
  # three of four is a silent divergence, not an error. The objective is
  # J-dependent, so anova() refuses to compare fits that disagree on it.
  # Stamped on .fit below, NOT on .ret: `.ret$env` does not exist, so assigning
  # into it creates a plain LIST named `env` and the field is unreachable --
  # `fit$env` dispatches to nlmixr2est's nmObjGet.env, which returns the fit
  # environment. anova()'s refusal was dead code for exactly that reason.
  # The whole set is recorded, not just a length-1 one, so a mixed-resolution
  # fit is comparable only with another built the same way.
  .Jn <- sort(unique(unlist(lapply(studies, function(s) s[[".adm_strata_nodes"]]))))
  nlmixr2est::.nlmixr2FitUpdateParams(.ret)
  handle_ctl(.ctl, .ret)
  if (exists("control", .ui)) rm(list = "control", envir = .ui)
  .ret$control <- .admToFoceiControl(.ctl, .admCovSkip(cov, .ui))
  # Pre-fetch so nlmixr2CreateOutputFromUi does not compile the model a second
  # time (rxUiGet.foceiOptEnv finds it on .ret via the foceiEnv mechanism).
  .focei_model <- suppressMessages(tryCatch(.ui$foceiModel, error = function(e) NULL))
  if (!is.null(.focei_model)) .ret$model <- .focei_model

  .fit <- nlmixr2est::nlmixr2CreateOutputFromUi(
    .ui, data = .admDummyData(.ui, multi_out, studies),
    control = .ret$control,
    table = .ret$table, env = .ret, est = est)

  .fit$env$method <- est
  if (length(.Jn)) .fit$env$strataNodes <- .Jn
  # THE QUADRATURE RESOLUTION, stamped so anova() can refuse to difference two
  # objectives computed on different grids. adgh's objective is a quadrature
  # approximation and its VALUE moves with the node count, so a fit at
  # n_nodes = 5 and the same model at 7 are not evaluating the same number --
  # differencing them gives something that looks like a dOFV and is not. This
  # is .adghGrid's own rule (it refuses to change the point count mid-fit)
  # extended across fits.
  #
  # `[[ ]]`, not `$`: `$` partial-matches on a list, and a control gaining any
  # longer `n_nodes*` field would silently stamp the wrong number. NULL for the
  # three estimators that have no node grid, which makes the check a no-op
  # there rather than a special case.
  .fit$env$nNodes <- .ctl[["n_nodes"]]
  # THE MONTE CARLO RESOLUTION, stamped for the same reason: admc/adirmc's
  # objective is an average over `n_sim` draws, and that average moves with
  # `n_sim` exactly as the quadrature objective moves with `n_nodes`. NULL for
  # adfo/adgh, which makes the check below a no-op there.
  .fit$env$nSim <- .ctl[["n_sim"]]
  .admRestoreCovNames(.fit, cov_nms)
  .fit$env$studies <- studies
  .extra <- .ret[[extra_field]]
  .fit$env[[extra_field]] <- .extra
  # nlmixr2-style parameter history, so traceplot(fit) works natively.
  .admAttachParHist(.fit, .extra$all_traces, .extra$par_names, .ui)
  # Observed + predicted aggregate moments (E vector, V matrix) per study.
  .admAttachAggData(.fit, .extra, .ui)
  .old_cls <- class(.fit)
  .new_cls <- c("admFit", .old_cls)
  attr(.new_cls, ".foceiEnv") <- attr(.old_cls, ".foceiEnv")
  class(.fit) <- .new_cls

  .stats <- .admCalcObjStats(objective, length(ov$p0), studies)
  row.names(.stats$objDf) <- est
  .fit$env$logLik    <- .stats$ll
  .fit$env$nobs      <- .stats$nobs
  .fit$env$objDf     <- .stats$objDf
  .fit$env$OBJF      <- .stats$objDf$OBJF
  .fit$env$AIC       <- .stats$objDf$AIC
  .fit$env$BIC       <- .stats$objDf$BIC
  .fit$env$objective <- objective
  .fit$env$time      <- data.frame(
    optimize   = t_opt,
    covariance = t_cov,
    other      = 0,
    elapsed    = t_elapsed,
    row.names  = NULL
  )

  .fit
}
