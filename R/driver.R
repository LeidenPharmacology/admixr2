# Shared machinery for the four nlmixr2Est.* drivers.
#
# admc, adfo, adgh and adirmc differ in how they compute an objective and a
# gradient, and in almost nothing else. Everything that happens once the optimum
# is in hand -- assembling the fit environment, handing it to nlmixr2est, and
# repairing what nlmixr2est does to it on the way back -- was written out four
# times, token for token, differing only in an estimator name and a field name.

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
.admFinaliseFit <- function(.ret, .ui, .ctl, est, objective, ov, studies,
                            cov, cov_nms, multi_out, extra_field, handle_ctl,
                            t_opt, t_cov, t_elapsed) {
  nlmixr2est::.nlmixr2FitUpdateParams(.ret)
  handle_ctl(.ctl, .ret)
  if (exists("control", .ui)) rm(list = "control", envir = .ui)
  .ret$control <- .admToFoceiControl(.ctl, .admCovSkip(cov, .ui))
  # Pre-fetch so nlmixr2CreateOutputFromUi does not compile the model a second
  # time (rxUiGet.foceiOptEnv finds it on .ret via the foceiEnv mechanism).
  .focei_model <- suppressMessages(tryCatch(.ui$foceiModel, error = function(e) NULL))
  if (!is.null(.focei_model)) .ret$model <- .focei_model

  .fit <- nlmixr2est::nlmixr2CreateOutputFromUi(
    .ui, data = if (multi_out) admData(.admEndpointNames(.ui)) else admData(),
    control = .ret$control,
    table = .ret$table, env = .ret, est = est)

  .fit$env$method <- est
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
