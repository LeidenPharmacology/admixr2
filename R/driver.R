# Shared machinery for the four nlmixr2Est.* drivers.

# Apply per-fit settings from control to pinfo (forwarded by value to workers).
.admDriverPinfo <- function(.ui, .ctl) {
  pinfo <- .admParseIniDf(.ui$iniDf, .ui)
  pinfo$nDisplayProgress <- .ctl$nDisplayProgress %||% pinfo$nDisplayProgress
  pinfo$sigdig           <- .ctl$sigdig
  # Compiled-model cache path for workers (.admModelCacheFile).
  pinfo$sim_cache_file   <- tryCatch(.admModelCacheFile(.ui), error = function(e) NULL)
  # Residual-quadrature nodes travel pinfo -> arr -> .admResidApply/.admResidDeriv.
  pinfo$resid_nodes      <- .ctl$resid_nodes %||% .ADM_TBS_NODES
  # Covariate integration node count for adgh product grid.
  pinfo$cov_nodes        <- .ctl$cov_nodes %||% 7L
  pinfo$n_nodes          <- .ctl$n_nodes
  # Covariate integration method: "on", "off", or "sparse" (adgh only).
  pinfo$cov_integration  <- .ctl$cov_integration %||% "off"
  # Smolyak level for sparse covariate grid.
  pinfo$cov_sparse_level <- .ctl$cov_sparse_level %||% 3L
  pinfo
}


# Driver preamble: validate study names, materialise specs, and check identifiability.
.admDriverStudies <- function(.ui, .ctl, est) {
  studies <- .ctl$studies
  if (length(studies) == 0L)
    stop(est, "Control(studies=...) required", call. = FALSE)
  # Validate study names and ensure non-blank, unique names.
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
  # Materialise before shared validation reads study fields. The analysis
  # model's covariates go in: a source is cut into nodes only along the
  # directions this model can be moved by, which leaves the objective unchanged
  # -- see .admMaterialise().
  studies <- .admMaterialise(
    studies, analysis_covs = .admAllCovs(.ui, NULL))
  pinfo <- .admDriverPinfo(.ui, .ctl)
  .admWarnCovIdentifiability(.ui, pinfo, studies)
  list(studies = studies, pinfo = pinfo)
}

# Normalise studies, flatten to observation units, and derive model/joint flags.
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
# Caller invariants: snapshot `cov_nms` before nlmixr2CreateOutputFromUi
# (re-dimnames in place); compute `.admCovSkip` from this fit's cov
# (parFixedDf$SE aligns positionally); preserve `.foceiEnv` on class rewrite.
# `extra_field` ("admExtra"/"adirmcExtra") stays split -- it's user-visible on
# fit$env (plot.admFit resolves both), so unifying it would break the interface.
#
# Dummy frame for post-fit solve: fills model covariates with finite values
# (mean across studies, else 1) so rxode2 does not reject missing parameters.
.admDummyData <- function(.ui, multi_out, studies) {
  d  <- if (multi_out) admData(.admEndpointNames(.ui)) else admData()
  cv <- .admAllCovs(.ui, NULL)
  for (nm in cv) {
    # s$cov may be a list or named numeric vector
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
  # THE RESOLUTION THIS FIT WAS BUILT AT, per covariate THIS MODEL READS, for
  # anova()'s compatibility check. Stamped on .fit$env (not .ret, where env does
  # not yet exist as an environment).
  #
  # ONE NUMBER PER FIT WAS TOO COARSE, and it refused a comparison that is
  # exactly valid. A source is cut into nodes only along the covariates the
  # analysis model reads, so the NULL model of a nested pair -- which has
  # dropped the term -- legitimately has fewer nodes, or none. Measured: the
  # full fit at 10 studies against a null fit at 2 gives dOFV 409.162568, and
  # the same null refitted at the full resolution gives 409.162615, a difference
  # of 5e-05 -- while the old check refused the first pair outright because one
  # fit stamped `5` and the other stamped nothing.
  #
  # `"0"` for a covariate the model reads that nothing is conditional on --
  # integrated over the whole distribution, which is NOT `strata_nodes = 1L`
  # (one node pinned at the median). Stamping both `1` made anova() accept that
  # pair. And the whole MULTISET, `"3/9"`, since max() could not tell two
  # sources at 3 and 9 from both at 9; sorted, so study order cannot matter.
  .Jn <- local({
    .cv <- .admAllCovs(.ui)
    .cv <- .cv[vapply(.cv, function(cv) any(vapply(studies, function(s)
      cv %in% c(.admCovSpecNames(s[["cov_dist"]]),
                s[[".adm_strata_covs"]] %||% character(0)),
      logical(1))), logical(1))]
    if (!length(.cv)) return(character(0))
    vapply(.cv, function(cv) {
      j <- unlist(lapply(studies, function(s)
        if (cv %in% (s[[".adm_strata_covs"]] %||% character(0)))
          s[[".adm_strata_nodes"]] else NULL))
      if (!length(j)) "0"
      else paste(sort(unique(as.integer(j))), collapse = "/")
    }, character(1))
  })
  nlmixr2est::.nlmixr2FitUpdateParams(.ret)
  handle_ctl(.ctl, .ret)
  if (exists("control", .ui)) rm(list = "control", envir = .ui)
  .ret$control <- .admToFoceiControl(.ctl, .admCovSkip(cov, .ui))
  # Pre-fetch foceiModel to prevent re-compilation during output creation
  .focei_model <- suppressMessages(tryCatch(.ui$foceiModel, error = function(e) NULL))
  if (!is.null(.focei_model)) .ret$model <- .focei_model

  .fit <- nlmixr2est::nlmixr2CreateOutputFromUi(
    .ui, data = .admDummyData(.ui, multi_out, studies),
    control = .ret$control,
    table = .ret$table, env = .ret, est = est)

  .fit$env$method <- est
  if (length(.Jn)) .fit$env$strataNodes <- .Jn
  # Stamp quadrature (adgh) and MC (admc/adirmc) resolution so anova() rejects
  # comparisons across mismatched grids or sample sizes. `[[` avoids partial matching.
  .fit$env$nNodes <- .ctl[["n_nodes"]]
  .fit$env$nSim <- .ctl[["n_sim"]]
  .admRestoreCovNames(.fit, cov_nms)
  .fit$env$studies <- studies
  .extra <- .ret[[extra_field]]
  .fit$env[[extra_field]] <- .extra
  # nlmixr2-style parameter history for traceplot()
  .admAttachParHist(.fit, .extra$all_traces, .extra$par_names, .ui)
  # Observed + predicted aggregate moments (E vector, V matrix) per study
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
