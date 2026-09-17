# Recursion guard: print.nlmixr2FitCore may call print(x) on the full fit.
# Without a guard print.admFit would recurse infinitely.
.adm_print_guard <- new.env(parent = emptyenv())
.adm_print_guard$active <- FALSE

#' Print method for admFit objects
#'
#' Delegates to `print.nlmixr2FitCore` for the standard nlmixr2 coloured
#' output.  `admFit` class is kept on the object during the call so that
#' `head.admFit` intercepts any `head(fit)` calls that arise in the paged-
#' output path (R Markdown / notebooks), preventing the
#' `[.data.frame(.subset2(env, integer))` crash that occurs when an
#' environment-backed fit is subscripted like a plain list.
#'
#' @param x An `admFit` object.
#' @param ... Passed to `print.nlmixr2FitCore`.
#' @return `x`, invisibly.
#'
#' @examples
#' \donttest{
#' library(rxode2)
#' library(nlmixr2)
#'
#' data("examplomycin")
#' obs    <- examplomycin[examplomycin$EVID == 0, ]
#' obs    <- obs[order(obs$ID, obs$TIME), ]
#' times  <- sort(unique(obs$TIME))
#' ids    <- unique(obs$ID)
#' dv_mat <- do.call(rbind, lapply(ids, function(i) {
#'   sub <- obs[obs$ID == i, ]; sub$DV[order(sub$TIME)]
#' }))
#' E <- colMeans(dv_mat)
#' V <- cov.wt(dv_mat, method = "ML")$cov
#'
#' pk_model <- function() {
#'   ini({
#'     tcl <- log(5); tv <- log(30)
#'     prop.sd <- c(0, 0.2)
#'     eta.cl ~ 0.09; eta.v ~ 0.04
#'   })
#'   model({
#'     cl <- exp(tcl + eta.cl)
#'     v  <- exp(tv  + eta.v)
#'     d/dt(central) <- -(cl/v) * central
#'     cp <- central / v
#'     cp ~ prop(prop.sd)
#'   })
#' }
#'
#' fit <- nlmixr2(
#'   pk_model, admData(), est = "adfo",
#'   control = adfoControl(
#'     studies = list(study1 = list(E = E, V = V, n = length(ids),
#'                                  times = times, ev = et(amt = 100))),
#'     maxeval = 100L
#'   )
#' )
#' print(fit)
#' }
#'
#' @export
print.admFit <- function(x, ...) {
  if (.adm_print_guard$active) {
    cl <- class(x)
    class(x) <- cl[cl != "admFit"]
    on.exit(class(x) <- cl)
    tryCatch(print(x, ...), error = function(e) invisible(NULL))
    return(invisible(x))
  }
  .adm_print_guard$active <- TRUE
  on.exit(.adm_print_guard$active <- FALSE, add = TRUE)
  saved_cl <- class(x)
  on.exit(tryCatch(class(x) <- saved_cl, error = function(e) NULL), add = TRUE)

  # getS3method(), not get(..., envir = asNamespace("nlmixr2est")).
  # print.nlmixr2FitCore is NOT exported, so reaching into the namespace for it is
  # semantically a ::: call that merely evades R CMD check's syntactic scan -- and
  # it carries the same fragility the no-::: policy exists to avoid. It IS
  # registered as an S3 method, so the method lookup is the supported public route
  # to the same function.
  fn <- utils::getS3method("print", "nlmixr2FitCore")
  class(x) <- class(x)[class(x) != "nlmixr2FitData"]
  fn(x, ...)
  invisible(x)
}

# knitr auto-print handler. Registered as knit_print.admFit into knitr's
# namespace by .register_knit_print() in .onLoad (zzz.R). Not named
# knit_print.admFit to avoid roxygen2 S3-method detection and spurious
# NAMESPACE requirements.
.admKnitPrint <- function(x, ...) {
  fn <- utils::getS3method("print", "nlmixr2FitCore")
  saved_cl <- class(x)
  on.exit(tryCatch(class(x) <- saved_cl, error = function(e) NULL))
  class(x) <- class(x)[class(x) != "nlmixr2FitData"]
  fn(x, ...)
  invisible(NULL)
}

# Called by print.nlmixr2FitCore when it ends its console branch with
# print(head(x)). Unlike rmarkdown's print.paged_df, this method is not in
# nlmixr2est's import chain, so S3 dispatch falls through to admixr2's method
# table and finds us first -- which is the supported way to intercept this path,
# and the reason no namespace mutation is needed (see #58).
# Converts to a plain data frame before head() to avoid .subset2(env, integer).
#' @method head admFit
#' @export
head.admFit <- function(x, n = 6L, ...) {
  cl <- class(x)
  class(x) <- cl[cl != "admFit"]
  on.exit(class(x) <- cl)
  tryCatch(head(as.data.frame(x), n = n), error = function(e) data.frame())
}

# Intercepts utils::head(x, n) inside rmarkdown's paged_table_html when x is an
# admFit environment whose class was replaced to paged_df by .pagedPrint.
# rmarkdown and nlmixr2est do not define head.paged_df, so S3 dispatch falls
# through to the search path and finds this method first.
#' @method head paged_df
#' @export
head.paged_df <- function(x, n = 6L, ...) {
  if (is.environment(x)) return(data.frame())
  NextMethod()
}

## Shared display spec for parameter-trace rendering.
##
## Returns the per-parameter display names, optimizer-scale -> natural-scale
## back-transforms, and iniDf-driven facet order used by both the custom
## `plot(fit, which = "par")` panel and the nlmixr2 `traceplot()` bridge
## (`.admBuildParHistData`). Keeping this in one place ensures both renderings
## label and scale parameters identically.
##
## - struct thetas: back-transformed via `ui$muRefCurEval` transform
## - omega diagonal (`log(Omega_ii)`): `exp()` -> variance, labelled `V(eta)`
## - omega off-diagonal (raw `L[i,j]`): identity, labelled `eta_i,eta_j`
## - sigma (`log(sigma^2)`): `exp(v/2)` -> SD
##
## Returns `NULL` when `pinfo` or `par_names` is unavailable.
.admTraceDisplaySpec <- function(pinfo, par_names, iniDf = NULL) {
  if (is.null(pinfo) || is.null(par_names)) return(NULL)

  disp_nms <- setNames(par_names, par_names)
  for (k in base::which(pinfo$chol_diag)) {
    nm <- pinfo$omega_par_names[k]
    disp_nms[[nm]] <- paste0("V(", pinfo$eta_names[pinfo$chol_i[k]], ")")
  }
  for (k in base::which(!pinfo$chol_diag)) {
    nm <- pinfo$omega_par_names[k]
    disp_nms[[nm]] <- paste0(pinfo$eta_names[pinfo$chol_i[k]], ",",
                             pinfo$eta_names[pinfo$chol_j[k]])
  }

  struct_nms        <- names(pinfo$struct_transforms)
  omega_diag_nms    <- pinfo$omega_par_names[pinfo$chol_diag]
  omega_offdiag_nms <- pinfo$omega_par_names[!pinfo$chol_diag]
  # Residual parameters go through .admSigmaReportFn() so the trace is on the
  # same scale print(fit) reports -- every sigma_role, not a list of the ones
  # that existed when this was written.
  back_fns <- setNames(lapply(par_names, function(nm) {
    if (nm %in% struct_nms)             function(v) .admBackTransform(v, pinfo$struct_transforms[[nm]])
    else if (nm %in% omega_diag_nms)    exp
    else if (nm %in% omega_offdiag_nms) identity
    else .admSigmaReportFn(pinfo, nm)
  }), par_names)

  # Order par_names by iniDf row position so facets follow ini() block order.
  # Omega params are positioned at their lead eta's iniDf row; off-diagonals
  # get a +0.5 fractional offset so they appear just after their diagonal.
  param_order <- if (!is.null(iniDf)) {
    ini_nms <- iniDf$name
    pos <- vapply(par_names, function(nm) {
      if (nm %in% pinfo$omega_par_names) {
        k <- match(nm, pinfo$omega_par_names)
        p <- match(pinfo$eta_names[pinfo$chol_i[k]], ini_nms)
        if (is.na(p)) Inf else p + if (!pinfo$chol_diag[k]) 0.5 else 0.0
      } else {
        p <- match(nm, ini_nms)
        if (is.na(p)) Inf else as.double(p)
      }
    }, double(1))
    unname(disp_nms[par_names[order(pos)]])
  } else NULL

  list(disp_nms = disp_nms, back_fns = back_fns, param_order = param_order)
}

## Build a nlmixr2-style `parHistData` frame from collected optimizer traces.
##
## nlmixr2's `traceplot()` generic reads `fit$parHistStacked`, which
## `nmObjGet.parHistStacked` derives from `fit$env$parHistData` -- a wide data.frame
## with a `type` column (it keeps `type == "Unscaled"`), an `iter` column, and one
## column per parameter. Populating that slot is all `traceplot(fit)` needs; no S3
## registration, because admFit already inherits `nlmixr2FitCore`.
##
## Semantics chosen here:
## - single chain = the best restart (lowest final NLL); nlmixr2's shape stores one
##   value per parameter per iter, so multi-restart overlay is not expressible --
##   that stays in `plot(fit, which = "par")`.
## - natural scale under `"Unscaled"`, using the same back-transforms and display
##   names as the custom par panel (`.admTraceDisplaySpec`).
## - no burn-in marker: `parHist`'s class carries no `niter` attribute, so
##   `traceplot()` draws no vline.
##
## The `iter` axis indexes improving optimizer evaluations, not raw nloptr
## iterations. Returns `NULL` when no usable trace is available.
.admBuildParHistData <- function(all_traces, par_names, ui) {
  if (is.null(all_traces) || length(all_traces) == 0L || is.null(par_names))
    return(NULL)

  # Best restart = lowest final NLL (NA traces ignored).
  finals <- vapply(all_traces, function(tr) {
    nt <- tr$nll_trace
    if (is.null(nt) || length(nt) == 0L) NA_real_ else nt[length(nt)]
  }, double(1))
  if (all(is.na(finals))) return(NULL)
  best <- all_traces[[which.min(finals)]]

  pt <- best$par_trace
  if (is.null(pt) || nrow(pt) == 0L) return(NULL)
  df <- as.data.frame(pt)
  if (ncol(df) != length(par_names)) return(NULL)
  colnames(df) <- par_names

  pinfo <- tryCatch(.admParseIniDf(ui$iniDf, ui), error = function(e) NULL)
  iniDf <- tryCatch(ui$iniDf, error = function(e) NULL)
  spec  <- .admTraceDisplaySpec(pinfo, par_names, iniDf)
  # Without the display spec we cannot back-transform to natural scale; emit no
  # parHistData rather than a raw optimizer-scale trace mislabelled "Unscaled".
  if (is.null(spec)) return(NULL)

  cols <- lapply(par_names, function(nm) as.numeric(spec$back_fns[[nm]](df[[nm]])))
  names(cols) <- vapply(par_names, function(nm) spec$disp_nms[[nm]], character(1))

  # Follow iniDf facet order when available so traceplot panels match the
  # custom par panel.
  if (!is.null(spec$param_order))
    cols <- cols[spec$param_order]

  data.frame(type = "Unscaled",
             iter = seq_len(nrow(df)),
             cols,
             check.names = FALSE,
             stringsAsFactors = FALSE)
}

## Attach a nlmixr2-style `parHistData` slot to a freshly constructed admFit when
## a usable trace is available. Shared by the admc/adfo/adgh/adirmc estimators so the
## binding logic lives in one place. A `NULL` build result must not be bound --
## `env$x <- NULL` still satisfies `exists()` and would leave a stale slot that
## `nmObjGet.parHistStacked` treats as present -- so we guard on non-NULL.
## `fit$env` is an environment, so the assignment is in place.
.admAttachParHist <- function(fit, all_traces, par_names, ui) {
  ph <- .admBuildParHistData(all_traces, par_names, ui)
  if (!is.null(ph)) fit$env$parHistData <- ph
  invisible(fit)
}

## Observed and predicted aggregate moments per study.
##
## Runs one MC simulation per study at the fitted parameters -- using the same
## quasi-random sampling and residual-error (sigma) handling as the diagnostic
## mean/cov panels -- and returns, per study, the observed and predicted mean
## vector `E` and (co)variance matrix `V`. Shared by `plot.admFit()` (mean/cov
## panels) and `.admAttachAggData()` (the fit's `aggData` slot) so the two never
## disagree.
##
## Returns a named list, one entry per study. Each entry is `NULL` when the
## simulation model is unavailable or the study simulation failed, otherwise:
##   list(times = <numeric>, n = <int>,
##        obs  = list(E = <named numeric>, V = <matrix>),
##        pred = list(E = <named numeric>, V = <matrix>))
## The observation times label `E` (names) and `V` (dimnames). `warn = TRUE`
## emits the user-facing warnings used on the interactive plot path; the fit
## attachment path passes `warn = FALSE` so a non-simulable fit stays quiet.
.admAggData <- function(extra, ui, n_sim = NULL, seed = 1L, warn = TRUE) {
  studies   <- extra$studies
  n_sim     <- n_sim %||% extra$n_sim %||% 5000L
  omega     <- extra$omega
  n_eta     <- nrow(omega)
  L         <- extra$L %||% tryCatch(t(chol(omega)), error = function(e) NULL)
  sv        <- extra$sigma_var
  eta_nms   <- extra$eta_col_names %||% character(0)
  sig_nms   <- names(sv)

  empty <- setNames(vector("list", length(studies)), names(studies))
  rxMod <- tryCatch(ui$simulationModel, error = function(e) NULL)
  if (is.null(rxMod)) {
    if (warn)
      warning("plot.admFit: could not retrieve simulation model from fit object",
              call. = FALSE)
    return(empty)
  }
  # Detect the simulation output variable (e.g. "ipredSim" for linCmt models)
  # rather than assuming "cp" -- matches the detection used on the fit path.
  out_var <- tryCatch(.admOutputVar(ui), error = function(e) "cp")
  # Re-parse the ui so the plotted bands use exactly the fit's residual error
  # model (per-endpoint spec, error form, sigma roles). Falls back to a
  # name-based guess only when the ui cannot be parsed.
  #
  # suppressWarnings: this is a RE-parse of a ui that was already parsed (and
  # already warned about) at fit time, so .admBuildResidSpecs()'s advisory
  # warnings -- e.g. an estimated t(nu) being non-identifiable from aggregate
  # data -- would otherwise be re-emitted on every plot() as if they were new.
  pinfo_r <- tryCatch(suppressWarnings(.admParseIniDf(ui$iniDf, ui)),
                      error = function(e) NULL)
  if (is.null(pinfo_r))
    pinfo_r <- list(sigma_names   = sig_nms,
                    sigma_output  = rep(NA_character_, length(sv)),
                    sigma_is_prop  = as.list(grepl("prop",  sig_nms, ignore.case = TRUE)),
                    sigma_is_lnorm = as.list(grepl("lnorm", sig_nms, ignore.case = TRUE)))
  # No cov_map is rebuilt here: every study uses the plain Omega and carries its
  # covariates as per-row data or as a shifted eta column, neither of which
  # needs one. .admStudyCovRows() still needs n_eta, so that stays.
  if (is.null(pinfo_r$n_eta)) pinfo_r$n_eta <- n_eta
  # .admParseIniDf() carries no resid_nodes -- only the DRIVERS set it, from the
  # control. Restore the count the fit actually used, or the diagnostics rebuild
  # V_pred on the 81-node default and a fit run with resid_nodes = 31L (or 201L)
  # is diagnosed against a different model than it was fitted with.
  pinfo_r$resid_nodes <- extra$resid_nodes %||% .ADM_TBS_NODES
  sig_output <- pinfo_r$sigma_output

  .sim_study <- function(s) {
    ov <- s$output %||% out_var
    tryCatch(rxode2::rxLoad(rxMod), error = function(e) NULL)
    set.seed(seed)
    if (n_eta > 0 && !is.null(L)) {
      .samp <- extra$sampling %||% "sobol"
      z_s   <- switch(.samp,
        sobol  = qnorm(randtoolbox::sobol( n = n_sim, dim = n_eta)),
        halton = qnorm(randtoolbox::halton(n = n_sim, dim = n_eta)),
        torus  = qnorm(randtoolbox::torus( n = n_sim, dim = n_eta)),
        lhs    = qnorm(.lhsSample(n_sim, n_eta)),
        rnorm  = matrix(rnorm(n_sim * n_eta), nrow = n_sim),
        qnorm(randtoolbox::sobol(n = n_sim, dim = n_eta))
      )
      eta_mat <- z_s %*% t(L)
      colnames(eta_mat) <- eta_nms
    } else {
      eta_mat <- matrix(0, nrow = n_sim, ncol = 0)
    }
    # ... and on the general path every simulated subject carries its own
    # covariate value. Without this the solve succeeds AT THE COVARIATE MEAN and
    # the residual is composed onto a var_f with the covariate spread removed,
    # while the OBSERVED V in the same panel still carries it: predicted
    # covariances measured 29-35% low and off-diagonals 61% low on an audited
    # model, i.e. structured standardised residuals for a fit that is fine.
    # .admStudyCovRows() draws nothing unless the path is "rows", so any study
    # routed some other way solves at the covariate MEAN here even though it
    # declares a distribution -- every panel for it then describes a model with
    # no covariate spread at all. The estimators are entitled to their own
    # reduction; this draw is not, because its etas are ordinary draws from
    # Omega rather than whatever the reduction re-aimed them onto. The general
    # per-row representation is always valid, just slower, and cost does not
    # matter for one diagnostic draw, so force it for anything with a cov_dist.
    if (!is.null(s[["cov_dist"]])) s$.adm_cov_path <- "rows"
    s <- tryCatch(.admStudyCovRows(s, pinfo_r, nrow(eta_mat)),
                  error = function(e) s)
    # One residual placeholder per observed output (rxerr.<output>); a
    # multi-endpoint solve needs every endpoint's rxerr present. rxSolve
    # defaults everything else (CMT, hard-coded constants).
    rxerr_nms <- { so <- unique(sig_output[!is.na(sig_output)])
                   if (length(so)) paste0("rxerr.", so) else "rxerr.cp" }
    col_nms   <- c(names(extra$struct), eta_nms, sig_nms, rxerr_nms)
    params_df <- as.data.frame(matrix(0, nrow = n_sim, ncol = length(col_nms),
                                      dimnames = list(NULL, col_nms)))
    params_df[, rxerr_nms] <- 1
    tryCatch(
      # A joint unit's outputs share one set of etas and are stacked into one
      # vector; .admSimulate() solves a single output, so it returned the first
      # endpoint's trajectory for every row. Same shared-eta solve the estimators
      # use, for the same unit.
      # sigdig: the tolerance the FIT solved at (extra$sigdig), not rxode2's
      # default. Without it the predicted-mean, residual and predicted-covariance
      # panels are computed from a different integration than the objective was
      # minimised on, and the standardised-residual panel can show structure the
      # fit never saw -- an artefact that moves with sigdig and disappears when
      # it is NULL. NULL for a fit made before this field existed, which is
      # exactly the previous behaviour.
      if (isTRUE(s$is_joint))
        .admSimulateJoint(rxMod, extra$struct, sig_nms, eta_mat, s, params_df, 1L,
                          sigdig = extra$sigdig)
      else
        .admSimulate(rxMod, extra$struct, sig_nms, eta_mat, s, ov, params_df, 1L,
                     sigdig = extra$sigdig),
      error = function(e) {
        if (warn) warning("plot.admFit: simulation failed: ", e$message, call. = FALSE)
        NULL
      })
  }

  # Returns BOTH the residual-adjusted covariance and mean: lnorm rescales the
  # mean, and the predicted E must carry that scaling just as the NLL does.
  .add_sigma <- function(V, mu, ov = out_var, times = NULL, phi = NULL,
                         cp = NULL) {
    # beta: the precision is SOLVED and rides back on the simulated matrix. Every
    # estimator patches it in; this path did not, so after a perfectly ordinary
    # beta fit the predicted-covariance heatmap, the standardised-residual panels
    # and the +-1 SD ribbon were all NA, silently.
    arr <- .admUnitResidRows(pinfo_r, ov, sv, length(mu), phi = phi)
    # Without `times` + the structural covariance the off-diagonal forms (ar,
    # ordinal) were dropped, so the predicted-covariance diagnostic panel showed
    # an independent-residual V for exactly the models whose off-diagonal is the
    # point of fitting them.
    # Match the objective: a TBS endpoint composes at each DRAW, so what a
    # diagnostic draws is what the fit was actually scored against. Needs the
    # simulated matrix, which is why `cp` is threaded in; without it the panel
    # falls back to the expansion, which is the previous behaviour.
    if (!is.null(cp)) {
      .ex <- .admResidNodeMomentsTBS(cp, rep(1, nrow(cp)), arr, times)
      if (!is.null(.ex)) return(list(V = .ex$V, mu = .ex$E))
    }
    ap  <- .admResidApply(mu, diag(V), arr, times, V)
    list(V = .admApplyResidTail(V, ap), mu = ap$mu)
  }

  setNames(lapply(names(studies), function(nm) {
    s      <- studies[[nm]]
    cp_mat <- .sim_study(s)
    if (is.null(cp_mat)) return(NULL)
    mu     <- colMeans(cp_mat)
    # A JOINT (same-subject, multi-output) unit stacks several endpoints into one
    # mean vector, so a single `output` cannot describe its rows: passing one made
    # .admResidRows() build the whole array from the FIRST endpoint's spec, and the
    # diagnostic panels then showed a covariance the fit never used (plasma's
    # prop() applied to the brain rows, and so on). Route it through
    # .admJointResidual() -- the estimators' own per-row-output path -- rather than
    # reconstructing the residual here for a second time.
    res    <- if (isTRUE(s$is_joint))
      .admJointResidual(mu, crossprod(sweep(cp_mat, 2L, mu)) / nrow(cp_mat),
                        s, pinfo_r, sv)
    else
      .add_sigma(crossprod(sweep(cp_mat, 2L, mu)) / nrow(cp_mat), mu,
                 s$output %||% out_var, s$times, attr(cp_mat, "phi"), cp_mat)
    V_pred <- res$V; mu <- res$mu
    obs_E  <- as.numeric(s$E)
    obs_V  <- as.matrix(s$V)
    # Joint row labels repeat the times across endpoints, so label them by the
    # endpoint they belong to; a plain unit keeps the bare times it always had.
    tnm    <- if (isTRUE(s$is_joint))
      paste0(.admRowOutput(s, length(mu)), "@", .admRowTimes(s, length(mu)))
    else as.character(s$times)
    names(mu) <- names(obs_E) <- tnm
    dimnames(V_pred) <- dimnames(obs_V) <- list(tnm, tnm)
    list(times = s$times, n = s$n,
         obs  = list(E = obs_E, V = obs_V),
         pred = list(E = mu,    V = V_pred))
  }), names(studies))
}

## Attach an `aggData` slot (observed + predicted aggregate moments per study) to
## a freshly constructed admFit. Shared by the admc/adfo/adgh/adirmc estimators.
## Computed at the fitted parameters with the fit's own `n_sim` and a fixed seed
## so `fit$env$aggData` matches the default `plot(fit)` mean/cov panels. A failure
## to simulate must not break fit construction, so the whole thing is guarded and
## a `NULL`/all-NULL result simply leaves the slot unset.
.admAttachAggData <- function(fit, extra, ui, seed = 1L) {
  ad <- tryCatch(.admAggData(extra, ui, n_sim = extra$n_sim, seed = seed, warn = FALSE),
                 error = function(e) NULL)
  if (!is.null(ad) && any(!vapply(ad, is.null, logical(1))))
    fit$env$aggData <- ad
  invisible(fit)
}

## `n` distinct colours. Okabe-Ito while it lasts -- colour-blind safe, and
## already the package's choice for the restart traces -- and a continuous HCL
## ramp past it.
##
## NOT recycled. Okabe-Ito has seven colours, and `rep_len` gave the eighth
## source the same black as the first: two entries in one legend, identically
## coloured, with nothing saying they were different studies. A per-band
## design with ten sources is an ordinary thing to plot.
.admOkabeIto <- function(n, black = TRUE) {
  ok <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#0072B2", "#D55E00", "#CC79A7")
  if (!black) ok <- ok[-1L]
  if (n <= length(ok)) ok[seq_len(max(n, 0L))]
  else grDevices::hcl.colors(n, "Dark 3")
}

## The covariate values a study is CONDITIONED at, as a title fragment.
##
## `stratify` splits one source into strata named `<source>_s1`, `_s2`, ... and
## nothing in that name says which slice is which, so a banded fit printed four
## panels per source with no way to tell the men from the women. The stratum
## itself knows: its `cov_dist` entry for a stratified covariate is a `.point`
## spec and `cov` carries the value it was conditioned at. `at`/`by` land in the
## same shape with no `cov_dist` entry at all.
##
## Covariates the study MARGINALISES over are deliberately left out -- `cov`
## holds a filled-in mean for those too (.admCheckCovariates), and printing it
## would label a panel with a value no simulated subject was held at.
## Returns "" when nothing is conditioned.
.admStudyCovLabel <- function(s) {
  cov <- s[["cov"]]
  if (is.null(cov) || !length(cov)) return("")
  # Through .admCovStudySpec, so a covariate DROPPED from the design still
  # counts as the marginalised one it is. It keeps a `cov` value and has no
  # `cov_dist` entry, which read directly would label it as conditioned -- the
  # title would assert "solved at CRCL = 90" for a study that integrated over a
  # distribution centred there.
  fix <- Filter(function(cv)
    .admCovStudyKind(s, cv) == "conditional", names(cov))
  if (!length(fix)) return("")
  paste(vapply(fix, function(cv)
    paste0(cv, " = ", format(cov[[cv]], digits = 3)), character(1)),
    collapse = ", ")
}

## `Study 'x'` for a pooled study, `Study 'x' [SEX = 0]` for a conditioned one.
.admStudyTitle <- function(s, nm) {
  lab <- .admStudyCovLabel(s)
  if (nzchar(lab)) sprintf("Study '%s' [%s]", nm, lab)
  else             sprintf("Study '%s'", nm)
}

## ---- covariate panels -------------------------------------------------------
##
## WHY THESE ARE SEPARATE PANELS RATHER THAN A COLUMN IN THE OTHERS. The mean
## and cov grids are per study, in isolation, and with aggregate data a
## covariate effect is identified BETWEEN studies: the renal exponent in the
## covariates vignette comes entirely from three cohorts sitting at three
## different CRCL medians, none of which fitted a renal term. A diagnostic that
## never puts two studies in the same frame structurally cannot show it -- three
## sources can each fit beautifully on their own panel while the covariate
## relation tying them together is wrong.

## What a study holds covariate `cv` at, as a quantile of its own declared
## distribution. A conditioned/stratified covariate is a `.point` spec, so its
## value comes off `cov` and does not move with `u`. NA when the study says
## nothing about `cv`.
## The distribution a study declared for `cv`, whether or not the design kept
## it. A covariate the model never reads is dropped from `cov_dist` (it cannot
## move the prediction, so integrating over it is waste) but the source still
## described it, and that description is what the residual panel plots against.
.admCovStudySpec <- function(s, cv)
  s[["cov_dist"]][[cv]] %||% s[[".adm_cov_dropped"]][[cv]]

.admCovStudyQ <- function(s, cv, u = 0.5) {
  sp <- .admCovStudySpec(s, cv)
  if (!is.null(sp) && !isTRUE(sp[[".point"]]))
    return(tryCatch(as.numeric(.admCovQuantile(sp, u))[1L],
                    error = function(e) NA_real_))
  v <- s[["cov"]][[cv]]
  if (is.null(v)) NA_real_ else as.numeric(v)[1L]
}

## Where a study SITS on the covariate's axis -- the one number that stands for
## it in both panels.
##
## The median, except for a covariate marginalised over declared LEVELS, where
## it is the probability-weighted mean. A binary margin's median is a step
## function of its probabilities: `probs = c(0.45, 0.55)` and `c(0.55, 0.45)`
## are all but the same cohort, and .admCovQuantile() reads their medians as 1
## and 0 -- the full width of the axis apart. Worse, every study with an even
## split lands on the SAME level, so a covariate the sources genuinely differ
## on shows no between-study contrast at all and its residual facet is dropped.
## The mean moves with the mixture, which is what a marginalised source IS.
##
## Conditioned studies are untouched: a point spec has no mixture to average,
## and its centre is the value it was solved at.
.admCovSpecCentre <- function(sp) {
  if (is.null(sp)) return(NA_real_)
  if (!isTRUE(sp[[".point"]]) && !is.null(sp$values)) {
    v  <- as.numeric(sp$values)
    pr <- as.numeric(sp$probs %||% rep(1 / length(v), length(v)))
    if (length(pr) == length(v) && sum(pr) > 0) return(sum(v * pr) / sum(pr))
  }
  tryCatch(as.numeric(.admCovQuantile(sp, 0.5))[1L],
           error = function(e) NA_real_)
}

.admCovStudyCentre <- function(s, cv) {
  sp <- .admCovStudySpec(s, cv)
  if (!is.null(sp) && !isTRUE(sp[[".point"]])) {
    v <- .admCovSpecCentre(sp)
    if (is.finite(v)) return(v)
  }
  .admCovStudyQ(s, cv, 0.5)
}

## The declared levels of `cv` a study actually covers.
##
## A level is covered when the study conditions at it, or when its declared
## margin puts non-zero probability there. This is the SUPPORT, and it is what
## decides whether a level on a discrete axis is extrapolation -- not whether
## any study's centre happens to land on it, which for an even binary split is
## no level at all.
.admCovStudySupport <- function(s, cv, levels) {
  sp <- .admCovStudySpec(s, cv)
  if (!is.null(sp) && !isTRUE(sp[[".point"]]) && !is.null(sp$values)) {
    v  <- as.numeric(sp$values)
    pr <- as.numeric(sp$probs %||% rep(1 / length(v), length(v)))
    if (length(pr) != length(v)) pr <- rep(1 / length(v), length(v))
    return(levels[vapply(levels, function(l)
      any(abs(v - l) < 1e-8 & pr > 0), logical(1))])
  }
  # A continuous margin marginalised onto a level axis, or a point: covered
  # where its 2.5th-97.5th reaches.
  lo <- .admCovStudyQ(s, cv, 0.025); hi <- .admCovStudyQ(s, cv, 0.975)
  if (!is.finite(lo) || !is.finite(hi)) return(numeric(0))
  levels[levels >= lo - 1e-8 & levels <= hi + 1e-8]
}

## How a study enters covariate `cv`: CONDITIONAL or MARGINAL?
##
## This is a property of the SOURCE'S OWN published model, and the mechanism
## follows from it rather than the other way round.
##
## CONDITIONAL means that analyst ESTIMATED this effect. Their result can
## therefore be read at one value of the covariate -- by subgroup (`at`, `by`),
## or by banding the source into strata (`stratify`) -- and the study carries a
## point spec. It contributes a point of contrast against the other sources.
##
## MARGINAL means they did not. No single paper has a contrast to report, so
## splitting it on that covariate would MANUFACTURE one; admixr2 integrates
## over the population the paper enrolled instead, and the study carries a
## distribution. It constrains the effect through the mixture that induces.
##
## Both feed the same fit, and they are different kinds of evidence, so they
## have to look different. Drawing a conditioned source as a centre and a bar
## would show it as a degenerate version of a marginal one.
.admCovStudyKind <- function(s, cv) {
  sp <- .admCovStudySpec(s, cv)
  if (is.null(sp) || isTRUE(sp[[".point"]])) "conditional" else "marginal"
}

## The levels of a DISCRETE covariate, or NULL when it is continuous.
##
## Discrete means declared levels (`values`), or a covariate every study
## conditions at a point -- `stratify`, `at`, `by`. Both panels have to agree
## about this: one drawing `SEX` on a swept axis while the other treats it as
## two levels is the same kind of split that already mislabelled a dropped
## covariate once.
##
## `mid` is each study's median, so a handful of distinct conditioned values
## are read as levels. The cap stops a continuous covariate that happens to be
## conditioned in every study from becoming a 30-level factor.
.admCovLevels <- function(cv, studies, mid, max_lev = 8L) {
  lev <- unlist(lapply(studies, function(s) .admCovStudySpec(s, cv)$values))
  if (length(lev)) return(sort(unique(as.numeric(lev))))
  if (!all(vapply(studies, .admCovStudyKind, character(1),
                  cv = cv) == "conditional")) return(NULL)
  u <- sort(unique(mid[is.finite(mid)]))
  if (length(u) <= max_lev) u else NULL
}

## One colour per source, held by NAME across both covariate panels.
##
## The two panels are handed different name sets -- the effect panel draws one
## mark per SOURCE while the residual panel keeps the strata apart, so it sees
## `x` where the other sees `x_s1` and `x_s2` -- and an unnamed palette hands
## the same source different positions in the vector, so it comes out orange on
## one panel and blue on the other. Built once from the union of both.
## BLACK IS RESERVED for the estimated effect, which both panels draw in it.
## A source in black could not be told from the fit it is being compared
## against -- and on the level axis, where both are a line joining two points,
## they were indistinguishable.
.admCovPalette <- function(nms) {
  nms <- sort(unique(nms[!is.na(nms)]))
  stats::setNames(.admOkabeIto(length(nms), black = FALSE), nms)
}

## Point AREA is the study's sample size, on both covariate panels.
##
## Area rather than radius, via scale_size_area(), so the encoding is the one a
## reader actually decodes and zero maps to zero. This follows multinma's
## `weight_nodes`, which scales its network nodes by sample size.
##
## WITH a legend. The residual panel had `guide = "none"`, so area carried n and
## nothing on the figure said so -- a reader could see that one mark was bigger
## and had no way to learn what bigger meant. A single study, or a set that all
## report the same n, gets no legend: there is nothing to compare.
.admCovSizeScale <- function(n, max_size = 5.5) {
  n  <- n[is.finite(n)]
  br <- if (length(unique(n)) > 1L)
    unique(round(range(n))) else ggplot2::waiver()
  if (!length(n) || length(unique(n)) < 2L)
    ggplot2::scale_size_area(max_size = max_size * 0.75, guide = "none")
  else
    ggplot2::scale_size_area(
      max_size = max_size, breaks = br, name = "n",
      guide = ggplot2::guide_legend(order = 3L, override.aes =
                                      list(colour = "grey40", shape = 16L)))
}

## Scales and styling shared by the two covariate panels.
##
## Both key colour on the source and shape on marginal/conditioned, and both
## want the small grey subtitle. Held once so a styling change cannot land on
## one panel and not the other -- the legends are meant to be read across the
## pair, so a source that is orange in one and blue in the other is worse than
## either choice alone.
##
## The guides carry an explicit `order`. Left at the default ggplot2 sorts
## equal-priority guides by a hash of their contents, which is not stable across
## sessions -- the source legend and the shape legend swapped places between two
## runs of the same code, which is both a confusing figure and a snapshot test
## that fails at random. Source first, since it is the one a reader looks up.
.admCovPanelStyle <- function(pal, kinds)
  list(
    ggplot2::scale_colour_manual(
      values = pal, name = NULL,
      guide = ggplot2::guide_legend(order = 1L)),
    ggplot2::scale_shape_manual(
      values = c(marginal = 16L, conditional = 18L), name = NULL,
      breaks = intersect(c("marginal", "conditional"), kinds),
      guide = ggplot2::guide_legend(order = 2L)),
    ggplot2::theme_bw(),
    ggplot2::theme(plot.subtitle = ggplot2::element_text(
      size = 7, colour = "grey40", face = "plain")))

## An axis-breaks function that ticks a discrete facet at its LEVELS only.
##
## `scales = "free_x"` gives each facet its own scale but they share one breaks
## function, so it has to decide from the panel limits alone: a panel whose span
## is filled by known levels is the discrete one. Everything else falls back to
## base `pretty()`, which is what an untouched continuous axis would use.
##
## Without this a binary covariate is ticked at 0.25 and 0.75 -- values it does
## not have, and that the model was never asked about.
.admLevelBreaks <- function(levels) {
  lv <- sort(unique(levels))
  function(lims) {
    inside <- lv[lv >= lims[1L] & lv <= lims[2L]]
    # 0.8, not 0.95: the levels sit inside ggplot2's default 5%-a-side
    # expansion, and a shaded end level pushes the limits out further still.
    if (length(inside) >= 2L && diff(range(inside)) >= 0.8 * diff(lims))
      inside
    else pretty(lims)
  }
}

## The source a study came from: `stratify` splits one into `<source>_s1`,
## `_s2`, ... and the strata of one source are a PAIRED set, not independent
## points. That pairing is the evidence stratification creates, so the panels
## need to be able to recover it.
.admCovSource <- function(nm) sub("_s[0-9]+$", "", nm)

## Covariates worth a facet: any covariate any study describes.
##
## NOT restricted to the ones the model reads. A covariate the model omits is
## dropped from the design -- correctly, it cannot move the prediction -- and
## that is precisely the case the residual panel exists for: the analyst left a
## term out and wants to know whether it belonged. `.admCovEffectData()` returns
## NULL for such a covariate on its own, since there is no fitted effect to
## draw, so only the residual panel picks it up. A covariate no study describes
## at all cannot be placed on an axis and is excluded here.
.admCovPanelCovs <- function(ui, studies) {
  covs <- unique(c(tryCatch(ui$allCovs, error = function(e) character(0)),
                   unlist(lapply(studies, function(s)
                     c(.admCovSpecNames(s[["cov_dist"]]),
                       names(s[[".adm_cov_dropped"]] %||% list()),
                       names(s[["cov"]] %||% list()))), use.names = FALSE)))
  if (!length(covs)) return(character(0))
  covs[vapply(covs, function(cv) any(vapply(studies, function(s)
    !is.null(s[["cov"]][[cv]]) || !is.null(.admCovStudySpec(s, cv)),
    logical(1))), logical(1))]
}

## The n-weighted median across studies, per covariate -- where the OTHER
## covariates are held while one of them is swept. The alternative, each
## study's own values, puts the per-study marks off the curve for reasons that
## have nothing to do with the covariate being swept.
.admCovPooled <- function(covs, studies) {
  w <- vapply(studies, function(s) {
    n <- suppressWarnings(as.numeric(s[["n"]] %||% 1))
    if (length(n) != 1L || !is.finite(n) || n <= 0) 1 else n
  }, double(1))
  setNames(lapply(covs, function(cv) {
    v  <- vapply(studies, .admCovStudyCentre, double(1), cv = cv)
    ok <- is.finite(v)
    if (!any(ok)) 1 else sum(v[ok] * w[ok]) / sum(w[ok])
  }), covs)
}

## One panel per (parameter, covariate) the model ESTIMATES an effect for.
##
## The question the panel answers is "does the fitted covariate effect agree
## with the sources it was fitted from". So it carries exactly three things:
##
##   * the ESTIMATED EFFECT, as a dotted line across the covariate axis -- the
##     MBMA model's parameter at the fitted thetas.
##   * each SOURCE at its own parameter value, from that source's own published
##     model. This is the comparison: a source sitting off the dotted line is
##     one the meta-analysis does not reproduce.
##   * how far along the axis each source SPEAKS FOR, and that is where
##     conditional and marginal differ. A conditional source was banded or
##     reported at a value, so it covers a stated range: a solid line. A
##     marginal one reported no contrast and admixr2 integrates over the
##     population it enrolled: a whisker, because that is a distribution.
##
## A covariate whose coefficient the model does not estimate gets no panel.
## There is no fitted effect to agree or disagree with, and a facet drawn for
## one -- a fixed 0.75 allometric exponent, say -- invites a reader to check an
## agreement that was never in question.
##
## Returns NULL when the model estimates nothing for `cv`, when no line reads
## it, when no study gives it a finite value, or when nothing varies.
.admCovEffectData <- function(ui, cv, studies, struct, src = NULL,
                              n_grid = 120L, pad = 0.15) {
  ml <- .admModelLines(ui)
  if (is.null(ml)) return(NULL)
  hit <- .admLinesReading(ml, cv)
  if (!length(hit)) return(NULL)
  est <- tryCatch(.admCovCoefThetas(ui, cv, NULL), error = function(e) NULL)
  if (is.null(est) || !length(est)) return(NULL)

  mid <- vapply(studies, .admCovStudyCentre, double(1), cv = cv)
  lo  <- vapply(studies, .admCovStudyQ, double(1), cv = cv, u = 0.1)
  hi  <- vapply(studies, .admCovStudyQ, double(1), cv = cv, u = 0.9)
  # The tails as well as the body: a marginalised source is a DISTRIBUTION the
  # estimator integrates over, and 10th-90th alone draws it as if it stopped
  # there. The 2.5th-97.5th also decides where extrapolation starts.
  lo2 <- vapply(studies, .admCovStudyQ, double(1), cv = cv, u = 0.025)
  hi2 <- vapply(studies, .admCovStudyQ, double(1), cv = cv, u = 0.975)
  knd <- vapply(studies, .admCovStudyKind, character(1), cv = cv)
  nn  <- vapply(studies, function(z) {
    v <- suppressWarnings(as.numeric(z[["n"]] %||% NA_real_)[1L])
    if (length(v) != 1L || !is.finite(v) || v <= 0) NA_real_ else v
  }, double(1))
  if (!any(is.finite(mid))) return(NULL)

  # DISCRETE: declared levels, or a covariate every study conditions at a
  # point. Swept continuously it draws the model at SEX = 0.37, which is not a
  # patient and not a prediction anyone can act on; its axis is the levels.
  disc    <- .admCovLevels(cv, studies, mid)
  is_disc <- !is.null(disc) && length(disc) > 1L

  covered <- range(c(lo2, lo, mid, hi, hi2), na.rm = TRUE)
  span    <- diff(covered)
  # Every study at one value: no covariate effect to show, and a curve through
  # territory no source spoke to would be invention.
  if (!is.finite(span) || span <= 0) return(NULL)

  grid <- if (is_disc) disc
          else seq(covered[1L] - pad * span, covered[2L] + pad * span,
                   length.out = n_grid)

  # THE ESTIMATED EFFECT. Other covariates at the pooled centre: this is one
  # line, not one per level of something else. A panel that splits is a panel
  # answering a second question, and the effect is the slope of this one.
  base <- .admCovPooled(.admCovPanelCovs(ui, studies), studies)
  at   <- utils::modifyList(base, stats::setNames(list(grid), cv))
  vals <- .admEvalModelLines(ml, at, struct, keep = hit)
  vals <- Filter(function(z) length(z$value) == length(grid) &&
                   all(is.finite(z$value)) && diff(range(z$value)) > 0, vals)
  if (!length(vals)) return(NULL)

  curve <- do.call(rbind, lapply(vals, function(z)
    data.frame(cov = cv, param = z$name, x = grid, y = as.numeric(z$value),
               disc = is_disc, stringsAsFactors = FALSE)))
  params <- unique(curve$param)

  # EACH SOURCE AT ITS OWN VALUE, once per position it speaks for ON THIS AXIS.
  #
  # Per SOURCE, not per stratum. A source banded on SEX expands into two
  # studies, and evaluating each of them drew that source TWICE on the CRCL
  # panel -- two whiskers at the same renal value, differing only in a
  # covariate this facet is not about. Every other covariate is held at the
  # source's OWN centre instead, which collapses the strata back to one mark.
  #
  # On the axis the source IS conditional on, the positions are the values it
  # reported, so a source banded on SEX gets one mark at each level -- and the
  # gap between them is that paper's own effect, to read against the gap
  # between the fitted points.
  #
  # A source carrying no model of its own gets no mark. Falling back to the
  # fitted curve would put it exactly on the dotted line and read as agreement
  # with a claim it never made.
  by_src <- split(seq_along(studies), .admCovSource(names(studies)))
  .or <- function(v, i) if (is.finite(v[i])) v[i] else mid[i]

  mk <- do.call(rbind, lapply(names(by_src), function(sn) {
    ii <- by_src[[sn]]
    ii <- ii[is.finite(mid[ii])]
    if (!length(ii)) return(NULL)
    so <- src[[sn]]
    if (is.null(so) || is.null(so[["ui"]])) return(NULL)
    sml <- .admModelLines(so[["ui"]])
    if (is.null(sml)) return(NULL)
    at_own <- .admCovSourceAt(so, so[["ui"]])

    cond <- any(knd[ii] == "conditional")
    # Conditional: the distinct values it reported. Marginal: one position,
    # its centre -- the strata share the distribution this facet is about.
    pos <- if (cond) sort(unique(mid[ii])) else mid[ii][1L]
    ref <- if (cond) ii[match(pos, mid[ii])] else ii[1L]
    # `n` for a POSITION, not for a stratum. `stratify` divides a source's n
    # among its strata, so a source banded on sex contributes half its patients
    # at each sex level -- but on every other covariate's axis both strata land
    # on the same position and the mark speaks for the whole source. Summing per
    # position gets both right without a special case.
    pos_n <- vapply(pos, function(v)
      sum(nn[ii][abs(mid[ii] - v) < 1e-8], na.rm = TRUE), double(1))
    pos_n[!is.finite(pos_n) | pos_n <= 0] <- NA_real_


    do.call(rbind, lapply(seq_along(pos), function(k) {
      a <- utils::modifyList(at_own, stats::setNames(list(pos[k]), cv))
      v <- tryCatch(.admEvalModelLines(sml, a), error = function(e) list())
      v <- Filter(function(z) length(z$value) == 1L && is.finite(z$value), v)
      if (!length(v)) return(NULL)
      o <- stats::setNames(lapply(v, `[[`, "value"),
                           vapply(v, `[[`, "", "name"))
      i0 <- ref[k]
      do.call(rbind, lapply(intersect(params, names(o)), function(pp)
        data.frame(
          cov   = cv,
          param = pp,
          study = sn,
          kind  = if (cond) "conditional" else "marginal",
          x     = pos[k],
          xlo   = .or(lo,  i0), xhi  = .or(hi,  i0),
          xlo2  = .or(lo2, i0), xhi2 = .or(hi2, i0),
          y     = as.numeric(o[[pp]]),
          n     = pos_n[k],
          stringsAsFactors = FALSE)))
    }))
  }))

  # THE SOURCE'S OWN REGRESSION, over the range it covers.
  #
  # Only for a source CONDITIONAL on this covariate, because that is the source
  # that reported a relationship here: its model estimated the effect, which is
  # what let it be banded or read at a value in the first place. A marginal
  # source reported no contrast along this axis and has no line of its own --
  # its whisker says what it covered, and nothing about slope.
  #
  # This is the comparison the panel exists for. The dotted line is the
  # meta-analysis; a source's own line running at a different slope over a
  # range that source actually enrolled is a paper the fit does not reproduce,
  # and the two being parallel but offset is a different finding from the two
  # crossing.
  slines <- do.call(rbind, Filter(Negate(is.null),
    lapply(names(by_src), function(sn) {
      ii <- by_src[[sn]]
      ii <- ii[is.finite(mid[ii])]
      if (!length(ii) || !any(knd[ii] == "conditional")) return(NULL)
      so <- src[[sn]]
      if (is.null(so) || is.null(so[["ui"]])) return(NULL)
      sml <- .admModelLines(so[["ui"]])
      if (is.null(sml) || !length(.admLinesReading(sml, cv))) return(NULL)
      # On a level axis the source's own line joins the levels it reported, the
      # same way the dotted estimated effect does. On a continuous one it needs
      # the range it covers, and a source that declared none has no extent to
      # draw over -- a point value gets a diamond and no line.
      g <- if (is_disc) sort(unique(mid[ii])) else {
        rg <- .admCovSourceRange(so, cv)
        if (is.null(rg)) return(NULL)
        seq(rg[1L], rg[2L], length.out = 40L)
      }
      if (length(g) < 2L) return(NULL)
      a <- utils::modifyList(.admCovSourceAt(so, so[["ui"]]),
                             stats::setNames(list(g), cv))
      v <- tryCatch(.admEvalModelLines(sml, a), error = function(e) list())
      v <- Filter(function(z) length(z$value) == length(g) &&
                    all(is.finite(z$value)), v)
      if (!length(v)) return(NULL)
      do.call(rbind, lapply(Filter(function(z) z$name %in% params, v),
        function(z) data.frame(
          cov = cv, param = z$name, study = sn, x = g,
          y = as.numeric(z$value), stringsAsFactors = FALSE)))
    })))

  # EXTRAPOLATION, the same idea on both kinds of axis: grey marks where the
  # fit is speaking past its sources. On a continuous axis that is the padding
  # beyond the range they cover; on a discrete one it is a declared LEVEL no
  # study sits at. From the studies' SUPPORT, not their centres: a source
  # marginalising over an even binary split has a centre at 0.5 and sits on no
  # level, and reading centres would grey a level every source sampled.
  shade <- if (!is_disc) data.frame(
      cov  = cv,
      xmin = c(min(grid), covered[2L]),
      xmax = c(covered[1L], max(grid)),
      stringsAsFactors = FALSE)
    else {
      sup <- unique(unlist(lapply(studies, .admCovStudySupport,
                                  cv = cv, levels = grid)))
      un  <- setdiff(grid, sup)
      if (!length(un) || length(grid) < 2L) NULL else {
        # Clamped to the level range: a rect hanging past the outermost level
        # widens the limits enough that .admLevelBreaks() falls back to
        # pretty() and ticks a three-level factor at 0.5 and 1.5.
        w <- 0.4 * min(diff(sort(grid)))
        data.frame(cov = cv,
                   xmin = pmax(un - w, min(grid)),
                   xmax = pmin(un + w, max(grid)),
                   stringsAsFactors = FALSE)
      }
    }

  list(curve = curve, marks = mk, shade = shade, slines = slines)
}


## Between-study mean standardised residual against a covariate.
##
## One point per study: its mean z over observation times against the covariate
## value it sits at. A covariate form that is wrong -- an exponent where the
## truth is a threshold, a linear term where the truth is allometric -- shows up
## here as a TREND across studies, which is the thing no single study's panel
## can contain.
##
## Returns NULL when fewer than two studies carry a finite value, since a single
## point has no between-study contrast to read.
.admCovResidData <- function(cv, studies, agg) {
  df <- do.call(rbind, lapply(names(studies), function(nm) {
    s  <- studies[[nm]]
    ag <- agg[[nm]]
    if (is.null(ag)) return(NULL)
    x <- .admCovStudyCentre(s, cv)
    if (!is.finite(x)) return(NULL)
    n  <- as.numeric(s[["n"]] %||% NA_real_)
    se <- sqrt(diag(ag$pred$V) / n)
    z  <- (as.numeric(ag$obs$E) - as.numeric(ag$pred$E)) / se
    z  <- z[is.finite(z)]
    if (!length(z)) return(NULL)
    # The covariate range each study speaks for, so a point that sits at 62 on
    # the axis is not read as a study that only ever saw 62. A conditioned
    # study genuinely did see one value, and gets a zero-width span.
    xl <- .admCovStudyQ(s, cv, 0.1); xh <- .admCovStudyQ(s, cv, 0.9)
    data.frame(cov = cv, study = nm, source = .admCovSource(nm),
               kind = .admCovStudyKind(s, cv),
               x = x, xlo = if (is.finite(xl)) xl else x,
               xhi = if (is.finite(xh)) xh else x,
               z = mean(z), n = n,
               label = .admStudyCovLabel(s), stringsAsFactors = FALSE)
  }))
  if (is.null(df) || nrow(df) < 2L) return(NULL)
  # A covariate with no between-study contrast has nothing for this panel to
  # read, and the test is against the WITHIN-study spread rather than against
  # zero. Three cohorts drawing weight from the same distribution have medians
  # a few hundred grams apart -- sampling noise in the fitted margins, not
  # evidence -- and an exact-equality guard let that draw a facet, blow it up to
  # full panel width on a free x scale, and fit a trend through it.
  #
  # A TENTH of the typical 10th-90th is the bar. Two cohorts at 70 kg and 90 kg
  # with a 20% CV overlap heavily and still carry a real contrast, so the bar
  # cannot be set near the spread itself; what this has to catch is sources that
  # declared the SAME distribution and differ only in what a finite sample of it
  # estimated, which lands two orders of magnitude below.
  w <- stats::median(df$xhi - df$xlo, na.rm = TRUE)
  if (!is.finite(w)) w <- 0
  if (diff(range(df$x)) <= max(0, 0.1 * w)) return(NULL)
  # A discrete covariate is read as a CONTRAST, not a trend, and the two are
  # drawn differently: see the panel code.
  # Over ALL studies, not the ones that reached this data frame. A study with
  # no `aggData` entry drops out above, and deciding discreteness from what
  # survived lets this panel call a covariate discrete while the effect panel,
  # which sees every study, calls it continuous -- different axis ticks on the
  # two halves of one figure.
  df$disc <- !is.null(.admCovLevels(
    cv, studies, vapply(studies, .admCovStudyCentre, double(1), cv = cv)))
  # A source whose strata both survived is a PAIR -- the within-source contrast
  # stratifying on this covariate produced. A source appearing once has no pair
  # and its line would be a dot.
  df$paired <- df$source %in% names(which(table(df$source) > 1L))
  df
}

## ---- what each SOURCE itself claims -----------------------------------------
##
## Two study lists live on a fit and they are not the same object.
## `admExtra$studies` is the EXPANDED one the estimator integrates: strata
## already cut, stripped to `E`, `V`, `n`, `cov`, `cov_dist`. Every other panel
## reads it, correctly -- it is what was fitted.
##
## The estimator's control keeps the list the user actually passed, untouched:
## each source's OWN model at its own published values, the population it
## enrolled, and the `range` it declared. None of that survives expansion, and
## all of it is what a source's own covariate claim is made of.
.admFitSourceStudies <- function(fit) {
  e <- fit[["env"]]
  if (is.null(e)) return(NULL)
  # One name per estimator, as nmObjHandleControlObject.* assigns them.
  ctl <- e$adghControl %||% e$adfoControl %||% e$admControl %||% e$adirmcControl
  st  <- ctl$studies
  if (is.null(st) || !length(st)) NULL else unclass(st)
}

## The covariate range a source speaks for, most authoritative first: the range
## it DECLARED (`range =`, which is the enrolled range `stratify` truncates
## to), the range it actually enrolled if a population table was given, and
## failing both the body of its declared distribution.
##
## NULL when the source says nothing, which is a source with no line to draw
## rather than one to draw over a guessed range.
.admCovSourceRange <- function(s, cv) {
  r <- s[["range"]][[cv]]
  if (length(r) == 2L && all(is.finite(r)) && r[1L] < r[2L])
    return(sort(as.numeric(r)))
  p <- s[["population"]]
  # A raw baseline table: the range it literally enrolled.
  if (is.data.frame(p) && is.numeric(p[[cv]])) {
    rr <- range(p[[cv]], na.rm = TRUE)
    if (all(is.finite(rr)) && rr[1L] < rr[2L]) return(rr)
  }
  # Otherwise the declared margin. `population` is normally a covDist by the
  # time it reaches here -- admStudy() fits margins to the table and keeps
  # those, not the rows -- so it holds specs under the covariate names, exactly
  # like `cov_dist`. Reading only the data frame found a range for nobody.
  sp <- s[["cov_dist"]][[cv]] %||%
    (if (!is.null(p) && !is.data.frame(p)) p[[cv]])
  if (!is.null(sp)) {
    rr <- tryCatch(c(.admCovQuantile(sp, 0.025), .admCovQuantile(sp, 0.975)),
                   error = function(e) NULL)
    if (length(rr) == 2L && all(is.finite(rr)) && rr[1L] < rr[2L])
      return(as.numeric(rr))
  }
  NULL
}

## Where a source holds its OTHER covariates while one is swept: at its OWN
## patients' values, not the pooled ones. A paper's claim about weight is a
## claim about the people it enrolled, and holding them at the pooled centre
## would redraw that claim for somebody else's cohort.
.admCovSourceAt <- function(s, ui) {
  base <- .admCovNominal(ui, s[["cov_dist"]] %||% s[["population"]])
  p <- s[["population"]]
  if (is.data.frame(p)) {
    for (nm in intersect(names(base), names(p)))
      if (is.numeric(p[[nm]])) {
        m <- stats::median(p[[nm]], na.rm = TRUE)
        if (is.finite(m)) base[[nm]] <- m
      }
  } else {
    # THE SAME CENTRE the estimated effect holds its other covariates at.
    # .admCovNominal() takes the median, and the median of an even binary split
    # is the upper level -- so the dotted line sat at the pooled mean while
    # every source sat at SEX = 1, and each one looked displaced from a fit it
    # actually agreed with.
    for (nm in intersect(names(base), names(p))) {
      v <- .admCovSpecCentre(p[[nm]])
      if (is.finite(v)) base[[nm]] <- v
    }
  }
  # PINNED values last, and they are not optional. `at`/`by` take a covariate
  # out of `population` -- a population cannot also give a pinned covariate a
  # distribution -- so .admCovNominal() finds nothing for it and falls back to
  # 1. A source reported at CRCL = 62 then had its own model evaluated at
  # CRCL = 1, which put it off the bottom of every other panel.
  for (k in c("at", "by")) {
    v <- s[[k]]
    if (is.list(v))
      for (nm in intersect(names(base), names(v))) {
        x <- suppressWarnings(as.numeric(v[[nm]])[1L])
        if (is.finite(x)) base[[nm]] <- x
      }
  }
  base
}

## ---- covariate panel BUILDERS ----------------------------------------------
##
## Data frames in, a ggplot out. These were 200 lines inline in plot.admFit(),
## which meant the only way to reach them was to run a fit -- so every bug in
## them was a bug the test suite structurally could not see, and both rounds of
## review found rendering defects the suite had passed through. Out here they
## take synthetic study lists through .admCovEffectData()/.admCovResidData()
## and render in milliseconds, which is what makes a snapshot test possible.
##
## `pal` is shared and comes from .admCovPalette(): the two panels are read as
## a pair, so a source has to keep one colour across both.

## The fitted-effect panel. `eff` is the list .admCovEffectData() returns, one
## element per covariate. NULL when no covariate produced a curve.
.admCovEffectPanel <- function(eff, pal, src = NULL) {
  if (!length(eff)) return(NULL)
  curve_df <- do.call(rbind, lapply(eff, `[[`, "curve"))
  marks_df <- do.call(rbind, lapply(eff, `[[`, "marks"))
  shade_df <- do.call(rbind, lapply(eff, `[[`, "shade"))
  sline_df <- do.call(rbind, lapply(eff, `[[`, "slines"))
  if (is.null(curve_df) || !nrow(curve_df)) return(NULL)
  # Degenerate padding (a covariate covering the whole grid) would draw a
  # zero-width rect; harmless, but it puts a stray border on the panel.
  if (!is.null(shade_df))
    shade_df <- shade_df[shade_df$xmax > shade_df$xmin, , drop = FALSE]
  # Carry the shading onto the (cov, param) pairs that EXIST. Faceting on
  # ~cov+param over a layer keyed by `cov` alone makes ggplot2 conjure a panel
  # for every combination, so a parameter that never reads CRCL still got an
  # empty "CRCL, v" panel containing nothing but the grey rect.
  if (!is.null(shade_df) && nrow(shade_df))
    shade_df <- merge(shade_df, unique(curve_df[, c("cov", "param")]),
                      by = "cov")

  p_eff <- ggplot2::ggplot(curve_df, ggplot2::aes(x = x, y = y))
  if (!is.null(shade_df) && nrow(shade_df))
    p_eff <- p_eff + ggplot2::geom_rect(
      data = shade_df, inherit.aes = FALSE,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      fill = "grey85", alpha = 0.55)

  # THE ESTIMATED EFFECT: one dotted line, heavy enough to read as the
  # reference the sources are being compared against. On a level axis the line
  # joins the levels and they are marked with points -- the points are where
  # the model was actually asked, and the line is there so the effect reads the
  # same way on both kinds of axis. Dotted, and only between marked levels, so
  # it is a connector and not a claim about the space between them.
  p_eff <- p_eff + ggplot2::geom_line(
    data = curve_df, ggplot2::aes(group = param),
    colour = "black", linewidth = 1.1, linetype = "dotted")
  # An OPEN SQUARE, which is deliberately not in the source vocabulary. The
  # shape legend reads "circle = marginal, diamond = conditional" and those are
  # statements about a SOURCE; the fit is not a source, and marking its levels
  # with a filled black circle had it reading as a marginalised study.
  if (any(curve_df$disc))
    p_eff <- p_eff + ggplot2::geom_point(
      data = curve_df[curve_df$disc, , drop = FALSE],
      colour = "black", fill = "white", shape = 22L, stroke = 0.9,
      size = 3.2)

  # EACH SOURCE at its own published parameter value, over the stretch of the
  # axis it speaks for. CONDITIONAL: a solid line over the range it was banded
  # or reported across -- a stated extent. MARGINAL: a whisker, because what it
  # reported is a distribution, and drawing that as a line would claim the
  # source covers its tails as evenly as its middle.
  # A CONDITIONAL source's own regression, solid, over the range it covers.
  if (!is.null(sline_df) && nrow(sline_df))
    p_eff <- p_eff + ggplot2::geom_line(
      data = sline_df, inherit.aes = FALSE,
      ggplot2::aes(x = x, y = y, colour = study,
                   group = paste(study, param)),
      linewidth = 1.1)
  if (!is.null(marks_df) && nrow(marks_df)) {
    marg <- marks_df[marks_df$kind == "marginal", , drop = FALSE]
    if (nrow(marg))
      p_eff <- p_eff +
        ggplot2::geom_segment(
          data = marg, inherit.aes = FALSE,
          ggplot2::aes(x = xlo2, xend = xhi2, y = y, yend = y, colour = study),
          linewidth = 0.5, alpha = 0.7) +
        ggplot2::geom_segment(
          data = marg, inherit.aes = FALSE,
          ggplot2::aes(x = xlo, xend = xhi, y = y, yend = y, colour = study),
          linewidth = 1.6, alpha = 0.55)
    p_eff <- p_eff +
      ggplot2::geom_point(data = marks_df,
                          ggplot2::aes(colour = study, shape = kind,
                                       size = n)) +
      .admCovSizeScale(marks_df$n) +
      .admCovPanelStyle(pal, marks_df$kind)
  }

  # facet_wrap, not facet_grid: most parameters read only some covariates, and
  # a grid spends half the figure on empty (v, CRCL)-style panels.
  p_eff <- p_eff +
    ggplot2::scale_x_continuous(
      breaks = .admLevelBreaks(curve_df$x[curve_df$disc])) +
    # The strip says which parameter is on y against which covariate on x, so
    # the axis titles do not have to and the panel needs no key to read it.
    ggplot2::facet_wrap(~ cov + param, scales = "free",
                        labeller = function(d)
                          list(paste0(d$param, "  vs  ", d$cov))) +
    ggplot2::labs(
      title = "Estimated covariate effect against its sources",
      x = "Covariate value", y = "Parameter value",
      subtitle = paste("black dotted, open squares: the ESTIMATED effect",
                       " |  coloured: each source's own published model",
                       "\nCONDITIONAL: its own regression over the range it",
                       "covers  |  MARGINAL: a whisker, 10th-90th over",
                       "2.5th-97.5th"))
  # Styling comes from .admCovPanelStyle(), added with the marks above. A fit
  # with no marks at all still needs it.
  if (is.null(marks_df) || !nrow(marks_df))
    p_eff <- p_eff + .admCovPanelStyle(pal, character(0))
  p_eff
}


## The between-study residual panel. `res` is the list .admCovResidData()
## returns, one element per covariate. NULL when no covariate produced rows.
.admCovResidPanel <- function(res, pal) {
  if (!length(res)) return(NULL)
  res_df <- do.call(rbind, res)
  if (is.null(res_df) || !nrow(res_df)) return(NULL)
  p_cres <- ggplot2::ggplot(res_df, ggplot2::aes(x = x, y = z)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
    ggplot2::geom_hline(yintercept = c(-1.96, 1.96), linetype = "dashed",
                        colour = "grey60")
  # A TREND, but only where one means anything. Regressing z on a covariate
  # every study conditions at one of two levels fits a line through
  # territory that has no patients in it, and reports as a slope what is
  # really a difference between two groups. Continuous covariates only, and
  # only with more than two studies -- two points define a line through
  # themselves and say nothing.
  # PER COVARIATE, not across the figure: facet_wrap fits a separate `lm`
  # per panel, so a single covariate with 3+ studies used to add the layer
  # for all of them and draw a "trend" through another covariate's two
  # points -- exactly what the paragraph above says must not happen.
  cont_r <- res_df[!res_df$disc, , drop = FALSE]
  if (nrow(cont_r))
    cont_r <- cont_r[cont_r$cov %in%
                       names(which(table(cont_r$cov) > 2L)), , drop = FALSE]
  if (nrow(cont_r))
    p_cres <- p_cres + ggplot2::geom_smooth(
      data = cont_r, method = "lm", formula = y ~ x, se = FALSE,
      na.rm = TRUE, colour = "#2166AC", linewidth = 0.7,
      linetype = "longdash")
  # A CONTRAST, where a trend is meaningless. `stratify` splits one source
  # into strata that differ only in this covariate, so the pair is what the
  # banding bought: the same patients, the same study, one covariate moved.
  # Joining them shows each source's own within-source contrast, and
  # several sources tilting the same way is the mis-specification -- which a
  # regression over the pooled cloud of strata cannot show, because it
  # averages the pairs away.
  pair_r <- res_df[res_df$disc & res_df$paired, , drop = FALSE]
  # Neutral: the line joins two studies, so colouring it by `study` would
  # split one connector across two colours. It is a connector, not a series.
  if (nrow(pair_r))
    p_cres <- p_cres + ggplot2::geom_line(
      data = pair_r, ggplot2::aes(group = source),
      colour = "grey45", linewidth = 0.6, alpha = 0.7)
  # Studies go in a LEGEND rather than inline text: banded sources sit at
  # nearly the same z, and six strata printed in place overlap into a smear.
  # The 10th-90th span a marginalised study speaks for, so a point plotted
  # at its centre is not read as a study that only ever saw that value.
  # Conditioned studies have a zero-width span and drop out of the layer.
  marg_r <- res_df[res_df$kind == "marginal" & res_df$xhi > res_df$xlo, ,
                   drop = FALSE]
  if (nrow(marg_r))
    p_cres <- p_cres + ggplot2::geom_segment(
      data = marg_r, inherit.aes = FALSE,
      ggplot2::aes(x = xlo, xend = xhi, y = z, yend = z, colour = source),
      linewidth = 1.1, alpha = 0.45)
  p_cres +
    # Colour on the SOURCE, not the stratum, so a source keeps one colour
    # across both panels -- the effect panel draws one mark per source too.
    # Strata of one source are told apart by where they sit and by the grey
    # line joining them, which is the whole point of their being a pair.
    ggplot2::geom_point(ggplot2::aes(size = n, colour = source,
                                     shape = kind), alpha = 0.85) +
    .admCovSizeScale(res_df$n) +
    .admCovPanelStyle(pal, res_df$kind) +
    # Breaks from the CONDITIONAL rows only. A marginal source on a level axis
    # sits at the mean of the mixture it declared -- 0.5 for an even split --
    # and that is not a level, so feeding it in ticked a binary covariate at
    # 0.5047619.
    ggplot2::scale_x_continuous(
      breaks = .admLevelBreaks(
        res_df$x[res_df$disc & res_df$kind == "conditional"])) +
    ggplot2::facet_wrap(~ cov, scales = "free_x", nrow = 1L) +
    ggplot2::labs(
      title = "Between-study residual vs covariate",
      x = "Covariate value (point = centre, bar = 10th-90th)",
      y = "Mean standardised residual",
      subtitle = paste("a SLOPE across sources, or a consistent TILT in the",
                       "grey within-source pairs, is a mis-specified",
                       "covariate form",
                       "\narea = n  |  dashed blue: lm  |  +/-1.96 is",
                       "indicative, not a test: z is averaged over correlated",
                       "times"))
}

## ---- putting a banded source back together ---------------------------------
##
## `stratify` is a LIKELIHOOD device. A covariate a study marginalises over is
## not identified against a random effect on the same parameter -- its effect
## enters only through the mixture it induces, which is exactly what the random
## effect does -- so admixr2 bands the source into one stratum per level and the
## objective sums over them. The strata have to be the unit there.
##
## They are not the unit a READER recognises. `A_s1` is an index into an
## internal expansion; the paper is `A`, and "does the fit reproduce this
## source" is a question about the paper. So every panel asking it puts the
## strata back together first.
##
## The mean is the n-weighted mean of the strata. THE VARIANCE IS NOT. It is the
## law of total variance -- within PLUS BETWEEN -- and the between term is the
## covariate effect the banding created. Verified against the same source fitted
## unbanded: with both terms the collapsed SD matches to 1.000 at every time;
## with the within term alone it reads 0.80-0.88, so dropping it would draw a
## correctly specified fit as under-predicting the reported spread by a fifth.
.admMixMoments <- function(E_list, V_list, w) {
  w  <- w / sum(w)
  E  <- Reduce(`+`, Map(function(e, a) a * as.numeric(e), E_list, w))
  Vw <- Reduce(`+`, Map(function(v, a) a * as.matrix(v), V_list, w))
  Vb <- Reduce(`+`, Map(function(e, a) a * tcrossprod(as.numeric(e) - E),
                        E_list, w))
  list(E = E, V = Vw + Vb)
}

## Source-level `studies` and `agg`, for the panels that are about a paper.
##
## The covariate panels keep the STRATA -- the between-level contrast is the
## whole signal there, and collapsing would destroy the thing they exist to
## show. So this returns new lists rather than replacing the originals.
##
## A single-stratum source collapses to itself (the between term is zero) and
## keeps its `cov`, so an `at`-pinned source still titles with the value it was
## solved at. Only a genuinely banded source loses that, because it no longer
## sits at one value of the covariate it was banded on.
.admCollapseSources <- function(studies, agg) {
  if (!length(studies) || is.null(names(studies)))
    return(list(studies = studies, agg = agg))
  grp <- split(names(studies), .admCovSource(names(studies)))
  st2 <- list(); ag2 <- list()
  for (sn in names(grp)) {
    ks  <- grp[[sn]]
    # Strata whose simulation failed drop out and the weights renormalise over
    # what is left: a partial collapse is a worse answer than a whole one, but
    # it is a better answer than none, and the alternative is losing the source.
    has <- ks[vapply(ks, function(k) !is.null(agg[[k]]), logical(1))]
    if (!length(has)) { st2[[sn]] <- studies[[ks[1L]]]; next }
    nk <- vapply(has, function(k) {
      v <- suppressWarnings(as.numeric(studies[[k]][["n"]] %||% NA_real_)[1L])
      if (!is.finite(v) || v <= 0) 1 else v
    }, double(1))

    obs <- .admMixMoments(lapply(has, function(k) studies[[k]]$E),
                          lapply(has, function(k) studies[[k]]$V), nk)
    prd <- .admMixMoments(lapply(has, function(k) agg[[k]]$pred$E),
                          lapply(has, function(k) agg[[k]]$pred$V), nk)

    s0 <- studies[[has[1L]]]
    s0$E <- obs$E; s0$V <- obs$V; s0$n <- sum(nk)
    if (length(has) > 1L) {
      # A banded source sits at no single value of what it was banded on, so a
      # `[SEX = 0]` title would be a claim about half of it.
      s0[["cov"]] <- NULL; s0[["cov_dist"]] <- NULL
      s0[[".adm_cov_dropped"]] <- NULL
    }
    st2[[sn]] <- s0
    ag2[[sn]] <- utils::modifyList(agg[[has[1L]]],
                                   list(obs = obs, pred = prd))
  }
  list(studies = st2, agg = ag2[names(st2)])
}

#' Diagnostic plots for an admixr2 fit
#'
#' Generates up to five diagnostic panels:
#'
#' 1. `"mean"` -- Observed vs predicted mean per study (2x2 grid). Upper row:
#'    observed and predicted mean lines with +/-1 SD ribbon on a shared y scale
#'    (black throughout). Lower row: raw residual lollipop with +/-2 SE band and
#'    standardised residual z-scores with +/-1.96 reference lines.
#' 2. `"cov"` -- Observed vs predicted (co)variance heatmaps per study (2x2
#'    grid). Upper row shares a common colour scale (blue-white-red). Lower row
#'    uses distinct diverging scales: residual (red-white-green) and
#'    standardised residual (gold-white-purple). Significance stars overlaid on
#'    the standardised residual panel.
#' 3. `"covariate"` -- Two BETWEEN-study covariate panels, one facet per
#'    covariate the model reads. `covariate_effect` sweeps each covariate across
#'    the range the sources between them cover and draws every model parameter
#'    that reads it, at the fitted thetas, once per level of any *conditioned*
#'    covariate and with the continuous others at the pooled median -- so the
#'    swept covariate's effect is the slope and a conditioned one's is the gap
#'    between the lines; the region no source sampled is shaded as extrapolation,
#'    and a discrete covariate is drawn on its levels rather than swept through
#'    the values between them. `covariate_resid` plots each study's mean
#'    standardised residual against the covariate value it sits at, with an `lm`
#'    trend: a slope there is a mis-specified covariate form. Produced only for
#'    a fit whose studies declare covariates, so it is silently absent
#'    otherwise.
#' 4. `"nll"` -- NLL trace per restart over optimizer evaluations. Restarts
#'    coloured with the Okabe-Ito palette.
#' 5. `"par"` -- Parameter trace per restart on the natural scale (struct thetas
#'    back-transformed, sigma as SD, omega diagonal as variance labelled
#'    `V(eta.x)`). Facets ordered as in the model `ini()` block. Restarts
#'    coloured with the Okabe-Ito palette.
#'
#' @section Why the covariate panels are separate:
#' With aggregate data a covariate effect is identified BETWEEN studies: the
#' renal exponent in `vignette("covariates")` is recovered from three cohorts
#' sitting at three different creatinine-clearance medians, none of which fitted
#' a renal term at all. The `"mean"` and `"cov"` panels are per study and in
#' isolation, so every source can sit beautifully on its own panel while the
#' relation tying them together is wrong. That contrast is what `"covariate"`
#' draws.
#'
#' @section Marginal and conditioned sources:
#' Both panels distinguish how each study enters a covariate, because the two
#' are different statements about the source rather than degrees of the same
#' one. A **marginalised** covariate -- one the estimator integrates over a
#' declared distribution -- is drawn as a round point at the study's median with
#' a 10th-90th percentile bar and a 2.5th-97.5th whisker: the source constrains
#' the effect through the spread it induces, and the whole span is what it
#' speaks for. A **conditioned** covariate -- `stratify`, `at` or `by`, so the
#' model is solved at one value -- is drawn as a diamond at that value with no
#' spread, because it has none; its contribution is a point of contrast against
#' the other studies. A banded source therefore appears as several diamonds on
#' the covariate it was banded on, and as one merged marginal mark on every
#' other, where its strata coincide.
#'
#' @section Reading `covariate_resid`:
#' A **continuous** covariate is read as a trend: the dashed `lm` across studies,
#' where a slope is a mis-specified covariate form.
#'
#' A **discrete** one is read as a contrast instead, and is drawn that way. Its
#' axis is ticked at its levels and nowhere else, no regression is fitted
#' through it -- a line across the levels of a factor reports as a slope what is
#' a difference between groups -- and the strata `stratify` cut from one source
#' are joined, because that pairing is the evidence banding creates: the same
#' study, one covariate moved. Several sources tilting the same way is the
#' mis-specification. A regression over the pooled cloud of strata cannot show
#' it, since it averages the pairs away.
#'
#' With few sources, covariates whose study medians happen to move together
#' cannot be told apart here -- three cohorts whose weights and renal function
#' both decline will show a slope in both facets whichever one is
#' mis-specified. The panel localises the problem to a set of covariates, not to
#' one; separating them needs a source that breaks the pattern.
#'
#' @param x An `admFit` object returned by `nlmixr2()` with
#'   `est = "adfo"`, `est = "admc"`, `est = "adgh"`, or `est = "adirmc"`.
#' @param which Character vector selecting which panel types to produce.
#'   Any subset of `c("mean", "cov", "covariate", "nll", "par")`. Defaults to
#'   all five. Note `"cov"` (predicted vs observed covariance) and `"covariate"`
#'   (covariate effect) are different panels.
#' @param n_sim Number of MC samples for the final prediction. Defaults to the
#'   value used during fitting. Only used when `"mean"`, `"cov"` or
#'   `"covariate"` is in `which`.
#' @param seed Random seed for reproducibility.
#' @param ... Unused.
#'
#' @return A named list of ggplot2 objects, invisibly. Prints each selected
#'   top-level panel. For the `"mean"` and `"cov"` panels the returned list also
#'   contains each sub-panel individually so a single panel (or a few) can be
#'   extracted in code without reprinting the whole grid. Elements can be pulled
#'   out by name -- `plot(fit, which = "mean")$mean_study1_pred` or
#'   `plot(fit, which = "cov")$cov_study1_std_resid` -- or by position, with the
#'   combined 2x2 grid stored first per study
#'   (`plot(fit, which = "mean")[[1]]` is the full grid, `[1]` the length-1
#'   named sub-list). The sub-panel keys are `<type>_<study>_obs`, `_pred`,
#'   `_resid`, and `_std_resid`; the combined grid stays under `<type>_<study>`.
#'   The extra sub-panel keys are not printed on their own.
#'
#' @section Aggregate data slot:
#' Every admixr2 fit also carries the observed and predicted aggregate data in
#' `fit$env$aggData`, a named list with one entry per study. Each entry holds the
#' observation `times`, the study `n`, and two moment sets -- `obs` (from the
#' data) and `pred` (predicted at the fitted parameters) -- each a list with the
#' mean vector `E` and the (co)variance matrix `V`:
#' \preformatted{
#'   fit$env$aggData$study1$obs$E    # observed mean vector
#'   fit$env$aggData$study1$obs$V    # observed covariance matrix
#'   fit$env$aggData$study1$pred$E   # predicted mean vector
#'   fit$env$aggData$study1$pred$V   # predicted covariance matrix
#' }
#' The predicted moments are computed by one MC simulation at the fitted
#' parameters using the fit's own `n_sim` and a fixed seed, so they match the
#' default `plot(fit)` mean/cov panels. The slot is absent only when the fit
#' cannot be simulated (no simulation model available).
#'
#' @section nlmixr2 `traceplot()`:
#' admixr2 fits also plug into the nlmixr2 `traceplot()` generic. During fitting
#' the parameter iteration history of the best restart is stored on the fit in
#' the standard `parHistData` slot (natural scale), so `traceplot(fit)` produces
#' the familiar per-parameter, free-y facetted trace used elsewhere in the
#' nlmixr2 ecosystem. There is no burn-in marker (admixr2 records optimizer
#' evaluations, not SAEM iterations), and only the best restart is shown -- the
#' per-restart overlay and the NLL trace remain available via
#' `plot(fit, which = c("par", "nll"))`. The trace stores only improving
#' evaluations (steps that lowered the best NLL), so the `iter` axis indexes
#' those improvement steps rather than raw optimizer iterations.
#'
#' @examples
#' \donttest{
#' library(rxode2)
#' library(nlmixr2)
#'
#' data("examplomycin")
#' obs    <- examplomycin[examplomycin$EVID == 0, ]
#' obs    <- obs[order(obs$ID, obs$TIME), ]
#' times  <- sort(unique(obs$TIME))
#' ids    <- unique(obs$ID)
#' dv_mat <- do.call(rbind, lapply(ids, function(i) {
#'   sub <- obs[obs$ID == i, ]; sub$DV[order(sub$TIME)]
#' }))
#' E <- colMeans(dv_mat)
#' V <- cov.wt(dv_mat, method = "ML")$cov
#'
#' pk_model <- function() {
#'   ini({
#'     tcl <- log(5); tv <- log(30)
#'     prop.sd <- c(0, 0.2)
#'     eta.cl ~ 0.09; eta.v ~ 0.04
#'   })
#'   model({
#'     cl <- exp(tcl + eta.cl)
#'     v  <- exp(tv  + eta.v)
#'     d/dt(central) <- -(cl/v) * central
#'     cp <- central / v
#'     cp ~ prop(prop.sd)
#'   })
#' }
#'
#' fit <- nlmixr2(
#'   pk_model, admData(), est = "adfo",
#'   control = adfoControl(
#'     studies = list(study1 = list(E = E, V = V, n = length(ids),
#'                                  times = times, ev = et(amt = 100))),
#'     maxeval = 100L
#'   )
#' )
#' plot(fit)
#' }
#'
#' @export
plot.admFit <- function(x, which = c("mean", "cov", "covariate", "nll", "par"),
                        n_sim = NULL, seed = 1L, ...) {
  # "cov" and "covariate" are both exact choices, so pmatch takes each to
  # itself; they are different panels (residual (co)variance vs covariate
  # effect) and the near-collision is the price of not renaming "cov".
  which <- match.arg(which, c("mean", "cov", "covariate", "nll", "par"),
                     several.ok = TRUE)
  fit <- x
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 required for plot.admFit", call. = FALSE)

  extra <- fit$env$admExtra %||% fit$env$adirmcExtra %||%
    stop("No admExtra/adirmcExtra on fit object", call. = FALSE)

  studies  <- extra$studies
  n_sim    <- n_sim %||% extra$n_sim %||% 5000L

  # The covariate residual panel reads the same predicted moments the mean and
  # cov panels do, so asking for "covariate" alone still has to simulate.
  # Which covariates earn a facet -- read off the studies, so it costs nothing
  # and is known BEFORE deciding whether to simulate.
  cov_nms <- if ("covariate" %in% which)
    .admCovPanelCovs(fit$env$ui, studies) else character(0)
  # The covariate RESIDUAL panel reads the same predicted moments the mean and
  # cov panels do, so asking for "covariate" has to simulate -- but only when
  # there is a covariate to plot. A fit that declares none would otherwise pay
  # for a full n_sim simulation and then draw nothing with it.
  need_sim_local <- any(c("mean", "cov") %in% which) || length(cov_nms) > 0L
  # Observed + predicted aggregate moments per study (mean vector + cov matrix).
  # Reuse the fit's stored `aggData` slot when it matches the requested n_sim/seed
  # (avoids a redundant simulation); otherwise recompute via the shared helper.
  agg <- if (!need_sim_local) {
    setNames(vector("list", length(studies)), names(studies))
  } else {
    cached <- fit$env$aggData
    # The stored slot was built at n_sim = extra$n_sim and seed 1L. Compare
    # numerically (not via identical()) so a double n_sim -- e.g. plot(fit,
    # n_sim = 5000) against a stored 5000L -- still hits the cache.
    if (!is.null(cached) &&
        isTRUE(n_sim == (extra$n_sim %||% 5000L)) && isTRUE(seed == 1L))
      cached
    else
      .admAggData(extra, fit$env$ui, n_sim = n_sim, seed = seed, warn = TRUE)
  }

  plots <- list()
  # Keys to actually display. Individual mean/cov sub-panels are added to
  # `plots` for programmatic extraction but are not printed on their own -- only
  # the combined 2x2 grid (or, without patchwork, the sub-panel list) is shown.
  print_keys <- character(0)

  # Fail loudly rather than silently overwrite: a study whose name ends in a
  # reserved suffix (e.g. "s1_pred") would derive a key that collides with
  # another study's sub-panel key. Only possible with pathological names.
  .check_panel_keys <- function(keys) {
    dup <- keys[keys %in% names(plots)]
    if (length(dup))
      stop("plot.admFit: panel keys collide with existing entries (",
           paste(dup, collapse = ", "),
           "); rename the study to avoid the reserved suffixes ",
           "_obs/_pred/_resid/_std_resid.", call. = FALSE)
  }

  # -- Mean diagnostics: 2x2 grid (Obs | Pred / Residual | Standardised residual)
  # Obs/Pred: shared y scale; black mean line + point + \u00b11 SD ribbon (black, alpha 0.15).
  # The ribbon is sqrt(diag(V_obs)) on the ML scale -- what the fit consumed --
  # so a study declared `v_denom = "unbiased"` shows a ribbon slightly narrower
  # than the SD read off its source figure. That is deliberate: it is compared
  # against sqrt(diag(V_pred)), a population quantity, and pairing an (n-1)
  # observed SD with it would build a 1/(2n) mismatch into the diagnostic.
  # Residual: raw (E_obs - mu_pred) lollipop with \u00b12 SE band (SE = sqrt(V_pred[t,t]/n)).
  # Standardised residual: z[t] = (E_obs[t] - mu[t]) / sqrt(V_pred[t,t]/n) ~ N(0,1).
  # Stars: |z| > 1.96 (*), > 2.58 (**), > 3.29 (***). Requires patchwork for 2x2.
  # The mean and cov panels are about a PAPER, so they read the source-level
  # lists: a banded source is put back together first, by the mixture law. The
  # covariate panels below keep the strata, where the between-level contrast is
  # the signal. See .admCollapseSources().
  .src <- if (any(c("mean", "cov") %in% which))
    .admCollapseSources(studies, agg) else list(studies = studies, agg = agg)
  studies_src <- .src$studies
  agg_src     <- .src$agg

  if ("mean" %in% which) for (nm in names(studies_src)) {
    s   <- studies_src[[nm]]
    ag  <- agg_src[[nm]]
    if (is.null(ag)) next

    n_obs      <- s$n
    mu         <- ag$pred$E
    V_pred     <- ag$pred$V
    pred_sd    <- sqrt(diag(V_pred))
    obs_sd     <- sqrt(diag(s$V))
    resid_mean <- as.numeric(s$E) - mu
    se_mean    <- sqrt(diag(V_pred) / n_obs)
    z_mean     <- resid_mean / se_mean
    sig_mean   <- ifelse(abs(z_mean) > 3.29, "***",
                  ifelse(abs(z_mean) > 2.58, "**",
                  ifelse(abs(z_mean) > 1.96, "*", "")))

    df_obs  <- data.frame(time     = s$times,
                          obs_mean = as.numeric(s$E),
                          obs_lo   = as.numeric(s$E) - obs_sd,
                          obs_hi   = as.numeric(s$E) + obs_sd)
    df_pred <- data.frame(time      = s$times,
                          pred_mean = mu,
                          pred_lo   = mu - pred_sd,
                          pred_hi   = mu + pred_sd)
    df_res  <- data.frame(time  = s$times,
                          resid = resid_mean,
                          lo    = -2 * se_mean,
                          hi    =  2 * se_mean)
    df_z    <- data.frame(time    = s$times,
                          z       = z_mean,
                          z_label = sig_mean,
                          z_vjust = ifelse(z_mean >= 0, -0.5, 1.5))

    p_obs <- ggplot2::ggplot(df_obs, ggplot2::aes(x = time)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = obs_lo, ymax = obs_hi),
                           fill = "black", alpha = 0.15) +
      ggplot2::geom_line(ggplot2::aes(y = obs_mean), colour = "black", linewidth = 1) +
      ggplot2::geom_point(ggplot2::aes(y = obs_mean), colour = "black", size = 2) +
      ggplot2::labs(title = "Observed", x = NULL, y = "Concentration",
                    subtitle = "ribbon +/-1 SD  [sqrtdiag(V_obs)]  |  shared y scale") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 7, colour = "grey40",
                                                           face = "plain"))

    p_pred <- ggplot2::ggplot(df_pred, ggplot2::aes(x = time)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = pred_lo, ymax = pred_hi),
                           fill = "black", alpha = 0.15) +
      ggplot2::geom_line(ggplot2::aes(y = pred_mean), colour = "black", linewidth = 1) +
      ggplot2::geom_point(ggplot2::aes(y = pred_mean), colour = "black", size = 2) +
      ggplot2::labs(title = "Predicted", x = NULL, y = NULL,
                    subtitle = "ribbon +/-1 SD  [sqrtdiag(V_pred);  BSV + sigma]  |  shared y scale") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 7, colour = "grey40",
                                                           face = "plain"))

    p_res <- ggplot2::ggplot(df_res, ggplot2::aes(x = time)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi),
                           fill = "black", alpha = 0.15) +
      ggplot2::geom_segment(ggplot2::aes(y = resid, xend = time, yend = 0),
                            colour = "grey50") +
      ggplot2::geom_point(ggplot2::aes(y = resid), colour = "black", size = 3) +
      ggplot2::labs(title = "Residual", x = "Time",
                    y = "E_obs - mu_pred",
                    subtitle = "band +/-2 SE(mean)  [SE = sqrt(V_pred[t,t]/n)]") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 7, colour = "grey40",
                                                           face = "plain"))

    p_z <- ggplot2::ggplot(df_z, ggplot2::aes(x = time, y = z)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
      ggplot2::geom_hline(yintercept = c(-1.96, 1.96),
                          linetype = "dashed", colour = "grey60") +
      ggplot2::geom_segment(ggplot2::aes(xend = time, yend = 0), colour = "grey50") +
      ggplot2::geom_point(size = 3, colour = "black") +
      ggplot2::geom_text(ggplot2::aes(label = z_label, vjust = z_vjust),
                         size = 4, colour = "black", fontface = "bold") +
      # Extra vertical headroom so significance stars placed above/below the
      # extreme points (via z_vjust) are not clipped at the panel edge.
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.15)) +
      ggplot2::labs(title = "Standardised residual",
                    x = "Time", y = "z-score",
                    subtitle = "dashed +/-1.96  |  z = resid/SE  |  *p<.05 **p<.01 ***p<.001") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 7, colour = "grey40",
                                                           face = "plain"))

    # Combined 2x2 grid first so positional extraction returns the whole panel
    # (plot(fit, which = "mean")[[1]]); the individual sub-panels follow under
    # their own keys for one-at-a-time extraction, e.g. $mean_study1_pred.
    .check_panel_keys(paste0("mean_", nm, c("", "_obs", "_pred", "_resid", "_std_resid")))
    if (requireNamespace("patchwork", quietly = TRUE)) {
      plots[[paste0("mean_", nm)]] <- (p_obs | p_pred) / (p_res | p_z) +
        patchwork::plot_annotation(
          title = paste(.admStudyTitle(s, nm), "-- Mean diagnostics"))
    } else {
      plots[[paste0("mean_", nm)]] <- list(obs = p_obs, pred = p_pred,
                                            resid = p_res, std_resid = p_z)
    }
    plots[[paste0("mean_", nm, "_obs")]]       <- p_obs
    plots[[paste0("mean_", nm, "_pred")]]      <- p_pred
    plots[[paste0("mean_", nm, "_resid")]]     <- p_res
    plots[[paste0("mean_", nm, "_std_resid")]] <- p_z
    print_keys <- c(print_keys, paste0("mean_", nm))
  }

  # -- Covariance heatmaps: 2x2 grid ----------------------------------------
  # Top: Observed | Predicted -- shared colour scale (blue-white-red, cov_lim).
  # Bottom: Residual (red-white-green) | Standardised residual (gold-white-purple).
  # Standardised residual SEs (asymptotic MVN): diag sqrt(2*V^2/(n-1)), off-diag sqrt((V_ii*V_jj+V_ij^2)/(n-1)).
  # Stars: |z| > 1.96 (*), > 2.58 (**), > 3.29 (***). No multiple-testing correction.
  # Requires patchwork for 2x2; falls back to 4 separate plots.
  if ("cov" %in% which) {
  .mat_df <- function(mat, times, lower_only = FALSE) {
    n <- length(times)
    if (lower_only) mat[upper.tri(mat)] <- NA_real_
    df <- data.frame(t_row = factor(rep(times, times = n), levels = rev(times)),
                     t_col = factor(rep(times, each  = n), levels = times),
                     value = as.vector(mat),
                     stringsAsFactors = FALSE)
    if (lower_only) df <- df[!is.na(df$value), ]
    df
  }
  .heat_tile <- function(df, lim, fill_name, low = "#2166AC", high = "#D6604D") {
    ggplot2::ggplot(df, ggplot2::aes(x = t_col, y = t_row, fill = value)) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
      ggplot2::scale_fill_gradient2(low = low, mid = "white", high = high,
                                    midpoint = 0, limits = c(-lim, lim), name = fill_name) +
      ggplot2::labs(x = NULL, y = NULL) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1),
                     plot.subtitle = ggplot2::element_text(size = 7, colour = "grey40",
                                                           face = "plain"))
  }

  for (nm in names(studies_src)) {
    s   <- studies_src[[nm]]
    ag  <- agg_src[[nm]]
    if (is.null(ag)) next

    n_obs  <- s$n
    mu     <- ag$pred$E
    V_pred <- ag$pred$V
    v_diag <- diag(V_pred)
    n_t    <- length(s$times)
    times  <- s$times

    resid_mat <- s$V - V_pred
    z_mat <- matrix(NA_real_, n_t, n_t)
    for (i in seq_len(n_t))
      for (j in seq_len(n_t)) {
        se_ij <- if (i == j)
          sqrt(2 * v_diag[i]^2 / (n_obs - 1L))
        else
          sqrt((v_diag[i] * v_diag[j] + V_pred[i, j]^2) / (n_obs - 1L))
        z_mat[i, j] <- resid_mat[i, j] / se_ij
      }

    sig_mat <- ifelse(abs(z_mat) > 3.29, "***",
               ifelse(abs(z_mat) > 2.58, "**",
               ifelse(abs(z_mat) > 1.96, "*", "")))

    cov_lim <- max(abs(c(s$V, V_pred)), na.rm = TRUE)
    res_lim <- max(abs(resid_mat), 1e-6, na.rm = TRUE)
    z_lim   <- max(abs(z_mat), 1.96, na.rm = TRUE)

    df_z         <- .mat_df(z_mat, times, lower_only = TRUE)
    df_z$z_label <- as.vector(sig_mat)[as.vector(!upper.tri(z_mat))]

    p_obs  <- .heat_tile(.mat_df(s$V,       times, lower_only = TRUE), cov_lim, "Cov") +
      ggplot2::ggtitle("Observed", subtitle = "sample (co)variance from data  [shared scale]")
    p_pred <- .heat_tile(.mat_df(V_pred,    times, lower_only = TRUE), cov_lim, "Cov") +
      ggplot2::ggtitle("Predicted", subtitle = "MC cov + sigma  [shared scale]")
    p_res  <- .heat_tile(.mat_df(resid_mat, times, lower_only = TRUE), res_lim, "DeltaCov",
                         low = "#4C0690", high = "#048590") +
      ggplot2::ggtitle("Residual", subtitle = "V_obs - V_pred")
    p_z    <- .heat_tile(df_z, z_lim, "  z", low = "#C101AC", high = "#D2D214") +
      ggplot2::ggtitle("Standardised residual",
                       subtitle = "z = DeltaCov/SE  |  *p<.05 **p<.01 ***p<.001") +
      ggplot2::geom_text(ggplot2::aes(label = z_label), size = 4, colour = "black", fontface = "bold")

    # Combined 2x2 grid first so positional extraction returns the whole panel
    # (plot(fit, which = "cov")[[1]]); the individual sub-panels follow under
    # their own keys for one-at-a-time extraction, e.g. $cov_study1_std_resid.
    .check_panel_keys(paste0("cov_", nm, c("", "_obs", "_pred", "_resid", "_std_resid")))
    if (requireNamespace("patchwork", quietly = TRUE)) {
      plots[[paste0("cov_", nm)]] <- (p_obs | p_pred) / (p_res | p_z) +
        patchwork::plot_annotation(
          title = paste(.admStudyTitle(s, nm), "-- Covariance diagnostics"))
    } else {
      plots[[paste0("cov_", nm)]] <- list(obs = p_obs, pred = p_pred,
                                           resid = p_res, std_resid = p_z)
    }
    plots[[paste0("cov_", nm, "_obs")]]       <- p_obs
    plots[[paste0("cov_", nm, "_pred")]]      <- p_pred
    plots[[paste0("cov_", nm, "_resid")]]     <- p_res
    plots[[paste0("cov_", nm, "_std_resid")]] <- p_z
    print_keys <- c(print_keys, paste0("cov_", nm))
  }
  } # end if ("cov" %in% which)

  # -- Covariate panels ------------------------------------------------------
  # Both are BETWEEN-study plots; see the block comment on .admCovEffectData().
  # Silently absent for a fit with no covariates, which is why "covariate" can
  # sit in the default `which` without changing what a non-covariate fit prints.
  if (length(cov_nms)) {
    # The ORIGINAL studies, off the estimator's control: each source's own
    # published model, which is what puts a source at its own parameter value
    # rather than on the fitted line. See .admFitSourceStudies().
    src_st <- tryCatch(.admFitSourceStudies(fit), error = function(e) NULL)
    eff <- Filter(Negate(is.null), lapply(cov_nms, function(cv)
      tryCatch(.admCovEffectData(fit$env$ui, cv, studies, extra$struct, src_st),
               error = function(e) NULL)))
    res <- Filter(Negate(is.null), lapply(cov_nms, function(cv)
      tryCatch(.admCovResidData(cv, studies, agg), error = function(e) NULL)))

    # Each source's OWN published claim, read off the UNEXPANDED studies the
    # control kept. Restricted to the facets the pooled curve draws, so a
    # source cannot conjure a panel of its own.
    # Both panels' colour scale, from the union of the names they use. See
    # .admCovPalette(): the sets differ, and the figure is read as a pair.
    cov_pal <- .admCovPalette(c(
      unlist(lapply(eff, function(z) z$marks$study), use.names = FALSE),
      unlist(lapply(res, `[[`, "source"), use.names = FALSE)))

    p_eff <- .admCovEffectPanel(eff, cov_pal)
    if (!is.null(p_eff)) {
      plots[["covariate_effect"]] <- p_eff
      print_keys <- c(print_keys, "covariate_effect")
    }
    p_cres <- .admCovResidPanel(res, cov_pal)
    if (!is.null(p_cres)) {
      plots[["covariate_resid"]] <- p_cres
      print_keys <- c(print_keys, "covariate_resid")
    }
  }

  # -- NLL trace per restart -------------------------------------------------
  all_traces <- extra$all_traces
  par_names  <- extra$par_names

  if (any(c("nll", "par") %in% which) && !is.null(all_traces) && length(all_traces) > 0) {

  pinfo_pt       <- tryCatch(.admParseIniDf(fit$env$ui$iniDf, fit$env$ui), error = function(e) NULL)
  ini_df_pt      <- tryCatch(fit$env$ui$iniDf, error = function(e) NULL)
  spec_pt        <- .admTraceDisplaySpec(pinfo_pt, par_names, ini_df_pt)
  disp_nms       <- spec_pt$disp_nms
  back_fns       <- spec_pt$back_fns
  pt_param_order <- spec_pt$param_order

  if ("nll" %in% which) {
    df_nll <- do.call(rbind, lapply(all_traces, function(tr) {
      nt <- tr$nll_trace
      if (is.null(nt) || length(nt) == 0) return(NULL)
      data.frame(eval    = seq_along(nt),
                 nll     = nt,
                 restart = factor(tr$restart_id))
    }))

    if (!is.null(df_nll) && nrow(df_nll) > 0) {
      plots[["nll_trace"]] <-
        ggplot2::ggplot(df_nll, ggplot2::aes(x = eval, y = nll, colour = restart)) +
        ggplot2::geom_line() +
        ggplot2::scale_colour_manual(
          values = setNames(.admOkabeIto(nlevels(df_nll$restart)),
                            levels(df_nll$restart)),
          name = "Restart") +
        ggplot2::labs(title = "NLL trace per restart",
                      subtitle = "each line = one restart; lower = better",
                      x = "Iteration", y = "-2LL") +
        ggplot2::theme_bw() +
        ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 7, colour = "grey40",
                                                             face = "plain"))
      print_keys <- c(print_keys, "nll_trace")
    }
  } # end if ("nll" %in% which)

    # -- Parameter trace per restart -----------------------------------------
  if ("par" %in% which && !is.null(par_names)) {
      df_par <- do.call(rbind, lapply(all_traces, function(tr) {
        pt <- tr$par_trace
        if (is.null(pt) || nrow(pt) == 0) return(NULL)
        df <- as.data.frame(pt)
        colnames(df) <- par_names
        df$iter    <- seq_len(nrow(df))
        df$restart <- factor(tr$restart_id)
        do.call(rbind, lapply(par_names, function(pnm) {
          fn  <- if (!is.null(back_fns)) back_fns[[pnm]] else identity
          data.frame(iter    = df$iter,
                     restart = df$restart,
                     param   = if (!is.null(disp_nms)) disp_nms[[pnm]] else pnm,
                     value   = fn(df[[pnm]]))
        }))
      }))

      if (!is.null(df_par) && nrow(df_par) > 0) {
        if (!is.null(pt_param_order))
          df_par$param <- factor(df_par$param, levels = pt_param_order)
        plots[["par_trace"]] <-
          ggplot2::ggplot(df_par, ggplot2::aes(x = iter, y = value, colour = restart)) +
          ggplot2::geom_line() +
          ggplot2::scale_colour_manual(
            values = setNames(.admOkabeIto(nlevels(df_par$restart)),
                              levels(df_par$restart)),
            name = "Restart") +
          ggplot2::facet_wrap(~param, scales = "free_y") +
          ggplot2::labs(title = "Parameter trace per restart",
                        subtitle = "Natural scale; struct = back-transformed, sigma = SD, V(eta) = variance, off-diagonal = Cholesky L[i,j]",
                        x = "Iteration", y = "Value") +
          ggplot2::theme_bw() +
          ggplot2::theme(strip.text    = ggplot2::element_text(size = 8),
                         plot.subtitle = ggplot2::element_text(size = 8, colour = "grey40",
                                                               face = "plain"))
        print_keys <- c(print_keys, "par_trace")
      }
    }
  }

  # Print only the top-level panels (combined mean/cov grids, nll, par). The
  # individual mean/cov sub-panels remain in `plots` for programmatic extraction
  # but are not printed here to avoid duplicating the combined grid output.
  for (key in print_keys) {
    p <- plots[[key]]
    if (is.list(p) && !inherits(p, "gg"))
      for (pp in p) print(pp)
    else
      print(p)
  }
  invisible(plots)
}

