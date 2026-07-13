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

  # Distinguish explicit print(fit) from RStudio auto-print (chunk evaluates
  # bare `fit`). RStudio's auto-print wrapper is an anonymous function; explicit
  # print() appears as "print" in sys.calls(). Inject the rmarkdown stub only
  # for explicit calls so auto-print keeps the paged-table path.
  .stk <- sys.calls()
  .n   <- length(.stk)
  .parent_fn <- if (.n >= 2L)
    tryCatch(deparse(.stk[[.n - 1L]][[1L]])[1L], error = function(e) "")
  else ""

  if (identical(.parent_fn, "print") && isNamespaceLoaded("rmarkdown")) {
    .rm_ns <- asNamespace("rmarkdown")
    .orig_paged <- tryCatch(
      get("print.paged_df", envir = .rm_ns, inherits = FALSE),
      error = function(e) NULL
    )
    if (!is.null(.orig_paged)) {
      utils::assignInNamespace("print.paged_df",
                               function(x, ...) { cat(" "); invisible(x) },
                               ns = "rmarkdown")
      on.exit(
        utils::assignInNamespace("print.paged_df", .orig_paged, ns = "rmarkdown"),
        add = TRUE
      )
    }
  }

  fn <- get("print.nlmixr2FitCore", envir = asNamespace("nlmixr2est"), inherits = FALSE)
  class(x) <- class(x)[class(x) != "nlmixr2FitData"]
  fn(x, ...)
  invisible(x)
}

# knitr auto-print handler. Registered as knit_print.admFit into knitr's
# namespace by .register_knit_print() in .onLoad (zzz.R). Not named
# knit_print.admFit to avoid roxygen2 S3-method detection and spurious
# NAMESPACE requirements. No rmarkdown stub injection -> .pagedPrint detects
# the knitr side-channel and produces paged tables.
.admKnitPrint <- function(x, ...) {
  fn <- get("print.nlmixr2FitCore", envir = asNamespace("nlmixr2est"), inherits = FALSE)
  saved_cl <- class(x)
  on.exit(tryCatch(class(x) <- saved_cl, error = function(e) NULL))
  class(x) <- class(x)[class(x) != "nlmixr2FitData"]
  fn(x, ...)
  invisible(NULL)
}

# Called by print.nlmixr2FitCore's paged path when it does head(fit).
# Unlike print.paged_df, this method is not in nlmixr2est's import chain, so
# S3 dispatch falls through to admixr2's method table and finds us first.
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
  sigma_pow_nms <- pinfo$sigma_names[.admSigmaRole(pinfo) == "pow_exp"]
  back_fns <- setNames(lapply(par_names, function(nm) {
    if (nm %in% struct_nms)             function(v) .admBackTransform(v, pinfo$struct_transforms[[nm]])
    else if (nm %in% omega_diag_nms)    exp
    else if (nm %in% omega_offdiag_nms) identity
    else if (nm %in% sigma_pow_nms)     identity   # pow() exponent: already natural
    else function(v) exp(v / 2)   # sigma: log(sigma^2) -> SD
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
## nlmixr2's `traceplot()` generic (`traceplot.nlmixr2FitCore` in nlmixr2plot)
## reads `fit$parHistStacked`, which `nmObjGet.parHistStacked` derives from
## `fit$env$parHistData` -- a wide data.frame with a `type` column (it keeps
## `type == "Unscaled"`), an `iter` column, and one column per parameter.
## Populating that slot is all that is required for `traceplot(fit)` to work on
## an admixr2 fit; no S3 method registration is needed because admFit already
## inherits `nlmixr2FitCore`.
##
## Semantics chosen here:
## - single chain = the best restart (lowest final NLL); nlmixr2's shape stores
##   one value per parameter per iter, so multi-restart overlay is not
##   expressible -- that stays in `plot(fit, which = "par")`.
## - natural scale under `"Unscaled"`, using the same back-transforms and
##   display names as the custom par panel (`.admTraceDisplaySpec`).
## - no burn-in marker: we leave `parHist`'s class without a `niter` attribute,
##   so `traceplot()` draws no vline (the trace records improving nloptr
##   evaluations, not SAEM iterations).
##
## Note the `iter` axis indexes improving optimizer evaluations (only steps that
## lowered the best NLL are stored), not raw nloptr iterations.
##
## Returns `NULL` when no usable trace is available.
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
  pinfo_r <- tryCatch(.admParseIniDf(ui$iniDf, ui), error = function(e) NULL)
  if (is.null(pinfo_r))
    pinfo_r <- list(sigma_names   = sig_nms,
                    sigma_output  = rep(NA_character_, length(sv)),
                    sigma_is_prop  = as.list(grepl("prop",  sig_nms, ignore.case = TRUE)),
                    sigma_is_lnorm = as.list(grepl("lnorm", sig_nms, ignore.case = TRUE)))
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
      .admSimulate(rxMod, extra$struct, sig_nms, eta_mat, s, ov, params_df, 1L),
      error = function(e) {
        if (warn) warning("plot.admFit: simulation failed: ", e$message, call. = FALSE)
        NULL
      })
  }

  # Returns BOTH the residual-adjusted covariance and mean: lnorm rescales the
  # mean, and the predicted E must carry that scaling just as the NLL does.
  .add_sigma <- function(V, mu, ov = out_var) {
    arr <- .admResidRows(pinfo_r, ov, sv, length(mu))
    ap  <- .admResidApply(mu, diag(V), arr)
    diag(V) <- ap$dv
    list(V = V, mu = ap$mu)
  }

  setNames(lapply(names(studies), function(nm) {
    s      <- studies[[nm]]
    cp_mat <- .sim_study(s)
    if (is.null(cp_mat)) return(NULL)
    mu     <- colMeans(cp_mat)
    res    <- .add_sigma(crossprod(sweep(cp_mat, 2L, mu)) / nrow(cp_mat), mu,
                         s$output %||% out_var)
    V_pred <- res$V; mu <- res$mu
    obs_E  <- as.numeric(s$E)
    obs_V  <- as.matrix(s$V)
    tnm    <- as.character(s$times)
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

#' Diagnostic plots for an admixr2 fit
#'
#' Generates up to four diagnostic panels:
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
#' 3. `"nll"` -- NLL trace per restart over optimizer evaluations. Restarts
#'    coloured with the Okabe-Ito palette.
#' 4. `"par"` -- Parameter trace per restart on the natural scale (struct thetas
#'    back-transformed, sigma as SD, omega diagonal as variance labelled
#'    `V(eta.x)`). Facets ordered as in the model `ini()` block. Restarts
#'    coloured with the Okabe-Ito palette.
#'
#' @param x An `admFit` object returned by `nlmixr2()` with
#'   `est = "adfo"`, `est = "admc"`, `est = "adgh"`, or `est = "adirmc"`.
#' @param which Character vector selecting which panel types to produce.
#'   Any subset of `c("mean", "cov", "nll", "par")`. Defaults to all four.
#' @param n_sim Number of MC samples for the final prediction. Defaults to the
#'   value used during fitting. Only used when `"mean"` or `"cov"` is in `which`.
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
plot.admFit <- function(x, which = c("mean", "cov", "nll", "par"),
                        n_sim = NULL, seed = 1L, ...) {
  which <- match.arg(which, c("mean", "cov", "nll", "par"), several.ok = TRUE)
  fit <- x
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 required for plot.admFit", call. = FALSE)

  extra <- fit$env$admExtra %||% fit$env$adirmcExtra %||%
    stop("No admExtra/adirmcExtra on fit object", call. = FALSE)

  studies  <- extra$studies
  n_sim    <- n_sim %||% extra$n_sim %||% 5000L

  need_sim_local <- any(c("mean", "cov") %in% which)
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
  # Residual: raw (E_obs - mu_pred) lollipop with \u00b12 SE band (SE = sqrt(V_pred[t,t]/n)).
  # Standardised residual: z[t] = (E_obs[t] - mu[t]) / sqrt(V_pred[t,t]/n) ~ N(0,1).
  # Stars: |z| > 1.96 (*), > 2.58 (**), > 3.29 (***). Requires patchwork for 2x2.
  if ("mean" %in% which) for (nm in names(studies)) {
    s   <- studies[[nm]]
    ag  <- agg[[nm]]
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
          title = sprintf("Study '%s' -- Mean diagnostics", nm))
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

  for (nm in names(studies)) {
    s   <- studies[[nm]]
    ag  <- agg[[nm]]
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
          title = sprintf("Study '%s' -- Covariance diagnostics", nm))
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
          values = setNames(
            rep_len(c("#000000", "#E69F00", "#56B4E9", "#009E73",
                      "#0072B2", "#D55E00", "#CC79A7"), nlevels(df_nll$restart)),
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
            values = setNames(
              rep_len(c("#000000", "#E69F00", "#56B4E9", "#009E73",
                        "#0072B2", "#D55E00", "#CC79A7"), nlevels(df_par$restart)),
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

