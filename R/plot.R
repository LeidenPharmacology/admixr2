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
  back_fns <- setNames(lapply(par_names, function(nm) {
    if (nm %in% struct_nms)             function(v) .admBackTransform(v, pinfo$struct_transforms[[nm]])
    else if (nm %in% omega_diag_nms)    exp
    else if (nm %in% omega_offdiag_nms) identity
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
##   so `traceplot()` draws no vline (we have nloptr evals, not SAEM iterations).
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

  cols <- lapply(par_names, function(nm) {
    fn <- if (!is.null(spec)) spec$back_fns[[nm]] else identity
    as.numeric(fn(df[[nm]]))
  })
  disp <- if (!is.null(spec))
    vapply(par_names, function(nm) spec$disp_nms[[nm]], character(1)) else par_names
  names(cols) <- disp

  # Follow iniDf facet order when available so traceplot panels match the
  # custom par panel.
  if (!is.null(spec) && !is.null(spec$param_order))
    cols <- cols[spec$param_order]

  data.frame(type = "Unscaled",
             iter = seq_len(nrow(df)),
             cols,
             check.names = FALSE,
             stringsAsFactors = FALSE)
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
#'   `est = "adfo"`, `est = "admc"`, or `est = "adirmc"`.
#' @param which Character vector selecting which panel types to produce.
#'   Any subset of `c("mean", "cov", "nll", "par")`. Defaults to all four.
#' @param n_sim Number of MC samples for the final prediction. Defaults to the
#'   value used during fitting. Only used when `"mean"` or `"cov"` is in `which`.
#' @param seed Random seed for reproducibility.
#' @param ... Unused.
#'
#' @return A named list of ggplot2 objects, invisibly. Prints each selected plot.
#'
#' @section nlmixr2 `traceplot()`:
#' admixr2 fits also plug into the nlmixr2 `traceplot()` generic. During fitting
#' the parameter iteration history of the best restart is stored on the fit in
#' the standard `parHistData` slot (natural scale), so `traceplot(fit)` produces
#' the familiar per-parameter, free-y facetted trace used elsewhere in the
#' nlmixr2 ecosystem. There is no burn-in marker (admixr2 records nloptr
#' evaluations, not SAEM iterations), and only the best restart is shown -- the
#' per-restart overlay and the NLL trace remain available via
#' `plot(fit, which = c("par", "nll"))`.
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
  omega    <- extra$omega
  n_eta    <- nrow(omega)
  L        <- extra$L %||% tryCatch(t(chol(omega)), error = function(e) NULL)
  sv        <- extra$sigma_var
  sig_prop  <- extra$sigma_is_prop  %||% grepl("prop",  names(sv), ignore.case = TRUE)
  sig_lnorm <- extra$sigma_is_lnorm %||% grepl("lnorm", names(sv), ignore.case = TRUE)
  eta_nms  <- extra$eta_col_names %||% character(0)
  sig_nms  <- names(sv)

  need_sim_local <- any(c("mean", "cov") %in% which)
  rxMod <- if (need_sim_local)
    tryCatch(fit$env$ui$simulationModel, error = function(e) NULL)
  else NULL
  rxMod_ok <- !is.null(rxMod)
  if (need_sim_local && !rxMod_ok)
    warning("plot.admFit: could not retrieve simulation model from fit object", call. = FALSE)

  .sim_study <- function(s) {
    if (!rxMod_ok) return(NULL)
    tryCatch(rxode2::rxLoad(rxMod), error = function(e) NULL)
    set.seed(seed)
    if (n_eta > 0 && !is.null(L)) {
      .samp   <- extra$sampling %||% "sobol"
      z_s     <- switch(.samp,
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
    col_nms   <- c(names(extra$struct), eta_nms, sig_nms, "rxerr.cp")
    params_df <- as.data.frame(matrix(0, nrow = n_sim, ncol = length(col_nms),
                                      dimnames = list(NULL, col_nms)))
    params_df[["rxerr.cp"]] <- 1
    tryCatch(
      .admSimulate(rxMod, extra$struct, sig_nms, eta_mat, s, "cp", params_df, 1L),
      error = function(e) { warning("plot.admFit: simulation failed: ", e$message, call. = FALSE); NULL })
  }

  .add_sigma <- function(V, mu) {
    for (k in seq_along(sv)) {
      s <- sv[[k]]
      if (sig_lnorm[k]) {
        mu_adj <- mu * exp(s / 2)
        diag(V) <- diag(V) + mu_adj^2 * (exp(s) - 1)
      } else if (sig_prop[k]) {
        diag(V) <- diag(V) + s * mu^2
      } else {
        diag(V) <- diag(V) + s
      }
    }
    V
  }

  # Simulate once per study; only when mean or cov panels are requested
  sims <- if (need_sim_local)
    setNames(lapply(names(studies), function(nm) .sim_study(studies[[nm]])),
             names(studies))
  else
    setNames(vector("list", length(studies)), names(studies))

  plots <- list()

  # -- Mean diagnostics: 2x2 grid (Obs | Pred / Residual | Standardised residual)
  # Obs/Pred: shared y scale; black mean line + point + \u00b11 SD ribbon (black, alpha 0.15).
  # Residual: raw (E_obs - mu_pred) lollipop with \u00b12 SE band (SE = sqrt(V_pred[t,t]/n)).
  # Standardised residual: z[t] = (E_obs[t] - mu[t]) / sqrt(V_pred[t,t]/n) ~ N(0,1).
  # Stars: |z| > 1.96 (*), > 2.58 (**), > 3.29 (***). Requires patchwork for 2x2.
  if ("mean" %in% which) for (nm in names(studies)) {
    s      <- studies[[nm]]
    cp_mat <- sims[[nm]]
    if (is.null(cp_mat)) next

    n_obs      <- s$n
    mu         <- colMeans(cp_mat)
    V_pred     <- .add_sigma(crossprod(sweep(cp_mat, 2L, mu)) / nrow(cp_mat), mu)
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
      ggplot2::labs(title = "Standardised residual",
                    x = "Time", y = "z-score",
                    subtitle = "dashed +/-1.96  |  z = resid/SE  |  *p<.05 **p<.01 ***p<.001") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 7, colour = "grey40",
                                                           face = "plain"))

    if (requireNamespace("patchwork", quietly = TRUE)) {
      plots[[paste0("mean_", nm)]] <- (p_obs | p_pred) / (p_res | p_z) +
        patchwork::plot_annotation(
          title = sprintf("Study '%s' -- Mean diagnostics", nm))
    } else {
      plots[[paste0("mean_", nm)]] <- list(obs = p_obs, pred = p_pred,
                                            resid = p_res, std_resid = p_z)
    }
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
    s      <- studies[[nm]]
    cp_mat <- sims[[nm]]
    if (is.null(cp_mat)) next

    n_obs  <- s$n
    mu     <- colMeans(cp_mat)
    V_pred <- .add_sigma(crossprod(sweep(cp_mat, 2L, mu)) / nrow(cp_mat), mu)
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

    if (requireNamespace("patchwork", quietly = TRUE)) {
      plots[[paste0("cov_", nm)]] <- (p_obs | p_pred) / (p_res | p_z) +
        patchwork::plot_annotation(
          title = sprintf("Study '%s' -- Covariance diagnostics", nm))
    } else {
      plots[[paste0("cov_", nm)]] <- list(obs = p_obs, pred = p_pred,
                                           resid = p_res, std_resid = p_z)
    }
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
      }
    }
  }

  for (p in plots) {
    if (is.list(p) && !inherits(p, "gg"))
      for (pp in p) print(pp)
    else
      print(p)
  }
  invisible(plots)
}

