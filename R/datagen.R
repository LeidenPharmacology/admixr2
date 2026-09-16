#' Control parameters for [datagen()]
#'
#' @param method Moment approximation used to generate `E` and `V`:
#'   `"mc"` (default) draws Monte Carlo samples over the IIV distribution, as in
#'   `est = "admc"`; `"fo"` uses the deterministic First-Order expansion
#'   (`mu = f(theta, 0)`, `V = J Omega J' + Sigma`), matching `est = "adfo"`;
#'   `"gh"` uses deterministic Gauss-Hermite quadrature over the random-effects
#'   prior, matching `est = "adgh"` -- unbiased at any IIV magnitude and
#'   noise-free.  Use `"fo"` or `"gh"` for design evaluation where the
#'   data-generating and data-analytic models must coincide.
#' @param n_sim Number of Monte Carlo samples used to approximate population
#'   moments. Ignored when `method = "fo"` or `"gh"`.
#' @param n_nodes Number of Gauss-Hermite nodes per eta dimension for
#'   `method = "gh"` (default 5). Total nodes = `n_nodes^n_eta`. Ignored for
#'   `"mc"` and `"fo"`.
#' @param resid_nodes Gauss-Hermite nodes used to integrate the RESIDUAL for a
#'   transform-both-sides endpoint (`boxCox`, `yeoJohnson`, `logitNorm`,
#'   `probitNorm`), where `y = g(h(f) + sigma*eps)` has no closed-form mean and
#'   variance. Ignored by every other error model. Default 81 -- the same default
#'   the four estimator controls use, so `datagen()` and the fit it feeds agree
#'   unless you deliberately change one of them. See [admControl()] for the
#'   measured convergence.
#' @param sampling Quasi-random sampling method: `"sobol"` (default),
#'   `"halton"`, `"torus"`, `"lhs"`, or `"rnorm"`. Ignored when `method = "fo"`
#'   or `"gh"`.
#' @param cov_nodes Gauss-Hermite nodes per covariate when `method = "gh"`
#'   integrates a study's `cov_dist` (default 7). Total covariate points are
#'   `cov_nodes^p` for `p` covariates. Ignored by `"mc"`, which draws a covariate
#'   value per simulated subject instead, and by `"fo"`, which cannot integrate a
#'   covariate at all.
#' @param seed Integer seed.  Applied before stochastic methods
#'   (`"rnorm"`, `"lhs"`). Ignored when `method = "fo"` or `"gh"`.
#' @param cores Number of `rxSolve` threads.
#' @param return_samples Include the raw `n_sim x length(times)`
#'   prediction matrix as `$samples` in each study's output. No effect when
#'   `method = "fo"` or `"gh"` (those methods draw no samples).
#'
#' @return A list of class `"datagenControl"`.
#' @seealso [datagen()]
#' @examples
#' ctrl <- datagenControl(n_sim = 2000L)
#' ctrl$sampling  # "sobol"
#'
#' # Deterministic FO moments for design evaluation:
#' datagenControl(method = "fo")$method  # "fo"
#'
#' # GH quadrature moments (unbiased, noise-free):
#' datagenControl(method = "gh", n_nodes = 5L)$n_nodes
#' @export
datagenControl <- function(
  method         = c("mc", "fo", "gh"),
  n_sim          = 5000L,
  n_nodes        = 5L,
  sampling       = c("sobol", "halton", "torus", "lhs", "rnorm"),
  seed           = 12345L,
  cores          = 1L,
  return_samples = FALSE,
  # Keep at end of signature to preserve positional argument compatibility.
  resid_nodes    = 81L,
  cov_nodes      = 7L) {
  method   <- match.arg(method)
  sampling <- match.arg(sampling)
  checkmate::assertIntegerish(n_sim,    lower = 1L, len = 1L)
  checkmate::assertIntegerish(n_nodes,  lower = 1L, len = 1L)
  checkmate::assertIntegerish(resid_nodes, lower = 5L, len = 1L)
  checkmate::assertIntegerish(cov_nodes, lower = 1L, len = 1L)
  checkmate::assertIntegerish(seed,                 len = 1L)
  checkmate::assertIntegerish(cores,    lower = 1L, len = 1L)
  checkmate::assertFlag(return_samples)
  structure(
    list(
      method         = method,
      cov_nodes      = as.integer(cov_nodes),
      n_sim          = as.integer(n_sim),
      n_nodes        = as.integer(n_nodes),
      resid_nodes    = as.integer(resid_nodes),
      sampling       = sampling,
      seed           = as.integer(seed),
      cores          = as.integer(cores),
      return_samples = return_samples
    ),
    class = "datagenControl"
  )
}


#' Generate aggregate study data from (possibly different) pharmacometric models
#'
#' Generates population mean vectors (`E`) and covariance matrices
#' (`V`) for each study by integrating over the IIV distribution -- either by
#' Monte Carlo (the default) or by a deterministic First-Order expansion
#' (`method = "fo"`, see [datagenControl()]).  Each study may specify its own PK/PD model (as would be the
#' case when digitising data from several published studies, each fit with a
#' different structural model).  True parameter values are taken from the
#' `ini()` block of each study's model.  Each element of the returned list
#' is ready to supply directly to `admControl(studies = ...)`.
#'
#' @param studies A named list of study specifications.  Each element is a list
#'   with:
#'   \describe{
#'     \item{`model`}{An nlmixr2-style model function with `ini()` and
#'       `model()` blocks.  Serves as the data-generating model for this
#'       study.  May differ between studies.  Can be omitted if a top-level
#'       default is supplied via the `model` argument.}
#'     \item{`times`}{Numeric vector of observation times.}
#'     \item{`ev`}{A dosing event table created with `rxode2::et()`.}
#'     \item{`n`}{(Optional) integer sample size; stored as metadata and
#'       used when supplying the result to `admControl()`.}
#'     \item{`cov_dist`}{(Optional) the covariate distribution this study's
#'       subjects span --- see [covDraw()] for the grammar. The generated
#'       `E`/`V` are MARGINAL over it, which is what a publication reports.
#'       Needs `datagenControl(method = "mc")` or `"gh"`; `"fo"` integrates
#'       over the random effects only and is refused. Prefer `"gh"`, which
#'       adds no Monte Carlo noise to data that is meant to BE the reference.}
#'     \item{`stratify`}{(Optional) `TRUE` to stratify on every covariate this
#'       study's OWN data-generating model conditions on, marginalising the
#'       rest --- the split is read from the model, so it cannot disagree with
#'       it. Two sources sharing one `cov_dist` therefore stratify differently,
#'       each according to what it fitted. A character vector names the
#'       covariates explicitly instead. The study is expanded into
#'       one ordinary study per covariate stratum, named `<study>_s1`,
#'       `<study>_s2`, ..., each pinned at its own covariate value, carrying its
#'       own effective size `n_k` (the quadrature weight times `n`, summing to
#'       `n`), and marginalising the remaining covariates over their
#'       distribution CONDITIONAL on that stratum. Use it when the source
#'       reports --- or, being a published model, can report --- summaries by
#'       covariate subgroup. Stratify only on what the source actually fitted:
#'       a source with no term in a covariate has no contrast in it to give, and
#'       nodes that vary it manufacture a null one. See [covStrata()].}
#'     \item{`strata_nodes`}{(Optional) strata per stratified covariate
#'       (default 5); a discrete covariate is cut at its levels instead. Each
#'       stratum costs a solve and more of them do not buy accuracy: a matched
#'       one-covariate fit recovers 0.7000 / 0.7002 / 0.7005 at 3 / 4 / 10
#'       strata against a true 0.700. See [covStrata()].}
#'     \item{`observations`}{(Optional) a named list to generate data for several
#'       observed outputs (multi-compartment). Each entry gives one output's
#'       `output` (model prediction variable, e.g. `"cp"`), `times`, and
#'       optionally `ev`/`n` (inherited from the study otherwise). When present,
#'       the study result carries a matching `observations` list of per-output
#'       `E`/`V`, ready to pass straight to `admControl(studies = ...)`.}
#'   }
#' @param model Optional default model function used for any study that does not
#'   supply its own `model` element.  At least one of `model` or each
#'   study's `model` must be non-`NULL`.
#' @param control A [datagenControl()] object.
#' @return A named list with one element per study.  Each element contains:
#'   \describe{
#'     \item{`E`}{Population mean vector at `times`.}
#'     \item{`V`}{Population covariance matrix
#'       (`length(times)` x `length(times)`; ML denominator `n_sim` for
#'       `method = "mc"`, the analytical FO covariance for `method = "fo"`,
#'       or the GH weighted covariance for `method = "gh"`).
#'       The diagonal carries the model's residual-error variance; to generate
#'       residual-free (IIV-only) moments, omit the error term from the model.}
#'     \item{`n`}{Sample size (`NA_integer_` if not supplied).}
#'     \item{`times`}{Observation times.}
#'     \item{`ev`}{Dosing event table.}
#'     \item{`samples`}{Raw `n_sim x length(times)` prediction matrix
#'       (only when `control$return_samples = TRUE`).}
#'   }
#'
#' @details
#' With `control = datagenControl(method = "mc")` (the default) population
#' moments are computed via the same Monte Carlo engine as `est = "admc"`:
#' \deqn{E_t = \bar{f}_s(\hat\theta_s, \eta_i, t)}
#' \deqn{V_{ts} = \widehat{\mathrm{Cov}}_\eta[f_{s,t}, f_{s,s'}] + \Sigma_s}
#' where \eqn{f_s} and \eqn{\hat\theta_s} are the model and initial estimates
#' from the `ini()` block of study \eqn{s}, the sample covariance uses the
#' ML denominator `n_sim`, and \eqn{\Sigma_s} is diagonal with entries
#' determined by that study model's residual error type (additive, proportional,
#' or log-normal).
#'
#' With `method = "fo"` the moments are instead the deterministic First-Order
#' expansion used by `est = "adfo"`:
#' \deqn{E = f_s(\hat\theta_s, 0)}
#' \deqn{V = J \Omega_s J^\top + \Sigma_s, \quad J_{tj} = \partial f_{s,t}/\partial \eta_j |_{\eta = 0}}
#' with the Jacobian \eqn{J} obtained from the sensitivity model (or finite
#' differences if that is unavailable). This is the natural choice for design
#' evaluation and optimal design: the moments are fast and reproducible, and
#' because the data-generating and data-analytic models coincide, the FO Hessian
#' of the log-likelihood (the expected information matrix) is evaluated at the
#' true maximum rather than at a point that is not an MLE of the generated data.
#' Note `est = "adfo"` always adds \eqn{\Sigma} to its predicted covariance, so
#' for a consistent FIM keep the residual error in the generating model; omit it
#' only when residual-free (IIV-only) moments are genuinely what you want.
#'
#' With `method = "gh"` the moments are computed by deterministic
#' Gauss-Hermite quadrature over the random-effects prior \eqn{\eta \sim N(0, \Omega)}:
#' \deqn{E = \sum_q w_q f(\hat\theta, \eta_q), \quad V = \sum_q w_q (f_q - E)(f_q - E)^\top + \Sigma}
#' where \eqn{(\eta_q, w_q)} are the Cholesky-scaled tensor-product GH nodes and
#' weights. Unlike FO this is unbiased at any IIV magnitude; unlike MC the result
#' is noise-free and exactly reproducible. Matching the moments of `est = "adgh"`
#' makes `method = "gh"` the natural choice for optimal design with that estimator.
#'
#' Models are compiled and cached on first use (keyed by model expression
#' digest), so repeated calls or multiple studies sharing the same model incur
#' only a single compilation.
#'
#' @seealso [datagenControl()], [admControl()]
#' @examples
#' \donttest{
#' library(rxode2)
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
#' study_data <- datagen(
#'   studies = list(
#'     study1 = list(times = c(1, 2, 4, 8, 12, 24),
#'                   ev = rxode2::et(amt = 100), n = 200L)
#'   ),
#'   model   = pk_model,
#'   control = datagenControl(n_sim = 2000L)
#' )
#'
#' # E and V plug directly into admControl(studies = ...)
#' round(study_data$study1$E, 2)
#' }
#' @export
datagen <- function(studies, model = NULL, control = datagenControl()) {
  checkmate::assertList(studies, min.len = 1L)
  if (!inherits(control, "datagenControl"))
    stop("`control` must be created via `datagenControl()`", call. = FALSE)

  # Free rxode2 model registry on exit so repeated runs stay bounded.
  on.exit({ gc(FALSE); rxode2::rxUnloadAll() }, add = TRUE)

  # Ensure studies are named
  study_names <- names(studies) %||% paste0("study", seq_along(studies))

  # Expand stratified studies into per-stratum studies before downstream processing.
  .ex <- .admExpandStrata(studies, study_names, model)
  studies <- .ex$studies; study_names <- .ex$names

  # Validate study specs and resolve per-study model
  study_models <- vector("list", length(studies))
  for (i in seq_along(studies)) {
    nm <- study_names[[i]]
    s  <- studies[[i]]
    # Use `[[` to avoid partial matching on lists.
    m  <- s[["model"]] %||% model
    if (is.null(m))
      stop(sprintf(
        "Study '%s' has no `model` and no top-level default was supplied.", nm),
        call. = FALSE)
    # Accepts rxUi or function (rxode2() is idempotent on rxUi).
    if (!is.function(m) && !inherits(m, "rxUi"))
      stop(sprintf("Study '%s': `model` must be a function or an rxUi.", nm),
           call. = FALSE)
    if (!is.null(s$covariate))
      stop("datagen(): `covariate` (node-quadrature generation) was removed. ",
           "Give the study a `cov_dist` instead -- ONE aggregate (E, V), ",
           "marginal over the covariate distribution, which is what a ",
           "publication reports. For summaries BY covariate stratum, generate ",
           "one ordinary study per stratum with its own `cov` and `n`.",
           call. = FALSE)
    # Multi-output combined with covariate quadrature is unsupported.
    if (!is.null(s$observations)) {
      if (!is.list(s$observations) || length(s$observations) == 0L)
        stop(sprintf("Study '%s': `observations` must be a non-empty list.", nm),
             call. = FALSE)
      for (k in seq_along(s$observations)) {
        o <- s$observations[[k]]
        if (is.null(o$times %||% s$times))
          stop(sprintf("Study '%s' observation %d is missing `times`.", nm, k),
               call. = FALSE)
        if (is.null(o$ev %||% s$ev))
          stop(sprintf("Study '%s' observation %d is missing `ev`.", nm, k),
               call. = FALSE)
      }
    } else {
      if (is.null(s$times))
        stop(sprintf("Study '%s' is missing `times`.", nm), call. = FALSE)
      if (is.null(s$ev))
        stop(sprintf("Study '%s' is missing `ev`.", nm), call. = FALSE)
    }
    study_models[[i]] <- m
  }

  # --- simulation loop ---
  # Models are cached by digest in .admLoadModel().

  if (control$sampling %in% c("rnorm", "lhs")) set.seed(control$seed)

  results <- vector("list", length(studies))
  for (i in seq_along(studies)) {
    s   <- studies[[i]]
    mdl <- study_models[[i]]
    nm  <- study_names[[i]]

    # Parse this study's model
    ui      <- rxode2::rxode2(mdl)
    pinfo   <- .admParseIniDf(ui$iniDf, ui)
    # Forward residual and covariate quadrature node counts via pinfo.
    pinfo$resid_nodes <- control$resid_nodes %||% .ADM_TBS_NODES
    pinfo$cov_nodes   <- control$cov_nodes %||% 7L
    out_var <- .admOutputVar(ui)
    # Model-implied moments lack parameter uncertainty needed for sampling laws.
    pars    <- .admUnpack(.admBuildOptVec(pinfo)$p0, pinfo)

    # Mixed continuous/count endpoints are unsupported; see .admCheckMixedEndpoints.
    .admCheckMixedEndpoints(ui)

    # FO cannot evaluate beta precision phi = b1 + b2 at eta = 0.
    if (control$method == "fo" && !is.null(.admBetaPair(ui)))
      stop("datagen(method = 'fo') does not support a beta() endpoint: the beta ",
           "precision is derived from the solved shapes, which the FO ",
           "linearisation has no path to. Use method = 'mc' or 'gh'.",
           call. = FALSE)

    # Ordinal endpoints require cross-category joint covariance, unsupported here.
    if (any(as.character(tryCatch(ui$predDf$distribution,
                                  error = function(e) character(0))) %in%
            c("ordinal", "dordinal")))
      stop("datagen() does not support an ordinal endpoint: its categories form ",
           "ONE joint\n  observation whose covariance carries the -p_j*p_k term ",
           "between categories at the\n  same time, and datagen() derives each ",
           "observed output separately. Build the\n  study from simulated ",
           "category counts instead, one observation block per category.",
           call. = FALSE)

    # Load sens model before .admLoadModel() to preserve cached compilation order.
    sensModel <- if (control$method == "fo" && pinfo$n_eta > 0L) {
      sm <- tryCatch(.admLoadSensModel(ui), error = function(e) NULL)
      if (is.null(sm))
        warning(sprintf(
          "datagen(method = 'fo'): sensitivity model unavailable for study '%s'; using finite differences for the Jacobian.",
          nm), call. = FALSE)
      sm
    } else NULL

    rxMod   <- .admLoadModel(ui)

    # Resolve observed compartments (observations list or single implicit output).
    obs_specs <- if (!is.null(s$observations)) {
      onm <- names(s$observations) %||% paste0("obs", seq_along(s$observations))
      lapply(seq_along(s$observations), function(k) {
        o <- s$observations[[k]]
        list(name   = onm[k],
             output = o$output %||% out_var,
             times  = o$times  %||% s$times,
             ev     = o$ev     %||% s$ev,
             n      = o$n      %||% s$n)
      })
    } else {
      list(list(name = NULL, output = out_var, times = s$times,
                ev = s$ev, n = s$n))
    }
    # Several observed outputs -> tag each observation's records with its output
    # compartment (nlmixr2's multi-endpoint simulation model routes by cmt).
    is_multi <- length(unique(vapply(obs_specs, `[[`, character(1),
                                     "output"))) > 1L

    grid <- if (control$method == "gh")
      .adghNodeGrid(control$n_nodes, pinfo$n_eta) else NULL

    # Moments (mu, V) for one observed compartment via the chosen method.
    # A study carrying `cov_dist` is generated MARGINAL over that distribution:
    # each simulated subject gets its own covariate value, exactly as the
    # estimator's general path does, so datagen() and the fit integrate the
    # covariate identically rather than by two constructions that could drift.
    # This is the ADM idiom -- a published model plus a study DESIGN (its dosing,
    # its sampling times, its population) produces that study's aggregate data.
    cov_rows_of <- function(n) {
      if (is.null(s[["cov_dist"]])) return(NULL)
      # `fo` linearises in the random effects around a single solve and has no
      # covariate integral at all. `gh` does: .adghGrid crosses the covariate
      # grid with the eta grid, which is the same construction the estimator
      # uses, so it needs no samples and adds no Monte Carlo noise to data that
      # is supposed to BE the reference.
      if (identical(control$method, "fo"))
        stop(sprintf(paste("Study '%s': `cov_dist` needs datagenControl(method =",
                           "\"mc\") or \"gh\"; the fo moment path integrates over",
                           "the random effects only."), nm), call. = FALSE)
      if (!identical(control$method, "mc")) return(NULL)
      .admCovRowsFor(s[["cov_dist"]], n, pinfo$n_eta)
    }
    # Fall back to fixed covariate value if cov_dist is absent.
    cov_ref_of <- function() {
      if (!is.null(s[["cov"]])) return(s[["cov"]])
      if (is.null(s[["cov_dist"]])) return(NULL)
      # Exclude covariance metadata before computing means.
      .cd <- s[["cov_dist"]][setdiff(names(s[["cov_dist"]]),
                                     .ADM_COV_META)]
      stats::setNames(lapply(.cd, .admCovMeanOf), names(.cd))
    }

    compute_moments <- function(spec) {
      ov  <- spec$output
      n_t <- length(spec$times)
      arr <- .admResidRows(pinfo, ov, pars$sigma_var, n_t)
      evf <- if (is_multi) spec$ev |> rxode2::et(spec$times, cmt = ov)
             else          spec$ev |> rxode2::et(spec$times)
      # Pass out_pair for beta endpoints (phi = b1 + b2) and fixed cov values.
      study_tmp <- list(ev_full = evf, times = spec$times,
                        out_pair = .admBetaPair(ui),
                        cov = cov_ref_of(),
                        cov_rows = cov_rows_of(control$n_sim))
      # GH integrates cov_dist over its quadrature grid, not at the mean.
      if (control$method == "gh")
        study_tmp$cov_dist <- s[["cov_dist"]]

      if (control$method == "gh") {
        m <- .adghMoments(pars, pinfo, study_tmp, rxMod, ov, grid, control$cores)
        list(mu = m$E, V = m$V, cp_mat = NULL)
      } else if (control$method == "fo") {
        params_mat <- .admMakeParamsList(1L, pinfo, 1L)[[1L]]
        mj <- .adfoGetMuJ(pars, pinfo, study_tmp, sensModel, rxMod, ov,
                          params_mat, control$cores)
        vp <- .adfoVpred(mj$mu, mj$J, pars$L, arr, n_t, pinfo$n_eta,
                           study_tmp$times)
        list(mu = vp$mu_sigma, V = vp$V, cp_mat = NULL)
      } else {
        z_list      <- .admMakeZ(control$n_sim, pinfo, 1L, control$sampling)
        params_list <- .admMakeParamsList(control$n_sim, pinfo, 1L)
        if (pinfo$n_eta > 0L) {
          z <- z_list[[1L]]
          if (!is.matrix(z)) z <- matrix(z, ncol = 1L)  # sobol dim=1 edge case
          eta_mat <- z %*% t(pars$L)
          colnames(eta_mat) <- pinfo$eta_col_names
        } else {
          eta_mat <- matrix(0, control$n_sim, 0L)
        }
        # Integrate at default solver tolerances without truncation.
        cp_mat <- .admSimulate(rxMod, pars$struct, pinfo$sigma_names, eta_mat,
                               study_tmp, ov, params_list[[1L]], control$cores)
        # beta precision (SOLVED) rides back on cp_mat; fold it into the row array
        arr  <- .admUnitResidRows(pinfo, ov, pars$sigma_var, n_t,
                                  phi = attr(cp_mat, "phi"))
        mu   <- colMeans(cp_mat)
        cp_c <- sweep(cp_mat, 2L, mu)
        V    <- crossprod(cp_c) / control$n_sim
        # Apply residual error structure; TBS endpoints integrate per draw.
        .ex <- .admResidNodeMomentsTBS(cp_mat, rep(1, nrow(cp_mat)), arr,
                                       study_tmp$times)
        if (!is.null(.ex))
          return(list(mu = .ex$E, V = .ex$V, cp_mat = cp_mat))
        ap   <- .admResidApply(mu, diag(V), arr, study_tmp$times, V)
        list(mu = ap$mu, V = .admApplyResidTail(V, ap), cp_mat = cp_mat)
      }
    }

    one_result <- function(spec) {
      m     <- compute_moments(spec)
      t_lbl <- as.character(spec$times)
      mu <- m$mu; V <- m$V
      names(mu) <- t_lbl; dimnames(V) <- list(t_lbl, t_lbl)
      r <- list(E = mu, V = V, n = spec$n %||% NA_integer_,
                times = spec$times, ev = spec$ev)
      # Retain cov_dist so downstream estimators marginalise over the same population.
      if (!is.null(s[["cov_dist"]])) { r$cov_dist <- s[["cov_dist"]] }
      # Datagen uses ML denominator n_sim for V.
      r$v_denom <- "ml"
      if (!is.null(cov_ref_of()))     { r$cov      <- cov_ref_of() }
      # Mark as model source (.adm_src) so standard errors can be refused.
      r$.adm_src <- TRUE
      # Carry stratum resolution for downstream consistency checks (e.g. anova).
      if (!is.null(s[[".adm_strata_nodes"]]))
        r$.adm_strata_nodes <- s[[".adm_strata_nodes"]]
      if (!is.null(spec$output)) r$output <- spec$output
      if (control$return_samples && !is.null(m$cp_mat)) r$samples <- m$cp_mat
      r
    }

    if (!is.null(s$observations)) {
      obs_res <- lapply(obs_specs, one_result)
      names(obs_res) <- vapply(obs_specs, `[[`, character(1), "name")
      results[[i]] <- list(observations = obs_res, n = s$n %||% NA_integer_,
                           .adm_src = TRUE)
    } else {
      results[[i]] <- one_result(obs_specs[[1L]])
    }
  }

  # setNames preserves duplicates so name collisions can be diagnosed below.
  if (anyDuplicated(unlist(study_names)))
    stop("admixr2: datagen() produced duplicate study name(s) ",
         paste(sQuote(unique(unlist(study_names)[
           duplicated(unlist(study_names))])), collapse = ", "),
         ". Studies are matched by name downstream, so a duplicate would ",
         "silently drop all but one. Rename the study, or the banded source ",
         "whose strata collide with it.", call. = FALSE)
  stats::setNames(results, unlist(study_names))
}

# Strips .adm_src so simulated data can be treated as an observed study with SE.
.admDatagenSim <- function(...) {
  lapply(datagen(...), function(u) {
    u[[".adm_src"]] <- NULL
    if (is.list(u[["observations"]]))
      u$observations <- lapply(u$observations, function(o) {
        o[[".adm_src"]] <- NULL
        o
      })
    u
  })
}

# --- Who may be given a standard error --------------------------------------

.admStudyHasModelSource <- function(study) {
  if (!is.list(study)) return(FALSE)
  if (inherits(study, "admStudy")) return(!is.null(study$ui))
  if (isTRUE(study[[".adm_src"]])) return(TRUE)
  observations <- study[["observations"]]
  is.list(observations) && any(vapply(
    observations, function(x) is.list(x) && isTRUE(x[[".adm_src"]]), logical(1)
  ))
}

.admHasModelSource <- function(studies) {
  if (inherits(studies, "admStudies")) studies <- unclass(studies)
  if (!is.list(studies) || !length(studies)) return(FALSE)
  any(vapply(studies, .admStudyHasModelSource, logical(1)))
}

# Warn when multiple sources are pooled but sample size `n` is missing or invalid.
.admWarnSourceWeight <- function(studies) {
  # Require positive sample size when pooling multiple sources.
  if (!is.list(studies) || !length(studies)) return(invisible(NULL))
  # Guard against malformed study elements.
  .is_src <- vapply(studies, .admStudyHasModelSource, logical(1))
  if (length(studies) > 1L && any(.is_src)) {
    study_names <- names(studies) %||% paste0("study", seq_along(studies))
    bad_n <- study_names[.is_src][vapply(studies[.is_src], function(s) {
      nn <- as.numeric(s[["n"]] %||% NA_real_)
      !is.finite(nn) || nn <= 0 }, logical(1))]
    if (length(bad_n))
      warning("admixr2: model source", if (length(bad_n) > 1L) "s " else " ",
              paste(sQuote(bad_n), collapse = ", "), " ",
              if (length(bad_n) > 1L) "have" else "has",
              " no usable `n`, and this fit combines several sources. Set `n` ",
              "to the true sample size of the dataset used to develop the ",
              "source model; it determines that study's contribution to the ",
              "pooled point estimate.", call. = FALSE)
  }
  invisible(NULL)
}

.admResolveCovMethod <- function(covMethod, studies, explicit) {
  .admWarnSourceWeight(studies)
  if (!.admHasModelSource(studies) || identical(covMethod, "none"))
    return(covMethod)
  if (isTRUE(explicit))
    stop("admixr2: covMethod = ", dQuote(covMethod), " was requested, but a ",
         "study contributes model-implied moments rather than observed sample ",
         "summaries. `n` is the source study's true sample size, but the source ",
         "parameter uncertainty needed for a sampling law is unavailable. ",
         "Pass covMethod = \"none\".", call. = FALSE)
  message("admixr2: a study contributes as a published MODEL, so covMethod ",
          "has been set to \"none\". `n` remains the source study's sample ",
          "size, but its parameter uncertainty is unavailable, so no sampling ",
          "law can be built for a standard error.")
  "none"
}

