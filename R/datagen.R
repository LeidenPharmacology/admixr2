#' Control parameters for [datagen()]
#'
#' @param n_sim Number of Monte Carlo samples used to approximate population
#'   moments.
#' @param sampling Quasi-random sampling method: `"sobol"` (default),
#'   `"halton"`, `"torus"`, `"lhs"`, or `"rnorm"`.
#' @param seed Integer seed.  Applied before stochastic methods
#'   (`"rnorm"`, `"lhs"`).
#' @param cores Number of `rxSolve` threads.
#' @param add_residual_error Add residual-error variance to the diagonal of
#'   `V` (`TRUE` by default), matching the admixr2 NLL convention.
#' @param return_samples Include the raw `n_sim x length(times)`
#'   prediction matrix as `$samples` in each study's output.
#'
#' @return A list of class `"datagenControl"`.
#' @seealso [datagen()]
#' @examples
#' ctrl <- datagenControl(n_sim = 2000L)
#' ctrl$sampling  # "sobol"
#' @export
datagenControl <- function(
  n_sim              = 5000L,
  sampling           = c("sobol", "halton", "torus", "lhs", "rnorm"),
  seed               = 12345L,
  cores              = 1L,
  add_residual_error = TRUE,
  return_samples     = FALSE
) {
  sampling <- match.arg(sampling)
  checkmate::assertIntegerish(n_sim, lower = 1L, len = 1L)
  checkmate::assertIntegerish(seed,              len = 1L)
  checkmate::assertIntegerish(cores, lower = 1L, len = 1L)
  checkmate::assertFlag(add_residual_error)
  checkmate::assertFlag(return_samples)
  structure(
    list(
      n_sim              = as.integer(n_sim),
      sampling           = sampling,
      seed               = as.integer(seed),
      cores              = as.integer(cores),
      add_residual_error = add_residual_error,
      return_samples     = return_samples
    ),
    class = "datagenControl"
  )
}


#' Generate aggregate study data from (possibly different) pharmacometric models
#'
#' Simulates population mean vectors (`E`) and covariance matrices
#' (`V`) for each study using Monte Carlo integration over the IIV
#' distribution.  Each study may specify its own PK/PD model (as would be the
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
#'   }
#' @param model Optional default model function used for any study that does not
#'   supply its own `model` element.  At least one of `model` or each
#'   study's `model` must be non-`NULL`.
#' @param control A [datagenControl()] object.
#'
#' @return A named list with one element per study.  Each element contains:
#'   \describe{
#'     \item{`E`}{Population mean vector at `times`.}
#'     \item{`V`}{Population covariance matrix
#'       (`length(times)` x `length(times)`, ML denominator
#'       `n_sim`).  Residual error is added to the diagonal when
#'       `control$add_residual_error = TRUE`.}
#'     \item{`n`}{Sample size (`NA_integer_` if not supplied).}
#'     \item{`times`}{Observation times.}
#'     \item{`ev`}{Dosing event table.}
#'     \item{`samples`}{Raw `n_sim x length(times)` prediction matrix
#'       (only when `control$return_samples = TRUE`).}
#'   }
#'
#' @details
#' Population moments are computed via the same Monte Carlo engine as
#' `est = "admc"`:
#' \deqn{E_t = \bar{f}_s(\hat\theta_s, \eta_i, t)}
#' \deqn{V_{ts} = \widehat{\mathrm{Cov}}_\eta[f_{s,t}, f_{s,s'}] + \Sigma_s}
#' where \eqn{f_s} and \eqn{\hat\theta_s} are the model and initial estimates
#' from the `ini()` block of study \eqn{s}, the sample covariance uses the
#' ML denominator `n_sim`, and \eqn{\Sigma_s} is diagonal with entries
#' determined by that study model's residual error type (additive, proportional,
#' or log-normal).
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

  # Ensure studies are named
  study_names <- names(studies) %||% paste0("study", seq_along(studies))

  # Validate study specs and resolve per-study model
  study_models <- vector("list", length(studies))
  for (i in seq_along(studies)) {
    nm <- study_names[[i]]
    s  <- studies[[i]]
    m  <- s$model %||% model
    if (is.null(m))
      stop(sprintf(
        "Study '%s' has no `model` and no top-level default was supplied.", nm),
        call. = FALSE)
    if (!is.function(m))
      stop(sprintf("Study '%s': `model` must be a function.", nm), call. = FALSE)
    if (is.null(s$times))
      stop(sprintf("Study '%s' is missing `times`.", nm), call. = FALSE)
    if (is.null(s$ev))
      stop(sprintf("Study '%s' is missing `ev`.", nm), call. = FALSE)
    study_models[[i]] <- m
  }

  # --- simulation loop ---
  # .admLoadModel() is cached by model digest, so identical models across
  # studies compile only once.

  if (control$sampling %in% c("rnorm", "lhs")) set.seed(control$seed)

  results <- vector("list", length(studies))
  for (i in seq_along(studies)) {
    s   <- studies[[i]]
    mdl <- study_models[[i]]

    # Parse this study's model
    ui      <- rxode2::rxode2(mdl)
    pinfo   <- .admParseIniDf(ui$iniDf, ui)
    out_var <- .admOutputVar(ui)
    pars    <- .admUnpack(.admBuildOptVec(pinfo)$p0, pinfo)
    rxMod   <- .admLoadModel(ui)

    # Quasi-random draws (n_sim x n_eta) for this study
    z_list      <- .admMakeZ(control$n_sim, pinfo, 1L, control$sampling)
    params_list <- .admMakeParamsList(control$n_sim, pinfo, 1L)

    # Correlated random effects: eta_mat = z L'
    if (pinfo$n_eta > 0L) {
      z <- z_list[[1L]]
      if (!is.matrix(z)) z <- matrix(z, ncol = 1L)  # sobol dim=1 edge case
      eta_mat <- z %*% t(pars$L)
      colnames(eta_mat) <- pinfo$eta_col_names
    } else {
      eta_mat <- matrix(0, control$n_sim, 0L)
    }

    # Forward simulation — n_sim x n_times matrix
    study_tmp <- list(ev_full = s$ev |> rxode2::et(s$times), times = s$times)
    cp_mat <- .admSimulate(
      rxMod, pars$struct, pinfo$sigma_names,
      eta_mat, study_tmp, out_var, params_list[[1L]], control$cores
    )

    # Population mean and IIV covariance (ML denominator, matches nll_cov_cpp)
    mu   <- colMeans(cp_mat)
    cp_c <- sweep(cp_mat, 2L, mu)
    V    <- crossprod(cp_c) / control$n_sim

    # Diagonal residual error (mirrors nll_cov_cpp sigma_type dispatch)
    if (control$add_residual_error && length(pars$sigma_var) > 0L) {
      sv <- unname(pars$sigma_var[[1L]])
      if (isTRUE(pinfo$sigma_is_lnorm[[1L]])) {
        mu <- mu * exp(sv / 2)                          # lnorm mean correction
        diag(V) <- diag(V) + mu^2 * (exp(sv) - 1)
      } else if (isTRUE(pinfo$sigma_is_prop[[1L]])) {
        diag(V) <- diag(V) + sv * mu^2
      } else {
        diag(V) <- diag(V) + sv
      }
    }

    t_lbl     <- as.character(s$times)
    names(mu) <- t_lbl
    dimnames(V) <- list(t_lbl, t_lbl)

    result_i <- list(
      E     = mu,
      V     = V,
      n     = s$n %||% NA_integer_,
      times = s$times,
      ev    = s$ev
    )
    if (control$return_samples) result_i$samples <- cp_mat
    results[[i]] <- result_i
  }

  setNames(results, study_names)
}
