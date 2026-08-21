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
#' @param cov_nodes Gauss-Hermite nodes per covariate when `method = "gh"`
#'   integrates a study's `cov_dist` (default 7). Total covariate points are
#'   `cov_nodes^p` for `p` covariates. Ignored by `"mc"`, which draws a covariate
#'   value per simulated subject instead, and by `"fo"`, which cannot integrate a
#'   covariate at all.
#'   `"halton"`, `"torus"`, `"lhs"`, or `"rnorm"`. Ignored when `method = "fo"`
#'   or `"gh"`.
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
  # LAST on purpose: inserting an argument mid-signature silently rebinds every
  # positional call (datagenControl("mc", 2000L, 7L) used to set n_nodes = 7).
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
#'     \item{`model_cov`}{(Optional) the source model's OWN parameter
#'       covariance, as that model's analysis reported it. Supplying it is what
#'       makes a standard error meaningful for a study generated from a
#'       published model: such a study is not a sample, its `E`/`V` are exact
#'       functions of the source's parameters, so the only uncertainty in the
#'       chain is the source's own and the fit propagates it by the delta
#'       method. Without it `n` is read as a sample size instead, and the
#'       reported SE then falls as `1/sqrt(n)` --- a factor you choose by
#'       typing a number, not a property of the evidence.
#'       A **named** matrix, whose dimnames are the source model's `ini()`
#'       parameter names. The names are load-bearing twice over: they say which
#'       parameter each row is, and they fix the SCALE, which is the scale that
#'       model's `ini()` uses --- log for a theta written `tcl <- log(5)`, an SD
#'       for a residual, a variance for an omega. That is also the scale
#'       nlmixr2 prints estimates on.
#'       It must cover every parameter the source ESTIMATES. A parameter left
#'       out contributes no uncertainty at all, which asserts the source knew it
#'       exactly and reports an SE that is too small, so an incomplete matrix
#'       warns and no standard error is reported. A `fix()`ed parameter is an
#'       assertion and must NOT appear. A diagonal --- SEs from a paper's
#'       \%RSE, with no correlations --- is a valid fallback: it is exact at the
#'       source's own covariate reference and degrades only as you extrapolate
#'       away from it.}
#'     \item{`cov_dist`}{(Optional) the covariate distribution this study's
#'       subjects span --- see [covDraw()] for the grammar. The generated
#'       `E`/`V` are MARGINAL over it, which is what a publication reports.
#'       Needs `datagenControl(method = "mc")`.}
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

  # Per-study model loading populates rxode2's global model registry; free it with
  # rxode2's own idiom on exit so repeated datagen() runs stay bounded.
  on.exit({ gc(FALSE); rxode2::rxUnloadAll() }, add = TRUE)

  # Ensure studies are named
  study_names <- names(studies) %||% paste0("study", seq_along(studies))

  # A study declaring `stratify` is expanded HERE, into one ordinary study per
  # covariate stratum, before anything else looks at it. Everything downstream
  # -- generation, and then the fit -- then sees plain studies and needs no
  # knowledge of where they came from.
  .ex <- .admExpandStrata(studies, study_names, model)
  studies <- .ex$studies; study_names <- .ex$names

  # Validate study specs and resolve per-study model
  study_models <- vector("list", length(studies))
  for (i in seq_along(studies)) {
    nm <- study_names[[i]]
    s  <- studies[[i]]
    # `[[ ]]`, NOT `$`: a study also carries `model_cov`, and `$` PARTIAL-MATCHES
    # on lists -- `s$model` silently returned the covariance MATRIX the moment
    # that field was added, and every study then failed as "must be a function".
    m  <- s[["model"]] %||% model
    if (is.null(m))
      stop(sprintf(
        "Study '%s' has no `model` and no top-level default was supplied.", nm),
        call. = FALSE)
    # An rxUi is accepted alongside a function because rxode2::rxode2() is
    # idempotent on one, and the model-source Jacobian re-generates the blocks
    # at PERTURBED parameter values -- which is a modified ui, not a function.
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
    # Covariate marginalisation expands ONE study into one sub-study per
    # quadrature node, each carrying a single E/V. A multi-output study is
    # already a list of per-output blocks, and nothing downstream defines what
    # the product of the two should be -- so refuse it rather than emit a shape
    # no estimator reads.
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
  # .admLoadModel() is cached by model digest, so identical models across
  # studies compile only once.

  if (control$sampling %in% c("rnorm", "lhs")) set.seed(control$seed)

  results <- vector("list", length(studies))
  for (i in seq_along(studies)) {
    s   <- studies[[i]]
    mdl <- study_models[[i]]
    nm  <- study_names[[i]]

    # Parse this study's model
    ui      <- rxode2::rxode2(mdl)
    pinfo   <- .admParseIniDf(ui$iniDf, ui)
    # Residual-quadrature nodes travel on pinfo -> arr -> .admResidApply(), the
    # same route the estimators use, so a study generated here and the fit that
    # consumes it integrate the residual identically.
    pinfo$resid_nodes <- control$resid_nodes %||% .ADM_TBS_NODES
    # Reaches .admCovGrid through .adghGrid; without it the control argument is
    # inert and the grid silently uses its own default.
    pinfo$cov_nodes   <- control$cov_nodes %||% 7L
    out_var <- .admOutputVar(ui)
    .src_cov  <- .admSrcCov(s[["model_cov"]] %||% attr(model, "model_cov"), ui, nm)
    .src_prov <- list(id    = s[[".adm_src_id"]] %||% nm,
                      model = mdl,
                      theta = .admSrcTheta(ui),
                      cov   = .src_cov$cov,
                      par   = .src_cov$par,
                      missing = .src_cov$missing,
                      control = control)
    pars    <- .admUnpack(.admBuildOptVec(pinfo)$p0, pinfo)

    # A model mixing a continuous endpoint with a COUNT one is refused here for
    # exactly the reason the estimators refuse it (.admCheckMixedEndpoints is
    # called by nlmixr2Est.adfo/.adgh/.admc): a count endpoint's output is the
    # DISTRIBUTION'S ARGUMENT (`y ~ pois(lam)` is read through `lam`), a model
    # variable rather than a compartment, so the `cmt = ov` tagging this function
    # applies below when `is_multi` matches no observation record. Generating such
    # a model either errored from inside .admSimulate() with no mention of the
    # endpoint, or recycled into an E/V the user then fed straight back into a fit
    # -- while the SAME model fitted directly was refused with an explanation.
    .admCheckMixedEndpoints(ui)

    # method = "fo" has no path to a beta endpoint's precision: .adfoVpred builds
    # V from J Omega J' + Sigma at eta = 0 and never sees the solved b1 + b2, so
    # it would emit a V whose diagonal is NA. This is the same refusal
    # nlmixr2Est.adfo() makes for the same reason -- said here rather than left to
    # produce NAs, because datagen() has no fit to fail afterwards.
    if (control$method == "fo" && !is.null(.admBetaPair(ui)))
      stop("datagen(method = 'fo') does not support a beta() endpoint: the beta ",
           "precision is derived from the solved shapes, which the FO ",
           "linearisation has no path to. Use method = 'mc' or 'gh'.",
           call. = FALSE)

    # An ordinal endpoint is a JOINT observation: its categories are one stacked
    # vector whose covariance carries the -p_j*p_k term between categories at the
    # same time. datagen() computes moments one observation spec at a time, each
    # with its own `arr` and its own rows, so that cross-category block cannot be
    # formed here at all -- the study it emitted would be scored against a
    # covariance missing exactly the multinomial structure an ordinal model exists
    # to capture, and .admCheckOrdinal() would then refuse it on the way back in.
    # Refuse at the point of generation instead of emitting something unusable.
    if (any(as.character(tryCatch(ui$predDf$distribution,
                                  error = function(e) character(0))) %in%
            c("ordinal", "dordinal")))
      stop("datagen() does not support an ordinal endpoint: its categories form ",
           "ONE joint\n  observation whose covariance carries the -p_j*p_k term ",
           "between categories at the\n  same time, and datagen() derives each ",
           "observed output separately. Build the\n  study from simulated ",
           "category counts instead, one observation block per category.",
           call. = FALSE)

    # FO needs the sensitivity model for the Jacobian df/d(eta)|_0. Load it
    # before .admLoadModel() to respect the compilation-ordering invariant
    # (.admLoadModel() poisons the cached inner model on the first-compile path).
    sensModel <- if (control$method == "fo" && pinfo$n_eta > 0L) {
      sm <- tryCatch(.admLoadSensModel(ui), error = function(e) NULL)
      if (is.null(sm))
        warning(sprintf(
          "datagen(method = 'fo'): sensitivity model unavailable for study '%s'; using finite differences for the Jacobian.",
          nm), call. = FALSE)
      sm
    } else NULL

    rxMod   <- .admLoadModel(ui)

    # Resolve the observed compartments for this study. A study may carry an
    # `observations` list (one entry per observed output, each with its own
    # output/times/ev/n); a legacy study describes a single implicit observation.
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
    # A study may fix a covariate VALUE (`cov`) instead of, or as well as, a
    # distribution. Deriving only from `cov_dist` left a study with a fixed
    # covariate unable to solve at all -- the model reads the covariate and
    # nothing supplies it.
    cov_ref_of <- function() {
      if (!is.null(s[["cov"]])) return(s[["cov"]])
      if (is.null(s[["cov_dist"]])) return(NULL)
      # `rho`, `Sigma` and `joint` are metadata and a sampler, not covariate
      # specs -- .admCovMeanOf() has nothing to compute from a function.
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
      # A beta endpoint's prediction is DERIVED from two solved columns and its
      # precision phi = b1 + b2 comes back with them -- the same pair the
      # estimators put on every study. Without it .admOutputVar() resolves to the
      # first shape parameter, so datagen() returned an `E` that was a shape (an
      # arbitrary positive number, not a probability) and a `V` whose diagonal was
      # entirely NA, with no error and no warning.
      # `cov` rides on the study exactly as it does on the fit path, so the
      # solve paths pick it up through the same channel.
      study_tmp <- list(ev_full = evf, times = spec$times,
                        out_pair = .admBetaPair(ui),
                        cov = cov_ref_of(),
                        cov_rows = cov_rows_of(control$n_sim))
      # `gh` integrates the covariate on its own grid, so it needs the
      # DISTRIBUTION, not just the reference value. Passing only `cov` left it
      # solving at the covariate mean -- the ecological plug-in, generating data
      # for a population that does not exist. Measured against the mc path on a
      # lognormal covariate: 2.1e-02 on the mean and 2.9e-01 on the covariance.
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
        # No sigdig: datagen() generates the TRUTH, so it integrates at rxode2's
        # own tolerances regardless of what any later fit asks for. That matches
        # the default fit exactly (the estimator controls default sigdig = NULL);
        # a fit that opts into a looser sigdig is then measured against a
        # reference tighter than itself, which is the right way round -- the
        # alternative attributes the solver gap to the estimator.
        cp_mat <- .admSimulate(rxMod, pars$struct, pinfo$sigma_names, eta_mat,
                               study_tmp, ov, params_list[[1L]], control$cores)
        # beta precision (SOLVED) rides back on cp_mat; fold it into the row array
        arr  <- .admUnitResidRows(pinfo, ov, pars$sigma_var, n_t,
                                  phi = attr(cp_mat, "phi"))
        # This output's residual error only. `times` + the structural covariance are
        # needed by the off-diagonal forms (ar, ordinal); without them datagen()
        # emitted a V that contradicted the model it was handed -- and disagreed
        # with its own method = "gh" branch, which went through .adghMoments and
        # did include them.
        m <- .admResidSampleMoments(cp_mat, arr, study_tmp$times)
        list(mu = m$mu, V = m$V, cp_mat = cp_mat)
      }
    }

    one_result <- function(spec) {
      m     <- compute_moments(spec)
      t_lbl <- as.character(spec$times)
      mu <- m$mu; V <- m$V
      names(mu) <- t_lbl; dimnames(V) <- list(t_lbl, t_lbl)
      r <- list(E = mu, V = V, n = spec$n %||% NA_integer_,
                times = spec$times, ev = spec$ev)
      # Carry the covariate distribution onto the result so the generated study
      # is directly fittable: the estimator must marginalise over the same
      # population the data were generated for, and making the caller restate it
      # is a way for the two to disagree.
      if (!is.null(s[["cov_dist"]])) { r$cov_dist <- s[["cov_dist"]] }
      # Self-describing: datagen builds V with the ML denominator, so say so
      # rather than leaving the consumer to rely on the default meaning the same
      # thing. A generated study can then be mixed with a digitised one that
      # declares "unbiased" and both are converted correctly.
      r$v_denom <- "ml"
      if (!is.null(cov_ref_of()))     { r$cov      <- cov_ref_of() }
      # Provenance for the standard error. See .admSrcCov(): a generated block
      # is not a sample, so its contribution to the covariance runs through the
      # SOURCE's own uncertainty rather than through `n`. `id` groups the strata
      # of one banded source, which must count as ONE contribution however fine
      # the banding.
      r$.adm_src <- .src_prov
      # The stratum resolution, carried for the same reason and by the same
      # route. .admExpandStrata() stamps it on the study, but one_result() builds
      # its output from an explicit field list, so it was being DROPPED here --
      # which silently disabled everything downstream that reads it:
      # .admFinaliseFit() never recorded `strataNodes`, so anova()'s refusal to
      # compare two fits built at different resolutions could not fire on any
      # generated study, which is the normal path.
      if (!is.null(s[[".adm_strata_nodes"]]))
        r$.adm_strata_nodes <- s[[".adm_strata_nodes"]]
      if (!is.null(spec$output)) r$output <- spec$output
      if (control$return_samples && !is.null(m$cp_mat)) r$samples <- m$cp_mat
      r
    }

    if (!is.null(s$observations)) {
      obs_res <- lapply(obs_specs, one_result)
      names(obs_res) <- vapply(obs_specs, `[[`, character(1), "name")
      results[[i]] <- list(observations = obs_res, n = s$n %||% NA_integer_)
    } else {
      results[[i]] <- one_result(obs_specs[[1L]])
    }
  }

  out <- list()
  for (i in seq_along(results)) out[[study_names[[i]]]] <- results[[i]]
  out
}

# NOTE: these live at the END of the file deliberately. Inserting them above
# datagenControl() put them BETWEEN that function and its roxygen block, so
# the block bound to the helper instead -- datagenControl stopped being
# exported and an internal was exported and documented in its place. Nothing
# follows here, so nothing can be orphaned.
# =============================================================================
# Model sources: the provenance the standard error needs
# =============================================================================
#
# A study generated from a published MODEL is not a sample. Its (E, V) are exact
# functions of that model's parameters, so the only random thing in the whole
# chain is `theta_src_hat` -- the estimate the source published -- and the
# covariance of our fit is the delta method through it:
#
#     Var(theta_hat) = G C_src G' ,      G = d theta_hat / d theta_src
#
# `n` is NOT a precision statement for such a study. It divides straight out of
# a lone source's estimating equation (measured: the estimate is 0.75000 at
# n = 100, 400, 1600 and 6400, unchanged), and sets only the RELATIVE WEIGHT
# against other sources. Reading `n` as precision is what makes the reported SE
# fall as exactly 1/sqrt(n) -- a factor the analyst chooses by typing a number.
#
# So the covariance is declared HERE, in the datagen block, beside the model it
# belongs to. Nothing has to be restated at fit time, and a generated study
# carries everything its own standard error needs.
#
#   datagen(list(trial1 = list(times = ..., ev = ..., n = 240,
#                              model     = published_mod,
#                              model_cov = C)),      # <- the source's own
#           control = datagenControl(method = "gh"))
#
# SCALE. `C_src` is on the scale of the SOURCE MODEL'S OWN `ini()` block, which
# is the scale that gets perturbed to form the Jacobian: log for a theta written
# `tcl <- log(5)`, natural for one written `bwt <- 0.75`. This is checked rather
# than documented -- the dimnames must be `iniDf` parameter names, so a matrix
# built against the wrong parameterisation is refused instead of silently
# rescaling every reported interval.
.admSrcCov <- function(cov, ui, nm) {
  if (is.null(cov)) return(NULL)
  bad <- function(...) stop("admixr2: study '", nm, "': `model_cov` ", ...,
                            call. = FALSE)
  cov <- as.matrix(cov)
  if (nrow(cov) != ncol(cov)) bad("must be square; got ", nrow(cov), " x ",
                                  ncol(cov), ".")
  rn <- rownames(cov) %||% colnames(cov)
  if (is.null(rn))
    bad("must carry the parameter NAMES as dimnames -- they are what says which ",
        "parameter each row is, and on which scale. Use the names from the ",
        "source model's `ini()` block, e.g. dimnames(C) <- list(c(\"tcl\", ",
        "\"bwt\"), c(\"tcl\", \"bwt\")).")
  if (!is.null(colnames(cov)) && !identical(rownames(cov), colnames(cov)))
    bad("has different row and column names, so it does not describe one ",
        "parameter set.")
  ini <- tryCatch(ui$iniDf, error = function(e) NULL)
  if (is.null(ini)) bad("cannot be checked: the source model would not parse.")
  # ACCEPT A FIT'S COVARIANCE AS IT COMES. The obvious thing to pass is the
  # source fit's own `$cov`, and nlmixr2 names its omega rows on the REPORTING
  # convention -- `om.eta.cl` for the variance of `eta.cl`, `cov.a.b` for an
  # off-diagonal -- while an `ini()` block calls that row `eta.cl`. Refusing
  # `om.` would mean every user renaming a matrix by hand, so translate it.
  # The scale already agrees: both are the variance.
  .om <- grepl("^om[.]", rn)
  if (any(.om)) {
    cand <- sub("^om[.]", "", rn[.om])
    okm  <- cand %in% ini$name[!is.na(ini$neta1) & ini$neta1 == ini$neta2]
    if (any(okm)) rn[which(.om)[okm]] <- cand[okm]
  }
  .cv <- grepl("^cov[.]", rn)
  if (any(.cv))
    bad("names ", paste(sQuote(rn[.cv]), collapse = ", "),
        ", which are OFF-DIAGONAL omega entries. Those are not yet mapped from ",
        "the reporting names onto `ini()` rows; drop them and their rows, ",
        "keeping the `om.` diagonal, or rename them to the names the source ",
        "model's `ini()` uses.")
  dimnames(cov) <- list(rn, rn)
  # Only ESTIMATED parameters carry uncertainty. A fix()ed one is an assertion
  # -- the source claims to know it -- so it contributes no variance, and naming
  # it is a sign the matrix came from somewhere other than that model's fit.
  est <- ini$name[!ini$fix]
  unknown <- setdiff(rn, ini$name)
  if (length(unknown))
    bad("names ", paste(sQuote(unknown), collapse = ", "),
        ", which the source model's `ini()` does not declare. Declared: ",
        paste(sQuote(ini$name), collapse = ", "), ".")
  fixed <- intersect(rn, ini$name[ini$fix])
  if (length(fixed))
    bad("names ", paste(sQuote(fixed), collapse = ", "),
        ", which the source model fix()es. A fixed parameter is an ASSERTION, ",
        "so it carries no uncertainty to propagate -- drop ",
        if (length(fixed) == 1L) "it" else "them", " from `model_cov`.")
  if (!isTRUE(all.equal(unname(cov), unname(t(cov)), tolerance = 1e-8)))
    bad("is not symmetric, so it is not a covariance matrix.")
  ev <- tryCatch(eigen(cov, symmetric = TRUE, only.values = TRUE)$values,
                 error = function(e) NULL)
  if (is.null(ev) || min(ev) < -1e-10 * max(1, max(ev)))
    bad("is not positive semi-definite (smallest eigenvalue ",
        sprintf("%.3g", if (is.null(ev)) NA_real_ else min(ev)),
        "), so it describes no distribution. A matrix rebuilt from published ",
        "SEs and correlations can fail this if the correlations were rounded; ",
        "supplying SEs only, as a DIAGONAL, is a valid fallback.")
  miss <- setdiff(est, rn)
  # AN INCOMPLETE C_src IS NOT A PARTIAL ANSWER. A parameter the source
  # ESTIMATED but reported no covariance for contributes zero to G C_src G',
  # which asserts the source knew it exactly and makes the reported standard
  # error too SMALL -- the dangerous direction, and nothing about the matrix
  # looks wrong. Said HERE, at generation, because a warning raised later from
  # inside CalcCov does not reach the user: the nlmixr2est stack swallows it,
  # which is why the drivers report a missing covariance from their own frame.
  if (length(miss))
    warning("admixr2: study '", nm, "': `model_cov` covers ",
            paste(sQuote(rn), collapse = ", "), " but the source model also ",
            "ESTIMATES ", paste(sQuote(miss), collapse = ", "),
            ". A parameter with no covariance contributes none, which asserts ",
            "the source knew it exactly and reports a standard error that is ",
            "too small, so no standard error will be reported at all.
",
            "  Supply the missing ", if (length(miss) == 1L) "row" else "rows",
            " -- a DIAGONAL entry from the paper's %RSE is a valid fallback -- ",
            "or fix() ", if (length(miss) == 1L) "it" else "them",
            " in the source model if it asserted ",
            if (length(miss) == 1L) "it" else "them", " rather than ",
            "estimating ", if (length(miss) == 1L) "it." else "them.",
            call. = FALSE)
  list(cov = cov, par = rn, missing = miss)
}

# Which parameters does a source model actually ESTIMATE, and at what values?
#
# Read at generation time rather than from the function text, because a model
# routinely reads its values from variables -- `ini({ tcl <- log(CLp) })` is the
# idiom for supplying a published model's numbers -- so the function alone does
# not say what was generated.
.admSrcTheta <- function(ui) {
  ini <- tryCatch(ui$iniDf, error = function(e) NULL)
  if (is.null(ini)) return(NULL)
  keep <- !ini$fix
  stats::setNames(ini$est[keep], ini$name[keep])
}
