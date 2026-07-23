# -- Control object -------------------------------------------------------------

#' Control settings for the ADM estimator
#'
#' Constructs a control object for `est = "admc"`, the Monte Carlo aggregate
#' data modelling estimator.
#'
#' @param studies Named list of study specifications. Each element is a list with:
#'   - `E` -- observed mean vector
#'   - `V` -- observed covariance matrix or variance vector (auto-detected)
#'   - `n` -- sample size
#'   - `times` -- numeric vector of observation times
#'   - `ev` -- `rxode2::et()` dosing event table
#'   - `method` -- `"cov"` or `"var"` (optional; auto-detected from `V`)
#'
#'   **Multi-compartment (multiple observed outputs).** To fit several observed
#'   compartments simultaneously (e.g. plasma and brain/CSF), give the study an
#'   `observations` list instead of top-level `E`/`V`/`times`. Each entry is one
#'   observed output with its own `output` (the model prediction variable, e.g.
#'   `"cp"` or `"cCSF"`), `times`, `E`, `V` and -- for independent fits -- `ev`
#'   and `n`. Pass the endpoint names to [admData()], e.g.
#'   `admData(c("cp", "cCSF"))`, so nlmixr2 recognises every endpoint. There are
#'   two modes:
#'
#'   * *Independent* -- each observed output has its own `n`/`ev` (separate
#'     experiments / subjects, e.g. a plasma study and a brain study combined for
#'     meta-analysis). The outputs are independent likelihood blocks and the
#'     aggregate `-2LL` is their sum.
#'   * *Joint (same subjects)* -- the outputs are measured on the SAME subjects.
#'     Give the study a shared `n` and `ev`, and a joint covariance either as a
#'     study-level full matrix `V` (blocks in `observations` order) or as
#'     per-output marginal `V` plus a `cross` list of cross-covariance blocks
#'     keyed `"outA:outB"` (each `length(times_A)` x `length(times_B)`; omitted
#'     pairs are zero). The compartments are then scored by a single MVN over the
#'     stacked vector with shared random effects. `est = "adirmc"` does not
#'     support multiple observed outputs; use `"admc"`, `"adfo"` or `"adgh"`.
#'
#'   **Long format (one row per endpoint/time).** As an alternative to the
#'   `observations` list, a study may carry a `data` frame that keys each observed
#'   summary by endpoint, the way nlmixr2 keys observations by `DVID`/`CMT`. The
#'   frame needs an endpoint column (`DVID`, `CMT` or `output`), a time column
#'   (`TIME`), a mean column (`E`) and -- unless a joint `V` is given -- a
#'   variance column (`V`) or an SD column (`SD`). It is normalised into exactly
#'   the same units as the `observations` form, so the two are interchangeable:
#'
#'   ```
#'   # independent blocks: per-row variances; optional per-endpoint `n` column
#'   # and per-endpoint `ev` (a list of event tables keyed by endpoint)
#'   list(n = 60L, ev = ev,
#'        data = data.frame(DVID = c("cp", "cp", "cCSF"), TIME = c(1, 2, 2),
#'                          E = c(9.1, 7.4, 2.2), V = c(1.2, 0.9, 0.1)))
#'
#'   # joint (same subjects): ONE stacked covariance whose rows/cols align with
#'   # the rows of `data` -- no `cross` blocks to assemble by hand
#'   list(n = 60L, ev = ev, data = data.frame(DVID = ..., TIME = ..., E = ...),
#'        V = V_joint)
#'   ```
#'
#'   A study-level `V` (or an explicit `joint = TRUE`) marks the endpoints as
#'   same-subject; without one, each endpoint is an independent likelihood block.
#'   Endpoints are stacked in the order they first appear in `data`.
#' @param resid_nodes Gauss-Hermite nodes used to integrate the RESIDUAL for a
#'   transform-both-sides endpoint (`boxCox`, `yeoJohnson`, `logitNorm`,
#'   `probitNorm`), where `y = g(h(f) + sigma*eps)` has no closed-form mean and
#'   variance. Ignored by every other error model, which has closed forms. Default
#'   81. Measured worst-case relative error against an independent quadrature, over
#'   all four transforms and residual SD of 0.5, 1, 2 and 3: n = 15 gives 5.7e-2,
#'   31 gives 4.5e-3, 81 gives 5.0e-5. The error is dominated by large residual SD;
#'   at SD <= 1, n = 31 already gives 1e-7 or better.
#'
#'   This is an ACCURACY dial, not a speed one. The quadrature is linear in
#'   `resid_nodes` in isolation (~50 us at 15, 300 us at 81 for an 8-row study) but
#'   negligible beside the ODE solve: a full NLL evaluation measured 0.750 s per 60
#'   evaluations at BOTH 31 and 81 nodes. Raise it if you have a saturating endpoint
#'   with a large residual SD; there is little to gain by lowering it.
#' @param n_sim Number of Monte Carlo samples per NLL evaluation.
#' @param sampling Sampling method for eta draws: `"sobol"` (Sobol, default),
#'   `"halton"` (Halton), `"torus"` (Kronecker/torus), `"lhs"` (Latin hypercube),
#'   or `"rnorm"` (iid normal).
#' @param algorithm nloptr algorithm string, or `NULL` (default) to pick the
#'   default that matches `grad`: `"NLOPT_LD_LBFGS"` with a gradient,
#'   `"NLOPT_LN_BOBYQA"` when `grad = "none"`. Any algorithm reported by
#'   [nloptr::nloptr.print.options()] is accepted (e.g. `"NLOPT_LD_MMA"`,
#'   `"NLOPT_LN_NELDERMEAD"`). An explicit algorithm is reconciled with `grad`:
#'   when `grad = "none"` a gradient-based algorithm (`NLOPT_LD_*` /
#'   `NLOPT_GD_*`) falls back to `"NLOPT_LN_BOBYQA"`; when a gradient is
#'   requested a derivative-free algorithm (`NLOPT_LN_*` / `NLOPT_GN_*`) turns
#'   the gradient off. Both emit a message.
#' @param maxeval Maximum number of optimizer function evaluations.
#' @param ftol_rel Relative function-value tolerance for convergence.
#' @param print Print progress every this many evaluations (0 = silent).
#' @param seed Random seed for reproducibility.
#' @param cores Number of OpenMP threads for `rxSolve()`. Defaults to
#'   `rxode2::rxCores()`. `rxSolve()` parallelises over subjects, so this is the
#'   main speed lever for the MC estimators; when `workers > 1` it is a *total*
#'   budget, split across the workers.
#' @param nDisplayProgress Passed to `rxSolve()`: the solver shows its text
#'   progress bar only once a single solve exceeds this many subjects. The
#'   default (`.Machine$integer.max`) keeps the bar off, which is what you want
#'   for scripts, vignettes and logs; lower it (e.g. `1000L`) to see solver
#'   progress during long interactive fits.
#' @param grad Gradient mode: `"sens"` (sensitivity equations, default), `"fd"`
#'   (forward finite differences), `"cfd"` (central finite differences), or
#'   `"none"` (derivative-free). A warning is issued when `"sens"` is requested
#'   but the sensitivity model is unavailable; the estimator then falls back to
#'   forward finite differences.
#' @param grad_h Step size for finite-difference gradient evaluation during
#'   optimization (used by `grad = "fd"` or `"cfd"`). The default 1e-4 is near
#'   the optimal balance between truncation error (grows with `h`) and MC noise
#'   amplification (grows as `1/h`) for forward FD. Central FD (`"cfd"`) has a
#'   slightly wider optimum around 1e-3, but 1e-4 works well for both.
#' @param cov_h Inner FD step for the gradient-based Hessian (only used when
#'   `covMethod = "r"` and `grad != "none"`). Each gradient evaluation has MC
#'   noise of order `sigma / cov_h`; the Hessian divides that noise by the outer
#'   step, giving total noise `sigma / (cov_h * cov_h_outer * |p|)`. `cov_h = 1e-3`
#'   balances truncation error and noise amplification. Increase to `1e-2` if the
#'   Hessian is non-positive definite.
#' @param cov_h_outer Outer step scale for the numerical Hessian. The actual step
#'   for parameter `p` is `max(|p|, 0.1) * cov_h_outer`. Applied to both the
#'   gradient-FD Hessian (`grad != "none"`) and the NLL-FD Hessian
#'   (`grad = "none"`). Default `eps^(1/5)` (~2.5e-3) is larger than the
#'   textbook `eps^(1/4)` to account for MC noise in NLL and gradient evaluations;
#'   empirically it matches the analytical (sensitivity-equation) Hessian ground
#'   truth. Increase (e.g. to `5e-3` or `1e-2`) if the Hessian is non-positive
#'   definite.
#' @param grad_bounds Box-constraint half-width when using gradients.
#' @param covMethod Covariance method: `"r"` (numerical Hessian over the
#'   structural, residual-error and omega parameters) or `"none"`. Omega is
#'   included because excluding it also biases the STRUCTURAL standard errors
#'   downward -- a theta carrying an eta is correlated with that eta's variance.
#'   If the weakly-identified omega Cholesky makes the Hessian non-positive
#'   definite, the structural + residual sub-block is reported with a warning.
#'
#'   All three blocks are reported on the scale the ESTIMATES are printed on, as
#'   `nlmixr2est` does: structural thetas on the log/optimizer scale, residual
#'   error as an SD, and omega as the variance/covariance entries (named
#'   `om.<eta>` and `cov.<eta_i>.<eta_j>`). The omega block is rotated by the
#'   full Jacobian of Omega with respect to the log-Cholesky, which is not
#'   diagonal once omega is correlated.
#' @param cov_n_sim Number of MC samples for the covariance (Hessian) step.
#'   More samples reduce MC noise in NLL evaluations. The NLL-based Hessian
#'   (`grad = "none"`) uses a central second difference of the NLL with the
#'   same Sobol sequence (CRN) at every perturbed point, so noise largely
#'   cancels and `cov_n_sim = 10000` (default) is sufficient for most models.
#' @param n_restarts Number of optimization restarts. Runs in parallel when
#'   `workers > 1`.
#' @param restart_sd Standard deviation of structural theta perturbations for
#'   restart initialisation.
#' @param workers Number of parallel workers for multi-restart. `1` (default)
#'   runs restarts sequentially. Values `> 1` run the restarts on a pool of
#'   background R processes (mirai daemons), which behaves the same way on every
#'   platform. Requires the `mirai` package. Workers are stopped automatically
#'   after the restart phase so all cores are available for the Hessian step; if
#'   a fit is interrupted, `admStopWorkers()` cleans up any survivors.
#' @param rxControl `rxode2::rxControl()` object. Created automatically when `NULL`.
#' @param addProp How combined additive+proportional error is parameterised in
#'   the nlmixr2 output tables: `"combined2"` (default, variance form) or
#'   `"combined1"` (SD form). Has no effect on admixr2's own estimation; passed
#'   to `nlmixr2est::foceiControl()` for the table/output machinery only.
#' @param calcTables,compress,ci,sigdig,sigdigTable,optExpression,sumProd,literalFix
#'   Passed to `nlmixr2est::foceiControl()` for the table/output machinery.
#' @param returnAdmr If `TRUE`, return a plain list instead of a full nlmixr2
#'   fit object (useful for debugging).
#' @param ... Additional arguments (none allowed; triggers an error).
#'
#' @return An object of class `admControl`.
#'
#' @examples
#' # Minimal control object -- inspect defaults
#' ctl <- admControl()
#' ctl$n_sim
#' ctl$algorithm
#'
#' # Override key settings without fitting
#' ctl2 <- admControl(
#'   n_sim    = 2000L,
#'   maxeval  = 300L,
#'   grad     = "fd",
#'   seed     = 42L
#' )
#'
#' \donttest{
#' library(rxode2)
#' library(nlmixr2)
#'
#' data("examplomycin")
#' obs   <- examplomycin[examplomycin$EVID == 0, ]
#' obs   <- obs[order(obs$ID, obs$TIME), ]
#' times <- sort(unique(obs$TIME))
#' ids   <- unique(obs$ID)
#' dv_mat <- do.call(rbind, lapply(ids, function(i) {
#'   sub <- obs[obs$ID == i, ]; sub$DV[order(sub$TIME)]
#' }))
#' E <- colMeans(dv_mat)
#' V <- cov.wt(dv_mat, method = "ML")$cov
#'
#' pk_model <- function() {
#'   ini({
#'     tcl <- log(5);  tv1 <- log(12); tv2 <- log(25)
#'     tq  <- log(12); tka <- log(1.2)
#'     prop.sd <- c(0, 0.2)
#'     eta.cl ~ 0.09; eta.v1 ~ 0.09; eta.v2 ~ 0.09
#'     eta.q  ~ 0.09; eta.ka ~ 0.09
#'   })
#'   model({
#'     cl <- exp(tcl + eta.cl); v1 <- exp(tv1 + eta.v1)
#'     v2 <- exp(tv2 + eta.v2); q  <- exp(tq  + eta.q)
#'     ka <- exp(tka + eta.ka)
#'     d/dt(depot)      <- -ka * depot
#'     d/dt(central)    <- ka * depot - (cl/v1 + q/v1) * central + (q/v2) * peripheral
#'     d/dt(peripheral) <- (q/v1) * central - (q/v2) * peripheral
#'     cp <- central / v1
#'     cp ~ prop(prop.sd)
#'   })
#' }
#'
#' fit <- nlmixr2(
#'   pk_model, admData(), est = "admc",
#'   control = admControl(
#'     studies  = list(study1 = list(E = E, V = V, n = length(ids),
#'                                   times = times, ev = et(amt = 100))),
#'     n_sim    = 1000L,
#'     maxeval  = 200L
#'   )
#' )
#' print(fit)
#' }
#'
#' @export
admControl <- function(
    studies    = list(),
    n_sim      = 5000L,
    sampling   = c("sobol", "halton", "torus", "lhs", "rnorm"),
    algorithm  = NULL,
    maxeval    = 500L,
    ftol_rel   = .Machine$double.eps^2,
    print      = 10L,
    seed       = 12345L,
    cores      = rxode2::rxCores(),
    nDisplayProgress = .Machine$integer.max,
    grad        = c("sens", "fd", "cfd", "none"),
    grad_h      = 1e-4,
    cov_h       = 1e-3,
    cov_h_outer = .Machine$double.eps^(1/5),
    grad_bounds = 5,
    covMethod   = c("r", "none"),
    cov_n_sim   = 10000L,
    n_restarts  = 1L,
    restart_sd  = 0.5,
    workers     = 1L,
    rxControl     = NULL,
    calcTables    = FALSE,
    compress      = TRUE,
    ci            = 0.95,
    sigdig        = 4,
    sigdigTable   = NULL,
    addProp       = c("combined2", "combined1"),
    optExpression = TRUE,
    sumProd       = FALSE,
    literalFix    = TRUE,
    returnAdmr    = FALSE,
    # LAST on purpose: inserting an argument mid-signature silently rebinds every
    # positional call -- admControl(studies, 20000L) used to set n_sim = 20000.
    resid_nodes = 81L,
    ...) {

  .xtra <- list(...)
  if (length(.xtra) > 0)
    stop("admControl: unused argument(s): ",
         paste(paste0("'", names(.xtra), "'"), collapse = ", "), call. = FALSE)

  addProp  <- match.arg(addProp)
  grad     <- match.arg(grad)
  sampling <- match.arg(sampling)

  checkmate::assertList(studies)
  checkmate::assertIntegerish(n_sim,   lower = 1L, len = 1, .var.name = "n_sim")
  # A residual quadrature needs a real grid. .adghNodes1() refuses m < 1, but it
  # accepts 1..4 happily and returns a rule that integrates nothing usefully --
  # the measured error at 5 nodes is already 3.3e-1. Refuse here, where the
  # message can name the argument, rather than silently scoring a wrong NLL.
  checkmate::assertIntegerish(resid_nodes, lower = 5L, len = 1,
                              .var.name = "resid_nodes")
  checkmate::assertIntegerish(maxeval, lower = 1L, len = 1, .var.name = "maxeval")
  checkmate::assertNumeric(ftol_rel,   lower = 0,  len = 1, .var.name = "ftol_rel")
  checkmate::assertIntegerish(print,   lower = 0L, len = 1, .var.name = "print")
  checkmate::assertIntegerish(seed,                len = 1, .var.name = "seed")
  checkmate::assertIntegerish(cores,   lower = 1L, len = 1, .var.name = "cores")
  checkmate::assertIntegerish(nDisplayProgress, lower = 1L, len = 1,
                              .var.name = "nDisplayProgress")
  checkmate::assertNumeric(grad_h,      lower = 0,  len = 1, .var.name = "grad_h")
  checkmate::assertNumeric(cov_h,       lower = 0, len = 1, .var.name = "cov_h")
  checkmate::assertNumeric(cov_h_outer, lower = 0, len = 1, .var.name = "cov_h_outer")
  checkmate::assertNumeric(grad_bounds, lower = 0,  len = 1, .var.name = "grad_bounds")
  covMethod <- match.arg(covMethod)
  checkmate::assertIntegerish(cov_n_sim,   lower = 1L, len = 1, .var.name = "cov_n_sim")
  checkmate::assertIntegerish(n_restarts,  lower = 1L, len = 1, .var.name = "n_restarts")
  checkmate::assertNumeric(restart_sd,     lower = 0,  len = 1, .var.name = "restart_sd")
  checkmate::assertIntegerish(workers,     lower = 1L, len = 1, .var.name = "workers")
  if (workers > 1L && cores < workers)
    message(sprintf(
      "admControl: cores (%d) < workers (%d) -- each worker will request 1 rxSolve thread.",
      as.integer(cores), as.integer(workers)
    ))
  checkmate::assertNumeric(ci,         lower = 0, upper = 1, len = 1, .var.name = "ci")
  checkmate::assertIntegerish(sigdig,  lower = 1L, len = 1, .var.name = "sigdig")
  checkmate::assertLogical(returnAdmr,             len = 1, .var.name = "returnAdmr")

  .algo     <- .admResolveAlgorithm(algorithm, grad,
                                    .var.name = "admControl: algorithm")
  algorithm <- .algo$algorithm
  grad      <- .algo$grad

  if (is.null(rxControl))   rxControl   <- rxode2::rxControl(sigdig = sigdig)
  if (is.null(sigdigTable)) sigdigTable <- max(round(sigdig), 3L)

  .ret <- list(
    studies       = studies,
    resid_nodes   = as.integer(resid_nodes),
    n_sim         = as.integer(n_sim),
    sampling      = sampling,
    algorithm     = algorithm,
    maxeval       = as.integer(maxeval),
    ftol_rel      = ftol_rel,
    print         = as.integer(print),
    seed          = as.integer(seed),
    cores         = as.integer(cores),
    nDisplayProgress = as.integer(nDisplayProgress),
    grad          = grad,
    grad_h        = grad_h,
    cov_h         = cov_h,
    cov_h_outer   = cov_h_outer,
    grad_bounds   = grad_bounds,
    covMethod     = covMethod,
    cov_n_sim     = as.integer(cov_n_sim),
    n_restarts    = as.integer(n_restarts),
    restart_sd    = restart_sd,
    workers       = as.integer(workers),
    rxControl     = rxControl,
    calcTables    = calcTables,
    compress      = compress,
    ci            = ci,
    sigdig        = sigdig,
    sigdigTable   = as.integer(sigdigTable),
    addProp       = addProp,
    optExpression = optExpression,
    sumProd       = sumProd,
    literalFix    = literalFix,
    returnAdmr    = returnAdmr
  )
  class(.ret) <- "admControl"
  .ret
}

# -- nlmixr2 S3 hooks -----------------------------------------------------------

#' @noRd
getValidNlmixrCtl.admc <- function(control) {
  if (inherits(control, "admControl")) return(control)
  .ctl <- control[[1]]
  if (inherits(.ctl, "admControl")) return(.ctl)
  if (is.list(.ctl) && "studies" %in% names(.ctl))
    return(do.call(admControl, .ctl[intersect(names(.ctl), names(formals(admControl)))]))
  admControl()
}

#' @noRd
nmObjHandleControlObject.admControl <- function(control, env) {
  assign("admControl", control, envir = env)
}

#' @noRd
nmObjGetControl.admc <- function(x, ...) {
  .env <- x[[1]]
  for (.nm in c("admControl", "control")) {
    if (exists(.nm, .env)) {
      .ctl <- get(.nm, .env)
      if (inherits(.ctl, "admControl")) return(.ctl)
    }
  }
  stop("cannot find admc control object", call. = FALSE)
}

# -- Aggregate -2LL ------------------------------------------------------------

.admNLL <- function(p, pinfo, studies, z_list, rxMod, output_var,
                    params_list, cores) {
  pars <- tryCatch(.admUnpack(p, pinfo), error = function(e) NULL)
  if (is.null(pars)) return(Inf)
  if (pinfo$n_eta > 0 && any(diag(pars$omega) <= 0)) return(Inf)

  nll2     <- 0
  for (i in seq_along(studies)) {
    s       <- studies[[i]]
    ov      <- s$output %||% output_var
    z       <- z_list[[i]]
    eta_mat <- z %*% t(pars$L)
    colnames(eta_mat) <- pinfo$eta_col_names

    # Joint (same-subject) unit: one shared-eta solve produces every output;
    # score the stacked vector with a single MVN over the joint covariance.
    if (isTRUE(s$is_joint)) {
      cp_mat <- tryCatch(
        .admSimulateJoint(rxMod, pars$struct, pinfo$sigma_names, eta_mat, s,
                          params_list[[i]], cores, pinfo$nDisplayProgress),
        error = function(e) NULL)
      if (is.null(cp_mat) || anyNA(cp_mat)) return(Inf)
      mu_struct <- colMeans(cp_mat)
      V_pred    <- crossprod(sweep(cp_mat, 2L, mu_struct)) / nrow(cp_mat)
      jr        <- .admJointResidual(mu_struct, V_pred, s, pinfo, pars$sigma_var)
      nll2      <- nll2 + nll_cov_cpp(as.numeric(s$E), s$V, jr$mu, jr$V, s$n)
      if (!is.finite(nll2)) return(Inf)
      next
    }

    cp_mat <- tryCatch(
      .admSimulate(rxMod, pars$struct, pinfo$sigma_names, eta_mat, s,
                   ov, params_list[[i]], cores, pinfo$nDisplayProgress),
      error = function(e) NULL)
    if (is.null(cp_mat) || anyNA(cp_mat)) return(Inf)

    # Restrict the residual error to this output's sigma(s) so a multi-output fit
    # does not apply another compartment's error model here (no-op single-output).
    ar <- .admUnitResidRows(pinfo, ov, pars$sigma_var, length(s$times),
                            phi = attr(cp_mat, "phi"))
    if (.admResidCppOK(ar)) {
      # forms the fused C++ kernels implement (combined1/2, lnorm, no ar)
      if (identical(s$method, "var")) {
        nll2 <- nll2 + nll_var_from_samples_cpp(cp_mat, as.numeric(s$E), s$v_diag,
                                                s$n, ar$form, ar$a2, ar$b2, ar$cc)
      } else {
        nll2 <- nll2 + nll_cov_from_samples_cpp(cp_mat, as.numeric(s$E), s$V,
                                                s$n, ar$form, ar$a2, ar$b2, ar$cc)
      }
    } else {
      # everything else: assemble the moments in R (form-agnostic), then score
      # with the plain kernels. See .admResidCppOK() for why.
      mu_s <- colMeans(cp_mat)
      cpc  <- sweep(cp_mat, 2L, mu_s)
      Vs   <- crossprod(cpc) / nrow(cp_mat)
      ap   <- .admResidApply(mu_s, diag(Vs), ar, s$times, Vs)
      if (identical(s$method, "var")) {
        nll2 <- nll2 + nll_var_cpp(as.numeric(s$E), s$v_diag, ap$mu, ap$dv, s$n)
      } else {
        Vp <- .admApplyResidTail(Vs, ap)
        nll2 <- nll2 + nll_cov_cpp(as.numeric(s$E), s$V, ap$mu, Vp, s$n)
      }
    }
    if (!is.finite(nll2)) return(Inf)
  }
  nll2
}

# Plain forward/central FD of the aggregate NLL. Used when a fit contains a
# joint (same-subject) unit, whose stacked-MVN gradient is not covered by the
# per-unit analytical decomposition in .admGrad. The fixed z_list makes this a
# common-random-number FD, stable despite MC noise.
.admNLLGradFD <- function(p, pinfo, studies, z_list, rxMod, output_var,
                          params_list, cores, h, use_central = FALSE) {
  f0 <- if (use_central) NA_real_ else
    .admNLL(p, pinfo, studies, z_list, rxMod, output_var, params_list, cores)
  vapply(seq_along(p), function(k) {
    pp <- p; pp[k] <- p[k] + h
    fp <- .admNLL(pp, pinfo, studies, z_list, rxMod, output_var, params_list, cores)
    if (use_central) {
      pm <- p; pm[k] <- p[k] - h
      fm <- .admNLL(pm, pinfo, studies, z_list, rxMod, output_var, params_list, cores)
      (fp - fm) / (2 * h)
    } else {
      (fp - f0) / h
    }
  }, double(1))
}

# -- Gradient (forward / central FD + sensitivity) -----------------------------

.admGrad <- function(p, pinfo, studies, z_list, rxMod, output_var,
                     params_list, cores, h, sensModel = NULL,
                     use_central = FALSE) {
  pars <- tryCatch(.admUnpack(p, pinfo), error = function(e) NULL)
  if (is.null(pars)) return(rep(NA_real_, length(p)))

  n_s   <- length(pinfo$struct_names)
  n_e   <- length(pinfo$sigma_names)
  n_eta <- pinfo$n_eta
  n_sim <- nrow(z_list[[1]])

  eta_col_names <- pinfo$eta_col_names

  grad <- numeric(length(p))
  names(grad) <- names(p)

  for (si in seq_along(studies)) {
    s   <- studies[[si]]
    ov  <- s$output %||% output_var
    z   <- z_list[[si]]
    pdf <- params_list[[si]]

    eta_mat           <- z %*% t(pars$L)
    colnames(eta_mat) <- eta_col_names

    unpaired_k <- which(vapply(pinfo$struct_names, function(nm)
      is.null(pinfo$struct_has_eta) || !isTRUE(pinfo$struct_has_eta[nm]), logical(1)))
    n_unp <- length(unpaired_k)

    # --- Joint (same-subject) analytical gradient ----------------------------
    # Stacked-MVN gradient over all outputs: eta/omega/sigma analytical on the
    # joint covariance (shared-eta sensitivities per output), paired struct
    # thetas via the eta path, unpaired struct thetas via CRN-FD. Requires the
    # sens model + etas; the driver routes joint fits without sens to full FD.
    if (isTRUE(s$is_joint)) {
      if (n_eta == 0L || is.null(sensModel)) return(rep(NA_real_, length(p)))
      js <- .admSimulateJointSens(sensModel, pars$struct, pinfo$sigma_names,
                                  eta_mat, s, cores, pinfo$nDisplayProgress,
                                  pars$sigma_var)
      if (is.null(js) || anyNA(js$cp_mat)) return(rep(NA_real_, length(p)))
      cp_mat <- js$cp_mat; dpred_list <- js$dpred_list
      n_t    <- s$n_total
      mu_struct <- colMeans(cp_mat)
      V_struct  <- crossprod(sweep(cp_mat, 2L, mu_struct)) / n_sim
      var_f     <- diag(V_struct)                    # Var_eta(f), pre-residual
      jr <- .admJointResidual(mu_struct, V_struct, s, pinfo, pars$sigma_var)
      mu <- jr$mu; V <- jr$V
      cp_c  <- sweep(cp_mat, 2L, mu_struct)
      r     <- as.numeric(s$E) - mu
      cholV <- tryCatch(chol(V), error = function(e) NULL)
      if (is.null(cholV)) return(rep(NA_real_, length(p)))
      invV         <- chol2inv(cholV)
      dNLL_dmu     <- s$n * as.numeric(-2 * invV %*% r)
      dNLL_dV      <- s$n * (invV - invV %*% (s$V + tcrossprod(r)) %*% invV)
      dNLL_dV_diag <- diag(dNLL_dV)

      # sigma mu-path coupling and sigma gradient, each on its output's rows.
      arr <- .admResidRows(pinfo, .admRowOutput(s, n_t), pars$sigma_var, n_t)
      .rt_j <- .admRowTimes(s, n_t)
      .dres <- .admResidDeriv(mu_struct, var_f, arr, pinfo)   # once, reused below
      sigma_mu_scale <- .admResidMuCoupling(mu_struct, arr, pinfo,
                                            dNLL_dV_diag, dNLL_dmu, var_f,
                                            dNLL_dV, V_struct, .rt_j, deriv = .dres)
      # V_pred -> V_struct chain (see the single-output branch below)
      vchain     <- .admResidVChain(mu_struct, var_f, arr, pinfo,
                                    .admRowTimes(s, length(mu_struct)), deriv = .dres)
      dNLL_dV_s  <- dNLL_dV * vchain
      diag(dNLL_dV_s) <- diag(dNLL_dV_s) +
        dNLL_dmu * (attr(vchain, "dmu_dv0") %||% numeric(n_t))

      if (n_eta > 0L) {
        eta_rows_df  <- pinfo$eta_rows_df
        D_mat        <- do.call(cbind, dpred_list)
        z_diag_scale <- sweep(z, 2L, diag(pars$L) / 2, "*")
        go <- adm_grad_eta_omega_cpp(
          cp_c, D_mat, z_diag_scale, z, dNLL_dV_s, dNLL_dmu, sigma_mu_scale,
          as.integer(eta_rows_df$neta1), as.integer(eta_rows_df$neta2),
          n_t, n_eta)
        for (j in seq_len(n_eta))
          if (!is.null(pinfo$struct_eta_idx) && !is.na(pinfo$struct_eta_idx[j]))
            grad[pinfo$struct_eta_idx[j]] <- grad[pinfo$struct_eta_idx[j]] + go$eta_grad[j]
        k_om <- n_s + n_e
        for (r_idx in seq_len(nrow(eta_rows_df))) {
          k_om <- k_om + 1L
          grad[k_om] <- grad[k_om] + go$omega_grad[r_idx]
        }
      }

      grad[n_s + seq_len(n_e)] <- grad[n_s + seq_len(n_e)] +
        .admSigmaGrad(mu_struct, arr, pinfo, dNLL_dV_diag, dNLL_dmu, var_f,
                    dNLL_dV, .rt_j, V_struct, deriv = .dres)

      # Unpaired struct thetas. The augmented sens model carries d(pred)/d(theta)
      # directly (js$dtheta_list): it enters exactly like an eta direction, so the
      # same partial kernel the FD path feeds serves it -- but exactly, with no
      # step size. Without those columns (plain sens model), CRN-FD of the joint
      # contribution at fixed eta_mat, as before.
      if (n_unp > 0L) {
        eff_dmu_j <- dNLL_dmu + sigma_mu_scale
        if (!is.null(js$dtheta_list)) {
          for (bi in seq_len(n_unp)) {
            k_s   <- unpaired_k[bi]
            dpred <- js$dtheta_list[[pinfo$struct_names[k_s]]]
            grad[k_s] <- grad[k_s] +
              adm_grad_partial_cpp(cp_c, dpred, dNLL_dV_s, eff_dmu_j, 1 / n_sim)
          }
        } else {
          nll0 <- nll_cov_cpp(as.numeric(s$E), s$V, mu, V, s$n)
          for (bi in seq_len(n_unp)) {
            k_s <- unpaired_k[bi]
            pp  <- p; pp[k_s] <- p[k_s] + h
            pars_p <- .admUnpack(pp, pinfo)
            cp_p   <- .admSimulateJoint(rxMod, pars_p$struct, pinfo$sigma_names,
                                        eta_mat, s, pdf, cores, pinfo$nDisplayProgress)
            mus_p  <- colMeans(cp_p)
            jr_p   <- .admJointResidual(mus_p,
                                        crossprod(sweep(cp_p, 2L, mus_p)) / n_sim,
                                        s, pinfo, pars_p$sigma_var)
            nllp <- nll_cov_cpp(as.numeric(s$E), s$V, jr_p$mu, jr_p$V, s$n)
            grad[k_s] <- grad[k_s] + (nllp - nll0) / h
          }
        }
      }
      next
    }

    use_sens <- !is.null(sensModel) && n_eta > 0
    batched_hi <- NULL; batched_lo <- NULL
    theta_sens <- NULL
    if (use_sens) {
      sens_out <- .admSimulateSens(sensModel, pars$struct, pinfo$sigma_names,
                                   eta_mat, s, cores, pinfo$nDisplayProgress,
                                   pars$sigma_var)
      if (is.null(sens_out) || anyNA(sens_out$cp_mat)) {
        use_sens <- FALSE
      } else {
        cp_mat     <- sens_out$cp_mat
        dpred_list <- sens_out$dpred_list
        # d(pred)/d(theta) for the unpaired thetas (augmented sens model only);
        # NULL -> the FD block below handles them, as before.
        theta_sens <- sens_out$dtheta_list
      }
    }
    if (!use_sens) {
      n_t       <- length(s$times)
      col_nms   <- colnames(pdf)
      n_cols    <- length(col_nms)
      n_fwd_eta <- if (n_eta > 0L) (if (use_central) 2L * n_eta else n_eta) else 0L
      n_fwd_unp <- if (n_unp > 0L) (if (use_central) 2L * n_unp else n_unp) else 0L
      n_runs    <- 1L + n_fwd_eta + n_fwd_unp

      pdf_big <- matrix(0, nrow = n_runs * n_sim, ncol = n_cols,
                        dimnames = list(NULL, col_nms))
      pdf_big[, grep("^rxerr", colnames(pdf_big), value = TRUE)] <- 1L
      for (nm in names(pars$struct)) pdf_big[, nm] <- pars$struct[nm]
      for (nm in pinfo$sigma_names)  pdf_big[, nm] <- 0
      if (n_eta > 0L) pdf_big[seq_len(n_sim), eta_col_names] <- eta_mat

      # eta perturbation rows
      if (n_eta > 0L) {
        if (use_central) {
          for (j in seq_len(n_eta)) {
            rows_hi <- n_sim * (2L*j - 1L) + seq_len(n_sim)
            rows_lo <- n_sim * (2L*j)      + seq_len(n_sim)
            eta_hi  <- eta_mat; eta_hi[, j] <- eta_hi[, j] + h
            eta_lo  <- eta_mat; eta_lo[, j] <- eta_lo[, j] - h
            pdf_big[rows_hi, eta_col_names] <- eta_hi
            pdf_big[rows_lo, eta_col_names] <- eta_lo
          }
        } else {
          for (j in seq_len(n_eta)) {
            rows_hi <- n_sim * j + seq_len(n_sim)
            eta_hi  <- eta_mat; eta_hi[, j] <- eta_hi[, j] + h
            pdf_big[rows_hi, eta_col_names] <- eta_hi
          }
        }
      }

      # struct perturbation rows appended after eta block
      if (n_unp > 0L) {
        if (use_central) {
          for (bi in seq_len(n_unp)) {
            rows_hi <- n_sim * (n_fwd_eta + 2L*bi - 1L) + seq_len(n_sim)
            rows_lo <- n_sim * (n_fwd_eta + 2L*bi)      + seq_len(n_sim)
            nm_u    <- pinfo$struct_names[unpaired_k[bi]]
            if (n_eta > 0L) {
              pdf_big[rows_hi, eta_col_names] <- eta_mat
              pdf_big[rows_lo, eta_col_names] <- eta_mat
            }
            pdf_big[rows_hi, nm_u] <- pars$struct[nm_u] + h
            pdf_big[rows_lo, nm_u] <- pars$struct[nm_u] - h
          }
        } else {
          for (bi in seq_len(n_unp)) {
            rows_hi <- n_sim * (n_fwd_eta + bi) + seq_len(n_sim)
            nm_u    <- pinfo$struct_names[unpaired_k[bi]]
            if (n_eta > 0L) pdf_big[rows_hi, eta_col_names] <- eta_mat
            pdf_big[rows_hi, nm_u] <- pars$struct[nm_u] + h
          }
        }
      }

      out_b  <- rxode2::rxSolve(rxMod, params = as.data.frame(pdf_big),
                                 events = s$ev_full, cores = cores,
                                 nDisplayProgress = pinfo$nDisplayProgress)
      keep_b <- out_b[["time"]] %in% s$times
      # beta: the prediction is DERIVED from two solved columns, mu = b1/(b1+b2),
      # and the precision phi = b1 + b2 is solved rather than fitted. Reading
      # s$output alone gave this path the raw first shape parameter (~40, not a
      # probability) and left arr$phi at NA, so the FD gradient scored a different
      # function than .admNLL and every entry of V_pred came back NA.
      # Inlined, matching .admSimulate -- see the dev-mode daemon note in simulate.R.
      .phi_b <- NULL
      vals_b <- if (!is.null(s$out_pair)) {
        .b1 <- out_b[[s$out_pair[[1L]]]][keep_b]; .b2 <- out_b[[s$out_pair[[2L]]]][keep_b]
        .phi_b <- .b1 + .b2
        .b1 / { .d <- .phi_b; .d[.d == 0] <- .Machine$double.eps; .d }
      } else {
        .v <- out_b[[ov]][keep_b]
        if (is.null(.v)) out_b[["ipredSim"]][keep_b] else .v
      }

      cp_mat <- matrix(vals_b[seq_len(n_sim * n_t)],
                       nrow = n_sim, ncol = n_t, byrow = TRUE)
      if (!is.null(.phi_b))
        attr(cp_mat, "phi") <- matrix(.phi_b[seq_len(n_sim * n_t)], nrow = n_sim,
                                      ncol = n_t, byrow = TRUE)[1L, ]
      if (anyNA(cp_mat)) return(rep(NA_real_, length(p)))

      dpred_list <- if (n_eta > 0L) {
        if (use_central) {
          lapply(seq_len(n_eta), function(j) {
            off_hi <- n_sim * n_t * (2L*j - 1L)
            off_lo <- n_sim * n_t * (2L*j)
            (matrix(vals_b[off_hi + seq_len(n_sim * n_t)], nrow = n_sim, ncol = n_t, byrow = TRUE) -
             matrix(vals_b[off_lo + seq_len(n_sim * n_t)], nrow = n_sim, ncol = n_t, byrow = TRUE)) / (2 * h)
          })
        } else {
          lapply(seq_len(n_eta), function(j) {
            off_hi <- n_sim * n_t * j
            (matrix(vals_b[off_hi + seq_len(n_sim * n_t)], nrow = n_sim, ncol = n_t, byrow = TRUE) -
             cp_mat) / h
          })
        }
      } else list()

      if (n_unp > 0L) {
        batched_hi <- lapply(seq_len(n_unp), function(bi) {
          off <- n_sim * n_t * (n_fwd_eta + if (use_central) 2L*bi - 1L else bi)
          matrix(vals_b[off + seq_len(n_sim * n_t)], nrow = n_sim, ncol = n_t, byrow = TRUE)
        })
        if (use_central) batched_lo <- lapply(seq_len(n_unp), function(bi) {
          off <- n_sim * n_t * (n_fwd_eta + 2L*bi)
          matrix(vals_b[off + seq_len(n_sim * n_t)], nrow = n_sim, ncol = n_t, byrow = TRUE)
        })
      }
    }

    n_t       <- length(s$times)
    # beta: the precision phi is SOLVED, not fitted -- .admSimulate returns it as
    # an attribute on cp_mat (see there). Everything else leaves arr$phi NA.
    arr       <- .admUnitResidRows(pinfo, ov, pars$sigma_var, n_t,
                                   phi = attr(cp_mat, "phi"))
    mu_struct <- colMeans(cp_mat)
    cp_c      <- sweep(cp_mat, 2L, mu_struct)

    is_var <- identical(s$method, "var")
    cov_f  <- NULL                                   # only the cov branch has one
    if (is_var) {
      var_f <- adm_col_sq_sum_cpp(cp_c) / n_sim      # Var_eta(f), pre-residual
      ap <- .admResidApply(mu_struct, var_f, arr)
      mu <- ap$mu; pv <- ap$dv
      r  <- as.numeric(s$E) - mu
      dNLL_dmu     <- s$n * as.numeric(-2 * r / pv)
      dNLL_dV_diag <- s$n * (1 / pv - s$v_diag / pv^2 - r^2 / pv^2)
    } else {
      cov_f <- crossprod(cp_c) / n_sim               # STRUCTURAL Cov_eta(f); keep it,
      V  <- cov_f                                    # the ms/sigma chain needs it
      var_f <- diag(V)                               # Var_eta(f), pre-residual
      ap <- .admResidApply(mu_struct, var_f, arr, s$times, cov_f)
      mu <- ap$mu
      V  <- .admApplyResidTail(V, ap)
      r  <- as.numeric(s$E) - mu
      cholV <- tryCatch(chol(V), error = function(e) NULL)
      if (is.null(cholV)) return(rep(NA_real_, length(p)))
      invV         <- chol2inv(cholV)
      dNLL_dmu     <- s$n * as.numeric(-2 * invV %*% r)
      dNLL_dV      <- s$n * (invV - invV %*% (s$V + tcrossprod(r)) %*% invV)
      dNLL_dV_diag <- diag(dNLL_dV)
    }

    # sigma_mu_scale: how the residual couples a change in mu_struct into the
    # objective (lnorm mean scaling + the residual variance's dependence on mu).
    # Computed once per study and reused across every gradient term.
    .dres <- .admResidDeriv(mu_struct, var_f, arr, pinfo)   # once, reused across terms
    sigma_mu_scale <- .admResidMuCoupling(mu_struct, arr, pinfo,
                                          dNLL_dV_diag, dNLL_dmu, var_f,
                                          if (is_var) NULL else dNLL_dV, cov_f,
                                          s$times, deriv = .dres)
    eff_dmu <- dNLL_dmu + sigma_mu_scale
    inv_n <- 1 / n_sim

    # V-path chain. dNLL_dV/dNLL_dV_diag are derivatives w.r.t. V_PRED, but every
    # kernel below differentiates the STRUCTURAL covariance Cov_eta(f). Since
    # V_pred = ms(x)ms o Cov_eta(f) + diag(E[v]) -- and E[v] itself depends on
    # diag(Cov_eta(f)) whenever the residual is f-dependent -- the two differ by
    # .admResidVChain(): ms_i*ms_j off the diagonal, dv_dv0_i on it. Without this
    # the eta/omega/theta gradients silently assume d(V_pred)/d(V_struct) = I,
    # which is true only for purely additive error (chain == identity there, so
    # add() models are bit-identical).
    vchain <- .admResidVChain(mu_struct, var_f, arr, pinfo, s$times, deriv = .dres)
    # TBS only: mu depends on Var_eta(f), so the mean contributes to the same
    # d(var_f)/d(param) the kernels already chain. Zero for every other form.
    .dmv <- attr(vchain, "dmu_dv0") %||% numeric(n_t)
    dNLL_dV_diag_s <- dNLL_dV_diag * diag(vchain) + dNLL_dmu * .dmv
    if (!is_var) {
      dNLL_dV_s <- dNLL_dV * vchain
      diag(dNLL_dV_s) <- diag(dNLL_dV_s) + dNLL_dmu * .dmv
    }

    # Eta + omega gradient: one C++ call; var variant avoids n_txn_t intermediates.
    if (n_eta > 0L) {
      eta_rows_df  <- pinfo$eta_rows_df
      D_mat        <- do.call(cbind, dpred_list)
      z_diag_scale <- sweep(z, 2L, diag(pars$L) / 2, "*")
      neta1 <- as.integer(eta_rows_df$neta1)
      neta2 <- as.integer(eta_rows_df$neta2)
      go <- if (is_var)
        adm_grad_eta_omega_var_cpp(
          cp_c, D_mat, z_diag_scale, z,
          dNLL_dV_diag_s, dNLL_dmu, sigma_mu_scale,
          neta1, neta2,
          n_t, n_eta)
      else
        adm_grad_eta_omega_cpp(
          cp_c, D_mat, z_diag_scale, z,
          dNLL_dV_s, dNLL_dmu, sigma_mu_scale,
          neta1, neta2,
          n_t, n_eta)
      for (j in seq_len(n_eta)) {
        if (!is.null(pinfo$struct_eta_idx) && !is.na(pinfo$struct_eta_idx[j]))
          grad[pinfo$struct_eta_idx[j]] <- grad[pinfo$struct_eta_idx[j]] + go$eta_grad[j]
      }
      k_om <- n_s + n_e
      for (r_idx in seq_len(nrow(eta_rows_df))) {
        k_om <- k_om + 1L
        grad[k_om] <- grad[k_om] + go$omega_grad[r_idx]
      }
    }

    # Unpaired struct theta gradient.
    # theta_sens path: the augmented sens model returned d(pred)/d(theta) in the
    #   SAME solve as the eta sensitivities -- exact, and it removes the extra
    #   rxSolve the FD path needs (an rxSolve costs ~11 ms before it integrates
    #   anything). An unpaired theta enters mu and V exactly like an eta, so it
    #   feeds the same partial kernel; only the derivative source changes.
    # !use_sens path: batched_hi/batched_lo already extracted from the single big rxSolve above.
    # use_sens without theta columns (plain sens model): separate rxSolve for struct perturbations.
    if (n_unp > 0L) {
      if (!is.null(theta_sens)) {
        for (bi in seq_len(n_unp)) {
          k_s   <- unpaired_k[bi]
          dpred <- theta_sens[[pinfo$struct_names[k_s]]]
          grad[k_s] <- grad[k_s] +
            if (is_var)
              adm_grad_partial_var_cpp(cp_c, dpred, dNLL_dV_diag_s, eff_dmu, inv_n)
            else
              adm_grad_partial_cpp(cp_c, dpred, dNLL_dV_s, eff_dmu, inv_n)
        }
      } else if (!is.null(batched_hi)) {
        for (bi in seq_len(n_unp)) {
          k_s     <- unpaired_k[bi]
          cp_hi_s <- batched_hi[[bi]]
          dpred <- if (use_central && !is.null(batched_lo))
            (cp_hi_s - batched_lo[[bi]]) / (2 * h)
          else
            (cp_hi_s - cp_mat) / h
          grad[k_s] <- grad[k_s] +
            if (is_var)
              adm_grad_partial_var_cpp(cp_c, dpred, dNLL_dV_diag_s, eff_dmu, inv_n)
            else
              adm_grad_partial_cpp(cp_c, dpred, dNLL_dV_s, eff_dmu, inv_n)
        }
      } else {
        col_nms <- colnames(pdf)
        n_cols  <- length(col_nms)
        pdf_hi  <- matrix(0, nrow = n_unp * n_sim, ncol = n_cols,
                          dimnames = list(NULL, col_nms))
        pdf_hi[, grep("^rxerr", colnames(pdf_hi), value = TRUE)] <- 1L
        for (nm in names(pars$struct))       pdf_hi[, nm]              <- pars$struct[nm]
        for (j  in seq_along(eta_col_names)) pdf_hi[, eta_col_names[j]] <- rep(eta_mat[, j], n_unp)
        for (nm in pinfo$sigma_names)        pdf_hi[, nm]              <- 0
        for (bi in seq_len(n_unp)) {
          rows <- (bi - 1L) * n_sim + seq_len(n_sim)
          nm   <- pinfo$struct_names[unpaired_k[bi]]
          pdf_hi[rows, nm] <- pars$struct[nm] + h
        }
        out_hi  <- rxode2::rxSolve(rxMod, params = as.data.frame(pdf_hi),
                                    events = s$ev_full, cores = cores,
                                    nDisplayProgress = pinfo$nDisplayProgress)
        keep_hi <- out_hi[["time"]] %in% s$times
        vals_hi <- out_hi[[ov]][keep_hi]
        if (is.null(vals_hi)) vals_hi <- out_hi[["ipredSim"]][keep_hi]

        if (use_central) {
          pdf_lo <- pdf_hi
          for (bi in seq_len(n_unp)) {
            rows <- (bi - 1L) * n_sim + seq_len(n_sim)
            nm   <- pinfo$struct_names[unpaired_k[bi]]
            pdf_lo[rows, nm] <- pars$struct[nm] - h
          }
          out_lo  <- rxode2::rxSolve(rxMod, params = as.data.frame(pdf_lo),
                                      events = s$ev_full, cores = cores,
                                      nDisplayProgress = pinfo$nDisplayProgress)
          keep_lo <- out_lo[["time"]] %in% s$times
          vals_lo <- out_lo[[ov]][keep_lo]
          if (is.null(vals_lo)) vals_lo <- out_lo[["ipredSim"]][keep_lo]
        }

        for (bi in seq_len(n_unp)) {
          k_s     <- unpaired_k[bi]
          idx     <- (bi - 1L) * n_sim * n_t + seq_len(n_sim * n_t)
          cp_hi_s <- matrix(vals_hi[idx], nrow = n_sim, ncol = n_t, byrow = TRUE)
          dpred <- if (use_central) {
            cp_lo_s <- matrix(vals_lo[idx], nrow = n_sim, ncol = n_t, byrow = TRUE)
            (cp_hi_s - cp_lo_s) / (2 * h)
          } else {
            (cp_hi_s - cp_mat) / h
          }
          grad[k_s] <- grad[k_s] +
            if (is_var)
              adm_grad_partial_var_cpp(cp_c, dpred, dNLL_dV_diag_s, eff_dmu, inv_n)
            else
              adm_grad_partial_cpp(cp_c, dpred, dNLL_dV_s, eff_dmu, inv_n)
        }
      }
    }

    grad[n_s + seq_len(n_e)] <- grad[n_s + seq_len(n_e)] +
      .admSigmaGrad(mu_struct, arr, pinfo, dNLL_dV_diag, dNLL_dmu, var_f,
                    if (is_var) NULL else dNLL_dV, s$times, cov_f, deriv = .dres)
  }

  grad
}

# -- Batched NLL evaluation ----------------------------------------------------
# Evaluates NLL for a list of parameter vectors via one rxSolve call per
# study per chunk (instead of one call per vector). Reduces rxSolve call
# overhead from O(n_configs) to O(ceil(n_configs / chunk_size)).
# chunk_size controls peak memory: n_chunk * n_sim rows per rxSolve call.
.admNLLBatch <- function(p_list, pinfo, studies, z_list, rxMod, output_var,
                          params_list, cores, chunk_size = 30L) {
  n_c      <- length(p_list)
  if (n_c == 0L) return(numeric(0))
  # Joint (same-subject) units are scored by a stacked MVN that the batched
  # single-output solve below does not build; evaluate them one config at a time.
  if (any(vapply(studies, function(u) isTRUE(u$is_joint), logical(1))))
    return(vapply(p_list, function(p)
      .admNLL(p, pinfo, studies, z_list, rxMod, output_var, params_list, cores),
      double(1)))
  n_sim   <- nrow(z_list[[1L]])
  col_nms <- colnames(params_list[[1L]])
  n_cols  <- length(col_nms)

  pars_list <- vector("list", n_c)
  valid     <- logical(n_c)
  for (ci in seq_len(n_c)) {
    pars <- tryCatch(.admUnpack(p_list[[ci]], pinfo), error = function(e) NULL)
    if (!is.null(pars) && (pinfo$n_eta == 0L || all(diag(pars$omega) > 0))) {
      pars_list[[ci]] <- pars; valid[ci] <- TRUE
    }
  }

  nlls   <- rep(0.0, n_c)
  finite <- valid

  chunks <- split(seq_len(n_c), ceiling(seq_len(n_c) / chunk_size))

  for (si in seq_along(studies)) {
    s <- studies[[si]]
    ov <- s$output %||% output_var
    z <- z_list[[si]]
    n_t <- length(s$times)

    for (chunk in chunks) {
      n_chunk <- length(chunk)
      pdf_mat <- matrix(0, nrow = n_chunk * n_sim, ncol = n_cols,
                        dimnames = list(NULL, col_nms))
      pdf_mat[, grep("^rxerr", colnames(pdf_mat), value = TRUE)] <- 1L

      for (cii in seq_along(chunk)) {
        ci <- chunk[cii]
        if (!valid[ci] || !finite[ci]) next
        pars <- pars_list[[ci]]
        rows <- (cii - 1L) * n_sim + seq_len(n_sim)
        for (nm in pinfo$struct_names) pdf_mat[rows, nm] <- pars$struct[nm]
        if (pinfo$n_eta > 0L) {
          eta_mat <- z %*% t(pars$L)
          pdf_mat[rows, pinfo$eta_col_names] <- eta_mat
        }
      }

      out <- tryCatch(
        rxode2::rxSolve(rxMod, params = as.data.frame(pdf_mat),
                        events = s$ev_full, cores = cores,
                        nDisplayProgress = pinfo$nDisplayProgress),
        error = function(e) NULL)
      if (is.null(out)) { for (ci in chunk) finite[ci] <- FALSE; next }

      keep <- out[["time"]] %in% s$times
      # beta: the prediction is DERIVED from two solved columns and the precision
      # phi = b1 + b2 comes back with it. Reading s$output alone handed this path
      # the raw first shape parameter -- not a probability -- and left phi NULL, so
      # THIS function, which is the covMethod = "r" objective evaluator, scored a
      # different model from the one that was fitted. Same shape as .admSimulate;
      # inlined for the same reason (see the daemon note in simulate.R).
      .phi_all <- NULL
      vals <- if (!is.null(s$out_pair)) {
        .b1 <- out[[s$out_pair[[1L]]]][keep]; .b2 <- out[[s$out_pair[[2L]]]][keep]
        .phi_all <- .b1 + .b2
        .b1 / { .d <- .phi_all; .d[.d == 0] <- .Machine$double.eps; .d }
      } else {
        .v <- out[[ov]][keep]
        if (is.null(.v)) out[["ipredSim"]][keep] else .v
      }

      for (cii in seq_along(chunk)) {
        ci <- chunk[cii]
        if (!valid[ci] || !finite[ci]) next
        pars <- pars_list[[ci]]
        idx  <- (cii - 1L) * n_sim * n_t + seq_len(n_sim * n_t)
        cp   <- matrix(vals[idx], nrow = n_sim, ncol = n_t, byrow = TRUE)
        if (anyNA(cp)) { finite[ci] <- FALSE; next }
        # phi is eta-independent, so this configuration's first simulated row is
        # representative of its whole block
        .ph <- if (is.null(.phi_all)) attr(cp, "phi") else
          matrix(.phi_all[idx], nrow = n_sim, ncol = n_t, byrow = TRUE)[1L, ]
        ar <- .admUnitResidRows(pinfo, ov, pars$sigma_var, n_t, phi = .ph)
        # SAME gate as .admNLL(). The fused kernels implement forms 0/1/2 only and
        # have no off-diagonal channel, so a TBS/count/beta/ordinal/ar model fell
        # into adm_apply_residual's `else` branch and was scored as combined2.
        # This function IS the post-fit Hessian evaluator (covMethod = "r"), so
        # without the gate every standard error and RSE for those models came from
        # a different objective than the one that was fitted. Measured on a boxCox
        # model: .admNLL 190.28 vs .admNLLBatch 49.46.
        nll_ci <- if (.admResidCppOK(ar)) {
          if (identical(s$method, "var"))
            nll_var_from_samples_cpp(cp, as.numeric(s$E), s$v_diag,
                                     s$n, ar$form, ar$a2, ar$b2, ar$cc)
          else
            nll_cov_from_samples_cpp(cp, as.numeric(s$E), s$V,
                                     s$n, ar$form, ar$a2, ar$b2, ar$cc)
        } else {
          mu_s <- colMeans(cp)
          cpc  <- sweep(cp, 2L, mu_s)
          Vs   <- crossprod(cpc) / nrow(cp)
          ap   <- .admResidApply(mu_s, diag(Vs), ar, s$times, Vs)
          if (identical(s$method, "var")) {
            nll_var_cpp(as.numeric(s$E), s$v_diag, ap$mu, ap$dv, s$n)
          } else {
            Vp <- .admApplyResidTail(Vs, ap)
            nll_cov_cpp(as.numeric(s$E), s$V, ap$mu, Vp, s$n)
          }
        }
        if (is.finite(nll_ci)) nlls[ci] <- nlls[ci] + nll_ci
        else                    finite[ci] <- FALSE
      }
    }
  }

  nlls[!valid | !finite] <- Inf
  nlls
}

# -- Batched gradient evaluation -----------------------------------------------
# Evaluates the gradient for a list of parameter vectors via batched rxSolve
# calls (one per study). Returns a (n_c x np) gradient matrix.
# Used by .admCalcCov (use_grad=TRUE) to compute the Hessian via forward FD
# of the gradient -- all np+1 configs packed into a single rxSolve call per study.
.admGradBatch <- function(p_list, pinfo, studies, z_list, rxMod, output_var,
                           params_list, cores, h, sensModel = NULL,
                           use_central = FALSE) {
  n_c   <- length(p_list)
  if (n_c == 0L) return(matrix(0, 0, length(p_list[[1L]])))

  # Joint (same-subject) units are not handled by the batched path -- its stacked
  # matrices are shaped for a single output and it errors with "non-conformable
  # arguments". .admNLLBatch() already falls back per config for these; mirror that
  # here rather than relying on the driver's `!any_joint` guard staying in place.
  if (any(vapply(studies, function(u) isTRUE(u$is_joint), logical(1))))
    return(t(vapply(p_list, function(.p)
      .admGrad(.p, pinfo, studies, z_list, rxMod, output_var, params_list, cores,
               h, sensModel),
      numeric(length(p_list[[1L]])))))

  np            <- length(p_list[[1L]])
  n_eta         <- pinfo$n_eta
  n_s           <- length(pinfo$struct_names)
  n_e           <- length(pinfo$sigma_names)
  eta_col_names <- pinfo$eta_col_names

  unpaired_k <- which(vapply(pinfo$struct_names, function(nm)
    is.null(pinfo$struct_has_eta) || !isTRUE(pinfo$struct_has_eta[nm]), logical(1)))
  n_unp <- length(unpaired_k)

  pars_list <- lapply(p_list, function(p)
    tryCatch(.admUnpack(p, pinfo), error = function(e) NULL))
  valid <- !vapply(pars_list, is.null, logical(1))

  grad_acc <- matrix(0, nrow = n_c, ncol = np,
                     dimnames = list(NULL, names(p_list[[1L]])))

  for (si in seq_along(studies)) {
    s     <- studies[[si]]
    z     <- z_list[[si]]
    n_sim <- nrow(z)
    n_t   <- length(s$times)
    col_nms <- colnames(params_list[[si]])
    n_cols  <- length(col_nms)

    # This unit's output. .admGrad restricts residual error to this output's own
    # sigma(s); .admGradBatch did not -- it read `output_var` (the model-level
    # default), so on a multi-output fit the gradient-FD Hessian applied every
    # endpoint's residual error to every block. The residual spec is now looked up
    # per output (.admResidRows below), so the two paths cannot drift apart again.
    ovb <- s$output %||% output_var

    eta_mats <- lapply(seq_len(n_c), function(ci) {
      if (!valid[ci]) return(NULL)
      em <- z %*% t(pars_list[[ci]]$L)
      colnames(em) <- eta_col_names
      em
    })

    cp_mats     <- vector("list", n_c)
    dpred_lists <- vector("list", n_c)
    dtheta_lists <- vector("list", n_c)

    # --- Sens model path -------------------------------------------------------
    use_sens <- !is.null(sensModel) && n_eta > 0L
    if (use_sens) {
      rmap      <- sensModel$rename_map
      all_src   <- c(pinfo$struct_names, pinfo$sigma_names, eta_col_names)
      inner_nms <- rmap[all_src]; inner_nms <- inner_nms[!is.na(inner_nms)]

      inner_df <- as.data.frame(matrix(0, nrow = n_c * n_sim,
                                        ncol = length(inner_nms),
                                        dimnames = list(NULL, unname(inner_nms))))
      # An ESTIMATED boxCox/yeoJohnson lambda is a SIGMA name, so the zero-fill of
      # inner_df hands the solve lambda = 0 -- a plain log transform -- while the
      # back-transform below would invert with pred_tbs$lam, the model's STARTING
      # lambda, held constant across every configuration. Two different transforms:
      # the same mismatch .admSimulateSens documents as making the sens gradient
      # ~60x wrong for boxCox and NaN for yeoJohnson. This function IS the
      # covMethod = "r" Hessian evaluator, so it would put that error into every
      # reported standard error -- and lambda's own row would be insensitive to
      # lambda. Each configuration carries its OWN lambda, so it is written per
      # block of rows and inverted with that same number below.
      # Inlined, not factored out -- see the dev-mode daemon note in simulate.R.
      .tb0    <- sensModel$pred_tbs
      .lam_nm <- if (is.null(.tb0)) NA_character_ else .tb0$lam_name %||% NA_character_
      .lam_ci <- rep(if (is.null(.tb0)) NA_real_ else .tb0$lam, n_c)
      .lam_mapped <- if (!is.na(.lam_nm)) rmap[.lam_nm] else NA_character_

      for (ci in seq_len(n_c)) {
        if (!valid[ci]) next
        rows <- (ci - 1L) * n_sim + seq_len(n_sim)
        pars <- pars_list[[ci]]; eta <- eta_mats[[ci]]
        for (nm in pinfo$struct_names) {
          mapped <- rmap[nm]; if (!is.na(mapped)) inner_df[rows, mapped] <- pars$struct[nm]
        }
        for (j in seq_along(eta_col_names)) {
          mapped <- rmap[eta_col_names[j]]
          if (!is.na(mapped)) inner_df[rows, mapped] <- eta[, j]
        }
        if (!is.na(.lam_nm) && .lam_nm %in% names(pars$sigma_var)) {
          .lam_ci[ci] <- unname(pars$sigma_var[[.lam_nm]])
          if (!is.na(.lam_mapped) && .lam_mapped %in% names(inner_df))
            inner_df[rows, .lam_mapped] <- .lam_ci[ci]
        }
      }
      # fixed thetas: constants the loops above never write (see the note at the
      # top of simulate.R -- inlined on purpose; a new helper would be missing in a
      # dev-mode daemon, which cannot ADD bindings to the installed namespace)
      for (nm in names(sensModel$fixed_theta))
        inner_df[[nm]] <- rep(unname(sensModel$fixed_theta[[nm]]), nrow(inner_df))
      # do.call + sensModel$solve_args: DDE sensitivity solves are forced onto pure
      # dop853 (see .admLoadSensModel); NULL, hence a no-op, for every other model.
      out <- tryCatch(
        suppressWarnings(
          do.call(rxode2::rxSolve,
                  c(list(sensModel$mod, params = inner_df,
                         events = s$ev_full, cores = cores,
                         nDisplayProgress = pinfo$nDisplayProgress),
                    sensModel$solve_args))),
        error = function(e) NULL)
      if (is.null(out) || !all(sensModel$sens_cols %in% names(out))) {
        use_sens <- FALSE
      } else {
        keep      <- out[["time"]] %in% s$times
        vals_pred <- out[["rx_pred_"]][keep]
        vals_sens <- lapply(sensModel$sens_cols, function(col) out[[col]][keep])
        # d(pred)/d(theta) columns (augmented sens model); NULL keeps the FD path
        tsc       <- sensModel$theta_sens_cols
        vals_th   <- if (!is.null(tsc) && all(tsc %in% names(out)))
          lapply(tsc, function(col) out[[col]][keep]) else NULL
        # lnorm endpoint: rx_pred_ is log(f). Back-transform here too (chain rule),
        # or this batched path would feed .admGradBatch a log-scale prediction while
        # the NLL scores the natural scale. Inlined -- see the note in simulate.R.
        .tb <- sensModel$pred_tbs; .plog <- !is.null(.tb)
        for (ci in seq_len(n_c)) {
          if (!valid[ci]) next
          idx <- (ci - 1L) * n_sim * n_t + seq_len(n_sim * n_t)
          cpm <- matrix(vals_pred[idx], nrow = n_sim, ncol = n_t, byrow = TRUE)
          # this configuration's lambda -- the one written into the solve above
          .lm <- if (.plog) .lam_ci[ci] else NA_real_
          .gp <- if (.plog) .admTBSid(cpm, .lm, .tb$yj, .tb$lo, .tb$hi) else NULL
          if (.plog) cpm <- .admTBSi(cpm, .lm, .tb$yj, .tb$lo, .tb$hi)
          cp_mats[[ci]]     <- cpm
          dpred_lists[[ci]] <- lapply(vals_sens, function(vs) {
            D <- matrix(vs[idx], nrow = n_sim, ncol = n_t, byrow = TRUE)
            if (.plog) D * .gp else D
          })
          if (!is.null(vals_th))
            dtheta_lists[[ci]] <- stats::setNames(
              lapply(vals_th, function(vs) {
                D <- matrix(vs[idx], nrow = n_sim, ncol = n_t, byrow = TRUE)
                if (.plog) D * .gp else D
              }),
              names(tsc))
        }
      }
    }

    # --- Forward/central FD / no-eta fallback ----------------------------------
    if (!use_sens) {
      if (n_eta == 0L) {
        pdf_mat <- matrix(0, nrow = n_c * n_sim, ncol = n_cols,
                          dimnames = list(NULL, col_nms))
        pdf_mat[, grep("^rxerr", colnames(pdf_mat), value = TRUE)] <- 1L
        for (ci in seq_len(n_c)) {
          if (!valid[ci]) next
          rows <- (ci - 1L) * n_sim + seq_len(n_sim)
          pars <- pars_list[[ci]]
          for (nm in names(pars$struct)) pdf_mat[rows, nm] <- pars$struct[nm]
        }
        out <- tryCatch(rxode2::rxSolve(rxMod, params = as.data.frame(pdf_mat),
                                         events = s$ev_full, cores = cores,
                                         nDisplayProgress = pinfo$nDisplayProgress),
                        error = function(e) NULL)
        if (is.null(out)) { valid[] <- FALSE } else {
          keep <- out[["time"]] %in% s$times
          vals <- out[[ovb]][keep]
          if (is.null(vals)) vals <- out[["ipredSim"]][keep]
          for (ci in seq_len(n_c)) {
            if (!valid[ci]) next
            idx <- (ci - 1L) * n_sim * n_t + seq_len(n_sim * n_t)
            cp_mats[[ci]]     <- matrix(vals[idx], nrow = n_sim, ncol = n_t, byrow = TRUE)
            dpred_lists[[ci]] <- list()
          }
        }
      } else {
        # central FD: [base, eta1_hi, eta1_lo, ..., etaN_hi, etaN_lo] per config.
        # forward FD: [base, eta1_hi, eta2_hi, ..., etaN_hi] per config.
        n_blk   <- if (use_central) 1L + 2L * n_eta else 1L + n_eta
        pdf_mat <- matrix(0, nrow = n_c * n_blk * n_sim, ncol = n_cols,
                          dimnames = list(NULL, col_nms))
        pdf_mat[, grep("^rxerr", colnames(pdf_mat), value = TRUE)] <- 1L
        for (ci in seq_len(n_c)) {
          if (!valid[ci]) next
          pars     <- pars_list[[ci]]; eta <- eta_mats[[ci]]
          cfg_base <- (ci - 1L) * n_blk * n_sim
          rows_b   <- cfg_base + seq_len(n_sim)
          for (nm in names(pars$struct)) pdf_mat[rows_b, nm] <- pars$struct[nm]
          pdf_mat[rows_b, eta_col_names] <- eta
          if (use_central) {
            for (j in seq_len(n_eta)) {
              rows_hi <- cfg_base + n_sim * (2L*j - 1L) + seq_len(n_sim)
              rows_lo <- cfg_base + n_sim * (2L*j)      + seq_len(n_sim)
              eta_hi  <- eta; eta_hi[, j] <- eta_hi[, j] + h
              eta_lo  <- eta; eta_lo[, j] <- eta_lo[, j] - h
              for (nm in names(pars$struct)) {
                pdf_mat[rows_hi, nm] <- pars$struct[nm]
                pdf_mat[rows_lo, nm] <- pars$struct[nm]
              }
              pdf_mat[rows_hi, eta_col_names] <- eta_hi
              pdf_mat[rows_lo, eta_col_names] <- eta_lo
            }
          } else {
            for (j in seq_len(n_eta)) {
              rows_hi <- cfg_base + n_sim * j + seq_len(n_sim)
              eta_hi  <- eta; eta_hi[, j] <- eta_hi[, j] + h
              for (nm in names(pars$struct)) pdf_mat[rows_hi, nm] <- pars$struct[nm]
              pdf_mat[rows_hi, eta_col_names] <- eta_hi
            }
          }
        }
        out <- tryCatch(rxode2::rxSolve(rxMod, params = as.data.frame(pdf_mat),
                                         events = s$ev_full, cores = cores,
                                         nDisplayProgress = pinfo$nDisplayProgress),
                        error = function(e) NULL)
        if (is.null(out)) { valid[] <- FALSE } else {
          keep <- out[["time"]] %in% s$times
          vals <- out[[ovb]][keep]
          if (is.null(vals)) vals <- out[["ipredSim"]][keep]
          for (ci in seq_len(n_c)) {
            if (!valid[ci]) next
            cfg_out_base <- (ci - 1L) * n_blk * n_sim * n_t
            cp_mats[[ci]] <- matrix(vals[cfg_out_base + seq_len(n_sim * n_t)],
                                    nrow = n_sim, ncol = n_t, byrow = TRUE)
            dpred_lists[[ci]] <- if (use_central) {
              lapply(seq_len(n_eta), function(j) {
                off_hi <- cfg_out_base + n_sim * n_t * (2L*j - 1L)
                off_lo <- cfg_out_base + n_sim * n_t * (2L*j)
                (matrix(vals[off_hi + seq_len(n_sim * n_t)], nrow = n_sim, ncol = n_t, byrow = TRUE) -
                 matrix(vals[off_lo + seq_len(n_sim * n_t)], nrow = n_sim, ncol = n_t, byrow = TRUE)) / (2 * h)
              })
            } else {
              lapply(seq_len(n_eta), function(j) {
                off_hi <- cfg_out_base + n_sim * n_t * j
                (matrix(vals[off_hi + seq_len(n_sim * n_t)], nrow = n_sim, ncol = n_t, byrow = TRUE) -
                 cp_mats[[ci]]) / h
              })
            }
          }
        }
      }
    }

    # --- Unpaired theta FD across all valid configs ---------------------------
    # All (ci x bi) configs stacked into one rxSolve call per pass.
    # central FD: two passes (hi + lo); forward FD: hi only, base reused.
    cp_hi_store <- vector("list", n_c)
    cp_lo_store <- if (use_central) vector("list", n_c) else NULL
    # The augmented sens model already returned d(pred)/d(theta) for every config
    # in the solve above -- skip the FD passes entirely (they are the expensive
    # part of this function: one extra rxSolve over all n_c x n_unp configs, two
    # with central differences).
    have_theta_sens <- n_unp > 0L &&
      any(!vapply(dtheta_lists, is.null, logical(1)))
    if (n_unp > 0L && !have_theta_sens) {
      for (ci in seq_len(n_c)) {
        cp_hi_store[[ci]] <- vector("list", n_unp)
        if (use_central) cp_lo_store[[ci]] <- vector("list", n_unp)
      }
      cu_idx <- expand.grid(bi = seq_len(n_unp), ci = seq_len(n_c))
      n_cu   <- nrow(cu_idx)
      pdf_hi <- matrix(0, nrow = n_cu * n_sim, ncol = n_cols,
                       dimnames = list(NULL, col_nms))
      pdf_hi[, grep("^rxerr", colnames(pdf_hi), value = TRUE)] <- 1L
      for (cuki in seq_len(n_cu)) {
        ci <- cu_idx$ci[cuki]; bi <- cu_idx$bi[cuki]
        if (!valid[ci] || is.null(cp_mats[[ci]])) next
        rows <- (cuki - 1L) * n_sim + seq_len(n_sim)
        pars <- pars_list[[ci]]; eta <- eta_mats[[ci]]
        for (nm in names(pars$struct))       pdf_hi[rows, nm]               <- pars$struct[nm]
        for (j  in seq_along(eta_col_names)) pdf_hi[rows, eta_col_names[j]] <- eta[, j]
        nm_u <- pinfo$struct_names[unpaired_k[bi]]
        pdf_hi[rows, nm_u] <- pars$struct[nm_u] + h
      }
      out_hi <- tryCatch(rxode2::rxSolve(rxMod, params = as.data.frame(pdf_hi),
                                          events = s$ev_full, cores = cores,
                                          nDisplayProgress = pinfo$nDisplayProgress),
                         error = function(e) NULL)
      if (!is.null(out_hi)) {
        vals_hi <- out_hi[[ovb]][out_hi[["time"]] %in% s$times]
        if (is.null(vals_hi)) vals_hi <- out_hi[["ipredSim"]][out_hi[["time"]] %in% s$times]
        for (cuki in seq_len(n_cu)) {
          ci <- cu_idx$ci[cuki]; bi <- cu_idx$bi[cuki]
          if (!valid[ci] || is.null(cp_mats[[ci]])) next
          idx <- (cuki - 1L) * n_sim * n_t + seq_len(n_sim * n_t)
          cp_hi_store[[ci]][[bi]] <- matrix(vals_hi[idx], nrow = n_sim, ncol = n_t, byrow = TRUE)
        }
      }

      if (use_central) {
        pdf_lo <- pdf_hi
        for (cuki in seq_len(n_cu)) {
          ci <- cu_idx$ci[cuki]; bi <- cu_idx$bi[cuki]
          if (!valid[ci] || is.null(cp_mats[[ci]])) next
          rows <- (cuki - 1L) * n_sim + seq_len(n_sim)
          pars <- pars_list[[ci]]
          nm_u <- pinfo$struct_names[unpaired_k[bi]]
          pdf_lo[rows, nm_u] <- pars$struct[nm_u] - h
        }
        out_lo <- tryCatch(rxode2::rxSolve(rxMod, params = as.data.frame(pdf_lo),
                                            events = s$ev_full, cores = cores,
                                            nDisplayProgress = pinfo$nDisplayProgress),
                           error = function(e) NULL)
        if (!is.null(out_lo)) {
          vals_lo <- out_lo[[ovb]][out_lo[["time"]] %in% s$times]
          if (is.null(vals_lo)) vals_lo <- out_lo[["ipredSim"]][out_lo[["time"]] %in% s$times]
          for (cuki in seq_len(n_cu)) {
            ci <- cu_idx$ci[cuki]; bi <- cu_idx$bi[cuki]
            if (!valid[ci] || is.null(cp_mats[[ci]])) next
            idx <- (cuki - 1L) * n_sim * n_t + seq_len(n_sim * n_t)
            cp_lo_store[[ci]][[bi]] <- matrix(vals_lo[idx], nrow = n_sim, ncol = n_t, byrow = TRUE)
          }
        }
      }
    }

    # --- Gradient formula per config ------------------------------------------
    for (ci in seq_len(n_c)) {
      if (!valid[ci] || is.null(cp_mats[[ci]])) next
      cp_mat <- cp_mats[[ci]]
      if (anyNA(cp_mat)) { valid[ci] <- FALSE; next }
      dpred_list <- dpred_lists[[ci]]
      pars       <- pars_list[[ci]]
      eta_mat    <- eta_mats[[ci]]

      arr       <- .admResidRows(pinfo, ovb, pars$sigma_var, n_t)
      mu_struct <- colMeans(cp_mat)
      cp_c      <- sweep(cp_mat, 2L, mu_struct)

      is_var <- identical(s$method, "var")
      cov_f  <- NULL
      if (is_var) {
        var_f <- adm_col_sq_sum_cpp(cp_c) / n_sim
        ap <- .admResidApply(mu_struct, var_f, arr)
        mu <- ap$mu; pv <- ap$dv
        r  <- as.numeric(s$E) - mu
        dNLL_dmu     <- s$n * as.numeric(-2 * r / pv)
        dNLL_dV_diag <- s$n * (1 / pv - s$v_diag / pv^2 - r^2 / pv^2)
      } else {
        cov_f <- crossprod(cp_c) / n_sim               # keep the STRUCTURAL cov
        V  <- cov_f
        var_f <- diag(V)
        ap <- .admResidApply(mu_struct, var_f, arr, s$times, cov_f)
        mu <- ap$mu
        V  <- .admApplyResidTail(V, ap)
        r  <- as.numeric(s$E) - mu
        cholV <- tryCatch(chol(V), error = function(e) NULL)
        if (is.null(cholV)) { valid[ci] <- FALSE; next }
        invV         <- chol2inv(cholV)
        dNLL_dmu     <- s$n * as.numeric(-2 * invV %*% r)
        dNLL_dV      <- s$n * (invV - invV %*% (s$V + tcrossprod(r)) %*% invV)
        dNLL_dV_diag <- diag(dNLL_dV)
      }

      .dres <- .admResidDeriv(mu_struct, var_f, arr, pinfo)   # once, reused across terms
      sigma_mu_scale <- .admResidMuCoupling(mu_struct, arr, pinfo,
                                            dNLL_dV_diag, dNLL_dmu, var_f,
                                            if (is_var) NULL else dNLL_dV, cov_f,
                                            s$times, deriv = .dres)
      eff_dmu <- dNLL_dmu + sigma_mu_scale
      inv_n <- 1 / n_sim
      # V_pred -> V_struct chain (see .admGrad)
      vchain <- .admResidVChain(mu_struct, var_f, arr, pinfo, s$times, deriv = .dres)
      .dmv <- attr(vchain, "dmu_dv0") %||% numeric(n_t)
      dNLL_dV_diag_s <- dNLL_dV_diag * diag(vchain) + dNLL_dmu * .dmv
      if (!is_var) {
        dNLL_dV_s <- dNLL_dV * vchain
        diag(dNLL_dV_s) <- diag(dNLL_dV_s) + dNLL_dmu * .dmv
      }

      if (n_eta > 0L) {
        D_mat        <- do.call(cbind, dpred_list)
        eta_rows_df  <- pinfo$eta_rows_df
        z_diag_scale <- sweep(z, 2L, diag(pars$L) / 2, "*")
        neta1 <- as.integer(eta_rows_df$neta1)
        neta2 <- as.integer(eta_rows_df$neta2)
        go <- if (is_var)
          adm_grad_eta_omega_var_cpp(
            cp_c, D_mat, z_diag_scale, z,
            dNLL_dV_diag_s, dNLL_dmu, sigma_mu_scale,
            neta1, neta2,
            n_t, n_eta)
        else
          adm_grad_eta_omega_cpp(
            cp_c, D_mat, z_diag_scale, z,
            dNLL_dV_s, dNLL_dmu, sigma_mu_scale,
            neta1, neta2,
            n_t, n_eta)
        for (j in seq_len(n_eta)) {
          if (!is.null(pinfo$struct_eta_idx) && !is.na(pinfo$struct_eta_idx[j]))
            grad_acc[ci, pinfo$struct_eta_idx[j]] <- grad_acc[ci, pinfo$struct_eta_idx[j]] + go$eta_grad[j]
        }
        k_om <- n_s + n_e
        for (r_idx in seq_len(nrow(eta_rows_df))) {
          k_om <- k_om + 1L
          grad_acc[ci, k_om] <- grad_acc[ci, k_om] + go$omega_grad[r_idx]
        }
      }

      dth_ci <- dtheta_lists[[ci]]
      for (bi in seq_len(n_unp)) {
        k_s   <- unpaired_k[bi]
        dpred <- if (!is.null(dth_ci)) {
          dth_ci[[pinfo$struct_names[k_s]]]            # exact, from the sens solve
        } else {
          cp_hi_s <- cp_hi_store[[ci]][[bi]]
          if (is.null(cp_hi_s)) NULL
          else if (use_central && !is.null(cp_lo_store[[ci]][[bi]]))
            (cp_hi_s - cp_lo_store[[ci]][[bi]]) / (2 * h)
          else (cp_hi_s - cp_mat) / h
        }
        if (is.null(dpred)) next
        grad_acc[ci, k_s] <- grad_acc[ci, k_s] +
          if (is_var)
            adm_grad_partial_var_cpp(cp_c, dpred, dNLL_dV_diag_s, eff_dmu, inv_n)
          else
            adm_grad_partial_cpp(cp_c, dpred, dNLL_dV_s, eff_dmu, inv_n)
      }

      # .admSigmaGrad returns a full-length n_e vector indexed by GLOBAL sigma
      # position, zero for sigmas belonging to other endpoints. That is what keeps
      # a second endpoint's sigma gradient out of the first endpoint's slot.
      grad_acc[ci, n_s + seq_len(n_e)] <- grad_acc[ci, n_s + seq_len(n_e)] +
        .admSigmaGrad(mu_struct, arr, pinfo, dNLL_dV_diag, dNLL_dmu, var_f,
                    if (is_var) NULL else dNLL_dV, s$times, cov_f, deriv = .dres)
    }
  }

  grad_acc[!valid, ] <- NA_real_
  grad_acc
}

# -- Post-fit covariance (numerical Hessian, R method: 2*H^-1) -----------------

.admCalcCov <- function(p_hat, pinfo, studies, z_list, rxMod, output_var,
                        params_list, cores, cov_n_sim = NULL,
                        use_grad = FALSE, grad_h = 1e-4, cov_h = 1e-3,
                        cov_h_outer = .Machine$double.eps^(1/5),
                        sensModel = NULL, use_central = FALSE,
                        sampling = "sobol") {
  np    <- length(p_hat)
  nms   <- names(p_hat)

  # Hessian over struct + sigma + omega (falls back to struct+sigma if not PD).
  # Matches nlmixr2 FOCEI: omega entries are in the optimizer but skipped for cov.
  n_s     <- length(pinfo$struct_names)
  n_e     <- length(pinfo$sigma_names)
  n_o     <- length(pinfo$omega_par)
  # The Hessian spans struct + sigma + OMEGA -- see .adghCalcCov() for the
  # measurement. Excluding omega made the STRUCTURAL SEs too small (reported SE /
  # empirical sampling SD over simulated datasets went from 0.67 to 1.17 on prop
  # and 0.67 to 1.06 on lnorm when omega was put back), because a theta carrying
  # an eta is correlated with that eta's variance. Falls back to the struct+sigma
  # sub-block if the weakly-identified omega Cholesky makes the full H indefinite.
  n_sub   <- n_s + n_e
  cov_idx <- seq_len(n_sub + n_o)
  np_cov  <- length(cov_idx)
  nms_cov  <- nms[cov_idx]

  if (!is.null(cov_n_sim) && cov_n_sim != nrow(z_list[[1]])) {
    z_list      <- .admMakeZ(cov_n_sim, pinfo, length(studies), sampling)
    params_list <- .admMakeParamsList(cov_n_sim, pinfo, length(studies))
  }

  nll_fn <- function(p)
    suppressMessages(.admNLL(p, pinfo, studies, z_list, rxMod, output_var, params_list, cores))

  nll0 <- nll_fn(p_hat)
  if (!is.finite(nll0)) {
    warning("admCalcCov: NLL not finite at p_hat -- covariance not computed")
    return(NULL)
  }

  if (use_grad) {
    h_fwd    <- pmax(abs(p_hat[cov_idx]), 0.1) * cov_h_outer
    # Inner step: larger than grad_h to reduce gradient noise amplification.
    # Hessian FD divides by h_outer, so gradient noise is scaled up by 1/h_outer.
    h_inner  <- cov_h
    # np_cov+1 param vectors: p_hat followed by np_cov forward-perturbed versions
    # (only struct+sigma entries perturbed; omega stays fixed at p_hat).
    p_list <- c(list(p_hat), lapply(seq_len(np_cov), function(jj) {
      ph <- p_hat; ph[cov_idx[jj]] <- ph[cov_idx[jj]] + h_fwd[jj]; ph
    }))
    grads <- .admGradBatch(p_list, pinfo, studies, z_list, rxMod, output_var,
                            params_list, cores, h_inner, sensModel,
                            use_central = use_central)
    g0 <- grads[1L, cov_idx]
    H  <- matrix(0, np_cov, np_cov, dimnames = list(nms_cov, nms_cov))
    for (jj in seq_len(np_cov)) {
      gj     <- grads[jj + 1L, cov_idx]
      H[, jj] <- if (anyNA(gj)) 0 else (gj - g0) / h_fwd[jj]
    }
    H <- (H + t(H)) / 2
  } else {
    h_gill  <- pmax(abs(p_hat[cov_idx]), 0.1) * cov_h_outer
    n_off   <- np_cov * (np_cov - 1L) / 2L

    # Perturb only struct+sigma entries; omega stays fixed at p_hat.
    diag_p <- vector("list", 2L * np_cov)
    for (k in seq_len(np_cov)) {
      ki <- cov_idx[k]
      ph <- p_hat; ph[ki] <- ph[ki] + h_gill[k]; diag_p[[2L*k - 1L]] <- ph
      pl <- p_hat; pl[ki] <- pl[ki] - h_gill[k]; diag_p[[2L*k]]      <- pl
    }
    off_p  <- vector("list", 4L * n_off)
    off_ij <- matrix(0L, n_off, 2L)
    oci    <- 0L
    for (i in seq_len(np_cov - 1L)) {
      for (j in seq(i + 1L, np_cov)) {
        oci <- oci + 1L; off_ij[oci, ] <- c(i, j)
        ii <- cov_idx[i]; ji <- cov_idx[j]
        hi <- h_gill[i];  hj <- h_gill[j]
        p_pp <- p_hat; p_pp[ii] <- p_pp[ii] + hi; p_pp[ji] <- p_pp[ji] + hj
        p_pm <- p_hat; p_pm[ii] <- p_pm[ii] + hi; p_pm[ji] <- p_pm[ji] - hj
        p_mp <- p_hat; p_mp[ii] <- p_mp[ii] - hi; p_mp[ji] <- p_mp[ji] + hj
        p_mm <- p_hat; p_mm[ii] <- p_mm[ii] - hi; p_mm[ji] <- p_mm[ji] - hj
        off_p[[(oci-1L)*4L + 1L]] <- p_pp; off_p[[(oci-1L)*4L + 2L]] <- p_pm
        off_p[[(oci-1L)*4L + 3L]] <- p_mp; off_p[[(oci-1L)*4L + 4L]] <- p_mm
      }
    }
    all_p   <- c(diag_p, off_p)
    nll_all <- .admNLLBatch(all_p, pinfo, studies, z_list, rxMod, output_var,
                             params_list, cores)

    H <- matrix(0, np_cov, np_cov, dimnames = list(nms_cov, nms_cov))
    for (k in seq_len(np_cov)) {
      hk <- h_gill[k]
      H[k, k] <- (nll_all[2L*k - 1L] - 2*nll0 + nll_all[2L*k]) / hk^2
    }
    for (oci in seq_len(n_off)) {
      i <- off_ij[oci, 1L]; j <- off_ij[oci, 2L]
      hi <- h_gill[i]; hj <- h_gill[j]
      base <- 2L * np_cov + (oci - 1L) * 4L
      H[i, j] <- H[j, i] <-
        (nll_all[base + 1L] - nll_all[base + 2L] -
         nll_all[base + 3L] + nll_all[base + 4L]) / (4 * hi * hj)
    }
  }

  if (!all(is.finite(H))) {
    warning("admCalcCov: Hessian has non-finite entries -- covariance not computed")
    return(NULL)
  }

  eig_dec <- tryCatch(eigen(H, symmetric = TRUE), error = function(e) NULL)
  H_eigs  <- if (!is.null(eig_dec)) eig_dec$values else rep(NA_real_, np_cov)

  # If the weakly-identified omega Cholesky makes the full Hessian indefinite,
  # drop back to the struct+sigma sub-block rather than reporting nothing.
  .red <- .admReduceNpdOmega(H, H_eigs, eig_dec, nms_cov, n_o, n_sub)
  if (.red$reduced)
    warning("admCalcCov: the full Hessian including omega was not positive ",
            "definite; reporting structural and sigma standard errors only.",
            call. = FALSE)
  H <- .red$H; nms_cov <- .red$nms_cov; np_cov <- .red$np_cov
  eig_dec <- .red$eig_dec; H_eigs <- .red$H_eigs

  if (!is.null(eig_dec) && min(H_eigs) < 0) {
    hint <- if (use_grad)
      sprintf("Try increasing cov_h_outer (currently %.3e) or cov_h (currently %.3e) in admControl(), e.g. cov_h_outer = %.3e.",
              cov_h_outer, cov_h, cov_h_outer * 4)
    else
      sprintf("Try increasing cov_h_outer (currently %.3e) in admControl(), e.g. cov_h_outer = %.3e.",
              cov_h_outer, cov_h_outer * 4)
    warning(sprintf(
      "admCalcCov: Hessian not positive definite (min eigenvalue %.3e). Covariance not computed. %s",
      min(H_eigs), hint), call. = FALSE)
    return(NULL)
  }

  inv_method <- "chol"
  Hinv <- tryCatch(
    chol2inv(chol(H)),
    error = function(e) {
      inv_method <<- "solve"
      tryCatch(
        solve(H),
        error = function(e2) {
          inv_method <<- "sqrtm"
          tryCatch({
            if (!requireNamespace("expm", quietly = TRUE))
              stop("expm package needed for sqrtm fallback")
            solve(expm::sqrtm(H %*% t(H)))
          }, error = function(e3) { inv_method <<- "failed"; NULL })
        }
      )
    }
  )
  if (is.null(Hinv)) {
    warning("admCalcCov: Hessian inversion failed -- covariance not computed")
    return(NULL)
  }

  cov_full <- (2 * Hinv + t(2 * Hinv)) / 2
  dimnames(cov_full) <- list(nms_cov, nms_cov)
  # Rotate onto the reported scale (residual delta factors + omega Jacobian). One
  # shared implementation for all three estimators -- see .admScaleReportedCov().
  .admScaleReportedCov(cov_full, p_hat, pinfo, n_s, n_e, n_o, n_sub)
}

# -- Restart worker ------------------------------------------------------------

.admRestartWorker <- function(restart_id, p_init, ui_lstExpr, pinfo,
                              ov_lower, ov_upper, scale_c = NULL, studies, n_sim,
                              seed, algorithm, ftol_rel, maxeval,
                              use_grad, grad_h, grad_bounds,
                              output_var = "cp",
                              sampling = "sobol",
                              use_central = FALSE,
                              print_progress = TRUE, print = 10L,
                              cores = NULL, no_lock = FALSE,
                              sens_cache_file = NULL, sens_cols = NULL, sens_rename = NULL,
                              rxMod_direct = NULL, sensModel_direct = NULL) {
  library(admixr2)

  # Dev mode: patch the installed namespace with any dev functions in .GlobalEnv.
  # A daemon is patched by .admDaemonRestart() before it gets here; this covers a
  # direct (sequential) call. tryCatch guards against the installed package
  # predating this function (run devtools::install() once).
  tryCatch(.admPatchDevNamespace(), error = function(e) NULL)

  m <- .admWorkerLoadModels(ui_lstExpr, rxMod_direct, cores,
                            sens_cache_file, sens_cols, sens_rename, sensModel_direct,
                            pinfo)

  set.seed(seed)
  z_list      <- .admMakeZ(n_sim, pinfo, length(studies), sampling)
  set.seed(seed + restart_id)
  params_list <- .admMakeParamsList(n_sim, pinfo, length(studies))

  nll_fn <- function(p)
    .admNLL(p, pinfo, studies, z_list, m$rxMod, output_var, params_list, m$cores_w)

  # Derived here (not a worker argument) so the worker signature stays stable for
  # the dev-mode worker path. A joint unit's stacked-MVN gradient is analytical
  # when the sens model is available, else FD of the aggregate NLL.
  .any_joint <- any(vapply(studies, function(u) isTRUE(u$is_joint), logical(1)))
  grad_fn <- if (.any_joint && is.null(m$sensModel))
    function(p)
      .admNLLGradFD(p, pinfo, studies, z_list, m$rxMod, output_var,
                    params_list, m$cores_w, grad_h, use_central = use_central)
  else
    function(p)
      .admGrad(p, pinfo, studies, z_list, m$rxMod, output_var,
               params_list, m$cores_w, grad_h, m$sensModel, use_central = use_central)

  # Lock only when the model was loaded from cache in this process; never across
  # workers (each holds its own independent model instance -> no_lock).
  lock_rxMod <- if (is.null(rxMod_direct) && !no_lock) m$rxMod else NULL

  .admScaledOptimize(restart_id, p_init, ov_lower, ov_upper, scale_c,
                     use_grad, grad_bounds, algorithm, ftol_rel, maxeval,
                     nll_fn, grad_fn, pinfo, print_progress, print,
                     lock_rxMod = lock_rxMod)
}

# -- Multi-restart orchestration -----------------------------------------------

# Parallel restarts run on a pool of mirai daemons: background R processes that
# behave identically on every platform, so there is exactly one worker code path
# (no fork/PSOCK split). Daemons never share the parent's memory, so everything
# a restart needs is serialised to it; compiled model DLLs cannot cross that
# boundary and are reloaded from the qs2 cache inside the daemon.
#
# The pool lives on its own mirai compute profile ("admixr2") so that starting
# and stopping it never disturbs daemons the user set up for their own code.
.adm_compute <- "admixr2"

.adm_worker_env <- new.env(parent = emptyenv())
.adm_worker_env$n <- 0L   # number of daemons currently running

# Clean up on R exit (onexit = TRUE) or when env is GC'd
reg.finalizer(.adm_worker_env, function(e) {
  if (e$n > 0L)
    tryCatch(mirai::daemons(0L, .compute = .adm_compute), error = function(err) NULL)
}, onexit = TRUE)

#' Stop parallel workers
#'
#' Stops any worker processes (mirai daemons) started by a parallel-restart fit
#' (`admControl(workers = N)`). Workers are stopped automatically after the
#' restart phase completes, so this function is only needed if a fit was
#' interrupted before cleanup could run.
#'
#' @return `NULL`, invisibly.
#'
#' @examples
#' # Safe to call at any time; no-op if no workers are running
#' admStopWorkers()
#'
#' @export
admStopWorkers <- function() {
  n <- .admStopDaemons()
  if (n == 0L)
    message("No admixr2 workers running.")
  else
    message(sprintf("%d admixr2 worker(s) stopped.", n))
  invisible(NULL)
}

# Silent counterpart of admStopWorkers() for internal use: returns the number of
# daemons that were shut down (0 when the pool was already empty).
.admStopDaemons <- function() {
  n <- .adm_worker_env$n
  if (n == 0L) return(invisible(0L))
  tryCatch(mirai::daemons(0L, .compute = .adm_compute), error = function(e) NULL)
  .adm_worker_env$n <- 0L
  invisible(n)
}

# Start the daemon pool for a multi-restart fit. The caller stops it again via
# .admStopDaemons() once the restarts are done, so all cores are free for the
# Hessian step, and registers an on.exit() so an interrupted fit leaves no
# orphaned processes behind.
.admSetupDaemons <- function(.ctl, n_r) {
  if (.ctl$workers <= 1L) return(invisible(NULL))
  if (!requireNamespace("mirai", quietly = TRUE))
    stop("Package 'mirai' required for workers > 1. Install it with install.packages('mirai').",
         call. = FALSE)
  if (.ctl$workers > n_r)
    warning(sprintf(
      "admControl(workers=%d) exceeds n_restarts=%d: only %d worker process(es) will be started",
      .ctl$workers, n_r, n_r
    ), call. = FALSE)

  n_w <- min(.ctl$workers, n_r)
  # Reset any pool left behind by an interrupted fit before starting a fresh one.
  .admStopDaemons()
  message(sprintf("  Starting %d worker(s)", n_w))
  mirai::daemons(n_w, .compute = .adm_compute)
  .adm_worker_env$n <- n_w
  invisible(n_w)
}

# Entry point evaluated inside a daemon: attach admixr2, patch in any dev-mode
# functions serialised from the parent (empty list in installed mode), then run
# one restart. Keeping this here -- rather than adding arguments to the restart
# workers -- is what lets the workers keep stable signatures: a daemon resolves
# them from the *installed* namespace, so a new worker argument would throw
# `unused argument` before the patched dev version could ever be reached.
.admDaemonRestart <- function(r, worker_fn_name, fn_list, inits, all_args,
                              cores_vec, effective_workers) {
  library(admixr2)
  .adm_ns <- asNamespace("admixr2")
  for (.nm in names(fn_list))
    tryCatch(utils::assignInNamespace(.nm, fn_list[[.nm]], ns = .adm_ns),
             error = function(e) NULL)
  wfn <- get(worker_fn_name, envir = .adm_ns, inherits = FALSE)
  args <- all_args
  args$cores <- cores_vec[[(r - 1L) %% effective_workers + 1L]]
  do.call(wfn, c(list(restart_id = r, p_init = inits[[r]]), args))
}

# Patch the installed admixr2 namespace with any dev-mode functions found in
# .GlobalEnv. Retained for direct (non-daemon) worker calls; daemons are patched
# by .admDaemonRestart() from the serialised fn_list instead. No-op when
# .GlobalEnv has no matching dev functions (installed mode).
.admPatchDevNamespace <- function() {
  .adm_dev_nms <- ls(envir = .GlobalEnv, all.names = TRUE,
                     pattern = "^\\.(adm|adfo|adirmc|softmax|logdmvnorm)")
  if (length(.adm_dev_nms) == 0L) return(invisible(NULL))
  .adm_ns <- asNamespace("admixr2")
  for (.nm in .adm_dev_nms) {
    .fn <- get(.nm, envir = .GlobalEnv, inherits = FALSE)
    if (is.function(.fn))
      tryCatch(utils::assignInNamespace(.nm, .fn, ns = .adm_ns),
               error = function(e) NULL)
  }
  invisible(length(.adm_dev_nms))
}

# -- Shared restart-worker helpers ---------------------------------------------
#
# The four estimator restart workers (.adfoRestartWorker, .admRestartWorker,
# .adghRestartWorker, .adirmcRestartWorker) shared two near-identical blocks:
# model/sens-model loading from cache, and -- for the three single-nloptr
# estimators (adfo/admc/adgh) -- the scaled box-constrained optimisation loop
# with progress tracking. Both are factored out here so each worker only
# supplies its estimator-specific NLL/gradient closures.

# Resolve worker cores, load the simulation model (direct or from qs2 cache),
# and (optionally) the sensitivity model. Returns a list(cores_w, rxMod,
# sensModel). adirmc passes no sens_* args -> sensModel is NULL (unused).
.admWorkerLoadModels <- function(ui_lstExpr, rxMod_direct = NULL, cores = NULL,
                                 sens_cache_file = NULL, sens_cols = NULL,
                                 sens_rename = NULL, sensModel_direct = NULL,
                                 pinfo = NULL) {
  cores_w <- if (!is.null(cores)) {
    cores
  } else if (!is.null(rxMod_direct)) {
    max(1L, parallel::detectCores() - 1L)
  } else {
    1L
  }

  if (!is.null(rxMod_direct)) {
    rxMod <- rxMod_direct
  } else {
    .cacheFile <- file.path(rxode2::rxTempDir(),
                            paste0("adm-sim-", digest::digest(ui_lstExpr), ".qs2"))
    rxMod <- qs2::qs_read(.cacheFile)
    rxode2::rxLoad(rxMod)
  }

  # The cache file holds the full sens result list (type/mod/sens_cols/...), so
  # the object to rxLoad() is m$mod; rxLoad(m) errors and lands in the NULL
  # branch, which is what silently forced workers onto grad = "fd".
  sensModel <- if (!is.null(sensModel_direct)) {
    sensModel_direct
  } else if (!is.null(sens_cache_file) && file.exists(sens_cache_file)) {
    tryCatch({
      m <- qs2::qs_read(sens_cache_file)
      rxode2::rxLoad(m$mod)
      # PREFER the parent's values over whatever is in the file. The worker cannot
      # re-derive these (it has no ui), so a cache written by an older admixr2 --
      # with the position-indexed rename_map, which puts a theta's value in the
      # wrong THETA[k] slot -- would otherwise be used verbatim and the parallel
      # fit would silently disagree with the sequential one. (The cache key now
      # carries a schema tag too, so such a file is no longer even a hit; this is
      # the belt to that pair of braces.)
      if (!is.null(sens_cols))   m$sens_cols  <- sens_cols
      if (!is.null(sens_rename)) m$rename_map <- sens_rename
      # m$theta_sens_cols and m$fixed_theta are NOT overwritten here -- the parent
      # does not thread them through (adding a worker-load argument would trip the
      # dev-mode stale-daemon `unused argument` trap; see .admRestartWorker's note).
      # Their staleness is guarded solely by the "dirs-jump+fixed-theta" schema tag
      # in the sens cache key: any change to how those fields are DERIVED must bump
      # that tag, or a stale file would become a false hit and a worker could fill
      # the wrong constant into a fixed theta's THETA[k] column, silently diverging
      # from the sequential fit. Overwriting sens_cols/rename_map above is the belt;
      # the schema tag is the braces for these two.
      #
      # pred_tbs MUST be re-derived, exactly as .admLoadSensModel() re-derives it
      # on the parent's cache-hit path. The cache key digests ui$lstExpr -- the
      # model({}) block only -- while lambda's starting value and its fix() status
      # live in ini({}), so `lam <- fix(0.5)` and `lam <- 0.5` COLLIDE on one key.
      # The parent overwrites the field; the worker did not, so a parallel restart
      # could invert the transform with a different lambda from the sequential fit
      # -- the same silent parent/worker divergence the block above exists to
      # prevent, and invisible because the NLL stays bit-identical.
      # Derived from `pinfo`, which the worker already holds, rather than from a
      # new worker ARGUMENT -- see .admRestartWorker's note on stale daemons.
      if (!is.null(pinfo) && !is.null(m$pred_tbs)) {
        .sp <- .admResidSpecs(pinfo)
        if (length(.sp)) {
          .s1  <- .sp[[1L]]
          .nat <- tryCatch(.admSigmaNat(pinfo$sigma_init, pinfo),
                           error = function(e) NULL)
          .est <- !is.null(.s1$k_lam) && !is.na(.s1$k_lam)
          m$pred_tbs <- list(
            lam      = if (.est && !is.null(.nat)) unname(.nat[[.s1$k_lam]])
                       else if (is.finite(.s1$lam_fixed %||% NA_real_)) .s1$lam_fixed
                       else m$pred_tbs$lam,
            yj       = .s1$yj %||% m$pred_tbs$yj,
            lam_name = if (.est) pinfo$sigma_names[[.s1$k_lam]] else NA_character_,
            lo       = if (is.finite(.s1$tr_lo %||% NA_real_)) .s1$tr_lo else 0,
            hi       = if (is.finite(.s1$tr_hi %||% NA_real_)) .s1$tr_hi else 1)
        }
      }
      m
    }, error = function(e) NULL)
  } else {
    NULL
  }

  list(cores_w = cores_w, rxMod = rxMod, sensModel = sensModel)
}

# Scaled, box-constrained single-nloptr optimisation with NLL/par trace
# tracking and live progress rows. Shared by the adfo/admc/adgh workers; each
# supplies `nll_fn(p)` and (optionally) `grad_fn(p)`. `lock_rxMod` non-NULL
# wraps the optimisation in rxLock/rxUnlock (used by the MC estimator when it
# loads the model from cache). Returns the standard restart-result list.
.admScaledOptimize <- function(restart_id, p_init, ov_lower, ov_upper, scale_c,
                               use_grad, grad_bounds, algorithm, ftol_rel, maxeval,
                               nll_fn, grad_fn, pinfo, print_progress, print,
                               lock_rxMod = NULL) {
  .iter      <- 0L
  .best_nll  <- Inf
  .nll_trace <- numeric(0)
  .par_trace <- NULL

  eval_f <- function(p) {
    .iter <<- .iter + 1L
    val <- nll_fn(p)
    if (is.finite(val) && val < .best_nll) {
      .best_nll  <<- val
      .nll_trace <<- c(.nll_trace, val)
      .par_trace <<- rbind(.par_trace, p)
    }
    if (print_progress && print > 0L && .iter %% print == 0L) {
      row <- .admProgressRow(sprintf("%04d", .iter), val, p, pinfo)
      if (!is.null(row)) message(row)
    }
    val
  }

  eval_grad_f <- if (use_grad) grad_fn else NULL

  lb <- if (use_grad) pmax(ov_lower, p_init - grad_bounds) else ov_lower
  ub <- if (use_grad) pmin(ov_upper, p_init + grad_bounds) else ov_upper

  sc    <- if (!is.null(scale_c)) scale_c else rep(1.0, length(p_init))
  p_sc  <- p_init / sc
  lb_sc <- lb / sc; ub_sc <- ub / sc
  eval_f_sc    <- function(p_s) eval_f(p_s * sc)
  eval_grad_sc <- if (!is.null(eval_grad_f)) function(p_s) eval_grad_f(p_s * sc) * sc else NULL

  # rxLock is a system-wide named mutex -- only used when the model was loaded
  # from cache in this process (lock_rxMod set by caller), never across parallel
  # workers (each has its own independent model instance).
  if (!is.null(lock_rxMod)) {
    tryCatch(rxode2::rxLock(lock_rxMod), error = function(e) NULL)
    on.exit(tryCatch(rxode2::rxUnlock(lock_rxMod), error = function(e) NULL), add = TRUE)
  }

  t0 <- proc.time()
  opt <- tryCatch(
    nloptr::nloptr(
      x0 = p_sc, eval_f = eval_f_sc,
      eval_grad_f = eval_grad_sc,
      lb = lb_sc, ub = ub_sc,
      opts = list(algorithm = algorithm, ftol_rel = ftol_rel, maxeval = maxeval)
    ),
    error = function(e) list(objective = Inf, solution = NULL,
                             message = conditionMessage(e))
  )
  list(restart_id = restart_id,
       objective  = opt$objective,
       solution   = if (!is.null(opt$solution)) opt$solution * sc else p_init,
       n_iter     = .iter,
       nll_trace  = .nll_trace,
       par_trace  = .par_trace,
       elapsed    = as.numeric((proc.time() - t0)["elapsed"]),
       message    = opt$message)
}

.admRunRestarts <- function(worker_fn, p0, ov, pinfo, .ctl, ui, studies,
                            extra_args = list()) {
  n_r <- .ctl$n_restarts
  set.seed(.ctl$seed)
  n_struct <- length(pinfo$struct_names)
  inits <- lapply(seq_len(n_r), function(r) {
    if (r == 1L) return(p0)
    p_new <- p0
    p_new[seq_len(n_struct)] <- p0[seq_len(n_struct)] +
      rnorm(n_struct, sd = .ctl$restart_sd)
    p_new
  })

  ui_lstExpr <- ui$lstExpr
  ov_lower   <- ov$lower
  ov_upper   <- ov$upper
  scale_c    <- ov$scale_c
  seed_val   <- .ctl$seed

  base_args <- list(
    ui_lstExpr = ui_lstExpr, pinfo = pinfo,
    ov_lower = ov_lower, ov_upper = ov_upper, scale_c = scale_c,
    studies = studies, n_sim = .ctl$n_sim, seed = seed_val
  )
  all_args <- c(base_args, extra_args)

  run_one <- function(r) {
    do.call(worker_fn, c(list(restart_id = r, p_init = inits[[r]]), all_args))
  }

  .worker_fn_name <- deparse(substitute(worker_fn))
  pkg_env    <- tryCatch(asNamespace("admixr2"), error = function(e) globalenv())
  pkg_locked <- tryCatch(environmentIsLocked(asNamespace("admixr2")), error = function(e) FALSE)
  if (pkg_locked) {
    .fn_list <- list()
  } else {
    .fn_names <- ls(pkg_env, all.names = TRUE)
    .fn_names <- .fn_names[grepl("^\\.(adm|adfo|adirmc|adgh|softmax|logdmvnorm)", .fn_names)]
    .fn_list  <- setNames(lapply(.fn_names, get, envir = pkg_env), .fn_names)
    .fn_list[[.worker_fn_name]] <- worker_fn
  }

  use_parallel <- n_r > 1L &&
    .ctl$workers > 1L &&
    requireNamespace("mirai", quietly = TRUE) &&
    .adm_worker_env$n > 1L

  .restart_msg <- function(r, res) {
    sec <- if (!is.null(res$elapsed)) res$elapsed else NA_real_
    if (isTRUE(res$final_row_printed)) {
      return(.admProgressTimingRow(sec, pinfo))
    }
    row <- .admProgressRow(sprintf("%04d \u2713", res$n_iter), res$objective, res$solution, pinfo)
    if (!is.null(row))
      return(paste0(row, "\n", .admProgressTimingRow(sec, pinfo)))
    .admProgressTimingRow(sec, pinfo)
  }

  if (use_parallel) {
    effective_workers <- .adm_worker_env$n
    base_tpw          <- max(1L, floor(.ctl$cores / effective_workers))
    remainder         <- max(0L, .ctl$cores - base_tpw * effective_workers)
    cores_vec         <- c(rep(base_tpw + 1L, remainder), rep(base_tpw, effective_workers - remainder))
    tpw_label         <- if (remainder > 0L)
      sprintf("%d-%d", base_tpw, base_tpw + 1L) else as.character(base_tpw)
    n_batches         <- ceiling(n_r / effective_workers)
    batch_label       <- if (n_batches > 1L) sprintf(", %d sequential batch(es)", n_batches) else ""

    message(sprintf("  Running %d restarts in parallel (%d workers, %s threads/worker%s)",
                    n_r, effective_workers, tpw_label, batch_label))
    message(.admProgressHeader(pinfo, bottom = FALSE))

    # Dev-mode stale-install guard. When admixr2 is loaded with devtools::load_all()
    # the parent runs the dev source, but the worker daemons `library(admixr2)` the
    # INSTALLED package (a daemon cannot ADD new bindings to a locked installed
    # namespace, so .fn_list is empty and no dev patch is applied). If the installed
    # package is older than the loaded source -- e.g. it predates a function the
    # parent's gradient path now uses -- the daemons silently compute a DIFFERENT
    # objective than the sequential path, with no error. Warn once so this cannot be
    # mistaken for a real numerical difference.
    #
    # Detect dev mode by the `.__DEVTOOLS__` marker devtools::load_all() stamps
    # into the namespace, NOT environmentIsLocked() -- an INSTALLED namespace is
    # locked too, so pkg_locked is TRUE in production (where daemons and parent run
    # the same code and there is nothing to warn about). Reading the marker
    # directly avoids a dependency on pkgload (which R CMD check would flag as an
    # undeclared `::` import). In production the marker is absent and this never
    # fires.
    .dev_loaded <- isTRUE(tryCatch(
      exists(".__DEVTOOLS__", envir = asNamespace("admixr2"), inherits = FALSE),
      error = function(e) FALSE))
    if (.dev_loaded &&
        !exists("dev_daemon_stale", envir = .adm_warn_env, inherits = FALSE)) {
      assign("dev_daemon_stale", TRUE, envir = .adm_warn_env)
      warning("admixr2: parallel restarts under devtools::load_all() -- worker ",
              "daemons load the INSTALLED admixr2, not your loaded source. If it ",
              "is stale, parallel results silently diverge from sequential. Run ",
              "devtools::install() (once) before comparing parallel vs sequential.",
              call. = FALSE)
    }

    # Daemons are separate processes on every platform: compiled DLLs cannot be
    # serialised, so the worker reloads them from the qs2 cache, and rxEt event
    # tables (~130 MB each) are stripped to plain data frames before sending.
    all_args_par <- all_args
    all_args_par$no_lock <- TRUE
    if ("print_progress"   %in% names(all_args_par)) all_args_par$print_progress   <- FALSE
    if ("rxMod_direct"     %in% names(all_args_par)) all_args_par$rxMod_direct     <- NULL
    if ("sensModel_direct" %in% names(all_args_par)) all_args_par$sensModel_direct <- NULL
    if (!is.null(extra_args$sensModel_direct)) {
      # Take the cache path recorded by .admLoadSensModel() when it wrote the
      # file. Re-deriving it from digest(ui$foceiModel$inner) misses: that access
      # returns a different object than the one digested at save time, so the
      # workers silently fell back to grad = "fd".
      .sm  <- extra_args$sensModel_direct
      .scf <- .sm$cache_file
      if (!is.null(.scf) && file.exists(.scf)) {
        all_args_par$sens_cache_file <- .scf
        all_args_par$sens_cols       <- .sm$sens_cols
        all_args_par$sens_rename     <- .sm$rename_map
      } else {
        message("  [parallel] Sens model cache unavailable; workers will use grad = 'fd'.")
      }
    }
    all_args_par$studies <- lapply(all_args_par$studies, function(s) {
      for (.ev_field in c("ev", "ev_full"))
        if (!is.null(s[[.ev_field]])) s[[.ev_field]] <- as.data.frame(s[[.ev_field]])
      s
    })

    # Batch loop -- print after each batch of effective_workers restarts.
    # .fn_list is empty in installed mode and holds the dev-mode functions under
    # devtools::load_all(); .admDaemonRestart() patches them into the daemon's
    # namespace before dispatching.
    batches <- split(seq_len(n_r), ceiling(seq_len(n_r) / effective_workers))
    results <- vector("list", n_r)

    .map_args <- list(worker_fn_name    = .worker_fn_name,
                      fn_list           = .fn_list,
                      inits             = inits,
                      all_args          = all_args_par,
                      cores_vec         = cores_vec,
                      effective_workers = effective_workers)

    for (.batch in batches) {
      .br <- mirai::mirai_map(.batch, .admDaemonRestart,
                              .args = .map_args, .compute = .adm_compute)[]
      for (i in seq_along(.batch)) {
        res <- .br[[i]]
        if (inherits(res, "errorValue") || inherits(res, "miraiError"))
          stop(sprintf("admixr2: parallel restart %d failed: %s",
                       .batch[[i]], paste(as.character(res), collapse = " ")),
               call. = FALSE)
        results[[.batch[[i]]]] <- res
        message(.restart_msg(.batch[[i]], res))
      }
    }
  } else {
    if (n_r > 1L) {
      if (!requireNamespace("mirai", quietly = TRUE)) {
        message(sprintf("  Running %d restarts sequentially (install mirai for parallel)",
                        n_r))
      } else {
        message(sprintf(paste0("  Running %d restarts sequentially",
                               " (set workers=%d in admControl() for parallel)"),
                        n_r, n_r))
      }
      message(.admProgressHeader(pinfo, bottom = FALSE))
    }
    results <- lapply(seq_len(n_r), function(r) {
      message(.admProgressRestart(r, n_r, pinfo))
      res <- run_one(r)
      message(.restart_msg(r, res))
      res
    })
  }

  nlls <- vapply(results, function(r) r$objective, double(1))
  best <- which.min(nlls)
  best_result <- results[[best]]
  best_result$all_traces <- lapply(results, function(r)
    list(restart_id = r$restart_id,
         nll_trace  = r$nll_trace,
         par_trace  = r$par_trace))
  best_result
}

# -- Main estimation entry point -----------------------------------------------

#' Fit an aggregate data model via Monte Carlo (admc estimator)
#'
#' Called automatically by `nlmixr2(model, admData(), est = "admc", control = admControl(...))`.
#' Not typically called directly.
#'
#' @param env nlmixr2 environment containing `ui` and `control`.
#' @param ... Unused.
#'
#' @return An `admFit` nlmixr2 fit object.
#'
#' @method nlmixr2Est admc
#' @importFrom nlmixr2est nlmixr2Est
#' @export
nlmixr2Est.admc <- function(env, ...) {
  .ui  <- env$ui
  .ctl <- env$control

  if (!inherits(.ctl, "admControl")) .ctl <- getValidNlmixrCtl.admc(.ctl)
  if (!inherits(.ctl, "admControl"))
    stop("Could not recover admControl", call. = FALSE)
  assign("control", .ctl, envir = .ui)

  studies <- .ctl$studies
  if (length(studies) == 0L)
    stop("admControl(studies=...) required", call. = FALSE)
  if (is.null(names(studies)))
    names(studies) <- paste0("study", seq_along(studies))

  pinfo      <- .admParseIniDf(.ui$iniDf, .ui)
  pinfo$nDisplayProgress <- .ctl$nDisplayProgress %||% pinfo$nDisplayProgress
  # Residual-quadrature nodes travel on pinfo -> arr -> .admResidApply/.admResidDeriv.
  pinfo$resid_nodes      <- .ctl$resid_nodes %||% .ADM_TBS_NODES
  output_var <- .admOutputVar(.ui)

  for (nm in names(studies))
    studies[[nm]] <- .admNormaliseStudy(studies[[nm]], nm, output_var)
  # Flatten to observation units (independent single-output blocks, or one joint
  # unit per same-subject study) and attach ev_full. multi_out is model-level
  # (the model has >1 endpoint) so a multi-endpoint model always tags obs by cmt.
  studies    <- .admFlattenStudies(studies)
  multi_out  <- length(.admOutputVars(.ui)) > 1L
  any_joint  <- any(vapply(studies, function(u) isTRUE(u$is_joint), logical(1)))
  studies    <- .admBuildEvFull(studies, tag_cmt = multi_out)

  .admCheckAR(pinfo, studies)
  .admCheckOrdinal(pinfo, studies)
  .admCheckMixedEndpoints(.ui)

  # A beta endpoint's prediction is derived from TWO solved columns; the pair
  # travels on each study so the solve paths can combine them (see .admSimulate).
  .bpair <- .admBetaPair(.ui)
  if (!is.null(.bpair)) {
    studies <- lapply(studies, function(u) { u$out_pair <- .bpair; u })
    # ... and it is fitted DERIVATIVE-FREE. beta's conditional variance is
    # mu(1-mu)/(1+phi) with phi = b1 + b2 SOLVED from the structural model, so a
    # structural theta reaches the objective through phi as well as through mu.
    # Every gradient path here chains through mu only (dpred), which makes the
    # analytic/FD-of-the-prediction gradient a gradient of the wrong function --
    # it holds phi fixed at the value it had before the perturbation. BOBYQA
    # differences the objective itself, where phi moves with the thetas as it
    # should. .adfoNLL/.adirmcNLL refuse beta outright for the related reason
    # that they have no phi at all.
    if (.ctl$grad != "none") {
      message("admControl: a beta() endpoint is fitted derivative-free ",
              "(grad = \"none\"): its precision is solved from the structural ",
              "model, and the gradient paths carry only d(prediction)/d(theta).")
      .ctl$grad      <- "none"
      .ctl$algorithm <- .admDefaultAlgorithm("none")
    }
  }

  want_grad    <- .ctl$grad != "none"
  want_sens    <- .ctl$grad == "sens"
  want_central <- .ctl$grad == "cfd"
  # Multi-compartment gradient. Independent blocks use the per-output analytical/
  # sens gradient. Joint (same-subject) fits are scored by a stacked MVN: with a
  # sens model + etas the joint gradient is analytical (.admGrad joint branch);
  # otherwise it falls back to the common-random-number FD of the aggregate NLL
  # (joint_fd, computed below once the sens model has been resolved).

  .unpaired <- if (!is.null(pinfo$struct_has_eta))
    names(pinfo$struct_has_eta)[!pinfo$struct_has_eta] else character(0)


  # ORDERING INVARIANT: .admLoadSensModel() must run before .admLoadModel().
  # .admLoadModel() calls rxode2::rxode2(ui) which triggers nlmixr2est's foceiModel
  # compilation via its FD path for linCmt, caching inner=NULL. Calling this first
  # ensures the foceiModel (inner != NULL) is compiled and cached before that happens.
  sensModel <- if (want_sens) {
    sm <- tryCatch(.admLoadSensModel(.ui), error = function(e) NULL)
    if (is.null(sm))
      warning("admControl(grad='sens'): sensitivity model unavailable -- falling back to forward FD")
    else if (isTRUE(sm$is_lincmt))
      warning("admControl(grad='sens'): linCmt sensitivity model detected; grad='fd' is typically faster for linCmt models -- consider switching to admControl(grad='fd')")
    sm
  } else NULL

  # Unpaired (non-mu-referenced) struct thetas. The sens model carries an explicit
  # THETA_j_ direction per unpaired theta (.admBuildThetaSens), so it returns d(pred)/d(theta)
  # for them and the whole gradient stays analytical -- including a model with NO
  # mu-referenced theta at all, which used to force a full FD gradient. If the
  # augmentation was unavailable (theta_sens_cols NULL) the old behaviour stands:
  # sens for the paired thetas + FD for the unpaired, or full FD when none is paired.
  if (!any_joint && pinfo$n_eta > 0L && length(.unpaired) > 0L) {
    .theta_sens <- want_sens && !is.null(sensModel) &&
      !is.null(sensModel$theta_sens_cols) &&
      all(.unpaired %in% names(sensModel$theta_sens_cols))
    if (.theta_sens) {
      message(sprintf(
        "admc: struct theta(s) without mu-referencing: %s. Sens model carries their sensitivities (no FD).",
        paste(.unpaired, collapse = ", ")))
    } else if (want_sens && all(!pinfo$struct_has_eta)) {
      message(sprintf("admc: no mu-referenced struct thetas (%s); falling back to full forward FD.",
                      paste(.unpaired, collapse = ", ")))
      want_sens <- FALSE
      sensModel <- NULL
    } else {
      message(sprintf(
        "admc: struct theta(s) without mu-referencing: %s. %s",
        paste(.unpaired, collapse = ", "),
        if (want_sens) "Sens model for paired thetas; FD for unpaired."
        else "FD gradient for these parameters."
      ))
    }
  }

  rxMod <- .admLoadModel(.ui)
  rxode2::rxLock(rxMod)
  # Reclaim compiled models with rxode2's own idiom -- the gc(); rxUnloadAll()
  # nlmixr2est runs per fit -- so a session of many fits does not accumulate models
  # (and RSS) unbounded. rxUnloadAll() keeps the last getOption("rxode2.dontUnload",
  # 10) models, and this fit registers only ~6 (sim + sens + the four foceiModel
  # companions), so the model driving the returned fit stays loaded for nlmixr2's
  # post-fit output/table solve; only OLDER models (from earlier fits) are freed --
  # exactly the semantics nlmixr2est's own rxUnloadAll() at fit start has.
  on.exit({ rxode2::rxUnlock(rxMod); rxode2::rxSolveFree(); gc(FALSE); rxode2::rxUnloadAll() },
          add = TRUE)

  set.seed(.ctl$seed)
  z_list      <- .admMakeZ(.ctl$n_sim, pinfo, length(studies), .ctl$sampling)
  params_list <- .admMakeParamsList(.ctl$n_sim, pinfo, length(studies))

  ov     <- .admBuildOptVec(pinfo)
  .iter  <- 0L
  cores  <- .ctl$cores
  grad_h <- .ctl$grad_h

  .nll_trace <- numeric(0)
  .par_trace <- NULL
  .best_nll  <- Inf

  eval_f <- function(p) {
    .iter <<- .iter + 1L
    val <- .admNLL(p, pinfo, studies, z_list, rxMod, output_var, params_list, cores)
    if (is.finite(val) && val < .best_nll) {
      .best_nll  <<- val
      .nll_trace <<- c(.nll_trace, val)
      .par_trace <<- rbind(.par_trace, p)
    }
    if (.ctl$print > 0L && .iter %% .ctl$print == 0L) {
      row <- .admProgressRow(sprintf("%04d", .iter), val, p, pinfo)
      if (!is.null(row)) message(row)
    }
    val
  }

  # Joint fits use FD only when no sens model is available; otherwise .admGrad's
  # joint branch computes the analytical stacked-MVN gradient.
  joint_fd <- any_joint && is.null(sensModel)
  eval_grad_f <- if (!want_grad) NULL
    else if (joint_fd)
      function(p) .admNLLGradFD(p, pinfo, studies, z_list, rxMod, output_var,
                                params_list, cores, grad_h,
                                use_central = want_central)
    else function(p) .admGrad(p, pinfo, studies, z_list, rxMod, output_var,
                              params_list, cores, grad_h, sensModel,
                              use_central = want_central)

  grad_label <- if (!want_grad) "none"
  else if (joint_fd) (if (want_central) "central FD (joint)" else "forward FD (joint)")
  else if (any_joint) "Sens (joint)"
  else if (!is.null(sensModel))
    if (pinfo$has_kappa) "Sens+FD" else "Sens"
  else if (want_central) "central FD" else "forward FD"
  message("=== admixr2: Aggregate Data Modeling (MC) ===")
  message(sprintf("  Obs units: %d | MC samples: %d | Params: %d | Cores: %d | Grad: %s | Restarts: %d",
                  length(studies), .ctl$n_sim, length(ov$p0), cores,
                  grad_label, .ctl$n_restarts))
  t0 <- proc.time()

  lb <- if (want_grad) pmax(ov$lower, ov$p0 - .ctl$grad_bounds) else ov$lower
  ub <- if (want_grad) pmin(ov$upper, ov$p0 + .ctl$grad_bounds) else ov$upper

  sc           <- ov$scale_c
  p0_sc        <- ov$p0 / sc
  lb_sc        <- lb    / sc
  ub_sc        <- ub    / sc
  eval_f_sc    <- function(p_s) eval_f(p_s * sc)
  eval_grad_sc <- if (!is.null(eval_grad_f)) {
    function(p_s) eval_grad_f(p_s * sc) * sc
  } else NULL

  if (.ctl$n_restarts == 1L) {
    message(.admProgressHeader(pinfo))
    opt_raw <- nlmixr2est::nlmixrWithTiming("admc", {
      nloptr::nloptr(x0 = p0_sc, eval_f = eval_f_sc,
                     eval_grad_f = eval_grad_sc,
                     lb = lb_sc, ub = ub_sc,
                     opts = list(algorithm = .ctl$algorithm,
                                 ftol_rel  = .ctl$ftol_rel,
                                 maxeval   = .ctl$maxeval))
    })
    opt <- list(objective = opt_raw$objective,
                solution  = opt_raw$solution * sc,
                message   = opt_raw$message)
    if (.ctl$print > 0L) {
      row <- .admProgressRow(sprintf("%04d \u2713", .iter), opt$objective, opt$solution, pinfo)
      if (!is.null(row)) message(paste0(row, "\n",
        .admProgressTimingRow((proc.time() - t0)["elapsed"], pinfo)))
    }
    opt$all_traces <- list(list(restart_id = 1L,
                                nll_trace  = .nll_trace,
                                par_trace  = .par_trace))
  } else {
    .admSetupDaemons(.ctl, .ctl$n_restarts)
    on.exit(.admStopDaemons(), add = TRUE)
    opt <- .admRunRestarts(
      worker_fn  = .admRestartWorker,
      p0         = ov$p0, ov = ov, pinfo = pinfo,
      .ctl       = .ctl, ui = .ui, studies = studies,
      extra_args = list(algorithm    = .ctl$algorithm,
                        ftol_rel     = .ctl$ftol_rel,
                        maxeval      = .ctl$maxeval,
                        use_grad     = want_grad,
                        grad_h       = .ctl$grad_h,
                        grad_bounds  = .ctl$grad_bounds,
                        output_var   = output_var,
                        sampling     = .ctl$sampling,
                        use_central  = want_central,
                        print_progress   = TRUE,
                        print            = .ctl$print,
                        cores            = .ctl$cores,
                        rxMod_direct     = rxMod,
                        sensModel_direct = sensModel)
    )
    .admStopDaemons()
    .iter <- opt$n_iter
  }

  t_opt     <- (proc.time() - t0)["elapsed"]
  final     <- .admUnpack(opt$solution, pinfo)
  fullTheta <- .admFullTheta(final, pinfo)

  p_hat  <- setNames(opt$solution, names(ov$p0))
  t0_cov <- proc.time()
  .cov <- if (.ctl$covMethod == "r") {
    # Multi-output / joint fits use the NLL-FD Hessian (via .admNLLBatch); the
    # grad-FD Hessian relies on the single-output analytical grad batch.
    use_grad_cov <- want_grad && !multi_out && !any_joint
    use_cent_cov <- want_central
    # struct + sigma + OMEGA: the Hessian spans all three, so the evaluation
    # count must too.
    np_cov <- length(pinfo$struct_names) + length(pinfo$sigma_names) +
              length(pinfo$omega_par)
    n_evals <- if (use_grad_cov) {
      np_cov + 1L
    } else {
      n_off <- np_cov * (np_cov - 1L) / 2L
      np_cov * 2L + n_off * 4L + 1L
    }
    evals_label <- if (use_grad_cov) "gradient evaluations" else "NLL evaluations"
    hess_label  <- if (!use_grad_cov) "" else if (!is.null(sensModel))
      ", Sens-Hessian" else if (use_cent_cov) ", cFD-Hessian" else ", FD-Hessian"
    message(sprintf("  Computing covariance (R method%s, %d %s)",
                    hess_label, n_evals, evals_label))
    tryCatch(
      .admCalcCov(p_hat, pinfo, studies, z_list, rxMod, output_var,
                  params_list, cores, cov_n_sim = .ctl$cov_n_sim,
                  use_grad = use_grad_cov, grad_h = .ctl$grad_h,
                  cov_h = .ctl$cov_h, cov_h_outer = .ctl$cov_h_outer,
                  sensModel = sensModel, use_central = use_cent_cov,
                  sampling = .ctl$sampling),
      error = function(e) { warning("admCalcCov failed: ", conditionMessage(e)); NULL })
  } else NULL
  # A NULL covariance used to be completely silent: no warning reached the user,
  # `warnings()` was empty, covMethod came back "" and every SE was NA with no
  # indication why. Say so once, from the driver, where it cannot be swallowed.
  if (isTRUE(.ctl$covMethod == "r") && is.null(.cov))
    warning("covariance could not be computed (the Hessian was singular or ",
            "non-finite); standard errors are unavailable for this fit.",
            call. = FALSE)
  # iniDf order first (nlmixr2est maps SEs positionally), then snapshot the names
  # BEFORE nlmixr2est sees it -- .admCovThetaOrder()/.admRestoreCovNames().
  .cov      <- .admCovThetaOrder(.cov, .ui)
  .cov_nms  <- .admCovNames(.cov)
  t_cov     <- (proc.time() - t0_cov)["elapsed"]
  t_elapsed <- t_opt + t_cov

  if (.ctl$returnAdmr) {
    return(list(objective = opt$objective, fullTheta = fullTheta,
                struct = final$struct, sigma_var = final$sigma_var,
                omega = final$omega, L = final$L, nloptr = opt,
                cov = .cov))
  }

  .ret            <- new.env(parent = emptyenv())
  .ret$table      <- env$table
  .ret$ui         <- .ui
  .ret$fullTheta  <- fullTheta
  .ret$objective  <- opt$objective
  .ret$est        <- "admc"
  .ret$ofvType    <- "admc"
  .ret$adjObf     <- FALSE
  .ret$covMethod  <- if (!is.null(.cov)) "r" else ""
  .ret$cov        <- .cov
  .ret$message    <- opt$message
  .ret$extra      <- ""
  .ret$origData   <- studies

  .ret$admExtra <- list(struct        = final$struct,
                        sigma_var     = final$sigma_var,
                        sigma_is_prop  = pinfo$sigma_is_prop,
                        sigma_is_lnorm = pinfo$sigma_is_lnorm,
                        omega         = final$omega,
                        L             = final$L,
                        eta_col_names = pinfo$eta_col_names,
                        par_names     = names(ov$p0),
                        npar          = length(ov$p0),
                        nloptr        = opt,
                        nll_trace     = .nll_trace,
                        par_trace     = .par_trace,
                        all_traces    = opt$all_traces,
                        n_iter        = .iter,
                        time          = t_elapsed,
                        t_opt         = t_opt,
                        t_cov         = t_cov,
                        studies       = studies,
                        n_sim         = .ctl$n_sim,
                        sampling      = .ctl$sampling)

  nlmixr2est::.nlmixr2FitUpdateParams(.ret)
  nmObjHandleControlObject.admControl(.ctl, .ret)
  if (exists("control", .ui)) rm(list = "control", envir = .ui)
  .ret$control <- .admToFoceiControl(.ctl, .admCovSkip(.cov, .ui))
  .focei_model <- suppressMessages(tryCatch(.ui$foceiModel, error = function(e) NULL))
  if (!is.null(.focei_model)) .ret$model <- .focei_model

  .fit <- nlmixr2est::nlmixr2CreateOutputFromUi(
    .ui, data = if (multi_out) admData(.admEndpointNames(.ui)) else admData(),
    control = .ret$control,
    table = .ret$table, env = .ret, est = "admc")

  .fit$env$method   <- "admc"
  .admRestoreCovNames(.fit, .cov_nms)
  .fit$env$studies  <- studies
  .fit$env$admExtra <- .ret$admExtra
  # Populate nlmixr2-style parameter history so traceplot(fit) works natively.
  .admAttachParHist(.fit, .ret$admExtra$all_traces, .ret$admExtra$par_names, .ui)
  # Store observed + predicted aggregate moments (E vector, V matrix) per study.
  .admAttachAggData(.fit, .ret$admExtra, .ui)
  .old_cls <- class(.fit)
  .new_cls <- c("admFit", .old_cls)
  attr(.new_cls, ".foceiEnv") <- attr(.old_cls, ".foceiEnv")
  class(.fit) <- .new_cls

  .stats <- .admCalcObjStats(opt$objective, length(ov$p0), studies)
  row.names(.stats$objDf) <- "admc"
  .fit$env$logLik    <- .stats$ll
  .fit$env$nobs      <- .stats$nobs
  .fit$env$objDf     <- .stats$objDf
  .fit$env$OBJF      <- .stats$objDf$OBJF
  .fit$env$AIC       <- .stats$objDf$AIC
  .fit$env$BIC       <- .stats$objDf$BIC
  .fit$env$objective <- opt$objective
  .fit$env$time     <- data.frame(
    optimize   = t_opt,
    covariance = t_cov,
    other      = 0,
    elapsed    = t_elapsed,
    row.names  = NULL
  )

  .fit
}
