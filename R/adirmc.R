# -- Control object -------------------------------------------------------------

#' Control settings for the IRMC estimator
#'
#' Constructs a control object for `est = "adirmc"`, the Iterative Reweighting
#' Monte Carlo estimator.
#'
#' Multi-compartment fits (a study `observations` list with several observed
#' outputs) are **not** supported by `adirmc`; use `est = "admc"`, `"adfo"`, or
#' `"adgh"` for those. Single-output studies are fit as usual.
#'
#' @inheritParams admControl
#' @param grad Gradient mode for the inner optimiser: `"analytical"` (default,
#'   closed-form weight-path gradient), `"none"` (derivative-free BOBYQA), or
#'   `"fd"` (finite differences). Note: `"sens"` and `"cfd"` are not available
#'   for the IRMC estimator.
#' @param kappa_method Kappa correction method for models with non-mu-referenced
#'   struct thetas: `"exact"` (default, re-evaluates population prediction `f(theta, 0)`
#'   via rxSolve at each inner step), `"linearized"` (precomputes `J = df/d(theta)`
#'   once per outer iteration using `f(theta, 0)` as baseline — zero rxSolve per inner step),
#'   or `"linearized_gh"` (same linear approximation but baseline and Jacobian use
#'   Gauss-Hermite quadrature `E_GH[f(theta, eta)]` instead of `f(theta, 0)` — more
#'   accurate baseline at any IIV magnitude, still zero rxSolve per inner step).
#' @param kappa_n_nodes Number of GH nodes per eta dimension for
#'   `kappa_method = "linearized_gh"` (default 5). Total quadrature points =
#'   `kappa_n_nodes^n_eta`. Ignored for other kappa methods.
#' @param outer_iter Maximum inner optimiser iterations per phase.
#' @param omega_expansion Inflate proposal Omega by this factor (>= 1).
#' @param phases Numeric vector of box-constraint half-widths, one per phase.
#'   Phases progressively tighten the search region.
#' @param convcrit Convergence criterion: phase ends when `|approx - exact| < convcrit`.
#' @param max_worse Stop a phase after this many consecutive worsening iterations.
#'
#' @return An object of class `adirmcControl`.
#'
#' @examples
#' # Inspect defaults
#' ctl <- adirmcControl()
#' ctl$phases
#' ctl$omega_expansion
#'
#' # Tighter phases, more restarts
#' ctl2 <- adirmcControl(
#'   n_sim           = 1000L,
#'   omega_expansion = 1.5,
#'   phases          = c(2, 1, 0.5, 0.01),
#'   n_restarts      = 3L
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
#' V <- diag(diag(cov.wt(dv_mat, method = "ML")$cov))
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
#'   pk_model, admData(), est = "adirmc",
#'   control = adirmcControl(
#'     studies = list(study1 = list(E = E, V = V, n = length(ids),
#'                                  times = times, ev = et(amt = 100))),
#'     n_sim   = 500L
#'   )
#' )
#' print(fit)
#' }
#'
#' @export
adirmcControl <- function(
    studies         = list(),
    n_sim           = 2500L,
    outer_iter      = 50L,
    sampling        = c("sobol", "halton", "torus", "lhs", "rnorm"),
    algorithm       = NULL,
    maxeval         = 5000L,
    ftol_rel        = .Machine$double.eps,
    print           = 1L,
    omega_expansion = 1.0,
    seed            = 12345L,
    cores           = 1L,
    nDisplayProgress = .Machine$integer.max,
    grad            = c("analytical", "none", "fd"),
    kappa_method    = c("exact", "linearized", "linearized_gh"),
    kappa_n_nodes   = 5L,
    grad_h          = 1e-4,
    cov_h           = 1e-3,
    cov_h_outer     = .Machine$double.eps^(1/5),
    phases          = c(2, 1, 0.5, 0.01),
    convcrit        = 1e-5,
    max_worse       = 5L,
    covMethod       = c("r", "none"),
    cov_n_sim       = 10000L,
    n_restarts      = 1L,
    restart_sd      = 0.2,
    workers         = 1L,
    rxControl       = NULL,
    calcTables      = FALSE,
    compress        = TRUE,
    ci              = 0.95,
    sigdig          = 4,
    sigdigTable     = NULL,
    addProp         = c("combined2", "combined1"),
    optExpression   = TRUE,
    sumProd         = FALSE,
    literalFix      = TRUE,
    returnAdmr      = FALSE,
    ...) {

  .xtra <- list(...)
  if (length(.xtra) > 0)
    stop("adirmcControl: unused argument(s): ",
         paste(paste0("'", names(.xtra), "'"), collapse = ", "), call. = FALSE)

  addProp      <- match.arg(addProp)
  grad         <- match.arg(grad)
  kappa_method <- match.arg(kappa_method)
  sampling     <- match.arg(sampling)

  checkmate::assertList(studies)
  checkmate::assertIntegerish(n_sim,        lower = 1L,  len = 1)
  checkmate::assertIntegerish(outer_iter,   lower = 1L,  len = 1)
  checkmate::assertIntegerish(maxeval,      lower = 1L,  len = 1)
  checkmate::assertNumeric(ftol_rel,        lower = 0,   len = 1)
  checkmate::assertIntegerish(print,        lower = 0L,  len = 1)
  checkmate::assertNumeric(omega_expansion, lower = 1,   len = 1)
  checkmate::assertIntegerish(seed,                      len = 1)
  checkmate::assertIntegerish(cores,        lower = 1L,  len = 1)
  checkmate::assertIntegerish(nDisplayProgress, lower = 1L, len = 1,
                              .var.name = "nDisplayProgress")
  checkmate::assertNumeric(grad_h,          lower = 0,   len = 1)
  checkmate::assertNumeric(cov_h,       lower = 0, len = 1, .var.name = "cov_h")
  checkmate::assertNumeric(cov_h_outer, lower = 0, len = 1, .var.name = "cov_h_outer")
  checkmate::assertNumeric(phases,          lower = 0,   min.len = 1)
  checkmate::assertNumeric(convcrit,        lower = 0,   len = 1)
  checkmate::assertIntegerish(max_worse,    lower = 1L,  len = 1)
  checkmate::assertIntegerish(kappa_n_nodes, lower = 1L, len = 1)
  covMethod <- match.arg(covMethod)
  checkmate::assertIntegerish(cov_n_sim,    lower = 1L,  len = 1)
  checkmate::assertIntegerish(n_restarts,   lower = 1L,  len = 1)
  checkmate::assertNumeric(restart_sd,      lower = 0,   len = 1)
  checkmate::assertIntegerish(workers,      lower = 1L,  len = 1)

  .algo     <- .admResolveAlgorithm(algorithm, grad,
                                    .var.name = "adirmcControl: algorithm")
  algorithm <- .algo$algorithm
  grad      <- .algo$grad

  if (is.null(rxControl))   rxControl   <- rxode2::rxControl(sigdig = sigdig)
  if (is.null(sigdigTable)) sigdigTable <- max(round(sigdig), 3L)

  .ret <- list(
    studies         = studies,
    n_sim           = as.integer(n_sim),
    outer_iter      = as.integer(outer_iter),
    sampling        = sampling,
    algorithm       = algorithm,
    maxeval         = as.integer(maxeval),
    ftol_rel        = ftol_rel,
    print           = as.integer(print),
    omega_expansion = omega_expansion,
    seed            = as.integer(seed),
    cores           = as.integer(cores),
    nDisplayProgress = as.integer(nDisplayProgress),
    grad            = grad,
    kappa_method    = kappa_method,
    kappa_n_nodes   = as.integer(kappa_n_nodes),
    grad_h          = grad_h,
    cov_h           = cov_h,
    cov_h_outer     = cov_h_outer,
    phases          = phases,
    convcrit        = convcrit,
    max_worse       = as.integer(max_worse),
    covMethod       = covMethod,
    cov_n_sim       = as.integer(cov_n_sim),
    n_restarts      = as.integer(n_restarts),
    restart_sd      = restart_sd,
    workers         = as.integer(workers),
    rxControl       = rxControl,
    calcTables      = calcTables,
    compress        = compress,
    ci              = ci,
    sigdig          = sigdig,
    sigdigTable     = as.integer(sigdigTable),
    addProp         = addProp,
    optExpression   = optExpression,
    sumProd         = sumProd,
    literalFix      = literalFix,
    returnAdmr      = returnAdmr
  )
  class(.ret) <- "adirmcControl"
  .ret
}

# -- nlmixr2 S3 hooks -----------------------------------------------------------

#' @noRd
getValidNlmixrCtl.adirmc <- function(control) {
  if (inherits(control, "adirmcControl")) return(control)
  .ctl <- control[[1]]
  if (inherits(.ctl, "adirmcControl")) return(.ctl)
  if (inherits(.ctl, "admControl") || inherits(control, "admControl")) {
    src    <- if (inherits(control, "admControl")) control else .ctl
    shared <- intersect(names(src), names(formals(adirmcControl)))
    return(do.call(adirmcControl, src[shared]))
  }
  if (is.list(.ctl) && "studies" %in% names(.ctl))
    return(do.call(adirmcControl,
                   .ctl[intersect(names(.ctl), names(formals(adirmcControl)))]))
  stop("est='adirmc' requires adirmcControl(...), not admControl(...)", call. = FALSE)
}

#' @noRd
nmObjHandleControlObject.adirmcControl <- function(control, env) {
  assign("adirmcControl", control, envir = env)
}

#' @noRd
nmObjGetControl.adirmc <- function(x, ...) {
  .env <- x[[1]]
  for (.nm in c("adirmcControl", "control")) {
    if (exists(.nm, .env)) {
      .ctl <- get(.nm, .env)
      if (inherits(.ctl, "adirmcControl")) return(.ctl)
    }
  }
  stop("cannot find irmc control object", call. = FALSE)
}

# -- MVN utilities -------------------------------------------------------------

.softmax <- function(lw) softmax_cpp(lw)

# -- Analytical gradient of IRMC inner NLL -------------------------------------
# Proposals are fixed -> inner NLL is deterministic -> gradient is exact.
# has_kappa=FALSE (all params mu-referenced, weight-shift):
#   struct(paired) -> mean_new -> w -> mu_w, V_w -> -2LL  (analytical via d_logback_dp)
#   omega          -> log_new -> w -> mu_w, V_w -> -2LL   (analytical via S kernel)
#   sigma          -> V_w -> -2LL                         (direct, analytical)
# has_kappa=TRUE (any unpaired param, combined weight-shift + kappa):
#   struct(paired) -> mean_new -> w -> mu_w  (weight-shift, analytical via d_logback_dp)
#   struct(single) -> kappa_jac (linearized) or kappa_fn_batch CFD (exact) -> mu_w
#   omega          -> log_new -> w -> mu_w, V_w -> -2LL   (analytical via S kernel)
#   sigma          -> V_w -> -2LL                         (direct, analytical)

# Fused NLL + analytical gradient in one study pass.
# Called by the memo-cached eval_f / eval_grad_f pair in .adirmcPhaseLoop for
# grad_mode = "analytical".  Returns list(nll = scalar, grad = numeric(np)).
# Eliminates the redundant recomputation of IS weights, weighted mean/cov,
# sigma corrections, and kappa that would otherwise happen when the optimiser
# calls eval_f(p) and eval_grad_f(p) sequentially at the same p.
.adirmcNLLAndGrad <- function(p, pinfo, studies_snap, proposals) {
  pars <- tryCatch(.admUnpack(p, pinfo), error = function(e) NULL)
  if (is.null(pars)) return(list(nll = Inf, grad = rep(NA_real_, length(p))))
  if (pinfo$n_eta > 0 && any(diag(pars$omega) <= 0))
    return(list(nll = Inf, grad = rep(NA_real_, length(p))))

  h     <- 1e-5
  n_s   <- length(pinfo$struct_names)
  n_e   <- length(pinfo$sigma_names)
  n_eta <- pinfo$n_eta
  np    <- length(p)
  grad  <- numeric(np); names(grad) <- names(p)
  nll2  <- 0.0

  if (n_eta > 0L) {
    cholO <- tryCatch(chol(pars$omega), error = function(e) NULL)
    if (is.null(cholO)) return(list(nll = Inf, grad = rep(NA_real_, np)))
    invO  <- chol2inv(cholO)
  } else {
    invO <- matrix(0, 0, 0)
  }

  for (si in seq_along(studies_snap)) {
    prop  <- proposals[[si]]
    s     <- studies_snap[[si]]
    F     <- prop$rawpreds
    bi    <- prop$bi

    paired_nms <- names(prop$log_origbeta)
    mean_new_paired <- if (length(paired_nms) > 0L)
      compute_mean_new_cpp(pars$struct[paired_nms], prop$log_origbeta,
                           prop$paired_types, prop$paired_lows, prop$paired_his)
    else numeric(0)

    n_eta_bi <- ncol(bi)
    mean_new <- rep(0.0, n_eta_bi)
    if (length(mean_new_paired) > 0L)
      mean_new[prop$paired_eta_pos] <- mean_new_paired

    L_cur   <- if (!is.null(pars$L)) pars$L else matrix(0, 0, 0)
    log_new <- logdmvnorm_batch_cpp(bi, mean_new, L_cur)
    lw      <- log_new - prop$log_prop
    w       <- .softmax(lw)

    wmc <- weighted_meancov_cpp(F, w)
    mu  <- as.numeric(wmc$mu)
    V   <- wmc$V

    mu_struct <- mu
    for (k in seq_along(pars$sigma_var))
      if (prop$sigma_type[k] == 2L)
        mu <- mu * exp(pars$sigma_var[k] / 2)
    for (k in seq_along(pars$sigma_var)) {
      sv <- pars$sigma_var[k]
      if (prop$sigma_type[k] == 1L)
        diag(V) <- diag(V) + sv * mu_struct^2
      else if (prop$sigma_type[k] == 2L)
        diag(V) <- diag(V) + mu^2 * (exp(sv) - 1)
      else
        diag(V) <- diag(V) + sv
    }

    mu_sigma <- mu
    kappa_delta <- if (!is.null(prop$mu_pop))
      prop$kappa_fn(pars$struct) - prop$mu_pop
    else numeric(0)
    if (length(kappa_delta) > 0L) mu <- mu + kappa_delta

    obs_E <- s$E
    r     <- obs_E - mu
    if (identical(s$method, "var")) {
      v_pred       <- diag(V)
      nll2         <- nll2 + s$n * sum(log(v_pred) + s$v_diag / v_pred + r^2 / v_pred)
      if (!is.finite(nll2)) return(list(nll = Inf, grad = rep(NA_real_, np)))
      dNLL_dmu     <- s$n * as.numeric(-2 * r / v_pred)
      dNLL_dV_diag <- s$n * (1 / v_pred - s$v_diag / v_pred^2 - r^2 / v_pred^2)
    } else {
      cholV <- tryCatch(chol(V), error = function(e) NULL)
      if (is.null(cholV)) return(list(nll = Inf, grad = rep(NA_real_, np)))
      invV         <- chol2inv(cholV)
      logdet       <- 2 * sum(log(diag(cholV)))
      nll2         <- nll2 + s$n * (logdet + sum(s$V * invV) + as.numeric(t(r) %*% invV %*% r))
      if (!is.finite(nll2)) return(list(nll = Inf, grad = rep(NA_real_, np)))
      dNLL_dmu     <- s$n * as.numeric(-2 * invV %*% r)
      dNLL_dV      <- s$n * (invV - invV %*% (s$V + tcrossprod(r)) %*% invV)
      dNLL_dV_diag <- diag(dNLL_dV)
    }

    # eff_dNLL_dmu folds in sigma V sensitivity w.r.t. mu_struct.
    # prop:  dV_diag/dmu_struct = 2*sv*mu_struct -> +2*sv*mu_sigma*dNLL_dV_diag
    # lnorm: also scales residual path by exp(sv/2): +(exp(sv/2)-1)*dNLL_dmu
    #        plus V path: +2*exp(sv/2)*mu_sigma*(exp(sv)-1)*dNLL_dV_diag
    eff_dNLL_dmu <- dNLL_dmu
    for (k in seq_along(pars$sigma_var)) {
      sv <- pars$sigma_var[k]
      if (prop$sigma_type[k] == 1L) {
        eff_dNLL_dmu <- eff_dNLL_dmu + 2 * sv * mu_sigma * dNLL_dV_diag
      } else if (prop$sigma_type[k] == 2L) {
        eff_dNLL_dmu <- eff_dNLL_dmu +
          (exp(sv / 2) - 1) * dNLL_dmu +
          2 * exp(sv / 2) * mu_sigma * (exp(sv) - 1) * dNLL_dV_diag
      }
    }

    d_mat <- sweep(bi, 2L, mean_new)

    gk <- if (identical(s$method, "var"))
      irmc_grad_kernel_var_cpp(F, w, mu, d_mat, invO, eff_dNLL_dmu, dNLL_dV_diag)
    else
      irmc_grad_kernel_cpp(F, w, mu, d_mat, invO, eff_dNLL_dmu, dNLL_dV)
    dNLL_dmean_new <- as.numeric(gk$dNLL_dmean_new)

    for (k in seq_along(paired_nms)) {
      nm  <- paired_nms[k]
      .tr <- pinfo$struct_transforms[[nm]]
      d_logback_dp <- if (is.null(.tr) || .tr$curEval %in% c("exp", "log")) {
        1.0
      } else if (.tr$curEval %in% c("expit", "logit")) {
        p_val <- pars$struct[[nm]]
        s_val <- 1 / (1 + exp(-p_val))
        back  <- .tr$low + (.tr$hi - .tr$low) * s_val
        (.tr$hi - .tr$low) * s_val * (1 - s_val) / back
      } else if (.tr$curEval %in% c("probitInv", "probit")) {
        p_val <- pars$struct[[nm]]
        back  <- .tr$low + (.tr$hi - .tr$low) * pnorm(p_val)
        (.tr$hi - .tr$low) * dnorm(p_val) / back
      } else {
        (.admLogBackTransform(pars$struct[[nm]] + 1e-7, .tr) -
           .admLogBackTransform(pars$struct[[nm]] - 1e-7, .tr)) / (2e-7)
      }
      k_struct <- match(nm, pinfo$struct_names)
      eta_pos  <- if (k <= length(prop$paired_eta_pos)) prop$paired_eta_pos[k] else k
      grad[k_struct] <- grad[k_struct] + dNLL_dmean_new[eta_pos] * d_logback_dp
    }

    if (!is.null(prop$kappa_jac)) {
      for (i_kb in seq_along(prop$kappa_grad_idxs)) {
        ki <- prop$kappa_grad_idxs[i_kb]
        grad[ki] <- grad[ki] + sum(eff_dNLL_dmu * prop$kappa_jac[i_kb, ])
      }
    } else if (!is.null(prop$mu_pop) && !is.null(prop$kappa_fn_batch)) {
      n_kb        <- length(prop$kappa_grad_idxs)
      struct_list <- vector("list", 2L * n_kb)
      for (i_kb in seq_along(prop$kappa_grad_idxs)) {
        ki    <- prop$kappa_grad_idxs[i_kb]
        nm_kb <- pinfo$struct_names[ki]
        sh <- pars$struct; sh[nm_kb] <- sh[nm_kb] + h
        sl <- pars$struct; sl[nm_kb] <- sl[nm_kb] - h
        struct_list[[2L * i_kb - 1L]] <- sh
        struct_list[[2L * i_kb]]      <- sl
      }
      kappa_batch <- prop$kappa_fn_batch(struct_list)
      for (i_kb in seq_along(prop$kappa_grad_idxs)) {
        ki     <- prop$kappa_grad_idxs[i_kb]
        dkappa <- (kappa_batch[2L * i_kb - 1L, ] - kappa_batch[2L * i_kb, ]) / (2 * h)
        grad[ki] <- grad[ki] + sum(eff_dNLL_dmu * dkappa)
      }
    }

    if (n_eta > 0L) {
      GL   <- invO %*% gk$S %*% (invO %*% pars$L)
      k_om <- n_s + n_e
      for (r_idx in seq_len(nrow(pinfo$eta_rows_df))) {
        k_om <- k_om + 1L
        ei   <- pinfo$eta_rows_df$neta1[r_idx]
        ej   <- pinfo$eta_rows_df$neta2[r_idx]
        ri   <- max(ei, ej); ci <- min(ei, ej)
        grad[k_om] <- grad[k_om] + if (ei == ej)
          GL[ri, ci] * pars$L[ri, ci] / 2
        else
          GL[ri, ci]
      }
    }

    k_sig <- n_s
    for (k in seq_along(pars$sigma_var)) {
      k_sig <- k_sig + 1L
      sv <- pars$sigma_var[k]
      if (prop$sigma_type[k] == 1L) {
        grad[k_sig] <- grad[k_sig] + sum(dNLL_dV_diag * sv * mu^2)
      } else if (prop$sigma_type[k] == 2L) {
        grad[k_sig] <- grad[k_sig] + sv * (
          sum(dNLL_dV_diag * mu^2 * (2 * exp(sv) - 1)) +
          sum(dNLL_dmu * mu) / 2
        )
      } else {
        grad[k_sig] <- grad[k_sig] + sum(dNLL_dV_diag) * sv
      }
    }
  }

  list(nll = nll2, grad = grad)
}

.adirmcInnerGrad <- function(p, pinfo, studies_snap, proposals)
  .adirmcNLLAndGrad(p, pinfo, studies_snap, proposals)$grad

# -- Inner NLL -----------------------------------------------------------------
# Deterministic given fixed proposals. Sequence:
# 1. IS log-weights: log p(bi|Omega_new, mean_new) - log p(bi|Omega_prop) -> softmax
#    mean_new encodes paired theta shift: log(back_fn(p_new)) - log(back_fn(p_orig))
# 2. Weighted mean + covariance (ML, no n-1)
# 3. Sigma added to V diagonal using pre-kappa mu (sigma_prop scales with pre-kappa mu^2)
# 4. Kappa correction: mu += kappa_fn(struct) - mu_pop
#    kappa_fn = f(theta,0) + (1/2) sum_k omega_k_prop * d2f/deta_k2(theta,0)  [2nd-order delta]
# 5. -2LL = n * (log|V| + tr(V_obs * V^{-1}) + r' V^{-1} r)

.adirmcInnerNLL <- function(pars, prop, s) {
  mean_new_paired <- if (length(prop$log_origbeta) > 0L)
    compute_mean_new_cpp(pars$struct[names(prop$log_origbeta)],
                         prop$log_origbeta, prop$paired_types,
                         prop$paired_lows,  prop$paired_his)
  else numeric(0)

  # Expand to full n_eta vector: 0 for unpaired etas, shift for paired etas.
  n_eta_bi <- ncol(prop$bi)
  mean_new <- rep(0.0, n_eta_bi)
  if (length(mean_new_paired) > 0L)
    mean_new[prop$paired_eta_pos] <- mean_new_paired

  L_omega <- if (!is.null(pars$L)) pars$L else matrix(0, 0, 0)

  kappa_delta <- if (!is.null(prop$mu_pop))
    prop$kappa_fn(pars$struct) - prop$mu_pop
  else numeric(0)

  irmc_inner_nll_cpp(
    prop$rawpreds, prop$bi, mean_new, L_omega, prop$log_prop,
    s$E, s$V, s$n,
    pars$sigma_var, prop$sigma_type,
    kappa_delta,
    as.integer(identical(s$method, "var"))
  )
}

# -- Proposal generation -------------------------------------------------------

# has_kappa  -- whether any struct theta is non-mu-referenced (single_beta)
# use_grad   -- additionally build kappa_fn_batch for exact kappa (needed only for analytical gradient)
# Approach: IS weight-shift for paired (mu-referenced) thetas always; kappa only for single_betas.
# When has_kappa=TRUE, mu_pop computed by appending an eta=0 row to the main rxSolve batch.
.adirmcProposal <- function(rxMod, struct_theta, sigma_names, sigma_is_prop, sigma_is_lnorm,
                          omega, omega_expansion,
                          study, z, output_var, params_df, cores,
                          eta_col_names, has_kappa = FALSE,
                          kappa_method = "exact",
                          kappa_n_nodes = 5L,
                          struct_transforms = NULL, struct_eta_idx = NULL,
                          use_grad = TRUE, ndp = .Machine$integer.max) {
  Omega_prop <- omega * omega_expansion
  L_prop     <- tryCatch(t(chol(Omega_prop)), error = function(e) {
    message("adirmc: proposal:chol(Omega_prop) failed: ", conditionMessage(e))
    NULL
  })
  if (is.null(L_prop)) return(NULL)

  bi           <- z %*% t(L_prop)
  colnames(bi) <- eta_col_names

  for (nm in names(struct_theta)) params_df[, nm]          <- struct_theta[nm]
  if (ncol(bi) > 0L)              params_df[, eta_col_names] <- bi
  for (nm in sigma_names)         params_df[, nm]            <- 0

  # Build kappa rows and append to main batch (avoids a separate rxSolve).
  # "linearized":  1 f0 row + n_s FD rows at eta=0 (J precomputed once per outer iter).
  # "exact":        1 f0 row only at eta=0 (kappa_fn re-solves per inner step).
  # Identify single_betas: struct thetas not in any mu-referenced eta pairing.
  # Paired struct thetas get IS weight-shift; single_betas get kappa correction.
  n_eta_prop       <- ncol(bi)
  paired_s_idxs    <- if (!is.null(struct_eta_idx) && length(struct_eta_idx) == n_eta_prop)
    unique(struct_eta_idx[!is.na(struct_eta_idx)]) else seq_along(struct_theta)
  single_beta_idxs <- setdiff(seq_along(struct_theta), paired_s_idxs)
  single_beta_nms  <- names(struct_theta)[single_beta_idxs]
  has_single_betas <- length(single_beta_idxs) > 0L

  params_df_1        <- NULL
  h_jac              <- NULL
  n_s_single         <- 0L
  n_gh               <- 1L
  W_gh               <- 1.0
  struct0_opt_single <- NULL
  if (has_kappa && has_single_betas) {
    params_df_1 <- params_df[1L, , drop = FALSE]
    if (length(eta_col_names) > 0L) params_df_1[1L, eta_col_names] <- 0

    if (kappa_method == "linearized") {
      n_s_single         <- length(single_beta_nms)
      struct0_opt_single <- unlist(struct_theta)[single_beta_nms]
      h_jac              <- rep(1e-4, n_s_single)

      params_fd <- params_df_1[rep(1L, n_s_single), , drop = FALSE]
      for (j in seq_len(n_s_single))
        params_fd[j, single_beta_nms[j]] <- struct0_opt_single[j] + h_jac[j]

      params_df_ext <- rbind(params_df, params_df_1, params_fd)

    } else if (kappa_method == "linearized_gh") {
      n_s_single         <- length(single_beta_nms)
      struct0_opt_single <- unlist(struct_theta)[single_beta_nms]
      h_jac              <- rep(1e-4, n_s_single)

      gh_grid <- .adghNodeGrid(kappa_n_nodes, n_eta_prop)
      n_gh    <- length(gh_grid$W)
      W_gh    <- gh_grid$W
      L_base  <- tryCatch(t(chol(omega)), error = function(e) diag(max(1L, nrow(omega))))
      eta_gh  <- if (n_eta_prop > 0L) gh_grid$X %*% t(L_base) else matrix(0, n_gh, 0L)

      # Baseline: n_gh rows with GH eta, theta at current values
      params_gh_base <- params_df[rep(1L, n_gh), , drop = FALSE]
      for (nm in names(struct_theta)) params_gh_base[, nm] <- struct_theta[nm]
      if (n_eta_prop > 0L) params_gh_base[, eta_col_names] <- eta_gh
      for (nm in sigma_names) params_gh_base[, nm] <- 0

      # FD rows: for each single_beta, n_gh rows with that param perturbed
      params_gh_fd <- params_df[rep(1L, n_s_single * n_gh), , drop = FALSE]
      for (j in seq_len(n_s_single)) {
        row_ids <- (j - 1L) * n_gh + seq_len(n_gh)
        for (nm in names(struct_theta)) params_gh_fd[row_ids, nm] <- struct_theta[nm]
        params_gh_fd[row_ids, single_beta_nms[j]] <- struct0_opt_single[j] + h_jac[j]
        if (n_eta_prop > 0L) params_gh_fd[row_ids, eta_col_names] <- eta_gh
        for (nm in sigma_names) params_gh_fd[row_ids, nm] <- 0
      }

      params_df_ext <- rbind(params_df, params_gh_base, params_gh_fd)

    } else {
      params_df_ext <- rbind(params_df, params_df_1)
    }
  } else {
    params_df_ext <- params_df
  }

  out      <- rxode2::rxSolve(rxMod, params = as.data.frame(params_df_ext),
                              events = study$ev_full,
                              cores = cores,
                              nDisplayProgress = ndp)
  keep     <- out[["time"]] %in% study$times
  obs_vals <- out[[output_var]][keep]
  if (is.null(obs_vals)) obs_vals <- out[["ipredSim"]][keep]

  n_sim      <- nrow(z)
  n_times    <- length(study$times)
  n_kappa_rows <- if (!has_kappa || !has_single_betas) 0L else
    if (kappa_method == "linearized") 1L + n_s_single else
    if (kappa_method == "linearized_gh") n_gh * (1L + n_s_single) else
    1L
  n_expected   <- (n_sim + n_kappa_rows) * n_times
  if (length(obs_vals) != n_expected) {
    message(sprintf("adirmc: proposal:rxSolve returned %d values, expected %d",
                    length(obs_vals), n_expected))
    return(NULL)
  }

  rawpreds <- matrix(obs_vals[seq_len(n_sim * n_times)],
                     nrow = n_sim, ncol = n_times, byrow = TRUE)

  mu_pop         <- NULL
  kappa_fn       <- NULL
  kappa_jac      <- NULL
  kappa_fn_batch <- NULL
  if (has_kappa && has_single_betas) {

    if (kappa_method == "linearized_gh") {
      # GH-averaged baseline: mu_pop = sum_q W_q f(theta, eta_q)
      gh_base_vals <- matrix(obs_vals[n_sim * n_times + seq_len(n_gh * n_times)],
                             nrow = n_gh, ncol = n_times, byrow = TRUE)
      mu_pop <- drop(crossprod(W_gh, gh_base_vals))

      # GH-averaged FD Jacobian: each single_beta has n_gh rows in the batch
      kappa_jac <- matrix(0, n_s_single, n_times)
      for (j in seq_len(n_s_single)) {
        fd_offset <- n_sim * n_times + n_gh * n_times + (j - 1L) * n_gh * n_times
        fd_gh_j   <- matrix(obs_vals[fd_offset + seq_len(n_gh * n_times)],
                             nrow = n_gh, ncol = n_times, byrow = TRUE)
        kappa_jac[j, ] <- (drop(crossprod(W_gh, fd_gh_j)) - mu_pop) / h_jac[j]
      }

      kappa_fn <- local({
        .kappa_jac <- kappa_jac
        .mu_pop    <- mu_pop
        .struct0   <- struct0_opt_single
        .names     <- single_beta_nms
        function(struct_cand) {
          delta <- unlist(struct_cand)[.names] - .struct0
          .mu_pop + drop(delta %*% .kappa_jac)
        }
      })

    } else {
      mu_pop <- obs_vals[n_sim * n_times + seq_len(n_times)]
    }

    if (kappa_method == "linearized") {
      # Linearised kappa: vary single_betas only; J has n_s_single rows.
      fd_start  <- n_sim * n_times + n_times
      fd_vals   <- obs_vals[fd_start + seq_len(n_s_single * n_times)]
      fd_mat    <- matrix(fd_vals, nrow = n_s_single, ncol = n_times, byrow = TRUE)
      kappa_jac <- sweep(sweep(fd_mat, 2L, mu_pop, "-"), 1L, h_jac, "/")

      kappa_fn <- local({
        .kappa_jac <- kappa_jac
        .mu_pop    <- mu_pop
        .struct0   <- struct0_opt_single
        .names     <- single_beta_nms
        function(struct_cand) {
          delta <- unlist(struct_cand)[.names] - .struct0
          .mu_pop + drop(delta %*% .kappa_jac)
        }
      })

    } else if (kappa_method != "linearized_gh") {
      # Exact kappa: re-evaluate f(theta_single, 0) per inner step; paired cols stay at orig.
      # Matches admr: beta_use <- origbeta; beta_use[single_betas] <- pnew[single_betas].
      kappa_fn <- local({
        .params_base <- params_df_1
        .single_nms  <- single_beta_nms
        .cache_key   <- NULL
        .cache_val   <- NULL
        function(struct_cand) {
          if (!is.null(.cache_key) && identical(struct_cand, .cache_key))
            return(.cache_val)
          params_cand <- .params_base[1L, , drop = FALSE]
          for (nm in .single_nms)
            params_cand[1L, nm] <- struct_cand[[nm]]
          out_c  <- rxode2::rxSolve(rxMod, params = as.data.frame(params_cand),
                                    events = study$ev_full, cores = 1L,
                                    nDisplayProgress = ndp)
          keep_c  <- out_c[["time"]] %in% study$times
          vals_c1 <- out_c[[output_var]][keep_c]
          if (is.null(vals_c1)) vals_c1 <- out_c[["ipredSim"]][keep_c]
          val <- as.numeric(vals_c1)
          .cache_key <<- struct_cand
          .cache_val <<- val
          val
        }
      })

      if (use_grad) {
        kappa_fn_batch <- local({
          .params_base <- params_df_1
          .single_nms  <- single_beta_nms
          function(struct_list) {
            n_cands    <- length(struct_list)
            params_bat <- .params_base[rep(1L, n_cands), , drop = FALSE]
            for (ci in seq_len(n_cands)) {
              sc <- struct_list[[ci]]
              for (nm in .single_nms)
                params_bat[ci, nm] <- sc[[nm]]
            }
            out_b  <- rxode2::rxSolve(rxMod, params = as.data.frame(params_bat),
                                      events = study$ev_full, cores = 1L,
                                      nDisplayProgress = ndp)
            keep_b  <- out_b[["time"]] %in% study$times
            vals_b1 <- out_b[[output_var]][keep_b]
            if (is.null(vals_b1)) vals_b1 <- out_b[["ipredSim"]][keep_b]
            matrix(vals_b1, nrow = n_cands, ncol = n_times, byrow = TRUE)
          }
        })
      }
    }
  }

  # Weight-shift IS fields: always computed for paired (mu-referenced) struct thetas.
  # paired_eta_pos / log_origbeta are empty when no etas are mu-referenced (all single_betas).
  # Single_betas are corrected via kappa, not weight-shift.
  valid_idx_prop    <- if (!is.null(struct_eta_idx) && length(struct_eta_idx) == n_eta_prop)
    struct_eta_idx[!is.na(struct_eta_idx)] else seq_len(n_eta_prop)
  paired_struct_nms <- names(struct_theta)[valid_idx_prop]

  paired_eta_pos <- if (!is.null(struct_eta_idx) && length(struct_eta_idx) == n_eta_prop)
    which(!is.na(struct_eta_idx))
  else seq_len(n_eta_prop)

  .type_code <- function(nm) {
    cv <- if (!is.null(struct_transforms) && !is.null(struct_transforms[[nm]]))
      struct_transforms[[nm]]$curEval else "exp"
    switch(cv, exp = 0L, log = 0L, expit = 1L, logit = 1L,
           probitInv = 2L, probit = 2L, 3L)
  }
  .bound <- function(nm, field) {
    if (!is.null(struct_transforms) && !is.null(struct_transforms[[nm]]))
      struct_transforms[[nm]][[field]] else NA_real_
  }

  log_origbeta <- vapply(paired_struct_nms, function(nm) {
    .admLogBackTransform(struct_theta[[nm]], struct_transforms[[nm]])
  }, double(1))
  paired_types <- vapply(paired_struct_nms, .type_code, integer(1))
  paired_lows  <- vapply(paired_struct_nms, .bound, double(1), field = "low")
  paired_his   <- vapply(paired_struct_nms, .bound, double(1), field = "hi")

  list(rawpreds         = rawpreds,
       bi               = bi,
       paired_eta_pos   = paired_eta_pos,
       log_origbeta     = log_origbeta,
       paired_types     = paired_types,
       paired_lows      = paired_lows,
       paired_his       = paired_his,
       log_prop         = logdmvnorm_batch_cpp(bi, rep(0, ncol(bi)), L_prop),
       sigma_type       = ifelse(sigma_is_lnorm, 2L, ifelse(sigma_is_prop, 1L, 0L)),
       mu_pop           = mu_pop,
       kappa_fn         = kappa_fn,
       kappa_jac        = kappa_jac,
       kappa_fn_batch   = kappa_fn_batch,
       kappa_grad_idxs  = single_beta_idxs)
}

# -- Outer IRMC NLL ------------------------------------------------------------

.adirmcNLL <- function(p, pinfo, studies_snap, proposals) {
  pars <- tryCatch(.admUnpack(p, pinfo), error = function(e) NULL)
  if (is.null(pars)) return(Inf)
  if (pinfo$n_eta > 0 && any(diag(pars$omega) <= 0)) return(Inf)

  nll2 <- 0
  for (i in seq_along(studies_snap)) {
    nll2 <- nll2 + .adirmcInnerNLL(pars, proposals[[i]], studies_snap[[i]])
    if (!is.finite(nll2)) return(Inf)
  }
  nll2
}

# -- Shared outer phase loop ---------------------------------------------------

# Runs all IRMC phases for one restart. draw_proposals_inner / draw_proposals_exact
# are closures that already capture rxMod, studies, z_list, etc.
# Returns list(best_p, best_nll, global_iter, nll_trace, par_trace).
.adirmcPhaseLoop <- function(
  best_p, best_nll,
  phases, outer_iter, convcrit, max_worse, print_every, grad_mode,
  algorithm, ftol_rel, maxeval,
  ov_lower, ov_upper,
  pinfo, studies,
  draw_proposals_inner,
  draw_proposals_exact,
  print_progress = TRUE
) {
  phase_names      <- c("Wide", "Focused", "Fine-tuning", "Precision")
  global_iter      <- 0L
  nll_trace        <- numeric(0)
  par_trace        <- NULL
  last_opt_message <- ""

  for (ph_idx in seq_along(phases)) {
    ph_step <- phases[ph_idx]; worse_counter <- 0L; p_cur <- best_p
    if (print_progress)
      message(.admProgressPhase(ph_idx, phase_names[min(ph_idx, 4L)], ph_step, pinfo))

    for (iter_idx in seq_len(outer_iter)) {
      global_iter <- global_iter + 1L
      proposals   <- draw_proposals_inner(p_cur)
      if (is.null(proposals)) {
        if (print_progress) message("proposal generation failed -- skipping iter")
        next
      }
      if (grad_mode == "analytical") {
        # Memo cache: LBFGS always calls eval_f(p) then eval_grad_f(p) at the same p.
        # Fusing into one pass eliminates redundant recomputation of IS weights,
        # weighted mean/cov, sigma corrections, and kappa (~50% of inner work).
        .cache_p <- NULL; .cache_res <- NULL
        eval_fused <- function(p) {
          if (!identical(p, .cache_p)) {
            .cache_res <<- .adirmcNLLAndGrad(p, pinfo, studies, proposals)
            .cache_p   <<- p
          }
          .cache_res
        }
        eval_f          <- function(p) eval_fused(p)$nll
        eval_grad_inner <- function(p) eval_fused(p)$grad
      } else {
        eval_f <- function(p) .adirmcNLL(p, pinfo, studies, proposals)
        eval_grad_inner <- switch(grad_mode,
          fd = local({
            h <- 1e-6
            function(p) {
              f0 <- eval_f(p)
              g  <- numeric(length(p))
              for (k in seq_along(p)) { p_h <- p; p_h[k] <- p_h[k] + h; g[k] <- (eval_f(p_h) - f0) / h }
              g
            }
          }),
          NULL)
      }
      # After .admResolveAlgorithm: grad_mode == "none" <=> derivative-free
      # algorithm; grad_mode != "none" <=> gradient-based algorithm. Either way
      # the user's chosen algorithm matches the available gradient, so honour it
      # (the tryCatch below falls back to BOBYQA if the inner solve errors).
      algorithm_inner <- algorithm
      lb_inner <- pmax(ov_lower, p_cur - ph_step)
      ub_inner <- pmin(ov_upper, p_cur + ph_step)
      opt <- tryCatch(
        nloptr::nloptr(x0 = p_cur, eval_f = eval_f, eval_grad_f = eval_grad_inner,
                       lb = lb_inner, ub = ub_inner,
                       opts = list(algorithm = algorithm_inner,
                                   ftol_rel = ftol_rel, maxeval = maxeval)),
        error = function(e)
          nloptr::nloptr(x0 = p_cur, eval_f = eval_f, eval_grad_f = NULL,
                         lb = lb_inner, ub = ub_inner,
                         opts = list(algorithm = "NLOPT_LN_BOBYQA",
                                     ftol_rel = ftol_rel, maxeval = maxeval)))
      last_opt_message <- opt$message
      p_new <- opt$solution; nll_approx <- opt$objective
      props_exact <- draw_proposals_exact(p_new)
      nll_exact   <- if (!is.null(props_exact))
        .adirmcNLL(p_new, pinfo, studies, props_exact) else Inf
      if (nll_exact < best_nll) { best_nll <- nll_exact; best_p <- p_new }
      if (is.finite(nll_exact)) {
        nll_trace <- c(nll_trace, nll_exact)
        par_trace <- rbind(par_trace, p_new)
      }
      converged <- is.finite(nll_exact) && is.finite(nll_approx) &&
        abs(nll_exact - nll_approx) < convcrit
      if (converged) {
        if (print_progress) {
          row <- .admProgressRow(sprintf("%04d \u2713", global_iter), nll_exact, p_new, pinfo)
          if (!is.null(row)) message(row)
        }
        p_cur <- p_new; break
      }
      if (print_progress && print_every > 0L && iter_idx %% print_every == 0L) {
        row <- .admProgressRow(sprintf("%04d", global_iter), nll_exact, p_new, pinfo)
        if (!is.null(row)) message(row)
      }
      if (nll_exact > best_nll) {
        worse_counter <- worse_counter + 1L
        if (worse_counter >= max_worse) {
          if (print_progress) {
            row <- .admProgressRow(sprintf("%04d st", global_iter), nll_exact, p_new, pinfo)
            if (!is.null(row)) message(row)
          }
          p_cur <- best_p; break
        }
      } else { worse_counter <- 0L }
      p_cur <- p_new
    }
  }
  list(best_p = best_p, best_nll = best_nll, global_iter = global_iter,
       nll_trace = nll_trace, par_trace = par_trace,
       last_opt_message = last_opt_message)
}

# -- Restart worker ------------------------------------------------------------

.adirmcRestartWorker <- function(restart_id, p_init, ui_lstExpr, pinfo,
                               ov_lower, ov_upper, scale_c = NULL, studies, n_sim, seed,
                               phases, outer_iter, maxeval, ftol_rel,
                               algorithm, omega_expansion, convcrit,
                               max_worse, grad_mode = "none",
                               output_var = "cp",
                               kappa_method = "exact",
                               kappa_n_nodes = 5L,
                               sampling = "sobol",
                               print_progress = TRUE, print = 1L,
                               cores = NULL, no_lock = FALSE,
                               rxMod_direct = NULL) {
  library(admixr2)

  # Dev mode (PSOCK workers): patch installed namespace with dev functions from
  # .GlobalEnv (serialised there by furrr globals). tryCatch guards against
  # the installed package predating this function (run devtools::install() once).
  tryCatch(.admPatchDevNamespace(), error = function(e) NULL)

  # adirmc has no sensitivity model (analytical inner gradient) -> no sens_* args.
  m       <- .admWorkerLoadModels(ui_lstExpr, rxMod_direct, cores)
  cores_w <- m$cores_w
  rxMod   <- m$rxMod

  set.seed(seed)
  z_list      <- .admMakeZ(n_sim, pinfo, length(studies), sampling)
  set.seed(seed + restart_id)
  params_list <- .admMakeParamsList(n_sim, pinfo, length(studies))

  p_cur    <- p_init
  best_nll <- Inf
  best_p   <- p_cur
  t0_worker <- proc.time()
  needs_lock <- is.null(rxMod_direct) && !no_lock
  if (needs_lock) {
    tryCatch(rxode2::rxLock(rxMod), error = function(e) NULL)
    on.exit(tryCatch(rxode2::rxUnlock(rxMod), error = function(e) NULL), add = TRUE)
  }

  .draw_proposals_inner <- function(p) {
    pars  <- .admUnpack(p, pinfo)
    props <- lapply(seq_along(studies), function(si)
      .adirmcProposal(rxMod, pars$struct, pinfo$sigma_names,
                    pinfo$sigma_is_prop, pinfo$sigma_is_lnorm,
                    pars$omega, omega_expansion, studies[[si]],
                    z_list[[si]], output_var, params_list[[si]], cores_w,
                    pinfo$eta_col_names,
                    has_kappa         = pinfo$has_kappa,
                    kappa_method      = kappa_method,
                    kappa_n_nodes     = kappa_n_nodes,
                    struct_transforms = pinfo$struct_transforms,
                    struct_eta_idx    = pinfo$struct_eta_idx,
                    use_grad          = grad_mode == "analytical",
                    ndp               = pinfo$nDisplayProgress))
    if (any(vapply(props, is.null, logical(1)))) NULL else props
  }

  .draw_proposals_exact <- function(p) {
    pars  <- .admUnpack(p, pinfo)
    props <- lapply(seq_along(studies), function(si)
      .adirmcProposal(rxMod, pars$struct, pinfo$sigma_names,
                    pinfo$sigma_is_prop, pinfo$sigma_is_lnorm,
                    pars$omega, omega_expansion, studies[[si]],
                    z_list[[si]], output_var, params_list[[si]], cores_w,
                    pinfo$eta_col_names,
                    has_kappa         = pinfo$has_kappa,
                    kappa_method      = kappa_method,
                    kappa_n_nodes     = kappa_n_nodes,
                    struct_transforms = pinfo$struct_transforms,
                    struct_eta_idx    = pinfo$struct_eta_idx,
                    use_grad          = FALSE,
                    ndp               = pinfo$nDisplayProgress))
    if (any(vapply(props, is.null, logical(1)))) NULL else props
  }

  pl <- .adirmcPhaseLoop(
    best_p = best_p, best_nll = best_nll,
    phases = phases, outer_iter = outer_iter,
    convcrit = convcrit, max_worse = max_worse,
    print_every = print, grad_mode = grad_mode,
    algorithm = algorithm, ftol_rel = ftol_rel, maxeval = maxeval,
    ov_lower = ov_lower, ov_upper = ov_upper,
    pinfo = pinfo, studies = studies,
    draw_proposals_inner = .draw_proposals_inner,
    draw_proposals_exact = .draw_proposals_exact,
    print_progress = print_progress
  )

  list(restart_id = restart_id, objective = pl$best_nll,
       solution = pl$best_p, n_iter = pl$global_iter,
       nll_trace = pl$nll_trace, par_trace = pl$par_trace,
       elapsed = as.numeric((proc.time() - t0_worker)["elapsed"]),
       final_row_printed = TRUE)
}

# -- Main estimation entry point -----------------------------------------------

#' Fit an aggregate data model via Iterative Reweighting MC (adirmc estimator)
#'
#' Called automatically by
#' `nlmixr2(model, admData(), est = "adirmc", control = adirmcControl(...))`.
#' Not typically called directly.
#'
#' @param env nlmixr2 environment containing `ui` and `control`.
#' @param ... Unused.
#'
#' @return An `admFit` nlmixr2 fit object.
#'
#' @method nlmixr2Est adirmc
#' @export
nlmixr2Est.adirmc <- function(env, ...) {
  .ui  <- env$ui
  .ctl <- env$control

  if (!inherits(.ctl, "adirmcControl")) .ctl <- getValidNlmixrCtl.adirmc(.ctl)
  if (!inherits(.ctl, "adirmcControl"))
    stop("Could not recover adirmcControl", call. = FALSE)
  assign("control", .ctl, envir = .ui)

  studies <- .ctl$studies
  if (length(studies) == 0L)
    stop("adirmcControl(studies=...) required", call. = FALSE)
  if (is.null(names(studies)))
    names(studies) <- paste0("study", seq_along(studies))

  pinfo      <- .admParseIniDf(.ui$iniDf, .ui)
  pinfo$nDisplayProgress <- .ctl$nDisplayProgress
  output_var <- .admOutputVar(.ui)

  for (nm in names(studies))
    studies[[nm]] <- .admNormaliseStudy(studies[[nm]], nm, output_var)
  studies    <- .admFlattenStudies(studies)
  .adm_multi_out <- length(.admOutputVars(.ui)) > 1L
  .adm_joint     <- any(vapply(studies, function(u) isTRUE(u$is_joint), logical(1)))
  if (.adm_multi_out || .adm_joint)
    stop("adirmc does not yet support multiple observed outputs (multi-compartment observations). ",
         "Use est = 'admc', 'adfo', or 'adgh' for multi-output fits.", call. = FALSE)
  studies    <- .admBuildEvFull(studies)

  if (pinfo$has_kappa) {
    .unpaired <- names(pinfo$struct_has_eta)[!pinfo$struct_has_eta]
    message(sprintf(
      "adirmc: struct theta(s) not mu-referenced: %s. Using %s kappa correction.",
      paste(.unpaired, collapse = ", "), .ctl$kappa_method
    ))
  }

  # ORDERING INVARIANT: .admLoadSensModel() must run before .admLoadModel().
  # See model.R for rationale (linCmt foceiModel FD-path caching inner=NULL).
  sensModel <- if (.ctl$grad == "analytical")
    tryCatch(.admLoadSensModel(.ui), error = function(e) NULL)
  else NULL

  rxMod <- .admLoadModel(.ui)
  rxode2::rxLock(rxMod)
  on.exit({ rxode2::rxUnlock(rxMod); rxode2::rxSolveFree() }, add = TRUE)

  studies_snap <- studies

  set.seed(.ctl$seed)
  z_list      <- .admMakeZ(.ctl$n_sim, pinfo, length(studies), .ctl$sampling)
  params_list <- .admMakeParamsList(.ctl$n_sim, pinfo, length(studies))

  ov    <- .admBuildOptVec(pinfo)
  p_cur <- ov$p0
  cores <- .ctl$cores

  irmc_grad_label <- if (.ctl$grad == "none") "none" else {
    cov_label <- if (!is.null(sensModel)) "+Sens-Hessian" else "+FD-Hessian"
    grad_inner_label <- if (.ctl$grad == "fd") "forward FD" else "analytic"
    paste0(grad_inner_label, if (.ctl$covMethod == "r") cov_label else "")
  }
  message("=== admixr2: Aggregate Data Modeling (IR-MC) ===")
  message(sprintf("  Studies: %d | MC samples: %d | Phases: %d | Iters/phase: %d | Expansion: %.2f | Grad: %s | Restarts: %d",
                  length(studies), .ctl$n_sim, length(.ctl$phases),
                  .ctl$outer_iter, .ctl$omega_expansion,
                  irmc_grad_label, .ctl$n_restarts))
  t0 <- proc.time()

  best_nll  <- Inf
  best_p    <- p_cur
  n_iter    <- 0L
  nll_trace <- numeric(0)
  par_trace <- NULL

  grad_mode_inner <- if (.ctl$grad == "none") "none" else if (.ctl$grad == "fd") "fd" else "analytical"

  .draw_proposals_inner <- function(p) {
    pars <- .admUnpack(p, pinfo)
    props <- lapply(seq_along(studies_snap), function(si)
      .adirmcProposal(rxMod, pars$struct, pinfo$sigma_names,
                    pinfo$sigma_is_prop, pinfo$sigma_is_lnorm,
                    pars$omega, .ctl$omega_expansion,
                    studies_snap[[si]], z_list[[si]], output_var,
                    params_list[[si]], cores, pinfo$eta_col_names,
                    has_kappa         = pinfo$has_kappa,
                    kappa_method      = .ctl$kappa_method,
                    kappa_n_nodes     = .ctl$kappa_n_nodes,
                    struct_transforms = pinfo$struct_transforms,
                    struct_eta_idx    = pinfo$struct_eta_idx,
                    use_grad          = grad_mode_inner == "analytical",
                    ndp               = pinfo$nDisplayProgress))
    if (any(vapply(props, is.null, logical(1)))) NULL else props
  }

  .draw_proposals_exact <- function(p) {
    pars <- .admUnpack(p, pinfo)
    props <- lapply(seq_along(studies_snap), function(si)
      .adirmcProposal(rxMod, pars$struct, pinfo$sigma_names,
                    pinfo$sigma_is_prop, pinfo$sigma_is_lnorm,
                    pars$omega, .ctl$omega_expansion,
                    studies_snap[[si]], z_list[[si]], output_var,
                    params_list[[si]], cores, pinfo$eta_col_names,
                    has_kappa         = pinfo$has_kappa,
                    kappa_method      = .ctl$kappa_method,
                    kappa_n_nodes     = .ctl$kappa_n_nodes,
                    struct_transforms = pinfo$struct_transforms,
                    struct_eta_idx    = pinfo$struct_eta_idx,
                    use_grad          = FALSE,
                    ndp               = pinfo$nDisplayProgress))
    if (any(vapply(props, is.null, logical(1)))) NULL else props
  }

  if (.ctl$n_restarts == 1L) {
    message(.admProgressHeader(pinfo, bottom = FALSE))
    pl <- .adirmcPhaseLoop(
      best_p = best_p, best_nll = best_nll,
      phases = .ctl$phases, outer_iter = .ctl$outer_iter,
      convcrit = .ctl$convcrit, max_worse = .ctl$max_worse,
      print_every = .ctl$print, grad_mode = grad_mode_inner,
      algorithm = .ctl$algorithm, ftol_rel = .ctl$ftol_rel, maxeval = .ctl$maxeval,
      ov_lower = ov$lower, ov_upper = ov$upper,
      pinfo = pinfo, studies = studies_snap,
      draw_proposals_inner = .draw_proposals_inner,
      draw_proposals_exact = .draw_proposals_exact
    )
    best_p    <- pl$best_p
    best_nll  <- pl$best_nll
    nll_trace <- pl$nll_trace
    par_trace <- pl$par_trace
    n_iter    <- pl$global_iter
    message(.admProgressTimingRow((proc.time() - t0)["elapsed"], pinfo))
  } else {
    .irmc_old_plan <- .admSetupParallelPlan(.ctl, .ctl$n_restarts)
    if (!is.null(.irmc_old_plan)) on.exit(future::plan(.irmc_old_plan), add = TRUE)
    opt_restart <- .admRunRestarts(
      worker_fn  = .adirmcRestartWorker,
      p0         = p_cur, ov = ov, pinfo = pinfo,
      .ctl       = .ctl, ui = .ui, studies = studies_snap,
      extra_args = list(
        phases          = .ctl$phases,
        outer_iter      = .ctl$outer_iter,
        maxeval         = .ctl$maxeval,
        ftol_rel        = .ctl$ftol_rel,
        algorithm       = .ctl$algorithm,
        omega_expansion = .ctl$omega_expansion,
        convcrit        = .ctl$convcrit,
        max_worse       = .ctl$max_worse,
        grad_mode       = grad_mode_inner,
        output_var      = output_var,
        kappa_method    = .ctl$kappa_method,
        kappa_n_nodes   = .ctl$kappa_n_nodes,
        sampling        = .ctl$sampling,
        print_progress  = .ctl$print > 0L,
        print           = .ctl$print,
        cores           = .ctl$cores,
        rxMod_direct    = rxMod
      )
    )
    admStopWorkers()
    best_nll  <- opt_restart$objective
    best_p    <- opt_restart$solution
    n_iter    <- opt_restart$n_iter
  }

  t_opt     <- (proc.time() - t0)["elapsed"]
  final     <- .admUnpack(best_p, pinfo)
  fullTheta <- .admFullTheta(final, pinfo)

  all_traces <- if (.ctl$n_restarts == 1L) {
    list(list(restart_id = 1L, nll_trace = nll_trace, par_trace = par_trace))
  } else {
    opt_restart$all_traces
  }

  p_hat_irmc <- setNames(best_p, names(ov$p0))
  t0_cov <- proc.time()
  .cov <- if (.ctl$covMethod == "r") {
    np_cov       <- length(pinfo$struct_names) + length(pinfo$sigma_names)
    use_grad_cov <- .ctl$grad != "none"
    n_evals <- if (use_grad_cov) {
      np_cov + 1L
    } else {
      n_off <- np_cov * (np_cov - 1L) / 2L
      np_cov * 2L + n_off * 4L + 1L
    }
    evals_label <- if (use_grad_cov) "gradient evaluations" else "NLL evaluations"
    hess_label  <- if (!use_grad_cov) "" else if (!is.null(sensModel))
      ", Sens-Hessian" else ", FD-Hessian"
    message(sprintf("  Computing covariance (R method, MC NLL%s, %d %s)",
                    hess_label, n_evals, evals_label))
    tryCatch(
      .admCalcCov(p_hat_irmc, pinfo, studies_snap, z_list, rxMod, output_var,
                  params_list, cores, cov_n_sim = .ctl$cov_n_sim,
                  use_grad = use_grad_cov, grad_h = .ctl$grad_h,
                  cov_h = .ctl$cov_h, cov_h_outer = .ctl$cov_h_outer,
                  sensModel = sensModel, sampling = .ctl$sampling),
      error = function(e) {
        warning("admCalcCov (adirmc) failed: ", conditionMessage(e))
        NULL
      })
  } else NULL
  t_cov     <- (proc.time() - t0_cov)["elapsed"]
  t_elapsed <- t_opt + t_cov

  if (.ctl$returnAdmr)
    return(list(objective = best_nll, fullTheta = fullTheta,
                struct = final$struct, sigma_var = final$sigma_var,
                omega = final$omega, nll_trace = nll_trace,
                cov = .cov))

  .ret           <- new.env(parent = emptyenv())
  .ret$table     <- env$table
  .ret$ui        <- .ui
  .ret$fullTheta <- fullTheta
  .ret$objective <- best_nll
  .ret$est       <- "adirmc"
  .ret$ofvType   <- "adirmc"
  .ret$adjObf    <- FALSE
  .ret$covMethod <- if (!is.null(.cov)) "r" else ""
  .ret$cov       <- .cov
  .ret$message   <- if (.ctl$n_restarts > 1L) opt_restart$message else pl$last_opt_message
  .ret$extra     <- ""
  .ret$origData  <- studies
  .ret$adirmcExtra <- list(struct         = final$struct,
                         sigma_var      = final$sigma_var,
                         sigma_is_prop  = pinfo$sigma_is_prop,
                         sigma_is_lnorm = pinfo$sigma_is_lnorm,
                         omega          = final$omega,
                         L              = final$L,
                         eta_col_names  = pinfo$eta_col_names,
                         par_names     = names(ov$p0),
                         npar          = length(ov$p0),
                         nll_trace     = nll_trace,
                         par_trace     = par_trace,
                         all_traces    = all_traces,
                         n_iter        = n_iter,
                         time          = t_elapsed,
                         t_opt         = t_opt,
                         t_cov         = t_cov,
                         studies       = studies,
                         n_sim         = .ctl$n_sim)

  nlmixr2est::.nlmixr2FitUpdateParams(.ret)
  nmObjHandleControlObject.adirmcControl(.ctl, .ret)
  if (exists("control", .ui)) rm(list = "control", envir = .ui)
  .ret$control <- .admToFoceiControl(.ctl)
  .focei_model <- suppressMessages(tryCatch(.ui$foceiModel, error = function(e) NULL))
  if (!is.null(.focei_model)) .ret$model <- .focei_model

  .fit <- nlmixr2est::nlmixr2CreateOutputFromUi(
    .ui, data = admData(), control = .ret$control,
    table = .ret$table, env = .ret, est = "adirmc")

  .fit$env$method    <- "adirmc"
  .fit$env$studies   <- studies
  .fit$env$adirmcExtra <- .ret$adirmcExtra
  # Populate nlmixr2-style parameter history so traceplot(fit) works natively.
  .admAttachParHist(.fit, .ret$adirmcExtra$all_traces, .ret$adirmcExtra$par_names, .ui)
  .old_cls <- class(.fit)
  .new_cls <- c("admFit", .old_cls)
  attr(.new_cls, ".foceiEnv") <- attr(.old_cls, ".foceiEnv")
  class(.fit) <- .new_cls

  .stats <- .admCalcObjStats(best_nll, length(ov$p0), studies)
  row.names(.stats$objDf) <- "adirmc"
  .fit$env$logLik    <- .stats$ll
  .fit$env$nobs      <- .stats$nobs
  .fit$env$objDf     <- .stats$objDf
  .fit$env$OBJF      <- .stats$objDf$OBJF
  .fit$env$AIC       <- .stats$objDf$AIC
  .fit$env$BIC       <- .stats$objDf$BIC
  .fit$env$objective <- best_nll
  .fit$env$time      <- data.frame(
    optimize   = t_opt,
    covariance = t_cov,
    other      = 0,
    elapsed    = t_elapsed,
    row.names  = NULL
  )

  .fit
}
