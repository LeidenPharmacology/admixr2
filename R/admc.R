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
#' @param n_sim Number of Monte Carlo samples per NLL evaluation.
#' @param sampling Sampling method for eta draws: `"sobol"` (Sobol, default),
#'   `"halton"` (Halton), `"torus"` (Kronecker/torus), `"lhs"` (Latin hypercube),
#'   or `"rnorm"` (iid normal).
#' @param algorithm nloptr algorithm string. Automatically switched to
#'   `"NLOPT_LD_LBFGS"` when `grad != "none"`.
#' @param maxeval Maximum number of optimizer function evaluations.
#' @param ftol_rel Relative function-value tolerance for convergence.
#' @param print Print progress every this many evaluations (0 = silent).
#' @param seed Random seed for reproducibility.
#' @param cores Number of OpenMP threads for `rxSolve()`.
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
#' @param covMethod Covariance method: `"r"` (numerical Hessian for structural
#'   and residual-error parameters only; omega/IIV SEs are not computed,
#'   consistent with nlmixr2 FOCEI) or `"none"`.
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
#'   runs restarts sequentially. Values `> 1` use a PSOCK cluster on Windows
#'   and fork workers on Unix/macOS. Workers are stopped automatically after
#'   the restart phase so all cores are available for the Hessian step.
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
    algorithm  = "NLOPT_LN_BOBYQA",
    maxeval    = 500L,
    ftol_rel   = .Machine$double.eps^2,
    print      = 10L,
    seed       = 12345L,
    cores      = 1L,
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
  checkmate::assertString(algorithm,                         .var.name = "algorithm")
  checkmate::assertIntegerish(maxeval, lower = 1L, len = 1, .var.name = "maxeval")
  checkmate::assertNumeric(ftol_rel,   lower = 0,  len = 1, .var.name = "ftol_rel")
  checkmate::assertIntegerish(print,   lower = 0L, len = 1, .var.name = "print")
  checkmate::assertIntegerish(seed,                len = 1, .var.name = "seed")
  checkmate::assertIntegerish(cores,   lower = 1L, len = 1, .var.name = "cores")
  checkmate::assertNumeric(grad_h,      lower = 0,  len = 1, .var.name = "grad_h")
  checkmate::assertNumeric(cov_h,       lower = 0, len = 1, .var.name = "cov_h")
  checkmate::assertNumeric(cov_h_outer, lower = 0, len = 1, .var.name = "cov_h_outer")
  checkmate::assertNumeric(grad_bounds, lower = 0,  len = 1, .var.name = "grad_bounds")
  covMethod <- match.arg(covMethod)
  checkmate::assertIntegerish(cov_n_sim,   lower = 1L, len = 1, .var.name = "cov_n_sim")
  checkmate::assertIntegerish(n_restarts,  lower = 1L, len = 1, .var.name = "n_restarts")
  checkmate::assertNumeric(restart_sd,     lower = 0,  len = 1, .var.name = "restart_sd")
  checkmate::assertIntegerish(workers,     lower = 1L, len = 1, .var.name = "workers")
  checkmate::assertNumeric(ci,         lower = 0, upper = 1, len = 1, .var.name = "ci")
  checkmate::assertIntegerish(sigdig,  lower = 1L, len = 1, .var.name = "sigdig")
  checkmate::assertLogical(returnAdmr,             len = 1, .var.name = "returnAdmr")

  if (grad != "none" && algorithm == "NLOPT_LN_BOBYQA")
    algorithm <- "NLOPT_LD_LBFGS"

  if (is.null(rxControl))   rxControl   <- rxode2::rxControl(sigdig = sigdig)
  if (is.null(sigdigTable)) sigdigTable <- max(round(sigdig), 3L)

  .ret <- list(
    studies       = studies,
    n_sim         = as.integer(n_sim),
    sampling      = sampling,
    algorithm     = algorithm,
    maxeval       = as.integer(maxeval),
    ftol_rel      = ftol_rel,
    print         = as.integer(print),
    seed          = as.integer(seed),
    cores         = as.integer(cores),
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
  sig_type <- ifelse(pinfo$sigma_is_lnorm, 2L, ifelse(pinfo$sigma_is_prop, 1L, 0L))
  for (i in seq_along(studies)) {
    s       <- studies[[i]]
    z       <- z_list[[i]]
    eta_mat <- z %*% t(pars$L)
    colnames(eta_mat) <- pinfo$eta_col_names
    cp_mat <- tryCatch(
      .admSimulate(rxMod, pars$struct, pinfo$sigma_names, eta_mat, s,
                   output_var, params_list[[i]], cores),
      error = function(e) NULL)
    if (is.null(cp_mat) || anyNA(cp_mat)) return(Inf)

    if (identical(s$method, "var")) {
      nll2 <- nll2 + nll_var_from_samples_cpp(cp_mat, as.numeric(s$E), s$v_diag,
                                               s$n, pars$sigma_var, sig_type)
    } else {
      nll2 <- nll2 + nll_cov_from_samples_cpp(cp_mat, as.numeric(s$E), s$V,
                                               s$n, pars$sigma_var, sig_type)
    }
    if (!is.finite(nll2)) return(Inf)
  }
  nll2
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
    z   <- z_list[[si]]
    pdf <- params_list[[si]]

    eta_mat           <- z %*% t(pars$L)
    colnames(eta_mat) <- eta_col_names

    unpaired_k <- which(vapply(pinfo$struct_names, function(nm)
      is.null(pinfo$struct_has_eta) || !isTRUE(pinfo$struct_has_eta[nm]), logical(1)))
    n_unp <- length(unpaired_k)

    use_sens <- !is.null(sensModel) && n_eta > 0
    batched_hi <- NULL; batched_lo <- NULL
    if (use_sens) {
      sens_out <- .admSimulateSens(sensModel, pars$struct, pinfo$sigma_names,
                                   eta_mat, s, cores)
      if (is.null(sens_out) || anyNA(sens_out$cp_mat)) {
        use_sens <- FALSE
      } else {
        cp_mat     <- sens_out$cp_mat
        dpred_list <- sens_out$dpred_list
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
      pdf_big[, "rxerr.cp"] <- 1L
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

      extra_b <- setdiff(rxMod$params, colnames(pdf_big))
      if (length(extra_b) > 0L)
        pdf_big <- cbind(pdf_big, matrix(0, nrow(pdf_big), length(extra_b),
                                         dimnames = list(NULL, extra_b)))
      out_b  <- rxode2::rxSolve(rxMod, params = as.data.frame(pdf_big),
                                 events = s$ev_full, cores = cores,
                                 nDisplayProgress = .Machine$integer.max)
      keep_b <- out_b[["time"]] %in% s$times
      vals_b <- out_b[[output_var]][keep_b]
      if (is.null(vals_b)) vals_b <- out_b[["ipredSim"]][keep_b]

      cp_mat <- matrix(vals_b[seq_len(n_sim * n_t)],
                       nrow = n_sim, ncol = n_t, byrow = TRUE)
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

    mu_sim <- colMeans(cp_mat)
    mu     <- mu_sim
    for (k in seq_along(pars$sigma_var))
      if (pinfo$sigma_is_lnorm[k])
        mu <- mu * exp(pars$sigma_var[k] / 2)
    cp_c <- sweep(cp_mat, 2L, mu_sim)
    r    <- as.numeric(s$E) - mu

    is_var <- identical(s$method, "var")
    if (is_var) {
      pv <- adm_col_sq_sum_cpp(cp_c) / n_sim
      for (k in seq_along(pars$sigma_var)) {
        sv <- pars$sigma_var[k]
        if (pinfo$sigma_is_prop[k])
          pv <- pv + sv * mu_sim^2
        else if (pinfo$sigma_is_lnorm[k])
          pv <- pv + mu^2 * (exp(sv) - 1)
        else
          pv <- pv + sv
      }
      dNLL_dmu     <- s$n * as.numeric(-2 * r / pv)
      dNLL_dV_diag <- s$n * (1 / pv - s$v_diag / pv^2 - r^2 / pv^2)
    } else {
      V <- crossprod(cp_c) / (n_sim - 1L)
      for (k in seq_along(pars$sigma_var)) {
        sv <- pars$sigma_var[k]
        if (pinfo$sigma_is_prop[k])
          diag(V) <- diag(V) + sv * mu_sim^2
        else if (pinfo$sigma_is_lnorm[k])
          diag(V) <- diag(V) + mu^2 * (exp(sv) - 1)
        else
          diag(V) <- diag(V) + sv
      }
      cholV <- tryCatch(chol(V), error = function(e) NULL)
      if (is.null(cholV)) return(rep(NA_real_, length(p)))
      invV         <- chol2inv(cholV)
      dNLL_dmu     <- s$n * as.numeric(-2 * invV %*% r)
      dNLL_dV      <- s$n * (invV - invV %*% (s$V + tcrossprod(r)) %*% invV)
      dNLL_dV_diag <- diag(dNLL_dV)
    }

    n_t   <- length(s$times)

    # sigma_mu_scale: error-V sensitivity w.r.t. mu_sim, reused across all gradient terms.
    # For lnorm, also folds in the residual scaling by exp(sv/2): (exp(sv/2)-1)*dNLL_dmu.
    # Computed once per study (depends only on pars, dNLL_dV_diag, dNLL_dmu, mu_sim, mu).
    sigma_mu_scale <- numeric(n_t)
    for (k in seq_along(pars$sigma_var)) {
      sv <- pars$sigma_var[k]
      if (pinfo$sigma_is_prop[k]) {
        sigma_mu_scale <- sigma_mu_scale + 2 * sv * dNLL_dV_diag * mu_sim
      } else if (pinfo$sigma_is_lnorm[k]) {
        sigma_mu_scale <- sigma_mu_scale +
          (exp(sv / 2) - 1) * dNLL_dmu +
          2 * exp(sv / 2) * mu * (exp(sv) - 1) * dNLL_dV_diag
      }
    }
    eff_dmu <- dNLL_dmu + sigma_mu_scale
    inv_nm1 <- 1 / (n_sim - 1L)

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
          dNLL_dV_diag, dNLL_dmu, sigma_mu_scale,
          neta1, neta2,
          n_t, n_eta)
      else
        adm_grad_eta_omega_cpp(
          cp_c, D_mat, z_diag_scale, z,
          dNLL_dV, dNLL_dmu, sigma_mu_scale,
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

    # Unpaired struct theta gradient via FD.
    # !use_sens path: batched_hi/batched_lo already extracted from the single big rxSolve above.
    # use_sens path: sens model handled etas; run a separate rxSolve for struct perturbations.
    if (n_unp > 0L) {
      if (!is.null(batched_hi)) {
        for (bi in seq_len(n_unp)) {
          k_s     <- unpaired_k[bi]
          cp_hi_s <- batched_hi[[bi]]
          dpred <- if (use_central && !is.null(batched_lo))
            (cp_hi_s - batched_lo[[bi]]) / (2 * h)
          else
            (cp_hi_s - cp_mat) / h
          grad[k_s] <- grad[k_s] +
            if (is_var)
              adm_grad_partial_var_cpp(cp_c, dpred, dNLL_dV_diag, eff_dmu, inv_nm1)
            else
              adm_grad_partial_cpp(cp_c, dpred, dNLL_dV, eff_dmu, inv_nm1)
        }
      } else {
        col_nms <- colnames(pdf)
        n_cols  <- length(col_nms)
        pdf_hi  <- matrix(0, nrow = n_unp * n_sim, ncol = n_cols,
                          dimnames = list(NULL, col_nms))
        pdf_hi[, "rxerr.cp"] <- 1L
        for (nm in names(pars$struct))       pdf_hi[, nm]              <- pars$struct[nm]
        for (j  in seq_along(eta_col_names)) pdf_hi[, eta_col_names[j]] <- rep(eta_mat[, j], n_unp)
        for (nm in pinfo$sigma_names)        pdf_hi[, nm]              <- 0
        for (bi in seq_len(n_unp)) {
          rows <- (bi - 1L) * n_sim + seq_len(n_sim)
          nm   <- pinfo$struct_names[unpaired_k[bi]]
          pdf_hi[rows, nm] <- pars$struct[nm] + h
        }
        extra_hi <- setdiff(rxMod$params, colnames(pdf_hi))
        if (length(extra_hi) > 0L)
          pdf_hi <- cbind(pdf_hi, matrix(0, nrow(pdf_hi), length(extra_hi),
                                          dimnames = list(NULL, extra_hi)))
        out_hi  <- rxode2::rxSolve(rxMod, params = as.data.frame(pdf_hi),
                                    events = s$ev_full, cores = cores,
                                    nDisplayProgress = .Machine$integer.max)
        keep_hi <- out_hi[["time"]] %in% s$times
        vals_hi <- out_hi[[output_var]][keep_hi]
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
                                      nDisplayProgress = .Machine$integer.max)
          keep_lo <- out_lo[["time"]] %in% s$times
          vals_lo <- out_lo[[output_var]][keep_lo]
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
              adm_grad_partial_var_cpp(cp_c, dpred, dNLL_dV_diag, eff_dmu, inv_nm1)
            else
              adm_grad_partial_cpp(cp_c, dpred, dNLL_dV, eff_dmu, inv_nm1)
        }
      }
    }

    k_sig <- n_s
    for (k in seq_along(pars$sigma_var)) {
      k_sig <- k_sig + 1L
      sv <- pars$sigma_var[k]
      if (pinfo$sigma_is_prop[k]) {
        grad[k_sig] <- grad[k_sig] + sum(dNLL_dV_diag * sv * mu_sim^2)
      } else if (pinfo$sigma_is_lnorm[k]) {
        grad[k_sig] <- grad[k_sig] + sv * (
          sum(dNLL_dV_diag * mu^2 * (2 * exp(sv) - 1)) +
          sum(dNLL_dmu * mu) / 2
        )
      } else {
        grad[k_sig] <- grad[k_sig] + sum(dNLL_dV_diag) * sv
      }
    }
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
  sig_type <- ifelse(pinfo$sigma_is_lnorm, 2L, ifelse(pinfo$sigma_is_prop, 1L, 0L))
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
    z <- z_list[[si]]
    n_t <- length(s$times)

    for (chunk in chunks) {
      n_chunk <- length(chunk)
      pdf_mat <- matrix(0, nrow = n_chunk * n_sim, ncol = n_cols,
                        dimnames = list(NULL, col_nms))
      pdf_mat[, "rxerr.cp"] <- 1L

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

      extra_m <- setdiff(rxMod$params, colnames(pdf_mat))
      if (length(extra_m) > 0L)
        pdf_mat <- cbind(pdf_mat, matrix(0, nrow(pdf_mat), length(extra_m),
                                         dimnames = list(NULL, extra_m)))
      out <- tryCatch(
        rxode2::rxSolve(rxMod, params = as.data.frame(pdf_mat),
                        events = s$ev_full, cores = cores,
                        nDisplayProgress = 1000L),
        error = function(e) NULL)
      if (is.null(out)) { for (ci in chunk) finite[ci] <- FALSE; next }

      keep <- out[["time"]] %in% s$times
      vals <- out[[output_var]][keep]
      if (is.null(vals)) vals <- out[["ipredSim"]][keep]

      for (cii in seq_along(chunk)) {
        ci <- chunk[cii]
        if (!valid[ci] || !finite[ci]) next
        pars <- pars_list[[ci]]
        idx  <- (cii - 1L) * n_sim * n_t + seq_len(n_sim * n_t)
        cp   <- matrix(vals[idx], nrow = n_sim, ncol = n_t, byrow = TRUE)
        if (anyNA(cp)) { finite[ci] <- FALSE; next }
        nll_ci <- if (identical(s$method, "var"))
          nll_var_from_samples_cpp(cp, as.numeric(s$E), s$v_diag,
                                    s$n, pars$sigma_var, sig_type)
        else
          nll_cov_from_samples_cpp(cp, as.numeric(s$E), s$V,
                                    s$n, pars$sigma_var, sig_type)
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

    eta_mats <- lapply(seq_len(n_c), function(ci) {
      if (!valid[ci]) return(NULL)
      em <- z %*% t(pars_list[[ci]]$L)
      colnames(em) <- eta_col_names
      em
    })

    cp_mats     <- vector("list", n_c)
    dpred_lists <- vector("list", n_c)

    # --- Sens model path -------------------------------------------------------
    use_sens <- !is.null(sensModel) && n_eta > 0L
    if (use_sens) {
      rmap      <- sensModel$rename_map
      all_src   <- c(pinfo$struct_names, pinfo$sigma_names, eta_col_names)
      inner_nms <- rmap[all_src]; inner_nms <- inner_nms[!is.na(inner_nms)]

      inner_df <- as.data.frame(matrix(0, nrow = n_c * n_sim,
                                        ncol = length(inner_nms),
                                        dimnames = list(NULL, unname(inner_nms))))
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
      }
      out <- tryCatch(
        suppressWarnings(
          rxode2::rxSolve(sensModel$mod, params = inner_df,
                          events = s$ev_full, cores = cores,
                          nDisplayProgress = 1000L)),
        error = function(e) NULL)
      if (is.null(out) || !all(sensModel$sens_cols %in% names(out))) {
        use_sens <- FALSE
      } else {
        keep      <- out[["time"]] %in% s$times
        vals_pred <- out[["rx_pred_"]][keep]
        vals_sens <- lapply(sensModel$sens_cols, function(col) out[[col]][keep])
        for (ci in seq_len(n_c)) {
          if (!valid[ci]) next
          idx <- (ci - 1L) * n_sim * n_t + seq_len(n_sim * n_t)
          cp_mats[[ci]]     <- matrix(vals_pred[idx], nrow = n_sim, ncol = n_t, byrow = TRUE)
          dpred_lists[[ci]] <- lapply(vals_sens, function(vs)
            matrix(vs[idx], nrow = n_sim, ncol = n_t, byrow = TRUE))
        }
      }
    }

    # --- Forward/central FD / no-eta fallback ----------------------------------
    if (!use_sens) {
      if (n_eta == 0L) {
        pdf_mat <- matrix(0, nrow = n_c * n_sim, ncol = n_cols,
                          dimnames = list(NULL, col_nms))
        pdf_mat[, "rxerr.cp"] <- 1L
        for (ci in seq_len(n_c)) {
          if (!valid[ci]) next
          rows <- (ci - 1L) * n_sim + seq_len(n_sim)
          pars <- pars_list[[ci]]
          for (nm in names(pars$struct)) pdf_mat[rows, nm] <- pars$struct[nm]
        }
        extra_m <- setdiff(rxMod$params, colnames(pdf_mat))
        if (length(extra_m) > 0L)
          pdf_mat <- cbind(pdf_mat, matrix(0, nrow(pdf_mat), length(extra_m),
                                           dimnames = list(NULL, extra_m)))
        out <- tryCatch(rxode2::rxSolve(rxMod, params = as.data.frame(pdf_mat),
                                         events = s$ev_full, cores = cores,
                                         nDisplayProgress = 1000L),
                        error = function(e) NULL)
        if (is.null(out)) { valid[] <- FALSE } else {
          keep <- out[["time"]] %in% s$times
          vals <- out[[output_var]][keep]
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
        pdf_mat[, "rxerr.cp"] <- 1L
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
        extra_m <- setdiff(rxMod$params, colnames(pdf_mat))
        if (length(extra_m) > 0L)
          pdf_mat <- cbind(pdf_mat, matrix(0, nrow(pdf_mat), length(extra_m),
                                           dimnames = list(NULL, extra_m)))
        out <- tryCatch(rxode2::rxSolve(rxMod, params = as.data.frame(pdf_mat),
                                         events = s$ev_full, cores = cores,
                                         nDisplayProgress = 1000L),
                        error = function(e) NULL)
        if (is.null(out)) { valid[] <- FALSE } else {
          keep <- out[["time"]] %in% s$times
          vals <- out[[output_var]][keep]
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
    if (n_unp > 0L) {
      for (ci in seq_len(n_c)) {
        cp_hi_store[[ci]] <- vector("list", n_unp)
        if (use_central) cp_lo_store[[ci]] <- vector("list", n_unp)
      }
      cu_idx <- expand.grid(bi = seq_len(n_unp), ci = seq_len(n_c))
      n_cu   <- nrow(cu_idx)
      pdf_hi <- matrix(0, nrow = n_cu * n_sim, ncol = n_cols,
                       dimnames = list(NULL, col_nms))
      pdf_hi[, "rxerr.cp"] <- 1L
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
      extra_hi <- setdiff(rxMod$params, colnames(pdf_hi))
      if (length(extra_hi) > 0L)
        pdf_hi <- cbind(pdf_hi, matrix(0, nrow(pdf_hi), length(extra_hi),
                                        dimnames = list(NULL, extra_hi)))
      out_hi <- tryCatch(rxode2::rxSolve(rxMod, params = as.data.frame(pdf_hi),
                                          events = s$ev_full, cores = cores,
                                          nDisplayProgress = 1000L),
                         error = function(e) NULL)
      if (!is.null(out_hi)) {
        vals_hi <- out_hi[[output_var]][out_hi[["time"]] %in% s$times]
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
                                            nDisplayProgress = 1000L),
                           error = function(e) NULL)
        if (!is.null(out_lo)) {
          vals_lo <- out_lo[[output_var]][out_lo[["time"]] %in% s$times]
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

      mu_sim <- colMeans(cp_mat)
      mu     <- mu_sim
      for (k in seq_along(pars$sigma_var))
        if (pinfo$sigma_is_lnorm[k])
          mu <- mu * exp(pars$sigma_var[k] / 2)
      cp_c <- sweep(cp_mat, 2L, mu_sim)
      r    <- as.numeric(s$E) - mu

      is_var <- identical(s$method, "var")
      if (is_var) {
        pv <- adm_col_sq_sum_cpp(cp_c) / n_sim
        for (k in seq_along(pars$sigma_var)) {
          sv <- pars$sigma_var[k]
          if (pinfo$sigma_is_prop[k])
            pv <- pv + sv * mu_sim^2
          else if (pinfo$sigma_is_lnorm[k])
            pv <- pv + mu^2 * (exp(sv) - 1)
          else
            pv <- pv + sv
        }
        dNLL_dmu     <- s$n * as.numeric(-2 * r / pv)
        dNLL_dV_diag <- s$n * (1 / pv - s$v_diag / pv^2 - r^2 / pv^2)
      } else {
        V <- crossprod(cp_c) / (n_sim - 1L)
        for (k in seq_along(pars$sigma_var)) {
          sv <- pars$sigma_var[k]
          if (pinfo$sigma_is_prop[k])
            diag(V) <- diag(V) + sv * mu_sim^2
          else if (pinfo$sigma_is_lnorm[k])
            diag(V) <- diag(V) + mu^2 * (exp(sv) - 1)
          else
            diag(V) <- diag(V) + sv
        }
        cholV <- tryCatch(chol(V), error = function(e) NULL)
        if (is.null(cholV)) { valid[ci] <- FALSE; next }
        invV         <- chol2inv(cholV)
        dNLL_dmu     <- s$n * as.numeric(-2 * invV %*% r)
        dNLL_dV      <- s$n * (invV - invV %*% (s$V + tcrossprod(r)) %*% invV)
        dNLL_dV_diag <- diag(dNLL_dV)
      }

      sigma_mu_scale <- numeric(n_t)
      for (k in seq_along(pars$sigma_var)) {
        sv <- pars$sigma_var[k]
        if (pinfo$sigma_is_prop[k]) {
          sigma_mu_scale <- sigma_mu_scale + 2 * sv * dNLL_dV_diag * mu_sim
        } else if (pinfo$sigma_is_lnorm[k]) {
          sigma_mu_scale <- sigma_mu_scale +
            (exp(sv / 2) - 1) * dNLL_dmu +
            2 * exp(sv / 2) * mu * (exp(sv) - 1) * dNLL_dV_diag
        }
      }
      eff_dmu <- dNLL_dmu + sigma_mu_scale
      inv_nm1 <- 1 / (n_sim - 1L)

      if (n_eta > 0L) {
        D_mat        <- do.call(cbind, dpred_list)
        eta_rows_df  <- pinfo$eta_rows_df
        z_diag_scale <- sweep(z, 2L, diag(pars$L) / 2, "*")
        neta1 <- as.integer(eta_rows_df$neta1)
        neta2 <- as.integer(eta_rows_df$neta2)
        go <- if (is_var)
          adm_grad_eta_omega_var_cpp(
            cp_c, D_mat, z_diag_scale, z,
            dNLL_dV_diag, dNLL_dmu, sigma_mu_scale,
            neta1, neta2,
            n_t, n_eta)
        else
          adm_grad_eta_omega_cpp(
            cp_c, D_mat, z_diag_scale, z,
            dNLL_dV, dNLL_dmu, sigma_mu_scale,
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

      for (bi in seq_len(n_unp)) {
        k_s     <- unpaired_k[bi]
        cp_hi_s <- cp_hi_store[[ci]][[bi]]
        if (is.null(cp_hi_s)) next
        dpred <- if (use_central && !is.null(cp_lo_store[[ci]][[bi]])) {
          (cp_hi_s - cp_lo_store[[ci]][[bi]]) / (2 * h)
        } else {
          (cp_hi_s - cp_mat) / h
        }
        grad_acc[ci, k_s] <- grad_acc[ci, k_s] +
          if (is_var)
            adm_grad_partial_var_cpp(cp_c, dpred, dNLL_dV_diag, eff_dmu, inv_nm1)
          else
            adm_grad_partial_cpp(cp_c, dpred, dNLL_dV, eff_dmu, inv_nm1)
      }

      k_sig <- n_s
      for (k in seq_along(pars$sigma_var)) {
        k_sig <- k_sig + 1L
        sv <- pars$sigma_var[k]
        if (pinfo$sigma_is_prop[k]) {
          grad_acc[ci, k_sig] <- grad_acc[ci, k_sig] + sum(dNLL_dV_diag * sv * mu_sim^2)
        } else if (pinfo$sigma_is_lnorm[k]) {
          grad_acc[ci, k_sig] <- grad_acc[ci, k_sig] + sv * (
            sum(dNLL_dV_diag * mu^2 * (2 * exp(sv) - 1)) +
            sum(dNLL_dmu * mu) / 2
          )
        } else {
          grad_acc[ci, k_sig] <- grad_acc[ci, k_sig] + sum(dNLL_dV_diag) * sv
        }
      }
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

  # Hessian restricted to struct + sigma only (omega Cholesky excluded).
  # Matches nlmixr2 FOCEI: omega entries are in the optimizer but skipped for cov.
  n_s      <- length(pinfo$struct_names)
  n_e      <- length(pinfo$sigma_names)
  cov_idx  <- seq_len(n_s + n_e)
  np_cov   <- length(cov_idx)
  nms_cov  <- nms[cov_idx]
  message("  Note: covMethod='r' computes covariance for structural and sigma ",
          "parameters only; omega (IIV) SEs are not computed (matching nlmixr2 ",
          "FOCEI behavior).")

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
  struct_nms <- pinfo$struct_names
  cov_full[struct_nms, struct_nms, drop = FALSE]
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

  # Dev mode (load_all): furrr serializes updated functions into worker .GlobalEnv.
  # Patch the installed namespace so all downstream calls use the dev versions.
  .adm_dev_nms <- ls(envir = .GlobalEnv, all.names = TRUE,
                     pattern = "^\\.(adm|adfo|adirmc|softmax|logdmvnorm)")
  if (length(.adm_dev_nms) > 0L) {
    .adm_ns <- asNamespace("admixr2")
    for (.adm_nm in .adm_dev_nms) {
      .adm_fn <- get(.adm_nm, envir = .GlobalEnv, inherits = FALSE)
      if (is.function(.adm_fn))
        tryCatch(utils::assignInNamespace(.adm_nm, .adm_fn, ns = .adm_ns),
                 error = function(e) NULL)
    }
    rm(.adm_dev_nms, .adm_ns, .adm_nm, .adm_fn)
  }

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

  sensModel <- if (!is.null(sensModel_direct)) {
    sensModel_direct
  } else if (!is.null(sens_cache_file) && file.exists(sens_cache_file)) {
    .smod <- tryCatch({ m <- qs2::qs_read(sens_cache_file); rxode2::rxLoad(m); m },
                      error = function(e) NULL)
    if (!is.null(.smod)) list(type = "ode", mod = .smod,
                               sens_cols = sens_cols, rename_map = sens_rename)
    else NULL
  } else {
    NULL
  }

  set.seed(seed)
  z_list      <- .admMakeZ(n_sim, pinfo, length(studies), sampling)
  set.seed(seed + restart_id)
  params_list <- .admMakeParamsList(n_sim, pinfo, length(studies))

  .iter      <- 0L
  .best_nll  <- Inf
  .nll_trace <- numeric(0)
  .par_trace <- NULL
  eval_f <- function(p) {
    .iter <<- .iter + 1L
    val <- .admNLL(p, pinfo, studies, z_list, rxMod, output_var, params_list, cores_w)
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
  eval_grad_f <- if (use_grad) {
    function(p) .admGrad(p, pinfo, studies, z_list, rxMod, output_var,
                         params_list, cores_w, grad_h, sensModel,
                         use_central = use_central)
  } else NULL

  lb <- if (use_grad) pmax(ov_lower, p_init - grad_bounds) else ov_lower
  ub <- if (use_grad) pmin(ov_upper, p_init + grad_bounds) else ov_upper

  sc <- if (!is.null(scale_c)) scale_c else rep(1.0, length(p_init))
  p_sc  <- p_init / sc
  lb_sc <- lb     / sc
  ub_sc <- ub     / sc
  eval_f_sc    <- function(p_s) eval_f(p_s * sc)
  eval_grad_sc <- if (!is.null(eval_grad_f)) {
    function(p_s) eval_grad_f(p_s * sc) * sc
  } else NULL

  # rxLock is a system-wide named mutex -- do NOT use across PSOCK worker processes,
  # each of which has its own independent model instance.
  needs_lock <- is.null(rxMod_direct) && !no_lock
  t0 <- proc.time()
  if (needs_lock) tryCatch(rxode2::rxLock(rxMod), error = function(e) NULL)
  opt <- tryCatch(
    nloptr::nloptr(
      x0 = p_sc, eval_f = eval_f_sc,
      eval_grad_f = eval_grad_sc,
      lb = lb_sc, ub = ub_sc,
      opts = list(algorithm = algorithm, ftol_rel = ftol_rel, maxeval = maxeval)
    ),
    finally = if (needs_lock)
      tryCatch(rxode2::rxUnlock(rxMod), error = function(e) NULL)
  )
  list(restart_id = restart_id, objective = opt$objective,
       solution = opt$solution * sc, n_iter = .iter,
       nll_trace = .nll_trace, par_trace = .par_trace,
       elapsed = as.numeric((proc.time() - t0)["elapsed"]),
       message = opt$message)
}

# -- Multi-restart orchestration -----------------------------------------------

# Persistent PSOCK worker pool -- created once, reused across fits in the same
# session. Fork workers (Unix/macOS) are always ephemeral and not stored here.
# The future::plan is set per-fit and restored afterwards (cheap: no worker
# spawning), so the cluster plan never leaks into unrelated nlmixr2() calls.
.adm_worker_env <- new.env(parent = emptyenv())
.adm_worker_env$cluster <- NULL   # parallel::makeCluster result
.adm_worker_env$n       <- 0L     # number of workers in cluster
# Clean up on R exit (onexit = TRUE) or when env is GC'd
reg.finalizer(.adm_worker_env, function(e) {
  if (!is.null(e$cluster))
    tryCatch(parallel::stopCluster(e$cluster), error = function(err) NULL)
}, onexit = TRUE)

#' Stop PSOCK workers
#'
#' Stops any PSOCK worker processes started by a parallel-restart fit
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
  if (is.null(.adm_worker_env$cluster)) {
    message("No persistent admixr2 workers running.")
    return(invisible(NULL))
  }
  n <- .adm_worker_env$n
  tryCatch(parallel::stopCluster(.adm_worker_env$cluster), error = function(e) NULL)
  .adm_worker_env$cluster <- NULL
  .adm_worker_env$n       <- 0L
  message(sprintf("%d admixr2 worker(s) stopped.", n))
  invisible(NULL)
}

# Set up parallel plan for restarts. Returns old_plan for on.exit restore.
# Workers are stopped by the caller (via admStopWorkers()) after restarts
# complete so all cores are free for the Hessian step.
.admSetupParallelPlan <- function(.ctl, n_r) {
  if (.ctl$workers <= 1L) return(invisible(NULL))
  if (!requireNamespace("future", quietly = TRUE))
    stop("Package 'future' required for workers > 1. Install it with install.packages('future').",
         call. = FALSE)
  if (!requireNamespace("furrr", quietly = TRUE))
    stop("Package 'furrr' required for workers > 1. Install it with install.packages('furrr').",
         call. = FALSE)
  if (.ctl$workers > n_r)
    warning(sprintf(
      "admControl(workers=%d) exceeds n_restarts=%d: only %d worker process(es) will be started",
      .ctl$workers, n_r, n_r
    ), call. = FALSE)

  old_plan <- future::plan()

  # Fork (Unix/macOS): ephemeral workers, cleaned up automatically.
  if (future::supportsMulticore()) {
    future::plan(future::multicore, workers = min(.ctl$workers, n_r))
    return(old_plan)
  }

  # PSOCK (Windows): always create a fresh cluster; stopped after restarts.
  n_w <- .ctl$workers
  if (!is.null(.adm_worker_env$cluster))
    tryCatch(parallel::stopCluster(.adm_worker_env$cluster), error = function(e) NULL)
  message(sprintf("  Starting %d PSOCK worker(s)", n_w))
  .adm_worker_env$cluster <- parallel::makeCluster(n_w)
  .adm_worker_env$n       <- n_w

  # Use workers= (not cluster=) -- future::cluster ignores cluster= and falls
  # back to detectCores() without it, making nbrOfWorkers() wrong.
  future::plan(future::cluster, workers = .adm_worker_env$cluster)
  old_plan
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
    .fn_names <- .fn_names[grepl("^\\.(adm|adfo|adirmc|softmax|logdmvnorm)", .fn_names)]
    .fn_list  <- setNames(lapply(.fn_names, get, envir = pkg_env), .fn_names)
    .fn_list[[.worker_fn_name]] <- worker_fn
  }

  use_parallel <- n_r > 1L &&
    .ctl$workers > 1L &&
    requireNamespace("furrr",  quietly = TRUE) &&
    requireNamespace("future", quietly = TRUE) &&
    future::nbrOfWorkers() > 1L

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
    n_workers         <- future::nbrOfWorkers()
    effective_workers <- min(n_workers, n_r)
    base_tpw          <- max(1L, floor(.ctl$cores / effective_workers))
    remainder         <- .ctl$cores - base_tpw * effective_workers
    cores_vec         <- c(rep(base_tpw + 1L, remainder), rep(base_tpw, effective_workers - remainder))
    tpw_label         <- if (remainder > 0L)
      sprintf("%d-%d", base_tpw, base_tpw + 1L) else as.character(base_tpw)
    n_batches         <- ceiling(n_r / effective_workers)
    batch_label       <- if (n_batches > 1L) sprintf(", %d sequential batch(es)", n_batches) else ""

    # use_fork when multicore is supported (Unix/macOS outside RStudio).
    use_fork <- .ctl$workers > 1L && future::supportsMulticore()

    message(sprintf("  Running %d restarts in parallel (%d workers, %s threads/worker%s%s)",
                    n_r, effective_workers, tpw_label, batch_label,
                    if (.ctl$workers > 1L) sprintf(", %s", if (use_fork) "fork" else "PSOCK") else ""))
    message(.admProgressHeader(pinfo, bottom = FALSE))

    # Prepare shared parallel args (common to all paths below)
    if (use_fork) {
      all_args_par <- all_args
      all_args_par$no_lock <- TRUE
      if ("print_progress" %in% names(all_args_par)) all_args_par$print_progress <- FALSE
    } else {
      # PSOCK: compiled DLLs cannot serialize; reload from qs2 cache.
      all_args_par <- all_args
      all_args_par$cores   <- base_tpw
      all_args_par$no_lock <- TRUE
      if ("print_progress"   %in% names(all_args_par)) all_args_par$print_progress   <- FALSE
      if ("rxMod_direct"     %in% names(all_args_par)) all_args_par$rxMod_direct     <- NULL
      if ("sensModel_direct" %in% names(all_args_par)) all_args_par$sensModel_direct <- NULL
      if (!is.null(extra_args$sensModel_direct)) {
        .sm    <- extra_args$sensModel_direct
        .inner <- tryCatch(ui$foceiModel$inner, error = function(e) NULL)
        if (!is.null(.inner)) {
          .scf <- file.path(rxode2::rxTempDir(),
                            paste0("adm-sens-", digest::digest(.inner), ".qs2"))
          if (file.exists(.scf)) {
            all_args_par$sens_cache_file <- .scf
            all_args_par$sens_cols       <- .sm$sens_cols
            all_args_par$sens_rename     <- .sm$rename_map
          }
        }
      }
      all_args_par$studies <- lapply(all_args_par$studies, function(s) {
        for (.ev_field in c("ev", "ev_full"))
          if (!is.null(s[[.ev_field]])) s[[.ev_field]] <- as.data.frame(s[[.ev_field]])
        s
      })
    }

    # Batch loop -- print after each batch of effective_workers restarts.
    batches <- split(seq_len(n_r), ceiling(seq_len(n_r) / effective_workers))
    results <- vector("list", n_r)

    if (use_fork) {
      for (.batch in batches) {
        .br <- furrr::future_map(
          .batch, function(r) {
            args <- all_args_par
            args$cores <- cores_vec[[(r - 1L) %% length(cores_vec) + 1L]]
            do.call(worker_fn, c(list(restart_id = r, p_init = inits[[r]]), args))
          },
          .options = furrr::furrr_options(seed = NULL, globals = FALSE)
        )
        for (i in seq_along(.batch)) {
          results[[.batch[[i]]]] <- .br[[i]]
          message(.restart_msg(.batch[[i]], .br[[i]]))
        }
      }
    } else if (pkg_locked) {
      .par_lambda <- function(r) {
        wfn  <- get(.worker_fn_name, envir = asNamespace("admixr2"), inherits = FALSE)
        args <- all_args_par
        args$cores <- cores_vec[[(r - 1L) %% length(cores_vec) + 1L]]
        do.call(wfn, c(list(restart_id = r, p_init = inits[[r]]), args))
      }
      .par_lambda_env <- new.env(parent = baseenv())
      .par_lambda_env$.worker_fn_name <- .worker_fn_name
      .par_lambda_env$inits            <- inits
      .par_lambda_env$all_args_par     <- all_args_par
      .par_lambda_env$cores_vec        <- cores_vec
      environment(.par_lambda) <- .par_lambda_env
      .furrr_opts <- furrr::furrr_options(seed = NULL, packages = "admixr2", globals = FALSE)
      for (.batch in batches) {
        .br <- furrr::future_map(.batch, .par_lambda, .options = .furrr_opts)
        for (i in seq_along(.batch)) {
          results[[.batch[[i]]]] <- .br[[i]]
          message(.restart_msg(.batch[[i]], .br[[i]]))
        }
      }
    } else {
      .furrr_opts <- furrr::furrr_options(
        seed = NULL,
        globals = c(.fn_list, list(inits = inits, all_args_par = all_args_par,
                                   cores_vec = cores_vec, worker_fn = worker_fn))
      )
      for (.batch in batches) {
        .br <- furrr::future_map(
          .batch, function(r) {
            args <- all_args_par
            args$cores <- cores_vec[[(r - 1L) %% length(cores_vec) + 1L]]
            do.call(worker_fn, c(list(restart_id = r, p_init = inits[[r]]), args))
          },
          .options = .furrr_opts
        )
        for (i in seq_along(.batch)) {
          results[[.batch[[i]]]] <- .br[[i]]
          message(.restart_msg(.batch[[i]], .br[[i]]))
        }
      }
    }
  } else {
    if (n_r > 1L) {
      if (!requireNamespace("furrr", quietly = TRUE) ||
          !requireNamespace("future", quietly = TRUE)) {
        message(sprintf("  Running %d restarts sequentially (install furrr+future for parallel)",
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
  for (nm in names(studies))
    studies[[nm]] <- .admNormaliseStudy(studies[[nm]], nm)

  pinfo      <- .admParseIniDf(.ui$iniDf, .ui)
  output_var <- .admOutputVar(.ui)

  want_grad    <- .ctl$grad != "none"
  want_sens    <- .ctl$grad == "sens"
  want_central <- .ctl$grad == "cfd"

  if (pinfo$n_eta > 0L && any(!pinfo$struct_has_eta)) {
    .unpaired <- names(pinfo$struct_has_eta)[!pinfo$struct_has_eta]
    if (want_sens && all(!pinfo$struct_has_eta)) {
      message(sprintf("admc: no mu-referenced struct thetas (%s); falling back to full forward FD.",
                      paste(.unpaired, collapse = ", ")))
      want_sens <- FALSE
    } else {
      message(sprintf(
        "admc: struct theta(s) without mu-referencing: %s. %s",
        paste(.unpaired, collapse = ", "),
        if (want_sens) "Sens model for paired thetas; FD for unpaired."
        else "FD gradient for these parameters."
      ))
    }
  }

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

  rxMod <- .admLoadModel(.ui)
  rxode2::rxLock(rxMod)
  on.exit({ rxode2::rxUnlock(rxMod); rxode2::rxSolveFree() }, add = TRUE)

  for (nm in names(studies)) {
    s  <- studies[[nm]]
    ev <- if (!is.null(s$ev)) s$ev else rxode2::et(amt = 100)
    studies[[nm]]$ev_full <- ev |> rxode2::et(s$times)
  }

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

  eval_grad_f <- if (want_grad) {
    function(p) .admGrad(p, pinfo, studies, z_list, rxMod, output_var,
                         params_list, cores, grad_h, sensModel,
                         use_central = want_central)
  } else NULL

  grad_label <- if (!want_grad) "none" else if (!is.null(sensModel))
    if (pinfo$has_kappa) "Sens+FD" else "Sens"
  else if (want_central) "central FD" else "forward FD"
  message("=== admixr2: Aggregate Data Modeling (MC) ===")
  message(sprintf("  Studies: %d | MC samples: %d | Params: %d | Cores: %d | Grad: %s | Restarts: %d",
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
    .adm_old_plan <- .admSetupParallelPlan(.ctl, .ctl$n_restarts)
    if (!is.null(.adm_old_plan)) on.exit(future::plan(.adm_old_plan), add = TRUE)
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
    admStopWorkers()
    .iter <- opt$n_iter
  }

  t_opt     <- (proc.time() - t0)["elapsed"]
  final     <- .admUnpack(opt$solution, pinfo)
  fullTheta <- .admFullTheta(final, pinfo)

  p_hat  <- setNames(opt$solution, names(ov$p0))
  t0_cov <- proc.time()
  .cov <- if (.ctl$covMethod == "r") {
    use_grad_cov <- want_grad
    use_cent_cov <- want_central
    np_cov <- length(pinfo$struct_names) + length(pinfo$sigma_names)
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
  .ret$control <- .admToFoceiControl(.ctl)
  .focei_model <- suppressMessages(tryCatch(.ui$foceiModel, error = function(e) NULL))
  if (!is.null(.focei_model)) .ret$model <- .focei_model

  .fit <- nlmixr2est::nlmixr2CreateOutputFromUi(
    .ui, data = admData(), control = .ret$control,
    table = .ret$table, env = .ret, est = "admc")

  .fit$env$method   <- "admc"
  .fit$env$studies  <- studies
  .fit$env$admExtra <- .ret$admExtra
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
