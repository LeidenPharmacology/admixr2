# -- FO (First-Order) aggregate data estimator ---------------------------------
# Approximates the aggregate NLL analytically:
#   mu_pred = f(theta, 0)        population prediction at eta = 0
#   J[t, j] = df/d(eta_j)|_0    Jacobian via sensitivity equations or FD
#   V_pred  = J Omega J^T + sigma_contribution
#
# One rxSolve call (via sens model) per study per NLL evaluation, vs n_sim for
# MC. Substantially faster per iteration; suitable for initial estimates or
# models where the FO approximation is adequate.

# -- Internal helpers ----------------------------------------------------------

# Get mu = f(theta, 0) and J = df/d(eta)|_0 for one study.
# Returns list(mu, J): mu is length-n_t, J is n_t x n_eta.
# Prefers sens model (one pass); falls back to rxMod + forward FD.
.adfoGetMuJ <- function(pars, pinfo, s, sensModel, rxMod, output_var, params_mat, cores) {
  n_t   <- length(s$times)
  n_eta <- pinfo$n_eta

  if (n_eta > 0L && !is.null(sensModel)) {
    eta0 <- matrix(0, 1L, n_eta, dimnames = list(NULL, pinfo$eta_col_names))
    res  <- .admSimulateSens(sensModel, pars$struct, pinfo$sigma_names, eta0, s, cores,
                             pinfo$nDisplayProgress)
    if (!is.null(res)) {
      mu <- as.numeric(res$cp_mat)
      J  <- do.call(cbind, lapply(res$dpred_list, as.numeric))
      return(list(mu = mu, J = J))
    }
  }

  # FD fallback (or n_eta == 0)
  if (n_eta > 0L) {
    eta0 <- matrix(0, 1L, n_eta, dimnames = list(NULL, pinfo$eta_col_names))
    mu   <- as.numeric(.admSimulate(rxMod, pars$struct, pinfo$sigma_names,
                                     eta0, s, output_var, params_mat, cores,
                                     pinfo$nDisplayProgress))
    eps  <- 1e-6
    J    <- matrix(0, n_t, n_eta)
    for (j in seq_len(n_eta)) {
      eta_p    <- eta0; eta_p[1L, j] <- eps
      mu_p     <- as.numeric(.admSimulate(rxMod, pars$struct, pinfo$sigma_names,
                                           eta_p, s, output_var, params_mat, cores,
                                           pinfo$nDisplayProgress))
      J[, j]   <- (mu_p - mu) / eps
    }
    list(mu = mu, J = J)
  } else {
    eta_e <- matrix(numeric(0), 1L, 0L)
    mu    <- as.numeric(.admSimulate(rxMod, pars$struct, pinfo$sigma_names,
                                      eta_e, s, output_var, params_mat, cores,
                                      pinfo$nDisplayProgress))
    list(mu = mu, J = matrix(0, n_t, 0))
  }
}

# Build V_pred and mu_sigma (lnorm-corrected mean) for one study.
# Returns list(V, mu_sigma, JL): JL = J %*% L cached for gradient reuse.
# V = tcrossprod(JL) + sigma contributions (additive, proportional, or lognormal).
.adfoVpred <- function(mu_pred, J, L, sigma_var, sigma_is_prop, sigma_is_lnorm,
                        n_t, n_eta) {
  mu_sigma <- mu_pred
  for (i in seq_along(sigma_var))
    if (sigma_is_lnorm[[i]]) mu_sigma <- mu_sigma * exp(sigma_var[[i]] / 2)

  JL <- if (n_eta > 0L) J %*% L else matrix(0, n_t, 0)
  V  <- if (n_eta > 0L) tcrossprod(JL) else matrix(0, n_t, n_t)
  d_V <- diag(V)
  for (i in seq_along(sigma_var)) {
    sv <- sigma_var[[i]]
    if (sigma_is_prop[[i]]) {
      d_V <- d_V + sv * mu_pred^2
    } else if (sigma_is_lnorm[[i]]) {
      d_V <- d_V + mu_sigma^2 * (exp(sv) - 1)
    } else {
      d_V <- d_V + sv
    }
  }
  diag(V) <- d_V
  list(V = V, mu_sigma = mu_sigma, JL = JL)
}

# -- NLL -----------------------------------------------------------------------

# FO population mean + Jacobian for a JOINT (same-subject) unit: mu = f(theta,0)
# stacked across outputs, J = d(stacked)/d(eta)|_0 (n_total x n_eta). Prefers the
# sensitivity model (analytical J, one solve per output at eta = 0) so the
# struct-theta FD of the joint NLL is not corrupted by an inner FD Jacobian;
# falls back to finite differences when no sens model is available.
.adfoGetMuJJoint <- function(pars, pinfo, unit, sensModel, rxMod, params_mat, cores) {
  n_eta   <- pinfo$n_eta
  n_total <- unit$n_total
  eta0 <- matrix(0, 1L, n_eta, dimnames = list(NULL, pinfo$eta_col_names))
  if (n_eta > 0L && !is.null(sensModel)) {
    res <- .admSimulateJointSens(sensModel, pars$struct, pinfo$sigma_names,
                                 eta0, unit, cores, pinfo$nDisplayProgress)
    if (!is.null(res))
      return(list(mu = as.numeric(res$cp_mat),
                  J  = do.call(cbind, lapply(res$dpred_list, as.numeric))))
  }
  mu <- as.numeric(.admSimulateJoint(rxMod, pars$struct, pinfo$sigma_names,
                                     eta0, unit, params_mat, cores,
                                     pinfo$nDisplayProgress))
  if (n_eta == 0L) return(list(mu = mu, J = matrix(0, n_total, 0)))
  eps <- 1e-6; J <- matrix(0, n_total, n_eta)
  for (j in seq_len(n_eta)) {
    etap <- eta0; etap[1L, j] <- eps
    mup  <- as.numeric(.admSimulateJoint(rxMod, pars$struct, pinfo$sigma_names,
                                         etap, unit, params_mat, cores,
                                         pinfo$nDisplayProgress))
    J[, j] <- (mup - mu) / eps
  }
  list(mu = mu, J = J)
}

#' @noRd
.adfoNLL <- function(p, pinfo, studies, sensModel, rxMod, output_var,
                      params_list, cores) {
  pars  <- tryCatch(.admUnpack(p, pinfo), error = function(e) NULL)
  if (is.null(pars)) return(Inf)
  total <- 0

  for (s_idx in seq_along(studies)) {
    s   <- studies[[s_idx]]

    # Joint (same-subject) unit: stacked FO mean/Jacobian -> full V = J Omega J'
    # (cross blocks from the shared J) + per-output residual; single MVN.
    if (isTRUE(s$is_joint)) {
      mj <- .adfoGetMuJJoint(pars, pinfo, s, sensModel, rxMod, params_list[[s_idx]], cores)
      V  <- if (pinfo$n_eta > 0L) tcrossprod(mj$J %*% pars$L) else
        matrix(0, s$n_total, s$n_total)
      jr <- .admJointResidual(mj$mu, V, s, pinfo, pars$sigma_var)
      nll_s <- nll_cov_cpp(s$E, s$V, jr$mu, jr$V, s$n)
      if (!is.finite(nll_s)) return(Inf)
      total <- total + nll_s
      next
    }

    ov  <- s$output %||% output_var
    n_t <- length(s$times)
    sel <- .admSigmaSel(pinfo, ov)

    mj  <- .adfoGetMuJ(pars, pinfo, s, sensModel, rxMod, ov,
                        params_list[[s_idx]], cores)
    vp  <- .adfoVpred(mj$mu, mj$J, pars$L, pars$sigma_var[sel],
                       pinfo$sigma_is_prop[sel], pinfo$sigma_is_lnorm[sel],
                       n_t, pinfo$n_eta)

    nll_s <- if (s$method == "var") {
      nll_var_cpp(s$E, s$v_diag, vp$mu_sigma, diag(vp$V), s$n)
    } else {
      nll_cov_cpp(s$E, s$V, vp$mu_sigma, vp$V, s$n)
    }
    if (!is.finite(nll_s)) return(Inf)
    total <- total + nll_s
  }
  total
}

# -- Gradient ------------------------------------------------------------------

# Gradient of FO NLL w.r.t. optimizer parameter vector p.
#
# Omega/sigma: analytical (chain rule through V_pred = J*Omega*J^T + sigma).
# Struct thetas: forward FD of full NLL -- required because J also depends on
#   theta (V-path: d(J*Omega*J^T)/d(theta) needs second-order sensitivities
#   unavailable from rxode2).
#
# (mu, J) are computed ONCE per study and shared between:
#   - the FD baseline NLL for struct theta FD
#   - the analytical omega/sigma gradient
# This avoids the redundant extra rxSolve that would occur if .adfoNLL() were
# called separately for the baseline and then .admSimulateSens() called again
# for the analytical gradient.
#' @noRd
.adfoGrad <- function(p, pinfo, studies, sensModel, rxMod, output_var,
                       params_list, cores, grad_h = 1e-4) {
  pars <- tryCatch(.admUnpack(p, pinfo), error = function(e) NULL)
  if (is.null(pars)) return(rep(NA_real_, length(p)))

  n_s   <- length(pinfo$struct_names)
  n_e   <- length(pinfo$sigma_names)
  n_eta <- pinfo$n_eta
  n_o   <- length(pinfo$omega_par)
  grad  <- numeric(length(p))
  names(grad) <- names(p)

  # --- Pass 1: compute (mu, J) and NLL per study at current p ---------------
  # Results cached for struct theta FD baseline and omega/sigma gradient.
  nll_0    <- 0
  muj_cache <- vector("list", length(studies))

  for (s_idx in seq_along(studies)) {
    s   <- studies[[s_idx]]

    # Joint (same-subject) unit: stacked FO mean/Jacobian; V = J Omega J' + per-
    # output residual. Analytical omega/sigma in Pass 3; struct thetas are FD of
    # the joint-aware .adfoNLL in Pass 2 (unchanged).
    if (isTRUE(s$is_joint)) {
      mj <- .adfoGetMuJJoint(pars, pinfo, s, sensModel, rxMod, params_list[[s_idx]], cores)
      JL <- if (n_eta > 0L) mj$J %*% pars$L else matrix(0, s$n_total, 0)
      Vs <- if (n_eta > 0L) tcrossprod(JL) else matrix(0, s$n_total, s$n_total)
      jr <- .admJointResidual(mj$mu, Vs, s, pinfo, pars$sigma_var)
      nll_s <- nll_cov_cpp(s$E, s$V, jr$mu, jr$V, s$n)
      if (!is.finite(nll_s)) { nll_0 <- Inf; break }
      nll_0 <- nll_0 + nll_s
      muj_cache[[s_idx]] <- list(is_joint = TRUE, mu = mj$mu, J = mj$J, JL = JL,
                                 V = jr$V, mu_sigma = jr$mu, n_t = s$n_total)
      next
    }

    ov  <- s$output %||% output_var
    sel <- .admSigmaSel(pinfo, ov)
    n_t <- length(s$times)

    mj  <- .adfoGetMuJ(pars, pinfo, s, sensModel, rxMod, ov,
                        params_list[[s_idx]], cores)
    vp  <- .adfoVpred(mj$mu, mj$J, pars$L, pars$sigma_var[sel],
                       pinfo$sigma_is_prop[sel], pinfo$sigma_is_lnorm[sel], n_t, n_eta)

    nll_s <- if (s$method == "var") {
      nll_var_cpp(s$E, s$v_diag, vp$mu_sigma, diag(vp$V), s$n)
    } else {
      nll_cov_cpp(s$E, s$V, vp$mu_sigma, vp$V, s$n)
    }
    if (!is.finite(nll_s)) { nll_0 <- Inf; break }
    nll_0 <- nll_0 + nll_s
    muj_cache[[s_idx]] <- list(mu = mj$mu, J = mj$J, JL = vp$JL, vp = vp,
                               n_t = n_t, sel = sel)
  }

  # --- Pass 2: struct theta forward FD (reuses cached nll_0 as baseline) -------
  if (n_s > 0L && is.finite(nll_0)) {
    for (k in seq_len(n_s)) {
      hk  <- pmax(abs(p[k]), 0.1) * grad_h
      p_p <- p; p_p[k] <- p[k] + hk
      nll_p <- .adfoNLL(p_p, pinfo, studies, sensModel, rxMod, output_var,
                         params_list, cores)
      grad[k] <- (nll_p - nll_0) / hk
    }
  }

  # --- Pass 3: analytical omega and sigma gradient (uses cached mu/J/vp) ----
  for (s_idx in seq_along(studies)) {
    mc <- muj_cache[[s_idx]]
    if (is.null(mc)) next
    s      <- studies[[s_idx]]

    # Joint (same-subject) analytical omega/sigma on the stacked V = J Omega J'.
    if (isTRUE(mc$is_joint)) {
      mu_pred <- mc$mu; J <- mc$J; JL <- mc$JL; V_pred <- mc$V
      mu_sigma <- mc$mu_sigma
      r    <- as.numeric(s$E) - mu_sigma
      invV <- tryCatch(chol2inv(chol(V_pred)),
                       error = function(e) tryCatch(solve(V_pred), error = function(e2) NULL))
      if (is.null(invV)) next
      dNLL_dV      <- s$n * (invV - invV %*% (s$V + tcrossprod(r)) %*% invV)
      dNLL_dV_diag <- diag(dNLL_dV)
      if (n_eta > 0L && n_o > 0L) {
        ML <- crossprod(J, dNLL_dV %*% JL)
        for (r_idx in seq_along(pinfo$omega_par)) {
          i <- pinfo$chol_i[r_idx]; jj <- pinfo$chol_j[r_idx]
          grad[n_s + n_e + r_idx] <- grad[n_s + n_e + r_idx] +
            if (pinfo$chol_diag[r_idx]) as.numeric(ML[i, i]) * pars$L[i, i]
            else 2 * as.numeric(ML[i, jj])
        }
      }
      dNLL_dmu_full <- drop(-2 * s$n * invV %*% r)
      for (blk in s$blocks) {
        rows <- blk$rows
        for (k in which(.admSigmaSel(pinfo, blk$output))) {
          k_sig <- n_s + k; sv <- pars$sigma_var[[k]]
          if (pinfo$sigma_is_prop[[k]])
            grad[k_sig] <- grad[k_sig] + sum(dNLL_dV_diag[rows] * mu_pred[rows]^2) * sv
          else if (pinfo$sigma_is_lnorm[[k]])
            grad[k_sig] <- grad[k_sig] + sv * (
              sum(dNLL_dV_diag[rows] * mu_sigma[rows]^2 * (2 * exp(sv) - 1)) +
              sum(dNLL_dmu_full[rows] * mu_sigma[rows]) / 2)
          else
            grad[k_sig] <- grad[k_sig] + sum(dNLL_dV_diag[rows]) * sv
        }
      }
      next
    }

    sel    <- mc$sel
    mu_pred <- mc$mu; J <- mc$J; vp <- mc$vp; n_t <- mc$n_t
    V_pred  <- vp$V;  mu_sigma <- vp$mu_sigma
    r       <- as.numeric(s$E) - mu_sigma
    is_var  <- identical(s$method, "var")

    if (is_var) {
      v_pred       <- diag(V_pred)
      dNLL_dv_pred <-      s$n * (1/v_pred - (s$v_diag + r^2) / v_pred^2)
      dNLL_dV_diag <- dNLL_dv_pred
    } else {
      invV <- tryCatch(
        chol2inv(chol(V_pred)),
        error = function(e) tryCatch(solve(V_pred), error = function(e2) NULL))
      if (is.null(invV)) next
      V_obs_hat    <- s$V + tcrossprod(r)
      dNLL_dV      <- s$n * (invV - invV %*% V_obs_hat %*% invV)
      dNLL_dV_diag <- diag(dNLL_dV)
    }

    # Omega gradient: d(-2LL)/d(L[i,j]) via ML = t(J) dNLL_dV JL
    # (JL cached from Pass 1 -- avoids materialising M = t(J) dNLL_dV J separately)
    # Diagonal p = log(Omega_ii) = 2*log(L_ii):
    #   d(NLL)/d(L_ii) = 2*(ML)[i,i]  (factor-of-2 from symmetric d(LL^T)/d(L_ii))
    #   d(L_ii)/dp     = L_ii/2        (chain rule through L_ii = exp(p/2))
    #   => d(NLL)/dp   = (ML)[i,i] * L_ii   (the two factors of 2 cancel)
    # Off-diagonal p = L_ij:    chain rule gives 2*(ML)[i,j]
    if (n_eta > 0L && n_o > 0L) {
      JL <- mc$JL
      ML <- if (is_var) crossprod(J, JL * dNLL_dv_pred) else crossprod(J, dNLL_dV %*% JL)
      d  <- pinfo$chol_diag
      for (r_idx in seq_along(pinfo$omega_par)) {
        i  <- pinfo$chol_i[r_idx]; jj <- pinfo$chol_j[r_idx]
        if (d[r_idx]) {
          grad[n_s + n_e + r_idx] <- grad[n_s + n_e + r_idx] +
            as.numeric(ML[i, i]) * pars$L[i, i]
        } else {
          grad[n_s + n_e + r_idx] <- grad[n_s + n_e + r_idx] +
            2 * as.numeric(ML[i, jj])
        }
      }
    }

    # Sigma gradient: d(-2LL)/d(p_sigma) = d(-2LL)/d(sv) * sv. Only this output's
    # residual-error parameters contribute (no-op for single-output fits).
    for (i in which(sel)) {
      k_sig <- n_s + i
      sv <- pars$sigma_var[[i]]
      if (pinfo$sigma_is_prop[[i]]) {
        grad[k_sig] <- grad[k_sig] + sum(dNLL_dV_diag * mu_pred^2) * sv
      } else if (pinfo$sigma_is_lnorm[[i]]) {
        dNLL_dmu     <- if (is_var) -2 * s$n * r / diag(V_pred) else
          drop(-2 * s$n * invV %*% r)
        grad[k_sig] <- grad[k_sig] + sv * (
          sum(dNLL_dV_diag * mu_sigma^2 * (2 * exp(sv) - 1)) +
          sum(dNLL_dmu * mu_sigma) / 2
        )
      } else {
        grad[k_sig] <- grad[k_sig] + sum(dNLL_dV_diag) * sv
      }
    }
  }

  grad
}

# Pure finite-difference gradient of FO NLL w.r.t. all optimizer parameters.
# Used for grad = "fd" (forward FD) and grad = "cfd" (central FD).
# No analytical chain rule -- every parameter differentiated by perturbing the NLL.
# The sens model is still used inside .adfoNLL() for efficient J computation.
#' @noRd
.adfoFDGrad <- function(p, pinfo, studies, sensModel, rxMod, output_var,
                         params_list, cores, grad_h = 1e-4, use_central = FALSE) {
  g <- numeric(length(p)); names(g) <- names(p)
  if (use_central) {
    for (k in seq_along(p)) {
      hk <- pmax(abs(p[k]), 0.1) * grad_h
      pp <- p; pp[k] <- p[k] + hk
      pm <- p; pm[k] <- p[k] - hk
      g[k] <- (.adfoNLL(pp, pinfo, studies, sensModel, rxMod, output_var, params_list, cores) -
               .adfoNLL(pm, pinfo, studies, sensModel, rxMod, output_var, params_list, cores)) / (2 * hk)
    }
  } else {
    nll0 <- .adfoNLL(p, pinfo, studies, sensModel, rxMod, output_var, params_list, cores)
    for (k in seq_along(p)) {
      hk <- pmax(abs(p[k]), 0.1) * grad_h
      ph <- p; ph[k] <- p[k] + hk
      g[k] <- (.adfoNLL(ph, pinfo, studies, sensModel, rxMod, output_var, params_list, cores) - nll0) / hk
    }
  }
  g
}

# -- Covariance ----------------------------------------------------------------

# Post-fit covariance via numerical Hessian.
# use_grad=TRUE: forward FD of gradient (np_cov+1 grad evals -- faster).
# use_grad=FALSE: full NLL-FD quadratic form (1+2*np_cov+4*n_off NLL evals).
# Hessian restricted to struct + sigma parameters (omega Cholesky excluded).
.adfoCalcCov <- function(p_hat, pinfo, studies, sensModel, rxMod, output_var,
                          params_list, cores,
                          use_grad = FALSE, grad_h = 1e-4, cov_h = 1e-3,
                          cov_h_outer = .Machine$double.eps^(1/5)) {
  n_s     <- length(pinfo$struct_names)
  n_e     <- length(pinfo$sigma_names)
  cov_idx <- seq_len(n_s + n_e)
  np_cov  <- length(cov_idx)
  nms_cov <- names(p_hat)[cov_idx]
  message("  Note: covMethod='r' computes covariance for structural and sigma ",
          "parameters only; omega (IIV) SEs are not computed (matching nlmixr2 ",
          "FOCEI behavior).")

  nll_fn <- function(p)
    suppressMessages(.adfoNLL(p, pinfo, studies, sensModel, rxMod, output_var,
                               params_list, cores))
  grad_fn <- function(p)
    suppressMessages(.adfoGrad(p, pinfo, studies, sensModel, rxMod, output_var,
                                params_list, cores, grad_h = cov_h))

  nll0 <- nll_fn(p_hat)
  if (!is.finite(nll0)) {
    warning("adfoCalcCov: NLL not finite at p_hat -- covariance not computed",
            call. = FALSE)
    return(NULL)
  }

  H <- matrix(0, np_cov, np_cov, dimnames = list(nms_cov, nms_cov))

  if (use_grad) {
    h_fwd <- pmax(abs(p_hat[cov_idx]), 0.1) * cov_h_outer
    g0    <- grad_fn(p_hat)[cov_idx]
    for (jj in seq_len(np_cov)) {
      ph        <- p_hat; ph[cov_idx[jj]] <- ph[cov_idx[jj]] + h_fwd[jj]
      gj        <- grad_fn(ph)[cov_idx]
      H[, jj]   <- if (anyNA(gj)) 0 else (gj - g0) / h_fwd[jj]
    }
    H <- (H + t(H)) / 2
  } else {
    h_gill <- pmax(abs(p_hat[cov_idx]), 0.1) * cov_h_outer
    for (k in seq_len(np_cov)) {
      ki  <- cov_idx[k]; hk  <- h_gill[k]
      p_p <- p_hat; p_p[ki] <- p_p[ki] + hk
      p_m <- p_hat; p_m[ki] <- p_m[ki] - hk
      H[k, k] <- (nll_fn(p_p) - 2 * nll0 + nll_fn(p_m)) / hk^2
    }
    for (i in seq_len(np_cov - 1L)) {
      for (j in seq(i + 1L, np_cov)) {
        ii <- cov_idx[i]; ji <- cov_idx[j]
        hi <- h_gill[i];  hj <- h_gill[j]
        p_pp <- p_hat; p_pp[ii] <- p_pp[ii] + hi; p_pp[ji] <- p_pp[ji] + hj
        p_pm <- p_hat; p_pm[ii] <- p_pm[ii] + hi; p_pm[ji] <- p_pm[ji] - hj
        p_mp <- p_hat; p_mp[ii] <- p_mp[ii] - hi; p_mp[ji] <- p_mp[ji] + hj
        p_mm <- p_hat; p_mm[ii] <- p_mm[ii] - hi; p_mm[ji] <- p_mm[ji] - hj
        H[i, j] <- H[j, i] <-
          (nll_fn(p_pp) - nll_fn(p_pm) - nll_fn(p_mp) + nll_fn(p_mm)) / (4 * hi * hj)
      }
    }
  }

  if (!all(is.finite(H))) {
    warning("adfoCalcCov: Hessian has non-finite entries -- covariance not computed",
            call. = FALSE)
    return(NULL)
  }


  eig_dec <- tryCatch(eigen(H, symmetric = TRUE), error = function(e) NULL)
  H_eigs  <- if (!is.null(eig_dec)) eig_dec$values else rep(NA_real_, np_cov)

  if (!is.null(eig_dec) && min(H_eigs) < 0) {
    warning(sprintf(
      "adfoCalcCov: Hessian not positive definite (min eigenvalue %.3e). Covariance not computed. Try increasing cov_h_outer (currently %.3e), e.g. cov_h_outer = %.3e.",
      min(H_eigs), cov_h_outer, cov_h_outer * 4), call. = FALSE)
    return(NULL)
  }

  Hinv <- tryCatch(
    chol2inv(chol(H)),
    error = function(e) tryCatch(solve(H), error = function(e2) NULL))
  if (is.null(Hinv)) {
    warning("adfoCalcCov: Hessian inversion failed -- covariance not computed",
            call. = FALSE)
    return(NULL)
  }

  cov_full <- (2 * Hinv + t(2 * Hinv)) / 2
  dimnames(cov_full) <- list(nms_cov, nms_cov)
  cov_full[pinfo$struct_names, pinfo$struct_names, drop = FALSE]
}

# -- Restart worker ------------------------------------------------------------

# Self-contained FO optimization run (one restart); serializable for furrr workers.
# Signature mirrors .admRestartWorker: same base_args from .admRunRestarts().
# n_sim, sampling accepted for interface compatibility but not used.
# use_pure_fd: TRUE for grad="fd"/"cfd"; dispatches to .adfoFDGrad() instead of .adfoGrad().
.adfoRestartWorker <- function(restart_id, p_init, ui_lstExpr, pinfo,
                                ov_lower, ov_upper, scale_c = NULL, studies, n_sim,
                                seed, algorithm, ftol_rel, maxeval,
                                use_grad, grad_h, grad_bounds,
                                output_var = "cp",
                                sampling = "sobol",
                                use_central = FALSE,
                                use_pure_fd = FALSE,
                                print_progress = TRUE, print = 10L,
                                cores = NULL, no_lock = FALSE,
                                sens_cache_file = NULL, sens_cols = NULL,
                                sens_rename = NULL,
                                rxMod_direct = NULL, sensModel_direct = NULL) {
  library(admixr2)

  # Dev mode (PSOCK workers): patch installed namespace with dev functions from
  # .GlobalEnv (serialised there by furrr globals). tryCatch guards against
  # the installed package predating this function (run devtools::install() once).
  tryCatch(.admPatchDevNamespace(), error = function(e) NULL)

  m <- .admWorkerLoadModels(ui_lstExpr, rxMod_direct, cores,
                            sens_cache_file, sens_cols, sens_rename, sensModel_direct)

  params_list <- .admMakeParamsList(1L, pinfo, length(studies))

  set.seed(seed + restart_id)

  nll_fn <- function(p)
    .adfoNLL(p, pinfo, studies, m$sensModel, m$rxMod, output_var, params_list, m$cores_w)

  grad_fn <- if (use_pure_fd) {
    function(p) .adfoFDGrad(p, pinfo, studies, m$sensModel, m$rxMod, output_var,
                            params_list, m$cores_w, grad_h, use_central)
  } else {
    function(p) .adfoGrad(p, pinfo, studies, m$sensModel, m$rxMod, output_var,
                          params_list, m$cores_w, grad_h)
  }

  # adfo loads its own model in-process and does not lock (single-nloptr path).
  .admScaledOptimize(restart_id, p_init, ov_lower, ov_upper, scale_c,
                     use_grad, grad_bounds, algorithm, ftol_rel, maxeval,
                     nll_fn, grad_fn, pinfo, print_progress, print,
                     lock_rxMod = NULL)
}

# -- Control object ------------------------------------------------------------

#' Control settings for the FO (First-Order) estimator
#'
#' Creates a control object for `nlmixr2(est = "adfo")`.  The FO estimator
#' linearises model predictions at \eqn{\eta = 0}: it is faster than the MC
#' estimator but less accurate for models with large IIV or strongly
#' non-linear individual predictions.
#'
#' @param studies Named list of study specifications (same format as
#'   [admControl()]: `E`, `V`, `n`, `times`, `ev`, optional `method`; or an
#'   `observations` list for multi-compartment fits -- see [admControl()]).
#' @param grad Gradient mode.  `"none"` (default) uses derivative-free BOBYQA;
#'   `"analytical"` uses the closed-form FO gradient (requires sensitivity
#'   equations); `"fd"` uses forward finite differences of the full NLL;
#'   `"cfd"` uses central finite differences for struct theta gradient
#'   (more accurate than `"fd"`, roughly twice as many NLL evaluations per step).
#' @param algorithm nloptr algorithm, or `NULL` (default) to pick the default
#'   that matches `grad`: `"NLOPT_LD_LBFGS"` with a gradient, `"NLOPT_LN_BOBYQA"`
#'   when `grad = "none"`. Any algorithm reported by
#'   [nloptr::nloptr.print.options()] is accepted. An explicit algorithm is
#'   reconciled with `grad`: when `grad = "none"` a gradient-based algorithm
#'   (`NLOPT_LD_*` / `NLOPT_GD_*`) falls back to `"NLOPT_LN_BOBYQA"`; when a
#'   gradient is requested a derivative-free algorithm (`NLOPT_LN_*` /
#'   `NLOPT_GN_*`) turns the gradient off. Both emit a message.
#' @param maxeval Maximum function evaluations (default 500).
#' @param ftol_rel Relative tolerance (default `sqrt(.Machine$double.eps)`).
#' @param print Print-frequency for live progress (0 = silent).
#' @param seed Random seed (used for restarts).
#' @param cores OpenMP threads for `rxSolve()` (default 1).
#' @param nDisplayProgress Passed to `rxSolve()`: show the solver's text
#'   progress bar only once a single solve exceeds this many subjects. The
#'   default (`.Machine$integer.max`) keeps it off for clean script/vignette
#'   output; lower it (e.g. `1000L`) to see progress during long fits.
#' @param grad_h Finite-difference step for unpaired struct theta gradient and
#'   FD Jacobian.
#' @param grad_bounds Box-constraint half-width when using gradients.
#' @param cov_h_outer Outer step scale for NLL-FD Hessian.
#' @param covMethod `"r"` computes covariance via numerical Hessian for
#'   structural and residual-error parameters only (omega/IIV SEs are not
#'   computed, consistent with nlmixr2 FOCEI); `"none"` skips it.
#' @param n_restarts Number of optimizer restarts (1 = no multi-start).
#' @param restart_sd Standard deviation for random perturbations of initial
#'   struct thetas at each restart (> 1).
#' @param workers Number of parallel PSOCK/fork workers for multi-restart
#'   (default 1 = sequential).
#' @param cov_h Inner FD step for the gradient-based Hessian (only used when
#'   `covMethod = "r"` and `grad != "none"`). Default 1e-3.
#' @param rxControl `rxode2::rxControl()` object. Created automatically when `NULL`.
#' @param calcTables,compress,ci,sigdig,sigdigTable,optExpression,sumProd,literalFix
#'   Passed to `nlmixr2est::foceiControl()` for the table/output machinery.
#' @param addProp How combined additive+proportional error is parameterised in
#'   the nlmixr2 output tables: `"combined2"` (default, variance form) or
#'   `"combined1"` (SD form). Has no effect on admixr2's own estimation.
#' @param returnAdmr If `TRUE`, return a plain list instead of the full
#'   nlmixr2 fit object.
#' @param ... Unused arguments (trigger an error).
#'
#' @return An `adfoControl` object (a named list).
#'
#' @seealso [admControl()], [adirmcControl()]
#'
#' @examples
#' # Inspect defaults
#' ctl <- adfoControl()
#' ctl$grad
#' ctl$maxeval
#'
#' # Analytical gradient, more evaluations
#' ctl2 <- adfoControl(grad = "analytical", maxeval = 1000L)
#'
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
adfoControl <- function(
    studies    = list(),
    grad        = c("none", "analytical", "fd", "cfd"),
    algorithm  = NULL,
    maxeval    = 500L,
    ftol_rel   = .Machine$double.eps^(1/2),
    print      = 10L,
    seed       = 12345L,
    cores      = 1L,
    nDisplayProgress = .Machine$integer.max,
    grad_h      = 1e-4,
    grad_bounds = 5,
    cov_h       = 1e-3,
    cov_h_outer = .Machine$double.eps^(1/5),
    covMethod   = c("r", "none"),
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
    stop("adfoControl: unused argument(s): ",
         paste(paste0("'", names(.xtra), "'"), collapse = ", "), call. = FALSE)

  addProp  <- match.arg(addProp)
  grad     <- match.arg(grad)
  covMethod <- match.arg(covMethod)

  checkmate::assertList(studies)
  checkmate::assertIntegerish(maxeval,    lower = 1L, len = 1)
  checkmate::assertNumeric(ftol_rel,      lower = 0,  len = 1)
  checkmate::assertIntegerish(print,      lower = 0L, len = 1)
  checkmate::assertIntegerish(seed,                   len = 1)
  checkmate::assertIntegerish(cores,      lower = 1L, len = 1)
  checkmate::assertIntegerish(nDisplayProgress, lower = 1L, len = 1,
                              .var.name = "nDisplayProgress")
  checkmate::assertNumeric(grad_h,        lower = 0,  len = 1)
  checkmate::assertNumeric(grad_bounds,   lower = 0,  len = 1)
  checkmate::assertNumeric(cov_h,         lower = 0,  len = 1)
  checkmate::assertNumeric(cov_h_outer,   lower = 0,  len = 1)
  checkmate::assertIntegerish(n_restarts, lower = 1L, len = 1)
  checkmate::assertNumeric(restart_sd,    lower = 0,  len = 1)
  checkmate::assertIntegerish(workers,    lower = 1L, len = 1)
  checkmate::assertNumeric(ci, lower = 0, upper = 1,  len = 1)
  checkmate::assertIntegerish(sigdig,     lower = 1L, len = 1)
  checkmate::assertLogical(returnAdmr,                len = 1)

  .algo     <- .admResolveAlgorithm(algorithm, grad,
                                    .var.name = "adfoControl: algorithm")
  algorithm <- .algo$algorithm
  grad      <- .algo$grad

  if (is.null(rxControl))   rxControl   <- rxode2::rxControl(sigdig = sigdig)
  if (is.null(sigdigTable)) sigdigTable <- max(round(sigdig), 3L)

  .ret <- list(
    studies       = studies,
    n_sim         = 1L,        # interface compatibility with .admRunRestarts()
    sampling      = "sobol",   # idem
    grad          = grad,
    algorithm     = algorithm,
    maxeval       = as.integer(maxeval),
    ftol_rel      = ftol_rel,
    print         = as.integer(print),
    seed          = as.integer(seed),
    cores         = as.integer(cores),
    nDisplayProgress = as.integer(nDisplayProgress),
    grad_h        = grad_h,
    grad_bounds   = grad_bounds,
    cov_h         = cov_h,
    cov_h_outer   = cov_h_outer,
    covMethod     = covMethod,
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
  class(.ret) <- "adfoControl"
  .ret
}

# -- nlmixr2 S3 hooks ----------------------------------------------------------

#' @noRd
getValidNlmixrCtl.adfo <- function(control) {
  if (inherits(control, "adfoControl")) return(control)
  .ctl <- control[[1]]
  if (inherits(.ctl, "adfoControl")) return(.ctl)
  if (is.list(.ctl) && "studies" %in% names(.ctl))
    return(do.call(adfoControl, .ctl[intersect(names(.ctl), names(formals(adfoControl)))]))
  if (is.list(control) && length(names(control)) > 0L)
    return(do.call(adfoControl, control[intersect(names(control), names(formals(adfoControl)))]))
  adfoControl()
}

#' @noRd
nmObjHandleControlObject.adfoControl <- function(control, env) {
  assign("adfoControl", control, envir = env)
}

#' @noRd
nmObjGetControl.adfo <- function(x, ...) {
  .env <- x[[1]]
  for (.nm in c("adfoControl", "control")) {
    if (exists(.nm, .env)) {
      .ctl <- get(.nm, .env)
      if (inherits(.ctl, "adfoControl")) return(.ctl)
    }
  }
  stop("cannot find adfo control object", call. = FALSE)
}

# -- Main estimation entry point -----------------------------------------------

#' Fit an aggregate data model via First-Order (FO) approximation
#'
#' Called automatically by `nlmixr2(model, admData(), est = "adfo",
#' control = adfoControl(...))`.  Not typically called directly.
#'
#' @param env nlmixr2 environment containing `ui` and `control`.
#' @param ... Unused.
#'
#' @return An `admFit` nlmixr2 fit object.
#'
#' @method nlmixr2Est adfo
#' @importFrom nlmixr2est nlmixr2Est
#' @export
nlmixr2Est.adfo <- function(env, ...) {
  .ui  <- env$ui
  .ctl <- env$control

  if (!inherits(.ctl, "adfoControl")) .ctl <- getValidNlmixrCtl.adfo(.ctl)
  if (!inherits(.ctl, "adfoControl"))
    stop("Could not recover adfoControl", call. = FALSE)
  assign("control", .ctl, envir = .ui)

  studies <- .ctl$studies
  if (length(studies) == 0L)
    stop("adfoControl(studies=...) required", call. = FALSE)
  if (is.null(names(studies)))
    names(studies) <- paste0("study", seq_along(studies))

  pinfo      <- .admParseIniDf(.ui$iniDf, .ui)
  pinfo$nDisplayProgress <- .ctl$nDisplayProgress %||% pinfo$nDisplayProgress
  output_var <- .admOutputVar(.ui)

  for (nm in names(studies))
    studies[[nm]] <- .admNormaliseStudy(studies[[nm]], nm, output_var)
  studies    <- .admFlattenStudies(studies)
  multi_out  <- length(.admOutputVars(.ui)) > 1L
  any_joint  <- any(vapply(studies, function(u) isTRUE(u$is_joint), logical(1)))
  studies    <- .admBuildEvFull(studies, tag_cmt = multi_out)

  want_grad   <- .ctl$grad != "none"
  want_sens   <- .ctl$grad == "analytical"
  use_central <- .ctl$grad == "cfd"
  use_pure_fd <- .ctl$grad %in% c("fd", "cfd")
  # Joint (same-subject) fits keep the analytical FO gradient: .adfoGrad's joint
  # branch applies the omega/sigma chain rule to the stacked V = J Omega J' + res
  # (struct thetas via FD of the joint-aware .adfoNLL, as in the single-output
  # path). grad = "fd"/"cfd" still use the pure FD gradient.
  # Multi-compartment fits use the same per-unit gradient as single-output: the
  # analytical omega/sigma chain rule and struct-theta FD apply per observed
  # output (each is an independent block with its own residual error).

  if (pinfo$n_eta > 0L && !is.null(pinfo$struct_has_eta) && any(!pinfo$struct_has_eta)) {
    .unpaired <- names(pinfo$struct_has_eta)[!pinfo$struct_has_eta]
    message(sprintf("adfo: struct theta(s) without mu-referencing: %s. FD for these parameters.",
                    paste(.unpaired, collapse = ", ")))
  }

  # ORDERING INVARIANT: .admLoadSensModel() must run before .admLoadModel().
  sensModel <- if (want_sens) {
    sm <- tryCatch(.admLoadSensModel(.ui), error = function(e) NULL)
    if (is.null(sm))
      warning("adfoControl(grad='analytical'): sensitivity model unavailable -- falling back to forward FD")
    sm
  } else NULL

  rxMod <- .admLoadModel(.ui)
  rxode2::rxLock(rxMod)
  on.exit({ rxode2::rxUnlock(rxMod); rxode2::rxSolveFree() }, add = TRUE)

  params_list <- .admMakeParamsList(1L, pinfo, length(studies))

  ov    <- .admBuildOptVec(pinfo)
  .iter <- 0L
  cores <- .ctl$cores

  .nll_trace <- numeric(0)
  .par_trace <- NULL
  .best_nll  <- Inf

  eval_f <- function(p) {
    .iter <<- .iter + 1L
    val <- .adfoNLL(p, pinfo, studies, sensModel, rxMod, output_var, params_list, cores)
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

  eval_grad_f <- if (!want_grad) {
    NULL
  } else if (use_pure_fd) {
    function(p) .adfoFDGrad(p, pinfo, studies, sensModel, rxMod, output_var,
                              params_list, cores, .ctl$grad_h, use_central)
  } else {
    function(p) .adfoGrad(p, pinfo, studies, sensModel, rxMod, output_var,
                           params_list, cores, .ctl$grad_h)
  }

  grad_label <- if (!want_grad) "none" else if (!is.null(sensModel)) "Analytical"
    else if (use_central) "CFD" else "FD"
  message("=== admixr2: Aggregate Data Modeling (FO) ===")
  message(sprintf("  Obs units: %d | Params: %d | Cores: %d | Grad: %s | Restarts: %d",
                  length(studies), length(ov$p0), cores, grad_label, .ctl$n_restarts))
  t0 <- proc.time()

  lb <- if (want_grad) pmax(ov$lower, ov$p0 - .ctl$grad_bounds) else ov$lower
  ub <- if (want_grad) pmin(ov$upper, ov$p0 + .ctl$grad_bounds) else ov$upper

  sc           <- ov$scale_c
  p0_sc        <- ov$p0 / sc
  lb_sc        <- lb  / sc
  ub_sc        <- ub  / sc
  eval_f_sc    <- function(p_s) eval_f(p_s * sc)
  eval_grad_sc <- if (!is.null(eval_grad_f)) {
    function(p_s) eval_grad_f(p_s * sc) * sc
  } else NULL

  if (.ctl$n_restarts == 1L) {
    message(.admProgressHeader(pinfo))
    opt_raw <- nlmixr2est::nlmixrWithTiming("adfo", {
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
      worker_fn  = .adfoRestartWorker,
      p0         = ov$p0, ov = ov, pinfo = pinfo,
      .ctl       = .ctl, ui = .ui, studies = studies,
      extra_args = list(algorithm       = .ctl$algorithm,
                        ftol_rel        = .ctl$ftol_rel,
                        maxeval         = .ctl$maxeval,
                        use_grad        = want_grad,
                        use_central     = use_central,
                        use_pure_fd     = use_pure_fd,
                        grad_h          = .ctl$grad_h,
                        grad_bounds     = .ctl$grad_bounds,
                        output_var      = output_var,
                        print_progress  = TRUE,
                        print           = .ctl$print,
                        cores           = .ctl$cores,
                        rxMod_direct    = rxMod,
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
    np_cov     <- length(pinfo$struct_names) + length(pinfo$sigma_names)
    # Joint fits use the NLL-FD Hessian (the grad-FD Hessian calls the
    # single-output analytical .adfoGrad, which does not handle joint units).
    use_grad_cov <- want_grad && !any_joint
    n_evals <- if (use_grad_cov) {
      np_cov + 1L
    } else {
      n_off <- np_cov * (np_cov - 1L) / 2L
      1L + 2L * np_cov + 4L * n_off
    }
    evals_label <- if (use_grad_cov) "gradient evaluations" else "NLL evaluations"
    hess_label  <- if (!use_grad_cov) "" else if (want_sens) ", Analytical-Hessian" else ", FD-Hessian"
    message(sprintf("  Computing covariance (R method%s, %d %s)", hess_label, n_evals, evals_label))
    tryCatch(
      .adfoCalcCov(p_hat, pinfo, studies, sensModel, rxMod, output_var,
                   params_list, cores, use_grad = use_grad_cov,
                   grad_h = .ctl$grad_h, cov_h = .ctl$cov_h,
                   cov_h_outer = .ctl$cov_h_outer),
      error = function(e) { warning("adfoCalcCov failed: ", conditionMessage(e)); NULL })
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
  .ret$est        <- "adfo"
  .ret$ofvType    <- "adfo"
  .ret$adjObf     <- FALSE
  .ret$covMethod  <- if (!is.null(.cov)) "r" else ""
  .ret$cov        <- .cov
  .ret$message    <- opt$message
  .ret$extra      <- ""
  .ret$origData   <- studies

  .ret$admExtra <- list(struct         = final$struct,
                        sigma_var      = final$sigma_var,
                        sigma_is_prop  = pinfo$sigma_is_prop,
                        sigma_is_lnorm = pinfo$sigma_is_lnorm,
                        omega          = final$omega,
                        L              = final$L,
                        eta_col_names  = pinfo$eta_col_names,
                        par_names      = names(ov$p0),
                        npar           = length(ov$p0),
                        nloptr         = opt,
                        nll_trace      = .nll_trace,
                        par_trace      = .par_trace,
                        all_traces     = opt$all_traces,
                        n_iter         = .iter,
                        time           = t_elapsed,
                        t_opt          = t_opt,
                        t_cov          = t_cov,
                        studies        = studies,
                        n_sim          = 5000L,  # for diagnostic plots; FO itself is n_sim=1
                        sampling       = "sobol")

  nlmixr2est::.nlmixr2FitUpdateParams(.ret)
  nmObjHandleControlObject.adfoControl(.ctl, .ret)
  if (exists("control", .ui)) rm(list = "control", envir = .ui)
  .ret$control <- .admToFoceiControl(.ctl)
  .focei_model <- suppressMessages(tryCatch(.ui$foceiModel, error = function(e) NULL))
  if (!is.null(.focei_model)) .ret$model <- .focei_model

  .fit <- nlmixr2est::nlmixr2CreateOutputFromUi(
    .ui, data = if (multi_out) admData(.admOutputVars(.ui)) else admData(),
    control = .ret$control,
    table = .ret$table, env = .ret, est = "adfo")

  .fit$env$method   <- "adfo"
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
  row.names(.stats$objDf) <- "adfo"
  .fit$env$logLik    <- .stats$ll
  .fit$env$nobs      <- .stats$nobs
  .fit$env$objDf     <- .stats$objDf
  .fit$env$OBJF      <- .stats$objDf$OBJF
  .fit$env$AIC       <- .stats$objDf$AIC
  .fit$env$BIC       <- .stats$objDf$BIC
  .fit$env$objective <- opt$objective
  .fit$env$time <- data.frame(
    optimize   = t_opt,
    covariance = t_cov,
    other      = 0,
    elapsed    = t_elapsed,
    row.names  = NULL
  )

  .fit
}
