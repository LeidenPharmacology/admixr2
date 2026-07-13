# -- adgh: aggregate Gauss-Hermite quadrature estimator -------------------------
# Computes population moments E[f] and Cov[f] for eta ~ N(0, Omega) by
# deterministic Gauss-Hermite quadrature over the random-effects distribution,
# then plugs them into the same aggregate MVN -2LL as adfo/admc.
#
# Structurally this is admc with a small fixed deterministic node grid in place
# of n_sim random draws, and quadrature weights in place of uniform ones.
# The objective is noise-free -> clean gradient/Hessian, fast reproducible opt.
#
# The measure is prior N(0, Omega) (prior-predictive population moments, not
# data-conditional posterior), so plain (non-adaptive) GH is exactly right.

# -- Node grid -----------------------------------------------------------------

# Probabilists' GH nodes/weights for E_{N(0,1)}[g] = sum_i w_i g(x_i).
# Golub-Welsch via symmetric tridiagonal eigendecomposition. No external deps.
# sum(w) = 1, sum(w * x^2) = 1.
.adghNodes1 <- function(m) {
  if (m < 1L) stop("n_nodes must be >= 1")
  if (m == 1L) return(list(x = 0, w = 1))
  i <- seq_len(m - 1L)
  J <- matrix(0, m, m)
  J[cbind(i, i + 1L)] <- sqrt(i)
  J[cbind(i + 1L, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = (e$vectors[1L, ])^2)
}

# Tensor-product GH grid for n_eta dimensions.
# Returns X (n_node x n_eta standard-normal nodes) and W (length n_node weights).
.adghNodeGrid <- function(m, n_eta) {
  if (n_eta == 0L) return(list(X = matrix(0, 1L, 0L), W = 1))
  g <- .adghNodes1(m)
  X <- as.matrix(expand.grid(rep(list(g$x), n_eta)))
  W <- as.numeric(apply(expand.grid(rep(list(g$w), n_eta)), 1L, prod))
  dimnames(X) <- NULL
  list(X = X, W = W)
}

# -- Moments -------------------------------------------------------------------

# Population moments (E, V) for one study via GH quadrature.
# One batched .admSimulate over the node grid; weighted ML mean/cov; residual
# error added to the diagonal exactly as adfo/admc.
# Quadrature grid (eta nodes + weights) for the current Omega. Shared by the
# single and batched moment paths.
.adghGrid <- function(pars, pinfo, grid) {
  if (pinfo$n_eta > 0L) {
    eta <- grid$X %*% t(pars$L)
    colnames(eta) <- pinfo$eta_col_names
    list(eta = eta, W = grid$W)
  } else {
    list(eta = matrix(0, 1L, 0L), W = 1)
  }
}

# Weighted moments + residual error from an already-solved quadrature matrix.
# Split out of .adghMoments so the solve and the assembly can be driven
# independently: the assembly depends on sigma, but the SOLVE does not (sigma is
# zeroed into it and re-added analytically here), so a set of configurations
# that share a solve can each be assembled cheaply.
.adghMomentsFromCp <- function(cp, W, pars, pinfo, out_var) {
  mu  <- as.numeric(crossprod(W, cp))
  cpc <- sweep(cp, 2L, mu)
  V   <- crossprod(cpc, W * cpc)

  # Restrict residual error to this output's sigma(s) (no-op single-output).
  mu_sigma <- mu
  for (i in which(.admSigmaSel(pinfo, out_var))) {
    sv <- pars$sigma_var[[i]]
    if (isTRUE(pinfo$sigma_is_lnorm[[i]])) {
      mu_sigma <- mu_sigma * exp(sv / 2)
      diag(V)  <- diag(V) + mu_sigma^2 * (exp(sv) - 1)
    } else if (isTRUE(pinfo$sigma_is_prop[[i]])) {
      diag(V) <- diag(V) + sv * mu^2
    } else {
      diag(V) <- diag(V) + sv
    }
  }
  list(E = mu_sigma, V = V)
}

.adghMoments <- function(pars, pinfo, study, rxMod, out_var, grid, cores) {
  g  <- .adghGrid(pars, pinfo, grid)
  pm <- .admMakeParamsList(nrow(g$eta), pinfo, 1L)[[1L]]
  cp <- .admSimulate(rxMod, pars$struct, pinfo$sigma_names, g$eta, study,
                     out_var, pm, cores, pinfo$nDisplayProgress)
  .adghMomentsFromCp(cp, g$W, pars, pinfo, out_var)
}

# Moments for a SET of structural-theta configurations in ONE rxSolve.
# The node grid and Omega are identical across configurations -- only the
# structural thetas move -- so the n_cfg quadrature solves stack into one call
# of n_cfg * n_node subjects instead of n_cfg calls of n_node.
.adghMomentsBatch <- function(struct_mat, pars, pinfo, study, rxMod, out_var, grid, cores) {
  g     <- .adghGrid(pars, pinfo, grid)
  Q     <- nrow(g$eta)
  n_cfg <- nrow(struct_mat)

  sm_big  <- struct_mat[rep(seq_len(n_cfg), each = Q), , drop = FALSE]
  eta_big <- g$eta[rep(seq_len(Q), times = n_cfg), , drop = FALSE]
  colnames(eta_big) <- colnames(g$eta)

  pm     <- .admMakeParamsList(n_cfg * Q, pinfo, 1L)[[1L]]
  cp_all <- .admSimulateRows(rxMod, sm_big, pinfo$sigma_names, eta_big, study,
                             out_var, pm, cores, pinfo$nDisplayProgress)

  lapply(seq_len(n_cfg), function(k)
    .adghMomentsFromCp(cp_all[(k - 1L) * Q + seq_len(Q), , drop = FALSE],
                       g$W, pars, pinfo, out_var))
}

# GH-quadrature joint moments for a same-subject unit: one shared-eta node grid
# gives every output; stacked weighted mean/cov (full, cross blocks included) +
# per-output residual. Returns list(E = mu_sigma, V).
.adghMomentsJoint <- function(pars, pinfo, unit, rxMod, grid, cores) {
  n_eta <- pinfo$n_eta
  if (n_eta > 0L) {
    eta <- grid$X %*% t(pars$L); colnames(eta) <- pinfo$eta_col_names; W <- grid$W
  } else { eta <- matrix(0, 1L, 0L); W <- 1 }
  pm <- .admMakeParamsList(nrow(eta), pinfo, 1L)[[1L]]
  cp <- .admSimulateJoint(rxMod, pars$struct, pinfo$sigma_names, eta, unit, pm, cores,
                          pinfo$nDisplayProgress)
  mu  <- as.numeric(crossprod(W, cp))
  cpc <- sweep(cp, 2L, mu)
  V   <- crossprod(cpc, W * cpc)
  jr  <- .admJointResidual(mu, V, unit, pinfo, pars$sigma_var)
  list(E = jr$mu, V = jr$V)
}

# -- NLL -----------------------------------------------------------------------

#' @noRd
.adghNLL <- function(p, pinfo, studies, rxMod, out_var, grid, cores) {
  pars <- tryCatch(.admUnpack(p, pinfo), error = function(e) NULL)
  if (is.null(pars)) return(Inf)
  total <- 0
  for (s in studies) {
    if (isTRUE(s$is_joint)) {
      m   <- .adghMomentsJoint(pars, pinfo, s, rxMod, grid, cores)
      nll <- nll_cov_cpp(s$E, s$V, m$E, m$V, s$n)
    } else {
      m <- .adghMoments(pars, pinfo, s, rxMod, s$output %||% out_var, grid, cores)
      nll <- if (identical(s$method, "var"))
        nll_var_cpp(s$E, s$v_diag, m$E, diag(m$V), s$n)
      else
        nll_cov_cpp(s$E, s$V, m$E, m$V, s$n)
    }
    if (!is.finite(nll)) return(Inf)
    total <- total + nll
  }
  total
}

# -- Analytic gradient ---------------------------------------------------------

# Analytic gradient of the GH NLL w.r.t. optimizer vector p.
# One batched sensitivity solve per study over the node grid; closed-form
# contractions for struct thetas, omega Cholesky, sigma (add/prop/lnorm).
# Unpaired struct thetas: forward FD of .adghNLL (like admc).
#
# Var-method studies use a diagonal derivative path; cov-method uses full B.
# For lnorm sigma: Jl scaled by exp(sv/2) for mean path; analytical sigma grad.
.adghGrad <- function(p, pinfo, studies, sensModel, rxMod, out_var, grid, cores,
                       grad_h = 1e-4) {
  pars  <- .admUnpack(p, pinfo)
  L     <- pars$L
  n_eta <- pinfo$n_eta
  n_s   <- length(pinfo$struct_names)
  n_e   <- length(pinfo$sigma_names)
  X     <- grid$X
  W     <- grid$W
  grad  <- numeric(length(p)); names(grad) <- names(p)

  # Which struct thetas are unpaired (no mu-referencing eta)?
  # struct_has_eta is struct-indexed (length n_s); struct_eta_idx is eta-indexed
  # (length n_eta) so is.na() on it never flags unpaired struct thetas.
  unpaired_k <- if (!is.null(pinfo$struct_has_eta))
    which(!pinfo$struct_has_eta) else integer(0)

  for (s in studies) {
    eta <- X %*% t(L)
    colnames(eta) <- pinfo$eta_col_names

    # --- Joint (same-subject) analytical quadrature gradient -----------------
    # Stacked weighted moments over all outputs (shared-eta node grid); paired
    # struct + omega + sigma analytical on the joint covariance, per output rows.
    # Unpaired struct thetas use the FD block after this loop (.adghNLL is joint).
    if (isTRUE(s$is_joint)) {
      js <- .admSimulateJointSens(sensModel, pars$struct, pinfo$sigma_names, eta, s, cores,
                                  pinfo$nDisplayProgress)
      if (is.null(js)) next
      f  <- js$cp_mat; Jl <- js$dpred_list
      mu  <- as.numeric(crossprod(W, f))
      cpc <- sweep(f, 2L, mu)
      # per-row lnorm mean scaling; residual-adjusted mean/cov via shared helper
      ls_vec <- rep(1, s$n_total)
      for (blk in s$blocks) for (k in which(.admSigmaSel(pinfo, blk$output)))
        if (isTRUE(pinfo$sigma_is_lnorm[[k]]))
          ls_vec[blk$rows] <- ls_vec[blk$rows] * exp(pars$sigma_var[[k]] / 2)
      jr <- .admJointResidual(mu, crossprod(cpc, W * cpc), s, pinfo, pars$sigma_var)
      mu_sigma <- jr$mu; V <- jr$V
      r  <- as.numeric(s$E) - mu_sigma
      G  <- tryCatch(chol2inv(chol(V)),
                     error = function(e) tryCatch(solve(V), error = function(e2) NULL))
      if (is.null(G)) next
      B    <- s$n * (G - G %*% (s$V + tcrossprod(r)) %*% G)
      dNLL_dmu_sig <- as.numeric(-2 * s$n * (G %*% r))
      Bdiag <- diag(B); Bt <- cpc %*% B
      contrib_j <- function(gmat) {              # gmat = d(mu_sigma)/dpsi rows
        dmu <- as.numeric(crossprod(W, gmat))
        sum(dNLL_dmu_sig * dmu) + 2 * sum(W * rowSums(gmat * Bt))
      }
      sig_V_extra <- function(dmu_raw) {          # dmu_raw = d(mu)/dpsi (pre-lnorm)
        acc <- 0
        for (blk in s$blocks) { rows <- blk$rows
          for (k in which(.admSigmaSel(pinfo, blk$output))) {
            sv <- pars$sigma_var[[k]]
            if (isTRUE(pinfo$sigma_is_prop[[k]]))
              acc <- acc + sum(Bdiag[rows] * 2 * sv * mu[rows] * dmu_raw[rows])
            else if (isTRUE(pinfo$sigma_is_lnorm[[k]]))
              acc <- acc + sum(Bdiag[rows] * 2 * (exp(sv) - 1) * ls_vec[rows] *
                                 mu_sigma[rows] * dmu_raw[rows])
          }
        }
        acc
      }
      # paired struct thetas
      for (k in seq_len(n_s)) {
        ei <- which(pinfo$struct_eta_idx == k); if (length(ei) == 0L) next; ei <- ei[[1L]]
        gmat    <- sweep(Jl[[ei]], 2L, ls_vec, "*")
        dmu_raw <- as.numeric(crossprod(W, Jl[[ei]]))
        grad[k] <- grad[k] + contrib_j(gmat) + sig_V_extra(dmu_raw)
      }
      # omega Cholesky
      if (n_eta > 0L) for (rr in seq_along(pinfo$omega_par)) {
        i <- pinfo$chol_i[rr]; j <- pinfo$chol_j[rr]
        base    <- Jl[[i]] * X[, j]
        gmat    <- sweep(base, 2L, ls_vec, "*")
        dmu_raw <- as.numeric(crossprod(W, base))
        dL  <- contrib_j(gmat) + sig_V_extra(dmu_raw)
        pos <- n_s + n_e + rr
        grad[pos] <- grad[pos] + if (pinfo$chol_diag[rr]) dL * L[i, i] / 2 else dL
      }
      # sigma, per output rows
      for (blk in s$blocks) { rows <- blk$rows
        for (k in which(.admSigmaSel(pinfo, blk$output))) {
          sv <- pars$sigma_var[[k]]; k_sig <- n_s + k
          if (isTRUE(pinfo$sigma_is_lnorm[[k]]))
            grad[k_sig] <- grad[k_sig] + (sum(dNLL_dmu_sig[rows] * mu_sigma[rows] / 2) +
              sum(Bdiag[rows] * mu_sigma[rows]^2 * (2 * exp(sv) - 1))) * sv
          else if (isTRUE(pinfo$sigma_is_prop[[k]]))
            grad[k_sig] <- grad[k_sig] + sum(Bdiag[rows] * mu[rows]^2) * sv
          else
            grad[k_sig] <- grad[k_sig] + sum(Bdiag[rows]) * sv
        }
      }
      next
    }

    ov <- s$output %||% out_var
    sk <- which(.admSigmaSel(pinfo, ov))   # this output's sigma indices

    res <- .admSimulateSens(sensModel, pars$struct, pinfo$sigma_names, eta, s, cores,
                            pinfo$nDisplayProgress)
    # .admSimulateSens returns NULL when the solve fails. The joint branch above
    # guards this; the non-joint branch did not, and errored mid-gradient on
    # `res$cp_mat` instead of degrading like every other path.
    if (is.null(res)) next
    f   <- res$cp_mat     # Q x n_t
    Jl  <- res$dpred_list # list n_eta of Q x n_t

    mu  <- as.numeric(crossprod(W, f))
    cpc <- sweep(f, 2L, mu)
    V   <- crossprod(cpc, W * cpc)

    # Sigma contributions to V (and lnorm scaling of mean) -- this output only
    mu_sigma <- mu
    sv_vec   <- pars$sigma_var
    for (i in sk) {
      sv <- sv_vec[[i]]
      if (isTRUE(pinfo$sigma_is_lnorm[[i]])) {
        mu_sigma <- mu_sigma * exp(sv / 2)
        diag(V)  <- diag(V) + mu_sigma^2 * (exp(sv) - 1)
      } else if (isTRUE(pinfo$sigma_is_prop[[i]])) {
        diag(V) <- diag(V) + sv * mu^2
      } else {
        diag(V) <- diag(V) + sv
      }
    }

    r <- as.numeric(s$E) - mu_sigma

    is_var <- identical(s$method, "var")

    if (is_var) {
      # ------ Var method: diagonal derivative path ----------------------------
      V_diag        <- diag(V)
      dNLL_dmu_sig  <- as.numeric(-2 * s$n * r / V_diag)  # d(NLL)/d(mu_sigma)
      dNLL_dV_diag  <- s$n * (1/V_diag - (s$v_diag + r^2) / V_diag^2)

      contrib <- function(gmat_mu) {
        # gmat_mu: Q x n_t, derivative of mu (not mu_sigma) w.r.t. psi
        dmu     <- as.numeric(crossprod(W, gmat_mu))
        dV_diag <- 2 * colSums(W * cpc * gmat_mu)
        sum(dNLL_dmu_sig * dmu) + sum(dNLL_dV_diag * dV_diag)
      }

      contrib_sigma <- function(i_sig) {
        sv       <- sv_vec[[i_sig]]
        is_lnorm <- isTRUE(pinfo$sigma_is_lnorm[[i_sig]])
        is_prop  <- isTRUE(pinfo$sigma_is_prop[[i_sig]])
        if (is_lnorm) {
          # d(NLL)/d(sv) via mu_sigma and V_diag
          d_mu_sig <- sum(dNLL_dmu_sig * mu_sigma / 2)
          d_V_diag <- sum(dNLL_dV_diag * mu_sigma^2 * (2 * exp(sv) - 1))
          (d_mu_sig + d_V_diag) * sv
        } else if (is_prop) {
          sum(dNLL_dV_diag * mu^2) * sv
        } else {
          sum(dNLL_dV_diag) * sv
        }
      }

    } else {
      # ------ Cov method: full-matrix derivative path -------------------------
      G    <- tryCatch(chol2inv(chol(V)),
                       error = function(e) tryCatch(solve(V), error = function(e2) NULL))
      if (is.null(G)) next
      Vhat      <- s$V + tcrossprod(r)
      B         <- s$n * (G - G %*% Vhat %*% G)
      dNLL_dmu_sig <- as.numeric(-2 * s$n * (G %*% r))  # d(NLL)/d(mu_sigma)
      Bdiag     <- diag(B)
      Bt        <- cpc %*% B  # Q x n_t; row q = cpc[q,] B

      contrib <- function(gmat_mu) {
        # gmat_mu: derivative of mu (not mu_sigma) w.r.t. psi
        dmu      <- as.numeric(crossprod(W, gmat_mu))
        term_mu  <- sum(dNLL_dmu_sig * dmu)
        term_cov <- 2 * sum(W * rowSums(gmat_mu * Bt))
        term_mu + term_cov
      }

      contrib_sigma <- function(i_sig) {
        sv       <- sv_vec[[i_sig]]
        is_lnorm <- isTRUE(pinfo$sigma_is_lnorm[[i_sig]])
        is_prop  <- isTRUE(pinfo$sigma_is_prop[[i_sig]])
        if (is_lnorm) {
          d_mu_sig <- sum(dNLL_dmu_sig * mu_sigma / 2)
          d_V_diag <- sum(Bdiag * mu_sigma^2 * (2 * exp(sv) - 1))
          (d_mu_sig + d_V_diag) * sv
        } else if (is_prop) {
          sum(Bdiag * mu^2) * sv
        } else {
          sum(Bdiag) * sv
        }
      }
    }

    # For lnorm: sens Jacobians give d(mu)/d(eta); need d(mu_sigma)/d(eta).
    # Scale factor: exp(sv/2) for each lnorm sigma (multiplicative on mu).
    lnorm_scale <- 1
    for (i in sk)
      if (isTRUE(pinfo$sigma_is_lnorm[[i]])) lnorm_scale <- lnorm_scale * exp(sv_vec[[i]] / 2)

    # Precompute the d(NLL)/d(V_diag) vector used for sigma V-correction terms.
    # For prop:  extra_t = 2*sv*mu_t*d(mu_t)/dpsi;  for lnorm: 2*(exp(sv)-1)*scale*mu_sigma_t*d(mu_t)/dpsi
    # These couple through mu so must be evaluated per-parameter inside the loops.
    Bvec <- if (!is_var) Bdiag else dNLL_dV_diag  # length n_t

    .sigma_V_extra <- function(dmu_raw) {
      if (length(sk) == 0L) return(0)
      acc <- 0
      for (i_sig in sk) {
        sv <- sv_vec[[i_sig]]
        if (isTRUE(pinfo$sigma_is_prop[[i_sig]]))
          acc <- acc + sum(Bvec * 2 * sv * mu * dmu_raw)
        else if (isTRUE(pinfo$sigma_is_lnorm[[i_sig]]))
          acc <- acc + sum(Bvec * 2 * (exp(sv) - 1) * lnorm_scale * mu_sigma * dmu_raw)
      }
      acc
    }

    # Struct thetas (paired with etas; unpaired handled by FD below).
    # struct_eta_idx is eta-indexed (value = struct paired with each eta), so the
    # eta for struct k is which(struct_eta_idx == k).
    for (k in seq_len(n_s)) {
      if (!is.null(pinfo$struct_has_eta) && !pinfo$struct_has_eta[k]) next  # unpaired
      ei <- which(pinfo$struct_eta_idx == k)[1L]  # struct k -> its eta dim
      if (is.na(ei)) next  # nocov -- defensive; ei always found when struct_has_eta[k]
      gmat    <- if (lnorm_scale != 1) Jl[[ei]] * lnorm_scale else Jl[[ei]]
      dmu_raw <- as.numeric(crossprod(W, Jl[[ei]]))  # d(mu_t)/d(psi) before lnorm scaling
      grad[k] <- grad[k] + contrib(gmat) + .sigma_V_extra(dmu_raw)
    }

    # Omega Cholesky L: d(eta[q,])/d(L_ij) = x[q,j] * e_i (unit vector eta dim i)
    # So d(f[q,])/d(L_ij) = Jl[[i]][q,] * X[q,j]
    # Chain: L_ii stored as log(Omega_ii) -> d(L_ii)/dp = L_ii/2.
    if (n_eta > 0L) for (rr in seq_along(pinfo$omega_par)) {
      i <- pinfo$chol_i[rr]; j <- pinfo$chol_j[rr]
      gmat    <- if (lnorm_scale != 1) Jl[[i]] * X[, j] * lnorm_scale else Jl[[i]] * X[, j]
      dmu_raw <- as.numeric(crossprod(W, Jl[[i]] * X[, j]))
      dL      <- contrib(gmat) + .sigma_V_extra(dmu_raw)
      pos <- n_s + n_e + rr
      grad[pos] <- grad[pos] + if (pinfo$chol_diag[rr]) dL * L[i, i] / 2 else dL
    }

    # Sigma: d(NLL)/d(sv) * d(sv)/d(p) where sv = exp(p) so d(sv)/d(p) = sv.
    # Only this output's residual-error parameters (no-op for single-output).
    for (i in sk) {
      sig_g <- contrib_sigma(i)
      grad[n_s + i] <- grad[n_s + i] + sig_g
    }
  }

  # Unpaired struct thetas: forward FD of .adghNLL.
  #
  # The baseline and every perturbed configuration differ only in their
  # structural thetas -- same node grid, same Omega, same sigma -- so they all
  # share one solve per study (.adghMomentsBatch), with the baseline carried as
  # configuration 1. That also removes the separate baseline .adghNLL(p) pass,
  # which re-solved every study to recompute a value this loop already needs.
  #
  # Joint units keep the per-configuration path (their solve is per output block).
  if (length(unpaired_k) > 0L) {
    n_u <- length(unpaired_k)
    hs  <- pmax(abs(p[unpaired_k]), 0.1) * grad_h
    p_pert <- lapply(seq_len(n_u), function(i) {
      pp <- p; pp[unpaired_k[i]] <- p[unpaired_k[i]] + hs[i]; pp
    })

    if (any(vapply(studies, function(u) isTRUE(u$is_joint), logical(1)))) {
      nll0 <- .adghNLL(p, pinfo, studies, rxMod, out_var, grid, cores)
      for (i in seq_len(n_u))
        grad[unpaired_k[i]] <-
          (.adghNLL(p_pert[[i]], pinfo, studies, rxMod, out_var, grid, cores) - nll0) / hs[i]
    } else {
      # configuration 1 = baseline, 1 + i = unpaired theta i perturbed
      struct_mat <- do.call(rbind, c(
        list(pars$struct),
        lapply(p_pert, function(pp) .admUnpack(pp, pinfo)$struct)))
      colnames(struct_mat) <- names(pars$struct)

      nll_cfg <- numeric(n_u + 1L)
      for (s in studies) {
        ovs <- s$output %||% out_var
        ms  <- .adghMomentsBatch(struct_mat, pars, pinfo, s, rxMod, ovs, grid, cores)
        for (cfg in seq_len(n_u + 1L)) {
          if (!is.finite(nll_cfg[cfg])) next
          m     <- ms[[cfg]]
          nll_c <- if (identical(s$method, "var"))
            nll_var_cpp(s$E, s$v_diag, m$E, diag(m$V), s$n)
          else
            nll_cov_cpp(s$E, s$V, m$E, m$V, s$n)
          nll_cfg[cfg] <- if (is.finite(nll_c)) nll_cfg[cfg] + nll_c else Inf
        }
      }
      grad[unpaired_k] <- (nll_cfg[-1L] - nll_cfg[1L]) / hs
    }
  }

  grad
}

# -- FD gradient ---------------------------------------------------------------

.adghFDGrad <- function(p, pinfo, studies, rxMod, out_var, grid, cores,
                          grad_h = 1e-4, use_central = FALSE) {
  g <- numeric(length(p)); names(g) <- names(p)
  if (use_central) {
    for (k in seq_along(p)) {
      hk <- pmax(abs(p[k]), 0.1) * grad_h
      pp <- p; pp[k] <- p[k] + hk
      pm <- p; pm[k] <- p[k] - hk
      g[k] <- (.adghNLL(pp, pinfo, studies, rxMod, out_var, grid, cores) -
               .adghNLL(pm, pinfo, studies, rxMod, out_var, grid, cores)) / (2 * hk)
    }
  } else {
    nll0 <- .adghNLL(p, pinfo, studies, rxMod, out_var, grid, cores)
    for (k in seq_along(p)) {
      hk <- pmax(abs(p[k]), 0.1) * grad_h
      ph <- p; ph[k] <- p[k] + hk
      g[k] <- (.adghNLL(ph, pinfo, studies, rxMod, out_var, grid, cores) - nll0) / hk
    }
  }
  g
}

# -- Covariance ----------------------------------------------------------------

# Post-fit covariance via numerical Hessian (struct + sigma params only).
# Noise-free GH surface -> use tighter eps^(1/4) default step vs admc's eps^(1/5).
# use_grad=TRUE: forward FD of gradient (np+1 grad evals).
# use_grad=FALSE: full NLL-FD quadratic form (1+2*np+4*n_off NLL evals).
.adghCalcCov <- function(p_hat, pinfo, studies, sensModel, rxMod, out_var,
                           grid, cores,
                           use_grad = TRUE, grad_h = 1e-3,
                           cov_h_outer = .Machine$double.eps^(1/4)) {
  n_s     <- length(pinfo$struct_names)
  n_e     <- length(pinfo$sigma_names)
  cov_idx <- seq_len(n_s + n_e)
  np_cov  <- length(cov_idx)
  nms_cov <- names(p_hat)[cov_idx]
  message("  Note: covMethod='r' computes covariance for structural and sigma ",
          "parameters only; omega (IIV) SEs are not computed (matching nlmixr2 ",
          "FOCEI behavior).")

  nll_fn  <- function(p)
    suppressMessages(.adghNLL(p, pinfo, studies, rxMod, out_var, grid, cores))
  grad_fn <- function(p)
    suppressMessages(.adghGrad(p, pinfo, studies, sensModel, rxMod, out_var,
                                grid, cores, grad_h = grad_h))

  nll0 <- nll_fn(p_hat)
  if (!is.finite(nll0)) {
    warning("adghCalcCov: NLL not finite at p_hat -- covariance not computed",
            call. = FALSE)
    return(NULL)
  }

  H <- matrix(0, np_cov, np_cov, dimnames = list(nms_cov, nms_cov))

  if (use_grad) {
    h_fwd <- pmax(abs(p_hat[cov_idx]), 0.1) * cov_h_outer
    g0    <- grad_fn(p_hat)[cov_idx]
    for (jj in seq_len(np_cov)) {
      ph       <- p_hat; ph[cov_idx[jj]] <- ph[cov_idx[jj]] + h_fwd[jj]
      gj       <- grad_fn(ph)[cov_idx]
      H[, jj]  <- if (anyNA(gj)) 0 else (gj - g0) / h_fwd[jj]
    }
    H <- (H + t(H)) / 2
  } else {
    h_gill <- pmax(abs(p_hat[cov_idx]), 0.1) * cov_h_outer
    for (k in seq_len(np_cov)) {
      ki <- cov_idx[k]; hk <- h_gill[k]
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
    warning("adghCalcCov: Hessian has non-finite entries -- covariance not computed",
            call. = FALSE)
    return(NULL)
  }

  eig_dec <- tryCatch(eigen(H, symmetric = TRUE), error = function(e) NULL)
  H_eigs  <- if (!is.null(eig_dec)) eig_dec$values else rep(NA_real_, np_cov)

  if (!is.null(eig_dec) && min(H_eigs) < 0) {
    warning(sprintf(
      "adghCalcCov: Hessian not positive definite (min eigenvalue %.3e). Covariance not computed. Try increasing cov_h_outer (currently %.3e), e.g. cov_h_outer = %.3e.",
      min(H_eigs), cov_h_outer, cov_h_outer * 4), call. = FALSE)
    return(NULL)
  }

  Hinv <- tryCatch(
    chol2inv(chol(H)),
    error = function(e) tryCatch(solve(H), error = function(e2) NULL))
  if (is.null(Hinv)) {
    warning("adghCalcCov: Hessian inversion failed -- covariance not computed",
            call. = FALSE)
    return(NULL)
  }

  cov_full <- (2 * Hinv + t(2 * Hinv)) / 2
  dimnames(cov_full) <- list(nms_cov, nms_cov)
  cov_full[pinfo$struct_names, pinfo$struct_names, drop = FALSE]
}

# -- Restart worker ------------------------------------------------------------

# Self-contained GH optimization run (one restart); serializable to a worker.
# Signature mirrors .adfoRestartWorker: same base_args from .admRunRestarts().
# n_sim, sampling accepted for interface compatibility but not used.
.adghRestartWorker <- function(restart_id, p_init, ui_lstExpr, pinfo,
                                ov_lower, ov_upper, scale_c = NULL, studies, n_sim,
                                seed, n_nodes, algorithm, ftol_rel, maxeval,
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
  tryCatch(.admPatchDevNamespace(), error = function(e) NULL)

  m <- .admWorkerLoadModels(ui_lstExpr, rxMod_direct, cores,
                            sens_cache_file, sens_cols, sens_rename, sensModel_direct)

  grid <- .adghNodeGrid(n_nodes, pinfo$n_eta)
  set.seed(seed + restart_id)

  nll_fn <- function(p)
    .adghNLL(p, pinfo, studies, m$rxMod, output_var, grid, m$cores_w)

  grad_fn <- if (use_pure_fd) {
    function(p) .adghFDGrad(p, pinfo, studies, m$rxMod, output_var, grid, m$cores_w,
                            grad_h, use_central)
  } else {
    function(p) .adghGrad(p, pinfo, studies, m$sensModel, m$rxMod, output_var,
                          grid, m$cores_w, grad_h)
  }

  # adgh loads its own model in-process and does not lock (single-nloptr path).
  .admScaledOptimize(restart_id, p_init, ov_lower, ov_upper, scale_c,
                     use_grad, grad_bounds, algorithm, ftol_rel, maxeval,
                     nll_fn, grad_fn, pinfo, print_progress, print,
                     lock_rxMod = NULL)
}

# -- Control object ------------------------------------------------------------

#' Control settings for the Gauss-Hermite (GH) quadrature estimator
#'
#' Creates a control object for `nlmixr2(est = "adgh")`. The GH estimator
#' integrates model predictions against the random-effects prior
#' \eqn{\eta \sim N(0, \Omega)} using a deterministic tensor-product
#' Gauss-Hermite quadrature grid. It is unbiased at any IIV magnitude (unlike
#' FO), noise-free (unlike MC), and much faster than MC for models with up to
#' ~4 etas.
#'
#' @param studies Named list of study specifications (same format as
#'   [admControl()]: `E`, `V`, `n`, `times`, `ev`, optional `method`; or an
#'   `observations` list for multi-compartment fits -- see [admControl()]).
#' @param n_nodes Number of quadrature nodes per eta dimension (default 5).
#'   Total nodes = `n_nodes^n_eta`. `n_nodes = 5` achieves near-exact covariance
#'   moments for IIV SD up to ~0.5; `n_nodes = 7` extends coverage to SD ~0.7.
#'   For models with >= 5 etas the node count grows steeply; consider reducing
#'   `n_nodes` or using a different estimator.
#' @param grad Gradient mode. `"analytical"` (default) uses closed-form
#'   contractions through the sensitivity equations -- cheapest and exact.
#'   `"fd"` uses forward finite differences; `"cfd"` uses central FD.
#'   `"none"` uses derivative-free BOBYQA.
#' @param algorithm nloptr algorithm. Automatically coerced to
#'   `"NLOPT_LD_LBFGS"` when `grad != "none"`.
#' @param maxeval Maximum function evaluations (default 500).
#' @param ftol_rel Relative tolerance (default `sqrt(.Machine$double.eps)`).
#' @param print Print-frequency for live progress (0 = silent).
#' @param seed Random seed (used for restarts).
#' @param cores OpenMP threads for `rxSolve()`. Defaults to
#'   `rxode2::rxCores()`. When `workers > 1` it is a *total* budget, split
#'   across the workers.
#' @param nDisplayProgress Passed to `rxSolve()`: show the solver's text
#'   progress bar only once a single solve exceeds this many subjects. The
#'   default (`.Machine$integer.max`) keeps it off for clean script/vignette
#'   output; lower it (e.g. `1000L`) to see progress during long fits.
#' @param grad_h Finite-difference step for unpaired struct theta gradient and
#'   FD Jacobian fallback.
#' @param grad_bounds Box-constraint half-width when using gradients.
#' @param cov_h Inner FD step for the gradient-based Hessian (only used when
#'   `covMethod = "r"` and `grad != "none"`).
#' @param cov_h_outer Outer step scale for numerical Hessian. Default
#'   `eps^(1/4)` (tighter than admc's `eps^(1/5)` because the GH surface is
#'   noise-free).
#' @param covMethod `"r"` computes covariance via numerical Hessian for
#'   structural and residual-error parameters only; `"none"` skips it.
#' @param n_restarts Number of optimizer restarts (1 = no multi-start).
#' @param restart_sd SD of random perturbations of initial struct thetas at
#'   each restart.
#' @param workers Number of parallel workers (mirai daemons) for multi-restart
#'   (default 1 = sequential). Requires the `mirai` package.
#' @param rxControl `rxode2::rxControl()` object. Created automatically when `NULL`.
#' @param calcTables,compress,ci,sigdig,sigdigTable,optExpression,sumProd,literalFix
#'   Passed to `nlmixr2est::foceiControl()` for the table/output machinery.
#' @param addProp How combined additive+proportional error is parameterised in
#'   the nlmixr2 output tables: `"combined2"` (default) or `"combined1"`.
#' @param returnAdmr If `TRUE`, return a plain list instead of the full
#'   nlmixr2 fit object.
#' @param ... Unused arguments (trigger an error).
#'
#' @return An `adghControl` object (a named list).
#'
#' @seealso [admControl()], [adfoControl()], [adirmcControl()]
#'
#' @examples
#' ctl <- adghControl()
#' ctl$n_nodes
#' ctl$grad
#'
#' # More nodes for large IIV, analytical gradient
#' ctl2 <- adghControl(n_nodes = 7L, grad = "analytical", maxeval = 300L)
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
#'   pk_model, admData(), est = "adgh",
#'   control = adghControl(
#'     studies = list(study1 = list(E = E, V = V, n = length(ids),
#'                                  times = times, ev = et(amt = 100)))
#'   )
#' )
#' }
#'
#' @export
adghControl <- function(
    studies     = list(),
    n_nodes     = 5L,
    grad        = c("analytical", "fd", "cfd", "none"),
    algorithm   = "NLOPT_LN_BOBYQA",
    maxeval     = 500L,
    ftol_rel    = .Machine$double.eps^(1/2),
    print       = 10L,
    seed        = 12345L,
    cores       = rxode2::rxCores(),
    nDisplayProgress = .Machine$integer.max,
    grad_h      = 1e-4,
    grad_bounds = 5,
    cov_h       = 1e-3,
    cov_h_outer = .Machine$double.eps^(1/4),
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
  if (length(.xtra) > 0L)
    stop("adghControl: unused argument(s): ",
         paste(paste0("'", names(.xtra), "'"), collapse = ", "), call. = FALSE)

  addProp   <- match.arg(addProp)
  grad      <- match.arg(grad)
  covMethod <- match.arg(covMethod)

  checkmate::assertList(studies)
  checkmate::assertIntegerish(n_nodes,     lower = 1L, len = 1)
  checkmate::assertString(algorithm)
  checkmate::assertIntegerish(maxeval,     lower = 1L, len = 1)
  checkmate::assertNumeric(ftol_rel,       lower = 0,  len = 1)
  checkmate::assertIntegerish(print,       lower = 0L, len = 1)
  checkmate::assertIntegerish(seed,                    len = 1)
  checkmate::assertIntegerish(cores,       lower = 1L, len = 1)
  checkmate::assertIntegerish(nDisplayProgress, lower = 1L, len = 1,
                              .var.name = "nDisplayProgress")
  checkmate::assertNumeric(grad_h,         lower = 0,  len = 1)
  checkmate::assertNumeric(grad_bounds,    lower = 0,  len = 1)
  checkmate::assertNumeric(cov_h,          lower = 0,  len = 1)
  checkmate::assertNumeric(cov_h_outer,    lower = 0,  len = 1)
  checkmate::assertIntegerish(n_restarts,  lower = 1L, len = 1)
  checkmate::assertNumeric(restart_sd,     lower = 0,  len = 1)
  checkmate::assertIntegerish(workers,     lower = 1L, len = 1)
  checkmate::assertNumeric(ci, lower = 0, upper = 1,   len = 1)
  checkmate::assertIntegerish(sigdig,      lower = 1L, len = 1)
  checkmate::assertLogical(returnAdmr,                 len = 1)

  if (grad != "none" && algorithm == "NLOPT_LN_BOBYQA")
    algorithm <- "NLOPT_LD_LBFGS"

  if (is.null(rxControl))   rxControl   <- rxode2::rxControl(sigdig = sigdig)
  if (is.null(sigdigTable)) sigdigTable <- max(round(sigdig), 3L)

  .ret <- list(
    studies       = studies,
    n_nodes       = as.integer(n_nodes),
    n_sim         = 1L,       # interface compat with .admRunRestarts()
    sampling      = "sobol",  # idem
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
  class(.ret) <- "adghControl"
  .ret
}

# -- nlmixr2 S3 hooks ----------------------------------------------------------

#' @noRd
getValidNlmixrCtl.adgh <- function(control) {
  if (inherits(control, "adghControl")) return(control)
  .ctl <- control[[1]]
  if (inherits(.ctl, "adghControl")) return(.ctl)
  if (is.list(.ctl) && "studies" %in% names(.ctl))
    return(do.call(adghControl, .ctl[intersect(names(.ctl), names(formals(adghControl)))]))
  if (is.list(control) && length(names(control)) > 0L)
    return(do.call(adghControl, control[intersect(names(control), names(formals(adghControl)))]))
  adghControl()
}

#' @noRd
nmObjHandleControlObject.adghControl <- function(control, env) {
  assign("adghControl", control, envir = env)
}

#' @noRd
nmObjGetControl.adgh <- function(x, ...) {
  .env <- x[[1]]
  for (.nm in c("adghControl", "control")) {
    if (exists(.nm, .env)) {
      .ctl <- get(.nm, .env)
      if (inherits(.ctl, "adghControl")) return(.ctl)
    }
  }
  stop("cannot find adgh control object", call. = FALSE)
}

# -- Main estimation entry point -----------------------------------------------

#' Fit an aggregate data model via Gauss-Hermite quadrature
#'
#' Called automatically by `nlmixr2(model, admData(), est = "adgh",
#' control = adghControl(...))`. Not typically called directly.
#'
#' @param env nlmixr2 environment containing `ui` and `control`.
#' @param ... Unused.
#'
#' @return An `admFit` nlmixr2 fit object.
#'
#' @method nlmixr2Est adgh
#' @importFrom nlmixr2est nlmixr2Est
#' @export
nlmixr2Est.adgh <- function(env, ...) {
  .ui  <- env$ui
  .ctl <- env$control

  if (!inherits(.ctl, "adghControl")) .ctl <- getValidNlmixrCtl.adgh(.ctl)
  if (!inherits(.ctl, "adghControl"))
    stop("Could not recover adghControl", call. = FALSE)
  assign("control", .ctl, envir = .ui)

  studies <- .ctl$studies
  if (length(studies) == 0L)
    stop("adghControl(studies=...) required", call. = FALSE)
  if (is.null(names(studies)))
    names(studies) <- paste0("study", seq_along(studies))

  pinfo      <- .admParseIniDf(.ui$iniDf, .ui)
  pinfo$nDisplayProgress <- .ctl$nDisplayProgress %||% pinfo$nDisplayProgress
  output_var <- .admOutputVar(.ui)
  n_nodes    <- .ctl$n_nodes

  for (nm in names(studies))
    studies[[nm]] <- .admNormaliseStudy(studies[[nm]], nm, output_var)
  studies    <- .admFlattenStudies(studies)
  multi_out  <- length(.admOutputVars(.ui)) > 1L
  any_joint  <- any(vapply(studies, function(u) isTRUE(u$is_joint), logical(1)))
  studies    <- .admBuildEvFull(studies, tag_cmt = multi_out)

  want_grad    <- .ctl$grad != "none"
  want_sens    <- .ctl$grad == "analytical"
  use_central  <- .ctl$grad == "cfd"
  use_pure_fd  <- .ctl$grad %in% c("fd", "cfd")
  # Joint (same-subject) fits keep the analytical quadrature gradient: .adghGrad's
  # joint branch computes the stacked-MVN gradient from shared-eta per-output
  # sensitivities (grad = "analytical"). grad = "fd"/"cfd" use the FD gradient.

  if (pinfo$n_eta > 0L) {
    n_total <- n_nodes^pinfo$n_eta
    if (n_total > 5000L)
      message(sprintf(
        "adgh: n_nodes=%d x n_eta=%d = %d nodes. Consider reducing n_nodes or using est='admc'.",
        n_nodes, pinfo$n_eta, n_total))
  }

  if (!is.null(pinfo$struct_has_eta) && any(!pinfo$struct_has_eta)) {
    .unpaired <- names(pinfo$struct_has_eta)[!pinfo$struct_has_eta]
    message(sprintf("adgh: struct theta(s) without mu-referencing: %s. FD for these parameters.",
                    paste(.unpaired, collapse = ", ")))
  }

  # ORDERING INVARIANT: .admLoadSensModel() before .admLoadModel().
  sensModel <- if (want_sens) {
    sm <- tryCatch(.admLoadSensModel(.ui), error = function(e) NULL)
    if (is.null(sm)) {
      warning("adghControl(grad='analytical'): sensitivity model unavailable -- falling back to FD")
      want_sens   <- FALSE
      use_pure_fd <- TRUE
    }
    sm
  } else NULL

  rxMod <- .admLoadModel(.ui)
  rxode2::rxLock(rxMod)
  on.exit({ rxode2::rxUnlock(rxMod); rxode2::rxSolveFree() }, add = TRUE)

  # Node grid: fixed in standard-normal space; L applied per-eval in .adghMoments.
  grid  <- .adghNodeGrid(n_nodes, pinfo$n_eta)

  ov    <- .admBuildOptVec(pinfo)
  cores <- .ctl$cores
  .iter <- 0L

  .nll_trace <- numeric(0)
  .par_trace <- NULL
  .best_nll  <- Inf

  eval_f <- function(p) {
    .iter <<- .iter + 1L
    val <- .adghNLL(p, pinfo, studies, rxMod, output_var, grid, cores)
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
    function(p) .adghFDGrad(p, pinfo, studies, rxMod, output_var, grid, cores,
                              .ctl$grad_h, use_central)
  } else {
    function(p) .adghGrad(p, pinfo, studies, sensModel, rxMod, output_var,
                           grid, cores, .ctl$grad_h)
  }

  grad_label <- if (!want_grad) "none"
                else if (!is.null(sensModel)) "Analytical"
                else if (use_central) "CFD"
                else "FD"
  n_total_nodes <- if (pinfo$n_eta > 0L) n_nodes^pinfo$n_eta else 1L
  message("=== admixr2: Aggregate Data Modeling (GH) ===")
  message(sprintf("  Obs units: %d | Params: %d | Nodes: %d^%d=%d | Cores: %d | Grad: %s | Restarts: %d",
                  length(studies), length(ov$p0),
                  n_nodes, pinfo$n_eta, n_total_nodes,
                  cores, grad_label, .ctl$n_restarts))
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
    opt_raw <- nlmixr2est::nlmixrWithTiming("adgh", {
      nloptr::nloptr(x0 = p0_sc, eval_f = eval_f_sc,
                     eval_grad_f = eval_grad_sc,
                     lb = lb_sc, ub = ub_sc,
                     opts = list(algorithm = .ctl$algorithm,
                                 ftol_rel  = .ctl$ftol_rel,
                                 maxeval   = .ctl$maxeval))
    })
    opt <- list(objective  = opt_raw$objective,
                solution   = opt_raw$solution * sc,
                message    = opt_raw$message,
                all_traces = list(list(restart_id = 1L,
                                       nll_trace  = .nll_trace,
                                       par_trace  = .par_trace)))
    if (.ctl$print > 0L) {
      row <- .admProgressRow(sprintf("%04d \u2713", .iter), opt$objective, opt$solution, pinfo)
      if (!is.null(row)) message(paste0(row, "\n",
        .admProgressTimingRow((proc.time() - t0)["elapsed"], pinfo)))
    }
  } else {
    .admSetupDaemons(.ctl, .ctl$n_restarts)
    on.exit(.admStopDaemons(), add = TRUE)
    opt <- .admRunRestarts(
      worker_fn  = .adghRestartWorker,
      p0         = ov$p0, ov = ov, pinfo = pinfo,
      .ctl       = .ctl, ui = .ui, studies = studies,
      extra_args = list(
        n_nodes          = n_nodes,
        algorithm        = .ctl$algorithm,
        ftol_rel         = .ctl$ftol_rel,
        maxeval          = .ctl$maxeval,
        use_grad         = want_grad,
        use_central      = use_central,
        use_pure_fd      = use_pure_fd,
        grad_h           = .ctl$grad_h,
        grad_bounds      = .ctl$grad_bounds,
        output_var       = output_var,
        print_progress   = TRUE,
        print            = .ctl$print,
        cores            = .ctl$cores,
        rxMod_direct     = rxMod,
        sensModel_direct = sensModel
      )
    )
    .admStopDaemons()
    .iter <- opt$n_iter
  }

  t_opt  <- (proc.time() - t0)["elapsed"]
  final  <- .admUnpack(opt$solution, pinfo)
  fullTheta <- .admFullTheta(final, pinfo)
  p_hat  <- setNames(opt$solution, names(ov$p0))

  t0_cov <- proc.time()
  .cov <- if (.ctl$covMethod == "r") {
    np_cov    <- length(pinfo$struct_names) + length(pinfo$sigma_names)
    use_grad_cov <- want_grad && !is.null(sensModel)
    n_evals   <- if (use_grad_cov) np_cov + 1L
                 else { n_off <- np_cov * (np_cov - 1L) / 2L; 1L + 2L * np_cov + 4L * n_off }
    evals_lbl <- if (use_grad_cov) "gradient evaluations" else "NLL evaluations"
    hess_lbl  <- if (!use_grad_cov) "" else if (!is.null(sensModel)) ", Analytical-Hessian" else ", FD-Hessian"
    message(sprintf("  Computing covariance (R method%s, %d %s)", hess_lbl, n_evals, evals_lbl))
    tryCatch(
      .adghCalcCov(p_hat, pinfo, studies, sensModel, rxMod, output_var, grid, cores,
                   use_grad    = use_grad_cov,
                   grad_h      = .ctl$cov_h,
                   cov_h_outer = .ctl$cov_h_outer),
      error = function(e) { warning("adghCalcCov failed: ", conditionMessage(e)); NULL })
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
  .ret$est        <- "adgh"
  .ret$ofvType    <- "adgh"
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
                        n_nodes        = n_nodes,
                        n_sim          = 5000L,
                        sampling       = "sobol",
                        n_gh           = n_total_nodes)

  nlmixr2est::.nlmixr2FitUpdateParams(.ret)
  nmObjHandleControlObject.adghControl(.ctl, .ret)
  if (exists("control", .ui)) rm(list = "control", envir = .ui)
  .ret$control <- .admToFoceiControl(.ctl)
  .focei_model <- suppressMessages(tryCatch(.ui$foceiModel, error = function(e) NULL))
  if (!is.null(.focei_model)) .ret$model <- .focei_model

  .fit <- nlmixr2est::nlmixr2CreateOutputFromUi(
    .ui, data = if (multi_out) admData(.admOutputVars(.ui)) else admData(),
    control = .ret$control,
    table = .ret$table, env = .ret, est = "adgh")

  .fit$env$method   <- "adgh"
  .fit$env$studies  <- studies
  .fit$env$admExtra <- .ret$admExtra
  .admAttachParHist(.fit, .ret$admExtra$all_traces, .ret$admExtra$par_names, .ui)
  # Store observed + predicted aggregate moments (E vector, V matrix) per study.
  .admAttachAggData(.fit, .ret$admExtra, .ui)
  .old_cls <- class(.fit)
  .new_cls <- c("admFit", .old_cls)
  attr(.new_cls, ".foceiEnv") <- attr(.old_cls, ".foceiEnv")
  class(.fit) <- .new_cls

  .stats <- .admCalcObjStats(opt$objective, length(ov$p0), studies)
  row.names(.stats$objDf) <- "adgh"
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
