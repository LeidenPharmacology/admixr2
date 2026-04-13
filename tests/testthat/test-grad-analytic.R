# Gradient correctness tests — Tier 1 (no rxode2, no rxSolve).
#
# Uses an analytic 1-cmt bolus model to generate cp_mat and its eta-Jacobian
# without rxSolve, then verifies omega/sigma gradients against CRN FD of
# nll_cov_cpp using the same z draws.
#
# Covers the full chain:
#   Omega -> L -> eta = L*z -> cp_mat (analytic) -> nll_cov_cpp -> gradient

# ---- Analytic 1-cmt helpers --------------------------------------------------

.cp1cmt <- function(tcl, tv, eta_cl, eta_v, dose = 100, times = c(0.5, 1, 2, 4)) {
  cl <- exp(tcl + eta_cl); v <- exp(tv + eta_v)
  (dose / v) * exp(-(cl / v) * times)
}

# Partial derivatives w.r.t. etas (analytic chain rule)
.dcp_dcl <- function(tcl, tv, eta_cl, eta_v, dose = 100, times = c(0.5, 1, 2, 4)) {
  cl <- exp(tcl + eta_cl); v <- exp(tv + eta_v)
  cp <- (dose / v) * exp(-(cl / v) * times)
  cp * (-times * cl / v)
}
.dcp_dv <- function(tcl, tv, eta_cl, eta_v, dose = 100, times = c(0.5, 1, 2, 4)) {
  cl <- exp(tcl + eta_cl); v <- exp(tv + eta_v)
  cp <- (dose / v) * exp(-(cl / v) * times)
  cp * (times * cl / v - 1)
}

# Build cp_mat (n_sim x n_t) and D_mat (n_sim x (n_t * n_eta)) from z and L
.cpmat_and_D <- function(z, L, tcl0, tv0, dose = 100, times = c(0.5, 1, 2, 4)) {
  eta    <- z %*% t(L)
  cp_mat <- t(apply(eta, 1, function(e) .cp1cmt(tcl0, tv0, e[1], e[2], dose, times)))
  D_cl   <- t(apply(eta, 1, function(e) .dcp_dcl(tcl0, tv0, e[1], e[2], dose, times)))
  D_v    <- t(apply(eta, 1, function(e) .dcp_dv(tcl0, tv0, e[1], e[2], dose, times)))
  list(cp_mat = cp_mat, D_mat = cbind(D_cl, D_v))
}

# NLL via nll_cov_cpp and its derivatives dNLL/dV, dNLL/dmu
.nll_derivs <- function(cp_mat, E_obs, V_obs, n_obs, sigma_var = 0) {
  n_sim  <- nrow(cp_mat)
  mu     <- colMeans(cp_mat)
  cp_c   <- sweep(cp_mat, 2L, mu)
  V_pred <- crossprod(cp_c) / (n_sim - 1L)
  diag(V_pred) <- diag(V_pred) + sigma_var
  cholV  <- chol(V_pred)
  invV   <- chol2inv(cholV)
  r      <- E_obs - mu
  list(
    nll      = admixr2:::nll_cov_cpp(E_obs, V_obs, mu, V_pred, n_obs),
    mu       = mu,
    cp_c     = cp_c,
    dNLL_dmu = n_obs * as.numeric(-2 * invV %*% r),
    dNLL_dV  = n_obs * (invV - invV %*% (V_obs + tcrossprod(r)) %*% invV)
  )
}

# Perturb one omega optimizer parameter in L.
# Diagonal: optimizer p = 2*log(L_ii), so L_ii' = exp((p + delta)/2).
# Off-diagonal: optimizer p = L_ij directly.
.perturb_L <- function(L, row, col, delta) {
  L_new <- L
  if (row == col) {
    L_new[row, col] <- exp((2 * log(L[row, col]) + delta) / 2)
  } else {
    L_new[row, col] <- L_new[row, col] + delta
  }
  L_new
}

# ==============================================================================

test_that("MC omega gradient vs CRN FD of nll_cov_cpp (no rxode2)", {
  set.seed(77)
  n_sim <- 600L
  times <- c(0.5, 1, 2, 4)
  n_t   <- length(times)
  tcl0  <- log(5); tv0 <- log(20); dose <- 100

  # Slightly off-truth omega so gradients are non-trivial
  Omega <- matrix(c(0.09, 0.018, 0.018, 0.04), 2, 2)
  L     <- t(chol(Omega))
  z     <- matrix(rnorm(n_sim * 2), n_sim, 2)

  E_obs <- .one_cmt_mean(5, 20, dose, times) * 0.9   # perturb from truth
  V_obs <- diag((.3 * .one_cmt_mean(5, 20, dose, times))^2)
  n_obs <- 200L

  res <- .cpmat_and_D(z, L, tcl0, tv0, dose, times)
  drv <- .nll_derivs(res$cp_mat, E_obs, V_obs, n_obs)

  # Analytical gradient: correct z_diag_scale (factor-of-2 fix)
  neta1        <- c(1L, 2L, 2L)
  neta2        <- c(1L, 2L, 1L)
  z_diag_scale <- sweep(z, 2L, diag(L) / 2, "*")

  g_ana <- admixr2:::adm_grad_eta_omega_cpp(
    drv$cp_c, res$D_mat, z_diag_scale, z,
    drv$dNLL_dV, drv$dNLL_dmu, rep(0.0, n_t),
    neta1, neta2, n_t, 2L
  )

  # CRN FD reference: same z, perturb each omega optimizer param
  h <- 1e-4
  omega_entries <- list(c(1L, 1L), c(2L, 2L), c(2L, 1L))

  g_fd <- vapply(omega_entries, function(rc) {
    L_h <- .perturb_L(L, rc[1], rc[2], +h)
    L_l <- .perturb_L(L, rc[1], rc[2], -h)
    nll_h <- .nll_derivs(.cpmat_and_D(z, L_h, tcl0, tv0, dose, times)$cp_mat,
                         E_obs, V_obs, n_obs)$nll
    nll_l <- .nll_derivs(.cpmat_and_D(z, L_l, tcl0, tv0, dose, times)$cp_mat,
                         E_obs, V_obs, n_obs)$nll
    (nll_h - nll_l) / (2 * h)
  }, double(1))

  ratio <- g_ana$omega_grad / g_fd
  bad   <- which(abs(ratio - 1) > 0.05)
  expect_equal(length(bad), 0L,
    info = paste("omega_grad ratio > 5% for entries:", paste(bad, collapse = ", "),
                 "| ratios:", paste(round(ratio, 3), collapse = ", ")))
})

test_that("MC sigma gradient vs FD of nll_cov_cpp (no rxode2)", {
  set.seed(88)
  n_sim <- 600L
  times <- c(0.5, 1, 2, 4)
  tcl0  <- log(5); tv0 <- log(20); dose <- 100

  Omega <- matrix(c(0.09, 0, 0, 0.04), 2, 2)
  L     <- t(chol(Omega))
  z     <- matrix(rnorm(n_sim * 2), n_sim, 2)

  E_obs    <- .one_cmt_mean(5, 20, dose, times) * 0.9
  V_obs    <- diag((.3 * .one_cmt_mean(5, 20, dose, times))^2)
  n_obs    <- 200L
  sigma_var <- 0.01  # additive variance

  res <- .cpmat_and_D(z, L, tcl0, tv0, dose, times)
  drv <- .nll_derivs(res$cp_mat, E_obs, V_obs, n_obs, sigma_var)

  # Analytical sigma gradient on log scale: d(NLL)/d(log sigma_var) = sigma_var * sum(dNLL/dV_diag)
  g_ana_sigma <- sigma_var * sum(diag(drv$dNLL_dV))

  # FD on log(sigma_var) scale (optimizer encoding: sigma stored as log(sigma_var))
  h    <- 1e-4
  sv_h <- sigma_var * exp(+h)
  sv_l <- sigma_var * exp(-h)
  nll_h <- .nll_derivs(res$cp_mat, E_obs, V_obs, n_obs, sv_h)$nll
  nll_l <- .nll_derivs(res$cp_mat, E_obs, V_obs, n_obs, sv_l)$nll
  g_fd_sigma <- (nll_h - nll_l) / (2 * h)

  ratio <- g_ana_sigma / g_fd_sigma
  expect_equal(ratio, 1, tolerance = 0.01,
    info = sprintf("sigma gradient ratio = %.4f (expected ~1)", ratio))
})

test_that("omega diagonal gradient: z_diag_scale gives half the gradient of raw eta_mat", {
  # Regression for factor-of-2 bug documented in CLAUDE.md.
  # Passing eta_mat = z * L_ii (wrong) instead of z * L_ii/2 (correct)
  # doubles the diagonal omega gradient — verified here with analytic cp_mat.
  set.seed(55)
  n_sim <- 400L
  times <- c(0.5, 1, 2, 4)
  n_t   <- length(times)
  tcl0  <- log(5); tv0 <- log(20)

  L  <- matrix(c(0.3, 0, 0.06, 0.2), 2, 2)  # arbitrary lower-triangular
  z  <- matrix(rnorm(n_sim * 2), n_sim, 2)

  res <- .cpmat_and_D(z, L, tcl0, tv0)
  drv <- .nll_derivs(res$cp_mat,
                      E_obs = .one_cmt_mean(5, 20, 100, times),
                      V_obs = diag((.3 * .one_cmt_mean(5, 20, 100, times))^2),
                      n_obs = 150L)

  neta1          <- c(1L, 2L, 2L); neta2 <- c(1L, 2L, 1L)
  eta_mat_wrong  <- sweep(z, 2L, diag(L),     "*")   # bug: z * L_ii
  z_diag_correct <- sweep(z, 2L, diag(L) / 2, "*")   # fix: z * L_ii/2

  g_wrong <- admixr2:::adm_grad_eta_omega_cpp(
    drv$cp_c, res$D_mat, eta_mat_wrong, z,
    drv$dNLL_dV, drv$dNLL_dmu, rep(0.0, n_t),
    neta1, neta2, n_t, 2L)
  g_correct <- admixr2:::adm_grad_eta_omega_cpp(
    drv$cp_c, res$D_mat, z_diag_correct, z,
    drv$dNLL_dV, drv$dNLL_dmu, rep(0.0, n_t),
    neta1, neta2, n_t, 2L)

  # First two entries are diagonal (neta1 == neta2): wrong must be 2x correct
  for (k in c(1L, 2L))
    expect_equal(g_wrong$omega_grad[k], 2 * g_correct$omega_grad[k], tolerance = 1e-10,
      info = paste("Diagonal omega entry", k))
})
