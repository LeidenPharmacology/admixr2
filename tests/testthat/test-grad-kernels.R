test_that("adm_col_sq_sum_cpp agrees with colSums(m^2) for general matrix", {
  set.seed(20)
  m <- matrix(rnorm(60), 12, 5)
  expect_equal(admixr2:::adm_col_sq_sum_cpp(m), colSums(m^2), tolerance = 1e-12)
})

test_that("adm_col_sq_sum_cpp: single column", {
  m <- matrix(c(1.0, 2.0, 3.0), 3, 1)
  expect_equal(admixr2:::adm_col_sq_sum_cpp(m), sum(m^2), tolerance = 1e-12)
})

test_that("adm_col_sq_sum_cpp: single row", {
  m <- matrix(c(1.0, 2.0, 3.0, 4.0), 1, 4)
  expect_equal(admixr2:::adm_col_sq_sum_cpp(m), as.numeric(m^2), tolerance = 1e-12)
})

# ---- Central regression test: z_diag_scale vs eta_mat factor-of-2 bug -------
# CLAUDE.md documents: passing raw eta_mat (= z * L_ii) as the eta_mat argument
# gives a diagonal omega gradient that is 2x too large.
# Correct call uses z_diag_scale = sweep(z, 2L, diag(L)/2, "*") = z * L_ii/2.
# neta1/neta2 are integer vectors from eta_rows_df (all eta rows, incl. off-diag).

test_that("adm_grad_eta_omega_cpp: z_diag_scale (correct) gives half the gradient of eta_mat (wrong)", {
  set.seed(30)
  n_sim <- 500L; n_t <- 2L
  pinfo  <- admixr2:::.admParseIniDf(make_inidf_1eta())
  L_ii   <- sqrt(0.09)   # = 0.3
  L      <- matrix(L_ii)
  z      <- matrix(rnorm(n_sim), n_sim, 1)
  eta_mat <- z * L_ii   # wrong: what the buggy code would pass

  cp_c   <- matrix(rnorm(n_sim * n_t), n_sim, n_t)
  D_mat  <- matrix(rnorm(n_sim * n_t), n_sim, n_t)
  dNLL_dV     <- diag(c(0.5, 0.3))
  dNLL_dmu    <- c(0.1, -0.2)
  sigma_mu_scale <- c(0.0, 0.0)
  neta1 <- as.integer(pinfo$eta_rows_df$neta1)
  neta2 <- as.integer(pinfo$eta_rows_df$neta2)

  go_wrong <- admixr2:::adm_grad_eta_omega_cpp(
    cp_c, D_mat, eta_mat, z,
    dNLL_dV, dNLL_dmu, sigma_mu_scale,
    neta1, neta2, n_t, pinfo$n_eta
  )

  z_diag_scale <- z * (L_ii / 2)
  go_correct <- admixr2:::adm_grad_eta_omega_cpp(
    cp_c, D_mat, z_diag_scale, z,
    dNLL_dV, dNLL_dmu, sigma_mu_scale,
    neta1, neta2, n_t, pinfo$n_eta
  )

  # Diagonal omega gradient: wrong is exactly 2x correct
  expect_equal(go_wrong$omega_grad[1], 2 * go_correct$omega_grad[1], tolerance = 1e-10)
})

test_that("adm_grad_eta_omega_cpp: 2-eta model returns 3 omega_grad entries", {
  set.seed(31)
  n_sim <- 200L; n_t <- 3L
  pinfo <- admixr2:::.admParseIniDf(make_inidf_2eta())
  L     <- t(chol(Omega_2x2))
  z     <- matrix(rnorm(n_sim * pinfo$n_eta), n_sim, pinfo$n_eta)
  z_diag_scale <- sweep(z, 2L, diag(L) / 2, "*")

  cp_c  <- matrix(rnorm(n_sim * n_t), n_sim, n_t)
  n_col <- n_t * pinfo$n_eta
  D_mat <- matrix(rnorm(n_sim * n_col), n_sim, n_col)
  A     <- diag(c(0.5, 0.3, 0.2))
  dNLL_dV     <- A
  dNLL_dmu    <- rnorm(n_t)
  sigma_mu_scale <- rep(0.0, n_t)
  neta1 <- as.integer(pinfo$eta_rows_df$neta1)
  neta2 <- as.integer(pinfo$eta_rows_df$neta2)

  g <- admixr2:::adm_grad_eta_omega_cpp(
    cp_c, D_mat, z_diag_scale, z,
    dNLL_dV, dNLL_dmu, sigma_mu_scale,
    neta1, neta2, n_t, pinfo$n_eta
  )

  expect_true(all(is.finite(g$eta_grad)))
  expect_true(all(is.finite(g$omega_grad)))
  # 2 diagonal rows + 1 off-diagonal row in eta_rows_df -> length 3
  expect_equal(length(g$omega_grad), nrow(pinfo$eta_rows_df))
})

test_that("adm_grad_eta_omega_var_cpp: factor-of-2 bug same as cov version", {
  set.seed(32)
  n_sim <- 200L; n_t <- 3L
  pinfo  <- admixr2:::.admParseIniDf(make_inidf_1eta())
  L_ii   <- sqrt(0.09)
  z      <- matrix(rnorm(n_sim), n_sim, 1)
  eta_mat <- z * L_ii         # wrong (factor-of-2)
  z_diag_scale <- z * (L_ii / 2)  # correct

  cp_c  <- matrix(rnorm(n_sim * n_t), n_sim, n_t)
  D_mat <- matrix(rnorm(n_sim * n_t), n_sim, n_t)
  dNLL_dV_diag   <- c(0.5, 0.3, 0.2)
  dNLL_dmu       <- c(0.1, -0.2, 0.05)
  sigma_mu_scale <- rep(0.0, n_t)
  neta1 <- as.integer(pinfo$eta_rows_df$neta1)
  neta2 <- as.integer(pinfo$eta_rows_df$neta2)

  gw <- admixr2:::adm_grad_eta_omega_var_cpp(
    cp_c, D_mat, eta_mat, z,
    dNLL_dV_diag, dNLL_dmu, sigma_mu_scale,
    neta1, neta2, n_t, pinfo$n_eta
  )
  gc <- admixr2:::adm_grad_eta_omega_var_cpp(
    cp_c, D_mat, z_diag_scale, z,
    dNLL_dV_diag, dNLL_dmu, sigma_mu_scale,
    neta1, neta2, n_t, pinfo$n_eta
  )

  expect_equal(gw$omega_grad[1], 2 * gc$omega_grad[1], tolerance = 1e-10)
  expect_true(all(is.finite(gc$eta_grad)))
})

test_that("adm_grad_partial_cpp vs adm_grad_partial_var_cpp agree for diagonal dNLL_dV", {
  set.seed(33)
  n_sim <- 100L; n_t <- 2L
  cp_c    <- matrix(rnorm(n_sim * n_t), n_sim, n_t)
  dpred   <- matrix(rnorm(n_sim * n_t), n_sim, n_t)
  dV_diag <- c(0.5, 0.3)
  dV      <- diag(dV_diag)
  eff_dmu <- c(0.1, -0.2)
  inv_nm1 <- 1.0 / (n_sim - 1L)

  g_cov <- admixr2:::adm_grad_partial_cpp(cp_c, dpred, dV, eff_dmu, inv_nm1)
  g_var <- admixr2:::adm_grad_partial_var_cpp(cp_c, dpred, dV_diag, eff_dmu, inv_nm1)

  expect_equal(g_cov, g_var, tolerance = 1e-12)
})

test_that("irmc_grad_kernel_cpp: dNLL_dlw sums to zero (softmax Jacobian property)", {
  # d/dlw of (sum w_i f(x_i)) where sum w_i = 1: the row-sum of the Jacobian = 0.
  # This holds for any choice of F, mu, dNLL_dV, dNLL_dmu, d_mat, invO.
  set.seed(34)
  n_sim <- 50L; n_t <- 2L; n_eta <- 1L
  F_mat <- matrix(rnorm(n_sim * n_t), n_sim, n_t)
  lw    <- rnorm(n_sim)
  w     <- softmax_r(lw)
  mu    <- colSums(w * F_mat)
  # d_mat: demeaned eta samples (bi - mean_new), here just zeros (mean_new=0, bi~N(0,1))
  bi    <- matrix(rnorm(n_sim * n_eta), n_sim, n_eta)
  d_mat <- bi  # sweep(bi, 2, 0) = bi
  invO  <- matrix(1 / 0.09)  # inverse of 1x1 Omega = 0.09
  dNLL_dmu <- c(0.1, -0.2)
  dNLL_dV  <- diag(c(0.5, 0.3))
  eff_dNLL_dmu <- dNLL_dmu

  gk <- admixr2:::irmc_grad_kernel_cpp(F_mat, w, mu, d_mat, invO, eff_dNLL_dmu, dNLL_dV)

  expect_equal(sum(gk$dNLL_dlw), 0, tolerance = 1e-10)
  expect_true(all(is.finite(gk$dNLL_dlw)))
})

test_that("irmc_grad_kernel_cpp vs irmc_grad_kernel_var_cpp agree for diagonal dNLL_dV", {
  set.seed(35)
  n_sim <- 60L; n_t <- 3L; n_eta <- 1L
  F_mat <- matrix(rnorm(n_sim * n_t), n_sim, n_t)
  lw    <- rnorm(n_sim)
  w     <- softmax_r(lw)
  mu    <- colSums(w * F_mat)
  bi    <- matrix(rnorm(n_sim * n_eta), n_sim, n_eta)
  d_mat <- bi
  invO  <- matrix(1 / 0.09)
  eff_dNLL_dmu <- rnorm(n_t)
  dNLL_dV_diag <- c(0.4, 0.2, 0.3)
  dNLL_dV      <- diag(dNLL_dV_diag)

  gk_cov <- admixr2:::irmc_grad_kernel_cpp(F_mat, w, mu, d_mat, invO, eff_dNLL_dmu, dNLL_dV)
  gk_var <- admixr2:::irmc_grad_kernel_var_cpp(F_mat, w, mu, d_mat, invO, eff_dNLL_dmu, dNLL_dV_diag)

  expect_equal(gk_cov$dNLL_dlw,       gk_var$dNLL_dlw,       tolerance = 1e-10)
  expect_equal(gk_cov$dNLL_dmean_new, gk_var$dNLL_dmean_new, tolerance = 1e-10)
})
