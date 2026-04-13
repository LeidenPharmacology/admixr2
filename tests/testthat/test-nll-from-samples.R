# sigma_type codes: 0=additive, 1=proportional, 2=lognormal

test_that("nll_cov_from_samples_cpp agrees with manual two-step (no sigma)", {
  set.seed(40)
  n_sim <- 200L; n_t <- 3L; n <- 50L
  cp_mat <- make_cp_mat(n_sim, n_t)$mat
  E_obs  <- c(1.0, 2.0, 1.5)
  A      <- matrix(c(0.2, 0.05, 0, 0.05, 0.3, 0.05, 0, 0.05, 0.2), 3, 3)
  V_obs  <- A %*% t(A)

  # Manual two-step
  cp_c   <- sweep(cp_mat, 2L, colMeans(cp_mat))
  V_pred <- crossprod(cp_c) / n_sim
  E_pred <- colMeans(cp_mat)
  r_nll  <- nll_cov_r(E_obs, V_obs, E_pred, V_pred, n)

  cpp_nll <- admixr2:::nll_cov_from_samples_cpp(
    cp_mat, E_obs, V_obs, n, sigma_var = 0.0, sigma_type = 0L
  )
  expect_equal(cpp_nll, r_nll, tolerance = 1e-8)
})

test_that("nll_var_from_samples_cpp agrees with manual two-step (no sigma)", {
  set.seed(41)
  n_sim <- 200L; n_t <- 3L; n <- 50L
  cp_mat <- make_cp_mat(n_sim, n_t)$mat
  E_obs  <- c(1.0, 2.0, 1.5)
  v_obs  <- c(0.01, 0.02, 0.015)

  # Manual two-step
  cp_c   <- sweep(cp_mat, 2L, colMeans(cp_mat))
  v_pred <- admixr2:::adm_col_sq_sum_cpp(cp_c) / n_sim
  E_pred <- colMeans(cp_mat)
  r_nll  <- nll_var_r(E_obs, v_obs, E_pred, v_pred, n)

  cpp_nll <- admixr2:::nll_var_from_samples_cpp(
    cp_mat, E_obs, v_obs, n, sigma_var = 0.0, sigma_type = 0L
  )
  expect_equal(cpp_nll, r_nll, tolerance = 1e-8)
})

test_that("nll_cov_from_samples_cpp: additive sigma adds sigma_var to V diagonal", {
  set.seed(42)
  n_sim <- 200L; n_t <- 2L; n <- 50L
  cp_mat <- make_cp_mat(n_sim, n_t, seed = 42)$mat[, 1:2]
  E_obs  <- c(1.0, 2.0)
  V_obs  <- diag(c(0.01, 0.02))
  sv     <- 0.05

  cp_c   <- sweep(cp_mat, 2L, colMeans(cp_mat))
  V_pred <- crossprod(cp_c) / n_sim
  E_pred <- colMeans(cp_mat)
  V_pred_sv <- V_pred; diag(V_pred_sv) <- diag(V_pred_sv) + sv
  expected <- nll_cov_r(E_obs, V_obs, E_pred, V_pred_sv, n)

  cpp_nll <- admixr2:::nll_cov_from_samples_cpp(
    cp_mat, E_obs, V_obs, n, sigma_var = sv, sigma_type = 0L
  )
  expect_equal(cpp_nll, expected, tolerance = 1e-8)
})

test_that("nll_cov_from_samples_cpp: proportional sigma adds sigma_var*mu^2 to V diagonal", {
  set.seed(43)
  n_sim <- 200L; n_t <- 2L; n <- 50L
  cp_mat <- make_cp_mat(n_sim, n_t, seed = 43)$mat[, 1:2]
  E_obs  <- c(1.0, 2.0)
  V_obs  <- diag(c(0.01, 0.02))
  sv     <- 0.05

  cp_c   <- sweep(cp_mat, 2L, colMeans(cp_mat))
  V_pred <- crossprod(cp_c) / n_sim
  mu     <- colMeans(cp_mat)
  V_pred_sv <- V_pred; diag(V_pred_sv) <- diag(V_pred_sv) + sv * mu^2
  expected <- nll_cov_r(E_obs, V_obs, mu, V_pred_sv, n)

  cpp_nll <- admixr2:::nll_cov_from_samples_cpp(
    cp_mat, E_obs, V_obs, n, sigma_var = sv, sigma_type = 1L
  )
  expect_equal(cpp_nll, expected, tolerance = 1e-8)
})

test_that("nll_var_from_samples_cpp: additive sigma agrees with manual", {
  set.seed(44)
  n_sim <- 200L; n_t <- 2L; n <- 50L
  cp_mat <- make_cp_mat(n_sim, n_t, seed = 44)$mat[, 1:2]
  E_obs  <- c(1.0, 2.0)
  v_obs  <- c(0.01, 0.02)
  sv     <- 0.04

  cp_c   <- sweep(cp_mat, 2L, colMeans(cp_mat))
  v_pred <- admixr2:::adm_col_sq_sum_cpp(cp_c) / n_sim
  mu     <- colMeans(cp_mat)
  v_pred_sv <- v_pred + sv
  expected <- nll_var_r(E_obs, v_obs, mu, v_pred_sv, n)

  cpp_nll <- admixr2:::nll_var_from_samples_cpp(
    cp_mat, E_obs, v_obs, n, sigma_var = sv, sigma_type = 0L
  )
  expect_equal(cpp_nll, expected, tolerance = 1e-8)
})

test_that("nll_var_from_samples_cpp: proportional sigma adds sigma_var*mu^2 to v_pred", {
  set.seed(45)
  n_sim <- 200L; n_t <- 2L; n <- 50L
  cp_mat <- make_cp_mat(n_sim, n_t, seed = 45)$mat[, 1:2]
  E_obs  <- c(1.0, 2.0)
  v_obs  <- c(0.01, 0.02)
  sv     <- 0.05

  cp_c   <- sweep(cp_mat, 2L, colMeans(cp_mat))
  v_pred <- admixr2:::adm_col_sq_sum_cpp(cp_c) / n_sim
  mu     <- colMeans(cp_mat)
  v_pred_sv <- v_pred + sv * mu^2   # proportional: += sigma_var * mu^2
  expected  <- nll_var_r(E_obs, v_obs, mu, v_pred_sv, n)

  cpp_nll <- admixr2:::nll_var_from_samples_cpp(
    cp_mat, E_obs, v_obs, n, sigma_var = sv, sigma_type = 1L
  )
  expect_equal(cpp_nll, expected, tolerance = 1e-8)
})

test_that("nll_cov_from_samples_cpp returns Inf when predicted V is not PD", {
  # All identical rows -> zero variance -> non-PD V_pred
  n_sim <- 10L; n_t <- 2L
  cp_mat <- matrix(1.0, n_sim, n_t)
  E_obs  <- c(1.0, 1.0)
  V_obs  <- diag(2)
  result <- admixr2:::nll_cov_from_samples_cpp(cp_mat, E_obs, V_obs, 20L, 0.0, 0L)
  expect_true(is.infinite(result))
})

test_that("nll_cov_from_samples_cpp: lognormal sigma (sigma_type=2) agrees with manual", {
  # lnorm: mu *= exp(sv/2), V_diag += mu_adj^2 * (exp(sv) - 1)
  set.seed(46)
  n_sim <- 200L; n_t <- 2L; n <- 50L
  cp_mat <- make_cp_mat(n_sim, n_t, seed = 46)$mat[, 1:2]
  E_obs  <- c(1.0, 2.0)
  V_obs  <- diag(c(0.01, 0.02))
  sv     <- 0.05

  mu     <- colMeans(cp_mat)
  cp_c   <- sweep(cp_mat, 2L, mu)
  V_pred <- crossprod(cp_c) / n_sim
  mu_adj <- mu * exp(sv / 2)
  V_adj  <- V_pred
  diag(V_adj) <- diag(V_adj) + mu_adj^2 * (exp(sv) - 1)
  expected <- nll_cov_r(E_obs, V_obs, mu_adj, V_adj, n)

  cpp_nll <- admixr2:::nll_cov_from_samples_cpp(cp_mat, E_obs, V_obs, n, sv, 2L)
  expect_equal(cpp_nll, expected, tolerance = 1e-8)
})

test_that("nll_var_from_samples_cpp: lognormal sigma (sigma_type=2) agrees with manual", {
  set.seed(47)
  n_sim <- 200L; n_t <- 2L; n <- 50L
  cp_mat <- make_cp_mat(n_sim, n_t, seed = 47)$mat[, 1:2]
  E_obs  <- c(1.0, 2.0)
  v_obs  <- c(0.01, 0.02)
  sv     <- 0.04

  mu     <- colMeans(cp_mat)
  cp_c   <- sweep(cp_mat, 2L, mu)
  v_pred <- admixr2:::adm_col_sq_sum_cpp(cp_c) / n_sim
  mu_adj <- mu * exp(sv / 2)
  v_adj  <- v_pred + mu_adj^2 * (exp(sv) - 1)
  expected <- nll_var_r(E_obs, v_obs, mu_adj, v_adj, n)

  cpp_nll <- admixr2:::nll_var_from_samples_cpp(cp_mat, E_obs, v_obs, n, sv, 2L)
  expect_equal(cpp_nll, expected, tolerance = 1e-8)
})
