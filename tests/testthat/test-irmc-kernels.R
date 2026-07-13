test_that("compute_mean_new_cpp: exp/log transform (type 0) gives p - log_orig", {
  p        <- c(log(5), log(10))
  log_orig <- c(log(4), log(8))
  n        <- length(p)
  # transform_type/low/hi must be vectors of length n — scalar causes OOB read in C++
  result   <- admixr2:::compute_mean_new_cpp(p, log_orig,
                                             rep(0L,  n),
                                             rep(0.0, n),
                                             rep(1.0, n))
  expected <- p - log_orig
  expect_equal(result, expected, tolerance = 1e-12)
})

test_that("compute_mean_new_cpp: no change (p == log_orig) gives mean_new = 0", {
  log_orig <- c(log(5), log(3))
  n        <- length(log_orig)
  result   <- admixr2:::compute_mean_new_cpp(log_orig, log_orig,
                                             rep(0L,  n),
                                             rep(0.0, n),
                                             rep(1.0, n))
  expect_equal(result, c(0.0, 0.0), tolerance = 1e-12)
})

test_that("compute_mean_new_cpp: expit transform (type 1) gives log(lo+(hi-lo)*sigmoid(p)) - log_orig", {
  lo <- 0.1; hi <- 5.0
  p  <- c(0.5, -1.0)
  log_orig <- log(c(1.0, 0.5))
  s        <- 1 / (1 + exp(-p))
  expected <- log(lo + (hi - lo) * s) - log_orig
  result   <- admixr2:::compute_mean_new_cpp(p, log_orig,
                                             as.integer(rep(1L, 2)),
                                             rep(lo, 2), rep(hi, 2))
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("irmc_inner_nll_cpp: uniform log-weights give same NLL as weighted_meancov + nll_cov_cpp", {
  set.seed(50)
  n_sim <- 300L; n_t <- 2L; n <- 50L

  rawpreds <- matrix(rnorm(n_sim * n_t, mean = 2), n_sim, n_t)
  bi_mat   <- matrix(0.0, n_sim, 1L)
  mean_new <- 0.0
  L_omega  <- matrix(1.0)          # 1x1 omega
  log_prop <- rep(0.0, n_sim)      # uniform weights
  E_obs    <- c(2.0, 2.1)
  A        <- matrix(c(0.15, 0.02, 0.02, 0.12), 2, 2)
  V_obs    <- A %*% t(A)

  irmc_nll <- .irmc_nll(
    rawpreds, bi_mat, mean_new, L_omega, log_prop,
    E_obs, V_obs, n = n, sigma_var = 0.0, sigma_type = 0L,
    kappa_delta = numeric(0), use_var = 0L
  )

  # With uniform weights: weighted_meancov gives ML mean/cov
  w   <- rep(1.0 / n_sim, n_sim)
  wmc <- admixr2:::weighted_meancov_cpp(rawpreds, w)
  direct_nll <- nll_cov_r(E_obs, V_obs, wmc$mu, wmc$V, n)

  expect_equal(irmc_nll, direct_nll, tolerance = 1e-6)
})

test_that("irmc_inner_nll_cpp: use_var=1 vs use_var=0 agree for diagonal V_obs", {
  set.seed(51)
  n_sim <- 200L; n_t <- 2L; n <- 40L

  rawpreds <- matrix(rnorm(n_sim * n_t, mean = 1.5), n_sim, n_t)
  bi_mat   <- matrix(0.0, n_sim, 1L)
  mean_new <- 0.0; L_omega <- matrix(1.0)
  log_prop <- rnorm(n_sim)
  E_obs    <- c(1.5, 1.6)
  v_obs    <- c(0.04, 0.05)
  V_obs    <- diag(v_obs)

  nll_cov <- .irmc_nll(
    rawpreds, bi_mat, mean_new, L_omega, log_prop,
    E_obs, V_obs, n, 0.0, 0L, numeric(0), use_var = 0L
  )
  nll_var <- .irmc_nll(
    rawpreds, bi_mat, mean_new, L_omega, log_prop,
    E_obs, V_obs, n, 0.0, 0L, numeric(0), use_var = 1L
  )
  # Cov and var paths differ in numerical approach (Cholesky vs elementwise) for
  # small n_sim -> allow a modest tolerance.
  expect_equal(nll_cov, nll_var, tolerance = 1e-2)
})

test_that("irmc_inner_nll_cpp: proportional sigma (sigma_type=1) adds sigma_var*mu^2, not sigma_var", {
  set.seed(55)
  n_sim <- 300L; n_t <- 2L; n <- 50L
  rawpreds <- matrix(rnorm(n_sim * n_t, mean = 3), n_sim, n_t)  # mu != 1 so prop != add
  bi_mat   <- matrix(0.0, n_sim, 1L)
  log_prop <- rep(0.0, n_sim)
  E_obs    <- c(3.0, 3.1)
  V_obs    <- diag(c(0.09, 0.10))
  sv <- 0.05

  nll_add  <- .irmc_nll(rawpreds, bi_mat, 0.0, matrix(1.0), log_prop,
                                            E_obs, V_obs, n, sv, 0L, numeric(0), 0L)
  nll_prop <- .irmc_nll(rawpreds, bi_mat, 0.0, matrix(1.0), log_prop,
                                            E_obs, V_obs, n, sv, 1L, numeric(0), 0L)

  # With mu ~3, proportional adds 0.05*9=0.45 vs additive adds 0.05: clearly different
  expect_true(is.finite(nll_add))
  expect_true(is.finite(nll_prop))
  expect_false(isTRUE(all.equal(nll_add, nll_prop, tolerance = 1e-4)))
})

test_that("irmc_inner_nll_cpp: non-empty kappa_delta shifts predicted mean and changes NLL", {
  set.seed(56)
  n_sim <- 200L; n_t <- 2L; n <- 50L
  rawpreds <- matrix(rnorm(n_sim * n_t, mean = 2), n_sim, n_t)
  bi_mat   <- matrix(0.0, n_sim, 1L)
  log_prop <- rep(0.0, n_sim)
  E_obs    <- c(2.0, 2.1)
  V_obs    <- diag(c(0.04, 0.05))

  nll_no_kappa <- .irmc_nll(rawpreds, bi_mat, 0.0, matrix(1.0), log_prop,
                                               E_obs, V_obs, n, 0.0, 0L, numeric(0), 0L)
  nll_kappa    <- .irmc_nll(rawpreds, bi_mat, 0.0, matrix(1.0), log_prop,
                                               E_obs, V_obs, n, 0.0, 0L, c(0.5, -0.3), 0L)

  expect_true(is.finite(nll_no_kappa))
  expect_true(is.finite(nll_kappa))
  expect_false(isTRUE(all.equal(nll_no_kappa, nll_kappa, tolerance = 1e-6)))
})

test_that("irmc_inner_nll_cpp / logdmvnorm_batch_cpp: zero mean_new (no_ws invariant) is finite", {
  # has_kappa=TRUE models always pass mean_new = rep(0, n_eta) — no weight-shift.
  # Confirm irmc_inner_nll_cpp and logdmvnorm_batch_cpp handle this without crash.
  set.seed(99)
  n_sim <- 50L; n_t <- 2L; n_eta <- 2L; n <- 40L

  bi_mat   <- matrix(rnorm(n_sim * n_eta), n_sim, n_eta)
  rawpreds <- matrix(rnorm(n_sim * n_t, mean = 2), n_sim, n_t)
  log_prop <- rep(0.0, n_sim)
  L_prop   <- diag(2)  # identity Cholesky

  mean_new_full <- rep(0.0, n_eta)
  L_omega       <- diag(2) * 0.9

  nll <- .irmc_nll(
    rawpreds, bi_mat, mean_new_full, L_omega, log_prop,
    c(2.0, 2.1), diag(c(0.04, 0.05)), n,
    0.0, 0L, numeric(0), 0L
  )
  expect_true(is.finite(nll))

  # Also confirm logdmvnorm_batch_cpp accepts full-length zero mean without crash
  ld <- admixr2:::logdmvnorm_batch_cpp(bi_mat, mean_new_full, L_omega)
  expect_equal(length(ld), n_sim)
  expect_true(all(is.finite(ld)))
})

test_that("irmc_inner_nll_cpp: additive sigma shifts V diagonal", {
  set.seed(52)
  n_sim <- 200L; n_t <- 2L; n <- 50L
  rawpreds <- matrix(rnorm(n_sim * n_t, mean = 1), n_sim, n_t)
  bi_mat   <- matrix(0.0, n_sim, 1L)
  log_prop <- rep(0.0, n_sim)
  E_obs    <- c(1.0, 1.1)
  V_obs    <- diag(c(0.01, 0.02))
  sv <- 0.05

  nll_0  <- .irmc_nll(rawpreds, bi_mat, 0.0, matrix(1.0), log_prop,
                                         E_obs, V_obs, n, 0.0, 0L, numeric(0), 0L)
  nll_sv <- .irmc_nll(rawpreds, bi_mat, 0.0, matrix(1.0), log_prop,
                                         E_obs, V_obs, n, sv, 0L, numeric(0), 0L)

  # With sigma, V_pred diagonal is larger -> lower sensitivity to residual, different NLL
  expect_false(isTRUE(all.equal(nll_0, nll_sv, tolerance = 1e-6)))
  expect_true(is.finite(nll_sv))
})

test_that("irmc_inner_nll_cpp: lognormal sigma (sigma_type=2) is finite and differs from additive", {
  # lnorm: mu *= exp(sv/2), V_diag += mu_adj^2 * (exp(sv) - 1)
  # With mu ~3, lnorm adds 3^2*(exp(0.05)-1) ~0.46 vs additive 0.05 -- clearly different.
  set.seed(57)
  n_sim <- 300L; n_t <- 2L; n <- 50L
  rawpreds <- matrix(rnorm(n_sim * n_t, mean = 3), n_sim, n_t)
  bi_mat   <- matrix(0.0, n_sim, 1L)
  log_prop <- rep(0.0, n_sim)
  E_obs    <- c(3.0, 3.1)
  V_obs    <- diag(c(0.09, 0.10))
  sv       <- 0.05

  nll_add   <- .irmc_nll(rawpreds, bi_mat, 0.0, matrix(1.0), log_prop,
                                            E_obs, V_obs, n, sv, 0L, numeric(0), 0L)
  nll_lnorm <- .irmc_nll(rawpreds, bi_mat, 0.0, matrix(1.0), log_prop,
                                            E_obs, V_obs, n, sv, 2L, numeric(0), 0L)

  expect_true(is.finite(nll_lnorm))
  expect_false(isTRUE(all.equal(nll_add, nll_lnorm, tolerance = 1e-4)))
})
