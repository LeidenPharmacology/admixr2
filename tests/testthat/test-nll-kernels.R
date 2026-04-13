test_that("nll_cov_cpp agrees with pure-R reference", {
  set.seed(1)
  n_t  <- 3L
  n    <- 50L
  E_pred <- c(1.0, 2.0, 3.0)
  E_obs  <- c(1.1, 1.9, 3.05)
  # Build a valid PD covariance matrix
  A <- matrix(c(0.3, 0.1, 0.0,
                0.1, 0.4, 0.1,
                0.0, 0.1, 0.3), 3, 3)
  V_pred <- A %*% t(A)
  V_obs  <- V_pred * 0.8

  cpp_val <- admixr2:::nll_cov_cpp(E_obs, V_obs, E_pred, V_pred, n)
  r_val   <- nll_cov_r(E_obs, V_obs, E_pred, V_pred, n)
  expect_equal(cpp_val, r_val, tolerance = 1e-10)
})

test_that("nll_cov_cpp: when E_obs==E_pred and V_obs==V_pred, result = n*(logdet + p)", {
  # tr(V^-1 V) = p (dimension), and quadratic term r'V^-1 r = 0
  set.seed(2)
  n_t <- 2L; n <- 30L
  A <- matrix(c(0.2, 0.05, 0.05, 0.3), 2, 2)
  V <- A %*% t(A)
  E <- c(1.0, 2.0)

  logdet <- 2 * sum(log(diag(chol(V))))
  expected <- n * (logdet + n_t)

  cpp_val <- admixr2:::nll_cov_cpp(E, V, E, V, n)
  expect_equal(cpp_val, expected, tolerance = 1e-10)
})

test_that("nll_cov_cpp returns Inf for non-PD V_pred", {
  E   <- c(1.0, 2.0)
  V_bad <- matrix(c(-1, 0, 0, 1), 2, 2)  # not PD
  V_obs <- diag(2)
  result <- admixr2:::nll_cov_cpp(E, V_obs, E, V_bad, n = 10L)
  expect_true(is.infinite(result) && result > 0)
})

test_that("nll_var_cpp agrees with pure-R reference", {
  set.seed(3)
  n_t  <- 4L; n <- 80L
  v_pred <- c(0.1, 0.2, 0.15, 0.3)
  v_obs  <- v_pred * 0.9
  E_pred <- c(1.0, 2.0, 1.5, 3.0)
  E_obs  <- E_pred + c(0.05, -0.1, 0.02, 0.0)

  cpp_val <- admixr2:::nll_var_cpp(E_obs, v_obs, E_pred, v_pred, n)
  r_val   <- nll_var_r(E_obs, v_obs, E_pred, v_pred, n)
  expect_equal(cpp_val, r_val, tolerance = 1e-10)
})

test_that("nll_var_cpp: when E_obs==E_pred and v_obs==v_pred, result = n*sum(log(v)+1)", {
  n <- 20L
  v <- c(0.1, 0.2, 0.3)
  E <- c(1.0, 2.0, 3.0)
  expected <- n * sum(log(v) + 1)
  cpp_val  <- admixr2:::nll_var_cpp(E, v, E, v, n)
  expect_equal(cpp_val, expected, tolerance = 1e-10)
})

test_that("nll_var_cpp is always >= its minimum over v_pred", {
  set.seed(4)
  n <- 50L
  v_obs  <- c(0.1, 0.2, 0.3)
  E_obs  <- c(1.0, 2.0, 1.5)
  E_pred <- c(1.05, 1.9, 1.55)
  r2     <- (E_obs - E_pred)^2

  # Analytic optimum: v_opt[t] = v_obs[t] + r2[t]
  v_opt   <- v_obs + r2
  nll_opt <- admixr2:::nll_var_cpp(E_obs, v_obs, E_pred, v_opt, n)

  for (v_other in list(v_obs, v_opt * 1.5, v_opt * 0.7 + 0.01)) {
    nll_other <- admixr2:::nll_var_cpp(E_obs, v_obs, E_pred, v_other, n)
    expect_gte(nll_other, nll_opt - 1e-8)
  }
})

test_that("softmax_cpp: outputs sum to 1 and are in (0, 1)", {
  set.seed(5)
  lw <- rnorm(100)
  w  <- admixr2:::softmax_cpp(lw)
  expect_equal(sum(w), 1.0, tolerance = 1e-12)
  expect_true(all(w > 0))
  expect_true(all(w < 1))
})

test_that("softmax_cpp: permuting input permutes output identically", {
  set.seed(6)
  lw   <- rnorm(20)
  perm <- sample(20)
  w    <- admixr2:::softmax_cpp(lw)
  wp   <- admixr2:::softmax_cpp(lw[perm])
  expect_equal(wp, w[perm], tolerance = 1e-14)
})

test_that("softmax_cpp: agrees with pure-R reference", {
  set.seed(7)
  lw <- c(-1.0, 0.5, 2.0, -3.0, 1.2)
  expect_equal(admixr2:::softmax_cpp(lw), softmax_r(lw), tolerance = 1e-12)
})

test_that("softmax_cpp: non-finite inputs clipped to max finite value", {
  lw <- c(1.0, -Inf, 2.0, Inf)
  # should not error and should return valid weights
  w <- admixr2:::softmax_cpp(lw)
  expect_equal(sum(w), 1.0, tolerance = 1e-12)
  expect_true(all(is.finite(w)))
})

test_that("weighted_meancov_cpp: uniform weights give ML mean and biased cov", {
  set.seed(8)
  n <- 50L; d <- 3L
  F_mat <- matrix(rnorm(n * d), n, d)
  w     <- rep(1.0 / n, n)

  res <- admixr2:::weighted_meancov_cpp(F_mat, w)

  expect_equal(res$mu, colMeans(F_mat), tolerance = 1e-12)
  # biased ML cov (no n-1 divisor) = (n-1)/n * sample_cov
  expected_V <- cov(F_mat) * (n - 1) / n
  expect_equal(res$V, expected_V, tolerance = 1e-10)
})

test_that("weighted_meancov_cpp: V is symmetric and PSD", {
  set.seed(9)
  n <- 30L; d <- 4L
  F_mat <- matrix(rnorm(n * d), n, d)
  w     <- abs(rnorm(n)); w <- w / sum(w)

  res <- admixr2:::weighted_meancov_cpp(F_mat, w)

  expect_equal(res$V, t(res$V), tolerance = 1e-12)
  ev <- eigen(res$V, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(ev >= -1e-10))
})

test_that("weighted_meancov_cpp: all weight on one row gives that row as mean and V=0", {
  n <- 10L; d <- 2L
  F_mat <- matrix(1:(n * d), n, d) * 0.1
  w     <- c(1.0, rep(0.0, n - 1L))

  res <- admixr2:::weighted_meancov_cpp(F_mat, w)

  expect_equal(res$mu, F_mat[1, ], tolerance = 1e-12)
  expect_equal(res$V, matrix(0.0, d, d), tolerance = 1e-12)
})

test_that("logdmvnorm_batch_cpp: standard normal at 0 gives -d/2*log(2*pi)", {
  d  <- 3L
  bi <- matrix(0.0, 1L, d)
  mu <- rep(0.0, d)
  L  <- diag(d)
  ld <- admixr2:::logdmvnorm_batch_cpp(bi, mu, L)
  expected <- -d / 2 * log(2 * pi)
  expect_equal(ld, expected, tolerance = 1e-10)
})

test_that("logdmvnorm_batch_cpp: log-density at mean equals standard normal log-density at 0", {
  set.seed(10)
  d  <- 3L
  mu <- c(1.0, -0.5, 2.0)
  Om <- Omega_2x2[1:2, 1:2]
  Om3 <- diag(3); Om3[1:2, 1:2] <- Om
  L  <- t(chol(Om3))
  bi <- matrix(mu, nrow = 1)

  ld      <- admixr2:::logdmvnorm_batch_cpp(bi, mu, L)
  ld_ref  <- -d / 2 * log(2 * pi) - sum(log(diag(L)))  # at the mode, z=0
  expect_equal(ld, ld_ref, tolerance = 1e-10)
})

# ---- helper-analytical.R self-validation ------------------------------------
# Ensures pure-R reference implementations agree with C++ so that test failures
# from them indicate a real disagreement, not a bug in the helper itself.

test_that("helper: nll_cov_r matches nll_cov_cpp on known 3-time case", {
  E_obs  <- c(1.0, 2.0, 1.5)
  V_obs  <- diag(c(0.01, 0.02, 0.015))
  E_pred <- c(1.1, 1.9, 1.6)
  V_pred <- diag(c(0.02, 0.03, 0.025))
  n      <- 100L
  expect_equal(nll_cov_r(E_obs, V_obs, E_pred, V_pred, n),
               admixr2:::nll_cov_cpp(E_obs, V_obs, E_pred, V_pred, n),
               tolerance = 1e-10)
})

test_that("helper: nll_var_r matches nll_var_cpp on known 3-time case", {
  E_obs  <- c(1.0, 2.0, 1.5)
  v_obs  <- c(0.01, 0.02, 0.015)
  E_pred <- c(1.1, 1.9, 1.6)
  v_pred <- c(0.02, 0.03, 0.025)
  n      <- 100L
  expect_equal(nll_var_r(E_obs, v_obs, E_pred, v_pred, n),
               admixr2:::nll_var_cpp(E_obs, v_obs, E_pred, v_pred, n),
               tolerance = 1e-10)
})

test_that("helper: logdmvnorm_r matches logdmvnorm_batch_cpp on 3 samples", {
  set.seed(321)
  L  <- t(chol(matrix(c(0.09, 0.018, 0.018, 0.04), 2, 2)))
  bi <- matrix(rnorm(6), nrow = 3, ncol = 2)
  expect_equal(logdmvnorm_r(bi, c(0, 0), L),
               admixr2:::logdmvnorm_batch_cpp(bi, c(0, 0), L),
               tolerance = 1e-10)
})

test_that("logdmvnorm_batch_cpp: agrees with pure-R reference (batch)", {
  set.seed(11)
  d    <- 2L; n_sim <- 20L
  Omega <- Omega_2x2
  L     <- t(chol(Omega))
  mu    <- c(0.5, -0.3)
  bi    <- matrix(rnorm(n_sim * d), n_sim, d)

  cpp_ld <- admixr2:::logdmvnorm_batch_cpp(bi, mu, L)
  r_ld   <- logdmvnorm_r(bi, mu, L)
  expect_equal(cpp_ld, r_ld, tolerance = 1e-10)
})

test_that("logdmvnorm_batch_cpp: agrees with mnormt::dmnorm", {
  skip_if_not_installed("mnormt")
  set.seed(42)
  d     <- 2L; n_sim <- 50L
  Omega <- Omega_2x2
  L     <- t(chol(Omega))
  mu    <- c(0.3, -0.1)
  bi    <- matrix(rnorm(n_sim * d), n_sim, d)

  cpp_ld    <- admixr2:::logdmvnorm_batch_cpp(bi, mu, L)
  mnormt_ld <- mnormt::dmnorm(bi, mean = mu, varcov = Omega, log = TRUE)
  expect_equal(cpp_ld, mnormt_ld, tolerance = 1e-12)
})
