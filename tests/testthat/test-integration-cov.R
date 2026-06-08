skip_if_not_installed("rxode2")
skip_on_cran()

# Setup in helper-integration.R. .admCalcCov() uses .admNLLBatch() (use_grad=FALSE)
# or .admGradBatch() (use_grad=TRUE) to compute a numerical Hessian, then
# returns 2*H^-1 restricted to struct + sigma parameters (omega Cholesky excluded).
# All rxSolve calls are batched: one call per study for all perturbed configs.
#
# .int_cov_setup() evaluates at p_cov where sigma_sd = 1 (not the true 0.1).
# At the true params, sigma contributes only 0.01 variance vs ~1.7 from IIV —
# making H[sigma,sigma] near-zero and non-PD regardless of step size or n_sim.
# Shifting sigma to sigma_sd = 1 makes it identifiable for structural tests.

# ---- NLL-FD Hessian (use_grad = FALSE, default) ------------------------------

test_that("admCalcCov NLL-FD: result is matrix of correct dimensions", {
  env_cov <- .int_cov_setup()
  result  <- env_cov$result_nll
  expect_false(is.null(result), info = "NLL-FD Hessian should be PD at sigma_sd=1")
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(env_cov$n_struct, env_cov$n_struct))
  expect_equal(rownames(result), env_cov$struct_names)
})

test_that("admCalcCov NLL-FD: result is symmetric", {
  env_cov <- .int_cov_setup()
  result  <- env_cov$result_nll
  expect_false(is.null(result), info = "NLL-FD Hessian should be PD at sigma_sd=1")
  expect_equal(result, t(result), tolerance = 1e-10)
})

test_that("admCalcCov NLL-FD: result is positive definite", {
  env_cov <- .int_cov_setup()
  result  <- env_cov$result_nll
  expect_false(is.null(result), info = "NLL-FD Hessian should be PD at sigma_sd=1")
  eigs <- eigen(result, only.values = TRUE)$values
  expect_true(all(eigs > 0))
})

test_that("admCalcCov NLL-FD: omega params excluded from returned matrix", {
  env_cov <- .int_cov_setup()
  result  <- env_cov$result_nll
  expect_false(is.null(result), info = "NLL-FD Hessian should be PD at sigma_sd=1")
  expect_false(any(rownames(result) %in% env_cov$env$pinfo$omega_par_names))
})

test_that("admCalcCov: covMethod='r' path reports omega SE exclusion", {
  env_cov <- .int_cov_setup()
  expect_message(
    suppressWarnings(admixr2:::.admCalcCov(
      env_cov$p_cov, env_cov$env$pinfo, env_cov$env$studies, env_cov$env$z_list,
      env_cov$env$rxMod, env_cov$env$output_var, env_cov$env$params_list, 1L,
      cov_n_sim = 5000L, use_grad = FALSE
    )),
    "omega \\(IIV\\) SEs are not computed"
  )
})

# ---- Grad-FD Hessian (use_grad = TRUE) ---------------------------------------

test_that("admCalcCov grad-FD: result dimensions match NLL-FD path", {
  env_cov  <- .int_cov_setup()
  res_nll  <- env_cov$result_nll
  res_grad <- env_cov$result_grad
  expect_false(is.null(res_nll),  info = "NLL-FD Hessian should be PD at sigma_sd=1")
  expect_false(is.null(res_grad), info = "grad-FD Hessian should be PD at sigma_sd=1")
  expect_equal(dim(res_grad), dim(res_nll))
})

# ---- Error path: non-finite NLL at p_hat -------------------------------------

test_that("admCalcCov: non-finite NLL at p_hat warns and returns NULL", {
  env <- .int_grad_setup()

  p_nonpd <- env$vec$p0
  p_nonpd[env$pinfo$omega_par_names[env$pinfo$chol_diag][1L]] <- -1e10

  result <- NULL
  expect_warning(
    { result <- admixr2:::.admCalcCov(
        p_nonpd, env$pinfo, env$studies, env$z_list,
        env$rxMod, env$output_var, env$params_list, 1L,
        cov_n_sim = NULL, use_grad = FALSE
      ) },
    "NLL not finite"
  )
  expect_null(result)
})
