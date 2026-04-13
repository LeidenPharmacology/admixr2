skip_if_not_installed("rxode2")
skip_on_cran()

# Setup lives in helper-integration.R. All rxSolve calls happen inside
# .int_grad_setup(); these tests only inspect pre-computed or cheap results.

# ---- .admSimulateSens: basic structure ----------------------------------------

test_that("admSimulateSens: returns list with cp_mat and dpred_list", {
  env  <- .int_grad_setup()
  if (is.null(env$sensModel)) skip("sens model unavailable")
  pars <- admixr2:::.admUnpack(env$vec$p0, env$pinfo)

  eta_mat <- matrix(0, nrow = 4L, ncol = env$pinfo$n_eta,
                    dimnames = list(NULL, env$pinfo$eta_col_names))

  res <- admixr2:::.admSimulateSens(
    env$sensModel, pars$struct, env$pinfo$sigma_names,
    eta_mat, env$studies[[1]], cores = 1L
  )

  expect_type(res, "list")
  expect_named(res, c("cp_mat", "dpred_list"))
})

test_that("admSimulateSens: cp_mat has correct dimensions", {
  env  <- .int_grad_setup()
  if (is.null(env$sensModel)) skip("sens model unavailable")
  pars <- admixr2:::.admUnpack(env$vec$p0, env$pinfo)

  n_sim <- 8L
  eta_mat <- matrix(0, nrow = n_sim, ncol = env$pinfo$n_eta,
                    dimnames = list(NULL, env$pinfo$eta_col_names))

  res <- admixr2:::.admSimulateSens(
    env$sensModel, pars$struct, env$pinfo$sigma_names,
    eta_mat, env$studies[[1]], cores = 1L
  )

  expect_true(is.matrix(res$cp_mat))
  expect_equal(dim(res$cp_mat), c(n_sim, length(env$times)))
})

test_that("admSimulateSens: dpred_list has n_eta matrices of correct dimensions", {
  env  <- .int_grad_setup()
  if (is.null(env$sensModel)) skip("sens model unavailable")
  pars <- admixr2:::.admUnpack(env$vec$p0, env$pinfo)

  n_sim <- 8L
  eta_mat <- matrix(0, nrow = n_sim, ncol = env$pinfo$n_eta,
                    dimnames = list(NULL, env$pinfo$eta_col_names))

  res <- admixr2:::.admSimulateSens(
    env$sensModel, pars$struct, env$pinfo$sigma_names,
    eta_mat, env$studies[[1]], cores = 1L
  )

  expect_length(res$dpred_list, env$pinfo$n_eta)
  for (j in seq_len(env$pinfo$n_eta))
    expect_equal(dim(res$dpred_list[[j]]), c(n_sim, length(env$times)))
})

test_that("admSimulateSens: all values finite", {
  env  <- .int_grad_setup()
  if (is.null(env$sensModel)) skip("sens model unavailable")
  pars <- admixr2:::.admUnpack(env$vec$p0, env$pinfo)

  n_sim <- 8L
  eta_mat <- matrix(0, nrow = n_sim, ncol = env$pinfo$n_eta,
                    dimnames = list(NULL, env$pinfo$eta_col_names))

  res <- admixr2:::.admSimulateSens(
    env$sensModel, pars$struct, env$pinfo$sigma_names,
    eta_mat, env$studies[[1]], cores = 1L
  )

  expect_true(all(is.finite(res$cp_mat)))
  for (j in seq_len(env$pinfo$n_eta))
    expect_true(all(is.finite(res$dpred_list[[j]])))
})

# ---- .admSimulateSens: correctness at eta = 0 --------------------------------

test_that("admSimulateSens: cp_mat at eta=0 matches analytic values", {
  env  <- .int_grad_setup()
  if (is.null(env$sensModel)) skip("sens model unavailable")
  pars <- admixr2:::.admUnpack(env$vec$p0, env$pinfo)

  n_sim <- 6L
  eta_mat <- matrix(0, nrow = n_sim, ncol = env$pinfo$n_eta,
                    dimnames = list(NULL, env$pinfo$eta_col_names))

  res <- admixr2:::.admSimulateSens(
    env$sensModel, pars$struct, env$pinfo$sigma_names,
    eta_mat, env$studies[[1]], cores = 1L
  )

  # All rows identical at eta=0
  for (i in seq_len(n_sim))
    expect_equal(res$cp_mat[i, ], res$cp_mat[1, ], tolerance = 1e-10)

  # Match analytic mean
  expect_equal(res$cp_mat[1, ], env$E_true, tolerance = 1e-4)
})

test_that("admSimulateSens: cp_mat agrees with admSimulate at non-zero eta", {
  env  <- .int_grad_setup()
  if (is.null(env$sensModel)) skip("sens model unavailable")
  pars <- admixr2:::.admUnpack(env$vec$p0, env$pinfo)

  set.seed(1L)
  n_sim   <- 10L
  eta_mat <- matrix(rnorm(n_sim * env$pinfo$n_eta, sd = 0.1),
                    nrow = n_sim, ncol = env$pinfo$n_eta,
                    dimnames = list(NULL, env$pinfo$eta_col_names))

  res_sens <- admixr2:::.admSimulateSens(
    env$sensModel, pars$struct, env$pinfo$sigma_names,
    eta_mat, env$studies[[1]], cores = 1L
  )

  params_mat <- admixr2:::.admMakeParamsList(n_sim, env$pinfo, 1L)[[1]]
  cp_sim <- admixr2:::.admSimulate(
    env$rxMod, pars$struct, env$pinfo$sigma_names,
    eta_mat, env$studies[[1]], env$output_var, params_mat, cores = 1L
  )

  expect_equal(res_sens$cp_mat, cp_sim, tolerance = 1e-4)
})

test_that("admSimulateSens: dpred_list sensitivities are non-zero at eta=0", {
  env  <- .int_grad_setup()
  if (is.null(env$sensModel)) skip("sens model unavailable")
  pars <- admixr2:::.admUnpack(env$vec$p0, env$pinfo)

  n_sim <- 4L
  eta_mat <- matrix(0, nrow = n_sim, ncol = env$pinfo$n_eta,
                    dimnames = list(NULL, env$pinfo$eta_col_names))

  res <- admixr2:::.admSimulateSens(
    env$sensModel, pars$struct, env$pinfo$sigma_names,
    eta_mat, env$studies[[1]], cores = 1L
  )

  # d(pred)/d(eta_j) at eta=0 should be non-trivially non-zero for a
  # mu-referenced lognormal model: d(C)/d(eta_cl) = -C * t * cl/v (< 0),
  # d(C)/d(eta_v) = -C (< 0).
  for (j in seq_len(env$pinfo$n_eta))
    expect_true(any(abs(res$dpred_list[[j]]) > 1e-6))
})

# ---- .admSimulateSens: n_sim consistency -------------------------------------

test_that("admSimulateSens: output dimensions scale with n_sim", {
  env  <- .int_grad_setup()
  if (is.null(env$sensModel)) skip("sens model unavailable")
  pars <- admixr2:::.admUnpack(env$vec$p0, env$pinfo)

  for (n_sim in c(1L, 5L, 20L)) {
    eta_mat <- matrix(0, nrow = n_sim, ncol = env$pinfo$n_eta,
                      dimnames = list(NULL, env$pinfo$eta_col_names))
    res <- admixr2:::.admSimulateSens(
      env$sensModel, pars$struct, env$pinfo$sigma_names,
      eta_mat, env$studies[[1]], cores = 1L
    )
    expect_equal(nrow(res$cp_mat), n_sim)
    expect_equal(nrow(res$dpred_list[[1]]), n_sim)
  }
})
