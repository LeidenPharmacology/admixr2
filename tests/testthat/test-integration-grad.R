skip_if_not_installed("rxode2")
skip_on_cran()

# Setup functions (.int_grad_setup, .int_lincmt_setup) live in helper-integration.R
# so they are accessible to all test files.
# All rxSolve calls happen there; tests only check pre-computed values.

# ---- MC gradient tests -------------------------------------------------------

test_that("admGrad: finite at initial parameters", {
  env <- .int_grad_setup()
  expect_true(all(is.finite(env$g_ana_p0)))
})

test_that("admGrad vs CFD of admNLL: ratio within 5% (ML denominator, same z_list)", {
  env   <- .int_grad_setup()
  ratio <- env$g_ana_p0 / env$g_fd_p0
  ok    <- is.finite(ratio) & abs(env$g_fd_p0) > 1e-6
  ratio_ok <- ratio[ok]
  expect_true(all(abs(ratio_ok - 1) < 0.05),
    info = sprintf("max ratio deviation: %.4f (param: %s)",
                   max(abs(ratio_ok - 1)),
                   names(ratio_ok)[which.max(abs(ratio_ok - 1))]))
})

test_that("ODE: .admLoadSensModel() returns non-NULL", {
  env <- .int_grad_setup()
  expect_false(is.null(env$sensModel),
    info = "ODE sensitivity model should be available (ODE sens equations always present)")
})

test_that("admLoadSensModel: foceiModel pinned to .adm_pin_env after load", {
  env     <- .int_grad_setup()
  pin_key <- paste0("focei_", digest::digest(env$ui$lstExpr))
  expect_true(exists(pin_key, envir = admixr2:::.adm_pin_env, inherits = FALSE),
    info = "foceiModel companions must be pinned in .adm_pin_env to prevent Windows GC heap corruption")
  pinned <- get(pin_key, envir = admixr2:::.adm_pin_env, inherits = FALSE)
  expect_false(is.null(pinned$inner),
    info = "pinned foceiModel must have non-NULL $inner (sens equations)")
})

test_that("admLoadSensModel: in-memory cache returns identical result on repeat call", {
  env   <- .int_grad_setup()
  sens1 <- admixr2:::.admLoadSensModel(env$ui)
  sens2 <- admixr2:::.admLoadSensModel(env$ui)
  expect_true(identical(sens1, sens2),
    info = "second call should return same object from .adm_pin_env cache, not re-read from disk")
})

test_that("admLoadSensModel: sens result cached in .adm_pin_env", {
  env      <- .int_grad_setup()
  sens_key <- paste0("sens_", digest::digest(env$ui$lstExpr))
  expect_true(exists(sens_key, envir = admixr2:::.adm_pin_env, inherits = FALSE),
    info = "sens model result must be cached in .adm_pin_env for in-memory cache to work")
})

# ---- linCmt sensitivity model tests -----------------------------------------

test_that("linCmt: .admNLL() is finite (simulationModel multi-subject fix)", {
  env <- .int_lincmt_setup()
  expect_true(is.finite(env$nll_p0))
})

test_that("linCmt: .admLoadSensModel() returns non-NULL", {
  env <- .int_lincmt_setup()
  expect_false(is.null(env$sensModel),
    info = "linCmt sensitivity model should be available (inner has linCmtB jacobians)")
})

test_that("linCmt: admLoadSensModel pins foceiModel to .adm_pin_env", {
  env     <- .int_lincmt_setup()
  pin_key <- paste0("focei_", digest::digest(env$ui$lstExpr))
  expect_true(exists(pin_key, envir = admixr2:::.adm_pin_env, inherits = FALSE),
    info = "linCmt foceiModel companions must be pinned (different digest from ODE model)")
})

# ---- IRMC inner gradient tests -----------------------------------------------

test_that("irmcInnerGrad: finite at initial parameters", {
  env <- .int_grad_setup()
  if (!env$proposals_ok) skip("proposal draw failed")
  expect_true(all(is.finite(env$g_irmc_ana)))
})

test_that("irmcInnerGrad vs FD of irmcNLL: ratio within 5% (proposals fixed, no extra rxSolve)", {
  env <- .int_grad_setup()
  if (!env$proposals_ok) skip("proposal draw failed")

  ratio <- env$g_irmc_ana / env$g_irmc_fd
  bad   <- names(ratio)[abs(ratio - 1) > 0.05]
  expect_equal(length(bad), 0L,
    info = paste("Params with |ratio - 1| > 0.05:",
                 paste(sprintf("%s=%.4f", bad, ratio[bad]), collapse = ", ")))
})

# ---- admSimulate unit test ---------------------------------------------------

test_that("admSimulate: returns n_sim x n_times matrix with expected analytic values at eta=0", {
  env  <- .int_grad_setup()
  pars <- admixr2:::.admUnpack(env$vec$p0, env$pinfo)

  n_rows     <- 6L
  eta_mat    <- matrix(0, nrow = n_rows, ncol = env$pinfo$n_eta,
                       dimnames = list(NULL, env$pinfo$eta_col_names))
  params_mat <- admixr2:::.admMakeParamsList(n_rows, env$pinfo, 1L)[[1]]

  cp_mat <- admixr2:::.admSimulate(
    env$rxMod, pars$struct, env$pinfo$sigma_names,
    eta_mat, env$studies[[1]], env$output_var, params_mat, cores = 1L
  )

  expect_true(is.matrix(cp_mat))
  expect_equal(dim(cp_mat), c(n_rows, length(env$times)))
  expect_true(all(is.finite(cp_mat)))
  for (i in seq_len(n_rows))
    expect_equal(cp_mat[i, ], cp_mat[1, ], tolerance = 1e-10)
  expect_equal(cp_mat[1, ], env$E_true, tolerance = 1e-4)
})

# ---- admNLLBatch consistency -------------------------------------------------

test_that("admNLLBatch: two identical p0 entries return equal finite NLL values", {
  env <- .int_grad_setup()
  p0  <- env$vec$p0

  result <- admixr2:::.admNLLBatch(
    list(p0, p0), env$pinfo, env$studies, env$z_list,
    env$rxMod, env$output_var, env$params_list, 1L
  )

  expect_length(result, 2L)
  expect_true(all(is.finite(result)))
  expect_gt(result[[1]], 0)
  expect_equal(result[[1]], result[[2]], tolerance = 1e-10)
})

# ---- admGradBatch consistency ------------------------------------------------

test_that("admGradBatch: single-element batch returns (1 x np) matrix of finite values", {
  env <- .int_grad_setup()
  p0  <- env$vec$p0

  result <- admixr2:::.admGradBatch(
    list(p0), env$pinfo, env$studies, env$z_list,
    env$rxMod, env$output_var, env$params_list, 1L, h = 1e-3
  )

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(1L, length(p0)))
  expect_true(all(is.finite(result)))
})

test_that("admGradBatch: multi-element batch returns (n x np) matrix with named columns", {
  env <- .int_grad_setup()
  p0  <- env$vec$p0

  result <- admixr2:::.admGradBatch(
    list(p0, p0, p0), env$pinfo, env$studies, env$z_list,
    env$rxMod, env$output_var, env$params_list, 1L, h = 1e-3
  )

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(3L, length(p0)))
  expect_equal(colnames(result), names(p0))
  expect_true(all(is.finite(result)))
})

test_that("admGradBatch: rows for identical p0 entries are equal", {
  env <- .int_grad_setup()
  p0  <- env$vec$p0

  result <- admixr2:::.admGradBatch(
    list(p0, p0), env$pinfo, env$studies, env$z_list,
    env$rxMod, env$output_var, env$params_list, 1L, h = 1e-3
  )

  expect_equal(result[1L, ], result[2L, ], tolerance = 1e-10)
})

# ---- admNLLBatch additional coverage -----------------------------------------

test_that("admNLLBatch: empty list returns numeric(0)", {
  env <- .int_grad_setup()
  result <- admixr2:::.admNLLBatch(
    list(), env$pinfo, env$studies, env$z_list,
    env$rxMod, env$output_var, env$params_list, 1L
  )
  expect_equal(result, numeric(0))
})

test_that("admNLLBatch: non-PD omega returns Inf for that entry", {
  env <- .int_grad_setup()
  p0  <- env$vec$p0

  p_nonpd <- p0
  p_nonpd[env$pinfo$omega_par_names[env$pinfo$chol_diag][1L]] <- -1e10

  result <- admixr2:::.admNLLBatch(
    list(p0, p_nonpd, p0), env$pinfo, env$studies, env$z_list,
    env$rxMod, env$output_var, env$params_list, 1L
  )

  expect_length(result, 3L)
  expect_true(is.finite(result[[1L]]))
  expect_equal(result[[2L]], Inf)
  expect_true(is.finite(result[[3L]]))
})
