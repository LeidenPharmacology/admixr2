# FO integration tests — Tier 2 (requires rxode2, slow)
# Skip on CRAN; each test guards with skip_on_cran() + skip_if_not_installed().
#
# Tests .adfoNLL() and .adfoGrad() against FD with the 1-cmt ODE model.
# Setup function .int_adfo_setup() lives in helper-integration.R.

# ==============================================================================

test_that("FO NLL is finite and positive at truth", {
  skip_on_cran()
  env <- .int_adfo_setup()
  expect_true(is.finite(env$nll_p0))
  expect_gt(env$nll_p0, 0)
})

test_that("FO NLL is finite and positive at perturbed params", {
  skip_on_cran()
  env <- .int_adfo_setup()
  expect_true(is.finite(env$nll_p_bad))
  expect_gt(env$nll_p_bad, 0)
})

test_that("FO NLL increases when struct theta perturbed away from truth", {
  skip_on_cran()
  env <- .int_adfo_setup()
  # NLL at true theta should be lower than at perturbed theta
  expect_lt(env$nll_p0, env$nll_p_bad)
})

test_that("FO analytical gradient is finite at p0", {
  skip_on_cran()
  env <- .int_adfo_setup()
  expect_true(all(is.finite(env$g_ana)))
})

test_that("FO analytical gradient direction agrees with FD (ratio within 5% for all params)", {
  skip_on_cran()
  env <- .int_adfo_setup()

  ratio <- env$g_ana / env$g_fd
  # Exclude near-zero FD gradients to avoid 0/0 issues
  ok <- abs(env$g_fd) > 1e-6
  if (sum(ok) == 0) skip("All FD gradients near-zero at p0")

  ratio_ok <- ratio[ok]
  expect_true(all(abs(ratio_ok - 1) < 0.05),
    info = sprintf("max ratio deviation: %.4f (param: %s)",
                   max(abs(ratio_ok - 1)),
                   names(ratio_ok)[which.max(abs(ratio_ok - 1))]))
})

test_that("FO analytical diagonal omega gradients agree with FD (ratio within 5%)", {
  skip_on_cran()
  env <- .int_adfo_setup()

  n_s <- length(env$pinfo$struct_names)
  n_e <- length(env$pinfo$sigma_names)
  diag_idx <- which(env$pinfo$chol_diag) + n_s + n_e

  ratio <- env$g_ana[diag_idx] / env$g_fd[diag_idx]
  ok <- abs(env$g_fd[diag_idx]) > 1e-6
  if (sum(ok) == 0) skip("All diagonal omega FD gradients near-zero at p0")

  ratio_ok <- ratio[ok]
  expect_true(all(abs(ratio_ok - 1) < 0.05),
    info = sprintf("max diagonal omega ratio deviation: %.4f (param: %s)",
                   max(abs(ratio_ok - 1)),
                   names(ratio_ok)[which.max(abs(ratio_ok - 1))]))
})

test_that("FO analytical sigma gradient agrees with FD (ratio within 5%)", {
  skip_on_cran()
  env <- .int_adfo_setup()

  n_s    <- length(env$pinfo$struct_names)
  n_e    <- length(env$pinfo$sigma_names)
  sig_idx <- n_s + seq_len(n_e)

  ratio <- env$g_ana[sig_idx] / env$g_fd[sig_idx]
  ok    <- abs(env$g_fd[sig_idx]) > 1e-6
  if (sum(ok) == 0L) skip("All sigma FD gradients near-zero at p0")

  ratio_ok <- ratio[ok]
  expect_true(all(abs(ratio_ok - 1) < 0.05),
    info = sprintf("max sigma ratio deviation: %.4f (param: %s)",
                   max(abs(ratio_ok - 1)),
                   names(ratio_ok)[which.max(abs(ratio_ok - 1))]))
})

test_that("FO analytical gradient matches FD for unpaired struct theta", {
  skip_on_cran()
  env <- .int_adfo_kappa_setup()
  fd_zero_tol <- 1e-6

  expect_false(env$pinfo$struct_has_eta[["tsc"]])

  ratio <- env$g_ana["tsc"] / env$g_fd["tsc"]
  if (!is.finite(env$g_fd["tsc"]) || abs(env$g_fd["tsc"]) <= fd_zero_tol)
    skip("Unpaired struct theta FD gradient near-zero at p0")

  expect_true(abs(ratio - 1) < 0.05,
    info = sprintf("unpaired struct theta ratio: %.4f", ratio))
})


test_that("FO .adfoNLL() produces no R-level warnings at true params", {
  skip_on_cran()
  env <- .int_adfo_setup()

  warns <- character(0)
  withCallingHandlers(
    admixr2:::.adfoNLL(env$p0, env$pinfo, env$studies, env$sensModel,
                        env$rxMod, env$output_var, env$params_list, 1L),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(warns, 0)
})
