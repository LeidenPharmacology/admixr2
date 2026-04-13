skip_if_not_installed("rxode2")
skip_on_cran()

test_that("NLL evaluates to finite at initial parameters", {
  env <- .int_setup()
  expect_true(is.finite(env$nll_p0))
  expect_gt(env$nll_p0, 0)
})

test_that("NLL at true params < NLL at substantially perturbed params", {
  env <- .int_setup()
  expect_lt(env$nll_p0, env$nll_p_bad)
})

test_that("NLL with non-PD perturbed omega returns Inf", {
  env <- .int_setup()
  expect_true(is.infinite(env$nll_p_nonpd) && env$nll_p_nonpd > 0)
})

test_that("NLL evaluation produces no R-level message/warning output", {
  env <- .int_setup()
  expect_length(env$nll_warnings, 0L)
})
