# ---- admControl --------------------------------------------------------------

test_that("admControl() returns correct class and key defaults", {
  ctl <- admControl()
  expect_s3_class(ctl, "admControl")
  expect_equal(ctl$n_sim,    5000L)
  expect_equal(ctl$sampling, "sobol")
  expect_equal(ctl$grad,     "sens")
  expect_equal(ctl$maxeval,  500L)
  expect_equal(ctl$seed,     12345L)
  expect_equal(ctl$n_restarts, 1L)
  expect_equal(ctl$workers,    1L)
  expect_equal(ctl$covMethod,  "r")
})

test_that("admControl(): grad != 'none' + BOBYQA switches to LBFGS", {
  ctl <- admControl(grad = "fd")
  expect_equal(ctl$algorithm, "NLOPT_LD_LBFGS")
})

test_that("admControl(): grad != 'none' + explicit non-BOBYQA algorithm kept", {
  ctl <- admControl(grad = "fd", algorithm = "NLOPT_LD_SLSQP")
  expect_equal(ctl$algorithm, "NLOPT_LD_SLSQP")
})

test_that("admControl(): grad = 'none' keeps BOBYQA", {
  ctl <- admControl(grad = "none")
  expect_equal(ctl$algorithm, "NLOPT_LN_BOBYQA")
})

test_that("admControl(): n_sim stored as integer", {
  ctl <- admControl(n_sim = 1000)
  expect_type(ctl$n_sim, "integer")
  expect_equal(ctl$n_sim, 1000L)
})

test_that("admControl(): negative n_sim errors via checkmate", {
  expect_error(admControl(n_sim = -1), regexp = "n_sim")
})

test_that("admControl(): zero n_sim errors", {
  expect_error(admControl(n_sim = 0), regexp = "n_sim")
})

test_that("admControl(): unknown argument errors with informative message", {
  expect_error(admControl(nonexistent_arg = 42), regexp = "unused argument")
})

test_that("admControl(): sampling match.arg works for abbreviation", {
  ctl <- admControl(sampling = "rnorm")
  expect_equal(ctl$sampling, "rnorm")
})

test_that("admControl(): invalid sampling errors", {
  expect_error(admControl(sampling = "invalid_method"))
})

test_that("admControl(): cov_h_outer default is eps^(1/5)", {
  ctl <- admControl()
  expect_equal(ctl$cov_h_outer, .Machine$double.eps^(1/5), tolerance = 1e-15)
})

test_that("admControl(): studies stored as-is", {
  st <- list(s1 = list(E = 1, V = 1, n = 10L, times = 1))
  ctl <- admControl(studies = st)
  expect_identical(ctl$studies, st)
})

test_that("admControl(): returnAdmr stored correctly", {
  ctl <- admControl(returnAdmr = TRUE)
  expect_true(ctl$returnAdmr)
})

# ---- adirmcControl ----------------------------------------------------------

test_that("adirmcControl() returns correct class and key defaults", {
  ctl <- adirmcControl()
  expect_s3_class(ctl, "adirmcControl")
  expect_equal(ctl$n_sim,           2500L)
  expect_equal(ctl$grad,            "analytical")
  expect_equal(ctl$phases,          c(2, 1, 0.5, 0.01))
  expect_equal(ctl$outer_iter,      50L)
  expect_equal(ctl$omega_expansion, 1.0)
  expect_equal(ctl$convcrit,        1e-5)
  expect_equal(ctl$max_worse,       5L)
  expect_equal(ctl$kappa_method,    "exact")
})

test_that("adirmcControl(): omega_expansion < 1 errors", {
  expect_error(adirmcControl(omega_expansion = 0.5), regexp = "omega_expansion")
})

test_that("adirmcControl(): omega_expansion = 1 is valid (boundary)", {
  expect_no_error(adirmcControl(omega_expansion = 1.0))
})

test_that("adirmcControl(): grad != 'none' + BOBYQA switches to LBFGS", {
  ctl <- adirmcControl(grad = "fd")
  expect_equal(ctl$algorithm, "NLOPT_LD_LBFGS")
})

test_that("adirmcControl(): outer_iter stored as integer", {
  ctl <- adirmcControl(outer_iter = 100)
  expect_type(ctl$outer_iter, "integer")
  expect_equal(ctl$outer_iter, 100L)
})

test_that("adirmcControl(): phases must be positive", {
  expect_error(adirmcControl(phases = c(2, -1)), regexp = "phases")
})

test_that("adirmcControl(): unknown argument errors", {
  expect_error(adirmcControl(bad_arg = 1), regexp = "unused argument")
})

test_that("adirmcControl(): cov_h and cov_h_outer stored", {
  ctl <- adirmcControl(cov_h = 0.01, cov_h_outer = 0.005)
  expect_equal(ctl$cov_h,       0.01)
  expect_equal(ctl$cov_h_outer, 0.005)
})

test_that("adirmcControl(): kappa_method accepts exact", {
  ctl <- adirmcControl(kappa_method = "exact")
  expect_equal(ctl$kappa_method, "exact")
})

test_that("adirmcControl(): invalid kappa_method errors", {
  expect_error(adirmcControl(kappa_method = "second-order"))
  expect_error(adirmcControl(kappa_method = "third-order"))
})

