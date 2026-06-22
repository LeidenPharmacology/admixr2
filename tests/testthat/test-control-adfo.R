# ---- adfoControl -------------------------------------------------------------

test_that("adfoControl() returns correct class and key defaults", {
  ctl <- adfoControl()
  expect_s3_class(ctl, "adfoControl")
  expect_equal(ctl$grad,       "none")
  expect_equal(ctl$maxeval,    500L)
  expect_equal(ctl$seed,       12345L)
  expect_equal(ctl$cores,      1L)
  expect_equal(ctl$n_restarts, 1L)
  expect_equal(ctl$workers,    1L)
  expect_equal(ctl$covMethod,  "r")
  expect_equal(ctl$returnAdmr, FALSE)
})

test_that("adfoControl(): grad = 'analytical' defaults to LBFGS", {
  ctl <- adfoControl(grad = "analytical")
  expect_equal(ctl$algorithm, "NLOPT_LD_LBFGS")
})

test_that("adfoControl(): grad = 'fd' defaults to LBFGS", {
  ctl <- adfoControl(grad = "fd")
  expect_equal(ctl$algorithm, "NLOPT_LD_LBFGS")
})

test_that("adfoControl(): grad = 'cfd' defaults to LBFGS", {
  ctl <- adfoControl(grad = "cfd")
  expect_equal(ctl$algorithm, "NLOPT_LD_LBFGS")
})

test_that("adfoControl(): grad = 'none' defaults to BOBYQA", {
  ctl <- adfoControl(grad = "none")
  expect_equal(ctl$algorithm, "NLOPT_LN_BOBYQA")
})

test_that("adfoControl(): grad != 'none' + explicit gradient algorithm kept", {
  ctl <- adfoControl(grad = "fd", algorithm = "NLOPT_LD_SLSQP")
  expect_equal(ctl$algorithm, "NLOPT_LD_SLSQP")
})

test_that("adfoControl(): MMA selectable with a gradient method", {
  ctl <- adfoControl(algorithm = "NLOPT_LD_MMA", grad = "analytical")
  expect_equal(ctl$algorithm, "NLOPT_LD_MMA")
  expect_equal(ctl$grad, "analytical")
})

test_that("adfoControl(): gradient algorithm + grad 'none' falls back to BOBYQA", {
  ctl <- suppressMessages(adfoControl(algorithm = "NLOPT_LD_MMA", grad = "none"))
  expect_equal(ctl$algorithm, "NLOPT_LN_BOBYQA")
  expect_equal(ctl$grad, "none")
})

test_that("adfoControl(): derivative-free algorithm drops the gradient", {
  ctl <- suppressMessages(
    adfoControl(algorithm = "NLOPT_LN_NELDERMEAD", grad = "analytical"))
  expect_equal(ctl$algorithm, "NLOPT_LN_NELDERMEAD")
  expect_equal(ctl$grad, "none")
})

test_that("adfoControl(): invalid algorithm errors", {
  expect_error(adfoControl(algorithm = "NLOPT_LD_NOTREAL"),
               regexp = "not a valid nloptr algorithm")
})

test_that("adfoControl(): internal n_sim is 1L for .admRunRestarts() compat", {
  ctl <- adfoControl()
  expect_type(ctl$n_sim, "integer")
  expect_equal(ctl$n_sim, 1L)
})

test_that("adfoControl(): internal sampling is 'sobol' for .admRunRestarts() compat", {
  ctl <- adfoControl()
  expect_equal(ctl$sampling, "sobol")
})

test_that("adfoControl(): maxeval stored as integer", {
  ctl <- adfoControl(maxeval = 200)
  expect_type(ctl$maxeval, "integer")
  expect_equal(ctl$maxeval, 200L)
})

test_that("adfoControl(): workers stored as integer", {
  ctl <- adfoControl(workers = 2)
  expect_type(ctl$workers, "integer")
  expect_equal(ctl$workers, 2L)
})

test_that("adfoControl(): n_restarts stored as integer", {
  ctl <- adfoControl(n_restarts = 3)
  expect_type(ctl$n_restarts, "integer")
  expect_equal(ctl$n_restarts, 3L)
})

test_that("adfoControl(): cov_h_outer default is eps^(1/5)", {
  ctl <- adfoControl()
  expect_equal(ctl$cov_h_outer, .Machine$double.eps^(1/5), tolerance = 1e-15)
})

test_that("adfoControl(): studies stored as-is", {
  st <- list(s1 = list(E = 1, V = 1, n = 10L, times = 1))
  ctl <- adfoControl(studies = st)
  expect_identical(ctl$studies, st)
})

test_that("adfoControl(): returnAdmr stored correctly", {
  ctl <- adfoControl(returnAdmr = TRUE)
  expect_true(ctl$returnAdmr)
})

test_that("adfoControl(): unknown argument errors with informative message", {
  expect_error(adfoControl(nonexistent_arg = 42), regexp = "unused argument")
})

test_that("adfoControl(): invalid grad errors", {
  expect_error(adfoControl(grad = "newton"))
})

test_that("adfoControl(): zero cores errors via checkmate", {
  expect_error(adfoControl(cores = 0), regexp = "cores")
})

test_that("adfoControl(): negative maxeval errors via checkmate", {
  expect_error(adfoControl(maxeval = 0), regexp = "maxeval")
})

test_that("adfoControl(): invalid covMethod errors", {
  expect_error(adfoControl(covMethod = "hessian"))
})

test_that("adfoControl(): covMethod 'none' accepted", {
  ctl <- adfoControl(covMethod = "none")
  expect_equal(ctl$covMethod, "none")
})

test_that("adfoControl(): grad_h, grad_bounds, cov_h stored", {
  ctl <- adfoControl(grad_h = 1e-5, grad_bounds = 3, cov_h = 5e-4)
  expect_equal(ctl$grad_h,      1e-5)
  expect_equal(ctl$grad_bounds, 3)
  expect_equal(ctl$cov_h,       5e-4)
})

test_that("adfoControl(): getValidNlmixrCtl.adfo passes through adfoControl", {
  ctl <- adfoControl(maxeval = 99L)
  out <- admixr2:::getValidNlmixrCtl.adfo(ctl)
  expect_s3_class(out, "adfoControl")
  expect_equal(out$maxeval, 99L)
})

test_that("adfoControl(): getValidNlmixrCtl.adfo returns default for NULL", {
  out <- admixr2:::getValidNlmixrCtl.adfo(NULL)
  expect_s3_class(out, "adfoControl")
})

test_that("adfoControl(): getValidNlmixrCtl.adfo converts plain list", {
  out <- admixr2:::getValidNlmixrCtl.adfo(list(maxeval = 77L))
  expect_s3_class(out, "adfoControl")
  expect_equal(out$maxeval, 77L)
})

test_that("adfoControl(): getValidNlmixrCtl.adfo extracts nested adfoControl", {
  ctl <- adfoControl(maxeval = 66L)
  out <- admixr2:::getValidNlmixrCtl.adfo(list(ctl))
  expect_s3_class(out, "adfoControl")
  expect_equal(out$maxeval, 66L)
})

test_that("adfoControl(): getValidNlmixrCtl.adfo builds from nested list with studies key", {
  plain <- list(list(studies = list(), maxeval = 55L))
  out   <- admixr2:::getValidNlmixrCtl.adfo(plain)
  expect_s3_class(out, "adfoControl")
})

test_that("adfoControl(): getValidNlmixrCtl.adfo falls back to default for unrecognised input", {
  # list(NULL): control[[1]] is NULL → not adfoControl, not list-with-studies,
  # no named keys on outer list → falls through to adfoControl() default.
  out <- admixr2:::getValidNlmixrCtl.adfo(list(NULL))
  expect_s3_class(out, "adfoControl")
})

test_that("nmObjHandleControlObject.adfoControl: assigns control to env", {
  ctl <- adfoControl(maxeval = 22L)
  e   <- new.env(parent = emptyenv())
  admixr2:::nmObjHandleControlObject.adfoControl(ctl, e)
  expect_true(exists("adfoControl", envir = e, inherits = FALSE))
  expect_equal(get("adfoControl", envir = e)$maxeval, 22L)
})

test_that("nmObjGetControl.adfo: retrieves control from env", {
  ctl <- adfoControl(maxeval = 33L)
  e   <- new.env(parent = emptyenv())
  e$adfoControl <- ctl
  out <- admixr2:::nmObjGetControl.adfo(list(e))
  expect_s3_class(out, "adfoControl")
  expect_equal(out$maxeval, 33L)
})

test_that("nmObjGetControl.adfo: falls back to 'control' slot", {
  ctl <- adfoControl(maxeval = 44L)
  e   <- new.env(parent = emptyenv())
  e$control <- ctl
  out <- admixr2:::nmObjGetControl.adfo(list(e))
  expect_s3_class(out, "adfoControl")
  expect_equal(out$maxeval, 44L)
})

test_that("nmObjGetControl.adfo: errors when no control found", {
  e <- new.env(parent = emptyenv())
  expect_error(
    admixr2:::nmObjGetControl.adfo(list(e)),
    regexp = "cannot find adfo control"
  )
})
