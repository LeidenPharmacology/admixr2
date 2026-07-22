# ---- resid_nodes -------------------------------------------------------------
# The residual quadrature's node count. Only a transform-both-sides endpoint uses
# it, but it is accepted (and validated) by all four controls so a user does not
# have to know which estimator routes through which code path.

test_that("resid_nodes: all four controls default to 81 and store an integer", {
  for (ctl in list(admControl(), adfoControl(), adghControl(), adirmcControl()))
    expect_identical(ctl$resid_nodes, 81L)
  # A double is coerced, like every other count parameter.
  expect_identical(admControl(resid_nodes = 31)$resid_nodes,    31L)
  expect_identical(adghControl(resid_nodes = 31)$resid_nodes,   31L)
  expect_identical(adfoControl(resid_nodes = 31)$resid_nodes,   31L)
  expect_identical(adirmcControl(resid_nodes = 31)$resid_nodes, 31L)
})

test_that("resid_nodes: a degenerate node count is refused, not passed through", {
  # 0 nodes would reach .adghNodes1() and divide by an empty weight vector, and a
  # handful is meaningless as a quadrature -- both would otherwise surface as an
  # obscure error thousands of lines from the control that caused them.
  for (f in list(admControl, adfoControl, adghControl, adirmcControl)) {
    expect_error(f(resid_nodes = 0L))
    expect_error(f(resid_nodes = -5L))
    expect_error(f(resid_nodes = c(31L, 81L)))
    expect_error(f(resid_nodes = "many"))
  }
})

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

test_that("admControl(): grad != 'none' defaults to LBFGS", {
  ctl <- admControl(grad = "fd")
  expect_equal(ctl$algorithm, "NLOPT_LD_LBFGS")
})

test_that("admControl(): grad != 'none' + explicit gradient algorithm kept", {
  ctl <- admControl(grad = "fd", algorithm = "NLOPT_LD_SLSQP")
  expect_equal(ctl$algorithm, "NLOPT_LD_SLSQP")
})

test_that("admControl(): grad = 'none' defaults to BOBYQA", {
  ctl <- admControl(grad = "none")
  expect_equal(ctl$algorithm, "NLOPT_LN_BOBYQA")
})

test_that("admControl(): default algorithm derived from grad", {
  expect_equal(admControl()$algorithm,             "NLOPT_LD_LBFGS")  # grad "sens"
  expect_equal(admControl(grad = "none")$algorithm, "NLOPT_LN_BOBYQA")
})

test_that("admControl(): MMA selectable with the default gradient", {
  ctl <- admControl(algorithm = "NLOPT_LD_MMA")
  expect_equal(ctl$algorithm, "NLOPT_LD_MMA")
  expect_equal(ctl$grad, "sens")
})

test_that("admControl(): gradient algorithm + grad 'none' falls back to BOBYQA", {
  ctl <- suppressMessages(admControl(algorithm = "NLOPT_LD_MMA", grad = "none"))
  expect_equal(ctl$algorithm, "NLOPT_LN_BOBYQA")
  expect_equal(ctl$grad, "none")
  expect_message(admControl(algorithm = "NLOPT_LD_MMA", grad = "none"),
                 regexp = "grad = 'none'")
})

test_that("admControl(): gradient algorithm + explicit grad kept", {
  ctl <- admControl(algorithm = "NLOPT_LD_MMA", grad = "fd")
  expect_equal(ctl$algorithm, "NLOPT_LD_MMA")
  expect_equal(ctl$grad, "fd")
})

test_that("admControl(): derivative-free algorithm drops the gradient", {
  ctl <- suppressMessages(
    admControl(algorithm = "NLOPT_LN_NELDERMEAD", grad = "sens"))
  expect_equal(ctl$algorithm, "NLOPT_LN_NELDERMEAD")
  expect_equal(ctl$grad, "none")
  expect_message(admControl(algorithm = "NLOPT_LN_NELDERMEAD", grad = "sens"),
                 regexp = "derivative-free")
})

test_that("admControl(): invalid algorithm errors", {
  expect_error(admControl(algorithm = "NLOPT_LD_NOTREAL"),
               regexp = "not a valid nloptr algorithm")
})

test_that("admControl(): AUGLAG meta-algorithm warns", {
  expect_warning(admControl(algorithm = "NLOPT_LD_AUGLAG", grad = "fd"),
                 regexp = "local_opts")
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

test_that("admControl(): workers stored as integer", {
  ctl <- admControl(workers = 2)
  expect_type(ctl$workers, "integer")
  expect_equal(ctl$workers, 2L)
})

test_that("admControl(): n_restarts stored as integer", {
  ctl <- admControl(n_restarts = 3)
  expect_type(ctl$n_restarts, "integer")
  expect_equal(ctl$n_restarts, 3L)
})

test_that("admControl(): cov_n_sim stored as integer", {
  ctl <- admControl(cov_n_sim = 20000)
  expect_type(ctl$cov_n_sim, "integer")
  expect_equal(ctl$cov_n_sim, 20000L)
})

test_that("admControl(): maxeval stored as integer", {
  ctl <- admControl(maxeval = 300)
  expect_type(ctl$maxeval, "integer")
  expect_equal(ctl$maxeval, 300L)
})

test_that("admControl(): cores stored as integer", {
  ctl <- admControl(cores = 2)
  expect_type(ctl$cores, "integer")
  expect_equal(ctl$cores, 2L)
})

test_that("admControl(): cores = 0 errors via checkmate", {
  expect_error(admControl(cores = 0), regexp = "cores")
})

test_that("admControl(): maxeval = 0 errors via checkmate", {
  expect_error(admControl(maxeval = 0), regexp = "maxeval")
})

test_that("admControl(): workers = 0 errors via checkmate", {
  expect_error(admControl(workers = 0), regexp = "workers")
})

test_that("admControl(): n_restarts = 0 errors via checkmate", {
  expect_error(admControl(n_restarts = 0), regexp = "n_restarts")
})

test_that("admControl(): cov_n_sim = 0 errors via checkmate", {
  expect_error(admControl(cov_n_sim = 0), regexp = "cov_n_sim")
})

test_that("admControl(): grad_h, grad_bounds, cov_h stored", {
  ctl <- admControl(grad_h = 1e-5, grad_bounds = 3, cov_h = 5e-4)
  expect_equal(ctl$grad_h,      1e-5)
  expect_equal(ctl$grad_bounds, 3)
  expect_equal(ctl$cov_h,       5e-4)
})

test_that("admControl(): covMethod 'none' accepted", {
  ctl <- admControl(covMethod = "none")
  expect_equal(ctl$covMethod, "none")
})

test_that("admControl(): grad = 'cfd' accepted", {
  ctl <- admControl(grad = "cfd")
  expect_equal(ctl$grad, "cfd")
  expect_equal(ctl$algorithm, "NLOPT_LD_LBFGS")
})

test_that("admControl(): addProp 'combined1' accepted", {
  ctl <- admControl(addProp = "combined1")
  expect_equal(ctl$addProp, "combined1")
})

test_that("admControl(): print stored as integer", {
  ctl <- admControl(print = 5)
  expect_type(ctl$print, "integer")
  expect_equal(ctl$print, 5L)
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

test_that("adirmcControl(): grad != 'none' defaults to LBFGS", {
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

test_that("adirmcControl(): kappa_method 'linearized' accepted", {
  ctl <- adirmcControl(kappa_method = "linearized")
  expect_equal(ctl$kappa_method, "linearized")
})

test_that("adirmcControl(): max_worse stored as integer", {
  ctl <- adirmcControl(max_worse = 10)
  expect_type(ctl$max_worse, "integer")
  expect_equal(ctl$max_worse, 10L)
})

test_that("adirmcControl(): convcrit stored", {
  ctl <- adirmcControl(convcrit = 1e-4)
  expect_equal(ctl$convcrit, 1e-4)
})

test_that("adirmcControl(): n_restarts and workers stored", {
  ctl <- adirmcControl(n_restarts = 3L, workers = 2L)
  expect_equal(ctl$n_restarts, 3L)
  expect_equal(ctl$workers, 2L)
})

test_that("adirmcControl(): omega_expansion stored", {
  ctl <- adirmcControl(omega_expansion = 2.0)
  expect_equal(ctl$omega_expansion, 2.0)
})

test_that("adirmcControl(): n_sim = 0 errors", {
  expect_error(adirmcControl(n_sim = 0))
})

test_that("adirmcControl(): outer_iter = 0 errors", {
  expect_error(adirmcControl(outer_iter = 0))
})

test_that("adirmcControl(): max_worse = 0 errors", {
  expect_error(adirmcControl(max_worse = 0))
})

test_that("adirmcControl(): workers = 0 errors", {
  expect_error(adirmcControl(workers = 0))
})

test_that("adirmcControl(): n_restarts = 0 errors", {
  expect_error(adirmcControl(n_restarts = 0))
})

test_that("adirmcControl(): cores = 0 errors", {
  expect_error(adirmcControl(cores = 0))
})

test_that("adirmcControl(): covMethod 'none' accepted", {
  ctl <- adirmcControl(covMethod = "none")
  expect_equal(ctl$covMethod, "none")
})

test_that("adirmcControl(): returnAdmr stored", {
  ctl <- adirmcControl(returnAdmr = TRUE)
  expect_true(ctl$returnAdmr)
})

test_that("adirmcControl(): sampling stored", {
  ctl <- adirmcControl(sampling = "halton")
  expect_equal(ctl$sampling, "halton")
})

test_that("adirmcControl(): grad = 'none' defaults to BOBYQA", {
  ctl <- adirmcControl(grad = "none")
  expect_equal(ctl$algorithm, "NLOPT_LN_BOBYQA")
})

test_that("adirmcControl(): MMA selectable with the default gradient", {
  ctl <- adirmcControl(algorithm = "NLOPT_LD_MMA")
  expect_equal(ctl$algorithm, "NLOPT_LD_MMA")
  expect_equal(ctl$grad, "analytical")
})

test_that("adirmcControl(): gradient algorithm + grad 'none' falls back to BOBYQA", {
  ctl <- suppressMessages(adirmcControl(algorithm = "NLOPT_LD_MMA", grad = "none"))
  expect_equal(ctl$algorithm, "NLOPT_LN_BOBYQA")
  expect_equal(ctl$grad, "none")
})

test_that("adirmcControl(): derivative-free algorithm drops the gradient", {
  ctl <- suppressMessages(
    adirmcControl(algorithm = "NLOPT_LN_SBPLX", grad = "analytical"))
  expect_equal(ctl$algorithm, "NLOPT_LN_SBPLX")
  expect_equal(ctl$grad, "none")
})

test_that("adirmcControl(): invalid algorithm errors", {
  expect_error(adirmcControl(algorithm = "NLOPT_LD_NOTREAL"),
               regexp = "not a valid nloptr algorithm")
})

# ---- admc S3 dispatch helpers ------------------------------------------------

test_that("getValidNlmixrCtl.admc: passes through admControl unchanged", {
  ctl <- admControl(maxeval = 77L)
  out <- admixr2:::getValidNlmixrCtl.admc(ctl)
  expect_s3_class(out, "admControl")
  expect_equal(out$maxeval, 77L)
})

test_that("getValidNlmixrCtl.admc: returns default for NULL (wrapped in list)", {
  out <- admixr2:::getValidNlmixrCtl.admc(list(NULL))
  expect_s3_class(out, "admControl")
})

test_that("getValidNlmixrCtl.admc: extracts nested admControl", {
  ctl <- admControl(maxeval = 55L)
  out <- admixr2:::getValidNlmixrCtl.admc(list(ctl))
  expect_s3_class(out, "admControl")
  expect_equal(out$maxeval, 55L)
})

test_that("getValidNlmixrCtl.admc: builds from plain list with studies key", {
  plain <- list(list(studies = list(), n_sim = 200L))
  out   <- admixr2:::getValidNlmixrCtl.admc(plain)
  expect_s3_class(out, "admControl")
})

test_that("nmObjHandleControlObject.admControl: assigns control to env", {
  ctl <- admControl(maxeval = 33L)
  e   <- new.env(parent = emptyenv())
  admixr2:::nmObjHandleControlObject.admControl(ctl, e)
  expect_true(exists("admControl", envir = e, inherits = FALSE))
  expect_equal(get("admControl", envir = e)$maxeval, 33L)
})

test_that("nmObjGetControl.admc: retrieves control from env", {
  ctl <- admControl(maxeval = 44L)
  e   <- new.env(parent = emptyenv())
  e$admControl <- ctl
  out <- admixr2:::nmObjGetControl.admc(list(e))
  expect_s3_class(out, "admControl")
  expect_equal(out$maxeval, 44L)
})

test_that("nmObjGetControl.admc: falls back to 'control' slot", {
  ctl <- admControl(maxeval = 88L)
  e   <- new.env(parent = emptyenv())
  e$control <- ctl
  out <- admixr2:::nmObjGetControl.admc(list(e))
  expect_s3_class(out, "admControl")
  expect_equal(out$maxeval, 88L)
})

test_that("nmObjGetControl.admc: errors when no control found", {
  e <- new.env(parent = emptyenv())
  expect_error(
    admixr2:::nmObjGetControl.admc(list(e)),
    regexp = "cannot find admc control"
  )
})

# ---- adirmc S3 dispatch helpers -----------------------------------------------

test_that("getValidNlmixrCtl.adirmc: passes through adirmcControl unchanged", {
  ctl <- adirmcControl(outer_iter = 20L)
  out <- admixr2:::getValidNlmixrCtl.adirmc(ctl)
  expect_s3_class(out, "adirmcControl")
  expect_equal(out$outer_iter, 20L)
})

test_that("getValidNlmixrCtl.adirmc: extracts nested adirmcControl", {
  ctl <- adirmcControl(outer_iter = 25L)
  out <- admixr2:::getValidNlmixrCtl.adirmc(list(ctl))
  expect_s3_class(out, "adirmcControl")
  expect_equal(out$outer_iter, 25L)
})

test_that("getValidNlmixrCtl.adirmc: converts admControl to adirmcControl", {
  # grad must be compatible with adirmcControl ("analytical"/"none"/"fd", not "sens")
  amc <- admControl(n_sim = 300L, grad = "none")
  out <- admixr2:::getValidNlmixrCtl.adirmc(amc)
  expect_s3_class(out, "adirmcControl")
})

test_that("getValidNlmixrCtl.adirmc: converts nested admControl to adirmcControl", {
  amc <- admControl(n_sim = 300L, grad = "none")
  out <- admixr2:::getValidNlmixrCtl.adirmc(list(amc))
  expect_s3_class(out, "adirmcControl")
})

test_that("getValidNlmixrCtl.adirmc: builds from list with studies key", {
  plain <- list(list(studies = list(), n_sim = 100L))
  out   <- admixr2:::getValidNlmixrCtl.adirmc(plain)
  expect_s3_class(out, "adirmcControl")
})

test_that("getValidNlmixrCtl.adirmc: errors for unrecognised plain list", {
  bad <- list(list(foo = 1))
  expect_error(
    admixr2:::getValidNlmixrCtl.adirmc(bad),
    regexp = "adirmcControl"
  )
})

test_that("nmObjHandleControlObject.adirmcControl: assigns control to env", {
  ctl <- adirmcControl(outer_iter = 15L)
  e   <- new.env(parent = emptyenv())
  admixr2:::nmObjHandleControlObject.adirmcControl(ctl, e)
  expect_true(exists("adirmcControl", envir = e, inherits = FALSE))
  expect_equal(get("adirmcControl", envir = e)$outer_iter, 15L)
})

test_that("nmObjGetControl.adirmc: retrieves control from env", {
  ctl <- adirmcControl(outer_iter = 30L)
  e   <- new.env(parent = emptyenv())
  e$adirmcControl <- ctl
  out <- admixr2:::nmObjGetControl.adirmc(list(e))
  expect_s3_class(out, "adirmcControl")
  expect_equal(out$outer_iter, 30L)
})

test_that("nmObjGetControl.adirmc: falls back to 'control' slot", {
  ctl <- adirmcControl(outer_iter = 20L)
  e   <- new.env(parent = emptyenv())
  e$control <- ctl
  out <- admixr2:::nmObjGetControl.adirmc(list(e))
  expect_s3_class(out, "adirmcControl")
  expect_equal(out$outer_iter, 20L)
})

test_that("nmObjGetControl.adirmc: errors when no control found", {
  e <- new.env(parent = emptyenv())
  expect_error(
    admixr2:::nmObjGetControl.adirmc(list(e)),
    regexp = "cannot find irmc control"
  )
})

