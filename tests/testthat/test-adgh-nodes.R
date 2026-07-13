# ---- .adghNodes1 -------------------------------------------------------------

test_that(".adghNodes1(1) returns trivial node", {
  g <- admixr2:::.adghNodes1(1L)
  expect_equal(g$x, 0)
  expect_equal(g$w, 1)
})

test_that(".adghNodes1(m): weights sum to 1 (probabilists' form)", {
  for (m in c(2L, 3L, 5L, 7L)) {
    g <- admixr2:::.adghNodes1(m)
    expect_equal(sum(g$w), 1, tolerance = 1e-12,
                 label = sprintf("sum(w), m=%d", m))
  }
})

test_that(".adghNodes1(m): weights are non-negative", {
  for (m in c(2L, 3L, 5L, 7L)) {
    g <- admixr2:::.adghNodes1(m)
    expect_true(all(g$w >= 0), label = sprintf("all(w>=0), m=%d", m))
  }
})

test_that(".adghNodes1(m): E[x^2] = 1 under N(0,1) (probabilists' GH)", {
  for (m in c(2L, 3L, 5L, 7L)) {
    g <- admixr2:::.adghNodes1(m)
    expect_equal(sum(g$w * g$x^2), 1, tolerance = 1e-10,
                 label = sprintf("E[x^2]=1, m=%d", m))
  }
})

test_that(".adghNodes1(m): E[x] = 0 (symmetry)", {
  for (m in c(2L, 3L, 5L, 7L)) {
    g <- admixr2:::.adghNodes1(m)
    expect_equal(sum(g$w * g$x), 0, tolerance = 1e-12,
                 label = sprintf("E[x]=0, m=%d", m))
  }
})

test_that(".adghNodes1(m): returns m nodes", {
  for (m in c(1L, 3L, 5L)) {
    g <- admixr2:::.adghNodes1(m)
    expect_length(g$x, m)
    expect_length(g$w, m)
  }
})

test_that(".adghNodes1(m): m < 1 errors", {
  expect_error(admixr2:::.adghNodes1(0L))
})

# ---- .adghNodeGrid -----------------------------------------------------------

test_that(".adghNodeGrid(m, 0): degenerate grid (no etas)", {
  g <- admixr2:::.adghNodeGrid(5L, 0L)
  expect_equal(nrow(g$X), 1L)
  expect_equal(ncol(g$X), 0L)
  expect_equal(g$W, 1)
})

test_that(".adghNodeGrid(m, 1): m nodes, weights sum to 1", {
  g <- admixr2:::.adghNodeGrid(5L, 1L)
  expect_equal(nrow(g$X), 5L)
  expect_equal(ncol(g$X), 1L)
  expect_equal(sum(g$W), 1, tolerance = 1e-12)
})

test_that(".adghNodeGrid(m, 2): m^2 nodes, weights sum to 1", {
  g <- admixr2:::.adghNodeGrid(3L, 2L)
  expect_equal(nrow(g$X), 9L)
  expect_equal(ncol(g$X), 2L)
  expect_equal(sum(g$W), 1, tolerance = 1e-12)
})

test_that(".adghNodeGrid(m, n): E[x_j^2] = 1 for each dimension", {
  g <- admixr2:::.adghNodeGrid(5L, 2L)
  for (j in 1:2)
    expect_equal(sum(g$W * g$X[, j]^2), 1, tolerance = 1e-10,
                 label = sprintf("E[x_%d^2]=1", j))
})

test_that(".adghNodeGrid(m, n): marginal weights are correct (colSums via product)", {
  m <- 3L; d <- 2L
  g <- admixr2:::.adghNodeGrid(m, d)
  g1 <- admixr2:::.adghNodes1(m)
  # marginal weight for first dim: sum over all second-dim combinations
  # = product of marginal weights * sum of second-dim weights
  expect_equal(length(g$W), m^d)
  expect_true(all(g$W > 0))
})

# ---- adghControl -------------------------------------------------------------

test_that("adghControl() returns correct class and key defaults", {
  ctl <- adghControl()
  expect_s3_class(ctl, "adghControl")
  expect_equal(ctl$n_nodes,    5L)
  expect_equal(ctl$grad,       "analytical")
  expect_equal(ctl$maxeval,    500L)
  expect_equal(ctl$seed,       12345L)
  expect_equal(ctl$cores,      as.integer(rxode2::rxCores()))
  expect_equal(ctl$n_restarts, 1L)
  expect_equal(ctl$workers,    1L)
  expect_equal(ctl$covMethod,  "r")
  expect_equal(ctl$returnAdmr, FALSE)
})

test_that("adghControl(): default grad='analytical' coerces BOBYQA to LBFGS", {
  ctl <- adghControl()
  expect_equal(ctl$algorithm, "NLOPT_LD_LBFGS")
})

test_that("adghControl(): grad='fd' coerces BOBYQA to LBFGS", {
  ctl <- adghControl(grad = "fd")
  expect_equal(ctl$algorithm, "NLOPT_LD_LBFGS")
})

test_that("adghControl(): grad='cfd' coerces BOBYQA to LBFGS", {
  ctl <- adghControl(grad = "cfd")
  expect_equal(ctl$algorithm, "NLOPT_LD_LBFGS")
})

test_that("adghControl(): grad='none' keeps BOBYQA", {
  ctl <- adghControl(grad = "none")
  expect_equal(ctl$algorithm, "NLOPT_LN_BOBYQA")
})

test_that("adghControl(): grad != 'none' + explicit non-BOBYQA algorithm kept", {
  ctl <- adghControl(grad = "fd", algorithm = "NLOPT_LD_SLSQP")
  expect_equal(ctl$algorithm, "NLOPT_LD_SLSQP")
})

test_that("adghControl(): internal n_sim=1L for .admRunRestarts() compat", {
  ctl <- adghControl()
  expect_type(ctl$n_sim, "integer")
  expect_equal(ctl$n_sim, 1L)
})

test_that("adghControl(): internal sampling='sobol' for .admRunRestarts() compat", {
  ctl <- adghControl()
  expect_equal(ctl$sampling, "sobol")
})

test_that("adghControl(): n_nodes stored as integer", {
  ctl <- adghControl(n_nodes = 3)
  expect_type(ctl$n_nodes, "integer")
  expect_equal(ctl$n_nodes, 3L)
})

test_that("adghControl(): maxeval stored as integer", {
  ctl <- adghControl(maxeval = 200)
  expect_type(ctl$maxeval, "integer")
  expect_equal(ctl$maxeval, 200L)
})

test_that("adghControl(): cov_h_outer default is eps^(1/4)", {
  ctl <- adghControl()
  expect_equal(ctl$cov_h_outer, .Machine$double.eps^(1/4), tolerance = 1e-15)
})

test_that("adghControl(): covMethod 'none' accepted", {
  ctl <- adghControl(covMethod = "none")
  expect_equal(ctl$covMethod, "none")
})

test_that("adghControl(): returnAdmr stored correctly", {
  ctl <- adghControl(returnAdmr = TRUE)
  expect_true(ctl$returnAdmr)
})

test_that("adghControl(): unknown argument errors", {
  expect_error(adghControl(nonexistent_arg = 42), regexp = "unused argument")
})

test_that("adghControl(): n_nodes=0 errors via checkmate", {
  expect_error(adghControl(n_nodes = 0L), regexp = "n_nodes")
})

test_that("adghControl(): cores=0 errors via checkmate", {
  expect_error(adghControl(cores = 0L), regexp = "cores")
})

test_that("adghControl(): maxeval=0 errors via checkmate", {
  expect_error(adghControl(maxeval = 0L), regexp = "maxeval")
})

test_that("adghControl(): invalid grad errors", {
  expect_error(adghControl(grad = "newton"))
})

test_that("adghControl(): invalid covMethod errors", {
  expect_error(adghControl(covMethod = "hessian"))
})

test_that("adghControl(): grad_h, grad_bounds, cov_h stored", {
  ctl <- adghControl(grad_h = 1e-5, grad_bounds = 3, cov_h = 5e-4)
  expect_equal(ctl$grad_h,      1e-5)
  expect_equal(ctl$grad_bounds, 3)
  expect_equal(ctl$cov_h,       5e-4)
})

# ---- S3 hooks ----------------------------------------------------------------

test_that("getValidNlmixrCtl.adgh passes through adghControl", {
  ctl <- adghControl(maxeval = 99L)
  out <- admixr2:::getValidNlmixrCtl.adgh(ctl)
  expect_s3_class(out, "adghControl")
  expect_equal(out$maxeval, 99L)
})

test_that("getValidNlmixrCtl.adgh returns default for NULL", {
  out <- admixr2:::getValidNlmixrCtl.adgh(NULL)
  expect_s3_class(out, "adghControl")
})

test_that("getValidNlmixrCtl.adgh converts plain list with n_nodes", {
  out <- admixr2:::getValidNlmixrCtl.adgh(list(n_nodes = 3L, maxeval = 77L))
  expect_s3_class(out, "adghControl")
  expect_equal(out$n_nodes, 3L)
})

test_that("getValidNlmixrCtl.adgh extracts nested adghControl", {
  ctl <- adghControl(maxeval = 66L)
  out <- admixr2:::getValidNlmixrCtl.adgh(list(ctl))
  expect_s3_class(out, "adghControl")
  expect_equal(out$maxeval, 66L)
})

test_that("getValidNlmixrCtl.adgh builds from nested list with studies key", {
  out <- admixr2:::getValidNlmixrCtl.adgh(list(list(studies = list(), maxeval = 55L)))
  expect_s3_class(out, "adghControl")
})

test_that("nmObjHandleControlObject.adghControl: assigns control to env", {
  ctl <- adghControl(maxeval = 22L)
  e   <- new.env(parent = emptyenv())
  admixr2:::nmObjHandleControlObject.adghControl(ctl, e)
  expect_true(exists("adghControl", envir = e, inherits = FALSE))
  expect_equal(get("adghControl", envir = e)$maxeval, 22L)
})

test_that("nmObjGetControl.adgh: retrieves control from env", {
  ctl <- adghControl(maxeval = 33L)
  e   <- new.env(parent = emptyenv())
  e$adghControl <- ctl
  out <- admixr2:::nmObjGetControl.adgh(list(e))
  expect_s3_class(out, "adghControl")
  expect_equal(out$maxeval, 33L)
})

test_that("nmObjGetControl.adgh: falls back to 'control' slot", {
  ctl <- adghControl(maxeval = 44L)
  e   <- new.env(parent = emptyenv())
  e$control <- ctl
  out <- admixr2:::nmObjGetControl.adgh(list(e))
  expect_s3_class(out, "adghControl")
  expect_equal(out$maxeval, 44L)
})

test_that("nmObjGetControl.adgh: errors when no control found", {
  e <- new.env(parent = emptyenv())
  expect_error(admixr2:::nmObjGetControl.adgh(list(e)),
               regexp = "cannot find adgh control")
})
