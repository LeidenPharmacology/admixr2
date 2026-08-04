# Gill (1983) step selection, end to end.
#
# The fixed `pmax(abs(p), 0.1) * h` scale is one guess about how much noise the
# objective carries, applied identically to every parameter -- and it is the
# guess behind the "Hessian not positive definite ... try increasing
# cov_h_outer" warning. Gill83 probes the objective and balances condition
# against truncation error per parameter. What must hold: the covariance is
# still SPD, the standard errors agree with the fixed-step ones to within the
# accuracy either can claim, and `gill = FALSE` changes nothing at all.
skip_on_cran()
skip_if_not_installed("rxode2")
skip_if_not_installed("nlmixr2est")

.gill_fit <- function(gill, maxeval = 40L, ...) {
  env <- .int_grad_setup()
  suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
    one_cmt_fn, admData(), est = "adfo",
    control = adfoControl(studies = env$studies, maxeval = maxeval, seed = 1L,
                          covMethod = "r", gill = gill, ...))))
}

test_that("gill = TRUE gives an SPD covariance and SEs consistent with the fixed step", {
  f0 <- .gill_fit(FALSE)
  f1 <- .gill_fit(TRUE)

  expect_s3_class(f1, "admFit")
  # This model builds the order-2 sensitivity block, so NO structural theta is
  # finite-differenced and .fd_idx is empty: the gradient probe never runs and
  # gill reaches only the covariance step. The objective is therefore still
  # bit-identical -- but that is a property of THIS model, not of the flag, which
  # is why the fd/cfd case below is tested separately.
  expect_identical(f1$objective, f0$objective)

  skip_if(is.null(f0$cov) || is.null(f1$cov), "covariance unavailable")
  expect_true(all(eigen(f1$cov, symmetric = TRUE, only.values = TRUE)$values > 0))
  expect_identical(dimnames(f1$cov), dimnames(f0$cov))

  se0 <- sqrt(diag(f0$cov)); se1 <- sqrt(diag(f1$cov))
  expect_true(all(is.finite(se1) & se1 > 0))
  # Two different-but-valid step choices for the same Hessian: they should agree
  # to a few percent, not to machine precision. A tight tolerance here would be
  # asserting that the step does not matter, which is the opposite of the point.
  expect_lt(max(abs(se1 - se0) / se0), 0.15)
})

test_that("gill = FALSE leaves the covariance exactly as it was", {
  # The default must not move anyone's standard errors.
  a <- .gill_fit(FALSE)
  b <- .gill_fit(FALSE)
  skip_if(is.null(a$cov) || is.null(b$cov), "covariance unavailable")
  expect_equal(a$cov, b$cov)
})

test_that("gill reaches the optimizer's FD gradient under grad = 'fd'", {
  # grad = "fd" finite-differences EVERY parameter through the whole NLL, so the
  # gradient probe does run here. The fit is allowed to move -- a different step
  # is a different gradient, and a different iterate sequence -- but it must not
  # land somewhere WORSE, which is the only thing that makes a step choice
  # legitimate.
  #
  # maxeval is deliberately generous. At maxeval = 40 the two runs stop at
  # different points along their descent (measured 713.3 vs 715.6, the Gill run
  # ahead), and asserting they agree there would be asserting that the step makes
  # no difference -- the opposite of the point.
  g0 <- .gill_fit(FALSE, grad = "fd", maxeval = 300L)
  g1 <- .gill_fit(TRUE,  grad = "fd", maxeval = 300L)

  expect_s3_class(g1, "admFit")
  expect_true(is.finite(g1$objective))
  # Converged to the same optimum, by differently-stepped gradients.
  expect_equal(g1$objective, g0$objective, tolerance = 1e-2)
  th0 <- g0$parFixedDf$Estimate; th1 <- g1$parFixedDf$Estimate
  expect_equal(th1, th0, tolerance = 5e-2)
})

test_that("the Gill step vector really is per-parameter and reaches the sites", {
  # Two regressions in one: (a) the acceptance rule must not reject the ordinary
  # "High Grad Error" return -- doing so silently made gill a no-op everywhere;
  # (b) .admGH must index a vector and pass a scalar straight through, which is
  # what lets a vector reach every FD site with no signature change.
  f <- function(p) sum((p - c(1, 2, 3))^2)
  p <- c(0.5, 1.5, 2.5)

  h <- admixr2:::.admGillGradH(f, p, c(1L, 3L), 1e-4, scaled = TRUE)
  expect_length(h, 3L)
  expect_identical(h[2], 1e-4)                 # not requested -> untouched
  expect_true(all(h[c(1, 3)] != 1e-4))         # requested -> actually measured
  expect_true(all(is.finite(h) & h > 0))
  # scaled = TRUE means the site's own pmax(abs(p), 0.1) * h reproduces Gill's
  # absolute step; so the two requested parameters differ in proportion to |p|.
  expect_equal(h[1] / h[3], pmax(abs(p[3]), 0.1) / pmax(abs(p[1]), 0.1),
               tolerance = 1e-6)

  expect_identical(admixr2:::.admGH(1e-4, c(1L, 3L)), 1e-4)
  expect_identical(admixr2:::.admGH(c(1, 2, 3, 4), c(1L, 3L)), c(1, 3))
  # An empty index set must skip the probe entirely and change nothing.
  expect_true(all(admixr2:::.admGillGradH(f, p, integer(0), 1e-4) == 1e-4))
})

test_that(".admGrad with a uniform step vector reproduces the scalar exactly", {
  # The real risk in teaching .admGrad/.admGradBatch to take a per-parameter
  # vector is not the arithmetic, it is the INDEXING: every site had to be keyed
  # to what it actually differences (unpaired_k[bi] for a struct theta, k_s in
  # the joint branch, .admGH0() in eta space). Feed it a vector whose entries are
  # all the old constant: any site that indexes the wrong way, or that reaches an
  # eta perturbation through .admGH() instead of .admGH0(), still returns the
  # same number here -- but a site that dropped its division, or recycles, does
  # not. Bit-identical is the bar; a uniform vector IS the scalar.
  env <- .int_grad_setup()
  h   <- 1e-3
  p0  <- env$vec$p0

  g_scalar <- admixr2:::.admGrad(p0, env$pinfo, env$studies, env$z_list, env$rxMod,
                                 env$output_var, env$params_list, cores = 1L,
                                 h = h, sensModel = env$sensModel)
  h_vec <- rep(h, length(p0))
  attr(h_vec, "h0") <- h
  g_vector <- admixr2:::.admGrad(p0, env$pinfo, env$studies, env$z_list, env$rxMod,
                                 env$output_var, env$params_list, cores = 1L,
                                 h = h_vec, sensModel = env$sensModel)
  expect_identical(g_vector, g_scalar)

  # ... and with the sens columns hidden, which is the path that actually
  # finite-differences the unpaired thetas and the etas.
  sm <- env$sensModel
  if (!is.null(sm)) sm$theta_sens_cols <- NULL
  f_scalar <- admixr2:::.admGrad(p0, env$pinfo, env$studies, env$z_list, env$rxMod,
                                 env$output_var, env$params_list, cores = 1L,
                                 h = h, sensModel = sm)
  f_vector <- admixr2:::.admGrad(p0, env$pinfo, env$studies, env$z_list, env$rxMod,
                                 env$output_var, env$params_list, cores = 1L,
                                 h = h_vec, sensModel = sm)
  expect_identical(f_vector, f_scalar)
})

test_that(".admGrad rejects a step vector of the wrong length", {
  # .admGrad accepts one step or one per parameter, and NOTHING else. A vector of
  # any other length would recycle into the eta-space perturbations
  # (`eta_hi[, j] + h` over n_sim rows) and the n_sim x n_t divisions without a
  # warning, handing every draw a different perturbation and returning a
  # plausible, wrong gradient. The guard fires before any solve, so no fixture.
  wrong <- tryCatch({
    admixr2:::.admGrad(rep(0, 4), list(), list(), list(), NULL, "cp",
                       list(), 1L, h = rep(1e-4, 3))   # 3 steps, 4 parameters
    ""
  }, error = function(e) conditionMessage(e))
  expect_match(wrong, "one step or one per parameter")

  # The two accepted shapes must get PAST the guard. Both still fail afterwards
  # on the empty pinfo, so assert on WHICH failure rather than on success.
  for (h in list(1e-4, rep(1e-4, 4))) {
    msg <- tryCatch({
      admixr2:::.admGrad(rep(0, 4), list(), list(), list(), NULL, "cp",
                         list(), 1L, h = h)
      ""
    }, error = function(e) conditionMessage(e))
    expect_false(grepl("one step or one per parameter", msg, fixed = TRUE))
  }
})

test_that("gill reaches adirmc's inner FD gradient", {
  # adirmc's inner FD used a hard-coded 1e-6 and ignored grad_h entirely, so this
  # is the only coverage of the path that now honours it. The step is measured on
  # the FIRST outer iteration and reused, because proposals are redrawn each
  # iteration -- so the probe must fire exactly once, not per iteration.
  env <- .int_grad_setup()
  # The counter lives in an env spliced into the tracer expression: a tracer runs
  # in the TRACED function's frame, so `fired <<- fired + 1L` cannot see a local
  # of this test and errors with "object 'fired' not found".
  cnt <- new.env(parent = emptyenv()); cnt$n <- 0L
  # assign()/get(), not `$<-`: bquote splices the environment in as a literal
  # object, and `<env>$n <- ...` is then "target of assignment expands to
  # non-language object".
  trace(admixr2:::.admGillGradH,
        tracer = bquote(assign("n", get("n", envir = .(cnt)) + 1L, envir = .(cnt))),
        print = FALSE)
  on.exit(untrace(admixr2:::.admGillGradH), add = TRUE)

  run <- function(gill) suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
    one_cmt_fn, admData(), est = "adirmc",
    control = adirmcControl(studies = env$studies, seed = 1L, grad = "fd",
                            n_sim = 300L, phases = c(1, 0.5), outer_iter = 2L,
                            maxeval = 20L, covMethod = "none", gill = gill))))

  g0 <- run(FALSE)
  expect_identical(cnt$n, 0L)          # nothing probed when the flag is off
  g1 <- run(TRUE)
  expect_gt(cnt$n, 0L)                 # ... and probed when it is on

  expect_true(is.finite(g0$objective))
  expect_true(is.finite(g1$objective))
  expect_equal(g1$objective, g0$objective, tolerance = 1e-2)
})

test_that("gill reaches admc's gradient when a struct theta is finite-differenced", {
  # admc probes only parameters it actually steps in PARAMETER space. With the
  # sensitivity model supplying a column for every theta -- the default -- that
  # set is empty and the probe is skipped, which is why no other test here
  # reaches this branch. grad = "fd" turns the sens model off, and
  # one_cmt_kappa_fn has a theta with no mu-referencing eta, so the set is not
  # empty and the measurement runs.
  env <- .int_grad_setup()
  cnt <- new.env(parent = emptyenv()); cnt$n <- 0L
  trace(admixr2:::.admGillGradH,
        tracer = bquote(assign("n", get("n", envir = .(cnt)) + 1L, envir = .(cnt))),
        print = FALSE)
  on.exit(untrace(admixr2:::.admGillGradH), add = TRUE)

  run <- function(gill) suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
    one_cmt_kappa_fn, admData(), est = "admc",
    control = admControl(studies = env$studies, seed = 1L, grad = "fd",
                         n_sim = 300L, maxeval = 15L, covMethod = "none",
                         gill = gill))))

  a <- run(FALSE)
  expect_identical(cnt$n, 0L)
  b <- run(TRUE)
  expect_gt(cnt$n, 0L)

  expect_true(is.finite(a$objective))
  expect_true(is.finite(b$objective))
  # Same objective surface, differently-stepped gradients at a small maxeval:
  # they need only be in the same place, not identical.
  expect_equal(b$objective, a$objective, tolerance = 5e-2)
})
