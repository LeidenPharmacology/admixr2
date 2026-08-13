# Shi (2021) step selection, end to end.
#
# The fixed `pmax(abs(p), 0.1) * h` scale is one guess about how much noise the
# objective carries, applied identically to every parameter -- and it is the
# guess behind the "Hessian not positive definite ... try increasing
# cov_h_outer" warning. Shi21 probes the objective and balances truncation
# against noise per parameter, and since 0.4.1 it is the unconditional default
# with that fixed scale demoted to a fallback. There is no flag to toggle, so
# what these tests assert is that the measured steps REACH each site and that
# the fits they produce are sane: SPD covariance, finite SEs, convergence.
skip_on_cran()
skip_if_not_installed("rxode2")
skip_if_not_installed("nlmixr2est")

.shi_fit <- function(maxeval = 40L, ...) {
  env <- .int_grad_setup()
  suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
    one_cmt_fn, admData(), est = "adfo",
    control = adfoControl(studies = env$studies, maxeval = maxeval, seed = 1L,
                          covMethod = "r", ...))))
}

test_that("the measured covariance step gives an SPD covariance and finite SEs", {
  f1 <- .shi_fit()
  expect_s3_class(f1, "admFit")
  expect_true(is.finite(f1$objective))
  skip_if(is.null(f1$cov), "covariance unavailable")
  expect_true(all(eigen(f1$cov, symmetric = TRUE, only.values = TRUE)$values > 0))
  se1 <- sqrt(diag(f1$cov))
  expect_true(all(is.finite(se1) & se1 > 0))
})

test_that("the step measurement is deterministic", {
  # Shi21 probes with fixed offsets and ECnoise with a fixed grid -- no RNG
  # anywhere -- so two identical fits must agree bit-for-bit. If this ever fails
  # something stochastic has entered the step choice, which would make every
  # covariance irreproducible.
  a <- .shi_fit()
  b <- .shi_fit()
  skip_if(is.null(a$cov) || is.null(b$cov), "covariance unavailable")
  expect_identical(a$cov, b$cov)
  expect_identical(a$objective, b$objective)
})

test_that("grad = 'fd' converges to the same optimum as the analytic gradient", {
  # grad = "fd" finite-differences EVERY parameter through the whole NLL, so the
  # gradient step probe does run here. The FD fit is a different iterate
  # sequence from the analytic one, but it must reach the same place -- that is
  # the only thing that makes a step choice legitimate. maxeval is deliberately
  # generous so both are compared AT convergence rather than mid-descent.
  ga <- .shi_fit(grad = "analytical", maxeval = 300L)
  gf <- .shi_fit(grad = "fd",         maxeval = 300L)

  expect_s3_class(gf, "admFit")
  expect_true(is.finite(gf$objective))
  expect_equal(gf$objective, ga$objective, tolerance = 1e-2)
  expect_equal(gf$parFixedDf$Estimate, ga$parFixedDf$Estimate, tolerance = 5e-2)
})

test_that("the measured step vector really is per-parameter and reaches the sites", {
  # .admGH must index a vector and pass a scalar straight through, which is what
  # lets a measured vector reach every FD site with no signature change.
  # Deliberately NOT a quadratic. A quadratic has no third derivative, so the
  # probe cannot measure a step and correctly falls back for every parameter --
  # which would make the per-parameter assertions below vacuously true.
  f <- function(p) sum(exp(p))
  p <- c(0.5, 1.5, 2.5)

  h <- admixr2:::.admShi21GradH(f, p, c(1L, 3L), 1e-4, scaled = TRUE)
  expect_length(h, 3L)
  expect_identical(h[2], 1e-4)                 # not requested -> untouched
  expect_true(all(h[c(1, 3)] != 1e-4))         # requested -> actually measured
  expect_true(all(is.finite(h) & h > 0))
  # `scaled` is the whole contract: the stored value is the measured ABSOLUTE
  # step divided by whatever the call site will multiply it back by, so that
  # `pmax(abs(p), 0.1) * h` reproduces the measurement exactly. Assert that
  # directly. (It is NOT that the entries differ in proportion to |p| -- that
  # held for Gill83, whose step keys off f'', but Shi21's keys off f''', so the
  # ratio carries the curvature too and the old assertion was measuring nothing.)
  h_unscaled <- admixr2:::.admShi21GradH(f, p, c(1L, 3L), 1e-4, scaled = FALSE)
  expect_equal(h[c(1, 3)] * pmax(abs(p[c(1, 3)]), 0.1), h_unscaled[c(1, 3)],
               tolerance = 1e-8)

  expect_identical(admixr2:::.admGH(1e-4, c(1L, 3L)), 1e-4)
  expect_identical(admixr2:::.admGH(c(1, 2, 3, 4), c(1L, 3L)), c(1, 3))
  # An empty index set must skip the probe entirely and change nothing.
  expect_true(all(admixr2:::.admShi21GradH(f, p, integer(0), 1e-4) == 1e-4))
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

test_that("the measured step reaches adirmc's inner FD gradient, once", {
  # adirmc's inner FD used a hard-coded 1e-6 and ignored grad_h entirely, so this
  # is the only coverage of the path that now measures its step.
  env <- .int_grad_setup()
  # The counter lives in an env spliced into the tracer expression: a tracer runs
  # in the TRACED function's frame, so `fired <<- fired + 1L` cannot see a local
  # of this test and errors with "object 'fired' not found".
  cnt <- new.env(parent = emptyenv()); cnt$n <- 0L
  # assign()/get(), not `$<-`: bquote splices the environment in as a literal
  # object, and `<env>$n <- ...` is then "target of assignment expands to
  # non-language object".
  trace(admixr2:::.admShi21GradH,
        tracer = bquote(assign("n", get("n", envir = .(cnt)) + 1L, envir = .(cnt))),
        print = FALSE)
  on.exit(untrace(admixr2:::.admShi21GradH), add = TRUE)

  g1 <- suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
    one_cmt_fn, admData(), est = "adirmc",
    control = adirmcControl(studies = env$studies, seed = 1L, grad = "fd",
                            n_sim = 300L, phases = c(1, 0.5), outer_iter = 2L,
                            maxeval = 20L, covMethod = "none"))))

  # Fired -- but ONCE, not once per outer iteration. Proposals are redrawn every
  # iteration, so re-probing would cost 10 x n_par extra inner NLL evaluations
  # each time; the measurement is taken on the first and reused.
  expect_identical(cnt$n, 1L)
  expect_true(is.finite(g1$objective))
})

test_that("the measured step reaches admc's gradient for a finite-differenced theta", {
  # admc probes only parameters it actually steps in PARAMETER space. With the
  # sensitivity model supplying a column for every theta -- the default -- that
  # set is empty and the probe is skipped entirely, which is why no other test
  # here reaches this branch. grad = "fd" turns the sens model off, and
  # one_cmt_kappa_fn has a theta with no mu-referencing eta, so the set is not
  # empty and the measurement runs.
  env <- .int_grad_setup()
  cnt <- new.env(parent = emptyenv()); cnt$n <- 0L
  trace(admixr2:::.admShi21GradH,
        tracer = bquote(assign("n", get("n", envir = .(cnt)) + 1L, envir = .(cnt))),
        print = FALSE)
  on.exit(untrace(admixr2:::.admShi21GradH), add = TRUE)

  b <- suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
    one_cmt_kappa_fn, admData(), est = "admc",
    control = admControl(studies = env$studies, seed = 1L, grad = "fd",
                         n_sim = 300L, maxeval = 15L, covMethod = "none"))))

  expect_gt(cnt$n, 0L)
  expect_true(is.finite(b$objective))
})
