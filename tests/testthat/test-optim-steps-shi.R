# Shi (2021) adaptive central-difference intervals, and the noise estimate they
# need. Tier 1 -- analytic functions only, no rxode2.
#
# The oracle test is the important one: nlmixr2est ships this algorithm as
# `shi21CentralWrap`, unexported, so admixr2 reimplements it (see the note in
# optim-steps.R). Two independent implementations of the same published
# procedure should agree, and that is the only cheap way to know the
# reimplementation is faithful rather than merely plausible.

test_that(".admEcNoise recovers a known injected noise level", {
  p <- c(1.6, 3.0, -2.4)
  for (nz in c(1e-12, 1e-9, 1e-6)) {
    set.seed(1)
    fn <- function(x) sum(exp(x)) + nz * (2 * stats::runif(1) - 1)
    e <- admixr2:::.admEcNoise(fn, p, 1L)
    expect_true(is.finite(e))
    # Moré-Wild is an order-of-magnitude estimator; h* scales as eps_f^(1/3), so
    # even a factor of 2 here moves the step by only 26%.
    expect_gt(e / nz, 0.2)
    expect_lt(e / nz, 5)
  }
})

test_that(".admEcNoise on a noiseless function returns the rounding floor", {
  p <- c(1.6, 3.0, -2.4)
  f <- function(x) sum(exp(x))
  e <- admixr2:::.admEcNoise(f, p, 1L)
  expect_true(is.finite(e))
  # |f| * eps is the smallest noise a double-precision evaluation can have.
  expect_lt(e, 1e-12)
  expect_gt(e, 0)
})

test_that(".admShi21Central lands near the theoretical optimum", {
  # h* = (3 eps_f / |f'''|)^(1/3), and for sum(exp(x)) the third derivative in
  # coordinate i is exactly exp(x_i) -- so the target is known in closed form.
  p <- c(1.6, 3.0, -2.4)
  f <- function(x) sum(exp(x))
  for (eps_f in c(1e-4, 1e-7, 1e-11)) {
    for (k in seq_along(p)) {
      r <- admixr2:::.admShi21Central(f, p, k, eps_f)
      hstar <- (3 * eps_f / exp(p[k]))^(1/3)
      expect_gt(r$h / hstar, 0.5)
      expect_lt(r$h / hstar, 2)
      # The derivative tolerance must scale with the noise, not be a constant:
      # at the optimum both error terms are O(eps_f^(2/3)), so eps_f = 1e-4
      # cannot do better than ~1e-3 no matter how good the step is. A flat 1e-3
      # here failed on exactly that row -- while the upstream oracle was LESS
      # accurate on it than this implementation.
      expect_equal(r$gr, exp(p[k]), tolerance = 3 * eps_f^(2/3))
    }
  }
})

test_that(".admShi21Central terminates instead of limit-cycling", {
  # Regression: the first design iterated h toward a fixed point of the h*
  # formula. It could not converge -- at the optimum the third difference is
  # only ~2x its own noise, so any floor high enough to trust it rejects h*
  # itself -- and it burned every iteration, 4 evaluations each, returning the
  # unrefined starting guess. Cap the budget and assert we stay well inside it.
  p <- c(1.6, 3.0, -2.4)
  n <- 0L
  f <- function(x) { n <<- n + 1L; sum(exp(x)) }
  r <- admixr2:::.admShi21Central(f, p, 1L, 1e-7)
  expect_true(is.finite(r$h) && r$h > 0)
  expect_lt(n, 30L)     # maxiter * 4 would be 40 and would mean no convergence
})

test_that(".admShi21Central beats a fixed step on a smooth objective", {
  p <- c(1.6, 3.0, -2.4)
  f <- function(x) sum(exp(x))
  eps <- admixr2:::.admEcNoise(f, p, 1L)
  fwd <- function(h, k) { q <- p; q[k] <- q[k] + h; (f(q) - f(p)) / h }
  for (k in seq_along(p)) {
    shi <- admixr2:::.admShi21Central(f, p, k, eps)$gr
    ref <- exp(p[k])
    expect_lt(abs(shi - ref) / ref, abs(fwd(1e-6, k) - ref) / ref)
  }
})

test_that(".admShi21Steps returns one positive finite step per requested index", {
  p <- c(1.6, 3.0, -2.4)
  f <- function(x) sum(exp(x))
  h <- admixr2:::.admShi21Steps(f, p, idx = c(1L, 3L), fallback = c(1e-5, 1e-5))
  expect_length(h, 2L)
  expect_true(all(is.finite(h)))
  expect_true(all(h > 0))
})

test_that(".admShi21Steps falls back when the noise level cannot be estimated", {
  p <- c(1.6, 3.0)
  # Constant in every direction: no noise information, and no step is better
  # than another.
  flat <- function(x) 1
  expect_warning(
    h <- admixr2:::.admShi21Steps(flat, p, idx = 1:2, fallback = c(1e-5, 2e-5),
                                  .var.name = "unit"),
    "noise level")
  expect_identical(h, c(1e-5, 2e-5))
})

test_that("the reimplementation agrees with nlmixr2est's shi21CentralWrap", {
  skip_if_not_installed("nlmixr2est")
  ora <- tryCatch(get("shi21CentralWrap", asNamespace("nlmixr2est")),
                  error = function(e) NULL)
  skip_if(is.null(ora), "nlmixr2est does not expose shi21CentralWrap")
  p <- c(1.6, 3.0, -2.4)
  f <- function(x) sum(exp(x))
  f0 <- f(p)
  for (eps_f in c(1e-4, 1e-7, 1e-9, 1e-11)) {
    for (k in seq_along(p)) {
      mine <- admixr2:::.admShi21Central(f, p, k, eps_f)
      theirs <- ora(f, p, f0, k, eps_f)
      # Same procedure, different stopping details: agreement to a factor of 2
      # on the interval is what "faithful" means here, not bit-identity. The
      # error curve is flat enough near the optimum that a factor of 2 in h is
      # a factor of ~1.25 in the derivative error.
      expect_gt(mine$h / theirs$h[1], 0.5)
      expect_lt(mine$h / theirs$h[1], 2)
      # Both derivatives carry O(eps_f^(2/3)) error of their own, so they can
      # only be asked to agree to that -- neither is the truth here.
      expect_equal(mine$gr, as.numeric(theirs$gr[1]),
                   tolerance = 10 * eps_f^(2/3))
    }
  }
})

test_that(".admHessSteps uses the SECOND-difference exponent, not the first", {
  # Regression: .admHessSteps first delegated straight to .admShi21Steps, which
  # returns the first-derivative optimum (~eps_f^(1/3)). A Hessian takes a second
  # difference, whose optimum scales as eps_f^(1/4) -- about 10x larger at a
  # machine-precision objective -- and whose noise term is 4*eps_f/h^2, so too
  # fine a step is exactly what makes a marginal Hessian indefinite. Measured 6x
  # to 385x worse before this was separated.
  p <- c(1.6, 3.0, -2.4)
  f <- function(x) sum(exp(x))
  h_hess <- admixr2:::.admHessSteps(f, p, 1L, cov_h_outer = 1e-4)
  h_grad <- admixr2:::.admShi21Steps(f, p, 1L, fallback = 1e-5)
  expect_true(is.finite(h_hess) && h_hess > 0)
  expect_gt(h_hess, h_grad)

  # ... and it really is better for a second difference.
  d2 <- function(h) { q1 <- p; q1[1] <- q1[1] + h; q2 <- p; q2[1] <- q2[1] - h
                      (f(q1) - 2 * f(p) + f(q2)) / h^2 }
  ref <- exp(p[1])
  expect_lt(abs(d2(h_hess) - ref) / ref, abs(d2(h_grad) - ref) / ref)
})

test_that(".admHessSteps falls back to the fixed scale when noise cannot be measured", {
  p <- c(0.5, -1.2)
  flat <- function(x) 1
  expect_warning(h <- admixr2:::.admHessSteps(flat, p, 1:2, cov_h_outer = 1e-3,
                                              .var.name = "unit"),
                 "noise level")
  expect_equal(h, pmax(abs(p), 0.1) * 1e-3)
})

test_that("cov_h_outer still scales the Hessian step", {
  # The documented escape hatch for an indefinite Hessian is "increase
  # cov_h_outer". Under a fallback-only design that argument would be inert
  # whenever the noise estimate succeeds -- which is almost always -- so it
  # multiplies the measured step instead. A user who raises it by 100 must get a
  # step 100x larger, or the advice is a lie.
  p <- c(1.6, 3.0, -2.4)
  f <- function(x) sum(exp(x))
  h1 <- admixr2:::.admHessSteps(f, p, 1L, cov_h_outer = .Machine$double.eps^(1/5))
  h2 <- admixr2:::.admHessSteps(f, p, 1L, cov_h_outer = .Machine$double.eps^(1/5) * 100)
  expect_equal(h2 / h1, 100, tolerance = 1e-8)
})

test_that("a locally QUADRATIC objective falls back instead of returning the cap", {
  # f''' is exactly zero for a quadratic, so the third difference never clears
  # the noise floor, there is no h* to compute, and the probe exits holding the
  # grown h -- which is the cap, a step the size of the parameter. Returning that
  # is worse than not measuring at all: a converged NLL is near-quadratic, so
  # this is the regime fits actually end in, and a step that size perturbs omega
  # out of positive-definiteness or sigma to Inf.
  p  <- c(0.5, 1.5, 2.5)
  fq <- function(x) sum((x - c(1, 2, 3))^2)
  r <- admixr2:::.admShi21Central(fq, p, 1L, 1e-12)
  expect_false(isTRUE(r$measured))

  fb <- c(1e-5, 2e-5, 3e-5)
  h  <- admixr2:::.admShi21Steps(fq, p, 1:3, fallback = fb)
  expect_identical(h, fb)

  # A cubic HAS a third derivative, so the same call measures rather than falls
  # back -- the guard must not swallow the ordinary case.
  fc <- function(x) sum((x - c(1, 2, 3))^3) + sum(x^2)
  expect_true(isTRUE(admixr2:::.admShi21Central(fc, p, 1L, 1e-12)$measured))
})
