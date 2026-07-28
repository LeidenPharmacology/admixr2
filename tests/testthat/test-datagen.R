# ---- datagenControl() unit tests (no rxode2) -----------------------------------

test_that("datagenControl() returns correct class and defaults", {
  ctl <- datagenControl()
  expect_s3_class(ctl, "datagenControl")
  expect_equal(ctl$method,             "mc")
  expect_equal(ctl$n_sim,              5000L)
  expect_equal(ctl$sampling,           "sobol")
  expect_equal(ctl$seed,               12345L)
  expect_equal(ctl$cores,              1L)
  expect_false(ctl$return_samples)
  # Must match the estimator controls' default, or a study generated here and the
  # fit that consumes it would integrate a transform-both-sides residual at
  # different accuracies -- a discrepancy that looks like model misspecification.
  expect_identical(ctl$resid_nodes,     81L)
  expect_identical(admControl()$resid_nodes, ctl$resid_nodes)
})

test_that("datagenControl(): resid_nodes is coerced and validated", {
  expect_identical(datagenControl(resid_nodes = 31)$resid_nodes, 31L)
  expect_error(datagenControl(resid_nodes = 0L))
  expect_error(datagenControl(resid_nodes = c(31L, 81L)))
})

test_that("datagenControl(): method = 'fo' accepted", {
  ctl <- datagenControl(method = "fo")
  expect_equal(ctl$method, "fo")
})

test_that("datagenControl(): invalid method errors", {
  expect_error(datagenControl(method = "laplace"))
})

test_that("datagenControl(): n_sim stored as integer", {
  ctl <- datagenControl(n_sim = 1000)
  expect_type(ctl$n_sim, "integer")
  expect_equal(ctl$n_sim, 1000L)
})

test_that("datagenControl(): seed stored as integer", {
  ctl <- datagenControl(seed = 42)
  expect_type(ctl$seed, "integer")
  expect_equal(ctl$seed, 42L)
})

test_that("datagenControl(): cores stored as integer", {
  ctl <- datagenControl(cores = 2)
  expect_type(ctl$cores, "integer")
  expect_equal(ctl$cores, 2L)
})

test_that("datagenControl(): sampling match.arg works", {
  for (s in c("sobol", "halton", "torus", "lhs", "rnorm")) {
    ctl <- datagenControl(sampling = s)
    expect_equal(ctl$sampling, s)
  }
})

test_that("datagenControl(): invalid sampling errors", {
  expect_error(datagenControl(sampling = "invalid"))
})

test_that("datagenControl(): n_sim = 0 errors via checkmate", {
  expect_error(datagenControl(n_sim = 0), regexp = "n_sim")
})

test_that("datagenControl(): negative n_sim errors", {
  expect_error(datagenControl(n_sim = -10), regexp = "n_sim")
})

test_that("datagenControl(): return_samples = TRUE accepted", {
  ctl <- datagenControl(return_samples = TRUE)
  expect_true(ctl$return_samples)
})

# ---- datagen() validation errors (no rxode2 needed) ---------------------------

test_that("datagen(): errors when no model supplied and no per-study model", {
  studies <- list(s1 = list(times = 1:3, ev = list(), n = 100L))
  expect_error(
    datagen(studies, model = NULL),
    regexp = "no.*model.*no top-level default"
  )
})

test_that("datagen(): errors when per-study model is not a function", {
  studies <- list(s1 = list(model = "not_a_function",
                            times = 1:3, ev = list(), n = 100L))
  expect_error(
    datagen(studies),
    regexp = "must be a function"
  )
})

test_that("datagen(): errors when study missing times", {
  fake_model <- function() {}
  studies <- list(s1 = list(model = fake_model, ev = list(), n = 100L))
  expect_error(
    datagen(studies),
    regexp = "missing.*times"
  )
})

test_that("datagen(): errors when study missing ev", {
  fake_model <- function() {}
  studies <- list(s1 = list(model = fake_model, times = 1:3, n = 100L))
  expect_error(
    datagen(studies),
    regexp = "missing.*ev"
  )
})

test_that("datagen(): errors when studies is empty", {
  expect_error(datagen(list()))
})

test_that("datagen(): errors when control is not datagenControl", {
  fake_model <- function() {}
  studies <- list(s1 = list(model = fake_model, times = 1:3, ev = list()))
  expect_error(
    datagen(studies, control = list()),
    regexp = "datagenControl"
  )
})

# ---- datagen() integration test (requires rxode2) --------------------------------

test_that("datagen(): produces E and V with correct dimensions for single study", {
  skip_if_not_installed("rxode2")
  skip_on_cran()

  pk_model <- function() {
    ini({
      tcl <- log(5)
      tv  <- log(30)
      prop.sd <- c(0, 0.2)
      eta.cl ~ 0.09
      eta.v  ~ 0.04
    })
    model({
      cl <- exp(tcl + eta.cl)
      v  <- exp(tv  + eta.v)
      d/dt(central) <- -(cl/v) * central
      cp <- central / v
      cp ~ prop(prop.sd)
    })
  }

  times   <- c(1, 2, 4, 8)
  studies <- list(
    s1 = list(times = times, ev = rxode2::et(amt = 100), n = 200L)
  )

  out <- datagen(studies, model = pk_model,
                 control = datagenControl(n_sim = 500L))

  expect_named(out, "s1")
  s <- out$s1
  expect_length(s$E, length(times))
  expect_equal(dim(s$V), c(length(times), length(times)))
  expect_true(all(is.finite(s$E)))
  expect_true(all(is.finite(s$V)))
  expect_equal(s$n, 200L)
  expect_equal(s$times, times)
})

test_that("datagen(): return_samples=TRUE includes samples matrix", {
  skip_if_not_installed("rxode2")
  skip_on_cran()

  pk_model <- function() {
    ini({
      tcl <- log(5)
      tv  <- log(30)
      prop.sd <- c(0, 0.2)
      eta.cl ~ 0.09
      eta.v  ~ 0.04
    })
    model({
      cl <- exp(tcl + eta.cl)
      v  <- exp(tv  + eta.v)
      d/dt(central) <- -(cl/v) * central
      cp <- central / v
      cp ~ prop(prop.sd)
    })
  }

  times   <- c(1, 4)
  n_sim   <- 300L
  studies <- list(s1 = list(times = times, ev = rxode2::et(amt = 100)))

  out <- datagen(studies, model = pk_model,
                 control = datagenControl(n_sim = n_sim, return_samples = TRUE))

  expect_true(!is.null(out$s1$samples))
  expect_equal(dim(out$s1$samples), c(n_sim, length(times)))
})

test_that("datagen(): a model without a residual-error term yields IIV-only V", {
  skip_if_not_installed("rxode2")
  skip_on_cran()

  # With residual error.
  pk_err <- function() {
    ini({
      tcl <- log(5)
      tv  <- log(30)
      prop.sd <- c(0, 0.2)
      eta.cl ~ 0.09
      eta.v  ~ 0.04
    })
    model({
      cl <- exp(tcl + eta.cl)
      v  <- exp(tv  + eta.v)
      d/dt(central) <- -(cl/v) * central
      cp <- central / v
      cp ~ prop(prop.sd)
    })
  }
  # Same structural model, no error term -> V should carry IIV only.
  pk_noerr <- function() {
    ini({
      tcl <- log(5)
      tv  <- log(30)
      eta.cl ~ 0.09
      eta.v  ~ 0.04
    })
    model({
      cl <- exp(tcl + eta.cl)
      v  <- exp(tv  + eta.v)
      d/dt(central) <- -(cl/v) * central
      cp <- central / v
    })
  }

  times   <- c(1, 4)
  studies <- list(s1 = list(times = times, ev = rxode2::et(amt = 100)))

  v_err   <- datagen(studies, model = pk_err,
                     control = datagenControl(n_sim = 500L))$s1$V
  v_noerr <- datagen(studies, model = pk_noerr,
                     control = datagenControl(n_sim = 500L))$s1$V

  # Proportional error inflates only the diagonal; off-diagonal IIV structure
  # matches up to MC noise.
  expect_true(all(diag(v_err) > diag(v_noerr)))
})

test_that("datagen(): multi-study with per-study model", {
  skip_if_not_installed("rxode2")
  skip_on_cran()

  pk_model <- function() {
    ini({
      tcl <- log(5)
      tv  <- log(30)
      prop.sd <- c(0, 0.2)
      eta.cl ~ 0.09
      eta.v  ~ 0.04
    })
    model({
      cl <- exp(tcl + eta.cl)
      v  <- exp(tv  + eta.v)
      d/dt(central) <- -(cl/v) * central
      cp <- central / v
      cp ~ prop(prop.sd)
    })
  }

  studies <- list(
    s1 = list(times = c(1, 2), ev = rxode2::et(amt = 100), n = 100L),
    s2 = list(times = c(1, 4, 8), ev = rxode2::et(amt = 200), n = 50L)
  )

  out <- datagen(studies, model = pk_model,
                 control = datagenControl(n_sim = 300L))

  expect_named(out, c("s1", "s2"))
  expect_length(out$s1$E, 2L)
  expect_length(out$s2$E, 3L)
})

test_that("datagen(method='fo'): produces finite E and V with correct dims", {
  skip_if_not_installed("rxode2")
  skip_on_cran()

  pk_model <- function() {
    ini({
      tcl <- log(5)
      tv  <- log(30)
      prop.sd <- c(0, 0.2)
      eta.cl ~ 0.09
      eta.v  ~ 0.04
    })
    model({
      cl <- exp(tcl + eta.cl)
      v  <- exp(tv  + eta.v)
      d/dt(central) <- -(cl/v) * central
      cp <- central / v
      cp ~ prop(prop.sd)
    })
  }

  times   <- c(1, 2, 4, 8)
  studies <- list(s1 = list(times = times, ev = rxode2::et(amt = 100), n = 200L))

  out <- datagen(studies, model = pk_model,
                 control = datagenControl(method = "fo"))

  s <- out$s1
  expect_length(s$E, length(times))
  expect_equal(dim(s$V), c(length(times), length(times)))
  expect_true(all(is.finite(s$E)))
  expect_true(all(is.finite(s$V)))
  # V is a valid covariance: symmetric, positive diagonal
  expect_equal(s$V, t(s$V))
  expect_true(all(diag(s$V) > 0))
  expect_equal(s$n, 200L)
})

test_that("datagen(method='fo'): E matches the closed-form f(theta, 0)", {
  skip_if_not_installed("rxode2")
  skip_on_cran()

  # 1-cmt bolus, eta = 0: cl = 5, v = 30, dose 100 -> cp(t) = (100/30) exp(-t/6).
  # Proportional error leaves the FO mean uncorrected, so E should equal it.
  pk_model <- function() {
    ini({
      tcl <- log(5)
      tv  <- log(30)
      prop.sd <- c(0, 0.2)
      eta.cl ~ 0.09
      eta.v  ~ 0.04
    })
    model({
      cl <- exp(tcl + eta.cl)
      v  <- exp(tv  + eta.v)
      d/dt(central) <- -(cl/v) * central
      cp <- central / v
      cp ~ prop(prop.sd)
    })
  }

  times   <- c(1, 2, 4, 8)
  studies <- list(s1 = list(times = times, ev = rxode2::et(amt = 100)))

  out <- datagen(studies, model = pk_model,
                 control = datagenControl(method = "fo"))

  expected_E <- (100 / 30) * exp(-times / 6)
  expect_equal(unname(out$s1$E), expected_E, tolerance = 1e-4)
})

test_that("datagen(method='fo'): deterministic (seed/n_sim irrelevant)", {
  skip_if_not_installed("rxode2")
  skip_on_cran()

  pk_model <- function() {
    ini({
      tcl <- log(5)
      tv  <- log(30)
      prop.sd <- c(0, 0.2)
      eta.cl ~ 0.09
      eta.v  ~ 0.04
    })
    model({
      cl <- exp(tcl + eta.cl)
      v  <- exp(tv  + eta.v)
      d/dt(central) <- -(cl/v) * central
      cp <- central / v
      cp ~ prop(prop.sd)
    })
  }

  studies <- list(s1 = list(times = c(1, 4), ev = rxode2::et(amt = 100)))

  a <- datagen(studies, model = pk_model,
               control = datagenControl(method = "fo", seed = 1L))
  b <- datagen(studies, model = pk_model,
               control = datagenControl(method = "fo", seed = 99L))

  expect_identical(a$s1$E, b$s1$E)
  expect_identical(a$s1$V, b$s1$V)
})

test_that("datagen(): unnamed studies get auto-names", {
  skip_if_not_installed("rxode2")
  skip_on_cran()

  pk_model <- function() {
    ini({
      tcl <- log(5)
      tv  <- log(30)
      prop.sd <- c(0, 0.2)
      eta.cl ~ 0.09
      eta.v  ~ 0.04
    })
    model({
      cl <- exp(tcl + eta.cl)
      v  <- exp(tv  + eta.v)
      d/dt(central) <- -(cl/v) * central
      cp <- central / v
      cp ~ prop(prop.sd)
    })
  }

  studies <- list(
    list(times = c(1, 2), ev = rxode2::et(amt = 100))
  )

  out <- datagen(studies, model = pk_model,
                 control = datagenControl(n_sim = 300L))

  expect_equal(names(out), "study1")
})

test_that("datagen() refuses the endpoints it cannot generate, rather than emitting NAs", {
  skip_on_cran(); skip_if_not_installed("rxode2")

  # ordinal: its categories are ONE joint observation whose covariance carries the
  # -p_j*p_k term between categories at the same time. datagen() derives each
  # observed output separately, so that block cannot be formed here at all.
  ord <- function() {
    ini({ tcl <- log(5); tv <- log(20); eta.cl ~ 0.09 })
    model({ cl <- exp(tcl + eta.cl); v <- exp(tv)
            d/dt(central) <- -(cl / v) * central
            cp <- central / v
            p1 <- 1 / (1 + exp(-(cp - 1)))
            p2 <- 1 - p1
            y ~ c(p1, p2) })
  }
  expect_error(
    datagen(list(s = list(times = c(1, 2), ev = rxode2::et(amt = 100), n = 50L)),
            model = ord, control = datagenControl(n_sim = 100L, seed = 1L)),
    "ordinal")

  # beta under method = "fo": .adfoVpred never sees the solved b1 + b2, so the
  # diagonal of the emitted V would be entirely NA.
  bet <- function() {
    ini({ tem <- log(0.6); tphi <- log(30); eta.em ~ 0.05 })
    model({ em <- exp(tem + eta.em); phi <- exp(tphi)
            mu <- em / (1 + em); b1 <- mu * phi; b2 <- (1 - mu) * phi
            y ~ beta(b1, b2) })
  }
  expect_error(
    datagen(list(s = list(times = c(1, 2), ev = rxode2::et(amt = 0), n = 50L)),
            model = bet, control = datagenControl(method = "fo", seed = 1L)),
    "beta")
})
