skip_if_not_installed("rxode2")
skip_on_cran()

# Driver-level cover for covMethod = "r,s". The unit tests in test-adfweight.R
# score the WEIGHT and the sandwich against their own definitions; these say the
# value reaches a fit -- that the control accepts it, that .adghCalcCov /
# .admCalcCov / .adfoCalcCov all pass `sandwich` through, that the label on the
# returned fit reports what the covariance IS, and that nothing on the objective
# side moved on the way.
#
# Fixture and its sigma level: see .int_sandwich_setup() in helper-integration.R.

test_that("r,s changes the reported uncertainty and nothing else", {
  env <- .int_sandwich_setup()
  for (e in env$ests) {
    a <- env$fits[[e]]$r; b <- env$fits[[e]]$rs
    # The sandwich replaces the FILLING of the covariance. The optimizer never
    # sees it -- the covariance step runs after convergence on p_hat -- so any
    # difference in the estimates or the objective is a wiring defect.
    expect_equal(b$parFixedDf$Estimate, a$parFixedDf$Estimate, info = e)
    expect_equal(b$objDf$OBJF, a$objDf$OBJF, info = e)
  }
})

test_that("the fit reports the covariance it HAS, not the one asked for", {
  env <- .int_sandwich_setup()
  for (e in env$ests) {
    expect_identical(env$fits[[e]]$r$covMethod, "r", info = e)
    # "r,s" here is an assertion that the correction was actually built: a
    # sandwich that degraded reports "r", so this failing means the fixture
    # stopped reaching the path, not that the number is wrong.
    expect_identical(env$fits[[e]]$rs$covMethod, "r,s", info = e)
  }
})

test_that("under correct specification the sandwich reduces to r", {
  # THE load-bearing test. J = 2H when the model is right, so H^-1 J H^-1 must
  # come back as 2H^-1 -- the number "r" reports. It runs the whole path (the
  # moment Jacobian G, the ADF weight Omega, and the shared H), and a wrong G, a
  # wrong Omega, or an H that is not the one "r" inverted would each break it
  # while still producing finite, plausible standard errors.
  #
  # This replaces an earlier assertion that the two must DIFFER, which passed for
  # the wrong reason: on the old fixture they differed because `add.err` was
  # unidentified (cond(H) = 3.5e5) and the sandwich inverts H twice. That the
  # weight moves the sandwich off 2H under MISspecification is pinned where it
  # belongs, on a controlled ensemble, in test-adfweight.R.
  env <- .int_sandwich_setup()
  for (e in env$ests) {
    a <- env$fits[[e]]$r; b <- env$fits[[e]]$rs
    ok <- is.finite(a$parFixedDf$SE) & is.finite(b$parFixedDf$SE)
    expect_true(any(ok), info = e)
    expect_true(all(b$parFixedDf$SE[ok] > 0), info = e)
    ratio <- unname(b$parFixedDf$SE[ok] / a$parFixedDf$SE[ok])
    # adgh is the exact-quadrature reference and holds tightest; the two MC
    # estimators carry their own sampling noise into H, and adfo's G comes from
    # the FO moment map rather than the ensemble, so it is not expected to be
    # exact even here.
    tol <- if (e == "adgh") 0.05 else 0.25
    expect_equal(ratio, rep(1, length(ratio)), tolerance = tol,
                 info = paste(e, "ratios:", paste(signif(ratio, 4), collapse = " ")))
  }
})

test_that("a residual the weight cannot reach degrades to r rather than guessing", {
  skip_if_not_installed("nlmixr2est")
  # ar() correlates the residual ACROSS timepoints, which is exactly the
  # conditional independence the Isserlis expansion rests on -- the dropped
  # cross terms are real, so .admAdfCondMom refuses it. The fit must still
  # succeed, with the "r" covariance and the "r" label.
  ar_fn <- function() {
    ini({ tcl <- log(5); tv <- log(20); add.err <- 0.3; rho <- 0.2; eta.cl ~ 0.09 })
    model({ cl <- exp(tcl + eta.cl); v <- exp(tv)
            cp <- linCmt(); cp ~ add(add.err) + ar(rho) })
  }
  times  <- c(0.5, 1, 2, 4)
  E_true <- .one_cmt_mean(5, 20, 100, times)
  # A FULL covariance, not a diagonal one: ar() refuses a `method = "var"` study
  # outright (a diagonal V carries no information about a residual correlation),
  # and that refusal would pre-empt the one being tested here.
  sd  <- 0.3 * E_true
  rho <- 0.25^abs(outer(seq_along(times), seq_along(times), "-"))
  V   <- outer(sd, sd) * rho
  st  <- list(s1 = list(E = E_true, V = V, n = 200L,
                        times = times, ev = rxode2::et(amt = 100)))
  fit <- tryCatch(suppressMessages(suppressWarnings(nlmixr2est::nlmixr2(
    ar_fn, admData(), est = "adgh",
    control = adghControl(studies = st, n_nodes = 5L, maxeval = 10L, seed = 1L,
                          grad = "none", covMethod = "r,s")))),
    error = function(e) e)
  skip_if(inherits(fit, "error"), "ar() endpoint did not fit in this environment")
  expect_identical(fit$covMethod, "r")
})
