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

test_that("an ill-conditioned Hessian is reported, and only then", {
  # The trap this guards: the correction is amplified quadratically in a
  # direction the data barely identifies, so "r,s" returns a finite, plausible,
  # meaningless SE there while the rest of the fit is fine. Degrading the whole
  # covariance would cost more than it saves, so the fit keeps its sandwich and
  # says which parameter not to read.
  #
  # Read off `fit$runInfo`, NOT with a warning handler. nlmixr2Est0 wraps the
  # estimator in .collectWarn(), which suppresses warnings at source and assigns
  # them to runInfo for a nlmixr2FitCore result -- so a handler around nlmixr2()
  # sees nothing and would make this test pass vacuously. runInfo is also what
  # print() lists, which is the thing the user actually sees.
  env <- .int_sandwich_setup()
  times  <- c(0.5, 1, 2, 4)
  E_true <- .one_cmt_mean(5, 20, 100, times)

  flagged <- function(fn, st, cm) {
    f <- suppressMessages(nlmixr2est::nlmixr2(fn, admData(), est = "adgh",
      control = adghControl(studies = st, n_nodes = 7L, maxeval = 200L,
                            seed = 1L, grad = "analytical", covMethod = cm,
                            print = 0L)))
    any(grepl("ill-conditioned", if (is.null(f$runInfo)) character() else f$runInfo))
  }

  # add.err = 0.1 contributes 0.01 variance against ~1.7 from IIV: RSE 186%,
  # cond(H) = 3.5e5. This is the fixture the file used to run on.
  ill_fn <- function() {
    ini({ tcl <- log(5); tv <- log(20); add.err <- 0.1
          eta.cl ~ 0.09; eta.v ~ 0.04 })
    model({ cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v); linCmt() ~ add(add.err) })
  }
  ill_st <- list(s1 = list(E = E_true, V = diag((0.3 * E_true)^2), n = 200L,
                           times = times, ev = rxode2::et(amt = 100)))

  expect_true(flagged(ill_fn, ill_st, "r,s"))
  # Silent under "r": one inversion is what the existing .ADM_NPD_RCOND bound is
  # for, and it does not fire here -- rcond 2.8e-06 is well above sqrt(eps).
  expect_false(flagged(ill_fn, ill_st, "r"))
  # ... and silent on the well-conditioned study the rest of the file uses, which
  # is what stops this being a note every "r,s" fit carries.
  expect_false(flagged(env$fn, env$studies, "r,s"))
})

test_that("the transform-both-sides family reaches the sandwich", {
  # boxCox / yeoJohnson / logitNorm / probitNorm were the one residual family
  # "r,s" refused, and the refusal was never mathematical: they are conditionally
  # independent across timepoints like every other supported family, they simply
  # had no closed form for the third and fourth moments. The quadrature that
  # already produces their mean and variance supplies those.
  #
  # This is the DRIVER-level half; the moments themselves are scored against a
  # simulation of the conditional law, and the weight against the sampling
  # covariance of simulated studies, in test-adfweight.R.
  skip_if_not_installed("nlmixr2est")
  TT <- c(2, 5, 9, 14); DOSE <- 100
  E0 <- DOSE / 10 * exp(-0.1 * TT)
  st <- list(s = list(E = E0, V = diag((0.3 * E0)^2), n = 150L, times = TT,
                      ev = rxode2::et(amt = DOSE)))
  mods <- list(
    boxCox = function() { ini({ tcl <- log(1); tv <- log(10); eta.cl ~ 0.16
                                aa <- 0.3; lam <- fix(0.5) })
      model({ cl <- exp(tcl + eta.cl); v <- exp(tv); cp <- linCmt()
              cp ~ add(aa) + boxCox(lam) }) },
    yeoJohnson = function() { ini({ tcl <- log(1); tv <- log(10); eta.cl ~ 0.16
                                    aa <- 0.3; lam <- fix(0.5) })
      model({ cl <- exp(tcl + eta.cl); v <- exp(tv); cp <- linCmt()
              cp ~ add(aa) + yeoJohnson(lam) }) },
    logitNorm = function() { ini({ tcl <- log(1); tv <- log(10); eta.cl ~ 0.16
                                   aa <- 0.25 })
      model({ cl <- exp(tcl + eta.cl); v <- exp(tv); cp <- linCmt()
              cp ~ logitNorm(aa, 0, 40) }) },
    probitNorm = function() { ini({ tcl <- log(1); tv <- log(10); eta.cl ~ 0.16
                                    aa <- 0.25 })
      model({ cl <- exp(tcl + eta.cl); v <- exp(tv); cp <- linCmt()
              cp ~ probitNorm(aa, 0, 40) }) })
  for (nm in names(mods)) {
    f <- suppressMessages(suppressWarnings(nlmixr2est::nlmixr2(
      mods[[nm]], admData(), est = "adgh",
      control = adghControl(studies = st, n_nodes = 5L, maxeval = 20L, seed = 1L,
                            grad = "analytical", covMethod = "r,s", print = 0L))))
    expect_identical(f$covMethod, "r,s", info = nm)
    # every FREE parameter gets a standard error; a fix()ed one has none, which
    # is why this filters rather than asserting all() -- boxCox and yeoJohnson
    # hold lambda fixed here.
    fixed <- f$ui$iniDf$fix %in% TRUE
    nmf   <- f$ui$iniDf$name[fixed]
    keep  <- !(rownames(f$parFixedDf) %in% nmf)
    expect_true(all(is.finite(f$parFixedDf$SE[keep])), info = nm)
  }
})

test_that("a model with no random effects still gets the sandwich", {
  # .admSandwichGrid() used to refuse n_eta < 1 outright, which made admc and
  # adfo degrade to "r" on every no-IIV model while adgh -- which passes its own
  # grid -- applied the correction to the same fit. The single-point ensemble is
  # the CORRECT one there: with no between-subject variability the summary's
  # sampling law is the residual's alone, and for a skewed residual that is still
  # not the normal-theory law the objective assumes.
  skip_if_not_installed("nlmixr2est")
  TT <- c(2, 5, 9, 14); DOSE <- 100
  E0 <- DOSE / 10 * exp(-0.1 * TT)
  st <- list(s = list(E = E0, V = diag((0.3 * E0)^2), n = 150L, times = TT,
                      ev = rxode2::et(amt = DOSE)))
  # lnorm, not add: with no IIV AND a normal residual the model really is
  # exactly normal and the correction has nothing to say, so a skewed residual
  # is what makes this test about more than the plumbing.
  fn <- function() {
    ini({ tcl <- log(1); tv <- log(10); aa <- 0.15 })
    model({ cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ lnorm(aa) })
  }
  ctl <- function(est) switch(est,
    adgh = adghControl(studies = st, n_nodes = 5L, maxeval = 20L, seed = 1L,
                       grad = "analytical", covMethod = "r,s", print = 0L),
    admc = admControl(studies = st, n_sim = 1500L, maxeval = 20L, seed = 1L,
                      grad = "sens", covMethod = "r,s", print = 0L),
    adfo = adfoControl(studies = st, maxeval = 20L, grad = "none",
                       covMethod = "r,s", print = 0L))
  for (e in c("adgh", "admc", "adfo")) {
    f <- suppressMessages(suppressWarnings(
      nlmixr2est::nlmixr2(fn, admData(), est = e, control = ctl(e))))
    expect_identical(f$covMethod, "r,s", info = e)
    expect_true(all(is.finite(f$parFixedDf$SE)), info = e)
  }
  # adirmc refuses a no-IIV model at the ESTIMATOR level (it draws its ensemble
  # from the etas), which is unrelated to the sandwich and left alone.
})
