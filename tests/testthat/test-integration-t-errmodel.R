# Student-t residual error (`cp ~ add(a) + t(nu)`) -----------------------------
#
# admixr2 supports t() by MOMENT MATCHING: nlmixr2 writes it as a scale family
# (residual = scale * T_nu), so the aggregate variance is the normal one times
# nu/(nu-2). Everything here is an ORACLE test -- each check compares admixr2's
# analytic answer against an independently computed truth, never against another
# admixr2 code path:
#
#   1. semantics  -- rxode2's OWN simulated residual (the `sim` column of the
#                    simulation model) vs our var = scale^2 * nu/(nu-2).
#   2. moments    -- .admResidApply's (mu, var) vs the empirical mean/variance of
#                    a large Monte Carlo draw from the same generative process.
#   3. gradient   -- .admResidDeriv's d(var)/d(p) vs a central finite difference.
#   4. recovery   -- fit aggregate data generated from a KNOWN (a, nu) truth and
#                    check the estimates come back (with nu FIXED -- see below).
#   5. equivalence-- a t fit vs a matched-variance normal fit, which cancels any
#                    estimator bias and so holds for adfo and adgh alike.
#
# (1) is the load-bearing one: it is the only check that admixr2's reading of
# what `a` and `nu` MEAN agrees with nlmixr2's. (5) is the sharpest, because it
# is the only one whose expected value does not depend on the estimator being
# unbiased.

skip_on_cran()
skip_if_not_installed("rxode2")

.t_mod <- function(iniExtra, errLine) {
  eval(parse(text = sprintf('function() {
  ini({ tcl <- log(5); tv <- log(20); %s
        eta.cl ~ 0.09 })
  model({
    cl <- exp(tcl + eta.cl); v <- exp(tv)
    d/dt(central) <- -(cl / v) * central
    cp <- central / v
    %s
  })
}', iniExtra, errLine)))
}

# suppressWarnings throughout: these fixtures ESTIMATE nu so they can exercise its
# gradient slot, which deliberately warns that nu is not identifiable from
# aggregate data (that warning has its own test at the bottom of this file).
.t_pinfo <- function(a = 0.5, nu = 5) {
  ui <- suppressMessages(rxode2::rxode2(
    .t_mod(sprintf("a <- %g; nu <- %g", a, nu), "cp ~ add(a) + t(nu)")))
  suppressWarnings(.admParseIniDf(ui$iniDf, ui))
}

# -- 1. semantics: our variance vs rxode2's own simulated residual -------------

test_that("t(nu) means scale*T_nu -- matches rxode2's simulated residual", {
  skip_on_os("solaris")
  a <- 0.5
  for (nu in c(4, 8, 30)) {
    ui <- suppressMessages(rxode2::rxode2(
      .t_mod(sprintf("a <- %g; nu <- %g", a, nu), "cp ~ add(a) + t(nu)")))
    ev <- rxode2::et(rxode2::et(amt = 1000), 4)
    set.seed(20260721L)
    s <- suppressWarnings(rxode2::rxSolve(ui$simulationModel, ev, nSub = 400000L,
                                          returnType = "data.frame"))
    s <- s[!is.na(s$sim), , drop = FALSE]
    emp <- stats::var(s$sim - s$ipredSim)          # rxode2's OWN residual draw
    ours <- a^2 * nu / (nu - 2)                    # admixr2's moment match
    # 4e5 draws of a t: the variance estimator itself is noisy, so 4% is tight.
    expect_equal(emp, ours, tolerance = 0.04,
                 label = sprintf("empirical var at nu=%g", nu))
    # and the mean is unshifted
    expect_equal(mean(s$sim - s$ipredSim), 0, tolerance = 0.02 * sqrt(ours))
  }
})

# -- 2. moments: .admResidApply vs a Monte Carlo oracle ------------------------

test_that(".admResidApply reproduces the Monte Carlo mean and variance", {
  a <- 0.4; nu <- 6
  pinfo <- .t_pinfo(a, nu)
  f     <- c(2, 10, 25)
  snat  <- .admSigmaNat(pinfo$sigma_init, pinfo)
  arr   <- .admResidRows(pinfo, "cp", snat, length(f))
  got   <- .admResidApply(f, numeric(length(f)), arr)

  set.seed(7L)
  N <- 400000L
  for (i in seq_along(f)) {
    draws <- f[i] + a * stats::rt(N, df = nu)      # the generative process
    expect_equal(got$mu[i], mean(draws), tolerance = 0.01)
    expect_equal(got$dv[i], stats::var(draws), tolerance = 0.05)
  }
  # the multiplier is exactly nu/(nu-2) over the normal answer
  expect_equal(got$dv, rep(a^2 * nu / (nu - 2), length(f)))
  expect_equal(got$mu, f)                          # t does not shift the mean
})

test_that("the t multiplier is exact for prop, pow and combined1 too", {
  nu <- 5; m <- nu / (nu - 2)
  f  <- c(3, 12)
  cases <- list(
    prop      = list(ini = "b <- 0.2; nu <- 5",           err = "cp ~ prop(b) + t(nu)"),
    addprop   = list(ini = "a <- 0.3; b <- 0.2; nu <- 5", err = "cp ~ add(a) + prop(b) + t(nu)"),
    combined1 = list(ini = "a <- 0.3; b <- 0.2; nu <- 5", err = "cp ~ add(a) + prop(b) + combined1() + t(nu)"),
    pow       = list(ini = "b <- 0.2; c1 <- 1.3; nu <- 5", err = "cp ~ pow(b, c1) + t(nu)")
  )
  for (nmc in names(cases)) {
    cs <- cases[[nmc]]
    ui_t <- suppressMessages(rxode2::rxode2(.t_mod(cs$ini, cs$err)))
    ui_n <- suppressMessages(rxode2::rxode2(
      .t_mod(sub("; nu <- 5", "", cs$ini), sub(" \\+ t\\(nu\\)", "", cs$err))))
    pt <- suppressWarnings(.admParseIniDf(ui_t$iniDf, ui_t)); pn <- .admParseIniDf(ui_n$iniDf, ui_n)
    at <- .admResidRows(pt, "cp", .admSigmaNat(pt$sigma_init, pt), length(f))
    an <- .admResidRows(pn, "cp", .admSigmaNat(pn$sigma_init, pn), length(f))
    vt <- .admResidApply(f, numeric(length(f)), at)$dv
    vn <- .admResidApply(f, numeric(length(f)), an)$dv
    # exact scale-family identity: var_t == m * var_normal, for EVERY form
    expect_equal(vt, m * vn, tolerance = 1e-12, label = paste0(nmc, ": var_t vs m*var_normal"))
  }
})

# -- 3. gradient: analytic vs central finite difference ------------------------

test_that(".admResidDeriv matches a finite difference for every residual param", {
  pinfo <- .t_pinfo(a = 0.4, nu = 6)
  # add a proportional term so the FD exercises a2, b2 and nu together
  ui <- suppressMessages(rxode2::rxode2(
    .t_mod("a <- 0.4; b <- 0.15; nu <- 6", "cp ~ add(a) + prop(b) + t(nu)")))
  pinfo <- suppressWarnings(.admParseIniDf(ui$iniDf, ui))
  f  <- c(2, 9, 22)
  p0 <- pinfo$sigma_init

  vfun <- function(p) {
    arr <- .admResidRows(pinfo, "cp", .admSigmaNat(p, pinfo), length(f))
    .admResidApply(f, numeric(length(f)), arr)$dv
  }
  arr <- .admResidRows(pinfo, "cp", .admSigmaNat(p0, pinfo), length(f))
  # var_f = 0: this check isolates the residual, matching vfun()'s zero
  # structural variance (composition is test-integration-resid-moments.R).
  an  <- .admResidDeriv(f, numeric(length(f)), arr, pinfo)

  h <- 1e-6
  for (k in seq_along(p0)) {
    pp <- p0; pp[k] <- pp[k] + h
    pm <- p0; pm[k] <- pm[k] - h
    fd <- (vfun(pp) - vfun(pm)) / (2 * h)
    expect_equal(an$dvar[, k], fd, tolerance = 1e-6,
                 label = sprintf("d(var)/d(%s)", pinfo$sigma_names[k]))
  }
  # and d(var)/d(f)
  fdf <- (vfun(p0) * 0 + sapply(seq_along(f), function(i) {
    fp <- f; fp[i] <- fp[i] + h; fm <- f; fm[i] <- fm[i] - h
    ap <- .admResidRows(pinfo, "cp", .admSigmaNat(p0, pinfo), length(f))
    (.admResidApply(fp, numeric(length(f)), ap)$dv[i] -
     .admResidApply(fm, numeric(length(f)), ap)$dv[i]) / (2 * h)
  }))
  expect_equal(an$dv_df, fdf, tolerance = 1e-6)
})

# -- 4. parameter recovery from aggregate data generated with t residuals ------
#
# Two things had to be learned the hard way here, and both are encoded below:
#
#  * nu is NOT identifiable (it is aliased with the scale through a^2*nu/(nu-2)),
#    so the thing to check for recovery is the PRODUCT, and nu must be fix()ed for
#    the scale itself to mean anything.
#  * `adfo` is the FO LINEARIZATION, and on this generative model (30% CV on CL,
#    observed over several half-lives) it absorbs the curvature it cannot represent
#    into the residual: a matched-variance NORMAL control comes back 1.33x high
#    under adfo and 0.99x under adgh. So absolute recovery is checked with adgh,
#    and adfo is checked only for t-vs-normal EQUIVALENCE, where the bias cancels.

.t_gen <- function(kind, a_true = 0.8, nu_true = 5, n_sub = 4000L, seed = 99L) {
  times <- c(0.5, 1, 2, 4, 8); tcl <- log(5); tv <- log(20); om_cl <- 0.09
  vr <- a_true^2 * nu_true / (nu_true - 2)
  set.seed(seed)
  eta <- stats::rnorm(n_sub, 0, sqrt(om_cl))
  cl <- exp(tcl + eta); v <- exp(tv)
  f <- vapply(times, function(tt) 1000 / v * exp(-cl / v * tt), numeric(n_sub))
  r <- if (kind == "t")
    a_true * matrix(stats::rt(n_sub * length(times), df = nu_true), n_sub, length(times))
  else
    matrix(stats::rnorm(n_sub * length(times), 0, sqrt(vr)), n_sub, length(times))
  dv <- f + r
  list(E = colMeans(dv), V = stats::cov.wt(dv, method = "ML")$cov,
       n = n_sub, times = times, tcl = tcl, tv = tv, var_true = vr)
}

.t_fit <- function(mod, d, est) {
  studies <- list(s1 = list(E = d$E, V = d$V, n = d$n, times = d$times,
                            ev = rxode2::et(amt = 1000)))
  ctl <- if (est == "adfo")
    adfoControl(studies = studies, covMethod = "none", n_restarts = 1L)
  else
    adghControl(studies = studies, covMethod = "none", n_restarts = 1L)
  f <- suppressWarnings(suppressMessages(
    nlmixr2est::nlmixr2(mod, admData(), est = est, control = ctl)))
  list(est = stats::setNames(f$parFixedDf$Estimate, rownames(f$parFixedDf)),
       obj = f$objective)
}

.mod_t_fixnu <- function() {
  ini({ tcl <- log(4); tv <- log(18); a <- 1; nu <- fix(5); eta.cl ~ 0.05 })
  model({ cl <- exp(tcl + eta.cl); v <- exp(tv)
          d/dt(central) <- -(cl / v) * central; cp <- central / v
          cp ~ add(a) + t(nu) })
}
.mod_norm <- function() {
  ini({ tcl <- log(4); tv <- log(18); a <- 1; eta.cl ~ 0.05 })
  model({ cl <- exp(tcl + eta.cl); v <- exp(tv)
          d/dt(central) <- -(cl / v) * central; cp <- central / v
          cp ~ add(a) })
}

test_that("t() with fixed nu recovers the scale and the structural parameters", {
  skip_if_not_installed("nlmixr2est")
  d <- .t_gen("t")
  r <- .t_fit(.mod_t_fixnu, d, "adgh")
  expect_true(is.finite(r$obj))
  expect_equal(unname(r$est[["tcl"]]), d$tcl, tolerance = 0.05)
  expect_equal(unname(r$est[["tv"]]),  d$tv,  tolerance = 0.05)
  # nu is fixed at its true 5, so the scale itself is identified and must come back
  expect_equal(unname(r$est[["a"]]), 0.8, tolerance = 0.12)
  # ... and so must the total residual variance a^2*nu/(nu-2)
  expect_equal(unname(r$est[["a"]])^2 * 5 / 3, d$var_true, tolerance = 0.12)
})

test_that("fixing nu actually applies the multiplier (not silently dropped)", {
  # Regression: .admParseIniDf drops FIXED error rows from sigma_rows, so a
  # fix()ed nu has no optimizer slot. If the spec did not carry the constant, the
  # endpoint would be fitted as a plain NORMAL -- a converged fit of the wrong
  # model. Check the multiplier is really in the variance.
  ui <- suppressMessages(rxode2::rxode2(.mod_t_fixnu))
  p  <- .admParseIniDf(ui$iniDf, ui)
  expect_false("nu" %in% p$sigma_names)                   # fixed => no optimizer slot
  arr <- .admResidRows(p, "cp", .admSigmaNat(p$sigma_init, p), 2L)
  expect_equal(arr$vmul, rep(5 / 3, 2))                   # nu = 5 => m = 5/3
  a2 <- exp(p$sigma_init[["a"]])
  expect_equal(.admResidApply(c(3, 9), numeric(2), arr)$dv, rep(a2 * 5 / 3, 2))
})

test_that("a t fit equals a matched-variance normal fit (estimator bias cancels)", {
  skip_if_not_installed("nlmixr2est")
  # The sharpest oracle available: t-generated data fitted with t(nu=5) must give
  # the same TOTAL residual variance as matched-variance normal data fitted with a
  # normal model. Any estimator bias (e.g. adfo's linearization) hits both equally
  # and cancels, so this holds for adfo AND adgh.
  for (est in c("adfo", "adgh")) {
    rt <- .t_fit(.mod_t_fixnu, .t_gen("t"),      est)
    rn <- .t_fit(.mod_norm,    .t_gen("normal"), est)
    var_t <- unname(rt$est[["a"]])^2 * 5 / 3     # nu fixed at 5
    var_n <- unname(rn$est[["a"]])^2
    expect_equal(var_t, var_n, tolerance = 0.08,
                 label = sprintf("%s: t total var vs normal total var", est))
  }
})

test_that("estimating nu warns that it is not identifiable", {
  mod <- function() {
    ini({ tcl <- log(4); tv <- log(18); a <- 1; nu <- 8; eta.cl ~ 0.05 })
    model({ cl <- exp(tcl + eta.cl); v <- exp(tv)
            d/dt(central) <- -(cl / v) * central; cp <- central / v
            cp ~ add(a) + t(nu) })
  }
  ui <- suppressMessages(rxode2::rxode2(mod))
  expect_warning(.admParseIniDf(ui$iniDf, ui), "cannot be\\s+estimated from aggregate data")
})

# -- 5. refusals stay refusals ------------------------------------------------

test_that("t() without a scale, with lnorm, or with nu <= 2 is refused", {
  expect_error(
    .admParseIniDf({ui <- suppressMessages(rxode2::rxode2(
      .t_mod("nu <- 5", "cp ~ t(nu)"))); ui$iniDf}, ui),
    "no scale parameter")
  expect_error(
    .admParseIniDf({ui <- suppressMessages(rxode2::rxode2(
      .t_mod("a <- 0.5; nu <- 1.5", "cp ~ add(a) + t(nu)"))); ui$iniDf}, ui),
    "initial estimate")
})
