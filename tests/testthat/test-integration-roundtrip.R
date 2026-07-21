# THE ROUND TRIP: simulate from a fitted model and you get back what was fitted.
#
# admixr2 fits a study mean E and covariance V. The property a user relies on is
# that those are the moments of the model they now hold: take the parameters,
# simulate individuals with rxode2, aggregate, and recover the same E and V.
#
# This is a stronger statement than "the conditional moments are right" (which the
# oracle tests pin down with eta = 0) and than "the parameters are recovered"
# (simulation-estimation). It composes the residual model, the law of total
# variance and the eta integration all at once, and it checks the OFF-DIAGONALS --
# which is where lnorm's ms(x)ms scaling and ordinal's cross-category term live.
#
# rxode2's own `sim` column is the reference: never a hand-rolled residual, so the
# test cannot agree with admixr2 by sharing a mistake.
#
# NOT covered here: ar(). rxode2's ar() SIMULATION is not stationary when a dose
# record precedes the first observation (a zero-amount dose reproduces it, and
# nlmixr2's own focei cannot recover rho from it either), while its ESTIMATION
# lines imply the stationary process admixr2 scores. See .admBuildResidSpecs().

.rt_cases <- list(
  add        = list("a <- 0.5",                  "cp ~ add(a)"),
  prop       = list("b <- 0.15",                 "cp ~ prop(b)"),
  combined2  = list("a <- 0.3; b <- 0.10",       "cp ~ add(a) + prop(b)"),
  combined1  = list("a <- 0.3; b <- 0.10",       "cp ~ add(a) + prop(b) + combined1()"),
  lnorm      = list("a <- 0.20",                 "cp ~ lnorm(a)"),
  t_fix      = list("a <- 0.5; nu <- fix(6)",    "cp ~ add(a) + t(nu)"),
  boxCox     = list("a <- 0.3; lam <- fix(0.5)", "cp ~ add(a) + boxCox(lam)"),
  bc_prop    = list("a <- 0.3; b <- 0.1; lam <- fix(0.5)",
                    "cp ~ add(a) + prop(b) + boxCox(lam)"),
  logitNorm  = list("a <- 0.25",                 "cp ~ logitNorm(a, 0, 12)")
)

test_that("predicted (E, V) equals simulate-from-the-model-and-aggregate", {
  skip_on_cran()
  skip_if_not_installed("rxode2")
  skip_if_not_installed("lotri")

  OM <- 0.09; N <- 40000L
  times <- c(0.5, 2, 6, 12)
  ev <- rxode2::et(amt = 100)

  for (nm in names(.rt_cases)) {
    cs  <- .rt_cases[[nm]]
    src <- sprintf('function() { ini({ tcl <- log(1); tv <- log(20); %s
                                       eta.cl ~ %g })
                     model({ cl <- exp(tcl + eta.cl); v <- exp(tv)
                             d/dt(central) <- -(cl/v)*central; cp <- central/v
                             %s }) }', cs[[1L]], OM, cs[[2L]])
    ui    <- suppressMessages(rxode2::rxode2(eval(parse(text = src))))
    pinfo <- suppressWarnings(.admParseIniDf(ui$iniDf, ui))
    rxMod <- .admLoadModel(ui)

    s <- .admNormaliseStudy(list(E = rep(0, length(times)), V = diag(length(times)),
                                 n = 100L, times = times, ev = ev), "s")
    s$ev_full <- rxode2::et(s$ev, s$times)
    pars <- .admUnpack(.admBuildOptVec(pinfo)$p0, pinfo)
    grid <- .adghNodeGrid(15L, pinfo$n_eta)
    m    <- .adghMoments(pars, pinfo, s, rxMod, .admOutputVar(ui), grid, 1L)

    set.seed(9)
    sim <- suppressWarnings(rxode2::rxSolve(
      ui$simulationModel, rxode2::et(ev, times), nSub = N,
      omega = lotri::lotri(eta.cl ~ OM), sigma = lotri::lotri(rxerr.cp ~ 1),
      returnType = "data.frame", nDisplayProgress = .Machine$integer.max))
    k  <- sim$time %in% times
    Y  <- matrix(sim$sim[k], N, length(times), byrow = TRUE)
    Ee <- colMeans(Y); Ve <- stats::cov.wt(Y, method = "ML")$cov

    # Judged against the simulation's OWN Monte-Carlo error, so the tolerance is
    # a statistical statement rather than a number chosen to make it pass.
    sdv <- apply(Y, 2, stats::sd)
    seE <- sdv / sqrt(N)
    seV <- outer(sdv, sdv) * sqrt(2 / N)
    od  <- upper.tri(Ve)

    expect_lt(max(abs(m$E - Ee) / seE),                    6, label = paste(nm, "mean"))
    expect_lt(max(abs(diag(m$V) - diag(Ve)) / diag(seV)),  6, label = paste(nm, "var"))
    expect_lt(max(abs((m$V - Ve)[od]) / seV[od]),          6, label = paste(nm, "cov"))
  }
})

test_that("the round trip holds for a depot model, where f == 0 at t = 0", {
  skip_on_cran()
  skip_if_not_installed("rxode2")
  skip_if_not_installed("lotri")
  # A structural prediction of exactly 0 is where the moment expansion has a pole
  # for pow(c < 1), and where the C++/R kernels disagreed until adm_mom_f was
  # guarded. Keep an observation AT t = 0 so the zero is actually exercised.
  OM <- 0.09; N <- 40000L
  times <- c(0, 0.5, 2, 6, 12)
  ev <- rxode2::et(amt = 100, cmt = "depot")
  for (cs in list(list("a <- 0.3; b <- 0.1", "cp ~ add(a) + prop(b)"),
                  list("a <- 0.3; b <- 0.15; c1 <- 0.4", "cp ~ add(a) + pow(b, c1)"))) {
    src <- sprintf('function() { ini({ tka <- log(1.2); tcl <- log(1); tv <- log(20)
                                       %s; eta.cl ~ %g })
                     model({ ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
                             d/dt(depot) <- -ka*depot
                             d/dt(central) <- ka*depot - (cl/v)*central
                             cp <- central/v
                             %s }) }', cs[[1L]], OM, cs[[2L]])
    ui    <- suppressMessages(rxode2::rxode2(eval(parse(text = src))))
    pinfo <- suppressWarnings(.admParseIniDf(ui$iniDf, ui))
    rxMod <- .admLoadModel(ui)
    s <- .admNormaliseStudy(list(E = rep(0, length(times)), V = diag(length(times)),
                                 n = 100L, times = times, ev = ev), "s")
    s$ev_full <- rxode2::et(s$ev, s$times)
    pars <- .admUnpack(.admBuildOptVec(pinfo)$p0, pinfo)
    m <- .adghMoments(pars, pinfo, s, rxMod, .admOutputVar(ui),
                      .adghNodeGrid(15L, pinfo$n_eta), 1L)
    expect_true(all(is.finite(m$E)))
    expect_true(all(is.finite(m$V)))
    expect_gte(min(diag(m$V)), 0)

    set.seed(9)
    sim <- suppressWarnings(rxode2::rxSolve(
      ui$simulationModel, rxode2::et(ev, times), nSub = N,
      omega = lotri::lotri(eta.cl ~ OM), sigma = lotri::lotri(rxerr.cp ~ 1),
      returnType = "data.frame", nDisplayProgress = .Machine$integer.max))
    k  <- sim$time %in% times
    Y  <- matrix(sim$sim[k], N, length(times), byrow = TRUE)
    sdv <- apply(Y, 2, stats::sd)
    expect_lt(max(abs(m$E - colMeans(Y)) / (sdv / sqrt(N))), 6)
    expect_lt(max(abs(diag(m$V) - diag(stats::cov.wt(Y, method = "ML")$cov)) /
                    (sdv^2 * sqrt(2 / N))), 6)
  }
})
