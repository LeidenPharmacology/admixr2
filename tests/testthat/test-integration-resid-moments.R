# Aggregate residual moments and their gradients ------------------------------
#
# admixr2 composes the residual with the eta-integral by the LAW OF TOTAL VARIANCE:
#
#   mu_pred = ms * E[f]
#   V_pred  = ms^2 (x) Cov_eta(f) + diag( E_eta[Var(y|eta)] )
#
# It used to compute `Cov_eta(f) + Sigma(mu)` instead -- evaluating the residual at
# the population mean rather than averaging it over individual predictions. That is
# exact only for ADDITIVE error; for prop/pow it understates diag(V_pred) by
# b^2*Var_eta(f), and for lnorm it additionally drops an exp(s) factor from every
# OFF-diagonal. Measured against individual-level simulation, the old formulas
# carried fixed 15-20% biases that did not shrink with sample size.
#
# Three layers of check, each against an INDEPENDENT reference:
#   1. closed-form moments vs a large individual-level Monte Carlo
#   2. .admResidApply / .admResidDeriv vs closed form and central FD
#   3. each estimator's ANALYTIC gradient vs central FD of its own NLL
#
# (3) is the one that matters most: it is what was missing. `lnorm` appeared in no
# gradient test at all, which let two independent defects survive --
#   * the sensitivity model returns rx_pred_ = log(f) for an lnorm endpoint while
#     .admSimulate returns f, so .admGrad differentiated a different function; and
#   * the d(V_pred)/d(V_struct) chain was assumed to be the identity.

skip_on_cran()
skip_if_not_installed("rxode2")

.rm_mk <- function(ini, err) suppressMessages(rxode2::rxode2(eval(parse(text = sprintf(
  'function() { ini({ tcl <- log(5); tv <- log(20); %s
                      eta.cl ~ 0.09; eta.v ~ 0.04 })
     model({ cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
             d/dt(central) <- -(cl/v)*central; cp <- central/v
             %s }) }', ini, err)))))

.rm_cases <- list(
  add   = c("a <- 0.5",                     "cp ~ add(a)"),
  prop  = c("b <- 0.15",                    "cp ~ prop(b)"),
  pow   = c("b <- 0.15; c1 <- 1.3",         "cp ~ pow(b, c1)"),
  comb2 = c("a <- 0.5; b <- 0.15",          "cp ~ add(a) + prop(b)"),
  comb1 = c("a <- 0.5; b <- 0.15",          "cp ~ add(a) + prop(b) + combined1()"),
  lnorm = c("a <- 0.3",                     "cp ~ lnorm(a)"),
  propt = c("b <- 0.15; nu <- fix(5)",      "cp ~ prop(b) + t(nu)")
)

.rm_pinfo <- function(nm) {
  cs <- .rm_cases[[nm]]
  ui <- .rm_mk(cs[1], cs[2])
  suppressWarnings(.admParseIniDf(ui$iniDf, ui))
}

# -- 1. the moment formulas vs an individual-level simulation ------------------

test_that("V_pred is the law of total variance, not Sigma at the mean", {
  skip_on_os("solaris")
  set.seed(202L)
  times <- c(1, 2, 4, 8); N <- 400000L
  tcl <- log(5); tv <- log(20); om <- 0.16          # 40% CV, deliberately large
  eta <- stats::rnorm(N, 0, sqrt(om))
  f   <- vapply(times, function(tt) 1000/exp(tv)*exp(-exp(tcl+eta)/exp(tv)*tt),
                numeric(N))
  mu_f <- colMeans(f); Cf <- crossprod(sweep(f, 2L, mu_f))/N
  Z <- function() matrix(stats::rnorm(N*length(times)), N, length(times))

  gen <- list(
    prop  = function() f + 0.15*f*Z(),
    lnorm = function() f * exp(sqrt(0.09)*Z())
  )
  for (nm in names(gen)) {
    p   <- .rm_pinfo(nm)
    arr <- .admResidRows(p, "cp", .admSigmaNat(p$sigma_init, p), length(times))
    ap  <- .admResidApply(mu_f, diag(Cf), arr)
    V   <- Cf * tcrossprod(ap$ms); diag(V) <- ap$dv

    y  <- gen[[nm]]()
    Ve <- crossprod(sweep(y, 2L, colMeans(y)))/N
    off <- upper.tri(Ve)
    expect_equal(as.numeric(ap$mu), colMeans(y), tolerance = 0.01,
                 label = paste0(nm, ": mu_pred"))
    expect_equal(diag(V), diag(Ve), tolerance = 0.05,
                 label = paste0(nm, ": diag(V_pred)"))
    expect_equal(V[off], Ve[off], tolerance = 0.05,
                 label = paste0(nm, ": offdiag(V_pred)"))
  }
})

# -- 2. .admResidApply / .admResidDeriv ---------------------------------------

test_that(".admResidApply matches the closed form for every residual form", {
  mu0 <- c(4, 12, 30); v0 <- c(1.5, 6.0, 20.0)
  Ef  <- function(k) mu0^k + k*(k-1)/2 * mu0^(k-2) * v0
  for (nm in names(.rm_cases)) {
    p   <- .rm_pinfo(nm)
    arr <- .admResidRows(p, "cp", .admSigmaNat(p$sigma_init, p), length(mu0))
    ap  <- .admResidApply(mu0, v0, arr)
    a2 <- arr$a2; b2 <- arr$b2; cc <- arr$cc
    if (arr$form[1] == 2L) {                    # lnorm
      ms <- exp(a2/2); ev <- Ef(2) * exp(a2) * (exp(a2) - 1)
    } else if (arr$form[1] == 1L) {              # combined1
      ms <- rep(1, 3); ev <- a2 + 2*sqrt(a2*b2)*Ef(cc) + b2*Ef(2*cc)
    } else {                                     # combined2
      ms <- rep(1, 3); ev <- a2 + b2*Ef(2*cc)
    }
    expect_equal(as.numeric(ap$mu), ms*mu0, label = paste0(nm, ": mu"))
    expect_equal(as.numeric(ap$dv), ms^2*v0 + ev, label = paste0(nm, ": dv"))
  }
})

test_that("additive error is bit-identical to v0 + a^2", {
  p   <- .rm_pinfo("add")
  arr <- .admResidRows(p, "cp", .admSigmaNat(p$sigma_init, p), 3L)
  v0  <- c(1.5, 6.0, 20.0)
  expect_identical(.admResidApply(c(4, 12, 30), v0, arr)$dv, v0 + 0.25)
})

test_that(".admResidDeriv matches central FD, including dv_dv0", {
  mu0 <- c(4, 12, 30); v0 <- c(1.5, 6.0, 20.0); h <- 1e-6
  for (nm in names(.rm_cases)) {
    p  <- .rm_pinfo(nm)
    p0 <- p$sigma_init
    ap <- function(ps, m, vv)
      .admResidApply(m, vv, .admResidRows(p, "cp", .admSigmaNat(ps, p), length(m)))
    d  <- .admResidDeriv(mu0, v0, .admResidRows(p, "cp", .admSigmaNat(p0, p), 3L), p)

    for (k in seq_along(p0)) {
      pp <- p0; pp[k] <- pp[k]+h; pm <- p0; pm[k] <- pm[k]-h
      expect_equal(d$dvar[, k], (ap(pp,mu0,v0)$dv - ap(pm,mu0,v0)$dv)/(2*h),
                   tolerance = 1e-5, label = sprintf("%s d(var)/d(%s)", nm, names(p0)[k]))
      expect_equal(d$dmu[, k], (ap(pp,mu0,v0)$mu - ap(pm,mu0,v0)$mu)/(2*h),
                   tolerance = 1e-5, label = sprintf("%s d(mu)/d(%s)", nm, names(p0)[k]))
    }
    fd_df <- vapply(seq_along(mu0), function(i) { mp<-mu0; mp[i]<-mp[i]+h; mm<-mu0; mm[i]<-mm[i]-h
      (ap(p0,mp,v0)$dv[i] - ap(p0,mm,v0)$dv[i])/(2*h) }, numeric(1))
    fd_v0 <- vapply(seq_along(v0), function(i) { vp<-v0; vp[i]<-vp[i]+h; vm<-v0; vm[i]<-vm[i]-h
      (ap(p0,mu0,vp)$dv[i] - ap(p0,mu0,vm)$dv[i])/(2*h) }, numeric(1))
    expect_equal(d$dv_df,  fd_df, tolerance = 1e-5, label = paste0(nm, ": dv_df"))
    expect_equal(d$dv_dv0, fd_v0, tolerance = 1e-5, label = paste0(nm, ": dv_dv0"))
  }
})

# -- 3. each estimator's analytic gradient vs central FD of its own NLL --------

test_that("analytic gradients match FD for add, prop and lnorm", {
  skip_if_not_installed("nlmixr2est")
  times <- c(0.5, 1, 2, 4)
  for (nm in c("add", "prop", "lnorm")) {
    cs    <- .rm_cases[[nm]]
    ui    <- .rm_mk(cs[1], cs[2])
    pinfo <- suppressWarnings(.admParseIniDf(ui$iniDf, ui))
    sens  <- suppressMessages(tryCatch(.admLoadSensModel(ui), error = function(e) NULL))
    rxMod <- .admLoadModel(ui)
    E <- 100/20*exp(-(5/20)*times)
    s <- .admNormaliseStudy(list(E = E, V = diag((0.3*E)^2), n = 200L,
                                 times = times, ev = rxode2::et(amt = 100)), "s")
    s$ev_full <- rxode2::et(s$ev, s$times)
    st <- list(s = s)
    z  <- .admMakeZ(400L, pinfo, 1L, "sobol")
    pm <- .admMakeParamsList(400L, pinfo, 1L)
    grid <- .adghNodeGrid(5L, pinfo$n_eta)
    p0 <- .admBuildOptVec(pinfo)$p0 * 1.03 + 0.01     # off the optimum

    runs <- list(
      adfo = list(n = function(p) .adfoNLL(p, pinfo, st, sens, rxMod, "cp", pm, 1L),
                  g = function(p) .adfoGrad(p, pinfo, st, sens, rxMod, "cp", pm, 1L, 1e-5)),
      adgh = list(n = function(p) .adghNLL(p, pinfo, st, rxMod, "cp", grid, 1L),
                  g = function(p) .adghGrad(p, pinfo, st, sens, rxMod, "cp", grid, 1L, 1e-5)),
      admc = list(n = function(p) .admNLL(p, pinfo, st, z, rxMod, "cp", pm, 1L),
                  g = function(p) .admGrad(p, pinfo, st, z, rxMod, "cp", pm, 1L, 1e-5, sens))
    )
    for (est in names(runs)) {
      an <- runs[[est]]$g(p0)
      h  <- 1e-5
      fd <- vapply(seq_along(p0), function(k) {
        pp <- p0; pp[k] <- pp[k]+h; pmn <- p0; pmn[k] <- pmn[k]-h
        (runs[[est]]$n(pp) - runs[[est]]$n(pmn))/(2*h) }, numeric(1))
      # adfo finite-differences its own struct thetas, so it is FD-vs-FD and
      # noisier than the exact-sensitivity estimators.
      tol <- if (est == "adfo") 5e-3 else 1e-3
      expect_lt(max(abs(an - fd)/pmax(abs(fd), 1e-6)), tol,
                label = sprintf("%s/%s max rel gradient error", nm, est))
    }
  }
})
