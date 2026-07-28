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
  propt = c("b <- 0.15; nu <- fix(5)",      "cp ~ prop(b) + t(nu)"),
  # The forms added when the error-model factory was completed. CLAUDE.md's rule
  # is that a new residual form extends THIS file; they were landing with the
  # closed-form layer alone and no analytic-gradient-vs-FD check, which is exactly
  # what let the lnorm sensitivity defect survive for so long.
  bc     = c("a <- 0.3; lam <- 0.4",        "cp ~ add(a) + boxCox(lam)"),
  yj     = c("a <- 0.3; lam <- 0.4",        "cp ~ add(a) + yeoJohnson(lam)"),
  logitn = c("a <- 0.3",                    "cp ~ logitNorm(a, 0, 60)"),
  probn  = c("a <- 0.3",                    "cp ~ probitNorm(a, 0, 60)"),
  pois   = c("",                            "y ~ pois(cp)"),
  nbin   = c("k <- 4",                      "y ~ nbinomMu(k, cp)")
)

# The closed forms written out in the first test below cover the three ANALYTIC
# variance shapes (combined2, combined1, lnorm). The transform-both-sides and
# count families are integrated/derived differently -- they are checked against
# rxode2's own rx_r_ in test-errmodel-rxode2-oracle.R and, here, through their
# derivatives and each estimator's gradient.
.rm_closed <- c("add", "prop", "pow", "comb2", "comb1", "lnorm", "propt")

.rm_pinfo <- function(nm) {
  cs <- .rm_cases[[nm]]
  ui <- .rm_mk(cs[1], cs[2])
  suppressWarnings(admixr2:::.admParseIniDf(ui$iniDf, ui))
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
    arr <- admixr2:::.admResidRows(p, "cp", admixr2:::.admSigmaNat(p$sigma_init, p), length(times))
    ap  <- admixr2:::.admResidApply(mu_f, diag(Cf), arr)
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
  for (nm in .rm_closed) {
    p   <- .rm_pinfo(nm)
    arr <- admixr2:::.admResidRows(p, "cp", admixr2:::.admSigmaNat(p$sigma_init, p), length(mu0))
    ap  <- admixr2:::.admResidApply(mu0, v0, arr)
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
  arr <- admixr2:::.admResidRows(p, "cp", admixr2:::.admSigmaNat(p$sigma_init, p), 3L)
  v0  <- c(1.5, 6.0, 20.0)
  expect_identical(admixr2:::.admResidApply(c(4, 12, 30), v0, arr)$dv, v0 + 0.25)
})

test_that(".admResidDeriv matches central FD, including dv_dv0", {
  mu0 <- c(4, 12, 30); v0 <- c(1.5, 6.0, 20.0); h <- 1e-6
  for (nm in names(.rm_cases)) {
    p  <- .rm_pinfo(nm)
    p0 <- p$sigma_init
    ap <- function(ps, m, vv)
      admixr2:::.admResidApply(m, vv, admixr2:::.admResidRows(p, "cp", admixr2:::.admSigmaNat(ps, p), length(m)))
    d  <- admixr2:::.admResidDeriv(mu0, v0, admixr2:::.admResidRows(p, "cp", admixr2:::.admSigmaNat(p0, p), 3L), p)

    for (k in seq_along(p0)) {
      pp <- p0; pp[k] <- pp[k]+h; pm <- p0; pm[k] <- pm[k]-h
      expect_equal(d$dvar[, k], (ap(pp,mu0,v0)$dv - ap(pm,mu0,v0)$dv)/(2*h),
                   tolerance = 1e-5, label = sprintf("%s d(var)/d(%s)", nm, names(p0)[k]))
      expect_equal(d$dmu[, k], (ap(pp,mu0,v0)$mu - ap(pm,mu0,v0)$mu)/(2*h),
                   tolerance = 1e-5, label = sprintf("%s d(mu)/d(%s)", nm, names(p0)[k]))
    }
    # A LARGER step for the prediction/variance directions. The transform-both-
    # sides forms evaluate their moments by Gauss-Hermite quadrature, so
    # differencing them at 1e-6 sits on the cancellation floor: the reference, not
    # the derivative, is what is wrong there. Verified by Richardson -- the FD
    # converges onto the analytic value as the step grows (max abs difference
    # 3e-5 at h = 1e-6, 1e-7 at h = 1e-4, 4e-8 at h = 1e-3).
    h_f <- 1e-4
    fd_df <- vapply(seq_along(mu0), function(i) { mp<-mu0; mp[i]<-mp[i]+h_f; mm<-mu0; mm[i]<-mm[i]-h_f
      (ap(p0,mp,v0)$dv[i] - ap(p0,mm,v0)$dv[i])/(2*h_f) }, numeric(1))
    fd_v0 <- vapply(seq_along(v0), function(i) { vp<-v0; vp[i]<-vp[i]+h_f; vm<-v0; vm[i]<-vm[i]-h_f
      (ap(p0,mu0,vp)$dv[i] - ap(p0,mu0,vm)$dv[i])/(2*h_f) }, numeric(1))
    expect_equal(d$dv_df,  fd_df, tolerance = 1e-5, label = paste0(nm, ": dv_df"))
    expect_equal(d$dv_dv0, fd_v0, tolerance = 1e-5, label = paste0(nm, ": dv_dv0"))
  }
})

# -- 3. each estimator's analytic gradient vs central FD of its own NLL --------

test_that("analytic gradients match FD for every residual form", {
  skip_if_not_installed("nlmixr2est")
  times <- c(0.5, 1, 2, 4)
  # EVERY form in .rm_cases, not the three this test was written with. The new
  # families each brought their own dv_df / dv_dv0 / dms / rmat gradient path, and
  # a residual whose gradient disagrees with its own objective is precisely the
  # defect class this file exists to catch.
  for (nm in names(.rm_cases)) {
    cs    <- .rm_cases[[nm]]
    ui    <- .rm_mk(cs[1], cs[2])
    pinfo <- suppressWarnings(admixr2:::.admParseIniDf(ui$iniDf, ui))
    sens  <- suppressMessages(tryCatch(admixr2:::.admLoadSensModel(ui), error = function(e) NULL))
    rxMod <- admixr2:::.admLoadModel(ui)
    E <- 100/20*exp(-(5/20)*times)
    s <- admixr2:::.admNormaliseStudy(list(E = E, V = diag((0.3*E)^2), n = 200L,
                                           times = times, ev = rxode2::et(amt = 100)), "s")
    s$ev_full <- rxode2::et(s$ev, s$times)
    st <- list(s = s)
    z  <- admixr2:::.admMakeZ(400L, pinfo, 1L, "sobol")
    pm <- admixr2:::.admMakeParamsList(400L, pinfo, 1L)
    grid <- admixr2:::.adghNodeGrid(5L, pinfo$n_eta)
    p0 <- admixr2:::.admBuildOptVec(pinfo)$p0 * 1.03 + 0.01     # off the optimum

    runs <- list(
      adfo = list(n = function(p) admixr2:::.adfoNLL(p, pinfo, st, sens, rxMod, "cp", pm, 1L),
                  g = function(p) admixr2:::.adfoGrad(p, pinfo, st, sens, rxMod, "cp", pm, 1L, 1e-5)),
      adgh = list(n = function(p) admixr2:::.adghNLL(p, pinfo, st, rxMod, "cp", grid, 1L),
                  g = function(p) admixr2:::.adghGrad(p, pinfo, st, sens, rxMod, "cp", grid, 1L, 1e-5)),
      admc = list(n = function(p) admixr2:::.admNLL(p, pinfo, st, z, rxMod, "cp", pm, 1L),
                  g = function(p) admixr2:::.admGrad(p, pinfo, st, z, rxMod, "cp", pm, 1L, 1e-5, sens))
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

# A residual whose ONLY parameter is fix()ed is prediction-dependent all the same:
# `.all_fixed_resid` keeps the spec alive with n_sig == 0, and Var(y|eta) still
# moves with f (b^2*E[f^2] for a fixed prop). .admResidDeriv used to early-return
# dv_df = 0 / dv_dv0 = 1 for every non-count form when n_sig == 0 -- dropping that
# dependence from BOTH the struct-theta and omega gradients, so the optimizer
# descended a gradient the objective did not follow. Guard: the analytic gradient
# of the struct thetas and omega must still match a central FD of the NLL.
test_that("a fix()ed prediction-dependent residual keeps its gradient chain", {
  skip_if_not_installed("nlmixr2est")
  times <- c(0.5, 1, 2, 4)
  fixed_cases <- list(
    prop  = c("b <- fix(0.2)",              "cp ~ prop(b)"),
    lnorm = c("a <- fix(0.3)",              "cp ~ lnorm(a)"),
    comb  = c("a <- fix(0.4); b <- fix(0.2)", "cp ~ add(a) + prop(b)"))
  for (nm in names(fixed_cases)) {
    cs    <- fixed_cases[[nm]]
    ui    <- .rm_mk(cs[1], cs[2])
    pinfo <- suppressWarnings(admixr2:::.admParseIniDf(ui$iniDf, ui))
    # the premise: no estimated residual parameter, yet the spec is kept
    expect_length(pinfo$sigma_names, 0L)
    sens  <- suppressMessages(tryCatch(admixr2:::.admLoadSensModel(ui), error = function(e) NULL))
    rxMod <- admixr2:::.admLoadModel(ui)
    E <- 100/20*exp(-(5/20)*times)
    s <- admixr2:::.admNormaliseStudy(list(E = E, V = diag((0.3*E)^2), n = 200L,
                                           times = times, ev = rxode2::et(amt = 100)), "s")
    s$ev_full <- rxode2::et(s$ev, s$times)
    st <- list(s = s)
    z    <- admixr2:::.admMakeZ(400L, pinfo, 1L, "sobol")
    pm   <- admixr2:::.admMakeParamsList(400L, pinfo, 1L)
    grid <- admixr2:::.adghNodeGrid(5L, pinfo$n_eta)
    p0   <- admixr2:::.admBuildOptVec(pinfo)$p0 * 1.03 + 0.01
    runs <- list(
      adfo = list(n = function(p) admixr2:::.adfoNLL(p, pinfo, st, sens, rxMod, "cp", pm, 1L),
                  g = function(p) admixr2:::.adfoGrad(p, pinfo, st, sens, rxMod, "cp", pm, 1L, 1e-5)),
      adgh = list(n = function(p) admixr2:::.adghNLL(p, pinfo, st, rxMod, "cp", grid, 1L),
                  g = function(p) admixr2:::.adghGrad(p, pinfo, st, sens, rxMod, "cp", grid, 1L, 1e-5)),
      admc = list(n = function(p) admixr2:::.admNLL(p, pinfo, st, z, rxMod, "cp", pm, 1L),
                  g = function(p) admixr2:::.admGrad(p, pinfo, st, z, rxMod, "cp", pm, 1L, 1e-5, sens)))
    for (est in names(runs)) {
      an <- runs[[est]]$g(p0)
      h  <- 1e-5
      fd <- vapply(seq_along(p0), function(k) {
        pp <- p0; pp[k] <- pp[k]+h; pmn <- p0; pmn[k] <- pmn[k]-h
        (runs[[est]]$n(pp) - runs[[est]]$n(pmn))/(2*h) }, numeric(1))
      tol <- if (est == "adfo") 5e-3 else 1e-3
      expect_lt(max(abs(an - fd)/pmax(abs(fd), 1e-6)), tol,
                label = sprintf("fixed-%s/%s max rel gradient error", nm, est))
    }
  }
})

# yeoJohnson at a large lambda over a LOW prediction pushes the +-12 SD quadrature
# tail nodes past the inverse transform's support -> +-Inf/NaN, and m2 - m*m then
# becomes Inf - Inf = NaN, turning the whole moment, the NLL AND the analytic
# gradient NaN. nloptr tolerates a NaN OBJECTIVE (it rejects the step) but errors on
# a NaN GRADIENT, so the DEFAULT adgh fit crashed mid-optimization. Those tail nodes
# carry GH weight ~1e-30, so .admTBSMomentsD drops the non-finite ones; the moments
# and gradient must stay finite at a parameter the optimizer legitimately reaches.
test_that("TBS moments stay finite when tail nodes overflow (yeoJohnson, large lambda)", {
  ui    <- .rm_mk("a <- 0.15; lam <- 1.3", "cp ~ add(a) + yeoJohnson(lam)")
  pinfo <- suppressWarnings(admixr2:::.admParseIniDf(ui$iniDf, ui))
  sig   <- admixr2:::.admSigmaNat(c(a = 2 * log(0.3), lam = 2.3), pinfo)
  mu    <- c(3, 2, 1, 0.5, 0.2); vf <- (0.3 * mu)^2 * 0.5
  arr   <- admixr2:::.admResidRows(pinfo, rep("cp", length(mu)), sig, length(mu))
  ap    <- admixr2:::.admResidApply(mu, vf, arr)
  d     <- admixr2:::.admResidDeriv(mu, vf, arr, pinfo)
  expect_true(all(is.finite(ap$mu)),  label = "TBS moment mean finite")
  expect_true(all(is.finite(ap$dv)),  label = "TBS moment variance finite")
  expect_true(all(is.finite(d$dvar)) && all(is.finite(d$dmu)) &&
              all(is.finite(d$dv_df)) && all(is.finite(d$dv_dv0)),
              label = "TBS derivative finite")
})

# End to end: the DEFAULT adgh gradient is analytical, and it used to die on the
# above. A yeoJohnson fit must now run to completion under it.
test_that("adgh fits a yeoJohnson endpoint under its default (analytical) gradient", {
  skip_if_not_installed("nlmixr2est")
  set.seed(11L)
  tm <- c(0.25, 0.5, 1, 2, 4, 8, 12); n <- 400L
  eta <- rnorm(n, 0, sqrt(0.09)); cl <- exp(log(3) + eta); v <- exp(log(30))
  f  <- outer(cl, tm, function(cl_, t_) 100 / v * exp(-(cl_ / v) * t_))
  yt <- admixr2:::.admTBS(f, 1.3, 1L, 0, 1) + 0.15 * matrix(rnorm(length(f)), n, length(tm))
  y  <- admixr2:::.admTBSi(yt, 1.3, 1L, 0, 1)
  study <- list(E = colMeans(y), V = cov.wt(y, method = "ML")$cov, n = n,
                times = tm, ev = rxode2::et(amt = 100))
  mod <- function() {
    ini({ tcl <- log(3); tv <- log(30); a <- 0.15; lam <- 1.3; eta.cl ~ 0.09 })
    model({ cl <- exp(tcl + eta.cl); v <- exp(tv)
            d/dt(central) <- -(cl/v)*central; cp <- central/v
            cp ~ add(a) + yeoJohnson(lam) }) }
  fit <- suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
    mod, admData(), est = "adgh", control = adghControl(studies = list(s = study),
                                                        covMethod = "none"))))
  expect_s3_class(fit, "admFit")
  expect_true(is.finite(fit$objective))
  expect_equal(unname(fit$parFixedDf["lam", "Estimate"]), 1.3, tolerance = 0.25)
})
