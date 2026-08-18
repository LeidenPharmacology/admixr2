# Covariate marginalisation, against an INDEPENDENT reference.
#
# Every check here compares the package against exact nested Gauss-Hermite
# quadrature computed in plain R on an analytic 1-cmt solution, or compares two
# independent code paths against each other. Nothing is pinned against its own
# output: the two bugs this feature shipped during development (a discarded
# return value that disabled the whole thing, and a covariate grid that dropped
# the distribution's tails) both produced finite, plausible numbers, so
# "the fit ran and the estimates look sane" catches neither.

skip_on_cran()
skip_if_not_installed("rxode2")

.cov_TCL <- log(1.0); .cov_TV <- log(10); .cov_OM <- 0.30
.cov_ADD <- 0.30; .cov_DOSE <- 100; .cov_TIMES <- c(0.5, 1, 2, 4, 8)

.cov_conc <- function(cl)
  .cov_DOSE / exp(.cov_TV) * exp(outer(-cl / exp(.cov_TV), .cov_TIMES))

# probabilists' Gauss-Hermite (Golub-Welsch), for the reference only
.cov_gh <- function(k) {
  i <- seq_len(k - 1L)
  J <- matrix(0, k, k); J[cbind(i, i + 1L)] <- sqrt(i); J[cbind(i + 1L, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE); o <- order(e$values)
  list(x = e$values[o], w = e$vectors[1L, o]^2)
}

# Exact structural moments over (WT, eta); residual NOT added.
.cov_ref <- function(cl_of, meanlog, sdlog, k = 40L) {
  q <- .cov_gh(k); m1 <- 0; M2 <- 0
  for (ia in seq_along(q$x)) {
    wt <- exp(meanlog + sdlog * q$x[ia])
    for (ib in seq_along(q$x)) {
      w <- q$w[ia] * q$w[ib]
      Y <- as.numeric(.cov_conc(cl_of(wt, .cov_OM * q$x[ib])))
      m1 <- m1 + w * Y; M2 <- M2 + w * outer(Y, Y)
    }
  }
  list(E = m1, V = M2 - outer(m1, m1))
}

.cov_allometric <- function() {
  ini({tcl <- log(1.0); tv <- log(10); tcov <- 0.75
       eta.cl ~ 0.09; add.err <- 0.3})
  model({cl <- exp(tcl + eta.cl) * (WT / 70)^tcov
         v  <- exp(tv)
         d/dt(centr) <- -cl / v * centr
         cp <- centr / v
         cp ~ add(add.err)})
}
.cov_linear <- function() {
  ini({tcl <- log(1.0); tv <- log(10); tcov <- 0.75
       eta.cl ~ 0.09; add.err <- 0.3})
  model({cl <- exp(tcl + tcov * WT + eta.cl)
         v  <- exp(tv)
         d/dt(centr) <- -cl / v * centr
         cp <- centr / v
         cp ~ add(add.err)})
}

# Build the internal state the NLL works from, for one study.
.cov_setup <- function(fn, cov_dist, cov_ref, E, V) {
  ui    <- suppressMessages(rxode2::rxode2(fn))
  pinfo <- admixr2:::.admParseIniDf(ui$iniDf, ui)
  pinfo$nDisplayProgress <- .Machine$integer.max
  s <- list(E = E, V = V, n = 300L, times = .cov_TIMES,
            ev = rxode2::et(amt = .cov_DOSE),
            cov = cov_ref, cov_dist = cov_dist)
  st <- admixr2:::.admFlattenStudies(
          list(s1 = admixr2:::.admNormaliseStudy(s, "s1", "cp")))
  st <- admixr2:::.admBuildEvFull(st)
  st <- admixr2:::.admCheckCovariates(ui, pinfo, st, "none")
  ov <- admixr2:::.admBuildOptVec(pinfo)
  list(ui = ui, pinfo = pinfo, st = st, ov = ov,
       rxMod = admixr2:::.admLoadModel(ui),
       pars = admixr2:::.admUnpack(ov$p0, pinfo))
}

# Predicted STRUCTURAL moments through the package's own covariate path.
.cov_pred <- function(d, n_sim = 6000L) {
  z <- admixr2:::.admMakeZ(n_sim, d$pinfo, 1L, "sobol")[[1L]]
  if (!is.matrix(z)) z <- matrix(z, ncol = 1L)
  eta <- z %*% t(admixr2:::.admStudyL(d$pars, d$pinfo, d$st[[1L]]))
  colnames(eta) <- d$pinfo$eta_col_names
  # Dispatch on the SAME field .admNLL dispatches on. An earlier version of this
  # helper handled only "uq" and silently held the covariate at its reference on
  # the "rows" path -- i.e. it reproduced the very bug the tests exist to catch.
  su <- d$st[[1L]]
  if (identical(su$.adm_cov_path, "rows")) {
    su <- admixr2:::.admStudyCovRows(su, d$pinfo, nrow(eta))
    expect_false(is.null(su$cov_rows))
  }
  pm <- admixr2:::.admMakeParamsList(n_sim, d$pinfo, 1L)[[1L]]
  cp <- admixr2:::.admSimulate(d$rxMod, d$pars$struct, d$pinfo$sigma_names, eta,
                               su, "cp", pm, 1L,
                               .Machine$integer.max, NULL)
  mu <- colMeans(cp)
  list(E = mu, V = crossprod(sweep(cp, 2L, mu)) / n_sim)
}

test_that(".admShiftDelta measures Delta(a) from the model, exactly", {
  # The replacement for the retired .admCovDelta, which measured the same shift
  # with an extra rxSolve per objective evaluation. .admShiftDelta evaluates the
  # model's OWN parameter assignment in R instead, so it costs no solve at all
  # -- and it is the quantity .admShiftVerify() then checks against the compiled
  # model, which is what makes the shift path admissible where uq was not.
  cd <- list(WT = list(meanlog = log(72), sdlog = 0.28))
  d  <- .cov_setup(.cov_allometric, cd, list(WT = 72),
                   rep(1, length(.cov_TIMES)), diag(length(.cov_TIMES)))
  sp <- admixr2:::.admShiftSpec(d$ui, "WT", d$pinfo$eta_col_names)
  expect_false(is.null(sp))
  expect_identical(sp$param, "cl")
  expect_identical(sp$eta, "eta.cl")

  a  <- exp(log(72) + 0.28 * admixr2:::.adghNodes1(15L)$x)
  X  <- matrix(a, ncol = 1L, dimnames = list(NULL, "WT"))
  D  <- admixr2:::.admShiftDelta(sp, admixr2:::.admShiftStruct(d$pinfo),
                                 X, list(WT = 72))
  # the model is cl = exp(tcl+eta)*(WT/70)^tcov, so Delta(a) = tcov*log(a/72)
  # measured against the reference value of 72
  expect_equal(as.numeric(D), 0.75 * log(a / 72), tolerance = 1e-12)
})

test_that("the general path reproduces exact nested quadrature for an ALLOMETRIC effect", {
  ml <- log(72); sl <- 0.28
  ref <- .cov_ref(function(wt, eta) exp(.cov_TCL + eta) * (wt / 70)^0.75, ml, sl)
  Vo  <- ref$V; diag(Vo) <- diag(Vo) + .cov_ADD^2
  d <- .cov_setup(.cov_allometric, list(WT = list(meanlog = ml, sdlog = sl)),
                  list(WT = exp(ml)), ref$E, Vo)
  expect_identical(d$st[[1L]]$.adm_cov_path, "rows")
  p <- .cov_pred(d)
  expect_lt(max(abs(p$E / ref$E - 1)), 5e-3)
  expect_lt(max(abs(p$V / ref$V - 1)), 3e-2)
})

test_that("the general path matches the retired collapse's closed form", {
  # `cl <- exp(tcl + tcov*WT + eta.cl)` with a NORMAL WT is the one and only
  # configuration the retired "collapse" path covered. Its closed form is
  #   Omega* = Omega + tcov^2 * Var(WT),
  # solved at the covariate mean -- so the general path, which knows none of
  # that, must reproduce exactly those moments. This is the check that the
  # capability was not lost with the code.
  mu <- 0.20; sd <- 0.35
  d <- .cov_setup(.cov_linear, list(WT = list(mu = mu, sd = sd)),
                  list(WT = mu), rep(1, length(.cov_TIMES)),
                  diag(length(.cov_TIMES)))
  expect_identical(d$st[[1L]]$.adm_cov_path, "rows")
  p <- .cov_pred(d, n_sim = 12000L)

  # the closed form, computed here rather than by the package
  om2 <- .cov_OM^2 + 0.75^2 * sd^2
  ref <- .cov_ref(function(wt, eta) exp(.cov_TCL + 0.75 * mu + eta),
                  0, 0)                      # covariate spread folded into eta
  q   <- .cov_gh(60L); m1 <- 0; M2 <- 0
  for (i in seq_along(q$x)) {
    Y  <- as.numeric(.cov_conc(exp(.cov_TCL + 0.75 * mu + sqrt(om2) * q$x[i])))
    m1 <- m1 + q$w[i] * Y; M2 <- M2 + q$w[i] * outer(Y, Y)
  }
  Vc <- M2 - outer(m1, m1)
  expect_lt(max(abs(p$E / m1 - 1)), 5e-3)
  expect_lt(max(abs(p$V / Vc - 1)), 3e-2)
})

test_that("a declared covariate distribution is never silently ignored", {
  # If the covariate path is disabled, the predicted covariance must CHANGE.
  # A discarded return value once turned the whole feature into a no-op while
  # every reported number stayed plausible; this is the check that catches it.
  ml <- log(72); sl <- 0.28
  ref <- .cov_ref(function(wt, eta) exp(.cov_TCL + eta) * (wt / 70)^0.75, ml, sl)
  Vo  <- ref$V; diag(Vo) <- diag(Vo) + .cov_ADD^2
  d <- .cov_setup(.cov_allometric, list(WT = list(meanlog = ml, sdlog = sl)),
                  list(WT = exp(ml)), ref$E, Vo)
  p_on <- .cov_pred(d)
  d$st[[1L]]$.adm_cov_path <- NULL        # covariate held at its mean
  p_off <- .cov_pred(d)
  expect_gt(max(abs(diag(p_on$V) / diag(p_off$V) - 1)), 0.05)
})

# ---- the general path: structures the shift assumption cannot represent ------

.cov_ref2 <- function(cl_of, v_of, meanlog, sdlog, k = 40L) {
  q <- .cov_gh(k); m1 <- 0; M2 <- 0
  for (ia in seq_along(q$x)) {
    wt <- exp(meanlog + sdlog * q$x[ia])
    for (ib in seq_along(q$x)) {
      w  <- q$w[ia] * q$w[ib]
      cl <- cl_of(wt, .cov_OM * q$x[ib]); v <- v_of(wt)
      Y  <- as.numeric(.cov_DOSE / v * exp(-cl / v * .cov_TIMES))
      m1 <- m1 + w * Y; M2 <- M2 + w * outer(Y, Y)
    }
  }
  list(E = m1, V = M2 - outer(m1, m1))
}

.cov_both <- function() {
  ini({tcl <- log(1.0); tv <- log(10); pcl <- 0.75; pv <- 1.0
       eta.cl ~ 0.09; add.err <- 0.3})
  model({cl <- exp(tcl + eta.cl) * (WT / 70)^pcl
         v  <- exp(tv) * (WT / 70)^pv
         d/dt(centr) <- -cl / v * centr
         cp <- centr / v
         cp ~ add(add.err)})
}
.cov_noeta <- function() {
  ini({tcl <- log(1.0); tv <- log(10); pv <- 1.0
       eta.cl ~ 0.09; add.err <- 0.3})
  model({cl <- exp(tcl + eta.cl)
         v  <- exp(tv) * (WT / 70)^pv
         d/dt(centr) <- -cl / v * centr
         cp <- centr / v
         cp ~ add(add.err)})
}

test_that("a covariate on SEVERAL parameters goes through the general path", {
  ml <- log(72); sl <- 0.28
  r  <- .cov_ref2(function(wt, e) exp(.cov_TCL + e) * (wt / 70)^0.75,
                  function(wt) exp(.cov_TV) * (wt / 70)^1.0, ml, sl)
  Vo <- r$V; diag(Vo) <- diag(Vo) + .cov_ADD^2
  d  <- .cov_setup(.cov_both, list(WT = list(meanlog = ml, sdlog = sl)),
                   list(WT = exp(ml)), r$E, Vo)
  # WT appears twice, so neither the collapse nor u-quantile can represent it:
  # u-quantile would freeze the effect on `v` at the reference value (measured
  # V 99.5% wrong before the general path existed).
  expect_identical(d$st[[1L]]$.adm_cov_path, "rows")
  p <- .cov_pred(d, n_sim = 12000L)
  expect_lt(max(abs(p$E / r$E - 1)), 5e-3)
  expect_lt(max(abs(p$V / r$V - 1)), 3e-2)
})

test_that("a covariate on a parameter with NO random effect is supported", {
  ml <- log(72); sl <- 0.28
  r  <- .cov_ref2(function(wt, e) exp(.cov_TCL + e),
                  function(wt) exp(.cov_TV) * (wt / 70)^1.0, ml, sl)
  Vo <- r$V; diag(Vo) <- diag(Vo) + .cov_ADD^2
  d  <- .cov_setup(.cov_noeta, list(WT = list(meanlog = ml, sdlog = sl)),
                   list(WT = exp(ml)), r$E, Vo)
  # there is no eta on `v` to carry a shift at all -- only the general path can
  expect_identical(d$st[[1L]]$.adm_cov_path, "rows")
  p <- .cov_pred(d, n_sim = 12000L)
  expect_lt(max(abs(p$E / r$E - 1)), 5e-3)
  expect_lt(max(abs(p$V / r$V - 1)), 3e-2)
})

test_that("the general path is the default route for every covariate form", {
  # Neither "collapse" (a bare theta*COV product, a NORMAL covariate, the study
  # solved at the covariate mean, and grad = "none") nor "uq" (four conditions
  # inferred from the model TEXT, each measured to be silently wrong when
  # assumed and false) is routed to any more. cov_integration = "shift"/"auto"
  # is the only route off "rows", and it is admitted numerically.
  ml <- log(72); sl <- 0.28
  E0 <- rep(1, length(.cov_TIMES)); V0 <- diag(length(.cov_TIMES))
  expect_identical(
    .cov_setup(.cov_linear, list(WT = list(mu = 0.2, sd = 0.35)),
               list(WT = 0.2), E0, V0)$st[[1L]]$.adm_cov_path, "rows")
  expect_identical(
    .cov_setup(.cov_allometric, list(WT = list(meanlog = ml, sdlog = sl)),
               list(WT = exp(ml)), E0, V0)$st[[1L]]$.adm_cov_path, "rows")
  expect_identical(
    .cov_setup(.cov_linear, list(WT = list(values = c(0, 1), probs = c(.6, .4))),
               list(WT = 0.4), E0, V0)$st[[1L]]$.adm_cov_path, "rows")
})

# ---- gradients on the general path -------------------------------------------

test_that("the general path supports ANALYTIC gradients (vs central FD)", {
  # The covariate is DATA here -- a per-row params column -- so the existing
  # sensitivity directions already differentiate the function the NLL evaluates:
  # a covariate coefficient is an unpaired struct theta with its own THETA_j
  # direction, and the eta draws are untouched. Nothing new is derived, but the
  # FINITE-DIFFERENCE frames in .admGrad build their params matrix by hand and
  # had to be given the covariate columns explicitly (tiled per block, or the
  # difference stops being common-random-numbers).
  ml <- log(72); sl <- 0.28
  r  <- .cov_ref2(function(wt, e) exp(.cov_TCL + e) * (wt / 70)^0.75,
                  function(wt) exp(.cov_TV) * (wt / 70)^1.0, ml, sl)
  Vo <- r$V; diag(Vo) <- diag(Vo) + .cov_ADD^2
  st0 <- list(s1 = list(E = r$E, V = Vo, n = 300L, times = .cov_TIMES,
                        ev = rxode2::et(amt = .cov_DOSE),
                        cov_dist = list(WT = list(meanlog = ml, sdlog = sl))))
  ui <- suppressMessages(rxode2::rxode2(.cov_both))
  ovar <- admixr2:::.admOutputVar(ui)

  for (g in c("sens", "fd")) {
    ctl   <- admControl(studies = st0, grad = g, n_sim = 3000L, print = 0L,
                        covMethod = "none")
    pinfo <- admixr2:::.admDriverPinfo(ui, ctl)
    ov    <- admixr2:::.admBuildOptVec(pinfo)
    u     <- admixr2:::.admDriverUnits(st0, ui, ovar)
    stu   <- admixr2:::.admCheckCovariates(ui, pinfo, u$studies, g)
    expect_identical(stu[[1L]]$.adm_cov_path, "rows")   # gradients allowed here

    zl <- admixr2:::.admMakeZ(3000L, pinfo, 1L, "sobol")
    pl <- admixr2:::.admMakeParamsList(3000L, pinfo, 1L)
    rx <- admixr2:::.admLoadModel(ui)
    sm <- if (g == "sens") tryCatch(admixr2:::.admLoadSensModel(ui),
                                    error = function(e) NULL) else NULL
    f  <- function(pp) admixr2:::.admNLL(pp, pinfo, stu, zl, rx, ovar, pl, 1L)
    ga <- admixr2:::.admGrad(ov$p0, pinfo, stu, zl, rx, ovar, pl, 1L, 1e-4,
                             sensModel = sm)
    h  <- 1e-5
    gf <- vapply(seq_along(ov$p0), function(k) {
      a <- b <- ov$p0; a[k] <- a[k] + h; b[k] <- b[k] - h
      (f(a) - f(b)) / (2 * h)
    }, numeric(1))
    expect_true(all(is.finite(ga)))
    # both covariate coefficients must be right, not just the ones with an eta
    expect_lt(max(abs(ga - gf) / pmax(abs(gf), 1e-8)), 2e-2)
  }
})

test_that("`grad` does not change the covariate path, end to end", {
  E0 <- rep(1, length(.cov_TIMES)); V0 <- diag(length(.cov_TIMES))
  d_none <- .cov_setup(.cov_linear, list(WT = list(mu = 0.2, sd = 0.35)),
                       list(WT = 0.2), E0, V0)                   # grad = "none"
  expect_identical(d_none$st[[1L]]$.adm_cov_path, "rows")
  ui <- d_none$ui; pinfo <- d_none$pinfo
  st <- admixr2:::.admFlattenStudies(list(s1 = admixr2:::.admNormaliseStudy(
    list(E = E0, V = V0, n = 300L, times = .cov_TIMES,
         ev = rxode2::et(amt = .cov_DOSE), cov = list(WT = 0.2),
         cov_dist = list(WT = list(mu = 0.2, sd = 0.35))), "s1", "cp")))
  for (g in c("sens", "fd"))
    expect_identical(
      admixr2:::.admCheckCovariates(ui, pinfo, st, g)[[1L]]$.adm_cov_path, "rows")
})

# ---- adgh: deterministic product grid ---------------------------------------

test_that("adgh marginalises a covariate by a PRODUCT GRID, not Monte Carlo", {
  # adgh's analogue of admc's per-row draws is a product grid over the covariate
  # quadrature and the eta grid -- still ONE rxSolve, but deterministic, so adgh
  # keeps its noise-free objective. Accuracy is correspondingly ~1e-6 rather
  # than admc's ~1e-4/1e-3 at a comparable cost.
  ml <- log(72); sl <- 0.28
  r  <- .cov_ref2(function(wt, e) exp(.cov_TCL + e) * (wt / 70)^0.75,
                  function(wt) exp(.cov_TV) * (wt / 70)^1.0, ml, sl)
  Vo <- r$V; diag(Vo) <- diag(Vo) + .cov_ADD^2
  st0 <- list(s = list(E = r$E, V = Vo, n = 300L, times = .cov_TIMES,
                       ev = rxode2::et(amt = .cov_DOSE),
                       cov_dist = list(WT = list(meanlog = ml, sdlog = sl))))
  ui   <- suppressMessages(rxode2::rxode2(.cov_both))
  ovar <- admixr2:::.admOutputVar(ui)
  for (g in c("none", "analytical")) {
    ctl   <- adghControl(studies = st0, grad = g, n_nodes = 9L, print = 0L,
                         covMethod = "none")
    pinfo <- admixr2:::.admDriverPinfo(ui, ctl)
    u     <- admixr2:::.admDriverUnits(st0, ui, ovar)
    stu   <- admixr2:::.admCheckCovariates(ui, pinfo, u$studies, g, "adgh")
    expect_identical(stu[[1L]]$.adm_cov_path, "rows")
    prs  <- admixr2:::.admUnpack(admixr2:::.admBuildOptVec(pinfo)$p0, pinfo)
    grid <- admixr2:::.adghNodeGrid(9L, pinfo$n_eta)
    mm   <- admixr2:::.adghMoments(prs, pinfo, stu[[1L]],
                                   admixr2:::.admLoadModel(ui), ovar, grid, 1L)
    expect_lt(max(abs(mm$E / r$E - 1)), 1e-4)
    # diag(0.09) would be a 0x0 matrix -- the diag(scalar) trap
    expect_lt(max(abs((mm$V - diag(.cov_ADD^2, length(.cov_TIMES))) / r$V - 1)), 1e-4)
  }
})

test_that("estimators without a covariate path REFUSE cov_dist", {
  # The dangerous outcome is silence: every study also carries a covariate VALUE,
  # so an unwired estimator does not fail -- it solves at the covariate mean and
  # reports a fit whose omega has absorbed the covariate spread.
  st <- list(a = list(cov_dist = list(WT = list(mu = 0, sd = 1))),
             b = list())
  for (est in c("adfo", "adirmc"))
    expect_error(admixr2:::.admRefuseCovariates(st, est),
                 "does not support covariate marginalisation")
  expect_silent(admixr2:::.admRefuseCovariates(list(a = list(), b = list()), "adfo"))
})

test_that("adgh's ANALYTIC gradient carries the covariate product grid", {
  # .adghGradNLL used to build its quadrature from pars$L directly, so it never
  # saw the covariate grid .adghNLL evaluates on -- it differentiated a different
  # function than the objective, during the FIT. It now derives the grid per
  # study through .adghGrid(), the same helper the objective uses.
  #
  # Evaluated AWAY from the optimum on purpose: the reference data are generated
  # at the true parameters, so at p0 every component is ~0 and the comparison
  # could not tell a correct gradient from a broken one.
  ml <- log(72); sl <- 0.28
  r  <- .cov_ref2(function(wt, e) exp(.cov_TCL + e) * (wt / 70)^0.75,
                  function(wt) exp(.cov_TV) * (wt / 70)^1.0, ml, sl)
  Vo <- r$V; diag(Vo) <- diag(Vo) + .cov_ADD^2
  st0 <- list(s = list(E = r$E, V = Vo, n = 300L, times = .cov_TIMES,
                       ev = rxode2::et(amt = .cov_DOSE),
                       cov_dist = list(WT = list(meanlog = ml, sdlog = sl))))
  ui    <- suppressMessages(rxode2::rxode2(.cov_both))
  ovar  <- admixr2:::.admOutputVar(ui)
  ctl   <- adghControl(studies = st0, grad = "analytical", n_nodes = 7L,
                       print = 0L, covMethod = "none")
  pinfo <- admixr2:::.admDriverPinfo(ui, ctl)
  u     <- admixr2:::.admDriverUnits(st0, ui, ovar)
  stu   <- admixr2:::.admCheckCovariates(ui, pinfo, u$studies, "analytical", "adgh")
  expect_identical(stu[[1L]]$.adm_cov_path, "rows")

  sm   <- tryCatch(admixr2:::.admLoadSensModel(ui), error = function(e) NULL)
  rx   <- admixr2:::.admLoadModel(ui)
  ov   <- admixr2:::.admBuildOptVec(pinfo)
  grid <- admixr2:::.adghNodeGrid(7L, pinfo$n_eta)
  f    <- function(pp) admixr2:::.adghNLL(pp, pinfo, stu, rx, ovar, grid, 1L)
  p1   <- ov$p0 + c(0.15, -0.10, 0.20, -0.18, 0.25, 0.30)[seq_along(ov$p0)]
  ga   <- admixr2:::.adghGrad(p1, pinfo, stu, sm, rx, ovar, grid, 1L, 1e-4)
  h    <- 1e-5
  gf   <- vapply(seq_along(p1), function(k) {
    a <- b <- p1; a[k] <- a[k] + h; b[k] <- b[k] - h; (f(a) - f(b)) / (2 * h)
  }, numeric(1))
  expect_gt(max(abs(gf)), 100)                       # the test has real signal
  expect_lt(max(abs(ga - gf) / pmax(abs(gf), 1e-6)), 1e-3)
})

# -- A covariate that scales the DOSE, not just the parameters -----------------
#
# `f(centr) <- WT` (a mg/kg dose) is the PAGE case study's shape: the same WT
# draw has to reach the dosing modifier AND cl/vp/q in the same solve. Nothing
# else covers a covariate entering a dosing modifier, and the failure mode is
# quiet -- a WT that reaches the parameters but not f() gives a fit that
# converges to plausible numbers for a dose it was not given.
#
# Reference is INDEPENDENT: subjects drawn with their own WT and etas, solved by
# plain rxode2, reduced to (E, V). No admixr2 covariate machinery in it.

test_that("a covariate scaling the dose is marginalised correctly", {
  WT <- list(meanlog = log(18.4), sdlog = 0.45)
  TT <- c(0.5, 1, 2, 4, 6)
  EV <- rxode2::et(amt = 17.5, dur = 1)
  PR <- 0.12

  mk <- function(dose_scaled) {
    if (dose_scaled)
      function() {
        ini({lcl <- log(0.16); lvc <- log(3.86); clwt <- 0.97
             eta.cl ~ 0.062; prop.err <- 0.12})
        model({cl <- exp(lcl + eta.cl) * WT^clwt; vc <- exp(lvc)
               f(centr) <- WT
               d/dt(centr) <- -cl / vc * centr
               cp <- centr / vc; cp ~ prop(prop.err)})
      }
    else
      function() {
        ini({lcl <- log(0.16); lvc <- log(3.86); clwt <- 0.97
             eta.cl ~ 0.062; prop.err <- 0.12})
        model({cl <- exp(lcl + eta.cl) * WT^clwt; vc <- exp(lvc)
               d/dt(centr) <- -cl / vc * centr
               cp <- centr / vc; cp ~ prop(prop.err)})
      }
  }

  ref <- function(fn, n = 40000L) {
    ui <- suppressMessages(rxode2::rxode2(fn))
    rx <- suppressMessages(rxode2::rxode2(ui$simulationModel))
    set.seed(11L)
    p  <- data.frame(lcl = log(0.16), lvc = log(3.86), clwt = 0.97,
                     prop.err = PR, rxerr.cp = 0,
                     eta.cl = sqrt(0.062) * stats::rnorm(n),
                     WT = exp(WT$meanlog + WT$sdlog * stats::rnorm(n)))
    out <- rxode2::rxSolve(rx, params = p, events = EV |> rxode2::et(TT),
                           cores = 1L, nDisplayProgress = .Machine$integer.max)
    cpm <- matrix(out[["cp"]][out[["time"]] %in% TT], nrow = n, byrow = TRUE)
    mu  <- colMeans(cpm); V <- crossprod(sweep(cpm, 2L, mu)) / n
    diag(V) <- diag(V) + PR^2 * (mu^2 + diag(V))   # law of total variance
    list(E = mu, V = V)
  }

  adm <- function(fn) {
    ui    <- suppressMessages(rxode2::rxode2(fn))
    ovar  <- admixr2:::.admOutputVar(ui)
    ctl   <- adghControl(studies = list(A = list(
               E = rep(1, length(TT)), V = diag(length(TT)) + 0.01, n = 20L,
               times = TT, ev = EV, cov_dist = list(WT = WT))),
               grad = "analytical", n_nodes = 9L, print = 0L,
               covMethod = "none")
    pinfo <- admixr2:::.admDriverPinfo(ui, ctl)
    u  <- admixr2:::.admCheckCovariates(
            ui, pinfo, admixr2:::.admDriverUnits(ctl$studies, ui, ovar)$studies,
            "analytical", "adgh")
    pars <- admixr2:::.admUnpack(admixr2:::.admBuildOptVec(pinfo)$p0, pinfo)
    list(m = admixr2:::.adghMoments(pars, pinfo, u[[1L]], ui |>
                                      admixr2:::.admLoadModel(), ovar,
                                    admixr2:::.adghNodeGrid(9L, pinfo$n_eta), 1L),
         path = u[[1L]]$.adm_cov_path)
  }

  a <- adm(mk(TRUE)); r <- ref(mk(TRUE))
  # per-subject covariates: the same draw feeds f() and cl in one solve
  expect_equal(a$path, "rows")
  expect_lt(max(abs(a$m$E - r$E) / abs(r$E)), 0.01)
  expect_lt(max(abs(a$m$V - r$V)) / max(abs(r$V)), 0.02)

  # the check must be SENSITIVE to the dose term, or it proves nothing: dropping
  # f(centr) <- WT moves the mean by more than an order of magnitude
  r0 <- ref(mk(FALSE))
  expect_gt(max(abs(r$E - r0$E) / r0$E), 5)
  a0 <- adm(mk(FALSE))
  expect_lt(max(abs(a0$m$E - r0$E) / abs(r0$E)), 0.01)
})

# -- The node-study guard must fire THROUGH A DRIVER --------------------------
#
# test-covariate.R exercises .admRefuseNodeStudies() on a raw study list, which
# is not the shape any driver hands it. That gap hid a real defect: the drivers
# called the guard AFTER .admDriverUnits(), which strips the top-level `weight`
# and `cov_method` fields it reads, so the guard never fired and a node study
# list FITTED in all four estimators -- at exactly twice the correct objective
# (720.715 against 360.358), weights silently ignored. Only an end-to-end test
# can catch that, so this one goes through nlmixr2().

test_that("a node study list is refused by every estimator, end to end", {
  skip_if_not_installed("nlmixr2est")
  m <- function() {
    ini({tcl <- log(1); tv <- log(10); add.err <- 0.3; eta.cl ~ 0.09})
    model({cl <- exp(tcl + eta.cl); v <- exp(tv)
           d/dt(centr) <- -cl / v * centr; cp <- centr / v; cp ~ add(add.err)})
  }
  TT <- c(0.5, 1, 2, 4, 8); EV <- rxode2::et(amt = 100)
  g  <- datagen(list(s = list(model = m, times = TT, ev = EV, n = 100L)),
                control = datagenControl(n_sim = 2000L, seed = 3L))
  nodes <- list(n1 = c(g$s, list(weight = 0.5)), n2 = c(g$s, list(weight = 0.5)))

  for (est in c("admc", "adgh", "adfo", "adirmc")) {
    ctl <- switch(est, admc = admControl, adgh = adghControl,
                  adfo = adfoControl, adirmc = adirmcControl)
    expect_error(
      suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
        m, admData(), est = est,
        control = ctl(studies = nodes, print = 0L, covMethod = "none",
                      maxeval = 3L)))),
      "was removed", info = est)
  }

  # a multi-output study can carry the field per OBSERVATION rather than at the
  # top level; the guard has to look there too
  nodes2 <- list(s = list(n = 100L, ev = EV, observations = list(
    o1 = list(E = g$s$E, V = g$s$V, times = TT, output = "cp", weight = 0.5))))
  expect_error(
    suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
      m, admData(), est = "adgh",
      control = adghControl(studies = nodes2, print = 0L, covMethod = "none",
                            maxeval = 3L)))),
    "was removed")

  # ... and an ordinary study list must still fit
  expect_no_error(
    suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
      m, admData(), est = "adgh",
      control = adghControl(studies = g, print = 0L, covMethod = "none",
                            maxeval = 3L)))))
})

test_that("a DEPENDENT covariate distribution supports analytic gradients", {
  # The `joint` sampler is the newest covariate path and the one a copula or an
  # R-vine arrives through. It reaches the solve as per-row params columns, the
  # same as the independent path, so the sensitivity directions still
  # differentiate the function the NLL evaluates -- but the FD frames in
  # .admGrad/.adghGrad tile those columns by hand, and a sampler consumes TWO
  # uniform columns rather than one. Nothing else pins that, and a joint sampler
  # dropped from a tiled frame would be a silently wrong gradient under a
  # finite, plausible objective.
  skip_if_not_installed("rxode2")
  ml <- log(72); sl <- 0.28; mc <- log(90); sc <- 0.30; rho <- 0.6
  # A Gaussian copula written the way a user would: consume the supplied
  # uniforms, return one column per declared covariate.
  # NOTE the clamp AFTER pnorm as well as before qnorm. The quadrature grid's
  # tail nodes reach |z| ~ 6.4, a copula's mixing step scales that by up to
  # ~1.4, and pnorm() of the result rounds to exactly 1 -- after which
  # qlnorm(1) is Inf. A per-subject sample never gets that far out, so this
  # only bites the deterministic paths.
  cl <- function(x) pmin(pmax(x, 1e-12), 1 - 1e-12)
  jf <- function(u) {
    z  <- stats::qnorm(cl(u))
    z2 <- rho * z[, 1] + sqrt(1 - rho^2) * z[, 2]
    cbind(WT   = stats::qlnorm(cl(stats::pnorm(z[, 1])), ml, sl),
          CRCL = stats::qlnorm(cl(stats::pnorm(z2)),     mc, sc))
  }
  cd <- list(WT   = list(quantile = function(u) stats::qlnorm(u, ml, sl)),
             CRCL = list(quantile = function(u) stats::qlnorm(u, mc, sc)),
             joint = jf)
  mod <- function() {
    ini({tcl <- log(3.2); tv <- log(21); bwt <- 0.75; bcr <- 0.45
         eta.cl ~ 0.09; add.err <- 0.6})
    model({cl <- exp(tcl + eta.cl) * (WT / 70)^bwt * (CRCL / 95)^bcr
           v  <- exp(tv)
           d/dt(centr) <- -cl / v * centr
           cp <- centr / v
           cp ~ add(add.err)})
  }
  ui   <- suppressMessages(rxode2::rxode2(mod))
  ovar <- admixr2:::.admOutputVar(ui)
  tms  <- c(1, 2, 4, 8)
  st0  <- list(s1 = list(E = rep(1.5, length(tms)), V = diag(length(tms)),
                         n = 200L, times = tms,
                         ev = rxode2::et(amt = 100), cov_dist = cd))

  ## ---- admc: per-subject draws --------------------------------------------
  ctl   <- admControl(studies = st0, grad = "sens", n_sim = 4000L, print = 0L,
                      covMethod = "none")
  pinfo <- admixr2:::.admDriverPinfo(ui, ctl)
  ov    <- admixr2:::.admBuildOptVec(pinfo)
  u     <- admixr2:::.admDriverUnits(st0, ui, ovar)
  stu   <- admixr2:::.admCheckCovariates(ui, pinfo, u$studies, "sens")
  expect_identical(stu[[1L]]$.adm_cov_path, "rows")
  zl <- admixr2:::.admMakeZ(4000L, pinfo, 1L, "sobol")
  pl <- admixr2:::.admMakeParamsList(4000L, pinfo, 1L)
  rx <- admixr2:::.admLoadModel(ui)
  sm <- tryCatch(admixr2:::.admLoadSensModel(ui), error = function(e) NULL)
  f  <- function(pp) admixr2:::.admNLL(pp, pinfo, stu, zl, rx, ovar, pl, 1L)
  ga <- admixr2:::.admGrad(ov$p0, pinfo, stu, zl, rx, ovar, pl, 1L, 1e-4,
                           sensModel = sm)
  h  <- 1e-5
  gf <- vapply(seq_along(ov$p0), function(k) {
    a <- b <- ov$p0; a[k] <- a[k] + h; b[k] <- b[k] - h
    (f(a) - f(b)) / (2 * h)
  }, numeric(1))
  expect_true(all(is.finite(ga)))
  expect_lt(max(abs(ga - gf) / pmax(abs(gf), 1e-8)), 2e-2)

  ## ---- adgh now INTEGRATES a dependent joint, on the u-space grid ---------
  # It used to refuse one, on the grounds that a product grid assumes
  # independence. That is true of a grid over covariate MARGINS and false of a
  # grid over the sampler's UNIFORMS: a copula maps independent uniforms to
  # dependent values, so the product rule is exact there whatever the
  # dependence. The check that matters is that adgh's moments match the
  # per-subject draws admc uses -- the two estimators must see one distribution.
  fit <- suppressMessages(nlmixr2est::nlmixr2(
    mod, admData(), est = "adgh",
    control = adghControl(studies = st0, print = 0L, covMethod = "none",
                          maxeval = 2L)))
  expect_s3_class(fit, "admFit")
  expect_true(is.finite(fit$objective))

  ctlg  <- adghControl(studies = st0, grad = "analytical", print = 0L,
                       covMethod = "none", n_nodes = 7L, cov_nodes = 15L)
  ping  <- admixr2:::.admDriverPinfo(ui, ctlg)
  stug  <- admixr2:::.admCheckCovariates(
             ui, ping, admixr2:::.admDriverUnits(st0, ui, ovar)$studies,
             "analytical", "adgh")
  gridg <- admixr2:::.adghNodeGrid(7L, ping$n_eta)
  parsg <- admixr2:::.admUnpack(admixr2:::.admBuildOptVec(ping)$p0, ping)
  mg    <- admixr2:::.adghMoments(parsg, ping, stug[[1L]], rx, ovar, gridg, 1L)

  # the same moments from the per-subject path admc drives
  crow <- admixr2:::.admStudyCovRows(stu[[1L]], pinfo, 60000L)
  zz   <- admixr2:::.admMakeZ(60000L, pinfo, 1L, "sobol")[[1L]]
  etam <- zz %*% t(parsg$L); colnames(etam) <- pinfo$eta_col_names
  pmm  <- admixr2:::.admMakeParamsList(60000L, pinfo, 1L)[[1L]]
  cpm  <- admixr2:::.admSimulate(rx, parsg$struct, pinfo$sigma_names, etam,
                                 crow, ovar, pmm, 1L,
                                 .Machine$integer.max, NULL)
  mum  <- colMeans(cpm); Vm <- crossprod(sweep(cpm, 2L, mum)) / nrow(cpm)
  arrm <- admixr2:::.admUnitResidRows(ping, ovar, parsg$sigma_var, length(mum))
  apm  <- admixr2:::.admResidApply(mum, diag(Vm), arrm, crow$times, Vm)
  Em   <- apm$mu; Vmc <- admixr2:::.admApplyResidTail(Vm, apm)

  expect_lt(max(abs(mg$E - Em) / abs(Em)), 5e-3)
  expect_lt(max(abs(mg$V - Vmc)) / max(abs(Vmc)), 5e-3)

  # ... and the dependence must actually MOVE the moments, or agreeing here
  # would prove nothing (an arm that silently dropped `cor` would also pass).
  st_ind <- st0                                   # NOTE: the study is `s1`
  st_ind$s1$cov_dist <- st0$s1$cov_dist[
    setdiff(names(st0$s1$cov_dist), c("cor", "rho", "Sigma", "joint"))]
  expect_null(st_ind$s1$cov_dist$joint)           # or the contrast is vacuous
  stui <- admixr2:::.admCheckCovariates(
            ui, ping, admixr2:::.admDriverUnits(st_ind, ui, ovar)$studies,
            "analytical", "adgh")
  mi <- admixr2:::.adghMoments(parsg, ping, stui[[1L]], rx, ovar, gridg, 1L)
  expect_gt(max(abs(mg$V - mi$V)) / max(abs(mi$V)), 1e-3)
})

# -- Second-order Taylor covariate integration (cov_integration = "taylor") ----
#
# The Tier-1 file checks the expansion itself against exact nested quadrature on
# an analytic solution. What is only testable here is the PIPELINE: that the
# design points reach rxSolve as per-row covariates, that the residual and the
# NLL are formed from the expanded moments, that the analytic gradient
# differentiates the function the objective evaluates, and that asking for
# "quadrature" leaves every number exactly where it was.

.tay_setup <- function(ci, grad = "analytical", n_nodes = 7L, ml, sl,
                       hfrac = 0.5, E, V) {
  st0 <- list(s = list(E = E, V = V, n = 300L, times = .cov_TIMES,
                       ev = rxode2::et(amt = .cov_DOSE),
                       cov_dist = list(WT = list(meanlog = ml, sdlog = sl))))
  ui    <- suppressMessages(rxode2::rxode2(.cov_both))
  ovar  <- admixr2:::.admOutputVar(ui)
  ctl   <- adghControl(studies = st0, grad = grad, n_nodes = n_nodes, print = 0L,
                       covMethod = "none", cov_integration = ci,
                       cov_taylor_h = hfrac)
  pinfo <- admixr2:::.admDriverPinfo(ui, ctl)
  u     <- admixr2:::.admDriverUnits(st0, ui, ovar)
  stu   <- admixr2:::.admCheckCovariates(ui, pinfo, u$studies, grad, "adgh")
  list(ui = ui, ovar = ovar, pinfo = pinfo, stu = stu,
       ov = admixr2:::.admBuildOptVec(pinfo),
       grid = admixr2:::.adghNodeGrid(n_nodes, pinfo$n_eta),
       rxMod = admixr2:::.admLoadModel(ui))
}

test_that("cov_integration = 'taylor' expands the marginal moments, in 1 + 2p points", {
  ml <- log(72); sl <- 0.28
  cl_of <- function(wt, e) exp(.cov_TCL + e) * (wt / 70)^0.75
  v_of  <- function(wt)    exp(.cov_TV) * (wt / 70)^1.0
  r  <- .cov_ref2(cl_of, v_of, ml, sl)
  # the ecological plug-in: the same model solved AT the covariate mean, which
  # is what a fit with a point `cov` and no marginalisation would report
  pg <- .cov_ref2(cl_of, v_of, log(exp(ml + sl^2 / 2)), 1e-8)
  Vo <- r$V; diag(Vo) <- diag(Vo) + .cov_ADD^2

  d  <- .tay_setup("taylor", ml = ml, sl = sl, E = r$E, V = Vo)
  expect_identical(d$stu[[1L]]$.adm_cov_path, "rows")
  prs <- admixr2:::.admUnpack(d$ov$p0, d$pinfo)
  mm  <- admixr2:::.adghMoments(prs, d$pinfo, d$stu[[1L]], d$rxMod, d$ovar,
                                d$grid, 1L)
  Vs  <- mm$V - diag(.cov_ADD^2, length(.cov_TIMES))
  eE  <- max(abs(mm$E / r$E - 1)); eV <- max(abs(Vs / r$V - 1))
  # Measured 1.9e-03 / 1.0e-01. This is a demanding regime for the expansion:
  # WT enters v with NO random effect at all, so on that channel the covariate
  # is the only source of variability and the effective ratio is unbounded.
  expect_lt(eE, 5e-3)
  expect_lt(eV, 2e-1)
  # It must still be a large improvement on the plug-in, which is the thing it
  # is an alternative to -- getting that backwards is the whole risk.
  expect_lt(eV, max(abs(pg$V / r$V - 1)) / 5)
  expect_lt(eE, max(abs(pg$E / r$E - 1)) / 5)

  # 1 + 2p design points, not cov_nodes^p: the params frame the solve sees has
  # 3 * n_nodes rows where quadrature would have had 11 * n_nodes.
  g <- admixr2:::.adghGrid(prs, d$pinfo, d$grid, d$stu[[1L]])
  expect_identical(nrow(g$eta), 3L * nrow(d$grid$X))
  expect_identical(nrow(unique(g$cov_rows)), 3L)
  dq <- .tay_setup("quadrature", ml = ml, sl = sl, E = r$E, V = Vo)
  gq <- admixr2:::.adghGrid(prs, dq$pinfo, dq$grid, dq$stu[[1L]])
  # read the default rather than hard-coding it: cov_nodes is a tuning
  # default and pinning its VALUE here made this assertion stale the moment
  # it moved.  What matters is that quadrature costs cov_nodes^p and the
  # expansion costs 1 + 2p.
  expect_identical(nrow(gq$eta),
                   as.integer(dq$pinfo$cov_nodes) * nrow(dq$grid$X))
  expect_gt(nrow(gq$eta), nrow(g$eta))
})

test_that("the Taylor path carries an ANALYTIC gradient (vs central FD)", {
  # The rank-p term sum_j v_j g'_j g'_j' is quadratic in the conditional means,
  # so it contributes a derivative the weighted-crossproduct contraction does
  # not. Evaluated AWAY from the optimum, where every component is large.
  ml <- log(72); sl <- 0.28
  r  <- .cov_ref2(function(wt, e) exp(.cov_TCL + e) * (wt / 70)^0.75,
                  function(wt) exp(.cov_TV) * (wt / 70)^1.0, ml, sl)
  Vo <- r$V; diag(Vo) <- diag(Vo) + .cov_ADD^2
  d  <- .tay_setup("taylor", ml = ml, sl = sl, E = r$E, V = Vo)
  sm <- tryCatch(admixr2:::.admLoadSensModel(d$ui), error = function(e) NULL)
  f  <- function(pp)
    admixr2:::.adghNLL(pp, d$pinfo, d$stu, d$rxMod, d$ovar, d$grid, 1L)
  p1 <- d$ov$p0 + c(0.15, -0.10, 0.20, -0.18, 0.25, 0.30)[seq_along(d$ov$p0)]
  ga <- admixr2:::.adghGrad(p1, d$pinfo, d$stu, sm, d$rxMod, d$ovar, d$grid,
                            1L, 1e-4)
  h  <- 1e-5
  gf <- vapply(seq_along(p1), function(k) {
    a <- b <- p1; a[k] <- a[k] + h; b[k] <- b[k] - h; (f(a) - f(b)) / (2 * h)
  }, numeric(1))
  expect_gt(max(abs(gf)), 100)                       # the test has real signal
  expect_lt(max(abs(ga - gf) / pmax(abs(gf), 1e-6)), 1e-3)
})

test_that("cov_integration = 'quadrature' is the default and changes nothing", {
  ml <- log(72); sl <- 0.28
  r  <- .cov_ref2(function(wt, e) exp(.cov_TCL + e) * (wt / 70)^0.75,
                  function(wt) exp(.cov_TV) * (wt / 70)^1.0, ml, sl)
  Vo <- r$V; diag(Vo) <- diag(Vo) + .cov_ADD^2
  st0 <- list(s = list(E = r$E, V = Vo, n = 300L, times = .cov_TIMES,
                       ev = rxode2::et(amt = .cov_DOSE),
                       cov_dist = list(WT = list(meanlog = ml, sdlog = sl))))
  ui   <- suppressMessages(rxode2::rxode2(.cov_both))
  ovar <- admixr2:::.admOutputVar(ui)
  rxM  <- admixr2:::.admLoadModel(ui)
  sm   <- tryCatch(admixr2:::.admLoadSensModel(ui), error = function(e) NULL)
  one <- function(ctl) {
    pinfo <- admixr2:::.admDriverPinfo(ui, ctl)
    u     <- admixr2:::.admDriverUnits(st0, ui, ovar)
    stu   <- admixr2:::.admCheckCovariates(ui, pinfo, u$studies, "analytical", "adgh")
    ov    <- admixr2:::.admBuildOptVec(pinfo)
    grid  <- admixr2:::.adghNodeGrid(7L, pinfo$n_eta)
    p1    <- ov$p0 + c(0.15, -0.10, 0.20, -0.18, 0.25, 0.30)[seq_along(ov$p0)]
    list(nll = admixr2:::.adghNLL(p1, pinfo, stu, rxM, ovar, grid, 1L),
         grad = admixr2:::.adghGrad(p1, pinfo, stu, sm, rxM, ovar, grid, 1L, 1e-4))
  }
  base <- one(adghControl(studies = st0, grad = "analytical", n_nodes = 7L,
                          print = 0L, covMethod = "none"))
  same <- one(adghControl(studies = st0, grad = "analytical", n_nodes = 7L,
                          print = 0L, covMethod = "none",
                          cov_integration = "quadrature", cov_taylor_h = 0.9))
  # BIT-identical, not merely close: the quadrature path must not move by an ulp
  expect_identical(same$nll,  base$nll)
  expect_identical(same$grad, base$grad)
})

test_that("the Taylor path ENUMERATES a discrete covariate, per study", {
  # A discrete covariate is not expanded, it is enumerated: its levels and
  # probabilities ARE the integration rule, exactly, so a study whose covariate
  # is two-point costs 2 design points and reproduces the enumeration.  The
  # continuous study alongside it is unaffected and still costs 1 + 2p.
  .st <- function(cd) list(E = rep(1, length(.cov_TIMES)),
                           V = diag(length(.cov_TIMES)), n = 100L,
                           times = .cov_TIMES,
                           ev = rxode2::et(amt = .cov_DOSE), cov_dist = cd)
  st0 <- list(s1 = .st(list(WT = list(meanlog = log(72), sdlog = 0.28))),
              s2 = .st(list(WT = list(values = c(60, 85), probs = c(0.4, 0.6)))))
  ui   <- suppressMessages(rxode2::rxode2(.cov_both))
  ovar <- admixr2:::.admOutputVar(ui)
  u    <- admixr2:::.admDriverUnits(st0, ui, ovar)

  ctl <- adghControl(studies = st0, grad = "analytical", print = 0L,
                     covMethod = "none", cov_integration = "taylor")
  stu <- admixr2:::.admCheckCovariates(ui, admixr2:::.admDriverPinfo(ui, ctl),
                                       u$studies, "analytical", "adgh")
  td2 <- stu$s2$.adm_cov_taylor
  expect_identical(td2$n_pt, 2L)
  expect_identical(td2$n_cell, 2L)
  expect_identical(td2$n_cpt, 1L)              # no cubature points at all
  expect_equal(td2$c, c(0.4, 0.6))             # the level probabilities exactly
  expect_equal(as.numeric(td2$X[, "WT"]), c(60, 85))
  expect_identical(stu$s1$.adm_cov_taylor$n_pt, 3L)

  # and the very same studies go through on the quadrature route
  ctlq <- adghControl(studies = st0, grad = "analytical", print = 0L,
                      covMethod = "none")
  expect_no_error(
    admixr2:::.admCheckCovariates(ui, admixr2:::.admDriverPinfo(ui, ctlq),
                                  u$studies, "analytical", "adgh"))

  # What IS still refused: a discrete covariate DEPENDENT on a continuous one.
  # A level is then a truncation of the latent normal rather than a point, so
  # the continuous conditional differs cell by cell and expanding it about the
  # marginal mean would be the wrong expansion in every cell.
  Rz <- matrix(c(1, 0.3, 0.3, 1), 2L, 2L,
               dimnames = list(c("WT", "SEX"), c("WT", "SEX")))
  expect_error(
    admixr2:::.admCovTaylorDesign(
      list(WT = list(meanlog = log(72), sdlog = 0.28),
           SEX = list(values = c(0, 1), probs = c(0.5, 0.5)), latentR = Rz)),
    "DEPENDENT")
})

# =============================================================================
# Shift path: the covariate leaves the solver
# =============================================================================

.shift_mod <- function() {
  ini({tcl <- log(4); tv <- log(30); b1 <- 0.75
       eta.cl ~ 0.09; eta.v ~ 0.04; a <- 0.1})
  model({cl <- exp(tcl + eta.cl) * (WT/70)^b1; v <- exp(tv + eta.v)
         d/dt(centr) <- -cl/v*centr; cp <- centr/v; cp ~ add(a)})
}
.shift_setup <- function(integ, grad = "none", cd = NULL, mod = .shift_mod) {
  ui <- suppressMessages(rxode2::rxode2(mod))
  ov <- admixr2:::.admOutputVar(ui)
  rx <- admixr2:::.admLoadModel(ui)
  sM <- if (identical(grad, "analytical"))
    tryCatch(admixr2:::.admLoadSensModel(ui), error = function(e) NULL) else NULL
  if (is.null(cd)) cd <- covDist(WT = c(mean = 78, sd = 18), dist = "lnorm")
  st <- list(s = list(E = c(9.9, 9.2, 8.1, 6.4, 5.2, 4.3, 2.9, 2.0, 0.9),
                      V = diag(0.4, 9L) + 0.05, n = 200L,
                      times = c(0.25,0.5,1,2,3,4,6,8,12),
                      ev = rxode2::et(amt = 500), cov_dist = cd))
  ctl <- adghControl(studies = st, grad = grad, print = 0L, covMethod = "none",
                     n_nodes = 7L, cov_nodes = 7L, cov_integration = integ)
  pin <- admixr2:::.admDriverPinfo(ui, ctl)
  stu <- admixr2:::.admCheckCovariates(
    ui, pin, admixr2:::.admDriverUnits(st, ui, ov)$studies, grad, "adgh")
  list(pin = pin, stu = stu, rx = rx, sM = sM, ov = ov,
       g = admixr2:::.adghNodeGrid(7L, pin$n_eta),
       p0 = admixr2:::.admBuildOptVec(pin)$p0)
}

test_that("the shift path reproduces the product grid at a fraction of the rows", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  Sq <- .shift_setup("quadrature"); Ss <- .shift_setup("shift")
  expect_identical(Ss$stu[[1L]]$.adm_cov_path, "shift")
  pq <- admixr2:::.admUnpack(Sq$p0, Sq$pin)
  ps <- admixr2:::.admUnpack(Ss$p0, Ss$pin)
  mq <- admixr2:::.adghMoments(pq, Sq$pin, Sq$stu[[1L]], Sq$rx, Sq$ov, Sq$g, 1L)
  ms <- admixr2:::.adghMoments(ps, Ss$pin, Ss$stu[[1L]], Ss$rx, Ss$ov, Ss$g, 1L)
  expect_equal(ms$E, mq$E, tolerance = 1e-5)
  expect_equal(ms$V, mq$V, tolerance = 1e-4)
  # ... and it costs far fewer solve rows, a count that does NOT grow with the
  # number of covariates (the product grid's does, as n_cov^p)
  nq <- nrow(admixr2:::.adghGrid(pq, Sq$pin, Sq$g, Sq$stu[[1L]])$eta)
  ns <- nrow(admixr2:::.adghGrid(ps, Ss$pin, Ss$g, Ss$stu[[1L]])$eta)
  expect_lt(ns, nq / 3)
})

test_that("the shift path carries an ANALYTIC gradient (vs central FD)", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  S <- .shift_setup("shift", grad = "analytical")
  skip_if(is.null(S$sM), "no sensitivity model")
  f  <- function(p) admixr2:::.adghNLL(p, S$pin, S$stu, S$rx, S$ov, S$g, 1L)
  ga <- admixr2:::.adghGrad(S$p0, S$pin, S$stu, S$sM, S$rx, S$ov, S$g, 1L, 1e-4)
  h  <- 1e-5
  gf <- vapply(seq_along(S$p0), function(j) {
    a <- b <- S$p0; a[j] <- a[j] + h; b[j] <- b[j] - h
    (f(a) - f(b)) / (2 * h) }, 0)
  # the covariate coefficient moves the u nodes through Delta, which is the
  # term the eta-column chain alone would miss entirely
  expect_lt(max(abs(ga - gf) / pmax(abs(gf), 1e-4)), 2e-3)
})

test_that("the shift path is REFUSED when its identity does not hold", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  # an ADDITIVE covariate effect: no shift of eta reproduces it
  expect_error(.shift_setup("shift", mod = function() {
    ini({tcl <- log(4); tv <- log(30); b1 <- 0.02; eta.cl ~ 0.09; a <- 0.1})
    model({cl <- exp(tcl + eta.cl) + b1*(WT - 70); v <- exp(tv)
           d/dt(centr) <- -cl/v*centr; cp <- centr/v; cp ~ add(a)}) }),
    "shift identity")
  # a covariate on a parameter with NO random effect: nothing to shift
  expect_error(.shift_setup("shift", mod = function() {
    ini({tcl <- log(4); tv <- log(30); b2 <- 1.0; eta.cl ~ 0.09; a <- 0.1})
    model({cl <- exp(tcl + eta.cl); v <- exp(tv)*(WT/70)^b2
           d/dt(centr) <- -cl/v*centr; cp <- centr/v; cp ~ add(a)}) }),
    "exactly one random effect")
})

test_that('cov_integration = "auto" takes the shift exactly where it verifies', {
  skip_on_cran(); skip_if_not_installed("rxode2")
  Sa <- .shift_setup("auto")
  Ss <- .shift_setup("shift")
  expect_identical(Sa$stu[[1L]]$.adm_cov_path, "shift")
  # it must build the SAME shift the explicit setting does, not a variant
  expect_equal(Sa$stu[[1L]]$.adm_cov_shift, Ss$stu[[1L]]$.adm_cov_shift)

  # ... and where the identity fails it FALLS BACK, never errors. "auto" is a
  # speed lever and the fallback is the more accurate path, so a refusal must
  # cost solve rows and nothing else. The reason is recorded on the study as
  # well as messaged -- a message is easy to lose in a fit's output.
  expect_message(
    Sb <- .shift_setup("auto", mod = function() {
      ini({tcl <- log(4); tv <- log(30); b1 <- 0.02; eta.cl ~ 0.09; a <- 0.1})
      model({cl <- exp(tcl + eta.cl) + b1*(WT - 70); v <- exp(tv)
             d/dt(centr) <- -cl/v*centr; cp <- centr/v; cp ~ add(a)}) }),
    "shift identity")
  expect_identical(Sb$stu[[1L]]$.adm_cov_path, "rows")
  expect_match(Sb$stu[[1L]]$.adm_cov_shift_why, "shift identity")
  # the default is unchanged: "quadrature" never tries the shift
  expect_identical(.shift_setup("quadrature")$stu[[1L]]$.adm_cov_path, "rows")
})

test_that("a shift study's FD gradient is not batched onto a stale node grid", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  # .adghMomentsBatch scores every perturbed configuration against ONE node grid
  # and ONE weight vector -- true once "collapse" is gone, EXCEPT on the shift
  # path, where Delta moves the u nodes and their weights with the structural
  # thetas. Batching there differences the objective on the unperturbed grid: a
  # finite, plausible, wrong gradient. .adghGrad routes shift studies to the
  # per-configuration path instead.
  S  <- .shift_setup("shift", grad = "analytical")
  skip_if(is.null(S$sM), "no sensitivity model")
  # Drop the theta directions to force the FD block: .admThetaSens() returns
  # NULL without them, which is exactly what a model that cannot build an
  # augmented sens model does. (Passing sensModel = NULL is NOT the same thing
  # -- .admSimulateSens dereferences it unconditionally.)
  sM <- S$sM; sM$theta_sens_cols <- NULL
  f  <- function(p) admixr2:::.adghNLL(p, S$pin, S$stu, S$rx, S$ov, S$g, 1L)
  ga <- admixr2:::.adghGrad(S$p0, S$pin, S$stu, sM, S$rx, S$ov, S$g, 1L, 1e-4)
  h  <- 1e-5
  gf <- vapply(seq_along(S$p0), function(j) {
    a <- b <- S$p0; a[j] <- a[j] + h; b[j] <- b[j] - h
    (f(a) - f(b)) / (2 * h) }, 0)
  expect_true(all(is.finite(ga)))
  expect_lt(max(abs(ga - gf) / pmax(abs(gf), 1e-4)), 1e-2)
})
