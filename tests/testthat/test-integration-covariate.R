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
  pinfo$cov_map <- admixr2:::.admCovMap(ui)
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
  if (identical(su$.adm_cov_path, "uq")) {
    eta <- admixr2:::.admCovUQEta(eta, z, d$pars, d$pinfo, su, d$rxMod, 1L)
    expect_false(is.null(eta))
    colnames(eta) <- d$pinfo$eta_col_names
  } else if (identical(su$.adm_cov_path, "rows")) {
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

test_that(".admCovDelta measures Delta(a) from the model, exactly", {
  cd <- list(WT = list(meanlog = log(72), sdlog = 0.28))
  d  <- .cov_setup(.cov_allometric, cd, list(WT = 72),
                   rep(1, length(.cov_TIMES)), diag(length(.cov_TIMES)))
  dl <- admixr2:::.admCovDelta(d$rxMod, d$pars, d$pinfo, d$st[[1L]], "WT", "cl", 1L)
  expect_false(is.null(dl))
  # the model is cl = exp(tcl+eta)*(WT/70)^tcov, so Delta(a) = tcov*log(a/72)
  # measured against the study's own reference value of 72
  expect_equal(dl$delta, 0.75 * log(dl$a / 72), tolerance = 1e-9)
  expect_equal(sum(dl$w), 1)
})

test_that("u-quantile reproduces exact nested quadrature for an ALLOMETRIC effect", {
  ml <- log(72); sl <- 0.28
  ref <- .cov_ref(function(wt, eta) exp(.cov_TCL + eta) * (wt / 70)^0.75, ml, sl)
  Vo  <- ref$V; diag(Vo) <- diag(Vo) + .cov_ADD^2
  d <- .cov_setup(.cov_allometric, list(WT = list(meanlog = ml, sdlog = sl)),
                  list(WT = exp(ml)), ref$E, Vo)
  # this model is NOT a bare theta*WT product -> u-quantile, not the collapse
  expect_false(isTRUE(d$st[[1L]]$.adm_cov_collapse))
  expect_false(is.null(d$st[[1L]]$.adm_cov_uq))
  p <- .cov_pred(d)
  expect_lt(max(abs(p$E / ref$E - 1)), 5e-3)
  expect_lt(max(abs(p$V / ref$V - 1)), 3e-2)
})

test_that("collapse and u-quantile agree on a LINEAR model, where both are valid", {
  # Two independent code paths -- closed-form Omega inflation vs a measured
  # Delta and an inverse-CDF draw -- must land on the same moments.
  mu <- 0.20; sd <- 0.35
  d_uq <- .cov_setup(.cov_linear, list(WT = list(mu = mu, sd = sd)),
                     list(WT = mu), rep(1, length(.cov_TIMES)),
                     diag(length(.cov_TIMES)))
  expect_true(isTRUE(d_uq$st[[1L]]$.adm_cov_collapse))   # collapse claimed
  p_col <- .cov_pred(d_uq)

  # force the u-quantile route on the same model/study
  d_uq$st[[1L]]$.adm_cov_collapse <- FALSE
  pe <- admixr2:::.admCovParamEta(d_uq$ui, "WT", d_uq$pinfo$eta_col_names)
  expect_equal(pe$param, "cl")
  d_uq$st[[1L]]$.adm_cov_uq <- list(list(cov = "WT", param = pe$param,
                                         eta = pe$eta, idx = 1L))
  p_uq <- .cov_pred(d_uq)

  expect_lt(max(abs(p_uq$E / p_col$E - 1)), 5e-3)
  expect_lt(max(abs(p_uq$V / p_col$V - 1)), 3e-2)
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
  d$st[[1L]]$.adm_cov_uq <- NULL          # covariate held at its mean
  d$st[[1L]]$.adm_cov_path <- NULL
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

test_that("routing prefers the most efficient VALID path", {
  ml <- log(72); sl <- 0.28
  E0 <- rep(1, length(.cov_TIMES)); V0 <- diag(length(.cov_TIMES))
  # bare theta*WT + normal -> closed form
  expect_identical(
    .cov_setup(.cov_linear, list(WT = list(mu = 0.2, sd = 0.35)),
               list(WT = 0.2), E0, V0)$st[[1L]]$.adm_cov_path, "collapse")
  # allometric, single parameter with an eta -> measured Delta
  expect_identical(
    .cov_setup(.cov_allometric, list(WT = list(meanlog = ml, sdlog = sl)),
               list(WT = exp(ml)), E0, V0)$st[[1L]]$.adm_cov_path, "uq")
  # bare theta*WT but NOT normal -> the closed form does not apply
  expect_identical(
    .cov_setup(.cov_linear, list(WT = list(values = c(0, 1), probs = c(.6, .4))),
               list(WT = 0.4), E0, V0)$st[[1L]]$.adm_cov_path, "uq")
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

test_that("a gradient mode selects the general path, and none keeps the collapse", {
  E0 <- rep(1, length(.cov_TIMES)); V0 <- diag(length(.cov_TIMES))
  d_none <- .cov_setup(.cov_linear, list(WT = list(mu = 0.2, sd = 0.35)),
                       list(WT = 0.2), E0, V0)
  expect_identical(d_none$st[[1L]]$.adm_cov_path, "collapse")
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
    stu   <- admixr2:::.admCheckCovariates(ui, pinfo, u$studies, g)
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
