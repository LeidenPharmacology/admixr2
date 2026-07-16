# Tier 2: the sens model's OWN columns, checked directly against finite
# differences of the same model's rx_pred_.
#
# WHY THIS FILE EXISTS. Every other gradient test compares .admGrad against a
# finite difference of .admNLL. That cannot see a broken sensitivity column:
#   * a missing THETA column silently falls back to FD, which is also correct;
#   * a ZERO eta column makes the analytic and FD gradients agree on being wrong.
# A parameter entering f()/lag()/rate()/dur() had exactly that -- an identically
# zero column -- and the whole suite stayed green. So: check the derivative the
# model actually returns, against the derivative of the model's own prediction.
#
# Helper + models: .int_sens_col_errs() in helper-integration.R.

skip_on_cran()
skip_if_not_installed("rxode2")
# This file compiles ~14 distinct rxode2 models (ODE, linCmt, every dosing
# modifier, transit, delay, multi-endpoint, ...) plus a sensitivity model for
# each. On the R-devel CI canary that compilation load hangs the runner, while it
# runs in <5 min on every release R. The columns are fully covered on ubuntu
# release, Windows, macOS and the compat-oldstack integration job, so skip the
# file on R-devel rather than let the non-blocking canary hang.
skip_if(grepl("Under development", R.version$status),
        "heavy model-compilation suite skipped on R-devel (covered on release)")

TOL   <- 1e-4                                   # FD reference is good to ~1e-6
times <- c(0.5, 1, 2, 4, 8, 12, 24)
bolus <- rxode2::et(amt = 100) |> rxode2::et(times)

expect_cols_exact <- function(errs, expected_cols) {
  expect_setequal(names(errs), expected_cols)
  for (nm in names(errs))
    expect_lt(errs[[nm]], TOL, label = sprintf("max rel err of column %s", nm))
  # a silently-zero column would give exactly 1.0; make that impossible to miss
  expect_true(all(errs < 0.5))
}

test_that("ODE: eta columns and an eta-less theta's column are exact", {
  errs <- .int_sens_col_errs(one_cmt_kappa_fn, bolus, times)
  expect_cols_exact(errs, c("eta:eta.cl", "eta:eta.v", "theta:tsc"))
})

test_that("linCmt: eta columns and an eta-less theta's column are exact", {
  # linCmt has no state sensitivities; these columns come from rxode2's linCmtB
  # chain rule. (This is the case nlmixr2est's augmented outer model cannot build.)
  errs <- .int_sens_col_errs(one_cmt_lincmt_kappa_fn, bolus, times)
  expect_cols_exact(errs, c("eta:eta.cl", "eta:eta.v", "theta:tsc"))
})

test_that("multi-endpoint linCmt: per-endpoint routing + all columns are exact", {
  # A linCmt PK endpoint (cp) + an algebraic PD endpoint (eff). The endpoint
  # pseudo-compartments are numbered after linCmt's implicit `central` (rxState),
  # NOT after the empty rxStateOde -- getting that wrong mis-routes the
  # CMT-conditional rx_pred_/rx_f1_ columns. nlmixr2est's inner model returns NULL
  # for this model, so this is the ONLY analytical sensitivity path for it. Check
  # every eta + THETA column, per endpoint, against a finite difference of the
  # emitted model's own rx_pred_.
  ui  <- suppressMessages(rxode2::rxode2(two_endpoint_lincmt_fn))
  unp <- admixr2:::.admUnpairedThetas(ui)
  expect_setequal(unp, c("tv", "tec50"))

  sm <- suppressMessages(admixr2:::.admBuildThetaSens(ui, unp))
  # The whole point of this test is that the emitter does NOT bail here.
  expect_false(is.null(sm))
  expect_setequal(names(sm$theta_sens_cols), c("tv", "tec50"))

  th_rows <- ui$iniDf[!is.na(ui$iniDf$ntheta), ]
  th <- stats::setNames(th_rows$est, paste0("THETA[", th_rows$ntheta, "]"))
  et <- c(`ETA[1]` = 0.1)

  sol <- function(params, cmt) {
    ev <- rxode2::et(amt = 100) |> rxode2::et(times, cmt = cmt)
    d  <- rxode2::rxSolve(sm$mod, params = params, events = ev,
                          returnType = "data.frame", addDosing = FALSE,
                          atol = 1e-12, rtol = 1e-12)
    d[d$time > 0, , drop = FALSE]
  }
  scal <- function(x, y) max(abs(x - y)) / max(max(abs(y)), 1e-8)

  # Routing: cp is a concentration (>> 1); eff = cp/(ec50+cp) lives in (0, 1).
  expect_gt(mean(sol(c(th, et), "cp")$rx_pred_),  1)
  expect_lt(max(sol(c(th, et), "eff")$rx_pred_),  1)

  for (cmt in c("cp", "eff")) {
    base <- sol(c(th, et), cmt)
    hh   <- 1e-6
    # eta column
    ep <- et; ep[1] <- ep[1] + hh; em <- et; em[1] <- em[1] - hh
    fd <- (sol(c(th, ep), cmt)$rx_pred_ - sol(c(th, em), cmt)$rx_pred_) / (2 * hh)
    expect_lt(scal(base[[sm$sens_cols[1]]], fd), TOL,
              label = sprintf("eta col, endpoint %s", cmt))
    # THETA columns
    for (tnm in names(sm$theta_sens_cols)) {
      k  <- th_rows$ntheta[match(tnm, th_rows$name)]; pn <- paste0("THETA[", k, "]")
      tp <- th; tp[pn] <- tp[pn] + hh; tm <- th; tm[pn] <- tm[pn] - hh
      fd <- (sol(c(tp, et), cmt)$rx_pred_ - sol(c(tm, et), cmt)$rx_pred_) / (2 * hh)
      expect_lt(scal(base[[sm$theta_sens_cols[[tnm]]]], fd), TOL,
                label = sprintf("THETA col %s, endpoint %s", tnm, cmt))
    }
  }
})

test_that("shared eta: tcl and tv each get their OWN exact column", {
  # eta.cl drives both cl and v, so d(pred)/d(eta.cl) collects a path through v
  # that d(pred)/d(tcl) does not have -- reusing the eta column would be wrong.
  errs <- .int_sens_col_errs(one_cmt_shared_eta_fn, bolus, times)
  expect_cols_exact(errs, c("eta:eta.cl", "eta:eta.ka", "theta:tcl", "theta:tv"))
})

test_that("covariate, expit-bounded, and parameter-dependent-IC columns are exact", {
  # One rich model carries all three independent eta-less theta types, so the whole
  # set is verified in a single compile (see one_cmt_feat_fn).
  errs <- .int_sens_col_errs(one_cmt_feat_fn, bolus, times)
  expect_cols_exact(errs, c("eta:eta.cl", "eta:eta.v",
                            "theta:bwt", "theta:tfr", "theta:tinit"))
})

# ---- Dosing modifiers: the regression this file was written for --------------
# These are ZERO without eventSens = "jump" (nlmixr2est's inner model carries no
# event/dose-parameter sensitivities; FOCEI finite-differences them separately).
#
# Jump support is version-dependent: rxode2 5.1.2 has no lag() jumps, 5.1.3 does.
# .admJumpCovers() feature-detects it and REFUSES the sens model when a modifier
# the model uses cannot be differentiated -- so the estimators fall back to FD
# instead of silently using a zero column. The tests below skip in that case.

test_that("a dosing modifier rxode2 cannot differentiate refuses the sens model", {
  # Contract: sens model or FD, never a silently-zero column. Both branches are
  # correct; which one runs depends on the installed rxode2.
  ui <- suppressMessages(rxode2::rxode2(one_cmt_dose_fn))
  sm <- suppressMessages(admixr2:::.admLoadSensModel(ui))

  expect_setequal(admixr2:::.admDoseMods(ui$loadPruneSens), c("f", "lag"))

  if (is.null(sm) || !identical(sm$type, "dirs")) {
    # refused: this rxode2 cannot differentiate a modifier one of our directions
    # feeds, so there is no theta column for tlag and the estimators FD it. The
    # point of the contract is that we never hand back a zero column instead.
    expect_null(sm$theta_sens_cols)
  } else {
    # accepted -> the jump derivatives must be there for every modifier a
    # direction of ours actually feeds
    expect_true(admixr2:::.admJumpCovers(sm$mod, ui$loadPruneSens, sm$dirs))
    d0 <- rxode2::rxSolve(sm$mod,
      params = c(`THETA[1]` = log(5), `THETA[2]` = log(20), `THETA[3]` = log(1),
                 `THETA[4]` = 0.5, `THETA[5]` = log(0.3), `THETA[6]` = 0.1,
                 `ETA[1]` = 0.2, `ETA[2]` = 0.1),
      events = rxode2::et(amt = 100) |> rxode2::et(c(0.5, 1, 2, 4)),
      returnType = "data.frame", addDosing = FALSE)
    d0 <- d0[d0$time > 0, ]
    expect_false(all(d0[[sm$sens_cols[2]]] == 0))               # eta.f in f()
    expect_false(all(d0[[sm$theta_sens_cols[["tlag"]]]] == 0))  # tlag in alag()
  }
})

test_that("bioavailability f() and lag time: eta AND theta columns are exact", {
  errs <- .int_sens_col_errs(one_cmt_dose_fn, bolus, times)
  # eta.f drives f(depot); tlag drives alag(depot); tka/tv are ordinary eta-less
  expect_cols_exact(errs, c("eta:eta.cl", "eta:eta.f",
                            "theta:tv", "theta:tka", "theta:tlag"))
})

test_that("modelled infusion rate(): theta column is exact", {
  # rate = -1 in the event table selects the modelled rate
  ev <- rxode2::et(amt = 100, rate = -1) |> rxode2::et(times)
  errs <- .int_sens_col_errs(one_cmt_rate_fn, ev, times)
  expect_cols_exact(errs, c("eta:eta.cl", "eta:eta.v", "theta:trate"))
})

test_that("modelled infusion dur(): theta column is exact", {
  # rate = -2 in the event table selects the modelled duration
  ev <- rxode2::et(amt = 100, rate = -2) |> rxode2::et(times)
  errs <- .int_sens_col_errs(one_cmt_dur_fn, ev, times)
  expect_cols_exact(errs, c("eta:eta.cl", "eta:eta.v", "theta:tdur"))
})

# ---- Cross-check against nlmixr2est's own implementation ---------------------

test_that("our emitted model agrees with nlmixr2est's inner model (eta columns + pred)", {
  # Two INDEPENDENT implementations of the same derivative: admixr2 emits its own
  # direction-set model (.admBuildThetaSens), nlmixr2est builds the focei inner
  # model. Their prediction and their eta sensitivity columns must coincide.
  # (There is no counterpart for the THETA columns -- nlmixr2est's inner model
  # does not emit them, which is the whole point of the emitter.)
  times <- c(0.5, 1, 2, 4, 8, 12, 24)
  ev    <- rxode2::et(amt = 100) |> rxode2::et(times)

  for (nm in c("one_cmt_kappa_fn", "one_cmt_lincmt_kappa_fn", "one_cmt_dose_fn",
               "one_cmt_feat_fn")) {
    ui <- suppressMessages(rxode2::rxode2(get(nm)))
    unp <- admixr2:::.admUnpairedThetas(ui)
    n_eta <- sum(!is.na(ui$iniDf$neta1) & ui$iniDf$neta1 == ui$iniDf$neta2 & !ui$iniDf$fix)

    ours   <- suppressMessages(admixr2:::.admBuildThetaSens(ui, unp))
    theirs <- suppressMessages(admixr2:::.admSensFromInner(ui, NULL, n_eta, tempfile()))
    # Either side may legitimately refuse on an rxode2 that cannot differentiate a
    # dosing modifier it needs (5.1.2 has no lag() jumps) -- the comparison only
    # means something when both exist.
    if (is.null(ours) || is.null(theirs)) next

    th_rows <- ui$iniDf[!is.na(ui$iniDf$ntheta), ]
    th <- stats::setNames(th_rows$est, paste0("THETA[", th_rows$ntheta, "]"))
    et <- stats::setNames(rep(0.1, n_eta), paste0("ETA[", seq_len(n_eta), "]"))
    sol <- function(mod) {
      need <- setdiff(rxode2::rxModelVars(mod)$params, c(names(th), names(et)))
      cv   <- if (length(need)) stats::setNames(rep(70, length(need)), need) else NULL
      d <- rxode2::rxSolve(mod, params = c(th, et, cv), events = ev,
                           returnType = "data.frame", addDosing = FALSE,
                           atol = 1e-12, rtol = 1e-12)
      d[d$time > 0, , drop = FALSE]
    }
    a <- sol(ours$mod); b <- sol(theirs$mod)
    scal <- function(x, y) max(abs(x - y)) / max(max(abs(y)), 1e-8)
    expect_lt(scal(a$rx_pred_, b$rx_pred_), 1e-6)
    for (j in seq_len(n_eta))
      expect_lt(scal(a[[ours$sens_cols[j]]], b[[theirs$sens_cols[j]]]), 1e-6)
    # This loop builds a fresh model each iteration; reclaim with rxode2's own
    # idiom so the registry does not accumulate all five at once.
    gc(FALSE); rxode2::rxUnloadAll()
  }
})

test_that("delay() and transit() models get exact columns", {
  # .rxSens() accumulates the delay-sensitivity augmentation itself, and a
  # non-constant delay's pre-history (s$..pastLines) is emitted too.
  times <- c(0.5, 1, 2, 4, 8, 12, 24)
  ev    <- rxode2::et(amt = 100) |> rxode2::et(times)

  errs <- .int_sens_col_errs(one_cmt_delay_fn, ev, times)
  expect_cols_exact(errs, c("eta:eta.cl", "eta:eta.v", "theta:tka"))

  errs <- .int_sens_col_errs(one_cmt_transit_fn, ev, times)
  expect_cols_exact(errs, c("eta:eta.cl", "eta:eta.v",
                            "theta:tka", "theta:tn", "theta:tmtt"))
})

test_that("d(pred)/d(lag) is undefined AT the dose boundary (documented, not a bug)", {
  # At an observation coinciding EXACTLY with the lagged dose time the derivative
  # DOES NOT EXIST. The analytic jump returns a one-sided value and a central
  # difference straddles the discontinuity, so the two disagree -- but the exact
  # values are implementation-defined and rxode2-version-dependent (older builds
  # returned the post-jump value with FD ~ analytic/2; the 5.1.3 r-universe build
  # returns 0 at the pre-dose side). So assert only what is actually guaranteed:
  # the column is EXACT everywhere off the boundary, and analytic and FD DISAGREE
  # at it. Keep lag/rate values off the observation grid in every other test.
  lag_on_grid_fn <- function() {
    ini({tcl <- log(5); tv <- log(20); tka <- log(1); tlag <- log(2); add.err <- 0.1
         eta.cl ~ 0.09; eta.v ~ 0.04})
    model({
      cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v); ka <- exp(tka)
      d/dt(depot)   <- -ka * depot
      alag(depot)   <- exp(tlag)         # lag = 2 h, and t = 2 IS an observation
      d/dt(central) <-  ka * depot - (cl / v) * central
      cp <- central / v
      cp ~ add(add.err)
    })
  }
  ui <- suppressMessages(rxode2::rxode2(lag_on_grid_fn))
  sm <- suppressMessages(admixr2:::.admLoadSensModel(ui))
  skip_if(is.null(sm) || is.null(sm$theta_sens_cols),
          "rxode2 build has no dosing-modifier (jump) sensitivities")
  th <- stats::setNames(c(log(5), log(20), log(1), log(2), 0.1),
                        paste0("THETA[", 1:5, "]"))
  n_all <- 2L + length(sm$theta_sens_cols)
  et <- stats::setNames(rep(0, n_all), paste0("ETA[", seq_len(n_all), "]"))
  et[1:2] <- 0.1
  sol <- function(th) {
    d <- rxode2::rxSolve(sm$mod, params = c(th, et), events = bolus,
                         returnType = "data.frame", addDosing = FALSE,
                         atol = 1e-12, rtol = 1e-12)
    d[d$time > 0, , drop = FALSE]
  }
  h  <- 1e-5
  d0 <- sol(th)
  tp <- th; tp[4] <- tp[4] + h
  tm <- th; tm[4] <- tm[4] - h
  fd  <- (sol(tp)$rx_pred_ - sol(tm)$rx_pred_) / (2 * h)
  ana <- d0[[sm$theta_sens_cols[["tlag"]]]]
  at_boundary <- d0$time == 2

  # the FD reference and the analytic column DISAGREE at the boundary (the
  # derivative is undefined there) -- the specific values are version-dependent
  expect_gt(abs(ana[at_boundary] - fd[at_boundary]), 1e-2)
  # ... and the column is EXACT everywhere else (the real guarantee)
  expect_lt(max(abs(ana[!at_boundary] - fd[!at_boundary]) /
                pmax(abs(fd[!at_boundary]), 1e-8)), 1e-4)
})
