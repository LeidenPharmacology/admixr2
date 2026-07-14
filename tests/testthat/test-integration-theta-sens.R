# Tier 2: structural-theta sensitivities for UNPAIRED (non-mu-referenced) thetas.
#
# The sens model is augmented with one dummy eta per unpaired theta
# (.admBuildSensUi), so d(pred)/d(theta) comes out of the SAME solve as the eta
# sensitivities and no longer needs a finite-difference rxSolve. Setup in
# helper-integration.R (.int_theta_sens_setup): ODE + linCmt, each with the
# unpaired theta `tsc`.

skip_on_cran()
skip_if_not_installed("rxode2")

test_that("sens model carries a direction per eta plus one per unpaired theta", {
  env <- .int_theta_sens_setup()
  for (case in c("ode", "lin")) {
    e  <- env[[case]]
    sm <- e$sensModel
    expect_false(is.null(sm), info = case)
    expect_equal(e$unpaired, "tsc", info = case)
    expect_equal(sm$type, "dirs", info = case)   # our emitter, not the inner fallback

    # one direction per eta (mu-referenced thetas reuse theirs for free) plus one
    # for the unpaired theta -- tsc is THETA[3] in this model
    expect_equal(sm$dirs, c("ETA_1_", "ETA_2_", "THETA_3_"), info = case)
    expect_equal(sm$sens_cols, c("rx_f1_ETA_1_", "rx_f1_ETA_2_"), info = case)
    expect_equal(sm$theta_sens_cols, c(tsc = "rx_f1_THETA_3_"), info = case)
    expect_length(sm$sens_cols, e$pinfo$n_eta)
  }
})

test_that("linCmt sens model carries a theta column (linCmtB chain rule)", {
  env <- .int_theta_sens_setup()
  expect_true(isTRUE(env$lin$sensModel$is_lincmt))
  expect_false(is.null(env$lin$sensModel$theta_sens_cols))
})

test_that(".admGrad theta-sens vs CFD of .admNLL: ratio within 5% (ODE + linCmt)", {
  env <- .int_theta_sens_setup()
  for (case in c("ode", "lin")) {
    e     <- env[[case]]
    ratio <- e$g_ana / e$g_fd
    ok    <- is.finite(ratio) & abs(e$g_fd) > 1e-6
    expect_true(all(abs(ratio[ok] - 1) < 0.05),
      info = sprintf("%s: max ratio deviation %.4f (param %s)", case,
                     max(abs(ratio[ok] - 1)),
                     names(ratio[ok])[which.max(abs(ratio[ok] - 1))]))
  }
})

test_that("the unpaired theta's gradient is MORE accurate than the FD path it replaces", {
  # Reference: central FD of the NLL (accurate to ~1e-5 here). g_plain is the
  # path this feature replaces, at the production default (forward FD).
  env <- .int_theta_sens_setup()
  for (case in c("ode", "lin")) {
    e   <- env[[case]]
    k   <- which(names(e$p0) == "tsc")
    ref <- e$g_fd[[k]]
    err_sens  <- abs(e$g_ana[[k]]   - ref) / abs(ref)
    err_plain <- abs(e$g_plain[[k]] - ref) / abs(ref)
    expect_lt(err_sens, 1e-3)
    expect_lt(err_sens, err_plain)
  }
})

test_that("FD fallback still works when the sens model has no theta columns", {
  env <- .int_theta_sens_setup()
  for (case in c("ode", "lin")) {
    e     <- env[[case]]
    ratio <- e$g_plain / e$g_fd
    ok    <- is.finite(ratio) & abs(e$g_fd) > 1e-6
    expect_true(all(abs(ratio[ok] - 1) < 0.05),
      info = sprintf("%s: FD-fallback max ratio deviation %.4f", case,
                     max(abs(ratio[ok] - 1))))
  }
})

test_that(".adghGrad uses the theta columns and matches CFD of .adghNLL", {
  env <- .int_theta_sens_setup()
  for (case in c("ode", "lin")) {
    e     <- env[[case]]
    ratio <- e$g_adgh / e$g_adgh_fd
    ok    <- is.finite(ratio) & abs(e$g_adgh_fd) > 1e-6
    expect_true(all(abs(ratio[ok] - 1) < 0.05),
      info = sprintf("%s: adgh max ratio deviation %.4f (param %s)", case,
                     max(abs(ratio[ok] - 1)),
                     names(ratio[ok])[which.max(abs(ratio[ok] - 1))]))
  }
})

test_that("sensitivities of DOSING-MODIFIER parameters are not silently zero", {
  # Regression: nlmixr2est's inner model carries no event/dose-parameter
  # sensitivities (FOCEI finite-differences those separately via predNoLhs), so
  # reading its columns directly gave a sensitivity of exactly ZERO for any
  # parameter entering f() / lag() / rate() / dur(). .admSensFromUi now recompiles
  # the inner model with eventSens = "jump" (rxode2's analytic variational jumps
  # at dose times), which is what nlmixr2est also does for its augmented model.
  #
  # Covers BOTH directions: eta.f (a real eta -- this was wrong before the fix)
  # and tlag (an eta-less theta -- its dummy-eta column would be wrong too).
  # Model: one_cmt_dose_fn (helper-integration.R). Its lag (0.3 h) is deliberately
  # OFF the observation grid -- d(pred)/d(lag) does not exist at an observation
  # coinciding exactly with the lagged dose time (the analytic jump gives the
  # one-sided derivative; a central difference averages across the discontinuity
  # and returns exactly half), so an FD reference is meaningless there.
  ui <- suppressMessages(rxode2::rxode2(one_cmt_dose_fn))
  sm <- suppressMessages(admixr2:::.admLoadSensModel(ui))
  expect_false(is.null(sm))
  expect_true("tlag" %in% names(sm$theta_sens_cols))

  times <- c(0.5, 1, 2, 4, 8, 12, 24)
  ev    <- rxode2::et(amt = 100) |> rxode2::et(times)
  th    <- setNames(c(log(5), log(20), log(1), 0.5, log(0.3), 0.1),
                    paste0("THETA[", 1:6, "]"))
  et    <- setNames(c(0.2, 0.1), paste0("ETA[", 1:2, "]"))   # eta.cl, eta.f

  sol <- function(th, et) {
    d <- rxode2::rxSolve(sm$mod, params = c(th, et), events = ev,
                         returnType = "data.frame", addDosing = FALSE,
                         atol = 1e-10, rtol = 1e-10)
    d[d$time > 0, , drop = FALSE]
  }
  H  <- 1e-5
  d0 <- sol(th, et)

  # eta.f (real eta in f(depot)) -- was identically zero before eventSens="jump"
  ep <- et; ep[2] <- ep[2] + H
  em <- et; em[2] <- em[2] - H
  fd_etaf  <- (sol(th, ep)$rx_pred_ - sol(th, em)$rx_pred_) / (2 * H)
  ana_etaf <- d0[[sm$sens_cols[2]]]
  expect_false(all(ana_etaf == 0))
  expect_lt(max(abs(ana_etaf - fd_etaf) / pmax(abs(fd_etaf), 1e-8)), 1e-4)

  # tlag (eta-less theta in alag(depot)) via its own THETA direction
  k <- 5L                                   # THETA[5] = tlag
  tp <- th; tp[k] <- tp[k] + H
  tm <- th; tm[k] <- tm[k] - H
  fd_tlag  <- (sol(tp, et)$rx_pred_ - sol(tm, et)$rx_pred_) / (2 * H)
  ana_tlag <- d0[[sm$theta_sens_cols[["tlag"]]]]
  expect_false(all(ana_tlag == 0))
  expect_lt(max(abs(ana_tlag - fd_tlag) / pmax(abs(fd_tlag), 1e-8)), 1e-4)
})

test_that("the sens model reproduces the simulation model's prediction", {
  # The sens model's rx_pred_ is what .admGrad uses as cp_mat, so it must be the
  # same prediction the NLL is built from. Tolerance 1e-5, not exact: rxode2's
  # adaptive solver takes different steps for a system carrying sensitivity
  # compartments (the same effect CLAUDE.md documents for a changed output grid).
  # That is pre-existing -- .admGrad has always taken cp_mat from the sens solve.
  env   <- .int_theta_sens_setup()
  e     <- env$ode
  pinfo <- e$pinfo
  ui    <- e$ui

  rxMod <- admixr2:::.admLoadModel(ui)
  rxode2::rxLoad(rxMod)

  times <- c(0.5, 1, 2, 4)
  E_t   <- .one_cmt_mean(5, 20, 100, times)
  study <- admixr2:::.admNormaliseStudy(
    list(E = E_t, V = diag((0.3 * E_t)^2), n = 200L, times = times,
         ev = rxode2::et(amt = 100)), "s")
  study$ev_full <- study$ev |> rxode2::et(study$times)

  pars    <- admixr2:::.admUnpack(e$p0, pinfo)
  n_sim   <- 32L
  z       <- admixr2:::.admMakeZ(n_sim, pinfo, 1L, "sobol")[[1]]
  eta_mat <- z %*% t(pars$L)
  colnames(eta_mat) <- pinfo$eta_col_names
  plist   <- admixr2:::.admMakeParamsList(n_sim, pinfo, 1L)

  cp_sim <- admixr2:::.admSimulate(rxMod, pars$struct, pinfo$sigma_names, eta_mat,
                                   study, "cp", plist[[1]], 1L)
  out    <- admixr2:::.admSimulateSens(e$sensModel, pars$struct, pinfo$sigma_names,
                                       eta_mat, study, 1L)
  expect_equal(out$cp_mat, cp_sim, tolerance = 1e-5)

  expect_false(is.null(out$dtheta_list))
  expect_equal(names(out$dtheta_list), "tsc")
  expect_equal(dim(out$dtheta_list$tsc), c(n_sim, length(times)))

  # hiding the theta columns is what the nlmixr2est-inner fallback looks like ->
  # the estimators finite-difference those thetas
  plain <- e$sensModel; plain$theta_sens_cols <- NULL
  plain_out <- admixr2:::.admSimulateSens(plain, pars$struct, pinfo$sigma_names,
                                          eta_mat, study, 1L)
  expect_null(plain_out$dtheta_list)
})
