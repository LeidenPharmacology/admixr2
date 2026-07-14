# Tier 2: structural-theta sensitivities for UNPAIRED (non-mu-referenced) thetas.
#
# The sens model is augmented with one dummy eta per unpaired theta
# (.admBuildSensUi), so d(pred)/d(theta) comes out of the SAME solve as the eta
# sensitivities and no longer needs a finite-difference rxSolve. Setup in
# helper-integration.R (.int_theta_sens_setup): ODE + linCmt, each with the
# unpaired theta `tsc`.

skip_on_cran()
skip_if_not_installed("rxode2")

test_that("sens model is augmented with a dummy eta per unpaired theta", {
  env <- .int_theta_sens_setup()
  for (case in c("ode", "lin")) {
    e  <- env[[case]]
    sm <- e$sensModel
    expect_false(is.null(sm), info = case)
    expect_equal(e$unpaired, "tsc", info = case)

    # theta columns present, named by the unpaired theta
    expect_false(is.null(sm$theta_sens_cols), info = case)
    expect_equal(names(sm$theta_sens_cols), "tsc", info = case)

    # existing consumers unaffected: sens_cols is still exactly the real etas
    expect_length(sm$sens_cols, e$pinfo$n_eta)

    # the dummy eta is the (n_eta + 1)-th and must be pinned at 0 by the solve paths
    expect_equal(sm$dummy_eta_inner, paste0("ETA[", e$pinfo$n_eta + 1L, "]"), info = case)
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
  # ETA[1] eta.cl, ETA[2] eta.f, ETA[3..] dummy etas (tka, tlag) -- pinned at 0
  n_all <- 2L + length(sm$theta_sens_cols)
  et    <- setNames(rep(0, n_all), paste0("ETA[", seq_len(n_all), "]"))
  et[1:2] <- c(0.2, 0.1)

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

  # tlag (eta-less theta in alag(depot)) via its dummy-eta column
  k <- 5L                                   # THETA[5] = tlag
  tp <- th; tp[k] <- tp[k] + H
  tm <- th; tm[k] <- tm[k] - H
  fd_tlag  <- (sol(tp, et)$rx_pred_ - sol(tm, et)$rx_pred_) / (2 * H)
  ana_tlag <- d0[[sm$theta_sens_cols[["tlag"]]]]
  expect_false(all(ana_tlag == 0))
  expect_lt(max(abs(ana_tlag - fd_tlag) / pmax(abs(fd_tlag), 1e-8)), 1e-4)
})

test_that("the dummy eta does not change the prediction", {
  # At eta_dummy = 0 the augmented model must reproduce the PLAIN sens model's
  # prediction -- otherwise every NLL/moment computed from it is wrong.
  #
  # Compared against the plain SENS model, not the simulation model: the two
  # differ by ~1e-6 anyway because rxode2's adaptive solver takes different steps
  # for a system with sensitivity compartments (same effect CLAUDE.md documents
  # for a changed output grid). That is pre-existing -- .admGrad has always taken
  # cp_mat from the sens solve -- and is not what this test is about.
  env   <- .int_theta_sens_setup()
  e     <- env$ode
  pinfo <- e$pinfo
  ui    <- e$ui

  plain <- suppressMessages(admixr2:::.admSensFromUi(ui, ui, character(0)))
  expect_false(is.null(plain))
  expect_null(plain$theta_sens_cols)          # plain model has no theta columns

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

  aug_out   <- admixr2:::.admSimulateSens(e$sensModel, pars$struct, pinfo$sigma_names,
                                          eta_mat, study, 1L)
  plain_out <- admixr2:::.admSimulateSens(plain, pars$struct, pinfo$sigma_names,
                                          eta_mat, study, 1L)
  expect_equal(aug_out$cp_mat, plain_out$cp_mat, tolerance = 1e-8)
  # the real etas' sensitivities are unchanged by the augmentation too
  expect_equal(aug_out$dpred_list, plain_out$dpred_list, tolerance = 1e-8)

  expect_false(is.null(aug_out$dtheta_list))
  expect_equal(names(aug_out$dtheta_list), "tsc")
  expect_equal(dim(aug_out$dtheta_list$tsc), c(n_sim, length(times)))
  expect_null(plain_out$dtheta_list)           # plain model -> FD fallback
})
