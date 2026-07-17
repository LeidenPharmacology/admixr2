# Tier 1: the sens-model direction set (.admUnpairedThetas, .admMuRefPairs) and
# the emitted model (.admBuildThetaSens -- needs rxode2).

test_that(".admUnpairedThetas mirrors struct_has_eta", {
  ini <- data.frame(
    ntheta = c(1L, 2L, 3L, 4L, NA, NA),
    neta1  = c(NA, NA, NA, NA, 1L, 2L),
    neta2  = c(NA, NA, NA, NA, 1L, 2L),
    name   = c("tcl", "tv", "tka", "add.err", "eta.cl", "eta.v"),
    est    = c(1, 3, 0, 0.1, 0.09, 0.04),
    fix    = FALSE,
    err    = c(NA, NA, NA, "add", NA, NA),
    stringsAsFactors = FALSE)
  # tcl/tv mu-referenced; tka is not; add.err is a sigma (never a direction)
  ui <- list(iniDf = ini,
             muRefDataFrame = data.frame(theta = c("tcl", "tv"),
                                         eta = c("eta.cl", "eta.v"),
                                         stringsAsFactors = FALSE))
  expect_equal(admixr2:::.admUnpairedThetas(ui), "tka")

  # no mu-referencing at all (e.g. cl <- exp(tcl) * exp(eta.cl)): every struct
  # theta is unpaired, and each gets its own direction
  ui0 <- list(iniDf = ini, muRefDataFrame = data.frame())
  expect_equal(admixr2:::.admUnpairedThetas(ui0), c("tcl", "tv", "tka"))

  # a FIXED theta is not estimated -> no direction
  ini_fix <- ini; ini_fix$fix[3] <- TRUE
  expect_equal(admixr2:::.admUnpairedThetas(list(iniDf = ini_fix,
                                                 muRefDataFrame = data.frame())),
               c("tcl", "tv"))
})

test_that(".admMuRefPairs drops a SHARED eta (its theta cannot reuse the eta column)", {
  ini <- data.frame(
    ntheta = c(1L, 2L, 3L, NA, NA),
    neta1  = c(NA, NA, NA, 1L, 2L),
    neta2  = c(NA, NA, NA, 1L, 2L),
    name   = c("tcl", "tv", "add.err", "eta.cl", "eta.ka"),
    est    = c(1, 3, 0.1, 0.09, 0.04),
    fix    = FALSE,
    lower  = -Inf, upper = Inf,
    err    = c(NA, NA, "add", NA, NA),
    stringsAsFactors = FALSE)
  mrd <- data.frame(theta = c("tcl", "tv"), eta = c("eta.cl", "eta.ka"),
                    stringsAsFactors = FALSE)

  # eta.cl used once -> tcl may reuse its column; eta.ka used once -> tv likewise
  ui_ok <- list(iniDf = ini, muRefDataFrame = mrd,
                lstExpr = list(quote(cl <- exp(tcl + eta.cl)),
                               quote(ka <- exp(tv + eta.ka))))
  expect_equal(admixr2:::.admMuRefPairs(ui_ok)$theta, c("tcl", "tv"))
  expect_equal(admixr2:::.admUnpairedThetas(ui_ok), character(0))

  # eta.cl now ALSO appears in v: d(pred)/d(eta.cl) picks up a path through v that
  # d(pred)/d(tcl) does not have, so the reuse is invalid -> tcl becomes unpaired
  # and gets its own (dummy-eta) direction.
  ui_shared <- list(iniDf = ini, muRefDataFrame = mrd,
                    lstExpr = list(quote(cl <- exp(tcl + eta.cl)),
                                   quote(v  <- exp(tv + eta.cl)),
                                   quote(ka <- exp(tv + eta.ka))))
  expect_equal(admixr2:::.admMuRefPairs(ui_shared)$theta, "tv")
  expect_equal(admixr2:::.admUnpairedThetas(ui_shared), "tcl")

  # no model text (mock ui) -> the guard cannot run; muRefDataFrame is trusted,
  # exactly as before this feature existed
  ui_blind <- list(iniDf = ini, muRefDataFrame = mrd)
  expect_equal(admixr2:::.admMuRefPairs(ui_blind)$theta, c("tcl", "tv"))
  expect_equal(admixr2:::.admUnpairedThetas(ui_blind), character(0))

  # information exists and says NOTHING is paired (non-mu-referenced model) ->
  # a zero-row frame, NOT NULL. NULL would send struct_eta_idx to its identity
  # fallback and double-count the eta path against the theta columns.
  ui_none <- list(iniDf = ini,
                  muRefDataFrame = mrd[0, , drop = FALSE],
                  lstExpr = list(quote(cl <- exp(tcl) * exp(eta.cl))))
  expect_false(is.null(admixr2:::.admMuRefPairs(ui_none)))
  expect_equal(nrow(admixr2:::.admMuRefPairs(ui_none)), 0L)
  expect_equal(admixr2:::.admUnpairedThetas(ui_none), c("tcl", "tv"))
  pinfo_none <- admixr2:::.admParseIniDf(ini, ui_none)
  expect_true(all(!pinfo_none$struct_has_eta))
  expect_true(all(is.na(pinfo_none$struct_eta_idx)))   # no identity fallback
})

test_that("struct_has_eta honours the shared-eta guard", {
  ini <- data.frame(
    ntheta = c(1L, 2L, NA),
    neta1  = c(NA, NA, 1L),
    neta2  = c(NA, NA, 1L),
    name   = c("tcl", "tv", "eta.cl"),
    est    = c(1, 3, 0.09),
    fix    = FALSE,
    lower  = -Inf, upper = Inf,
    err    = c(NA, NA, NA),
    stringsAsFactors = FALSE)
  ui <- list(iniDf = ini,
             muRefDataFrame = data.frame(theta = "tcl", eta = "eta.cl",
                                         stringsAsFactors = FALSE),
             lstExpr = list(quote(cl <- exp(tcl + eta.cl)),
                            quote(v  <- exp(tv + eta.cl))))   # eta.cl shared
  pinfo <- admixr2:::.admParseIniDf(ini, ui)
  # tcl must NOT be treated as mu-referenced: the eta-reuse identity fails
  expect_false(pinfo$struct_has_eta[["tcl"]])
  expect_true(is.na(pinfo$struct_eta_idx[[1L]]))
})

test_that(".admBuildThetaSens emits a direction per eta plus one per unpaired theta", {
  skip_on_cran()
  skip_if_not_installed("rxode2")

  ui <- suppressMessages(rxode2::rxode2(one_cmt_kappa_fn))   # tsc = THETA[3], no eta
  sm <- suppressMessages(admixr2:::.admBuildThetaSens(ui, "tsc"))
  expect_false(is.null(sm))

  # mu-referenced tcl/tv reuse their etas' directions -- only tsc gets its own
  expect_equal(sm$dirs, c("ETA_1_", "ETA_2_", "THETA_3_"))
  expect_equal(sm$sens_cols, c("rx_f1_ETA_1_", "rx_f1_ETA_2_"))
  expect_equal(sm$theta_sens_cols, c(tsc = "rx_f1_THETA_3_"))

  # the emitted model really carries those columns, plus the prediction
  lhs <- sm$mod$lhs
  expect_true(all(c("rx_pred_", "rx_f1_ETA_1_", "rx_f1_ETA_2_", "rx_f1_THETA_3_") %in% lhs))

  # ... and one variational compartment per direction per state
  st <- rxode2::rxModelVars(sm$mod)$state
  expect_true("rx__sens_central_BY_THETA_3___" %in% st)
  expect_true("rx__sens_central_BY_ETA_1___" %in% st)

  # no unpaired thetas -> eta directions only, no theta columns
  sm0 <- suppressMessages(admixr2:::.admBuildThetaSens(ui, character(0)))
  expect_equal(sm0$dirs, c("ETA_1_", "ETA_2_"))
  expect_null(sm0$theta_sens_cols)

  # a theta that is not in the model -> NULL, and the caller falls back to
  # nlmixr2est's inner model + FD
  expect_null(suppressMessages(admixr2:::.admBuildThetaSens(ui, "not_a_theta")))
})

test_that(".admBuildThetaSens works for linCmt (no ODE states -> linCmtB chain rule)", {
  skip_on_cran()
  skip_if_not_installed("rxode2")

  ui <- suppressMessages(rxode2::rxode2(one_cmt_lincmt_kappa_fn))
  # rxStateOde() is EMPTY for a linCmt model: there is nothing for .rxSens() to
  # augment (it would error). The emitter therefore drops the state-sensitivity
  # sum and D(pred, dir) alone resolves through rxode2's linCmtB derivative rules.
  expect_length(rxode2::rxStateOde(ui$loadPruneSens), 0L)

  sm <- suppressMessages(admixr2:::.admBuildThetaSens(ui, "tsc"))
  expect_false(is.null(sm))
  expect_equal(sm$theta_sens_cols, c(tsc = "rx_f1_THETA_3_"))

  mv <- rxode2::rxModelVars(sm$mod)
  expect_true(any(grepl("linCmtB", mv$model, fixed = TRUE)))
  # no variational compartments of OUR making. (rxode2 adds its own linCmtB
  # pseudo-compartments, rx__sens_central_BY_p1/v1 -- the same ones nlmixr2est's
  # inner model gets -- so assert on the direction-named ones specifically.)
  expect_false(any(grepl("_BY_(ETA|THETA)_", mv$state)))
})
