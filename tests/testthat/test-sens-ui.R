# Tier 1: the sens-model augmentation helpers (.admSubstSym, .admUnpairedThetas,
# .admBuildSensUi). Only .admBuildSensUi needs rxode2.

test_that(".admSubstSym substitutes a symbol everywhere it appears", {
  map <- list(tka = quote((tka + eta.admSens.tka)))

  # simple parameter line
  expect_equal(admixr2:::.admSubstSym(quote(ka <- exp(tka)), map),
               quote(ka <- exp((tka + eta.admSens.tka))))

  # every occurrence, including several in one expression (total derivative)
  expect_equal(admixr2:::.admSubstSym(quote(y <- tka * exp(tka)), map),
               quote(y <- (tka + eta.admSens.tka) * exp((tka + eta.admSens.tka))))

  # nested calls and d/dt() lines
  expect_equal(admixr2:::.admSubstSym(quote(d/dt(depot) <- -exp(tka) * depot), map),
               quote(d/dt(depot) <- -exp((tka + eta.admSens.tka)) * depot))

  # a call HEAD is never substituted (a theta is never a function name)
  expect_equal(admixr2:::.admSubstSym(quote(tka(x)), list(tka = quote(zzz))),
               quote(tka(x)))

  # untouched symbols and literals survive
  expect_equal(admixr2:::.admSubstSym(quote(v <- exp(tv) + 2), map),
               quote(v <- exp(tv) + 2))
})

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

test_that(".admBuildSensUi adds one dummy eta per unpaired theta", {
  skip_on_cran()
  skip_if_not_installed("rxode2")

  ui  <- suppressMessages(rxode2::rxode2(one_cmt_kappa_fn))
  aug <- suppressMessages(admixr2:::.admBuildSensUi(ui, "tsc"))
  expect_false(is.null(aug))

  ini0 <- ui$iniDf
  ini1 <- aug$iniDf
  eta0 <- ini0[!is.na(ini0$neta1), ]
  eta1 <- ini1[!is.na(ini1$neta1), ]

  # one extra eta, appended AFTER the real ones (so ETA[k] indices of the real
  # etas -- and every existing sens column -- are unchanged)
  expect_equal(nrow(eta1), nrow(eta0) + 1L)
  expect_equal(eta1$name[seq_len(nrow(eta0))], eta0$name)
  expect_equal(eta1$name[nrow(eta1)], "eta.admSens.tsc")

  # theta block untouched: same names, same order, same ntheta numbering
  th0 <- ini0[is.na(ini0$neta1), c("ntheta", "name")]
  th1 <- ini1[is.na(ini1$neta1), c("ntheta", "name")]
  expect_equal(th1, th0)

  # a theta that is not in the model cannot be substituted -> NULL, and the
  # caller falls back to the plain sens model + FD
  expect_null(suppressWarnings(admixr2:::.admBuildSensUi(ui, "not_a_theta")))
  # nothing to do -> NULL
  expect_null(admixr2:::.admBuildSensUi(ui, character(0)))
})
