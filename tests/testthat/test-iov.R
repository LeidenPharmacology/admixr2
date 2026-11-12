# Inter-occasion variability: the bits that need no solve.

test_that(".admOccStarts starts an occasion at the dose that opens it", {
  tt <- c(2, 6, 12, 26, 50, 98, 102, 110)
  lab <- c(1, 1, 1, 1, 1, 2, 2, 2)
  dose <- seq(0, 144, 24)
  # two doses (72, 96) fall between the occasions: ambiguous, so say which
  expect_error(admixr2:::.admOccStarts(tt, lab, c(1, 2), dose), "occasion_starts")
  st <- admixr2:::.admOccStarts(tt, lab, c(1, 2), dose, starts = 96)
  expect_equal(st[2L], 96)
  expect_lt(st[1L], 0)                           # first level open on the left
  expect_error(admixr2:::.admOccStarts(tt, lab, c(1, 2), dose, starts = 40), "must fall after")
  expect_error(admixr2:::.admOccStarts(tt, lab, c(1, 2), dose, starts = c(96, 97)), "one time per")
  # exactly one dose between the occasions: it opens the new one
  expect_equal(admixr2:::.admOccStarts(c(2, 26, 98), c(1, 1, 2), c(1, 2), c(0, 96))[2L], 96)
  # the reviewer's sparse case is refused rather than guessed
  expect_error(admixr2:::.admOccStarts(c(2, 150), c(1, 2), c(1, 2), seq(0, 144, 24)), "could have opened")
  # no et() dose table: the first observation of the occasion
  expect_equal(admixr2:::.admOccStarts(tt, lab, c(1, 2), NULL)[2L], 98)
  # a dose inside the earlier occasion's observations does not open the next one
  expect_equal(admixr2:::.admOccStarts(c(2, 6, 10, 14), c(1, 1, 2, 2), c(1, 2),
                                       c(0, 4))[2L], 10)
})

test_that(".admOccStarts gives an unobserved level an empty interval", {
  st <- admixr2:::.admOccStarts(c(2, 6, 98), c(1, 1, 3), c(1, 2, 3), NULL)
  expect_equal(st[2L], st[3L])                   # level 2: [t2, t3) is empty
  # a study that starts at level 2 never uses level 1
  st2 <- admixr2:::.admOccStarts(c(2, 6), c(2, 2), c(1, 2), NULL)
  expect_lt(st2[2L], 0)
  # no labels: one occasion for ever
  expect_true(all(admixr2:::.admOccStarts(c(2, 6), NULL, c(1, 2, 3))[-1L] >= 1e10))
})

test_that(".admOccStarts refuses labels that go back or are unknown", {
  expect_error(admixr2:::.admOccStarts(c(1, 2, 3), c(1, 2, 1), c(1, 2)),
               "earlier occasion")
  expect_error(admixr2:::.admOccStarts(c(1, 2), c(1, 5), c(1, 2)), "occasion levels")
})

test_that(".admEtDoseTimes expands addl and ii", {
  skip_if_not_installed("rxode2")
  ev <- rxode2::et(amt = 100, ii = 24, addl = 2)
  expect_equal(admixr2:::.admEtDoseTimes(ev), c(0, 24, 48))
  expect_null(admixr2:::.admEtDoseTimes(100))
})

test_that("admData(occasions =) adds the occasion levels", {
  d <- admData(occasions = c(occ = 3))
  expect_true("occ" %in% names(d))
  expect_setequal(unique(d$occ), 1:3)
  expect_error(admData(occasions = c(occ = 1)), ">= 2")
  expect_error(admData(occasions = 3), "named")
})

test_that("admStudy(occasions =) wants one label per time", {
  expect_error(admStudy(E = c(1, 2), sd = c(1, 1), n = 10, times = c(1, 2),
                        dose = 100, occasions = 1), "one numeric occasion label")
  s <- admStudy(E = c(1, 2), sd = c(1, 1), n = 10, times = c(1, 2),
                dose = 100, occasions = c(1, 2))
  expect_equal(s$occasions, c(1, 2))
})

test_that("a fixed diagonal eta stays in the model but out of the optimizer", {
  skip_if_not_installed("rxode2")
  m <- function() {
    ini({ tcl <- log(2.7); eta.cl ~ 0.09; eta.fx ~ fix(1); add.sd <- 0.3 })
    model({ cl <- exp(tcl + eta.cl + 0.5 * eta.fx); v <- 30
            d/dt(centr) = -cl / v * centr; cp = centr / v; cp ~ add(add.sd) })
  }
  ui <- suppressWarnings(suppressMessages(rxode2::rxode2(m)))
  pi <- admixr2:::.admParseIniDf(ui$iniDf, ui)
  expect_equal(pi$n_eta, 2L)
  expect_equal(pi$eta_col_names, c("eta.cl", "eta.fx"))
  expect_equal(length(pi$omega_par), 1L)         # only eta.cl is optimized
  pars <- admixr2:::.admUnpack(admixr2:::.admBuildOptVec(pi)$p0, pi)
  expect_equal(unname(pars$omega), diag(c(0.09, 1)))
})

test_that("an eta not named eta.* keeps its own name as the solve column", {
  skip_if_not_installed("rxode2")
  m <- function() {
    ini({ tcl <- log(2.7); b1 ~ 0.09; add.sd <- 0.3 })
    model({ cl <- exp(tcl + b1); v <- 30
            d/dt(centr) = -cl / v * centr; cp = centr / v; cp ~ add(add.sd) })
  }
  ui <- suppressMessages(rxode2::rxode2(m))
  expect_equal(admixr2:::.admParseIniDf(ui$iniDf, ui)$eta_col_names, "b1")
})

test_that("an IOV model's occasion variable and levels are read off the rewrite", {
  skip_if_not_installed("rxode2")
  m <- function() {
    ini({ tcl <- log(2.7); eta.cl ~ 0.09; add.sd <- 0.3 })
    model({ rx.iov <- 0.2 * (eta.cl * (occ == 1) + eta.cl * (occ == 3))
            cl <- exp(tcl + rx.iov); v <- 30
            d/dt(centr) = -cl / v * centr; cp = centr / v; cp ~ add(add.sd) })
  }
  ui <- suppressMessages(rxode2::rxode2(m))
  o <- admixr2:::.admIovOcc(ui)
  expect_equal(o$var, "occ")
  expect_equal(o$levels, c(1, 3))
  m2 <- function() {
    ini({ tcl <- log(2.7); eta.cl ~ 0.09; add.sd <- 0.3 })
    model({ cl <- exp(tcl + eta.cl); v <- 30
            d/dt(centr) = -cl / v * centr; cp = centr / v; cp ~ add(add.sd) })
  }
  expect_null(admixr2:::.admIovOcc(suppressMessages(rxode2::rxode2(m2))))
})

test_that("IOV is flagged as unidentified unless a study separates the occasions", {
  s_occ  <- list(occasions = c(1, 1, 2, 2), times = 1:4, V = diag(4))                 # full V
  s_sd   <- list(occasions = c(1, 1, 2, 2), times = 1:4, V = rep(1, 4))               # SDs only
  s_one  <- list(occasions = c(1, 1, 1, 1), times = 1:4, V = diag(4))                 # one occasion
  s_none <- list(times = 1:4, V = diag(4))                                            # no occasions
  s_mod  <- list(occasions = c(1, 1, 2, 2), ui = "a model source")
  expect_equal(unname(admixr2:::.admIovInformative(
    list(s_occ, s_sd, s_one, s_none, s_mod))), c(TRUE, FALSE, FALSE, FALSE, TRUE))
  occ <- list(var = "occ", levels = c(1, 2))
  expect_warning(admixr2:::.admIovPad(list(a = s_sd, b = s_none), occ),
                 "only their sum is")
  expect_no_warning(admixr2:::.admIovPad(list(a = s_occ, b = s_none), occ))
})
