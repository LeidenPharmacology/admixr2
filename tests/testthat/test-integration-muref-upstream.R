# =============================================================================
# Upstream watch: muRefCurEval and the "mu referencing 3.0" style
# =============================================================================
#
# rxode2 (<= 5.1.4, at least) loses the TRANSFORM of a mu-referenced parameter
# when the mu-referenced expression uses a variable computed earlier in the
# model -- `wt70 <- log(WT/70); cl <- exp(tcl + tcov*wt70 + eta.cl)` -- even
# though muRefDataFrame still identifies the pair correctly. admixr2 reads
# muRefCurEval to back-transform structural thetas for REPORTING, so on such a
# model a log-scale estimate is printed untransformed.
#
# This file is written so that NO upstream version breaks it:
#
#   * the premise -- that the two spellings are the same model -- is asserted
#     unconditionally, because it must hold on every version;
#   * the version-dependent part BRANCHES. While upstream is broken the test
#     records that and skips with a message; once upstream reports the
#     transform, the branch flips and the back-transform is asserted to be
#     right. Either way the file passes, and the skip message is the signal
#     that the workaround is still needed.
#
# A standalone reprex with no admixr2 dependency, for filing upstream, is in
# validation/rxode2-muref-curEval-reprex.R.

.mu3_inline <- function() {
  ini({tcl <- log(4); tv <- log(30); tcov <- 0.75; eta.cl ~ 0.09; a <- 0.1})
  model({cl <- exp(tcl + tcov * log(WT / 70) + eta.cl); v <- exp(tv)
         d/dt(centr) <- -cl / v * centr; cp <- centr / v; cp ~ add(a)})
}
.mu3_precomputed <- function() {
  ini({tcl <- log(4); tv <- log(30); tcov <- 0.75; eta.cl ~ 0.09; a <- 0.1})
  model({wt70 <- log(WT / 70)
         cl <- exp(tcl + tcov * wt70 + eta.cl); v <- exp(tv)
         d/dt(centr) <- -cl / v * centr; cp <- centr / v; cp ~ add(a)})
}

test_that("the two mu-3.0 spellings are the SAME model (the premise)", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  ui1 <- suppressMessages(rxode2::rxode2(.mu3_inline))
  ui2 <- suppressMessages(rxode2::rxode2(.mu3_precomputed))
  ev   <- rxode2::et(rxode2::et(amt = 500), c(0.5, 1, 2, 4, 8, 12))
  pars <- c(tcl = log(4), tv = log(30), tcov = 0.75, eta.cl = 0.2, WT = 92)
  sim <- function(ui) {
    m <- rxode2::rxode2(ui$simulationModel)
    p <- pars[intersect(names(pars), m$params)]
    for (q in setdiff(m$params, names(p))) p[q] <- 0
    s <- rxode2::rxSolve(m, params = p, events = ev, returnType = "data.frame",
                         addDosing = FALSE, atol = 1e-12, rtol = 1e-12)
    s$cp[!is.na(s$cp)]
  }
  # bit-identical, not merely close: any divergence here means the comparison
  # below is no longer about metadata
  expect_identical(sim(ui1), sim(ui2))
  # ... and BOTH are recognised as mu-referenced, so only the transform is at
  # issue
  for (ui in list(ui1, ui2)) {
    md <- ui$muRefDataFrame
    expect_true(NROW(md) >= 1L)
    expect_true("tcl" %in% md$theta)
  }
})

test_that("a mu-referenced theta back-transforms correctly, or upstream is why not", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  ui1 <- suppressMessages(rxode2::rxode2(.mu3_inline))
  ui2 <- suppressMessages(rxode2::rxode2(.mu3_precomputed))
  ce  <- function(ui) { d <- ui$muRefCurEval
    as.character(d$curEval[d$parameter == "tcl"])[1L] }
  tr  <- function(ui) {
    admixr2:::.admParseIniDf(ui$iniDf, ui)$struct_transforms[["tcl"]] }

  # the inline spelling must always work -- if THIS regresses it is admixr2's
  # problem, not upstream's
  expect_identical(ce(ui1), "exp")
  expect_equal(admixr2:::.admBackTransform(log(4), tr(ui1)), 4)

  # THE REPORTED VALUE MUST BE RIGHT EITHER WAY, and that is what is asserted --
  # the suite used to pin log(4) = 1.386 as the answer for a parameter whose
  # value is 4, on the ordinary mu-3.0 spelling, with the corrective branch
  # unreachable. A green light over a live, user-visible reporting bug.
  # .admParseIniDf() now reads the transform off ui$lstExpr whenever upstream
  # reports "" (.admCurEvalFromModel), so this holds on every rxode2 in play.
  expect_identical(admixr2:::.admCurEvalFromModel(ui2, "tcl")$curEval, "exp")
  # A COVARIATE COEFFICIENT IS NOT WRAPPED BY THE TRANSFORM. rxode2 blanks
  # curEval for every theta in the expression, coefficients included, and
  # `tcov` in exp(tcl + tcov*wt70 + eta.cl) is 0.75, not exp(0.75). The
  # fallback requires a mu-reference AND a bare additive term, so it declines.
  expect_identical(admixr2:::.admCurEvalFromModel(ui2, "tcov")$curEval, "")
  .tr2 <- admixr2:::.admParseIniDf(ui2$iniDf, ui2)$struct_transforms[["tcov"]]
  expect_equal(admixr2:::.admBackTransform(0.75, .tr2), 0.75)
  expect_equal(admixr2:::.admBackTransform(log(4), tr(ui2)), 4)
  # Upstream's own answer is still recorded, so the day it starts reporting the
  # transform the fallback is visibly no longer load-bearing.
  expect_true(ce(ui2) %in% c("", "exp"))
})
