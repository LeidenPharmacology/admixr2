# .admCovCols() appends each study's covariates in that study's own name order,
# so the positional rbind behind the batched solve gave the second study the
# first study's column names and solved it at swapped covariate values. Every
# dimension stayed valid, so no fallback saw it: an ordinary two-study GH fit
# moved from 5659.727 to 6243.142 with nothing to report.

test_that("batching aligns covariate columns by name, not position", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  f <- function() {
    ini({ tcl <- log(3.2); tv <- log(21); bwt <- .6; bage <- .2
          eta.cl ~ .09; add.err <- .6 })
    model({ cl <- exp(tcl + eta.cl)*(WT/70)^bwt*(AGE/40)^bage
            v <- exp(tv); cp <- linCmt(); cp ~ add(add.err) }) }
  ui  <- suppressMessages(rxode2::rxode2(f))
  pin <- admixr2:::.admDriverPinfo(ui, adghControl(studies = list(), print = 0L))
  rx  <- admixr2:::.admLoadModel(ui); rxode2::rxLoad(rx)
  mk <- function(cv, nm) admixr2:::.admNormaliseStudy(
    list(E = c(1, 1, 1), V = diag(3), n = 100L, times = c(1, 4, 8),
         ev = rxode2::et(amt = 100), cov = cv), nm, "cp")
  # SAME covariates, other order
  st <- admixr2:::.admBuildEvFull(admixr2:::.admFlattenStudies(list(
    a = mk(list(WT = 70, AGE = 40), "a"),
    b = mk(list(AGE = 60, WT = 90), "b"))))
  gr <- admixr2:::.adghNodeGrid(5L, pin$n_eta)
  p  <- admixr2:::.admBuildOptVec(pin)$p0

  # the batched objective against the per-study one it replaced
  got <- suppressWarnings(admixr2:::.adghNLL(p, pin, st, rx, "cp", gr, 1L))
  ref <- testthat::with_mocked_bindings(
    suppressWarnings(admixr2:::.adghNLL(p, pin, st, rx, "cp", gr, 1L)),
    .admSimulateMany = function(...) NULL, .package = "admixr2")
  expect_equal(got, ref, tolerance = 1e-12)

  # and the gradient, analytical against finite differences
  sens <- admixr2:::.admLoadSensModel(ui)
  expect_equal(
    admixr2:::.adghGradNLL(p, pin, st, sens, rx, "cp", gr, 1L)$grad,
    admixr2:::.adghFDGrad(p, pin, st, rx, "cp", gr, 1L, 1e-4),
    tolerance = 1e-5)
})

test_that(".admRbindParams binds by name and declines what it cannot reconcile", {
  # differing column SETS: no safe bind, so the caller solves per study
  expect_null(admixr2:::.admRbindParams(list(
    matrix(1, 2, 2, dimnames = list(NULL, c("WT", "AGE"))),
    matrix(1, 2, 2, dimnames = list(NULL, c("WT", "CRCL"))))))
  expect_null(admixr2:::.admRbindParams(list(
    matrix(1, 2, 2, dimnames = list(NULL, c("WT", "AGE"))),
    matrix(1, 2, 2))))
  # same names in another order are bound BY NAME
  al <- admixr2:::.admRbindParams(list(
    matrix(c(1, 2, 3, 4), 2, 2, dimnames = list(NULL, c("WT", "AGE"))),
    matrix(c(5, 6, 9, 8), 2, 2, dimnames = list(NULL, c("AGE", "WT")))))
  expect_equal(al[["WT"]],  c(1, 2, 9, 8))
  expect_equal(al[["AGE"]], c(3, 4, 5, 6))
})
