test_that("correlated covariate collapse differentiates the scored objective", {
  skip_on_cran()
  skip_if_not_installed("rxode2")
  skip_if_not_installed("randtoolbox")
  f <- function() {
    ini({ b1 <- .2; b2 <- .1; a <- .2; eta.a ~ .09 })
    model({ x <- exp(eta.a)
            cv <- exp(b1 * W1 + b2 * W2)
            cp <- x + cv; cp ~ add(a) })
  }
  ui <- suppressMessages(rxode2::rxode2(f))
  pin <- admixr2:::.admParseIniDf(ui$iniDf, ui)
  pin$cov_integration <- "on"; pin$cov_nodes <- 11L; pin$n_nodes <- 5L
  pin$nDisplayProgress <- .Machine$integer.max
  cd <- list(W1 = list(mu = 0, sd = 1),
             W2 = list(mu = 0, sd = 1), cor = .8)
  s <- list(E = 2.2, V = matrix(.3), n = 100L, times = 1,
            ev = rxode2::et(amt = 0), cov_dist = cd)
  st <- admixr2:::.admBuildEvFull(admixr2:::.admFlattenStudies(
    list(s = admixr2:::.admNormaliseStudy(s, "s", "cp"))))
  st <- admixr2:::.admCheckCovariates(ui, pin, st, "adgh")
  expect_false(is.null(st[[1L]]$.adm_cov_collapse))
  expect_null(st[[1L]]$.adm_cov_joint)

  sens <- admixr2:::.admLoadSensModel(ui)
  rx <- admixr2:::.admLoadModel(ui)
  rxode2::rxLoad(rx)
  p <- admixr2:::.admBuildOptVec(pin)$p0
  grid <- admixr2:::.adghNodeGrid(5L, pin$n_eta)
  got <- admixr2:::.adghGradNLL(p, pin, st, sens, rx, "cp", grid, 1L)$grad
  ref <- admixr2:::.adghFDGrad(p, pin, st, rx, "cp", grid, 1L, 1e-4)
  expect_equal(got, ref, tolerance = 1e-6)

  # A FAILED SENSITIVITY SOLVE must reach the finite-difference fallback. The
  # pre-pass stored its results with `res[[i]] <- NULL`, which DELETES the slot
  # rather than emptying it: with one study that left a zero-length `res` and
  # .adghGradNLL() subscripted past the end instead of seeing the NULL it tests
  # for. With several studies it silently shifted every later result by one.
  testthat::local_mocked_bindings(
    .admSimulateSensMany = function(...) list(NULL), .package = "admixr2")
  prs <- admixr2:::.admUnpack(p, pin)
  pre <- admixr2:::.adghGradPre(prs, pin, st, sens, grid, prs$L, 1L, p)
  expect_length(pre$res, length(st))
  expect_null(pre$res[[1L]])
  expect_equal(admixr2:::.adghGradNLL(p, pin, st, sens, rx, "cp", grid, 1L)$grad,
               ref, tolerance = 1e-6)
})
