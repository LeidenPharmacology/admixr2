test_that("Vector V expanded to diagonal matrix with method = 'var'", {
  s  <- list(E = c(1.0, 2.0), V = c(0.1, 0.2), n = 50L, times = c(1, 2))
  ns <- admixr2:::.admNormaliseStudy(s, "s1")

  expect_equal(ns$method, "var")
  expect_equal(ns$V, diag(c(0.1, 0.2)))
  expect_equal(ns$v_diag, c(0.1, 0.2))
})

test_that("Diagonal matrix auto-detected as method = 'var'", {
  s  <- list(E = c(1.0, 2.0), V = diag(c(0.1, 0.2)), n = 50L, times = c(1, 2))
  ns <- admixr2:::.admNormaliseStudy(s, "s2")
  expect_equal(ns$method, "var")
  expect_equal(ns$v_diag, c(0.1, 0.2))
})

test_that("Full matrix auto-detected as method = 'cov'", {
  V  <- matrix(c(0.1, 0.02, 0.02, 0.2), 2, 2)
  s  <- list(E = c(1.0, 2.0), V = V, n = 50L, times = c(1, 2))
  ns <- admixr2:::.admNormaliseStudy(s, "s3")
  expect_equal(ns$method, "cov")
  expect_null(ns$v_diag)
})

test_that("v_diag only set for 'var' studies", {
  V  <- matrix(c(0.1, 0.02, 0.02, 0.2), 2, 2)
  s  <- list(E = c(1.0, 2.0), V = V, n = 50L, times = c(1, 2))
  ns <- admixr2:::.admNormaliseStudy(s, "s4")
  expect_null(ns$v_diag)
})

test_that("Vector V + method='cov' warns and coerces to 'var'", {
  s <- list(E = c(1.0, 2.0), V = c(0.1, 0.2), n = 50L, times = c(1, 2),
            method = "cov")
  expect_warning(
    ns <- admixr2:::.admNormaliseStudy(s, "s5"),
    regexp = "method='var'"
  )
  expect_equal(ns$method, "var")
})

test_that("Non-diagonal V + method='var' warns about off-diagonal entries", {
  V <- matrix(c(0.1, 0.02, 0.02, 0.2), 2, 2)
  s <- list(E = c(1.0, 2.0), V = V, n = 50L, times = c(1, 2), method = "var")
  expect_warning(
    admixr2:::.admNormaliseStudy(s, "s6"),
    regexp = "off-diagonal"
  )
})

test_that("Missing E stops with informative message", {
  s <- list(V = diag(2), n = 50L, times = c(1, 2))
  expect_error(admixr2:::.admNormaliseStudy(s, "study_x"), regexp = "missing 'E'")
})

test_that("Missing V stops with informative message", {
  s <- list(E = c(1, 2), n = 50L, times = c(1, 2))
  expect_error(admixr2:::.admNormaliseStudy(s, "study_x"), regexp = "missing 'V'")
})

test_that("Missing n stops with informative message", {
  s <- list(E = c(1, 2), V = diag(2), times = c(1, 2))
  expect_error(admixr2:::.admNormaliseStudy(s, "study_x"), regexp = "missing 'n'")
})

test_that("Missing times stops with informative message", {
  s <- list(E = c(1, 2), V = diag(2), n = 50L)
  expect_error(admixr2:::.admNormaliseStudy(s, "study_x"), regexp = "missing 'times'")
})

test_that("Explicit method='cov' on full matrix is respected", {
  V  <- matrix(c(0.1, 0.02, 0.02, 0.2), 2, 2)
  s  <- list(E = c(1.0, 2.0), V = V, n = 50L, times = c(1, 2), method = "cov")
  ns <- admixr2:::.admNormaliseStudy(s, "s7")
  expect_equal(ns$method, "cov")
})

test_that("Study name appears in error message", {
  s <- list(V = diag(2), n = 50L, times = c(1, 2))
  expect_error(admixr2:::.admNormaliseStudy(s, "my_study"), regexp = "my_study")
})

