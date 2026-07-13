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

# -- long-format studies (data frame keyed by DVID / CMT) ----------------------

# One row per observed endpoint/time: the nlmixr2 way of keying observations.
.long_df <- function()
  data.frame(DVID = c("cp", "cp", "cp", "cb", "cb"),
             TIME = c(0.5, 1, 2, 1, 2),
             E    = c(8.8, 7.8, 6.9, 3.0, 3.4),
             V    = c(1.21, 0.81, 0.64, 0.09, 0.09))

test_that("Long format without a joint V gives one independent unit per endpoint", {
  ns <- admixr2:::.admNormaliseStudy(
    list(n = 60L, ev = "EV", data = .long_df()), "lit")

  expect_named(ns$observations, c("cp", "cb"))
  expect_false(isTRUE(ns$joint))
  expect_equal(ns$observations$cp$times, c(0.5, 1, 2))
  expect_equal(ns$observations$cb$E, c(3.0, 3.4))
  expect_equal(ns$observations$cb$output, "cb")
  # variance column -> diagonal V, method "var", n/ev inherited from the study
  expect_equal(ns$observations$cp$method, "var")
  expect_equal(ns$observations$cp$v_diag, c(1.21, 0.81, 0.64))
  expect_equal(ns$observations$cb$n, 60L)
  expect_equal(ns$observations$cb$ev, "EV")
})

test_that("Long format matches the equivalent observations spec", {
  long <- admixr2:::.admNormaliseStudy(
    list(n = 60L, ev = "EV", data = .long_df()), "lit")
  obs <- admixr2:::.admNormaliseStudy(list(n = 60L, ev = "EV", observations = list(
    cp = list(output = "cp", times = c(0.5, 1, 2), E = c(8.8, 7.8, 6.9),
              V = c(1.21, 0.81, 0.64)),
    cb = list(output = "cb", times = c(1, 2), E = c(3.0, 3.4), V = c(0.09, 0.09)))),
    "lit")

  for (k in c("cp", "cb"))
    expect_equal(long$observations[[k]][c("output", "times", "E", "V", "n", "method", "v_diag")],
                 obs$observations[[k]][c("output", "times", "E", "V", "n", "method", "v_diag")])
})

test_that("Long format + study-level V builds one joint (same-subject) unit", {
  V  <- diag(c(1.21, 0.81, 0.64, 0.09, 0.09))
  V[1, 4] <- V[4, 1] <- 0.15                    # plasma(t=0.5) <-> brain(t=1)
  ns <- admixr2:::.admNormaliseStudy(
    list(n = 60L, ev = "EV", data = .long_df(), V = V), "lit")

  expect_length(ns$observations, 1L)
  u <- ns$observations[[1L]]
  expect_true(isTRUE(u$is_joint))
  expect_equal(u$E, c(8.8, 7.8, 6.9, 3.0, 3.4))   # stacked, blocks in data order
  expect_equal(u$V, unname(V))                    # rows aligned to the data rows
  expect_equal(u$row_output, c(1L, 1L, 1L, 2L, 2L))
  expect_equal(vapply(u$blocks, `[[`, character(1), "output"), c("cp", "cb"))
  expect_equal(u$n, 60L)
})

test_that("Long format keeps E and V aligned when data rows are unordered", {
  d    <- .long_df()
  V    <- diag(c(1.21, 0.81, 0.64, 0.09, 0.09)); V[1, 4] <- V[4, 1] <- 0.15
  ord  <- c(4, 1, 5, 3, 2)                        # scramble the rows
  ns   <- admixr2:::.admNormaliseStudy(
    list(n = 60L, ev = "EV", data = d[ord, ], V = V[ord, ord]), "lit")
  u    <- ns$observations[[1L]]

  # blocks follow first appearance (cb now leads); E and V permute together
  stack <- c(4, 5, 1, 2, 3)
  expect_equal(vapply(u$blocks, `[[`, character(1), "output"), c("cb", "cp"))
  expect_equal(u$E, d$E[stack])
  expect_equal(u$V, unname(V[stack, stack]))
})

test_that("Long format accepts an SD column and per-endpoint n / ev", {
  d  <- data.frame(CMT = c("cp", "cp", "cb"), time = c(1, 2, 1),
                   mean = c(8.8, 7.8, 3.0), SD = c(1.1, 0.9, 0.3),
                   n = c(60, 60, 12))
  ns <- admixr2:::.admNormaliseStudy(
    list(ev = list(cp = "EVp", cb = "EVb"), data = d), "lit")

  expect_equal(ns$observations$cp$v_diag, c(1.21, 0.81))   # SD^2
  expect_equal(ns$observations$cp$n, 60)
  expect_equal(ns$observations$cb$n, 12)                   # separate experiment
  expect_equal(ns$observations$cb$ev, "EVb")
})

test_that("Long-format input errors are informative", {
  d <- .long_df()
  expect_error(admixr2:::.admNormaliseStudy(
    list(n = 1L, ev = "EV", data = d[, c("DVID", "TIME")]), "lit"),
    regexp = "mean column")
  expect_error(admixr2:::.admNormaliseStudy(
    list(n = 1L, ev = "EV", data = d[, c("DVID", "TIME", "E")]), "lit"),
    regexp = "variance column")
  expect_error(admixr2:::.admNormaliseStudy(
    list(n = 1L, ev = "EV", data = rbind(d, d[1, ])), "lit"),
    regexp = "duplicate endpoint/time")
  expect_error(admixr2:::.admNormaliseStudy(
    list(n = 1L, ev = "EV", data = d, V = diag(4)), "lit"),
    regexp = "must be 5 x 5")
  # a joint study is one experiment: it cannot carry per-endpoint n or ev
  expect_error(admixr2:::.admNormaliseStudy(
    list(ev = "EV", data = transform(d, n = ifelse(DVID == "cp", 60, 12)),
         V = diag(5)), "lit"),
    regexp = "one shared `n`")
  expect_error(admixr2:::.admNormaliseStudy(
    list(n = 1L, ev = list(cp = "a", cb = "b"), data = d, V = diag(5)), "lit"),
    regexp = "shares one `ev`")
})

