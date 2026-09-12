# anova() on nested admFits: the ORDINARY likelihood-ratio test.
#
# Everything here is about the arithmetic and the refusals. There is no
# reweighting to test -- that is the point of the plain form -- so what is left
# is that the test statistic is the objective difference, that Df counts the
# added parameters, and that the comparisons which are NOT likelihood ratios
# are refused rather than reported.

# A minimal stand-in for an admFit: anova() reads only `objective`, `env$method`,
# `env$nNodes`, `env$AIC`/`env$BIC` and the parameter names off admExtra.
.lrt_fit <- function(par_names, objective, method = "adgh", nNodes = 5L,
                     strataNodes = NULL,
                     AIC = NA_real_, BIC = NA_real_) {
  e <- new.env(parent = emptyenv())
  e$admExtra <- list(par_names = par_names)
  e$method   <- method
  e$nNodes   <- nNodes
  e$strataNodes <- strataNodes
  e$AIC      <- AIC
  e$BIC      <- BIC
  structure(list(env = e, objective = objective), class = "admFit")
}

test_that("the statistic is the objective difference on the added parameters", {
  full <- .lrt_fit(c("tcl", "tv", "b1"), 100)
  red  <- .lrt_fit(c("tcl", "tv"),       106)
  a <- anova(full, red)
  expect_s3_class(a, "anova.admFit")
  # smallest model first, the way anova.lm orders its table
  expect_identical(a$Npar, c(2L, 3L))
  expect_equal(a$dOFV[2L], 6)
  expect_identical(a$Df[2L], 1L)
  expect_equal(a$p[2L], stats::pchisq(6, 1, lower.tail = FALSE))
  # the first row is the reference model and carries no test
  expect_true(is.na(a$dOFV[1L]))
  expect_true(is.na(a$p[1L]))
})

test_that("Df counts EVERY added parameter, not just one", {
  full <- .lrt_fit(c("tcl", "tv", "b1", "b2", "b3"), 40)
  red  <- .lrt_fit(c("tcl", "tv"),                   55)
  a <- anova(full, red)
  expect_identical(a$Df[2L], 3L)
  expect_equal(a$p[2L], stats::pchisq(15, 3, lower.tail = FALSE))
})

test_that("a NEGATIVE dOFV is reported, not clamped to zero", {
  # The larger model cannot fit worse at its own optimum, so this says one of
  # the two did not converge. Rounding it to zero would hide that.
  full <- .lrt_fit(c("tcl", "tv", "b1"), 110)
  red  <- .lrt_fit(c("tcl", "tv"),       106)
  a <- anova(full, red)
  expect_lt(a$dOFV[2L], 0)
  expect_true(is.na(a$p[2L]))
})

test_that("a non-nested pair is refused, not approximated around", {
  a <- .lrt_fit(c("tcl", "tv", "b1"), 100)
  b <- .lrt_fit(c("tcl", "tv", "b2"), 101)
  expect_error(anova(a, b), "not nested")
  # and so is a pair with nothing to test
  c1 <- .lrt_fit(c("tcl", "tv"), 100)
  c2 <- .lrt_fit(c("tcl", "tv"), 101)
  expect_error(anova(c1, c2), "nothing to test")
})

test_that("fits from DIFFERENT estimators are refused", {
  # Each scores its own approximation to the same likelihood -- FO-linearised,
  # quadrature, Monte Carlo -- so the difference of two objectives is not a
  # likelihood ratio. This returned a perfectly finite p before it was caught.
  full <- .lrt_fit(c("tcl", "tv", "b1"), 100, method = "adgh")
  red  <- .lrt_fit(c("tcl", "tv"),       106, method = "adfo")
  expect_error(anova(full, red), "different estimators")
})

test_that("fits on different node counts are refused", {
  # The objective moves with the grid, which is why .adghGrid refuses to change
  # the point count mid-fit in the first place.
  full <- .lrt_fit(c("tcl", "tv", "b1"), 100, nNodes = 9L)
  red  <- .lrt_fit(c("tcl", "tv"),       106, nNodes = 5L)
  expect_error(anova(full, red), "node counts")
})

test_that("fits on different strata grids are refused", {
  full <- .lrt_fit(c("tcl", "tv", "b1"), 100, strataNodes = 9L)
  red  <- .lrt_fit(c("tcl", "tv"),       106, strataNodes = 5L)
  expect_error(anova(full, red), "stratum node counts")
})

test_that("anova() needs a pair, and needs admFits", {
  full <- .lrt_fit(c("tcl", "tv", "b1"), 100)
  expect_error(anova(full), "nothing to compare")
  expect_error(anova(full, 1), "must be an admFit")
})

test_that("printing keeps `Test` as written", {
  # stats::print.anova routes through printCoefmat, which calls data.matrix()
  # -- and that turns a CHARACTER column into its factor CODES, so "1 vs 2"
  # would print as 1.
  full <- .lrt_fit(c("tcl", "tv", "b1"), 100)
  red  <- .lrt_fit(c("tcl", "tv"),       106)
  out <- utils::capture.output(print(anova(full, red)))
  expect_true(any(grepl("Likelihood-ratio test", out)))
  expect_true(any(grepl("1 vs 2", out, fixed = TRUE)))
})
