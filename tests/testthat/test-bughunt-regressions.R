test_that("covDist rejects empty and duplicate covariates", {
  expect_error(covDist(data.frame(covariate = character(), mean = numeric(),
                                  sd = numeric())), "at least one")
  expect_error(covDist(data.frame(covariate = NA_character_, mean = 1, sd = 1)),
               "non-empty name")
  expect_error(covDist(A = c(mean = 0, sd = 1),
                       A = c(mean = 10, sd = 1)), "duplicate.*A")
})

test_that("transcription slips in a covariate table are caught", {
  # a range is not two proportions: this used to become a binary covariate
  expect_error(covDist(WT = c(min = 50, max = 110)), "looks like a range")
  expect_error(covDist(WT = c(Lower = 50, Upper = 110)), "looks like a range")
  # counts and percentages are relative weights, and stay accepted
  expect_equal(covDist(SEX = c(male = 55, female = 45))$SEX$probs, c(.55, .45))
  expect_silent(covDist(G = c(a = 1, b = 1)))
  expect_silent(covDist(DOSE = c(low = .4, high = .6)))
  # rounding is fine, a missing level is not
  expect_silent(covDist(G = c(a = .33, b = .33, c = .33)))
  expect_warning(covDist(G = c(a = .6, b = .3)), "sum to 0.9")
  # a cv is a percent
  expect_warning(admPopulation(WT = c(mean = 80, cv = 0.25)), "percent")
  expect_silent(admPopulation(WT = c(mean = 80, cv = 25)))
  expect_error(admPopulation(C = c(median = 200, iqr = c(62, 118))),
               "outside its own iqr")
})

test_that("lognormal moment matching remains finite at extreme scales", {
  x <- covDist(A = c(mean = 1e200, sd = 1e199), dist = "lnorm")$A
  expect_true(all(is.finite(unlist(x[c("meanlog", "sdlog")]))))
  expect_equal(x$meanlog, log(1e200) - log1p(0.1^2) / 2)
  expect_equal(x$sdlog, sqrt(log1p(0.1^2)))
})

test_that("study normalization rejects invalid size, variance, and endpoint", {
  study <- list(n = 10, times = 1, E = 2, V = 1)
  expect_error(admixr2:::.admNormaliseStudy(within(study, n <- 0), "s"),
               "finite positive")
  expect_error(admixr2:::.admNormaliseStudy(within(study, V <- -1), "s"),
               "non-negative")
  expect_error(admixr2:::.admNormaliseStudy(
    list(n = 10, times = 1:2, E = 1:2, V = matrix(c(1, 2, 2, 1), 2)), "s"),
    "positive-semidefinite")
  expect_error(admixr2:::.admNormaliseStudy(within(study, V <- NA_real_), "s"),
               "finite")
  expect_error(admixr2:::.admNormaliseStudy(
    list(n = 10, times = 1:2, E = 1:2, V = matrix(c(1, NA, NA, 1), 2)), "s"),
    "finite")
  expect_error(admixr2:::.admExpandLongStudy(
    list(n = 10, data = data.frame(DVID = factor(c("cp", NA)),
                                   TIME = c(1, 2), E = c(1, 99), V = c(1, 1))),
    "s"), "endpoint values")
})

test_that("joint covariance follows sorted observation times", {
  study <- list(
    n = 10, ev = "EV", V = diag(c(.4, .1, .2, .8)),
    observations = list(
      plasma = list(times = c(4, 1, 2), E = c(30, 10, 20)),
      csf = list(times = 8, E = 1)
    )
  )
  unit <- admixr2:::.admBuildJointUnit(study, "s", NULL)
  expect_equal(unit$E, c(10, 20, 30, 1))
  expect_equal(diag(unit$V), c(.1, .2, .4, .8))
  study$n <- 0
  expect_error(admixr2:::.admBuildJointUnit(study, "s", NULL), "finite positive")
})

test_that("fractional stratum sizes are summed before BIC", {
  stats <- admixr2:::.admCalcObjStats(
    100, 2,
    list(a = list(n = .6, times = 1:2), b = list(n = .4, times = 1:2))
  )
  expect_equal(stats$nobs, 2)
  expect_true(is.finite(stats$objDf$BIC))
})
