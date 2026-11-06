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

test_that("a plain study is put in time order, E and V with it", {
  V <- matrix(c(4, 1, 0, 1, 5, 2, 0, 2, 6), 3L)
  u <- admixr2:::.admNormaliseObs(
    list(E = c(30, 10, 20), V = V, n = 10, times = c(4, 1, 2)), "s")
  expect_equal(u$times, c(1, 2, 4))
  expect_equal(u$E, c(10, 20, 30))
  expect_equal(u$V, V[c(2, 3, 1), c(2, 3, 1)])
  s <- admixr2:::.admNormaliseStudy(
    list(E = c(2, 1), V = c(.2, .1), n = 10, times = c(2, 1)), "s")
  expect_equal(s$times, 1:2)
  expect_equal(s$v_diag, c(.1, .2))
})

test_that("a missing mean or time is refused rather than fitted", {
  ok <- list(E = c(1, 2), V = c(1, 1), n = 10, times = c(1, 2))
  expect_error(admixr2:::.admNormaliseObs(utils::modifyList(ok, list(E = c(1, NA))), "s"),
               "E must be finite")
  expect_error(admixr2:::.admNormaliseObs(utils::modifyList(ok, list(times = c(1, NA))), "s"),
               "must be finite numbers")
  expect_error(admixr2:::.admBuildJointUnit(
    list(n = 10, ev = "EV", V = diag(2),
         observations = list(a = list(times = 1, E = NA), b = list(times = 1, E = 1))),
    "s", NULL), "must be finite")
})

test_that("a lognormal margin keeps the centre the paper reported", {
  sp <- function(...) admixr2:::.admPopSpec(c(...), "WT", "lnorm")
  x <- sp(median = 92, sd = 30)
  expect_equal(exp(x$meanlog), 92)
  t2 <- exp(x$sdlog^2)
  expect_equal(92 * sqrt(t2 * (t2 - 1)), 30)
  x <- sp(median = 92, cv = 30)
  expect_equal(exp(x$meanlog), 92)
  expect_equal(sqrt(exp(x$sdlog^2) - 1), .3)
  expect_equal(exp(suppressMessages(sp(median = 92, range = c(40, 180)))$meanlog), 92)
  x <- suppressMessages(sp(mean = 92, iqr = c(70, 115)))
  expect_equal(exp(x$meanlog + x$sdlog^2 / 2), 92)
  x <- sp(mean = 92, sd = 30)
  expect_equal(exp(x$meanlog + x$sdlog^2 / 2), 92)
})

test_that("a margin that describes no distribution is refused where it is written", {
  expect_error(admPopulation(WT = c(mean = 72, sd = -16)), "POSITIVE")
  expect_error(admPopulation(WT = c(mean = 72, sd = 0)), "POSITIVE")
  expect_error(admPopulation(WT = c(mean = -72, sd = 16)), "lognormal")
  expect_error(admPopulation(WT = c(mean = NA, sd = 16)), "non-finite")
  expect_error(admPopulation(WT = c(meanlog = 4, sdlog = -.3)), "POSITIVE")
  expect_error(admPopulation(WT = c(median = 70, iqr = c(70, 70))), "POSITIVE")
  expect_error(admPopulation(WT = c(median = 5, iqr = c(0, 9))), "iqr")
  expect_error(covDist(WT = c(mu = 70, sd = -5)), "POSITIVE")
  expect_error(covDist(WT = c(mean = NA_real_, sd = 5)), "POSITIVE")
  expect_error(admPopulation(G = list(values = 0:2, probs = c(.5, .5))), "one non-negative")
  expect_error(admPopulation(G = list(values = 0:1, probs = c(-.5, 1.5))), "one non-negative")
  z <- admPopulation(Z = c(median = 0, iqr = c(-1, 1)), dist = "normal")
  expect_equal(z$Z$mu, 0)
})

test_that("a named correlation matrix is matched by name or refused", {
  R <- diag(3); dimnames(R) <- rep(list(c("WT", "CRCL", "AGE")), 2L)
  R["WT", "AGE"] <- R["AGE", "WT"] <- .7
  expect_error(covDist(WT = c(mean = 72, sd = 16), AGEE = c(mean = 40, sd = 9),
                       CRCL = c(mean = 90, sd = 20), cor = R), "no row for")
  ok <- covDist(WT = c(mean = 72, sd = 16), CRCL = c(mean = 90, sd = 20),
                AGE = c(mean = 40, sd = 9), cor = R[c(3, 1, 2), c(3, 1, 2)])
  expect_equal(ok$latentR[1L, 3L], .7)
})

test_that("an event table must be et()'s and dose a compartment the model has", {
  ui <- list(mv0 = list(state = c("depot", "central")))
  chk <- function(ev, ui_ = ui)
    admixr2:::.admCheckEv(list(list(label = "s", ev = ev)), ui_)
  et <- rxode2::et
  expect_error(chk(et(amt = 100, cmt = "nope")), "does not have")
  expect_error(chk(et(amt = 100, cmt = 7L)), "does not have")
  expect_silent(chk(et(amt = 100, cmt = "central")))
  expect_silent(chk(et(amt = 100, cmt = 2L)))
  expect_silent(chk(et(amt = 100)))
  expect_silent(chk(et(amt = 100, cmt = "central") |> et(1, cmt = "cp")))
  expect_silent(chk(NULL))
  expect_silent(chk(et(amt = 100, cmt = "nope"), list()))
  expect_error(chk(data.frame(TIME = 0, AMT = 100, EVID = 1)), "rxode2::et")
})

test_that("summary() of a fit is the fit", {
  fit <- structure(list(), class = "admFit")
  expect_identical(summary(fit), fit)
})

test_that("slips the review of the table checks turned up", {
  # a median is checked against its IQR even when a mean is the centre used
  x <- suppressMessages(admPopulation(X = c(mean = 140, median = 90, iqr = c(62, 118))))
  expect_gt(x$X$sdlog, 0)
  expect_error(admPopulation(X = c(mean = 90, median = 140, iqr = c(62, 118))),
               "median 140")
  # lo/hi are level labels unless the values cannot be proportions
  expect_equal(covDist(DOSE = c(lo = .4, hi = .6))$DOSE$probs, c(.4, .6))
  expect_error(covDist(WT = c(lo = 50, hi = 110)), "looks like a range")
  # a percentage table with a level missing
  expect_warning(covDist(SEX = c(male = 55, female = 40)), "sum to 95")
  expect_silent(covDist(SEX = c(male = 55, female = 45)))
  expect_silent(covDist(G = c(a = 1, b = 2, c = 1)))
  # half-named matrices are matched by name too
  n3 <- c(WT = 72, CRCL = 90, AGE = 40)
  R <- diag(3); R[1L, 3L] <- R[3L, 1L] <- .7
  mk <- function(M) covDist(AGE = c(mean = 40, sd = 9), WT = c(mean = 72, sd = 16),
                            CRCL = c(mean = 90, sd = 20), cor = M)
  dimnames(R) <- list(NULL, names(n3))
  expect_equal(mk(R)$latentR[1L, 2L], .7)
  dimnames(R) <- list(names(n3), NULL)
  expect_equal(mk(R)$latentR[1L, 2L], .7)
  dimnames(R) <- list(NULL, c("WT", "CRCL", "AGEE"))
  expect_error(mk(R), "no row for")
  # one warning per covariate, however many times its margin is resolved
  df <- data.frame(WT = c(60, 70, 80, 90), AGE = c(30, 40, 35, 50))
  n <- 0L
  withCallingHandlers(admPopulation(data = df, WT = c(mean = 80, cv = .25)),
    warning = function(w) { n <<- n + 1L; invokeRestart("muffleWarning") })
  expect_equal(n, 1L)
})

test_that("admStudy() refuses what it would otherwise drop or misread", {
  expect_error(admStudy(E = 1:2, sd = 1:2, n = 1:2, times = 1:2, dose = 100), "positive `n`")
  expect_error(admStudy(E = 1:2, sd = 1:2, n = 20, times = 1:2, dose = -5), "`dose`")
  expect_error(admStudy(E = 1:2, sd = 1:2, n = 20, times = 1:2, dose = c(1, 2)), "`dose`")
  expect_error(admStudy(E = 1:2, sd = 1:2, n = 20, times = 1:2, dose = "100"), "`dose`")
  expect_error(admStudy(E = 1:2, sd = 1:2, n = 20, times = c(1, NA), dose = 100), "`times`")
  expect_s3_class(admStudy(E = 1:2, sd = 1:2, n = 20, times = 1:2, dose = 0), "admStudy")
  expect_error(admPopulation(data = data.frame(WT = 70)), "two subjects")
})
