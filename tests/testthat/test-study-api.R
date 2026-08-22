# What a user has is a PAPER -- a parameter table with %RSE, a baseline
# demographics table and a design -- not an nlmixr2 fit. These pin the
# transcription, especially the conversions that are easy to do wrong by hand.

.sa_model <- function() {
  ini({ tcl <- log(5); tv <- log(50); bsex <- 0.15
        eta.cl ~ 0.09; add.err <- 0.08 })
  model({ cl <- exp(tcl + eta.cl) * (WT/70)^0.75 * exp(bsex * SEX)
          v  <- exp(tv) * (WT/70); cp <- linCmt(); cp ~ add(add.err) })
}

test_that("admPopulation reads a baseline table in the forms papers print", {
  skip_if_not_installed("randtoolbox")
  X <- covDraw(admPopulation(WT = c(mean = 75, sd = 16)), n = 40000L)
  expect_equal(mean(X[, "WT"]), 75, tolerance = 0.01)
  expect_equal(stats::sd(X[, "WT"]), 16, tolerance = 0.02)
  # a cv is a percent, as printed
  X <- covDraw(admPopulation(WT = c(mean = 75, cv = 20)), n = 40000L)
  expect_equal(stats::sd(X[, "WT"]), 15, tolerance = 0.02)
  # median + IQR, consistent with the lognormal, reproduces the quartiles
  X <- covDraw(admPopulation(CRCL = c(median = 92, iqr = c(70, 121))), n = 40000L)
  expect_equal(unname(stats::quantile(X[, "CRCL"], c(.25, .5, .75))),
               c(70, 92, 121), tolerance = 0.01)
  # a binary covariate is a proportion
  X <- covDraw(admPopulation(SEX = c(male = 0.55)), n = 20000L)
  expect_equal(mean(X[, "SEX"] == 1), 0.55, tolerance = 0.02)
})

test_that("a median and IQR inconsistent with the shape is SAID, not swallowed", {
  # Real quartiles are rarely symmetric about the median on any scale, so only
  # two of the three numbers can be honoured. Reproducing neither quartile
  # silently is the failure this avoids.
  expect_message(admPopulation(CRCL = c(median = 92, iqr = c(62, 118))),
                 "not consistent")
  expect_silent(admPopulation(CRCL = c(median = 92, iqr = c(70, 121))))
  # a min-max is a rule of thumb, and says so
  expect_message(admPopulation(WT = c(median = 75, range = c(45, 110))),
                 "rule of thumb")
})

test_that("cor names PAIRS, so an independent covariate needs no padding", {
  skip_if_not_installed("randtoolbox")
  p <- admPopulation(WT = c(mean = 75, sd = 16), CRCL = c(mean = 92, sd = 24),
                     SEX = c(male = 0.55), cor = c(WT.CRCL = 0.45))
  X <- covDraw(p, n = 40000L)
  expect_equal(stats::cor(log(X[, "WT"]), log(X[, "CRCL"])), 0.45,
               tolerance = 0.02)
  expect_equal(stats::cor(X[, "WT"], X[, "SEX"]), 0, tolerance = 0.02)
  expect_error(admPopulation(WT = c(mean = 75, sd = 16),
                             CRCL = c(mean = 92, sd = 24), cor = c(WT = 0.4)),
               "does not name two")
  # a discrete covariate cannot be latently correlated -- a level would be a
  # truncation of the latent normal rather than a point
  expect_error(admPopulation(WT = c(mean = 75, sd = 16), SEX = c(male = 0.55),
                             cor = c(WT.SEX = 0.3)), "DISCRETE")
})

test_that("%RSE converts on the scale each parameter actually lives on", {
  skip_if_not_installed("rxode2")
  # THE conversion to get wrong. A paper reporting CL = 5.2 with 4.1% RSE means
  # SE(CL) = 0.041 * 5.2, so SE(log CL) = 0.041 -- the estimate does NOT appear.
  # An ordinary coefficient takes the usual |estimate| * RSE/100. Doing both the
  # same way is a factor of log(5.2) = 1.65 on the clearance, silently.
  s <- admStudy(model = .sa_model,
                est = c(tcl = log(5.2), bsex = 0.17, eta.cl = 0.11,
                        add.err = 0.09),
                rse = c(tcl = 4.1, bsex = 31, eta.cl = 18, add.err = 9),
                n = 240, dose = 200, times = c(1, 4, 12), label = "s")
  expect_equal(s$se[["tcl"]], 0.041, tolerance = 1e-9)          # NOT 0.041*1.649
  expect_equal(s$se[["bsex"]], 0.31 * 0.17, tolerance = 1e-9)
  expect_equal(s$se[["eta.cl"]], 0.18 * 0.11, tolerance = 1e-9)
  expect_equal(s$se[["add.err"]], 0.09 * 0.09, tolerance = 1e-9)
  # and it lands in a covariance keyed by the ini() names
  expect_setequal(rownames(s$cov), c("tcl", "bsex", "eta.cl", "add.err"))
  expect_equal(sqrt(s$cov["tcl", "tcl"]), 0.041, tolerance = 1e-9)
  # a relative error on a zero estimate has no meaning
  expect_error(admStudy(model = .sa_model, est = c(bsex = 0),
                        rse = c(bsex = 20), n = 10, dose = 1, times = 1),
               "RELATIVE")
})

test_that("the paper's estimates go INTO the model, not into globals", {
  skip_if_not_installed("rxode2")
  s <- admStudy(model = .sa_model, est = c(tcl = log(7.5), bsex = 0.22),
                n = 100, dose = 200, times = c(1, 4), label = "s")
  expect_equal(s$ui$iniDf$est[s$ui$iniDf$name == "tcl"], log(7.5),
               tolerance = 1e-12)
  expect_equal(s$ui$iniDf$est[s$ui$iniDf$name == "bsex"], 0.22)
  # untouched parameters keep the model's own value
  expect_equal(s$ui$iniDf$est[s$ui$iniDf$name == "tv"], log(50))
  expect_error(admStudy(model = .sa_model, est = c(nope = 1), n = 10,
                        dose = 1, times = 1), "does not declare")
})

test_that("a study is one currency or the other, and says which", {
  skip_if_not_installed("rxode2")
  expect_error(admStudy(model = .sa_model, E = 1:3, n = 10, dose = 1,
                        times = 1:3), "BOTH")
  expect_error(admStudy(n = 10, dose = 1, times = 1:3), "either a .model.")
  expect_error(admStudy(E = 1:3, n = 10, dose = 1, times = 1:3), "spread")
  expect_error(admStudy(E = 1:3, sd = 1:3, dose = 1, times = 1:3), "positive .n.")
  expect_error(admStudy(E = 1:2, sd = 1:2, n = 5, dose = 1, times = 1:3),
               "2 values but .times. has 3")
})

test_that("a SEM is scaled back to a per-subject spread by sqrt(n)", {
  # V is a per-SUBJECT covariance. Reading a standard error of the mean as an SD
  # understates the spread by sqrt(n), which is the commonest error in
  # digitising a published figure.
  s <- admStudy(E = c(9, 7, 5), sem = c(0.2, 0.15, 0.1), n = 100, dose = 200,
                times = c(1, 4, 12), label = "s")
  expect_equal(sqrt(s$V), c(0.2, 0.15, 0.1) * 10, tolerance = 1e-12)
  expect_error(admStudy(E = c(9, 7), sd = c(1, 1), sem = c(1, 1), n = 10,
                        dose = 1, times = c(1, 2)), "both")
})

test_that("admStudies names studies from the objects they were built as", {
  skip_if_not_installed("rxode2")
  a <- admStudy(E = 1:3, sd = 1:3, n = 10, dose = 1, times = 1:3)
  b <- admStudy(E = 1:3, sd = 1:3, n = 20, dose = 1, times = 1:3)
  expect_named(admStudies(a, b), c("a", "b"))
  expect_named(admStudies(first = a, b), c("first", "b"))
  expect_error(admStudies(a, list(1)), "not one")
  expect_error(admStudies(x = a, x = b), "unique")
})

test_that("materialising a spec is deferred, and matches datagen directly", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  # Lazy on purpose: building a study solves nothing. What it eventually
  # produces must be exactly what the equivalent datagen() call produces.
  s <- admStudy(model = .sa_model, est = c(tcl = log(5.2)),
                n = 240, dose = 200, times = c(1, 4, 12),
                population = admPopulation(WT = c(mean = 75, sd = 16),
                                           SEX = c(male = 0.55)),
                label = "s")
  got <- admixr2:::.admMaterialise(admStudies(s))
  ui <- suppressMessages(rxode2::rxode2(.sa_model))
  d <- ui$iniDf; d$est[d$name == "tcl"] <- log(5.2); ui$iniDf <- d
  want <- suppressWarnings(suppressMessages(datagen(
    list(s = list(times = c(1, 4, 12), ev = rxode2::et(amt = 200), n = 240,
                  cov_dist = admPopulation(WT = c(mean = 75, sd = 16),
                                           SEX = c(male = 0.55)))),
    model = ui, control = datagenControl(method = "gh"))))
  expect_equal(got$s$E, want$s$E, tolerance = 1e-10)
  expect_equal(got$s$V, want$s$V, tolerance = 1e-10)
})

test_that("`by` expands into one study per level, splitting n by the proportion", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  # A paper reporting results separately by sex really did report each
  # subgroup, so each becomes an ordinary study with sex PINNED -- not a
  # stratified one.
  s <- admStudy(model = .sa_model, n = 200, dose = 200, times = c(1, 4, 12),
                population = admPopulation(WT = c(mean = 75, sd = 16),
                                           SEX = c(male = 0.6)),
                by = "SEX", label = "trial")
  got <- admixr2:::.admMaterialise(admStudies(s))
  expect_length(got, 2L)
  expect_equal(sum(vapply(got, function(g) g$n, 0)), 200)
  expect_setequal(vapply(got, function(g) g[["cov"]][["SEX"]], 0), c(0, 1))
  # 60% male -> the SEX = 1 study carries 120 subjects
  n1 <- got[[which(vapply(got, function(g) g[["cov"]][["SEX"]], 0) == 1)]]$n
  expect_equal(n1, 120)
  # and SEX is no longer marginalised inside either study
  expect_false("SEX" %in% admixr2:::.admCovSpecNames(got[[1L]]$cov_dist))
})
