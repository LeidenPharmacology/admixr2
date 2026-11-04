# Published-study transcription and conversion checks.

# `.sa_model` asserts WT's exponent and estimates `bsex`; `.sa_fitwt` ESTIMATES
# the exponent, so WT is conditional for it. The pair is what separates "reads a
# covariate" from "estimated a coefficient for one", which is the rule
# conditioning is derived by, so both belong at file scope rather than being
# retyped inside the tests that need the contrast.
.sa_fitwt <- function() {
  ini({ tcl <- log(5); tv <- log(50); bwt <- 0.75
        eta.cl ~ 0.09; add.err <- 0.08 })
  model({ cl <- exp(tcl + eta.cl) * (WT/70)^bwt
          v <- exp(tv); cp <- linCmt(); cp ~ add(add.err) })
}

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
  # Inconsistent median/IQR inputs must be reported.
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
  expect_error(admPopulation(WT = c(mean = 75, sd = 16),
                             CRCL = c(mean = 92, sd = 24),
                             cor = c(WT.CRCL = 0.2, CRCL.WT = 0.8)),
               "same covariate pair")
  # Discrete margins cannot use the continuous latent correlation model.
  expect_error(admPopulation(WT = c(mean = 75, sd = 16), SEX = c(male = 0.55),
                             cor = c(WT.SEX = 0.3)), "DISCRETE")
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
  expect_error(admStudy(model = .sa_model, est = c(tcl = 1, tcl = 2), n = 10,
                        dose = 1, times = 1), "must be unique")
  expect_error(admStudy(model = .sa_model, est = c(tcl = Inf), n = 10,
                        dose = 1, times = 1), "finite numeric")

  ui <- suppressMessages(rxode2::rxode2(.sa_model))
  original <- ui$iniDf$est[ui$iniDf$name == "tcl"]
  a <- admStudy(model = ui, est = c(tcl = log(6)), n = 10,
                dose = 1, times = 1)
  b <- admStudy(model = ui, est = c(tcl = log(7)), n = 10,
                dose = 1, times = 1)
  expect_equal(ui$iniDf$est[ui$iniDf$name == "tcl"], original)
  expect_equal(a$ui$iniDf$est[a$ui$iniDf$name == "tcl"], log(6))
  expect_equal(b$ui$iniDf$est[b$ui$iniDf$name == "tcl"], log(7))
})

test_that("a study is one currency or the other, and says which", {
  skip_if_not_installed("rxode2")
  expect_error(admStudy(model = .sa_model, E = 1:3, n = 10, dose = 1,
                        times = 1:3), "BOTH")
  expect_error(admStudy(model = .sa_model, V = 1, n = 10, dose = 1,
                        times = 1), "BOTH")
  expect_error(admStudy(E = 1, sd = 1, est = c(tcl = 1), n = 10, dose = 1,
                        times = 1), "cannot have `est`")
  expect_error(admStudy(n = 10, dose = 1, times = 1:3), "either a .model.")
  expect_error(admStudy(E = 1:3, n = 10, dose = 1, times = 1:3), "spread")
  expect_error(admStudy(E = 1:3, sd = 1:3, dose = 1, times = 1:3), "positive .n.")
  expect_error(admStudy(E = 1:2, sd = 1:2, n = 5, dose = 1, times = 1:3),
               "2 values but .times. has 3")
})

test_that("a SEM is scaled back to a per-subject spread by sqrt(n)", {
  # V is per-subject, so SEM covariance is multiplied by n.
  s <- admStudy(E = c(9, 7, 5), sem = c(0.2, 0.15, 0.1), n = 100, dose = 200,
                times = c(1, 4, 12), label = "s")
  expect_equal(sqrt(s$V), c(0.2, 0.15, 0.1) * 10, tolerance = 1e-12)
  expect_error(admStudy(E = c(9, 7), sd = c(1, 1), sem = c(1, 1), n = 10,
                        dose = 1, times = c(1, 2)), "both")
  expect_error(admStudy(E = c(9, 7), sd = c(-1, 1), n = 10,
                        dose = 1, times = c(1, 2)), "non-negative")
  expect_error(admStudy(E = c(9, 7), sem = c(1, Inf), n = 10,
                        dose = 1, times = c(1, 2)), "finite")
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
  # Materialisation must match the equivalent datagen() call.
  s <- admStudy(model = .sa_model, est = c(tcl = log(5.2)),
                n = 240, dose = 200, times = c(1, 4, 12),
                population = admPopulation(WT = c(mean = 75, sd = 16),
                                           SEX = c(male = 0.55)),
                label = "s")
  got <- suppressWarnings(suppressMessages(
    admixr2:::.admMaterialise(admStudies(s))))
  ui <- suppressMessages(rxode2::rxode2(.sa_model))
  d <- ui$iniDf; d$est[d$name == "tcl"] <- log(5.2); ui$iniDf <- d
  # The same call written out by hand, conditioned the way the derivation does:
  # `.sa_model` estimates the SEX coefficient and asserts the weight exponent.
  # BOTH SIDES MUST NAME THE SAME STRATA. Compared on `$s` alone this passed
  # while both sides were NULL -- materialise conditional, datagen conditional too, and
  # neither had a study called `s` to compare.
  want <- suppressWarnings(suppressMessages(datagen(
    list(s = list(times = c(1, 4, 12), ev = rxode2::et(amt = 200), n = 240,
                  stratify = "SEX",
                  cov_dist = admPopulation(WT = c(mean = 75, sd = 16),
                                           SEX = c(male = 0.55)))),
    model = ui, control = datagenControl(method = "gh"))))
  expect_setequal(names(got), names(want))
  expect_setequal(names(got), c("s_s1", "s_s2"))
  for (k in names(want)) {
    expect_equal(got[[k]]$E, want[[k]]$E, tolerance = 1e-10)
    expect_equal(got[[k]]$V, want[[k]]$V, tolerance = 1e-10)
    expect_equal(got[[k]]$n, want[[k]]$n, tolerance = 1e-10)
  }
})

test_that("`by` expands into one study per level, splitting n by the proportion", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  # Reported subgroups become ordinary pinned studies.
  s <- admStudy(model = .sa_model, n = 200, dose = 200, times = c(1, 4, 12),
                population = admPopulation(WT = c(mean = 75, sd = 16),
                                           SEX = c(male = 0.6)),
                by = "SEX", label = "trial")
  got <- admixr2:::.admMaterialise(admStudies(s))
  # `.sa_model` asserts the weight exponent and estimates only the SEX
  # coefficient, which `by` has already pinned one level per study -- so there
  # is nothing left to condition on and each level is one study. What `by` is
  # responsible for is the LEVELS and the n split, and the sums below hold
  # however many strata each level is cut into.
  sx <- vapply(got, function(g) g[["cov"]][["SEX"]], 0)
  expect_setequal(unique(sx), c(0, 1))
  expect_equal(sum(vapply(got, function(g) g$n, 0)), 200)
  # 60% male -> the SEX = 1 level carries 120 subjects across its strata
  expect_equal(sum(vapply(got[sx == 1], function(g) g$n, 0)), 120)
  expect_equal(sum(vapply(got[sx == 0], function(g) g$n, 0)), 80)
  # and SEX is no longer marginalised inside any of them
  expect_false(any(vapply(got, function(g)
    "SEX" %in% admixr2:::.admCovSpecNames(g$cov_dist), logical(1))))
})

test_that("materialising a subgroup refuses a colliding study name", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  # A model reading ONLY the `by` covariate, so the levels are the whole
  # expansion and the names are exactly `trial_SEX0`/`trial_SEX1` -- a model
  # estimating a second covariate's coefficient would cut each level further
  # and give the names a stratum suffix this test is not about.
  sexonly <- function() {
    ini({ tcl <- log(5); tv <- log(50); bsex <- 0.2
          add.err <- 0.1; eta.cl ~ 0.09 })
    model({ cl <- exp(tcl + eta.cl) * exp(bsex * SEX)
            v  <- exp(tv); cp <- linCmt(); cp ~ add(add.err) })
  }
  s <- admStudy(model = sexonly, n = 200, dose = 200, times = c(1, 4, 12),
                population = admPopulation(SEX = c(male = 0.6)),
                by = "SEX", label = "trial")
  other <- list(E = 1, V = 1, n = 10, times = 1, ev = rxode2::et(amt = 1))
  expect_error(admixr2:::.admMaterialise(list(trial_SEX0 = other, trial = s)),
               "duplicate name 'trial_SEX0'")
})

test_that("a model source withdraws the standard error, explicit or not", {
  skip_if_not_installed("rxode2")
  src <- admStudy(model = .sa_model, n = 60, dose = 100, times = c(1, 4, 12),
                  population = admPopulation(WT = c(mean = 75, sd = 16),
                                             SEX = c(male = 0.5)))
  dat <- admStudy(E = c(1.6, 0.9, 0.3), sd = c(.4, .25, .1), n = 48,
                  dose = 100, times = c(1, 4, 12))
  # The model, not its uncertainty fields, marks a model source.
  expect_true(admixr2:::.admHasModelSource(admStudies(a = src, b = dat)))
  expect_false(admixr2:::.admHasModelSource(admStudies(b = dat)))

  st <- admStudies(a = src, b = dat)
  expect_equal(suppressMessages(
    admixr2:::.admResolveCovMethod("r", st, FALSE)), "none")
  # Explicit covariance methods are refused for model sources.
  expect_error(admixr2:::.admResolveCovMethod("r", st, TRUE), "model-implied")
  expect_error(admixr2:::.admResolveCovMethod("r,s", st, TRUE), "model-implied")
  # "none" is what the refusal asks for, so it passes through either way
  expect_equal(admixr2:::.admResolveCovMethod("none", st, TRUE), "none")
  expect_equal(admixr2:::.admResolveCovMethod("none", st, FALSE), "none")
  # and a fit with no model source is untouched
  expect_equal(admixr2:::.admResolveCovMethod("r", admStudies(b = dat), FALSE), "r")
  expect_equal(admixr2:::.admResolveCovMethod("r,s", admStudies(b = dat), TRUE),
               "r,s")
  # it reaches the control objects, which is where a fit reads it
  expect_equal(suppressMessages(adghControl(studies = st))$covMethod, "none")
  expect_error(adghControl(studies = st, covMethod = "r"), "model-implied")
})

test_that("model-source provenance is recursive through observations", {
  nested <- list(pub = list(observations = list(
    cp = list(E = 1, V = 1, n = 20, times = 1, .adm_src = TRUE))))
  expect_true(admixr2:::.admHasModelSource(nested))
  expect_error(admixr2:::.admResolveCovMethod("r", nested, TRUE),
               "model-implied")
})

.sa_cohort <- function(n = 6000L, rho = 0.5, p_male = 0.42, seed = 7L) {
  set.seed(seed)
  z1 <- stats::rnorm(n)
  z2 <- rho * z1 + sqrt(1 - rho^2) * stats::rnorm(n)
  data.frame(WT   = exp(log(78) + 0.20 * z1),
             CRCL = exp(log(85) + 0.35 * z2),
             SEX  = stats::rbinom(n, 1L, p_male))
}

test_that("admPopulation(data=) reproduces the hand-written table", {
  skip_if_not_installed("randtoolbox")
  coh <- .sa_cohort()
  # The copula uses latent-scale correlation.
  byhand <- covDist(
    WT   = list(meanlog = mean(log(coh$WT)),   sdlog = stats::sd(log(coh$WT))),
    CRCL = list(meanlog = mean(log(coh$CRCL)), sdlog = stats::sd(log(coh$CRCL))),
    SEX  = list(values = c(0, 1),
                probs  = c(1 - mean(coh$SEX), mean(coh$SEX))),
    cor  = matrix(c(1, stats::cor(log(coh$WT), log(coh$CRCL)), 0,
                    stats::cor(log(coh$WT), log(coh$CRCL)), 1, 0,
                    0, 0, 1), 3L, byrow = TRUE,
                  dimnames = list(c("WT", "CRCL", "SEX"),
                                  c("WT", "CRCL", "SEX"))))
  X <- covDraw(admPopulation(data = coh), n = 60000L)
  Y <- covDraw(byhand,                    n = 60000L)
  # margins match the COHORT on the natural scale, which is the contract the
  # typed-out form has too
  expect_equal(mean(X[, "WT"]),   mean(coh$WT),   tolerance = 0.01)
  expect_equal(stats::sd(X[, "WT"]), stats::sd(coh$WT), tolerance = 0.02)
  expect_equal(mean(X[, "SEX"]),  mean(coh$SEX),  tolerance = 0.02)
  expect_equal(stats::cor(log(X[, "WT"]), log(X[, "CRCL"])),
               stats::cor(log(Y[, "WT"]), log(Y[, "CRCL"])), tolerance = 0.02)
})

test_that("data-derived dependence follows overridden margin definitions", {
  set.seed(42)
  z <- stats::rnorm(5000)
  d <- data.frame(WT = 75 + 10 * z,
                  CRCL = exp(4 + .2 * (.8 * z + .6 * stats::rnorm(5000))))

  p <- admPopulation(data = d, dist = "normal",
                     CRCL = c(meanlog = mean(log(d$CRCL)),
                              sdlog = stats::sd(log(d$CRCL))))
  expect_equal(p[["latentR"]][1L, 2L], stats::cor(d$WT, log(d$CRCL)),
               tolerance = 1e-12)

  p <- admPopulation(data = d, dist = "normal",
                     WT = c(mean = mean(d$WT), sd = stats::sd(d$WT)),
                     CRCL = c(mean = mean(d$CRCL), sd = stats::sd(d$CRCL)))
  expect_equal(p[["latentR"]][1L, 2L], stats::cor(d$WT, d$CRCL),
               tolerance = 1e-12)
})

test_that("character binary columns report dependence that is dropped", {
  z <- seq(-3, 3, length.out = 500)
  d <- data.frame(WT = exp(4 + .2 * z),
                  SEX = ifelse(z > 0, "male", "female"))
  expect_message(admPopulation(data = d), "DROPPED")
})

test_that("overriding a cohort margin keeps its data-derived correlations", {
  skip_if_not_installed("randtoolbox")
  coh <- .sa_cohort()
  p <- admPopulation(data = coh,
                     WT = c(mean = mean(coh$WT), sd = stats::sd(coh$WT)))
  X <- covDraw(p, n = 40000L)
  expect_equal(stats::cor(log(X[, "WT"]), log(X[, "CRCL"])),
               stats::cor(log(coh$WT), log(coh$CRCL)), tolerance = 0.02)
})

test_that("at cannot silently compete with a population margin", {
  p <- admPopulation(SEX = c(male = 0.6))
  expect_error(admStudy(E = 1, sd = 1, n = 20, dose = 1, times = 1,
                        population = p, at = list(SEX = 1)),
               "also gives it a distribution")
})

test_that("at requires named finite scalar values", {
  args <- list(E = 1, sd = 1, n = 20, dose = 1, times = 1)
  expect_error(do.call(admStudy, c(args, list(at = list(1)))),
               "unique, non-empty covariate names")
  expect_error(do.call(admStudy, c(args, list(at = list(SEX = c(0, 1))))),
               "one finite number")
  expect_error(do.call(admStudy, c(args, list(at = list(SEX = Inf)))),
               "one finite number")
})

test_that("digitised profiles refuse source-only expansion arguments", {
  args <- list(E = 1, sd = 1, n = 20, dose = 1, times = 1,
               population = admPopulation(SEX = c(male = 0.6)))
  expect_error(do.call(admStudy, c(args, list(by = "SEX"))), "digitised data")
  # `stratify` is not an argument at all any more -- conditioning is derived from
  # the source's model, and digitised data has none to derive from.
  expect_error(do.call(admStudy, c(args, list(stratify = "SEX"))),
               "unused argument")
})

test_that("a stated margin beats the data, and a dropped association is said", {
  skip_if_not_installed("randtoolbox")
  coh <- .sa_cohort()
  p <- admPopulation(SEX = c(male = 0.60), data = coh)
  expect_equal(mean(covDraw(p, n = 40000L)[, "SEX"]), 0.60, tolerance = 0.02)
  # a stated correlation is not silently replaced by the sample one either
  p2 <- admPopulation(data = coh, cor = c(WT.CRCL = 0.05))
  X  <- covDraw(p2, n = 60000L)
  expect_equal(stats::cor(log(X[, "WT"]), log(X[, "CRCL"])), 0.05,
               tolerance = 0.02)
  # A pair is unordered, so reversing its spelling must still override the
  # data-derived WT.CRCL entry rather than leaving both entries to be applied.
  p3 <- admPopulation(data = coh, cor = c(CRCL.WT = 0.05))
  X  <- covDraw(p3, n = 60000L)
  expect_equal(stats::cor(log(X[, "WT"]), log(X[, "CRCL"])), 0.05,
               tolerance = 0.02)
  # Unsupported discrete-continuous dependence is reported when dropped.
  coh2 <- coh; coh2$SEX <- as.integer(coh$WT > stats::median(coh$WT))
  expect_message(admPopulation(data = coh2), "being DROPPED")
  expect_silent(invisible(admPopulation(data = coh)))
})

test_that("conditioning is derived from the model with nothing declared", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  # NOTHING IS SET. Whether a covariate is conditional or marginal is a
  # property of the source's own model, so admixr2 derives it -- and it has to
  # work off the parsed rxUi admStudy() stores rather than the original
  # function.
  s <- admStudy(model = .sa_model, n = 200, dose = 200, times = c(1, 4, 12),
                population = admPopulation(WT = c(mean = 75, sd = 16),
                                           SEX = c(male = 0.55)),
                label = "trial")
  got <- suppressWarnings(suppressMessages(
    admixr2:::.admMaterialise(admStudies(s))))
  # .sa_model estimates bsex and reads WT/SEX, so SEX is conditional and WT is not:
  # a fixed allometric exponent carries no fitted effect to recover.
  expect_gt(length(got), 1L)
  expect_true(all(vapply(got, function(g) !is.null(g[["cov"]][["SEX"]]),
                         logical(1))))
  expect_equal(sum(vapply(got, function(g) g$n, 0)), 200)
})

test_that("`by` names one declared discrete population margin", {
  pop <- admPopulation(SEX = c(male = 0.6), WT = c(mean = 75, sd = 16))
  args <- list(model = .sa_model, n = 20, dose = 1, times = 1,
               population = pop)
  expect_error(do.call(admStudy, c(args, list(by = TRUE))), "one non-empty")
  expect_error(do.call(admStudy, c(args, list(by = c("SEX", "WT")))),
               "one non-empty")
  expect_error(do.call(admStudy, c(args, list(by = "CRCL"))), "does not declare")
  expect_error(do.call(admStudy, c(args, list(by = "WT"))), "discrete")
})

test_that("a cohort data frame IS a population, and n comes from it", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  coh <- .sa_cohort(n = 300L)
  coh$ID <- seq_len(nrow(coh))          # a real cohort carries data columns too
  s <- suppressMessages(
    admStudy(model = .sa_model, population = coh, dose = 200,
             times = c(1, 4, 12), label = "trial"))
  expect_equal(s$n, 300)
  # Keep covariates absent from the source model for cross-study inference.
  expect_setequal(admixr2:::.admCovSpecNames(s$population),
                  c("WT", "CRCL", "SEX"))
  # ...but a reserved data column is never a covariate
  expect_false("ID" %in% admixr2:::.admCovSpecNames(s$population))
  expect_message(admStudy(model = .sa_model, population = coh, dose = 200,
                          times = c(1, 4, 12)), "'ID' is a data column")
  # an explicit n still wins
  s2 <- suppressMessages(
    admStudy(model = .sa_model, population = coh, n = 120, dose = 200,
             times = c(1, 4, 12)))
  expect_equal(s2$n, 120)
})

test_that("print.admStudies flags a covariate no source can identify", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  # Report covariates marginalised identically across all studies.
  coh <- .sa_cohort(n = 300L)
  s <- suppressMessages(
    admStudy(model = .sa_model, population = coh, dose = 200,
             times = c(1, 4, 12), label = "a"))
  out <- paste(utils::capture.output(print(admStudies(a = s))), collapse = " ")
  # `.sa_model` ESTIMATES the SEX coefficient, so SEX is conditional. WT is read at a
  # fixed exponent, and CRCL is declared by the cohort and never mentioned by
  # the model: both are marginal, for the two different reasons print()
  # distinguishes, and neither is identified by any source here.
  expect_match(out, "SEX +conditional")
  expect_match(out, "WT +marginal")
  expect_match(out, "CRCL +marginal")
  expect_match(out, "marginal in every source")
  # that a published model carries no standard error is the other silent one:
  # said at transcription time, where the user can still change the design
  bare <- suppressMessages(
    admStudy(model = .sa_model, population = coh, dose = 200,
             times = c(1, 4, 12), label = "b"))
  out2 <- paste(utils::capture.output(print(admStudies(b = bare))),
                collapse = " ")
  expect_match(out2, "No standard error is reported")
})

test_that("a covariate constant within a study is pinned, not described", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  # Reject constant cohort columns at ingestion.
  coh <- .sa_cohort(n = 200L)
  coh$CRCL <- 62                       # one renal value for the whole study
  expect_error(
    suppressMessages(admStudy(model = .sa_model, population = coh, dose = 200,
                              times = c(1, 4, 12))),
    "is CONSTANT at 62")
  expect_error(
    suppressMessages(admStudy(model = .sa_model, population = coh, dose = 200,
                              times = c(1, 4, 12))),
    "at = list(CRCL = 62)", fixed = TRUE)
  # ...and the route it names actually works
  s <- suppressMessages(
    admStudy(model = .sa_model, population = coh[, c("WT", "SEX")],
             at = list(CRCL = 62), n = 200, dose = 200, times = c(1, 4, 12)))
  expect_equal(s$at$CRCL, 62)
  expect_false("CRCL" %in% admixr2:::.admCovSpecNames(s$population))
})

test_that("the derivation conditions on what the source model ESTIMATED", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  # ESTIMATED, not merely read. `.sa_model` reads WT at a literal ^0.75 and SEX
  # at an estimated coefficient: the fixed exponent is an assumption, so that
  # source holds no evidence about weight and conditioning on it would buy strata
  # and no evidence. Marginal covers both that and a covariate the population
  # describes and the model never mentions.
  s <- suppressMessages(
    admStudy(model = .sa_model, n = 200, dose = 200, times = c(1, 4, 12),
             population = admPopulation(WT = c(mean = 75, sd = 16),
                                        SEX = c(male = 0.55)),
             label = "trial"))
  got <- suppressWarnings(suppressMessages(
    admixr2:::.admMaterialise(admStudies(s))))
  expect_setequal(unlist(lapply(got, function(g) names(g[["cov"]]))), "SEX")
  # WT is still described, and integrated over.
  expect_true(all(vapply(got, function(g)
    !is.null(g$cov_dist[["WT"]]) && !isTRUE(g$cov_dist[["WT"]][[".point"]]),
    logical(1))))
  expect_equal(sum(vapply(got, function(g) g$n, 0)), 200)
})

test_that("admStudy resolves v_denom from the currency the study is written in", {
  ev <- rxode2::et(amt = 100)
  tt <- c(1, 2, 4); EE <- c(2, 1.5, 1)

  # Published sd/sem values use the n-1 denominator.
  s_sd <- admStudy(E = EE, sd = c(.4, .3, .2), n = 60L, times = tt, ev = ev)
  expect_identical(s_sd$v_denom, "unbiased")
  s_sem <- admStudy(E = EE, sem = c(.4, .3, .2) / sqrt(60), n = 60L,
                    times = tt, ev = ev)
  expect_identical(s_sem$v_denom, "unbiased")

  # Supplied covariance matrices default to ML.
  s_V <- admStudy(E = EE, V = diag(c(.16, .09, .04)), n = 60L, times = tt, ev = ev)
  expect_identical(s_V$v_denom, "ml")

  # given both, the V is what is used, so the V's convention is what applies
  s_both <- admStudy(E = EE, V = diag(c(.16, .09, .04)), sd = c(.4, .3, .2),
                     n = 60L, times = tt, ev = ev)
  expect_identical(s_both$v_denom, "ml")

  # explicit always wins
  expect_identical(admStudy(E = EE, sd = c(.4, .3, .2), n = 60L, times = tt,
                            ev = ev, v_denom = "ml")$v_denom, "ml")
  expect_error(admStudy(E = EE, sd = c(.4, .3, .2), n = 60L, times = tt,
                        ev = ev, v_denom = "n-1"), "must be")
})

test_that("the resolved denominator reaches the conversion, and is shown", {
  ev <- rxode2::et(amt = 100)
  tt <- c(1, 2, 4)
  n  <- 60L
  sdv <- c(.4, .3, .2)
  s <- admStudy(E = c(2, 1.5, 1), sd = sdv, n = n, times = tt, ev = ev)

  # .admVDenom is what applies it: unbiased V scaled by (n-1)/n
  vv <- function(V) if (is.matrix(V)) diag(V) else as.numeric(V)
  got <- admixr2:::.admVDenom(unclass(s), "s")
  expect_equal(vv(got$V), sdv^2 * (n - 1) / n, tolerance = 1e-12)
  # ...and it is idempotent, so normalising twice cannot convert twice
  expect_identical(got$v_denom, "ml")
  expect_equal(vv(admixr2:::.admVDenom(got, "s")$V),
               sdv^2 * (n - 1) / n, tolerance = 1e-12)

  # Materialisation preserves the declared denominator.
  materialised <- admixr2:::.admMaterialise(admStudies(s))$s
  expect_identical(materialised$v_denom, "unbiased")
  expect_equal(vv(admixr2:::.admVDenom(materialised, "s")$V),
               sdv^2 * (n - 1) / n, tolerance = 1e-12)

  # the resolved convention is visible rather than silent
  expect_output(print(s), "unbiased denominator")
  expect_output(print(admStudy(E = c(2, 1.5, 1), V = diag(sdv^2), n = n,
                               times = tt, ev = ev)), "ml denominator")
})

test_that("a `by` level keeps the correlations among the margins it retains", {
  # Dropping a margin must rebuild positional latentR.
  cd <- covDist(SEX = list(values = c(0, 1), probs = c(.45, .55)),
                WT = c(mean = 75, sd = 16), CRCL = c(mean = 90, sd = 25),
                cor = matrix(c(1, 0, 0, 0, 1, .5, 0, .5, 1), 3L, 3L,
                             dimnames = list(c("SEX", "WT", "CRCL"),
                                             c("SEX", "WT", "CRCL"))))
  out <- admixr2:::.admCovDropMargin(cd, "SEX")
  expect_identical(admixr2:::.admCovSpecNames(out), c("WT", "CRCL"))
  R <- out[["latentR"]]
  expect_identical(dim(R), c(2L, 2L))
  expect_equal(R[1L, 2L], 0.5)
  # Dropping either member of a correlated pair would silently change the
  # retained population distribution.
  expect_error(admixr2:::.admCovDropMargin(cd, "CRCL"),
               "cannot materialise correlated 'CRCL' subgroups")
})

test_that("`population` does not accept a caller's own `joint` yet", {
  skip_if_not_installed("rxode2")
  # THE SECOND DOOR. covDist() refuses it, but `population` also takes a plain
  # list and the canon lets an existing sampler through untouched.
  raw <- list(WT  = list(quantile = function(u) stats::qlnorm(u, log(70), .25)),
              AGE = list(quantile = function(u) stats::qgamma(u, 9, 0.3)),
              joint = function(u) cbind(WT = stats::qlnorm(u[, 1], log(70), .25),
                                        AGE = stats::qgamma(u[, 2], 9, 0.3)))
  expect_error(admStudy(model = .sa_model, population = raw, n = 100,
                        dose = 200, times = c(1, 4)),
               "own `joint` sampler")
  # The sampler admixr2 BUILDS from `cor` arrives here on every correlated
  # population and must not be caught by that: `jointOwn` is the difference.
  expect_s3_class(
    admStudy(model = .sa_model, n = 100, dose = 200, times = c(1, 4),
             population = admPopulation(WT = c(mean = 75, sd = 16),
                                        CRCL = c(mean = 92, sd = 24),
                                        cor = c(WT.CRCL = 0.45))),
    "admStudy")
})

test_that("by refuses to discard an opaque joint sampler", {
  cd <- with_opaque_joint(
    covDist(SEX = list(values = c(0, 1), probs = c(.5, .5)),
            WT = c(mean = 75, sd = 16), CRCL = c(mean = 90, sd = 25)),
    function(u) cbind(SEX = as.integer(u[, 1] > .5),
                      WT = stats::qlnorm(u[, 2], log(75), .2),
                      CRCL = stats::qlnorm(u[, 2], log(90), .2)))
  expect_error(admixr2:::.admCovDropMargin(cd, "SEX"), "user-supplied `joint`")
})

test_that("a population is canonicalised on build, so `by` sees the same object", {
  # A plain list never went through the canon, so `cor` was still spelled `cor`
  # and .admCovDropMargin() -- which carries only `latentR` across -- dropped
  # the dependence without a word.
  raw <- list(WT = list(mean = 70, sd = 15, dist = "lnorm"),
              CRCL = list(mean = 90, sd = 22, dist = "lnorm"),
              SEX = list(values = c(0, 1), probs = c(.45, .55)),
              cor = matrix(c(1, .6, 0, .6, 1, 0, 0, 0, 1), 3L, 3L,
                           dimnames = rep(list(c("WT", "CRCL", "SEX")), 2L)))
  expect_equal(admixr2:::.admCovDropMargin(raw, "SEX")[["latentR"]][1L, 2L], 0.6)

  s <- admStudy(model = .sa_model, n = 100, dose = 200, times = c(1, 4),
                population = raw, by = "SEX", label = "s")
  expect_false(is.null(s$population[["latentR"]]))
  expect_null(s$population[["cor"]])
  expect_error(admStudy(model = .sa_model, n = 10, dose = 1, times = 1,
                        population = list(WT = list(mean = 70, sd = 15,
                                                    dist = "weibull"))),
               "not a valid covariate specification")
})

test_that("print() does not misreport a `by =` source", {
  skip_if_not_installed("rxode2")
  # `by` pins its covariate one level per study, so it is out of the
  # conditional set deliberately -- and subtracting that set alone then called
  # SEX "read at an asserted coefficient" for a model that estimates `bsex`,
  # under a "reported by SEX" line saying otherwise.
  pop <- admPopulation(WT = c(mean = 70, sd = 15), SEX = c(male = .55))
  out <- capture.output(print(admStudy(
    model = .sa_model, n = 200, dose = 200, times = c(1, 4),
    population = pop, by = "SEX")))
  expect_true(any(grepl("reported by SEX", out)))
  expect_false(any(grepl("conditional on  nothing", out)))
  # WT genuinely IS asserted in `.sa_model` -- a fixed allometric exponent --
  # so that line stays; SEX must not be on it.
  .asserted <- grep("asserted coefficient", out, value = TRUE)
  expect_length(.asserted, 1L)
  expect_match(.asserted, "WT")
  expect_false(grepl("SEX", .asserted))
})

test_that("a `range` keyed by covariate names must be a list", {
  skip_if_not_installed("rxode2")
  # Read as the short form, `c(WT = 52, CRCL = 118)` attached both numbers to
  # whichever single covariate was conditional, as its low and high -- in the
  # fit and, once the plot was made to agree, in both.
  pop <- admPopulation(WT = c(mean = 70, sd = 15), CRCL = c(mean = 90, sd = 25))
  expect_error(suppressMessages(admixr2:::.admMaterialise(admStudies(
    s = admStudy(model = .sa_fitwt, population = pop, n = 100, dose = 200,
                 times = c(1, 4), range = c(WT = 52, CRCL = 118))))),
    "is not a list")
  # `lo`/`hi` are not covariates, so the short form still stands.
  expect_silent(suppressMessages(admixr2:::.admMaterialise(admStudies(
    s = admStudy(model = .sa_fitwt, population = pop, n = 100, dose = 200,
                 times = c(1, 4), range = c(lo = 52, hi = 118))))))
})

test_that("no enrolled-range warning for a DISCRETE conditional covariate", {
  skip_if_not_installed("rxode2")
  # A discrete margin is enumerated at its declared LEVELS, so a `range` buys
  # nothing and the warning's premise is false -- yet it fired for every source
  # estimating a discrete effect, asking for a span for a two-level factor.
  pop <- admPopulation(WT = c(mean = 70, sd = 15), SEX = c(male = .55))
  st  <- admStudies(s = admStudy(model = .sa_model, n = 100, dose = 200,
                                 times = c(1, 4), population = pop))
  w <- NULL
  withCallingHandlers(
    suppressMessages(admixr2:::.admMaterialise(st)),
    warning = function(z) { w <<- c(w, conditionMessage(z))
                            invokeRestart("muffleWarning") })
  expect_length(w, 0L)

  # A CONTINUOUS one with no range still warns, which is the case it is for.
  w2 <- NULL
  withCallingHandlers(
    suppressMessages(admixr2:::.admMaterialise(admStudies(
      s = admStudy(model = .sa_fitwt, n = 100, dose = 200, times = c(1, 4),
                   population = admPopulation(WT = c(mean = 70, sd = 15)))))),
    warning = function(z) { w2 <<- c(w2, conditionMessage(z))
                            invokeRestart("muffleWarning") })
  expect_length(w2, 1L)
  expect_match(w2, "FULL declared distribution")
})

test_that("conditioning is not a user option any more", {
  skip_if_not_installed("rxode2")
  pop <- admPopulation(WT = c(mean = 70, sd = 15), SEX = c(male = .55))
  # `stratify` is gone: whether a covariate is conditional or marginal follows
  # from whether the source's own model uses it.
  for (arg in list(list(stratify = TRUE), list(stratify = FALSE),
                   list(stratify = "SEX")))
    expect_error(
      do.call(admStudy, c(list(model = .sa_model, n = 100, dose = 200,
                               times = c(1, 4), population = pop), arg)),
      "unused argument")

  # ESTIMATED, not merely READ. `.sa_model` carries WT at a FIXED allometric
  # exponent and estimates `bsex`, so SEX is the only covariate it holds
  # evidence about: conditioning on WT would buy strata and no evidence, and credit
  # the source with information it never earned.
  sn <- admStudy(model = .sa_model, n = 100, dose = 200, times = c(1, 4),
                 population = pop, strata_nodes = 3L,
                 range = list(WT = c(50, 100)), label = "s")
  g <- admixr2:::.admMaterialise(admStudies(s = sn))
  expect_identical(length(g), 2L)          # the 2 SEX levels, and no WT nodes

  # ESTIMATE the exponent and WT becomes conditional, at the caller's
  # precision: `strata_nodes` is the one knob that is not a statement about
  # which covariates are conditional.
  .fit_wt <- function() {
    ini({ tcl <- log(5); tv <- log(50); bsex <- 0.15; bwt <- 0.75
          eta.cl ~ 0.09; add.err <- 0.08 })
    model({ cl <- exp(tcl + eta.cl) * (WT/70)^bwt * exp(bsex * SEX)
            v  <- exp(tv) * (WT/70); cp <- linCmt(); cp ~ add(add.err) })
  }
  g2 <- admixr2:::.admMaterialise(admStudies(s = admStudy(
    model = .fit_wt, n = 100, dose = 200, times = c(1, 4),
    population = pop, strata_nodes = 3L,
    range = list(WT = c(50, 100)), label = "s")))
  expect_identical(length(g2), 3L * 2L)    # 3 WT nodes x 2 SEX levels
  # THE NODES ARE THE CALLER'S, on every stratum: copied only inside the
  # "something was conditional" branch, a null fit took the default 9 while the full
  # fit took this 3, and a different resolution per study is what anova()
  # refuses to compare across.
  expect_true(all(vapply(g2, function(z) identical(z[[".adm_strata_nodes"]], 3L),
                         logical(1))))

  # And `range` reaches a MARGINAL covariate too. This model never mentions
  # CRCL, so CRCL is integrated over -- over the TRUNCATED margin, because the
  # enrolled span is true before admixr2 decides how to use the covariate.
  pop2 <- admPopulation(WT = c(mean = 70, sd = 15), SEX = c(male = .55),
                        CRCL = c(mean = 90, sd = 25))
  sp <- admixr2:::.admMaterialise(admStudies(s = admStudy(
    model = .sa_model, n = 100, dose = 200, times = c(1, 4),
    population = pop2, strata_nodes = 1L,
    range = list(CRCL = c(60, 120)), label = "s")))[[1L]]
  cr <- sp$cov_dist[["CRCL"]]
  expect_gte(admixr2:::.admCovQuantile(cr, 0.001), 60)
  expect_lte(admixr2:::.admCovQuantile(cr, 0.999), 120)

  # An unnamed `range` still needs exactly one conditional covariate to belong to,
  # and `.fit_wt` is conditional on two.
  expect_error(admixr2:::.admMaterialise(admStudies(s = admStudy(
    model = .fit_wt, n = 100, dose = 200, times = c(1, 4),
    population = pop, range = c(50, 100), label = "s"))),
    "does not say which covariate")
  # With exactly one it is accepted, and belongs to that covariate -- SEX here,
  # whose only two levels the range keeps.
  g3 <- suppressMessages(admixr2:::.admMaterialise(admStudies(s = admStudy(
    model = .sa_model, n = 100, dose = 200, times = c(1, 4),
    population = pop, range = c(0, 1), label = "s"))))
  expect_identical(length(g3), 2L)

  # WHICH COVARIATE AN UNNAMED `range` BELONGS TO IS A QUESTION ABOUT THE
  # SOURCE, asked before the analysis model has any say. Resolved against the
  # NARROWED set instead, the same `studies` object worked under a model that
  # reads the covariate and died under the null that drops it -- which is the
  # nested pair this is all for.
  .one_cov <- admStudies(s = admStudy(
    model = .sa_model, n = 100, dose = 200, times = c(1, 4),
    population = pop, range = c(0, 1), label = "s"))
  expect_identical(length(suppressMessages(
    admixr2:::.admMaterialise(.one_cov, analysis_covs = "SEX"))), 2L)
  # The null model reads no covariate: no nodes, and no error either.
  expect_identical(length(suppressMessages(
    admixr2:::.admMaterialise(.one_cov, analysis_covs = character(0)))), 1L)

  # ...and with TWO conditional covariates, narrowing to one used to make the
  # unnamed range unambiguous by accident and attach it to the survivor.
  expect_error(suppressMessages(admixr2:::.admMaterialise(
    admStudies(s = admStudy(model = .fit_wt, n = 100, dose = 200,
                            times = c(1, 4), population = pop,
                            range = c(50, 100), label = "s")),
    analysis_covs = "WT")),
    "does not say which covariate")

  # And the derivation resolves to the covariates themselves, so a study
  # prints the set it found rather than a flag it was handed -- and says why
  # the one it did not is marginal, which is the part a reader cannot see
  # from the model.
  out <- capture.output(print(admStudy(model = .sa_model, n = 100, dose = 200,
                                       times = c(1, 4), population = pop,
                                       label = "s")))
  expect_match(paste(out, collapse = "\n"), "conditional on +SEX")
  expect_match(paste(out, collapse = "\n"), "WT read at an asserted coefficient")
})

test_that("nodes are cut only along what the ANALYSIS model reads", {
  skip_if_not_installed("rxode2")
  # A source is conditional on what its OWN model estimated; which of those
  # admixr2 has to cut into nodes is narrowed to what the analysis model reads.
  # The narrowing is EXACT, not a saving with a cost: if a model's prediction
  # does not move across a source's nodes, the mixture those nodes collapse to
  # is a sufficient statistic for it. Measured at identical parameters, on a
  # source conditional on CRCL and WT with an analysis model reading WT only:
  #
  #   both conditional        50 studies   OFV -3041.72602426
  #   only WT conditional     10 studies   OFV -3041.72606427   diff -0.00004
  #   only CRCL conditional   10 studies   OFV -3115.32302414   diff -73.59700
  #
  # so collapsing the unread direction is free to four decimals and collapsing
  # a read one is not. Against a fully marginal reference the invariance is
  # -0.158 at 3 nodes, +0.00003 at 5 and -0.00000001 at 9 -- it is the
  # quadrature converging, not an approximation being tolerated.
  pop <- admPopulation(WT = c(mean = 70, sd = 15), SEX = c(male = .55))
  .src <- function() {
    ini({ tcl <- log(5); tv <- log(50); bsex <- 0.15; bwt <- 0.75
          eta.cl ~ 0.09; add.err <- 0.08 })
    model({ cl <- exp(tcl + eta.cl) * (WT/70)^bwt * exp(bsex * SEX)
            v  <- exp(tv); cp <- linCmt(); cp ~ add(add.err) })
  }
  st <- admStudies(s = admStudy(model = .src, n = 100, dose = 200,
                                times = c(1, 4), population = pop,
                                strata_nodes = 3L, label = "s"))
  # An analysis model reading both: 3 WT nodes x 2 SEX levels.
  g_both <- suppressWarnings(suppressMessages(
    admixr2:::.admMaterialise(st, analysis_covs = c("WT", "SEX"))))
  expect_identical(length(g_both), 3L * 2L)
  expect_true(all(vapply(g_both, function(z)
    isTRUE(z$cov_dist[["WT"]][[".point"]]), logical(1))))

  # Reading SEX only: the WT nodes collapse, the SEX levels stay, and WT is
  # integrated over the distribution the source declared.
  g_sex <- suppressWarnings(suppressMessages(
    admixr2:::.admMaterialise(st, analysis_covs = "SEX")))
  expect_identical(length(g_sex), 2L)
  expect_false(any(vapply(g_sex, function(z)
    isTRUE(z$cov_dist[["WT"]][[".point"]]), logical(1))))

  # Reading neither -- the null model of a nested pair -- leaves one study,
  # marginal over both, and `n` is conserved whichever way it is cut.
  g_none <- suppressWarnings(suppressMessages(
    admixr2:::.admMaterialise(st, analysis_covs = character(0))))
  expect_identical(length(g_none), 1L)
  for (g in list(g_both, g_sex, g_none))
    expect_equal(sum(vapply(g, function(z) z$n, 0)), 100)
})

test_that("admPopulation guards the data-frame route it advertises", {
  set.seed(4)
  d <- data.frame(WT = rlnorm(60, log(70), .2), CRCL = rlnorm(60, log(90), .25))
  # a discrete margin typed in for a covariate the listing does not report
  p <- admPopulation(SEX = c(male = .55), data = d)
  expect_setequal(admixr2:::.admCovSpecNames(p), c("SEX", "WT", "CRCL"))
  # NA is a hole in the margin, not a level -- sort() used to swallow it
  d2 <- d; d2$SEX <- factor(c(rep("m", 30), rep("f", 29), NA))
  expect_error(admPopulation(data = d2), "missing or non-finite")
  # a matrix `cor` is checked exactly as the named-vector form is
  M <- matrix(c(1, .5, .5, 1), 2L, 2L,
              dimnames = rep(list(c("WT", "SEX")), 2L))
  expect_error(admPopulation(WT = c(mean = 70, sd = 15), SEX = c(male = .55),
                             cor = M), "DISCRETE")
  # positional dimnames are not guessable once the data appends its own columns
  expect_error(admPopulation(data = d, cor = matrix(c(1, .3, .3, 1), 2L, 2L)),
               "must have dimnames")
  expect_message(admPopulation(data = d,
                               cor = matrix(c(1, .3, .3, 1), 2L, 2L,
                                            dimnames = rep(list(c("WT", "CRCL")), 2L))),
                 "REPLACED rather than merged")
})

test_that("admPopulation refuses duplicate covariate names", {
  expect_error(admPopulation(WT = c(mean = 70, sd = 10),
                             WT = c(mean = 80, sd = 12)),
               "covariate names must be unique")
})
