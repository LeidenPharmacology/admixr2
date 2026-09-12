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

test_that("materialising a subgroup refuses a colliding study name", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  s <- admStudy(model = .sa_model, n = 200, dose = 200, times = c(1, 4, 12),
                population = admPopulation(WT = c(mean = 75, sd = 16),
                                           SEX = c(male = 0.6)),
                by = "SEX", label = "trial")
  other <- list(E = 1, V = 1, n = 10, times = 1, ev = rxode2::et(amt = 1))
  expect_error(admixr2:::.admMaterialise(list(trial_SEX0 = other, trial = s)),
               "duplicate name 'trial_SEX0'")
})

# ---------------------------------------------------------------------------
# The default that would otherwise be a silent wrong answer, and the check that
# has to happen at construction because the fit swallows it
# ---------------------------------------------------------------------------

test_that("a model source withdraws the standard error, explicit or not", {
  skip_if_not_installed("rxode2")
  src <- admStudy(model = .sa_model, n = 60, dose = 100, times = c(1, 4, 12),
                  population = admPopulation(WT = c(mean = 75, sd = 16),
                                             SEX = c(male = 0.5)))
  dat <- admStudy(E = c(1.6, 0.9, 0.3), sd = c(.4, .25, .1), n = 48,
                  dose = 100, times = c(1, 4, 12))
  # A published model is a model source whatever it does or does not report
  # about its own uncertainty -- the marker is the MODEL, because that is what
  # makes the study something other than a sample.
  expect_true(admixr2:::.admHasModelSource(admStudies(a = src, b = dat)))
  expect_false(admixr2:::.admHasModelSource(admStudies(b = dat)))

  st <- admStudies(a = src, b = dat)
  expect_equal(suppressMessages(
    admixr2:::.admResolveCovMethod("r", st, FALSE)), "none")
  # AN EXPLICIT covMethod IS REFUSED, NOT HONOURED. Honouring one would leave
  # an SE that tracks the `n` the analyst typed exactly one argument away.
  expect_error(admixr2:::.admResolveCovMethod("r", st, TRUE), "published MODEL")
  expect_error(admixr2:::.admResolveCovMethod("r,s", st, TRUE), "published MODEL")
  # "none" is what the refusal asks for, so it passes through either way
  expect_equal(admixr2:::.admResolveCovMethod("none", st, TRUE), "none")
  expect_equal(admixr2:::.admResolveCovMethod("none", st, FALSE), "none")
  # and a fit with no model source is untouched
  expect_equal(admixr2:::.admResolveCovMethod("r", admStudies(b = dat), FALSE), "r")
  expect_equal(admixr2:::.admResolveCovMethod("r,s", admStudies(b = dat), TRUE),
               "r,s")
  # it reaches the control objects, which is where a fit reads it
  expect_equal(suppressMessages(adghControl(studies = st))$covMethod, "none")
  expect_error(adghControl(studies = st, covMethod = "r"), "published MODEL")
})

test_that("model-source provenance is recursive through observations", {
  nested <- list(pub = list(observations = list(
    cp = list(E = 1, V = 1, n = 20, times = 1, .adm_src = TRUE))))
  expect_true(admixr2:::.admHasModelSource(nested))
  expect_error(admixr2:::.admResolveCovMethod("r", nested, TRUE),
               "published MODEL")
})

# ---------------------------------------------------------------------------
# Reading a baseline table off a cohort instead of transcribing it
# ---------------------------------------------------------------------------

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
  # THE STEP THIS EXISTS FOR: the copula runs on the LATENT scale, so the
  # correlation wanted is cor(log(WT), log(CRCL)) and not cor(WT, CRCL) --
  # close enough to look right by hand and wrong enough to matter.
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

test_that("digitised profiles refuse source-only expansion arguments", {
  args <- list(E = 1, sd = 1, n = 20, dose = 1, times = 1,
               population = admPopulation(SEX = c(male = 0.6)))
  expect_error(do.call(admStudy, c(args, list(by = "SEX"))), "digitised data")
  expect_error(do.call(admStudy, c(args, list(stratify = "SEX"))),
               "digitised data")
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
  # a discrete covariate associated with a continuous one CANNOT be represented
  # (a level would be a truncation of the latent normal), so it is dropped --
  # out loud, because nothing downstream would show it
  coh2 <- coh; coh2$SEX <- as.integer(coh$WT > stats::median(coh$WT))
  expect_message(admPopulation(data = coh2), "being DROPPED")
  expect_silent(invisible(admPopulation(data = coh)))
})

test_that("stratify = TRUE accepts the parsed model admStudy() hands down", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  # admStudy() parses at construction so the transcription can be checked, so
  # the ONLY thing .admExpandStrata() ever sees from that route is an rxUi.
  # Demanding a function rejected the whole admStudy() path -- and only that
  # path, since a raw datagen() spec passes the function through.
  s <- admStudy(model = .sa_model, n = 200, dose = 200, times = c(1, 4, 12),
                population = admPopulation(WT = c(mean = 75, sd = 16),
                                           SEX = c(male = 0.55)),
                stratify = TRUE, label = "trial")
  got <- suppressWarnings(suppressMessages(
    admixr2:::.admMaterialise(admStudies(s))))
  # .sa_model estimates bsex and reads WT/SEX, so TRUE bands on SEX only
  expect_gt(length(got), 1L)
  expect_true(all(vapply(got, function(g) !is.null(g[["cov"]][["SEX"]]),
                         logical(1))))
  expect_equal(sum(vapply(got, function(g) g$n, 0)), 200)
})

test_that("a cohort data frame IS a population, and n comes from it", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  coh <- .sa_cohort(n = 300L)
  coh$ID <- seq_len(nrow(coh))          # a real cohort carries data columns too
  s <- suppressMessages(
    admStudy(model = .sa_model, population = coh, dose = 200,
             times = c(1, 4, 12), label = "trial"))
  expect_equal(s$n, 300)
  # EVERY covariate survives, not just the ones this source's model reads.
  # .sa_model never mentions CRCL -- and a covariate no published model fitted
  # is exactly what a meta-analysis is for, so dropping it would delete the
  # evidence that identifies its effect and leave a fit that still converges.
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
  # THE FAILURE IT PRE-EMPTS: a covariate marginalised identically everywhere
  # is not identified, and converges anyway -- 0.019 objective units across its
  # whole range, settling at -0.059 against a truth of +0.150. Nothing after
  # the fit shows it, so it has to be said before.
  coh <- .sa_cohort(n = 300L)
  s <- suppressMessages(
    admStudy(model = .sa_model, population = coh, dose = 200,
             times = c(1, 4, 12), stratify = "SEX", label = "a"))
  out <- paste(utils::capture.output(print(admStudies(a = s))), collapse = " ")
  expect_match(out, "SEX +banded")
  expect_match(out, "WT +marginal")
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
  # Left to run, mean/sd gives sd = 0 -> sdlog = 0, and the refusal surfaced
  # much later inside datagen() as "not a supported distribution" -- far from
  # the column that caused it.
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

test_that("stratify = TRUE bands what the source ESTIMATED, not what it reads", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  # .sa_model fits bsex and holds weight at a literal ^0.75, so TRUE must band
  # SEX alone. Banding WT too would manufacture the null contrast that
  # .admExpandStrata exists to avoid.
  s <- suppressMessages(
    admStudy(model = .sa_model, n = 200, dose = 200, times = c(1, 4, 12),
             population = admPopulation(WT = c(mean = 75, sd = 16),
                                        SEX = c(male = 0.55)),
             stratify = TRUE, label = "trial"))
  got <- suppressWarnings(suppressMessages(
    admixr2:::.admMaterialise(admStudies(s))))
  expect_setequal(unlist(lapply(got, function(g) names(g[["cov"]]))), "SEX")
  expect_equal(sum(vapply(got, function(g) g$n, 0)), 200)
})

test_that("admStudy resolves v_denom from the currency the study is written in", {
  ev <- rxode2::et(amt = 100)
  tt <- c(1, 2, 4); EE <- c(2, 1.5, 1)

  # A PUBLISHED SPREAD IS THE n-1 ONE. admStudy() had no `v_denom` at all, so
  # every study it built defaulted to "ml" -- and `sd =` is precisely the
  # digitised-figure path the constructor exists for, so the one convention a
  # paper never uses was the one silently assumed.
  s_sd <- admStudy(E = EE, sd = c(.4, .3, .2), n = 60L, times = tt, ev = ev)
  expect_identical(s_sd$v_denom, "unbiased")
  s_sem <- admStudy(E = EE, sem = c(.4, .3, .2) / sqrt(60), n = 60L,
                    times = tt, ev = ev)
  expect_identical(s_sem$v_denom, "unbiased")

  # a covariance MATRIX is computed, not transcribed -- the docs tell you to
  # use the ML denominator, so guessing "unbiased" there would corrupt it
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

  # Materialisation is the production path; it must preserve the convention
  # for the subsequent study normalisation rather than defaulting it to ML.
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
  # `cdk[[by]] <- NULL` alone leaves latentR -- indexed POSITIONALLY, no
  # dimnames -- describing the ORIGINAL set, so .admCovCollapse read the
  # (SEX, WT) block, i.e. the identity, and a declared WT-CRCL correlation of
  # 0.5 became independence in every by-level study.
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
  # a single surviving margin needs no correlation and must not error
  expect_identical(admixr2:::.admCovSpecNames(
    admixr2:::.admCovDropMargin(cd, "CRCL")), c("SEX", "WT"))
})

test_that("by refuses to discard an opaque joint sampler", {
  cd <- covDist(SEX = list(values = c(0, 1), probs = c(.5, .5)),
                WT = c(mean = 75, sd = 16), CRCL = c(mean = 90, sd = 25),
                joint = function(u) cbind(SEX = as.integer(u[, 1] > .5),
                                           WT = stats::qlnorm(u[, 2], log(75), .2),
                                           CRCL = stats::qlnorm(u[, 2], log(90), .2)))
  expect_error(admixr2:::.admCovDropMargin(cd, "SEX"), "user-supplied `joint`")
})
