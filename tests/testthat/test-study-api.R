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

# ---------------------------------------------------------------------------
# The default that would otherwise be a silent wrong answer, and the check that
# has to happen at construction because the fit swallows it
# ---------------------------------------------------------------------------

test_that("covMethod upgrades to the sandwich only for a model source", {
  skip_if_not_installed("rxode2")
  src <- suppressWarnings(suppressMessages(
    admStudy(model = .sa_model, se = c(tcl = .04, tv = .03, bsex = .05,
                                       add.err = .004, eta.cl = .01),
             n = 60, dose = 100, times = c(1, 4, 12),
             population = admPopulation(WT = c(mean = 75, sd = 16),
                                        SEX = c(male = 0.5)))))
  dat <- admStudy(E = c(1.6, 0.9, 0.3), sd = c(.4, .25, .1), n = 48,
                  dose = 100, times = c(1, 4, 12))
  expect_true(admixr2:::.admHasModelSource(admStudies(a = src, b = dat)))
  expect_false(admixr2:::.admHasModelSource(admStudies(b = dat)))
  # a source with NO reported uncertainty is not one either -- there is nothing
  # for the sandwich to carry
  bare <- admStudy(model = .sa_model, n = 60, dose = 100, times = c(1, 4, 12),
                   population = admPopulation(WT = c(mean = 75, sd = 16),
                                              SEX = c(male = 0.5)))
  expect_false(admixr2:::.admHasModelSource(admStudies(c = bare)))

  st <- admStudies(a = src, b = dat)
  expect_equal(suppressMessages(
    admixr2:::.admResolveCovMethod("r", st, FALSE)), "r,s")
  # an EXPLICIT covMethod is honoured, including an explicit "r"
  expect_equal(admixr2:::.admResolveCovMethod("r", st, TRUE), "r")
  expect_equal(admixr2:::.admResolveCovMethod("none", st, FALSE), "none")
  expect_equal(admixr2:::.admResolveCovMethod("r", admStudies(b = dat), FALSE), "r")
  # and it reaches the control objects, which is where a fit reads it
  expect_equal(suppressMessages(adghControl(studies = st))$covMethod, "r,s")
  expect_equal(adghControl(studies = st, covMethod = "r")$covMethod, "r")
})

test_that("an incomplete source covariance is refused where the user is standing", {
  skip_if_not_installed("rxode2")
  # THE FAILURE THIS PINS: .admSandwichCov() refuses a source whose C_src does
  # not cover every parameter the source ESTIMATED, and refusing one source
  # refuses the sandwich for the whole fit. datagen() does warn -- but it runs
  # inside the nlmixr2est stack, which swallows it, so covMethod came back "r"
  # with the naive standard errors printed and nothing said anywhere.
  expect_warning(
    admStudy(model = .sa_model, rse = c(tcl = 4.1, tv = 6.0),
             n = 60, dose = 100, times = c(1, 4, 12),
             population = admPopulation(WT = c(mean = 75, sd = 16),
                                        SEX = c(male = 0.5))),
    "no standard error will be reported")
  # naming which rows are missing is the actionable half
  w <- tryCatch(suppressMessages(
    admStudy(model = .sa_model, rse = c(tcl = 4.1), n = 60, dose = 100,
             times = c(1, 4, 12),
             population = admPopulation(WT = c(mean = 75, sd = 16),
                                        SEX = c(male = 0.5)))),
    warning = conditionMessage)
  expect_match(w, "bsex")

  # a malformed matrix ERRORS at construction, naming the argument the user
  # actually typed rather than datagen()'s internal `model_cov`
  C <- diag(2); dimnames(C) <- list(c("tcl", "nope"), c("tcl", "nope"))
  expect_error(admStudy(model = .sa_model, cov = C, n = 60, dose = 100,
                        times = c(1, 4, 12)),
               "`cov` names 'nope'")

  # a fit's own naming (`om.eta.cl`) is accepted and mapped onto the ini() row,
  # so print() and the fit agree about which rows are covered
  nm <- c("tcl", "tv", "bsex", "add.err", "om.eta.cl")
  C2 <- diag(c(.04, .03, .05, .004, .01)^2); dimnames(C2) <- list(nm, nm)
  s <- admStudy(model = .sa_model, cov = C2, n = 60, dose = 100,
                times = c(1, 4, 12))
  expect_equal(rownames(s$cov),
               c("tcl", "tv", "bsex", "add.err", "eta.cl"))
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
  # a source with a model but no uncertainty is the other silent one
  bare <- suppressMessages(
    admStudy(model = .sa_model, population = coh, dose = 200,
             times = c(1, 4, 12), label = "b"))
  out2 <- paste(utils::capture.output(print(admStudies(b = bare))),
                collapse = " ")
  expect_match(out2, "no uncertainty")
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

test_that("a source covariance may use nlmixr2's off-diagonal omega names", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  # `fit$cov` from a correlated two-eta model names the off-diagonal
  # `cov.eta.cl.eta.v`, which cannot be split on dots because eta names contain
  # them. It maps onto the `ini()` row `(eta.cl,eta.v)`.
  m <- function() {
    ini({ tcl <- log(5); tv <- log(50); a <- 0.1
          eta.cl + eta.v ~ c(0.09, 0.01, 0.04) })
    model({ cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
            cp <- linCmt(); cp ~ add(a) })
  }
  u  <- suppressMessages(rxode2::rxode2(m))
  nm <- c("tcl", "tv", "a", "om.eta.cl", "cov.eta.cl.eta.v", "om.eta.v")
  C  <- diag(c(.05, .04, .004, .01, .005, .008)^2); dimnames(C) <- list(nm, nm)
  got <- admixr2:::.admSrcCov(C, u, "src", "cov")
  expect_setequal(rownames(got$cov),
                  c("tcl", "tv", "a", "eta.cl", "(eta.cl,eta.v)", "eta.v"))
  expect_length(got$missing, 0L)
  # an off-diagonal the model does not declare is still refused
  nm2 <- c("tcl", "tv", "a", "om.eta.cl", "cov.eta.cl.eta.zz", "om.eta.v")
  C2 <- diag(rep(.01, 6)); dimnames(C2) <- list(nm2, nm2)
  expect_error(admixr2:::.admSrcCov(C2, u, "src", "cov"), "OFF-DIAGONAL")
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

  # the resolved convention is visible rather than silent
  expect_output(print(s), "unbiased denominator")
  expect_output(print(admStudy(E = c(2, 1.5, 1), V = diag(sdv^2), n = n,
                               times = tt, ev = ev)), "ml denominator")
})
