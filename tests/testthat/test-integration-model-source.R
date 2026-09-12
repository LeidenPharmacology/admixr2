# A study generated from a published MODEL is not a sample. Its (E, V) are exact
# functions of theta_src, so `n` sets that study's RELATIVE WEIGHT against the
# others rather than its precision, and there is no sampling law to build a
# standard error from. Reading `n` as precision anyway makes the reported SE
# fall as exactly 1/sqrt(n) -- a factor the analyst chooses by typing a number.
#
# So admixr2 reports NO standard error for a fit that includes a model source,
# and refuses an explicit covMethod rather than honouring one. These tests pin
# that contract, and pin that an ordinary digitised study still gets its
# sandwich -- the failure mode of removing this feature is over-cutting.

TIMES_MS <- c(0.5, 1, 2, 4, 8, 12, 24)
DOSE_MS  <- 200

.ms_published <- function() {
  ini({ tcl <- log(5); tv <- log(50)
        eta.cl ~ 0.05; add.err <- 0.08 })
  model({ cl <- exp(tcl + eta.cl); v <- exp(tv)
          cp <- linCmt(); cp ~ add(add.err) })
}
.ms_fit <- function() {
  ini({ tcl <- log(4); tv <- log(45)
        eta.cl ~ 0.1; add.err <- 0.1 })
  model({ cl <- exp(tcl + eta.cl); v <- exp(tv)
          cp <- linCmt(); cp ~ add(add.err) })
}
.ms_gen <- function(n = 400, times = TIMES_MS) {
  suppressWarnings(suppressMessages(datagen(
    list(t1 = list(times = times, ev = rxode2::et(amt = DOSE_MS), n = n)),
    model = .ms_published,
    control = datagenControl(method = "gh", seed = 1L))))
}
# The same numbers, as a SIMULATED data study rather than a published one. This
# is the control arm: identical (E, V), no model-source marker, so the sandwich
# must still run. Anything that refuses BOTH has over-cut.
.ms_as_data <- function(n = 400, times = TIMES_MS) {
  suppressWarnings(suppressMessages(.admDatagenSim(
    list(t1 = list(times = times, ev = rxode2::et(amt = DOSE_MS), n = n)),
    model = .ms_published,
    control = datagenControl(method = "gh", seed = 1L))))
}
# The control object is built OUTSIDE the suppressMessages(), because the
# refusal is raised there and a test that cannot hear it proves nothing.
.ms_ctl <- function(g, ...) adghControl(studies = g, print = 0L, cores = 2L, ...)
.ms_run <- function(g, ..., mod = .ms_fit) {
  ctl <- .ms_ctl(g, ...)
  suppressMessages(nlmixr2est::nlmixr2(mod, admData(), est = "adgh",
                                       control = ctl))
}
# covMethod = "none" is reported back as an empty string, not the literal word
.ms_no_se <- function(f) {
  expect_false(f$covMethod %in% c("r", "r,s"))
  expect_true(all(is.na(f$parFixedDf[["SE"]])))
}

test_that("a model source fit reports NO standard error, at any n", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  # The point estimates are still combined -- only the uncertainty claim is
  # withdrawn. covMethod resolves to "none" on its own, with a message saying
  # so, rather than silently producing a number that tracks `n`.
  for (n in c(100, 1600)) {
    g <- .ms_gen(n)
    expect_message(.ms_ctl(g), "no sampling law")
    f <- .ms_run(g)
    expect_s3_class(f, "admFit")
    expect_true(all(is.finite(f$parFixedDf[["Estimate"]])),
                info = paste("n =", n))
    .ms_no_se(f)
  }
})

test_that("an EXPLICIT covMethod is refused, not honoured", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  # The whole defect is a plausible-looking SE that moves with a number the
  # analyst typed. Honouring an explicit covMethod would leave that number one
  # argument away, so both sandwich and naive forms are refused by name.
  g <- .ms_gen(400)
  for (cm in c("r", "r,s"))
    expect_error(.ms_ctl(g, covMethod = cm), "published MODEL", info = cm)
  # "none" is what the refusal asks for, so it must be accepted silently
  expect_silent(.ms_ctl(g, covMethod = "none"))
  .ms_no_se(.ms_run(g, covMethod = "none"))
})

test_that("the refusal is keyed on the model source, not on generated moments", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  # Byte-identical (E, V); only the marker differs. A digitised study reports on
  # real patients, so `n` IS its precision and the sandwich is correct for it.
  # This is the arm that fails if the removal cut too deep.
  d <- .ms_as_data(400)
  # Not expect_silent(): .ms_run() already suppresses messages/warnings, and a
  # full nlmixr2() fit also prints rxode2 compilation output on a cold cache --
  # under Config/testthat/parallel harmless but order-dependent, since it only
  # stayed quiet here because an earlier test warmed the cache first. What this
  # test pins is the SE values below, not console silence.
  f <- .ms_run(d, covMethod = "r,s")
  expect_identical(f$covMethod, "r,s")
  se <- stats::setNames(f$parFixedDf[["SE"]], rownames(f$parFixedDf))
  expect_true(all(is.finite(se[c("tcl", "tv", "add.err")])))
  expect_true(all(se[c("tcl", "tv", "add.err")] > 0))
})

test_that(".admDatagenSim is the same numbers without the published claim", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  # The two doors differ ONLY in the claim, never in the moments. If they ever
  # diverge numerically then the marker has started changing the fit rather than
  # changing what may be said about it, which is not what it is for.
  pub <- .ms_gen(400); sim <- .ms_as_data(400)
  expect_equal(pub$t1$E, sim$t1$E, tolerance = 1e-12)
  expect_equal(pub$t1$V, sim$t1$V, tolerance = 1e-12)
  expect_true(isTRUE(pub$t1[[".adm_src"]]))
  expect_null(sim$t1[[".adm_src"]])
})

test_that("a digitised study's SE still scales with its n", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  # The ordinary path is untouched: for real patients 1/sqrt(n) is the right
  # law, and it is only wrong for a study that never sampled anyone.
  se_at <- function(n) {
    f <- .ms_run(.ms_as_data(n), covMethod = "r")
    stats::setNames(f$parFixedDf[["SE"]], rownames(f$parFixedDf))[["tcl"]]
  }
  expect_gt(se_at(100) / se_at(1600), 1.5)
})
