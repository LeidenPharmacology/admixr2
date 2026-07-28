# beta() endpoints, end to end.
#
# A beta endpoint is the one residual family whose variance needs a quantity that
# is neither observed nor fitted: with `y ~ beta(b1, b2)` the prediction is the
# DERIVED mean mu = b1/(b1+b2) and the conditional variance is mu(1-mu)/(1+phi)
# with phi = b1 + b2 SOLVED from the structural model. Every path that turns a
# solve into a prediction matrix therefore has to combine the pair AND carry phi
# back, and each of them used to do neither:
#
#   .admNLLBatch   -- the covMethod = "r" objective evaluator: read b1 alone, so
#                     the Hessian came from a different function than the fit
#   .admGrad (FD)  -- same read, plus arr$phi left at NA -> an all-NA gradient
#   .admSimulateRows / .adghMomentsBatch -- phi dropped by the matrix subset
#   datagen()      -- emitted E = a shape parameter and a V that was entirely NA
#   plot.admFit    -- all-NA predicted covariance after a converged fit
#
# None of that raised an error anywhere, which is why this file exists.

.beta_model <- function() {
  ini({ tem  <- log(0.6)
        tphi <- log(30)
        eta.em ~ 0.05 })
  model({
    em  <- exp(tem + eta.em)
    phi <- exp(tphi)
    mu  <- em / (1 + em)                 # in (0, 1)
    b1  <- mu * phi
    b2  <- (1 - mu) * phi
    y ~ beta(b1, b2)
  })
}

.beta_cache <- NULL
.beta_setup <- function() {
  if (!is.null(.beta_cache)) return(.beta_cache)
  skip_on_cran(); skip_if_not_installed("rxode2"); skip_if_not_installed("nlmixr2est")
  gen <- datagen(
    studies = list(s1 = list(times = c(1, 2, 4), ev = rxode2::et(amt = 0),
                             n = 300L)),
    model = .beta_model, control = datagenControl(n_sim = 20000L, seed = 1L))
  fit <- suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
    .beta_model, admData(), est = "admc",
    control = admControl(studies = gen, n_sim = 3000L, maxeval = 150L,
                         seed = 7L, covMethod = "r"))))
  .beta_cache <<- list(gen = gen, fit = fit)
  .beta_cache
}

test_that("datagen() emits a beta study on the probability scale, not a shape", {
  s <- .beta_setup()$gen$s1
  # mu = expit(log(0.6)) = 0.375 up to the IIV correction -- a probability. The
  # bug returned b1 (~11), an arbitrary positive number, with no complaint.
  expect_true(all(s$E > 0 & s$E < 1))
  expect_equal(unname(s$E), rep(0.6 / 1.6, 3L), tolerance = 0.05)
  expect_true(all(is.finite(s$V)))
  expect_true(all(diag(s$V) > 0))
})

test_that("admc recovers the generating parameters of a beta endpoint", {
  fit <- .beta_setup()$fit
  expect_s3_class(fit, "admFit")
  expect_true(is.finite(fit$objective))
  ex <- fit$env$admExtra
  expect_equal(unname(ex$struct[["tem"]]),  log(0.6), tolerance = 0.05)
  expect_equal(unname(ex$struct[["tphi"]]), log(30),  tolerance = 0.05)
})

test_that("a beta fit is driven derivative-free, and says so", {
  # The analytic and FD-of-the-prediction gradients both chain through mu only,
  # holding phi at its pre-perturbation value -- but a structural theta moves phi
  # too. BOBYQA differences the objective, where phi moves with it.
  skip_on_cran(); skip_if_not_installed("rxode2")
  gen <- .beta_setup()$gen
  expect_message(
    suppressWarnings(nlmixr2est::nlmixr2(
      .beta_model, admData(), est = "admc",
      control = admControl(studies = gen, n_sim = 500L, maxeval = 2L,
                           seed = 7L, covMethod = "none", grad = "sens"))),
    "derivative-free")
})

test_that("the beta covariance is real, and on the reported scale", {
  fit <- .beta_setup()$fit
  skip_if(is.null(fit$cov), "covariance not computed")
  # .admNLLBatch is the evaluator behind covMethod = "r". Reading b1 instead of
  # b1/(b1+b2) made it score a different objective from the one minimised, which
  # shows up as SEs that are finite and plausible but wrong -- so check the
  # objective the Hessian is built on directly, not only that the SEs exist.
  expect_true(all(is.finite(fit$cov)))
  expect_true(all(diag(fit$cov) > 0))
  pf <- fit$parFixedDf
  for (nm in intersect(rownames(pf), rownames(fit$cov)))
    expect_equal(unname(pf[nm, "SE"]), sqrt(fit$cov[nm, nm]), tolerance = 1e-8,
                 info = nm)
})

test_that(".admNLLBatch scores the same beta objective as .admNLL", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  env   <- .beta_setup()
  ui    <- suppressMessages(rxode2::rxode2(.beta_model()))
  pinfo <- suppressWarnings(admixr2:::.admParseIniDf(ui$iniDf, ui))
  rxMod <- admixr2:::.admLoadModel(ui)
  st    <- admixr2:::.admBuildEvFull(list(admixr2:::.admNormaliseStudy(
    c(env$gen$s1, list(ev = rxode2::et(amt = 0))), "s1")))
  st[[1L]]$out_pair <- admixr2:::.admBetaPair(ui)
  z  <- admixr2:::.admMakeZ(400L, pinfo, 1L, "sobol")
  pm <- admixr2:::.admMakeParamsList(400L, pinfo, 1L)
  p0 <- admixr2:::.admBuildOptVec(pinfo)$p0

  n1 <- admixr2:::.admNLL(p0, pinfo, st, z, rxMod, admixr2:::.admOutputVar(ui), pm, 1L)
  nb <- admixr2:::.admNLLBatch(list(p0), pinfo, st, z, rxMod,
                               admixr2:::.admOutputVar(ui), pm, 1L)
  expect_true(is.finite(n1))
  expect_equal(nb, n1, tolerance = 1e-10)
})

test_that("plot.admFit gives a beta fit a real predicted covariance", {
  fit <- .beta_setup()$fit
  ad  <- fit$env$aggData
  skip_if(is.null(ad), "aggregate diagnostics unavailable")
  expect_true(all(is.finite(ad[[1L]]$pred$V)))
  expect_true(all(diag(ad[[1L]]$pred$V) > 0))
  expect_true(all(ad[[1L]]$pred$E > 0 & ad[[1L]]$pred$E < 1))
  skip_if_not_installed("ggplot2")
  p <- plot(fit, which = "cov")
  expect_true(length(p) > 0L)
})

test_that("datagen(method = 'fo') refuses a beta endpoint rather than emitting NAs", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  expect_error(
    datagen(studies = list(s1 = list(times = c(1, 2), ev = rxode2::et(amt = 0),
                                     n = 50L)),
            model = .beta_model,
            control = datagenControl(method = "fo", seed = 1L)),
    "beta")
})
