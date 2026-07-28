# Edge cases surfaced by a multi-agent probe of the estimators (2026-07).
#   - a no-IIV (zero-eta / population-only) model: admc crashed with "argument is
#     not a matrix" (t(NULL)) where adfo/adgh guard it, and its grad-FD covariance
#     had no eta gradient to batch, so it returned NA SEs. adirmc, which draws its
#     proposals FROM the random effects, returned a silent objective = Inf.
#   - a bounded transform-both-sides residual whose PREDICTION crosses the bound at
#     an optimizer trial point: .admTBS returns NaN, and `any(dms != 0)` (without
#     na.rm) is NA, so the default analytic gradient crashed the fit.
skip_on_cran()
skip_if_not_installed("rxode2")
skip_if_not_installed("nlmixr2est")

.edge_study <- function() {
  set.seed(1L); N <- 300L; tm <- c(0.5, 1, 2, 4, 8); cl <- exp(log(3)); v <- exp(log(20))
  f <- outer(rep(cl, N), tm, function(c_, t_) 100 / v * exp(-(c_ / v) * t_))
  y <- f + 0.3 * matrix(stats::rnorm(length(f)), N)      # data with NO IIV
  list(E = colMeans(y), V = cov.wt(y, method = "ML")$cov, n = N, times = tm,
       ev = rxode2::et(amt = 100))
}
.edge_model0 <- function() {          # no eta -> zero-eta / population-only
  ini({ tcl <- log(3); tv <- log(20); a <- 0.3 })
  model({ cl <- exp(tcl); v <- exp(tv)
          d/dt(central) <- -(cl/v)*central; cp <- central/v; cp ~ add(a) })
}

test_that("a no-IIV (zero-eta) model fits with finite SEs under admc/adfo/adgh", {
  st <- .edge_study()
  ests <- c(admc = NA, adfo = NA, adgh = NA)
  se_a <- numeric(0)
  for (es in names(ests)) {
    ctl <- switch(es,
      admc = admControl(studies = list(s = st), covMethod = "r"),
      adfo = adfoControl(studies = list(s = st), covMethod = "r"),
      adgh = adghControl(studies = list(s = st), covMethod = "r"))
    fit <- suppressWarnings(suppressMessages(
      nlmixr2est::nlmixr2(.edge_model0, admData(), est = es, control = ctl)))
    expect_s3_class(fit, "admFit")
    expect_true(is.finite(fit$objective), label = paste(es, "objective"))
    expect_equal(unname(fit$parFixedDf["tcl", "Estimate"]), log(3), tolerance = 0.05,
                 label = paste(es, "tcl"))
    expect_true(is.finite(fit$parFixedDf["a", "SE"]), label = paste(es, "SE(a)"))
    se_a[es] <- fit$parFixedDf["a", "SE"]
  }
  # the three estimators must AGREE on the standard error (they fit the same model)
  expect_lt(max(abs(se_a - stats::median(se_a)) / stats::median(se_a)), 0.05,
            label = "SE(a) agreement across estimators")
})

test_that("adirmc refuses a no-IIV model cleanly (it needs random effects to propose)", {
  st <- .edge_study()
  expect_error(
    nlmixr2est::nlmixr2(.edge_model0, admData(), est = "adirmc",
                        control = adirmcControl(studies = list(s = st))),
    "requires at least one random effect")
})

test_that("a bounded TBS residual crossing its bound still fits under the analytic grad", {
  # logitNorm with a large SD: the +-12 SD quadrature AND a low prediction push f
  # toward/through the bound at optimizer trial points -> NaN dms; without na.rm the
  # analytic gradient (the adgh default) threw "missing value where TRUE/FALSE needed".
  set.seed(1L); N <- 500L; tm <- c(0.25, 0.5, 1, 2, 4, 8, 12)
  eta <- stats::rnorm(N, 0, 0.3); cl <- exp(log(3) + eta); v <- 30
  f <- outer(cl, tm, function(c_, t_) 100 / v * exp(-(c_ / v) * t_))
  y <- admixr2:::.admTBSi(admixr2:::.admTBS(f, 1, 4L, 0, 8) +
                            1.0 * matrix(stats::rnorm(length(f)), N), 1, 4L, 0, 8)
  st <- list(E = colMeans(y), V = cov.wt(y, method = "ML")$cov, n = N, times = tm,
             ev = rxode2::et(amt = 100))
  mf <- function() {
    ini({ tcl <- log(3); tv <- log(30); a <- 1.0; eta.cl ~ 0.09 })
    model({ cl <- exp(tcl + eta.cl); v <- exp(tv)
            d/dt(central) <- -(cl/v)*central; cp <- central/v
            cp ~ logitNorm(a, 0, 8) }) }
  for (es in c("adgh", "admc")) {
    ctl <- switch(es, adgh = adghControl(studies = list(s = st), covMethod = "none"),
                      admc = admControl(studies = list(s = st), covMethod = "none"))
    fit <- suppressWarnings(suppressMessages(
      nlmixr2est::nlmixr2(mf, admData(), est = es, control = ctl)))
    expect_s3_class(fit, "admFit")
    expect_true(is.finite(fit$objective), label = paste(es, "objective"))
    expect_equal(unname(fit$parFixedDf["a", "Estimate"]), 1.0, tolerance = 0.2,
                 label = paste(es, "a"))
  }
})
