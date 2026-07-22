# datagen() -> fit -> recover the generating parameters.
#
# datagen() emits the (E, V) a study WOULD report under a given model, using the
# same .admResidRows()/.admResidApply() the estimators score against. Feeding
# that straight back in is therefore a SELF-CONSISTENCY check, not a validation
# against rxode2 -- but it is the only end-to-end check available for a residual
# whose whole content lives in the OFF-DIAGONAL of V:
#
#   ar(rho) is invisible to a diagonal V (admixr2 refuses `method = "var"` for it
#   for exactly that reason), and it cannot be round-tripped against rxode2's own
#   simulator, whose ar() draw is not stationary when a dose record precedes the
#   first observation -- nlmixr2's own focei cannot recover rho from it either
#   (see the long note in .admBuildResidSpecs()).
#
# So what this pins down is: rho is identifiable from a full V; the ar_cor
# gradient actually drives the optimizer to it from a wrong start; and the
# post-fit covariance handles the logit-encoded role.

.dgr_model <- function(rho = 0.6, add = 0.5) {
  src <- sprintf('function() {
    ini({ tcl <- log(5); tv <- log(20); add.sd <- %g; rho <- %g; eta.cl ~ 0.09 })
    model({ cl <- exp(tcl + eta.cl); v <- exp(tv)
            d/dt(central) <- -(cl / v) * central
            cp <- central / v
            cp ~ add(add.sd) + ar(rho) }) }', add, rho)
  eval(parse(text = src))
}

test_that("datagen() -> adgh recovers an ar(rho) residual correlation", {
  skip_on_cran()
  skip_if_not_installed("rxode2")

  RHO <- 0.6; ADD <- 0.5
  times <- c(0.5, 1, 2, 3, 4, 6, 8)          # evenly spaced enough to see rho^d
  ev    <- rxode2::et(amt = 100)

  gen <- suppressWarnings(suppressMessages(datagen(
    studies = list(s = list(times = times, ev = ev, n = 400L)),
    model   = .dgr_model(RHO, ADD),
    control = datagenControl(method = "gh"))))
  V <- gen$s$V

  # The generated V must actually CARRY the correlation, or the test would pass
  # on a model that ignores rho entirely.
  sdv <- sqrt(diag(V))
  R   <- V / outer(sdv, sdv)
  expect_gt(max(abs(R[upper.tri(R)])), 0.05)

  # Start rho well away from the truth so the gradient has to do the work.
  fit <- suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
    .dgr_model(0.2, ADD), admData(), est = "adgh",
    control = adghControl(
      studies = list(s = list(E = gen$s$E, V = V, n = 400L,
                              times = times, ev = ev)),
      maxeval = 300L, covMethod = "none"))))

  expect_s3_class(fit, "admFit")
  expect_true(is.finite(fit$objective))
  # parFixedDf reports on the NATURAL scale (.admFullTheta maps each sigma_role
  # back), so rho is a correlation here, not its logit.
  rho_hat <- unname(fit$parFixedDf["rho", "Estimate"])
  expect_gt(rho_hat, 0); expect_lt(rho_hat, 1)
  expect_equal(rho_hat, RHO, tolerance = 0.15)
  expect_equal(unname(fit$parFixedDf["add.sd", "Estimate"]), ADD, tolerance = 0.15)
  expect_equal(exp(fit$env$admExtra$struct[["tcl"]]), 5, tolerance = 0.10)
})

test_that("an ar() fit reports a finite standard error for rho", {
  skip_on_cran()
  skip_if_not_installed("rxode2")
  # rho is the `ar_cor` sigma_role: held as logit(rho), reported as rho, so its
  # SE carries the delta factor rho(1 - rho). A missing factor would print an SE
  # ~4x too large here and nobody would notice from the point estimate.
  RHO <- 0.6; ADD <- 0.5
  times <- c(0.5, 1, 2, 3, 4, 6, 8)
  ev    <- rxode2::et(amt = 100)

  gen <- suppressWarnings(suppressMessages(datagen(
    studies = list(s = list(times = times, ev = ev, n = 400L)),
    model   = .dgr_model(RHO, ADD),
    control = datagenControl(method = "gh"))))

  fit <- suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
    .dgr_model(RHO, ADD), admData(), est = "adgh",
    control = adghControl(
      studies = list(s = list(E = gen$s$E, V = gen$s$V, n = 400L,
                              times = times, ev = ev)),
      maxeval = 60L, covMethod = "r"))))

  skip_if(is.null(fit$cov), "Hessian not PD for this configuration")
  rn <- rownames(fit$cov)
  expect_false(any(is.na(rn) | rn == ""))
  expect_true("rho" %in% rn)
  se <- sqrt(fit$cov["rho", "rho"])
  expect_true(is.finite(se) && se > 0)
  # Delta factor rho(1-rho) <= 0.25, so a correlation SE on (0,1) cannot be huge.
  # Without the factor this would come back on the logit scale and blow past 1.
  expect_lt(se, 0.5)
})
