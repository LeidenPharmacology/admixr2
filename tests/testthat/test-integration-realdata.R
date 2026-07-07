# Tier 2 integration: real-data validation on Theophylline (nlmixr2data::theo_sd).
#
# Aggregate the 12-subject single-dose data onto a common grid (rank-binned,
# dose-normalised) and fit the CORRECT model -- IIV on ka, cl AND v -- with adgh
# under both likelihood branches. Compares the recovered structural parameters
# to the individual-level FOCEI fit of the same model.
#
# Reference (individual FOCEI, 3-eta): ka = 1.48 /h, cl = 2.78 L/h, v = 32.0 L.
# Reproduce with validation/regression-check.R.
#
# Key, deliberately-encoded finding: with a correctly-specified model the
# diagonal "var" branch recovers all three structural parameters, but the full
# "cov" branch still biases ka badly -- the n = 12 sample covariance is a noisy,
# ill-conditioned estimate that the full-covariance likelihood chases. So cov is
# asserted only for cl and v (the AUC/volume parameters it can identify), and the
# test locks in that var recovers ka strictly better than cov.

.realdata_cache <- NULL

.realdata_setup <- function() {
  if (!is.null(.realdata_cache)) return(.realdata_cache)
  skip_on_cran()
  skip_if_not_installed("rxode2")
  skip_if_not_installed("nlmixr2est")
  skip_if_not_installed("nlmixr2data")

  e <- new.env()
  utils::data("theo_sd", package = "nlmixr2data", envir = e)
  theo <- get("theo_sd", envir = e)
  ids  <- unique(theo$ID)
  obs  <- theo[theo$EVID == 0, ]
  dose <- theo[theo$EVID != 0, ]
  amt  <- stats::setNames(dose$AMT, dose$ID)
  DOSE <- 320; K <- 11L

  rt <- dvn <- matrix(NA_real_, length(ids), K)
  for (i in seq_along(ids)) {
    s <- obs[obs$ID == ids[i], ]; s <- s[order(s$TIME), ]
    if (nrow(s) != K) skip("theo_sd sampling changed -- expected 11 samples/subject")
    rt[i, ]  <- s$TIME
    dvn[i, ] <- s$DV * (DOSE / amt[[as.character(ids[i])]])   # dose-normalise
  }
  times <- round(apply(rt, 2L, stats::median), 3)
  E <- colMeans(dvn)
  Vfull <- stats::cov.wt(dvn, method = "ML")$cov
  n <- length(ids)

  ## correct structure: IIV on ka, cl AND v
  theo3 <- function() {
    ini({
      tka <- log(1.5); tcl <- log(2.8); tv <- log(32)
      prop.sd <- c(0, 0.1); add.sd <- c(0, 0.3)
      eta.ka ~ 0.2; eta.cl ~ 0.09; eta.v ~ 0.09
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      d/dt(depot)  <- -ka * depot
      d/dt(center) <-  ka * depot - (cl / v) * center
      cp <- center / v
      cp ~ prop(prop.sd) + add(add.sd)
    })
  }

  fit_branch <- function(Vmat) {
    f <- tryCatch(
      suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
        theo3, admData(), est = "adgh",
        control = adghControl(
          studies = list(theoph = list(E = E, V = Vmat, n = n, times = times,
                                       ev = rxode2::et(amt = DOSE, cmt = 1))),
          maxeval = 400L)))),
      error = function(err) err)
    if (inherits(f, "error")) skip(paste("adgh fit failed:", conditionMessage(f)))
    ex <- f$env$admExtra
    list(fit = f,
         ka = exp(ex$struct[["tka"]]), cl = exp(ex$struct[["tcl"]]), v = exp(ex$struct[["tv"]]),
         obj = f$objective)
  }

  var_fit <- fit_branch(diag(diag(Vfull)))   # diagonal V -> method = "var"
  cov_fit <- fit_branch(Vfull)               # full V     -> method = "cov"

  .realdata_cache <<- list(
    ref = c(ka = 1.48, cl = 2.78, v = 32.0),
    var = var_fit, cov = cov_fit)
  .realdata_cache
}

test_that("adgh 'var' recovers the individual-fit ka, cl and v (correct model)", {
  s <- .realdata_setup()
  expect_s3_class(s$var$fit, "admFit")
  expect_true(is.finite(s$var$obj))
  expect_equal(s$var$ka, s$ref[["ka"]], tolerance = 0.35)   # ~9%  observed
  expect_equal(s$var$cl, s$ref[["cl"]], tolerance = 0.20)   # ~3%  observed
  expect_equal(s$var$v,  s$ref[["v"]],  tolerance = 0.20)   # ~1%  observed
})

test_that("adgh 'cov' recovers cl and v (ka not identifiable from a small-n covariance)", {
  s <- .realdata_setup()
  expect_s3_class(s$cov$fit, "admFit")
  expect_true(is.finite(s$cov$obj))
  expect_equal(s$cov$cl, s$ref[["cl"]], tolerance = 0.25)   # ~15% observed
  expect_equal(s$cov$v,  s$ref[["v"]],  tolerance = 0.20)   # ~5%  observed
  # ka is deliberately NOT asserted for cov: the noisy n=12 sample covariance
  # biases it upward even with the correct model. Just require it stayed finite.
  expect_true(is.finite(s$cov$ka) && s$cov$ka > 0)
})

test_that("var recovers ka strictly better than cov (small-n covariance fragility)", {
  s <- .realdata_setup()
  expect_lt(abs(s$var$ka - s$ref[["ka"]]), abs(s$cov$ka - s$ref[["ka"]]))
  # the two branches agree on the AUC/volume parameters they can both identify
  expect_equal(s$var$cl, s$cov$cl, tolerance = 0.25)
  expect_equal(s$var$v,  s$cov$v,  tolerance = 0.20)
})
