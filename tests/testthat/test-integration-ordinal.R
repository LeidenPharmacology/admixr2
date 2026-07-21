# Tier 2 integration: ordinal endpoints (`y ~ c(p1, p2)`).
#
# An ordinal observation is the VECTOR of K-1 category indicators at each time, so
# a study contributes a stacked (n_times * (K-1)) mean and its covariance. admixr2
# expresses that as a JOINT same-subject unit with one observation block per
# category probability -- the same machinery multi-compartment fits use.
#
# The load-bearing fact, and what most of this file exists to pin down: for two
# categories observed at the SAME time the law of total covariance gives
#
#   Cov(1_j, 1_k) = E[Cov(1_j,1_k | eta)] + Cov_eta(p_j, p_k)
#                 = -E[p_j]E[p_k] - Cov_eta + Cov_eta
#                 = -E[p_j] E[p_k]
#
# i.e. the STRUCTURAL covariance cancels exactly. The off-diagonal must therefore
# REPLACE V_struct rather than add to it. Getting that wrong leaves the entry too
# large by exactly Cov_eta(p_j, p_k), which is invisible when omega is small --
# hence the deliberately large omega in the oracle below.

.ord_model <- function() {
  function() {
    ini({ tcl <- log(5); tv <- log(20); eta.cl ~ 0.36 })
    model({
      cl <- exp(tcl + eta.cl); v <- exp(tv)
      d/dt(central) <- -(cl / v) * central
      cp <- central / v
      # scaled so p1 + p2 <= 0.8: a valid ordinal spec needs sum(p) <= 1
      p1 <- 0.45 / (1 + exp(-(cp - 2)))
      p2 <- 0.35 / (1 + exp(-(cp - 4)))
      y ~ c(p1, p2)
    })
  }
}

.ord_cache <- NULL
.ord_setup <- function() {
  if (!is.null(.ord_cache)) return(.ord_cache)
  skip_on_cran(); skip_if_not_installed("rxode2"); skip_if_not_installed("nlmixr2est")
  ui    <- suppressMessages(rxode2::rxode2(.ord_model()))
  pinfo <- suppressWarnings(.admParseIniDf(ui$iniDf, ui))
  rx    <- .admLoadModel(ui)
  times <- c(0.5, 1, 2, 4)
  ev    <- rxode2::et(amt = 100)

  set.seed(4); N <- 3000L
  pdf <- data.frame(tcl = log(5), tv = log(20),
                    eta.cl = stats::rnorm(N, 0, 0.6), rxerr.y = 0)
  sol <- rxode2::rxSolve(rx, params = pdf, events = rxode2::et(ev, times), cores = 1L,
                         returnType = "data.frame",
                         nDisplayProgress = .Machine$integer.max)
  keep <- sol$time %in% times
  P1 <- matrix(sol$p1[keep], N, length(times), byrow = TRUE)
  P2 <- matrix(sol$p2[keep], N, length(times), byrow = TRUE)
  U  <- matrix(stats::runif(N * length(times)), N, length(times))
  Y  <- cbind((U < P1) * 1, (U >= P1 & U < P1 + P2) * 1)
  nt <- length(times)
  .ord_cache <<- list(
    ui = ui, pinfo = pinfo, times = times, ev = ev, nt = nt, N = N,
    E = colMeans(Y), V = stats::cov.wt(Y, method = "ML")$cov)
  .ord_cache
}

.ord_study <- function(s) list(
  observations = list(
    c1 = list(output = "p1", times = s$times, E = s$E[seq_len(s$nt)],
              n = s$N, ev = s$ev),
    c2 = list(output = "p2", times = s$times, E = s$E[s$nt + seq_len(s$nt)],
              n = s$N, ev = s$ev)),
  V = s$V, n = s$N, ev = s$ev)

test_that("the ordinal spec is registered under EVERY category probability", {
  s  <- .ord_setup()
  sp <- .admResidSpecs(s$pinfo)
  # Registering only the first probability left every other category with no spec
  # at all -- form 0 and zero residual variance.
  expect_setequal(names(sp), c("p1", "p2"))
  expect_true(all(vapply(sp, function(x) identical(x$form, .ADM_RESID_ORDINAL),
                         logical(1))))
})

test_that("an ordinal endpoint gets no sensitivity model (rx_pred_ is the likelihood)", {
  s <- .ord_setup()
  # rx_pred_ for `y ~ c(p1,p2)` is the ordinal LOG-LIKELIHOOD, not a category
  # probability, so its sensitivity columns differentiate the wrong function.
  # NULL here is what routes every estimator onto the finite-difference path.
  expect_null(suppressMessages(.admLoadSensModel(s$ui)))
})

test_that("ordinal aggregate moments match a real multinomial simulation", {
  s <- .ord_setup()
  rx <- .admLoadModel(s$ui)
  set.seed(11); N <- 200000L
  tt <- c(0.5, 2)
  pdf <- data.frame(tcl = log(5), tv = log(20),
                    eta.cl = stats::rnorm(N, 0, 0.8), rxerr.y = 0)
  sol <- rxode2::rxSolve(rx, params = pdf, events = rxode2::et(s$ev, tt), cores = 1L,
                         returnType = "data.frame",
                         nDisplayProgress = .Machine$integer.max)
  keep <- sol$time %in% tt
  P1 <- matrix(sol$p1[keep], N, 2L, byrow = TRUE)
  P2 <- matrix(sol$p2[keep], N, 2L, byrow = TRUE)
  U  <- matrix(stats::runif(N * 2L), N, 2L)
  Y  <- cbind((U < P1) * 1, (U >= P1 & U < P1 + P2) * 1)
  Eemp <- colMeans(Y); Vemp <- stats::cov.wt(Y, method = "ML")$cov

  F   <- cbind(P1, P2); mu <- colMeans(F)
  Vst <- crossprod(sweep(F, 2L, mu)) / N
  rowt <- c(tt, tt)
  arr <- .admResidRows(s$pinfo, c("p1", "p1", "p2", "p2"),
                       .admSigmaNat(s$pinfo$sigma_init, s$pinfo), 4L)
  ap  <- .admResidApply(mu, diag(Vst), arr, rowt, Vst)
  V   <- Vst; diag(V) <- ap$dv
  V   <- V + ap$rmat

  mc_se <- 1 / sqrt(N)
  expect_equal(unname(ap$mu), unname(Eemp), tolerance = 0.01)
  expect_lt(max(abs(V - Vemp)), 5 * mc_se)

  # ... and the version WITHOUT the cancellation is measurably worse, so this
  # test actually discriminates rather than passing on a loose tolerance.
  Vbad <- Vst; diag(Vbad) <- ap$dv
  same <- outer(rowt, rowt, "=="); diag(same) <- FALSE
  Vbad[same] <- (Vst - outer(mu, mu))[same]
  expect_gt(max(abs(Vbad - Vemp)), 3 * max(abs(V - Vemp)))
})

test_that("admc and adgh recover the truth from ordinal aggregate data", {
  s  <- .ord_setup()
  st <- .ord_study(s)
  for (est in c("admc", "adgh")) {
    ctl <- if (est == "admc")
      admControl(studies = list(s = st), n_sim = 600L, maxeval = 120L, grad = "fd")
    else adghControl(studies = list(s = st), maxeval = 120L)
    fit <- suppressWarnings(suppressMessages(
      nlmixr2est::nlmixr2(.ord_model(), admData(), est = est, control = ctl)))
    expect_s3_class(fit, "admFit")
    expect_true(is.finite(fit$objective))
    ex <- fit$env$admExtra
    expect_equal(exp(ex$struct[["tcl"]]),  5, tolerance = 0.10)
    expect_equal(exp(ex$struct[["tv"]]),  20, tolerance = 0.10)
    expect_equal(ex$omega[1, 1],        0.36, tolerance = 0.20)
  }
})

test_that("adfo fits an ordinal endpoint, with the expected FO bias", {
  # adfo supports joint units, so ordinal runs -- but FO linearises the model at
  # eta = 0, and a category probability is a LOGISTIC function of the prediction.
  # With omega = 0.36 that linearisation is poor, so the point estimates are
  # visibly biased (measured: cl ~ 4.35 against a truth of 5) while admc and adgh
  # -- which integrate over eta properly -- land within a few percent. That is a
  # property of the FO approximation on a saturating endpoint, not a defect, so
  # this test pins down that the fit RUNS and stays in a sane range rather than
  # pretending FO is accurate here. Prefer admc/adgh for ordinal endpoints.
  s   <- .ord_setup()
  fit <- suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
    .ord_model(), admData(), est = "adfo",
    control = adfoControl(studies = list(s = .ord_study(s)), maxeval = 120L))))
  expect_s3_class(fit, "admFit")
  expect_true(is.finite(fit$objective))
  ex <- fit$env$admExtra
  expect_equal(exp(ex$struct[["tcl"]]),  5, tolerance = 0.35)
  expect_equal(exp(ex$struct[["tv"]]),  20, tolerance = 0.25)
})

test_that("ordinal refuses a diagonal V, a non-joint study, and adirmc", {
  s <- .ord_setup()
  nt <- s$nt

  bad <- .ord_study(s); bad$V <- NULL
  bad$observations$c1$V <- diag(s$V[seq_len(nt), seq_len(nt)])
  bad$observations$c2$V <- diag(s$V[nt + seq_len(nt), nt + seq_len(nt)])
  expect_error(
    suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
      .ord_model(), admData(), est = "admc",
      control = admControl(studies = list(s = bad), n_sim = 100L, maxeval = 2L)))),
    "full observed covariance")

  flat <- list(E = s$E[seq_len(nt)], V = s$V[seq_len(nt), seq_len(nt)],
               n = s$N, times = s$times, ev = s$ev, output = "p1")
  expect_error(
    suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
      .ord_model(), admData(), est = "admc",
      control = admControl(studies = list(s = flat), n_sim = 100L, maxeval = 2L)))),
    "one observation block per category")

  expect_error(
    suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
      .ord_model(), admData(), est = "adirmc",
      control = adirmcControl(studies = list(s = .ord_study(s)), maxeval = 2L)))),
    # adirmc refuses on the residual FORM before it ever reaches the joint check:
    # its inner kernel builds the importance-weighted mean internally, so an
    # ordinal residual cannot be scored on that path at all.
    "supports only add/prop/pow/combined")
})
