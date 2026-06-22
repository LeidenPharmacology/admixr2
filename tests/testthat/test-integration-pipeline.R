# Tier 2 integration: end-to-end estimator pipeline.
# Exercises the full nlmixr2(model, admData(), est = ...) entry points
# (nlmixr2Est.admc / .adfo / .adirmc) that the other integration files
# deliberately avoid by calling internals directly. A minimal end-to-end check:
# a valid admFit is returned with a finite objective and near-truth estimates --
# not a convergence test. All fits are built once in .int_pipeline_setup().

# Loose band around the data-generating truth. Init values equal truth, so a
# short fit cannot drift far; a wide band keeps the test about "did the pipeline
# return something sensible", not "did it converge".
.expect_sensible <- function(fit, tcl_true, tv_true) {
  th <- fit$theta
  expect_true(all(is.finite(th)))
  expect_true(is.finite(fit$objective))
  expect_gt(fit$objective, 0)
  expect_equal(unname(th[["tcl"]]), tcl_true, tolerance = 0.5)
  expect_equal(unname(th[["tv"]]),  tv_true,  tolerance = 0.5)
  om <- fit$omega
  expect_true(is.matrix(om) && all(is.finite(om)))
}

# ---- admc --------------------------------------------------------------------

test_that("admc pipeline: returns a finite admFit with method 'admc'", {
  env <- .int_pipeline_setup()
  fit <- env$fit_admc
  expect_s3_class(fit, "admFit")
  expect_identical(fit$env$method, "admc")
  .expect_sensible(fit, env$tcl_true, env$tv_true)
})

test_that("admc pipeline: objective statistics (logLik/AIC/BIC) are finite", {
  env <- .int_pipeline_setup()
  fit <- env$fit_admc
  expect_true(is.finite(as.numeric(logLik(fit))))
  expect_true(is.finite(AIC(fit)))
  expect_true(is.finite(BIC(fit)))
})

# ---- adfo --------------------------------------------------------------------

test_that("adfo pipeline: returns a finite admFit with method 'adfo'", {
  env <- .int_pipeline_setup()
  fit <- env$fit_adfo
  expect_s3_class(fit, "admFit")
  expect_identical(fit$env$method, "adfo")
  .expect_sensible(fit, env$tcl_true, env$tv_true)
})

# ---- adirmc (exercises .adirmcPhaseLoop end-to-end) --------------------------

test_that("adirmc pipeline: returns a finite admFit with method 'adirmc'", {
  env <- .int_pipeline_setup()
  fit <- env$fit_adirmc
  expect_s3_class(fit, "admFit")
  expect_identical(fit$env$method, "adirmc")
  .expect_sensible(fit, env$tcl_true, env$tv_true)
})

# ---- adgh --------------------------------------------------------------------

test_that("adgh pipeline: returns a finite admFit with method 'adgh'", {
  env <- .int_pipeline_setup()
  fit <- env$fit_adgh
  expect_s3_class(fit, "admFit")
  expect_identical(fit$env$method, "adgh")
  .expect_sensible(fit, env$tcl_true, env$tv_true)
})

# ---- multi-restart (.admRunRestarts best-result selection) -------------------

test_that("admc multi-restart: records one trace per restart", {
  env <- .int_pipeline_setup()
  fit <- env$fit_restart
  expect_identical(fit$env$method, "admc")
  expect_length(fit$env$admExtra$all_traces, 2L)
  expect_true(is.finite(fit$objective))
})

test_that("admc multi-restart: reported objective matches the best trace", {
  env <- .int_pipeline_setup()
  fit <- env$fit_restart
  # Best-result selection: the returned objective is the minimum over restarts.
  best_per_restart <- vapply(fit$env$admExtra$all_traces,
                             function(tr) min(tr$nll_trace), numeric(1))
  expect_equal(fit$objective, min(best_per_restart), tolerance = 1e-6)
})

# ---- multi-study -------------------------------------------------------------

test_that("admc multi-study: fit carries both studies and a finite objective", {
  env <- .int_pipeline_setup()
  fit <- env$fit_multistudy
  expect_length(fit$env$studies, 2L)
  expect_identical(names(fit$env$studies), c("s1", "s2"))
  .expect_sensible(fit, env$tcl_true, env$tv_true)
})

# ---- covMethod = "r" (in-pipeline .admCalcCov) -------------------------------

test_that("admc covMethod='r': pipeline attaches a covariance matrix", {
  env <- .int_pipeline_setup()
  fit <- env$fit_cov
  expect_identical(fit$env$covMethod, "r")
  expect_false(is.null(fit$cov))
  expect_true(is.matrix(fit$cov))
  expect_equal(nrow(fit$cov), ncol(fit$cov))
  expect_true(all(is.finite(fit$cov)))
  # Symmetric positive-definite Hessian-derived covariance.
  expect_equal(fit$cov, t(fit$cov), tolerance = 1e-8)
  expect_true(all(eigen(fit$cov, symmetric = TRUE, only.values = TRUE)$values > 0))
})

test_that("admc covMethod='r': struct thetas get finite standard errors", {
  env <- .int_pipeline_setup()
  fit <- env$fit_cov
  pf <- fit$parFixedDf
  se_struct <- pf[c("tcl", "tv"), "SE"]
  expect_true(all(is.finite(se_struct)))
  expect_true(all(se_struct > 0))
})

# ---- plot / print on a real fit object ---------------------------------------

test_that("plot.admFit on a real fit: nll/par trace panels are ggplots", {
  skip_if_not_installed("ggplot2")
  env <- .int_pipeline_setup()
  # which = c("nll","par") skips MC simulation -> fast, no rxSolve.
  p <- plot(env$fit_restart, which = c("nll", "par"))
  expect_true(is.list(p))
  expect_true(all(c("nll_trace", "par_trace") %in% names(p)))
  expect_s3_class(p$nll_trace, "ggplot")
  expect_s3_class(p$par_trace, "ggplot")
})

test_that("print.admFit on a real fit: runs without error", {
  env <- .int_pipeline_setup()
  expect_output(print(env$fit_admc))
})

# ---- restart / covariance / multi-study for the non-MC estimators ------------
# Covers the per-estimator restart workers (.adghRestartWorker /
# .adfoRestartWorker / .adirmcRestartWorker), the in-pipeline covariance
# branches (.adghCalcCov / .adfoCalcCov / .adirmcCalcCov) and adirmc multi-study
# -- paths previously only exercised for admc at the pipeline level.

# A deterministic-surface (adgh/adfo) Hessian covariance should be SPD; for the
# MC-based adirmc it may legitimately fail to be PD on a short fit, so callers
# pass require_pd accordingly.
.expect_cov_ok <- function(fit, require_pd) {
  expect_s3_class(fit, "admFit")
  expect_true(is.finite(fit$objective))
  if (is.null(fit$cov)) {
    if (require_pd) fail("covariance was not computed")
    succeed()
    return(invisible())
  }
  expect_identical(fit$env$covMethod, "r")
  expect_true(is.matrix(fit$cov) && all(is.finite(fit$cov)))
  expect_equal(nrow(fit$cov), ncol(fit$cov))
  expect_equal(fit$cov, t(fit$cov), tolerance = 1e-8)
  if (require_pd)
    expect_true(all(eigen(fit$cov, symmetric = TRUE, only.values = TRUE)$values > 0))
}

# -- adgh ----------------------------------------------------------------------

test_that("adgh covMethod='r': pipeline computes an SPD covariance", {
  env <- .int_pipeline_setup()
  .expect_cov_ok(env$fit_adgh_cov, require_pd = TRUE)
})

test_that("adgh multi-restart: records one trace per restart, finite objective", {
  env <- .int_pipeline_setup()
  fit <- env$fit_adgh_restart
  expect_identical(fit$env$method, "adgh")
  expect_length(fit$env$admExtra$all_traces, 2L)
  expect_true(is.finite(fit$objective))
})

# -- adfo ----------------------------------------------------------------------

test_that("adfo multi-restart: records one trace per restart, finite objective", {
  env <- .int_pipeline_setup()
  fit <- env$fit_adfo_restart
  expect_identical(fit$env$method, "adfo")
  expect_length(fit$env$admExtra$all_traces, 2L)
  expect_true(is.finite(fit$objective))
})

test_that("adfo covMethod='r': pipeline computes an SPD covariance", {
  env <- .int_pipeline_setup()
  .expect_cov_ok(env$fit_adfo_cov, require_pd = TRUE)
})

# -- adirmc --------------------------------------------------------------------

test_that("adirmc multi-restart: records one trace per restart, finite objective", {
  env <- .int_pipeline_setup()
  fit <- env$fit_adirmc_restart
  expect_identical(fit$env$method, "adirmc")
  # adirmc stores traces in adirmcExtra (admc/adgh use admExtra).
  expect_length(fit$env$adirmcExtra$all_traces, 2L)
  expect_true(is.finite(fit$objective))
})

test_that("adirmc covMethod='r': pipeline attaches a finite covariance", {
  env <- .int_pipeline_setup()
  .expect_cov_ok(env$fit_adirmc_cov, require_pd = FALSE)
})

test_that("adirmc multi-study: fit carries both studies and a finite objective", {
  env <- .int_pipeline_setup()
  fit <- env$fit_adirmc_multistudy
  expect_length(fit$env$studies, 2L)
  expect_identical(names(fit$env$studies), c("s1", "s2"))
  .expect_sensible(fit, env$tcl_true, env$tv_true)
})
