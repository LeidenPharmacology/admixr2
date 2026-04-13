skip_if_not_installed("rxode2")
skip_if_not_installed("ggplot2")
skip_on_cran()

# Setup in helper-integration.R. .int_plot_setup() reuses the cached
# .int_grad_setup() result: real rxMod + real iniDf + true parameters.
# The "mean" panel runs a genuine rxSolve; "nll"/"par" traces are
# representative values that exercise the back-transform and display-name paths.

.pdf_plot_int <- function(code) {
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  on.exit({ grDevices::dev.off(); unlink(f) }, add = TRUE)
  force(code)
}

# ---- Basic structure ---------------------------------------------------------

test_that("plot.admFit real rxMod: returns named list", {
  env <- .int_plot_setup()
  out <- .pdf_plot_int(plot(env$fit, which = c("mean", "nll", "par"), n_sim = 50L))
  expect_type(out, "list")
  expect_gt(length(out), 0L)
})

# ---- Mean panel: real simulation --------------------------------------------

test_that("plot.admFit real rxMod: mean panel produced for the study", {
  env <- .int_plot_setup()
  out <- .pdf_plot_int(plot(env$fit, which = "mean", n_sim = 50L))
  expect_true(any(startsWith(names(out), "mean_")))
})

test_that("plot.admFit real rxMod: mean panel is gg or list of gg objects", {
  env <- .int_plot_setup()
  out <- .pdf_plot_int(plot(env$fit, which = "mean", n_sim = 50L))
  key <- grep("^mean_", names(out), value = TRUE)[1]
  p   <- out[[key]]
  expect_true(
    inherits(p, "gg") ||
      (is.list(p) && length(p) > 0 && all(vapply(p, inherits, logical(1), "gg")))
  )
})

test_that("plot.admFit real rxMod: mean panel produces finite predictions", {
  env  <- .int_plot_setup()
  out  <- .pdf_plot_int(plot(env$fit, which = "mean", n_sim = 50L))
  key  <- grep("^mean_", names(out), value = TRUE)[1]
  p    <- out[[key]]
  # pred_mean column of the predicted sub-panel should be finite
  pred_panel <- if (inherits(p, "gg")) p else p$pred
  expect_true(all(is.finite(pred_panel$data$pred_mean)))
})

# ---- NLL and parameter traces with real back-transform ----------------------

test_that("plot.admFit real rxMod: nll_trace is a ggplot object", {
  env <- .int_plot_setup()
  out <- .pdf_plot_int(plot(env$fit, which = "nll"))
  expect_s3_class(out$nll_trace, "gg")
})

test_that("plot.admFit real rxMod: par_trace uses iniDf-driven display names", {
  env <- .int_plot_setup()
  out <- .pdf_plot_int(plot(env$fit, which = "par"))
  params <- unique(as.character(out$par_trace$data$param))
  # Real iniDf → omega diagonal shown as V(eta.x)
  expect_true(any(startsWith(params, "V(")))
})

test_that("plot.admFit real rxMod: all par_trace values finite", {
  env <- .int_plot_setup()
  out <- .pdf_plot_int(plot(env$fit, which = "par"))
  expect_true(all(is.finite(out$par_trace$data$value)))
})
