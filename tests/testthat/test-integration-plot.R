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

# Combined grid key = mean_<study> with no _obs/_pred/_resid/_std_resid suffix.
.mean_combined_key <- function(nms)
  grep("_(obs|pred|resid)$", grep("^mean_", nms, value = TRUE),
       invert = TRUE, value = TRUE)[1]

test_that("plot.admFit real rxMod: mean panel is gg or list of gg objects", {
  env <- .int_plot_setup()
  out <- .pdf_plot_int(plot(env$fit, which = "mean", n_sim = 50L))
  p   <- out[[.mean_combined_key(names(out))]]
  expect_true(
    inherits(p, "gg") ||
      (is.list(p) && length(p) > 0 && all(vapply(p, inherits, logical(1), "gg")))
  )
})

test_that("plot.admFit real rxMod: mean panel produces finite predictions", {
  env  <- .int_plot_setup()
  out  <- .pdf_plot_int(plot(env$fit, which = "mean", n_sim = 50L))
  # Predicted sub-panel is now individually extractable.
  pred_panel <- out[[grep("^mean_.*_pred$", names(out), value = TRUE)[1]]]
  expect_s3_class(pred_panel, "gg")
  expect_true(all(is.finite(pred_panel$data$pred_mean)))
})

test_that("plot.admFit real rxMod: mean sub-panels are individually extractable", {
  env   <- .int_plot_setup()
  out   <- .pdf_plot_int(plot(env$fit, which = "mean", n_sim = 50L))
  study <- sub("^mean_", "", .mean_combined_key(names(out)))
  for (suf in c("obs", "pred", "resid", "std_resid")) {
    key <- paste0("mean_", study, "_", suf)
    expect_true(key %in% names(out))
    expect_s3_class(out[[key]], "gg")
  }
})

test_that("plot.admFit real rxMod: cov sub-panels are individually extractable", {
  env   <- .int_plot_setup()
  out   <- .pdf_plot_int(plot(env$fit, which = "cov", n_sim = 50L))
  cov_combined <- grep("_(obs|pred|resid)$", grep("^cov_", names(out), value = TRUE),
                       invert = TRUE, value = TRUE)[1]
  study <- sub("^cov_", "", cov_combined)
  for (suf in c("obs", "pred", "resid", "std_resid")) {
    key <- paste0("cov_", study, "_", suf)
    expect_true(key %in% names(out))
    expect_s3_class(out[[key]], "gg")
  }
})

# ---- Aggregate data (.admAggData) -------------------------------------------

test_that("plot.admFit real rxMod: .admAggData returns obs/pred E vector + V matrix", {
  env   <- .int_plot_setup()
  extra <- env$fit$env$adirmcExtra
  ui    <- env$fit$env$ui
  agg   <- admixr2:::.admAggData(extra, ui, n_sim = extra$n_sim, seed = 1L, warn = FALSE)
  expect_type(agg, "list")
  study <- names(extra$studies)[1]
  a     <- agg[[study]]
  expect_false(is.null(a))
  n_t   <- length(extra$studies[[study]]$times)
  # Observed side is taken straight from the study spec.
  expect_equal(as.numeric(a$obs$E), as.numeric(extra$studies[[study]]$E))
  expect_equal(unname(a$obs$V), unname(as.matrix(extra$studies[[study]]$V)))
  # Predicted side: mean vector length n_t, symmetric finite n_t x n_t cov.
  expect_length(a$pred$E, n_t)
  expect_equal(dim(a$pred$V), c(n_t, n_t))
  expect_true(all(is.finite(a$pred$E)))
  expect_true(all(is.finite(a$pred$V)))
  expect_equal(a$pred$V, t(a$pred$V))
})

test_that("plot.admFit real rxMod: aggData pred matches the plot mean panel", {
  env   <- .int_plot_setup()
  extra <- env$fit$env$adirmcExtra
  ui    <- env$fit$env$ui
  agg   <- admixr2:::.admAggData(extra, ui, n_sim = extra$n_sim, seed = 1L, warn = FALSE)
  study <- names(extra$studies)[1]
  out   <- .pdf_plot_int(plot(env$fit, which = "mean", n_sim = extra$n_sim, seed = 1L))
  pred_panel <- out[[grep("^mean_.*_pred$", names(out), value = TRUE)[1]]]
  expect_equal(as.numeric(pred_panel$data$pred_mean),
               as.numeric(agg[[study]]$pred$E))
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
