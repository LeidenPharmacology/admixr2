skip_if_not_installed("ggplot2")

# ---- helpers -----------------------------------------------------------------

# Minimal admFit-like object for trace-panel tests (no rxode2, no real fit).
.make_mock_fit <- function(n_restarts = 1L, n_iter = 5L,
                           par_names = c("tcl", "log_omega_cl")) {
  all_traces <- lapply(seq_len(n_restarts), function(r) {
    list(
      restart_id = r,
      nll_trace  = seq(100 + r, 80 + r, length.out = n_iter),
      par_trace  = matrix(seq_len(n_iter * length(par_names)),
                          nrow = n_iter, ncol = length(par_names))
    )
  })
  env <- new.env(parent = emptyenv())
  env$admExtra <- list(
    all_traces     = all_traces,
    par_names      = par_names,
    studies        = list(),
    n_sim          = 100L,
    omega          = diag(0.09, 1L),
    L              = matrix(sqrt(0.09), 1L, 1L),
    sigma_var      = c(prop.sd = 0.04),
    sigma_is_prop  = TRUE,
    sigma_is_lnorm = FALSE,
    eta_col_names  = "eta.cl",
    struct         = c(tcl = log(5)),
    sampling       = "sobol"
  )
  env$ui <- list(iniDf = NULL, simulationModel = NULL)
  structure(list(env = env), class = c("admFit", "list"))
}

# Open a temporary PDF device so ggplot2 rendering side effects land somewhere.
.pdf_wrap <- function(code) {
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  on.exit({ grDevices::dev.off(); unlink(f) }, add = TRUE)
  force(code)
}

# ---- nll_trace panel ---------------------------------------------------------

test_that("plot.admFit nll: returns list with nll_trace ggplot object", {
  fit <- .make_mock_fit()
  out <- .pdf_wrap(plot(fit, which = "nll"))
  expect_type(out, "list")
  expect_true("nll_trace" %in% names(out))
  expect_s3_class(out$nll_trace, "gg")
})

test_that("plot.admFit nll: ggplot data has eval, nll, restart columns", {
  fit <- .make_mock_fit(n_iter = 7L)
  out <- .pdf_wrap(plot(fit, which = "nll"))
  d <- out$nll_trace$data
  expect_named(d, c("eval", "nll", "restart"), ignore.order = TRUE)
  expect_equal(nrow(d), 7L)
})

test_that("plot.admFit nll: multi-restart produces correct colour factor levels", {
  fit <- .make_mock_fit(n_restarts = 3L, n_iter = 4L)
  out <- .pdf_wrap(plot(fit, which = "nll"))
  expect_length(levels(out$nll_trace$data$restart), 3L)
})

test_that("plot.admFit nll: empty all_traces returns list without nll_trace", {
  fit <- .make_mock_fit()
  fit$env$admExtra$all_traces <- list()
  expect_no_error(out <- .pdf_wrap(plot(fit, which = "nll")))
  expect_false("nll_trace" %in% names(out))
})

test_that("plot.admFit nll: NULL all_traces returns list without nll_trace", {
  fit <- .make_mock_fit()
  fit$env$admExtra$all_traces <- NULL
  expect_no_error(out <- .pdf_wrap(plot(fit, which = "nll")))
  expect_false("nll_trace" %in% names(out))
})

# ---- par_trace panel ---------------------------------------------------------

test_that("plot.admFit par: returns list with par_trace ggplot object", {
  fit <- .make_mock_fit()
  out <- .pdf_wrap(plot(fit, which = "par"))
  expect_true("par_trace" %in% names(out))
  expect_s3_class(out$par_trace, "gg")
})

test_that("plot.admFit par: ggplot data covers all par_names", {
  par_names <- c("tcl", "log_omega_cl", "log_sigma_prop")
  fit <- .make_mock_fit(par_names = par_names)
  out <- .pdf_wrap(plot(fit, which = "par"))
  params_in_data <- unique(as.character(out$par_trace$data$param))
  expect_setequal(params_in_data, par_names)
})

test_that("plot.admFit par: multi-restart produces correct colour factor levels", {
  fit <- .make_mock_fit(n_restarts = 2L)
  out <- .pdf_wrap(plot(fit, which = "par"))
  expect_length(levels(out$par_trace$data$restart), 2L)
})

test_that("plot.admFit par: ggplot data has iter, restart, param, value columns", {
  fit <- .make_mock_fit(n_iter = 6L, par_names = c("tcl", "log_omega_cl"))
  out <- .pdf_wrap(plot(fit, which = "par"))
  d <- out$par_trace$data
  expect_named(d, c("iter", "restart", "param", "value"), ignore.order = TRUE)
  # 2 params * 6 iters * 1 restart
  expect_equal(nrow(d), 12L)
})

# ---- combined panels ---------------------------------------------------------

test_that("plot.admFit c('nll','par'): returns both keys", {
  fit <- .make_mock_fit()
  out <- .pdf_wrap(plot(fit, which = c("nll", "par")))
  expect_true("nll_trace" %in% names(out))
  expect_true("par_trace"  %in% names(out))
})

test_that("plot.admFit c('nll','par'): no mean_ or cov_ keys produced", {
  fit <- .make_mock_fit()
  out <- .pdf_wrap(plot(fit, which = c("nll", "par")))
  expect_false(any(startsWith(names(out), "mean_")))
  expect_false(any(startsWith(names(out), "cov_")))
})

# ---- graceful degradation when simulation model unavailable ------------------

test_that("plot.admFit mean: warns and returns no mean_ keys when rxMod NULL", {
  fit <- .make_mock_fit()
  fit$env$admExtra$studies <- list(s1 = list(
    E = c(1, 2), V = diag(2), n = 10L, times = c(1, 2)
  ))
  out <- suppressWarnings(.pdf_wrap(plot(fit, which = "mean")))
  expect_warning(
    .pdf_wrap(plot(fit, which = "mean")),
    "could not retrieve simulation model"
  )
  expect_type(out, "list")
  expect_false("mean_s1" %in% names(out))
})

test_that("plot.admFit cov: warns and returns no cov_ keys when rxMod NULL", {
  fit <- .make_mock_fit()
  fit$env$admExtra$studies <- list(s1 = list(
    E = c(1, 2), V = diag(2), n = 10L, times = c(1, 2)
  ))
  out <- suppressWarnings(.pdf_wrap(plot(fit, which = "cov")))
  expect_warning(
    .pdf_wrap(plot(fit, which = "cov")),
    "could not retrieve simulation model"
  )
  expect_type(out, "list")
  expect_false("cov_s1" %in% names(out))
})

# ---- head crash guards -------------------------------------------------------

test_that("head.admFit returns data.frame without error", {
  fit <- .make_mock_fit()
  out <- suppressWarnings(head(fit, n = 3L))
  expect_s3_class(out, "data.frame")
})

test_that("head.paged_df returns empty data.frame for environment input", {
  e <- new.env(parent = emptyenv())
  class(e) <- "paged_df"
  out <- head(e)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 0L)
})
