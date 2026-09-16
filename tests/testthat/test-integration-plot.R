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

.int_plot_result <- NULL
.int_plot_all <- function() {
  if (!is.null(.int_plot_result)) return(.int_plot_result)
  env <- .int_plot_setup()
  .int_plot_result <<- .pdf_plot_int(
    plot(env$fit, which = c("mean", "cov", "nll", "par"), n_sim = 50L)
  )
  .int_plot_result
}

# ---- Basic structure ---------------------------------------------------------

test_that("plot.admFit real rxMod: returns named list", {
  out <- .int_plot_all()
  expect_type(out, "list")
  expect_gt(length(out), 0L)
})

# ---- Mean panel: real simulation --------------------------------------------

test_that("plot.admFit real rxMod: mean panel produced for the study", {
  out <- .int_plot_all()
  expect_true(any(startsWith(names(out), "mean_")))
})

# Combined grid key = mean_<study> with no _obs/_pred/_resid/_std_resid suffix.
.mean_combined_key <- function(nms)
  grep("_(obs|pred|resid)$", grep("^mean_", nms, value = TRUE),
       invert = TRUE, value = TRUE)[1]

test_that("plot.admFit real rxMod: mean panel is gg or list of gg objects", {
  out <- .int_plot_all()
  p   <- out[[.mean_combined_key(names(out))]]
  expect_true(
    inherits(p, "gg") ||
      (is.list(p) && length(p) > 0 && all(vapply(p, inherits, logical(1), "gg")))
  )
})

test_that("plot.admFit real rxMod: mean panel produces finite predictions", {
  out <- .int_plot_all()
  # Predicted sub-panel is now individually extractable.
  pred_panel <- out[[grep("^mean_.*_pred$", names(out), value = TRUE)[1]]]
  expect_s3_class(pred_panel, "gg")
  expect_true(all(is.finite(pred_panel$data$pred_mean)))
})

test_that("plot.admFit real rxMod: mean sub-panels are individually extractable", {
  out   <- .int_plot_all()
  study <- sub("^mean_", "", .mean_combined_key(names(out)))
  for (suf in c("obs", "pred", "resid", "std_resid")) {
    key <- paste0("mean_", study, "_", suf)
    expect_true(key %in% names(out))
    expect_s3_class(out[[key]], "gg")
  }
})

test_that("plot.admFit real rxMod: cov sub-panels are individually extractable", {
  out   <- .int_plot_all()
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
  out   <- .int_plot_all()
  pred_panel <- out[[grep("^mean_.*_pred$", names(out), value = TRUE)[1]]]
  expect_equal(as.numeric(pred_panel$data$pred_mean),
               as.numeric(agg[[study]]$pred$E))
})

# ---- NLL and parameter traces with real back-transform ----------------------

test_that("plot.admFit real rxMod: nll_trace is a ggplot object", {
  out <- .int_plot_all()
  expect_s3_class(out$nll_trace, "gg")
})

test_that("plot.admFit real rxMod: par_trace uses iniDf-driven display names", {
  out <- .int_plot_all()
  params <- unique(as.character(out$par_trace$data$param))
  # Real iniDf → omega diagonal shown as V(eta.x)
  expect_true(any(startsWith(params, "V(")))
})

test_that("plot.admFit real rxMod: all par_trace values finite", {
  out <- .int_plot_all()
  expect_true(all(is.finite(out$par_trace$data$value)))
})

test_that("plot.admFit default which: a fit with no covariates is unchanged", {
  # "covariate" joined the default `which`, and it must be a no-op for a fit
  # whose studies declare none -- otherwise every existing caller of plot(fit)
  # gains an empty panel.
  env <- .int_plot_setup()
  out <- .pdf_plot_int(plot(env$fit, n_sim = 50L))
  expect_true(any(startsWith(names(out), "mean_")))
  expect_true(any(startsWith(names(out), "cov_")))
  expect_false(any(startsWith(names(out), "covariate_")))
})

# ---- covariate fits ---------------------------------------------------------

# A banded source becomes several studies named `<source>_s1`, `_s2`, ..., and a
# marginalised covariate has to reach the diagnostic draw as per-subject values:
# held at its mean instead, the predicted V loses the covariate spread the
# observed V in the same panel still carries.

.int_cov_plot_result <- NULL
.int_cov_plot <- function() {
  if (!is.null(.int_cov_plot_result)) return(.int_cov_plot_result)
  skip_if_not_installed("nlmixr2est")
  skip_if_not_installed("patchwork")
  set.seed(1)
  pop <- data.frame(WT  = rlnorm(120L, log(70), 0.2),
                    SEX = rbinom(120L, 1L, 0.5))
  mfn <- function() {
    ini({
      tcl <- log(5); tv <- log(50); bwt <- 0.6; bsex <- 0.2
      add.err <- 0.1
      eta.cl ~ 0.1
    })
    model({
      cl <- exp(tcl + eta.cl) * (WT / 70)^bwt * exp(bsex * SEX)
      v  <- exp(tv)
      cp <- linCmt()
      cp ~ add(add.err)
    })
  }
  st <- admStudies(A = admStudy(model = mfn, population = pop, dose = 100,
                                times = c(0.5, 1, 2, 4, 8),
                                stratify = "SEX", label = "A"))
  nlmixr2 <- nlmixr2est::nlmixr2
  fit <- suppressMessages(suppressWarnings(
    nlmixr2(mfn, admData(), est = "adgh",
            control = adghControl(studies = st, print = 0L,
                                  n_restart = 1L, maxeval = 3L))))
  .int_cov_plot_result <<- list(
    fit  = fit,
    out  = .pdf_plot_int(plot(fit, which = c("mean", "cov"), n_sim = 200L)))
  .int_cov_plot_result
}

test_that("plot.admFit covariates: one panel set per stratum", {
  skip_if_not_installed("nlmixr2")
  env <- .int_cov_plot()
  expect_setequal(grep("^mean_.*_(obs|pred|resid)$", grep("^mean_", names(env$out), value = TRUE),
                       invert = TRUE, value = TRUE),
                  c("mean_A_s1", "mean_A_s2"))
})

test_that("plot.admFit covariates: stratum panels are titled by covariate value", {
  skip_if_not_installed("nlmixr2")
  skip_if_not_installed("patchwork")
  env    <- .int_cov_plot()
  titles <- vapply(c("mean_A_s1", "mean_A_s2"),
                   function(k) env$out[[k]]$patches$annotation$title, character(1))
  # One stratum per sex level, each naming the level it was conditioned at --
  # not the bare `_s1`/`_s2` index, which says nothing about which is which.
  expect_setequal(unname(titles),
                  c("Study 'A_s1' [SEX = 0] -- Mean diagnostics",
                    "Study 'A_s2' [SEX = 1] -- Mean diagnostics"))
})

test_that("plot.admFit covariates: both covariate panels are produced", {
  skip_if_not_installed("nlmixr2est")
  env <- .int_cov_plot()
  out <- .pdf_plot_int(plot(env$fit, which = "covariate", n_sim = 200L))
  expect_s3_class(out$covariate_effect, "gg")
  expect_s3_class(out$covariate_resid,  "gg")
  # WT is swept; SEX is banded, so its facet is the two levels and nothing
  # between them. ONE estimated-effect line per facet -- other covariates sit
  # at the pooled centre, and the line is the thing the sources are compared
  # against.
  d  <- out$covariate_effect$data
  wt <- d[d$cov == "WT" & d$param == "cl", , drop = FALSE]
  expect_gt(nrow(wt), 0L)
  expect_true(all(diff(wt$y) > 0))
  expect_equal(sort(unique(d$x[d$cov == "SEX"])), c(0, 1))
})

# The claim the residual panel's subtitle makes -- "a SLOPE is a mis-specified
# covariate form" -- checked against a fit that IS mis-specified, rather than
# against its own output. Three cohorts at three renal medians are generated
# from a model with a real renal effect; the analysis model then has that effect
# restricted to zero. Nothing else differs between the two fits.
.int_cov_mis_result <- NULL
.int_cov_mis <- function() {
  if (!is.null(.int_cov_mis_result)) return(.int_cov_mis_result)
  skip_if_not_installed("nlmixr2est")
  nlmixr2 <- nlmixr2est::nlmixr2
  set.seed(11)
  draw <- function(n, crclm) data.frame(
    WT = rlnorm(n, log(76), 0.198), CRCL = rlnorm(n, log(crclm), 0.05))
  true_fn <- function() {
    ini({ tcl <- log(5); tv <- log(50); bcrcl <- 0.6; add.err <- 0.08
          eta.cl ~ 0.05 })
    model({ cl <- exp(tcl + eta.cl) * (WT/70)^0.75 * (CRCL/90)^bcrcl
            v  <- exp(tv) * (WT/70); cp <- linCmt(); cp ~ add(add.err) })
  }
  # The restriction written the way it now reads: the term is simply GONE, so
  # the model never mentions CRCL. The studies still declare it -- they describe
  # who was enrolled -- so it comes off the design, and the panel has to keep
  # plotting against it anyway. That is the whole case this diagnostic is for.
  null_fn <- function() {
    ini({ tcl <- log(5); tv <- log(50); add.err <- 0.08; eta.cl ~ 0.05 })
    model({ cl <- exp(tcl + eta.cl) * (WT/70)^0.75
            v  <- exp(tv) * (WT/70); cp <- linCmt(); cp ~ add(add.err) })
  }
  mk <- function(co) admStudy(model = true_fn, population = co, dose = 200,
                              times = c(0.5, 1, 2, 4, 8, 12, 24))
  st <- admStudies(normal   = mk(draw(260L, 95)),
                   mild     = mk(draw(210L, 62)),
                   moderate = mk(draw(180L, 38)))
  # ONE `studies` object for both models -- what dropping an unread covariate
  # instead of refusing it is for.
  slope <- function(fn) {
    fit <- suppressMessages(suppressWarnings(
      nlmixr2(fn, admData(), est = "adgh",
              control = adghControl(studies = st, print = 0L, n_restart = 1L,
                                    maxeval = 40L))))
    d <- .pdf_plot_int(plot(fit, which = "covariate",
                            n_sim = 300L))$covariate_resid$data
    d[d$cov == "CRCL", ]
  }
  .int_cov_mis_result <<- list(ok = slope(true_fn), bad = slope(null_fn))
  .int_cov_mis_result
}

test_that("plot.admFit covariates: a dropped covariate is still plotted against", {
  skip_on_cran()
  env <- .int_cov_mis()
  # The null model never reads CRCL, so it is off the design -- and it is still
  # the covariate the analyst needs the residual plotted against. It also keeps
  # its declared SPREAD: the source marginalised over a distribution, and the
  # `cov` value the drop leaves behind is a single number that would have been
  # mislabelled as a conditioned one.
  expect_equal(nrow(env$bad), 3L)
  expect_true(all(env$bad$kind == "marginal"))
  expect_true(all(env$bad$xhi > env$bad$xlo))
})

test_that("plot.admFit covariates: a correct covariate form leaves no trend", {
  skip_on_cran()
  env <- .int_cov_mis()
  expect_lt(abs(unname(coef(stats::lm(z ~ x, env$ok))[2L])), 0.01)
  expect_lt(diff(range(env$ok$z)), 1)
})

test_that("plot.admFit covariates: a dropped covariate effect shows as a slope", {
  skip_on_cran()
  env <- .int_cov_mis()
  b_ok  <- unname(coef(stats::lm(z ~ x, env$ok))[2L])
  b_bad <- unname(coef(stats::lm(z ~ x, env$bad))[2L])
  # Ordered in the covariate, well past the +/-1.96 the panel draws, and orders
  # of magnitude steeper than the correctly specified fit on the same data.
  expect_gt(abs(b_bad), 50 * abs(b_ok))
  expect_gt(max(abs(env$bad$z)), 1.96)
  expect_equal(order(env$bad$x), order(-env$bad$z))
})

test_that("plot.admFit covariates: predicted V carries the marginalised spread", {
  skip_if_not_installed("nlmixr2")
  env   <- .int_cov_plot()
  extra <- env$fit$env$admExtra
  ui    <- env$fit$env$ui
  with_cov <- admixr2:::.admAggData(extra, ui, n_sim = 200L, seed = 1L, warn = FALSE)
  # The same fit with WT pinned at its mean instead of integrated over: this is
  # what the panels showed when the diagnostic draw kept the estimator's own
  # covariate reduction rather than forcing per-row values.
  flat  <- extra
  flat$studies <- lapply(extra$studies, function(s) { s$cov_dist <- NULL; s })
  no_cov <- admixr2:::.admAggData(flat, ui, n_sim = 200L, seed = 1L, warn = FALSE)
  nm <- names(extra$studies)[1]
  expect_gt(max(abs(diag(with_cov[[nm]]$pred$V) - diag(no_cov[[nm]]$pred$V))),
            1e-6)
})
