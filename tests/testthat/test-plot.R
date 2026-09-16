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

# ---- conditioned-study panel titles ------------------------------------------

# `stratify` names the strata `<source>_s1`, `_s2`, ... and the index alone does
# not say which covariate slice a panel is. The label is read off the stratum's
# own point specs, so it has to skip the covariates the study marginalises --
# those carry a filled-in mean in `cov` as well.

test_that(".admStudyCovLabel names only the point-conditioned covariates", {
  s <- list(cov      = list(SEX = 0, WT = 72.7),
            cov_dist = list(SEX = list(.point = TRUE, quantile = function(u) 0),
                            WT  = list(meanlog = log(70), sdlog = 0.2)))
  expect_equal(.admStudyCovLabel(s), "SEX = 0")
})

test_that(".admStudyCovLabel is empty for a fully marginalised study", {
  s <- list(cov      = list(WT = 72.7),
            cov_dist = list(WT = list(meanlog = log(70), sdlog = 0.2)))
  expect_equal(.admStudyCovLabel(s), "")
})

test_that(".admStudyCovLabel labels a covariate held at a value with no cov_dist", {
  # `at`/`by` conditioning: a value and no distribution to marginalise.
  expect_equal(.admStudyCovLabel(list(cov = list(CRCL = 30))), "CRCL = 30")
})

test_that(".admStudyCovLabel is empty when the study declares no covariates", {
  expect_equal(.admStudyCovLabel(list(E = 1, V = 1)), "")
})

test_that(".admStudyTitle appends the conditioning, and omits it otherwise", {
  s_fix <- list(cov = list(SEX = 1), cov_dist = list(SEX = list(.point = TRUE)))
  expect_equal(.admStudyTitle(s_fix, "A_s2"), "Study 'A_s2' [SEX = 1]")
  expect_equal(.admStudyTitle(list(), "B"),   "Study 'B'")
})

# ---- covariate panel data ----------------------------------------------------

.cov_studies <- function()
  list(
    lo = list(n = 100L, cov = list(WT = 70, SEX = 0),
              cov_dist = list(WT  = list(meanlog = log(70), sdlog = 0.2),
                              SEX = list(.point = TRUE))),
    hi = list(n = 300L, cov = list(WT = 90, SEX = 1),
              cov_dist = list(WT  = list(meanlog = log(90), sdlog = 0.2),
                              SEX = list(.point = TRUE))))

.cov_ui <- function() {
  fn <- function() {
    ini({
      tcl <- log(5); bwt <- 0.75; bsex <- 0.2
      add.err <- 0.1
      eta.cl ~ 0.1
    })
    model({
      cl <- exp(tcl + eta.cl) * (WT / 70)^bwt * exp(bsex * SEX)
      v  <- exp(log(30))
      cp <- linCmt()
      cp ~ add(add.err)
    })
  }
  suppressMessages(rxode2::rxode2(fn))
}

test_that(".admCovStudyQ reads a quantile, a point value, and an absence", {
  st <- .cov_studies()
  expect_equal(.admCovStudyQ(st$lo, "WT", 0.5), 70, tolerance = 1e-8)
  expect_gt(.admCovStudyQ(st$lo, "WT", 0.9), .admCovStudyQ(st$lo, "WT", 0.1))
  # A point spec does not move with the quantile -- it is a conditioned value.
  expect_equal(.admCovStudyQ(st$lo, "SEX", 0.1), 0)
  expect_equal(.admCovStudyQ(st$lo, "SEX", 0.9), 0)
  expect_true(is.na(.admCovStudyQ(st$lo, "AGE", 0.5)))
})

test_that(".admCovPanelCovs keeps covariates the model reads and a study describes", {
  skip_if_not_installed("rxode2")
  expect_setequal(.admCovPanelCovs(.cov_ui(), .cov_studies()), c("WT", "SEX"))
  # Described by no study: nothing to put on an axis.
  expect_equal(.admCovPanelCovs(.cov_ui(), list(a = list(n = 1L))), character(0))
})

test_that(".admCovStudyKind separates a marginalised covariate from a conditioned one", {
  st <- .cov_studies()
  expect_equal(.admCovStudyKind(st$lo, "WT"),  "marginal")
  expect_equal(.admCovStudyKind(st$lo, "SEX"), "conditional")
  # A value with no distribution behind it is conditioned, not marginalised.
  expect_equal(.admCovStudyKind(list(cov = list(AGE = 40)), "AGE"),
               "conditional")
})

test_that(".admCovEffectData marks marginal spread and conditioned points apart", {
  skip_if_not_installed("rxode2")
  st <- .cov_studies()
  mw <- .admCovEffectData(.cov_ui(), "WT", st,
                          list(tcl = log(5), bwt = 0.75, bsex = 0.2))$marks
  expect_true(all(mw$kind == "marginal"))
  # A distribution, drawn as one: median inside a 10th-90th bar inside a
  # 2.5th-97.5th whisker.
  expect_true(all(mw$xlo2 < mw$xlo & mw$xlo < mw$x &
                  mw$x < mw$xhi & mw$xhi < mw$xhi2))

  ms <- .admCovEffectData(.cov_ui(), "SEX", st,
                          list(tcl = log(5), bwt = 0.75, bsex = 0.2))$marks
  expect_true(all(ms$kind == "conditional"))
  # Conditioned at one value: no spread to draw at all.
  expect_true(all(ms$xlo == ms$x & ms$xhi == ms$x & ms$xlo2 == ms$x))
})

test_that(".admCovResidData carries the span each study speaks for", {
  agg <- list(
    lo = list(obs = list(E = c(1, 2)), pred = list(E = c(1.1, 2.1),
                                                   V = diag(c(0.01, 0.04)))),
    hi = list(obs = list(E = c(1, 2)), pred = list(E = c(0.9, 1.9),
                                                   V = diag(c(0.01, 0.04)))))
  st <- .cov_studies()
  wt <- .admCovResidData("WT", st, agg)
  expect_true(all(wt$kind == "marginal") && all(wt$xhi > wt$xlo))
  sx <- .admCovResidData("SEX", st, agg)
  expect_true(all(sx$kind == "conditional") && all(sx$xhi == sx$x))
})

test_that(".admCovPooled weights each study's median by n", {
  # 70 at n = 100 and 90 at n = 300 -> 85, not the unweighted 80.
  expect_equal(.admCovPooled("WT", .cov_studies())$WT, 85, tolerance = 1e-6)
})

test_that(".admMergeCovMarks merges coincident strata under their source name", {
  df <- data.frame(
    cov = "WT", param = "cl", study = c("a_s1", "a_s2", "b_s1"),
    x = c(70, 70, 90), xlo = 60, xhi = 80, y = c(5, 5, 6),
    stringsAsFactors = FALSE)
  out <- .admMergeCovMarks(df)
  expect_equal(nrow(out), 2L)
  expect_setequal(out$study, c("a", "b"))
})

test_that(".admMergeCovMarks keeps strata that genuinely differ apart", {
  # The axis they were banded ON: same source, different x, so no merge.
  df <- data.frame(
    cov = "SEX", param = "cl", study = c("a_s1", "a_s2"),
    x = c(0, 1), xlo = c(0, 1), xhi = c(0, 1), y = c(5, 6),
    stringsAsFactors = FALSE)
  expect_equal(nrow(.admMergeCovMarks(df)), 2L)
})

test_that(".admCovEffectData sweeps a continuous covariate and pads past it", {
  skip_if_not_installed("rxode2")
  d <- .admCovEffectData(.cov_ui(), "WT", .cov_studies(),
                         list(tcl = log(5), bwt = 0.75, bsex = 0.2))
  expect_equal(unique(d$curve$param), "cl")
  expect_equal(nrow(d$curve), 120L)
  # Allometric with a positive exponent: monotone increasing in weight.
  expect_true(all(diff(d$curve$y) > 0))
  # Two shaded regions, one past each end of the range the studies cover.
  expect_equal(nrow(d$shade), 2L)
  expect_lt(min(d$curve$x), min(vapply(.cov_studies(), .admCovStudyQ,
                                       double(1), cv = "WT", u = 0.1)))
})

test_that(".admCovEffectData puts a discrete covariate on its levels", {
  skip_if_not_installed("rxode2")
  d <- .admCovEffectData(.cov_ui(), "SEX", .cov_studies(),
                         list(tcl = log(5), bwt = 0.75, bsex = 0.2))
  # The levels ARE the support: no sweep through SEX = 0.37, no extrapolation
  # region past them.
  expect_equal(sort(unique(d$curve$x)), c(0, 1))
  expect_true(all(d$curve$disc))
  expect_null(d$shade)
})

test_that(".admCovEffectData returns NULL when there is nothing to sweep", {
  skip_if_not_installed("rxode2")
  ui <- .cov_ui()
  # A covariate the model never reads.
  expect_null(.admCovEffectData(ui, "AGE", .cov_studies(), list()))
  # Every study at the same single value: no effect the sources speak to.
  flat <- lapply(.cov_studies(), function(s) {
    s$cov_dist <- list(WT = list(.point = TRUE)); s$cov <- list(WT = 70); s
  })
  expect_null(.admCovEffectData(ui, "WT", flat, list(tcl = log(5), bwt = 0.75)))
})

test_that("plot.admFit covariate panel survives an all-discrete figure", {
  skip_if_not_installed("rxode2")
  # Every covariate banded: no study contributes a shade region, so the rect
  # layer's data is NULL and `nrow(NULL)` is length zero, not FALSE.
  fit <- .make_mock_fit()
  fit$env$ui <- .cov_ui()
  fit$env$admExtra$struct <- list(tcl = log(5), bwt = 0.75, bsex = 0.2)
  fit$env$admExtra$studies <- list(
    a = list(E = c(1, 2), V = diag(2), n = 50L, times = c(1, 2),
             cov = list(WT = 70, SEX = 0),
             cov_dist = list(WT  = list(.point = TRUE),
                             SEX = list(.point = TRUE))),
    b = list(E = c(1, 2), V = diag(2), n = 50L, times = c(1, 2),
             cov = list(WT = 90, SEX = 1),
             cov_dist = list(WT  = list(.point = TRUE),
                             SEX = list(.point = TRUE))))
  out <- suppressWarnings(.pdf_wrap(plot(fit, which = "covariate")))
  expect_s3_class(out$covariate_effect, "gg")
  d <- out$covariate_effect$data
  expect_gt(nrow(d), 0L)
  expect_true(all(d$disc))
  # Two conditioned levels per covariate, and nothing swept between them.
  expect_setequal(d$x[d$cov == "WT"],  c(70, 90))
  expect_setequal(d$x[d$cov == "SEX"], c(0, 1))
})

test_that(".admCovResidData needs two studies to have a contrast", {
  agg <- list(
    lo = list(obs = list(E = c(1, 2)), pred = list(E = c(1.1, 2.1),
                                                   V = diag(c(0.01, 0.04)))),
    hi = list(obs = list(E = c(1, 2)), pred = list(E = c(0.9, 1.9),
                                                   V = diag(c(0.01, 0.04)))))
  st <- .cov_studies()
  out <- .admCovResidData("WT", st, agg)
  expect_equal(nrow(out), 2L)
  # z = (obs - pred)/sqrt(diag(V)/n), averaged over times.
  expect_equal(out$z[out$study == "lo"],
               mean((c(1, 2) - c(1.1, 2.1)) / sqrt(c(0.01, 0.04) / 100)),
               tolerance = 1e-8)
  expect_null(.admCovResidData("WT", st["lo"], agg["lo"]))
})
