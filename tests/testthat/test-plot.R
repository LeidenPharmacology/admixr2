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

test_that(".admStudyCovLabel does not call a dropped covariate conditioned", {
  # Dropped from the design, so no `cov_dist` entry -- but the source
  # MARGINALISED over it. Read directly, the panel title would assert
  # "Study 'normal' [CRCL = 90]" for a study that never solved at 90.
  s <- list(cov      = list(WT = 76, CRCL = 90),
            cov_dist = list(WT = list(meanlog = log(76), sdlog = 0.2)),
            .adm_cov_dropped = list(CRCL = list(meanlog = log(90),
                                                sdlog = 0.25)))
  expect_equal(.admStudyCovLabel(s), "")
  expect_equal(.admStudyTitle(s, "normal"), "Study 'normal'")
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

## The ORIGINAL studies for `.cov_studies()`, keyed by SOURCE: each carries
## its own published model, which is what puts a source at its own parameter
## value on the effect panel rather than on the fitted line.
.cov_src <- function()
  stats::setNames(lapply(c("lo", "hi"), function(k) list(
    ui = .cov_ui(),
    range = list(WT = c(55, 105), SEX = c(0, 1)),
    population = list(WT  = list(meanlog = log(if (k == "lo") 70 else 90),
                                 sdlog = 0.2),
                      SEX = list(values = c(0, 1), probs = c(0.5, 0.5))))),
    c("lo", "hi"))

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

test_that(".admCovPanelCovs keeps a covariate the model never reads", {
  skip_if_not_installed("rxode2")
  # Dropped from the design because it cannot move the prediction, and retained
  # on the study. This is the omitted-term case the residual panel exists for,
  # so the facet has to survive the drop.
  st <- .cov_studies()
  st$lo$.adm_cov_dropped <- list(CRCL = list(meanlog = log(90), sdlog = 0.1))
  st$hi$.adm_cov_dropped <- list(CRCL = list(meanlog = log(40), sdlog = 0.1))
  expect_true("CRCL" %in% .admCovPanelCovs(.cov_ui(), st))
})

test_that(".admCovStudySpec and friends read a dropped covariate's distribution", {
  s <- list(cov = list(CRCL = 90),
            .adm_cov_dropped = list(CRCL = list(meanlog = log(90),
                                                sdlog = 0.25)))
  expect_false(is.null(.admCovStudySpec(s, "CRCL")))
  # The source MARGINALISED over it; the `cov` value the drop leaves behind is
  # a single number, and reading that alone would label it as conditioned and
  # throw the spread away.
  expect_equal(.admCovStudyKind(s, "CRCL"), "marginal")
  expect_gt(.admCovStudyQ(s, "CRCL", 0.9), .admCovStudyQ(s, "CRCL", 0.1))
  expect_equal(.admCovStudyQ(s, "CRCL", 0.5), 90, tolerance = 1e-6)
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
                          list(tcl = log(5), bwt = 0.75, bsex = 0.2),
                          .cov_src())$marks
  expect_true(all(mw$kind == "marginal"))
  # A distribution, drawn as one: centre inside a 10th-90th bar inside a
  # 2.5th-97.5th whisker.
  expect_true(all(mw$xlo2 < mw$xlo & mw$xlo < mw$x &
                  mw$x < mw$xhi & mw$xhi < mw$xhi2))
  # And at the SOURCE's own parameter value, not read off the fitted line.
  # `.cov_ui()` at its own ini: cl = 5 * (WT/70)^0.75 * exp(0.2 * SEX), with
  # SEX at the source's own centre -- the prob-weighted mean of its declared
  # 50/50 split, which is the SAME centre the estimated effect holds its other
  # covariates at. Taking the median instead put the source at SEX = 1 and made
  # it look displaced from a fit it agrees with.
  expect_equal(mw$y[mw$study == "lo"],
               5 * (70 / 70)^0.75 * exp(0.2 * 0.5), tolerance = 1e-6)

  ms <- .admCovEffectData(.cov_ui(), "SEX", st,
                          list(tcl = log(5), bwt = 0.75, bsex = 0.2),
                          .cov_src())$marks
  expect_true(all(ms$kind == "conditional"))
  # The solid line is a RANGE, and only a continuous covariate gives one. SEX
  # is levels: a stratum at SEX = 0 covers that level and not the other, so it
  # stays a point. Drawn across the declared 0-1 it would claim both.
  expect_true(all(ms$xlo == ms$x & ms$xhi == ms$x))
  expect_setequal(ms$x, c(0, 1))
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

test_that(".admCovEffectData sweeps a continuous covariate and pads past it", {
  skip_if_not_installed("rxode2")
  d <- .admCovEffectData(.cov_ui(), "WT", .cov_studies(),
                         list(tcl = log(5), bwt = 0.75, bsex = 0.2),
                         .cov_src())
  expect_equal(unique(d$curve$param), "cl")
  # ONE line: the estimated effect. Other covariates sit at the pooled centre.
  expect_equal(nrow(d$curve), 120L)
  # Allometric with a positive exponent: monotone increasing in weight.
  expect_true(all(diff(d$curve$y) > 0))
  # Two shaded regions, one past each end of the range the studies cover.
  expect_equal(nrow(d$shade), 2L)
  expect_lt(min(d$curve$x), min(vapply(.cov_studies(), .admCovStudyQ,
                                       double(1), cv = "WT", u = 0.1)))
})

test_that(".admCovEffectData draws no mark for a source with no model", {
  skip_if_not_installed("rxode2")
  # Falling back to the fitted curve would put the study exactly on the dotted
  # line and read as agreement with a claim it never made.
  d <- .admCovEffectData(.cov_ui(), "WT", .cov_studies(),
                         list(tcl = log(5), bwt = 0.75, bsex = 0.2), NULL)
  expect_false(is.null(d$curve))
  expect_null(d$marks)
})

test_that(".admCovEffectData puts a discrete covariate on its levels", {
  skip_if_not_installed("rxode2")
  d <- .admCovEffectData(.cov_ui(), "SEX", .cov_studies(),
                         list(tcl = log(5), bwt = 0.75, bsex = 0.2),
                         .cov_src())
  # The levels ARE the support: no sweep through SEX = 0.37, no extrapolation
  # region past them.
  expect_equal(sort(unique(d$curve$x)), c(0, 1))
  expect_true(all(d$curve$disc))
  expect_null(d$shade)
})

test_that(".admCovEffectData shades a declared level no study sits at", {
  skip_if_not_installed("rxode2")
  # `values` can declare a group nobody enrolled. The model predicts for it
  # happily, and drawn like the studied levels that prediction looks equally
  # earned. Grey means "extrapolating" on a discrete axis too.
  fn <- function() {
    ini({ tcl <- log(5); bg <- 0.2; add.err <- 0.1; eta.cl ~ 0.1 })
    model({ cl <- exp(tcl + eta.cl) * exp(bg * GRP)
            v <- exp(log(30)); cp <- linCmt(); cp ~ add(add.err) })
  }
  ui <- suppressMessages(rxode2::rxode2(fn))
  mk <- function(at, vals) list(
    n = 50L, cov = list(GRP = at),
    cov_dist = list(GRP = c(list(values = vals), list(.point = TRUE))))

  # Three declared levels, only two ever conditioned at.
  d <- .admCovEffectData(ui, "GRP", list(a = mk(0, c(0, 1, 2)),
                                         b = mk(1, c(0, 1, 2))),
                         list(tcl = log(5), bg = 0.2))
  expect_equal(sort(unique(d$curve$x)), c(0, 1, 2))
  expect_equal(nrow(d$shade), 1L)
  # Reaches level 2 and stops there. Hanging past the outermost level widens
  # the panel limits enough that .admLevelBreaks() falls back to pretty() and
  # ticks the factor at 0.5 and 1.5.
  expect_true(d$shade$xmin < 2)
  expect_equal(d$shade$xmax, 2)

  # Every declared level studied: nothing to warn about.
  d2 <- .admCovEffectData(ui, "GRP", list(a = mk(0, c(0, 1)),
                                          b = mk(1, c(0, 1))),
                          list(tcl = log(5), bg = 0.2))
  expect_null(d2$shade)
})

test_that(".admCovEffectData gives no facet where nothing is estimated", {
  skip_if_not_installed("rxode2")
  # A FIXED allometric exponent: the model reads WT and varies with it, but
  # estimates no coefficient for it. There is no fitted effect to agree or
  # disagree with, so the panel would invite a reader to check an agreement
  # that was never in question.
  fn <- function() {
    ini({ tcl <- log(5); add.err <- 0.1; eta.cl ~ 0.1 })
    model({ cl <- exp(tcl + eta.cl) * (WT / 70)^0.75
            v <- exp(log(30)); cp <- linCmt(); cp ~ add(add.err) })
  }
  ui <- suppressMessages(rxode2::rxode2(fn))
  expect_null(.admCovEffectData(ui, "WT", .cov_studies(), list(tcl = log(5))))
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

test_that("plot.admFit does not simulate for a fit with no covariates", {
  # "covariate" is in the default `which`. A fit declaring no covariates must
  # not pay for a full n_sim simulation to then draw nothing with it.
  fit <- .make_mock_fit()
  fit$env$admExtra$studies <- list(s1 = list(
    E = c(1, 2), V = diag(2), n = 10L, times = c(1, 2)))
  out <- .pdf_wrap(plot(fit, which = "covariate"))
  # No panels, and -- the point of the test -- no simulation attempted on the
  # way to producing none. That warning is the tell that .admAggData() ran.
  expect_length(out, 0L)
  expect_silent(.pdf_wrap(plot(fit, which = "covariate")))
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

test_that(".admCovResidData drops a covariate with no between-study contrast", {
  # Every study at the same value: nothing for the panel to read, and a facet
  # that would stack every point on one x and fit a rank-deficient `lm`.
  agg <- list(
    lo = list(obs = list(E = c(1, 2)), pred = list(E = c(1.1, 2.1),
                                                   V = diag(c(0.01, 0.04)))),
    hi = list(obs = list(E = c(1, 2)), pred = list(E = c(0.9, 1.9),
                                                   V = diag(c(0.01, 0.04)))))
  flat <- lapply(.cov_studies(), function(s) {
    s$cov_dist <- NULL; s$cov <- list(WT = 70); s
  })
  expect_null(.admCovResidData("WT", flat, agg))
})

test_that(".admCovLevels finds levels only where a covariate has them", {
  st <- .cov_studies()
  # Conditioned at a point in every study: the levels are those values.
  expect_equal(.admCovLevels("SEX", st, c(0, 1)), c(0, 1))
  # Marginalised over a lognormal: continuous, no levels.
  expect_null(.admCovLevels("WT", st, c(70, 90)))
  # Too many distinct conditioned values to be a factor.
  flat <- lapply(seq_len(10L), function(i)
    list(cov = list(AGE = i), n = 10L))
  names(flat) <- paste0("s", seq_len(10L))
  expect_null(.admCovLevels("AGE", flat, seq_len(10L)))
})

test_that(".admCovSource recovers the source a stratum came from", {
  expect_equal(.admCovSource(c("normal_s1", "normal_s2", "plain")),
               c("normal", "normal", "plain"))
})

test_that(".admLevelBreaks ticks a discrete panel at its levels only", {
  brk <- .admLevelBreaks(c(0, 1))
  # A panel spanned by the levels is the discrete one.
  expect_equal(brk(c(-0.05, 1.05)), c(0, 1))
  # A continuous panel falls back to pretty(), not to the levels.
  expect_equal(brk(c(35, 98)), pretty(c(35, 98)))
})

test_that(".admCovResidData marks discreteness and pairs a source's strata", {
  agg <- setNames(rep(list(list(obs = list(E = c(1, 2)),
                                pred = list(E = c(1.1, 2.1),
                                            V = diag(c(0.01, 0.04))))), 4L),
                  c("a_s1", "a_s2", "b_s1", "b_s2"))
  mk <- function(sex, wt) list(
    n = 50L, cov = list(SEX = sex, WT = wt),
    cov_dist = list(SEX = list(.point = TRUE),
                    WT  = list(meanlog = log(wt), sdlog = 0.2)))
  st <- list(a_s1 = mk(0, 70), a_s2 = mk(1, 70),
             b_s1 = mk(0, 90), b_s2 = mk(1, 90))

  sx <- .admCovResidData("SEX", st, agg)
  # Conditioned in every study, so a contrast between levels -- not a trend.
  expect_true(all(sx$disc))
  expect_true(all(sx$paired))
  expect_equal(sort(unique(sx$source)), c("a", "b"))

  wt <- .admCovResidData("WT", st, agg)
  # Marginalised: a continuous axis, where a regression does mean something.
  expect_false(any(wt$disc))
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

test_that(".admCovStudyCentre reads a marginalised level mix as its mean", {
  # The median of an even binary split is a step function of `probs`: a cohort
  # at c(0.45, 0.55) and one at c(0.55, 0.45) are all but the same, and their
  # medians are the full width of the axis apart. The mean is the mixture.
  sp <- list(values = c(0, 1), probs = c(0.7, 0.3))
  s1 <- list(cov_dist = list(SEX = sp), cov = list(SEX = 0.3))
  expect_equal(.admCovStudyCentre(s1, "SEX"), 0.3, tolerance = 1e-12)
  # A conditioned study has no mixture to average: its centre is its value.
  s2 <- list(cov_dist = list(SEX = list(.point = TRUE)), cov = list(SEX = 1))
  expect_equal(.admCovStudyCentre(s2, "SEX"), 1)
  # A continuous margin is unchanged -- still the median.
  s3 <- list(cov_dist = list(WT = list(mu = 70, sd = 10)), cov = list(WT = 70))
  expect_equal(.admCovStudyCentre(s3, "WT"), 70, tolerance = 1e-8)
})

test_that(".admCovStudySupport reads declared levels, not the centre", {
  # An even split sits on NEITHER level by median, yet covers both. Reading
  # centres would grey a level every source sampled as extrapolation.
  s1 <- list(cov_dist = list(SEX = list(values = c(0, 1), probs = c(0.5, 0.5))),
             cov = list(SEX = 0.5))
  expect_equal(.admCovStudySupport(s1, "SEX", c(0, 1)), c(0, 1))
  # A zero probability is a level the source declared and did not enrol.
  s2 <- list(cov_dist = list(SEX = list(values = c(0, 1, 2),
                                        probs = c(0.5, 0.5, 0))),
             cov = list(SEX = 0.5))
  expect_equal(.admCovStudySupport(s2, "SEX", c(0, 1, 2)), c(0, 1))
  # A conditioned study covers the one level it was solved at.
  s3 <- list(cov_dist = list(SEX = list(.point = TRUE)), cov = list(SEX = 1))
  expect_equal(.admCovStudySupport(s3, "SEX", c(0, 1)), 1)
})

test_that(".admCovPalette gives a source the same colour in both panels", {
  # The effect panel draws one mark per SOURCE while the residual panel keeps
  # the strata apart, so it sees `a` where the other sees `a_s1`/`a_s2`. An
  # unnamed palette then hands `b` a different position in each and it changes
  # colour across one figure.
  pal <- .admCovPalette(c("a", "b", "a_s1", "a_s2"))
  expect_equal(names(pal), c("a", "a_s1", "a_s2", "b"))
  expect_equal(length(unique(pal)), 4L)
  expect_true(!is.null(names(.admCovPalette(character(0)))) ||
                length(.admCovPalette(character(0))) == 0L)
})

test_that(".admCovSourceRange prefers what the source declared", {
  # The enrolled range, then the table it enrolled, then the body of the
  # margin it declared -- most authoritative first.
  s1 <- list(range = list(WT = c(50, 110)),
             population = list(WT = list(meanlog = log(70), sdlog = 0.2)))
  expect_equal(.admCovSourceRange(s1, "WT"), c(50, 110))
  # `population` is a covDist by the time it reaches the plot -- admStudy()
  # fits margins and keeps those, not the rows -- so the specs live under the
  # covariate names. Reading only a data frame found a range for nobody.
  s2 <- list(population = list(WT = list(meanlog = log(70), sdlog = 0.2)))
  r <- .admCovSourceRange(s2, "WT")
  expect_true(r[1L] < 70 && r[2L] > 70)
  # A raw baseline table: what it literally enrolled.
  s3 <- list(population = data.frame(WT = c(52, 61, 88)))
  expect_equal(.admCovSourceRange(s3, "WT"), c(52, 88))
  # A source that says nothing gets no line rather than a guessed range.
  expect_null(.admCovSourceRange(list(), "WT"))
})

.pan_src <- function(cond = "WT") {
  pop <- list(WT  = list(meanlog = log(75), sdlog = 0.2),
              SEX = list(values = c(0, 1), probs = c(0.5, 0.5)))
  cd <- list(WT  = list(meanlog = log(75), sdlog = 0.2),
             SEX = list(values = c(0, 1), probs = c(0.5, 0.5)))
  for (cv in cond) cd[[cv]] <- list(.point = TRUE)
  list(src = list(paper = list(ui = .cov_ui(),
                               range = list(WT = c(60, 90), SEX = c(0, 1)),
                               population = pop)),
       studies = list(paper_s1 = list(n = 100L, cov = list(WT = 75, SEX = 0),
                                      cov_dist = cd)))
}

test_that(".admFitSourceStudies is absent rather than fatal", {
  expect_null(.admFitSourceStudies(list()))
  expect_null(.admFitSourceStudies(list(env = new.env())))
})

test_that(".admCovEffectData draws a conditional source's own regression", {
  skip_if_not_installed("rxode2")
  # A source CONDITIONAL on the covariate reported a relationship along this
  # axis -- its model estimated the effect, which is what let it be banded or
  # read at a value. Its own line over the range it covers, against the dotted
  # estimated effect, is the comparison the panel exists for.
  #
  # One source BANDED into two strata, which is what gives it a contrast of its
  # own to draw. A source that reported a single level has no slope of its own,
  # whatever its model estimates, and correctly gets no line.
  band <- function(sex) list(
    n = 100L, cov = list(WT = 80, SEX = sex),
    cov_dist = list(WT  = list(meanlog = log(80), sdlog = 0.2),
                    SEX = list(.point = TRUE)))
  st  <- list(a_s1 = band(0), a_s2 = band(1))
  src <- list(a = list(
    ui = .cov_ui(), range = list(WT = c(60, 100)),
    population = list(WT  = list(meanlog = log(80), sdlog = 0.2),
                      SEX = list(values = c(0, 1), probs = c(0.5, 0.5)))))

  d <- .admCovEffectData(.cov_ui(), "SEX", st,
                         list(tcl = log(5), bwt = 0.75, bsex = 0.2), src)
  expect_false(is.null(d$slines))
  expect_equal(unique(d$slines$study), "a")
  # Its own line joins the levels it reported. `.cov_ui()` at its own ini:
  # exp(0.2) between them.
  y <- d$slines$y[order(d$slines$x)]
  expect_equal(y[2L] / y[1L], exp(0.2), tolerance = 1e-8)

  # A source reporting ONE level has no contrast of its own to draw.
  one <- .admCovEffectData(.cov_ui(), "SEX", .cov_studies(),
                           list(tcl = log(5), bwt = 0.75, bsex = 0.2),
                           .cov_src())
  expect_null(one$slines)

  # A MARGINAL source has no line either. It reported no contrast along this
  # axis; its whisker says what it covered and nothing about slope.
  w <- .admCovEffectData(.cov_ui(), "WT", .cov_studies(),
                         list(tcl = log(5), bwt = 0.75, bsex = 0.2),
                         .cov_src())
  expect_true(all(w$marks$kind == "marginal"))
  expect_null(w$slines)
})

test_that(".admCovEffectData needs a declared range for a continuous one", {
  skip_if_not_installed("rxode2")
  # A conditional CONTINUOUS covariate: the line needs an extent, and the only
  # thing that gives one is the range the source declared it enrolled. A point
  # value gets a diamond and no line.
  # TEN distinct conditioned values, which is past .admCovLevels()' level cap.
  # Fewer than that and a set of pinned values reads as a factor, and the
  # discrete branch joins the levels instead -- needing no declared range.
  pin <- function(wt) list(
    n = 100L, cov = list(WT = wt, SEX = 0),
    cov_dist = list(SEX = list(.point = TRUE), WT = list(.point = TRUE)))
  wts <- c(58, 63, 68, 73, 78, 83, 88, 93, 98, 103)
  st  <- stats::setNames(lapply(wts, pin), sprintf("a_s%d", seq_along(wts)))

  pop <- list(SEX = list(values = c(0, 1), probs = c(0.5, 0.5)))
  with_rng <- list(a = list(ui = .cov_ui(), range = list(WT = c(60, 100)),
                            population = pop))
  no_rng   <- list(a = list(ui = .cov_ui(), population = pop))
  th <- list(tcl = log(5), bwt = 0.75, bsex = 0.2)

  d <- .admCovEffectData(.cov_ui(), "WT", st, th, with_rng)
  expect_false(is.null(d$slines))
  expect_equal(range(d$slines$x), c(60, 100))

  expect_null(.admCovEffectData(.cov_ui(), "WT", st, th, no_rng)$slines)
})

test_that(".admCovPalette leaves black to the estimated effect", {
  # Both panels draw the fit in black. A source in black could not be told from
  # the thing it is being compared against -- and on a level axis, where both
  # are a line joining two points, they were indistinguishable.
  expect_false("#000000" %in% .admCovPalette(c("a", "b", "c")))
})

test_that(".admMixMoments is the law of total variance, not an average", {
  # Two equally weighted strata, means (1, 3) and (2, 5), identity variances.
  # Mean: (1.5, 4). Variance: I + the between term, which is
  # tcrossprod(c(0.5, 1)) -- and it is the between term that carries the
  # covariate effect the banding created.
  E <- list(c(1, 3), c(2, 5))
  V <- list(diag(2), diag(2))
  m <- .admMixMoments(E, V, c(1, 1))
  expect_equal(m$E, c(1.5, 4))
  expect_equal(m$V, diag(2) + tcrossprod(c(0.5, 1)))
  # Weights are normalised, so the scale of `n` cannot matter.
  expect_equal(.admMixMoments(E, V, c(50, 50))$V, m$V)
  # n-WEIGHTED: a stratum with more patients pulls the mean further.
  expect_equal(.admMixMoments(E, V, c(3, 1))$E, c(1.25, 3.5))
  # One stratum collapses to itself -- the between term is zero.
  one <- .admMixMoments(E[1], V[1], 1)
  expect_equal(one$E, c(1, 3))
  expect_equal(one$V, diag(2))
})

test_that(".admCollapseSources leaves an unbanded source alone", {
  # A single-stratum source keeps its `cov`, so an `at`-pinned study still
  # titles with the value it was solved at. Only a genuinely banded source
  # loses that, because it no longer sits at one value.
  st <- list(solo = list(E = c(1, 2), V = diag(2), n = 50,
                         times = c(1, 2), cov = list(CRCL = 62)))
  ag <- list(solo = list(obs  = list(E = c(1, 2), V = diag(2)),
                         pred = list(E = c(1, 2), V = diag(2))))
  out <- .admCollapseSources(st, ag)
  expect_equal(names(out$studies), "solo")
  expect_equal(out$studies$solo$cov, list(CRCL = 62))
  expect_equal(out$studies$solo$n, 50)
})

test_that(".admCollapseSources renormalises over the strata it has", {
  # A stratum whose simulation failed drops out. A partial collapse is a worse
  # answer than a whole one and a better answer than losing the source.
  st <- list(a_s1 = list(E = c(2), V = matrix(1), n = 30, times = 1),
             a_s2 = list(E = c(4), V = matrix(1), n = 10, times = 1),
             a_s3 = list(E = c(9), V = matrix(1), n = 60, times = 1))
  ag <- list(a_s1 = list(pred = list(E = c(2), V = matrix(1))),
             a_s2 = list(pred = list(E = c(4), V = matrix(1))),
             a_s3 = NULL)
  out <- .admCollapseSources(st, ag)
  expect_equal(names(out$studies), "a")
  # Weighted over s1 and s2 only: (30*2 + 10*4)/40 = 2.5.
  expect_equal(as.numeric(out$studies$a$E), 2.5)
  expect_equal(as.numeric(out$studies$a$n), 40)
})

test_that(".admVarParts names what V actually contains", {
  # `sqrt(diag(V))` is the total SD of one observation across subjects -- NOT
  # between-subject variability, which is only one of its terms. Labelling it
  # as BSV points a reader at omega for a misfit that may be all error model.
  plain <- list(a = list(E = 1, V = 1, n = 10))
  expect_equal(.admVarParts(plain, n_eta = 1L), "BSV + sigma")
  # No random effects: no BSV term to name.
  expect_equal(.admVarParts(plain, n_eta = 0L), "sigma")
  # A marginalised covariate puts its own spread into V as well.
  cov <- list(a = list(E = 1, V = 1, n = 10,
                       cov_dist = list(WT = list(mu = 70, sd = 10))))
  expect_equal(.admVarParts(cov, n_eta = 1L),
               "BSV + covariate spread + sigma")
  # A covariate the model stopped reading still shaped the reported V.
  drop <- list(a = list(E = 1, V = 1, n = 10,
                        .adm_cov_dropped = list(WT = list(mu = 70, sd = 10))))
  expect_equal(.admVarParts(drop, n_eta = 1L),
               "BSV + covariate spread + sigma")
})

test_that(".admCollapseSources carries the structural variance too", {
  # The predicted total and the pre-sigma part both collapse by the mixture
  # law. The structural part's BETWEEN term is the banded covariate's own
  # contribution: banding moved that covariate out of each stratum's spread
  # and into the spacing between them.
  st <- list(a_s1 = list(E = c(2), V = matrix(1), n = 50, times = 1),
             a_s2 = list(E = c(6), V = matrix(1), n = 50, times = 1))
  ag <- list(
    a_s1 = list(pred = list(E = c(2), V = matrix(2), V_struct = matrix(1))),
    a_s2 = list(pred = list(E = c(6), V = matrix(2), V_struct = matrix(1))))
  out <- .admCollapseSources(st, ag)
  # within 1 + between (4) = 5 on the structural part; total adds sigma.
  expect_equal(as.numeric(out$agg$a$pred$V_struct), 5)
  expect_equal(as.numeric(out$agg$a$pred$V), 6)
  # A transforming error model drops V_struct rather than rescale it, and the
  # collapse has to survive that.
  ag2 <- ag; ag2$a_s2$pred$V_struct <- NULL
  expect_null(.admCollapseSources(st, ag2)$agg$a$pred$V_struct)
})

test_that("admMoments returns the numbers the panels are drawn from", {
  fit <- .make_mock_fit()
  fit$env$admExtra$studies <- list(
    a_s1 = list(E = c(2, 1), V = diag(c(1, 1)), n = 50, times = c(1, 2)),
    a_s2 = list(E = c(6, 3), V = diag(c(1, 1)), n = 50, times = c(1, 2)))
  fit$env$aggData <- NULL
  # Stand in for the simulation: the accessor's job is the tidying, and a mock
  # fit has no model to solve.
  ag <- list(
    a_s1 = list(times = c(1, 2), n = 50,
                obs  = list(E = c(2, 1), V = diag(c(1, 1))),
                pred = list(E = c(2, 1), V = diag(c(4, 4)),
                            V_struct = diag(c(1, 1)))),
    a_s2 = list(times = c(1, 2), n = 50,
                obs  = list(E = c(6, 3), V = diag(c(1, 1))),
                pred = list(E = c(6, 3), V = diag(c(4, 4)),
                            V_struct = diag(c(1, 1)))))
  local_mocked_bindings(.admAggData = function(...) ag, .package = "admixr2")

  # BY SOURCE by default: the two strata are one paper.
  m <- admMoments(fit)
  expect_equal(unique(m$study), "a")
  expect_equal(nrow(m), 2L)
  expect_setequal(names(m), c("study", "source", "time", "n", "obs_mean",
                              "pred_mean", "obs_sd", "pred_sd", "struct_sd",
                              "z"))
  expect_equal(m$n, c(100, 100))
  # Collapsed by the mixture law: means (2, 6) and (1, 3) -> (4, 2).
  expect_equal(m$obs_mean, c(4, 2))
  # struct_sd is the PRE-SIGMA part, so strictly inside pred_sd.
  expect_true(all(m$struct_sd < m$pred_sd))

  # BY STRATUM keeps them apart, which is what the covariate panels need.
  ms <- admMoments(fit, by = "stratum")
  expect_setequal(ms$study, c("a_s1", "a_s2"))
  expect_equal(nrow(ms), 4L)
})

test_that("admMoments reports NA struct_sd for a transforming error model", {
  fit <- .make_mock_fit()
  fit$env$admExtra$studies <- list(
    a = list(E = c(2, 1), V = diag(c(1, 1)), n = 50, times = c(1, 2)))
  ag <- list(a = list(times = c(1, 2), n = 50,
                      obs  = list(E = c(2, 1), V = diag(c(1, 1))),
                      pred = list(E = c(2, 1), V = diag(c(4, 4)),
                                  V_struct = NULL)))
  local_mocked_bindings(.admAggData = function(...) ag, .package = "admixr2")
  # Not rescaled and not guessed: a transforming error model moves the mean, so
  # the pre-sigma variance is on a different scale and is reported as absent.
  expect_true(all(is.na(admMoments(fit)$struct_sd)))
})

test_that(".admCovResidData is one row per SOURCE per position on THIS axis", {
  # Banding on SEX splits every source in two, and on another covariate's axis
  # both halves land on the same value -- so a two-paper fit drew four points at
  # two positions, each pair differing only in a covariate that facet is not
  # about, and each carrying half its paper's n. The effect panel already marks
  # one position per source; this now matches it.
  mk <- function(crcl, sex) list(
    n = 50L, times = c(1, 2),
    cov = list(CRCL = crcl, SEX = sex),
    cov_dist = list(SEX = list(.point = TRUE),
                    CRCL = list(meanlog = log(crcl), sdlog = 0.2)))
  st <- list(a_s1 = mk(40, 0), a_s2 = mk(40, 1),
             b_s1 = mk(90, 0), b_s2 = mk(90, 1))
  ag <- stats::setNames(lapply(names(st), function(nm) list(
    obs  = list(E = c(1, 2)),
    pred = list(E = c(1.1, 2.1), V = diag(c(0.01, 0.04))))), names(st))

  # CRCL: both sex strata of a source share its renal value, so they combine.
  r <- .admCovResidData("CRCL", st, ag)
  expect_equal(nrow(r), 2L)
  expect_setequal(r$study, c("a", "b"))
  expect_true(all(r$n == 100))              # the whole source, not a stratum

  # SEX: the strata ARE the contrast and stay apart, at their own n.
  rs <- .admCovResidData("SEX", st, ag)
  expect_equal(nrow(rs), 4L)
  expect_setequal(rs$study, names(st))
  expect_true(all(rs$n == 50))
  expect_setequal(rs$x, c(0, 1))
})

test_that("a source banded into quadrature nodes is ONE mark, not one per node", {
  # THE DEFECT: banding is derived from the source's own model, so a source
  # whose model uses two continuous covariates is cut into `nodes^2` strata --
  # 81 by default, 162 with a sex split. Each was read as a position: 99 dots
  # for one paper on the renal axis, most carrying under two patients, and the
  # outermost node (three SDs out) set the axis, which then ran to 653 mL/min
  # with every paper below 150.
  zz    <- seq(-3.2, 3.2, length.out = 9)
  nodes <- exp(log(60) + 0.5 * zz)
  wt    <- stats::dnorm(zz) / sum(stats::dnorm(zz))   # as banding apportions n
  st <- c(
    # the banded source: one point spec per node, n split between them
    stats::setNames(lapply(seq_along(nodes), function(k) list(
      n = 210 * wt[k], times = c(1, 2), cov = list(CRCL = nodes[k]),
      cov_dist = list(CRCL = list(.point = TRUE)))),
      paste0("mild_s", seq_along(nodes))),
    # a marginal source, which keeps the distribution it declared
    list(normal = list(n = 260L, times = c(1, 2), cov = list(CRCL = 95),
                       cov_dist = list(CRCL = list(meanlog = log(95),
                                                   sdlog = 0.2)))))
  ag <- stats::setNames(lapply(names(st), function(nm) list(
    obs  = list(E = c(1, 2)),
    pred = list(E = c(1.1, 2.1), V = diag(c(0.01, 0.04))))), names(st))

  r <- .admCovResidData("CRCL", st, ag)
  expect_equal(nrow(r), 2L)
  expect_setequal(r$study, c("mild", "normal"))
  # The whole source's patients, at its centre -- not a ninth of them at a node.
  expect_equal(r$n[r$study == "mild"], 210)
  expect_lt(r$x[r$study == "mild"], max(nodes))
  # And the span comes from the SOURCE, so the noise guard has a width to
  # measure: a stratum is a point spec whose own 10th and 90th are that node.
  expect_gt(r$xhi[r$study == "normal"], r$xlo[r$study == "normal"])
})
