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

test_that(".admMergeCovMarks merges coincident strata under their source name", {
  df <- data.frame(
    cov = "WT", param = "cl", study = c("a_s1", "a_s2", "b_s1"),
    x = c(70, 70, 90), xlo = 60, xhi = 80, y = c(5, 5, 6),
    stringsAsFactors = FALSE)
  out <- .admMergeCovMarks(df)
  expect_equal(nrow(out), 2L)
  expect_setequal(out$study, c("a", "b"))
})

test_that(".admMergeCovMarks unions the spread of the marks it absorbs", {
  # The merged mark is labelled for both strata, so it has to draw both their
  # ranges. Keeping the first row's silently showed one stratum's coverage
  # under a label claiming the pair's.
  df <- data.frame(
    cov = "WT", param = "cl", study = c("a_s1", "a_s2"),
    kind = "marginal", x = 70, y = 5,
    xlo = c(60, 55), xhi = c(80, 92), xlo2 = c(50, 45), xhi2 = c(90, 99),
    stringsAsFactors = FALSE)
  out <- .admMergeCovMarks(df)
  expect_equal(nrow(out), 1L)
  expect_equal(c(out$xlo, out$xhi, out$xlo2, out$xhi2), c(55, 92, 45, 99))
})

test_that(".admMergeCovMarks takes the weaker claim when merged kinds disagree", {
  # `conditional` is the claim-less reading: a diamond and no bar. Calling the
  # merge marginal would hand a study solved at one value the 10th-90th spread
  # of whichever study it happened to land on -- a distribution out of a point.
  df <- data.frame(
    cov = "WT", param = "cl", study = c("a_s1", "a_s2"),
    kind = c("marginal", "conditional"), x = 70, y = 5,
    xlo = c(60, 70), xhi = c(80, 70), xlo2 = c(50, 70), xhi2 = c(90, 70),
    stringsAsFactors = FALSE)
  out <- .admMergeCovMarks(df)
  expect_equal(out$kind, "conditional")
  # And no borrowed spread with it.
  expect_equal(c(out$xlo, out$xhi, out$xlo2, out$xhi2), c(60, 80, 50, 90))
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
  # The effect panel's marks come through .admMergeCovMarks(), which collapses
  # `a_s1`/`a_s2` into `a`; the residual panel keeps them apart. An unnamed
  # palette then hands `b` a different position in each and it changes colour
  # across one figure.
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
