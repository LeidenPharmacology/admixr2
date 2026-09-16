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
  df <- data.frame(
    cov = "WT", param = "cl", study = c("a_s1", "a_s2"),
    kind = c("marginal", "conditional"), x = 70, y = 5,
    xlo = 60, xhi = 80, xlo2 = 50, xhi2 = 90, stringsAsFactors = FALSE)
  expect_equal(.admMergeCovMarks(df)$kind, "marginal")
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
  # SEX is conditioned at 0 in one study and 1 in the other, so the sweep is
  # drawn once per level rather than once at the pooled 0.5 -- a patient that
  # does not exist. 120 grid points per level.
  expect_setequal(unique(d$curve$level), c("SEX = 0", "SEX = 1"))
  expect_equal(nrow(d$curve), 240L)
  # Allometric with a positive exponent: monotone increasing in weight, within
  # each level. Across the concatenation it is not, and should not be.
  for (lv in unique(d$curve$level))
    expect_true(all(diff(d$curve$y[d$curve$level == lv]) > 0))
  # The conditioned covariate's own effect is the GAP between the levels.
  y0 <- d$curve$y[d$curve$level == "SEX = 0"]
  y1 <- d$curve$y[d$curve$level == "SEX = 1"]
  expect_true(all(y1 > y0))
  expect_equal(unique(round(y1 / y0, 8)), round(exp(0.2), 8))
  # Two shaded regions, one past each end of the range the studies cover.
  expect_equal(nrow(d$shade), 2L)
  expect_lt(min(d$curve$x), min(vapply(.cov_studies(), .admCovStudyQ,
                                       double(1), cv = "WT", u = 0.1)))
})

test_that(".admCovEffectData puts each study's mark on its OWN level's line", {
  skip_if_not_installed("rxode2")
  # `lo` is conditioned at SEX = 0, `hi` at SEX = 1. Each belongs on the line
  # for the level it was solved at; on a single pooled line both would sit at a
  # height no level of the model predicts.
  m <- .admCovEffectData(.cov_ui(), "WT", .cov_studies(),
                         list(tcl = log(5), bwt = 0.75, bsex = 0.2))$marks
  expect_equal(m$level[m$study == "lo"], "SEX = 0")
  expect_equal(m$level[m$study == "hi"], "SEX = 1")
  # Tolerance is for the INTERPOLATION, not the model: a mark's y is read off
  # the 120-point grid with approx(), so a convex curve lands a couple of parts
  # per million away from the closed form.
  expect_equal(m$y[m$study == "hi"] /
                 (5 * (m$x[m$study == "hi"] / 70)^0.75),
               exp(0.2), tolerance = 1e-4)
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
  expect_true(d$shade$xmin < 2 && d$shade$xmax > 2)

  # Every declared level studied: nothing to warn about.
  d2 <- .admCovEffectData(ui, "GRP", list(a = mk(0, c(0, 1)),
                                          b = mk(1, c(0, 1))),
                          list(tcl = log(5), bg = 0.2))
  expect_null(d2$shade)
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
