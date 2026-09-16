## Visual regression for the covariate panels.
##
## WHY SNAPSHOTS AND NOT MORE ASSERTIONS. Every defect both rounds of review on
## these panels turned up was a RENDERING defect -- a studied level greyed as
## extrapolation, a blank key stealing `solid` from every real level, two
## identical black dots with nothing to tell them apart, a trend line through
## another covariate's two points. The data frames underneath were all correct,
## and the suite passed through every one of them. What was missing was not
## another `expect_equal()` on a column; it was a test of what gets drawn.
##
## No fit is involved. The panel builders take the data frames
## .admCovEffectData()/.admCovResidData() produce, so synthetic studies render
## in milliseconds and deterministically -- which is the whole reason
## plot.admFit() was split in the first place.
##
## One case per bug that actually happened. A snapshot nobody can explain is a
## snapshot nobody will dare to update, so each has the defect it pins.

skip_if_no_vdiffr <- function() {
  skip_if_not_installed("vdiffr")
  skip_if_not_installed("rxode2")
  skip_if_not_installed("ggplot2")
}

## ---- fixtures ---------------------------------------------------------------

## cl reads WT and SEX; v reads neither.
.pan_ui <- function() {
  fn <- function() {
    ini({
      tcl <- log(5); bwt <- 0.75; bsex <- 0.2
      add.err <- 0.1
      eta.cl ~ 0.1
    })
    model({
      cl <- exp(tcl + eta.cl) * (WT / 70)^bwt * exp(bsex * SEX)
      v  <- exp(log(30)) * (WT / 70)
      cp <- linCmt()
      cp ~ add(add.err)
    })
  }
  suppressMessages(rxode2::rxode2(fn))
}

## cl reads SEX and REGION, both discrete.
.pan_ui_two_disc <- function() {
  fn <- function() {
    ini({
      tcl <- log(5); bsex <- 0.2; breg <- 0.35
      add.err <- 0.1
      eta.cl ~ 0.1
    })
    model({
      cl <- exp(tcl + eta.cl + bsex * SEX + breg * REGION)
      v  <- exp(log(30))
      cp <- linCmt()
      cp ~ add(add.err)
    })
  }
  suppressMessages(rxode2::rxode2(fn))
}

## cl reads CRCL and WT.
.pan_ui_renal <- function() {
  fn <- function() {
    ini({
      tcl <- log(5); bcrcl <- 0.6
      add.err <- 0.1
      eta.cl ~ 0.1
    })
    model({
      cl <- exp(tcl + eta.cl) * (WT / 70)^0.75 * (CRCL / 90)^bcrcl
      v  <- exp(log(30))
      cp <- linCmt()
      cp ~ add(add.err)
    })
  }
  suppressMessages(rxode2::rxode2(fn))
}

.pan_lnorm <- function(m, s = 0.2) list(meanlog = log(m), sdlog = s)

## The original studies, keyed by SOURCE: each carries its own published model,
## which is what puts a source at its own parameter value rather than on the
## fitted line.
.pan_src_of <- function(ui, studies, rng)
  stats::setNames(lapply(names(studies), function(nm)
    list(ui = ui, range = rng,
         population = studies[[nm]]$cov_dist)), names(studies))

.pan_eff <- function(ui, covs, studies, rng = list()) {
  src <- .pan_src_of(ui, studies, rng)
  eff <- Filter(Negate(is.null),
                lapply(covs, function(cv)
                  .admCovEffectData(ui, cv, studies, NULL, src)))
  pal <- .admCovPalette(unlist(lapply(eff, function(z) z$marks$study),
                               use.names = FALSE))
  .admCovEffectPanel(eff, pal)
}

## Aggregate data shaped so each study's mean standardised residual is exactly
## `z[[nm]]`: se is sqrt(1/n) on an identity V, so an offset of z*se gives z.
.pan_agg <- function(studies, z)
  stats::setNames(lapply(names(studies), function(nm) {
    se <- sqrt(1 / studies[[nm]]$n)
    list(obs  = list(E = rep(0, 3)),
         pred = list(E = rep(-z[[nm]] * se, 3), V = diag(3)))
  }), names(studies))

.pan_resid <- function(covs, studies, z) {
  agg <- .pan_agg(studies, z)
  res <- Filter(Negate(is.null),
                lapply(covs, function(cv) .admCovResidData(cv, studies, agg)))
  .admCovResidPanel(res, .admCovPalette(unlist(lapply(res, `[[`, "study"),
                                               use.names = FALSE)))
}

## ---- effect panel -----------------------------------------------------------

test_that("effect panel: a split facet and an unsplit one on one figure", {
  skip_if_no_vdiffr()
  # PINS: the blank linetype key. The WT facets split on the conditioned SEX,
  # the SEX facet has nothing to split on and carries level "". Pooled into one
  # scale, "" sorted to the front of an unnamed palette, took `solid` from
  # `SEX = 0` and added an unlabelled key to the legend.
  st <- list(
    lo = list(n = 100L, cov = list(WT = 70, SEX = 0),
              cov_dist = list(WT = .pan_lnorm(70), SEX = list(.point = TRUE))),
    hi = list(n = 300L, cov = list(WT = 90, SEX = 1),
              cov_dist = list(WT = .pan_lnorm(90), SEX = list(.point = TRUE))))
  vdiffr::expect_doppelganger(
    "effect-split-and-unsplit",
    .pan_eff(.pan_ui(), c("WT", "SEX"), st))
})

test_that("effect panel: a discrete sweep with a second covariate conditioned", {
  skip_if_no_vdiffr()
  # PINS: indistinguishable discrete points. Sweeping SEX while conditioning
  # REGION draws two points per level at different heights; with no `level`
  # aesthetic on the points, nothing said which was REGION = 1.
  st <- list(
    north = list(n = 120L, cov = list(SEX = 0.5, REGION = 0),
                 cov_dist = list(SEX = list(values = c(0, 1),
                                            probs = c(0.5, 0.5)),
                                 REGION = list(.point = TRUE))),
    south = list(n = 180L, cov = list(SEX = 0.5, REGION = 1),
                 cov_dist = list(SEX = list(values = c(0, 1),
                                            probs = c(0.5, 0.5)),
                                 REGION = list(.point = TRUE))))
  vdiffr::expect_doppelganger(
    "effect-discrete-two-conditioned",
    .pan_eff(.pan_ui_two_disc(), c("SEX", "REGION"), st))
})

test_that("effect panel: a declared level nobody enrolled is shaded", {
  skip_if_no_vdiffr()
  # PINS: extrapolation marking on a discrete axis. A `values` spec can declare
  # a group with zero probability; the model predicts for it happily, and drawn
  # like the studied levels that prediction looks equally earned.
  st <- list(
    a = list(n = 100L, cov = list(WT = 70, SEX = 0),
             cov_dist = list(WT = .pan_lnorm(70),
                             SEX = list(values = c(0, 1, 2),
                                        probs = c(1, 0, 0)))),
    b = list(n = 100L, cov = list(WT = 90, SEX = 1),
             cov_dist = list(WT = .pan_lnorm(90),
                             SEX = list(values = c(0, 1, 2),
                                        probs = c(0, 1, 0)))))
  vdiffr::expect_doppelganger(
    "effect-unstudied-level-shaded",
    .pan_eff(.pan_ui(), c("SEX"), st))
})

test_that("effect panel: an even binary mix shades neither level", {
  skip_if_no_vdiffr()
  # PINS: the median collapse. .admCovQuantile() reads the median of a 50/50
  # split as the UPPER level, so every source landed on 1 and level 0 -- which
  # both sources sampled -- was greyed as extrapolation. Also pins the split
  # rule: SEX is marginalised by everyone, so the WT facet gets ONE curve, not
  # a SEX = 0 and a SEX = 1 line that no source is on.
  mix <- list(values = c(0, 1), probs = c(0.5, 0.5))
  st <- list(
    light = list(n = 100L, cov = list(WT = 60, SEX = 0.5),
                 cov_dist = list(WT = .pan_lnorm(60), SEX = mix)),
    heavy = list(n = 100L, cov = list(WT = 95, SEX = 0.5),
                 cov_dist = list(WT = .pan_lnorm(95), SEX = mix)))
  vdiffr::expect_doppelganger(
    "effect-even-binary-mix",
    .pan_eff(.pan_ui(), c("WT", "SEX"), st))
})

test_that("effect panel: a conditional CONTINUOUS covariate", {
  skip_if_no_vdiffr()
  # Ten sources each reported at one renal value, which is what a per-band
  # summary table gives you. Past the level cap, so the axis is swept and the
  # curve is drawn across it -- but every source is a diamond with no bar,
  # because a conditioned source has no distribution to show.
  crcl <- c(28, 36, 45, 55, 66, 78, 90, 104, 118, 135)
  st <- stats::setNames(lapply(crcl, function(v)
    list(n = 140L, cov = list(WT = 76, CRCL = v),
         cov_dist = list(WT = .pan_lnorm(76), CRCL = list(.point = TRUE)))),
    paste0("band", seq_along(crcl)))
  vdiffr::expect_doppelganger(
    "effect-conditional-continuous",
    .pan_eff(.pan_ui_renal(), c("CRCL"), st))
})

## ---- residual panel ---------------------------------------------------------

test_that("residual panel: two studies on a continuous covariate get no trend", {
  skip_if_no_vdiffr()
  # PINS: the global `lm` guard. CRCL has four studies and earns a trend; WT
  # has two and must not get one -- two points define a line through themselves
  # and say nothing. The guard used to be `any(table(cov) > 2)` across the whole
  # figure, so WT got a dashed line as soon as CRCL qualified.
  # Only a and b report a weight, so WT has exactly two points; all four
  # report renal function.
  st <- list(
    a = list(n = 100L, cov = list(WT = 66, CRCL = 40),
             cov_dist = list(WT = .pan_lnorm(66), CRCL = .pan_lnorm(40, 0.1))),
    b = list(n = 150L, cov = list(WT = 94, CRCL = 60),
             cov_dist = list(WT = .pan_lnorm(94), CRCL = .pan_lnorm(60, 0.1))),
    c = list(n = 200L, cov = list(CRCL = 85),
             cov_dist = list(CRCL = .pan_lnorm(85, 0.1))),
    d = list(n = 250L, cov = list(CRCL = 110),
             cov_dist = list(CRCL = .pan_lnorm(110, 0.1))))
  z <- list(a = 2.4, b = 0.9, c = -0.6, d = -2.8)
  vdiffr::expect_doppelganger(
    "resid-trend-only-where-earned",
    .pan_resid(c("WT", "CRCL"), st, z))
})

test_that("residual panel: strata of one source are joined, not regressed", {
  skip_if_no_vdiffr()
  # PINS: the contrast reading on a discrete axis. `stratify` cuts one source
  # into strata differing only in SEX, so the PAIR is the evidence banding
  # bought. Several sources tilting the same way is the mis-specification; a
  # regression over the pooled cloud averages the pairs away.
  mk <- function(wt, sex) list(
    n = 120L, cov = list(WT = wt, SEX = sex),
    cov_dist = list(WT = .pan_lnorm(wt), SEX = list(.point = TRUE)))
  st <- list(alpha_s1 = mk(70, 0), alpha_s2 = mk(70, 1),
             beta_s1  = mk(88, 0), beta_s2  = mk(88, 1))
  z <- list(alpha_s1 = -1.4, alpha_s2 = 1.1,
            beta_s1  = -1.9, beta_s2  = 0.7)
  vdiffr::expect_doppelganger(
    "resid-paired-strata",
    .pan_resid(c("SEX"), st, z))
})
