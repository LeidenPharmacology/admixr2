## Assertions on the two covariate panels.
##
## WHY NOT SNAPSHOTS. These were seven vdiffr snapshots, and the case against
## them is one of their own: `resid-paired-strata` was recorded with its palette
## keyed on `study` while the panel colours on `source`, so the manual scale
## matched nothing and ggplot2 fell back to default colours. It warned rather
## than failing, the SVG recorded the fallback, and the test passed for its
## whole life -- because nobody reads an SVG diff. A snapshot catches CHANGE,
## not correctness, and every time one fired here the answer was "yes, that
## change was intended".
##
## So each test below asserts the thing its defect was: which facets got a
## trend line, whether a level nobody enrolled is shaded, whether the axis
## ticks on integers, whether a source has a span at all. Those read in a
## failure message, they do not move when ggplot2 changes its fonts, and they
## cannot be satisfied by a picture that happens to look the same.
##
## No fit is involved. The panel builders take the data frames
## .admCovEffectData()/.admCovResidData() produce, so synthetic studies render
## in milliseconds and deterministically -- which is the whole reason
## plot.admFit() was split in the first place.

skip_if_no_panels <- function() {
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
  # The data and the palette come back with the plot: a test that can only see
  # the rendered result is a test that has to look at a picture.
  list(eff = eff, pal = pal, p = .admCovEffectPanel(eff, pal))
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
  # From `source`, which is what the panel colours on -- plot.admFit() builds
  # it the same way.
  pal <- .admCovPalette(unlist(lapply(res, `[[`, "source"), use.names = FALSE))
  list(res = res, pal = pal, p = .admCovResidPanel(res, pal))
}

## ---- assertions on the built panels ----------------------------------------

## What a layer maps, so a test can say "nothing keys on linetype any more"
## rather than recording a picture in which nothing does.
.pan_aes <- function(p) unique(unlist(lapply(p$layers, function(l)
  names(l$mapping))))

.pan_geoms <- function(p) vapply(p$layers, function(l) class(l$geom)[1L], "")

## The layers carrying a statistical summary, and the facets they were given.
## The `lm` guard is a per-facet decision, so the question is always "which
## facets got one", never "is there a line somewhere on the figure".
.pan_smooth_covs <- function(p) {
  keep <- vapply(p$layers, function(l)
    inherits(l$stat, "StatSmooth") || inherits(l$geom, "GeomSmooth"),
    logical(1))
  unique(unlist(lapply(p$layers[keep], function(l)
    if (is.data.frame(l$data)) unique(l$data$cov))))
}

## The x-axis breaks ggplot2 settled on, which is where a level axis ticked at
## 0.5 would show up.
.pan_xbreaks <- function(p) {
  b <- ggplot2::ggplot_build(p)
  br <- lapply(b$layout$panel_params, function(pp) pp$x$breaks)
  lapply(br, function(v) v[!is.na(v)])
}

## ---- effect panel -----------------------------------------------------------

test_that("effect panel: a split facet and an unsplit one on one figure", {
  skip_if_no_panels()
  # PINS: the blank linetype key. The WT facet split on the conditioned SEX and
  # the SEX facet had nothing to split on, so it carried level "". Pooled into
  # one scale, "" sorted to the front of an unnamed palette, took `solid` from
  # `SEX = 0` and added an unlabelled key to the legend.
  st <- list(
    lo = list(n = 100L, cov = list(WT = 70, SEX = 0),
              cov_dist = list(WT = .pan_lnorm(70), SEX = list(.point = TRUE))),
    hi = list(n = 300L, cov = list(WT = 90, SEX = 1),
              cov_dist = list(WT = .pan_lnorm(90), SEX = list(.point = TRUE))))
  o <- .pan_eff(.pan_ui(), c("WT", "SEX"), st)

  # There is no level aesthetic left to carry a blank key: the effect is ONE
  # line per facet, and a source is told apart by colour and shape.
  expect_false("linetype" %in% .pan_aes(o$p))
  expect_setequal(intersect(c("colour", "shape"), .pan_aes(o$p)),
                  c("colour", "shape"))
  # One curve per (facet, parameter), not one per level of the other covariate.
  cu <- o$eff[[1L]]$curve
  expect_equal(anyDuplicated(cu[cu$param == "cl", "x"]), 0L)
  # WT is marginalised by both and SEX conditioned by both, so the two facets
  # carry different kinds -- which is the case the pooled scale broke on.
  expect_setequal(o$eff[[1L]]$marks$kind, "marginal")
  expect_setequal(o$eff[[2L]]$marks$kind, "conditional")
})

test_that("effect panel: a discrete sweep with a second covariate conditioned", {
  skip_if_no_panels()
  # PINS: indistinguishable discrete points. Sweeping SEX while conditioning
  # REGION draws a mark per source at each level, and with no aesthetic telling
  # the sources apart nothing said which was REGION = 1.
  st <- list(
    north = list(n = 120L, cov = list(SEX = 0.5, REGION = 0),
                 cov_dist = list(SEX = list(values = c(0, 1),
                                            probs = c(0.5, 0.5)),
                                 REGION = list(.point = TRUE))),
    south = list(n = 180L, cov = list(SEX = 0.5, REGION = 1),
                 cov_dist = list(SEX = list(values = c(0, 1),
                                            probs = c(0.5, 0.5)),
                                 REGION = list(.point = TRUE))))
  o <- .pan_eff(.pan_ui_two_disc(), c("SEX", "REGION"), st)
  # Two sources, two colours, and the colours are actually distinct.
  expect_length(unique(o$pal), 2L)
  expect_setequal(names(o$pal), c("north", "south"))
  # Both sources present on the REGION facet, each at its own level.
  rg <- o$eff[[2L]]$marks
  expect_setequal(rg$study, c("north", "south"))
  expect_equal(length(unique(rg$x)), 2L)
})

test_that("effect panel: a declared level nobody enrolled is shaded", {
  skip_if_no_panels()
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
  o  <- .pan_eff(.pan_ui(), c("SEX"), st)
  sh <- o$eff[[1L]]$shade
  # Level 2 and only level 2.
  expect_equal(nrow(sh), 1L)
  expect_true(sh$xmin < 2 && sh$xmax >= 2 - 1e-9)
  expect_true(sh$xmin > 1)
  # PINS the clamp: a rect hanging past the outermost level widened the limits
  # enough that .admLevelBreaks() fell back to pretty() and ticked a
  # three-level factor at 0.5 and 1.5.
  expect_lte(sh$xmax, 2)
  for (br in .pan_xbreaks(o$p))
    expect_equal(br, br[br == round(br)])
})

test_that("effect panel: an even binary mix shades neither level", {
  skip_if_no_panels()
  # PINS: the median collapse. .admCovQuantile() reads the median of a 50/50
  # split as the UPPER level, so every source landed on 1 and level 0 -- which
  # both sources sampled -- was greyed as extrapolation.
  mix <- list(values = c(0, 1), probs = c(0.5, 0.5))
  st <- list(
    light = list(n = 100L, cov = list(WT = 60, SEX = 0.5),
                 cov_dist = list(WT = .pan_lnorm(60), SEX = mix)),
    heavy = list(n = 100L, cov = list(WT = 95, SEX = 0.5),
                 cov_dist = list(WT = .pan_lnorm(95), SEX = mix)))
  o <- .pan_eff(.pan_ui(), c("WT", "SEX"), st)
  # Both sources sampled both levels, so nothing is extrapolation.
  expect_null(o$eff[[2L]]$shade)
  # Each source sits at the mean of the mixture it declared, not on a level.
  expect_setequal(o$eff[[2L]]$marks$x, 0.5)
  expect_setequal(o$eff[[2L]]$marks$kind, "marginal")
  # And the WT facet gets ONE curve, not a SEX = 0 and a SEX = 1 line that no
  # source is on.
  cu <- o$eff[[1L]]$curve
  expect_equal(anyDuplicated(cu[cu$param == "cl", "x"]), 0L)
})

test_that("effect panel: a conditional CONTINUOUS covariate", {
  skip_if_no_panels()
  # Ten sources each reported at one renal value, which is what a per-band
  # summary table gives you. Past the level cap, so the axis is swept -- but
  # every source is a diamond with no bar, because a conditioned source has no
  # distribution to show.
  crcl <- c(28, 36, 45, 55, 66, 78, 90, 104, 118, 135)
  st <- stats::setNames(lapply(crcl, function(v)
    list(n = 140L, cov = list(WT = 76, CRCL = v),
         cov_dist = list(WT = .pan_lnorm(76), CRCL = list(.point = TRUE)))),
    paste0("band", seq_along(crcl)))
  o  <- .pan_eff(.pan_ui_renal(), c("CRCL"), st)
  mk <- o$eff[[1L]]$marks
  expect_equal(nrow(mk), length(crcl))
  expect_setequal(mk$kind, "conditional")
  # Zero-width spans: what each reported is a value, not a distribution.
  expect_equal(mk$xlo, mk$x)
  expect_equal(mk$xhi, mk$x)
  # Swept, not ticked at ten levels.
  expect_false(any(o$eff[[1L]]$curve$disc))
  expect_gt(nrow(o$eff[[1L]]$curve), length(crcl))
})

## ---- residual panel ---------------------------------------------------------

test_that("residual panel: two studies on a continuous covariate get no trend", {
  skip_if_no_panels()
  # PINS: the global `lm` guard. CRCL has four studies and earns a trend; WT
  # has two and must not get one -- two points define a line through themselves
  # and say nothing. The guard used to be `any(table(cov) > 2)` across the whole
  # figure, so WT got a dashed line as soon as CRCL qualified.
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
  o <- .pan_resid(c("WT", "CRCL"), st, z)
  expect_setequal(.pan_smooth_covs(o$p), "CRCL")
})

test_that("residual panel: strata of one source are joined, not regressed", {
  skip_if_no_panels()
  # PINS: the contrast reading on a discrete axis. Banding cuts one source into
  # strata differing only in SEX, so the PAIR is the evidence it bought.
  # Several sources tilting the same way is the mis-specification; a regression
  # over the pooled cloud averages the pairs away.
  mk <- function(wt, sex) list(
    n = 120L, cov = list(WT = wt, SEX = sex),
    cov_dist = list(WT = .pan_lnorm(wt), SEX = list(.point = TRUE)))
  st <- list(alpha_s1 = mk(70, 0), alpha_s2 = mk(70, 1),
             beta_s1  = mk(88, 0), beta_s2  = mk(88, 1))
  z <- list(alpha_s1 = -1.4, alpha_s2 = 1.1,
            beta_s1  = -1.9, beta_s2  = 0.7)
  o <- .pan_resid(c("SEX"), st, z)
  # No trend line anywhere: a level axis is read as a contrast.
  expect_length(.pan_smooth_covs(o$p), 0L)
  # Each source appears as a pair, joined by a connector grouped on `source`.
  expect_true(all(o$res[[1L]]$paired))
  expect_equal(nrow(o$res[[1L]]), 4L)
  grp <- vapply(o$p$layers, function(l) {
    g <- l$mapping$group
    if (is.null(g)) "" else deparse(rlang::quo_get_expr(g))
  }, "")
  expect_true("source" %in% grp)
  # PINS the palette keying. Built from `study` it matched nothing, ggplot2
  # warned rather than failed, and the figure fell back to default colours.
  expect_setequal(names(o$pal), c("alpha", "beta"))
  expect_length(unique(o$pal), 2L)
})
