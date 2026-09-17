# Covariate diagnostic panels -- runnable demo.
#
#   Rscript inst/examples/covariate-panels.R
#
# from the package root. Writes three PNGs into the working directory and says
# where it put them.
#
# Three cohorts sampled at three different renal medians, all generated WITH a
# renal effect. Fit A contains the renal term; fit B (the null) deletes it, so
# the effect it is missing has to surface in the between-study residual.
#
# The figures carry both kinds of covariate, which is the distinction the
# panels are built around and it is a property of the SOURCE's own model:
#
#   CONDITIONAL -- the analyst estimated this effect, so the paper can be
#     read at one value of it. Every source here fitted sex, so SEX bands
#     (`stratify = "SEX"`) and its strata draw as diamonds.
#   MARGINAL -- the analyst never estimated it, or never reported a
#     breakdown by it. No contrast to band on, so admixr2 integrates over
#     the population enrolled instead, and the source draws as a round
#     point with a whisker.
#
# `mild` below is CONDITIONAL on renal function and the other two are
# MARGINAL on it, so one figure carries both readings of the same
# covariate: its own regression line across the band it reported, against
# two whiskers over the distributions their populations imply.
#
# Both are used for estimation. Splitting a source on a covariate its model
# never saw would manufacture a contrast that is not in the literature.

message("[1/6] loading admixr2 ...")
if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
  pkgload::load_all(".", quiet = TRUE)
} else {
  library(admixr2)
}
suppressMessages({
  library(rxode2); library(nlmixr2); library(ggplot2)
})

set.seed(11)
draw <- function(n, crclm) data.frame(
  WT   = rlnorm(n, log(76), 0.198),
  CRCL = rlnorm(n, log(crclm), 0.05),
  SEX  = rbinom(n, 1, 0.5))

with_renal <- function() {
  ini({
    tcl <- log(5); tv <- log(50); bcrcl <- 0.6; bsex <- 0.18
    add.err <- 0.08
    eta.cl ~ 0.05
  })
  model({
    cl <- exp(tcl + eta.cl) * (WT/70)^0.75 * (CRCL/90)^bcrcl * exp(bsex * SEX)
    v  <- exp(tv) * (WT/70)
    cp <- linCmt()
    cp ~ add(add.err)
  })
}

# The null model: the renal term is simply DELETED. The studies still declare
# CRCL -- that is still who was enrolled -- so admixr2 takes it off the design
# and the residual panel keeps plotting against it.
no_renal <- function() {
  ini({
    tcl <- log(5); tv <- log(50); bsex <- 0.18
    add.err <- 0.08
    eta.cl ~ 0.05
  })
  model({
    cl <- exp(tcl + eta.cl) * (WT/70)^0.75 * exp(bsex * SEX)
    v  <- exp(tv) * (WT/70)
    cp <- linCmt()
    cp ~ add(add.err)
  })
}

message("[2/6] building studies ...")
TIMES <- c(0.5, 1, 2, 4, 8, 12, 24)

# MARGINAL on renal function: the paper reports who it enrolled, admixr2
# integrates over that. Banded on sex, which every analyst did fit.
mk <- function(co) admStudy(model = with_renal, population = co, dose = 200,
                            times = TIMES, stratify = "SEX")

# CONDITIONAL on renal function: this paper reported its result BY renal band,
# so it is banded on CRCL as well as sex. `range` is the band it enrolled --
# without it admixr2 cuts the strata from the whole declared distribution and
# credits the source with evidence outside its own band.
#
# On the CRCL panel it therefore draws as ITS OWN REGRESSION LINE across
# 50-75 -- its published model, at its published values, over the range it
# covers -- against the other two sources' whiskers. That is the comparison the
# panel exists for: a source's own line at a different slope from the dotted
# estimated effect is a paper the meta-analysis does not reproduce.
#
# Its population is given as explicit MARGINS rather than a table. A table
# carries the observed covariate correlations, banding on a correlated
# covariate then builds a joint sampler for the stratum, and the null model
# cannot drop a covariate out of one -- see `?covDist`.
mk_band <- function(crclm, band)
  admStudy(model = with_renal,
           population = covDist(WT   = c(meanlog = log(76), sdlog = 0.198),
                                CRCL = c(meanlog = log(crclm), sdlog = 0.05),
                                SEX  = c(female = 0.5, male = 0.5)),
           n = 210L, dose = 200, times = TIMES,
           stratify = c("SEX", "CRCL"), strata_nodes = 1L,
           range = list(CRCL = band))

studies <- admStudies(normal   = mk(draw(260L, 95)),
                      mild     = mk_band(62, c(50, 75)),
                      moderate = mk(draw(180L, 38)))

run <- function(m) nlmixr2(m, admData(), est = "adgh",
  control = adghControl(studies = studies, print = 0L, n_restart = 1L,
                        maxeval = 40L))

message("[3/6] fitting the model WITH the renal term (~30s) ...")
p_ok <- plot(run(with_renal), which = "covariate", n_sim = 300L)

message("[4/6] fitting the NULL model, renal term deleted (~30s) ...")
p_null <- plot(run(no_renal), which = "covariate", n_sim = 300L)

message("[5/6] writing PNGs ...")
out <- c("covariate-effect.png", "covariate-resid-ok.png",
         "covariate-resid-null.png")
ggsave(out[1], p_ok$covariate_effect,  width = 10, height = 6, dpi = 110)
ggsave(out[2], p_ok$covariate_resid,   width = 10, height = 4, dpi = 110)
ggsave(out[3], p_null$covariate_resid, width = 10, height = 4, dpi = 110)

message("[6/6] done. Wrote:")
for (f in out) message("   ", normalizePath(f, winslash = "\\"))
message("\nLook at covariate-resid-null.png: the SEX facet joins each source's")
message("two strata in grey (a contrast, flat = correctly specified), while")
message("CRCL keeps a dashed lm and plunges -- the deleted renal term.")
message("\nAnd covariate-effect.png carries BOTH readings on one figure.")
message("cl vs CRCL: 'mild' reported BY renal band, so it is CONDITIONAL --")
message("its OWN regression line across 50-75, to compare against the dotted")
message("estimated effect. 'normal' and 'moderate' reported only who they")
message("enrolled, so they are MARGINAL -- a whisker over the distribution")
message("admixr2 integrates. cl vs SEX: every analyst")
message("banded on sex, so each source has one mark per level, and the gap")
message("between its own marks is that paper's own sex effect -- read it")
message("against the gap in the dotted estimated effect.")
