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
# covariate: a diamond with a solid line across the band it reported,
# against two whiskers over the distributions their populations imply.
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

# CONDITIONAL on renal function: this paper reported its result AT a renal
# value, the way a per-subgroup table gives it to you. `at` pins it, so CRCL
# leaves `population` -- a population cannot also give a pinned covariate a
# distribution.
#
# On the CRCL panel it therefore draws as a DIAMOND with no spread, against the
# other two sources' whiskers. That is the difference the panel exists to show:
# this paper reported a result for a renal value, those two only reported who
# they enrolled.
#
# It does not band on sex. A paper that published a renal breakdown need not
# have published a sex one, and you can only band on what was reported.
mk_at <- function(co, v)
  admStudy(model = with_renal, population = co[c("WT", "SEX")], dose = 200,
           times = TIMES, at = list(CRCL = v))

studies <- admStudies(normal   = mk(draw(260L, 95)),
                      mild     = mk_at(draw(210L, 62), 62),
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
message("cl vs CRCL: 'mild' reported AT a renal value, so it is CONDITIONAL --")
message("a diamond, no spread. 'normal' and 'moderate' reported only who they")
message("enrolled, so they are MARGINAL -- a whisker over the distribution")
message("admixr2 integrates. cl vs SEX: the two banded analysts")
message("banded on sex, so each source has one mark per level, and the gap")
message("between its own marks is that paper's own sex effect -- read it")
message("against the gap in the dotted estimated effect.")
