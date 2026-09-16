# Covariate diagnostic panels -- runnable demo.
#
#   Rscript inst/examples/covariate-panels.R
#
# from the package root. Writes four PNGs into the working directory and says
# where it put them.
#
# Three cohorts sampled at three different renal medians, all generated WITH a
# renal effect. Fit A contains the renal term; fit B (the null) deletes it, so
# the effect it is missing has to surface in the between-study residual.

message("[1/7] loading admixr2 ...")
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

message("[2/7] building studies ...")
mk <- function(co) admStudy(model = with_renal, population = co, dose = 200,
                            times = c(0.5, 1, 2, 4, 8, 12, 24),
                            stratify = "SEX")
studies <- admStudies(normal   = mk(draw(260L, 95)),
                      mild     = mk(draw(210L, 62)),
                      moderate = mk(draw(180L, 38)))

run <- function(m) nlmixr2(m, admData(), est = "adgh",
  control = adghControl(studies = studies, print = 0L, n_restart = 1L,
                        maxeval = 40L))

message("[3/7] fitting the model WITH the renal term (~30s) ...")
p_ok <- plot(run(with_renal), which = "covariate", n_sim = 300L)

message("[4/7] fitting the NULL model, renal term deleted (~30s) ...")
p_null <- plot(run(no_renal), which = "covariate", n_sim = 300L)

message("[5/7] writing PNGs ...")
out <- c("covariate-effect.png", "covariate-resid-ok.png",
         "covariate-resid-null.png")
ggsave(out[1], p_ok$covariate_effect,  width = 10, height = 6, dpi = 110)
ggsave(out[2], p_ok$covariate_resid,   width = 10, height = 4, dpi = 110)
ggsave(out[3], p_null$covariate_resid, width = 10, height = 4, dpi = 110)

# ---------------------------------------------------------------------------
# A CONDITIONAL CONTINUOUS covariate.
#
# Above, CRCL is marginalised: each source declares a distribution and the
# estimator integrates over it, so the source draws as a round median with a
# 10th-90th bar and a 2.5th-97.5th whisker.
#
# Here the same covariate is CONDITIONED instead -- ten sources, each reported
# at one renal value (`at = list(CRCL = ...)`), which is what a per-band summary
# table gives you (`population` therefore describes only WT and SEX -- `at`
# pins CRCL, and a population cannot also give it a distribution). It is still
# continuous: ten distinct values is past the level
# cap, so the axis is swept and the fitted curve is drawn across it. What
# changes is the source mark. A conditioned source has no distribution to show,
# so it is a DIAMOND with no bar at all, sitting on the curve at its own value.
#
# Read together: the curve is the model's claim about the whole range, and the
# ten diamonds are the only places anybody measured. Grey is the rest.
message("[6/7] fitting a CONDITIONAL continuous covariate (~40s) ...")
crcls <- c(28, 36, 45, 55, 66, 78, 90, 104, 118, 135)
cond <- do.call(admStudies, setNames(lapply(crcls, function(v)
  admStudy(model = with_renal, population = draw(140L, 90)[c("WT", "SEX")],
           dose = 200, times = c(0.5, 2, 8, 24), at = list(CRCL = v))),
  paste0("band", seq_along(crcls))))

fit_c <- nlmixr2(with_renal, admData(), est = "adgh",
  control = adghControl(studies = cond, print = 0L, n_restart = 1L,
                        maxeval = 30L))
p_cond <- plot(fit_c, which = "covariate", n_sim = 200L)
out <- c(out, "covariate-effect-conditional.png")
ggsave(out[4], p_cond$covariate_effect, width = 10, height = 6, dpi = 110)

message("[7/7] done. Wrote:")
for (f in out) message("   ", normalizePath(f, winslash = "\\"))
message("\nLook at covariate-resid-null.png: the SEX facet joins each source's")
message("two strata in grey (a contrast, flat = correctly specified), while")
message("CRCL keeps a dashed lm and plunges -- the deleted renal term.")
message("\nAnd covariate-effect-conditional.png: CRCL swept continuously, with")
message("ten DIAMONDS and no bars -- conditioned sources have no spread to draw")
message("-- against WT on the same figure, marginalised, round with bars.")
