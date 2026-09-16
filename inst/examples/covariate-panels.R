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
mk <- function(co) admStudy(model = with_renal, population = co, dose = 200,
                            times = c(0.5, 1, 2, 4, 8, 12, 24),
                            stratify = "SEX")
studies <- admStudies(normal   = mk(draw(260L, 95)),
                      mild     = mk(draw(210L, 62)),
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
