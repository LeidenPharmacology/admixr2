# END TO END, against the product grid the projection replaces.
#
# Node counts and matched moments do not validate a likelihood: the reduction
# can be exact in the INTEGRAL while the finite rules in the two coordinate
# systems still carry different quadrature errors, and the recombination
# preserves the polynomial moments in its basis rather than the nonlinear
# integrand the objective actually contains. So the check is the objective
# itself, across parameter points, ranges, correlations and both supported
# source ranks.
#
# The two paths are deliberately held to DIFFERENT tolerances, because they are
# not equally exact and saying so is the point: untruncated margins agree to
# ~1e-08, while a truncated margin goes through the cloud and its recombination
# and lands at ~3e-05.

.pc_src <- function() {
  ini({ tcl <- log(3.2); tv <- log(21); bwt <- .60; bcr <- .40; balb <- .30
        eta.cl ~ .09; add.err <- .35 })
  model({ cl <- exp(tcl + eta.cl)*(WT/70)^bwt*(CRCL/95)^bcr*(ALB/40)^balb
          v <- exp(tv); cp <- linCmt(); cp ~ add(add.err) })
}
# rank two: the covariates reach the model through cl AND v
.pc_src2 <- function() {
  ini({ tcl <- log(3.2); tv <- log(21); bwt <- .60; bcr <- .40; balb <- .30
        eta.cl ~ .09; add.err <- .35 })
  model({ cl <- exp(tcl + eta.cl)*(WT/70)^bwt*(CRCL/95)^bcr
          v  <- exp(tv)*(ALB/40)^balb; cp <- linCmt(); cp ~ add(add.err) })
}
.pc_ana <- function() {
  ini({ tcl <- log(3.2); tv <- log(21)
        bwt <- fix(.52); bcr <- fix(.55); balb <- fix(.18)
        eta.cl ~ .09; add.err <- .6 })
  model({ cl <- exp(tcl + eta.cl)*(WT/70)^bwt*(CRCL/95)^bcr*(ALB/40)^balb
          v <- exp(tv); cp <- linCmt(); cp ~ add(add.err) })
}

# Worst relative gap between the two designs over several parameter points.
.pc_gap <- function(src, pop, rng) {
  au  <- suppressMessages(rxode2::rxode2(.pc_ana))
  pin <- admixr2:::.admDriverPinfo(au, adghControl(studies = list(), print = 0L))
  rx  <- admixr2:::.admLoadModel(au); rxode2::rxLoad(rx)
  gr  <- admixr2:::.adghNodeGrid(7L, pin$n_eta)
  p0  <- admixr2:::.admBuildOptVec(pin)$p0
  mk  <- function() admStudies(a = admStudy(
    model = src, population = pop, n = 300L, dose = 200, times = c(1, 4),
    strata_nodes = 5L, range = rng))
  gp <- suppressMessages(suppressWarnings(admixr2:::.admMaterialise(mk())))
  gj <- suppressMessages(suppressWarnings(admixr2:::.admMaterialise(
    mk(), analysis_covs = c("WT", "CRCL", "ALB"), analysis_ui = au)))
  ep <- admixr2:::.admBuildEvFull(gp); ej <- admixr2:::.admBuildEvFull(gj)
  worst <- 0
  for (d in c(-0.3, 0, 0.3)) {
    p <- p0; p[1L] <- p[1L] + d
    a <- suppressWarnings(admixr2:::.adghNLL(p, pin, ep, rx, "cp", gr, 1L))
    b <- suppressWarnings(admixr2:::.adghNLL(p, pin, ej, rx, "cp", gr, 1L))
    expect_true(is.finite(a) && is.finite(b))
    worst <- max(worst, abs(a - b) / (abs(a) + 1))
  }
  list(worst = worst, np = length(gp), nj = length(gj))
}

.pc_pop <- function(...) admPopulation(
  WT = c(mean = 76, sd = 15), CRCL = c(mean = 92, sd = 22),
  ALB = c(mean = 40, sd = 5), ...)

test_that("an UNtruncated projection reproduces the product objective", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  for (cfg in list(
    list(nm = "independent",  pop = .pc_pop()),
    list(nm = "cor .45",      pop = .pc_pop(cor = c(WT.CRCL = 0.45))),
    list(nm = "cor .8/-.5",   pop = .pc_pop(cor = c(WT.CRCL = 0.8,
                                                    WT.ALB = -0.5))))) {
    g <- .pc_gap(.pc_src, cfg$pop, NULL)
    expect_lt(g$nj, g$np)
    expect_lt(g$worst, 1e-6)
  }
})

test_that("a TRUNCATED projection tracks the product objective, less exactly", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  for (rng in list(list(WT = c(55, 105), CRCL = c(60, 130), ALB = c(31, 49)),
                   list(WT = c(65, 90),  CRCL = c(75, 110), ALB = c(35, 45)))) {
    g <- .pc_gap(.pc_src, .pc_pop(cor = c(WT.CRCL = 0.45)), rng)
    expect_lt(g$nj, g$np)
    # NOT 1e-6: the cloud and its recombination each cost accuracy, and the
    # certificate bounds that cost rather than removing it
    expect_lt(g$worst, 1e-3)
  }
})

test_that("a rank-two source projects and still tracks the product objective", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  g <- .pc_gap(.pc_src2, .pc_pop(cor = c(WT.CRCL = 0.45)), NULL)
  expect_lt(g$worst, 1e-6)
})
