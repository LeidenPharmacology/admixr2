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
# not equally exact and saying so is the point. Measured against the product
# grid as strata_nodes goes 3, 4, 5:
#
#   untruncated   9.7e-06   2.9e-07   5.4e-08     keeps falling
#   truncated     1.5e-04   2.3e-05   3.1e-05     plateaus
#
# The untruncated reduction is exact in the integral, so only the finite tensor
# rules in the two coordinate systems differ and that gap closes. The truncated
# one goes through a cloud whose own quadrature error more node count cannot
# remove, which is why it settles instead.

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
.pc_gap <- function(src, pop, rng, sn = 4L, ana = .pc_ana,
                    covs = c("WT", "CRCL", "ALB")) {
  au  <- suppressMessages(rxode2::rxode2(ana))
  pin <- admixr2:::.admDriverPinfo(au, adghControl(studies = list(), print = 0L))
  rx  <- admixr2:::.admLoadModel(au); rxode2::rxLoad(rx)
  gr  <- admixr2:::.adghNodeGrid(5L, pin$n_eta)
  p0  <- admixr2:::.admBuildOptVec(pin)$p0
  mk  <- function() admStudies(a = admStudy(
    model = src, population = pop, n = 120L, dose = 200, times = c(1, 4),
    strata_nodes = sn, range = rng))
  gp <- suppressMessages(suppressWarnings(admixr2:::.admMaterialise(mk())))
  gj <- suppressMessages(suppressWarnings(admixr2:::.admMaterialise(
    mk(), analysis_covs = covs, analysis_ui = au)))
  ep <- admixr2:::.admBuildEvFull(gp); ej <- admixr2:::.admBuildEvFull(gj)
  worst <- 0
  for (d in c(-0.3, 0.3)) {
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

# FOUR covariates, so a rank-two source and a rank-one analysis model span
# three directions inside a four-dimensional grid and the reduction is real.
# With three covariates that sum reaches the full dimension and the design
# declines, which is what the earlier version of this test was measuring
# without saying so: it compared the product grid with itself.
.pc_pop4 <- function(...) admPopulation(
  WT = c(mean = 76, sd = 15), CRCL = c(mean = 92, sd = 22),
  ALB = c(mean = 40, sd = 5), AGE = c(mean = 55, sd = 12), ...)
.pc_src2 <- function() {
  ini({ tcl <- log(3.2); tv <- log(21); bwt <- .60; bcr <- .40
        balb <- .30; bage <- .20; eta.cl ~ .09; add.err <- .35 })
  model({ cl <- exp(tcl + eta.cl)*(WT/70)^bwt*(CRCL/95)^bcr
          v  <- exp(tv)*(ALB/40)^balb*(AGE/55)^bage
          cp <- linCmt(); cp ~ add(add.err) })
}
.pc_ana4 <- function() {
  ini({ tcl <- log(3.2); tv <- log(21)
        bwt <- fix(.52); bcr <- fix(.55); balb <- fix(.18); bage <- fix(.12)
        eta.cl ~ .09; add.err <- .6 })
  model({ cl <- exp(tcl + eta.cl)*(WT/70)^bwt*(CRCL/95)^bcr*
                (ALB/40)^balb*(AGE/55)^bage
          v <- exp(tv); cp <- linCmt(); cp ~ add(add.err) })
}
.pc_covs4 <- c("WT", "CRCL", "ALB", "AGE")

test_that("a rank-THREE span reduces a four-covariate grid and tracks it", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  ctl <- adghControl(studies = list(), print = 0L)
  su  <- suppressMessages(rxode2::rxode2(.pc_src2))
  au  <- suppressMessages(rxode2::rxode2(.pc_ana4))
  ds  <- admixr2:::.admCovDirections(su, admixr2:::.admDriverPinfo(su, ctl),
                                     .pc_pop4())
  da  <- admixr2:::.admCovReachable(au, admixr2:::.admDriverPinfo(au, ctl),
                                    .pc_pop4())
  expect_equal(ds$r, 2L)              # source reaches through cl AND v
  expect_equal(da$r, 1L)              # analysis loading is fixed
  expect_equal(ds$pc, 4L)

  g <- .pc_gap(.pc_src2, .pc_pop4(cor = c(WT.CRCL = 0.45)), NULL, sn = 3L,
               ana = .pc_ana4, covs = .pc_covs4)
  expect_equal(g$np, 3L^4L)
  expect_equal(g$nj, 3L^3L)           # the reduction actually happened
  expect_lt(g$worst, 1e-4)
})

test_that("a rank-THREE TRUNCATED span is declined, on its own criterion", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  # The cloud converges here (2.4e-04 to 1.4e-05 over j) but only 20 atoms are
  # retained at this resolution, and 20 atoms do not carry a three-dimensional
  # truncated measure: the recombination sits at ~3e-02 against a 1e-3
  # tolerance at every j. So the certificate refuses it, and refusing means the
  # product grid rather than a rule that is quietly 30x out.
  rng <- list(WT = c(55, 105), CRCL = c(60, 130),
              ALB = c(31, 49), AGE = c(35, 75))
  g <- .pc_gap(.pc_src2, .pc_pop4(cor = c(WT.CRCL = 0.45)), rng, sn = 3L,
               ana = .pc_ana4, covs = .pc_covs4)
  expect_equal(g$nj, g$np)            # declined
  expect_equal(g$worst, 0)            # and therefore identical, not merely close
})

test_that("more nodes close the untruncated gap and not the truncated one", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  # The two error sources are distinguishable by how they RESPOND to
  # strata_nodes, which is the claim the documentation makes. Untruncated, only
  # the finite rules in the two coordinate systems differ and the gap closes by
  # orders of magnitude. Truncated, the cloud's own quadrature error is there at
  # every node count and the gap settles instead of vanishing.
  pop <- .pc_pop(cor = c(WT.CRCL = 0.45))
  rng <- list(WT = c(55, 105), CRCL = c(60, 130), ALB = c(31, 49))
  u3 <- .pc_gap(.pc_src, pop, NULL, sn = 3L)$worst
  u5 <- .pc_gap(.pc_src, pop, NULL, sn = 5L)$worst
  t3 <- .pc_gap(.pc_src, pop, rng,  sn = 3L)$worst
  t5 <- .pc_gap(.pc_src, pop, rng,  sn = 5L)$worst
  expect_lt(u5, u3 / 20)          # closes fast
  expect_gt(t5, t3 / 20)          # does not
  expect_gt(t5, u5)               # and the floor is the cloud's, not the rule's
})

test_that("admission certifies the resolution the source ASKED for", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  # `.adm_strata_nodes` is stamped on the expanded studies after the design is
  # cut, so reading it at admission saw nothing and fell back to the default 9
  # whatever the user set. The cloud was then certified at nine nodes per
  # direction and .admCovStrata() recombined it at the requested number --
  # a different rule from the one that passed. At three nodes only six atoms
  # survive and the rule is ~5e-02 out, well past its own tolerance.
  #
  # Entered through the study spec, not by handing the helper a dotted field,
  # which is what hid this.
  ana <- suppressMessages(rxode2::rxode2(.pc_ana))
  rng <- list(WT = c(55, 105), CRCL = c(60, 130), ALB = c(31, 49))
  adm <- function(sn) admixr2:::.admStrataProj(
    list(cov_dist = .pc_pop(cor = c(WT.CRCL = 0.45)), cov_range = rng,
         strata_nodes = sn, .adm_ana_ui = ana), ana, c("WT", "CRCL", "ALB"))
  # too few retained atoms to certify, so no projection at all
  for (sn in c(1L, 2L, 3L)) expect_null(adm(sn))
  # and where it does certify, it is the requested resolution that was used
  sp <- adm(9L)
  expect_false(is.null(sp))
  expect_lt(sp$recomb_err, 1e-2)
  expect_equal(admixr2:::.admStrataNodes(list(strata_nodes = 3L)), 3L)
  expect_equal(admixr2:::.admStrataNodes(list()), admixr2:::.ADM_STRATA_NODES)

  # end to end: a request the certificate refuses keeps the product grid
  mk <- function() admStudies(a = admStudy(
    model = .pc_src, population = .pc_pop(cor = c(WT.CRCL = 0.45)), n = 120L,
    dose = 200, times = c(1, 4), strata_nodes = 3L, range = rng))
  g <- suppressMessages(suppressWarnings(admixr2:::.admMaterialise(
    mk(), analysis_covs = c("WT", "CRCL", "ALB"), analysis_ui = ana)))
  expect_equal(length(g), 3L^3L)
})

test_that("a branch that rewrites an exponent is declined, across the branch", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  skip_if_not_installed("Deriv")
  # `ba` is set unconditionally and then overwritten inside an `if`, so a
  # straight-line substitution certifies `.18` while the solver may run `1.18`.
  # The covariate is not in the branch, so .admCovInBranch() does not see it.
  # Crossing the branch AFTER materialisation is the only way to catch it.
  br <- function() {
    ini({ tcl <- log(3.2); tv <- log(21); b <- -1
          eta.cl ~ .09; add.err <- .6 })
    model({ ba <- .18
            if (b > 0) { ba <- 1.18 }
            cl <- exp(tcl + eta.cl)*(WT/70)^.52*(CRCL/95)^.55*(ALB/40)^ba
            v <- exp(tv); cp <- linCmt(); cp ~ add(add.err) }) }
  au  <- suppressMessages(rxode2::rxode2(br))
  expect_true(is.na(admixr2:::.admCovLoadingInvariant(
    au, c("WT", "CRCL", "ALB"), c("tcl", "tv", "b"))))

  pop <- .pc_pop(cor = c(WT.CRCL = 0.45))
  mk  <- function() admStudies(a = admStudy(
    model = .pc_src, population = pop, n = 120L, dose = 200,
    times = c(1, 4), strata_nodes = 3L))
  gp <- suppressMessages(suppressWarnings(admixr2:::.admMaterialise(mk())))
  gj <- suppressMessages(suppressWarnings(admixr2:::.admMaterialise(
    mk(), analysis_covs = c("WT", "CRCL", "ALB"), analysis_ui = au)))
  expect_equal(length(gj), length(gp))        # declined

  pin <- admixr2:::.admDriverPinfo(au, adghControl(studies = list(), print = 0L))
  rx  <- admixr2:::.admLoadModel(au); rxode2::rxLoad(rx)
  gr  <- admixr2:::.adghNodeGrid(5L, pin$n_eta)
  p0  <- admixr2:::.admBuildOptVec(pin)$p0
  ib  <- which(pin$struct_names == "b")
  expect_length(ib, 1L)
  ep <- admixr2:::.admBuildEvFull(gp); ej <- admixr2:::.admBuildEvFull(gj)
  for (bv in c(-1, 1)) {                      # either side of the branch
    p <- p0; p[ib] <- bv
    a <- suppressWarnings(admixr2:::.adghNLL(p, pin, ep, rx, "cp", gr, 1L))
    d <- suppressWarnings(admixr2:::.adghNLL(p, pin, ej, rx, "cp", gr, 1L))
    expect_true(is.finite(a))
    expect_equal(d, a, tolerance = 1e-12)
  }
})
