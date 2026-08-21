# A study generated from a published MODEL is not a sample. Its (E, V) are exact
# functions of theta_src, so the only random object in the chain is the estimate
# the source published and the covariance of our fit is the delta method through
# it. Reading `n` as precision instead makes the reported SE fall as exactly
# 1/sqrt(n) -- a factor the analyst chooses by typing a number.
#
# See validation/SPEC-se-paths.md for the full taxonomy.

TIMES_MS <- c(0.5, 1, 2, 4, 8, 12, 24)
DOSE_MS  <- 200

.ms_published <- function() {
  ini({ tcl <- log(5); tv <- log(50)
        eta.cl ~ 0.05; add.err <- 0.08 })
  model({ cl <- exp(tcl + eta.cl); v <- exp(tv)
          cp <- linCmt(); cp ~ add(add.err) })
}
.ms_fit <- function() {
  ini({ tcl <- log(4); tv <- log(45)
        eta.cl ~ 0.1; add.err <- 0.1 })
  model({ cl <- exp(tcl + eta.cl); v <- exp(tv)
          cp <- linCmt(); cp ~ add(add.err) })
}
# a FULL C_src over every parameter the source estimates, with correlation
.ms_C <- function() {
  nm <- c("tcl", "tv", "add.err", "eta.cl")
  sd <- c(0.080, 0.060, 0.004, 0.010)
  R  <- diag(4); R[1, 2] <- R[2, 1] <- -0.5; R[1, 4] <- R[4, 1] <- 0.3
  C  <- diag(sd) %*% R %*% diag(sd)
  dimnames(C) <- list(nm, nm)
  C
}
# the banded arm needs a model that actually READS the covariate it bands on
.ms_pub_wt <- function() {
  ini({ tcl <- log(5); tv <- log(50); bwt <- 0.75
        eta.cl ~ 0.05; add.err <- 0.08 })
  model({ cl <- exp(tcl + eta.cl) * (WT/70)^bwt; v <- exp(tv)
          cp <- linCmt(); cp ~ add(add.err) })
}
.ms_fit_wt <- function() {
  ini({ tcl <- log(4); tv <- log(45); bwt <- 0.4
        eta.cl ~ 0.1; add.err <- 0.1 })
  model({ cl <- exp(tcl + eta.cl) * (WT/70)^bwt; v <- exp(tv)
          cp <- linCmt(); cp ~ add(add.err) })
}
.ms_C_wt <- function() {
  nm <- c("tcl", "tv", "bwt", "add.err", "eta.cl")
  sd <- c(0.080, 0.060, 0.055, 0.004, 0.010)
  C  <- diag(sd^2); dimnames(C) <- list(nm, nm); C
}
.ms_gen <- function(n, C = .ms_C(), strat = NULL, J = NULL, times = TIMES_MS,
                    mod = .ms_published) {
  sp <- list(times = times, ev = rxode2::et(amt = DOSE_MS), n = n,
             model_cov = C)
  if (!is.null(strat)) {
    sp$cov_dist <- covDist(WT = c(mean = 78, sd = 16), dist = "lnorm")
    sp$stratify <- strat; sp$strata_nodes <- J
    sp$cov_range <- list(WT = c(50, 115))
  }
  suppressWarnings(suppressMessages(datagen(
    list(t1 = sp), model = mod,
    control = datagenControl(method = "gh", seed = 1L))))
}
.ms_run <- function(g, cm = "r,s", mod = .ms_fit) {
  suppressMessages(nlmixr2est::nlmixr2(mod, admData(), est = "adgh",
    control = adghControl(studies = g, print = 0L, cores = 2L, covMethod = cm)))
}

test_that("a lone model source reports the SOURCE's own uncertainty, at every n", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  # Where our model can reproduce the source's, the discrepancy is exactly zero
  # at theta = theta_src, so theta_hat = theta_src and G = I -- which makes
  # Var(theta_hat) = C_src EXACTLY. That is the whole claim, and it is what the
  # old path got wrong by a factor of sqrt(n).
  C  <- .ms_C()
  sd <- sqrt(diag(C))
  for (n in c(100, 1600)) {
    f  <- .ms_run(.ms_gen(n))
    pf <- f$parFixedDf
    se <- stats::setNames(pf[["SE"]], rownames(pf))
    for (k in c("tcl", "tv", "add.err"))
      expect_equal(se[[k]], sd[[k]], tolerance = 1e-4,
                   info = paste("n =", n, "parameter", k))
    expect_equal(sqrt(diag(f$cov))[["om.eta.cl"]], sd[["eta.cl"]],
                 tolerance = 1e-4, info = paste("n =", n, "omega"))
  }
})

test_that("the source's parameter CORRELATIONS propagate", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  # A scalar `n` carries no correlation at all. This is one of the things the
  # covariance route gives that no choice of `n` ever could.
  f  <- .ms_run(.ms_gen(400))
  cc <- stats::cov2cor(f$cov)
  expect_equal(cc["tcl", "tv"], -0.5, tolerance = 1e-3)
  expect_equal(cc["tcl", "om.eta.cl"], 0.3, tolerance = 1e-2)
})

test_that("a BANDED model source counts as one contribution, not J", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  # C_src is applied ONCE across the stacked strata. Per-stratum application
  # would give J independent copies, so raising the resolution would silently
  # buy confidence -- the covariance analogue of the sum(n_k) = n rule the
  # objective already keeps.
  se_at <- function(J) {
    g <- .ms_gen(400, C = .ms_C_wt(), strat = "WT", J = J, mod = .ms_pub_wt)
    f <- .ms_run(g, mod = .ms_fit_wt)
    stats::setNames(f$parFixedDf[["SE"]], rownames(f$parFixedDf))[["bwt"]]
  }
  s4 <- se_at(4L); s9 <- se_at(9L)
  expect_equal(s4, s9, tolerance = 1e-3)
  # ... and it is still the source's own, not something J-dependent
  expect_equal(s9, sqrt(.ms_C_wt()["bwt", "bwt"]), tolerance = 1e-3)
})

test_that("an INCOMPLETE model_cov is refused rather than understating the SE", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  # A parameter the source ESTIMATED but did not report a covariance for
  # contributes zero, which asserts the source knew it exactly and makes the
  # reported SE too SMALL -- the dangerous direction, and nothing about the
  # matrix looks wrong.
  C <- .ms_C()[c("tcl", "tv"), c("tcl", "tv"), drop = FALSE]
  # said at GENERATION time -- a warning raised later, from inside CalcCov, is
  # swallowed by the nlmixr2est stack (the drivers record the same lesson)
  expect_warning(
    g <- suppressMessages(datagen(
      list(t1 = list(times = TIMES_MS, ev = rxode2::et(amt = DOSE_MS),
                     n = 400, model_cov = C)),
      model = .ms_published, control = datagenControl(method = "gh", seed = 1L))),
    "ESTIMATES")
  # and the fit then refuses the sandwich rather than reporting a too-small SE
  f <- suppressWarnings(.ms_run(g))
  expect_s3_class(f, "admFit")
  expect_false(isTRUE(all.equal(
    stats::setNames(f$parFixedDf[["SE"]], rownames(f$parFixedDf))[["tcl"]],
    sqrt(C["tcl", "tcl"]), tolerance = 1e-3)))
})

test_that("the covariance route does not disturb a DATA source", {
  skip_on_cran(); skip_if_not_installed("rxode2")
  # No `model_cov` means no provenance means the existing data weight, byte for
  # byte. The change must be invisible to every fit that does not opt in.
  g  <- .ms_gen(400, C = NULL)
  a  <- .ms_run(g, cm = "r,s")
  b  <- .ms_run(g, cm = "r")
  expect_true(all(is.finite(sqrt(diag(a$cov)))))
  expect_true(all(is.finite(sqrt(diag(b$cov)))))
  # and a data source's SE DOES scale with n, which is correct for real patients
  g2 <- .ms_gen(1600, C = NULL)
  s1 <- stats::setNames(a$parFixedDf[["SE"]], rownames(a$parFixedDf))[["tcl"]]
  s2 <- .ms_run(g2, cm = "r")
  s2 <- stats::setNames(s2$parFixedDf[["SE"]], rownames(s2$parFixedDf))[["tcl"]]
  expect_gt(s1 / s2, 1.5)
})
