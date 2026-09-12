# Tier 1 -- covariate marginalisation helpers (no rxode2).

# gl/gh/taylor node quadrature was removed, and so were the "collapse" and "uq" paths (see the header of
# R/covariate.R). What remains is the general per-row / product-grid route, the shift route, and the refusal
# that keeps old node-style study lists from being fitted as something else.

# ---- routing -----------------------------------------------------------------

# # n_nodes is present because a real adgh pinfo always has it, and the covariate # ladder keys on it to tell
# adgh from the estimators that have no eta grid -- # driver.R hands those cov_integration = "on" by fallback,
# so the # string alone cannot make that distinction.
.cov_pinfo <- function(n_eta = 1L, cov_integration = NULL, n_nodes = 5L)
  list(n_eta = n_eta, n_nodes = n_nodes,
       eta_col_names = if (n_eta == 1L) "eta.cl" else c("eta.cl", "eta.v"),
       cov_integration = cov_integration,
       chol_diag = rep(TRUE, n_eta),
       chol_i = seq_len(n_eta), chol_j = seq_len(n_eta))

.cov_ui <- function(cov = "WT",
                    expr = list(quote(cl <- exp(tcl + tcov * WT + eta.cl))))
  list(allCovs = cov, lstExpr = expr)

test_that(".admCovCols adds ONLY declared covariates, never a blanket fill", {
  m <- matrix(1, 3L, 1L, dimnames = list(NULL, "tcl"))
  # `vb` is a hard-coded model constant and `lam` an estimated TBS lambda: both are model params, neither is a
  # covariate, so neither may be added.
  got <- admixr2:::.admCovCols(m, c("tcl", "wt", "vb", "lam"), list(wt = 70))
  expect_equal(colnames(got), c("tcl", "wt"))
  expect_equal(unname(got[, "wt"]), rep(70, 3L))

  # a covariate the model does not read is not added either
  expect_equal(colnames(admixr2:::.admCovCols(m, c("tcl", "vb"), list(wt = 70))),
               "tcl")
  # no covariates at all -> untouched
  expect_identical(admixr2:::.admCovCols(m, c("tcl", "wt"), NULL), m)
})

test_that(".admCheckCovariates accepts a supported mu-referenced covariate", {
  st <- list(a = list(cov_dist = list(WT = list(mu = 0, sd = 0.6))))
  expect_silent(admixr2:::.admCheckCovariates(.cov_ui(), .cov_pinfo(), st))
})

test_that(".admCheckCovariates is a no-op when no study declares cov_dist", {
  expect_silent(admixr2:::.admCheckCovariates(.cov_ui(), .cov_pinfo(),
                                              list(a = list())))
})

test_that(".admCheckCovariates caches the product grid for adgh", {
  st <- list(a = list(cov = list(WT = 0),
                      cov_dist = list(WT = list(mu = 0, sd = 0.6))))
  got <- admixr2:::.admCheckCovariates(.cov_ui(), .cov_pinfo(), st, "adgh")
  expect_equal(got$a$.adm_cov_grid,
               admixr2:::.admCovGrid(got$a$cov_dist, 7L))
})

test_that(".admCheckCovariates routes to the general path by default", {
  ok_st <- list(a = list(cov = list(WT = 0), cov_dist = list(WT = list(mu = 0, sd = 0.6))))
  path <- function(ui = .cov_ui(), pi = .cov_pinfo(), st = ok_st)
    admixr2:::.admCheckCovariates(ui, pi, st)$a$.adm_cov_path

  # the bare theta*COV product the retired collapse needed -- now "rows" like everything else. adgh reaches
  # the same moments through the shift path when cov_integration says so, and that route is admitted
  # numerically rather than by rxode2's muRefCovariateDataFrame (which recognises 1 of 8 realistic covariate
  # parameterisations).
  expect_identical(path(), "rows")
  expect_identical(path(st = list(a = list(cov = list(WT = 0),
      cov_dist = list(WT = list(values = c(0, 1), probs = c(.6, .4)))))), "rows")
  expect_identical(path(ui = .cov_ui(expr = list(
      quote(cl <- exp(tcl + tcov * WT + eta.cl)),
      quote(v  <- exp(tv) * (1 + 0.01 * WT))))), "rows")
  expect_identical(path(pi = .cov_pinfo(0L)), "rows")
  # no `cov` value given -> derived from the distribution
  expect_identical(path(st = list(a = list(
      cov_dist = list(WT = list(mu = 0, sd = 0.6))))), "rows")
})

test_that("a reduction that does not apply falls back quietly", {
  # A covariate on a parameter carrying no random effect. Nothing reduces it, and that is the ordinary case
  # rather than an event: no error, and nothing announced, because the caller asked for "on" and not for any
  # particular reduction. The study is integrated on the full grid.
  ui <- .cov_ui(expr = list(quote(v <- exp(tv) * WT^vwt)))
  st <- list(a = list(cov_dist = list(WT = list(mu = 70, sd = 10))))
  expect_silent(
    got <- admixr2:::.admCheckCovariates(
      ui, .cov_pinfo(cov_integration = "on"), st))
  expect_identical(got$a$.adm_cov_path, "rows")
  expect_null(got$a$.adm_cov_joint)
})

test_that(".admCovMeanOf gives the solve value each path needs", {
  expect_equal(admixr2:::.admCovMeanOf(list(mu = 70, sd = 10)), 70)
  expect_equal(admixr2:::.admCovMeanOf(list(meanlog = log(70), sdlog = 0.2)),
               exp(log(70) + 0.2^2 / 2))
  expect_equal(admixr2:::.admCovMeanOf(list(values = c(0, 1), probs = c(0.6, 0.4))), 0.4)
  expect_equal(admixr2:::.admCovMeanOf(list(quantile = function(p) qnorm(p, 5, 1))),
               5, tolerance = 1e-6)
})

test_that("a study omitting `cov` has it filled in from `cov_dist`", {
  st <- admixr2:::.admCheckCovariates(
    .cov_ui(), .cov_pinfo(),
    list(a = list(cov_dist = list(WT = list(mu = 70, sd = 10)))))
  expect_equal(st$a$cov$WT, 70)
})

test_that(".admCheckCovariates still errors on genuinely unsupportable input", {
  ok_st <- list(a = list(cov = list(WT = 0), cov_dist = list(WT = list(mu = 0, sd = 0.6))))

  # covariate the model never reads -- almost always a typo
  expect_error(
    admixr2:::.admCheckCovariates(.cov_ui(cov = "AGE"), .cov_pinfo(), ok_st),
    "which the model never reads")

  # distributions we cannot draw from
  for (spec in list(list(mu = 0), list(mu = 0, sd = 0), list(mu = 0, sd = NA_real_),
                    list(values = numeric(0))))
    expect_error(
      admixr2:::.admCheckCovariates(.cov_ui(), .cov_pinfo(),
                                    list(a = list(cov = list(WT = 0),
                                                  cov_dist = list(WT = spec)))),
      "not a supported distribution")
})

test_that("an unidentifiable covariate coefficient is warned about", {
  # One population, or several with the SAME covariate distribution, cannot
  # identify a coefficient that shares its argument with a random effect: the
  # likelihood is exactly flat along
  #   theta' = theta + (b-b')*mu_a ,  omega'^2 = omega^2 + (b^2-b'^2)*sd_a^2
  # (verified bit-identical across b = 0.40 .. 1.10). Only between-study
  # variation in the covariate DISTRIBUTION breaks it.
  ui <- .cov_ui()   # cl <- exp(tcl + tcov*WT + eta.cl)
  one  <- list(a = list(cov_dist = list(WT = list(mu = 0, sd = 0.6))))
  same <- list(a = list(cov_dist = list(WT = list(mu = 0, sd = 0.6))),
               b = list(cov_dist = list(WT = list(mu = 0, sd = 0.6))))
  diff_mu <- list(a = list(cov_dist = list(WT = list(mu = 0.0, sd = 0.6))),
                  b = list(cov_dist = list(WT = list(mu = 0.4, sd = 0.6))))
  diff_sd <- list(a = list(cov_dist = list(WT = list(mu = 0, sd = 0.6))),
                  b = list(cov_dist = list(WT = list(mu = 0, sd = 0.9))))

  expect_warning(admixr2:::.admWarnCovIdentifiability(ui, .cov_pinfo(), one),
                 "not identifiable")
  expect_warning(admixr2:::.admWarnCovIdentifiability(ui, .cov_pinfo(), same),
                 "not identifiable")
  # differing MEANS break the first equation, differing SPREADS the second
  expect_silent(admixr2:::.admWarnCovIdentifiability(ui, .cov_pinfo(), diff_mu))
  expect_silent(admixr2:::.admWarnCovIdentifiability(ui, .cov_pinfo(), diff_sd))

  # a covariate on a parameter with NO random effect is NOT on the ridge: there is no omega for its variance
  # to hide in, so one population is enough
  ui2 <- .cov_ui(expr = list(quote(v <- exp(tv) * WT^vwt)))
  expect_silent(admixr2:::.admWarnCovIdentifiability(ui2, .cov_pinfo(), one))
})

test_that("conditioned strata identify the coefficient within one source", {
  # A stratum is a covariate distribution at zero spread, and supplies its own equation in (theta, gamma)
  # exactly as a second study's distribution does. Reading `cov_dist` alone called this unidentified -- a
  # false alarm on the one design that identifies gamma without a between-source contrast.
  ui <- .cov_ui()   # cl <- exp(tcl + tcov*WT + eta.cl)
  d  <- list(WT = list(mu = 0, sd = 0.6))

  two_strata <- list(a = list(cov = list(WT = -0.4)),
                     b = list(cov = list(WT =  0.4)))
  expect_silent(admixr2:::.admWarnCovIdentifiability(ui, .cov_pinfo(), two_strata))

  # mixed: one marginalised source plus strata from another
  mixed <- list(a = list(cov_dist = d),
                b = list(cov = list(WT = -0.4)),
                c = list(cov = list(WT =  0.4)))
  expect_silent(admixr2:::.admWarnCovIdentifiability(ui, .cov_pinfo(), mixed))

  # a SINGLE stratum is still one equation in two unknowns -> still warns
  expect_warning(admixr2:::.admWarnCovIdentifiability(
    ui, .cov_pinfo(), list(a = list(cov = list(WT = 0.4)))), "not identifiable")
})

test_that("a covariate entering a second parameter is off the ridge", {
  # eq:ridge holds only while the model sees gamma*a + omega*b. A covariate read by a second assignment
  # restores the separate dependence on (a, b), so it is identified from one population and must not be warned
  # about.
  ui <- .cov_ui(expr = list(quote(v <- exp(tv + vwt * WT))))
  expect_null(admixr2:::.admCovParamEta(ui, "WT", "eta.cl"))
  expect_silent(admixr2:::.admWarnCovIdentifiability(
    ui, .cov_pinfo(),
    list(a = list(cov_dist = list(WT = list(mu = 0, sd = 0.6))))))
})

test_that(".admCovSdOf reports the spread the ridge depends on", {
  expect_equal(admixr2:::.admCovSdOf(list(mu = 70, sd = 10)), 10)
  expect_equal(admixr2:::.admCovSdOf(list(meanlog = log(70), sdlog = 0.2)), 0.2)
  expect_equal(admixr2:::.admCovSdOf(list(values = c(0, 1), probs = c(0.5, 0.5))), 0.5)
})

test_that("the friendly cov_dist grammar is EXACTLY the hand-written one", {
  # mean/sd on the natural scale, and `cor`, exist so a user can transcribe a baseline-characteristics table
  # directly. They must therefore produce the identical distribution to the canonical spelling -- a
  # convenience that shifted the covariate distribution even slightly would move every estimate while looking
  # like a formatting choice.
  ml <- log(72^2 / sqrt(16^2 + 72^2)); sl <- sqrt(log(1 + 16^2 / 72^2))
  friendly <- admixr2:::.admCovDistCanon(
    list(WT = list(mean = 72, sd = 16, dist = "lnorm")))
  expect_equal(friendly$WT$meanlog, ml)
  expect_equal(friendly$WT$sdlog,   sl)
  # the natural-scale moments are matched, which is the whole point
  q <- stats::qlnorm(stats::ppoints(200000), friendly$WT$meanlog,
                     friendly$WT$sdlog)
  expect_equal(mean(q), 72, tolerance = 1e-3)
  expect_equal(stats::sd(q), 16, tolerance = 5e-3)

  # `cor` must reproduce a hand-written Gaussian copula draw for draw
  rho <- 0.6
  mc  <- log(90^2 / sqrt(25^2 + 90^2)); sc <- sqrt(log(1 + 25^2 / 90^2))
  # lognormal margins here on purpose: weight and creatinine clearance are positive, right-skewed
  # physiological quantities, and a Gaussian copula on lognormal margins has correlation exactly rho on the
  # LOG scale
  cd  <- admixr2:::.admCovDistCanon(
    list(WT = list(mean = 72, sd = 16, dist = "lnorm"),
         CRCL = list(mean = 90, sd = 25, dist = "lnorm"), cor = rho))
  expect_true(is.function(cd$joint))
  expect_null(cd$cor)                       # consumed, so nothing else sees it
  u <- randtoolbox::sobol(4000L, dim = 2)
  u <- pmin(pmax(u, 1e-12), 1 - 1e-12)
  got <- cd$joint(u)
  z   <- stats::qnorm(u)
  z2  <- rho * z[, 1] + sqrt(1 - rho^2) * z[, 2]
  want <- cbind(WT   = stats::qlnorm(stats::pnorm(z[, 1]), ml, sl),
                CRCL = stats::qlnorm(stats::pnorm(z2),     mc, sc))
  expect_equal(unname(got), unname(want), tolerance = 1e-8)
  expect_identical(colnames(got), c("WT", "CRCL"))
  # a Gaussian copula on lognormal margins has correlation exactly rho on the LOG scale, which is the check
  # that `cor` means what a reader expects
  expect_equal(stats::cor(log(got[, 1]), log(got[, 2])), rho, tolerance = 1e-2)

  # idempotent: applying it twice must not re-expand anything
  expect_equal(admixr2:::.admCovDistCanon(cd)$WT, cd$WT)

  # an explicit canonical spelling always wins over the shorthand
  keep <- admixr2:::.admCovDistCanon(
    list(WT = list(mean = 999, meanlog = ml, sdlog = sl)))
  expect_equal(keep$WT$meanlog, ml)
})

test_that("the friendly grammar refuses what it cannot represent", {
  expect_error(admixr2:::.admCovDistCanon(list(WT = list(mean = 72))), "sd")
  # a negative mean is legal for the NORMAL default and impossible for lnorm
  expect_silent(admixr2:::.admCovDistCanon(list(WT = list(mean = -1, sd = 2))))
  expect_error(admixr2:::.admCovDistCanon(
    list(WT = list(mean = -1, sd = 2, dist = "lnorm"))), "lognormal")
  expect_error(admixr2:::.admCovDistCanon(list(WT = list(mean = 72, sd = 16),
                                               cor = 0.5)),
               "at least two")
  # a scalar cor with three covariates is ambiguous, not a broadcast
  three <- list(A = list(mean = 1, sd = .1), B = list(mean = 1, sd = .1),
                C = list(mean = 1, sd = .1), cor = 0.5)
  expect_error(admixr2:::.admCovDistCanon(three), "two covariates")
  # a correlation matrix that is not a correlation matrix
  bad <- list(A = list(mean = 1, sd = .1), B = list(mean = 1, sd = .1),
              cor = matrix(c(1, 1.4, 1.4, 1), 2, 2))
  expect_error(admixr2:::.admCovDistCanon(bad), "positive definite")
  # a named matrix given in the other order is REORDERED, not misapplied
  R <- matrix(c(1, 0.7, 0.7, 1), 2, 2,
              dimnames = list(c("B", "A"), c("B", "A")))
  ok <- admixr2:::.admCovDistCanon(
    list(A = list(mean = 10, sd = 1), B = list(mean = 100, sd = 20), cor = R))
  expect_true(is.function(ok$joint))
  expect_identical(colnames(ok$joint(matrix(c(.2, .8, .4, .6), 2, 2))),
                   c("A", "B"))
})

# Sparse-grid covariate integration (cov_integration = "sparse")

# The reference is EXACT nested Gauss-Hermite over (covariate, eta) on an analytic one-compartment solution,
# computed here in plain R. Nothing is pinned against the package's own output: the failure mode of a moment
# expansion is a finite, plausible, biased number, so a self-comparison catches nothing.

.tay_gh <- function(k) {
  i <- seq_len(k - 1L)
  J <- matrix(0, k, k); J[cbind(i, i + 1L)] <- sqrt(i); J[cbind(i + 1L, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE); o <- order(e$values)
  list(x = e$values[o], w = e$vectors[1L, o]^2)
}
.tay_TCL <- log(1); .tay_TV <- log(10); .tay_OM <- 0.30; .tay_DOSE <- 100
.tay_TIMES <- c(0.5, 1, 1.5, 2, 3, 4, 5, 6, 8)
.tay_QE <- .tay_gh(21L)

# f(a, eta) for a 1-cmt bolus with cl = exp(tcl + tcov*a + eta), v fixed
.tay_conc <- function(cl)
  .tay_DOSE / exp(.tay_TV) * exp(outer(-cl / exp(.tay_TV), .tay_TIMES))

# EXACT marginal structural moments over (a, eta), by nested quadrature
.tay_exact <- function(tcov, mu_a, sd_a, k = 21L) {
  QA <- .tay_gh(k); m1 <- 0; M2 <- 0
  for (ia in seq_along(QA$x)) {
    a <- mu_a + sd_a * QA$x[ia]
    Y <- .tay_conc(exp(.tay_TCL + tcov * a + .tay_OM * .tay_QE$x))
    for (ib in seq_along(.tay_QE$x)) {
      w <- QA$w[ia] * .tay_QE$w[ib]; y <- Y[ib, ]
      m1 <- m1 + w * y; M2 <- M2 + w * outer(y, y)
    }
  }
  list(E = m1, V = M2 - outer(m1, m1))
}

# Conditional moments at ONE covariate value -- the ecological plug-in, and the thing the expansion
# differences.
.tay_cond <- function(a, tcov) {
  Y  <- .tay_conc(exp(.tay_TCL + tcov * a + .tay_OM * .tay_QE$x))
  mu <- as.numeric(crossprod(.tay_QE$w, Y)); Yc <- sweep(Y, 2L, mu)
  list(E = mu, V = t(Yc) %*% (Yc * .tay_QE$w))
}

.tay_relerr <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-12))

test_that(".adghStructMoments is BIT-IDENTICAL to the pooled formulas without a design", {
  # Every route passes through this one two-line form now that the Taylor design's second branch is gone --
  # this is the whole of what it computes.
  set.seed(11)
  cp <- matrix(stats::rnorm(40L * 6L, 5, 2), 40L, 6L)
  W  <- stats::runif(40L); W <- W / sum(W)
  got <- admixr2:::.adghStructMoments(cp, W)
  mu  <- as.numeric(crossprod(W, cp))
  cpc <- sweep(cp, 2L, mu)
  expect_identical(got$mu,  mu)
  expect_identical(got$cpc, cpc)
  expect_identical(got$V,   crossprod(cpc, W * cpc))
  expect_null(got$dE)
})

# Dependent covariates on the DETERMINISTIC (quadrature / sparse) paths

test_that(".admCovGrid integrates a DEPENDENT distribution on the u-space grid", {
  # A copula maps INDEPENDENT uniforms to dependent values, so a product rule in u is exact whatever the
  # dependence -- the weights genuinely factorise there. adgh used to refuse a `joint` spec on the grounds
  # that its product grid assumed independence; that is true of a grid over covariate MARGINS and false of a
  # grid over u.
  cd <- admixr2:::.admCovDistCanon(
    list(WT  = list(mu = 70, sd = 10, dist = "normal"),
         AGE = list(mu = 50, sd = 12, dist = "normal"), cor = 0.7))
  for (nn in c(5L, 9L)) {
    g <- admixr2:::.admCovGrid(cd, nn)
    expect_identical(dim(g$X), c(nn * nn, 2L))
    expect_identical(colnames(g$X), c("WT", "AGE"))
    expect_equal(sum(g$W), 1)
    m  <- as.numeric(crossprod(g$W, g$X))
    Xc <- sweep(g$X, 2L, m); S <- crossprod(Xc, g$W * Xc)
    # normal margins + Gaussian copula: the closed form is exact, so this is a real check and not a comparison
    # against another approximation
    expect_equal(unname(m), c(70, 50), tolerance = 1e-6)
    expect_equal(unname(sqrt(diag(S))), c(10, 12), tolerance = 1e-6)
    expect_equal(S[1, 2] / sqrt(S[1, 1] * S[2, 2]), 0.7, tolerance = 1e-6)
  }
})

test_that("the grid and the per-subject sampler see the SAME distribution", {
  # admc draws rows, adgh builds a grid. If these disagreed, the two estimators would fit different data and
  # only a side-by-side run would reveal it.
  cd <- admixr2:::.admCovDistCanon(
    list(WT = list(mean = 72, sd = 16), CRCL = list(mean = 90, sd = 25),
         cor = 0.6))
  g  <- admixr2:::.admCovGrid(cd, 21L)
  mg <- as.numeric(crossprod(g$W, g$X))
  Xg <- sweep(g$X, 2L, mg); Sg <- crossprod(Xg, g$W * Xg)
  rr <- admixr2:::.admCovRowsFor(cd, 200000L, 0L)
  mr <- colMeans(rr); Sr <- crossprod(sweep(rr, 2L, mr)) / nrow(rr)
  expect_equal(unname(mg), unname(mr), tolerance = 1e-3)
  expect_equal(unname(sqrt(diag(Sg))), unname(sqrt(diag(Sr))), tolerance = 1e-3)
  expect_equal(Sg[1, 2] / sqrt(Sg[1, 1] * Sg[2, 2]),
               Sr[1, 2] / sqrt(Sr[1, 1] * Sr[2, 2]), tolerance = 1e-3)
})

# The joint covariance a cov_dist actually induces, measured rather than read.

# This used to be .admCovDistMoments in R/, but the sparse grid works from the margins and latentR directly
# and nothing in the package needed joint moments any more -- so it lives here, where its only two users are.
# Measuring is the point: these tests check that the three spellings of a correlation produce the SAME
# distribution, which reading the declaration back would not test at all.
.cd_sigma <- function(cov_dist, n = 8192L) {
  cd  <- admixr2:::.admCovDistCanon(cov_dist)
  nms <- admixr2:::.admCovSpecNames(cd)
  jf  <- cd[["joint"]]
  if (!is.function(jf)) {
    g <- admixr2:::.admCovGrid(cd, 21L)
    m <- as.numeric(crossprod(g$W, g$X))
    return(crossprod(sweep(g$X, 2L, m), g$W * sweep(g$X, 2L, m)))
  }
  u <- randtoolbox::sobol(n, dim = length(nms))
  if (!is.matrix(u)) u <- matrix(u, nrow = n)
  u <- pmin(pmax(u, .Machine$double.eps), 1 - .Machine$double.eps)
  colnames(u) <- nms
  X <- admixr2:::.admCovJointEval(jf, u, nms)
  S <- stats::cov(X) * (n - 1) / n
  dimnames(S) <- list(nms, nms)
  S
}

test_that("`cor`, `rho` and `Sigma` are ONE statement, honoured by every path", {
  # They used to diverge: `rho` built a copula for the retired collapse and a DIAGONAL grid for everything
  # else, so the correlation was present in one path and silently absent in the other. All three spellings are
  # still accepted -- published specs are written every which way -- so all three must still land on the same
  # latent correlation.
  base <- list(A = list(mu = 0, sd = 2, dist = "normal"),
               B = list(mu = 0, sd = 3, dist = "normal"))
  sig <- function(x) .cd_sigma(x)
  s_cor <- sig(c(base, list(cor   = 0.5)))
  expect_equal(sig(c(base, list(rho   = 0.5))), s_cor, tolerance = 0)
  expect_equal(sig(c(base, list(Sigma = matrix(c(4, 3, 3, 9), 2, 2)))), s_cor,
               tolerance = 0)
  expect_equal(unname(s_cor), matrix(c(4, 3, 3, 9), 2, 2), tolerance = 1e-2)
  # ... and all three canonicalise to the same latent correlation matrix, which is the single object every
  # consumer (the joint sampler, the product grid, the sparse grid, the shift grid) reads.
  lr <- function(x) admixr2:::.admCovDistCanon(x)[["latentR"]]
  expect_equal(lr(c(base, list(rho = 0.5))), lr(c(base, list(cor = 0.5))))
  expect_equal(lr(c(base, list(Sigma = matrix(c(4, 3, 3, 9), 2, 2)))),
               lr(c(base, list(cor = 0.5))))
})

test_that("metadata keys are never mistaken for covariates", {
  # each of these used to be a different internal error: `$ operator is invalid for atomic vectors` on rho,
  # `object of type 'closure' is not subsettable` on the joint sampler.
  cd <- list(A = list(mu = 0, sd = 1), B = list(mu = 0, sd = 1), rho = 0.3)
  expect_identical(admixr2:::.admCovSpecNames(cd), c("A", "B"))
  expect_silent(admixr2:::.admCovGrid(admixr2:::.admCovDistCanon(cd), 3L))
  expect_identical(
    colnames(admixr2:::.admCovGrid(admixr2:::.admCovDistCanon(cd), 3L)$X),
    c("A", "B"))
})

test_that("a NAMED correlation matrix is ordered to the declared covariates", {
  # The reorder is the only thing standing between a user writing their cor matrix in a different order and
  # admixr2 correlating the wrong PAIR -- which is finite, plausible and silent. It briefly became dead code
  # by sitting after a `dimnames(R) <- NULL`.
  R  <- matrix(c(1, 0.8, 0.8, 1), 2, 2,
               dimnames = list(c("B", "A"), c("B", "A")))
  cd <- admixr2:::.admCovDistCanon(
    list(A = list(mu = 0, sd = 1, dist = "normal"),
         B = list(mu = 0, sd = 5, dist = "normal"), cor = R))
  S <- .cd_sigma(cd)
  # variances must land on the covariate they were DECLARED for
  expect_equal(unname(S[1, 1]), 1, tolerance = 5e-3)
  expect_equal(unname(S[2, 2]), 25, tolerance = 5e-3)
  expect_equal(unname(S[1, 2] / sqrt(S[1, 1] * S[2, 2])), 0.8, tolerance = 5e-3)
})

# Covariate STRATA -- the per-covariate stratify/marginalise split

test_that("a `cor` spec conditions on exact POINTS, matching closed form", {
  # admixr2 built this copula, so the conditional is closed form on the latent scale -- no pool, no importance
  # weighting, nothing to collapse. It is also the cheapest route at any useful stratum count.
  RHO <- 0.7; MU1 <- 70; S1 <- 10; MU2 <- 90; S2 <- 20
  cd <- covDist(WT = c(mu = MU1, sd = S1), CRCL = c(mu = MU2, sd = S2),
                cor = RHO)
  expect_false(is.null(cd[["latentR"]]))    # the latent correlation, stashed
  st <- admixr2:::.admCovStrata(cd, stratify = "WT", n_nodes = 4L)
  expect_length(st, 4L)
  expect_equal(sum(vapply(st, `[[`, 0, "weight")), 1)
  for (s in st) {
    X <- covDraw(s$cov_dist, n = 25000L)
    # the closed-form Gaussian conditional, at this stratum's own WT
    expect_equal(mean(X[, "CRCL"]),
                 MU2 + RHO * S2 / S1 * (s$cov$WT - MU1), tolerance = 2e-2)
    expect_equal(sd(X[, "CRCL"]), S2 * sqrt(1 - RHO^2), tolerance = 0.12)
    # the stratified covariate is a POINT, so it has no spread within a stratum
    expect_lt(sd(X[, "WT"]), 1e-8)
  }
  # the conditional mean must MOVE across strata
  em <- vapply(st, function(s) mean(covDraw(s$cov_dist, n = 20000L)[, "CRCL"]), 0)
  expect_gt(diff(range(em)), S2)
})

test_that("without a density, strata fall back to BANDS that partition", {
  # A `joint` sampler with no density cannot be importance-weighted, so the strata are equiprobable bins of
  # the pool instead -- correct, but coarser.
  cl <- function(x) pmin(pmax(x, 1e-12), 1 - 1e-12)
  gauss <- function(u) {
    z <- stats::qnorm(cl(u))
    z2 <- 0.7 * z[, 1] + sqrt(1 - 0.49) * z[, 2]
    cbind(WT = stats::qnorm(cl(stats::pnorm(z[, 1])), 70, 10),
          CRCL = stats::qnorm(cl(stats::pnorm(z2)), 90, 20))
  }
  cd <- covDist(WT = list(quantile = function(u) stats::qnorm(u, 70, 10)),
                CRCL = list(quantile = function(u) stats::qnorm(u, 90, 20)),
                joint = gauss)
  expect_null(cd$density)
  st <- admixr2:::.admCovStrata(cd, stratify = "WT", n_nodes = 4L)
  expect_length(st, 4L)
  # equiprobable BANDS, not quadrature weights
  expect_equal(vapply(st, `[[`, 0, "weight"), rep(0.25, 4L), tolerance = 1e-6)
  # a band is a range, so it carries the stratified covariate too
  expect_true(all(vapply(st, function(s) "WT" %in% names(s$cov_dist), TRUE)))
  # and it still conditions: checked against brute force on the same bins
  big <- covDraw(cd, n = 300000L)
  qb <- stats::quantile(big[, "WT"], seq(0, 1, 0.25), names = FALSE)
  qb[1L] <- -Inf; qb[length(qb)] <- Inf
  b <- cut(big[, "WT"], qb, labels = FALSE)
  got <- vapply(st, function(s)
    mean(covDraw(s$cov_dist, n = 25000L)[, "CRCL"]), 0)
  want <- vapply(sort(unique(b)), function(k) mean(big[b == k, "CRCL"]), 0)
  expect_equal(sort(got), sort(want), tolerance = 2e-2)
})

test_that("with INDEPENDENT covariates the conditioning is a no-op", {
  cd <- covDist(WT = c(mu = 70, sd = 10),
                CRCL = c(mu = 90, sd = 20))
  st <- admixr2:::.admCovStrata(cd, stratify = "WT", n_nodes = 4L)
  em <- vapply(st, function(s) mean(covDraw(s$cov_dist, n = 30000L)[, "CRCL"]), 0)
  expect_equal(em, rep(90, 4L), tolerance = 2e-2)
})

test_that("an opaque sampler is conditioned correctly when the stratified
           covariate is NOT the head of its cascade", {
  # This is the case the u-space route got silently wrong: on a cascade running AGE -> WT -> CRCL, stratifying
  # on WT left AGE at its unconditional mean in every stratum, against a true 50.1 to 60.4. Banding an opaque
  # sampler has to bin its OUTPUT; pinning its input uniforms bands only the cascade head.

  # The cascade is built by hand rather than fitted, so this needs no copula package. The conditioning is
  # ASYMMETRIC -- a squared term -- so no correlation matrix reproduces it and the sampler really is opaque.
  ml <- c(log(72), log(90), log(54)); sl <- c(0.22, 0.26, 0.20)
  cl <- function(x) pmin(pmax(x, 1e-12), 1 - 1e-12)
  jf <- function(u) {
    u  <- cl(u)
    zA <- stats::qnorm(u[, 3L])
    # A PURE FUNCTION OF `u`, row by row: admixr2 calls this with a node grid as well as with a sample, so
    # anything estimated from the batch (an sd, a rank) would make the same u give different covariates in the
    # two.
    zW <- 0.62 * zA + 0.18 * (zA^2 - 1) +
            sqrt(1 - 0.62^2) * stats::qnorm(u[, 1L])
    zC <- 0.55 * zW + sqrt(1 - 0.55^2) * stats::qnorm(u[, 2L])
    o  <- cbind(stats::qlnorm(stats::pnorm(zW), ml[1L], sl[1L]),
                stats::qlnorm(stats::pnorm(zC), ml[2L], sl[2L]),
                stats::qlnorm(stats::pnorm(zA), ml[3L], sl[3L]))
    colnames(o) <- c("WT", "CRCL", "AGE"); o
  }
  cd <- covDist(WT = list(quantile = function(u) stats::qlnorm(u, ml[1], sl[1])),
                CRCL = list(quantile = function(u) stats::qlnorm(u, ml[2], sl[2])),
                AGE = list(quantile = function(u) stats::qlnorm(u, ml[3], sl[3])),
                joint = jf)
  st <- admixr2:::.admCovStrata(cd, stratify = "WT", n_nodes = 4L)
  got <- vapply(st, function(s) {
    X <- covDraw(s$cov_dist, n = 30000L)
    c(mean(X[, "WT"]), mean(X[, "CRCL"]), mean(X[, "AGE"])) }, numeric(3))
  got <- got[, order(got[1L, ]), drop = FALSE]
  big <- covDraw(cd, n = 300000L)
  qb <- stats::quantile(big[, "WT"], seq(0, 1, 0.25), names = FALSE)
  qb[1L] <- -Inf; qb[length(qb)] <- Inf
  b <- cut(big[, "WT"], qb, labels = FALSE)
  want <- vapply(sort(unique(b)), function(k)
    c(mean(big[b == k, "CRCL"]), mean(big[b == k, "AGE"])), numeric(2))
  expect_equal(got[2L, ], want[1L, ], tolerance = 2e-2)
  expect_equal(got[3L, ], want[2L, ], tolerance = 2e-2)
  # AGE is correlated with WT, so its conditional mean MUST move
  expect_gt(diff(range(got[3L, ])), 0.3 * sd(big[, "AGE"]))
})

test_that("a DISCRETE stratified covariate enumerates its levels exactly", {
  cd <- list(SEX = list(values = c(0, 1), probs = c(0.6, 0.4)),
             WT  = list(mu = 70, sd = 10))
  st <- admixr2:::.admCovStrata(cd, stratify = "SEX", n_nodes = 7L)
  expect_length(st, 2L)                     # levels, NOT n_nodes
  expect_equal(sort(vapply(st, function(s) s$cov$SEX, 0)), c(0, 1))
  expect_equal(sort(vapply(st, `[[`, 0, "weight")), c(0.4, 0.6),
               tolerance = 1e-2)
})

test_that(".admExpandStrata turns one study into ordinary per-stratum studies", {
  s <- list(src = list(n = 300L, times = c(1, 2), ev = "EV",
                       cov_dist = list(WT = list(mu = 70, sd = 10),
                                       CRCL = list(mu = 90, sd = 20)),
                       stratify = "WT", strata_nodes = 4L))
  ex <- admixr2:::.admExpandStrata(s, names(s))
  expect_length(ex$studies, 4L)
  expect_identical(ex$names, paste0("src_s", 1:4))
  # n_k are quadrature weights times n, and must still sum to n
  expect_equal(sum(vapply(ex$studies, `[[`, 0, "n")), 300)
  # every stratum is a PLAIN study -- the flags are consumed, not passed on
  expect_true(all(vapply(ex$studies, function(x) is.null(x$stratify), TRUE)))
  expect_true(all(vapply(ex$studies, function(x) is.null(x$strata_nodes), TRUE)))
  # and each carries its own pinned value plus the conditional remainder
  expect_true(all(vapply(ex$studies, function(x) !is.null(x$cov$WT), TRUE)))
  expect_true(all(vapply(ex$studies, function(x) "CRCL" %in% names(x$cov_dist), TRUE)))
  # a study without `stratify` passes through untouched
  s2 <- list(a = list(n = 10L), b = list(n = 20L))
  expect_identical(admixr2:::.admExpandStrata(s2, names(s2))$studies, s2)
})

test_that("stratify refuses what it cannot cut", {
  mk <- function(...) list(src = utils::modifyList(
    list(n = 300L, times = c(1, 2), stratify = "WT"), list(...)))
  expect_error(admixr2:::.admExpandStrata(mk(), "src"), "no `cov_dist`")
  expect_error(admixr2:::.admExpandStrata(
    mk(cov_dist = list(WT = list(mu = 0, sd = 1)), n = NULL), "src"),
    "no positive.*`n`")
  expect_error(admixr2:::.admExpandStrata(
    mk(cov_dist = list(WT = list(mu = 0, sd = 1)),
       observations = list(a = 1)), "src"), "observations")
  expect_error(admixr2:::.admCovStrata(
    list(WT = list(mu = 0, sd = 1)), stratify = "AGE"), "does not declare")
})

test_that("`stratify = TRUE` derives the split from the SOURCE model", {
  skip_if_not_installed("rxode2")
  # The model is the only thing that knows which covariates a source conditioned on. Making the user restate
  # it duplicates information already present, and getting it wrong is the fabricated null contrast.
  mA <- function() {
    ini({lcl <- log(3); lv <- log(20); bwt <- 0.75; eta.cl ~ 0.09; a <- 0.1})
    model({cl <- exp(lcl + eta.cl) * (WT / 70)^bwt; v <- exp(lv)
           d/dt(ce) <- -cl / v * ce; cp <- ce / v; cp ~ add(a)})
  }
  mB <- function() {
    ini({lcl <- log(3); lv <- log(20); bwt <- 0.75; bcr <- 0.4
         eta.cl ~ 0.09; a <- 0.1})
    model({cl <- exp(lcl + eta.cl) * (WT / 70)^bwt * (CRCL / 90)^bcr
           v <- exp(lv); d/dt(ce) <- -cl / v * ce; cp <- ce / v; cp ~ add(a)})
  }
  cd <- list(WT = list(mean = 72, sd = 15), CRCL = list(mean = 90, sd = 22))
  mk <- function(m) list(s = list(model = m, n = 200L, cov_dist = cd,
                                  stratify = TRUE, strata_nodes = 3L))
  # A fitted WT only -> 3 strata, CRCL left to be marginalised
  eA <- admixr2:::.admExpandStrata(mk(mA), "s")
  expect_length(eA$studies, 3L)
  expect_identical(names(eA$studies[[1L]]$cov), "WT")
  # a stratum is a range, so it carries ALL covariates, including its own
  expect_setequal(admixr2:::.admCovSpecNames(eA$studies[[1L]]$cov_dist),
                  c("WT", "CRCL"))
  # B fitted both -> 9 strata, nothing left over
  eB <- admixr2:::.admExpandStrata(mk(mB), "s")
  expect_length(eB$studies, 9L)
  expect_setequal(names(eB$studies[[1L]]$cov), c("WT", "CRCL"))
  expect_setequal(admixr2:::.admCovSpecNames(eB$studies[[1L]]$cov_dist),
                  c("WT", "CRCL"))
  # both partition the same n
  expect_equal(sum(vapply(eA$studies, `[[`, 0, "n")), 200)
  expect_equal(sum(vapply(eB$studies, `[[`, 0, "n")), 200)
  # a source whose model reads NONE of the declared covariates has nothing to stratify on, and saying so beats
  # silently generating a null contrast
  mN <- function() {
    ini({lcl <- log(3); lv <- log(20); eta.cl ~ 0.09; a <- 0.1})
    model({cl <- exp(lcl + eta.cl); v <- exp(lv)
           d/dt(ce) <- -cl / v * ce; cp <- ce / v; cp ~ add(a)})
  }
  expect_error(admixr2:::.admExpandStrata(mk(mN), "s"), "no contrast to stratify")
  # and TRUE without any model to read it from is refused, not guessed
  expect_error(admixr2:::.admExpandStrata(
    list(s = list(n = 10L, cov_dist = cd, stratify = TRUE)), "s"),
    "no `model` was supplied")
})

test_that("every covariate the ANALYSIS model reads must be described", {
  # A model covariate a study never mentions is held at whatever rxSolve defaults it to -- the ecological
  # plug-in wearing a fit's clothes.
  ui <- .cov_ui(cov = c("WT", "CRCL"), expr = list(
    quote(cl <- exp(tcl + tcov * WT + tcr * CRCL + eta.cl))))
  chk <- function(st) admixr2:::.admCheckCovariates(ui, .cov_pinfo(), st)
  expect_error(chk(list(a = list(cov_dist = list(WT = list(mu = 0, sd = 1))))),
               "does not describe covariate")
  # a fixed value is a legitimate description when it does not vary
  expect_silent(chk(list(a = list(cov = list(CRCL = 90),
                                  cov_dist = list(WT = list(mu = 0, sd = 1))))))
  # ... as is a distribution for both
  expect_silent(chk(list(a = list(cov_dist = list(WT = list(mu = 0, sd = 1),
                                                  CRCL = list(mu = 0, sd = 1))))))
  # a study that declares NO cov_dist has not opted in and is left alone
  expect_silent(chk(list(a = list(cov = list(WT = 0, CRCL = 0)))))
})

test_that("covDist() accepts what a baseline table reports", {
  cd <- covDist(WT = c(mean = 72, sd = 16), CRCL = c(mean = 90, sd = 25),
                cor = 0.6, dist = "lnorm")
  expect_s3_class(cd, "covDist")
  expect_identical(admixr2:::.admCovSpecNames(cd), c("WT", "CRCL"))
  X <- covDraw(cd, n = 20000L)
  expect_equal(unname(colMeans(X)), c(72, 90), tolerance = 1e-2)
  expect_equal(unname(apply(X, 2L, sd)), c(16, 25), tolerance = 1e-2)
  # a named matrix, a scalar and the older spellings must all agree
  R <- matrix(c(1, .6, .6, 1), 2, 2)
  expect_equal(covDraw(covDist(WT = c(mean = 72, sd = 16),
                               CRCL = c(mean = 90, sd = 25), cor = R,
                               dist = "lnorm"), n = 5000L),
               covDraw(cd, n = 5000L), tolerance = 1e-8)
  # ... and it drops straight into every consumer a plain list served
  expect_length(covStrata(cd, stratify = "WT", n_nodes = 3L), 3L)
  expect_identical(nrow(admixr2:::.admCovGrid(cd, 5L)$X), 25L)
})

test_that("covDist() reads a categorical covariate from its labels", {
  cd <- covDist(SEX = c(female = 0.55, male = 0.45))
  expect_equal(cd$SEX$values, c(0, 1))
  expect_equal(cd$SEX$probs, c(0.55, 0.45))
  expect_identical(cd$SEX$labels, c("female", "male"))
  expect_equal(mean(covDraw(cd, n = 20000L)[, "SEX"]), 0.45, tolerance = 1e-2)
})

test_that("covDist() transcribes a data.frame", {
  tbl <- data.frame(covariate = c("WT", "CRCL"), mean = c(72, 90),
                    sd = c(16, 25), dist = "lnorm")
  X <- covDraw(covDist(tbl), n = 20000L)
  expect_identical(colnames(X), c("WT", "CRCL"))
  expect_equal(unname(colMeans(X)), c(72, 90), tolerance = 1e-2)
  # the naming column may be called several things, and `dist` is honoured
  tbl2 <- data.frame(name = "WT", mean = 0, sd = 1, dist = "normal")
  expect_equal(mean(covDraw(covDist(tbl2), n = 20000L)[, "WT"]), 0,
               tolerance = 5e-2)
  expect_error(covDist(data.frame(mean = 1, sd = 1)), "column naming the")
  expect_error(covDist(data.frame(covariate = "WT", mu = 1)), "`mean` and `sd`")
})

test_that("covDist() refuses ambiguity, at construction, naming the covariate", {
  # c(0, 1) is a mean and an SD, or two levels -- and no fit can tell which
  expect_error(covDist(WT = c(72, 16)), "needs NAMED values")
  expect_error(covDist(WT = c(72, 16)), "'WT'")
  # legal under the normal default; impossible once lnorm is asked for
  expect_s3_class(covDist(WT = c(mean = -5, sd = 2)), "covDist")
  expect_error(covDist(WT = c(mean = -5, sd = 2), dist = "lnorm"),
               "lognormal margin cannot")
  expect_error(covDist(A = c(mean = 1, sd = 1), B = c(mean = 1, sd = 1),
                       cor = 1.4), "not positive definite")
  expect_error(covDist(A = c(mean = 1, sd = 1), B = c(mean = 1, sd = 1),
                       cor = diag(c(4, 1))), "unit diagonal")
  expect_error(covDist(A = c(mean = 1, sd = 1), B = c(mean = 1, sd = 1),
                       cor = matrix(c(1, 0.9, 0, 1), 2)), "not symmetric")
  expect_error(covDist(A = c(mean = 1, sd = 1), B = c(mean = 1, sd = 1),
                       cor = matrix(c(1, NA, NA, 1), 2)), "finite")
  expect_error(covDist(list(mean = 1, sd = 1)), "NAMED argument per covariate")
  # `dist` sets the default margin, and a per-covariate one still wins
  expect_equal(mean(covDraw(covDist(A = c(mean = 0, sd = 1), dist = "normal"),
                            n = 20000L)[, "A"]), 0, tolerance = 5e-2)
})

test_that("print.covDist reports what was DECLARED, not what it canonicalised to", {
  # `cor` becomes a joint sampler internally; printing "joint sampler" would hide the number the user typed
  out <- utils::capture.output(print(covDist(WT = c(mean = 72, sd = 16),
                                             CRCL = c(mean = 90, sd = 25),
                                             cor = 0.6)))
  expect_true(any(grepl("cor = 0.6", out, fixed = TRUE)))
  expect_false(any(grepl("joint sampler", out, fixed = TRUE)))
  expect_true(any(grepl("normal", out)))       # the default margin
  out2 <- utils::capture.output(print(covDist(SEX = c(f = 0.6, m = 0.4))))
  expect_true(any(grepl("categorical", out2)))
  expect_true(any(grepl("f=0", out2, fixed = TRUE)))
  out3 <- utils::capture.output(print(covDist(SEX = c(male = 0.55))))
  expect_true(any(grepl("not male=0", out3, fixed = TRUE)))
  expect_true(any(grepl("male=1", out3, fixed = TRUE)))
})

# Shift path -- the pieces that need no compiled model

test_that(".admShiftDelta is the log-scale shift, vectorised over the nodes", {
  spec <- list(eta = "eta.cl", link = "exp",
               rhs = list(quote(exp(tcl + eta.cl) * (WT / 70)^b1)))
  st <- list(tcl = log(4), b1 = 0.75)
  X  <- matrix(c(55, 70, 92, 120), 4L, 1L, dimnames = list(NULL, "WT"))
  D  <- admixr2:::.admShiftDelta(spec, st, X, list(WT = 70))
  # one COLUMN per affected random effect
  expect_equal(dim(D), c(4L, 1L))
  # for a multiplicative allometric term the shift is b1 * log(WT/ref) exactly
  expect_equal(as.numeric(D), 0.75 * (log(X[, 1L]) - log(70)))
  # ... and it is ONE evaluation over the whole node set: a per-row loop that lost the column names would
  # leave every node at the reference and return zeros, which is finite, plausible and wrong
  expect_false(all(D == 0))

  # m = 2: a covariate on two mu-referenced parameters gives two columns, and the shift for each is that
  # parameter's own coefficient
  sp2 <- list(eta = c("eta.cl", "eta.v"), link = c("exp", "exp"),
              rhs = list(quote(exp(tcl + eta.cl) * (WT / 70)^b1),
                         quote(exp(tv + eta.v) * (WT / 70)^c1)))
  D2 <- admixr2:::.admShiftDelta(sp2, list(tcl = log(4), tv = log(30),
                                           b1 = 0.75, c1 = 1.0),
                                 X, list(WT = 70))
  expect_equal(dim(D2), c(4L, 2L))
  expect_equal(D2[, 1L], 0.75 * (log(X[, 1L]) - log(70)))
  expect_equal(D2[, 2L], 1.00 * (log(X[, 1L]) - log(70)))
})

test_that("the affine certificate decides the Gaussian branch in 2-D", {
  # Delta = c + B z is exactly (jointly) normal because admixr2 builds every continuous covariate as
  # F^-1(Phi(z)) from a jointly normal z. A moment test on Delta cannot certify this above one dimension --
  # Cramer-Wold needs ALL projections, so finitely many can only fail to find a counterexample.
  gh <- admixr2:::.adghNodes1; nc <- 21L; g <- gh(nc)
  ix <- as.matrix(expand.grid(seq_len(nc), seq_len(nc)))
  z  <- cbind(g$x[ix[, 1]], g$x[ix[, 2]])
  W  <- g$w[ix[, 1]] * g$w[ix[, 2]]; W <- W / sum(W)
  aff <- list(correlated  = z %*% t(matrix(c(0.40, 0.15, 0, 0.30), 2, 2)),
              independent = cbind(0.4 * z[, 1], 0.3 * z[, 2]),
              lognorm_log = cbind(0.4 * z[, 1], 0.5 * log(exp(0.2 * z[, 2]))))
  non <- list(squared   = cbind(0.4 * z[, 1], 0.3 * z[, 2]^2),
              lognormal = cbind(0.4 * z[, 1], 0.3 * exp(0.2 * z[, 2])))
  for (nm in names(aff)) {
    expect_lt(admixr2:::.admShiftAffineResid(aff[[nm]], W, z), 1e-8)
    expect_true(admixr2:::.admShiftGaussOK(aff[[nm]], W, z, 2L))
  }
  for (nm in names(non)) {
    expect_gt(admixr2:::.admShiftAffineResid(non[[nm]], W, z), 1e-3)
    expect_false(admixr2:::.admShiftGaussOK(non[[nm]], W, z, 2L))
  }
  # WITHOUT the latent scores a vector shift must not take the branch, even though every margin here is
  # normal: that is the case a moment test would wave through without being able to justify it.
  expect_false(admixr2:::.admShiftGaussOK(aff$correlated, W, NULL, 2L))
})

test_that(".admCovGrid returns the latent scores behind X", {
  d <- list(WT = list(mu = 70, sd = 8), AGE = list(mu = 50, sd = 10))
  g <- admixr2:::.admCovGrid(d, 5L)
  expect_equal(nrow(g$z), nrow(g$X))
  expect_equal(ncol(g$z), 2L)
  # X is the image of z under F^-1(Phi(.)), so for a normal margin it is affine
  expect_lt(admixr2:::.admShiftAffineResid(g$X, g$W, g$z), 1e-10)
  # a discrete margin has no score of its own but does NOT void the others: the scores come back for the
  # continuous columns only, with X's stride, so a Delta that never reaches the discrete covariate still
  # certifies
  d2 <- list(WT = list(mu = 70, sd = 8), SEX = list(values = c(0, 1)))
  g2 <- admixr2:::.admCovGrid(d2, 7L)
  expect_equal(nrow(g2$z), nrow(g2$X))     # stride, not just presence
  expect_equal(ncol(g2$z), 1L)
  D_no  <- matrix(0.8 * g2$X[, "WT"] / 70, ncol = 1)
  D_yes <- matrix(0.8 * g2$X[, "WT"] / 70 + 0.3 * g2$X[, "SEX"], ncol = 1)
  expect_lt(admixr2:::.admShiftAffineResid(D_no,  g2$W, g2$z), 1e-8)
  expect_gt(admixr2:::.admShiftAffineResid(D_yes, g2$W, g2$z), 1e-3)
  # a correlated grid holding a discrete margin goes through the pool instead, which has no latent structure
  # at all
  d3 <- list(WT = list(mu = 70, sd = 8), CRCL = list(mu = 90, sd = 22),
             SEX = list(values = c(0, 1)),
             cor = matrix(c(1, .5, 0, .5, 1, 0, 0, 0, 1), 3, 3))
  expect_null(admixr2:::.admCovGrid(d3, 5L)$z)
})

test_that("lognormal and correlated covariates reach the Gaussian branch", {
  # The two forms this exists for, on real grids rather than synthetic Delta. Correlated margins certify
  # because the Gaussian copula maps X_k = F_k^-1(Phi((L z)_k)) and (L z)_k is linear in z.
  cg <- admixr2:::.admCovGrid; ar <- admixr2:::.admShiftAffineResid
  yes <- list(
    list(list(WT = list(meanlog = log(70), sdlog = 0.17)),
         function(X) 0.75 * log(X[, "WT"] / 70)),
    list(list(WT = list(mu = 70, sd = 8), CRCL = list(mu = 90, sd = 22),
              cor = 0.6),
         function(X) cbind(0.8 * X[, "WT"] / 70, 0.3 * X[, "CRCL"] / 90)),
    list(list(WT = list(meanlog = log(70), sdlog = 0.17),
              CRCL = list(meanlog = log(90), sdlog = 0.25), cor = 0.6),
         function(X) cbind(0.75 * log(X[, "WT"] / 70),
                           0.4 * log(X[, "CRCL"] / 90))),
    list(list(WT = list(meanlog = log(70), sdlog = 0.17),
              AGE = list(mu = 50, sd = 12), cor = 0.4),
         function(X) cbind(0.75 * log(X[, "WT"] / 70), 0.02 * X[, "AGE"])))
  for (cs in yes) {
    g <- cg(cs[[1]], 7L)
    expect_lt(ar(as.matrix(cs[[2]](g$X)), g$W, g$z), 1e-8)
  }
  # a lognormal covariate entering RAW is not affine in the score, and is refused
  g <- cg(list(WT = list(meanlog = log(70), sdlog = 0.17)), 7L)
  expect_gt(ar(matrix(0.8 * g$X[, "WT"], ncol = 1), g$W, g$z), 1e-3)
})

test_that("an ALL-discrete covariate set still takes the product grid", {
  skip_on_cran()
  skip_if_not_installed("rxode2")
  skip_if_not_installed("rxode2")
  # The shift's saving is eliminating the CONTINUOUS covariate dimension. The discrete levels cost K cells on
  # either route, and n_u is inflated above n_nodes to resolve the widened u, so a stratified shift comes out
  # with MORE rows than the grid (26 vs 18 for one binary covariate) -- and as an approximation where the grid
  # enumerates exactly.
  .f <- function() {
    ini({ tcl <- log(1); tv <- log(10); tsex <- 0.2
          eta.cl ~ 0.09; add.err <- 0.3 })
    model({ cl <- exp(tcl + tsex * SEX + eta.cl); v <- exp(tv)
            cp <- linCmt(); cp ~ add(add.err) })
  }
  ui <- suppressMessages(rxode2::rxode2(.f))
  st <- list(s = list(E = 1:3, V = diag(3), n = 10L, times = 1:3,
                      ev = rxode2::et(amt = 100),
                      cov_dist = list(SEX = list(values = c(0, 1)))))
  ctl <- adghControl(studies = st, grad = "analytical", n_nodes = 5L,
                     print = 0L, covMethod = "none", cov_integration = "on")
  pin <- admixr2:::.admDriverPinfo(ui, ctl)
  u   <- admixr2:::.admDriverUnits(st, ui, admixr2:::.admOutputVar(ui))
  out <- suppressMessages(admixr2:::.admCheckCovariates(ui, pin, u$studies))
  # An all-discrete covariate set has no continuous dimension to reduce, so the product grid -- which
  # enumerates the levels exactly -- is what runs.
  expect_identical(out[[1L]]$.adm_cov_path, "rows")
})

test_that("the identifiability warning canonicalises the user's shorthand", {
  # It is the one entry point that reads the RAW study list -- every driver calls it before normalising -- so
  # it is the one that meets `mean`/`sd` un-expanded. .admCovMeanOf has no `mean` branch, so the mean came
  # back NA, between-study variation was invisible, and the warning fired on exactly the contrast that
  # identifies the coefficient.
  ui <- .cov_ui()                       # cl <- exp(tcl + tcov*WT + eta.cl)
  pin <- .cov_pinfo()
  diff_mu <- list(a = list(cov_dist = list(WT = list(mean = -0.4, sd = 0.6))),
                  b = list(cov_dist = list(WT = list(mean =  0.4, sd = 0.6))))
  same_mu <- list(a = list(cov_dist = list(WT = list(mean = 0, sd = 0.6))),
                  b = list(cov_dist = list(WT = list(mean = 0, sd = 0.6))))
  expect_silent(admixr2:::.admWarnCovIdentifiability(ui, pin, diff_mu))
  expect_warning(admixr2:::.admWarnCovIdentifiability(ui, pin, same_mu),
                 "not identifiable")
})

test_that("bands are bounded by the range the source REPORTED", {
  # Strata cut over the analyst's full `cov_dist` evaluate a published model in covariate bands that study
  # never enrolled, and credit its coefficient as evidence there. The overstatement is exactly
  # var(declared)/var(enrolled) -- 3.43x for a source spanning +/-1 SD, 12.38x at +/-0.5 SD. Not bias: false
  # confidence, and it only bites once a second source disagrees.
  cd <- covDist(WT = c(mean = 78, sd = 16), dist = "lnorm")
  full <- suppressWarnings(covStrata(cd, "WT", n_nodes = 5L, n = 300))
  rng  <- covStrata(cd, "WT", n_nodes = 5L, n = 300,
                    cov_range = list(WT = c(60, 100)))
  spanof <- function(st) range(vapply(st, function(s) s$cov$WT, numeric(1)))
  expect_gt(spanof(full)[2L] - spanof(full)[1L],
            spanof(rng)[2L]  - spanof(rng)[1L])
  expect_gte(spanof(rng)[1L], 60)
  expect_lte(spanof(rng)[2L], 100)
  # the weights still partition the study: sum to 1, and n is preserved
  expect_equal(sum(vapply(rng, `[[`, 0, "weight")), 1, tolerance = 1e-10)
  expect_equal(sum(vapply(rng, `[[`, 0, "n")), 300, tolerance = 1e-8)
  # and the credited spread -- which IS the information claimed -- falls
  cvar <- function(st) {
    x <- vapply(st, function(s) s$cov$WT, numeric(1))
    w <- vapply(st, `[[`, 0, "weight")
    sum(w * (x - sum(w * x))^2)
  }
  expect_gt(cvar(full) / cvar(rng), 1.5)
  # absent, it WARNS -- the silent case is the leaky one
  expect_warning(covStrata(cd, "WT", n_nodes = 5L, n = 300), "cov_range")
  # a range that excludes the distribution is an error, not an empty result
  expect_error(covStrata(cd, "WT", cov_range = list(WT = c(1e5, 2e5))),
               "covers essentially none")
  expect_error(covStrata(cd, "WT", cov_range = list(ZZ = c(1, 2))),
               "does not declare")
})

test_that("a source that ASSERTED a covariate's coefficient is not banded", {
  skip_on_cran()
  skip_if_not_installed("rxode2")
  # `allCovs` reports what a model READS, not what it ESTIMATED. A model carrying (WT/70)^0.75 -- the
  # allometric convention -- reads WT while asserting its coefficient, so it holds no evidence about WT at all
  # and the claimed-to-earned information ratio is unbounded at every stratum count.
  skip_if_not_installed("rxode2")
  cd <- covDist(WT = c(mean = 78, sd = 16), dist = "lnorm")
  st <- list(a = list(times = c(1, 4, 12), ev = rxode2::et(amt = 100),
                      n = 100L, cov_dist = cd, stratify = TRUE,
                      cov_range = list(WT = c(60, 100))))
  asserted <- function() {
    ini({ tcl <- log(5); tv <- log(50); eta.cl ~ .09; add.err <- .3 })
    model({ cl <- exp(tcl + eta.cl) * (WT/70)^0.75
            v <- exp(tv); cp <- linCmt(); cp ~ add(add.err) })
  }
  expect_error(admixr2:::.admExpandStrata(st, names(st), asserted),
               "ASSERTS the coefficient")
  estimated <- function() {
    ini({ tcl <- log(5); tv <- log(50); bwt <- 0.75; eta.cl ~ .09
          add.err <- .3 })
    model({ cl <- exp(tcl + eta.cl) * (WT/70)^bwt
            v <- exp(tv); cp <- linCmt(); cp ~ add(add.err) })
  }
  ex <- admixr2:::.admExpandStrata(st, names(st), estimated)
  expect_gt(length(ex$studies), 1L)
})

test_that(".admCovCoefThetas separates an estimated coefficient from an asserted one", {
  skip_on_cran()
  skip_if_not_installed("rxode2")
  # Answered numerically, and on the LOG scale, because these models are multiplicative: raising a SCALE
  # parameter changes the covariate's absolute effect, so differencing the raw prediction flags `tcl` as WT's
  # coefficient, which it is not. What a scale parameter leaves alone is the PROPORTIONAL effect -- d(log
  # f)/d(log WT) is 0.75 whatever tcl is.
  skip_if_not_installed("rxode2")
  cd <- covDist(WT = c(mean = 78, sd = 16), CRCL = c(mean = 90, sd = 20),
                dist = "lnorm")
  co <- function(m, cov) admixr2:::.admCovCoefThetas(
    suppressMessages(rxode2::rxode2(m)), cov, cd)
  lit <- function() {
    ini({ tcl <- log(5); tv <- log(50); eta.cl ~ .09; add.err <- .3 })
    model({ cl <- exp(tcl + eta.cl) * (WT/70)^0.75
            v <- exp(tv); cp <- linCmt(); cp ~ add(add.err) })
  }
  est <- function() {
    ini({ tcl <- log(5); tv <- log(50); bwt <- 0.75; eta.cl ~ .09; add.err <- .3 })
    model({ cl <- exp(tcl + eta.cl) * (WT/70)^bwt
            v <- exp(tv); cp <- linCmt(); cp ~ add(add.err) })
  }
  fx <- function() {
    ini({ tcl <- log(5); tv <- log(50); bwt <- fix(0.75); eta.cl ~ .09
          add.err <- .3 })
    model({ cl <- exp(tcl + eta.cl) * (WT/70)^bwt
            v <- exp(tv); cp <- linCmt(); cp ~ add(add.err) })
  }
  part <- function() {
    ini({ tcl <- log(5); tv <- log(50); bcr <- 0.6; eta.cl ~ .09; add.err <- .3 })
    model({ cl <- exp(tcl + eta.cl) * (WT/70)^0.75 * (CRCL/90)^bcr
            v <- exp(tv); cp <- linCmt(); cp ~ add(add.err) })
  }
  expect_length(co(lit, "WT"), 0L)          # literal exponent: asserted
  expect_identical(co(est, "WT"), "bwt")    # estimated
  expect_length(co(fx,  "WT"), 0L)          # fix()ed theta: asserted
  expect_length(co(part, "WT"), 0L)         # partial set: WT asserted ...
  expect_identical(co(part, "CRCL"), "bcr") # ... CRCL estimated
})

test_that("strata_nodes is a convergence parameter, and fits carry it", {
  skip_on_cran()
  skip_if_not_installed("rxode2")
  # The J -> infinity limit is the well-defined object; finite J approximates it, and under misspecification
  # the answer can jump basins rather than drift (1.57 at J=4 -> 1.6e-05 at J=5). The objective is J-dependent
  # -- 441 units across J = 1 to 500 with the model correct -- so two fits at different resolutions are on
  # different scales.
  skip_if_not_installed("rxode2")
  expect_gt(admixr2:::.ADM_STRATA_NODES, 5L)
  cd <- covDist(WT = c(mean = 78, sd = 16), dist = "lnorm")
  st <- list(a = list(times = c(1, 4, 12), ev = rxode2::et(amt = 100),
                      n = 100L, cov_dist = cd, stratify = TRUE,
                      cov_range = list(WT = c(60, 100))))
  m <- function() {
    ini({ tcl <- log(5); tv <- log(50); bwt <- 0.75; eta.cl ~ .09; add.err <- .3 })
    model({ cl <- exp(tcl + eta.cl) * (WT/70)^bwt
            v <- exp(tv); cp <- linCmt(); cp ~ add(add.err) })
  }
  ex <- admixr2:::.admExpandStrata(st, names(st), m)
  expect_equal(length(ex$studies), admixr2:::.ADM_STRATA_NODES)
  # every generated study is stamped, so a fit can refuse a cross-J comparison
  expect_true(all(vapply(ex$studies, function(s)
    identical(s[[".adm_strata_nodes"]], admixr2:::.ADM_STRATA_NODES),
    logical(1))))
  st2 <- st; st2$a$strata_nodes <- 4L
  ex2 <- admixr2:::.admExpandStrata(st2, names(st2), m)
  expect_equal(ex2$studies[[1L]][[".adm_strata_nodes"]], 4L)
})

test_that("banding uses the quadrature rule that CONVERGES", {
  # The objective has to converge in the stratum count or it is not a likelihood -- AIC, BIC and any ratio are
  # otherwise on an arbitrary scale.

  # The pooled branch cuts equiprobable bins and evaluates at a representative point, which is a midpoint
  # rule: measured on an independent pair, it moved the objective 749 units between J = 3 and J = 25 and was
  # still moving. The Gauss-Hermite branch moved 0.2 units in total. Independence is a KNOWN latent structure,
  # so it takes the fast branch too -- it used to fall to the pooled one purely because no `cor` had been
  # supplied.
  skip_if_not_installed("randtoolbox")
  cd <- covDist(WT = c(mean = 78, sd = 16), CRCL = c(mean = 90, sd = 20),
                dist = "lnorm")
  expect_null(cd[["latentR"]])            # no dependence declared ...
  st <- covStrata(cd, "WT", n_nodes = 5L, n = 400,
                  cov_range = list(WT = c(55, 108)))
  # ... and the GH branch still ran: it hands back a per-stratum n_pool_cell, which the pooled branch does not
  # produce the same way
  expect_length(st, 5L)
  expect_true(all(vapply(st, function(s) !is.null(s$cov$WT), logical(1))))
  # an explicit identity `cor` must give the SAME strata -- same distribution, so the branch must not change
  # the answer
  I2 <- diag(2); dimnames(I2) <- list(c("WT", "CRCL"), c("WT", "CRCL"))
  cdI <- covDist(WT = c(mean = 78, sd = 16), CRCL = c(mean = 90, sd = 20),
                 dist = "lnorm", cor = I2)
  stI <- covStrata(cdI, "WT", n_nodes = 5L, n = 400,
                   cov_range = list(WT = c(55, 108)))
  expect_equal(vapply(st,  function(s) s$cov$WT, numeric(1)),
               vapply(stI, function(s) s$cov$WT, numeric(1)), tolerance = 1e-10)
  expect_equal(vapply(st,  `[[`, 0, "weight"),
               vapply(stI, `[[`, 0, "weight"), tolerance = 1e-10)
  # an OPAQUE joint cannot be conditioned, so it keeps the pooled route
  cdj <- covDist(WT = c(mean = 78, sd = 16), CRCL = c(mean = 90, sd = 20),
                 dist = "lnorm",
                 joint = function(u) cbind(WT = stats::qlnorm(u[, 1], log(78), .2),
                                           CRCL = stats::qlnorm(u[, 2], log(90), .2)))
  expect_length(suppressWarnings(covStrata(cdj, "WT", n_nodes = 4L, n = 100)), 4L)
})

test_that("a declared DISCRETE covariate no longer forces the pooled rule", {
  # Discreteness was never the obstacle -- the SAMPLER was. A margin latently independent of the rest has
  # chol(R)[, j] = e_j, so its level is monotone in its own uniform and enumerates exactly. Before this, one
  # declared sex or genotype put the whole study on equiprobable bins: 515 units of J-dependence across J = 5
  # to 50, still rising, against 0.009 after.
  cd <- covDist(WT = c(mean = 78, sd = 16), dist = "lnorm",
                SEX = list(values = c(0, 1), probs = c(0.45, 0.55)))
  st <- covStrata(cd, "WT", n_nodes = 5L, n = 400,
                  cov_range = list(WT = c(50, 115)))
  expect_length(st, 5L)
  expect_equal(sum(vapply(st, `[[`, 0, "weight")), 1)
  for (s in st) {
    # the quadrature route hands out the DECLARED spec, not a sample standing in for it: no manufactured
    # `joint`, and nothing pooled to report
    expect_null(s$cov_dist[["joint"]])
    expect_identical(s$cov_dist$SEX$probs, c(0.45, 0.55))
    expect_true(is.na(s$n_pool_cell))
    # WT is held at a POINT
    X <- covDraw(s$cov_dist, n = 4000L)
    expect_lt(stats::sd(X[, "WT"]), 1e-8)
    expect_equal(mean(X[, "SEX"] == 1), 0.55, tolerance = 5e-3)
  }
})

test_that("stratifying on a discrete covariate crosses levels with the GH grid", {
  cd <- covDist(WT = c(mean = 78, sd = 16), dist = "lnorm",
                SEX = list(values = c(0, 1), probs = c(0.45, 0.55)))
  st <- covStrata(cd, c("WT", "SEX"), n_nodes = 5L, n = 400,
                  cov_range = list(WT = c(50, 115)))
  expect_length(st, 10L)                       # 5 nodes x 2 levels
  expect_equal(sum(vapply(st, `[[`, 0, "weight")), 1)
  # the level probability factorises out exactly -- the discrete block is latently independent, so the cell
  # weight is w_k * p_level
  w1 <- sum(vapply(st, function(s) if (s$cov$SEX == 1) s$weight else 0, 0))
  expect_equal(w1, 0.55, tolerance = 1e-12)
  # every stratum is a point in BOTH covariates
  for (s in st) {
    X <- covDraw(s$cov_dist, n = 2000L)
    expect_lt(stats::sd(X[, "WT"]), 1e-8)
    expect_equal(unique(X[, "SEX"]), s$cov$SEX)
  }
  # and the quadrature nodes are the SAME ones the unstratified-sex design used
  st0 <- covStrata(cd, "WT", n_nodes = 5L, n = 400,
                   cov_range = list(WT = c(50, 115)))
  expect_equal(sort(unique(vapply(st, function(s) s$cov$WT, 0))),
               sort(vapply(st0, function(s) s$cov$WT, 0)), tolerance = 1e-10)
})

test_that(".admCovGrid enumerates a separable discrete margin under `cor`", {
  # A discrete margin correlated with a continuous one is a TRUNCATION of the latent, not a point, and still
  # falls to the pool. One that is latently independent does not, even though the OTHER two covariates are
  # dependent.
  R <- diag(3); R[1, 2] <- R[2, 1] <- 0.7
  dimnames(R) <- list(c("WT", "CRCL", "SEX"), c("WT", "CRCL", "SEX"))
  cd <- covDist(WT = c(mean = 78, sd = 16), CRCL = c(mean = 90, sd = 20),
                SEX = list(values = c(0, 1), probs = c(0.45, 0.55)),
                dist = "lnorm", cor = R)
  expect_identical(cd[["discExact"]], "SEX")
  g <- admixr2:::.admCovGrid(cd, 7L)
  expect_equal(sum(g$W), 1)
  expect_equal(nrow(g$X), 7L * 7L * 2L)        # GH^2 crossed with 2 levels
  # EXACT level probabilities -- the pooled route reported 0.477 for a declared 0.55, and that error does not
  # shrink with cov_nodes
  expect_equal(sum(g$W[g$X[, "SEX"] == 1]), 0.55, tolerance = 1e-12)
  # the dependence between the correlated pair survives the crossing
  mw <- sum(g$W * g$X[, "WT"]); mc <- sum(g$W * g$X[, "CRCL"])
  cv <- sum(g$W * (g$X[, "WT"] - mw) * (g$X[, "CRCL"] - mc)) /
    sqrt(sum(g$W * (g$X[, "WT"] - mw)^2) * sum(g$W * (g$X[, "CRCL"] - mc)^2))
  expect_gt(cv, 0.5)
  # ... and a discrete margin that IS correlated keeps the pool
  R2 <- R; R2[2, 3] <- R2[3, 2] <- 0.4
  cd2 <- covDist(WT = c(mean = 78, sd = 16), CRCL = c(mean = 90, sd = 20),
                 SEX = list(values = c(0, 1), probs = c(0.45, 0.55)),
                 dist = "lnorm", cor = R2)
  expect_length(cd2[["discExact"]], 0L)
  g2 <- admixr2:::.admCovGrid(cd2, 7L)
  expect_true(all(abs(g2$W - g2$W[1L]) < 1e-12))     # equal-weight pool
})

test_that("a marginalised discrete covariate with no contrast is FLAGGED", {
  skip_on_cran()
  skip_if_not_installed("rxode2")
  # It is not identified -- its effect enters only through the mixture it induces, which is what an eta on the
  # same parameter does, and with the same level distribution in every study there is no between-study
  # contrast either. Measured: the profile moves 0.019 units across the coefficient's whole range and the
  # optimizer settled at -0.059 against a truth of +0.150, at every resolution. A deterministic optimizer on a
  # flat ridge stops in the same place every run, so it reads as converged. Silence is the hazard.
  m <- function() {
    ini({ tcl <- log(4); tv <- log(45); bsex <- 0.05
          eta.cl ~ 0.1; add.err <- 0.1 })
    model({ cl <- exp(tcl + eta.cl) * exp(bsex * SEX)
            v  <- exp(tv); cp <- linCmt(); cp ~ add(add.err) })
  }
  ui <- suppressMessages(rxode2::rxode2(m))
  sp <- list(values = c(0, 1), probs = c(0.45, 0.55))
  st <- list(a = list(cov_dist = list(SEX = sp)),
             b = list(cov_dist = list(SEX = sp)))
  expect_warning(admixr2:::.admCovDiscContrast(ui, st, names(st)),
                 "DISCRETE covariate")
  # a DIFFERENT level distribution in the second study is a contrast
  st2 <- st; st2$b$cov_dist$SEX$probs <- c(0.8, 0.2)
  expect_silent(admixr2:::.admCovDiscContrast(ui, st2, names(st2)))
  # so is a study that pins it at a value
  st3 <- st; st3$b$cov_dist <- NULL; st3$b$cov <- list(SEX = 1)
  expect_silent(admixr2:::.admCovDiscContrast(ui, st3, "a"))
  # NO random effect on the parameter the covariate modulates -- the mixture is then the only thing putting
  # spread on it, so it IS identified from V and a warning would be a false positive
  m0 <- function() {
    ini({ tcl <- log(4); tv <- log(45); bsex <- 0.05
          eta.v ~ 0.1; add.err <- 0.1 })
    model({ cl <- exp(tcl) * exp(bsex * SEX)
            v  <- exp(tv + eta.v); cp <- linCmt(); cp ~ add(add.err) })
  }
  expect_false(admixr2:::.admCovMeetsEta(
    suppressMessages(rxode2::rxode2(m0)), "SEX"))
  expect_silent(admixr2:::.admCovDiscContrast(
    suppressMessages(rxode2::rxode2(m0)), st, names(st)))
  # ... and it is followed TRANSITIVELY through an intermediate assignment
  m1 <- function() {
    ini({ tcl <- log(4); tv <- log(45); bsex <- 0.05
          eta.cl ~ 0.1; add.err <- 0.1 })
    model({ cl0 <- exp(tcl + eta.cl); cl <- cl0 * exp(bsex * SEX)
            v  <- exp(tv); cp <- linCmt(); cp ~ add(add.err) })
  }
  expect_true(admixr2:::.admCovMeetsEta(
    suppressMessages(rxode2::rxode2(m1)), "SEX"))
  # and an ASSERTED coefficient carries no claim to identify, so it is silent
  m2 <- function() {
    ini({ tcl <- log(4); tv <- log(45); bsex <- fix(0.05)
          eta.cl ~ 0.1; add.err <- 0.1 })
    model({ cl <- exp(tcl + eta.cl) * exp(bsex * SEX)
            v  <- exp(tv); cp <- linCmt(); cp ~ add(add.err) })
  }
  expect_silent(admixr2:::.admCovDiscContrast(
    suppressMessages(rxode2::rxode2(m2)), st, names(st)))
})

test_that("a correlated conditional still enumerates an independent discrete margin", {
  # The stratified covariate is correlated with a marginalised one, so that PAIR needs a conditional pool --
  # there is no closed-form spec for a latent shifted off the origin. A third covariate, discrete and latently
  # independent, must not be swept into it: it keeps its exact spec and rides the grid at its levels. This is
  # the path where `discExact` and the degenerate stratified specs coexist, so it is also what pins the `disc
  # & discExact` mask.
  R <- diag(3); R[1L, 2L] <- R[2L, 1L] <- 0.7
  dimnames(R) <- list(c("WT", "CRCL", "SEX"), c("WT", "CRCL", "SEX"))
  cd <- covDist(WT = c(mean = 78, sd = 16), CRCL = c(mean = 90, sd = 20),
                SEX = list(values = c(0, 1), probs = c(0.45, 0.55)),
                dist = "lnorm", cor = R)
  st <- covStrata(cd, "WT", n_nodes = 4L, n = 400,
                  cov_range = list(WT = c(50, 115)))
  expect_length(st, 4L)
  for (s in st) {
    expect_true(is.function(s$cov_dist[["joint"]]))
    expect_true("SEX" %in% s$cov_dist[["discExact"]])
    X <- covDraw(s$cov_dist, n = 8000L)
    expect_lt(stats::sd(X[, "WT"]), 1e-8)            # stratified: a point
    expect_equal(mean(X[, "SEX"] == 1), 0.55, tolerance = 0.02)
    g <- admixr2:::.admCovGrid(s$cov_dist, 5L)
    expect_equal(sum(g$W), 1)
    expect_gt(nrow(g$X), 1L)
    # exact level weights off the grid, not an equal-weight pool
    expect_equal(sum(g$W[g$X[, "SEX"] == 1]), 0.55, tolerance = 1e-12)
    expect_false(all(abs(g$W - g$W[1L]) < 1e-12))
  }
  # ... and the conditional mean of the correlated covariate still MOVES
  em <- vapply(st, function(s) mean(covDraw(s$cov_dist, n = 20000L)[, "CRCL"]), 0)
  expect_gt(diff(range(em)), 10)
})

test_that("covStrata's cov_range truncation survives the second canon", {
  # The first canon builds `joint` as a closure over the UNTRUNCATED margins and consumes `cor` into
  # `latentR`; the second early-returns on an existing `joint`. So the truncation was applied to the specs and
  # then read straight past, and the pool sampled the full declared support -- bands cut over a range the
  # source never enrolled, which is the var(declared)/var(enrolled) overstatement the cov_range gate exists to
  # prevent (3.4x at +/-1 SD).
  cd <- covDist(WT = c(mean = 75, sd = 25), CRCL = c(mean = 90, sd = 30),
                cor = c(WT.CRCL = 0.5))
  st <- suppressWarnings(admixr2:::.admCovStrata(cd, "WT", n_nodes = 3L,
                                                 n_pool = 4096L,
                                                 cov_range = list(WT = c(60, 100))))
  expect_false(is.null(st))
  # every banded value sits inside the enrolled range
  wt <- unlist(lapply(st, function(b) b[["cov"]][["WT"]]), use.names = FALSE)
  expect_true(all(wt >= 60 - 1e-8 & wt <= 100 + 1e-8))
  # ... and the correlation was NOT lost when the derived fields were rebuilt
  cd2 <- st[[1L]][["cov_dist"]]
  expect_true(is.function(cd2[["joint"]]) || !is.null(cd2[["latentR"]]))
})

test_that("a continuous summary is not read as a set of levels", {
  # c(median = 92, iqr = c(62, 118)) and c(mean = 72, cv = 22) are the forms admPopulation() documents; both
  # fell through to the categorical branch and became a 3- and a 2-level covariate, so the quadrature
  # integrated over a covariate that does not exist, with no error.
  s1 <- admixr2:::.admCovSpecFromVec(c(median = 92, iqr = c(62, 118)), "CRCL")
  expect_null(s1$values)
  expect_equal(s1$mu, 92)
  expect_equal(s1$sd, (118 - 62) / 1.34898, tolerance = 1e-8)
  s2 <- admixr2:::.admCovSpecFromVec(c(mean = 72, cv = 22), "WT")
  expect_null(s2$values)
  expect_equal(s2$mu, 72); expect_equal(s2$sd, 72 * 0.22)
  # a genuine proportion table is still categorical
  s3 <- admixr2:::.admCovSpecFromVec(c(female = 0.45, male = 0.55), "SEX")
  expect_identical(s3$values, c(0, 1))
})

test_that("a discrete spec with mismatched probs/values is refused", {
  # admc's .admCovQuantile renormalises over 2 entries so level 2 is UNREACHABLE; adgh's .admCovNodesFor
  # returns 3 nodes against 2 weights, which .admCovGrid recycles to a UNIFORM 3-level covariate. Two
  # estimators, two different distributions, both finite and plausible.
  ui <- list(allCovs = "GRP")
  pin <- list(cov_integration = "on")
  st  <- list(s = list(times = 1, cov_dist = list(
    GRP = list(values = c(0, 1, 2), probs = c(0.5, 0.5)))))
  expect_error(admixr2:::.admCheckCovariates(ui, pin, st),
               "3 `values` but 2 `probs`|one probability per level")
})

# Sparse-grid covariate integration (cov_integration = "sparse")

test_that("level 2 at one covariate IS the 3-node product grid, exactly", {
  # A live property of the sparse rule, not a historical note: at one covariate Smolyak level 2 is
  # bit-identical to .admCovGrid at n_nodes = 3, same points and same weights. (It is also what showed the
  # retired "taylor" design to be a re-labelling of Gauss-Hermite rather than a method of its own, which is
  # why that path could be removed rather than kept beside this one.)
  cd <- covDist(WT = c(mean = 75, sd = 16))
  g2 <- admixr2:::.admCovSparseGrid(cd, 2L)
  gh <- admixr2:::.admCovGrid(cd, 3L)
  expect_identical(nrow(g2$X), 3L)
  expect_equal(sort(as.numeric(g2$X)), sort(as.numeric(gh$X)))
  expect_equal(sort(g2$W), sort(as.numeric(gh$W)))
})

test_that("the sparse grid weights sum to one at every level and dimension", {
  # sum(W) == 1 is the one property the combination technique must never lose; sum|W| is allowed to grow, and
  # IS the price the level buys accuracy with.
  for (p in 1:3) {
    sp <- stats::setNames(
      lapply(seq_len(p), function(k) list(mu = 70 + k, sd = 10 + k)),
      paste0("C", seq_len(p)))
    cd <- do.call(covDist, sp)
    for (lv in 2:4) {
      g <- admixr2:::.admCovSparseGrid(cd, lv)
      expect_equal(sum(g$W), 1, tolerance = 1e-12,
                   info = sprintf("p=%d level=%d", p, lv))
      expect_identical(colnames(g$X), names(sp))
      expect_true(all(is.finite(g$X)))
      # signed weights are expected from level 3 up, and only there
      if (lv == 2L && p < 4L)
        expect_true(all(g$W > 0), info = sprintf("p=%d", p))
    }
    # ... and |W| grows with the level once there is more than one dimension for the combination technique to
    # cancel across (at p = 1 a Smolyak rule IS the 1-D rule, all weights positive, so it stays at exactly 1).
    a <- sum(abs(admixr2:::.admCovSparseGrid(cd, 2L)$W))
    b <- sum(abs(admixr2:::.admCovSparseGrid(cd, 3L)$W))
    if (p == 1L) expect_equal(c(a, b), c(1, 1), tolerance = 1e-12)
    else         expect_gt(b, a)
  }
})

test_that("level 3 is much more accurate than level 2, and correlation does not cost it", {
  # The measurement the default rests on. Reference is a 21-node product grid; the integrand is allometric
  # plus a saturable term so the mixed derivatives a level-2 rule cannot see are genuinely present.
  gq <- function(n) admixr2:::.adghNodes1(n)
  ff <- function(A) exp(0.75 * log(A[, 1L]) + 0.4 * log(A[, 2L])) +
                    0.6 * A[, 1L]^2 / (1 + 0.5 * A[, 2L])
  err <- function(rho) {
    # LOGNORMAL margins: the integrand takes log(A), and a normal margin's tail nodes go negative, which is
    # NaN rather than an accuracy question.
    cd <- covDist(A = c(mean = 1, sd = 0.3), B = c(mean = 1, sd = 0.3),
                  cor = c(A.B = rho), dist = "lnorm")
    ref <- admixr2:::.admCovGrid(cd, 21L)
    tru <- sum(ref$W * ff(ref$X))
    vapply(2:3, function(lv) {
      g <- admixr2:::.admCovSparseGrid(cd, lv)
      abs(sum(g$W * ff(g$X)) - tru) / abs(tru)
    }, numeric(1))
  }
  e0 <- err(0); e5 <- err(0.5); e8 <- err(0.85)
  # Level 3 beats level 2 everywhere -- and by FAR more once the covariates are correlated, which is the whole
  # reason 3 is the default. Measured ratios: 43x at rho = 0, 4551x at 0.5, 12397x at 0.85.
  expect_lt(e0[2L], e0[1L] / 20)
  expect_lt(e5[2L], e5[1L] / 500)
  expect_lt(e8[2L], e8[1L] / 500)
  # The two levels respond to correlation in OPPOSITE directions: level 2 loses an order of magnitude to it,
  # level 3 does not. That asymmetry is the reason the retired Taylor design could not be sold as "good with
  # correlation".
  expect_gt(e5[1L], 10 * e0[1L])
  expect_lt(e8[2L], 1e-5)
})

test_that("the sparse grid enumerates a discrete covariate exactly", {
  # A sparse rule is a statement about the CONTINUOUS dimensions. A level probability is not an approximation
  # of anything, so it is crossed in whole.
  cd <- covDist(WT = c(mean = 75, sd = 16),
                SEX = list(values = c(0, 1), probs = c(0.45, 0.55)))
  g <- admixr2:::.admCovSparseGrid(cd, 3L)
  expect_setequal(colnames(g$X), c("WT", "SEX"))
  expect_setequal(unique(g$X[, "SEX"]), c(0, 1))
  expect_equal(as.numeric(tapply(g$W, g$X[, "SEX"], sum)), c(0.45, 0.55),
               tolerance = 1e-12)
  expect_equal(sum(g$W), 1, tolerance = 1e-12)
})

test_that("the sparse grid refuses an opaque joint sampler and a bad level", {
  # The rotation needs the latent correlation, and the canoniser early-returns on a user closure BEFORE
  # recording it -- so reading it as independent is exactly the silent approximation this file refuses
  # everywhere else.
  cd <- covDist(A = c(mean = 1, sd = 0.2), B = c(mean = 1, sd = 0.2),
                joint = function(u) {
                  z <- stats::qnorm(pmin(pmax(u, 1e-12), 1 - 1e-12))
                  out <- cbind(A = exp(0.2 * z[, 1L]),
                               B = exp(0.2 * (0.8 * z[, 1L] +
                                              sqrt(1 - 0.64) * z[, 2L])))
                  out
                })
  expect_error(admixr2:::.admCovSparseGrid(cd, 3L), "joint. sampler")
  mixed <- covDist(
    SEX = list(values = c(0, 1), probs = c(.5, .5)),
    WT = c(mean = 0, sd = 1),
    joint = function(u)
      cbind(SEX = as.numeric(u[, 1L] > .5), WT = stats::qnorm(u[, 1L])))
  expect_error(admixr2:::.admCovSparseGrid(mixed, 3L), "joint. sampler")
  cd2 <- covDist(WT = c(mean = 75, sd = 16))
  expect_error(admixr2:::.admCovSparseGrid(cd2, 1L), "between 2 and")
  expect_error(admixr2:::.admCovSparseGrid(cd2, 99L), "between 2 and")
})

test_that("adghControl carries the sparse settings and validates the level", {
  st <- list(s = list(E = 1, V = 1, n = 10, times = 1,
                      ev = rxode2::et(amt = 1)))
  d <- adghControl(studies = st)
  expect_identical(d$cov_integration, "on")   # reductions on by default
  expect_identical(d$cov_sparse_level, 3L)
  s <- adghControl(studies = st, cov_integration = "sparse",
                   cov_sparse_level = 2L)
  expect_identical(s$cov_integration, "sparse")
  expect_identical(s$cov_sparse_level, 2L)
  expect_error(adghControl(studies = st, cov_sparse_level = 1L))
  expect_error(adghControl(studies = st, cov_integration = "taylor"))
})

test_that("the sparse grid is built ONCE and cached on the study", {
  # It is a pure function of cov_dist and the level, both DATA, and .adghGrid runs inside the objective -- so
  # it must not be rebuilt per evaluation.
  skip_on_cran()
  skip_if_not_installed("rxode2")
  f <- function() {
    ini({tcl <- log(5); tv <- log(30); bwt <- 0.75; eta.cl ~ 0.09; a <- 0.1})
    model({cl <- exp(tcl + eta.cl) * (WT / 70)^bwt; v <- exp(tv)
           cp <- linCmt(); cp ~ add(a)})
  }
  ui  <- suppressMessages(rxode2::rxode2(f))
  st  <- list(s = list(E = c(9, 7, 5), V = diag(c(1, 1, 1)), n = 50L,
                       times = c(1, 4, 12), ev = rxode2::et(amt = 100),
                       cov_dist = list(WT = list(mu = 75, sd = 16))))
  ctl <- adghControl(studies = st, cov_integration = "sparse", print = 0L,
                     covMethod = "none")
  pin <- admixr2:::.admDriverPinfo(ui, ctl)
  out <- suppressMessages(admixr2:::.admCheckCovariates(ui, pin, st, "adgh"))
  sg  <- out$s[[".adm_cov_sparse"]]
  expect_false(is.null(sg))
  expect_equal(sum(sg$W), 1, tolerance = 1e-12)
  # numeric only, so it survives serialisation to a restart worker by value
  expect_false(any(vapply(sg, is.function, logical(1))))
})

test_that("covDist honours `dist` for every vocabulary, without partial matching", {
  # `$` PARTIAL-MATCHES on lists, and .admPopSpec returns meanlog/sdlog on the lnorm branch -- where BOTH
  # `$mean` and `$sd` match. A defensive `if (!is.null(sp$mean))` therefore fired on a LOGNORMAL spec and
  # rewrote it as list(mu = meanlog, sd = sdlog): covDist(WT = c(mean = 72, cv = 22), dist = "lnorm") came
  # back as a NORMAL margin centred at 4.25 kg, so every quadrature node sat near 4 instead of near 72 and
  # (WT/70)^0.75 was evaluated at 0.06. Finite, plausible, and nowhere near the declared cohort.
  for (v in list(c(mean = 72, cv = 22), c(mean = 72, sd = 16),
                 c(median = 70, iqr = c(60, 84)))) {
    n <- covDist(WT = v, dist = "normal")$WT
    l <- covDist(WT = v, dist = "lnorm")$WT
    expect_named(n, c("mu", "sd"), ignore.order = TRUE)
    expect_named(l, c("meanlog", "sdlog"), ignore.order = TRUE)
    # the lognormal is on the LOG scale, so its centre is log(natural centre)
    expect_lt(l$meanlog, 6)
    expect_gt(n$mu, 50)
  }
})

test_that("a single named proportion is a BINARY covariate, not a constant", {
  # `SEX = c(male = 0.55)` is how a baseline table prints it and what admPopulation() documents -- but
  # covDist() normalised it like a level vector, dividing 0.55 by itself: ONE level at probability 1, i.e. a
  # constant. The covariate then contributed no variation, its coefficient was not identified, and covDraw()
  # returned the same value for every subject. Found by a worked example, not by a test, which is why this one
  # exists.
  s <- covDist(SEX = c(male = 0.55))$SEX
  expect_equal(s$values, c(0, 1))
  expect_equal(s$probs,  c(0.45, 0.55))
  # the two parsers must agree on the form the docs teach most
  expect_equal(s[c("values", "probs")],
               admixr2:::.admPopSpec(c(male = 0.55), "SEX", "norm")[c("values", "probs")])
  # and the draws actually vary
  expect_setequal(unique(covDraw(covDist(SEX = c(male = 0.55)), n = 60)[, "SEX"]),
                  c(0, 1))
  # multi-level specs are untouched
  expect_equal(covDist(SEX = c(female = 0.45, male = 0.55))$SEX$probs, c(0.45, 0.55))
  expect_equal(covDist(G = c(a = 1, b = 2, c = 1))$G$probs, c(0.25, 0.5, 0.25))
  # ... and a LONE value that is not a proportion has no second number to be normalised against, so reading it
  # as one gave P(absent) = -69. The multi-level branch is normalised by its sum and so absorbed any scale
  # silently; this one cannot, and a negative quadrature weight is a wrong answer rather than an error.
  expect_error(covDist(WT = c(kg = 70)), "[[]0, 1[]]")
  expect_error(covDist(SEX = c(male = -0.2)), "negative")
})
