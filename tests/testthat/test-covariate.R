# Tier 1 -- covariate marginalisation helpers (no rxode2).
#
# gl/gh/taylor node quadrature was removed, and so were the "collapse" and "uq"
# paths (see the header of R/covariate.R). What remains is the general per-row /
# product-grid route, the shift route, and the refusal that keeps old node-style
# study lists from being fitted as something else.

# ---- routing -----------------------------------------------------------------

.cov_pinfo <- function(n_eta = 1L, cov_integration = NULL)
  list(n_eta = n_eta,
       eta_col_names = if (n_eta == 1L) "eta.cl" else c("eta.cl", "eta.v"),
       cov_integration = cov_integration,
       chol_diag = rep(TRUE, n_eta),
       chol_i = seq_len(n_eta), chol_j = seq_len(n_eta))

.cov_ui <- function(cov = "WT",
                    expr = list(quote(cl <- exp(tcl + tcov * WT + eta.cl))))
  list(allCovs = cov, lstExpr = expr)

test_that(".admCovCols adds ONLY declared covariates, never a blanket fill", {
  m <- matrix(1, 3L, 1L, dimnames = list(NULL, "tcl"))
  # `vb` is a hard-coded model constant and `lam` an estimated TBS lambda:
  # both are model params, neither is a covariate, so neither may be added.
  got <- admixr2:::.admCovCols(m, c("tcl", "wt", "vb", "lam"), list(wt = 70))
  expect_equal(colnames(got), c("tcl", "wt"))
  expect_equal(unname(got[, "wt"]), rep(70, 3L))

  # a covariate the model does not read is not added either
  expect_equal(colnames(admixr2:::.admCovCols(m, c("tcl", "vb"), list(wt = 70))),
               "tcl")
  # no covariates at all -> untouched
  expect_identical(admixr2:::.admCovCols(m, c("tcl", "wt"), NULL), m)
})

test_that(".admStudyL is the plain Cholesky for EVERY study", {
  # The "collapse" path -- chol(Omega + J Sigma_a J'), which folded a normal
  # covariate's variance into Omega -- was retired, so no study-specific Omega
  # remains. Pinned because R/admc.R and R/plot.R call this on every objective
  # evaluation: an inflated Omega reaching them again would be a plausible fit
  # of a different model, which is exactly how the collapse used to fail.
  pars <- list(L = t(chol(matrix(0.09, 1L, 1L))), struct = c(tcov = 0.75))
  expect_identical(admixr2:::.admStudyL(pars, .cov_pinfo(), list()), pars$L)
  expect_identical(
    admixr2:::.admStudyL(pars, .cov_pinfo(),
                         list(cov_dist = list(WT = list(mu = 0, sd = 0.6)))),
    pars$L)
})

test_that(".admCheckCovariates accepts a supported mu-referenced covariate", {
  st <- list(a = list(cov_dist = list(WT = list(mu = 0, sd = 0.6))))
  expect_silent(admixr2:::.admCheckCovariates(.cov_ui(), .cov_pinfo(), st, "none"))
})

test_that(".admCheckCovariates is a no-op when no study declares cov_dist", {
  expect_silent(admixr2:::.admCheckCovariates(.cov_ui(), .cov_pinfo(),
                                              list(a = list()), "sens"))
})

test_that("`grad` no longer enters the covariate path choice", {
  # It used to: "collapse" and "uq" had no derivative, so a gradient mode forced
  # "rows" and grad = "none" got the closed form. Both of those paths are gone
  # and both survivors are differentiable, so every mode must land identically.
  # This is what made "collapse" unreachable in a real fit -- every estimator
  # defaults to a gradient -- and the reason it was retired rather than repaired.
  ok_st <- list(a = list(cov = list(WT = 0),
                         cov_dist = list(WT = list(mu = 0, sd = 0.6))))
  for (g in c("none", "sens", "analytical", "fd", "cfd"))
    expect_identical(
      admixr2:::.admCheckCovariates(.cov_ui(), .cov_pinfo(), ok_st, g)$a$.adm_cov_path,
      "rows")
})

test_that(".admCheckCovariates routes to the general path by default", {
  ok_st <- list(a = list(cov = list(WT = 0), cov_dist = list(WT = list(mu = 0, sd = 0.6))))
  path <- function(ui = .cov_ui(), pi = .cov_pinfo(), st = ok_st)
    admixr2:::.admCheckCovariates(ui, pi, st, "none")$a$.adm_cov_path

  # the bare theta*COV product the retired collapse needed -- now "rows" like
  # everything else. adgh reaches the same moments through the shift path when
  # cov_integration says so, and that route is admitted numerically rather than
  # by rxode2's muRefCovariateDataFrame (which recognises 1 of 8 realistic
  # covariate parameterisations).
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

test_that('cov_integration = "auto" falls back where "shift" errors', {
  # `v <- exp(tv) * WT^vwt` has no random effect in the assignment carrying the
  # covariate, so there is no column to shift. .admShiftSpec() decides that from
  # the model text alone, which is why it is reachable without a compiled model;
  # everything subtler is decided by .admShiftVerify() against the solver.
  ui <- .cov_ui(expr = list(quote(v <- exp(tv) * WT^vwt)))
  st <- list(a = list(cov_dist = list(WT = list(mu = 70, sd = 10))))

  expect_error(
    admixr2:::.admCheckCovariates(ui, .cov_pinfo(cov_integration = "shift"),
                                  st, "none"),
    "exactly one random effect")

  # "auto" is a SPEED lever, and the fallback is the more accurate path, so a
  # refusal must never be an error -- only a message and the reason recorded.
  expect_message(
    got <- admixr2:::.admCheckCovariates(
      ui, .cov_pinfo(cov_integration = "auto"), st, "none"),
    "auto")
  expect_identical(got$a$.adm_cov_path, "rows")
  expect_match(got$a$.adm_cov_shift_why, "exactly one random effect")
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
    list(a = list(cov_dist = list(WT = list(mu = 70, sd = 10)))), "none")
  expect_equal(st$a$cov$WT, 70)
})

test_that(".admCheckCovariates still errors on genuinely unsupportable input", {
  ok_st <- list(a = list(cov = list(WT = 0), cov_dist = list(WT = list(mu = 0, sd = 0.6))))

  # covariate the model never reads -- almost always a typo
  expect_error(
    admixr2:::.admCheckCovariates(.cov_ui(cov = "AGE"), .cov_pinfo(), ok_st, "none"),
    "which the model never reads")

  # distributions we cannot draw from
  for (spec in list(list(mu = 0), list(mu = 0, sd = 0), list(mu = 0, sd = NA_real_),
                    list(values = numeric(0))))
    expect_error(
      admixr2:::.admCheckCovariates(.cov_ui(), .cov_pinfo(),
                                    list(a = list(cov = list(WT = 0),
                                                  cov_dist = list(WT = spec))), "none"),
      "not a supported distribution")
})

test_that(".admShiftNodes solves F_u(u) = Phi(z) exactly, and smoothly", {
  # The shift path's inversion, and the direct replacement for the retired
  # .admCovUQuantile: Newton on the exact mixture CDF of u = Delta(a) + eta,
  # never interpolation on a grid. A grid makes u piecewise-linear in the
  # parameters, which is the mismatch that stops an analytic gradient being
  # consistent with the objective it differentiates.
  g   <- admixr2:::.adghNodes1(32L)
  D   <- 0.75 * log(exp(log(72) + 0.28 * g$x) / 72)   # allometric Delta
  om  <- 0.30
  # 25, not more: beyond ~8 SD pnorm() saturates at exactly 1, so the target
  # probability carries no information and Newton leaves those nodes wherever it
  # found them (documented in .admShiftNodes -- their GH weight is ~1e-23, so
  # they cannot move a moment). .adghNodes1(41) reaches there; 25 does not.
  n_u <- 25L
  un  <- admixr2:::.admShiftNodes(D, g$w, om, n_u)
  tg  <- pnorm(admixr2:::.adghNodes1(n_u)$x)

  Fu <- as.numeric(pnorm(outer(un$u, D, "-") / om) %*% g$w)
  expect_lt(max(abs(Fu - tg)), 1e-10)         # solves its defining equation
  # monotone IN z. The nodes come back in .adghNodes1's own order, which is
  # descending -- asserting diff(u) > 0 instead pins the ordering convention of
  # a helper this function only borrows.
  expect_true(all(diff(un$u[order(tg)]) > 0))
  expect_true(all(is.finite(un$u)))
  expect_equal(sum(un$w), 1)

  # SHIFTING every Delta by eps must shift u by exactly eps: an exact identity,
  # and one a grid interpolant cannot reproduce without quantisation.
  eps <- 1e-7
  u2  <- admixr2:::.admShiftNodes(D + eps, g$w, om, n_u)$u
  expect_equal(mean((u2 - un$u) / eps), 1, tolerance = 1e-6)
  expect_lt(sd((u2 - un$u) / eps), 1e-6)

  # u's first two moments are known in closed form (E[Delta], Var(Delta)+om^2),
  # so the node set must reproduce them -- that is what the quadrature is for.
  mD <- sum(g$w * D); vD <- sum(g$w * D^2) - mD^2
  expect_equal(sum(un$w * un$u), mD, tolerance = 1e-6)
  expect_equal(sum(un$w * (un$u - mD)^2), vD + om^2, tolerance = 1e-4)

  # the vector driver must agree with the scalar one at m = 1
  m1 <- admixr2:::.admShiftNodesMulti(matrix(D, ncol = 1L), g$w, om, n_u)
  expect_equal(as.numeric(m1$u), un$u)
  expect_equal(m1$w, un$w)
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

  # a covariate on a parameter with NO random effect is NOT on the ridge: there
  # is no omega for its variance to hide in, so one population is enough
  ui2 <- .cov_ui(expr = list(quote(v <- exp(tv) * WT^vwt)))
  expect_silent(admixr2:::.admWarnCovIdentifiability(ui2, .cov_pinfo(), one))
})

test_that("conditioned strata identify the coefficient within one source", {
  # A stratum is a covariate distribution at zero spread, and supplies its own
  # equation in (theta, gamma) exactly as a second study's distribution does.
  # Reading `cov_dist` alone called this unidentified -- a false alarm on the
  # one design that identifies gamma without a between-source contrast.
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
  # eq:ridge holds only while the model sees gamma*a + omega*b. A covariate
  # read by a second assignment restores the separate dependence on (a, b), so
  # it is identified from one population and must not be warned about.
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

# -- The removed node route must fail loudly, not quietly become something else -

test_that(".admRefuseNodeStudies rejects the removed node inputs", {
  # `weight` was the node's combination coefficient. Ignoring it would fit the
  # nodes as an ordinary unweighted study list -- a different objective, and a
  # perfectly plausible number.
  w <- list(a = list(n = 10, weight = 2), b = list(n = 10, weight = -3))
  expect_error(admixr2:::.admRefuseNodeStudies(w), "was removed")
  expect_error(admixr2:::.admRefuseNodeStudies(w), "'a', 'b'")
  m <- list(a = list(n = 10, cov_method = "gh"))
  expect_error(admixr2:::.admRefuseNodeStudies(m), "cov_method")

  # ordinary studies pass, including an explicit weight of 1 and the marginal
  # method spelled out
  expect_null(admixr2:::.admRefuseNodeStudies(list(a = list(n = 10))))
  expect_null(admixr2:::.admRefuseNodeStudies(list(a = list(n = 10, weight = 1))))
  expect_null(admixr2:::.admRefuseNodeStudies(
    list(a = list(n = 10, cov_method = "marginal"))))
})

test_that("the message points at BOTH replacements, not just cov_dist", {
  # A user with stratified summaries must be told to pass them as ordinary
  # studies with their own n -- that route is the same likelihood the node
  # code computed, and it is the one they should reach for.
  msg <- tryCatch(admixr2:::.admRefuseNodeStudies(list(a = list(n = 1, weight = 2))),
                  error = conditionMessage)
  expect_match(msg, "cov_dist")
  # strata must carry their own DISTRIBUTION, not a point `cov`: plugging in the
  # stratum mean is the ecological plug-in and biases the coefficient upward
  # (+17% at 2 strata, +4.3% at 4). The message has to say so.
  expect_match(msg, "its own `n` and its own `cov_dist`")
  expect_match(msg, "Do NOT give a stratum a point `cov`")
})

test_that("the friendly cov_dist grammar is EXACTLY the hand-written one", {
  # mean/sd on the natural scale, and `cor`, exist so a user can transcribe a
  # baseline-characteristics table directly. They must therefore produce the
  # identical distribution to the canonical spelling -- a convenience that
  # shifted the covariate distribution even slightly would move every estimate
  # while looking like a formatting choice.
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
  # lognormal margins here on purpose: weight and creatinine clearance are
  # positive, right-skewed physiological quantities, and a Gaussian copula on
  # lognormal margins has correlation exactly rho on the LOG scale
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
  # a Gaussian copula on lognormal margins has correlation exactly rho on the
  # LOG scale, which is the check that `cor` means what a reader expects
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

# =============================================================================
# Second-order Taylor covariate integration (cov_integration = "taylor")
# =============================================================================
#
# The reference is EXACT nested Gauss-Hermite over (covariate, eta) on an
# analytic one-compartment solution, computed here in plain R. Nothing is pinned
# against the package's own output: the failure mode of a moment expansion is a
# finite, plausible, biased number, so a self-comparison catches nothing.

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

# Conditional moments at ONE covariate value -- the ecological plug-in, and the
# thing the expansion differences.
.tay_cond <- function(a, tcov) {
  Y  <- .tay_conc(exp(.tay_TCL + tcov * a + .tay_OM * .tay_QE$x))
  mu <- as.numeric(crossprod(.tay_QE$w, Y)); Yc <- sweep(Y, 2L, mu)
  list(E = mu, V = t(Yc) %*% (Yc * .tay_QE$w))
}

# The PACKAGE's own path: build the design, lay the solved matrix out exactly as
# .adghGrid() does (design point slowest, eta node fastest), assemble.
.tay_pkg <- function(tcov, mu_a, sd_a, hfrac = 1) {   # package default radius
  td  <- admixr2:::.admCovTaylorDesign(list(A = list(mu = mu_a, sd = sd_a)), hfrac)
  nq  <- length(.tay_QE$x)
  tay <- admixr2:::.admCovTaylorRows(td, .tay_QE$w, nq)
  W   <- as.numeric(outer(.tay_QE$w, td$c))
  a   <- rep(td$X[, 1L], each = nq)
  eta <- rep(.tay_OM * .tay_QE$x, times = td$n_pt)
  cp  <- .tay_conc(exp(.tay_TCL + tcov * a + eta))
  c(admixr2:::.adghStructMoments(cp, W, tay), list(design = td, rows = tay))
}

.tay_relerr <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-12))

test_that(".admCovVarOf reports Var(a) on the covariate's OWN scale", {
  # NOT .admCovSdOf()^2: that one is deliberately the log-scale spread for a
  # lognormal spec, which paired with the natural-scale mean would expand about
  # the right point with the wrong second moment.
  expect_equal(admixr2:::.admCovVarOf(list(mu = 3, sd = 2)), 4)
  ml <- log(72); sl <- 0.28
  expect_equal(admixr2:::.admCovVarOf(list(meanlog = ml, sdlog = sl)),
               (exp(sl^2) - 1) * exp(2 * ml + sl^2))
  expect_false(isTRUE(all.equal(
    admixr2:::.admCovVarOf(list(meanlog = ml, sdlog = sl)),
    admixr2:::.admCovSdOf(list(meanlog = ml, sdlog = sl))^2)))
  # discrete: a two-point 0/1 covariate with p = 0.25 has variance 0.1875
  expect_equal(admixr2:::.admCovVarOf(list(values = c(0, 1), probs = c(3, 1))),
               0.1875)
  # a quantile spec is integrated on the same 1024 midpoints .admCovMeanOf uses,
  # which truncates the tails -- so it is a touch under the true 4, by design
  expect_equal(admixr2:::.admCovVarOf(
    list(quantile = function(u) stats::qnorm(u, 5, 2))), 4, tolerance = 5e-3)
})

test_that(".admCovTaylorDesign lays out 1 + 2p points with weights summing to 1", {
  # hfrac is a MULTIPLIER on the moment-matched radius sqrt(3*lambda), so the
  # default 1 puts the design points where 3-point Gauss-Hermite does: at
  # +/- sqrt(3) latent SDs, weights 1/6 and centre 1 - d/3.
  h  <- sqrt(3)
  td <- admixr2:::.admCovTaylorDesign(
    list(A = list(mu = 2, sd = 0.4), B = list(mu = -1, sd = 3)), hfrac = 1)
  expect_identical(td$n_pt, 5L)                       # 1 + 2*2
  expect_identical(td$names, c("A", "B"))
  # centre, then (+h, -h) per covariate with every OTHER covariate at its mean
  expect_equal(td$X[1L, ], c(A = 2, B = -1))
  expect_equal(unname(td$X[td$ip, ]),
               rbind(c(2 + 0.4*h, -1), c(2, -1 + 3*h)))
  expect_equal(unname(td$X[td$im, ]),
               rbind(c(2 - 0.4*h, -1), c(2, -1 - 3*h)))
  # E_marg = sum_k c_k E_k is an average, so the coefficients sum to one
  expect_equal(sum(td$c), 1)
  # 1 - 2/h^2 = 1 - 2/3, then 1/(2h^2) = 1/6 each
  expect_equal(td$c, c(1 - 2/3, rep(1/6, 4L)))
  # a smaller radius is still a valid degree-3 rule, just not moment-matched
  td5 <- admixr2:::.admCovTaylorDesign(
    list(A = list(mu = 2, sd = 0.4), B = list(mu = -1, sd = 3)), hfrac = 0.5)
  expect_equal(sum(td5$c), 1)
  # `var` is on the LATENT scale the expansion runs in, so it is 1 per
  # direction for independent covariates -- the covariate-scale spread enters
  # through the design POINTS instead (mu +/- hfrac*sd, above)
  expect_equal(unname(td$var), c(1, 1))
})

test_that("the Taylor design handles DEPENDENT covariates in the eigenbasis", {
  # 1 + 2p points at ANY correlation: the two Sigma-weighted terms of the
  # expansion are directional, so differencing along the eigenvectors of Sigma
  # carries the off-diagonal at no extra cost. (This used to be refused.)
  cd <- list(A = list(mu = 0, sd = 2, dist = "normal"),
             B = list(mu = 0, sd = 3, dist = "normal"), cor = 0.5)
  td <- admixr2:::.admCovTaylorDesign(cd, hfrac = 0.5)
  expect_identical(td$n_pt, 5L)                       # unchanged by rho
  expect_equal(sum(td$c), 1)
  expect_equal(crossprod(td$dir), diag(1, 2), tolerance = 1e-12)
  # the design must REPRODUCE Sigma: sum_k lam_k v_k v_k' == Sigma
  rec <- Reduce(`+`, lapply(seq_along(td$var),
    function(k) td$var[k] * tcrossprod(td$dir[, k])))
  expect_equal(rec, unname(td$Sigma), tolerance = 1e-10)
  # off the coordinate axes -- both covariates move together along direction 1
  expect_true(all(abs(td$X[td$ip[1L], ] - td$mu) > 1e-8))
  # and the three spellings of dependence must agree exactly
  eig <- function(x) admixr2:::.admCovTaylorDesign(x, 0.5)$var
  base <- list(A = list(mu = 0, sd = 2, dist = "normal"),
               B = list(mu = 0, sd = 3, dist = "normal"))
  expect_equal(eig(c(base, list(rho = 0.5))), td$var, tolerance = 0)
  expect_equal(eig(c(base, list(Sigma = matrix(c(4, 3, 3, 9), 2, 2)))), td$var,
               tolerance = 0)
})

test_that("the Taylor design REFUSES what it cannot expand, rather than approximating", {
  # (near-)collinear covariates: h = hfrac*sqrt(lambda) collapses, so the
  # second difference along that direction is cancellation, not a derivative.
  # Checked on the CORRELATION matrix, so covariates that merely differ in
  # UNITS are not caught by it.
  expect_error(admixr2:::.admCovTaylorDesign(
    list(A = list(mu = 0, sd = 1, dist = "normal"),
         B = list(mu = 0, sd = 1, dist = "normal"), cor = 1 - 1e-14)),
    "collinear")
  expect_silent(admixr2:::.admCovTaylorDesign(
    list(A = list(mu = 0, sd = 1e4, dist = "normal"),
         B = list(mu = 0, sd = 1e-4, dist = "normal"), cor = 0.3)))
  # a discrete covariate is ENUMERATED, not expanded -- what is refused is a
  # discrete covariate DEPENDENT on a continuous one, where a level is a
  # truncation of the latent normal rather than a point, so the continuous
  # conditional differs cell by cell.
  Rz <- matrix(c(1, 0.3, 0.3, 1), 2L, 2L,
               dimnames = list(c("WT", "SEX"), c("WT", "SEX")))
  expect_error(admixr2:::.admCovTaylorDesign(
    list(WT = list(mu = 70, sd = 15), SEX = list(values = c(0, 1)),
         latentR = Rz)), "DEPENDENT")
  # a degenerate covariate makes the differencing step zero
  expect_error(admixr2:::.admCovTaylorDesign(
    list(A = list(quantile = function(u) rep(5, length(u))))), "no spread")
  # and the step itself has to be a usable number
  expect_error(admixr2:::.admCovTaylorDesign(list(A = list(mu = 0, sd = 1)), 0),
               "positive finite")
  # a diagonal Sigma is fine -- it is only the off-diagonal that is refused
  Sd <- diag(c(1, 4)); dimnames(Sd) <- list(c("A", "B"), c("A", "B"))
  expect_silent(admixr2:::.admCovTaylorDesign(
    list(A = list(mu = 0, sd = 1), B = list(mu = 0, sd = 2), Sigma = Sd)))
})

test_that("the Taylor design ENUMERATES a discrete covariate exactly", {
  # A step function has no curvature, but it does not need one: its levels and
  # probabilities ARE the rule.  All-discrete therefore reduces to the exact
  # enumeration -- no expansion, no design points off the levels.
  td <- admixr2:::.admCovTaylorDesign(
    list(SEX = list(values = c(0, 1), probs = c(0.55, 0.45))))
  expect_equal(td$n_pt, 2L)
  expect_equal(as.numeric(td$X[, "SEX"]), c(0, 1))
  expect_equal(td$c, c(0.55, 0.45))            # the level probabilities exactly
  expect_equal(td$n_cell, 2L)
  expect_equal(td$n_cpt, 1L)                   # no cubature points at all

  # mixed: L cells x (1 + 2 p_continuous), weights the product, and the signed
  # combination still sums to one.
  tm <- admixr2:::.admCovTaylorDesign(
    list(WT = list(mu = 70, sd = 15), SEX = list(values = c(0, 1),
                                                 probs = c(0.55, 0.45))))
  expect_equal(tm$n_pt, 6L)
  expect_equal(tm$n_cell, 2L)
  expect_equal(tm$n_cpt, 3L)
  expect_equal(sum(tm$c), 1)
  expect_equal(sort(unique(tm$X[, "SEX"])), c(0, 1))
  # the derivative directions are per CELL, so the rank term is p_l-weighted
  expect_equal(length(tm$ip), 2L)              # 2 cells x 1 continuous direction
  expect_equal(unname(tm$var), c(0.55, 0.45))  # p_l * lambda_j, lambda = 1

  # and the row map carries one extra column per cell: E_l - E, the exact
  # between-cell term of the law of total variance.  Each such column must have
  # zero total weight, or it would not be a centred contrast.
  tr <- admixr2:::.admCovTaylorRows(tm, c(0.25, 0.5, 0.25), 3L)
  # 2 first-difference + 2 second-difference + 2 between-cell columns
  expect_equal(ncol(tr$Dw), 6L)
  expect_equal(length(tr$var), 6L)
  # every column is a CENTRED contrast: a first difference (+1,-1), a second
  # difference (+1,-2,+1) and a cell mean minus the global mean all sum to zero
  expect_equal(unname(colSums(tr$Dw)), rep(0, 6L))
  # the second-difference block carries p_l * lambda^2 / 2
  expect_equal(unname(tr$var[3:4]), c(0.55, 0.45) / 2)
})

test_that("the Taylor moments match the EXACT marginal, and beat the plug-in", {
  # ratio = tcov * sd_a / omega -- the covariate-induced spread relative to the
  # IIV. Reproduces validation/taylor-corrected.R's accuracy table.
  mu_a <- 0; sd_a <- 0.40
  want <- list("0.1" = c(1e-6, 1e-4), "0.5" = c(1e-4, 1e-2), "1" = c(2e-3, 1e-1))
  for (ratio in c(0.1, 0.5, 1)) {
    tcov <- ratio * .tay_OM / sd_a
    ex   <- .tay_exact(tcov, mu_a, sd_a)
    ta   <- .tay_pkg(tcov, mu_a, sd_a)
    cd   <- .tay_cond(mu_a, tcov)                  # the ecological plug-in
    lim  <- want[[as.character(ratio)]]
    expect_lt(.tay_relerr(ta$mu, ex$E), lim[[1L]])
    expect_lt(.tay_relerr(ta$V,  ex$V), lim[[2L]])
    # ... and it is an IMPROVEMENT on solving at the covariate mean, which is
    # what "we already have a covariate value" would silently give you
    expect_lt(.tay_relerr(ta$V, ex$V), .tay_relerr(cd$V, ex$V) / 10)
  }
})

test_that("the rank-one sd^2 g' g'^T term is PRESENT in the Taylor covariance", {
  # This is the term that carries the covariate effect into V at all; without it
  # the expansion returns E_a[Cov_eta(f|a)], the average WITHIN-covariate
  # covariance, and Cov_a(g(a)) is simply missing. Same class of mistake as
  # writing V_struct + Sigma(mu) for the residual: finite, plausible, biased low.
  mu_a <- 0; sd_a <- 0.40
  for (ratio in c(0.5, 1)) {
    tcov <- ratio * .tay_OM / sd_a
    ex   <- .tay_exact(tcov, mu_a, sd_a)
    ta   <- .tay_pkg(tcov, mu_a, sd_a)
    # Remove exactly that term and nothing else. Measured: the covariance error
    # goes 5.5e-03 -> 2.0e-01 at ratio 0.5 (36x) and 5.3e-02 -> 5.0e-01 at
    # ratio 1 (9.5x), and the absolute floor says the term is a fifth of V
    # rather than a rounding correction.
    V_no <- ta$V - crossprod(ta$dE, ta$rows$var * ta$dE)
    expect_gt(.tay_relerr(V_no, ex$V), 8 * .tay_relerr(ta$V, ex$V))
    expect_gt(.tay_relerr(V_no, ex$V), 0.19)
    # Cov_a(g) is now TWO blocks: the rank-p first-difference term
    # sum_j lam_j g'_j g'_j' and the second-difference term
    # (1/2) sum_j lam_j^2 g''_j g''_j', which is the next order of the same
    # expansion and costs no extra design point. Split them to check each.
    i1 <- 1L; i2 <- 2L                       # one covariate -> one direction each
    rk <- crossprod(ta$dE[i1, , drop = FALSE],
                    ta$rows$var[i1] * ta$dE[i1, , drop = FALSE])
    rk2 <- crossprod(ta$dE[i2, , drop = FALSE],
                     ta$rows$var[i2] * ta$dE[i2, , drop = FALSE])
    expect_equal(qr(rk)$rank, 1L)
    expect_equal(qr(rk2)$rank, 1L)
    expect_gt(min(eigen(rk2, symmetric = TRUE, only.values = TRUE)$values), -1e-12)
    expect_gt(min(eigen(rk, symmetric = TRUE, only.values = TRUE)$values), -1e-12)
    # ... and it is exactly Var(z) * (dg/dz)(dg/dz)', with the derivative a
    # central difference IN THE LATENT VARIABLE. For a normal margin z and a
    # differ only by the linear map a = mu + sd*z, so dg/dz = sd * dg/da and
    # the rank term is the same object either way -- which is the check.
    h  <- sqrt(3) * sd_a                    # the moment-matched radius
    gp <- (.tay_cond(mu_a + h, tcov)$E - .tay_cond(mu_a - h, tcov)$E) / (2 * h)
    expect_equal(as.numeric(ta$dE[i1, ]), sd_a * gp)
    expect_equal(rk, sd_a^2 * outer(gp, gp))
  }
})

test_that("the Taylor moments are the second-order expansion, term for term", {
  # E = g(mu) + (sd^2/2) g''(mu), differenced at h = hfrac*sd -- checked against
  # the three conditional moment sets directly, not against the package's own
  # assembly, so a wrong stride or a mis-signed weight cannot pass.
  mu_a <- 0.2; sd_a <- 0.35; tcov <- 0.9; h <- sqrt(3) * sd_a
  m0 <- .tay_cond(mu_a, tcov); mp <- .tay_cond(mu_a + h, tcov)
  mm <- .tay_cond(mu_a - h, tcov)
  dE  <- (mp$E - mm$E) / (2 * h)                       # g'
  d2E <- (mp$E - 2 * m0$E + mm$E) / h^2                # g''
  want_E <- m0$E + 0.5 * sd_a^2 * d2E
  # Cov_a(g) to second order: the rank-p first-derivative term PLUS
  # (1/2) lam^2 g'' g'', which uses the same three design points.
  want_V <- m0$V + 0.5 * sd_a^2 * (mp$V - 2 * m0$V + mm$V) / h^2 +
    sd_a^2 * outer(dE, dE) + 0.5 * sd_a^4 * outer(d2E, d2E)
  ta <- .tay_pkg(tcov, mu_a, sd_a)
  expect_equal(ta$mu, want_E)
  expect_equal(ta$V,  want_V)
})

test_that(".adghStructMoments is BIT-IDENTICAL to the pooled formulas without a design", {
  # cov_integration = "quadrature" must not move by one ulp, and tay = NULL is
  # the only thing that path passes -- this is the whole of what it computes.
  set.seed(11)
  cp <- matrix(stats::rnorm(40L * 6L, 5, 2), 40L, 6L)
  W  <- stats::runif(40L); W <- W / sum(W)
  got <- admixr2:::.adghStructMoments(cp, W, NULL)
  mu  <- as.numeric(crossprod(W, cp))
  cpc <- sweep(cp, 2L, mu)
  expect_identical(got$mu,  mu)
  expect_identical(got$cpc, cpc)
  expect_identical(got$V,   crossprod(cpc, W * cpc))
  expect_null(got$dE)
})

test_that("adghControl defaults to quadrature and validates the Taylor step", {
  expect_identical(adghControl()$cov_integration, "quadrature")
  # cov_taylor_h is now a multiplier on the moment-matched radius, default 1
  expect_identical(adghControl()$cov_taylor_h, 1)
  expect_identical(adghControl(cov_integration = "shift")$cov_integration,
                   "shift")
  expect_identical(adghControl(cov_integration = "taylor")$cov_integration,
                   "taylor")
  expect_error(adghControl(cov_integration = "gh"))
  expect_error(adghControl(cov_taylor_h = 0))
  expect_error(adghControl(cov_taylor_h = c(0.2, 0.4)))
  # the option is adgh's alone, so a stray argument elsewhere must be an error
  # rather than a silently ignored request for a different integration
  expect_error(admControl(cov_integration = "taylor"), "unused argument")
})

# =============================================================================
# Dependent covariates on the DETERMINISTIC (quadrature / Taylor) paths
# =============================================================================

test_that(".admCovGrid integrates a DEPENDENT distribution on the u-space grid", {
  # A copula maps INDEPENDENT uniforms to dependent values, so a product rule
  # in u is exact whatever the dependence -- the weights genuinely factorise
  # there. adgh used to refuse a `joint` spec on the grounds that its product
  # grid assumed independence; that is true of a grid over covariate MARGINS
  # and false of a grid over u.
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
    # normal margins + Gaussian copula: the closed form is exact, so this is a
    # real check and not a comparison against another approximation
    expect_equal(unname(m), c(70, 50), tolerance = 1e-6)
    expect_equal(unname(sqrt(diag(S))), c(10, 12), tolerance = 1e-6)
    expect_equal(S[1, 2] / sqrt(S[1, 1] * S[2, 2]), 0.7, tolerance = 1e-6)
  }
})

test_that("the grid and the per-subject sampler see the SAME distribution", {
  # admc draws rows, adgh builds a grid. If these disagreed, the two estimators
  # would fit different data and only a side-by-side run would reveal it.
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

test_that("`cor`, `rho` and `Sigma` are ONE statement, honoured by every path", {
  # They used to diverge: `rho` built a copula for the retired collapse and a
  # DIAGONAL grid for everything else, so the correlation was present in one
  # path and silently absent in the other. All three spellings are still
  # accepted -- published specs are written every which way -- so all three must
  # still land on the same latent correlation.
  base <- list(A = list(mu = 0, sd = 2, dist = "normal"),
               B = list(mu = 0, sd = 3, dist = "normal"))
  sig <- function(x) admixr2:::.admCovDistMoments(x)$Sigma
  s_cor <- sig(c(base, list(cor   = 0.5)))
  expect_equal(sig(c(base, list(rho   = 0.5))), s_cor, tolerance = 0)
  expect_equal(sig(c(base, list(Sigma = matrix(c(4, 3, 3, 9), 2, 2)))), s_cor,
               tolerance = 0)
  expect_equal(unname(s_cor), matrix(c(4, 3, 3, 9), 2, 2), tolerance = 1e-2)
  # ... and all three canonicalise to the same latent correlation matrix, which
  # is the single object every consumer (the joint sampler, the product grid,
  # the Taylor design, the shift grid) reads.
  lr <- function(x) admixr2:::.admCovDistCanon(x)[["latentR"]]
  expect_equal(lr(c(base, list(rho = 0.5))), lr(c(base, list(cor = 0.5))))
  expect_equal(lr(c(base, list(Sigma = matrix(c(4, 3, 3, 9), 2, 2)))),
               lr(c(base, list(cor = 0.5))))
})

test_that(".admCovDistMoments is diagonal exactly when the covariates are independent", {
  ind <- list(A = list(mu = 1, sd = 2), B = list(meanlog = 0, sdlog = 0.3))
  mi  <- admixr2:::.admCovDistMoments(ind)
  expect_true(mi$diagonal)
  expect_equal(mi$Sigma[1, 2], 0, tolerance = 0)
  # the lognormal entry must be the NATURAL-scale variance, not sdlog^2
  expect_equal(mi$Sigma[2, 2], (exp(0.3^2) - 1) * exp(0.3^2), tolerance = 1e-12)
  dep <- admixr2:::.admCovDistCanon(c(ind, list(cor = 0.4)))
  expect_false(admixr2:::.admCovDistMoments(dep)$diagonal)
})

test_that("metadata keys are never mistaken for covariates", {
  # each of these used to be a different internal error: `$ operator is invalid
  # for atomic vectors` on rho, `object of type 'closure' is not subsettable`
  # on the joint sampler.
  cd <- list(A = list(mu = 0, sd = 1), B = list(mu = 0, sd = 1), rho = 0.3)
  expect_identical(admixr2:::.admCovSpecNames(cd), c("A", "B"))
  expect_silent(admixr2:::.admCovGrid(admixr2:::.admCovDistCanon(cd), 3L))
  expect_identical(
    colnames(admixr2:::.admCovGrid(admixr2:::.admCovDistCanon(cd), 3L)$X),
    c("A", "B"))
})

test_that("a NAMED correlation matrix is ordered to the declared covariates", {
  # The reorder is the only thing standing between a user writing their cor
  # matrix in a different order and admixr2 correlating the wrong PAIR -- which
  # is finite, plausible and silent. It briefly became dead code by sitting
  # after a `dimnames(R) <- NULL`.
  R  <- matrix(c(1, 0.8, 0.8, 1), 2, 2,
               dimnames = list(c("B", "A"), c("B", "A")))
  cd <- admixr2:::.admCovDistCanon(
    list(A = list(mu = 0, sd = 1, dist = "normal"),
         B = list(mu = 0, sd = 5, dist = "normal"), cor = R))
  S <- admixr2:::.admCovDistMoments(cd)$Sigma
  # variances must land on the covariate they were DECLARED for
  expect_equal(unname(S[1, 1]), 1, tolerance = 5e-3)
  expect_equal(unname(S[2, 2]), 25, tolerance = 5e-3)
  expect_equal(unname(S[1, 2] / sqrt(S[1, 1] * S[2, 2])), 0.8, tolerance = 5e-3)
})

test_that("the Taylor design is built ONCE and cached on the study", {
  # .adghGrid runs inside the objective. Rebuilding the design there re-ran an
  # 8192-row Sobol pass through the joint sampler on every evaluation (8.6 ms a
  # call, ~21 s across a fit) for a quantity that is DATA and cannot change.
  st <- list(s = list(cov = list(WT = 72, CRCL = 90),
                      cov_dist = admixr2:::.admCovDistCanon(
                        list(WT = list(mean = 72, sd = 16),
                             CRCL = list(mean = 90, sd = 25), cor = 0.6))))
  ui  <- list(allCovs = c("WT", "CRCL"),
              lstExpr = list(quote(cl <- exp(tcl + bw*WT + bc*CRCL + eta.cl))))
  pin <- list(n_eta = 1L, eta_col_names = "eta.cl",
              cov_integration = "taylor", cov_taylor_h = 0.5)
  out <- admixr2:::.admCheckCovariates(ui, pin, st, "sens")
  td  <- out$s$.adm_cov_taylor
  expect_false(is.null(td))
  expect_identical(td$n_pt, 5L)
  # it travels to a parallel-restart worker by value, so it must hold no
  # closures and no external pointers
  expect_true(all(!vapply(td, is.function, logical(1))))
  # and it must equal what the objective would have built for itself
  expect_equal(td, admixr2:::.admCovTaylorDesign(out$s$cov_dist, 0.5))
  # quadrature must NOT pay for a design it does not use
  pin$cov_integration <- "quadrature"
  expect_null(admixr2:::.admCheckCovariates(ui, pin, st, "sens")$s$.adm_cov_taylor)
})

# =============================================================================
# Covariate STRATA -- the per-covariate stratify/marginalise split
# =============================================================================

test_that("a `cor` spec conditions on exact POINTS, matching closed form", {
  # admixr2 built this copula, so the conditional is closed form on the latent
  # scale -- no pool, no importance weighting, nothing to collapse. It is also
  # the cheapest route at any useful stratum count.
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
  # A `joint` sampler with no density cannot be importance-weighted, so the
  # strata are equiprobable bins of the pool instead -- correct, but coarser.
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

test_that("a vine is conditioned correctly even when the stratified covariate
           is NOT first in its structure", {
  skip_if_not_installed("rvinecopulib")
  # This is the case the u-space route got silently wrong: on a vine whose
  # cascade runs AGE -> WT -> CRCL, stratifying on WT left AGE at its
  # unconditional 55.09 in every stratum, against a true 50.1 to 60.4.
  set.seed(11); n <- 3000L
  R <- matrix(c(1, .65, .30, .65, 1, .15, .30, .15, 1), 3, 3)
  U <- stats::pnorm(matrix(stats::rnorm(n * 3), n, 3) %*% chol(R))
  vc <- rvinecopulib::vinecop(U, family_set = "all")
  ml <- c(log(72), log(90), log(54)); sl <- c(0.22, 0.26, 0.20)
  cl <- function(x) pmin(pmax(x, 1e-12), 1 - 1e-12)
  jf <- function(u) {
    w <- rvinecopulib::inverse_rosenblatt(cl(u), vc)
    o <- vapply(1:3, function(j) stats::qlnorm(cl(w[, j]), ml[j], sl[j]),
                numeric(nrow(u)))
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
  # The model is the only thing that knows which covariates a source
  # conditioned on. Making the user restate it duplicates information already
  # present, and getting it wrong is the fabricated null contrast.
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
  # a source whose model reads NONE of the declared covariates has nothing to
  # stratify on, and saying so beats silently generating a null contrast
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
  # A model covariate a study never mentions is held at whatever rxSolve
  # defaults it to -- the ecological plug-in wearing a fit's clothes.
  ui <- .cov_ui(cov = c("WT", "CRCL"), expr = list(
    quote(cl <- exp(tcl + tcov * WT + tcr * CRCL + eta.cl))))
  chk <- function(st) admixr2:::.admCheckCovariates(ui, .cov_pinfo(), st, "none")
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
  expect_error(covDist(list(mean = 1, sd = 1)), "NAMED argument per covariate")
  # `dist` sets the default margin, and a per-covariate one still wins
  expect_equal(mean(covDraw(covDist(A = c(mean = 0, sd = 1), dist = "normal"),
                            n = 20000L)[, "A"]), 0, tolerance = 5e-2)
})

test_that("print.covDist reports what was DECLARED, not what it canonicalised to", {
  # `cor` becomes a joint sampler internally; printing "joint sampler" would
  # hide the number the user typed
  out <- utils::capture.output(print(covDist(WT = c(mean = 72, sd = 16),
                                             CRCL = c(mean = 90, sd = 25),
                                             cor = 0.6)))
  expect_true(any(grepl("cor = 0.6", out, fixed = TRUE)))
  expect_false(any(grepl("joint sampler", out, fixed = TRUE)))
  expect_true(any(grepl("normal", out)))       # the default margin
  out2 <- utils::capture.output(print(covDist(SEX = c(f = 0.6, m = 0.4))))
  expect_true(any(grepl("categorical", out2)))
  expect_true(any(grepl("f=0", out2, fixed = TRUE)))
})

# =============================================================================
# Shift path -- the pieces that need no compiled model
# =============================================================================

test_that(".admShiftNodes inverts the mixture CDF of u = Delta + eta", {
  # u's law is sum_j W_j N(Delta_j, om^2); the nodes must sit exactly at its
  # Gauss-Hermite probabilities, which is what makes n_u rows stand in for the
  # whole covariate x eta product grid.
  D <- c(-0.4, -0.1, 0.05, 0.3, 0.7); W <- c(.1, .2, .4, .2, .1); om <- 0.3
  for (n_u in c(7L, 15L, 31L)) {
    un <- admixr2:::.admShiftNodes(D, W, om, n_u)
    g  <- admixr2:::.adghNodes1(n_u)
    Fu <- vapply(un$u, function(t) sum(W * stats::pnorm((t - D) / om)), 0)
    expect_equal(Fu, stats::pnorm(g$x), tolerance = 1e-10)
    expect_equal(sum(un$w), 1)
    # a quantile map is monotone, so the nodes come back in the SAME order as
    # the probabilities that generated them (.adghNodes1 is not sorted
    # ascending). Only where the probability is RESOLVABLE: by 31 nodes the
    # outermost sits near 9.9 SD, where pnorm() saturates to exactly 1 and the
    # quantile is not determined. Those nodes carry weight ~1e-23 and so cannot
    # move a moment, but they are not ordered.
    keep <- abs(g$x) < 8
    expect_identical(order(un$u[keep]), order(g$x[keep]))
  }
  # a degenerate covariate (no spread) must give back the plain eta nodes
  un <- admixr2:::.admShiftNodes(rep(0, 3L), rep(1/3, 3L), om, 9L)
  expect_equal(un$u, om * admixr2:::.adghNodes1(9L)$x, tolerance = 1e-9)
})

test_that(".admShiftDu gives the chain factors as conditional expectations", {
  # du/dtheta = E[dDelta/dtheta | u] and du/dom = E[(u-Delta)/om | u], both
  # under the covariate quadrature that built F_u -- checked against a finite
  # difference of the NODES themselves.
  # spec carries one entry per AFFECTED random effect: eta, link and the
  # model's own right-hand side, so `rhs` is a list even at m = 1
  spec <- list(eta = "eta.cl", link = "exp",
               rhs = list(quote(exp(tcl + eta.cl) * (WT / 70)^b1)))
  st <- list(tcl = log(4), b1 = 0.75)
  X  <- matrix(c(55, 62, 70, 85, 100), 5L, 1L, dimnames = list(NULL, "WT"))
  W  <- c(.1, .2, .4, .2, .1); ar <- list(WT = 70); om <- 0.3; n_u <- 11L
  D  <- admixr2:::.admShiftDelta(spec, st, X, ar)[, 1L]
  u  <- admixr2:::.admShiftNodes(D, W, om, n_u)$u
  du <- admixr2:::.admShiftDu(spec, st, X, ar, D, W, om, u)

  h <- 1e-3
  fd <- function(k) {
    s1 <- st; s1[[k]] <- s1[[k]] + h; s2 <- st; s2[[k]] <- s2[[k]] - h
    (admixr2:::.admShiftNodes(admixr2:::.admShiftDelta(spec, s1, X, ar)[, 1L],
                              W, om, n_u)$u -
     admixr2:::.admShiftNodes(admixr2:::.admShiftDelta(spec, s2, X, ar)[, 1L],
                              W, om, n_u)$u) / (2 * h) }
  expect_equal(du$du_dtheta[, "b1"], fd("b1"), tolerance = 1e-5)
  # the TYPICAL VALUE cancels out of Delta (a difference of the same expression
  # at two covariate values), so it contributes nothing through u. Analytically
  # zero; what survives is the round-off of the finite difference that measures
  # dDelta/dtheta, which is ~1e-11 against b1's ~0.6.
  expect_lt(max(abs(du$du_dtheta[, "tcl"])), 1e-8)
  expect_gt(max(abs(du$du_dtheta[, "b1"])), 0.1)
  fo <- (admixr2:::.admShiftNodes(D, W, om + h, n_u)$u -
         admixr2:::.admShiftNodes(D, W, om - h, n_u)$u) / (2 * h)
  expect_equal(du$du_domega, fo, tolerance = 1e-5)
})

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
  # ... and it is ONE evaluation over the whole node set: a per-row loop that
  # lost the column names would leave every node at the reference and return
  # zeros, which is finite, plausible and wrong
  expect_false(all(D == 0))

  # m = 2: a covariate on two mu-referenced parameters gives two columns, and
  # the shift for each is that parameter's own coefficient
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

test_that("the Gaussian shift branch fires on exactly the qualifying laws", {
  # The condition is that Delta(a) is NORMAL, not that the covariate is. That
  # admits the allometric case for free: for lognormal WT, tcov*log(WT/70) is
  # affine in the latent normal score. Deliberately NOT read off the model text
  # or off rxode2's covariate frames -- on 5.1.4 those split the three common
  # spellings across muRefCovariateDataFrame, muRefExtra and
  # mu2RefCovariateReplaceDataFrame, so no single one of them sees this.
  gh <- admixr2:::.adghNodes1; g <- gh(31L)
  z  <- g$x; W <- g$w / sum(g$w)
  qual <- list(normal_linear = 0.8 * (70 + 6 * z),
               normal_scaled = 0.8 * (70 + 6 * z) / 70,
               lognormal_log = 0.8 * log(exp(log(70) + 0.17 * z) / 70),
               degenerate    = rep(0.5, length(z)))
  nope <- list(lognormal_linear = 0.8 * exp(log(70) + 0.17 * z),
               normal_square    = 0.8 * (70 + 6 * z)^2 / 100,
               normal_log       = 0.8 * log(70 + 6 * z))
  for (nm in names(qual))
    expect_lt(admixr2:::.admShiftGaussResid(qual[[nm]], W), 1e-8)
  for (nm in names(nope))
    expect_gt(admixr2:::.admShiftGaussResid(nope[[nm]], W), 1e-3)
  # and the node builder flags it
  expect_true(isTRUE(admixr2:::.admShiftNodes(qual$lognormal_log, W, 0.3, 20L)$gauss))
  expect_null(admixr2:::.admShiftNodes(nope$normal_square, W, 0.3, 20L)$gauss)
})

test_that("the Gaussian branch is exact where the mixture route only converges", {
  # u = Delta + eta is exactly normal here, so the closed form is the truth and
  # the quadrature over it is exact. The mixture route degrades as the covariate
  # outruns the random effect because n_u nodes resolve a widely separated
  # mixture worst -- that is the regime the branch exists for.
  gh <- admixr2:::.adghNodes1; g <- gh(31L); z <- g$x; W <- g$w / sum(g$w)
  om <- 0.30
  for (ratio in c(0.5, 2, 8)) {
    D  <- ratio * om * z                       # exactly normal, mean 0
    nd <- admixr2:::.admShiftNodes(D, W, om, 20L)
    expect_true(isTRUE(nd$gauss))
    gg  <- gh(20L); w <- gg$w / sum(gg$w)
    # E[exp(u)] = exp(om^2/2) * E_a[exp(Delta)] exactly, eta independent of a
    truth <- exp(om^2 / 2) * sum(W * exp(D))
    expect_equal(sum(w * exp(nd$u)), truth, tolerance = 1e-10)
  }
})

test_that("shift node derivatives match a difference where u is determined", {
  # Only where the CDF is not saturated: the outermost node's target sits
  # ~1e-14 from 1, u is undetermined across a plateau there, and its difference
  # quotient is noise rather than a reference. Interior nodes are the check.
  gh <- admixr2:::.adghNodes1; g <- gh(31L); z <- g$x; W <- g$w / sum(g$w)
  om <- 0.30
  D  <- 0.8 * (70 + 6 * z); D <- D - mean(D)   # Gaussian branch
  nd <- admixr2:::.admShiftNodes(D, W, om, 20L)
  mD <- sum(W * D); vD <- max(sum(W * D^2) - mD^2, 0); s <- sqrt(vD + om^2)
  an <- (om / s) * ((nd$u - mD) / s)           # closed form used by .admShiftDu
  h  <- 1e-6
  fd <- (admixr2:::.admShiftNodes(D, W, om + h, 20L)$u -
         admixr2:::.admShiftNodes(D, W, om - h, 20L)$u) / (2 * h)
  expect_equal(an, fd, tolerance = 1e-6)
})

test_that("shift nodes refuse non-finite input instead of erroring", {
  gh <- admixr2:::.adghNodes1; g <- gh(15L); W <- g$w / sum(g$w)
  D <- 0.5 * g$x; D[3] <- NaN
  # previously reached `if (max(abs(st)) < tol)` as a missing value
  expect_null(admixr2:::.admShiftNodes(D, W, 0.3, 20L))
  D2 <- 0.5 * g$x; D2[1] <- Inf
  expect_null(admixr2:::.admShiftNodes(D2, W, 0.3, 20L))
})

test_that("the affine certificate decides the Gaussian branch in 2-D", {
  # Delta = c + B z is exactly (jointly) normal because admixr2 builds every
  # continuous covariate as F^-1(Phi(z)) from a jointly normal z. A moment test
  # on Delta cannot certify this above one dimension -- Cramer-Wold needs ALL
  # projections, so finitely many can only fail to find a counterexample.
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
  # WITHOUT the latent scores a vector shift must not take the branch, even
  # though every margin here is normal: that is the case a moment test would
  # wave through without being able to justify it.
  expect_false(admixr2:::.admShiftGaussOK(aff$correlated, W, NULL, 2L))
})

test_that("the 2-D Gaussian branch is exact and skips the recursion", {
  gh <- admixr2:::.adghNodes1; nc <- 21L; g <- gh(nc)
  ix <- as.matrix(expand.grid(seq_len(nc), seq_len(nc)))
  z  <- cbind(g$x[ix[, 1]], g$x[ix[, 2]])
  W  <- g$w[ix[, 1]] * g$w[ix[, 2]]; W <- W / sum(W)
  A  <- matrix(c(0.40, 0.15, 0, 0.30), 2, 2)
  D  <- z %*% t(A); om <- c(0.30, 0.25)
  r  <- admixr2:::.admShiftNodesMulti(D, W, om, 12L, z = z)
  expect_true(isTRUE(r$gauss))
  expect_equal(nrow(r$u), 144L)
  expect_equal(sum(r$w), 1)
  # u ~ N(0, A A' + diag(om^2)); score a smooth integrand against the truth
  S <- A %*% t(A) + diag(om^2)
  cf <- c(1, 0.5)
  expect_equal(sum(r$w * exp(r$u %*% cf)),
               as.numeric(exp(0.5 * t(cf) %*% S %*% cf)), tolerance = 1e-10)
  # and the non-affine case still goes through the recursion and stays finite
  D2 <- cbind(0.4 * z[, 1], 0.3 * exp(0.2 * z[, 2]))
  r2 <- admixr2:::.admShiftNodesMulti(D2, W, om, 8L, z = z)
  expect_null(r2$gauss)
  expect_true(all(is.finite(r2$u)))
})

test_that("a refused inversion propagates out of the shift recursion", {
  # matrix(NULL, ncol = 1) is a 0-row matrix, so without propagation the caller
  # receives a well-shaped node set holding no nodes.
  gh <- admixr2:::.adghNodes1; nc <- 11L; g <- gh(nc)
  ix <- as.matrix(expand.grid(seq_len(nc), seq_len(nc)))
  z  <- cbind(g$x[ix[, 1]], g$x[ix[, 2]])
  W  <- g$w[ix[, 1]] * g$w[ix[, 2]]; W <- W / sum(W)
  D  <- cbind(0.4 * z[, 1], 0.3 * exp(0.2 * z[, 2]))
  D[5, 2] <- NaN                       # non-affine, so the recursion is used
  expect_null(admixr2:::.admShiftNodesMulti(D, W, c(0.3, 0.25), 8L, z = z))
})

test_that(".admCovGrid returns the latent scores behind X", {
  d <- list(WT = list(mu = 70, sd = 8), AGE = list(mu = 50, sd = 10))
  g <- admixr2:::.admCovGrid(d, 5L)
  expect_equal(nrow(g$z), nrow(g$X))
  expect_equal(ncol(g$z), 2L)
  # X is the image of z under F^-1(Phi(.)), so for a normal margin it is affine
  expect_lt(admixr2:::.admShiftAffineResid(g$X, g$W, g$z), 1e-10)
  # a discrete margin has no score, and one is enough to void the grid
  d2 <- list(WT = list(mu = 70, sd = 8), SEX = list(values = c(0, 1)))
  expect_null(admixr2:::.admCovGrid(d2, 5L)$z)
})
