# Tier 1 -- covariate marginalisation helpers (no rxode2).
#
# gl/gh/taylor node quadrature was removed; what remains is the marginal route
# (collapse, u-quantile, general per-row) plus the refusal that keeps old
# node-style study lists from being fitted as something else.

# ---- covariate collapse ------------------------------------------------------

.cov_pinfo <- function(n_eta = 1L, cmap = data.frame(covariate = "WT",
                                                     coef = "tcov",
                                                     eta = "eta.cl",
                                                     stringsAsFactors = FALSE))
  list(n_eta = n_eta,
       eta_col_names = if (n_eta == 1L) "eta.cl" else c("eta.cl", "eta.v"),
       cov_map = cmap)

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

test_that(".admCovInflateL folds theta_cov^2 * Sigma_a into Omega", {
  Om   <- matrix(0.09, 1L, 1L)
  pars <- list(L = t(chol(Om)), struct = c(tcl = 0, tcov = 0.75))
  s    <- list(cov_dist = list(WT = list(mu = 0, sd = 0.6)))
  L2   <- admixr2:::.admCovInflateL(pars, .cov_pinfo(), s)
  # 0.09 + 0.75^2 * 0.6^2 = 0.2925
  expect_equal(as.numeric(tcrossprod(L2)), 0.2925)
  # and it really is a Cholesky of that
  expect_equal(L2[lower.tri(L2)], numeric(0))
})

test_that(".admCovInflateL puts the variance on the RIGHT eta only", {
  Om   <- diag(c(0.09, 0.04))
  pars <- list(L = t(chol(Om)), struct = c(tcov = 0.5))
  cmap <- data.frame(covariate = "WT", coef = "tcov", eta = "eta.v",
                     stringsAsFactors = FALSE)
  s    <- list(cov_dist = list(WT = list(mu = 0, sd = 2)))
  L2   <- admixr2:::.admCovInflateL(pars, .cov_pinfo(2L, cmap), s)
  Om2  <- tcrossprod(L2)
  expect_equal(Om2[1, 1], 0.09)              # eta.cl untouched
  expect_equal(Om2[2, 2], 0.04 + 0.5^2 * 4)  # eta.v inflated
  expect_equal(Om2[1, 2], 0)
})

test_that(".admStudyL returns the plain Cholesky when no cov_dist is declared", {
  pars <- list(L = t(chol(matrix(0.09, 1L, 1L))), struct = c(tcov = 0.75))
  expect_identical(admixr2:::.admStudyL(pars, .cov_pinfo(), list()), pars$L)
})

test_that(".admCheckCovariates accepts a supported mu-referenced covariate", {
  st <- list(a = list(cov_dist = list(WT = list(mu = 0, sd = 0.6))))
  expect_silent(admixr2:::.admCheckCovariates(.cov_ui(), .cov_pinfo(), st, "none"))
})

test_that(".admCheckCovariates is a no-op when no study declares cov_dist", {
  expect_silent(admixr2:::.admCheckCovariates(.cov_ui(), .cov_pinfo(),
                                              list(a = list()), "sens"))
})

test_that("a gradient mode routes to the general path rather than erroring", {
  # Every estimator defaults to a gradient, so refusing here would make a
  # covariate model fail out of the box. Only "rows" carries a gradient.
  ok_st <- list(a = list(cov = list(WT = 0),
                         cov_dist = list(WT = list(mu = 0, sd = 0.6))))
  for (g in c("sens", "analytical", "fd", "cfd"))
    expect_identical(
      admixr2:::.admCheckCovariates(.cov_ui(), .cov_pinfo(), ok_st, g)$a$.adm_cov_path,
      "rows")
  # derivative-free still gets the exact closed form
  expect_identical(
    admixr2:::.admCheckCovariates(.cov_ui(), .cov_pinfo(), ok_st, "none")$a$.adm_cov_path,
    "collapse")
})

test_that(".admCheckCovariates ROUTES rather than refuses, most efficient first", {
  ok_st <- list(a = list(cov = list(WT = 0), cov_dist = list(WT = list(mu = 0, sd = 0.6))))
  path <- function(ui = .cov_ui(), pi = .cov_pinfo(), st = ok_st)
    admixr2:::.admCheckCovariates(ui, pi, st, "none")$a$.adm_cov_path

  # bare theta*COV + normal + single occurrence -> closed form
  expect_identical(path(), "collapse")
  # non-normal spec: the closed form does not apply. It used to fall to the
  # u-quantile path; that path is no longer routed to, because on a DISCRETE
  # covariate its unbracketed Newton solve silently returns garbage quantiles
  # (measured: 13-20% on the mean, 3-5x on the covariance). The general per-row
  # path assumes nothing and was exact on every configuration tested.
  expect_identical(path(st = list(a = list(cov = list(WT = 0),
      cov_dist = list(WT = list(values = c(0, 1), probs = c(.6, .4)))))), "rows")
  # covariate used in a SECOND place -> its whole effect no longer fits in one
  # eta column, so the general per-row path takes over
  expect_identical(path(ui = .cov_ui(expr = list(
      quote(cl <- exp(tcl + tcov * WT + eta.cl)),
      quote(v  <- exp(tv) * (1 + 0.01 * WT))))), "rows")
  # no eta to inflate or to carry a shift -> general path
  expect_identical(path(pi = .cov_pinfo(0L)), "rows")
  # no `cov` value given -> derived from the distribution, so the path is
  # unaffected (and the user cannot state a mean that contradicts cov_dist)
  expect_identical(path(st = list(a = list(
      cov_dist = list(WT = list(mu = 0, sd = 0.6))))), "collapse")
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

test_that(".admCovUQuantile solves F_u(u) = p exactly, and smoothly", {
  # Newton on the exact mixture CDF, not interpolation on a grid. The grid made
  # u piecewise-linear in the parameters, which is the mismatch that stops an
  # analytic gradient being consistent with the objective it differentiates.
  g  <- admixr2:::.adghNodes1(32L)
  dl <- list(delta = 0.75 * log(exp(log(72) + 0.28 * g$x) / 72), w = g$w)
  s  <- 0.30
  p  <- (seq_len(2000L) - 0.5) / 2000L
  u  <- admixr2:::.admCovUQuantile(dl, s, p)

  Fu <- as.numeric(pnorm(outer(u, dl$delta, "-") / s) %*% dl$w)
  expect_lt(max(abs(Fu - p)), 1e-10)          # solves its defining equation
  expect_true(all(diff(u) > 0))               # monotone
  expect_true(all(is.finite(u)))

  # SHIFTING every Delta by eps must shift u by exactly eps: an exact identity,
  # and one a grid interpolant cannot reproduce without quantisation.
  eps <- 1e-7
  u2  <- admixr2:::.admCovUQuantile(list(delta = dl$delta + eps, w = dl$w), s, p)
  expect_equal(mean((u2 - u) / eps), 1, tolerance = 1e-6)
  expect_lt(sd((u2 - u) / eps), 1e-6)

  # the implicit-function derivative wrt the eta SD matches finite differences,
  # which is what a future analytic gradient would rely on
  e2 <- 1e-6
  fd <- (admixr2:::.admCovUQuantile(dl, s + e2, p) -
         admixr2:::.admCovUQuantile(dl, s - e2, p)) / (2 * e2)
  z  <- outer(u, dl$delta, "-") / s
  fu <- as.numeric(dnorm(z) %*% dl$w) / s
  an <- -as.numeric((dnorm(z) * (-z / s)) %*% dl$w) / fu
  expect_lt(max(abs(an - fd) / pmax(abs(fd), 1e-8)), 1e-6)
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
  expect_match(msg, "its own `cov` and its own `n`")
})
