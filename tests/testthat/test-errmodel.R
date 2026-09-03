# Residual error models: variance forms and analytical sigma gradients.
#
# admixr2 used to support exactly three residual models -- add, prop, lnorm --
# hard-coded as an if/else chain repeated at ~8 sites in R and 3 in C++. Anything
# else (combined1, pow, addPow) was warned about once and then treated as
# ADDITIVE, so those models fitted silently and wrongly. errmodel.R replaces that
# chain with a single per-endpoint spec plus closed-form derivatives.
#
# These tests pin down, for every supported form:
#   1. the variance it puts on diag(V), against the closed form; and
#   2. that the analytical d(var)/d(p) and d(mu)/d(p) match finite differences.
#
# (2) is the part that matters: the gradient is what the optimizer follows, and a
# wrong analytical gradient produces a plausible-looking but wrong fit.

.em_pinfo <- function(names, role, spec) {
  list(sigma_names = names, sigma_role = role, resid = list(cp = spec))
}

.em_spec <- function(form, k_add = NA_integer_, k_prop = NA_integer_,
                     k_pow = NA_integer_) {
  list(output = "cp", form = form, k_add = k_add, k_prop = k_prop, k_pow = k_pow)
}

# Residual variance / mean at optimizer vector p, for prediction f.
.em_apply <- function(pinfo, p, f) {
  nat <- admixr2:::.admSigmaNat(p, pinfo)
  arr <- admixr2:::.admResidRows(pinfo, "cp", nat, length(f))
  ap  <- admixr2:::.admResidApply(f, rep(0, length(f)), arr)
  list(mu = ap$mu, var = ap$dv, arr = arr)
}

# Analytical vs central-FD derivatives w.r.t. each optimizer parameter.
.em_expect_grad_matches_fd <- function(pinfo, p, f, tol = 1e-5) {
  nat <- admixr2:::.admSigmaNat(p, pinfo)
  arr <- admixr2:::.admResidRows(pinfo, "cp", nat, length(f))
  # var_f = 0 matches .em_apply(), which passes a zero structural variance --
  # these unit tests exercise the residual in isolation, not its composition
  # with Cov_eta(f) (that is test-integration-resid-moments.R).
  d   <- admixr2:::.admResidDeriv(f, rep(0, length(f)), arr, pinfo)

  h <- 1e-6
  for (k in seq_along(p)) {
    hi <- p; hi[k] <- p[k] + h
    lo <- p; lo[k] <- p[k] - h
    fd_var <- (.em_apply(pinfo, hi, f)$var - .em_apply(pinfo, lo, f)$var) / (2 * h)
    fd_mu  <- (.em_apply(pinfo, hi, f)$mu  - .em_apply(pinfo, lo, f)$mu)  / (2 * h)
    expect_equal(d$dvar[, k], fd_var, tolerance = tol,
                 info = paste("d(var)/dp for param", k))
    expect_equal(d$dmu[, k],  fd_mu,  tolerance = tol,
                 info = paste("d(mu)/dp for param", k))
  }

  # d(var)/d(f): the V-path every estimator uses to chain a struct theta through
  # the residual. A prop/pow variance depends on the prediction, so this is not zero.
  fd_df <- vapply(seq_along(f), function(t) {
    fh <- f; fh[t] <- f[t] + 1e-7
    fl <- f; fl[t] <- f[t] - 1e-7
    (admixr2:::.admResidApply(fh, rep(0, length(f)), arr)$dv[t] -
     admixr2:::.admResidApply(fl, rep(0, length(f)), arr)$dv[t]) / 2e-7
  }, numeric(1))
  expect_equal(d$dv_df, fd_df, tolerance = tol, info = "d(var)/d(f)")
}

f_test <- c(0.5, 2.0, 7.5)
a  <- 0.30; b <- 0.20; cpow <- 1.35
pa <- 2 * log(a)      # additive param on log-variance scale
pb <- 2 * log(b)      # proportional/power param on log-variance scale

# ---- existing forms (must be unchanged) --------------------------------------

test_that("add: var = a^2, gradient matches FD", {
  pinfo <- .em_pinfo("add.sd", "var", .em_spec(0L, k_add = 1L))
  expect_equal(.em_apply(pinfo, pa, f_test)$var, rep(a^2, 3))
  .em_expect_grad_matches_fd(pinfo, pa, f_test)
})

test_that("prop: var = (b*f)^2, gradient matches FD", {
  pinfo <- .em_pinfo("prop.sd", "var", .em_spec(0L, k_prop = 1L))
  expect_equal(.em_apply(pinfo, pb, f_test)$var, (b * f_test)^2)
  .em_expect_grad_matches_fd(pinfo, pb, f_test)
})

test_that("lnorm: mu = f*exp(s/2), var = mu^2*(exp(s)-1), gradient matches FD", {
  pinfo <- .em_pinfo("ln.sd", "var", .em_spec(2L, k_add = 1L))
  sv <- exp(pa)
  r  <- .em_apply(pinfo, pa, f_test)
  expect_equal(r$mu,  f_test * exp(sv / 2))
  expect_equal(r$var, (f_test * exp(sv / 2))^2 * (exp(sv) - 1))
  .em_expect_grad_matches_fd(pinfo, pa, f_test)
})

# ---- combined ----------------------------------------------------------------

test_that("combined2 (add + prop): var = a^2 + b^2*f^2, gradient matches FD", {
  pinfo <- .em_pinfo(c("add.sd", "prop.sd"), c("var", "var"),
                     .em_spec(0L, k_add = 1L, k_prop = 2L))
  expect_equal(.em_apply(pinfo, c(pa, pb), f_test)$var, a^2 + b^2 * f_test^2)
  .em_expect_grad_matches_fd(pinfo, c(pa, pb), f_test)
})

test_that("combined1 (add + prop): var = (a + b*f)^2, gradient matches FD", {
  # The form admixr2 previously could not represent at all: the independent
  # per-sigma addition structurally produces combined2, so an explicit
  # combined1() was silently computed as a^2 + b^2*f^2 -- missing the 2ab*f
  # cross term entirely.
  pinfo <- .em_pinfo(c("add.sd", "prop.sd"), c("var", "var"),
                     .em_spec(1L, k_add = 1L, k_prop = 2L))
  expect_equal(.em_apply(pinfo, c(pa, pb), f_test)$var, (a + b * f_test)^2)
  .em_expect_grad_matches_fd(pinfo, c(pa, pb), f_test)
})

test_that("combined1 and combined2 genuinely differ (cross term is present)", {
  p  <- c(pa, pb)
  c2 <- .em_pinfo(c("add.sd", "prop.sd"), c("var", "var"),
                  .em_spec(0L, k_add = 1L, k_prop = 2L))
  c1 <- .em_pinfo(c("add.sd", "prop.sd"), c("var", "var"),
                  .em_spec(1L, k_add = 1L, k_prop = 2L))
  v2 <- .em_apply(c2, p, f_test)$var
  v1 <- .em_apply(c1, p, f_test)$var
  expect_equal(v1 - v2, 2 * a * b * f_test)   # the cross term
  expect_true(all(v1 > v2))
})

# ---- power -------------------------------------------------------------------

test_that("pow: var = (b*f^c)^2, gradient matches FD (incl. the exponent)", {
  pinfo <- .em_pinfo(c("pow.sd", "pow.c"), c("var", "pow_exp"),
                     .em_spec(0L, k_prop = 1L, k_pow = 2L))
  p <- c(pb, cpow)
  expect_equal(.em_apply(pinfo, p, f_test)$var, (b * f_test^cpow)^2)
  .em_expect_grad_matches_fd(pinfo, p, f_test)
})

test_that("addPow combined2: var = a^2 + b^2*f^(2c), gradient matches FD", {
  pinfo <- .em_pinfo(c("add.sd", "pow.sd", "pow.c"), c("var", "var", "pow_exp"),
                     .em_spec(0L, k_add = 1L, k_prop = 2L, k_pow = 3L))
  p <- c(pa, pb, cpow)
  expect_equal(.em_apply(pinfo, p, f_test)$var, a^2 + b^2 * f_test^(2 * cpow))
  .em_expect_grad_matches_fd(pinfo, p, f_test)
})

test_that("addPow combined1: var = (a + b*f^c)^2, gradient matches FD", {
  pinfo <- .em_pinfo(c("add.sd", "pow.sd", "pow.c"), c("var", "var", "pow_exp"),
                     .em_spec(1L, k_add = 1L, k_prop = 2L, k_pow = 3L))
  p <- c(pa, pb, cpow)
  expect_equal(.em_apply(pinfo, p, f_test)$var, (a + b * f_test^cpow)^2)
  .em_expect_grad_matches_fd(pinfo, p, f_test)
})

test_that("pow with exponent 1 reduces exactly to prop", {
  pw <- .em_pinfo(c("pow.sd", "pow.c"), c("var", "pow_exp"),
                  .em_spec(0L, k_prop = 1L, k_pow = 2L))
  pr <- .em_pinfo("prop.sd", "var", .em_spec(0L, k_prop = 1L))
  expect_identical(.em_apply(pw, c(pb, 1.0), f_test)$var,
                   .em_apply(pr, pb,         f_test)$var)
})

# ---- the pow exponent is not a variance --------------------------------------

test_that("pow exponent uses an identity transform, not log-variance", {
  pinfo <- .em_pinfo(c("pow.sd", "pow.c"), c("var", "pow_exp"),
                     .em_spec(0L, k_prop = 1L, k_pow = 2L))
  nat <- admixr2:::.admSigmaNat(c(pb, cpow), pinfo)
  expect_equal(unname(nat[1]), b^2)     # variance param: exp(p)
  expect_equal(unname(nat[2]), cpow)    # exponent: p itself, NOT exp(p)
})

test_that("a negative pow exponent is representable (exponents are unconstrained)", {
  # A log-variance encoding cannot represent this at all: exp(p) > 0 always.
  pinfo <- .em_pinfo(c("pow.sd", "pow.c"), c("var", "pow_exp"),
                     .em_spec(0L, k_prop = 1L, k_pow = 2L))
  p <- c(pb, -0.4)
  expect_equal(.em_apply(pinfo, p, f_test)$var, (b * f_test^(-0.4))^2)
  .em_expect_grad_matches_fd(pinfo, p, f_test)
})

# ---- the zero prediction -----------------------------------------------------
#
# f == 0 is not exotic: a depot model has it at t = 0, and it is where the moment
# expansion E[f^k] ~ mu^k + k(k-1)/2 * mu^(k-2) * Var(f) has a POLE for every
# k < 2 -- i.e. for every pow() exponent below 1, and for combined1's f^c term.
# Both the R path (.admMomF) and the C++ kernel (adm_mom_f) cap that correction
# term against the leading one so the result stays finite; these tests pin the
# behaviour down at f == 0 exactly, which no closed-form test above reaches.

test_that("pow with c < 1 stays finite at a zero prediction", {
  pinfo <- .em_pinfo(c("pow.sd", "pow.c"), c("var", "pow_exp"),
                     .em_spec(0L, k_prop = 1L, k_pow = 2L))
  f0 <- c(0, 0.5, 2)
  for (cc in c(0.2, 0.4, 0.5, 0.9)) {
    r <- .em_apply(pinfo, c(pb, cc), f0)
    expect_true(all(is.finite(r$var)), info = paste("c =", cc))
    expect_true(all(r$var >= 0),       info = paste("c =", cc))
    # At f = 0 with c > 0 the power variance is exactly 0 -- no residual scatter
    # where there is no signal. (Degenerate as a MODEL; the point is that it is
    # 0 and not NaN.)
    expect_equal(r$var[1L], 0, info = paste("c =", cc))
    # And the closed form still holds away from zero.
    expect_equal(r$var[-1L], (b * f0[-1L]^cc)^2, info = paste("c =", cc))
  }
})

test_that("derivatives at a zero prediction are finite for pow and combined1", {
  # A NaN here would poison the whole gradient vector, not just this row, so the
  # optimizer would fail at iteration 0 with no indication of which observation
  # caused it. Checked directly rather than through an FD, because a central
  # difference straddles f = 0 into f < 0 where f^c is undefined.
  f0 <- c(0, 1.5)
  specs <- list(
    pow       = list(.em_pinfo(c("pow.sd", "pow.c"), c("var", "pow_exp"),
                               .em_spec(0L, k_prop = 1L, k_pow = 2L)), c(pb, 0.4)),
    addPow2   = list(.em_pinfo(c("a.sd", "b.sd", "b.c"), c("var", "var", "pow_exp"),
                               .em_spec(0L, k_add = 1L, k_prop = 2L, k_pow = 3L)),
                     c(pa, pb, 0.4)),
    addPow1   = list(.em_pinfo(c("a.sd", "b.sd", "b.c"), c("var", "var", "pow_exp"),
                               .em_spec(1L, k_add = 1L, k_prop = 2L, k_pow = 3L)),
                     c(pa, pb, 0.4)),
    combined1 = list(.em_pinfo(c("a.sd", "b.sd"), c("var", "var"),
                               .em_spec(1L, k_add = 1L, k_prop = 2L)), c(pa, pb))
  )
  for (nm in names(specs)) {
    pinfo <- specs[[nm]][[1L]]; p <- specs[[nm]][[2L]]
    nat <- admixr2:::.admSigmaNat(p, pinfo)
    arr <- admixr2:::.admResidRows(pinfo, "cp", nat, length(f0))
    # A NON-zero structural variance is the case that actually exercises the cap:
    # with var_f = 0 the correction term is 0 whatever the pole does.
    d   <- admixr2:::.admResidDeriv(f0, c(0.25, 0.25), arr, pinfo)
    for (fld in c("dvar", "dmu", "dv_df", "dv_dv0"))
      expect_true(all(is.finite(d[[fld]])), info = paste(nm, fld))
    ap <- admixr2:::.admResidApply(f0, c(0.25, 0.25), arr)
    expect_true(all(is.finite(ap$dv)) && all(is.finite(ap$mu)), info = nm)
    expect_true(all(ap$dv >= 0), info = nm)
  }
})

test_that("the C++ moment helper agrees with the R one at a zero prediction", {
  # adm_mom_f (src/nll.cpp) and .admMomF (R/errmodel.R) are two implementations of
  # the same expansion, and the fused MC kernels use the C++ one while the NLL
  # composition uses the R one. If their poles are capped differently, admc's NLL
  # and its gradient describe different functions at f = 0.
  cp <- matrix(c(0, 0, 0.5, 2.0), 2L, 2L)   # column 1 is identically zero
  E  <- c(0.02, 1.2); V <- diag(c(0.05, 0.3))
  for (cc in c(0.4, 0.9, 1.0, 1.6)) {
    arr <- list(form = rep(0L, 2L), a2 = rep(0.09, 2L),
                b2 = rep(0.04, 2L), cc = rep(cc, 2L))
    nll <- admixr2:::nll_cov_from_samples_cpp(cp, E, V, 100L,
                                              arr$form, arr$a2, arr$b2, arr$cc)
    expect_true(is.finite(nll), info = paste("c =", cc))
  }
})

test_that(".admMomFd differentiates exactly what .admMomF evaluates", {
  # The two carry one expansion between them: .admMomF is what the NLL scores and
  # .admMomFd is what the gradient chains through. They used to cap the divergent
  # mu^(k-2) term by DIFFERENT rules -- .admMomF against the leading term, .admMomFd
  # by zeroing it past .ADM_MOM_CAP -- so for pow(b, c) with c < 1 near a zero
  # prediction the optimizer was handed the gradient of a different function.
  F <- admixr2:::.admMomF; Fd <- admixr2:::.admMomFd
  binds <- function(mu, v0, k) {
    lead <- mu^k; corr <- (k * (k - 1) / 2) * mu^(k - 2) * v0
    isTRUE((k - 2) < 0 && is.finite(lead) && is.finite(corr) &&
             abs(corr) > abs(lead))
  }
  grid <- expand.grid(mu = c(1e-6, 1e-3, 0.1, 1, 5, 40),
                      v0 = c(0, 1e-6, 1e-2, 1, 9),
                      k  = c(0.5, 1.5, 2, 2.4, 3))
  # the capped branch is the one that regressed, so make sure it is reached
  expect_gt(sum(mapply(binds, grid$mu, grid$v0, grid$k)), 10L)

  # A finite difference is a valid reference only when both probes sit on the same
  # side of the cap (the capped expression is genuinely discontinuous there) and
  # the difference clears the cancellation floor of the two values it subtracts.
  ok <- function(fp, fm) abs(fp - fm) > 1e-9 * max(abs(fp), abs(fm), 1e-300)
  for (i in seq_len(nrow(grid))) {
    mu <- grid$mu[i]; v0 <- grid$v0[i]; k <- grid$k[i]
    a  <- Fd(mu, v0, k)
    lbl <- sprintf("mu=%g v0=%g k=%g", mu, v0, k)
    expect_equal(a$m, F(mu, v0, k), tolerance = 1e-12, info = lbl)

    h <- 1e-6 * mu
    fp <- F(mu + h, v0, k); fm <- F(mu - h, v0, k)
    if (binds(mu + h, v0, k) == binds(mu - h, v0, k) && ok(fp, fm))
      expect_equal(a$dmu, (fp - fm) / (2 * h), tolerance = 1e-4, info = lbl)

    if (v0 > 0) {
      h <- 1e-6 * v0
      fp <- F(mu, v0 + h, k); fm <- F(mu, v0 - h, k)
      if (binds(mu, v0 + h, k) == binds(mu, v0 - h, k) && ok(fp, fm))
        expect_equal(a$dv0, (fp - fm) / (2 * h), tolerance = 1e-4, info = lbl)
    }

    h <- 1e-6
    fp <- F(mu, v0, k + h); fm <- F(mu, v0, k - h)
    if (binds(mu, v0, k + h) == binds(mu, v0, k - h) && ok(fp, fm))
      expect_equal(a$dk, (fp - fm) / (2 * h), tolerance = 1e-4, info = lbl)
  }
})

# ---- legacy fallback ---------------------------------------------------------

test_that("a pinfo with only the legacy flags still parses (no $resid spec)", {
  # Hand-built pinfo (unit tests, Tier-1 mock iniDf) carries sigma_is_prop /
  # sigma_is_lnorm but no per-endpoint spec. That subset is exactly combined2
  # with exponent 1, and must keep producing bit-identical arithmetic.
  legacy <- list(sigma_names    = c("add.sd", "prop.sd"),
                 sigma_is_prop  = list(FALSE, TRUE),
                 sigma_is_lnorm = list(FALSE, FALSE),
                 sigma_output   = c(NA_character_, NA_character_))
  nat <- admixr2:::.admSigmaNat(c(pa, pb), legacy)
  arr <- admixr2:::.admResidRows(legacy, NULL, nat, length(f_test))
  v   <- admixr2:::.admResidApply(f_test, rep(0, length(f_test)), arr)$dv

  old <- rep(0, length(f_test))
  old <- old + exp(pa)                  # additive, then
  old <- old + exp(pb) * f_test^2       # proportional -- the old per-sigma loop
  expect_identical(v, old)
})

# ---- the shared moment tail --------------------------------------------------
#
# .admResidMoments() / .admResidSampleMoments() / .admResidChain() collapse the
# four-step tail (.admResidApply -> .admApplyResidTail -> .admResidDeriv ->
# .admResidVChain) every estimator used to write out by hand. Each test below
# pins one of the traps that made the hand-assembly worth removing.

.mt_pinfo <- function(ar_fixed = NA_real_) {
  sp <- list(output = "cp", form = admixr2:::.ADM_RESID_COMBINED2,
             k_add = 1L, k_prop = 2L, k_pow = NA_integer_, ar_fixed = ar_fixed)
  list(sigma_names = c("add.sd", "prop.sd"), sigma_role = c("var", "var"),
       resid = list(cp = sp))
}

.mt_setup <- function(ar_fixed = NA_real_) {
  pinfo <- .mt_pinfo(ar_fixed)
  p     <- c(2 * log(0.4), 2 * log(0.15))
  nat   <- admixr2:::.admSigmaNat(p, pinfo)
  mu    <- c(4, 12, 30)
  times <- c(1, 2, 4)
  arr   <- admixr2:::.admResidRows(pinfo, "cp", nat, length(mu))
  cov_f <- matrix(c(1.5, 0.4, 0.2, 0.4, 6.0, 0.9, 0.2, 0.9, 20.0), 3, 3)
  list(pinfo = pinfo, arr = arr, mu = mu, times = times, cov_f = cov_f,
       var_f = diag(cov_f))
}

test_that(".admResidMoments reproduces the hand-written objective tail", {
  st <- .mt_setup()
  # cov branch: apply with times + structural covariance, then compose
  ap  <- admixr2:::.admResidApply(st$mu, st$var_f, st$arr, st$times, st$cov_f)
  ref <- admixr2:::.admApplyResidTail(st$cov_f, ap)
  got <- admixr2:::.admResidMoments(st$mu, st$var_f, st$arr, st$cov_f, st$times)
  expect_identical(got$mu, ap$mu)
  expect_identical(got$dv, ap$dv)
  expect_identical(got$V,  ref)

  # var branch: no structural covariance -> no times, no composed V
  apv <- admixr2:::.admResidApply(st$mu, st$var_f, st$arr)
  gv  <- admixr2:::.admResidMoments(st$mu, st$var_f, st$arr)
  expect_identical(gv$mu, apv$mu)
  expect_identical(gv$dv, apv$dv)
  expect_null(gv$V)
})

# THE trap CLAUDE.md records: arr$rho (and the ordinal cross term) key off
# `times`, so a diagonal-path caller that forwards `times` changes its own NLL.
# .admResidMoments() only forwards `times` alongside `cov_f`, so the diagonal
# path cannot turn the off-diagonal terms on even if a caller passes times.
test_that(".admResidMoments cannot leak `times` into the diagonal path", {
  st <- .mt_setup(ar_fixed = 0.5)
  expect_false(is.na(st$arr$rho[[1L]]))          # premise: ar() really is active

  # with a structural covariance, ar() DOES reach the off-diagonal
  cv <- admixr2:::.admResidMoments(st$mu, st$var_f, st$arr, st$cov_f, st$times)
  expect_false(isTRUE(all.equal(cv$V[1L, 2L], st$cov_f[1L, 2L])))

  # without one, passing times must change nothing at all
  a <- admixr2:::.admResidMoments(st$mu, st$var_f, st$arr)
  b <- admixr2:::.admResidMoments(st$mu, st$var_f, st$arr, NULL, st$times)
  expect_identical(a$mu, b$mu)
  expect_identical(a$dv, b$dv)
  expect_null(b$V)
  expect_null(b$rmat)
})

test_that(".admResidSampleMoments reproduces the MC sample-moment block", {
  st <- .mt_setup()
  set.seed(7L)
  cp <- matrix(stats::rnorm(200L * 3L, rep(st$mu, each = 200L), 1.2), 200L, 3L)
  mu_s <- colMeans(cp)
  cpc  <- sweep(cp, 2L, mu_s)
  Vs   <- crossprod(cpc) / nrow(cp)
  ap   <- admixr2:::.admResidApply(mu_s, diag(Vs), st$arr, st$times, Vs)
  got  <- admixr2:::.admResidSampleMoments(cp, st$arr, st$times)
  expect_identical(got$mu, ap$mu)
  expect_identical(got$dv, ap$dv)
  expect_identical(got$V,  admixr2:::.admApplyResidTail(Vs, ap))
  expect_identical(got$cov_f, Vs)
})

test_that(".admResidChain reproduces the hand-written gradient tail", {
  st <- .mt_setup()
  n_t     <- length(st$mu)
  dNLL_dV <- matrix(0.15, n_t, n_t) + diag(0.6, n_t)
  dmu     <- c(0.3, -0.2, 0.11)
  dvdiag  <- diag(dNLL_dV)

  d   <- admixr2:::.admResidDeriv(st$mu, st$var_f, st$arr, st$pinfo)
  M   <- admixr2:::.admResidVChain(st$mu, st$var_f, st$arr, st$pinfo, st$times,
                                   deriv = d)
  dmv <- attr(M, "dmu_dv0") %||% numeric(n_t)
  ref_dV <- dNLL_dV * M
  diag(ref_dV) <- diag(ref_dV) + dmu * dmv

  ch <- admixr2:::.admResidChain(st$mu, st$var_f, st$arr, st$pinfo, dmu, dvdiag,
                                 dNLL_dV, st$cov_f, st$times)
  expect_identical(ch$vchain,  M)
  expect_identical(ch$dmu_dv0, dmv)
  expect_identical(ch$dmu_df,  d$dmu_df)
  expect_identical(ch$dV_diag, dvdiag * diag(M) + dmu * dmv)
  expect_identical(ch$dV,      ref_dV)

  # the two contractions, which take cov_f/times in OPPOSITE orders
  expect_identical(
    ch$sigma_grad(),
    admixr2:::.admSigmaGrad(st$mu, st$arr, st$pinfo, dvdiag, dmu, st$var_f,
                            dNLL_dV, st$times, st$cov_f, deriv = d))
  expect_identical(
    ch$mu_coupling(),
    admixr2:::.admResidMuCoupling(st$mu, st$arr, st$pinfo, dvdiag, dmu,
                                  st$var_f, dNLL_dV, st$cov_f, st$times,
                                  deriv = d))
})

test_that(".admResidChain var branch drops dV and the off-diagonal terms", {
  st <- .mt_setup()
  dmu    <- c(0.3, -0.2, 0.11)
  dvdiag <- c(0.6, 0.5, 0.4)
  ch <- admixr2:::.admResidChain(st$mu, st$var_f, st$arr, st$pinfo, dmu, dvdiag,
                                 NULL, NULL, st$times)
  expect_null(ch$dV)
  d <- admixr2:::.admResidDeriv(st$mu, st$var_f, st$arr, st$pinfo)
  expect_identical(
    ch$sigma_grad(),
    admixr2:::.admSigmaGrad(st$mu, st$arr, st$pinfo, dvdiag, dmu, st$var_f,
                            NULL, st$times, NULL, deriv = d))
})
