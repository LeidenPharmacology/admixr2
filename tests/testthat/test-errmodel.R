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
