# FO NLL and gradient correctness — Tier 1 (no rxode2).
#
# Uses an analytic 1-cmt bolus model to supply mu_pred and J without rxSolve,
# then validates:
#   - .adfoVpred() math for additive / proportional / lnorm sigma
#   - FO NLL formula agrees with nll_cov_cpp / nll_var_cpp
#   - Omega Cholesky gradient (cov + var branches) vs numerical FD of NLL
#   - Sigma additive/proportional gradient vs numerical FD
#   - Struct theta mu-path formula (sum(dNLL_dmu * J[:,j])) vs mu-only FD
#     (J held fixed).  .adfoGrad() uses full FD for struct thetas because J
#     also depends on theta (V-path contribution), but the mu-path formula
#     is internally consistent and tested here as a formula check.
#
# All checks are against numerical FD of the same closed-form NLL expression
# used in .adfoNLL(), so errors here indicate formula bugs, not solver issues.

# ---- Analytic 1-cmt helpers (bolus, eta=0) -----------------------------------

.fo1_mu <- function(tcl, tv, dose = 100, times = c(0.5, 1, 2, 4)) {
  cl <- exp(tcl); v <- exp(tv)
  (dose / v) * exp(-(cl / v) * times)
}

# J[t,j]: d(cp)/d(eta_j)|_{eta=0}.  col1 = eta_cl, col2 = eta_v.
.fo1_J <- function(tcl, tv, dose = 100, times = c(0.5, 1, 2, 4)) {
  cl <- exp(tcl); v <- exp(tv)
  mu    <- (dose / v) * exp(-(cl / v) * times)
  J_cl  <- -mu * times * cl / v          # d/d(eta_cl) at eta=0
  J_v   <-  mu * (times * cl / v - 1)   # d/d(eta_v)  at eta=0
  cbind(J_cl, J_v)
}

# FO NLL (cov branch) for given Cholesky L and additive sigma_var.
.fo_nll_cov <- function(L, mu_pred, J, sv_add, E_obs, V_obs, n) {
  omega  <- L %*% t(L)
  V_pred <- J %*% omega %*% t(J)
  diag(V_pred) <- diag(V_pred) + sv_add
  admixr2:::nll_cov_cpp(E_obs, V_obs, mu_pred, V_pred, n)
}

# FO NLL (var branch) for given Cholesky L and additive sigma_var.
.fo_nll_var <- function(L, mu_pred, J, sv_add, E_obs, v_obs, n) {
  omega  <- L %*% t(L)
  V_pred <- J %*% omega %*% t(J)
  v_pred <- diag(V_pred) + sv_add
  admixr2:::nll_var_cpp(E_obs, v_obs, mu_pred, v_pred, n)
}

# FO NLL (cov, proportional sigma).
.fo_nll_cov_prop <- function(L, mu_pred, J, sv_prop, E_obs, V_obs, n) {
  omega  <- L %*% t(L)
  V_pred <- J %*% omega %*% t(J)
  diag(V_pred) <- diag(V_pred) + sv_prop * mu_pred^2
  admixr2:::nll_cov_cpp(E_obs, V_obs, mu_pred, V_pred, n)
}

# FO NLL (var, proportional sigma).
.fo_nll_var_prop <- function(L, mu_pred, J, sv_prop, E_obs, v_obs, n) {
  omega  <- L %*% t(L)
  V_pred <- J %*% omega %*% t(J)
  v_pred <- diag(V_pred) + sv_prop * mu_pred^2
  admixr2:::nll_var_cpp(E_obs, v_obs, mu_pred, v_pred, n)
}

# Perturb a diagonal L entry (p = 2*log(L_ii)) by delta in optimizer scale.
.perturb_L_diag <- function(L, i, delta) {
  L_new      <- L
  L_new[i,i] <- exp((2 * log(L[i,i]) + delta) / 2)
  L_new
}

# Perturb an off-diagonal L entry (p = L_ij) directly.
.perturb_L_off <- function(L, i, j, delta) {
  L_new      <- L
  L_new[i,j] <- L_new[i,j] + delta
  L_new
}

# Analytical gradient of FO cov-NLL w.r.t. L and sigma (additive).
# Returns list with omega and sigma components.
.fo_grad_analytic_cov <- function(L, mu_pred, J, sv_add, E_obs, V_obs, n) {
  omega  <- L %*% t(L)
  V_pred <- J %*% omega %*% t(J)
  diag(V_pred) <- diag(V_pred) + sv_add

  invV      <- chol2inv(chol(V_pred))
  r         <- E_obs - mu_pred
  V_obs_hat <- V_obs + tcrossprod(r)
  dNLL_dV   <- n * (invV - invV %*% V_obs_hat %*% invV)

  M  <- t(J) %*% dNLL_dV %*% J
  ML <- M %*% L

  list(
    grad_L11   = as.numeric(ML[1,1]) * L[1,1],   # diagonal: p = 2*log(L_ii)
    grad_L22   = as.numeric(ML[2,2]) * L[2,2],
    grad_L21   = 2 * as.numeric(ML[2,1]),         # off-diagonal: p = L_ij directly
    grad_sigma = as.numeric(sum(diag(dNLL_dV)) * sv_add)
  )
}

# Analytical gradient of FO var-NLL w.r.t. L and sigma (additive).
.fo_grad_analytic_var <- function(L, mu_pred, J, sv_add, E_obs, v_obs, n) {
  omega  <- L %*% t(L)
  V_pred <- J %*% omega %*% t(J)
  v_pred <- diag(V_pred) + sv_add

  r            <- E_obs - mu_pred
  dNLL_dv_pred <- n * (1/v_pred - (v_obs + r^2) / v_pred^2)

  M  <- t(J) %*% (J * dNLL_dv_pred)   # n_eta x n_eta
  ML <- M %*% L

  list(
    grad_L11   = as.numeric(ML[1,1]) * L[1,1],
    grad_L22   = as.numeric(ML[2,2]) * L[2,2],
    grad_L21   = 2 * as.numeric(ML[2,1]),
    grad_sigma = as.numeric(sum(dNLL_dv_pred) * sv_add)
  )
}

# ---- Setup: known analytic parameters ----------------------------------------

.fo_setup <- function() {
  tcl   <- log(5); tv <- log(20)
  times <- c(0.5, 1, 2, 4)
  n     <- 80L;   dose <- 100
  n_t   <- length(times)
  mu    <- .fo1_mu(tcl, tv, dose, times)
  J     <- .fo1_J(tcl, tv, dose, times)
  L     <- matrix(c(0.3, 0.1, 0, 0.2), 2, 2)  # lower triangular
  sv_add  <- 0.04
  sv_prop <- 0.05
  E_obs   <- mu * 1.05    # slightly off truth
  V_obs   <- diag(0.01, n_t) + 0.002
  diag(V_obs) <- 0.015
  v_obs   <- diag(V_obs)
  list(tcl=tcl, tv=tv, times=times, n=n, dose=dose, n_t=n_t,
       mu=mu, J=J, L=L, sv_add=sv_add, sv_prop=sv_prop,
       E_obs=E_obs, V_obs=V_obs, v_obs=v_obs)
}

# ==============================================================================
# .adfoVpred() tests
#
# .adfoVpred no longer classifies the error model itself: it takes the row array
# built by .admResidRows(). Building that array here from the legacy
# sigma_is_prop / sigma_is_lnorm flags also exercises the fallback path used when
# a pinfo carries no per-endpoint `resid` spec (hand-built pinfo / mock iniDf).
# ==============================================================================

.fo_arr <- function(sigma_var, is_prop, is_lnorm, n_t) {
  pinfo <- list(sigma_names    = paste0("sig", seq_along(sigma_var)),
                sigma_is_prop  = is_prop,
                sigma_is_lnorm = is_lnorm,
                sigma_output   = rep(NA_character_, length(sigma_var)))
  admixr2:::.admResidRows(pinfo, NULL, sigma_var, n_t)
}

test_that(".adfoVpred(): additive sigma -- V = J*Omega*J^T + sv*I", {
  s   <- .fo_setup()
  arr <- .fo_arr(s$sv_add, list(FALSE), list(FALSE), s$n_t)
  vp  <- admixr2:::.adfoVpred(s$mu, s$J, s$L, arr, n_t = s$n_t, n_eta = 2L)
  V_expected <- s$J %*% (s$L %*% t(s$L)) %*% t(s$J)
  diag(V_expected) <- diag(V_expected) + s$sv_add

  expect_equal(vp$V, V_expected, tolerance = 1e-12)
  expect_equal(vp$mu_sigma, s$mu)          # no lnorm correction
})

test_that(".adfoVpred(): proportional sigma -- V diag += sv*E[f^2]", {
  s   <- .fo_setup()
  arr <- .fo_arr(s$sv_prop, list(TRUE), list(FALSE), s$n_t)
  vp  <- admixr2:::.adfoVpred(s$mu, s$J, s$L, arr, n_t = s$n_t, n_eta = 2L)
  V_expected <- s$J %*% (s$L %*% t(s$L)) %*% t(s$J)
  # Law of total variance: E_eta[Var(y|eta)] = sv*E[f^2] = sv*(Var_eta(f) + mu^2),
  # NOT sv*mu^2. The extra sv*Var_eta(f) is the eta-eps interaction term that the
  # old "evaluate the residual at the population mean" formula dropped.
  var_f <- diag(V_expected)
  diag(V_expected) <- var_f + s$sv_prop * (var_f + s$mu^2)

  expect_equal(vp$V, V_expected, tolerance = 1e-12)
  expect_equal(vp$mu_sigma, s$mu)
})

test_that(".adfoVpred(): lnorm sigma -- mean scaled, WHOLE covariance scaled", {
  s   <- .fo_setup()
  arr <- .fo_arr(s$sv_add, list(FALSE), list(TRUE), s$n_t)
  vp  <- admixr2:::.adfoVpred(s$mu, s$J, s$L, arr, n_t = s$n_t, n_eta = 2L)
  es  <- exp(s$sv_add)
  mu_sigma_exp <- s$mu * exp(s$sv_add / 2)
  V_struct <- s$J %*% (s$L %*% t(s$L)) %*% t(s$J)
  # E[y|eta] = f*exp(s/2), so Var_eta(E[y|eta]) scales the ENTIRE covariance by
  # exp(s) -- off-diagonals included -- and E_eta[Var(y|eta)] uses E[f^2], not mu^2.
  V_expected <- V_struct * es
  diag(V_expected) <- diag(V_struct) * es +
    (diag(V_struct) + s$mu^2) * es * (es - 1)

  expect_equal(vp$mu_sigma, mu_sigma_exp, tolerance = 1e-12)
  expect_equal(vp$V, V_expected, tolerance = 1e-12)
  # the two algebraic forms of the diagonal agree
  expect_equal(diag(vp$V), diag(V_struct) * es^2 + s$mu^2 * es * (es - 1),
               tolerance = 1e-12)
})

test_that(".adfoVpred(): no etas -- V is pure sigma contribution, J is ignored", {
  s   <- .fo_setup()
  J_empty  <- matrix(0, s$n_t, 0)
  omega_e  <- matrix(numeric(0), 0, 0)
  arr <- .fo_arr(s$sv_add, list(FALSE), list(FALSE), s$n_t)
  vp  <- admixr2:::.adfoVpred(s$mu, J_empty, omega_e, arr, n_t = s$n_t, n_eta = 0L)
  V_expected <- matrix(0, s$n_t, s$n_t)
  diag(V_expected) <- s$sv_add

  expect_equal(vp$V, V_expected, tolerance = 1e-12)
  expect_equal(vp$mu_sigma, s$mu)
})

test_that(".adfoVpred(): two sigma terms -- contributions additive (combined2)", {
  s   <- .fo_setup()
  # Two sigmas on one endpoint: one additive, one proportional. The legacy flags
  # describe combined2, so E[Var(y|eta)] = a^2 + b^2*E[f^2]. Only the ADDITIVE
  # term is free of the prediction; the proportional one picks up Var_eta(f) via
  # E[f^2] = Var_eta(f) + mu^2 (the eta-eps interaction).
  arr <- .fo_arr(c(s$sv_add, s$sv_prop), list(FALSE, TRUE), list(FALSE, FALSE), s$n_t)
  vp  <- admixr2:::.adfoVpred(s$mu, s$J, s$L, arr, n_t = s$n_t, n_eta = 2L)
  V_exp <- s$J %*% (s$L %*% t(s$L)) %*% t(s$J)
  var_f <- diag(V_exp)
  diag(V_exp) <- var_f + s$sv_add + s$sv_prop * (var_f + s$mu^2)

  expect_equal(vp$V, V_exp, tolerance = 1e-12)
})

# ==============================================================================
# FO NLL formula
# ==============================================================================

test_that("FO NLL (cov) is finite at non-degenerate setup", {
  # -2LL for continuous distributions can be negative (log-det term dominates
  # when V_pred is small), so we only check finiteness.
  s   <- .fo_setup()
  nll <- .fo_nll_cov(s$L, s$mu, s$J, s$sv_add, s$E_obs, s$V_obs, s$n)
  expect_true(is.finite(nll))
})

test_that("FO NLL (var) is finite at non-degenerate setup", {
  s   <- .fo_setup()
  nll <- .fo_nll_var(s$L, s$mu, s$J, s$sv_add, s$E_obs, s$v_obs, s$n)
  expect_true(is.finite(nll))
})

test_that("FO NLL (cov) increases when L diagonal perturbed away from near-truth", {
  s     <- .fo_setup()
  # Inflate omega by 3x (move away from near-truth)
  L_bad <- s$L * 1.8
  nll0  <- .fo_nll_cov(s$L,     s$mu, s$J, s$sv_add, s$E_obs, s$V_obs, s$n)
  nll1  <- .fo_nll_cov(L_bad,   s$mu, s$J, s$sv_add, s$E_obs, s$V_obs, s$n)
  expect_true(nll1 > nll0)
})

test_that("FO NLL (var) = nll_var_cpp applied to diagonal of FO V", {
  s     <- .fo_setup()
  omega <- s$L %*% t(s$L)
  V_fo  <- s$J %*% omega %*% t(s$J)
  v_fo  <- diag(V_fo) + s$sv_add
  nll_direct <- admixr2:::nll_var_cpp(s$E_obs, s$v_obs, s$mu, v_fo, s$n)
  nll_fn     <- .fo_nll_var(s$L, s$mu, s$J, s$sv_add, s$E_obs, s$v_obs, s$n)
  expect_equal(nll_fn, nll_direct, tolerance = 1e-12)
})

# ==============================================================================
# Omega gradient (cov branch)
# ==============================================================================

test_that("FO Omega L[1,1] gradient (cov) matches FD", {
  s   <- .fo_setup()
  ag  <- .fo_grad_analytic_cov(s$L, s$mu, s$J, s$sv_add, s$E_obs, s$V_obs, s$n)

  delta <- 1e-5
  nll_p <- .fo_nll_cov(.perturb_L_diag(s$L, 1, +delta), s$mu, s$J, s$sv_add, s$E_obs, s$V_obs, s$n)
  nll_m <- .fo_nll_cov(.perturb_L_diag(s$L, 1, -delta), s$mu, s$J, s$sv_add, s$E_obs, s$V_obs, s$n)
  fd <- (nll_p - nll_m) / (2 * delta)

  expect_equal(ag$grad_L11, fd, tolerance = 1e-5)
})

test_that("FO Omega L[2,2] gradient (cov) matches FD", {
  s   <- .fo_setup()
  ag  <- .fo_grad_analytic_cov(s$L, s$mu, s$J, s$sv_add, s$E_obs, s$V_obs, s$n)

  delta <- 1e-5
  nll_p <- .fo_nll_cov(.perturb_L_diag(s$L, 2, +delta), s$mu, s$J, s$sv_add, s$E_obs, s$V_obs, s$n)
  nll_m <- .fo_nll_cov(.perturb_L_diag(s$L, 2, -delta), s$mu, s$J, s$sv_add, s$E_obs, s$V_obs, s$n)
  fd <- (nll_p - nll_m) / (2 * delta)

  expect_equal(ag$grad_L22, fd, tolerance = 1e-5)
})

test_that("FO Omega L[2,1] off-diagonal gradient (cov) matches FD", {
  s   <- .fo_setup()
  ag  <- .fo_grad_analytic_cov(s$L, s$mu, s$J, s$sv_add, s$E_obs, s$V_obs, s$n)

  delta <- 1e-5
  nll_p <- .fo_nll_cov(.perturb_L_off(s$L, 2, 1, +delta), s$mu, s$J, s$sv_add, s$E_obs, s$V_obs, s$n)
  nll_m <- .fo_nll_cov(.perturb_L_off(s$L, 2, 1, -delta), s$mu, s$J, s$sv_add, s$E_obs, s$V_obs, s$n)
  fd <- (nll_p - nll_m) / (2 * delta)

  expect_equal(ag$grad_L21, fd, tolerance = 1e-5)
})

# ==============================================================================
# Sigma gradient (cov branch)
# ==============================================================================

test_that("FO sigma additive gradient (cov) matches FD", {
  s   <- .fo_setup()
  ag  <- .fo_grad_analytic_cov(s$L, s$mu, s$J, s$sv_add, s$E_obs, s$V_obs, s$n)

  # Optimizer parameter for sigma: p_sigma = 2*log(sigma), sigma_var = exp(p_sigma).
  # delta in p_sigma space.
  delta    <- 1e-5
  sv_p     <- s$sv_add * exp(+delta)  # exp(p + delta) = sv * exp(delta) ≈ sv + delta*sv
  sv_m     <- s$sv_add * exp(-delta)
  nll_p    <- .fo_nll_cov(s$L, s$mu, s$J, sv_p, s$E_obs, s$V_obs, s$n)
  nll_m    <- .fo_nll_cov(s$L, s$mu, s$J, sv_m, s$E_obs, s$V_obs, s$n)
  fd       <- (nll_p - nll_m) / (2 * delta)

  expect_equal(ag$grad_sigma, fd, tolerance = 1e-5)
})

test_that("FO sigma proportional gradient (cov) matches FD", {
  s <- .fo_setup()

  # Analytical gradient for proportional sigma (cov branch)
  omega    <- s$L %*% t(s$L)
  V_pred   <- s$J %*% omega %*% t(s$J)
  diag(V_pred) <- diag(V_pred) + s$sv_prop * s$mu^2
  invV     <- chol2inv(chol(V_pred))
  r        <- s$E_obs - s$mu
  V_eff    <- s$V_obs + tcrossprod(r)
  dNLL_dV  <- s$n * (invV - invV %*% V_eff %*% invV)
  grad_sigma_prop_analytic <- sum(diag(dNLL_dV) * s$mu^2) * s$sv_prop

  delta <- 1e-5
  nll_p <- .fo_nll_cov_prop(s$L, s$mu, s$J, s$sv_prop * exp(+delta), s$E_obs, s$V_obs, s$n)
  nll_m <- .fo_nll_cov_prop(s$L, s$mu, s$J, s$sv_prop * exp(-delta), s$E_obs, s$V_obs, s$n)
  fd <- (nll_p - nll_m) / (2 * delta)

  expect_equal(grad_sigma_prop_analytic, fd, tolerance = 1e-5)
})

# ==============================================================================
# Omega gradient (var branch)
# ==============================================================================

test_that("FO Omega L[1,1] gradient (var) matches FD", {
  s   <- .fo_setup()
  ag  <- .fo_grad_analytic_var(s$L, s$mu, s$J, s$sv_add, s$E_obs, s$v_obs, s$n)

  delta <- 1e-5
  nll_p <- .fo_nll_var(.perturb_L_diag(s$L, 1, +delta), s$mu, s$J, s$sv_add, s$E_obs, s$v_obs, s$n)
  nll_m <- .fo_nll_var(.perturb_L_diag(s$L, 1, -delta), s$mu, s$J, s$sv_add, s$E_obs, s$v_obs, s$n)
  fd <- (nll_p - nll_m) / (2 * delta)

  expect_equal(ag$grad_L11, fd, tolerance = 1e-5)
})

test_that("FO Omega L[2,2] gradient (var) matches FD", {
  s   <- .fo_setup()
  ag  <- .fo_grad_analytic_var(s$L, s$mu, s$J, s$sv_add, s$E_obs, s$v_obs, s$n)

  delta <- 1e-5
  nll_p <- .fo_nll_var(.perturb_L_diag(s$L, 2, +delta), s$mu, s$J, s$sv_add, s$E_obs, s$v_obs, s$n)
  nll_m <- .fo_nll_var(.perturb_L_diag(s$L, 2, -delta), s$mu, s$J, s$sv_add, s$E_obs, s$v_obs, s$n)
  fd <- (nll_p - nll_m) / (2 * delta)

  expect_equal(ag$grad_L22, fd, tolerance = 1e-5)
})

test_that("FO Omega L[2,1] off-diagonal gradient (var) matches FD", {
  s   <- .fo_setup()
  ag  <- .fo_grad_analytic_var(s$L, s$mu, s$J, s$sv_add, s$E_obs, s$v_obs, s$n)

  delta <- 1e-5
  nll_p <- .fo_nll_var(.perturb_L_off(s$L, 2, 1, +delta), s$mu, s$J, s$sv_add, s$E_obs, s$v_obs, s$n)
  nll_m <- .fo_nll_var(.perturb_L_off(s$L, 2, 1, -delta), s$mu, s$J, s$sv_add, s$E_obs, s$v_obs, s$n)
  fd <- (nll_p - nll_m) / (2 * delta)

  expect_equal(ag$grad_L21, fd, tolerance = 1e-5)
})

# ==============================================================================
# Sigma gradient (var branch)
# ==============================================================================

test_that("FO sigma additive gradient (var) matches FD", {
  s   <- .fo_setup()
  ag  <- .fo_grad_analytic_var(s$L, s$mu, s$J, s$sv_add, s$E_obs, s$v_obs, s$n)

  delta <- 1e-5
  nll_p <- .fo_nll_var(s$L, s$mu, s$J, s$sv_add * exp(+delta), s$E_obs, s$v_obs, s$n)
  nll_m <- .fo_nll_var(s$L, s$mu, s$J, s$sv_add * exp(-delta), s$E_obs, s$v_obs, s$n)
  fd <- (nll_p - nll_m) / (2 * delta)

  expect_equal(ag$grad_sigma, fd, tolerance = 1e-5)
})

test_that("FO sigma proportional gradient (var) matches FD", {
  s   <- .fo_setup()

  omega  <- s$L %*% t(s$L)
  V_pred <- s$J %*% omega %*% t(s$J)
  v_pred <- diag(V_pred) + s$sv_prop * s$mu^2

  r            <- s$E_obs - s$mu
  dNLL_dv_pred <- s$n * (1/v_pred - (s$v_obs + r^2) / v_pred^2)
  grad_sigma_prop_analytic <- sum(dNLL_dv_pred * s$mu^2) * s$sv_prop

  delta <- 1e-5
  nll_p <- .fo_nll_var_prop(s$L, s$mu, s$J, s$sv_prop * exp(+delta), s$E_obs, s$v_obs, s$n)
  nll_m <- .fo_nll_var_prop(s$L, s$mu, s$J, s$sv_prop * exp(-delta), s$E_obs, s$v_obs, s$n)
  fd <- (nll_p - nll_m) / (2 * delta)

  expect_equal(grad_sigma_prop_analytic, fd, tolerance = 1e-5)
})

# ==============================================================================
# Struct theta mu-path formula
# NOTE: .adfoGrad() uses full FD for struct thetas (both mu-path and V-path).
# These tests verify that the mu-path formula sum(dNLL_dmu * J[:,j]) is
# internally correct when J is held fixed, as a formula regression check.
# ==============================================================================

test_that("FO struct theta mu-path gradient (cov, additive) agrees with mu-only FD", {
  s <- .fo_setup()

  omega  <- s$L %*% t(s$L)
  V_pred <- s$J %*% omega %*% t(s$J)
  diag(V_pred) <- diag(V_pred) + s$sv_add
  invV         <- chol2inv(chol(V_pred))
  r            <- s$E_obs - s$mu
  dNLL_dmu     <- drop(-2 * s$n * invV %*% r)
  grad_tcl_analytic <- sum(dNLL_dmu * s$J[, 1])   # J[:,1] = d/d(eta_cl)

  delta   <- 1e-5
  tcl_p   <- log(5) + delta
  tcl_m   <- log(5) - delta
  mu_p    <- .fo1_mu(tcl_p, log(20), s$dose, s$times)
  mu_m    <- .fo1_mu(tcl_m, log(20), s$dose, s$times)

  # Perturb mu only (J held fixed) — isolates the mu-path contribution.
  V_pred_p <- s$J %*% omega %*% t(s$J); diag(V_pred_p) <- diag(V_pred_p) + s$sv_add
  V_pred_m <- V_pred_p
  nll_p    <- admixr2:::nll_cov_cpp(s$E_obs, s$V_obs, mu_p, V_pred_p, s$n)
  nll_m    <- admixr2:::nll_cov_cpp(s$E_obs, s$V_obs, mu_m, V_pred_m, s$n)
  fd <- (nll_p - nll_m) / (2 * delta)

  expect_equal(grad_tcl_analytic, fd, tolerance = 1e-5)
})

test_that("FO struct theta mu-path gradient (var, additive) agrees with mu-only FD", {
  s <- .fo_setup()

  omega  <- s$L %*% t(s$L)
  V_pred <- s$J %*% omega %*% t(s$J)
  v_pred <- diag(V_pred) + s$sv_add
  r      <- s$E_obs - s$mu
  dNLL_dmu <- -2 * s$n * r / v_pred
  grad_tcl_analytic <- sum(dNLL_dmu * s$J[, 1])

  delta   <- 1e-5
  tcl_p   <- log(5) + delta
  tcl_m   <- log(5) - delta
  mu_p    <- .fo1_mu(tcl_p, log(20), s$dose, s$times)
  mu_m    <- .fo1_mu(tcl_m, log(20), s$dose, s$times)

  nll_p <- admixr2:::nll_var_cpp(s$E_obs, s$v_obs, mu_p, v_pred, s$n)
  nll_m <- admixr2:::nll_var_cpp(s$E_obs, s$v_obs, mu_m, v_pred, s$n)
  fd <- (nll_p - nll_m) / (2 * delta)

  expect_equal(grad_tcl_analytic, fd, tolerance = 1e-5)
})

# ==============================================================================
# Special case: zero IIV (J irrelevant, V_pred = sigma*I diagonal)
# ==============================================================================

test_that("FO cov and var NLL agree when n_eta=0 (V_pred purely diagonal)", {
  s      <- .fo_setup()
  # Zero etas: J = 0 matrix, omega empty -> V_pred = sv_add * I (diagonal)
  J_zero <- matrix(0, s$n_t, 0)
  V_pred_diag <- diag(s$sv_add, s$n_t)   # V_pred = sv*I

  nll_cov <- admixr2:::nll_cov_cpp(s$E_obs, V_pred_diag, s$mu, V_pred_diag, s$n)
  nll_var <- admixr2:::nll_var_cpp(s$E_obs, diag(V_pred_diag), s$mu, diag(V_pred_diag), s$n)

  expect_equal(nll_cov, nll_var, tolerance = 1e-8)
})
