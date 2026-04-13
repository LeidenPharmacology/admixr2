# Pure-R reference implementations for cross-validation in tests.

nll_cov_r <- function(E_obs, V_obs, E_pred, V_pred, n) {
  r    <- E_obs - E_pred
  L    <- chol(V_pred)
  logdet <- 2 * sum(log(diag(L)))
  invV   <- chol2inv(L)
  n * (logdet + sum(V_obs * invV) + as.numeric(t(r) %*% invV %*% r))
}

nll_var_r <- function(E_obs, v_obs, E_pred, v_pred, n) {
  r2 <- (E_obs - E_pred)^2
  n * sum(log(v_pred) + v_obs / v_pred + r2 / v_pred)
}

softmax_r <- function(lw) {
  lw_clipped <- ifelse(is.finite(lw), lw, max(lw[is.finite(lw)], na.rm = TRUE))
  w <- exp(lw_clipped - max(lw_clipped))
  w / sum(w)
}

logdmvnorm_r <- function(bi, mean_vec, L) {
  # bi: n_sim x d, L: lower-triangular Cholesky of Omega
  d   <- ncol(bi)
  off <- sweep(bi, 2, mean_vec, "-")
  # solve L z = off_i  for each row
  log_const <- -d / 2 * log(2 * pi) - sum(log(diag(L)))
  apply(off, 1, function(r) {
    z <- forwardsolve(L, r)
    log_const - 0.5 * sum(z^2)
  })
}
