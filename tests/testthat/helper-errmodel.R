# Shims for the C++ kernels' pre-errmodel calling convention.
#
# irmc_inner_nll_cpp / nll_{cov,var}_from_samples_cpp used to take a per-sigma
# (sigma_var, sigma_type) pair and classify the error model themselves. They now
# take the row-indexed residual arrays built by .admResidRows(), which makes them
# error-model agnostic (and is what lets them support combined1/pow at all).
#
# These shims map the old scalar encoding onto the new arrays so the kernel tests
# keep asserting exactly the same maths against the same hand-computed references.
#
# legacy sigma_type: 0 = additive, 1 = proportional, 2 = lognormal
.res_arr <- function(sigma_var, sigma_type, n_t) {
  list(form = rep(if (sigma_type == 2L) 2L else 0L, n_t),
       # add and lnorm both live in the a2 slot (additive variance / log-variance);
       # proportional lives in b2, with exponent 1.
       a2 = rep(if (sigma_type == 1L) 0 else sigma_var, n_t),
       b2 = rep(if (sigma_type == 1L) sigma_var else 0, n_t),
       cc = rep(1, n_t))
}

.nll_cov_fs <- function(cp_mat, E_obs, V_obs, n, sigma_var = 0, sigma_type = 0L) {
  a <- .res_arr(sigma_var, sigma_type, ncol(cp_mat))
  admixr2:::nll_cov_from_samples_cpp(cp_mat, E_obs, V_obs, n, a$form, a$a2, a$b2, a$cc)
}

.nll_var_fs <- function(cp_mat, E_obs, v_obs, n, sigma_var = 0, sigma_type = 0L) {
  a <- .res_arr(sigma_var, sigma_type, ncol(cp_mat))
  admixr2:::nll_var_from_samples_cpp(cp_mat, E_obs, v_obs, n, a$form, a$a2, a$b2, a$cc)
}

.irmc_nll <- function(rawpreds, bi_mat, mean_new, L_omega, log_prop,
                      E_obs, V_obs, n, sigma_var = 0, sigma_type = 0L,
                      kappa_delta = numeric(0), use_var = 0L) {
  a <- .res_arr(sigma_var, sigma_type, length(E_obs))
  admixr2:::irmc_inner_nll_cpp(rawpreds, bi_mat, mean_new, L_omega, log_prop,
                               E_obs, V_obs, n, a$form, a$a2, a$b2, a$cc,
                               kappa_delta, use_var)
}
