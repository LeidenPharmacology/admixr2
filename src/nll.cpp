// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <cmath>

using namespace Rcpp;
using namespace Eigen;

// ===========================================================================
// Static helpers (not exported)
// ===========================================================================

static inline bool chol_logdet(const MatrixXd& V, MatrixXd& L, double& log_det) {
  LLT<MatrixXd> llt(V);
  if (llt.info() != Success) return false;
  L       = llt.matrixL();
  log_det = 2.0 * L.diagonal().array().log().sum();
  return true;
}

// Batch MVN log-density: rows of bi under N(mean_vec, L L')
static inline VectorXd logdmvnorm_batch_impl(
    const MatrixXd& bi,
    const VectorXd& mean_vec,
    const MatrixXd& L)
{
  int n_sim = bi.rows(), n_eta = bi.cols();
  double log_det   = 2.0 * L.diagonal().array().log().sum();
  double log_const = -0.5 * (n_eta * std::log(2.0 * M_PI) + log_det);
  MatrixXd centered = bi.rowwise() - mean_vec.transpose();
  MatrixXd solved   = L.triangularView<Lower>().solve(centered.transpose());
  // solved: n_eta × n_sim; squaredNorm per column = Mahalanobis dist per sample.
  // Must transpose (1×n_sim) → (n_sim×1) before subtracting from column vector.
  VectorXd mahal = solved.colwise().squaredNorm().transpose();
  return VectorXd::Constant(n_sim, log_const) - 0.5 * mahal;
}

// Numerically stable softmax with non-finite clipping
static inline VectorXd softmax_impl(const VectorXd& lw) {
  int n = lw.size();
  double m = lw.maxCoeff();
  VectorXd w = (lw.array() - m).exp();
  double max_finite = 0.0;
  for (int i = 0; i < n; ++i)
    if (std::isfinite(w[i]) && w[i] > max_finite) max_finite = w[i];
  for (int i = 0; i < n; ++i)
    if (!std::isfinite(w[i])) w[i] = max_finite;
  double s = w.sum();
  if (s > 0.0) w /= s;
  return w;
}

// Weighted mean + covariance: mu = F'w, V = F_c' diag(w) F_c
static inline void weighted_meancov_impl(
    const MatrixXd& F,
    const VectorXd& w,
    VectorXd& mu,
    MatrixXd& V)
{
  mu        = F.transpose() * w;
  MatrixXd Fc = F.rowwise() - mu.transpose();
  V         = Fc.transpose() * w.asDiagonal() * Fc;
}

// Diagonal -2LL (var branch)
static inline double nll_var_impl(
    const VectorXd& E_obs, const VectorXd& v_obs,
    const VectorXd& mu,    const VectorXd& v_pred, double n)
{
  if ((v_pred.array() <= 0.0).any()) return R_PosInf;
  ArrayXd r2  = (E_obs - mu).array().square();
  double  val = (v_pred.array().log() +
                 v_obs.array() / v_pred.array() +
                 r2   / v_pred.array()).sum();
  return n * val;
}

// Cholesky-based -2LL (no explicit inverse)
static inline double nll_cov_impl(
    const VectorXd& E_obs,
    const MatrixXd& V_obs,
    const VectorXd& mu,
    const MatrixXd& V_pred,
    double n)
{
  MatrixXd L; double log_det;
  if (!chol_logdet(V_pred, L, log_det)) return R_PosInf;
  VectorXd r    = E_obs - mu;
  VectorXd Lr   = L.triangularView<Lower>().solve(r);
  double quad   = Lr.squaredNorm();
  // tr(V_obs * V_pred^{-1}) via explicit inverse (small n_times, fast)
  MatrixXd invV = L.triangularView<Lower>().solve(
                    MatrixXd::Identity(V_pred.rows(), V_pred.cols()));
  invV          = L.triangularView<Lower>().adjoint().solve(invV);
  double trace  = (invV.array() * V_obs.array()).sum();
  return n * (log_det + trace + quad);
}

// ===========================================================================
// Exported functions
// ===========================================================================

// ---------------------------------------------------------------------------
// nll_cov_cpp: Cholesky-based MVN -2LL (covariance branch)
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
double nll_cov_cpp(
    const Eigen::VectorXd& E_obs,
    const Eigen::MatrixXd& V_obs,
    const Eigen::VectorXd& E_pred,
    const Eigen::MatrixXd& V_pred,
    double n
) {
  return nll_cov_impl(E_obs, V_obs, E_pred, V_pred, n);
}

// ---------------------------------------------------------------------------
// nll_var_cpp: diagonal MVN -2LL (variance branch)
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
double nll_var_cpp(
    const Eigen::VectorXd& E_obs,
    const Eigen::VectorXd& v_obs,
    const Eigen::VectorXd& E_pred,
    const Eigen::VectorXd& v_pred,
    double n
) {
  if ((v_pred.array() <= 0.0).any()) return R_PosInf;
  ArrayXd r2  = (E_obs - E_pred).array().square();
  double  val = (v_pred.array().log() +
                 v_obs.array() / v_pred.array() +
                 r2 / v_pred.array()).sum();
  return n * val;
}

// ---------------------------------------------------------------------------
// logdmvnorm_batch_cpp: batch MVN log-density
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
Eigen::VectorXd logdmvnorm_batch_cpp(
    const Eigen::MatrixXd& bi_mat,
    const Eigen::VectorXd& mean_vec,
    const Eigen::MatrixXd& L
) {
  return logdmvnorm_batch_impl(bi_mat, mean_vec, L);
}

// ---------------------------------------------------------------------------
// softmax_cpp: numerically stable softmax
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
Eigen::VectorXd softmax_cpp(const Eigen::VectorXd& lw) {
  return softmax_impl(lw);
}

// ---------------------------------------------------------------------------
// weighted_meancov_cpp: weighted ML mean and covariance
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List weighted_meancov_cpp(
    const Eigen::MatrixXd& F,
    const Eigen::VectorXd& w
) {
  VectorXd mu; MatrixXd V;
  weighted_meancov_impl(F, w, mu, V);
  return Rcpp::List::create(Rcpp::Named("mu") = mu,
                             Rcpp::Named("V")  = V);
}

// ---------------------------------------------------------------------------
// compute_mean_new_cpp
//
// Computes mean_new[i] = log_back_fn(struct_paired[i]) - log_origbeta[i]
// without R closure dispatch.  transform_type codes:
//   0 = exp/log   -> log_back_fn(p) = p
//   1 = expit     -> log_back_fn(p) = log(lo + (hi-lo)*sigmoid(p))
//   2 = probitInv -> log_back_fn(p) = log(lo + (hi-lo)*pnorm(p))
//   3 = other     -> log_back_fn(p) = log(p)
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
Eigen::VectorXd compute_mean_new_cpp(
    const Eigen::VectorXd& struct_paired,
    const Eigen::VectorXd& log_origbeta,
    const Eigen::VectorXi& transform_type,
    const Eigen::VectorXd& low,
    const Eigen::VectorXd& hi
) {
  int n = struct_paired.size();
  Eigen::VectorXd result(n);
  for (int i = 0; i < n; ++i) {
    double p = struct_paired[i];
    double lb;
    switch (transform_type[i]) {
      case 1: {
        double s = 1.0 / (1.0 + std::exp(-p));
        lb = std::log(low[i] + (hi[i] - low[i]) * s);
        break;
      }
      case 2:
        lb = std::log(low[i] + (hi[i] - low[i]) * R::pnorm(p, 0.0, 1.0, 1, 0));
        break;
      case 3:
        lb = std::log(p);
        break;
      default:  // 0: exp/log
        lb = p;
        break;
    }
    result[i] = lb - log_origbeta[i];
  }
  return result;
}

// ---------------------------------------------------------------------------
// adm_apply_residual: add the residual error to a predicted mean/variance.
//
// The kernels below are error-model agnostic: R hands them four row-indexed
// arrays (built by .admResidRows) and they evaluate the residual against
// whatever mu they compute. That matters for the IRMC kernel in particular,
// whose mu is the importance-weighted mean and so is not known to R in advance.
//
//   res_form  0 = combined2, 1 = combined1, 2 = lnorm
//   res_a2    additive VARIANCE (a^2), or the lnorm log-variance
//   res_b2    proportional/power VARIANCE (b^2)
//   res_cc    power exponent (1 => proportional)
//
//   combined2   var = a^2 + b^2 * f^(2c)
//   combined1   var = (a + b * f^c)^2
//   lnorm       mu  = f*exp(s/2);  var = mu^2 * (exp(s) - 1)
//
// The c == 1 fast path and the (fe*fe) association keep add/prop/lnorm models
// bit-for-bit identical to the old per-sigma loop this replaced.
// ---------------------------------------------------------------------------
static inline void adm_apply_residual(
    Eigen::VectorXd& mu,
    Eigen::ArrayXd& dv,
    const Eigen::VectorXi& res_form,
    const Eigen::VectorXd& res_a2,
    const Eigen::VectorXd& res_b2,
    const Eigen::VectorXd& res_cc
) {
  const int n_t = static_cast<int>(mu.size());
  for (int t = 0; t < n_t; ++t) {
    const double f = mu[t];
    if (res_form[t] == 2) {                       // lnorm
      const double sv = res_a2[t];
      mu[t] = f * std::exp(sv / 2.0);
      dv[t] += mu[t] * mu[t] * (std::exp(sv) - 1.0);
    } else {
      const double c  = res_cc[t];
      const double fe = (c == 1.0) ? f : std::pow(f, c);
      if (res_form[t] == 1) {                     // combined1 (SD-additive)
        const double s = std::sqrt(res_a2[t]) + std::sqrt(res_b2[t]) * fe;
        dv[t] += s * s;
      } else {                                    // combined2 (variance-additive)
        dv[t] += res_a2[t] + res_b2[t] * (fe * fe);
      }
    }
  }
}

// ---------------------------------------------------------------------------
// irmc_inner_nll_cpp
//
// Fused IRMC inner NLL for one study — now includes sigma and kappa.
//
// Arguments:
//   rawpreds      n_sim x n_times  prediction matrix (fixed proposals)
//   bi_mat        n_sim x n_eta    random-effect samples
//   mean_new      n_eta            log(beta_new/beta_orig) shift
//   L_omega       n_eta x n_eta    lower Cholesky of current Omega
//   log_prop      n_sim            log-density under proposal distribution
//   E_obs         n_times          observed mean
//   V_obs         n_times x n_times observed covariance
//   n             study sample size
//   res_form/a2/b2/cc  n_times      residual error arrays (see adm_apply_residual)
//   kappa_delta   n_times or empty  kappa_fn(struct) - mu_pop (pre-computed in R);
//                                   pass length-0 vector when no kappa correction
//   use_var       int (0/1)        1 = diagonal-only NLL (var branch); 0 = full covariance
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
double irmc_inner_nll_cpp(
    const Eigen::MatrixXd& rawpreds,
    const Eigen::MatrixXd& bi_mat,
    const Eigen::VectorXd& mean_new,
    const Eigen::MatrixXd& L_omega,
    const Eigen::VectorXd& log_prop,
    const Eigen::VectorXd& E_obs,
    const Eigen::MatrixXd& V_obs,
    double n,
    const Eigen::VectorXi& res_form,
    const Eigen::VectorXd& res_a2,
    const Eigen::VectorXd& res_b2,
    const Eigen::VectorXd& res_cc,
    const Eigen::VectorXd& kappa_delta,
    int use_var
) {
  // 1. IS log-weights: log p(bi|Omega_new) - log p(bi|Omega_prop), then softmax.
  //    mean_new = log(beta_new/beta_orig) encodes how mu-referenced paired thetas changed.
  VectorXd log_new = logdmvnorm_batch_impl(bi_mat, mean_new, L_omega);
  VectorXd w       = softmax_impl(log_new - log_prop);

  // 2. Weighted mean + covariance (ML estimator, no n-1 divisor — weights encode distribution).
  VectorXd mu; MatrixXd V;
  weighted_meancov_impl(rawpreds, w, mu, V);

  // 3. Residual error BEFORE kappa: kappa_delta shifts mu deterministically (not
  //    part of the IS average), so the residual is added here while mu still
  //    reflects the IS-weighted structural prediction.
  ArrayXd dv = V.diagonal().array();
  adm_apply_residual(mu, dv, res_form, res_a2, res_b2, res_cc);
  V.diagonal().array() = dv;

  // 4. Kappa correction: shifts mu only, V unchanged. kappa_delta = kappa_fn(struct_cand) - mu_pop,
  //    where mu_pop is fixed per outer iteration and accounts for unpaired theta effects.
  if (kappa_delta.size() > 0)
    mu += kappa_delta;

  // 5. MVN -2LL
  if (use_var) {
    return nll_var_impl(E_obs, V_obs.diagonal(), mu, V.diagonal(), n);
  }
  return nll_cov_impl(E_obs, V_obs, mu, V, n);
}

// ---------------------------------------------------------------------------
// nll_cov_from_samples_cpp: fused NLL from raw prediction matrix (cov branch)
//
// Avoids R-side allocation of cp_c and V by doing all computation in C++.
// Equivalent to: mu=colMeans, V=crossprod(cp_c)/n_sim + sigma, nll_cov_impl.
// Uses ML denominator n_sim (not n_sim-1) consistent with the MVN likelihood.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
double nll_cov_from_samples_cpp(
    const Eigen::MatrixXd& cp_mat,
    const Eigen::VectorXd& E_obs,
    const Eigen::MatrixXd& V_obs,
    double n,
    const Eigen::VectorXi& res_form,
    const Eigen::VectorXd& res_a2,
    const Eigen::VectorXd& res_b2,
    const Eigen::VectorXd& res_cc
) {
  int n_sim = cp_mat.rows();
  VectorXd mu   = cp_mat.colwise().mean();
  MatrixXd cp_c = cp_mat.rowwise() - mu.transpose();  // centred around mu_sim
  MatrixXd V    = cp_c.transpose() * cp_c * (1.0 / n_sim);
  ArrayXd  dv   = V.diagonal().array();
  adm_apply_residual(mu, dv, res_form, res_a2, res_b2, res_cc);
  V.diagonal().array() = dv;
  return nll_cov_impl(E_obs, V_obs, mu, V, n);
}

// ---------------------------------------------------------------------------
// nll_var_from_samples_cpp: fused NLL from raw prediction matrix (var branch)
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
double nll_var_from_samples_cpp(
    const Eigen::MatrixXd& cp_mat,
    const Eigen::VectorXd& E_obs,
    const Eigen::VectorXd& v_obs,
    double n,
    const Eigen::VectorXi& res_form,
    const Eigen::VectorXd& res_a2,
    const Eigen::VectorXd& res_b2,
    const Eigen::VectorXd& res_cc
) {
  int n_sim = cp_mat.rows();
  VectorXd mu   = cp_mat.colwise().mean();
  MatrixXd cp_c = cp_mat.rowwise() - mu.transpose();  // centred around mu_sim
  ArrayXd  pv   = cp_c.array().square().colwise().sum() / n_sim;
  adm_apply_residual(mu, pv, res_form, res_a2, res_b2, res_cc);
  if ((pv <= 0.0).any()) return R_PosInf;
  ArrayXd r2  = (E_obs.array() - mu.array()).square();
  return n * (pv.log() + v_obs.array() / pv + r2 / pv).sum();
}

// ---------------------------------------------------------------------------
// adm_grad_partial_cpp
//
// Scalar gradient contribution for one FD direction (struct or eta-FD path).
// Replaces: dmu=colMeans(dpred); sum(eff_dmu*dmu) + 2*inv_nm1*sum(dNLL_dV*cp_c'dpred)
//
// Arguments:
//   cp_c     n_sim x n_t  centred predictions
//   dpred    n_sim x n_t  FD sensitivity  (cp_hi - cp_lo) / 2h
//   dNLL_dV  n_t  x n_t
//   eff_dmu  n_t           dNLL_dmu + sigma_mu_scale
//   inv_nm1  1 / n_sim  (ML denominator; named inv_nm1 for ABI stability — do NOT pass 1/(n-1))
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
double adm_grad_partial_cpp(
    const Eigen::MatrixXd& cp_c,
    const Eigen::MatrixXd& dpred,
    const Eigen::MatrixXd& dNLL_dV,
    const Eigen::VectorXd& eff_dmu,
    double inv_nm1
) {
  VectorXd dmu = dpred.colwise().mean();
  MatrixXd cpd = cp_c.transpose() * dpred;   // n_t x n_t
  return eff_dmu.dot(dmu) + 2.0 * inv_nm1 * (dNLL_dV.array() * cpd.array()).sum();
}

// ---------------------------------------------------------------------------
// adm_grad_eta_omega_cpp
//
// Gradient contributions for all eta (direct) and omega (Cholesky L) parameters.
// Replaces the R loops that build HD_traces via index tricks and the omega loop
// that allocated dpred_dpar = D_ei * scale per iteration.
//
// Arguments:
//   cp_c            n_sim x n_t    centred predictions
//   D_mat           n_sim x (n_eta*n_t)  sensitivity matrix (do.call(cbind, dpred_list))
//   eta_mat         n_sim x n_eta  diagonal-Cholesky scale for the omega gradient.
//                                  Callers MUST pass sweep(z, 2L, diag(L)/2, "*")
//                                  (= z[:,i] * L_ii/2), NOT the realised effects
//                                  z %*% t(L). The L_ii/2 factor is the
//                                  d(L_ii)/d(log Omega_ii) chain-rule term; passing
//                                  raw z %*% t(L) double-counts the diagonal (the
//                                  factor-of-2 omega-gradient bug). Only used for
//                                  diagonal entries (ei==ej); off-diagonal uses z.
//   z               n_sim x n_eta  standard draws
//   dNLL_dV         n_t  x n_t
//   dNLL_dmu        n_t
//   sigma_mu_scale  n_t            sum_k [sigma_prop_k * 2*dNLL_dV_diag*mu]
//   neta1, neta2    n_o  integer, 1-indexed, from eta_rows_df
//   n_t, n_eta      dimensions
//
// Returns List:
//   eta_grad    n_eta  gradient w.r.t. paired structural theta (eta position)
//   omega_grad  n_o    gradient w.r.t. Cholesky parameters (in eta_rows_df order)
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List adm_grad_eta_omega_cpp(
    const Eigen::MatrixXd& cp_c,
    const Eigen::MatrixXd& D_mat,
    const Eigen::MatrixXd& eta_mat,
    const Eigen::MatrixXd& z,
    const Eigen::MatrixXd& dNLL_dV,
    const Eigen::VectorXd& dNLL_dmu,
    const Eigen::VectorXd& sigma_mu_scale,
    const Eigen::VectorXi& neta1,
    const Eigen::VectorXi& neta2,
    int n_t,
    int n_eta
) {
  int n_sim  = cp_c.rows();
  int n_o    = neta1.size();
  double inv_n = 1.0 / n_sim;
  VectorXd eff_dmu = dNLL_dmu + sigma_mu_scale;   // n_t

  // Every trace term here contracts some row-scaled cp_c against dNLL_dV:
  //
  //   sum_{t1,t2} dNLL_dV[t1,t2] * (X' D)[t1,t2]
  //     = sum_i < (X * dNLL_dV)[i,:], D[i,:] >
  //
  // so hoisting the single gemm cpdV = cp_c * dNLL_dV (n_sim x n_t) out of both
  // loops turns each eta trace into an O(n_sim*n_t) row-dot reduction and each
  // omega trace into an O(n_sim) dot product against the scale vector -- the row
  // scaling factors straight out of the sum. Previously this ran (n_eta + n_o)
  // separate O(n_sim*n_t^2) gemms, each materialising its own n_sim x n_t
  // temporary (n_eta=5 => 5 + 15 = 20 gemms per study per gradient evaluation).
  MatrixXd cpdV(n_sim, n_t);
  cpdV.noalias() = cp_c * dNLL_dV;

  // rowdot.col(j)[i] = < cpdV[i,:], D_j[i,:] >  -- shared by both loops
  MatrixXd rowdot(n_sim, n_eta);
  VectorXd eta_grad(n_eta);
  for (int j = 0; j < n_eta; ++j) {
    auto     D_j   = D_mat.middleCols(j * n_t, n_t);   // view, no copy
    VectorXd dmu_j = D_j.colwise().mean();             // n_t
    rowdot.col(j)  = (cpdV.array() * D_j.array()).rowwise().sum();
    eta_grad[j]    = eff_dmu.dot(dmu_j) + 2.0 * inv_n * rowdot.col(j).sum();
  }

  // Omega gradient: for Cholesky entry (ei,ej),
  //   scale = eta_mat[:,ei] (diagonal) or z[:,ej] (off-diagonal)
  //   trace = sum_i scale[i] * rowdot[i,ei]   (the row scaling factors out)
  VectorXd omega_grad(n_o);
  for (int r = 0; r < n_o; ++r) {
    int ei = neta1[r] - 1;
    int ej = neta2[r] - 1;
    auto     D_ei    = D_mat.middleCols(ei * n_t, n_t);  // view, no copy
    VectorXd scale   = (ei == ej) ? eta_mat.col(ei) : z.col(ej);
    VectorXd dmu_om  = D_ei.transpose() * scale / n_sim; // n_t, no intermediate
    double   trace_o = scale.dot(rowdot.col(ei));        // O(n_sim), no gemm
    omega_grad[r] = eff_dmu.dot(dmu_om) + 2.0 * inv_n * trace_o;
  }

  return Rcpp::List::create(
    Rcpp::Named("eta_grad")   = eta_grad,
    Rcpp::Named("omega_grad") = omega_grad
  );
}

// ---------------------------------------------------------------------------
// adm_grad_partial_var_cpp: var-method version of adm_grad_partial_cpp.
// Takes dNLL_dV_diag as a vector; avoids n_t×n_t matrix and full gemm.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
double adm_grad_partial_var_cpp(
    const Eigen::MatrixXd& cp_c,
    const Eigen::MatrixXd& dpred,
    const Eigen::VectorXd& dNLL_dV_diag,
    const Eigen::VectorXd& eff_dmu,
    double inv_nm1
) {
  VectorXd dmu = dpred.colwise().mean();
  // diag(cp_c' * dpred)[t] = (cp_c .* dpred).colwise().sum()[t]
  VectorXd diag_cpd = (cp_c.array() * dpred.array()).matrix().colwise().sum().transpose();
  return eff_dmu.dot(dmu) + 2.0 * inv_nm1 * dNLL_dV_diag.dot(diag_cpd);
}

// ---------------------------------------------------------------------------
// adm_grad_eta_omega_var_cpp: var-method version of adm_grad_eta_omega_cpp.
// Takes dNLL_dV_diag as a vector; avoids n_t×n_t intermediates in trace terms.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List adm_grad_eta_omega_var_cpp(
    const Eigen::MatrixXd& cp_c,
    const Eigen::MatrixXd& D_mat,
    const Eigen::MatrixXd& eta_mat,
    const Eigen::MatrixXd& z,
    const Eigen::VectorXd& dNLL_dV_diag,
    const Eigen::VectorXd& dNLL_dmu,
    const Eigen::VectorXd& sigma_mu_scale,
    const Eigen::VectorXi& neta1,
    const Eigen::VectorXi& neta2,
    int n_t,
    int n_eta
) {
  int n_sim  = cp_c.rows();
  int n_o    = neta1.size();
  double inv_n = 1.0 / n_sim;
  VectorXd eff_dmu = dNLL_dmu + sigma_mu_scale;

  // Same hoist as the cov branch, with the diagonal dNLL_dV: scaling cp_c's
  // COLUMNS by dNLL_dV_diag once lets every trace collapse to a row-dot, and the
  // omega row-scaling again factors out to a plain dot product. This also removes
  // the n_sim x n_t `cp_c_s` temporary that was allocated per omega entry.
  MatrixXd cpdV = cp_c.array().rowwise() * dNLL_dV_diag.transpose().array();

  MatrixXd rowdot(n_sim, n_eta);
  VectorXd eta_grad(n_eta);
  for (int j = 0; j < n_eta; ++j) {
    auto     D_j   = D_mat.middleCols(j * n_t, n_t);
    VectorXd dmu_j = D_j.colwise().mean();
    rowdot.col(j)  = (cpdV.array() * D_j.array()).rowwise().sum();
    eta_grad[j]    = eff_dmu.dot(dmu_j) + 2.0 * inv_n * rowdot.col(j).sum();
  }

  VectorXd omega_grad(n_o);
  for (int r = 0; r < n_o; ++r) {
    int ei = neta1[r] - 1;
    int ej = neta2[r] - 1;
    auto     D_ei    = D_mat.middleCols(ei * n_t, n_t);
    VectorXd scale   = (ei == ej) ? eta_mat.col(ei) : z.col(ej);
    VectorXd dmu_om  = D_ei.transpose() * scale / n_sim;
    double   trace_o = scale.dot(rowdot.col(ei));
    omega_grad[r] = eff_dmu.dot(dmu_om) + 2.0 * inv_n * trace_o;
  }

  return Rcpp::List::create(
    Rcpp::Named("eta_grad")   = eta_grad,
    Rcpp::Named("omega_grad") = omega_grad
  );
}

// ---------------------------------------------------------------------------
// irmc_grad_kernel_cpp
//
// Computes weight-path gradient quantities for .adirmcInnerGrad.
// All inputs are pre-computed in R (invO, eff_dNLL_dmu, dNLL_dV).
//
// Arguments:
//   F             n_sim x n_times  raw predictions
//   w             n_sim            softmax weights
//   mu            n_times          pre-computed weighted mean (= F'w); avoids recompute
//   d_mat         n_sim x n_eta    bi - mean_new (centred proposals)
//   invO          n_eta x n_eta    Omega^{-1}
//   eff_dNLL_dmu  n_times          effective dNLL/dmu (sigma-prop corrected)
//   dNLL_dV       n_times x n_times dNLL/dV
//
// Returns List with:
//   dNLL_dw        n_sim            gradient w.r.t. unnormalised weights
//   dNLL_dlw       n_sim            gradient w.r.t. log-weights (softmax Jacobian applied)
//   S              n_eta x n_eta    d_mat' diag(dNLL_dlw) d_mat
//   dNLL_dmean_new n_eta            invO %*% d_mat' dNLL_dlw
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List irmc_grad_kernel_cpp(
    const Eigen::MatrixXd& F,
    const Eigen::VectorXd& w,
    const Eigen::VectorXd& mu,        // pre-computed from weighted_meancov_cpp — skips F'w matvec
    const Eigen::MatrixXd& d_mat,
    const Eigen::MatrixXd& invO,
    const Eigen::VectorXd& eff_dNLL_dmu,
    const Eigen::MatrixXd& dNLL_dV
) {
  int n_sim = F.rows();

  // F_c computed from pre-passed mu; avoids recomputing mu = F'w (O(n_sim * n_times) matvec)
  MatrixXd F_c = F.rowwise() - mu.transpose();
  VectorXd term1    = F * eff_dNLL_dmu;                  // n_sim (uncentered, gemv)
  MatrixXd FdV      = F_c * dNLL_dV;                     // n_sim x n_times (gemm)
  VectorXd term2    = (FdV.array() * F_c.array()).rowwise().sum(); // n_sim
  VectorXd dNLL_dw  = term1 + term2;

  // Softmax Jacobian: dNLL/dlw_i = w_i * (dNLL_dw_i - sum_j w_j dNLL_dw_j)
  double wdw        = w.dot(dNLL_dw);
  VectorXd dNLL_dlw = w.array() * (dNLL_dw.array() - wdw);

  // S = d_mat' diag(dNLL_dlw) d_mat  (n_eta x n_eta)
  MatrixXd d_scaled = d_mat.array().colwise() * dNLL_dlw.array(); // n_sim x n_eta
  MatrixXd S        = d_mat.transpose() * d_scaled;               // n_eta x n_eta

  // dNLL/dmean_new = invO %*% colSums(d_mat * dNLL_dlw)
  //                = invO %*% d_mat' dNLL_dlw
  VectorXd dNLL_dmean_new = invO * (d_mat.transpose() * dNLL_dlw);

  return Rcpp::List::create(
    Rcpp::Named("dNLL_dw")        = dNLL_dw,
    Rcpp::Named("dNLL_dlw")       = dNLL_dlw,
    Rcpp::Named("S")              = S,
    Rcpp::Named("dNLL_dmean_new") = dNLL_dmean_new
  );
}

// ---------------------------------------------------------------------------
// adm_col_sq_sum_cpp: column sums of squared entries without a temporary matrix.
// Replaces colSums(cp_c^2) in .admGrad() var branch — avoids allocating the
// n_sim x n_t squared matrix that R's `^` operator produces.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
Eigen::VectorXd adm_col_sq_sum_cpp(const Eigen::MatrixXd& m) {
  return m.colwise().squaredNorm().transpose();
}

// ---------------------------------------------------------------------------
// irmc_grad_kernel_var_cpp: var-branch version of irmc_grad_kernel_cpp.
//
// Takes dNLL_dV_diag (n_times vector) instead of a full n_times x n_times matrix.
// Avoids the R-side diag(dNLL_dV_diag) allocation and reduces the dominant
// F_c * dNLL_dV gemm from O(n_sim * n_t^2) to an O(n_sim * n_t) column-scaled gemv.
// For diagonal dNLL_dV: term2[i] = sum_t F_c[i,t]^2 * dNLL_dV_diag[t].
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List irmc_grad_kernel_var_cpp(
    const Eigen::MatrixXd& F,
    const Eigen::VectorXd& w,
    const Eigen::VectorXd& mu,
    const Eigen::MatrixXd& d_mat,
    const Eigen::MatrixXd& invO,
    const Eigen::VectorXd& eff_dNLL_dmu,
    const Eigen::VectorXd& dNLL_dV_diag
) {
  MatrixXd F_c      = F.rowwise() - mu.transpose();
  VectorXd term1    = F * eff_dNLL_dmu;
  VectorXd term2    = F_c.array().square().matrix() * dNLL_dV_diag;
  VectorXd dNLL_dw  = term1 + term2;

  double   wdw      = w.dot(dNLL_dw);
  VectorXd dNLL_dlw = w.array() * (dNLL_dw.array() - wdw);

  MatrixXd d_scaled       = d_mat.array().colwise() * dNLL_dlw.array();
  MatrixXd S              = d_mat.transpose() * d_scaled;
  VectorXd dNLL_dmean_new = invO * (d_mat.transpose() * dNLL_dlw);

  return Rcpp::List::create(
    Rcpp::Named("dNLL_dw")        = dNLL_dw,
    Rcpp::Named("dNLL_dlw")       = dNLL_dlw,
    Rcpp::Named("S")              = S,
    Rcpp::Named("dNLL_dmean_new") = dNLL_dmean_new
  );
}
