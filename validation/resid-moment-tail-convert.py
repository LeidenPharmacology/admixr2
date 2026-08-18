"""Convert the residual moment-tail call sites to the shared helpers.

The helpers live in R/errmodel.R (already applied):

    .admResidMoments(mu_struct, var_f, arr, cov_f = NULL, times = NULL,
                     compose = !is.null(cov_f))
    .admResidSampleMoments(cp, arr, times = NULL)
    .admResidChain(mu_struct, var_f, arr, pinfo, dNLL_dmu, dNLL_dV_diag,
                   dNLL_dV = NULL, cov_f = NULL, times = NULL, deriv = NULL)

This script rewrites the CALL SITES, which live in files errmodel.R's owner
does not own. Every replacement is a literal, asserted, single-occurrence swap,
so it fails loudly rather than silently half-applying if a target has moved.
Prefer it over resid-moment-tail-consumers.patch when line numbers have drifted.

    python validation/resid-moment-tail-convert.py .

Call sites targeted (40 replacements over 7 files):

  R/admc.R
    .admNLL            non-CppOK sample-moment block        -> .admResidSampleMoments
    .admNLLBatch       the same block, verbatim             -> .admResidSampleMoments
    .admGrad  (joint)  deriv/muCoupling/vchain/dV fold      -> .admResidChain
    .admGrad  (joint)  .admSigmaGrad                        -> ch$sigma_grad()
    .admGrad  (single) var + cov objective tails            -> .admResidMoments
    .admGrad  (single) deriv/muCoupling/vchain/dV fold      -> .admResidChain
    .admGrad  (single) .admSigmaGrad                        -> ch$sigma_grad()
    .admGradBatch      the same four, verbatim              -> .admResidMoments/.admResidChain

  R/adgh.R
    .adghMomentsFromCp objective tail                       -> .admResidMoments
    .adghGrad (joint)  deriv/vchain/Bsj fold + sigma        -> .admResidChain
    .adghGrad (single) objective tail                       -> .admResidMoments
    .adghGrad (single) vchain/dV fold per branch + sigma    -> .admResidChain

  R/adfo.R
    .adfoVpred         objective tail                       -> .admResidMoments
    .adfoGrad (joint)  deriv/vchain/ML fold + sigma         -> .admResidChain
    .adfoGrad (single) deriv/vchain/ML/.G/muCoupling/sigma  -> .admResidChain

  R/adirmc.R
    .adirmcInnerNLL    objective tail                       -> .admResidMoments
    .adirmcNLLAndGrad  deriv/muCoupling/vchain/dV + sigma   -> .admResidChain

  R/datagen.R   compute_moments() mc branch                 -> .admResidSampleMoments
  R/plot.R      .admAggData() moment assembly               -> .admResidMoments
  R/studies.R   .admJointResidual()                         -> .admResidMoments

Plus four dead locals the conversion leaves behind (adirmc vchain/.dmv,
adfo vchain/.dmv, adfo dNLL_dv_pred).

Verified bit-identical: validation/resid-moment-tail-bitidentity.R (2568
doubles: NLL + gradient for adfo/adgh/admc/batch/datagen over add, prop,
combined1, combined2, lnorm, pow, prop+t, boxCox, ar, pois x var/cov studies)
and validation/resid-tail-adirmc-probe.R (64 doubles). 0 mismatches.
"""
import sys, os

ROOT = sys.argv[1]

def sub(rel, old, new, count=1):
    p = os.path.join(ROOT, rel)
    s = open(p, encoding='utf-8', newline='').read()
    # this repo mixes LF and CRLF files; match in the file's own convention
    crlf = chr(13) + chr(10)
    if crlf in s:
        old = old.replace(chr(10), crlf)
        new = new.replace(chr(10), crlf)
    n = s.count(old)
    assert n == count, "%s: expected %d, found %d" % (rel, count, n) + chr(10) + old
    s = s.replace(old, new)
    open(p, 'w', encoding='utf-8', newline='').write(s)
    print("  ok %s: %s" % (rel, old.strip().splitlines()[0][:62]))

# ---------------------------------------------------------------- admc.R ----
print("admc.R")

sub('R/admc.R', """      mu_s <- colMeans(cp_mat)
      cpc  <- sweep(cp_mat, 2L, mu_s)
      Vs   <- crossprod(cpc) / nrow(cp_mat)
      ap   <- .admResidApply(mu_s, diag(Vs), ar, s$times, Vs)
      if (identical(s$method, "var")) {
        nll2 <- nll2 + nll_var_cpp(as.numeric(s$E), s$v_diag, ap$mu, ap$dv, s$n)
      } else {
        Vp <- .admApplyResidTail(Vs, ap)
        nll2 <- nll2 + nll_cov_cpp(as.numeric(s$E), s$V, ap$mu, Vp, s$n)
      }
""", """      m <- .admResidSampleMoments(cp_mat, ar, s$times)
      if (identical(s$method, "var")) {
        nll2 <- nll2 + nll_var_cpp(as.numeric(s$E), s$v_diag, m$mu, m$dv, s$n)
      } else {
        nll2 <- nll2 + nll_cov_cpp(as.numeric(s$E), s$V, m$mu, m$V, s$n)
      }
""")

sub('R/admc.R', """          mu_s <- colMeans(cp)
          cpc  <- sweep(cp, 2L, mu_s)
          Vs   <- crossprod(cpc) / nrow(cp)
          ap   <- .admResidApply(mu_s, diag(Vs), ar, s$times, Vs)
          if (identical(s$method, "var")) {
            nll_var_cpp(as.numeric(s$E), s$v_diag, ap$mu, ap$dv, s$n)
          } else {
            Vp <- .admApplyResidTail(Vs, ap)
            nll_cov_cpp(as.numeric(s$E), s$V, ap$mu, Vp, s$n)
          }
""", """          m <- .admResidSampleMoments(cp, ar, s$times)
          if (identical(s$method, "var")) {
            nll_var_cpp(as.numeric(s$E), s$v_diag, m$mu, m$dv, s$n)
          } else {
            nll_cov_cpp(as.numeric(s$E), s$V, m$mu, m$V, s$n)
          }
""")

# .admGrad joint branch
sub('R/admc.R', """      .dres <- .admResidDeriv(mu_struct, var_f, arr, pinfo)   # once, reused below
      sigma_mu_scale <- .admResidMuCoupling(mu_struct, arr, pinfo,
                                            dNLL_dV_diag, dNLL_dmu, var_f,
                                            dNLL_dV, V_struct, .rt_j, deriv = .dres)
      # V_pred -> V_struct chain (see the single-output branch below)
      vchain     <- .admResidVChain(mu_struct, var_f, arr, pinfo,
                                    .admRowTimes(s, length(mu_struct)), deriv = .dres)
      dNLL_dV_s  <- dNLL_dV * vchain
      diag(dNLL_dV_s) <- diag(dNLL_dV_s) +
        dNLL_dmu * (attr(vchain, "dmu_dv0") %||% numeric(n_t))
""", """      # ONE moment tail for this unit: .admResidDeriv once, the V_pred ->
      # V_struct chain, the TBS mean-from-covariance diagonal fold, and both
      # parameter contractions -- see .admResidChain().
      ch <- .admResidChain(mu_struct, var_f, arr, pinfo, dNLL_dmu, dNLL_dV_diag,
                           dNLL_dV, V_struct, .rt_j)
      sigma_mu_scale <- ch$mu_coupling()
      dNLL_dV_s      <- ch$dV
""")

sub('R/admc.R', """      grad[n_s + seq_len(n_e)] <- grad[n_s + seq_len(n_e)] +
        .admSigmaGrad(mu_struct, arr, pinfo, dNLL_dV_diag, dNLL_dmu, var_f,
                    dNLL_dV, .rt_j, V_struct, deriv = .dres)
""", """      grad[n_s + seq_len(n_e)] <- grad[n_s + seq_len(n_e)] + ch$sigma_grad()
""")

# .admGrad single-output objective tail
sub('R/admc.R', """      var_f <- adm_col_sq_sum_cpp(cp_c) / n_sim      # Var_eta(f), pre-residual
      ap <- .admResidApply(mu_struct, var_f, arr)
      mu <- ap$mu; pv <- ap$dv
""", """      var_f <- adm_col_sq_sum_cpp(cp_c) / n_sim      # Var_eta(f), pre-residual
      # No times/cov_f on the diagonal path -- .admResidMoments() enforces it.
      pm <- .admResidMoments(mu_struct, var_f, arr)
      mu <- pm$mu; pv <- pm$dv
""")

sub('R/admc.R', """      cov_f <- crossprod(cp_c) / n_sim               # STRUCTURAL Cov_eta(f); keep it,
      V  <- cov_f                                    # the ms/sigma chain needs it
      var_f <- diag(V)                               # Var_eta(f), pre-residual
      ap <- .admResidApply(mu_struct, var_f, arr, s$times, cov_f)
      mu <- ap$mu
      V  <- .admApplyResidTail(V, ap)
""", """      cov_f <- crossprod(cp_c) / n_sim               # STRUCTURAL Cov_eta(f); keep it,
      var_f <- diag(cov_f)                           # the ms/sigma chain needs it
      pm <- .admResidMoments(mu_struct, var_f, arr, cov_f, s$times)
      mu <- pm$mu
      V  <- pm$V
""")

# .admGrad single-output chain tail
sub('R/admc.R', """    .dres <- .admResidDeriv(mu_struct, var_f, arr, pinfo)   # once, reused across terms
    sigma_mu_scale <- .admResidMuCoupling(mu_struct, arr, pinfo,
                                          dNLL_dV_diag, dNLL_dmu, var_f,
                                          if (is_var) NULL else dNLL_dV, cov_f,
                                          s$times, deriv = .dres)
    eff_dmu <- dNLL_dmu + sigma_mu_scale
""", """    ch <- .admResidChain(mu_struct, var_f, arr, pinfo, dNLL_dmu, dNLL_dV_diag,
                         if (is_var) NULL else dNLL_dV, cov_f, s$times)
    sigma_mu_scale <- ch$mu_coupling()
    eff_dmu <- dNLL_dmu + sigma_mu_scale
""")

sub('R/admc.R', """    vchain <- .admResidVChain(mu_struct, var_f, arr, pinfo, s$times, deriv = .dres)
    # TBS only: mu depends on Var_eta(f), so the mean contributes to the same
    # d(var_f)/d(param) the kernels already chain. Zero for every other form.
    .dmv <- attr(vchain, "dmu_dv0") %||% numeric(n_t)
    dNLL_dV_diag_s <- dNLL_dV_diag * diag(vchain) + dNLL_dmu * .dmv
    if (!is_var) {
      dNLL_dV_s <- dNLL_dV * vchain
      diag(dNLL_dV_s) <- diag(dNLL_dV_s) + dNLL_dmu * .dmv
    }
""", """    # ch$dV_diag / ch$dV already carry the TBS mean-from-covariance term: mu
    # depends on Var_eta(f), which folds onto the diagonal. Zero for every other
    # form, so add() stays bit-identical.
    dNLL_dV_diag_s <- ch$dV_diag
    if (!is_var) dNLL_dV_s <- ch$dV
""")

sub('R/admc.R', """    grad[n_s + seq_len(n_e)] <- grad[n_s + seq_len(n_e)] +
      .admSigmaGrad(mu_struct, arr, pinfo, dNLL_dV_diag, dNLL_dmu, var_f,
                    if (is_var) NULL else dNLL_dV, s$times, cov_f, deriv = .dres)
  }

  grad
}""", """    grad[n_s + seq_len(n_e)] <- grad[n_s + seq_len(n_e)] + ch$sigma_grad()
  }

  grad
}""")

# .admGradBatch objective tail
sub('R/admc.R', """        var_f <- adm_col_sq_sum_cpp(cp_c) / n_sim
        ap <- .admResidApply(mu_struct, var_f, arr)
        mu <- ap$mu; pv <- ap$dv
""", """        var_f <- adm_col_sq_sum_cpp(cp_c) / n_sim
        pm <- .admResidMoments(mu_struct, var_f, arr)   # no times/cov_f: var path
        mu <- pm$mu; pv <- pm$dv
""")

sub('R/admc.R', """        cov_f <- crossprod(cp_c) / n_sim               # keep the STRUCTURAL cov
        V  <- cov_f
        var_f <- diag(V)
        ap <- .admResidApply(mu_struct, var_f, arr, s$times, cov_f)
        mu <- ap$mu
        V  <- .admApplyResidTail(V, ap)
""", """        cov_f <- crossprod(cp_c) / n_sim               # keep the STRUCTURAL cov
        var_f <- diag(cov_f)
        pm <- .admResidMoments(mu_struct, var_f, arr, cov_f, s$times)
        mu <- pm$mu
        V  <- pm$V
""")

# .admGradBatch chain tail
sub('R/admc.R', """      .dres <- .admResidDeriv(mu_struct, var_f, arr, pinfo)   # once, reused across terms
      sigma_mu_scale <- .admResidMuCoupling(mu_struct, arr, pinfo,
                                            dNLL_dV_diag, dNLL_dmu, var_f,
                                            if (is_var) NULL else dNLL_dV, cov_f,
                                            s$times, deriv = .dres)
      eff_dmu <- dNLL_dmu + sigma_mu_scale
      inv_n <- 1 / n_sim
      # V_pred -> V_struct chain (see .admGrad)
      vchain <- .admResidVChain(mu_struct, var_f, arr, pinfo, s$times, deriv = .dres)
      .dmv <- attr(vchain, "dmu_dv0") %||% numeric(n_t)
      dNLL_dV_diag_s <- dNLL_dV_diag * diag(vchain) + dNLL_dmu * .dmv
      if (!is_var) {
        dNLL_dV_s <- dNLL_dV * vchain
        diag(dNLL_dV_s) <- diag(dNLL_dV_s) + dNLL_dmu * .dmv
      }
""", """      # The same moment tail .admGrad() builds -- one object, not a second copy.
      ch <- .admResidChain(mu_struct, var_f, arr, pinfo, dNLL_dmu, dNLL_dV_diag,
                           if (is_var) NULL else dNLL_dV, cov_f, s$times)
      sigma_mu_scale <- ch$mu_coupling()
      eff_dmu <- dNLL_dmu + sigma_mu_scale
      inv_n <- 1 / n_sim
      dNLL_dV_diag_s <- ch$dV_diag
      if (!is_var) dNLL_dV_s <- ch$dV
""")

sub('R/admc.R', """        .admSigmaGrad(mu_struct, arr, pinfo, dNLL_dV_diag, dNLL_dmu, var_f,
                    if (is_var) NULL else dNLL_dV, s$times, cov_f, deriv = .dres)
    }
  }
""", """        ch$sigma_grad()
    }
  }
""")

# ---------------------------------------------------------------- adgh.R ----
print("adgh.R")

sub('R/adgh.R', """  ap  <- .admResidApply(mu, diag(V), arr, times)
  list(E = ap$mu, V = .admApplyResidTail(V, ap))
""", """  m <- .admResidMoments(mu, diag(V), arr, V, times)
  list(E = m$mu, V = m$V)
""")

# .adghGrad joint branch: the tail moves down to where the score vectors exist
sub('R/adgh.R', """      dres   <- .admResidDeriv(mu, var_f, arr, pinfo)
      ls_vec <- dres$dmu_df            # d(mu_pred)/df -- see the single-output branch

      vchain <- .admResidVChain(mu, var_f, arr, pinfo,
                                .admRowTimes(s, length(mu)), deriv = dres)
      jr <- .admJointResidual(mu, V_str, s, pinfo, pars$sigma_var)
""", """      jr <- .admJointResidual(mu, V_str, s, pinfo, pars$sigma_var)
""")

sub('R/adgh.R', """      Bdiag <- diag(B)
      Bsj   <- B * vchain
      diag(Bsj) <- diag(Bsj) +            # mean-from-covariance path (TBS only)
        dNLL_dmu_sig * (attr(vchain, "dmu_dv0") %||% numeric(length(mu)))
      Bt <- cpc %*% Bsj
""", """      Bdiag <- diag(B)
      # ONE moment tail for this unit: .admResidDeriv, the V_pred -> V_struct
      # chain, the TBS mean-from-covariance diagonal fold, the sigma contraction.
      ch     <- .admResidChain(mu, var_f, arr, pinfo, dNLL_dmu_sig, Bdiag, B,
                               V_str, .admRowTimes(s, length(mu)))
      dres   <- ch$deriv
      ls_vec <- ch$dmu_df              # d(mu_pred)/df -- see the single-output branch
      Bsj    <- ch$dV
      Bt <- cpc %*% Bsj
""")

sub('R/adgh.R', """        .admSigmaGrad(mu, arr, pinfo, Bdiag, dNLL_dmu_sig, var_f, B,
                      .admRowTimes(s, length(mu)), V_str, deriv = dres)
""", """        ch$sigma_grad()
""")

# .adghGrad single-output objective tail
sub('R/adgh.R', """    ap    <- .admResidApply(mu, var_f, arr, s$times, cov_f)
    V <- .admApplyResidTail(V, ap)
    mu_sigma <- ap$mu
""", """    pm <- .admResidMoments(mu, var_f, arr, cov_f, s$times)
    V  <- pm$V
    mu_sigma <- pm$mu
""")

# .adghGrad single-output chain tail: one chain per branch, sharing `dres`
sub('R/adgh.R', """    lnorm_scale <- dres$dmu_df
    vchain      <- .admResidVChain(mu, var_f, arr, pinfo, s$times, deriv = dres)
""", """    lnorm_scale <- dres$dmu_df
""")

sub('R/adgh.R', """      dNLL_dV_dg_s  <- dNLL_dV_diag * diag(vchain) +      # -> d(NLL)/d(var_f)
        dNLL_dmu_sig * (attr(vchain, "dmu_dv0") %||% numeric(length(mu)))
""", """      # The moment tail, once: the V_pred -> V_struct chain plus the TBS
      # mean-from-covariance diagonal. NULL dNLL_dV/cov_f: diagonal path.
      ch            <- .admResidChain(mu, var_f, arr, pinfo, dNLL_dmu_sig,
                                      dNLL_dV_diag, NULL, NULL, s$times,
                                      deriv = dres)
      dNLL_dV_dg_s  <- ch$dV_diag                        # -> d(NLL)/d(var_f)
""")

sub('R/adgh.R', """      Bdiag     <- diag(B)
      Bs        <- B * vchain
      diag(Bs)  <- diag(Bs) +             # mean-from-covariance path (TBS only)
        dNLL_dmu_sig * (attr(vchain, "dmu_dv0") %||% numeric(length(mu)))
      Bt        <- cpc %*% Bs             # Q x n_t; chained to V_struct
""", """      Bdiag     <- diag(B)
      ch        <- .admResidChain(mu, var_f, arr, pinfo, dNLL_dmu_sig, Bdiag, B,
                                  cov_f, s$times, deriv = dres)
      Bs        <- ch$dV                  # mean-from-covariance fold included
      Bt        <- cpc %*% Bs             # Q x n_t; chained to V_struct
""")

sub('R/adgh.R', """      .admSigmaGrad(mu, arr, pinfo, Bvec, dNLL_dmu_sig, var_f,
                    if (is_var) NULL else B, s$times,
                    if (is_var) NULL else cov_f, deriv = dres)
""", """      ch$sigma_grad()
""")

# ---------------------------------------------------------------- adfo.R ----
print("adfo.R")

sub('R/adfo.R', """  ap <- .admResidApply(mu_pred, diag(V), arr, times, V)
  V <- .admApplyResidTail(V, ap)
  list(V = V, mu_sigma = ap$mu, JL = JL, ms = ap$ms, var_f = diag(tcrossprod(JL)))
""", """  pm <- .admResidMoments(mu_pred, diag(V), arr, V, times)
  list(V = pm$V, mu_sigma = pm$mu, JL = JL, ms = pm$ms,
       var_f = diag(tcrossprod(JL)))
""")

sub('R/adfo.R', """      var_f_j  <- diag(tcrossprod(JL))
      .dres_j  <- .admResidDeriv(mu_pred, var_f_j, mc$arr, pinfo)   # once, reused below
      vchain_j <- .admResidVChain(mu_pred, var_f_j, mc$arr, pinfo,
                                  .admRowTimes(s, length(mu_pred)), deriv = .dres_j)
      dNLL_dmu_full <- drop(-2 * s$n * invV %*% r)
      if (n_eta > 0L && n_o > 0L) {
        .bj <- dNLL_dV * vchain_j
        diag(.bj) <- diag(.bj) +           # mean-from-covariance path (TBS only)
          dNLL_dmu_full * (attr(vchain_j, "dmu_dv0") %||% numeric(length(mu_pred)))
        ML <- crossprod(J, .bj %*% JL)
""", """      .cov_f_j <- tcrossprod(JL)              # STRUCTURAL J Omega J'
      var_f_j  <- diag(.cov_f_j)
      dNLL_dmu_full <- drop(-2 * s$n * invV %*% r)
      # ONE moment tail for this unit: .admResidDeriv, the V_pred -> V_struct
      # chain, the TBS mean-from-covariance diagonal fold, the sigma contraction.
      ch_j <- .admResidChain(mu_pred, var_f_j, mc$arr, pinfo, dNLL_dmu_full,
                             dNLL_dV_diag, dNLL_dV, .cov_f_j,
                             .admRowTimes(s, length(mu_pred)))
      if (n_eta > 0L && n_o > 0L) {
        ML <- crossprod(J, ch_j$dV %*% JL)
""")

sub('R/adfo.R', """        .admSigmaGrad(mu_pred, mc$arr, pinfo, dNLL_dV_diag, dNLL_dmu_full, var_f_j,
                      dNLL_dV, .admRowTimes(s, length(mu_pred)), tcrossprod(JL),
                      deriv = .dres_j)
""", """        ch_j$sigma_grad()
""")

sub('R/adfo.R', """    var_f  <- vp$var_f
    .dres  <- .admResidDeriv(mu_pred, var_f, mc$arr, pinfo)   # once, reused below
    vchain <- .admResidVChain(mu_pred, var_f, mc$arr, pinfo, s$times, deriv = .dres)
    # Needed by the omega block below as well as the sigma block, so compute it here.
    dNLL_dmu <- if (is_var) -2 * s$n * r / diag(V_pred) else
      drop(-2 * s$n * invV %*% r)
""", """    var_f  <- vp$var_f
    # Needed by the omega block below as well as the sigma block, so compute it here.
    dNLL_dmu <- if (is_var) -2 * s$n * r / diag(V_pred) else
      drop(-2 * s$n * invV %*% r)
    # cov_f: the STRUCTURAL covariance J Omega J'. NULL on the diagonal path, so
    # the off-diagonal ms terms stay skipped exactly as they were.
    .cov_f <- if (is_var) NULL else tcrossprod(mc$JL)
    ch <- .admResidChain(mu_pred, var_f, mc$arr, pinfo, dNLL_dmu, dNLL_dV_diag,
                         if (is_var) NULL else dNLL_dV, .cov_f, s$times)
    vchain <- ch$vchain
""")

sub('R/adfo.R', """    .dmv <- attr(vchain, "dmu_dv0") %||% numeric(length(mu_pred))
""", """    .dmv <- ch$dmu_dv0
""")

sub('R/adfo.R', """      ML <- if (is_var)
        crossprod(J, JL * (dNLL_dv_pred * diag(vchain) + dNLL_dmu * .dmv))
      else {
        .b <- dNLL_dV * vchain
        diag(.b) <- diag(.b) + dNLL_dmu * .dmv
        crossprod(J, .b %*% JL)
      }
""", """      ML <- if (is_var) crossprod(J, JL * ch$dV_diag)
            else            crossprod(J, ch$dV %*% JL)
""")

sub('R/adfo.R', """      .G  <- if (is_var) .JO * (dNLL_dv_pred * diag(vchain) + dNLL_dmu * .dmv)
             else {
               .bt <- dNLL_dV * vchain
               diag(.bt) <- diag(.bt) + dNLL_dmu * .dmv
               .bt %*% .JO
             }
""", """      .G  <- if (is_var) .JO * ch$dV_diag else ch$dV %*% .JO
""")

sub('R/adfo.R', """      .dmu_eff <- dNLL_dmu +
        .admResidMuCoupling(mu_pred, mc$arr, pinfo, dNLL_dV_diag, dNLL_dmu, var_f,
                            if (is_var) NULL else dNLL_dV,
                            if (is_var) NULL else tcrossprod(mc$JL), s$times,
                            deriv = .dres)
""", """      .dmu_eff <- dNLL_dmu + ch$mu_coupling()
""")

sub('R/adfo.R', """      .admSigmaGrad(mu_pred, mc$arr, pinfo, dNLL_dV_diag, dNLL_dmu, var_f,
                    if (is_var) NULL else dNLL_dV, s$times,
                    if (is_var) NULL else tcrossprod(mc$JL), deriv = .dres)
""", """      ch$sigma_grad()
""")

# -------------------------------------------------------------- adirmc.R ----
print("adirmc.R")

sub('R/adirmc.R', """    ap    <- .admResidApply(mu_struct, var_f, arr, s$times, cov_f)
""", """    ap    <- .admResidMoments(mu_struct, var_f, arr, cov_f, s$times)$ap
""")

sub('R/adirmc.R', """    .is_var_s <- identical(s$method, "var")
    .dres        <- .admResidDeriv(mu_struct, var_f, arr, pinfo)   # once, reused below
    eff_dNLL_dmu <- dNLL_dmu +
      .admResidMuCoupling(mu_struct, arr, pinfo, dNLL_dV_diag, dNLL_dmu, var_f,
                          if (.is_var_s) NULL else dNLL_dV,
                          if (.is_var_s) NULL else cov_f, s$times, deriv = .dres)
""", """    .is_var_s <- identical(s$method, "var")
    # ONE moment tail for this study: .admResidDeriv, the V_pred -> V_struct
    # chain, the TBS mean-from-covariance fold and both contractions.
    ch <- .admResidChain(mu_struct, var_f, arr, pinfo, dNLL_dmu, dNLL_dV_diag,
                         if (.is_var_s) NULL else dNLL_dV,
                         if (.is_var_s) NULL else cov_f, s$times)
    .dres        <- ch$deriv
    eff_dNLL_dmu <- dNLL_dmu + ch$mu_coupling()
""")

sub('R/adirmc.R', """    vchain         <- .admResidVChain(mu_struct, var_f, arr, pinfo, s$times, deriv = .dres)
    .dmv           <- attr(vchain, "dmu_dv0") %||% numeric(length(mu_struct))
    dNLL_dV_diag_s <- dNLL_dV_diag * diag(vchain) + dNLL_dmu * .dmv
""", """    vchain         <- ch$vchain
    .dmv           <- ch$dmu_dv0
    dNLL_dV_diag_s <- ch$dV_diag
""")

sub('R/adirmc.R', """                           { .b <- dNLL_dV * vchain
                             diag(.b) <- diag(.b) + dNLL_dmu * .dmv; .b })
""", """                           ch$dV)
""")

sub('R/adirmc.R', """      .admSigmaGrad(mu_struct, arr, pinfo, dNLL_dV_diag, dNLL_dmu, var_f,
                    if (.is_var_s) NULL else dNLL_dV, s$times,
                    if (.is_var_s) NULL else cov_f, deriv = .dres)
""", """      ch$sigma_grad()
""")

# ------------------------------------------------------------- datagen.R ----
print("datagen.R")

sub('R/datagen.R', """        mu   <- colMeans(cp_mat)
        cp_c <- sweep(cp_mat, 2L, mu)
        V    <- crossprod(cp_c) / control$n_sim
""", """""")

sub('R/datagen.R', """        ap   <- .admResidApply(mu, diag(V), arr, study_tmp$times, V)
        list(mu = ap$mu, V = .admApplyResidTail(V, ap), cp_mat = cp_mat)
""", """        m <- .admResidSampleMoments(cp_mat, arr, study_tmp$times)
        list(mu = m$mu, V = m$V, cp_mat = cp_mat)
""")

# ---------------------------------------------------------------- plot.R ----
print("plot.R")

sub('R/plot.R', """    ap  <- .admResidApply(mu, diag(V), arr, times, V)
    list(V = .admApplyResidTail(V, ap), mu = ap$mu)
""", """    m <- .admResidMoments(mu, diag(V), arr, V, times)
    list(V = m$V, mu = m$mu)
""")

# ------------------------------------------------------------- studies.R ----
print("studies.R")

sub('R/studies.R', """  ap  <- .admResidApply(mu_struct, diag(V_pred), arr, rt, V_pred)
  list(mu = ap$mu, V = .admApplyResidTail(V_pred, ap))
""", """  m <- .admResidMoments(mu_struct, diag(V_pred), arr, V_pred, rt)
  list(mu = m$mu, V = m$V)
""")

# dead locals left behind by the conversion
sub('R/adirmc.R', """    vchain         <- ch$vchain
    .dmv           <- ch$dmu_dv0
    dNLL_dV_diag_s <- ch$dV_diag
""", """    dNLL_dV_diag_s <- ch$dV_diag
""")

sub('R/adfo.R', """                         if (is_var) NULL else dNLL_dV, .cov_f, s$times)
    vchain <- ch$vchain
""", """                         if (is_var) NULL else dNLL_dV, .cov_f, s$times)
""")

sub('R/adfo.R', """    .dmv <- ch$dmu_dv0
""", """""")

sub('R/adfo.R', """      v_pred       <- diag(V_pred)
      dNLL_dv_pred <-      s$n * (1/v_pred - (s$v_diag + r^2) / v_pred^2)
      dNLL_dV_diag <- dNLL_dv_pred
""", """      v_pred       <- diag(V_pred)
      dNLL_dV_diag <- s$n * (1/v_pred - (s$v_diag + r^2) / v_pred^2)
""")
print("done")
