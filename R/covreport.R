# Covariance reporting: rotate optimizer-scale Hessian to reported parameter scale.
# Omega uses non-diagonal Jacobian (Omega = L L'); parFixedDf$SE aligns
# positionally with skipCov; omega dimnames are preserved before upstream output creation.


# Jacobian of the REPORTED omega entries w.r.t. the optimizer's log-Cholesky.
#
# The optimizer holds L (lower-triangular) with diagonal p = log(Omega_ii)
# and off-diagonal raw L_ij; Omega = L L', so:
#   d(Omega_ij)/d(L_ab) = delta_ia L_jb + delta_ja L_ib
#   d(L_ab)/d(p)        = L_aa/2 on diagonal (L_aa = exp(p/2)), 1 off it.
# Returns an (n_report x n_par) matrix transforming as J %*% cov %*% t(J).
# Not reused from upstream: nlmixr2est's R is natural-scale (needs FOCEI EBEs,
# which an aggregate fit has none of); rxOmegaVarCovDeriv() differentiates
# w.r.t. Omega, not the Cholesky parameter; rxSymInvCholCreate() is a Cholesky
# of Omega^-1 with different thetas. Shared: the free-entry (row, col)
# enumeration, pinned to rxOmegaVarCovDeriv()$elements in test-cov-reporting.R.
.admOmegaJacobian <- function(pinfo, L) {
  n_o <- length(pinfo$omega_par)
  if (n_o == 0L) return(NULL)
  ci <- pinfo$chol_i; cj <- pinfo$chol_j; cd <- pinfo$chol_diag
  J  <- matrix(0, n_o, n_o)
  for (r in seq_len(n_o)) {           # reported entry Omega[i, j]
    i <- ci[r]; j <- cj[r]
    for (k in seq_len(n_o)) {         # parameter p_k -> L[a, b]
      a <- ci[k]; b <- cj[k]
      dOm <- (if (i == a) L[j, b] else 0) + (if (j == a) L[i, b] else 0)
      J[r, k] <- dOm * (if (cd[k]) L[a, a] / 2 else 1)
    }
  }
  J
}

# Names for reported omega entries matching nlmixr2est convention (.foceiOmegaCovNames):
# `om.<eta>` on diagonal, `cov.<eta_i>.<eta_j>` off it.
.admOmegaReportNames <- function(pinfo) {
  et <- pinfo$eta_names %||% pinfo$eta_col_names
  if (is.null(et)) return(pinfo$omega_par_names)
  vapply(seq_along(pinfo$omega_par), function(r) {
    i <- pinfo$chol_i[r]; j <- pinfo$chol_j[r]
    if (i == j) paste0("om.", et[i]) else paste0("cov.", et[i], ".", et[j])
  }, character(1))
}

# Scale optimizer-scale covariance to reported scale:
#   - struct thetas: identity scale factor 1;
#   - residual params: natural scale via .admSigmaReportJac();
#   - omega: transformed via J = .admOmegaJacobian(L).
# If ill-conditioned (min(eig) < sqrt(eps) * max(|eig|)), falls back to struct
# + sigma (rows 1:n_sub). Trigger is CONDITIONING not sign: an unidentified
# omega's curvature is FD noise whose sign is an accident -- under Shi21's
# finer step a flat direction can return ~0 from the positive side (measured
# 8.66e-19 vs largest eigenvalue 1.1e-3), so a sign test misses it and an IIV
# at 1e-13 gets an SE of 9.3e-10. sqrt(eps) is the standard inversion
# tolerance; measured ratio here is 1.8e-12, so 1e-12 was too tight to fire.
.ADM_NPD_RCOND <- sqrt(.Machine$double.eps)

# Covariance method predicates for Hessian inversion and sandwich calculation.
# adirmc.R carries its own note that an `==` in place of `%in%` here is the
# same shape of bug as the r,s-vs-r fallback confusion that once cost it a
# silently wrong covariance path.
.admCovWantsHessian  <- function(covMethod) covMethod %in% c("r", "r,s")
.admCovWantsSandwich <- function(covMethod) identical(covMethod, "r,s")

.admReduceNpdOmega <- function(H, H_eigs, eig_dec, nms_cov, n_o, n_sub) {
  .singular <- function(e) {
    if (is.null(e) || !length(e) || any(!is.finite(e))) return(TRUE)
    mx <- max(abs(e))
    if (!is.finite(mx) || mx <= 0) return(TRUE)
    min(e) < .ADM_NPD_RCOND * mx
  }
  if (n_o > 0L && (is.null(eig_dec) || .singular(H_eigs))) {
    .sub <- seq_len(n_sub)
    .Hs  <- H[.sub, .sub, drop = FALSE]
    .es  <- tryCatch(eigen(.Hs, symmetric = TRUE), error = function(e) NULL)
    if (!is.null(.es) && min(.es$values) >= 0)
      return(list(H = .Hs, nms_cov = nms_cov[.sub], np_cov = length(.sub),
                  eig_dec = .es, H_eigs = .es$values, reduced = TRUE))
  }
  list(H = H, nms_cov = nms_cov, np_cov = length(nms_cov),
       eig_dec = eig_dec, H_eigs = H_eigs, reduced = FALSE)
}

.admScaleReportedCov <- function(cov_full, p_hat, pinfo, n_s, n_e, n_o, n_sub) {
  .jac <- rep(1, n_sub)
  if (n_e > 0L)
    .jac[n_s + seq_len(n_e)] <- .admSigmaReportJac(p_hat[n_s + seq_len(n_e)], pinfo)
  .keep <- seq_len(min(n_sub, nrow(cov_full)))
  .out  <- cov_full[.keep, .keep, drop = FALSE] * tcrossprod(.jac[.keep])

  if (n_o > 0L && nrow(cov_full) >= n_sub + n_o) {
    .oi <- n_sub + seq_len(n_o)
    .L_cov <- tryCatch(.admUnpack(p_hat, pinfo)$L, error = function(e) NULL)
    .Jo <- if (is.null(.L_cov)) NULL else
      tryCatch(.admOmegaJacobian(pinfo, .L_cov), error = function(e) NULL)
    if (!is.null(.Jo)) {
      .co <- .Jo %*% cov_full[.oi, .oi, drop = FALSE] %*% t(.Jo)
      .cx <- cov_full[.keep, .oi, drop = FALSE] * .jac[.keep]
      .cx <- .cx %*% t(.Jo)
      .nm <- .admOmegaReportNames(pinfo)
      .big <- rbind(cbind(.out, .cx), cbind(t(.cx), .co))
      dimnames(.big) <- list(c(rownames(.out), .nm), c(colnames(.out), .nm))
      .out <- .big
    }
  }
  .out
}

# Snapshot dimnames on fit$env$cov BEFORE handing it to nlmixr2est: it
# re-dimnames using only theta names, blanking appended omega entries -- and
# foceiFitCpp_ mutates the same SEXP in place, so reading back afterward is too late.
.admCovNames <- function(cov) if (is.matrix(cov)) rownames(cov) else NULL

# Reorder theta rows of covariance into iniDf order so nlmixr2est's positional
# SE assignment maps to the correct parameters. Omega block remains at end.
.admCovThetaOrder <- function(cov, ui) {
  if (!is.matrix(cov) || is.null(rownames(cov))) return(cov)
  .th <- .admThetaIniDf(ui)
  if (is.null(.th)) return(cov)
  nms  <- rownames(cov)
  # the thetas this matrix carries, in iniDf's order; everything else (the omega
  # block) keeps its position after them, which is where nlmixr2est stops looking
  want <- .th$name[.th$name %in% nms]
  if (!length(want)) return(cov)
  ord <- c(match(want, nms), which(!(nms %in% .th$name)))
  cov[ord, ord, drop = FALSE]
}

# iniDf's theta rows, in ntheta order. NULL when there is nothing to order by.
.admThetaIniDf <- function(ui) {
  iniDf <- tryCatch(ui$iniDf, error = function(e) NULL)
  if (is.null(iniDf) || is.null(iniDf$ntheta)) return(NULL)
  .th <- iniDf[!is.na(iniDf$ntheta), , drop = FALSE]
  if (!nrow(.th)) return(NULL)
  .th[order(.th$ntheta), , drop = FALSE]
}

# skipCov indicator (TRUE = theta missing an SE). nlmixr2est's own default is
# version-dependent (older versions also skip every residual-error theta,
# which admixr2 does have an SE for) -- deriving it from the matrix instead is
# what makes that SE print on the older host.
.admCovSkip <- function(cov, ui) {
  if (!is.matrix(cov) || is.null(rownames(cov))) return(NULL)
  .th <- .admThetaIniDf(ui)
  if (is.null(.th)) return(NULL)
  # Indexed by ntheta to match nlmixr2est length expectation.
  skip <- rep(TRUE, max(.th$ntheta))
  skip[.th$ntheta] <- !(.th$name %in% rownames(cov))
  if (all(skip)) return(NULL)
  skip
}

.admRestoreCovNames <- function(fit, nms) {
  if (is.null(nms)) return(invisible(NULL))
  e <- tryCatch(fit$env, error = function(e) NULL)
  if (is.null(e) || is.null(e$cov) || !is.matrix(e$cov)) return(invisible(NULL))
  if (nrow(e$cov) != length(nms)) return(invisible(NULL))
  dimnames(e$cov) <- list(nms, nms)
  invisible(NULL)
}
