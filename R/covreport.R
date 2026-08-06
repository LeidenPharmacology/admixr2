# Covariance reporting: rotating the OPTIMIZER-scale Hessian onto the scale the
# estimates are PRINTED on, and surviving what nlmixr2est does to the result.
#
# Split out of utils.R, which had grown into seven unrelated concerns. Contents
# and order are unchanged -- DESCRIPTION has no Collate field, so R sources these
# alphabetically and every object here is a function referenced only from other
# functions.
#
# The load-bearing facts live with the code they constrain: omega is the block
# whose Jacobian is NOT a per-row factor (Omega = L L'), parFixedDf$SE is filled
# POSITIONALLY against skipCov, and nlmixr2est's C++ blanks the omega dimnames
# in place -- so the names must be snapshot BEFORE nlmixr2CreateOutputFromUi.


# Jacobian of the REPORTED omega entries w.r.t. the optimizer's log-Cholesky.
#
# The optimizer holds L (lower-triangular) with the diagonal stored as
# p = log(Omega_ii) and the off-diagonal as the raw L_ij; .admFullTheta() reports
# the Omega entries themselves. Omega = L L', so
#
#   d(Omega)/d(L_ab)  = E_ab L' + L E_ab'   =>   d(Omega_ij)/d(L_ab)
#                                              = delta_ia L_jb + delta_ja L_ib
#   d(L_ab)/d(p)      = L_aa/2 on the diagonal (L_aa = exp(p/2)), 1 off it.
#
# This is NOT diagonal once omega is correlated, which is why a per-row delta
# factor -- fine for sigma -- cannot be used here. Returns an (n_report x n_par)
# matrix so the covariance transforms as J %*% cov %*% t(J).
#
# WHY NOT REUSE nlmixr2est/rxode2 -- checked, and the answer is structural:
#
#   * nlmixr2est never needs this Jacobian. `.foceiCalcRanalytic()` builds its R
#     matrix on NATURAL-scale directions from the start (via `.omegaVarCovDeriv()`
#     -> `rxode2::rxOmegaVarCovDeriv()`, the derivatives of Omega^-1 and
#     log|Omega| w.r.t. the Omega ELEMENTS), so its result is already on the
#     reported scale -- hence the variable is literally called `.covNat`. There is
#     no transform of ours to borrow. We cannot borrow the assembler either: it
#     needs FOCEI's EBEs, `etaObf` and individual-level `dataSav`, none of which
#     exist in an aggregate fit. admixr2's covariance is a NUMERICAL Hessian of the
#     OPTIMIZER-scale objective, so a change of variables is unavoidable.
#   * `rxode2::rxOmegaVarCovDeriv()` differentiates the wrong way round (w.r.t.
#     Omega's own entries, not w.r.t. a Cholesky parameter).
#   * `rxode2::rxSymInvCholCreate()` is a Cholesky of **Omega^-1** with
#     `diag.xform` in {sqrt, log, identity} -- a different parameterisation from
#     admixr2's (Cholesky of Omega, diagonal held as log(Omega_ii)), so its thetas
#     are not ours and its derivatives are not the ones we need.
#
# What IS shared with upstream, deliberately: the (row, col) ENUMERATION of the
# free Omega entries -- pinned to `rxOmegaVarCovDeriv()$elements` in
# test-cov-reporting.R -- and the resulting names (`.foceiOmegaCovNames`'s
# `om.<eta>` / `cov.<eta_i>.<eta_j>`). The delta-transform PATTERN is theirs too:
# `.postEstimationBoundedTransformJacobian()` does `env$cov <- J cov J'` for
# bounded structural thetas, with a diagonal J; the sigma factors in
# .admSigmaReportJac() are that same pattern, and omega is the case where J is
# not diagonal.
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

# Names for the reported omega entries, matching nlmixr2est's own convention
# exactly (.foceiOmegaCovNames): `om.<eta>` on the diagonal, `cov.<eta_i>.<eta_j>`
# off it, with i the LATER eta -- which is how admixr2 stores chol_i/chol_j too
# (they come from iniDf's neta1/neta2, i.e. the lower triangle).
.admOmegaReportNames <- function(pinfo) {
  et <- pinfo$eta_names %||% pinfo$eta_col_names
  if (is.null(et)) return(pinfo$omega_par_names)
  vapply(seq_along(pinfo$omega_par), function(r) {
    i <- pinfo$chol_i[r]; j <- pinfo$chol_j[r]
    if (i == j) paste0("om.", et[i]) else paste0("cov.", et[i], ".", et[j])
  }, character(1))
}

# Rotate an optimizer-scale covariance onto the scale the ESTIMATES are printed on.
#
# Shared by .adfoCalcCov / .adghCalcCov / .admCalcCov -- it was ~46 identical lines
# in each, the exact "three places to forget when a sigma_role is added" hazard
# CLAUDE.md cites as the reason .admSigmaReportJac was extracted. That extraction
# stopped at the per-row factor; this finishes the job for the surrounding block.
#
# `cov_full` is the (already symmetrised) inverse-Hessian in OPTIMIZER order --
# struct thetas (n_s), then residual params (n_e), then omega Cholesky (n_o), with
# n_sub = n_s + n_e. nlmixr2est prints `Estimate +- 1.96*SE`, so the covariance
# must match .admFullTheta()'s reported scale:
#   - struct thetas: optimizer (log) scale, factor 1;
#   - residual params: NATURAL scale (SD, nu, rho, ...) via the per-row delta
#     factor .admSigmaReportJac() (the derivative of .admSigmaNat(), beside it);
#   - omega: reported as VARIANCE/COVARIANCE entries while the optimizer holds the
#     log-Cholesky -- NOT a per-row factor (dense once omega is correlated), so it
#     takes the full .admOmegaJacobian() and transforms as J cov J', with the
#     struct/sigma-vs-omega cross-block rotated to match. Dropped silently if the
#     omega block was not inverted (the non-PD fallback) or the Jacobian is absent.
# When including the omega Cholesky makes the full Hessian indefinite, drop back to
# the struct+sigma sub-block (rows seq_len(n_sub)) rather than reporting nothing --
# the weakly-identified off-diagonal L is the usual culprit. Shared by all three
# CalcCov functions so the threshold ("full H not positive definite") and the
# retained rows cannot drift apart across estimators (they had: adfo/admc reduced
# eigen-first, adgh via a private invert-first .invert(), on the SAME min-eig < 0
# threshold -- one helper removes the three-way replicate). Returns the (possibly
# reduced) H with its names and eigendecomposition; `reduced` flags whether the
# fallback fired, so each caller can warn once with its own estimator label.
#
# The trigger is CONDITIONING, not sign. `min(H_eigs) < 0` alone was relying on
# an accident: an unidentified omega has no curvature, so the entry it produces
# is numerical junk, and whether that junk lands negative depends on the FD step
# rather than on anything about the model. It did land negative under the old
# fixed step, which is why the sign test appeared to work. Under Shi21's finer,
# better step the same flat direction returns ~0 from the POSITIVE side
# (measured 8.66e-19 against a largest eigenvalue of 1.1e-3, a condition number
# of ~1e16) -- so the sign test silently stopped firing and an IIV estimated at
# 1e-13 was reported with an SE of 9.3e-10, which reads as a precise estimate.
# A singularity test should test for singularity; this one catches the old case
# too, since a negative eigenvalue of that magnitude is also below the tolerance.
# sqrt(eps), the conventional numerical-singularity tolerance: the covariance is
# H inverted, so a reciprocal condition number below this means the reported SE
# carries no significant digits. Measured on the flat-omega case the ratio is
# 1.8e-12 (omega eigenvalues 2.0e-08 against a largest of 11290, i.e. a condition
# number of 5.6e11) -- so a tighter 1e-12 threshold, tried first, did NOT fire.
.ADM_NPD_RCOND <- sqrt(.Machine$double.eps)

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

# Put the dimnames back on fit$env$cov after nlmixr2est has been through it.
#
# nlmixr2est's C++ `foceiFitCpp_` re-dimnames whatever covariance it finds in the
# fit environment using its OWN parameter-name vector, which knows about the
# thetas only -- so the omega rows we append come back named "". This is a known
# shape upstream, not a bug of ours: nlmixr2est ships `.impmapNameCov()` to
# repair exactly these blanks for its importance-sampling estimator, reading the
# omega names off the model. We do the same thing, except we already hold the
# authoritative names (we built the block), so we restore them verbatim.
#
# `nms` must be SNAPSHOT with .admCovNames() before the matrix is handed to
# nlmixr2est: foceiFitCpp_ sets the dimnames attribute IN PLACE on the same SEXP
# (no R-level copy happens when a matrix is merely assigned into an environment),
# so by the time we get here the driver's own `.cov` has been blanked as well.
# Reading the names back off it would restore nothing.
#
# Guarded on the length matching: if a future nlmixr2est returns a covariance of
# a different shape, leaving it untouched is the safe outcome -- a wrongly
# labelled SE is far worse than an unlabelled one.
.admCovNames <- function(cov) if (is.matrix(cov)) rownames(cov) else NULL

# Put the theta rows of the covariance in iniDf's OWN order.
#
# nlmixr2est fills its `SE` column POSITIONALLY: it walks the thetas in iniDf
# order and takes the next entry of `sqrt(diag(cov))` for each one it did not
# skip. admixr2 builds the covariance in OPTIMIZER order -- every structural
# theta, then every residual parameter -- and those two orders agree only when
# the model happens to declare its residual parameters last.
#
#   ini({ a <- 0.1; tcl <- log(3); tv <- log(30) })   # residual declared FIRST
#
# printed `a` with tcl's SE, tcl with tv's and tv with a's: a silent rotation,
# every value finite and plausible. So this is not cosmetic ordering -- it is
# what makes the SE belong to the parameter it is printed beside.
#
# Rows that are not thetas (the appended omega block) keep their position at the
# end. Anything unrecognised is left alone: a covariance we cannot map is better
# reported in the order we built it than permuted on a guess.
#
# The other half of the contract is .admCovSkip(), which tells nlmixr2est WHICH
# thetas this matrix carries -- without it, nlmixr2est < 6.2.0 skips every
# residual-error theta (FOCEI computes its covariance without them) and so reads
# the residual's row as the first structural theta's standard error.
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

# The `skipCov` vector to hand nlmixr2est: TRUE for every theta this covariance
# does NOT carry a standard error for.
#
# Derived from the matrix itself rather than from a convention, so the two cannot
# drift apart. nlmixr2est's own default is version-dependent -- 6.2.0 skips only
# fixed thetas, earlier versions skip every residual-error theta as well, because
# FOCEI's covariance genuinely does not include them. admixr2's does, so saying so
# is what makes the residual SE print at all on the older host, and what stops the
# structural SEs being read off the wrong rows there.
#
# NULL (leave nlmixr2est's default alone) when the model has no thetas or the
# covariance is missing, and when we do not hold every theta nlmixr2est expects to
# find -- the same refusal as .admCovThetaOrder(), for the same reason.
.admCovSkip <- function(cov, ui) {
  if (!is.matrix(cov) || is.null(rownames(cov))) return(NULL)
  .th <- .admThetaIniDf(ui)
  if (is.null(.th)) return(NULL)
  # Indexed BY ntheta, not by row position: nlmixr2est checks the length against
  # max(ntheta) and, if it disagrees, silently substitutes a vector of its own.
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
