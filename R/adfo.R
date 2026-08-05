# -- FO (First-Order) aggregate data estimator ---------------------------------
# Approximates the aggregate NLL analytically:
#   mu_pred = f(theta, 0)        population prediction at eta = 0
#   J[t, j] = df/d(eta_j)|_0    Jacobian via sensitivity equations or FD
#   V_pred  = J Omega J^T + sigma_contribution
#
# One rxSolve call (via sens model) per study per NLL evaluation, vs n_sim for
# MC. Substantially faster per iteration; suitable for initial estimates or
# models where the FO approximation is adequate.

# -- Internal helpers ----------------------------------------------------------

# Get mu = f(theta, 0) and J = df/d(eta)|_0 for one study.
# Returns list(mu, J): mu is length-n_t, J is n_t x n_eta.
# Prefers sens model (one pass); falls back to rxMod + forward FD.
.adfoGetMuJ <- function(pars, pinfo, s, sensModel, rxMod, output_var, params_mat, cores) {
  sm  <- matrix(pars$struct, nrow = 1L,
                dimnames = list(NULL, names(pars$struct)))
  res <- .adfoGetMuJBatch(sm, pinfo, s, sensModel, rxMod, output_var, cores,
                          pars$sigma_var)
  if (is.null(res)) return(NULL)
  res[[1L]]
}

# (mu, J) for a whole SET of structural-theta configurations, in ONE rxSolve.
#
# struct_mat is n_cfg x n_struct (natural scale). Returns a length-n_cfg list of
# list(mu, J), one per row -- i.e. n_cfg independent FO linearisations.
#
# This is the primitive behind the FD batching. Previously each configuration
# cost its own rxSolve, and under FO every one of those solves carries a single
# subject: an ~11 ms fixed call cost to do ~0.015 ms of integration. Stacking the
# configurations as rows collapses them into one call and, as a bonus, finally
# gives rxSolve's OpenMP something to parallelise over (under FO it previously
# had exactly one subject, so `cores` did nothing at all).
#
# The FD fallback (no sens model) stacks (1 + n_eta) rows per configuration --
# base plus one perturbed eta per dimension -- so it too is a single solve.
.adfoGetMuJBatch <- function(struct_mat, pinfo, s, sensModel, rxMod, output_var, cores,
                             sigma_var = NULL) {
  n_cfg <- nrow(struct_mat)
  n_eta <- pinfo$n_eta
  n_t   <- length(s$times)

  if (n_eta > 0L && !is.null(sensModel)) {
    eta0 <- matrix(0, n_cfg, n_eta, dimnames = list(NULL, pinfo$eta_col_names))
    res  <- .admSimulateSensRows(sensModel, struct_mat, pinfo$sigma_names, eta0, s,
                                 cores, pinfo$nDisplayProgress, sigma_var,
                                 pinfo$sigdig)
    if (!is.null(res))
      return(lapply(seq_len(n_cfg), function(k)
        list(mu = res$cp_mat[k, ],
             J  = matrix(vapply(res$dpred_list, function(D) D[k, ], numeric(n_t)),
                         nrow = n_t, ncol = n_eta),
             # dJ[[dir]] = d(J)/d(dir), n_t x n_eta -- the second-order block, or
             # NULL from an order-1 sens model (then the struct thetas are FD'd).
             dJ = if (is.null(res$d2_list)) NULL else
               lapply(res$d2_list, function(bl)
                 matrix(vapply(bl, function(D) D[k, ], numeric(n_t)),
                        nrow = n_t, ncol = n_eta)),
             # d(mu)/d(theta) for the UNPAIRED thetas. adfo ignored these while its
             # struct gradient was FD (a paired theta reads its eta column out of J,
             # an unpaired one has no column in J at all); the analytic mean path
             # below needs them.
             dth = if (is.null(res$dtheta_list)) NULL else
               lapply(res$dtheta_list, function(D) D[k, ]))))
  }

  # FD fallback, or n_eta == 0. One solve either way.
  if (n_eta == 0L) {
    eta_e <- matrix(numeric(0), n_cfg, 0L)
    pm    <- .admMakeParamsList(n_cfg, pinfo, 1L)[[1L]]
    vals  <- .admSimulateRows(rxMod, struct_mat, pinfo$sigma_names, eta_e, s,
                              output_var, pm, cores, pinfo$nDisplayProgress,
                              pinfo$sigdig)
    return(lapply(seq_len(n_cfg), function(k)
      list(mu = vals[k, ], J = matrix(0, n_t, 0))))
  }

  eps    <- 1e-6
  n_blk  <- 1L + n_eta                      # base + one perturbed eta per dim
  n_row  <- n_cfg * n_blk
  sm_big <- struct_mat[rep(seq_len(n_cfg), each = n_blk), , drop = FALSE]
  eta_big <- matrix(0, n_row, n_eta, dimnames = list(NULL, pinfo$eta_col_names))
  for (j in seq_len(n_eta))
    eta_big[seq(from = 1L + j, to = n_row, by = n_blk), j] <- eps

  pm   <- .admMakeParamsList(n_row, pinfo, 1L)[[1L]]
  vals <- .admSimulateRows(rxMod, sm_big, pinfo$sigma_names, eta_big, s,
                           output_var, pm, cores, pinfo$nDisplayProgress,
                              pinfo$sigdig)

  lapply(seq_len(n_cfg), function(k) {
    base <- (k - 1L) * n_blk + 1L
    mu   <- vals[base, ]
    J    <- matrix(0, n_t, n_eta)
    for (j in seq_len(n_eta)) J[, j] <- (vals[base + j, ] - mu) / eps
    list(mu = mu, J = J)
  })
}

# Build V_pred and mu_sigma (lnorm-corrected mean) for one study.
# Returns list(V, mu_sigma, JL): JL = J %*% L cached for gradient reuse.
# V = tcrossprod(JL) + the residual contribution to diag(V). `arr` is the row
# array from .admResidRows(); it, not this function, knows the error model.
.adfoVpred <- function(mu_pred, J, L, arr, n_t, n_eta, times = NULL) {
  JL <- if (n_eta > 0L) J %*% L else matrix(0, n_t, 0)
  V  <- if (n_eta > 0L) tcrossprod(JL) else matrix(0, n_t, n_t)
  # `times` and the structural covariance are REQUIRED by any residual form that
  # reaches the off-diagonal (ar's rho^|dt|, ordinal's cross-category term). This
  # used to pass neither and never add ap$rmat, so adfo scored an INDEPENDENT
  # residual: its NLL was exactly invariant in rho, while .adfoGrad() -- which does
  # pass s$times to .admSigmaGrad -- returned a non-zero rho gradient. The optimiser
  # therefore walked a direction the objective could not move along, and adfo
  # reported a completely different objective from adgh/admc on identical data.
  ap <- .admResidApply(mu_pred, diag(V), arr, times, V)
  V <- .admApplyResidTail(V, ap)
  list(V = V, mu_sigma = ap$mu, JL = JL, ms = ap$ms, var_f = diag(tcrossprod(JL)))
}

# -- NLL -----------------------------------------------------------------------

# FO population mean + Jacobian for a JOINT (same-subject) unit: mu = f(theta,0)
# stacked across outputs, J = d(stacked)/d(eta)|_0 (n_total x n_eta). Prefers the
# sensitivity model (analytical J, one solve per output at eta = 0) so the
# struct-theta FD of the joint NLL is not corrupted by an inner FD Jacobian;
# falls back to finite differences when no sens model is available.
.adfoGetMuJJoint <- function(pars, pinfo, unit, sensModel, rxMod, params_mat, cores) {
  n_eta   <- pinfo$n_eta
  n_total <- unit$n_total
  eta0 <- matrix(0, 1L, n_eta, dimnames = list(NULL, pinfo$eta_col_names))
  if (n_eta > 0L && !is.null(sensModel)) {
    res <- .admSimulateJointSens(sensModel, pars$struct, pinfo$sigma_names,
                                 eta0, unit, cores, pinfo$nDisplayProgress,
                                 pars$sigma_var, pinfo$sigdig)
    if (!is.null(res))
      return(list(mu = as.numeric(res$cp_mat),
                  J  = do.call(cbind, lapply(res$dpred_list, as.numeric))))
  }
  mu <- as.numeric(.admSimulateJoint(rxMod, pars$struct, pinfo$sigma_names,
                                     eta0, unit, params_mat, cores,
                                     pinfo$nDisplayProgress, pinfo$sigdig))
  if (n_eta == 0L) return(list(mu = mu, J = matrix(0, n_total, 0)))
  eps <- 1e-6; J <- matrix(0, n_total, n_eta)
  for (j in seq_len(n_eta)) {
    etap <- eta0; etap[1L, j] <- eps
    mup  <- as.numeric(.admSimulateJoint(rxMod, pars$struct, pinfo$sigma_names,
                                         etap, unit, params_mat, cores,
                                         pinfo$nDisplayProgress, pinfo$sigdig))
    J[, j] <- (mup - mu) / eps
  }
  list(mu = mu, J = J)
}

# (mu, J) memo for the FO solve.
#
# The FO solve depends ONLY on the structural thetas: eta is pinned at 0 and the
# sigmas are zeroed into it, while Omega and sigma enter afterwards, analytically
# (V = J Omega J' + Sigma). Verified empirically: perturbing any omega or sigma
# parameter leaves mu and J bit-for-bit unchanged.
#
# So an FD direction that perturbs an omega or a sigma needs NO solve at all --
# it reuses the base (mu, J). That is what makes the FD gradient and the NLL-FD
# Hessian cheap: their solve count collapses to the number of DISTINCT structural
# vectors they visit, which for the Hessian is far smaller than the number of
# NLL evaluations.
#
# `cache` is an environment created per top-level call (never persisted across
# parameter points), so there is no staleness risk.
# `struct` is the KEY, not the parameters: callers append an estimated Box-Cox /
# Yeo-Johnson lambda to it, because that lambda reaches rx_pred_ and so the cached
# (mu, J) is no longer a function of the structural thetas alone. Everything else
# about FO's "omega/sigma cost no solve" property is unchanged.
.adfoMuJKey <- function(pars, sensModel) {
  .tbn <- if (is.null(sensModel)) NULL else sensModel$pred_tbs$lam_name
  if (!is.null(.tbn) && !is.na(.tbn) && .tbn %in% names(pars$sigma_var))
    c(pars$struct, unname(pars$sigma_var[[.tbn]]))
  else pars$struct
}

.adfoMuJMemo <- function(cache, s_idx, struct, solve_fn) {
  if (is.null(cache)) return(solve_fn())
  key <- paste0(s_idx, ":", paste(sprintf("%.17g", struct), collapse = ","))
  hit <- cache[[key]]
  if (!is.null(hit)) return(hit)
  val <- solve_fn()
  cache[[key]] <- val
  val
}

#' @noRd
.adfoNLL <- function(p, pinfo, studies, sensModel, rxMod, output_var,
                      params_list, cores, muj_cache = NULL) {
  pars  <- tryCatch(.admUnpack(p, pinfo), error = function(e) NULL)
  if (is.null(pars)) return(Inf)
  # A NON-FINITE parameter cannot produce a finite objective, but it CAN be handed
  # to rxSolve, which then integrates garbage: the covariance probe legitimately
  # perturbs a sigma to exp(1e5/2) = Inf (cov_h_outer is deliberately huge in the
  # guard tests), and lsoda answers with ~120k `intdy -- t = <denormal> illegal` /
  # `h too small` warnings before the caller's finite-check rejects the result
  # anyway. Rejecting it HERE costs one comparison, removes a guaranteed-useless
  # solve, and keeps the console readable. Written inline rather than as a shared
  # predicate: this runs inside mirai daemons, where assignInNamespace
  # can replace a binding but not ADD one, so a brand-new helper would be missing
  # (see the note at the top of simulate.R).
  if (!(all(is.finite(pars$struct)) && all(is.finite(pars$sigma_var)) &&
        (pinfo$n_eta == 0L || all(is.finite(pars$omega))))) return(Inf)
  total <- 0
  .mkey <- .adfoMuJKey(pars, sensModel)

  for (s_idx in seq_along(studies)) {
    s   <- studies[[s_idx]]

    # Joint (same-subject) unit: stacked FO mean/Jacobian -> full V = J Omega J'
    # (cross blocks from the shared J) + per-output residual; single MVN.
    if (isTRUE(s$is_joint)) {
      mj <- .adfoMuJMemo(muj_cache, s_idx, .mkey, function()
        .adfoGetMuJJoint(pars, pinfo, s, sensModel, rxMod, params_list[[s_idx]], cores))
      V  <- if (pinfo$n_eta > 0L) tcrossprod(mj$J %*% pars$L) else
        matrix(0, s$n_total, s$n_total)
      jr <- .admJointResidual(mj$mu, V, s, pinfo, pars$sigma_var)
      nll_s <- nll_cov_cpp(s$E, s$V, jr$mu, jr$V, s$n)
      if (!is.finite(nll_s)) return(Inf)
      total <- total + nll_s
      next
    }

    ov  <- s$output %||% output_var
    n_t <- length(s$times)
    arr <- .admResidRows(pinfo, ov, pars$sigma_var, n_t)

    mj  <- .adfoMuJMemo(muj_cache, s_idx, .mkey, function()
      .adfoGetMuJ(pars, pinfo, s, sensModel, rxMod, ov,
                  params_list[[s_idx]], cores))
    vp  <- .adfoVpred(mj$mu, mj$J, pars$L, arr, n_t, pinfo$n_eta, s$times)

    nll_s <- if (s$method == "var") {
      nll_var_cpp(s$E, s$v_diag, vp$mu_sigma, diag(vp$V), s$n)
    } else {
      nll_cov_cpp(s$E, s$V, vp$mu_sigma, vp$V, s$n)
    }
    if (!is.finite(nll_s)) return(Inf)
    total <- total + nll_s
  }
  total
}

# -- Gradient ------------------------------------------------------------------

# Gradient of FO NLL w.r.t. optimizer parameter vector p.
#
# Omega/sigma: analytical (chain rule through V_pred = J*Omega*J^T + sigma).
# Struct thetas: forward FD of full NLL -- required because J also depends on
#   theta (V-path: d(J*Omega*J^T)/d(theta) needs second-order sensitivities
#   unavailable from rxode2).
#
# (mu, J) are computed ONCE per study and shared between:
#   - the FD baseline NLL for struct theta FD
#   - the analytical omega/sigma gradient
# This avoids the redundant extra rxSolve that would occur if .adfoNLL() were
# called separately for the baseline and then .admSimulateSens() called again
# for the analytical gradient.
#' @noRd
.adfoGrad <- function(p, pinfo, studies, sensModel, rxMod, output_var,
                       params_list, cores, grad_h = 1e-4) {
  pars <- tryCatch(.admUnpack(p, pinfo), error = function(e) NULL)
  if (is.null(pars)) return(rep(NA_real_, length(p)))
  # Non-finite parameters never reach rxSolve: the covariance probe can perturb
  # a sigma to Inf, and the solver answers a poisoned parameter with tens of
  # thousands of `intdy`/`h too small` warnings before the caller discards the
  # result anyway. Same guard as the NLL entry points -- and it has to be HERE
  # too, because the gradient unpacks `p` itself rather than going through the
  # NLL. Inline, not a shared predicate (mirai daemons cannot gain a new
  # binding -- see the note at the top of simulate.R).
  if (!(all(is.finite(pars$struct)) && all(is.finite(pars$sigma_var)) &&
        (pinfo$n_eta == 0L || all(is.finite(pars$omega))))) return(rep(NA_real_, length(p)))

  n_s   <- length(pinfo$struct_names)
  n_e   <- length(pinfo$sigma_names)
  n_eta <- pinfo$n_eta
  n_o   <- length(pinfo$omega_par)
  grad  <- numeric(length(p))
  names(grad) <- names(p)

  # --- Pass 1: compute (mu, J) and NLL per study at current p ---------------
  # Results cached for struct theta FD baseline and omega/sigma gradient.
  nll_0    <- 0
  muj_cache <- vector("list", length(studies))

  for (s_idx in seq_along(studies)) {
    s   <- studies[[s_idx]]

    # Joint (same-subject) unit: stacked FO mean/Jacobian; V = J Omega J' + per-
    # output residual. Analytical omega/sigma in Pass 3; struct thetas are FD of
    # the joint-aware .adfoNLL in Pass 2 (unchanged).
    if (isTRUE(s$is_joint)) {
      mj <- .adfoGetMuJJoint(pars, pinfo, s, sensModel, rxMod, params_list[[s_idx]], cores)
      JL <- if (n_eta > 0L) mj$J %*% pars$L else matrix(0, s$n_total, 0)
      Vs <- if (n_eta > 0L) tcrossprod(JL) else matrix(0, s$n_total, s$n_total)
      jr <- .admJointResidual(mj$mu, Vs, s, pinfo, pars$sigma_var)
      nll_s <- nll_cov_cpp(s$E, s$V, jr$mu, jr$V, s$n)
      if (!is.finite(nll_s)) { nll_0 <- Inf; break }
      nll_0 <- nll_0 + nll_s
      muj_cache[[s_idx]] <- list(is_joint = TRUE, mu = mj$mu, J = mj$J, JL = JL,
                                 V = jr$V, mu_sigma = jr$mu, n_t = s$n_total,
                                 arr = .admResidRows(pinfo, .admRowOutput(s, s$n_total),
                                                     pars$sigma_var, s$n_total))
      next
    }

    ov  <- s$output %||% output_var
    sel <- .admSigmaSel(pinfo, ov)
    n_t <- length(s$times)
    arr <- .admResidRows(pinfo, ov, pars$sigma_var, n_t)

    mj  <- .adfoGetMuJ(pars, pinfo, s, sensModel, rxMod, ov,
                        params_list[[s_idx]], cores)
    vp  <- .adfoVpred(mj$mu, mj$J, pars$L, arr, n_t, n_eta, s$times)

    nll_s <- if (s$method == "var") {
      nll_var_cpp(s$E, s$v_diag, vp$mu_sigma, diag(vp$V), s$n)
    } else {
      nll_cov_cpp(s$E, s$V, vp$mu_sigma, vp$V, s$n)
    }
    if (!is.finite(nll_s)) { nll_0 <- Inf; break }
    nll_0 <- nll_0 + nll_s
    muj_cache[[s_idx]] <- list(mu = mj$mu, J = mj$J, JL = vp$JL, vp = vp,
                               n_t = n_t, sel = sel, arr = arr,
                               dJ = mj$dJ, dth = mj$dth)
  }

  # Can the structural thetas be differentiated ANALYTICALLY this call?
  #
  # V_pred = J Omega J' + resid(mu, diag(J Omega J')), so d(NLL)/d(theta) needs
  # dJ/d(theta) = d2f/(d eta d theta) -- the cross block an order-2 sens model
  # supplies. Requires every study to have returned it (a joint unit, a
  # transformed endpoint or an order-1 fallback yields NULL), plus the direction
  # map that says which column of dJ belongs to which theta. Anything missing and
  # Pass 2's finite differences run exactly as before -- this is an accuracy
  # upgrade with a fallback, not a new mode.
  use_d2 <- n_s > 0L && is.finite(nll_0) && n_eta > 0L &&
    !is.null(sensModel) && !is.null(sensModel$d2_cols) &&
    !any(vapply(studies, function(u) isTRUE(u$is_joint), logical(1))) &&
    length(muj_cache) == length(studies) &&
    all(vapply(muj_cache, function(mc) !is.null(mc) && !is.null(mc$dJ), logical(1)))
  if (use_d2) {
    # struct theta -> its direction. Paired: its eta's direction (theta and eta
    # enter the parameter identically, so d/d(theta) == d/d(eta)). Unpaired: its
    # own THETA_j_ direction, which also needs its first-order column in dth.
    .dir_of <- character(n_s)
    for (k in seq_len(n_s)) {
      nm <- pinfo$struct_names[k]
      ei <- which(pinfo$struct_eta_idx == k)
      # `[`, NOT `[[`: theta_dirs is a named CHARACTER vector, and `[[` with an
      # unmatched name on an atomic vector THROWS "subscript out of bounds"
      # rather than returning NULL -- so `%||% NA_character_` could never fire
      # and the intended FD degradation was unreachable. .adfoGrad is not wrapped
      # in a tryCatch here, so that error propagated out of eval_grad_f and
      # killed the whole nloptr run. Reachable whenever pinfo's unpaired set and
      # the cached model's differ (theta_dirs is NULL when nothing was unpaired
      # at build time, and `character(0)[["tka"]]` is the same error). Single-
      # bracket indexing yields NA for an unmatched name, which is what the
      # anyNA() check below is written to catch.
      .td <- sensModel$theta_dirs %||% character(0)
      .dir_of[k] <- if (length(ei) > 0L && !is.na(ei[[1L]]))
        sensModel$eta_dirs[ei[[1L]]] %||% NA_character_
      else unname(.td[nm]) %||% NA_character_
    }
    unp <- vapply(seq_len(n_s), function(k) {
      ei <- which(pinfo$struct_eta_idx == k); length(ei) == 0L || is.na(ei[[1L]])
    }, logical(1))
    use_d2 <- !anyNA(.dir_of) &&
      all(vapply(muj_cache, function(mc)
        all(.dir_of %in% names(mc$dJ)) &&
          (!any(unp) || (!is.null(mc$dth) &&
                         all(pinfo$struct_names[unp] %in% names(mc$dth)))),
        logical(1)))
  }

  # --- Pass 2: struct theta forward FD (reuses cached nll_0 as baseline) -------
  #
  # The n_s perturbed configurations differ ONLY in their structural thetas, and
  # every one of them re-solves the same study. They are therefore stacked into a
  # single rxSolve per study (.adfoGetMuJBatch) rather than driven through n_s
  # separate .adfoNLL calls. Same arithmetic, same FD steps -- just one call
  # instead of (1 + n_s) per study.
  #
  # Joint units keep the original per-configuration path: their sensitivity solve
  # is per output block, so batching them needs a per-ID event table (deferred --
  # a subtle bug there yields wrong-but-plausible gradients rather than an error).
  # Skipped entirely when the second-order block is available: the struct thetas
  # are then accumulated analytically in Pass 3, which already holds dNLL/dV.
  if (n_s > 0L && is.finite(nll_0) && !use_d2) {
    hs   <- pmax(abs(p[seq_len(n_s)]), 0.1) * .admGH(grad_h, seq_len(n_s))
    p_pert <- lapply(seq_len(n_s), function(k) { pp <- p; pp[k] <- p[k] + hs[k]; pp })

    any_joint <- any(vapply(studies, function(u) isTRUE(u$is_joint), logical(1)))

    if (any_joint) {
      nll_cfg <- vapply(p_pert, function(pp)
        .adfoNLL(pp, pinfo, studies, sensModel, rxMod, output_var, params_list, cores),
        numeric(1))
    } else {
      struct_mat <- do.call(rbind, lapply(p_pert, function(pp)
        .admUnpack(pp, pinfo)$struct))
      colnames(struct_mat) <- names(pars$struct)

      nll_cfg <- numeric(n_s)
      for (s_idx in seq_along(studies)) {
        s    <- studies[[s_idx]]
        ovs  <- s$output %||% output_var
        n_ts <- length(s$times)
        arrs <- .admResidRows(pinfo, ovs, pars$sigma_var, n_ts)

        mjs <- .adfoGetMuJBatch(struct_mat, pinfo, s, sensModel, rxMod, ovs, cores,
                                pars$sigma_var)
        if (is.null(mjs)) { nll_cfg[] <- Inf; break }

        for (k in seq_len(n_s)) {
          if (!is.finite(nll_cfg[k])) next
          vpk <- .adfoVpred(mjs[[k]]$mu, mjs[[k]]$J, pars$L, arrs, n_ts, n_eta, s$times)
          nll_k <- if (s$method == "var") {
            nll_var_cpp(s$E, s$v_diag, vpk$mu_sigma, diag(vpk$V), s$n)
          } else {
            nll_cov_cpp(s$E, s$V, vpk$mu_sigma, vpk$V, s$n)
          }
          nll_cfg[k] <- if (is.finite(nll_k)) nll_cfg[k] + nll_k else Inf
        }
      }
    }

    grad[seq_len(n_s)] <- (nll_cfg - nll_0) / hs
  }

  # --- Pass 3: analytical omega and sigma gradient (uses cached mu/J/vp) ----
  for (s_idx in seq_along(studies)) {
    mc <- muj_cache[[s_idx]]
    if (is.null(mc)) next
    s      <- studies[[s_idx]]

    # Joint (same-subject) analytical omega/sigma on the stacked V = J Omega J'.
    if (isTRUE(mc$is_joint)) {
      mu_pred <- mc$mu; J <- mc$J; JL <- mc$JL; V_pred <- mc$V
      mu_sigma <- mc$mu_sigma
      r    <- as.numeric(s$E) - mu_sigma
      invV <- tryCatch(chol2inv(chol(V_pred)),
                       error = function(e) tryCatch(solve(V_pred), error = function(e2) NULL))
      if (is.null(invV)) next
      dNLL_dV      <- s$n * (invV - invV %*% (s$V + tcrossprod(r)) %*% invV)
      dNLL_dV_diag <- diag(dNLL_dV)
      # V_pred -> V_struct chain (see the single-output branch below)
      var_f_j  <- diag(tcrossprod(JL))
      .dres_j  <- .admResidDeriv(mu_pred, var_f_j, mc$arr, pinfo)   # once, reused below
      vchain_j <- .admResidVChain(mu_pred, var_f_j, mc$arr, pinfo,
                                  .admRowTimes(s, length(mu_pred)), deriv = .dres_j)
      dNLL_dmu_full <- drop(-2 * s$n * invV %*% r)
      if (n_eta > 0L && n_o > 0L) {
        .bj <- dNLL_dV * vchain_j
        diag(.bj) <- diag(.bj) +           # mean-from-covariance path (TBS only)
          dNLL_dmu_full * (attr(vchain_j, "dmu_dv0") %||% numeric(length(mu_pred)))
        ML <- crossprod(J, .bj %*% JL)
        for (r_idx in seq_along(pinfo$omega_par)) {
          i <- pinfo$chol_i[r_idx]; jj <- pinfo$chol_j[r_idx]
          grad[n_s + n_e + r_idx] <- grad[n_s + n_e + r_idx] +
            if (pinfo$chol_diag[r_idx]) as.numeric(ML[i, i]) * pars$L[i, i]
            else 2 * as.numeric(ML[i, jj])
        }
      }
      grad[n_s + seq_len(n_e)] <- grad[n_s + seq_len(n_e)] +
        .admSigmaGrad(mu_pred, mc$arr, pinfo, dNLL_dV_diag, dNLL_dmu_full, var_f_j,
                      dNLL_dV, .admRowTimes(s, length(mu_pred)), tcrossprod(JL),
                      deriv = .dres_j)
      next
    }

    sel    <- mc$sel
    mu_pred <- mc$mu; J <- mc$J; vp <- mc$vp; n_t <- mc$n_t
    V_pred  <- vp$V;  mu_sigma <- vp$mu_sigma
    r       <- as.numeric(s$E) - mu_sigma
    is_var  <- identical(s$method, "var")

    if (is_var) {
      v_pred       <- diag(V_pred)
      dNLL_dv_pred <-      s$n * (1/v_pred - (s$v_diag + r^2) / v_pred^2)
      dNLL_dV_diag <- dNLL_dv_pred
    } else {
      invV <- tryCatch(
        chol2inv(chol(V_pred)),
        error = function(e) tryCatch(solve(V_pred), error = function(e2) NULL))
      if (is.null(invV)) next
      V_obs_hat    <- s$V + tcrossprod(r)
      dNLL_dV      <- s$n * (invV - invV %*% V_obs_hat %*% invV)
      dNLL_dV_diag <- diag(dNLL_dV)
    }

    # dNLL_dV is w.r.t. V_PRED, but ML below differentiates the STRUCTURAL
    # covariance J Omega J'. They differ by .admResidVChain() whenever the
    # residual is f-dependent (ms_i*ms_j off-diagonal, dv_dv0 on it); identity
    # for purely additive error, so add() models are unchanged.
    var_f  <- vp$var_f
    .dres  <- .admResidDeriv(mu_pred, var_f, mc$arr, pinfo)   # once, reused below
    vchain <- .admResidVChain(mu_pred, var_f, mc$arr, pinfo, s$times, deriv = .dres)
    # Needed by the omega block below as well as the sigma block, so compute it here.
    dNLL_dmu <- if (is_var) -2 * s$n * r / diag(V_pred) else
      drop(-2 * s$n * invV %*% r)

    # + dNLL_dmu * dmu_dv0 on the diagonal: for TBS the predicted MEAN depends on
    # Var_eta(f) = diag(J Omega J'), so omega -- and a structural theta -- reach
    # the objective through the mean too. Zero for every other residual form (see
    # .admResidVChain). Hoisted out of the omega block so that block and the
    # struct-theta block below share ONE value by construction: they are two
    # contractions of the same object and must not drift apart.
    .dmv <- attr(vchain, "dmu_dv0") %||% numeric(length(mu_pred))

    # Omega gradient: d(-2LL)/d(L[i,j]) via ML = t(J) dNLL_dV JL
    # (JL cached from Pass 1 -- avoids materialising M = t(J) dNLL_dV J separately)
    # Diagonal p = log(Omega_ii) = 2*log(L_ii):
    #   d(NLL)/d(L_ii) = 2*(ML)[i,i]  (factor-of-2 from symmetric d(LL^T)/d(L_ii))
    #   d(L_ii)/dp     = L_ii/2        (chain rule through L_ii = exp(p/2))
    #   => d(NLL)/dp   = (ML)[i,i] * L_ii   (the two factors of 2 cancel)
    # Off-diagonal p = L_ij:    chain rule gives 2*(ML)[i,j]
    if (n_eta > 0L && n_o > 0L) {
      JL <- mc$JL
      ML <- if (is_var)
        crossprod(J, JL * (dNLL_dv_pred * diag(vchain) + dNLL_dmu * .dmv))
      else {
        .b <- dNLL_dV * vchain
        diag(.b) <- diag(.b) + dNLL_dmu * .dmv
        crossprod(J, .b %*% JL)
      }
      d  <- pinfo$chol_diag
      for (r_idx in seq_along(pinfo$omega_par)) {
        i  <- pinfo$chol_i[r_idx]; jj <- pinfo$chol_j[r_idx]
        if (d[r_idx]) {
          grad[n_s + n_e + r_idx] <- grad[n_s + n_e + r_idx] +
            as.numeric(ML[i, i]) * pars$L[i, i]
        } else {
          grad[n_s + n_e + r_idx] <- grad[n_s + n_e + r_idx] +
            2 * as.numeric(ML[i, jj])
        }
      }
    }

    # Struct thetas, ANALYTICALLY (order-2 sensitivity model). Replaces Pass 2's
    # forward FD of the whole NLL -- which was adfo's last finite-difference
    # component, and the noisiest thing in the package (differencing a log-det and
    # a quadratic form).
    #
    # theta moves the objective two ways, because V_pred = J Omega J' + resid:
    #   mean:  d(NLL)/d(mu_struct) . d(mu)/d(theta)
    #   cov:   <dNLL/dV_struct, dJ Omega J' + J Omega dJ'>
    #        = 2 * sum(dJ * (B J Omega))          (B and Omega both symmetric)
    # B is the SAME matrix the omega block above contracts (dNLL_dV rotated to the
    # structural covariance by vchain, plus the TBS mean-from-variance term on the
    # diagonal), so the theta and omega paths cannot drift apart -- they are two
    # contractions of one object. d(mu)/d(theta) is the first-order column: the
    # eta's for a mu-referenced theta, its own THETA_j_ column for an unpaired one.
    if (use_d2) {
      .JO <- mc$JL %*% t(pars$L)                         # J Omega
      .G  <- if (is_var) .JO * (dNLL_dv_pred * diag(vchain) + dNLL_dmu * .dmv)
             else {
               .bt <- dNLL_dV * vchain
               diag(.bt) <- diag(.bt) + dNLL_dmu * .dmv
               .bt %*% .JO
             }
      # The residual makes the objective depend on mu beyond the plain r = E - mu
      # path (its variance is mu-dependent, and for TBS so is the mean scale).
      # .admResidMuCoupling returns exactly that EXCESS, so it adds to dNLL_dmu --
      # the same helper admc/adgh/adirmc use, rather than a fourth copy of it.
      .dmu_eff <- dNLL_dmu +
        .admResidMuCoupling(mu_pred, mc$arr, pinfo, dNLL_dV_diag, dNLL_dmu, var_f,
                            if (is_var) NULL else dNLL_dV,
                            if (is_var) NULL else tcrossprod(mc$JL), s$times,
                            deriv = .dres)
      for (k in seq_len(n_s)) {
        .dmu_k <- if (unp[k]) mc$dth[[pinfo$struct_names[k]]]
                  else mc$J[, which(pinfo$struct_eta_idx == k)[[1L]]]
        grad[k] <- grad[k] + sum(.dmu_eff * .dmu_k) +
          2 * sum(mc$dJ[[.dir_of[k]]] * .G)
      }
    }

    # Sigma gradient. Only this output's residual parameters contribute; rows
    # belonging to other endpoints carry a zero derivative, so summing over all
    # rows is a no-op for them.
    # cov_f: the STRUCTURAL covariance J Omega J'. Needed because for lnorm/TBS the
    # mean scale ms is itself a function of a residual parameter, so sigma reaches
    # every OFF-diagonal of V_pred (= ms_i ms_j cov_ij) and not just the diagonal.
    # Under FO the struct thetas are finite-differenced through the full NLL, so
    # only the sigma path needs this; omega's ms factor is already in vchain.
    grad[n_s + seq_len(n_e)] <- grad[n_s + seq_len(n_e)] +
      .admSigmaGrad(mu_pred, mc$arr, pinfo, dNLL_dV_diag, dNLL_dmu, var_f,
                    if (is_var) NULL else dNLL_dV, s$times,
                    if (is_var) NULL else tcrossprod(mc$JL), deriv = .dres)
  }

  grad
}

# Pure finite-difference gradient of FO NLL w.r.t. all optimizer parameters.
# Used for grad = "fd" (forward FD) and grad = "cfd" (central FD).
# No analytical chain rule -- every parameter differentiated by perturbing the NLL.
# The sens model is still used inside .adfoNLL() for efficient J computation.
#' @noRd
.adfoFDGrad <- function(p, pinfo, studies, sensModel, rxMod, output_var,
                         params_list, cores, grad_h = 1e-4, use_central = FALSE) {
  g <- numeric(length(p)); names(g) <- names(p)
  # Shared (mu, J) memo: the omega and sigma directions perturb nothing the solve
  # sees, so they all hit the base entry and cost no rxSolve at all.
  cache <- new.env(parent = emptyenv())
  if (use_central) {
    for (k in seq_along(p)) {
      hk <- pmax(abs(p[k]), 0.1) * .admGH(grad_h, k)
      pp <- p; pp[k] <- p[k] + hk
      pm <- p; pm[k] <- p[k] - hk
      g[k] <- (.adfoNLL(pp, pinfo, studies, sensModel, rxMod, output_var, params_list, cores, cache) -
               .adfoNLL(pm, pinfo, studies, sensModel, rxMod, output_var, params_list, cores, cache)) / (2 * hk)
    }
  } else {
    nll0 <- .adfoNLL(p, pinfo, studies, sensModel, rxMod, output_var, params_list, cores, cache)
    for (k in seq_along(p)) {
      hk <- pmax(abs(p[k]), 0.1) * .admGH(grad_h, k)
      ph <- p; ph[k] <- p[k] + hk
      g[k] <- (.adfoNLL(ph, pinfo, studies, sensModel, rxMod, output_var, params_list, cores, cache) - nll0) / hk
    }
  }
  g
}

# -- Covariance ----------------------------------------------------------------

# Post-fit covariance via numerical Hessian.
# use_grad=TRUE: forward FD of gradient (np_cov+1 grad evals -- faster).
# use_grad=FALSE: full NLL-FD quadratic form (1+2*np_cov+4*n_off NLL evals).
# Hessian over struct + sigma + omega (falls back to struct+sigma if not PD).
.adfoCalcCov <- function(p_hat, pinfo, studies, sensModel, rxMod, output_var,
                          params_list, cores,
                          use_grad = FALSE, grad_h = 1e-4, cov_h = 1e-3,
                          cov_h_outer = .Machine$double.eps^(1/5),
                          gill = FALSE) {
  n_s     <- length(pinfo$struct_names)
  n_e     <- length(pinfo$sigma_names)
  n_o     <- length(pinfo$omega_par)
  # The Hessian spans struct + sigma + OMEGA -- see .adghCalcCov() for the
  # measurement. Excluding omega made the STRUCTURAL SEs too small (reported SE /
  # empirical sampling SD over simulated datasets went from 0.67 to 1.17 on prop
  # and 0.67 to 1.06 on lnorm when omega was put back), because a theta carrying
  # an eta is correlated with that eta's variance. Falls back to the struct+sigma
  # sub-block if the weakly-identified omega Cholesky makes the full H indefinite.
  n_sub   <- n_s + n_e
  cov_idx <- seq_len(n_sub + n_o)
  np_cov  <- length(cov_idx)
  nms_cov <- names(p_hat)[cov_idx]

  # One (mu, J) memo for the whole Hessian. The Hessian perturbs struct AND sigma
  # (omega Cholesky is excluded), and every sigma direction -- including the
  # struct x sigma cross terms -- reuses the solve of its struct component. The
  # solve count therefore collapses to the number of distinct structural vectors,
  # which is far below the 1 + 2*np_cov + 4*n_off NLL evaluations performed.
  cache  <- new.env(parent = emptyenv())
  nll_fn <- function(p)
    suppressMessages(.adfoNLL(p, pinfo, studies, sensModel, rxMod, output_var,
                               params_list, cores, cache))
  grad_fn <- function(p)
    suppressMessages(.adfoGrad(p, pinfo, studies, sensModel, rxMod, output_var,
                                params_list, cores, grad_h = cov_h))

  nll0 <- nll_fn(p_hat)
  if (!is.finite(nll0)) {
    warning("adfoCalcCov: NLL not finite at p_hat -- covariance not computed",
            call. = FALSE)
    return(NULL)
  }

  H <- matrix(0, np_cov, np_cov, dimnames = list(nms_cov, nms_cov))

  if (use_grad) {
    h_fwd <- pmax(abs(p_hat[cov_idx]), 0.1) * cov_h_outer
    g0    <- grad_fn(p_hat)[cov_idx]
    for (jj in seq_len(np_cov)) {
      ph        <- p_hat; ph[cov_idx[jj]] <- ph[cov_idx[jj]] + h_fwd[jj]
      gj        <- grad_fn(ph)[cov_idx]
      H[, jj]   <- if (anyNA(gj)) 0 else (gj - g0) / h_fwd[jj]
    }
    H <- (H + t(H)) / 2
  } else {
    # Step selection. `pmax(abs(p), 0.1) * cov_h_outer` is a guess about how much
    # noise the objective carries, applied identically to every parameter -- and
    # it is the guess behind the "Hessian not positive definite ... try
    # increasing cov_h_outer" warning below. Gill83 measures instead: it probes
    # THIS objective and returns the step where condition error and truncation
    # error balance, per parameter. Exact fit here, since the function it probes
    # is the one being differenced.
    h_gill <- if (isTRUE(gill))
      .admGillSteps(nll_fn, p_hat, cov_idx,
                    fallback = pmax(abs(p_hat[cov_idx]), 0.1) * cov_h_outer,
                    .var.name = "adfoCalcCov")
    else pmax(abs(p_hat[cov_idx]), 0.1) * cov_h_outer
    for (k in seq_len(np_cov)) {
      ki  <- cov_idx[k]; hk  <- h_gill[k]
      p_p <- p_hat; p_p[ki] <- p_p[ki] + hk
      p_m <- p_hat; p_m[ki] <- p_m[ki] - hk
      H[k, k] <- (nll_fn(p_p) - 2 * nll0 + nll_fn(p_m)) / hk^2
    }
    for (i in seq_len(np_cov - 1L)) {
      for (j in seq(i + 1L, np_cov)) {
        ii <- cov_idx[i]; ji <- cov_idx[j]
        hi <- h_gill[i];  hj <- h_gill[j]
        p_pp <- p_hat; p_pp[ii] <- p_pp[ii] + hi; p_pp[ji] <- p_pp[ji] + hj
        p_pm <- p_hat; p_pm[ii] <- p_pm[ii] + hi; p_pm[ji] <- p_pm[ji] - hj
        p_mp <- p_hat; p_mp[ii] <- p_mp[ii] - hi; p_mp[ji] <- p_mp[ji] + hj
        p_mm <- p_hat; p_mm[ii] <- p_mm[ii] - hi; p_mm[ji] <- p_mm[ji] - hj
        H[i, j] <- H[j, i] <-
          (nll_fn(p_pp) - nll_fn(p_pm) - nll_fn(p_mp) + nll_fn(p_mm)) / (4 * hi * hj)
      }
    }
  }

  if (!all(is.finite(H))) {
    warning("adfoCalcCov: Hessian has non-finite entries -- covariance not computed",
            call. = FALSE)
    return(NULL)
  }


  eig_dec <- tryCatch(eigen(H, symmetric = TRUE), error = function(e) NULL)
  H_eigs  <- if (!is.null(eig_dec)) eig_dec$values else rep(NA_real_, np_cov)

  # If the weakly-identified omega Cholesky makes the full Hessian indefinite,
  # drop back to the struct+sigma sub-block rather than reporting nothing.
  .red <- .admReduceNpdOmega(H, H_eigs, eig_dec, nms_cov, n_o, n_sub)
  if (.red$reduced)
    warning("adfoCalcCov: the full Hessian including omega was not positive ",
            "definite; reporting structural and sigma standard errors only.",
            call. = FALSE)
  H <- .red$H; nms_cov <- .red$nms_cov; np_cov <- .red$np_cov
  eig_dec <- .red$eig_dec; H_eigs <- .red$H_eigs

  if (!is.null(eig_dec) && min(H_eigs) < 0) {
    warning(sprintf(
      "adfoCalcCov: Hessian not positive definite (min eigenvalue %.3e). Covariance not computed. Try increasing cov_h_outer (currently %.3e), e.g. cov_h_outer = %.3e.",
      min(H_eigs), cov_h_outer, cov_h_outer * 4), call. = FALSE)
    return(NULL)
  }

  Hinv <- tryCatch(
    chol2inv(chol(H)),
    error = function(e) tryCatch(solve(H), error = function(e2) NULL))
  if (is.null(Hinv)) {
    warning("adfoCalcCov: Hessian inversion failed -- covariance not computed",
            call. = FALSE)
    return(NULL)
  }

  cov_full <- (2 * Hinv + t(2 * Hinv)) / 2
  dimnames(cov_full) <- list(nms_cov, nms_cov)
  # Rotate onto the reported scale (residual delta factors + omega Jacobian). One
  # shared implementation for all three estimators -- see .admScaleReportedCov().
  .admScaleReportedCov(cov_full, p_hat, pinfo, n_s, n_e, n_o, n_sub)
}

# -- Restart worker ------------------------------------------------------------

# Self-contained FO optimization run (one restart); serializable to a worker.
# Signature mirrors .admRestartWorker: same base_args from .admRunRestarts().
# n_sim, sampling accepted for interface compatibility but not used.
# use_pure_fd: TRUE for grad="fd"/"cfd"; dispatches to .adfoFDGrad() instead of .adfoGrad().
.adfoRestartWorker <- function(restart_id, p_init, ui_lstExpr, pinfo,
                                ov_lower, ov_upper, scale_c = NULL, studies, n_sim,
                                seed, algorithm, ftol_rel, maxeval,
                                use_grad, grad_h, grad_bounds,
                                output_var = "cp",
                                sampling = "sobol",
                                use_central = FALSE,
                                use_pure_fd = FALSE,
                                print_progress = TRUE, print = 10L,
                                cores = NULL, no_lock = FALSE,
                                sens_cache_file = NULL, sens_cols = NULL,
                                sens_rename = NULL,
                                rxMod_direct = NULL, sensModel_direct = NULL) {
  library(admixr2)

  # Dev mode: patch the installed namespace with any dev functions in .GlobalEnv.
  # A daemon is patched by .admDaemonRestart() before it gets here; this covers a
  # direct (sequential) call. tryCatch guards against the installed package
  # predating this function (run devtools::install() once).
  tryCatch(.admPatchDevNamespace(), error = function(e) NULL)

  m <- .admWorkerLoadModels(ui_lstExpr, rxMod_direct, cores,
                            sens_cache_file, sens_cols, sens_rename, sensModel_direct,
                            pinfo)

  params_list <- .admMakeParamsList(1L, pinfo, length(studies))

  set.seed(seed + restart_id)

  nll_fn <- function(p)
    .adfoNLL(p, pinfo, studies, m$sensModel, m$rxMod, output_var, params_list, m$cores_w)

  grad_fn <- if (use_pure_fd) {
    function(p) .adfoFDGrad(p, pinfo, studies, m$sensModel, m$rxMod, output_var,
                            params_list, m$cores_w, grad_h, use_central)
  } else {
    function(p) .adfoGrad(p, pinfo, studies, m$sensModel, m$rxMod, output_var,
                          params_list, m$cores_w, grad_h)
  }

  # adfo loads its own model in-process and does not lock (single-nloptr path).
  .admScaledOptimize(restart_id, p_init, ov_lower, ov_upper, scale_c,
                     use_grad, grad_bounds, algorithm, ftol_rel, maxeval,
                     nll_fn, grad_fn, pinfo, print_progress, print,
                     lock_rxMod = NULL)
}

# -- Control object ------------------------------------------------------------

#' Control settings for the FO (First-Order) estimator
#'
#' Creates a control object for `nlmixr2(est = "adfo")`.  The FO estimator
#' linearises model predictions at \eqn{\eta = 0}: it is faster than the MC
#' estimator but less accurate for models with large IIV or strongly
#' non-linear individual predictions.
#'
#' @param studies Named list of study specifications (same format as
#'   [admControl()]: `E`, `V`, `n`, `times`, `ev`, optional `method`; or an
#'   `observations` list for multi-compartment fits -- see [admControl()]).
#' @param resid_nodes Gauss-Hermite nodes used to integrate the RESIDUAL for a
#'   transform-both-sides endpoint (`boxCox`, `yeoJohnson`, `logitNorm`,
#'   `probitNorm`), where `y = g(h(f) + sigma*eps)` has no closed-form mean and
#'   variance. Ignored by every other error model, which has closed forms. Default
#'   81. Measured worst-case relative error against an independent quadrature, over
#'   all four transforms and residual SD of 0.5, 1, 2 and 3: n = 15 gives 5.7e-2,
#'   31 gives 4.5e-3, 81 gives 5.0e-5. The error is dominated by large residual SD;
#'   at SD <= 1, n = 31 already gives 1e-7 or better.
#'
#'   This is an ACCURACY dial, not a speed one. The quadrature is linear in
#'   `resid_nodes` in isolation (~50 us at 15, 300 us at 81 for an 8-row study) but
#'   negligible beside the ODE solve: a full NLL evaluation measured 0.750 s per 60
#'   evaluations at BOTH 31 and 81 nodes. Raise it if you have a saturating endpoint
#'   with a large residual SD; there is little to gain by lowering it.
#' @param grad Gradient mode.  `"analytical"` (default) uses the closed-form FO
#'   gradient with LBFGS; `"none"` uses derivative-free BOBYQA; `"fd"` uses
#'   forward finite differences of the full NLL; `"cfd"` uses central finite
#'   differences for struct theta gradient (more accurate than `"fd"`, roughly
#'   twice as many NLL evaluations per step).
#'
#'   The default was `"none"` up to 0.4.0, because the structural thetas were
#'   finite-differenced through the whole NLL and the resulting gradient was too
#'   noisy for a quasi-Newton step to pay off. They are now differentiated
#'   analytically from a second-order sensitivity model (relative error ~1e-7
#'   against a central difference, where the finite-difference pass reached
#'   1e-2), so LBFGS on the exact gradient is the better default. A model that
#'   cannot build that sensitivity model falls back to the finite-difference
#'   gradient automatically, and `grad = "none"` remains available.
#' @param algorithm nloptr algorithm, or `NULL` (default) to pick the default
#'   that matches `grad`: `"NLOPT_LD_LBFGS"` with a gradient, `"NLOPT_LN_BOBYQA"`
#'   when `grad = "none"`. Any algorithm reported by
#'   [nloptr::nloptr.print.options()] is accepted. An explicit algorithm is
#'   reconciled with `grad`: when `grad = "none"` a gradient-based algorithm
#'   (`NLOPT_LD_*` / `NLOPT_GD_*`) falls back to `"NLOPT_LN_BOBYQA"`; when a
#'   gradient is requested a derivative-free algorithm (`NLOPT_LN_*` /
#'   `NLOPT_GN_*`) turns the gradient off. Both emit a message.
#' @param maxeval Maximum function evaluations (default 500).
#' @param ftol_rel Relative tolerance (default `sqrt(.Machine$double.eps)`).
#' @param print Print-frequency for live progress (0 = silent).
#' @param seed Random seed (used for restarts).
#' @param cores OpenMP threads for `rxSolve()`. Defaults to
#'   `rxode2::rxCores()`. When `workers > 1` it is a *total* budget, split
#'   across the workers.
#' @param nDisplayProgress Passed to `rxSolve()`: show the solver's text
#'   progress bar only once a single solve exceeds this many subjects. The
#'   default (`.Machine$integer.max`) keeps it off for clean script/vignette
#'   output; lower it (e.g. `1000L`) to see progress during long fits.
#' @param grad_h Finite-difference step for unpaired struct theta gradient and
#'   FD Jacobian.
#' @param grad_bounds Box-constraint half-width when using gradients: the fit is
#'   confined to `p0 +/- grad_bounds` on the optimizer scale, which for a
#'   log-scale parameter is a factor of `exp(grad_bounds)` (~148 at the default
#'   5). This bound is admixr2's, not the model's -- an unbounded parameter has
#'   no other -- and nloptr reports normal convergence at a box corner, so a
#'   warning is emitted if an estimate finishes on it.
#' @param cov_h_outer Outer step scale for NLL-FD Hessian.
#' @param gill Use the Gill, Murray, Saunders & Wright (1983) algorithm to choose
#'   the finite-difference steps, instead of the fixed `cov_h_outer` and `grad_h`
#'   scales. `FALSE` (the default) keeps the existing steps, so nothing moves
#'   unless you ask.
#'
#'   A fixed scale is one guess about how much noise the objective carries,
#'   applied identically to every parameter -- and it is the guess behind the
#'   "Hessian not positive definite ... try increasing `cov_h_outer`" warning.
#'   Gill83 probes the objective instead and returns, per parameter, the step at
#'   which condition error and truncation error balance, so a parameter the
#'   objective is flat in and one it is sharp in get different steps and the
#'   choice tracks the ODE tolerance rather than assuming one. It is the same
#'   algorithm FOCEI uses for its own steps (`nlmixr2est::nlmixr2Gill83()`), so
#'   the two packages agree on what a finite difference is. Worth trying first
#'   when the Hessian comes back indefinite.
#'
#'   Applied in both of the places adfo finite-differences the OBJECTIVE:
#'   the post-fit covariance Hessian, and the optimizer's gradient -- every
#'   parameter under `grad = "fd"`/`"cfd"`, or the structural thetas when the
#'   second-order sensitivity model is unavailable and they are differenced
#'   through the whole NLL. When that model IS available no structural theta is
#'   finite-differenced at all and the gradient probe is skipped.
#'
#'   Measured ONCE in each case -- for the gradient at the starting values, which
#'   is where FOCEI measures too (`numericGrad` runs Gill83 at its first
#'   evaluation and reuses the steps for the rest of the fit). Costs up to ten
#'   extra objective evaluations per finite-differenced parameter. A parameter
#'   Gill83 reports as flat or linear falls back to the fixed scale.
#'
#'   Because it reaches the gradient, `gill = TRUE` can change the FIT and not
#'   only the standard errors -- but only when the gradient has a
#'   finite-differenced parameter to begin with.
#' @param covMethod `"r"` computes covariance via a numerical Hessian over the
#'   structural, residual-error and omega parameters; `"none"` skips it. Omega is
#'   included because excluding it also biases the STRUCTURAL standard errors
#'   downward -- a theta carrying an eta is correlated with that eta's variance.
#'   If the weakly-identified omega Cholesky makes the Hessian non-positive
#'   definite, the structural + residual sub-block is reported with a warning.
#'
#'   All three blocks are reported on the scale the ESTIMATES are printed on, as
#'   `nlmixr2est` does: structural thetas on the log/optimizer scale, residual
#'   error as an SD, and omega as the variance/covariance entries (named
#'   `om.<eta>` and `cov.<eta_i>.<eta_j>`). The omega block is rotated by the
#'   full Jacobian of Omega with respect to the log-Cholesky, which is not
#'   diagonal once omega is correlated.
#'
#'   **An adfo standard error describes scatter, not accuracy.** FO linearises the
#'   model at eta = 0, and on a non-additive residual (or a saturating endpoint,
#'   or a large omega) the resulting point estimates carry a bias of several
#'   standard errors -- measured 5-20 SE, giving 0% coverage for a nominal 95%
#'   interval even where the SE itself matches the sampling SD. Use `adgh` or
#'   `admc` when the uncertainty matters.
#' @param n_restarts Number of optimizer restarts (1 = no multi-start).
#' @param restart_sd Standard deviation for random perturbations of initial
#'   struct thetas at each restart (> 1).
#' @param workers Number of parallel workers (mirai daemons) for multi-restart
#'   (default 1 = sequential). Requires the `mirai` package.
#' @param cov_h Inner FD step for the gradient-based Hessian (only used when
#'   `covMethod = "r"` and `grad != "none"`). Default 1e-3.
#' @param rxControl `rxode2::rxControl()` object. Created automatically when `NULL`.
#' @param sigdig Significant digits asked of the ODE solver, or `NULL` (the
#'   default) to leave rxode2's own solver tolerances alone. When set, it is
#'   passed to `rxode2::rxSolve()`'s own `sigdig` argument for every solve the
#'   estimator issues -- rxode2 owns the mapping to `atol`/`rtol` and has changed
#'   it between releases, which is why the digits, not the tolerances, are what
#'   travels -- and to `nlmixr2est::foceiControl()` for the post-fit tables.
#'
#'   It is a speed lever, and an opt-in one because it is not free. The
#'   estimators finite-difference the solve with steps of the same order:
#'   `grad_h` (1e-4), `cov_h` (1e-3) and `cov_h_outer` (~2.5e-3), while
#'   `sigdig = 4` maps to a relative tolerance of ~1e-4 on current rxode2.
#'   Differencing a solution whose own noise is 1e-4 with a 1e-4 step returns
#'   noise, and it surfaces as a moved objective and an indefinite covariance
#'   Hessian (every `SE` reported `NA`) rather than as an error. Most worthwhile
#'   where the gradient is fully analytic and nothing differences the solve --
#'   `adfoControl(grad = "analytical")` measured ~4.8x faster at `sigdig = 4`
#'   with standard errors unchanged to 4 significant figures. Elsewhere, compare
#'   the objective and the standard errors against `NULL` before relying on it.
#'   Table formatting is unaffected either way: `sigdigTable` defaults to 4
#'   regardless.
#' @param calcTables,compress,ci,sigdigTable,optExpression,sumProd,literalFix
#'   Passed to `nlmixr2est::foceiControl()` for the table/output machinery.
#' @param addProp How combined additive+proportional error is parameterised in
#'   the nlmixr2 output tables: `"combined2"` (default, variance form) or
#'   `"combined1"` (SD form). Has no effect on admixr2's own estimation.
#' @param returnAdmr If `TRUE`, return a plain list instead of the full
#'   nlmixr2 fit object.
#' @param ... Unused arguments (trigger an error).
#'
#' @return An `adfoControl` object (a named list).
#'
#' @section Installing memuse:
#' `rxode2::rxSolve()` estimates free RAM on every call. When the `memuse`
#' package is not installed its fallback ends up shelling out to `vm_stat`, a
#' macOS-only command, so on Windows and Linux every solve spawns a process that
#' can only fail. Because the FO estimator issues many small solves, this
#' overhead is measurable (roughly 17% of an FO gradient). Installing `memuse`
#' makes the fallback unreachable:
#' \preformatted{install.packages("memuse")}
#'
#' @seealso [admControl()], [adirmcControl()]
#'
#' @examples
#' # Inspect defaults
#' ctl <- adfoControl()
#' ctl$grad
#' ctl$maxeval
#'
#' # Analytical gradient, more evaluations
#' ctl2 <- adfoControl(grad = "analytical", maxeval = 1000L)
#'
#' \donttest{
#' library(rxode2)
#' library(nlmixr2)
#'
#' data("examplomycin")
#' obs    <- examplomycin[examplomycin$EVID == 0, ]
#' obs    <- obs[order(obs$ID, obs$TIME), ]
#' times  <- sort(unique(obs$TIME))
#' ids    <- unique(obs$ID)
#' dv_mat <- do.call(rbind, lapply(ids, function(i) {
#'   sub <- obs[obs$ID == i, ]; sub$DV[order(sub$TIME)]
#' }))
#' E <- colMeans(dv_mat)
#' V <- cov.wt(dv_mat, method = "ML")$cov
#'
#' pk_model <- function() {
#'   ini({
#'     tcl <- log(5); tv <- log(30)
#'     prop.sd <- c(0, 0.2)
#'     eta.cl ~ 0.09; eta.v ~ 0.04
#'   })
#'   model({
#'     cl <- exp(tcl + eta.cl)
#'     v  <- exp(tv  + eta.v)
#'     d/dt(central) <- -(cl/v) * central
#'     cp <- central / v
#'     cp ~ prop(prop.sd)
#'   })
#' }
#'
#' fit <- nlmixr2(
#'   pk_model, admData(), est = "adfo",
#'   control = adfoControl(
#'     studies = list(study1 = list(E = E, V = V, n = length(ids),
#'                                  times = times, ev = et(amt = 100))),
#'     maxeval = 100L
#'   )
#' )
#' print(fit)
#' }
#'
#' @export
adfoControl <- function(
    studies    = list(),
    grad        = c("analytical", "none", "fd", "cfd"),
    algorithm  = NULL,
    maxeval    = 500L,
    ftol_rel   = .Machine$double.eps^(1/2),
    print      = 10L,
    seed       = 12345L,
    cores      = rxode2::rxCores(),
    nDisplayProgress = .Machine$integer.max,
    grad_h      = 1e-4,
    grad_bounds = 5,
    cov_h       = 1e-3,
    cov_h_outer = .Machine$double.eps^(1/5),
    covMethod   = c("r", "none"),
    n_restarts  = 1L,
    restart_sd  = 0.5,
    workers     = 1L,
    rxControl     = NULL,
    calcTables    = FALSE,
    compress      = TRUE,
    ci            = 0.95,
    sigdig        = NULL,
    sigdigTable   = NULL,
    addProp       = c("combined2", "combined1"),
    optExpression = TRUE,
    sumProd       = FALSE,
    literalFix    = TRUE,
    returnAdmr    = FALSE,
    # LAST on purpose: inserting an argument mid-signature silently rebinds every
    # positional call -- adfoControl(studies, "fd") used to set grad = "fd".
    resid_nodes = 81L,
    # ... and this one after it, for the same reason.
    gill        = FALSE,
    ...) {

  .xtra <- list(...)
  if (length(.xtra) > 0)
    stop("adfoControl: unused argument(s): ",
         paste(paste0("'", names(.xtra), "'"), collapse = ", "), call. = FALSE)

  addProp  <- match.arg(addProp)
  # Whether the user NAMED `grad`, recorded before match.arg() erases the
  # distinction. It decides how loudly the driver reports a fall-back to finite
  # differences: an unavailable sensitivity model is routine for a defaulted
  # grad (a fixed-effects-only model, an ordinal endpoint) and merits a message,
  # but silently ignoring an EXPLICIT grad = "analytical" does not -- a message is
  # swallowed by suppressMessages(), a knitr chunk with message = FALSE, or any
  # stderr-capturing wrapper, leaving no durable record that the run used a
  # gradient the control asked it not to use. A warning survives to warnings() and
  # to options(warn = 2).
  .grad_explicit <- !missing(grad)
  grad     <- match.arg(grad)
  covMethod <- match.arg(covMethod)

  checkmate::assertList(studies)
  # A residual quadrature needs a real grid. .adghNodes1() refuses m < 1, but it
  # accepts 1..4 happily and returns a rule that integrates nothing usefully --
  # the measured error at 5 nodes is already 3.3e-1. Refuse here, where the
  # message can name the argument, rather than silently scoring a wrong NLL.
  checkmate::assertIntegerish(resid_nodes, lower = 5L, len = 1)
  checkmate::assertIntegerish(maxeval,    lower = 1L, len = 1)
  checkmate::assertNumeric(ftol_rel,      lower = 0,  len = 1)
  checkmate::assertIntegerish(print,      lower = 0L, len = 1)
  checkmate::assertIntegerish(seed,                   len = 1)
  checkmate::assertIntegerish(cores,      lower = 1L, len = 1)
  checkmate::assertIntegerish(nDisplayProgress, lower = 1L, len = 1,
                              .var.name = "nDisplayProgress")
  checkmate::assertNumeric(grad_h,        lower = 0,  len = 1)
  checkmate::assertNumeric(grad_bounds,   lower = 0,  len = 1)
  checkmate::assertNumeric(cov_h,         lower = 0,  len = 1)
  checkmate::assertNumeric(cov_h_outer,   lower = 0,  len = 1)
  checkmate::assertFlag(gill, .var.name = "gill")
  checkmate::assertIntegerish(n_restarts, lower = 1L, len = 1)
  checkmate::assertNumeric(restart_sd,    lower = 0,  len = 1)
  checkmate::assertIntegerish(workers,    lower = 1L, len = 1)
  checkmate::assertNumeric(ci, lower = 0, upper = 1,  len = 1)
  if (!is.null(sigdig))
  checkmate::assertIntegerish(sigdig,     lower = 1L, len = 1)
  checkmate::assertLogical(returnAdmr,                len = 1)

  .algo     <- .admResolveAlgorithm(algorithm, grad,
                                    .var.name = "adfoControl: algorithm")
  algorithm <- .algo$algorithm
  grad      <- .algo$grad

  # sigdig = NULL (the DEFAULT) means "leave rxode2's own solver defaults alone".
  # It is the one setting whose meaning does not move under an rxode2 upgrade,
  # and it is the default because a looser solve is not free: this release is
  # what first routed sigdig into the estimators' own rxSolve calls, and every
  # finite-difference step that consumes those solves (grad_h 1e-4, cov_h 1e-3,
  # cov_h_outer ~2.5e-3) is the same order as the tolerance sigdig = 4 asks for.
  # rxode2 5.1.5 maps sigdig = 4 to rtol = 1e-4 (5.1.4 mapped it to 5e-7 -- 200x
  # tighter for the same request), so differencing with a 1e-4 step differences
  # noise: a moved objective and an indefinite Hessian, not an error. Shipping it
  # on by default would have changed the numerics of every existing script
  # silently, for a knob that looked like table formatting before this release.
  #
  # NULL is also the only way back: the sigdig -> tolerance map is
  # one-dimensional while rxode2's defaults are not (atol 1e-8 vs rtol 1e-6), so
  # no sigdig value reproduces them. The tables still need a number, so they fall
  # back to 4 -- i.e. sigdigTable is unchanged whichever way sigdig is set.
  if (is.null(rxControl))   rxControl   <- if (is.null(sigdig))
    rxode2::rxControl() else rxode2::rxControl(sigdig = sigdig)
  if (is.null(sigdigTable)) sigdigTable <- if (is.null(sigdig)) 4L else
    max(round(sigdig), 3L)

  .ret <- list(
    studies       = studies,
    resid_nodes   = as.integer(resid_nodes),
    n_sim         = 1L,        # interface compatibility with .admRunRestarts()
    sampling      = "sobol",   # idem
    grad          = grad,
    algorithm     = algorithm,
    maxeval       = as.integer(maxeval),
    ftol_rel      = ftol_rel,
    print         = as.integer(print),
    seed          = as.integer(seed),
    cores         = as.integer(cores),
    nDisplayProgress = as.integer(nDisplayProgress),
    grad_h        = grad_h,
    grad_bounds   = grad_bounds,
    cov_h         = cov_h,
    cov_h_outer   = cov_h_outer,
    gill          = gill,
    grad_explicit = .grad_explicit,
    covMethod     = covMethod,
    n_restarts    = as.integer(n_restarts),
    restart_sd    = restart_sd,
    workers       = as.integer(workers),
    rxControl     = rxControl,
    calcTables    = calcTables,
    compress      = compress,
    ci            = ci,
    sigdig        = sigdig,
    sigdigTable   = as.integer(sigdigTable),
    addProp       = addProp,
    optExpression = optExpression,
    sumProd       = sumProd,
    literalFix    = literalFix,
    returnAdmr    = returnAdmr
  )
  class(.ret) <- "adfoControl"
  .ret
}

# -- nlmixr2 S3 hooks ----------------------------------------------------------

#' @noRd
getValidNlmixrCtl.adfo <- function(control) {
  if (inherits(control, "adfoControl")) return(control)
  .ctl <- control[[1]]
  if (inherits(.ctl, "adfoControl")) return(.ctl)
  if (is.list(.ctl) && "studies" %in% names(.ctl))
    return(do.call(adfoControl, .ctl[intersect(names(.ctl), names(formals(adfoControl)))]))
  if (is.list(control) && length(names(control)) > 0L)
    return(do.call(adfoControl, control[intersect(names(control), names(formals(adfoControl)))]))
  adfoControl()
}

#' @noRd
nmObjHandleControlObject.adfoControl <- function(control, env) {
  assign("adfoControl", control, envir = env)
}

#' @noRd
nmObjGetControl.adfo <- function(x, ...) {
  .env <- x[[1]]
  for (.nm in c("adfoControl", "control")) {
    if (exists(.nm, .env)) {
      .ctl <- get(.nm, .env)
      if (inherits(.ctl, "adfoControl")) return(.ctl)
    }
  }
  stop("cannot find adfo control object", call. = FALSE)
}

# -- Main estimation entry point -----------------------------------------------

#' Fit an aggregate data model via First-Order (FO) approximation
#'
#' Called automatically by `nlmixr2(model, admData(), est = "adfo",
#' control = adfoControl(...))`.  Not typically called directly.
#'
#' @param env nlmixr2 environment containing `ui` and `control`.
#' @param ... Unused.
#'
#' @return An `admFit` nlmixr2 fit object.
#'
#' @method nlmixr2Est adfo
#' @importFrom nlmixr2est nlmixr2Est
#' @export
nlmixr2Est.adfo <- function(env, ...) {
  .ui  <- env$ui
  .ctl <- env$control

  if (!inherits(.ctl, "adfoControl")) .ctl <- getValidNlmixrCtl.adfo(.ctl)
  if (!inherits(.ctl, "adfoControl"))
    stop("Could not recover adfoControl", call. = FALSE)
  assign("control", .ctl, envir = .ui)

  studies <- .ctl$studies
  if (length(studies) == 0L)
    stop("adfoControl(studies=...) required", call. = FALSE)
  if (is.null(names(studies)))
    names(studies) <- paste0("study", seq_along(studies))

  pinfo      <- .admDriverPinfo(.ui, .ctl)
  output_var <- .admOutputVar(.ui)
  # A beta endpoint's precision phi is SOLVED, not fitted: .admSimulate() returns
  # it as an attribute on cp_mat and admc/adgh patch it into the residual rows.
  # This estimator has no such path -- arr$phi stayed NA and .admResidApply()
  # produced NA for every row, i.e. an NaN objective with no explanation. Refuse
  # with one instead, and name the estimators that do support it.
  if (any(vapply(.admResidSpecs(pinfo), function(sp)
        identical(sp$form, .ADM_RESID_BETA), logical(1))))
    stop("est='adfo' does not support a beta() endpoint: the beta precision is derived ",
         "from the solved shapes, which only the admc and adgh paths recover.
",
         "  Use est = \"admc\" or est = \"adgh\" for this model.", call. = FALSE)


  # ev_full = FALSE: the three checks below run BETWEEN the flatten and the
  # ev_full build, and must keep doing so. They inspect per-unit fields
  # (`method`, `is_joint`, `blocks`) that only exist on flattened units -- checking
  # the un-flattened studies made every ordinal fit fail its own block-count
  # check -- and adgh/admc call them at the equivalent point.
  .u         <- .admDriverUnits(studies, .ui, output_var, ev_full = FALSE)
  studies    <- .u$studies
  multi_out  <- .u$multi_out
  any_joint  <- .u$any_joint
  .admCheckAR(pinfo, studies)
  .admCheckOrdinal(pinfo, studies)
  .admCheckMixedEndpoints(.ui)
  studies    <- .admBuildEvFull(studies, tag_cmt = multi_out)

  want_grad   <- .ctl$grad != "none"
  want_sens   <- .ctl$grad == "analytical"
  use_central <- .ctl$grad == "cfd"
  use_pure_fd <- .ctl$grad %in% c("fd", "cfd")
  # Joint (same-subject) fits keep the analytical FO gradient: .adfoGrad's joint
  # branch applies the omega/sigma chain rule to the stacked V = J Omega J' + res
  # (struct thetas via FD of the joint-aware .adfoNLL, as in the single-output
  # path). grad = "fd"/"cfd" still use the pure FD gradient.
  # Multi-compartment fits use the same per-unit gradient as single-output: the
  # analytical omega/sigma chain rule and struct-theta FD apply per observed
  # output (each is an independent block with its own residual error).


  # ORDERING INVARIANT: .admLoadSensModel() must run before .admLoadModel().
  #
  # order = 2: adfo is the ONE estimator that needs the cross second-order block
  # d2f/(d eta d theta), because its V_pred = J Omega J' depends on theta through
  # J. With it .adfoGrad differentiates the struct thetas analytically instead of
  # finite-differencing the whole NLL (measured against a central difference on a
  # 1-cmt oral model: 2e-07..2e-06 relative, against 8e-04..1e-02 for the FD pass
  # it replaces). Only asked for under grad = "analytical" -- the default
  # grad = "none" runs BOBYQA and would pay the extra compartments for nothing.
  # A model that cannot build it (transformed endpoint, joint unit, an rxode2
  # that refuses) transparently gets the order-1 model back and keeps the FD pass.
  sensModel <- if (want_sens) {
    # order 1 for a JOINT (same-subject) fit: `have_d2` below excludes joint
    # units -- Pass 3 is not implemented for them -- so asking for order 2 would
    # compile the cross block and integrate it on every solve for a result
    # nothing reads. Same waste that .admLoadSensModel() already avoids for a
    # transformed endpoint; this is the other case that reaches it.
    .ord <- if (any_joint) 1L else 2L
    sm <- tryCatch(.admLoadSensModel(.ui, order = .ord), error = function(e) NULL)
    if (is.null(sm)) {
      # How loudly depends on whether the user ASKED for this gradient.
      #
      # Defaulted: a message. NULL here is very often BY DESIGN, not a failure --
      # .admLoadSensModel returns it for a fixed-effects-only model (n_eta == 0),
      # an ordinal endpoint, and mixed transformed/untransformed endpoints. Those
      # models ran BOBYQA silently under the old grad = "none" default, and making
      # "analytical" the default turned each of them into a warning on every fit,
      # for a condition the user cannot act on. The fit is correct either way,
      # just finite-differenced, and grad_label below already says "FD".
      #
      # Explicit: a warning. Quietly giving someone who wrote grad = "analytical"
      # the gradient they wrote it to avoid is a different matter, and a message
      # leaves no record of it -- see .grad_explicit in adfoControl().
      .msg <- "adfo: no sensitivity model for this model -- gradient by forward FD."
      if (isTRUE(.ctl$grad_explicit))
        warning(sprintf("adfoControl(grad = '%s'): %s", .ctl$grad,
                        sub("^adfo: ", "", .msg)), call. = FALSE)
      else message(.msg)
      # ... and actually fall back. This only warned: want_sens stayed TRUE, so
      # .adfoGrad's ANALYTICAL path still ran on an FD Jacobian. adgh.R does this
      # correctly; adfo did not.
      want_sens <- FALSE
    }
    sm
  } else NULL

  # Whether the struct thetas are differentiated analytically at all. The order-2
  # cross block is what .adfoGrad's Pass 3 contracts; without it Pass 2's forward
  # FD of the whole NLL runs instead, for EVERY struct theta (paired and unpaired
  # alike -- adfo's struct gradient has always been the FD one). Also FALSE for a
  # joint fit, whose Pass 3 is not implemented.
  have_d2 <- want_sens && !is.null(sensModel$d2_cols) && !any_joint

  # Say so whenever they really are finite-differenced.
  #
  # The gate used to also require any(!pinfo$struct_has_eta), which meant a model
  # whose thetas are ALL mu-referenced -- the common cl <- exp(tcl + eta.cl)
  # style -- got no notice of any kind when the order-2 build silently fell back
  # to order 1 (a linCmt promotion that fails, an indLin bail-out, .rxSens
  # returning nothing). The returned object is a perfectly valid order-1 sens
  # model, so the is.null(sm) branch above does not fire either: the fit ran
  # LBFGS on the 8e-04..1e-02 forward FD the analytic pass exists to replace,
  # invisibly, surfacing only as slow or stalled convergence.
  if (pinfo$n_eta > 0L && want_sens && !have_d2) {
    message(if (any_joint)
      "adfo: joint (same-subject) study -- struct thetas by forward FD."
      else
      "adfo: second-order sensitivities unavailable -- struct thetas by forward FD.")
  }

  rxMod <- .admLoadModel(.ui)
  rxode2::rxLock(rxMod)
  # Free the models this fit registered with rxode2's own idiom (the same
  # gc(); rxUnloadAll() nlmixr2est runs), so many fits in a session stay bounded.
  on.exit({ rxode2::rxUnlock(rxMod); rxode2::rxSolveFree(); gc(FALSE); rxode2::rxUnloadAll() },
          add = TRUE)

  params_list <- .admMakeParamsList(1L, pinfo, length(studies))

  ov    <- .admBuildOptVec(pinfo)
  .iter <- 0L
  cores <- .ctl$cores

  .nll_trace <- numeric(0)
  .par_trace <- NULL
  .best_nll  <- Inf

  eval_f <- function(p) {
    .iter <<- .iter + 1L
    val <- .adfoNLL(p, pinfo, studies, sensModel, rxMod, output_var, params_list, cores)
    if (is.finite(val) && val < .best_nll) {
      .best_nll  <<- val
      .nll_trace <<- c(.nll_trace, val)
      .par_trace <<- rbind(.par_trace, p)
    }
    if (.ctl$print > 0L && .iter %% .ctl$print == 0L) {
      row <- .admProgressRow(sprintf("%04d", .iter), val, p, pinfo)
      if (!is.null(row)) message(row)
    }
    val
  }

  # gill = TRUE: measure the gradient's FD steps ONCE, here, and let every later
  # difference reuse them (the mechanism FOCEI's numericGrad uses at nF == 1).
  # Only for the parameters that are ACTUALLY finite-differenced -- under
  # grad = "fd"/"cfd" that is all of them, otherwise just the struct thetas of
  # Pass 2, and when Pass 3 has the order-2 block there are none at all and the
  # probe is skipped entirely. `.ctl$grad_h` stays a scalar unless this fires, so
  # gill = FALSE is bit-for-bit unchanged.
  .fd_idx <- if (!want_grad) integer(0)
    else if (use_pure_fd) seq_along(ov$p0)
    else if (!have_d2) seq_len(length(pinfo$struct_names))
    else integer(0)
  grad_h_v <- if (isTRUE(.ctl$gill) && length(.fd_idx))
    .admGillGradH(function(p) .adfoNLL(p, pinfo, studies, sensModel, rxMod,
                                       output_var, params_list, cores),
                  ov$p0, .fd_idx, .ctl$grad_h, scaled = TRUE,
                  .var.name = "adfo gradient")
  else .ctl$grad_h

  eval_grad_f <- if (!want_grad) {
    NULL
  } else if (use_pure_fd) {
    function(p) .adfoFDGrad(p, pinfo, studies, sensModel, rxMod, output_var,
                              params_list, cores, grad_h_v, use_central)
  } else {
    function(p) .adfoGrad(p, pinfo, studies, sensModel, rxMod, output_var,
                           params_list, cores, grad_h_v)
  }

  # "Analytical" is claimed only when the struct thetas really are analytic. With
  # an order-1 sens model the omega/sigma blocks are analytic but every struct
  # theta is finite-differenced, which is a materially different gradient to be
  # told you are running.
  grad_label <- if (!want_grad) "none"
    else if (have_d2) "Analytical"
    else if (!is.null(sensModel)) "Analytical (struct FD)"
    else if (use_central) "CFD" else "FD"
  message("=== admixr2: Aggregate Data Modeling (FO) ===")
  message(sprintf("  Obs units: %d | Params: %d | Cores: %d | Grad: %s | Restarts: %d",
                  length(studies), length(ov$p0), cores, grad_label, .ctl$n_restarts))
  t0 <- proc.time()

  lb <- if (want_grad) pmax(ov$lower, ov$p0 - .ctl$grad_bounds) else ov$lower
  ub <- if (want_grad) pmin(ov$upper, ov$p0 + .ctl$grad_bounds) else ov$upper

  sc           <- ov$scale_c
  p0_sc        <- ov$p0 / sc
  lb_sc        <- lb  / sc
  ub_sc        <- ub  / sc
  eval_f_sc    <- function(p_s) eval_f(p_s * sc)
  eval_grad_sc <- if (!is.null(eval_grad_f)) {
    function(p_s) eval_grad_f(p_s * sc) * sc
  } else NULL

  if (.ctl$n_restarts == 1L) {
    message(.admProgressHeader(pinfo))
    opt_raw <- nlmixr2est::nlmixrWithTiming("adfo", {
      nloptr::nloptr(x0 = p0_sc, eval_f = eval_f_sc,
                     eval_grad_f = eval_grad_sc,
                     lb = lb_sc, ub = ub_sc,
                     opts = list(algorithm = .ctl$algorithm,
                                 ftol_rel  = .ctl$ftol_rel,
                                 maxeval   = .ctl$maxeval))
    })
    opt <- list(objective = opt_raw$objective,
                solution  = opt_raw$solution * sc,
                message   = opt_raw$message)
    if (.ctl$print > 0L) {
      row <- .admProgressRow(sprintf("%04d \u2713", .iter), opt$objective, opt$solution, pinfo)
      if (!is.null(row)) message(paste0(row, "\n",
        .admProgressTimingRow((proc.time() - t0)["elapsed"], pinfo)))
    }
    opt$all_traces <- list(list(restart_id = 1L,
                                nll_trace  = .nll_trace,
                                par_trace  = .par_trace))
  } else {
    .admSetupDaemons(.ctl, .ctl$n_restarts)
    on.exit(.admStopDaemons(), add = TRUE)
    opt <- .admRunRestarts(
      worker_fn  = .adfoRestartWorker,
      p0         = ov$p0, ov = ov, pinfo = pinfo,
      .ctl       = .ctl, ui = .ui, studies = studies,
      extra_args = list(algorithm       = .ctl$algorithm,
                        ftol_rel        = .ctl$ftol_rel,
                        maxeval         = .ctl$maxeval,
                        use_grad        = want_grad,
                        use_central     = use_central,
                        use_pure_fd     = use_pure_fd,
                        # The measured steps, not the constant -- restarts must
                        # difference the same way the sequential path does.
                        grad_h          = grad_h_v,
                        grad_bounds     = .ctl$grad_bounds,
                        output_var      = output_var,
                        print_progress  = TRUE,
                        print           = .ctl$print,
                        cores           = .ctl$cores,
                        rxMod_direct    = rxMod,
                        sensModel_direct = sensModel)
    )
    .admStopDaemons()
    .iter <- opt$n_iter
  }

  t_opt     <- (proc.time() - t0)["elapsed"]
  # A gradient fit is confined to <box centre> +/- grad_bounds; say so if it
  # stopped there rather than at an interior optimum. adfo only acquired this
  # constraint in 0.4.1, when grad = "analytical" became the default (want_grad
  # gates it). The centre is the WINNING RESTART's own init, not p0 -- see
  # .admScaledOptimize()'s box_centre; the sequential path has no restart and is
  # centred on p0, hence the fallback.
  if (want_grad)
    .admWarnOnBounds(opt$solution, opt$box_centre %||% ov$p0, ov,
                     .ctl$grad_bounds, pinfo)
  final     <- .admUnpack(opt$solution, pinfo)
  fullTheta <- .admFullTheta(final, pinfo)

  p_hat  <- setNames(opt$solution, names(ov$p0))
  t0_cov <- proc.time()
  .cov <- if (.ctl$covMethod == "r") {
    # struct + sigma + OMEGA: the Hessian spans all three, so the evaluation
    # count must too.
    np_cov     <- length(pinfo$struct_names) + length(pinfo$sigma_names) +
                  length(pinfo$omega_par)
    # Joint fits use the NLL-FD Hessian (the grad-FD Hessian calls the
    # single-output analytical .adfoGrad, which does not handle joint units).
    #
    # `have_d2`, not `want_grad`: forward-differencing the gradient is only a
    # sound way to build a Hessian if that gradient is smooth, i.e. actually
    # analytic. Gated on want_grad alone it fired in two cases where it is not:
    #   * the sens model failed to load -- want_sens FALSE but want_grad still
    #     TRUE -- so the Hessian forward-differenced an already-finite-differenced
    #     gradient. .adfoGrad turns any anyNA column into zeros, giving a singular
    #     H, a NULL cov and "standard errors are unavailable for this fit".
    #   * a transformed endpoint, where the order-2 block is deliberately absent
    #     and Pass 2 already forward-differences the whole NLL at grad_h = 1e-4 --
    #     so the covariance FD at cov_h_outer ~ 2.4e-3 nested on top of it.
    # Both then take the Gill NLL-FD Hessian, which is what 0.4.0 used by default
    # and what the covariance calibration in the docs was measured on.
    use_grad_cov <- want_grad && have_d2
    n_evals <- if (use_grad_cov) {
      np_cov + 1L
    } else {
      n_off <- np_cov * (np_cov - 1L) / 2L
      1L + 2L * np_cov + 4L * n_off
    }
    evals_label <- if (use_grad_cov) "gradient evaluations" else "NLL evaluations"
    # use_grad_cov implies have_d2 implies want_sens, so the gradient this
    # differences is always the analytic one.
    hess_label  <- if (use_grad_cov) ", Analytical-Hessian" else ""
    message(sprintf("  Computing covariance (R method%s, %d %s)", hess_label, n_evals, evals_label))
    tryCatch(
      .adfoCalcCov(p_hat, pinfo, studies, sensModel, rxMod, output_var,
                   params_list, cores, use_grad = use_grad_cov,
                   grad_h = .ctl$grad_h, cov_h = .ctl$cov_h,
                   cov_h_outer = .ctl$cov_h_outer, gill = .ctl$gill),
      error = function(e) { warning("adfoCalcCov failed: ", conditionMessage(e)); NULL })
  } else NULL
  # A NULL covariance used to be completely silent: no warning reached the user,
  # `warnings()` was empty, covMethod came back "" and every SE was NA with no
  # indication why. Say so once, from the driver, where it cannot be swallowed.
  if (isTRUE(.ctl$covMethod == "r") && is.null(.cov))
    warning("covariance could not be computed (the Hessian was singular or ",
            "non-finite); standard errors are unavailable for this fit.",
            call. = FALSE)
  # iniDf order first (nlmixr2est maps SEs positionally), then snapshot the names
  # BEFORE nlmixr2est sees it -- .admCovThetaOrder()/.admRestoreCovNames().
  .cov      <- .admCovThetaOrder(.cov, .ui)
  .cov_nms  <- .admCovNames(.cov)
  t_cov     <- (proc.time() - t0_cov)["elapsed"]
  t_elapsed <- t_opt + t_cov

  if (.ctl$returnAdmr) {
    return(list(objective = opt$objective, fullTheta = fullTheta,
                struct = final$struct, sigma_var = final$sigma_var,
                omega = final$omega, L = final$L, nloptr = opt,
                cov = .cov))
  }

  .ret            <- new.env(parent = emptyenv())
  .ret$table      <- env$table
  .ret$ui         <- .ui
  .ret$fullTheta  <- fullTheta
  .ret$objective  <- opt$objective
  .ret$est        <- "adfo"
  .ret$ofvType    <- "adfo"
  .ret$adjObf     <- FALSE
  .ret$covMethod  <- if (!is.null(.cov)) "r" else ""
  .ret$cov        <- .cov
  .ret$message    <- opt$message
  .ret$extra      <- ""
  .ret$origData   <- studies

  .ret$admExtra <- list(struct         = final$struct,
                        sigma_var      = final$sigma_var,
                        sigma_is_prop  = pinfo$sigma_is_prop,
                        sigma_is_lnorm = pinfo$sigma_is_lnorm,
                        # the TBS residual quadrature the FIT used -- plot.admFit()
                        # re-parses the ui, and .admParseIniDf() carries no
                        # resid_nodes, so without this the diagnostics recompute
                        # V_pred at the 81-node default and disagree with the fit
                        resid_nodes    = pinfo$resid_nodes,
                        # ... and the solver tolerance the FIT used, for the
                        # same reason: plot.admFit() re-solves the model to build
                        # the predicted mean and covariance panels, and a fit run
                        # at a looser sigdig diagnosed against rxode2's stock
                        # tolerances shows standardised-residual structure the
                        # objective was never minimised on.
                        sigdig         = pinfo$sigdig,
                        omega          = final$omega,
                        L              = final$L,
                        eta_col_names  = pinfo$eta_col_names,
                        par_names      = names(ov$p0),
                        npar           = length(ov$p0),
                        nloptr         = opt,
                        nll_trace      = .nll_trace,
                        par_trace      = .par_trace,
                        all_traces     = opt$all_traces,
                        n_iter         = .iter,
                        time           = t_elapsed,
                        t_opt          = t_opt,
                        t_cov          = t_cov,
                        studies        = studies,
                        n_sim          = 5000L,  # for diagnostic plots; FO itself is n_sim=1
                        sampling       = "sobol")

  .admFinaliseFit(.ret, .ui, .ctl, est = "adfo", objective = opt$objective,
                  ov = ov, studies = studies, cov = .cov,
                  cov_nms = .cov_nms, multi_out = multi_out,
                  extra_field = "admExtra",
                  handle_ctl = nmObjHandleControlObject.adfoControl,
                  t_opt = t_opt, t_cov = t_cov, t_elapsed = t_elapsed)
}
