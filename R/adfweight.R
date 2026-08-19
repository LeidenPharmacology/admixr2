# =============================================================================
# ADF weight matrix -- the sampling law of the reported summary
# =============================================================================
#
# The aggregate objective scores (ybar, V) as if it were the log-likelihood of N
# iid draws from N(yt, Vt). That is exact only when each subject's OBSERVATION
# VECTOR is multivariate normal, and it is not: y_i = f(theta, a_i, b_i) + eps_i
# with f nonlinear in b_i, so the marginal is a mixture over (a_i, b_i) and is
# normal only where f is linear in both. The covariate is not what breaks it --
# nonlinearity in the random effect alone is enough.
#
# What that costs is NOT the point estimates. The score is
# -2 (dtau/dPsi)' W^-1 (t - tau) and E[t] = tau at the true Psi for ANY W
# (Gourieroux-Monfort-Trognon), so every fit is consistent whatever the weight.
# It costs the reported UNCERTAINTY, in two ways that differ in kind:
#
#   Cov(V_ij, V_kl)   assumed (V_ik V_jl + V_il V_jk)/N   true (mu4 - V V)/N
#   Cov(ybar, vech V) assumed 0                           true mu3/N
#
# The first is mis-sized by the excess kurtosis. The second is a zero where a
# real correlation of 0.3-0.6 sits: for a multivariate normal the sample mean and
# sample covariance are exactly independent, and for anything else they are not.
# A sample that comes out high also comes out more spread, and the current
# objective counts the two channels as independent evidence.
#
# So score t = (ybar, vech V) against its own asymptotic law instead. That is
# Browne's ADF estimator, with Omega computed FROM THE MODEL rather than
# estimated from the sample -- which is what removes ADF's small-sample failure,
# since the sample estimate of a fourth moment is what needs enormous N.
# tau(Psi) is unchanged, so the covariate machinery, the shift/absorption paths
# and the quadrature all carry over: this replaces the scoring, not the model.

# Conditional central moments of the residual at every node, per timepoint.
#
# The node ensemble carries f; given the node the residual is a draw from the
# endpoint's own distribution with that node's f as its mean parameter. What the
# weight needs from it is the 2nd, 3rd and 4th CENTRAL moments -- the expansion
# below is exact for any residual that is independent across timepoints given
# the node, not only a normal one. Returns three Q x m matrices matching `cp`,
# or NULL when the family cannot supply them.
#
# Under normality t3 = 0 and q4 = 3 d^2, and the two extra terms in the
# expansion vanish identically, so add/prop/pow/combined stay bit-for-bit what
# they were. lnorm does NOT: its conditional law is lognormal, and treating it
# as normal understated the fourth moment by exactly the lognormal excess
# kurtosis -- which is the very quantity the magnitude estimate for this whole
# correction is written in terms of.
#
# Dispatch is PER ROW. A single unit has one form today (multi-output studies
# are separate units and joint units are refused), but keying the whole matrix
# off `form[[1]]` is a trap that costs nothing to avoid.
.admAdfCondMom <- function(cp, arr) {
  Q <- nrow(cp); m <- ncol(cp)
  # ar() correlates the residual ACROSS timepoints, so the products below no
  # longer factor and every cross term the expansion drops is real.
  if (any(!is.na(arr$rho))) return(NULL)
  d <- t3 <- q4 <- matrix(NA_real_, Q, m)
  col <- function(x, j) if (length(x) == 1L) x else x[[j]]
  for (j in seq_len(m)) {
    f  <- cp[, j]
    a2 <- col(arr$a2, j); b2 <- col(arr$b2, j); cc <- col(arr$cc, j)
    vm <- col(arr$vmul, j); sz <- col(arr$csz, j); ph <- col(arr$phi, j)
    fm <- abs(f)
    v <- switch(as.character(arr$form[[j]]),
      "0" = a2 + b2 * fm^(2 * cc),                     # combined2
      "1" = (sqrt(a2) + sqrt(b2) * fm^cc)^2,           # combined1
      NULL)
    if (!is.null(v)) {
      d[, j] <- v
      if (isTRUE(all.equal(vm, 1))) {                  # normal
        t3[, j] <- 0; q4[, j] <- 3 * v^2
      } else {
        # Student-t, folded in as vmul = nu/(nu-2). Symmetric, so t3 = 0; the
        # fourth moment needs nu > 4 and there is nothing sensible to return
        # below that -- a t_3 residual has no finite kurtosis, so the sampling
        # law of the reported V does not have the variance the weight is made of.
        nu <- 2 * vm / (vm - 1)
        if (!is.finite(nu) || nu <= 4) return(NULL)
        t3[, j] <- 0; q4[, j] <- 3 * v^2 * (nu - 2) / (nu - 4)
      }
      next
    }
    switch(as.character(arr$form[[j]]),
      "2" = {                                          # lnorm, exact
        # y | node ~ LogNormal(log f, sv): M1 = f exp(sv/2), w = exp(sv)
        sv <- a2; w <- exp(sv); M1 <- f * exp(sv / 2)
        d[, j]  <- M1^2 * (w - 1)
        t3[, j] <- M1^3 * (w - 1)^2 * (w + 2)
        q4[, j] <- M1^4 * (w - 1)^2 * (w^4 + 2 * w^3 + 3 * w^2 - 3)
      },
      "4" = {                                          # pois(f)
        lam <- f
        d[, j] <- lam; t3[, j] <- lam; q4[, j] <- lam + 3 * lam^2
      },
      "5" = {                                          # binom(N, f), N constant
        if (!is.finite(sz)) return(NULL)
        pp <- f; v <- sz * pp * (1 - pp)
        d[, j]  <- v
        t3[, j] <- v * (1 - 2 * pp)
        q4[, j] <- v * (1 + (3 * sz - 6) * pp * (1 - pp))
      },
      "6" = {                                          # nbinomMu(size, f)
        if (!is.finite(sz)) return(NULL)
        mu <- f; k <- sz; v <- mu + mu^2 / k
        d[, j]  <- v
        t3[, j] <- mu * (1 + 3 * mu / k + 2 * mu^2 / k^2)
        q4[, j] <- mu * (1 + 7 * mu / k + 12 * mu^2 / k^2 + 6 * mu^3 / k^3) +
                   3 * v^2
      },
      "7" = {                                          # beta(mu, phi)
        if (!is.finite(ph)) return(NULL)
        al <- f * ph; be <- (1 - f) * ph; sm <- al + be
        v  <- al * be / (sm^2 * (sm + 1))
        sk <- 2 * (be - al) * sqrt(sm + 1) / ((sm + 2) * sqrt(al * be))
        ku <- 3 + 6 * ((al - be)^2 * (sm + 1) - al * be * (sm + 2)) /
                      (al * be * (sm + 2) * (sm + 3))
        d[, j] <- v; t3[, j] <- sk * v^1.5; q4[, j] <- ku * v^2
      },
      # form 3 (TBS) needs the residual quadrature carried to third and fourth
      # order and form 8 (ordinal) is a joint unit, refused upstream.
      return(NULL))
  }
  if (!all(is.finite(d)) || !all(is.finite(t3)) || !all(is.finite(q4)))
    return(NULL)
  list(d = d, t3 = t3, q4 = q4)
}

# Omega / N: the asymptotic covariance of (ybar, vech V).
#
# REFERENCE IMPLEMENTATION. .admAdfWeightFast computes the same matrix with the
# node contraction hoisted out of the q x q loop and is what runs; this one is
# the readable statement of the expansion and the oracle the fast path is pinned
# against. Keep them in step.
#
# `C` is the CENTRED conditional means (Q x m), `w` the node weights summing to
# one, `Dv` the conditional residual variances (Q x m), `N` the subjects. The
# mu3 and mu4 blocks are Isserlis/Wick expansions at the node: given the node the
# residual is normal with diagonal covariance Dv, so every odd pairing collapses
# and the even ones are sums of products of C and Dv.
#
# One weight per source, frozen at a first stage -- see .admAdfFreeze.
.admAdfWeight <- function(C, w, Dv, N, T3 = NULL, Q4 = NULL) {
  w  <- w / sum(w)
  m  <- ncol(C)
  ij <- which(lower.tri(diag(m), diag = TRUE), arr.ind = TRUE)
  q  <- nrow(ij)
  if (is.null(T3)) T3 <- matrix(0, nrow(C), m)
  if (is.null(Q4)) Q4 <- 3 * Dv^2
  dbar <- colSums(w * Dv)
  S    <- crossprod(C, w * C); diag(S) <- diag(S) + dbar
  dl   <- function(i, j) if (i == j) Dv[, i] else 0
  # E[e_i e_j e_k] under conditional independence: non-zero only when all three
  # indices coincide. Zero for a normal residual, which is why the terms guarded
  # by tl() below have no effect on add/prop/combined.
  tl   <- function(i, j, k) if (i == j && j == k) T3[, i] else 0
  # E[e^4] - 3 Var^2: the excess-kurtosis correction to the one quartic case the
  # three delta-delta products get wrong. Also identically zero under normality.
  ke   <- function(i, j, k, l)
    if (i == j && j == k && k == l) Q4[, i] - 3 * Dv[, i]^2 else 0
  W    <- matrix(0, m + q, m + q)
  W[seq_len(m), seq_len(m)] <- S / N
  for (b in seq_len(q)) {
    k <- ij[b, 1L]; l <- ij[b, 2L]
    for (i in seq_len(m)) {
      t <- C[, i] * C[, k] * C[, l] +
           C[, i] * dl(k, l) + C[, k] * dl(i, l) + C[, l] * dl(i, k) +
           tl(i, k, l)
      W[i, m + b] <- W[m + b, i] <- sum(w * t) / N
    }
  }
  for (a in seq_len(q)) {
    i <- ij[a, 1L]; j <- ij[a, 2L]
    for (b in seq_len(q)) {
      k <- ij[b, 1L]; l <- ij[b, 2L]
      t <- C[, i] * C[, j] * C[, k] * C[, l] +
           C[, i] * C[, j] * dl(k, l) + C[, i] * C[, k] * dl(j, l) +
           C[, i] * C[, l] * dl(j, k) + C[, j] * C[, k] * dl(i, l) +
           C[, j] * C[, l] * dl(i, k) + C[, k] * C[, l] * dl(i, j) +
           dl(i, j) * dl(k, l) + dl(i, k) * dl(j, l) + dl(i, l) * dl(j, k) +
           C[, i] * tl(j, k, l) + C[, j] * tl(i, k, l) +
           C[, k] * tl(i, j, l) + C[, l] * tl(i, j, k) +
           ke(i, j, k, l)
      W[m + a, m + b] <- (sum(w * t) - S[i, j] * S[k, l]) / N
    }
  }
  W
}

# The same weight, with the node contraction done once instead of q^2 times.
#
# Every term in the Wick expansion is a weighted sum over nodes of a product of
# at most four C columns and Dv columns, so each DISTINCT contraction can be
# formed as one crossprod and the q x q assembly reduces to indexing:
#
#   P[, a] = C_i C_j            T1 = P' w P          the four-C term
#   PD     = P' w Dv            A_j = C' (w Dv_j) C  one C-pair with one Dv
#   B      = Dv' w Dv           the two-Dv terms
#
# Cost goes from O(q^2 Q) to O(m^3 Q + q^2), which is what makes the ceiling the
# handoff quotes (m ~ 30, q = 465) reachable rather than theoretical.
.admAdfWeightFast <- function(C, w, Dv, N, T3 = NULL, Q4 = NULL) {
  w  <- w / sum(w)
  m  <- ncol(C)
  ij <- which(lower.tri(diag(m), diag = TRUE), arr.ind = TRUE)
  q  <- nrow(ij)
  I  <- ij[, 1L]; J <- ij[, 2L]
  dbar <- colSums(w * Dv)
  S    <- crossprod(C, w * C); diag(S) <- diag(S) + dbar
  # Third/fourth conditional moments. NULL means "normal", where both extra
  # contractions are zero and every line below that touches them is a no-op --
  # which is what keeps add/prop/combined bit-identical to the version that had
  # no concept of them.
  TW <- if (is.null(T3)) NULL else colSums(w * T3)              # m
  CT <- if (is.null(T3)) NULL else crossprod(C, w * T3)         # m x m
  KE <- if (is.null(Q4)) NULL else colSums(w * (Q4 - 3 * Dv^2)) # m

  P  <- C[, I, drop = FALSE] * C[, J, drop = FALSE]     # Q x q
  wP <- w * P
  T1 <- crossprod(P, wP)                                # q x q
  PD <- crossprod(P, w * Dv)                            # q x m
  B  <- crossprod(Dv, w * Dv)                           # m x m
  A  <- lapply(seq_len(m), function(jj) crossprod(C, (w * Dv[, jj]) * C))
  CD <- crossprod(C, w * Dv)                            # m x m, mu3 helper

  W <- matrix(0, m + q, m + q)
  W[seq_len(m), seq_len(m)] <- S / N

  # mu3: E[C_i C_k C_l] + C_i D_kl + C_k D_il + C_l D_ik
  CC <- vapply(seq_len(q), function(b)
    as.numeric(crossprod(C, wP[, b])), numeric(m))      # m x q
  M3 <- CC
  for (b in seq_len(q)) {
    k <- I[b]; l <- J[b]
    # C_i D_kl  (needs k == l), then C_k D_il and C_l D_ik, whose delta fixes
    # the ROW: the i they select is l and k respectively.
    if (k == l) M3[, b] <- M3[, b] + CD[, k]
    M3[l, b] <- M3[l, b] + CD[k, l]
    M3[k, b] <- M3[k, b] + CD[l, k]
    # E[e_i e_k e_l] survives only at i = k = l
    if (!is.null(TW) && k == l) M3[k, b] <- M3[k, b] + TW[k]
  }
  W[seq_len(m), m + seq_len(q)] <- M3 / N
  W[m + seq_len(q), seq_len(m)] <- t(M3) / N

  M4 <- T1
  eqA <- I == J
  for (b in seq_len(q)) {
    k <- I[b]; l <- J[b]
    if (k == l) M4[, b] <- M4[, b] + PD[, k]
    M4[eqA, b] <- M4[eqA, b] + PD[b, I[eqA]]
    for (a in seq_len(q)) {
      i <- I[a]; j <- J[a]
      v <- 0
      if (j == l) v <- v + A[[j]][i, k]
      if (j == k) v <- v + A[[j]][i, l]
      if (i == l) v <- v + A[[i]][j, k]
      if (i == k) v <- v + A[[i]][j, l]
      if (i == j && k == l) v <- v + B[i, k]
      if (i == k && j == l) v <- v + B[i, j]
      if (i == l && j == k) v <- v + B[i, j]
      # the four C x E[eee] terms, each alive only where its own three indices
      # coincide, and the one quartic case the delta-delta products get wrong
      if (!is.null(CT)) {
        if (j == k && k == l) v <- v + CT[i, j]
        if (i == k && k == l) v <- v + CT[j, i]
        if (i == j && j == l) v <- v + CT[k, i]
        if (i == j && j == k) v <- v + CT[l, i]
      }
      if (!is.null(KE) && i == j && j == k && k == l) v <- v + KE[i]
      M4[a, b] <- M4[a, b] + v
    }
  }
  W[m + seq_len(q), m + seq_len(q)] <- (M4 - tcrossprod(S[cbind(I, J)])) / N
  W
}

# The normal-theory weight from the same S -- what the current objective implies.
# Kept because the difference between the two IS the correction, and a control
# arm that reproduces the current fit is the cheapest proof the rest is wired up.
.admAdfWeightNormal <- function(S, N) {
  m  <- ncol(S)
  ij <- which(lower.tri(diag(m), diag = TRUE), arr.ind = TRUE)
  q  <- nrow(ij)
  W  <- matrix(0, m + q, m + q)
  W[seq_len(m), seq_len(m)] <- S / N
  for (a in seq_len(q)) {
    i <- ij[a, 1L]; j <- ij[a, 2L]
    for (b in seq_len(q)) {
      k <- ij[b, 1L]; l <- ij[b, 2L]
      W[m + a, m + b] <- (S[i, k] * S[j, l] + S[i, l] * S[j, k]) / N
    }
  }
  W
}

# The node-level quantities the weight needs, alongside the predicted moments.
# Same path .adghMoments takes, so the two cannot describe different node sets.
.admAdfParts <- function(pars, pinfo, study, rxMod, out_var, grid, cores) {
  g     <- .adghGrid(pars, pinfo, grid, study)
  study <- .adghStudyCov(study, g)
  pm    <- .admMakeParamsList(nrow(g$eta), pinfo, 1L)[[1L]]
  cp    <- .admSimulate(rxMod, pars$struct, pinfo$sigma_names, g$eta, study,
                        out_var, pm, cores, pinfo$nDisplayProgress, pinfo$sigdig)
  sm  <- .adghStructMoments(cp, g$W, g$taylor)
  arr <- .admUnitResidRows(pinfo, out_var, pars$sigma_var, length(sm$mu),
                           phi = attr(cp, "phi"))
  m   <- .admResidMoments(sm$mu, diag(sm$V), arr, sm$V, study$times)
  # C carries the residual's MEAN SCALING. On an lnorm endpoint the conditional
  # mean is f exp(s/2), not f, so an unscaled C makes the weight's own S differ
  # from the V the objective scores against -- the two would then describe
  # different laws. ms is 1 on every other family, so this is a no-op there.
  cm  <- .admAdfCondMom(cp, arr)
  list(E = m$mu, V = m$V, w = g$W / sum(g$W),
       C  = if (is.null(m$ms)) sm$cpc else sweep(sm$cpc, 2L, m$ms, "*"),
       Dv = cm$d, T3 = cm$t3, Q4 = cm$q4)
}

# Freeze one weight per study at a first-stage estimate.
#
# FROZEN is not an approximation, it is a condition of the estimator. A
# Psi-dependent W inside log|W| contributes score terms whose expectation is not
# zero, so the two-stage form is the one that stays consistent -- and the current
# objective is a fine first stage, being consistent for any weight.
.admAdfFreeze <- function(p, pinfo, studies, rxMod, out_var, grid, cores) {
  pars <- .admUnpack(p, pinfo)
  lapply(studies, function(s) {
    pt <- .admAdfParts(pars, pinfo, s, rxMod, s$output %||% out_var, grid, cores)
    if (is.null(pt$Dv)) return(NULL)          # residual outside the free family
    W  <- .admAdfWeightFast(pt$C, pt$w, pt$Dv, as.numeric(s$n), pt$T3, pt$Q4)
    ch <- tryCatch(chol(W), error = function(e) NULL)
    if (is.null(ch)) return(NULL)
    list(Wi = chol2inv(ch), ldet = 2 * sum(log(diag(ch))))
  })
}

# -2 log L for the summary vector scored against its own sampling law.
#
# tau's covariance block is (N-1)/N * Vt, not Vt: under the ML denominator
# E[V] = (N-1)/N Vt, and scoring an unaligned tau makes the estimator MORE
# biased than the one it replaces rather than less -- the O(1/N) term is required,
# not cosmetic. log|W| is constant once W is frozen and is carried only so the
# objective stays on a comparable scale.
.admAdfNLL <- function(p, pinfo, studies, rxMod, out_var, grid, cores, Wl) {
  pars <- tryCatch(.admUnpack(p, pinfo), error = function(e) NULL)
  if (is.null(pars) || !.admParsFinite(pars, pinfo)) return(Inf)
  tot <- 0
  for (i in seq_along(studies)) {
    s <- studies[[i]]; wi <- Wl[[i]]
    if (is.null(wi)) return(Inf)
    pt <- tryCatch(.admAdfParts(pars, pinfo, s, rxMod, s$output %||% out_var,
                                grid, cores), error = function(e) NULL)
    if (is.null(pt) || !all(is.finite(pt$E)) || !all(is.finite(pt$V))) return(Inf)
    N  <- as.numeric(s$n)
    lo <- lower.tri(pt$V, diag = TRUE)
    d  <- c(as.numeric(s$E) - pt$E,
            as.numeric(s$V[lo]) - as.numeric(((N - 1) / N * pt$V)[lo]))
    tot <- tot + as.numeric(crossprod(d, wi$Wi %*% d)) + wi$ldet
  }
  if (is.finite(tot)) tot else Inf
}

# =============================================================================
# The sandwich: covMethod = "r,s"
# =============================================================================
#
# Avar = H^-1 J H^-1,   H = d2F/dPsi dPsi' at the optimum,
#                       J = sum_s G_s Omega_s G_s',  G_s = d2F_s/(dPsi dt_s')
#
# H is the Hessian of the objective ACTUALLY minimised -- the same one
# covMethod = "r" reports -- and that is what guarantees the reduction: under
# correct specification J = 2H, so Avar collapses to 2 H^-1 by construction
# rather than by hope.
#
# AN EARLIER VERSION OF THIS USED THE GLS SURROGATE and was wrong. Eq. (1) is
# GLS on t with the normal-theory weight only ASYMPTOTICALLY: F is LINEAR in V
# and quadratic in ybar, while a GLS criterion is quadratic in both, so the two
# share a score and an expected information at t = tau and nowhere else. Using
# (G' Wn^-1 G)^-1 as the bread and dtau/dPsi as G therefore drifts by terms in
# (t - tau) -- measured on a badly-fitting fixture as 1.0006 on a well-determined
# structural theta rising to 1.46 on log(om^2), which is the
# (t - tau) . d2tau/dPsi2 signature exactly.
#
# The cross-derivative, from F = N( log|Vt| + tr(Vt^-1 V) + r' Vt^-1 r ):
#
#   dF/dybar = 2N Vt^-1 r
#     => d2F/(dPsi dybar') = 2N [ (dVt^-1/dPsi) r  -  Vt^-1 dyt/dPsi ]
#
#   dF/dV_ij = N (Vt^-1)_ij            (x2 for an off-diagonal vech entry)
#     => d2F/(dPsi dV_ij) = N d(Vt^-1)_ij/dPsi      (x2 off-diagonal)
#
# Note what the V block does NOT contain: any dtau/dPsi. That half of F is
# linear in V, so only d(Vt^-1)/dPsi = -Vt^-1 (dVt/dPsi) Vt^-1 survives.
.admScoreCross <- function(E, V, dE, dV, s, N) {
  m   <- length(E)
  isv <- identical(s$method, "var")
  # A `var` study is scored by nll_var_cpp, whose objective is the DIAGONAL one
  #   sum_i N( log Vt_ii + V_ii/Vt_ii + r_i^2/Vt_ii )
  # so dF/dV_ii = N / Vt_ii, NOT N (Vt^-1)_ii, and dF/dybar_i = 2N r_i / Vt_ii.
  # Those coincide only when Vt is diagonal. Using the full inverse for both
  # branches broke the information equality on exactly this branch: eigenvalues
  # of J(Wn)/2H came out 0.06 / 0.41 / 3.29 / 1782 instead of all ones, while
  # the cov branch was exact.
  Vi <- if (isv) diag(1 / diag(V), m) else
    tryCatch(chol2inv(chol(V)), error = function(e) NULL)
  if (is.null(Vi) || !all(is.finite(Vi))) return(NULL)
  r  <- as.numeric(s$E) - as.numeric(E)
  p  <- length(dV)
  ij <- which(lower.tri(diag(m), diag = TRUE), arr.ind = TRUE)
  if (isv) ij <- ij[ij[, 1L] == ij[, 2L], , drop = FALSE]
  G <- matrix(0, p, m + nrow(ij))
  for (k in seq_len(p)) {
    # on the var branch only the diagonal of dV reaches the objective
    dVk <- if (isv) diag(diag(dV[[k]]), m) else dV[[k]]
    dVi <- -Vi %*% dVk %*% Vi
    G[k, seq_len(m)] <- 2 * N * (as.numeric(dVi %*% r) -
                                 as.numeric(Vi %*% dE[, k]))
    dup <- ifelse(ij[, 1L] == ij[, 2L], 1, 2)      # vech duplication
    G[k, m + seq_len(nrow(ij))] <- N * dup * dVi[ij]
  }
  G
}

.admSandwich <- function(H, G, Om) {
  Hi <- tryCatch(solve(H), error = function(e) NULL)
  if (is.null(Hi)) return(NULL)
  p <- nrow(H); J <- matrix(0, p, p)
  for (i in seq_along(G)) J <- J + G[[i]] %*% Om[[i]] %*% t(G[[i]])
  list(cov = Hi %*% J %*% Hi, bread = 2 * Hi, J = J, H = H)
}

# The summary a study actually reports, stacked: (ybar, vech V) for a full
# covariance and (ybar, diag V) for a variance-only study.
#
# A `method = "var"` study is not a degenerate covariance study -- it reports
# fewer numbers, and its weight is the corresponding MARGINAL of the full one
# rather than a different derivation. Scoring covariances the fit never saw
# would invent information.
.admTauVec <- function(E, V, s) {
  if (identical(s$method, "var")) c(as.numeric(E), diag(V))
  else c(as.numeric(E), V[lower.tri(V, diag = TRUE)])
}

# The weight the OBJECTIVE implicitly uses -- the baseline the sandwich corrects
# away from.
#
# For a `cov` study that is the normal-theory covariance of (ybar, vech V). For a
# `var` study it is NOT the marginal of that: nll_var_cpp scores
#   sum_i N( log v_i + V_ii/v_i + r_i^2/v_i )
# which treats the m variances, and the m mean residuals, as INDEPENDENT. The
# true normal-theory marginal still has Cov(V_ii, V_jj) = 2 V_ij^2 / N. So the
# var branch's baseline is WORKING INDEPENDENCE, and that -- not kurtosis -- is
# the bulk of what is wrong with it.
#
# Using the marginal here instead left eigen(J/2H) at 0.38 .. 3.26 rather than
# all ones, i.e. it described a weight the objective does not use.
.admWorkingWeight <- function(V, N, method) {
  if (!identical(method, "var")) return(.admAdfWeightNormal(V, N))
  v <- diag(V); m <- length(v)
  W <- matrix(0, 2 * m, 2 * m)
  W[seq_len(m), seq_len(m)] <- diag(v / N, m)
  W[m + seq_len(m), m + seq_len(m)] <- diag(2 * v^2 / N, m)
  W
}

# Restrict a full (m + q) weight to the rows a `var` study reports.
.admWeightSel <- function(W, m, method) {
  if (!identical(method, "var")) return(W)
  ij  <- which(lower.tri(diag(m), diag = TRUE), arr.ind = TRUE)
  keep <- c(seq_len(m), m + which(ij[, 1L] == ij[, 2L]))
  W[keep, keep, drop = FALSE]
}

# d(yt)/dPsi and d(Vt)/dPsi per study, by central difference on the MOMENTS.
#
# These are analytic from what the gradient machinery already forms; this is the
# reference until that extraction is written. Differencing the MOMENTS rather
# than the objective keeps it well conditioned, and it runs once, post-fit.
.admMomentDeriv <- function(p_hat, pinfo, studies, rxMod, out_var, grid, cores,
                            h = 1e-5, mom_fn = NULL) {
  # `mom_fn` is the ESTIMATOR's moment map, and it is separate from the ensemble
  # the weight is built on for a reason: G describes the objective that was
  # minimised, Omega describes the law the data actually came from. For adgh the
  # two coincide and the default is used. For adfo they do NOT -- adfo predicts
  # V = J Omega J' + Sigma, whose implied individual law is exactly normal, so
  # scoring it against its own assumption would make the sandwich identically
  # 2H^-1 and say nothing. Passing adfo's moment map here and keeping the
  # quadrature ensemble for Omega is what makes the correction meaningful there:
  # it scores an FO fit against the model's true nonlinear law.
  mom <- mom_fn %||% function(pp) {
    pars <- .admUnpack(pp, pinfo)
    lapply(studies, function(s)
      .admAdfParts(pars, pinfo, s, rxMod, s$output %||% out_var, grid, cores))
  }
  p  <- length(p_hat)
  b  <- mom(p_hat)
  dE <- lapply(b, function(x) matrix(0, length(x$E), p))
  dV <- lapply(b, function(x) rep(list(matrix(0, nrow(x$V), ncol(x$V))), p))
  for (k in seq_len(p)) {
    a  <- p_hat; a[k]  <- a[k] + h
    cc <- p_hat; cc[k] <- cc[k] - h
    ma <- mom(a); mc <- mom(cc)
    for (i in seq_along(b)) {
      dE[[i]][, k] <- (ma[[i]]$E - mc[[i]]$E) / (2 * h)
      dV[[i]][[k]] <- (ma[[i]]$V - mc[[i]]$V) / (2 * h)
    }
  }
  list(E = lapply(b, `[[`, "E"), V = lapply(b, `[[`, "V"), dE = dE, dV = dV)
}

# Post-fit sandwich covariance for a fitted parameter vector.
#
# `H` is the Hessian of the objective at the optimum -- the SAME one
# covMethod = "r" reports -- passed in rather than rebuilt here, so "r,s" and
# "r" cannot disagree about the half they share.
#
# Returns NULL rather than guessing whenever an ingredient is unavailable: a
# residual outside the conditionally-normal family, a singular weight, a failed
# moment solve. The caller falls back to "r" and says so.
.admSandwichCov <- function(p_hat, pinfo, studies, rxMod, out_var, grid, cores,
                            H, md = NULL, keep = NULL, mom_fn = NULL) {
  # H is checked FIRST, and here rather than being left to solve() inside .admSandwich,
  # whose tryCatch is there for a singular matrix and cannot tell that apart from
  # an H that was never supplied. A caller that forgot the argument then gets a
  # silent NULL -- which reads as "the sandwich does not apply to this fit"
  # rather than "you called it wrong", and a calibration study built on it
  # reported numbers for replicates it had skipped entirely.
  if (missing(H) || !is.matrix(H) || nrow(H) != ncol(H) || !all(is.finite(H)))
    stop(".admSandwichCov: `H` must be a finite square Hessian of the objective ",
         "at the optimum -- the same one covMethod = \"r\" inverts.", call. = FALSE)
  pars <- tryCatch(.admUnpack(p_hat, pinfo), error = function(e) NULL)
  if (is.null(pars)) return(NULL)
  # A joint (same-subject, multi-output) unit stacks several outputs into one
  # covariance. .admAdfParts solves ONE output, so its node ensemble is not the
  # one behind that block and the Wick expansion would be taken over the wrong
  # conditional law. Refuse rather than return a plausible wrong weight.
  if (any(vapply(studies, function(s) isTRUE(s$is_joint), logical(1))))
    return(NULL)
  md <- md %||% tryCatch(
    .admMomentDeriv(p_hat, pinfo, studies, rxMod, out_var, grid, cores,
                    mom_fn = mom_fn),
    error = function(e) NULL)
  if (is.null(md)) return(NULL)
  G <- Om <- vector("list", length(studies))
  for (i in seq_along(studies)) {
    s  <- studies[[i]]
    pt <- tryCatch(.admAdfParts(pars, pinfo, s, rxMod, s$output %||% out_var,
                                grid, cores), error = function(e) NULL)
    if (is.null(pt) || is.null(pt$Dv)) return(NULL)
    N <- as.numeric(s$n); m <- length(pt$E)
    Om[[i]] <- .admWeightSel(
      .admAdfWeightFast(pt$C, pt$w, pt$Dv, N, pt$T3, pt$Q4), m, s$method)
    G[[i]]  <- .admScoreCross(md$E[[i]], md$V[[i]], md$dE[[i]], md$dV[[i]], s, N)
    if (is.null(G[[i]])) return(NULL)
    # H may have been reduced to the struct+sigma sub-block (.admReduceNpdOmega),
    # so G must lose the same rows -- G is parameter-indexed by row.
    if (!is.null(keep)) G[[i]] <- G[[i]][keep, , drop = FALSE]
  }
  if (!is.null(keep) && nrow(H) != length(keep)) return(NULL)
  .admSandwich(H, G, Om)
}

# -- adfo ----------------------------------------------------------------------

# A quadrature grid for the sandwich WEIGHT, sized to the number of etas.
#
# adfo carries no node ensemble -- that is the point of FO -- but the weight
# needs one, because Omega is a property of the model's true nonlinear law and
# not of the linearisation used to fit it. The grid is built once, post-fit, so
# a node count that would be extravagant inside an optimisation loop is cheap
# here; it is still capped, since the product grid is NQ^n_eta and a 5-eta model
# at 9 nodes would be 59049 subjects in one solve for no accuracy that matters.
.admSandwichGrid <- function(pinfo, max_nodes = 5000L) {
  n_eta <- pinfo$n_eta
  if (is.null(n_eta) || n_eta < 1L) return(NULL)
  nq <- 9L
  while (nq > 3L && nq^n_eta > max_nodes) nq <- nq - 2L
  .adghNodeGrid(nq, n_eta)
}

# adfo's own (E, V) per study -- the moment map the FO objective minimises.
#
# Returned in the shape .admMomentDeriv expects, so the only thing that changes
# between estimators is this function. `V` is the full predicted covariance even
# for a `method = "var"` study: .admScoreCross takes the diagonal itself, which
# keeps the branch logic in one place rather than two.
.admAdfoMomFn <- function(pinfo, studies, sensModel, rxMod, out_var,
                          params_list, cores) {
  function(pp) {
    pars <- .admUnpack(pp, pinfo)
    lapply(seq_along(studies), function(i) {
      s   <- studies[[i]]
      ov  <- s$output %||% out_var
      n_t <- length(s$times)
      arr <- .admResidRows(pinfo, ov, pars$sigma_var, n_t)
      mj  <- .adfoGetMuJ(pars, pinfo, s, sensModel, rxMod, ov,
                         params_list[[i]], cores)
      vp  <- .adfoVpred(mj$mu, mj$J, pars$L, arr, n_t, pinfo$n_eta, s$times)
      list(E = vp$mu_sigma, V = vp$V)
    })
  }
}
