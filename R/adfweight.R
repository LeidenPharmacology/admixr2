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

# Conditional residual variance Var(y | node) at every node, per timepoint.
#
# The node ensemble carries f; the residual is conditionally normal given the
# node for the add/prop/combined family, so its variance is the same row-indexed
# formula .admResidRows encodes, evaluated at each node's own f rather than at
# the marginal mean. Returns Q x m to match `cp`.
.admAdfCondVar <- function(cp, arr) {
  a2 <- matrix(arr$a2, nrow(cp), ncol(cp), byrow = TRUE)
  b2 <- matrix(arr$b2, nrow(cp), ncol(cp), byrow = TRUE)
  cc <- matrix(arr$cc, nrow(cp), ncol(cp), byrow = TRUE)
  fm <- abs(cp)
  out <- switch(as.character(arr$form[[1L]]),
    "0" = a2 + b2 * fm^(2 * cc),                       # combined2
    "1" = (sqrt(a2) + sqrt(b2) * fm^cc)^2,             # combined1
    "2" = {                                            # lnorm
      sv <- a2
      (cp * exp(sv / 2))^2 * (exp(sv) - 1)
    },
    NULL)
  out
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
.admAdfWeight <- function(C, w, Dv, N) {
  w  <- w / sum(w)
  m  <- ncol(C)
  ij <- which(lower.tri(diag(m), diag = TRUE), arr.ind = TRUE)
  q  <- nrow(ij)
  dbar <- colSums(w * Dv)
  S    <- crossprod(C, w * C); diag(S) <- diag(S) + dbar
  dl   <- function(i, j) if (i == j) Dv[, i] else 0
  W    <- matrix(0, m + q, m + q)
  W[seq_len(m), seq_len(m)] <- S / N
  for (b in seq_len(q)) {
    k <- ij[b, 1L]; l <- ij[b, 2L]
    for (i in seq_len(m)) {
      t <- C[, i] * C[, k] * C[, l] +
           C[, i] * dl(k, l) + C[, k] * dl(i, l) + C[, l] * dl(i, k)
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
           dl(i, j) * dl(k, l) + dl(i, k) * dl(j, l) + dl(i, l) * dl(j, k)
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
.admAdfWeightFast <- function(C, w, Dv, N) {
  w  <- w / sum(w)
  m  <- ncol(C)
  ij <- which(lower.tri(diag(m), diag = TRUE), arr.ind = TRUE)
  q  <- nrow(ij)
  I  <- ij[, 1L]; J <- ij[, 2L]
  dbar <- colSums(w * Dv)
  S    <- crossprod(C, w * C); diag(S) <- diag(S) + dbar

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
  list(E = m$mu, V = m$V, C = sm$cpc, w = g$W / sum(g$W),
       Dv = .admAdfCondVar(cp, arr))
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
    W  <- .admAdfWeightFast(pt$C, pt$w, pt$Dv, as.numeric(s$n))
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
