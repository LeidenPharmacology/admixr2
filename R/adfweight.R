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
  # The conditional MEAN, filled only by the forms whose mean is NONLINEAR in the
  # structural prediction -- which is TBS and nothing else. Everywhere else the
  # mean is ms * f for a constant ms (f for combined1/2 and pois, f exp(sv/2) for
  # lnorm, size * f for binom, mu for nbinomMu and beta), so centring the nodes
  # and scaling by ms is EXACT and .admAdfParts keeps doing that. Leaving m1 NULL
  # there is what makes this change bit-identical for every existing family
  # rather than a re-derivation of six conditional means that could disagree with
  # the objective's.
  m1 <- NULL
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
      "3" = {                                          # TBS, by quadrature
        # boxCox / yeoJohnson / logitNorm / probitNorm. Conditionally independent
        # across timepoints like every other form here, so the expansion applies;
        # what it needed was the third and fourth central moments, which
        # .admTBSCentral() takes off the SAME node set the objective's mean and
        # variance come from.
        #
        # The sd is derived by .admTBSSd(), the function .admTBSRow() uses, so the
        # law conditioned on here is the law the objective composes V_pred from
        # rather than a second reading of `ftr`/`c1`/`cc`.
        # Accessed exactly as .admResidDeriv() accesses them: single-bracket, and
        # `tbs_ftr`/`tbs_c1` guarded for NULL, which they are on a row that never
        # went through the TBS builder. col() would have errored on NULL[[j]].
        # t() folds nu/(nu-2) into a2/b2, which is exact for the combined forms
        # because only the VARIANCE of the residual enters there. It is not exact
        # here: this branch integrates g() over the conditional law, so a
        # t-distributed error is not an inflated-sd normal one, and the third and
        # fourth moments this returns would be a normal's. The combined branch
        # above computes the t kurtosis explicitly; there is no such closed form
        # once the transform is applied, so refuse and let the fit report "r".
        # The OBJECTIVE keeps composing it as an inflated-sd normal, which is
        # what it has always done -- see .admTBSRow().
        vmj <- col(arr$vmul, j)
        if (!is.null(vmj) && !isTRUE(all.equal(vmj, 1))) return(NULL)
        lam <- arr$lam[j]; yjc <- arr$yj[j]
        lo  <- arr$tlo[j]; hi  <- arr$thi[j]
        ftr <- !is.null(arr$tbs_ftr) && isTRUE(arr$tbs_ftr[j])
        c1  <- !is.null(arr$tbs_c1)  && isTRUE(arr$tbs_c1[j])
        if (length(lam) != 1L || length(yjc) != 1L ||
            !is.finite(lam) || !is.finite(yjc)) return(NULL)
        sdv <- .admTBSSd(f, a2, b2, cc, lam, yjc, lo, hi, ftr, c1)$sdv
        if (any(!is.finite(sdv))) return(NULL)
        tm <- .admTBSCentral(f, sdv, lam, yjc, lo, hi,
                             arr$nodes %||% .ADM_TBS_NODES)
        d[, j] <- tm$v; t3[, j] <- tm$mu3; q4[, j] <- tm$mu4
        # E[y | node] is a nonlinear function of f here, so the linearisation
        # ms * (f - fbar) that serves every other family understates the spread of
        # the conditional means and Var(ybar) with it. Carry the exact means.
        if (is.null(m1)) m1 <- matrix(NA_real_, Q, m)
        m1[, j] <- tm$m
      },
      # form 8 (ordinal) is a joint unit, refused upstream.
      return(NULL))
  }
  if (!all(is.finite(d)) || !all(is.finite(t3)) || !all(is.finite(q4)))
    return(NULL)
  # A unit mixing a TBS endpoint with a closed-form one would leave m1 partly NA.
  # Fall back to the linearisation for the whole unit rather than mix the two.
  if (!is.null(m1) && !all(is.finite(m1))) m1 <- NULL
  list(d = d, t3 = t3, q4 = q4, m1 = m1)
}

# Which timepoints the OBJECTIVE's residual variance is exact at.
#
# pow()/combined() compose E_eta[Var(y|eta)] as a2 + b2 E[f^2c] (plus a cross term
# for combined1), and .admMomF supplies E[f^k] from (mu, var_f) by a second-order
# delta expansion that is EXACT only at k = 1 and k = 2. Everywhere else it is a
# truncation -- reached by pow()/combined() with an exponent outside {0.5, 1}, and
# by nothing else.
#
# .admAdfCondMom, by contrast, evaluates b2 |f|^2c AT EVERY NODE, so its node
# average is the exact integral. Both are defensible readings of the same model
# and they are not equal: measured 9e-05 relative at c = 0.75 with omega = 0.16,
# rising to 2.5e-02 at c = 1.5 with omega = 1 (against 8e-16 at c = 0.5 and c = 1).
# See .admAdfAlignDv for why that gap has to be closed rather than tolerated.
.admAdfMomExact <- function(arr, m) {
  col <- function(x, j) if (length(x) == 1L) x else x[[j]]
  vapply(seq_len(m), function(j) {
    fm <- arr$form[[j]]
    if (fm != .ADM_RESID_COMBINED2 && fm != .ADM_RESID_COMBINED1) return(TRUE)
    a2 <- col(arr$a2, j); b2 <- col(arr$b2, j); cc <- col(arr$cc, j)
    # b2 == 0 is add(): there is no E[f^k] to get wrong.
    if (!is.finite(b2) || b2 == 0) return(TRUE)
    if (!is.finite(cc)) return(FALSE)
    # c == 1 (prop/add+prop) needs only E[f^2] = mu^2 + var_f, exact.
    if (isTRUE(all.equal(cc, 1))) return(TRUE)
    # c == 0.5 needs only E[f], also exact -- unless the form is combined1, whose
    # cross term 2 sqrt(a2 b2) E[f^c] then asks for E[f^0.5], which is not. With
    # a2 == 0 that term is absent and the row is exact again.
    isTRUE(all.equal(cc, 0.5)) &&
      (fm == .ADM_RESID_COMBINED2 || !is.finite(a2) || a2 == 0)
  }, logical(1))
}

# Make the weight's implied residual variance the one the objective composes.
#
# S is rebuilt from (C, Dv) by the law of total variance, and the whole expansion
# downstream is scaled by it -- so if S disagrees with the V_pred that G is the
# cross-derivative of, the weight and the objective describe different laws and
# J = 2H fails at a correctly-specified fit. On a pow(c = 1.5) model that showed
# up as a reported "r,s" correction which was nothing but .admMomF's truncation.
#
# THE OBJECTIVE IS WHAT OMEGA HAS TO DESCRIBE, as it is for TBS one branch up: G
# is a derivative of the criterion that was actually minimised, so reading the
# moments from a second, better composition breaks the information equality by
# exactly the difference between the two. The node integral being the more
# accurate of the pair does not make it the right one here.
#
# Rescale each affected column's CONDITIONAL variance by the single factor that
# puts its node average on the objective's E[Var(y|eta)]. The correction is not
# applied at all where .admMomF is exact, which keeps every other family
# bit-for-bit what it was.
#
# MULTIPLICATIVE, NOT A SHIFT. The two differ where it matters: b2 |f|^2c at the
# outer nodes of a pow() grid runs down to 4e-08 while the gap to close is 2e-04,
# so the constant that fixes the average drives those nodes' variance NEGATIVE and
# the whole correction has to be abandoned. A factor cannot -- it preserves sign
# and relative spread across nodes, which is what the third and fourth moments are
# built out of, and it is the natural reading of a variance anyway.
.admAdfAlignDv <- function(cm, w, ev, arr) {
  if (is.null(cm) || is.null(cm$d) || is.null(ev)) return(cm)
  ex <- .admAdfMomExact(arr, ncol(cm$d))
  if (all(ex)) return(cm)
  jj <- which(!ex)
  d0 <- cm$d[, jj, drop = FALSE]
  r  <- ev[jj] / as.numeric(crossprod(w, d0))
  # A non-positive or non-finite ratio means the expansion and the node integral
  # do not describe the same law at all -- .admMomF's near-zero cap can put ev
  # below zero where the node integral cannot. Drop Dv and let the caller degrade
  # to "r" rather than weight by it.
  if (!all(is.finite(r)) || any(r <= 0)) { cm$d <- NULL; return(cm) }
  # The combined branch of .admAdfCondMom leaves t3 = 0 and sets q4 = kappa d^2
  # for a CONSTANT kappa (3, or the Student-t kurtosis), so scaling q4 by r^2
  # carries the same kappa through and is all the fourth moment needs.
  cm$q4[, jj] <- sweep(cm$q4[, jj, drop = FALSE], 2L, r^2, "*")
  cm$d[, jj]  <- sweep(d0, 2L, r, "*")
  cm
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
#
# `var_only`: a `method = "var"` study only ever reports the mean block plus
# the DIAGONAL of the covariance-summary block, so building the full
# q = m(m+1)/2 vech and then discarding the off-diagonal q - m of it is pure
# waste -- and a quadratic one: T1/M3/M4 below are q x q, so at m = 200 the
# full path allocates multiple ~3 GB matrices for a result that ends up
# 2m x 2m. Restricting `ij` to the diagonal up front makes q = m and every
# line below -- written generically in terms of `ij`/`I`/`J`/`q` -- comes out
# already in the shape a `var` study needs, in the same ascending-i order the
# old post-hoc trim (.admWeightSel, now folded into this switch) produced.
.admAdfWeightFast <- function(C, w, Dv, N, T3 = NULL, Q4 = NULL,
                              var_only = FALSE) {
  w  <- w / sum(w)
  m  <- ncol(C)
  ij <- if (var_only) cbind(seq_len(m), seq_len(m))
        else which(lower.tri(diag(m), diag = TRUE), arr.ind = TRUE)
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
  # One crossprod, NOT a vapply over the columns of wP. At m == 1 there is a
  # single vech entry and vapply(..., numeric(1)) returns a VECTOR, so the
  # M3[, b] indexing below threw "incorrect number of dimensions" -- which the
  # drivers' tryCatch swallowed, so every SINGLE-TIMEPOINT study warned that the
  # sandwich could not be computed and silently reported "r". The reference
  # .admAdfWeight() never had the problem, and the pinning tests all use m >= 3.
  CC <- crossprod(C, wP)                                # m x q
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
#
# `g`: .adghGrid(pars, pinfo, grid) depends on pars/grid alone, not on `study`,
# so it is identical across every study in one .admSandwichCov() call -- a
# caller looping over studies with the same `pars`/`grid` can build it ONCE and
# pass it in, rather than have every call recompute the same small matmul.
# Optional and NULL by default so every other caller (the tests, the default
# .admMomentDeriv() mom closure) is unaffected.
.admAdfParts <- function(pars, pinfo, study, rxMod, out_var, grid, cores, g = NULL) {
  g     <- g %||% .adghGrid(pars, pinfo, grid)
  pm    <- .admMakeParamsList(nrow(g$eta), pinfo, 1L)[[1L]]
  cp    <- .admSimulate(rxMod, pars$struct, pinfo$sigma_names, g$eta, study,
                        out_var, pm, cores, pinfo$nDisplayProgress, pinfo$sigdig)
  sm  <- .adghStructMoments(cp, g$W)
  arr <- .admUnitResidRows(pinfo, out_var, pars$sigma_var, length(sm$mu),
                           phi = attr(cp, "phi"))
  # A TBS unit composes at the NODES, exactly as the objective does. G has to
  # describe the objective, so reading the moments from a different composition
  # than .adghNLL uses breaks the information equality by precisely the
  # expansion's truncation -- measured as eigen(J/2H) drifting to 0.987/1.007 on
  # logitNorm while boxCox, whose truncation is ~1e-06, still looked fine.
  m   <- .admResidNodeMomentsTBS(cp, g$W, arr, study$times)
  m   <- if (!is.null(m)) list(mu = m$E, V = m$V, ms = NULL)
         else .admResidMoments(sm$mu, diag(sm$V), arr, sm$V, study$times)
  # C carries the residual's MEAN SCALING. On an lnorm endpoint the conditional
  # mean is f exp(s/2), not f, so an unscaled C makes the weight's own S differ
  # from the V the objective scores against -- the two would then describe
  # different laws. ms is 1 on every other family, so this is a no-op there.
  #
  # A TBS endpoint returns the conditional means THEMSELVES (cm$m1), because
  # there the mean is a nonlinear function of f and ms * (f - fbar) is only its
  # linearisation. Centring is by the SAME normalised node weights the rest of
  # the weight uses, so Var over nodes of the conditional mean comes out exact.
  cm  <- .admAdfCondMom(cp, arr)
  wn  <- g$W / sum(g$W)
  # Put the weight's residual variance on the objective's composition wherever
  # .admMomF's expansion is a truncation (pow()/combined() outside c in {0.5, 1});
  # a no-op for every other family. `m$ev` is NULL on the TBS branch, which
  # composes at the nodes already and is exact by construction.
  cm  <- .admAdfAlignDv(cm, wn, m$ev, arr)
  Cm  <- if (!is.null(cm$m1)) sweep(cm$m1, 2L, as.numeric(crossprod(wn, cm$m1)))
         else if (is.null(m$ms)) sm$cpc
         else sweep(sm$cpc, 2L, m$ms, "*")
  list(E = m$mu, V = m$V, w = wn,
       C  = Cm, Dv = cm$d, T3 = cm$t3, Q4 = cm$q4)
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

.admSandwich <- function(H, G, Om, Hinv = NULL) {
  # `Hinv` lets the caller thread through the inverse it already computed for
  # the "r" leg (chol2inv(chol(H)), with its own fallbacks) instead of paying
  # for a second O(P^3) factorisation of the identical H here.
  Hi <- Hinv %||% tryCatch(solve(H), error = function(e) NULL)
  if (is.null(Hi)) return(NULL)
  p <- nrow(H); J <- matrix(0, p, p)
  for (i in seq_along(G)) J <- J + G[[i]] %*% Om[[i]] %*% t(G[[i]])
  list(cov = Hi %*% J %*% Hi, bread = 2 * Hi, J = J, H = H)
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

# d(yt)/dPsi and d(Vt)/dPsi per study, by central difference on the MOMENTS.
#
# These are analytic from what the gradient machinery already forms; this is the
# reference until that extraction is written. Differencing the MOMENTS rather
# than the objective keeps it well conditioned, and it runs once, post-fit.
#
# `h` is only the FALLBACK step now, used verbatim where it always was. When
# `nll_fn` -- the estimator's own objective, already in scope at every call site
# for the Hessian FD -- is supplied, the step is measured per parameter by
# .admShi21Steps() instead, the same convention every other finite difference in
# the package follows (NEWS.md: "not optional"). This differences the MOMENTS,
# not the objective, but the objective's noise level is still the right proxy
# for how finely `p` can be perturbed before the perturbation is swamped by
# solver noise -- and measuring against a fixed 1e-5 regardless of parameter
# scale is exactly the failure mode Shi21 exists to avoid.
.admMomentDeriv <- function(p_hat, pinfo, studies, rxMod, out_var, grid, cores,
                            h = 1e-5, mom_fn = NULL, nll_fn = NULL) {
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
  hv <- if (is.null(nll_fn)) rep(h, p) else
    tryCatch(.admShi21Steps(nll_fn, p_hat, seq_len(p), fallback = rep(h, p),
                            .var.name = "admMomentDeriv"),
             error = function(e) rep(h, p))
  b  <- mom(p_hat)
  dE <- lapply(b, function(x) matrix(0, length(x$E), p))
  dV <- lapply(b, function(x) rep(list(matrix(0, nrow(x$V), ncol(x$V))), p))
  for (k in seq_len(p)) {
    hk <- hv[[k]]
    a  <- p_hat; a[k]  <- a[k] + hk
    cc <- p_hat; cc[k] <- cc[k] - hk
    ma <- mom(a); mc <- mom(cc)
    for (i in seq_along(b)) {
      dE[[i]][, k] <- (ma[[i]]$E - mc[[i]]$E) / (2 * hk)
      dV[[i]][[k]] <- (ma[[i]]$V - mc[[i]]$V) / (2 * hk)
    }
  }
  list(E = lapply(b, `[[`, "E"), V = lapply(b, `[[`, "V"), dE = dE, dV = dV)
}

# The conditioning bound the SANDWICH needs, which is not the one a single
# inversion needs.
#
# covreport.R calls H singular below .ADM_NPD_RCOND = sqrt(eps): one inversion
# loses about kappa * eps of relative accuracy, so that is where 2H^-1 is deemed
# hopeless. The sandwich inverts H TWICE -- H^-1 J H^-1 -- so its effective
# conditioning is kappa^2, and it reaches that same bound already at
# kappa = eps^(-1/4), i.e. rcond = eps^(1/4) ~ 1.2e-04 (cond ~ 8200).
#
# Between the two thresholds is a band where "r" is usable and "r,s" is not, and
# nothing else in the package looks at it. It is not a rounding problem: what is
# amplified is the genuine gap between J and 2H in a direction the data barely
# identifies, which is exactly where the two differ most and mean least. Measured
# on a 1-cmt fixture whose residual SD contributed 0.01 variance against 1.7 from
# IIV (RSE 186%, cond(H) = 3.5e05), the reported residual SE moved by 0.115 and
# two omega entries by 0.59 and 1.55; the same model and design with the residual
# identified (cond(H) = 247) reproduced "r" to four decimals throughout.
#
# WARN RATHER THAN DEGRADE. The number is not garbage in the way a non-finite or
# non-PD one is, and the well-determined directions of the same fit are fine --
# so withholding the whole covariance would cost more than it saves. What the
# user cannot do is NOTICE, which is what this fixes: the warning names the
# direction, so a surprising SE can be read against the parameter that caused it.
.ADM_SANDWICH_RCOND <- .Machine$double.eps^(1/4)

# RETURNS the diagnosis, it does not raise it. A warning() from here would be
# swallowed: this runs inside the nlmixr2est stack, which is why the driver
# raises "covariance could not be computed" itself rather than letting
# .admCalcCov do it, and why an incomplete source covariance once cost a fit its
# sandwich in silence. The message travels out on an attribute and each driver
# emits it beside the covMethod label, where it reaches the user.
#
# `eig_dec`: every driver already has `eigen(H, symmetric = TRUE)` in hand --
# .admReduceNpdOmega() built it deciding whether to drop the omega block, on
# the SAME (possibly-reduced) H this function is handed -- so recomputing an
# O(P^3) eigendecomposition of the identical matrix here on every "r,s" fit is
# pure waste. Supplied, it is used as-is; NULL (a caller with no eigendecomp on
# hand, e.g. a test) falls back to computing it, unchanged from before.
.admSandwichCond <- function(H, nms = NULL, eig_dec = NULL) {
  e <- eig_dec %||% tryCatch(eigen(H, symmetric = TRUE), error = function(e) NULL)
  if (is.null(e)) return(NULL)
  ev <- e$values
  mx <- max(abs(ev))
  if (!is.finite(mx) || mx <= 0) return(NULL)
  rc <- min(abs(ev)) / mx
  if (!is.finite(rc) || rc >= .ADM_SANDWICH_RCOND) return(NULL)
  # The offending direction is the eigenvector of the smallest eigenvalue; report
  # the parameter loading on it most heavily. That is the actionable half -- a
  # condition number alone does not tell anyone which SE to distrust.
  k  <- which.min(abs(ev))
  ld <- abs(e$vectors[, k])
  nm <- (nms %||% rownames(H) %||% paste0("p", seq_along(ev)))[which.max(ld)]
  sprintf(
    paste0("covMethod = \"r,s\": the Hessian is ill-conditioned (rcond %.1e, ",
           "cond %.3g), and the sandwich inverts it twice where \"r\" inverts it ",
           "once -- so the correction is amplified quadratically in the ",
           "weakly-identified direction, which loads mainly on `%s`. Check that ",
           "parameter's relative standard error before reading its \"r,s\" value ",
           "as a finding; the well-determined parameters are unaffected."),
    rc, 1 / rc, nm)
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
                            H, md = NULL, keep = NULL, mom_fn = NULL,
                            sensModel = NULL, nms = NULL, Hinv = NULL,
                            nll_fn = NULL, eig_dec = NULL) {
  # H is checked FIRST, and here rather than being left to solve() inside .admSandwich,
  # whose tryCatch is there for a singular matrix and cannot tell that apart from
  # an H that was never supplied. A caller that forgot the argument then gets a
  # silent NULL -- which reads as "the sandwich does not apply to this fit"
  # rather than "you called it wrong", and a calibration study built on it
  # reported numbers for replicates it had skipped entirely.
  if (missing(H) || !is.matrix(H) || nrow(H) != ncol(H) || !all(is.finite(H)))
    stop(".admSandwichCov: `H` must be a finite square Hessian of the objective ",
         "at the optimum -- the same one covMethod = \"r\" inverts.", call. = FALSE)
  # Diagnosed here, where H is, and carried out on an attribute of the RESULT --
  # so a sandwich that goes on to degrade for an unrelated reason returns NULL
  # and takes the conditioning message with it, rather than warning about a
  # covariance the fit never reported.
  .cond <- .admSandwichCond(H, nms, eig_dec = eig_dec)
  pars <- tryCatch(.admUnpack(p_hat, pinfo), error = function(e) NULL)
  if (is.null(pars)) return(NULL)
  # A joint (same-subject, multi-output) unit stacks several outputs into one
  # covariance. .admAdfParts solves ONE output, so its node ensemble is not the
  # one behind that block and the Wick expansion would be taken over the wrong
  # conditional law. Refuse rather than return a plausible wrong weight.
  if (any(vapply(studies, function(s) isTRUE(s$is_joint), logical(1))))
    return(NULL)
  # Analytic first, finite differences as the fallback. .admMomentJac returns
  # NULL rather than an approximation for any path it does not cover, so the
  # switch is on availability, not on a tolerance.
  #
  # The two differ by ~7e-5 on an ODE model and NOT as a function of the FD step,
  # so that gap is not truncation: the analytic route reads its predictions from
  # the SENSITIVITY model and the FD route from the plain simulation model, and
  # two separately compiled models take different adaptive steps. The sens model
  # is the right one here -- covMethod = "r"'s Hessian is built from the
  # analytic gradient, which reads the same model, so G and H stay consistent.
  if (is.null(md) && is.null(mom_fn))
    md <- tryCatch(.admMomentJac(p_hat, pinfo, studies, sensModel, rxMod,
                                 out_var, grid, cores), error = function(e) NULL)
  md <- md %||% tryCatch(
    .admMomentDeriv(p_hat, pinfo, studies, rxMod, out_var, grid, cores,
                    mom_fn = mom_fn, nll_fn = nll_fn),
    error = function(e) NULL)
  if (is.null(md)) return(NULL)
  G <- Om <- vector("list", length(studies))
  # .adghGrid(pars, pinfo, grid) is invariant across this loop (same pars, same
  # grid for every study), so it is built once here rather than once per study
  # inside .admAdfParts().
  g_node <- .adghGrid(pars, pinfo, grid)
  for (i in seq_along(studies)) {
    s  <- studies[[i]]
    pt <- tryCatch(.admAdfParts(pars, pinfo, s, rxMod, s$output %||% out_var,
                                grid, cores, g = g_node), error = function(e) NULL)
    if (is.null(pt) || is.null(pt$Dv)) return(NULL)
    N <- as.numeric(s$n)
    Om[[i]] <- .admAdfWeightFast(pt$C, pt$w, pt$Dv, N, pt$T3, pt$Q4,
                                 var_only = identical(s$method, "var"))
    G[[i]]  <- .admScoreCross(md$E[[i]], md$V[[i]], md$dE[[i]], md$dV[[i]], s, N)
    if (is.null(G[[i]])) return(NULL)
    # H may have been reduced to the struct+sigma sub-block (.admReduceNpdOmega),
    # so G must lose the same rows -- G is parameter-indexed by row.
    if (!is.null(keep)) G[[i]] <- G[[i]][keep, , drop = FALSE]
  }
  if (!is.null(keep) && nrow(H) != length(keep)) return(NULL)
  out <- .admSandwich(H, G, Om, Hinv = Hinv)
  if (!is.null(out)) attr(out, "illcond") <- .cond
  out
}

# Why the sandwich DOES NOT APPLY to this model, or NULL when it might.
#
# These are refusals by construction, not failures: the Wick expansion needs a
# residual that is independent across timepoints given the node and has a finite
# fourth moment, and it needs one output per unit. An ar() residual, an ordinal
# endpoint, a joint unit and a t with nu <= 4 each break one of those, and no
# amount of re-conditioning or re-solving would change the answer.
#
# Telling them apart from a failure MATTERS because "r,s" is the default: every
# ar() and joint fit was emitting "the sandwich correction could not be
# computed" -- a sentence shaped like something went wrong -- and nlmixr2est
# carries that onto fit$runInfo, where it reads as a defect in the fit. For these
# models the reason below is reported instead, once, as a message.
#
# Deliberately limited to causes that are cheap and CERTAIN to detect up front. A
# refusal decided deeper in (a non-finite TBS sd, a singular weight, a failed
# moment solve) still comes back as NULL and still warns -- that one is a failure
# of something that was expected to work.
.admSandwichNA <- function(p_hat, pinfo, studies, out_var) {
  if (any(vapply(studies, function(s) isTRUE(s$is_joint), logical(1))))
    return(paste("the study stacks several outputs per subject (a joint unit),",
                 "whose conditional law the weight's node ensemble does not",
                 "describe"))
  # .admSandwichGrid() refuses rather than degrade to a coarser-than-floor grid
  # once nq^n_eta cannot be covered at nq = 3 (n_eta >= 8, see its own comment).
  # Checked here, up front, for the same reason as every other cause below: the
  # driver's tryCatch(..., error = NULL) around the grid builder cannot tell
  # this refusal apart from an actual failure, so without this check it fell
  # through to the generic "could not be computed" warning instead.
  if (is.null(.admSandwichGrid(pinfo)))
    return(sprintf(paste("the model has %d random effects, too many for the",
                         "sandwich's capped product quadrature grid (floored",
                         "at 3 nodes/eta, capped at 5000 nodes total) to cover"),
                   pinfo$n_eta %||% NA_integer_))
  sv <- tryCatch(.admUnpack(p_hat, pinfo)$sigma_var, error = function(e) NULL)
  if (is.null(sv)) return(NULL)
  for (s in studies) {
    n_t <- length(s$times %||% integer(0))
    if (n_t == 0L) next
    arr <- tryCatch(.admResidRows(pinfo, s$output %||% out_var, sv, n_t),
                    error = function(e) NULL)
    if (is.null(arr)) next
    if (!is.null(arr$rho) && any(!is.na(arr$rho)))
      return(paste("the residual is autocorrelated (ar()), so the products the",
                   "expansion factors across timepoints do not factor"))
    if (any(arr$form == .ADM_RESID_ORDINAL))
      return("the endpoint is ordinal, which is scored as a joint unit")
    vm <- arr$vmul %||% 1
    nu <- ifelse(is.finite(vm) & vm > 1, 2 * vm / (vm - 1), Inf)
    if (any(nu <= 4))
      return(paste("the residual is Student-t with nu <= 4, which has no finite",
                   "kurtosis -- the sampling law of the reported V then has no",
                   "variance for the weight to be made of"))
    # TBS's conditional mean is nonlinear in f, so .admAdfCondMom integrates the
    # residual's law over the node ensemble rather than folding nu/(nu-2) into a
    # variance the way the combined forms do -- and there is no closed form for
    # that integral under a t rather than a normal law, at any nu. See the TBS
    # branch of .admAdfCondMom, which refuses on exactly this condition.
    if (any(arr$form == .ADM_RESID_TBS & is.finite(vm) &
            abs(vm - 1) > sqrt(.Machine$double.eps)))
      return(paste("the endpoint composes a transform-both-sides law",
                   "(boxCox/yeoJohnson/logitNorm/probitNorm) with a Student-t",
                   "residual; the transformed residual's third and fourth",
                   "moments have no closed form under a t rather than a normal",
                   "law, at any nu -- the objective keeps composing it as an",
                   "inflated-sd normal and the weight cannot"))
  }
  NULL
}

# Validates a candidate sandwich result and folds it into the *CalcCov "r"
# baseline, or falls back to that baseline with a warning. `label` names the
# caller (e.g. "adghCalcCov") for the fallback warning.
#
# `na` is .admSandwichNA()'s reason, for a model the correction does not apply to
# at all: the fallback is then reported as a MESSAGE naming the reason, rather
# than as a warning about a computation that failed.
#
# Shared by adgh/admc/adfo's *CalcCov -- what differs between them is how `sw`
# is BUILT (adgh's own grid + sensModel; admc's quadrature grid + sensModel;
# adfo's quadrature grid + moment map), not what happens to it once built.
#
# The acceptance gate is a full PSD check, not a diagonal one. J = sum(G Om G')
# is only guaranteed PSD if every per-study Om is -- and .admAdfAlignDv's
# rescaling (pow()/combined() with an exponent outside {0.5, 1}) is not itself
# checked for that, so a positive diagonal does not imply a valid covariance.
# Reuses .ADM_NPD_RCOND (covreport.R), the same reciprocal-condition-number
# tolerance the Hessian's own singularity test uses, so "PSD enough to trust"
# means the same thing on both sides of the sandwich.
.admIsPsd <- function(m) {
  e <- tryCatch(eigen(m, symmetric = TRUE, only.values = TRUE)$values,
               error = function(e) NULL)
  if (is.null(e) || !length(e) || any(!is.finite(e))) return(FALSE)
  mx <- max(abs(e))
  if (!is.finite(mx) || mx <= 0) return(FALSE)
  min(e) >= -.ADM_NPD_RCOND * mx
}

.admApplySandwich <- function(sw, cov_r, label, na = NULL) {
  ok <- !is.null(sw) && all(is.finite(sw$cov)) && .admIsPsd(sw$cov)
  if (!ok) {
    if (!is.null(na))
      message(sprintf(paste("  covMethod = \"r,s\" does not apply to this model:",
                            "%s. Reporting the covMethod = \"r\" covariance."), na))
    else
      warning(sprintf(paste("%s: the sandwich correction could not be computed;",
                            "reporting the covMethod = \"r\" covariance instead."),
                      label), call. = FALSE)
    return(list(cov_full = cov_r, sw_used = FALSE, sw_cond = NULL))
  }
  list(cov_full = (sw$cov + t(sw$cov)) / 2, sw_used = TRUE, sw_cond = attr(sw, "illcond"))
}

# Wires the sandwich into a *CalcCov driver: refuse via .admSandwichNA(), build
# `sw` from `build()` inside a tryCatch, and fold the result into `cov_r` via
# .admApplySandwich(). Shared because this ~15-line block was copy-pasted
# near-verbatim across adfo/adgh/admc's *CalcCov -- what differs between them is
# how the sandwich's own ingredients (grid, sensModel or mom_fn) get built, which
# is exactly what `build` (a zero-arg closure over the caller's locals) captures.
# Returns list(cov_full, sw_used, sw_cond), ready to assign back at the call site.
.admWireSandwich <- function(p_hat, pinfo, studies, out_var, cov_r, label, build) {
  # A model the correction does not APPLY to (ar(), ordinal, a joint unit, a t
  # with nu <= 4, TBS + t() at any nu, 8+ random effects) is reported as such
  # and not attempted -- see .admSandwichNA.
  .sw_na <- tryCatch(.admSandwichNA(p_hat, pinfo, studies, out_var),
                     error = function(e) NULL)
  sw <- if (!is.null(.sw_na)) NULL else tryCatch(build(), error = function(e) NULL)
  .admApplySandwich(sw, cov_r, label, na = .sw_na)
}

# Finalises the covariance a driver reports: warns once if none could be
# computed, re-raises any ill-conditioning note the sandwich attached (as a
# warning, so nlmixr2est carries it onto fit$runInfo -- see the call sites),
# and returns the covMethod label the fit should report ("r,s" / "r" / "").
.admFinalizeCovLabel <- function(cov, want_cov) {
  if (isTRUE(want_cov) && is.null(cov))
    warning("covariance could not be computed (the Hessian was singular or ",
            "non-finite); standard errors are unavailable for this fit.",
            call. = FALSE)
  if (!is.null(sw_cond <- attr(cov, "sandwich_illcond")))
    warning(sw_cond, call. = FALSE)
  if (is.null(cov)) "" else if (isTRUE(attr(cov, "sandwich"))) "r,s" else "r"
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
#
# n_eta == 0 IS A GRID, not a refusal. .adghNodeGrid() returns the single-point
# ensemble there, which is the correct one: with no between-subject variability
# every subject shares the structural prediction, and the summary's sampling law
# is the residual's alone. That law is still not the normal-theory one the
# objective assumes -- lnorm, pois, binom, beta and the TBS family are all skewed
# or over-dispersed conditional on the prediction -- so the correction still has
# something to say. Returning NULL here made admc and adfo degrade to "r" on
# every no-IIV model while adgh, which passes its own grid, applied the sandwich
# to the same fit; the two disagreed for no reason but this line.
.admSandwichGrid <- function(pinfo, max_nodes = 5000L) {
  n_eta <- pinfo$n_eta
  if (is.null(n_eta) || n_eta < 0L) return(NULL)
  if (n_eta == 0L) return(.adghNodeGrid(1L, 0L))
  nq <- 9L
  while (nq > 3L && nq^n_eta > max_nodes) nq <- nq - 2L
  # THREE IS THE FLOOR, and past it the answer is NULL rather than a smaller
  # grid. The decrement used to run down to nq = 1, which for n_eta >= 8 it
  # reached (3^8 = 6561 > 5000): .adghNodes1(1) is the single node at eta = 0, so
  # C centres to identically zero, S collapses to diag(E[Var(y|eta)]), and the
  # weight then describes a model with NO between-subject variability. That
  # result is finite with a positive diagonal, so .admApplySandwich() accepts it
  # and the fit reports covMethod = "r,s" over silently wrong SEs -- the one
  # failure mode worse than not applying the correction at all. Both callers
  # stop() on a NULL grid inside a tryCatch and degrade to "r", which is the
  # honest answer for a model with more etas than a product grid can cover.
  if (nq^n_eta > max_nodes) return(NULL)
  .adghNodeGrid(nq, n_eta)
}

# adfo's own (E, V) per study -- the moment map the FO objective minimises.
#
# Returned in the shape .admMomentDeriv expects, so the only thing that changes
# between estimators is this function. `V` is the full predicted covariance even
# for a `method = "var"` study: .admScoreCross takes the diagonal itself, which
# keeps the branch logic in one place rather than two.
#
# .admMomentDeriv's central-difference loop calls this 2p times (once per
# +h/-h perturbation of every struct/sigma/omega parameter), each of which
# used to cost .adfoGetMuJ() -- and so its own rxSolve -- per study: the "each
# configuration cost its own rxSolve" anti-pattern .adfoGetMuJBatch exists to
# collapse, reintroduced here by a loop that .adfoGetMuJBatch itself cannot
# see into. But .adfoGetMuJ's result depends on pp ONLY through pars$struct
# and, for a TBS endpoint whose lambda is an estimated sigma parameter, that
# one sigma entry (.adfoMuJKey; see its own comment for why -- the ONLY place
# a residual parameter feeds the structural solve). Every sigma/omega
# direction this FD loop perturbs OTHER than that leaves both inputs at the
# base point, so memoizing on .adfoMuJKey via .adfoMuJMemo -- exactly what
# .adfoNLL()/.adfoGrad() already do for the same reason -- collapses those
# directions' solves to one, reused, rather than reissuing them.
.admAdfoMomFn <- function(pinfo, studies, sensModel, rxMod, out_var,
                          params_list, cores) {
  cache <- new.env(parent = emptyenv())
  function(pp) {
    pars <- .admUnpack(pp, pinfo)
    key  <- .adfoMuJKey(pars, sensModel)
    lapply(seq_along(studies), function(i) {
      s   <- studies[[i]]
      ov  <- s$output %||% out_var
      n_t <- length(s$times)
      arr <- .admResidRows(pinfo, ov, pars$sigma_var, n_t)
      mj  <- .adfoMuJMemo(cache, i, key, function()
        .adfoGetMuJ(pars, pinfo, s, sensModel, rxMod, ov,
                    params_list[[i]], cores))
      vp  <- .adfoVpred(mj$mu, mj$J, pars$L, arr, n_t, pinfo$n_eta, s$times)
      list(E = vp$mu_sigma, V = vp$V)
    })
  }
}

# -- Analytic moment Jacobian --------------------------------------------------

# d(yt)/dPsi and d(Vt)/dPsi per study, analytically, from one sensitivity solve.
#
# G = d2F/(dPsi dt') is closed form in these two (see .admScoreCross), so this is
# the only place a derivative is taken at all. Both moments are LINEAR in the raw
# sensitivity column `graw = d(f)/dPsi`:
#
#   d(mu_struct)/dPsi = W' graw
#   d(V_struct)/dPsi  = A + A',   A = cpc' diag(W) graw
#
# which means a parameter reached through SEVERAL paths is handled by summing
# its columns before this is applied, exactly as .adghGradNLL sums its
# `contrib()` calls. That linearity is why this does not need to know which path
# a parameter took.
#
# The residual composition is then applied FORWARD, from the same
# .admResidDeriv() partials the gradient chains BACKWARD:
#
#   dE      = dmu_df o dmu_s + dmu_dv0 o diag(dV_s) + dmu %*% dsig
#   dV_ij   = ms_i ms_j dV_s_ij + (dms_i ms_j + ms_i dms_j) cov_f_ij   (i != j)
#   dV_ii   = dv_dv0_i dV_s_ii + dv_df_i dmu_s_i + dvar_i . dsig
#
# Deriving this tail by hand is the seventh consumer of the residual row arrays,
# and CLAUDE.md is explicit that the moment tail is where the misses happen --
# so it is pinned against the finite-difference version (.admMomentDeriv) across
# error models rather than trusted. FD stays as that oracle.
#
# Returns NULL, not an approximation, whenever a path is not covered: no sens
# model, a joint unit, unpaired thetas without their own columns, or an
# ar()/ordinal residual whose rmat carries an off-diagonal this forward map does
# not model. The caller falls back to FD.
.admMomentJac <- function(p_hat, pinfo, studies, sensModel, rxMod, out_var, grid,
                          cores) {
  if (is.null(sensModel)) return(NULL)
  pars  <- tryCatch(.admUnpack(p_hat, pinfo), error = function(e) NULL)
  if (is.null(pars) || !.admParsFinite(pars, pinfo)) return(NULL)
  n_s   <- length(pinfo$struct_names)
  n_e   <- length(pinfo$sigma_names)
  n_eta <- pinfo$n_eta
  L     <- pars$L
  np    <- length(p_hat)
  unpaired_k <- if (!is.null(pinfo$struct_has_eta))
    which(!pinfo$struct_has_eta) else integer(0)

  Eo <- Vo <- dEo <- dVo <- vector("list", length(studies))

  # The node grid depends on Omega alone, so it is the same for every study --
  # the same grid .adghNLL()/.adghGradNLL() build once per evaluation.
  gS  <- .adghGrid(pars, pinfo, grid)
  X   <- grid$X
  W   <- gS$W
  eta <- gS$eta; colnames(eta) <- pinfo$eta_col_names

  for (si in seq_along(studies)) {
    s <- studies[[si]]
    if (isTRUE(s$is_joint)) return(NULL)
    ov <- s$output %||% out_var

    res <- .admSimulateSens(sensModel, pars$struct, pinfo$sigma_names, eta, s,
                            cores, pinfo$nDisplayProgress, pars$sigma_var,
                            pinfo$sigdig)
    if (is.null(res)) return(NULL)
    if (length(unpaired_k) > 0L && is.null(res$dtheta_list)) return(NULL)
    Jl <- res$dpred_list

    sm    <- .adghStructMoments(res$cp_mat, W)
    mu    <- sm$mu; cpc <- sm$cpc; cov_f <- sm$V; var_f <- diag(cov_f)
    m     <- length(mu)
    arr   <- .admResidRows(pinfo, ov, pars$sigma_var, m)

    # ---- TBS: the node-wise chain, matching the objective --------------------
    .npj <- .admTBSNodeParts(res$cp_mat, arr)
    if (!is.null(.npj)) {
      .agj <- .admTBSAggregate(.npj, W)
      wnj <- .agj$wn; Mcj <- .agj$Mc
      dEi <- matrix(0, m, np)
      dVi <- rep(list(matrix(0, m, m)), np)
      .cnj <- function(A, dvn) {
        de <- as.numeric(crossprod(wnj, A))
        Ac <- sweep(A, 2L, de)
        Bc <- crossprod(Ac, wnj * Mcj)
        dV <- Bc + t(Bc)
        diag(dV) <- diag(dV) + as.numeric(crossprod(wnj, dvn))
        list(dE = de, dV = dV)
      }
      .cgj <- function(graw) .cnj(.npj$dm * graw, .npj$dv * graw)
      for (k in seq_len(n_s)) {
        graw <- if (is.null(pinfo$struct_has_eta) || pinfo$struct_has_eta[k]) {
          ei <- which(pinfo$struct_eta_idx == k)[1L]
          if (is.na(ei)) NULL else Jl[[ei]]
        } else res$dtheta_list[[pinfo$struct_names[k]]]
        if (is.null(graw)) return(NULL)
        r <- .cgj(graw); dEi[, k] <- r$dE; dVi[[k]] <- r$dV
      }
      if (n_e > 0L) for (k in seq_len(n_e)) {
        hk <- max(abs(p_hat[n_s + k]), 1) * 1e-5
        mv <- lapply(c(1, -1), function(sgn) {
          pk <- p_hat; pk[n_s + k] <- pk[n_s + k] + sgn * hk
          pp <- tryCatch(.admUnpack(pk, pinfo), error = function(e) NULL)
          if (is.null(pp)) return(NULL)
          .admTBSNodeParts(res$cp_mat, .admResidRows(pinfo, ov, pp$sigma_var, m))
        })
        if (any(vapply(mv, is.null, TRUE))) return(NULL)
        r <- .cnj((mv[[1L]]$m - mv[[2L]]$m) / (2 * hk),
                  (mv[[1L]]$v - mv[[2L]]$v) / (2 * hk))
        dEi[, n_s + k] <- r$dE; dVi[[n_s + k]] <- r$dV
      }
      if (n_eta > 0L) for (rr in seq_along(pinfo$omega_par)) {
        i <- pinfo$chol_i[rr]; j <- pinfo$chol_j[rr]
        r  <- .cgj(Jl[[i]] * X[, j])
        sc <- if (pinfo$chol_diag[rr]) L[i, i] / 2 else 1
        pos <- n_s + n_e + rr
        dEi[, pos] <- r$dE * sc; dVi[[pos]] <- r$dV * sc
      }
      Eo[[si]] <- .agj$E; Vo[[si]] <- .agj$V
      dEo[[si]] <- dEi;   dVo[[si]] <- dVi
      next
    }
    # ---- end TBS ------------------------------------------------------------

    pmres <- .admResidMoments(mu, var_f, arr, cov_f, s$times)
    # ar()/ordinal put an off-diagonal residual term in rmat whose sigma and mu
    # paths this forward map does not carry. .admAdfCondMom refuses those
    # residuals anyway; refusing here too keeps the two boundaries identical.
    if (!is.null(pmres$rmat) && any(pmres$rmat != 0, na.rm = TRUE)) return(NULL)
    dres  <- .admResidDeriv(mu, var_f, arr, pinfo)

    # structural moments from a raw sensitivity column
    mom <- function(graw) {
      dmu <- as.numeric(crossprod(W, graw))
      A   <- crossprod(cpc, W * graw)
      list(dmu = dmu, dV = A + t(A))
    }
    # forward residual composition
    ms  <- dres$ms
    tailf <- function(dmu_s, dV_s, dsig = NULL) {
      dms <- dres$dms_df * dmu_s
      dE  <- dres$dmu_df * dmu_s + dres$dmu_dv0 * diag(dV_s)
      dgd <- dres$dv_dv0 * diag(dV_s) + dres$dv_df * dmu_s
      if (!is.null(dsig)) {
        dms <- dms + as.numeric(dres$dms %*% dsig)
        dE  <- dE  + as.numeric(dres$dmu %*% dsig)
        dgd <- dgd + as.numeric(dres$dvar %*% dsig)
      }
      dV <- outer(ms, ms) * dV_s +
            (outer(dms, ms) + outer(ms, dms)) * cov_f
      diag(dV) <- dgd
      list(dE = dE, dV = dV)
    }

    dEi <- matrix(0, m, np)
    dVi <- rep(list(matrix(0, m, m)), np)
    zero <- matrix(0, nrow(res$cp_mat), m)

    # --- structural thetas: every path that reaches f, summed --------------
    for (k in seq_len(n_s)) {
      graw <- zero
      if (is.null(pinfo$struct_has_eta) || pinfo$struct_has_eta[k]) {
        ei <- which(pinfo$struct_eta_idx == k)[1L]
        if (!is.na(ei)) graw <- graw + Jl[[ei]]
      } else {
        Dt <- res$dtheta_list[[pinfo$struct_names[k]]]
        if (is.null(Dt)) return(NULL)
        graw <- graw + Dt
      }
      mm <- mom(graw); tf <- tailf(mm$dmu, mm$dV)
      dEi[, k] <- tf$dE; dVi[[k]] <- tf$dV
    }

    # --- sigma: no path through f, only the residual composition -----------
    for (k in seq_len(n_e)) {
      ds <- numeric(n_e); ds[k] <- 1
      tf <- tailf(numeric(m), matrix(0, m, m), ds)
      dEi[, n_s + k] <- tf$dE; dVi[[n_s + k]] <- tf$dV
    }

    # --- omega Cholesky ----------------------------------------------------
    if (n_eta > 0L) for (rr in seq_along(pinfo$omega_par)) {
      i <- pinfo$chol_i[rr]; j <- pinfo$chol_j[rr]
      # d(eta[q,])/d(L_ij) = X[q,j] e_i, so d(f[q,])/d(L_ij) = Jl[[i]][q,] X[q,j]
      # -- the same base .adghGradNLL() contracts.
      base <- Jl[[i]] * X[, j]
      mm  <- mom(base); tf <- tailf(mm$dmu, mm$dV)
      sc  <- if (pinfo$chol_diag[rr]) L[i, i] / 2 else 1
      pos <- n_s + n_e + rr
      dEi[, pos] <- tf$dE * sc; dVi[[pos]] <- tf$dV * sc
    }

    Eo[[si]] <- pmres$mu; Vo[[si]] <- pmres$V
    dEo[[si]] <- dEi;     dVo[[si]] <- dVi
  }
  list(E = Eo, V = Vo, dE = dEo, dV = dVo)
}
