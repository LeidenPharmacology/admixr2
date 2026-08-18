## ============================================================================
## covariate-integration-comparison-v2.R -- CORRECTED re-run
## ============================================================================
##
## WHAT WAS WRONG IN THE ORIGINAL
## ------------------------------
## The original compared methods against THREE "targets":
##
##   (A) -2LL at the MARGINAL moments        <- this is CONSTRUCTION 1, correct
##   (B) -2 log E_a[L(a)]  (mixture)         <- marginal obs vs conditional preds,
##                                              combined as a mixture. A valid
##                                              likelihood only for a study whose
##                                              subjects ALL share one unknown a,
##                                              which is not how the data here (or
##                                              a published summary) are generated.
##   (C) E_a[-2 log L(a)]                    <- CONSTRUCTION 4, the bug: it scores
##                                              a MARGINAL observation against a
##                                              CONDITIONAL prediction at each node
##                                              and averages the resulting -2LLs.
##
## Every "method" in the original (taylor2/taylor4, gl-n, IS-reweight) was then
## scored on how well it approximated (C). They approximate the category error
## well; that says nothing about integrating a covariate out of an aggregate
## likelihood. Section 5's h-sweep is the differencing step OF (C), so its numbers
## are the step-size sensitivity of the wrong functional.
##
## WHAT THIS VERSION DOES
##   - target (A) is renamed `gh_marginal` and is the reference throughout;
##   - a genuine STRATIFIED arm (construction 2) is added: per-node observed
##     blocks with n_k = w_k*n scored against per-node predictions;
##   - `taylor_moments` (construction 3) replaces taylor2/taylor4 as the cheap
##     approximation, and its h-sweep replaces section 5's;
##   - `gh_moments-K` replaces `gl-K`: the same K nodes, but integrating the
##     MOMENTS and evaluating the likelihood once;
##   - IS reweighting is evaluated BOTH ways: the ESS mechanics (unchanged, they
##     do not depend on what is combined) and a corrected IS-moments arm;
##   - the `s`-collapse check (section 4) is UNCHANGED -- it is an identity
##     between the marginal moments and nested quadrature and never touched
##     construction 4;
##   - (C) is retained under the name `mismatched` purely to size the error.
## ============================================================================

## ---- model: 1-cmt IV bolus, CL mu-referenced -------------------------------
## CL = exp(tcl + th_cov * a + b) = exp(tcl + s),   s = th_cov * a + b
tcl <- 1.0; Vd <- 10; dose <- 100
times <- c(0.5, 1, 2, 4, 8)

fpred <- function(s) dose / Vd * exp(outer(-exp(tcl + s) / Vd, times))

gh <- function(n) {
  i <- seq_len(n - 1L)
  J <- matrix(0, n, n); J[cbind(i, i + 1L)] <- sqrt(i); J[cbind(i + 1L, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE); o <- order(e$values)
  list(x = e$values[o], w = e$vectors[1L, o]^2)
}
GH <- gh(80L)

moments_s <- function(m, sd, q = GH) {
  Fm <- fpred(m + sd * q$x)
  mu <- as.numeric(crossprod(q$w, Fm))
  Fc <- sweep(Fm, 2L, mu)
  list(mu = mu, V = t(Fc) %*% (Fc * q$w))
}
nll2 <- function(E_obs, V_obs, mu, Vp, n) {
  ch <- tryCatch(chol(Vp), error = function(e) NULL)
  if (is.null(ch)) return(Inf)
  iv <- chol2inv(ch); r <- E_obs - mu
  n * (2 * sum(log(diag(ch))) + sum(iv * V_obs) + as.numeric(t(r) %*% iv %*% r))
}

mu_a <- 0; sd_a <- 1; om <- 0.30; n_sub <- 200L; sig2 <- 0.01
th_true <- 0.75
sd_s <- function(tc) sqrt(tc^2 * sd_a^2 + om^2)

## moments CONDITIONAL on covariate value a (eta integrated out), + residual
mom_at <- function(a, tc) {
  m <- moments_s(tc * a, om); Vp <- m$V; diag(Vp) <- diag(Vp) + sig2
  list(E = m$mu, V = Vp)
}
## CONSTRUCTION 1 -- exact marginal moments (via the s-collapse; section 4 proves
## it equals nested quadrature to 1e-14)
mom_marg <- function(tc) {
  m <- moments_s(tc * mu_a, sd_s(tc)); Vp <- m$V; diag(Vp) <- diag(Vp) + sig2
  list(E = m$mu, V = Vp)
}
## CONSTRUCTION 1 with K covariate nodes, done as nested quadrature over MOMENTS
mom_marg_K <- function(tc, K) {
  q <- gh(K); E <- 0; M2 <- 0; Vc <- 0
  for (i in seq_along(q$x)) {
    a <- mu_a + sd_a * q$x[i]; m <- moments_s(tc * a, om)
    E  <- E  + q$w[i] * m$mu
    Vc <- Vc + q$w[i] * m$V
    M2 <- M2 + q$w[i] * outer(m$mu, m$mu)
  }
  V <- Vc + M2 - outer(E, E); diag(V) <- diag(V) + sig2
  list(E = E, V = V)
}
## CONSTRUCTION 3 -- 2nd-order expansion of the MOMENTS
mom_taylor <- function(tc, hfrac = 0.5) {
  h <- hfrac * sd_a
  m0 <- mom_at(mu_a, tc); mp <- mom_at(mu_a + h, tc); mm <- mom_at(mu_a - h, tc)
  dE  <- (mp$E - mm$E) / (2 * h)
  d2E <- (mp$E - 2 * m0$E + mm$E) / h^2
  d2V <- (mp$V - 2 * m0$V + mm$V) / h^2
  list(E = m0$E + 0.5 * sd_a^2 * d2E,
       V = m0$V + 0.5 * sd_a^2 * d2V + sd_a^2 * outer(dE, dE))
}

make_obs <- function(tc) mom_marg(tc)                       # pooled/marginal data
make_obs_nodes <- function(tc, K = 9L) {                    # per-node data
  q <- gh(K)
  list(q = q, blocks = lapply(mu_a + sd_a * q$x, mom_at, tc = tc))
}
obs   <- make_obs(th_true)
obsN  <- make_obs_nodes(th_true)

## ---- objectives -------------------------------------------------------------
obj_gh_marginal <- function(tc, obs) {                       # [1]
  m <- mom_marg(tc); nll2(obs$E, obs$V, m$E, m$V, n_sub) }
obj_taylor_moments <- function(tc, obs, hfrac = 0.5) {       # [3]
  m <- mom_taylor(tc, hfrac); nll2(obs$E, obs$V, m$E, m$V, n_sub) }
obj_stratified <- function(tc, obsN) {                       # [2]
  q <- obsN$q; aa <- mu_a + sd_a * q$x
  sum(vapply(seq_along(aa), function(k) {
    m <- mom_at(aa[k], tc)
    nll2(obsN$blocks[[k]]$E, obsN$blocks[[k]]$V, m$E, m$V, q$w[k] * n_sub) }, numeric(1))) }
nll_at_a <- function(a, tc, obs) {                           # marginal obs vs cond pred
  m <- mom_at(a, tc); nll2(obs$E, obs$V, m$E, m$V, n_sub) }
obj_mismatched <- function(tc, obs, q = GH)                  # [4] the bug (old target C)
  sum(q$w * vapply(mu_a + sd_a * q$x, nll_at_a, numeric(1), tc, obs))
obj_mixture <- function(tc, obs, q = GH) {                   # old target B
  l <- vapply(mu_a + sd_a * q$x, nll_at_a, numeric(1), tc, obs)
  mn <- min(l); mn - 2 * log(sum(q$w * exp(-(l - mn) / 2))) }

hr <- function(s) cat("\n", strrep("=", 76), "\n", s, "\n", strrep("=", 76), "\n", sep = "")

## ---------------------------------------------------------------------------
hr("1. The constructions, and their minimisers")
cat("Data generated at the TRUE marginal moments with theta_cov = 0.75.\n")
cat("(One population, mu_a = 0: tcov and omega are only jointly identified through\n")
cat(" tcov^2*sd_a^2 + omega^2, so each column below profiles the other at truth.)\n\n")
cat(sprintf("%10s %14s %14s %14s %14s %14s\n", "theta_cov",
            "[1] gh_marg", "[2] stratif", "[3] taylor", "[4] MISMATCH", "mixture"))
for (tc in c(0.60, 0.70, 0.75, 0.80, 0.90))
  cat(sprintf("%10.2f %14.3f %14.3f %14.3f %14.3f %14.3f\n", tc,
              obj_gh_marginal(tc, obs), obj_stratified(tc, obsN),
              obj_taylor_moments(tc, obs), obj_mismatched(tc, obs), obj_mixture(tc, obs)))

cat("\nMinimiser of each (this is the theta_cov you would report):\n")
FNS <- list("[1] gh_marginal   (construction 1)" = function(tc) obj_gh_marginal(tc, obs),
            "[2] stratified    (construction 2)" = function(tc) obj_stratified(tc, obsN),
            "[3] taylor_moments(construction 3)" = function(tc) obj_taylor_moments(tc, obs),
            "[4] MISMATCHED    (construction 4)" = function(tc) obj_mismatched(tc, obs),
            "    mixture       (old target B)  " = function(tc) obj_mixture(tc, obs))
for (nm in names(FNS)) {
  o <- optimize(FNS[[nm]], c(0.2, 2.0), tol = 1e-8)
  cat(sprintf("   %-36s -> theta_cov = %.4f   (truth 0.7500, bias %+.4f, %+.1f%%)\n",
              nm, o$minimum, o$minimum - th_true, 100*(o$minimum - th_true)/th_true))
}

cat("\nOmega recovery, profiling omega with theta_cov fixed at truth:\n")
OFNS <- list("[1] gh_marginal" = function(tc) obj_gh_marginal(tc, obs),
             "[2] stratified"  = function(tc) obj_stratified(tc, obsN),
             "[3] taylor"      = function(tc) obj_taylor_moments(tc, obs),
             "[4] MISMATCHED"  = function(tc) obj_mismatched(tc, obs))
for (nm in names(OFNS)) {
  f <- OFNS[[nm]]
  o <- optimize(function(o_try) { om <<- o_try; v <- f(th_true); om <<- 0.30; v },
                c(0.05, 1.2), tol = 1e-8)
  cat(sprintf("   %-16s -> omega = %.4f   (truth 0.3000, bias %+.4f)\n",
              nm, o$minimum, o$minimum - 0.30))
}
cat("\nNOTE: [2]'s per-node observed blocks are regenerated at each candidate omega\n")
cat("      only through the prediction; the obs blocks are fixed at truth, as for\n")
cat("      every other arm.\n")

## ---------------------------------------------------------------------------
hr("2. How well does each cheap method hit CONSTRUCTION 1 (the right target)?")
cat("Reference = [1] gh_marginal at 80 nodes (exact, via the s-collapse).\n\n")
cat(sprintf("%10s %14s %14s %14s %14s %14s\n", "theta_cov", "[1] exact",
            "gh_mom-5", "gh_mom-9", "[3] taylor", "[4] MISMATCH"))
for (tc in c(0.50, 0.75, 1.00, 1.50)) {
  ex <- obj_gh_marginal(tc, obs)
  m5 <- mom_marg_K(tc, 5L);  n5 <- nll2(obs$E, obs$V, m5$E, m5$V, n_sub)
  m9 <- mom_marg_K(tc, 9L);  n9 <- nll2(obs$E, obs$V, m9$E, m9$V, n_sub)
  cat(sprintf("%10.2f %14.3f %14.3f %14.3f %14.3f %14.3f\n", tc, ex, n5, n9,
              obj_taylor_moments(tc, obs), obj_mismatched(tc, obs)))
}
cat("\nabsolute error against [1]:\n")
cat(sprintf("%10s %14s %14s %14s %14s\n", "theta_cov", "gh_mom-5", "gh_mom-9",
            "[3] taylor", "[4] MISMATCH"))
for (tc in c(0.50, 0.75, 1.00, 1.50)) {
  ex <- obj_gh_marginal(tc, obs)
  m5 <- mom_marg_K(tc, 5L); m9 <- mom_marg_K(tc, 9L)
  cat(sprintf("%10.2f %14.3e %14.3e %14.3e %14.3e\n", tc,
              abs(nll2(obs$E,obs$V,m5$E,m5$V,n_sub) - ex),
              abs(nll2(obs$E,obs$V,m9$E,m9$V,n_sub) - ex),
              abs(obj_taylor_moments(tc, obs) - ex),
              abs(obj_mismatched(tc, obs) - ex)))
}

## ---------------------------------------------------------------------------
hr("3. IS reweighting: ESS mechanics (UNCHANGED) + a corrected IS-moments arm")
cat("The ESS numbers below do not depend on what is combined afterwards, so they\n")
cat("carry over from the original run untouched. What IS corrected is the arm that\n")
cat("uses the reweighted pool: it now produces per-node MOMENTS which are combined\n")
cat("into marginal moments (construction 1), instead of per-node -2LLs.\n\n")
m_is <- function(tc, obs, n_nodes = 9L, n_pool = 20000L) {
  b  <- qnorm((seq_len(n_pool) - 0.5) / n_pool) * om
  Fm <- fpred(b)
  q  <- gh(n_nodes); nodes <- mu_a + sd_a * q$x
  ess <- numeric(n_nodes); nl <- numeric(n_nodes)
  E <- 0; Vc <- 0; M2 <- 0
  for (i in seq_along(nodes)) {
    D <- tc * nodes[i]
    w <- dnorm(b, D, om) / dnorm(b, 0, om); w <- w / sum(w)
    ess[i] <- 1 / sum(w^2) / n_pool
    mu <- as.numeric(crossprod(w, Fm)); Fc <- sweep(Fm, 2L, mu)
    Vp <- t(Fc) %*% (Fc * w)
    Vd <- Vp; diag(Vd) <- diag(Vd) + sig2
    nl[i] <- nll2(obs$E, obs$V, mu, Vd, n_sub)     # the OLD (construction 4) use
    E  <- E  + q$w[i] * mu
    Vc <- Vc + q$w[i] * Vp
    M2 <- M2 + q$w[i] * outer(mu, mu)
  }
  Vm <- Vc + M2 - outer(E, E); diag(Vm) <- diag(Vm) + sig2
  list(mismatched = sum(q$w * nl), moments = nll2(obs$E, obs$V, E, Vm, n_sub), ess = ess)
}
cat(sprintf("%12s %10s %14s %14s %12s %12s\n",
            "shift/omega", "theta_cov", "[1] exact", "IS-moments", "abs err", "min ESS/N"))
for (tc in c(0.05, 0.10, 0.20, 0.35, 0.50, 0.75)) {
  r <- m_is(tc, obs); ex <- obj_gh_marginal(tc, obs)
  cat(sprintf("%12.2f %10.2f %14.3f %14.3f %12.4f %12.2e\n",
              tc * 2 * sd_a / om, tc, ex, r$moments, abs(r$moments - ex), min(r$ess)))
}
cat("\nFor one node at shift D the weight is exp(D'Omega^-1 b - D'Omega^-1 D/2),\n")
cat("so ESS/N decays as exp(-(D/omega)^2)  [UNCHANGED from the original run]:\n")
for (d in c(0.5, 1, 2, 3, 4))
  cat(sprintf("     D/omega = %.1f -> ESS/N = %.2e  (of 20000 draws: %.1f effective)\n",
              d, exp(-d^2), 20000 * exp(-d^2)))

## ---------------------------------------------------------------------------
hr("4. s-collapse: is the nested integral needed at all?  [UNCHANGED]")
cat("This section never used construction 4. It is an identity between the\n")
cat("marginal moments and nested 2-D quadrature, and is reproduced verbatim.\n\n")
cat("Under mu-referencing with a LINEAR covariate effect and a normal covariate,\n")
cat("   s = theta_cov*a + b ~ N(theta_cov*mu_a, theta_cov^2*sd_a^2 + omega^2)\n")
cat("EXACTLY.\n\n")
m1 <- moments_s(th_true * mu_a, sd_s(th_true))
q <- gh(60L); acc_mu <- 0; acc_M2 <- 0
for (ia in seq_along(q$x)) {
  mm <- moments_s(th_true * (mu_a + sd_a * q$x[ia]), om)
  acc_mu <- acc_mu + q$w[ia] * mm$mu
  acc_M2 <- acc_M2 + q$w[ia] * (mm$V + outer(mm$mu, mm$mu))
}
acc_V <- acc_M2 - outer(acc_mu, acc_mu)
cat(sprintf("   max |rel diff| vs nested 2-D quadrature:   E %.2e    V %.2e\n",
            max(abs(m1$mu / acc_mu - 1)), max(abs(m1$V / acc_V - 1))))

## ---------------------------------------------------------------------------
hr("5. The differencing step h for the CORRECTED (moment) Taylor expansion")
cat("The original section 5 swept h for the expansion of construction 4 -- the\n")
cat("step size of the wrong functional. Redone here for construction 3, with h\n")
cat("expressed as a fraction of sd_a (the original's h was in ABSOLUTE covariate\n")
cat("units, which is a separate defect and is why h = 2.0 landed arbitrarily).\n\n")
relerr <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-12))
ex_m <- mom_marg(th_true); ex_n <- obj_gh_marginal(th_true, obs)
cat(sprintf("%10s %12s %16s %16s %16s\n", "hfrac", "h (abs)", "rel err E", "rel err V", "NLL err vs [1]"))
for (hf in c(0.05, 0.1, 0.25, 0.5, 0.75, 1.0, 2.0)) {
  ta <- mom_taylor(th_true, hf)
  cat(sprintf("%10.2f %12.2f %16.3e %16.3e %16.4f\n", hf, hf * sd_a,
              relerr(ta$E, ex_m$E), relerr(ta$V, ex_m$V),
              obj_taylor_moments(th_true, obs, hf) - ex_n))
}
cat("\nFor comparison, the OLD (construction 4) h-sweep, error in the MISMATCHED\n")
cat("target -- retained only to identify what those registry numbers measured:\n")
m_taylor_old <- function(tc, obs, h = 2.0, ord = 2L) {
  g <- function(k) nll_at_a(mu_a + k * h, tc, obs)
  r <- g(0) + 0.5 * sd_a^2 * (g(1) - 2 * g(0) + g(-1)) / h^2
  if (ord >= 4L) r <- r + sd_a^4 / 8 * (g(2) - 4*g(1) + 6*g(0) - 4*g(-1) + g(-2)) / h^4
  r
}
tr <- obj_mismatched(th_true, obs)
cat(sprintf("%10s %14s %14s\n", "h", "old taylor2", "old taylor4"))
for (h in c(0.05, 0.2, 0.5, 1.0, 2.0, 3.0))
  cat(sprintf("%10.2f %14.5f %14.5f\n", h,
              m_taylor_old(th_true, obs, h = h, ord = 2L) - tr,
              m_taylor_old(th_true, obs, h = h, ord = 4L) - tr))

## ---------------------------------------------------------------------------
hr("6. Regime sweep -- where is each CORRECT method adequate?")
cat("Ratio = theta_cov*sd_a / omega. Data regenerated at the true marginal for\n")
cat("each ratio, then evaluated there. Relative errors are against [1].\n")
cat("'mismatch size' is |[4] - [1]| / |[1]|: how large the category error is.\n\n")
cat(sprintf("%7s %9s %13s %13s %13s %13s %13s\n",
            "ratio", "theta_cov", "gh_mom-9", "[3] taylor", "IS-moments",
            "[2] stratif", "mismatch size"))
for (r in c(0.1, 0.25, 0.5, 1.0, 2.0, 3.0)) {
  tc <- r * om
  o2 <- make_obs(tc); o2N <- make_obs_nodes(tc)
  t1 <- obj_gh_marginal(tc, o2)
  m9 <- mom_marg_K(tc, 9L); n9 <- nll2(o2$E, o2$V, m9$E, m9$V, n_sub)
  t3 <- obj_taylor_moments(tc, o2)
  is <- m_is(tc, o2)$moments
  st <- obj_stratified(tc, o2N)
  t4 <- obj_mismatched(tc, o2)
  cat(sprintf("%7.2f %9.3f %13.2e %13.2e %13.2e %13.2e %13.2e\n", r, tc,
              abs(n9 - t1)/abs(t1), abs(t3 - t1)/abs(t1), abs(is - t1)/abs(t1),
              abs(st - t1)/abs(t1), abs(t4 - t1)/abs(t1)))
}
cat("\n([2] stratified is a DIFFERENT likelihood on DIFFERENT data -- it is not an\n")
cat(" approximation to [1], so its column is a magnitude, not an error.)\n\n")
