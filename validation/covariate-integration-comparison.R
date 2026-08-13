## ============================================================================
## Covariate marginalisation for AGGREGATE data -- comparing the options
## ============================================================================
##
## Self-contained: no admixr2, no rxode2, no packages. Run with
##   Rscript validation/covariate-integration-comparison.R
##
## The model has a closed form, so "truth" here is exact quadrature, not a
## bigger Monte Carlo. Everything below is 1-D Gauss-Hermite at 80 nodes.
##
## Methods compared
##   taylor2 / taylor4  what .adm_combine_nll(method = "taylor") does today
##   gl-n               what .adm_combine_nll(gl/gh) does today: sum_i w_i NLL_i
##   IS-reweight        the one-simulation-pool importance-weighting proposal
##   s-collapse         exact, and free, under mu-referencing + normal covariate
##
## ...against THREE targets, because these methods do not all estimate the same
## quantity, and for aggregate data only one of them is the right one:
##   (A) -2LL at the MARGINAL moments E_{a,b}[f], Cov_{a,b}(f)
##   (B) -2 log E_a[ L(a) ]   mixture likelihood: study has ONE unknown a
##   (C) E_a[ -2 log L(a) ]   average of per-stratum -2LLs

## ---- model: 1-cmt IV bolus, CL mu-referenced -------------------------------
## CL = exp(tcl + th_cov * a + b) = exp(tcl + s),   s = th_cov * a + b
tcl <- 1.0; Vd <- 10; dose <- 100
times <- c(0.5, 1, 2, 4, 8)

fpred <- function(s)                        # length(s) x length(times)
  dose / Vd * exp(outer(-exp(tcl + s) / Vd, times))

## ---- Gauss-Hermite for N(0,1) (Golub-Welsch, probabilists') ----------------
gh <- function(n) {
  i <- seq_len(n - 1L)
  J <- matrix(0, n, n); J[cbind(i, i + 1L)] <- sqrt(i); J[cbind(i + 1L, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE); o <- order(e$values)
  list(x = e$values[o], w = e$vectors[1L, o]^2)          # weights sum to 1
}
GH <- gh(80L)

## moments of f under s ~ N(m, sd^2)
moments_s <- function(m, sd, q = GH) {
  Fm <- fpred(m + sd * q$x)
  mu <- as.numeric(crossprod(q$w, Fm))
  Fc <- sweep(Fm, 2L, mu)
  list(mu = mu, V = t(Fc) %*% (Fc * q$w))
}

## ---- the aggregate -2LL admixr2 minimises ----------------------------------
nll2 <- function(E_obs, V_obs, mu, Vp, n) {
  ch <- tryCatch(chol(Vp), error = function(e) NULL)
  if (is.null(ch)) return(Inf)
  iv <- chol2inv(ch); r <- E_obs - mu
  n * (2 * sum(log(diag(ch))) + sum(iv * V_obs) + as.numeric(t(r) %*% iv %*% r))
}

## ---- scenario --------------------------------------------------------------
## Covariate on the log scale, so sd_a = 1 makes "h" and "h in SDs" coincide
## (section 5 breaks that on purpose).
mu_a <- 0; sd_a <- 1; om <- 0.30; n_sub <- 200L; sig2 <- 0.01
th_true <- 0.75

sd_s <- function(tc) sqrt(tc^2 * sd_a^2 + om^2)

make_obs <- function(tc) {                  # data AT the true marginal moments
  m <- moments_s(tc * mu_a, sd_s(tc))
  Vp <- m$V; diag(Vp) <- diag(Vp) + sig2
  list(E = m$mu, V = Vp)
}
obs <- make_obs(th_true)

## ---- the three targets -----------------------------------------------------
target_A <- function(tc, obs) {             # marginal moments (aggregate-correct)
  m <- moments_s(tc * mu_a, sd_s(tc))
  Vp <- m$V; diag(Vp) <- diag(Vp) + sig2
  nll2(obs$E, obs$V, m$mu, Vp, n_sub)
}
nll_at_a <- function(a, tc, obs) {          # per-stratum -2LL, b integrated out
  m <- moments_s(tc * a, om)
  Vp <- m$V; diag(Vp) <- diag(Vp) + sig2
  nll2(obs$E, obs$V, m$mu, Vp, n_sub)
}
target_C <- function(tc, obs, q = GH)
  sum(q$w * vapply(mu_a + sd_a * q$x, nll_at_a, numeric(1), tc, obs))
target_B <- function(tc, obs, q = GH) {
  l <- vapply(mu_a + sd_a * q$x, nll_at_a, numeric(1), tc, obs)
  mn <- min(l)
  mn - 2 * log(sum(q$w * exp(-(l - mn) / 2)))
}

## ---- methods ---------------------------------------------------------------
## Taylor exactly as .adm_combine_nll: nodes at mu +- h (ABSOLUTE covariate
## units), combined with Sigma_cov[j,j] * second-difference / h^2.
m_taylor <- function(tc, obs, h = 2.0, ord = 2L) {
  g <- function(k) nll_at_a(mu_a + k * h, tc, obs)
  r <- g(0) + 0.5 * sd_a^2 * (g(1) - 2 * g(0) + g(-1)) / h^2
  if (ord >= 4L)
    r <- r + sd_a^4 / 8 * (g(2) - 4*g(1) + 6*g(0) - 4*g(-1) + g(-2)) / h^4
  r
}
m_gl <- function(tc, obs, n_nodes = 9L) {
  q <- gh(n_nodes)
  sum(q$w * vapply(mu_a + sd_a * q$x, nll_at_a, numeric(1), tc, obs))
}
## The proposal: ONE pool of eta draws with NO covariate shift, reweighted to
## each covariate node by  w = N(b; th*a_i, Omega) / N(b; 0, Omega).
m_is <- function(tc, obs, n_nodes = 9L, n_pool = 20000L) {
  b  <- qnorm((seq_len(n_pool) - 0.5) / n_pool) * om     # deterministic pool
  Fm <- fpred(b)
  q  <- gh(n_nodes); nodes <- mu_a + sd_a * q$x
  ess <- nl <- numeric(n_nodes)
  for (i in seq_along(nodes)) {
    D <- tc * nodes[i]
    w <- dnorm(b, D, om) / dnorm(b, 0, om); w <- w / sum(w)
    ess[i] <- 1 / sum(w^2) / n_pool
    mu <- as.numeric(crossprod(w, Fm)); Fc <- sweep(Fm, 2L, mu)
    Vp <- t(Fc) %*% (Fc * w); diag(Vp) <- diag(Vp) + sig2
    nl[i] <- nll2(obs$E, obs$V, mu, Vp, n_sub)
  }
  list(nll = sum(q$w * nl), ess = ess)
}

hr <- function(s) cat("\n", strrep("=", 76), "\n", s, "\n", strrep("=", 76), "\n", sep = "")

## ---------------------------------------------------------------------------
hr("1. The three targets are different functions, with different minimisers")
cat("Data generated at the TRUE marginal moments with theta_cov = 0.75.\n\n")
cat(sprintf("%10s %14s %14s %14s\n", "theta_cov", "(A) marginal", "(B) mixture", "(C) E[NLL]"))
for (tc in c(0.60, 0.70, 0.75, 0.80, 0.90))
  cat(sprintf("%10.2f %14.3f %14.3f %14.3f\n",
              tc, target_A(tc, obs), target_B(tc, obs), target_C(tc, obs)))
cat("\nMinimiser of each (this is the theta_cov you would report):\n")
for (nm in c("A", "B", "C")) {
  f <- get(paste0("target_", nm))
  o <- optimize(function(tc) f(tc, obs), c(0.2, 2.0), tol = 1e-8)
  cat(sprintf("   target %s -> theta_cov = %.4f   (truth 0.7500, bias %+.4f)\n",
              nm, o$minimum, o$minimum - th_true))
}
cat("\nOmega recovery, profiling omega with theta_cov fixed at truth:\n")
for (nm in c("A", "C")) {
  f <- get(paste0("target_", nm))
  o <- optimize(function(o_try) { om <<- o_try; v <- f(th_true, obs); om <<- 0.30; v },
                c(0.05, 1.2), tol = 1e-8)
  cat(sprintf("   target %s -> omega = %.4f   (truth 0.3000, bias %+.4f)\n",
              nm, o$minimum, o$minimum - 0.30))
}

## ---------------------------------------------------------------------------
hr("2. How well does each method hit target (C), the thing it aims at?")
cat(sprintf("%10s %12s %12s %12s %12s\n",
            "theta_cov", "truth (C)", "taylor2", "taylor4", "gl-9"))
for (tc in c(0.50, 0.75, 1.00, 1.50))
  cat(sprintf("%10.2f %12.3f %12.3f %12.3f %12.3f\n", tc, target_C(tc, obs),
              m_taylor(tc, obs, ord = 2L), m_taylor(tc, obs, ord = 4L), m_gl(tc, obs)))

## ---------------------------------------------------------------------------
hr("3. The IS-reweighting proposal: accuracy and effective sample size")
cat("One pool of 20000 eta draws, reweighted to 9 covariate nodes.\n\n")
cat(sprintf("%12s %10s %12s %12s %11s %12s\n",
            "shift/omega", "theta_cov", "truth (C)", "IS-reweight", "abs err", "min ESS/N"))
for (tc in c(0.05, 0.10, 0.20, 0.35, 0.50, 0.75)) {
  r <- m_is(tc, obs); tr <- target_C(tc, obs)
  cat(sprintf("%12.2f %10.2f %12.3f %12.3f %11.4f %12.2e\n",
              tc * 2 * sd_a / om, tc, tr, r$nll, abs(r$nll - tr), min(r$ess)))
}
cat("\nFor one node at shift D the weight is exp(D'Omega^-1 b - D'Omega^-1 D/2),\n")
cat("so ESS/N decays as exp(-(D/omega)^2):\n")
for (d in c(0.5, 1, 2, 3, 4))
  cat(sprintf("     D/omega = %.1f -> ESS/N = %.2e  (of 20000 draws: %.1f effective)\n",
              d, exp(-d^2), 20000 * exp(-d^2)))

## ---------------------------------------------------------------------------
hr("4. s-collapse: is the nested integral needed at all?")
cat("Under mu-referencing with a LINEAR covariate effect and a normal covariate,\n")
cat("   s = theta_cov*a + b ~ N(theta_cov*mu_a, theta_cov^2*sd_a^2 + omega^2)\n")
cat("EXACTLY. So the marginal moments need one 1-D integral, no reweighting,\n")
cat("no nodes, and no extra ODE solves over a plain no-covariate fit.\n\n")
m1 <- moments_s(th_true * mu_a, sd_s(th_true))
q <- gh(60L); acc_mu <- 0; acc_M2 <- 0
for (ia in seq_along(q$x)) {                         # brute-force nested 2-D
  mm <- moments_s(th_true * (mu_a + sd_a * q$x[ia]), om)
  acc_mu <- acc_mu + q$w[ia] * mm$mu
  acc_M2 <- acc_M2 + q$w[ia] * (mm$V + outer(mm$mu, mm$mu))
}
acc_V <- acc_M2 - outer(acc_mu, acc_mu)
cat(sprintf("   max |rel diff| vs nested 2-D quadrature:   E %.2e    V %.2e\n",
            max(abs(m1$mu / acc_mu - 1)), max(abs(m1$V / acc_V - 1))))
cat("   solve cost: 1 set of nodes, vs n_a x n_b nested.\n")

## ---------------------------------------------------------------------------
hr("5. The Taylor step h is in ABSOLUTE covariate units, not SDs")
cat("wt_nodes = mu +- h, combined as Sigma_cov[j,j] * d2NLL / h^2, with h = 2.0\n")
cat("by default -- so the same h is 0.2 SD for a covariate with sd = 10 and\n")
cat("4 SD for one with sd = 0.5. Error in target (C) at theta_cov = 0.75:\n\n")
cat(sprintf("%10s %12s %14s %14s\n", "h", "h in SDs", "taylor2 err", "taylor4 err"))
tr <- target_C(th_true, obs)
for (h in c(0.05, 0.2, 0.5, 1.0, 2.0, 3.0))
  cat(sprintf("%10.2f %12.2f %14.5f %14.5f\n", h, h / sd_a,
              m_taylor(th_true, obs, h = h, ord = 2L) - tr,
              m_taylor(th_true, obs, h = h, ord = 4L) - tr))
cat("\n(sd_a = 1 here, so column 2 = column 1. For a covariate on a different\n")
cat(" scale the default h = 2.0 lands somewhere else entirely on this curve.)\n\n")

## ---------------------------------------------------------------------------
hr("6. Regime sweep -- where is each method adequate?")
cat("Ratio = theta_cov*sd_a / omega: how big the covariate effect is next to IIV.\n")
cat("Data regenerated at the true marginal for each ratio, then evaluated there.\n")
cat("'target gap' is |C - A| / |A|: how much the CHOICE of target matters.\n\n")
cat(sprintf("%7s %9s %12s %12s %12s %12s\n",
            "ratio", "theta_cov", "taylor2 rel", "gl-9 rel", "IS rel", "target gap"))
for (r in c(0.1, 0.25, 0.5, 1.0, 2.0, 3.0)) {
  tc <- r * om
  o2 <- make_obs(tc)
  tC <- target_C(tc, o2); tA <- target_A(tc, o2)
  t2 <- m_taylor(tc, o2, ord = 2L)
  gl <- m_gl(tc, o2)
  is <- m_is(tc, o2)$nll
  cat(sprintf("%7.2f %9.3f %12.2e %12.2e %12.2e %12.2e\n", r, tc,
              abs(t2 - tC) / abs(tC), abs(gl - tC) / abs(tC),
              abs(is - tC) / abs(tC), abs(tC - tA) / abs(tA)))
}
cat("\n")
