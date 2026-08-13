## ============================================================================
## The u-distribution framing, tested across functional forms
## ============================================================================
##
##   Rscript validation/covariate-u-distribution.R
##
## Self-contained: no admixr2, no rxode2, no packages.
##
## Premise under test: the model sees the covariate and the random effect only
## through their SUM in the mu-referenced argument,
##
##     u = Delta(a) + eta,      eta ~ N(0, Omega)
##
## for ANY covariate effect Delta -- linear, allometric, exponential, Emax,
## categorical. Delta(a) is arithmetic, so it is free to evaluate; only f(u)
## costs an ODE solve. So the whole solve budget should be spent representing
## the distribution of u as well as possible.
##
## Methods compared at an EQUAL budget of N solved subjects:
##   u-quantile   build u's distribution cheaply, then solve at N quantiles of it
##   per-row      draw N (a, eta) pairs jointly (Sobol) and solve at those
##   IS-marginal  Pyry: pool N draws from N(0, Omega), reweight to the marginal
##   collapse     exact closed form -- only defined when Delta is linear
##   mean-only    ignore the covariate spread (the baseline failure mode)

## ---- model: only f(u) costs a solve ----------------------------------------
TCL <- 1.0; VD <- 10; DOSE <- 100; TIMES <- c(0.5, 1, 2, 4, 8)
fpred <- function(u) DOSE / VD * exp(outer(-exp(TCL + u) / VD, TIMES))
OM <- 0.30

gh <- function(k) {
  i <- seq_len(k - 1L)
  J <- matrix(0, k, k); J[cbind(i, i + 1L)] <- sqrt(i); J[cbind(i + 1L, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE); o <- order(e$values)
  list(x = e$values[o], w = e$vectors[1L, o]^2)
}

## ---- the covariate models --------------------------------------------------
## Each carries Delta(a), the covariate's distribution as (nodes, weights) for
## an exact reference, and a quantile function for sampling-based methods.
qnorm_ <- function(mu, sd) function(p) qnorm(p, mu, sd)
FORMS <- list(
  linear = list(
    lab = "linear      0.75*a,        a ~ N(0, 0.6)",
    D = function(a) 0.75 * a,
    nodes = function(k) { g <- gh(k); list(x = 0 + 0.6 * g$x, w = g$w) },
    q = qnorm_(0, 0.6), linear = TRUE, mu_a = 0, sd_a = 0.6, slope = 0.75),
  allometric = list(
    lab = "allometric  0.75*log(a/70), a ~ lnorm(log70, 0.20)",
    D = function(a) 0.75 * log(a / 70),
    nodes = function(k) { g <- gh(k); list(x = exp(log(70) + 0.20 * g$x), w = g$w) },
    q = function(p) qlnorm(p, log(70), 0.20), linear = FALSE),
  exponential = list(
    lab = "exponential 0.02*(a-70),   a ~ N(70, 15)",
    D = function(a) 0.02 * (a - 70),
    nodes = function(k) { g <- gh(k); list(x = 70 + 15 * g$x, w = g$w) },
    q = qnorm_(70, 15), linear = FALSE),
  emax = list(
    # NOT a normal covariate: a normal weight puts quadrature mass at a < 0,
    # where 1 + 0.8a/(50+a) turns negative and log() is NaN. The assumed
    # covariate distribution has to stay inside the model's domain -- a real
    # constraint, and an argument for validating the node range against the
    # model rather than trusting mu/sd.
    lab = "Emax        log(1+0.8a/(50+a)), a ~ lnorm(log70, 0.20)",
    D = function(a) log(1 + 0.8 * a / (50 + a)),
    nodes = function(k) { g <- gh(k); list(x = exp(log(70) + 0.20 * g$x), w = g$w) },
    q = function(p) qlnorm(p, log(70), 0.20), linear = FALSE),
  categorical = list(
    lab = "categorical 0.35*a,        a in {0,1}, p = (0.6, 0.4)",
    D = function(a) 0.35 * a,
    nodes = function(k) list(x = c(0, 1), w = c(0.6, 0.4)),
    q = function(p) as.numeric(p > 0.6), linear = FALSE, disc = TRUE)
)

## ---- weighted moments of f over a set of u values --------------------------
mom_u <- function(u, w = rep(1 / length(u), length(u))) {
  F <- fpred(u); w <- w / sum(w)
  mu <- as.numeric(crossprod(w, F)); Fc <- sweep(F, 2L, mu)
  list(mu = mu, V = t(Fc) %*% (Fc * w))
}

## ---- TRUTH: dense nested quadrature over (a, eta) --------------------------
truth <- function(fm, ka = 60L, kb = 60L) {
  na <- fm$nodes(ka); ge <- gh(kb)
  u <- as.numeric(outer(fm$D(na$x), OM * ge$x, "+"))
  w <- as.numeric(outer(na$w, ge$w))
  mom_u(u, w)
}

## ---- METHOD 1: u-quantile (the proposal) -----------------------------------
## Build u's distribution from a FINE grid -- pure arithmetic, no solves -- then
## spend the budget on N quantiles of it.
## u = Delta(a) + eta is a convolution, and its CDF is available in CLOSED FORM
## given a quadrature over a:
##
##     F_u(t) = E_a[ Phi( (t - Delta(a)) / omega ) ] = sum_j w_j Phi(...)
##
## No sampling, no lattice artefacts, and Phi carries the tails analytically --
## which is what a grid of summed quantiles cannot do (an earlier version of this
## built u by summing 400 covariate quantiles with 100 eta quantiles, which
## silently truncated u at +-2.6 sigma and capped the accuracy at ~1e-2).
##
## Gauss-Hermite is the right rule for the a-integral here, because here we ARE
## integrating rather than representing a distribution.
m_uquant <- function(fm, N, ka = 120L, ngrid = 40001L) {
  na <- if (isTRUE(fm$disc)) fm$nodes(0L) else fm$nodes(ka)
  d  <- fm$D(na$x); w <- na$w / sum(na$w)
  t  <- seq(min(d) - 9 * OM, max(d) + 9 * OM, length.out = ngrid)
  Fu <- as.numeric(pnorm(outer(t, d, "-") / OM) %*% w)
  mom_u(approx(Fu, t, xout = (seq_len(N) - 0.5) / N, rule = 2)$y)
}

## ---- METHOD 2: per-row joint sampling --------------------------------------
## Sobol in 2 dims; dim 1 -> eta, dim 2 -> covariate (distinct dimensions).
m_perrow <- function(fm, N) {
  s <- .sobol2(N)
  mom_u(fm$D(fm$q(s[, 2L])) + OM * qnorm(s[, 1L]))
}

## ---- METHOD 3: Pyry's reweighting, targeted at the MARGINAL ----------------
## Pool from N(0, Omega); the marginal target is the mixture over a, so
##   w(b) = [ sum_j p_j N(b; Delta(a_j), Omega) ] / N(b; 0, Omega).
m_is <- function(fm, N, ka = 60L) {
  b <- OM * qnorm((seq_len(N) - 0.5) / N)
  na <- fm$nodes(ka); d <- fm$D(na$x)
  num <- rowSums(vapply(seq_along(d), function(j) na$w[j] * dnorm(b, d[j], OM),
                        numeric(N)))
  w <- num / dnorm(b, 0, OM)
  list(m = mom_u(b, w), ess = (sum(w)^2 / sum(w^2)) / N)
}

## ---- METHOD 4: collapse (linear only) --------------------------------------
m_collapse <- function(fm, N) {
  if (!isTRUE(fm$linear)) return(NULL)
  sd_u <- sqrt(fm$slope^2 * fm$sd_a^2 + OM^2)
  mom_u(fm$slope * fm$mu_a + sd_u * qnorm((seq_len(N) - 0.5) / N))
}

## ---- METHOD 5: mean-only ---------------------------------------------------
m_mean <- function(fm, N, ka = 60L) {
  na <- fm$nodes(ka)
  mom_u(fm$D(sum(na$w * na$x)) + OM * qnorm((seq_len(N) - 0.5) / N))
}

## van der Corput / Sobol-lite: 2-D stratified low-discrepancy without packages
.sobol2 <- function(n) {
  vdc <- function(n, base) vapply(seq_len(n), function(i) {
    f <- 1; r <- 0
    while (i > 0) { f <- f / base; r <- r + f * (i %% base); i <- i %/% base }
    r
  }, numeric(1))
  cbind(pmin(pmax(vdc(n, 2L), 1e-9), 1 - 1e-9),
        pmin(pmax(vdc(n, 3L), 1e-9), 1 - 1e-9))
}

relerr <- function(a, b) max(abs(a$mu / b$mu - 1), abs(a$V / b$V - 1))

cat("\nRelative error in the aggregate moments (max over E and V),\n")
cat("at an equal budget of N solved subjects. Truth = dense nested quadrature.\n")
for (N in c(500L, 2000L)) {
  cat("\n", strrep("=", 92), "\n  N = ", N, " solved subjects\n", strrep("=", 92), "\n", sep = "")
  cat(sprintf("%-46s %10s %10s %10s %10s %10s\n", "covariate model",
              "u-quantile", "per-row", "IS(marg)", "collapse", "mean-only"))
  for (nm in names(FORMS)) {
    fm <- FORMS[[nm]]; tr <- truth(fm)
    isr <- m_is(fm, N); co <- m_collapse(fm, N)
    cat(sprintf("%-46s %10.2e %10.2e %10.2e %10s %10.2e\n", fm$lab,
                relerr(m_uquant(fm, N), tr),
                relerr(m_perrow(fm, N), tr),
                relerr(isr$m, tr),
                if (is.null(co)) "n/a" else sprintf("%.2e", relerr(co, tr)),
                relerr(m_mean(fm, N), tr)))
  }
}

cat("\n", strrep("=", 92), "\n  Effective sample size of the IS pool (N = 2000)\n",
    strrep("=", 92), "\n", sep = "")
for (nm in names(FORMS))
  cat(sprintf("%-46s ESS/N = %.4f   (%6.1f effective of 2000)\n",
              FORMS[[nm]]$lab, m_is(FORMS[[nm]], 2000L)$ess,
              2000 * m_is(FORMS[[nm]], 2000L)$ess))
cat("\n")
