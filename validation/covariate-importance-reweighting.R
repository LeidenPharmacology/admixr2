## Pyry's importance-reweighting idea, pushed one step further.
##
## His proposal: draw eta once from N(0, Omega), then for covariate value a_i
## reweight by  w = phi(b; theta_cov*a_i, Omega) / phi(b; 0, Omega)  so the SAME
## simulations serve every covariate value.
##
## One step further: `a` is integrated over too, so average the weight over a
## before touching the ensemble --
##     E_a[w(a,b)] = p_u(b) / phi(b; 0, Omega),   u = Delta(a) + eta
## -- which is ONE weight per draw, the ratio of the u-density to the eta
## density. No quantile inversion, no Newton solve, and a discrete covariate is
## just a mixture density.
##
## Checked against exact nested Gauss-Hermite on an analytic 1-cmt model:
##   1. accuracy for normal / lognormal / DISCRETE covariates
##   2. several covariates sharing one eta
##   3. WHERE IT DEGENERATES -- effective sample size as Delta's spread grows
##      relative to omega, which is exactly the failure Pyry anticipated
kv <- function(k, v) cat(sprintf("%s\t%s\n", k, paste(v, collapse = ",")))
TCL <- log(1); TV <- log(10); DOSE <- 100; TT <- c(0.5, 1, 2, 4, 8)
conc <- function(cl) DOSE / exp(TV) * exp(outer(-cl / exp(TV), TT))

gh <- function(k) {
  i <- seq_len(k - 1L)
  J <- matrix(0, k, k); J[cbind(i, i + 1L)] <- sqrt(i); J[cbind(i + 1L, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE); o <- order(e$values)
  list(x = e$values[o], w = e$vectors[1L, o]^2)
}
Q <- gh(60L)

## --- EXACT: nested quadrature over (a, eta) ---------------------------------
exact <- function(dlt, om) {           # dlt: list(delta = values, w = probs)
  m1 <- 0; M2 <- 0
  for (i in seq_along(dlt$delta)) for (j in seq_along(Q$x)) {
    w <- dlt$w[i] * Q$w[j]
    y <- as.numeric(conc(exp(TCL + dlt$delta[i] + om * Q$x[j])))
    m1 <- m1 + w * y; M2 <- M2 + w * outer(y, y)
  }
  list(E = m1, V = M2 - outer(m1, m1))
}

## --- PYRY-IS: ONE ensemble at eta ~ N(0, om^2), reweighted by p_u/p_eta ------
pyry_is <- function(dlt, om, n = 4000L) {
  z   <- qnorm((seq_len(n) - 0.5) / n)          # deterministic N(0,1) grid
  eta <- om * z
  Y   <- conc(exp(TCL + eta))                   # the ONE set of simulations
  ## p_u(eta) = sum_i p_i * phi(eta; delta_i, om);  p_eta(eta) = phi(eta; 0, om)
  pu  <- rowSums(vapply(seq_along(dlt$delta),
          function(i) dlt$w[i] * dnorm(eta, dlt$delta[i], om), numeric(n)))
  w   <- pu / dnorm(eta, 0, om)
  w   <- w / sum(w)
  E   <- as.numeric(crossprod(w, Y))
  Yc  <- sweep(Y, 2L, E)
  list(E = E, V = t(Yc) %*% (Yc * w),
       ess = 1 / sum(w^2) / n)                  # effective sample size fraction
}

rel <- function(a, b) c(E = max(abs(a$E - b$E) / abs(b$E)),
                        V = max(abs(a$V - b$V)) / max(abs(b$V)))

## --- 1. accuracy across covariate distributions -----------------------------
mk_norm <- function(bcov, mu, sd, k = 60L) {
  g <- gh(k); list(delta = bcov * (mu + sd * g$x), w = g$w)
}
mk_lnorm <- function(bcov, ml, sl, k = 60L) {
  g <- gh(k); list(delta = bcov * exp(ml + sl * g$x), w = g$w)
}
mk_disc <- function(bcov, vals, p) list(delta = bcov * vals, w = p / sum(p))

cat(sprintf("%-34s %10s %10s %8s\n", "covariate distribution", "E rel", "V rel", "ESS"))
cases <- list(
  "normal  b=0.75 mu=0 sd=0.5"   = mk_norm(0.75, 0, 0.5),
  "normal  b=0.75 mu=0 sd=1.0"   = mk_norm(0.75, 0, 1.0),
  "lognormal b=0.75 (log 1, .3)" = mk_lnorm(0.75, 0, 0.3),
  "DISCRETE two-point 0/1 (.65)" = mk_disc(0.75, c(0, 1), c(.65, .35)),
  "DISCRETE three-point"         = mk_disc(0.75, c(0, 1, 2), c(.2, .5, .3)))
for (nm in names(cases)) {
  d <- cases[[nm]]; om <- 0.30
  r <- rel(pyry_is(d, om), exact(d, om))
  kv(sprintf("%-34s", nm),
     sprintf("%10.3e %10.3e %7.1f%%", r[["E"]], r[["V"]],
             100 * pyry_is(d, om)$ess))
}

## --- 2. several covariates sharing ONE eta ----------------------------------
## Delta is a SUM, so the mixture is a convolution -- no special handling.
g1 <- gh(30L); g2 <- gh(30L)
d2 <- list(delta = as.numeric(outer(0.75 * (0 + 0.4 * g1$x), 0.30 * (0 + 0.6 * g2$x), "+")),
           w     = as.numeric(outer(g1$w, g2$w)))
r <- rel(pyry_is(d2, 0.30), exact(d2, 0.30))
kv("two_covariates_one_eta",
   sprintf("E rel %.3e  V rel %.3e  ESS %.1f%%", r[["E"]], r[["V"]],
           100 * pyry_is(d2, 0.30)$ess))

## --- 3. WHERE IT DEGENERATES ------------------------------------------------
## Pyry: "at the extremes of the distributions the weighted mean and var-cov of
## the random effects will not be distributed ~N(0,Omega)". The controlling ratio
## is the spread of Delta against omega.
cat("\n")
cat(sprintf("%-8s %-8s %10s %10s %8s\n", "sd(Delta)", "omega", "E rel", "V rel", "ESS"))
for (sdd in c(0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5)) {
  om <- 0.30
  d  <- mk_norm(1.0, 0, sdd)
  p  <- pyry_is(d, om); r <- rel(p, exact(d, om))
  cat(sprintf("%-8.2f %-8.2f %10.3e %10.3e %7.1f%%\n",
              sdd, om, r[["E"]], r[["V"]], 100 * p$ess))
}

## --- 4. THE FIX: use the collapse's inflated covariance as the PROPOSAL ------
## The weights degenerate because N(0, omega^2) is a poor proposal for u, whose
## spread is omega^2 + var(Delta). admixr2 already computes that inflated
## covariance -- it is the closed-form collapse, Omega + J Sigma_a J'. Using it
## as the IS proposal rather than as an answer makes the weights near-1 by
## construction, and the residual weight corrects exactly the part the collapse
## gets wrong (non-normal a, non-linear effect).
pyry_is_inflated <- function(dlt, om, n = 4000L) {
  s_prop <- sqrt(om^2 + sum(dlt$w * (dlt$delta - sum(dlt$w * dlt$delta))^2))
  mu_prop <- sum(dlt$w * dlt$delta)
  z   <- qnorm((seq_len(n) - 0.5) / n)
  eta <- mu_prop + s_prop * z                  # draws from the INFLATED normal
  Y   <- conc(exp(TCL + eta))
  pu  <- rowSums(vapply(seq_along(dlt$delta),
          function(i) dlt$w[i] * dnorm(eta, dlt$delta[i], om), numeric(n)))
  w   <- pu / dnorm(eta, mu_prop, s_prop)
  w   <- w / sum(w)
  E   <- as.numeric(crossprod(w, Y)); Yc <- sweep(Y, 2L, E)
  list(E = E, V = t(Yc) %*% (Yc * w), ess = 1 / sum(w^2) / n)
}
cat("\n--- same cases, inflated proposal (Omega + J Sigma_a J') ---\n")
cat(sprintf("%-34s %10s %10s %8s\n", "covariate distribution", "E rel", "V rel", "ESS"))
for (nm in names(cases)) {
  d <- cases[[nm]]; om <- 0.30
  p <- pyry_is_inflated(d, om); r <- rel(p, exact(d, om))
  cat(sprintf("%-34s %10.3e %10.3e %7.1f%%\n", nm, r[["E"]], r[["V"]], 100 * p$ess))
}
cat("\n--- and across the degeneracy sweep that broke the plain proposal ---\n")
cat(sprintf("%-8s %10s %10s %8s\n", "sd(Delta)", "E rel", "V rel", "ESS"))
for (sdd in c(0.1, 0.3, 0.5, 0.75, 1.0, 1.5, 3.0)) {
  d <- mk_norm(1.0, 0, sdd)
  p <- pyry_is_inflated(d, 0.30); r <- rel(p, exact(d, 0.30))
  cat(sprintf("%-8.2f %10.3e %10.3e %7.1f%%\n", sdd, r[["E"]], r[["V"]], 100 * p$ess))
}
