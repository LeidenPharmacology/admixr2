## Is there a FORMAL LIKELIHOOD for synthesising published models?
##
## The projection admixr2 runs is an M-estimator, not a likelihood: it weights
## each source by n_s, and n_s is a choice, so the Hessian reports whatever
## precision was asked for. The alternative is to treat the PUBLISHED ESTIMATES as
## the data:
##
##     theta_hat_s ~ N( g_s(phi), Sigma_s )
##
## phi      = the unified model's parameters
## g_s(phi) = the INDUCED map: what study s would have estimated, under its own
##            (smaller, misspecified) model and its own design and population, if
##            the truth were phi.  = argmin_psi E_{P_s} KL( p_phi || p_psi^(s) ),
##            i.e. the projection run in the OPPOSITE direction.
## Sigma_s  = the published uncertainty (RSEs), which is real information and not
##            a weighting choice -- that is the whole point.
##
## Setting: the omitted-covariate case, where g_s is available in closed form and
## was verified against the full numerical projection to 4 dp in
## validation/model-synthesis-omitted-covariate.R:
##
##     a_s = b_WT + b_CRCL * rho_s * sd_CRCL / sd_WT
##
## FOUR sources and TWO parameters, so the problem is OVER-identified -- which is
## exactly where weighting starts to matter and the two objectives separate. The
## n_s are deliberately NOT aligned with the SEs: the two small studies sampled
## richly and published the precise estimates.

B1 <- 0.75; B2 <- 0.40; S1 <- 0.30; S2 <- 0.35
RHO <- c(0.10, 0.35, 0.60, 0.80)
NS  <- c( 150,   40,  200,   30)     # patients   -- the projection's weight
SE  <- c(0.10, 0.03, 0.12, 0.04)     # published  -- the likelihood's weight
X   <- cbind(1, RHO * S2 / S1)       # d g_s / d phi   (linear here)
TRU <- c(b_WT = B1, b_CRCL = B2)
a_true <- as.numeric(X %*% TRU)

wls <- function(y, w) {              # returns estimate + the three variances
  W  <- diag(w); XtWX <- t(X) %*% W %*% X; A <- solve(XtWX)
  ph <- as.numeric(A %*% t(X) %*% W %*% y)
  S  <- diag(SE^2)
  list(phi = ph,
       ## the honest variance: propagate the SOURCES' uncertainty (delta method /
       ## sandwich). Valid for ANY weight w, including a chosen one.
       V_delta = A %*% t(X) %*% W %*% S %*% W %*% X %*% A,
       ## the naive variance: pretend the weighted objective is a log-likelihood
       ## and read its curvature, scaled by the residual spread.
       V_naive = {
         r <- y - X %*% ph
         s2 <- max(as.numeric(t(r) %*% W %*% r) / max(length(y) - 2L, 1L), 1e-12)
         s2 * A })
}

set.seed(20260815)
R <- 4000L
lab <- c("projection  (weights = n_s)", "likelihood  (weights = 1/SE^2)")
est <- array(0, c(R, 2, 2), dimnames = list(NULL, lab, names(TRU)))
cvD <- cvN <- array(FALSE, c(R, 2, 2), dimnames = dimnames(est))
for (r in seq_len(R)) {
  y <- a_true + stats::rnorm(length(a_true), 0, SE)      # the published numbers
  for (m in 1:2) {
    f <- wls(y, if (m == 1L) NS else 1 / SE^2)
    est[r, m, ] <- f$phi
    sdD <- sqrt(diag(f$V_delta)); sdN <- sqrt(diag(f$V_naive))
    cvD[r, m, ] <- abs(f$phi - TRU) < 1.96 * sdD
    cvN[r, m, ] <- abs(f$phi - TRU) < 1.96 * sdN
  }
}

cat(sprintf("truth: b_WT %.3f   b_CRCL %.3f      %d replicates, 4 sources, 2 parameters\n\n",
            B1, B2, R))
cat(sprintf("%-31s %8s %8s %8s %10s %10s\n", "", "param", "bias", "SD",
            "cover:delta", "cover:naive"))
for (m in 1:2) for (p in 1:2)
  cat(sprintf("%-31s %8s %+8.4f %8.4f %10.3f %10.3f\n",
              if (p == 1L) lab[m] else "", names(TRU)[p],
              mean(est[, m, p]) - TRU[p], stats::sd(est[, m, p]),
              mean(cvD[, m, p]), mean(cvN[, m, p])))

cat(sprintf("\nefficiency  SD(projection) / SD(likelihood):  b_WT %.2fx   b_CRCL %.2fx\n",
            stats::sd(est[,1,1])/stats::sd(est[,2,1]),
            stats::sd(est[,1,2])/stats::sd(est[,2,2])))

## Over-identification also buys a GOODNESS-OF-FIT TEST, which the projection has
## no analogue of: with 4 sources and 2 parameters there are 2 df of disagreement
## that a single unified structure has to explain.
q <- vapply(seq_len(200L), function(i) {
  y <- a_true + stats::rnorm(4L, 0, SE)
  f <- wls(y, 1 / SE^2); r <- y - X %*% f$phi
  as.numeric(t(r) %*% diag(1/SE^2) %*% r)
}, 0)
cat(sprintf("\nQ statistic under a TRUE common structure: mean %.2f (df = 2), P(Q > chisq_.95) = %.3f\n",
            mean(q), mean(q > stats::qchisq(0.95, 2))))
## ...and under a source that does NOT share the structure (source 3 shifted):
q2 <- vapply(seq_len(200L), function(i) {
  y <- a_true + stats::rnorm(4L, 0, SE); y[3] <- y[3] + 0.35
  f <- wls(y, 1 / SE^2); r <- y - X %*% f$phi
  as.numeric(t(r) %*% diag(1/SE^2) %*% r)
}, 0)
cat(sprintf("Q with one source off-structure:            mean %.2f,          rejected %.3f\n",
            mean(q2), mean(q2 > stats::qchisq(0.95, 2))))
