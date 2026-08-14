## Does the inflated proposal COMPLETELY fix the weight degeneracy?
##
## No. It fixes a VARIANCE mismatch, which is the dominant part, and it fixes it
## exactly: p_u has spread Omega + var(Delta), and that is what the proposal is
## given. What it does not fix is a SHAPE mismatch. p_u is a mixture,
##
##     p_u(x) = sum_k pw_k * phi(x; Delta_k, omega),
##
## and a variance-matched NORMAL is only a good cover for it when that mixture is
## roughly normal. Two situations where it is not:
##
##   1. a DISCRETE covariate with small omega -- the mixture is then a set of
##      well-separated narrow spikes, and a single wide normal covers the gaps
##      between them, where p_u is nearly zero, while under-covering the spikes;
##   2. a heavy-tailed or strongly skewed covariate distribution.
##
## The controlling quantity for (1) is the SEPARATION of the mixture components
## against omega, which is a different ratio from the sd(Delta)/omega that drove
## the plain-proposal failure. So this sweeps it.
##
## Run:  Rscript validation/covariate-reweighting-proposal-limits.R
TCL <- 0; TV <- log(10); DOSE <- 100; TT <- c(0.5, 1, 2, 4, 8)
conc <- function(cl) DOSE / exp(TV) * exp(outer(-cl / exp(TV), TT))
gh <- function(k) {
  i <- seq_len(k - 1L)
  J <- matrix(0, k, k); J[cbind(i, i + 1L)] <- sqrt(i); J[cbind(i + 1L, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE); o <- order(e$values)
  list(x = e$values[o], w = e$vectors[1L, o]^2)
}
QE <- gh(60L)
exact <- function(dl, om) {
  m1 <- 0; M2 <- 0
  for (i in seq_along(dl$delta)) for (j in seq_along(QE$x)) {
    w <- dl$w[i] * QE$w[j]
    y <- as.numeric(conc(exp(TCL + dl$delta[i] + om * QE$x[j])))
    m1 <- m1 + w * y; M2 <- M2 + w * outer(y, y)
  }
  list(E = m1, V = M2 - outer(m1, m1))
}
## proposal: "inflated" = variance-matched normal;  "mixture" = p_u itself
rw <- function(dl, om, kind = "inflated", n = 4000L) {
  mu_d <- sum(dl$w * dl$delta); sd_d <- sqrt(sum(dl$w * (dl$delta - mu_d)^2))
  z <- qnorm((seq_len(n) - 0.5) / n)
  if (kind == "inflated") {
    s_p <- sqrt(om^2 + sd_d^2); u <- mu_d + s_p * z
    pp  <- dnorm(u, mu_d, s_p)
  } else {
    ## stratified draw from the mixture itself: allocate the n draws across the
    ## components in proportion to their weight, each stratum a normal about its
    ## own Delta_k. This IS p_u, so the weights are identically 1.
    cnt <- pmax(1L, round(dl$w / sum(dl$w) * n)); u <- numeric(0)
    for (k in seq_along(dl$delta)) {
      zz <- qnorm((seq_len(cnt[k]) - 0.5) / cnt[k])
      u  <- c(u, dl$delta[k] + om * zz)
    }
    n  <- length(u)
    pp <- rowSums(vapply(seq_along(dl$delta),
            function(k) dl$w[k] * dnorm(u, dl$delta[k], om), numeric(n)))
  }
  pu <- rowSums(vapply(seq_along(dl$delta),
          function(k) dl$w[k] * dnorm(u, dl$delta[k], om), numeric(length(u))))
  w  <- pu / pp; w[!is.finite(w)] <- 0
  if (sum(w) <= 0) return(NULL)
  wn <- w / sum(w)
  Y  <- conc(exp(TCL + u))
  E  <- as.numeric(crossprod(wn, Y)); Yc <- sweep(Y, 2L, E)
  list(E = E, V = t(Yc) %*% (Yc * wn), ess = (sum(w)^2 / sum(w^2)) / length(w))
}
rel <- function(a, b) c(max(abs(a$E - b$E) / abs(b$E)),
                        max(abs(a$V - b$V)) / max(abs(b$V)))

cat("1. DISCRETE covariate: separation of the mixture components against omega.\n")
cat("   Delta in {0, 0.75}; omega shrinks, so the two spikes separate.\n\n")
cat(sprintf("%-8s %-10s %10s %10s %8s   %10s %8s\n",
            "omega", "sep/omega", "E rel", "V rel", "ESS",
            "V rel", "ESS"))
cat(sprintf("%-8s %-10s %10s %10s %8s   %10s %8s\n",
            "", "", "(inflated)", "(inflated)", "(infl)",
            "(mixture)", "(mix)"))
for (om in c(0.60, 0.30, 0.20, 0.12, 0.07, 0.04)) {
  dl <- list(delta = c(0, 0.75), w = c(0.65, 0.35))
  x  <- exact(dl, om)
  a  <- rw(dl, om, "inflated"); b <- rw(dl, om, "mixture")
  ra <- rel(a, x); rb <- rel(b, x)
  cat(sprintf("%-8.2f %-10.1f %10.2e %10.2e %7.1f%%   %10.2e %7.1f%%\n",
              om, 0.75 / om, ra[1], ra[2], 100 * a$ess, rb[2], 100 * b$ess))
}

cat("\n2. Skewed continuous covariate: a lognormal Delta with a growing effect.\n\n")
cat(sprintf("%-10s %10s %10s %8s   %10s %8s\n",
            "coef", "E rel", "V rel", "ESS", "V rel", "ESS"))
QA <- gh(60L)
for (b in c(0.25, 0.75, 1.5, 3.0)) {
  dl <- list(delta = b * exp(0 + 0.5 * QA$x), w = QA$w)
  om <- 0.30
  x  <- exact(dl, om); a <- rw(dl, om, "inflated"); m <- rw(dl, om, "mixture")
  ra <- rel(a, x); rm_ <- rel(m, x)
  cat(sprintf("%-10.2f %10.2e %10.2e %7.1f%%   %10.2e %7.1f%%\n",
              b, ra[1], ra[2], 100 * a$ess, rm_[2], 100 * m$ess))
}
cat("\nThe mixture proposal is p_u itself, so its weights are identically 1 and\n")
cat("its ESS is 100% by construction -- it is what the inflated normal is an\n")
cat("approximation to, and the gap between the two columns is the shape mismatch.\n")
