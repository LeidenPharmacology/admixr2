## THREE-WAY comparison, to separate two different things that have been
## conflated all along.
##
##   marginal :  -2LL at the covariate-marginal moments        (my construction)
##   gh       :  sum_j w_j * NLL(a_j)                          (exact E_a[NLL])
##   taylor   :  NLL(mu_a) + 0.5 sd^2 * d2NLL/da2              (2nd-order approx of gh)
##
## gh is the EXACT quadrature of what taylor approximates, so
##   gh - taylor    isolates taylor's TRUNCATION error
##   marginal - gh  isolates the difference of TARGET
## Conflating those two is what made the earlier comparisons unreadable.
##
## Data are the EXACT marginal moments at the true parameters -- no sampling
## noise -- so a displaced optimum is a property of the objective. THREE
## populations with differing covariate means and spreads, because one population
## does not identify the coefficient at all.

TCL <- log(1); TV <- log(10); TCOV <- 0.75; OM <- 0.30; ADD <- 0.30; DOSE <- 100
POPS <- list(c(mu = -0.45, sd = 0.30), c(mu = 0.00, sd = 0.55), c(mu = 0.50, sd = 0.35))
RICH   <- c(0.5, 1, 1.5, 2, 3, 4, 5, 6, 8)
SPARSE <- c(6, 8)

gh_nodes <- function(k) {
  i <- seq_len(k - 1L)
  J <- matrix(0, k, k); J[cbind(i, i + 1L)] <- sqrt(i); J[cbind(i + 1L, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE); o <- order(e$values)
  list(x = e$values[o], w = e$vectors[1L, o]^2)
}
QE <- gh_nodes(21L)   # eta integration
QA <- gh_nodes(21L)   # covariate integration
conc <- function(cl, tt) DOSE / exp(TV) * exp(outer(-cl / exp(TV), tt))

## moments conditional on covariate value a
mom_cond <- function(a, tcl, tcov, om, add, tt) {
  Y  <- conc(exp(tcl + tcov * a + om * QE$x), tt)
  mu <- as.numeric(crossprod(QE$w, Y)); Yc <- sweep(Y, 2L, mu)
  V  <- t(Yc) %*% (Yc * QE$w); diag(V) <- diag(V) + add^2
  list(E = mu, V = V)
}
## marginal moments over (a, eta) for one population
mom_marg <- function(tcl, tcov, om, add, tt, pop) {
  m1 <- 0; M2 <- 0
  for (ia in seq_along(QA$x)) {
    a <- pop[["mu"]] + pop[["sd"]] * QA$x[ia]
    Y <- conc(exp(tcl + tcov * a + om * QE$x), tt)
    for (ib in seq_along(QE$x)) {
      w <- QA$w[ia] * QE$w[ib]; y <- Y[ib, ]
      m1 <- m1 + w * y; M2 <- M2 + w * outer(y, y)
    }
  }
  V <- M2 - outer(m1, m1); diag(V) <- diag(V) + add^2
  list(E = m1, V = V)
}
nll2 <- function(obs, pred, n) {
  ch <- tryCatch(chol(pred$V), error = function(e) NULL)
  if (is.null(ch)) return(1e12)
  iv <- chol2inv(ch); r <- obs$E - pred$E
  n * (2*sum(log(diag(ch))) + sum(iv * obs$V) + as.numeric(t(r) %*% iv %*% r))
}

obj_marginal <- function(p, obs, n, tt)
  sum(vapply(seq_along(POPS), function(k)
    nll2(obs[[k]], mom_marg(p[1], p[2], exp(p[3]), exp(p[4]), tt, POPS[[k]]), n),
    numeric(1)))

obj_gh <- function(p, obs, n, tt) {
  tot <- 0
  for (k in seq_along(POPS)) {
    pop <- POPS[[k]]
    tot <- tot + sum(QA$w * vapply(pop[["mu"]] + pop[["sd"]] * QA$x, function(a)
      nll2(obs[[k]], mom_cond(a, p[1], p[2], exp(p[3]), exp(p[4]), tt), n), numeric(1)))
  }
  tot
}

obj_taylor <- function(p, obs, n, tt, hfrac = 0.5) {
  tot <- 0
  for (k in seq_along(POPS)) {
    pop <- POPS[[k]]; h <- hfrac * pop[["sd"]]
    f0 <- nll2(obs[[k]], mom_cond(pop[["mu"]],     p[1],p[2],exp(p[3]),exp(p[4]), tt), n)
    fp <- nll2(obs[[k]], mom_cond(pop[["mu"]] + h, p[1],p[2],exp(p[3]),exp(p[4]), tt), n)
    fm <- nll2(obs[[k]], mom_cond(pop[["mu"]] - h, p[1],p[2],exp(p[3]),exp(p[4]), tt), n)
    if (!all(is.finite(c(f0, fp, fm)))) return(1e12)
    tot <- tot + f0 + 0.5 * pop[["sd"]]^2 * (fp - 2*f0 + fm) / h^2
  }
  tot
}

LO <- 0.05; HI <- 2.20
argmin_tcov <- function(objfn, obs, n, tt) {
  prof <- function(tc) optim(c(TCL, log(OM), log(ADD)),
      function(v) objfn(c(v[1], tc, v[2], v[3]), obs, n, tt),
      method = "Nelder-Mead", control = list(maxit = 400, reltol = 1e-11))$value
  o <- optimize(prof, c(LO, HI), tol = 5e-4)
  list(at = o$minimum, boundary = (o$minimum < LO + 0.02 || o$minimum > HI - 0.02))
}

cat(sprintf("Truth: tcov = %.4f   (exact marginal data, 3 identifying populations)\n\n", TCOV))
cat(sprintf("%-7s %6s %10s %10s %10s | %10s %10s %10s\n",
            "design", "n", "marginal", "gh", "taylor",
            "marg bias", "gh bias", "tayl bias"))
for (dn in c("rich", "sparse")) {
  tt  <- if (dn == "rich") RICH else SPARSE
  obs <- lapply(POPS, function(pp) mom_marg(TCL, TCOV, OM, ADD, tt, pp))
  for (n in c(10L, 100L, 1000L)) {
    a1 <- argmin_tcov(obj_marginal, obs, n, tt)
    a2 <- argmin_tcov(obj_gh,       obs, n, tt)
    a3 <- argmin_tcov(obj_taylor,   obs, n, tt)
    flag <- paste0(ifelse(a1$boundary, "!", ""), ifelse(a2$boundary, "!", ""),
                   ifelse(a3$boundary, "!", ""))
    cat(sprintf("%-7s %6d %10.4f %10.4f %10.4f | %+10.4f %+10.4f %+10.4f %s\n",
                dn, n, a1$at, a2$at, a3$at,
                a1$at - TCOV, a2$at - TCOV, a3$at - TCOV, flag))
  }
}
cat("\n('!' marks an argmin at the search boundary -- not a converged answer)\n")

## Objective VALUES at truth: the two gaps, separated, per subject
cat("\nPer-subject objective at the true parameters (n divided out):\n")
cat(sprintf("%-7s %12s %12s %12s | %12s %12s\n", "design",
            "marginal/n", "gh/n", "taylor/n", "gh-marginal", "taylor-gh"))
for (dn in c("rich", "sparse")) {
  tt  <- if (dn == "rich") RICH else SPARSE
  obs <- lapply(POPS, function(pp) mom_marg(TCL, TCOV, OM, ADD, tt, pp))
  p0  <- c(TCL, TCOV, log(OM), log(ADD))
  for (n in c(10L, 1000L)) {
    m <- obj_marginal(p0, obs, n, tt)/n
    g <- obj_gh(p0, obs, n, tt)/n
    t_ <- obj_taylor(p0, obs, n, tt)/n
    cat(sprintf("%-7s %12.5f %12.5f %12.5f | %12.5f %12.5f   (n=%d)\n",
                dn, m, g, t_, g - m, t_ - g, n))
  }
}
