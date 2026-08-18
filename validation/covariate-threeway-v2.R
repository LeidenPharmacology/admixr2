## ============================================================================
## covariate-threeway-v2.R  -- CORRECTED re-run of covariate-threeway.R
## ============================================================================
##
## WHAT WAS WRONG IN THE ORIGINAL
## ------------------------------
## The original had three arms: `marginal`, `gh`, `taylor`.
##
##   `gh`     = sum_k w_k * nll2(obs_MARGINAL, pred_CONDITIONAL(a_k), n)
##   `taylor` = a second-order expansion of exactly that quantity
##
## Both score MARGINAL observations against CONDITIONAL predictions. That is a
## mismatch of TARGETS, not a quadrature method: the observed (E, V) already has
## the covariate integrated out, while every predicted block pretends the whole
## study sat at one covariate value. No amount of quadrature makes the two sides
## refer to the same population. Labelling it `gh` implied it was "Gauss-Hermite
## integration of the covariate", which it is not.
##
## The arm ACTUALLY named `marginal` in the original IS the correct Gauss-Hermite
## implementation -- it integrates the MOMENTS over the covariate with GH nodes
## (`mom_marg`, 21 nodes) and then evaluates the likelihood ONCE. It is renamed
## `gh_marginal` here so nothing reads it as a non-quadrature alternative.
##
## THE FOUR CONSTRUCTIONS
##   1 gh_marginal      E = E_a[g(a)], V = E_a[Vc(a)] + Cov_a(g(a)); ONE nll2.
##   2 stratified       per-node obs vs per-node pred (not testable here: this
##                      script's data are pooled; see covariate-matched-
##                      conditional-v2.R and covariate-node-retest-v2.R)
##   3 taylor_moments   2nd-order expansion of the MOMENTS, then ONE nll2
##                      (`mom_taylor`, validated in taylor-corrected.R)
##   4 mismatched       sum_k w_k nll2(obs_marginal, pred_cond(a_k)) -- the bug,
##                      retained ONLY as a display of the error, explicitly named
## ============================================================================

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
## CONSTRUCTION 1: marginal moments over (a, eta) -- the correct GH implementation
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
## CONSTRUCTION 3: 2nd-order expansion of the MOMENTS (reference impl from
## taylor-corrected.R, verbatim)
mom_taylor <- function(tcl, tcov, om, add, tt, pop, hfrac = 0.5) {
  mu <- pop[["mu"]]; s <- pop[["sd"]]; h <- hfrac * s
  m0 <- mom_cond(mu,     tcl, tcov, om, add, tt)
  mp <- mom_cond(mu + h, tcl, tcov, om, add, tt)
  mm <- mom_cond(mu - h, tcl, tcov, om, add, tt)
  dE  <- (mp$E - mm$E) / (2 * h)
  d2E <- (mp$E - 2 * m0$E + mm$E) / h^2
  d2V <- (mp$V - 2 * m0$V + mm$V) / h^2
  list(E = m0$E + 0.5 * s^2 * d2E,
       V = m0$V + 0.5 * s^2 * d2V + s^2 * outer(dE, dE))
}
nll2 <- function(obs, pred, n) {
  ch <- tryCatch(chol(pred$V), error = function(e) NULL)
  if (is.null(ch)) return(1e12)
  iv <- chol2inv(ch); r <- obs$E - pred$E
  n * (2*sum(log(diag(ch))) + sum(iv * obs$V) + as.numeric(t(r) %*% iv %*% r))
}

## ---- objectives -------------------------------------------------------------
obj_gh_marginal <- function(p, obs, n, tt)          # construction 1
  sum(vapply(seq_along(POPS), function(k)
    nll2(obs[[k]], mom_marg(p[1], p[2], exp(p[3]), exp(p[4]), tt, POPS[[k]]), n),
    numeric(1)))

obj_taylor_moments <- function(p, obs, n, tt, hfrac = 0.5)   # construction 3
  sum(vapply(seq_along(POPS), function(k)
    nll2(obs[[k]], mom_taylor(p[1], p[2], exp(p[3]), exp(p[4]), tt, POPS[[k]], hfrac), n),
    numeric(1)))

obj_mismatched <- function(p, obs, n, tt) {         # construction 4 -- THE BUG
  tot <- 0
  for (k in seq_along(POPS)) {
    pop <- POPS[[k]]
    tot <- tot + sum(QA$w * vapply(pop[["mu"]] + pop[["sd"]] * QA$x, function(a)
      nll2(obs[[k]], mom_cond(a, p[1], p[2], exp(p[3]), exp(p[4]), tt), n), numeric(1)))
  }
  tot
}
## the original `taylor` arm: a 2nd-order approx OF construction 4. Kept only so
## the report can say what the old number measured.
obj_mismatched_taylor <- function(p, obs, n, tt, hfrac = 0.5) {
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
cat("Arms: [1] gh_marginal = GH over the moments, one likelihood  (CORRECT)\n")
cat("      [3] taylor_moments = 2nd-order expansion of the moments (approx of [1])\n")
cat("      [4] mismatched = sum_k w_k NLL(marginal obs, conditional pred)  (THE BUG)\n")
cat("      [4t] mismatched_taylor = 2nd-order expansion of [4] (the old `taylor` arm)\n\n")
cat(sprintf("%-7s %6s %11s %11s %11s %11s | %10s %10s %10s %10s\n",
            "design", "n", "gh_marg[1]", "taylor[3]", "mismat[4]", "mism_tay[4t]",
            "bias[1]", "bias[3]", "bias[4]", "bias[4t]"))
for (dn in c("rich", "sparse")) {
  tt  <- if (dn == "rich") RICH else SPARSE
  obs <- lapply(POPS, function(pp) mom_marg(TCL, TCOV, OM, ADD, tt, pp))
  for (n in c(10L, 100L, 1000L)) {
    a1 <- argmin_tcov(obj_gh_marginal,       obs, n, tt)
    a3 <- argmin_tcov(obj_taylor_moments,    obs, n, tt)
    a4 <- argmin_tcov(obj_mismatched,        obs, n, tt)
    a5 <- argmin_tcov(obj_mismatched_taylor, obs, n, tt)
    flag <- paste0(ifelse(a1$boundary, "!", ""), ifelse(a3$boundary, "!", ""),
                   ifelse(a4$boundary, "!", ""), ifelse(a5$boundary, "!", ""))
    cat(sprintf("%-7s %6d %11.4f %11.4f %11.4f %11.4f | %+10.4f %+10.4f %+10.4f %+10.4f %s\n",
                dn, n, a1$at, a3$at, a4$at, a5$at,
                a1$at - TCOV, a3$at - TCOV, a4$at - TCOV, a5$at - TCOV, flag))
  }
}
cat("\n('!' marks an argmin at the search boundary -- not a converged answer)\n")

## Objective VALUES at truth, per subject. The meaningful gap is now
## [3] - [1] (the truncation error of the corrected expansion). [4] - [1] is a
## difference of TARGETS and is not a quadrature error at all.
cat("\nPer-subject objective at the true parameters (n divided out):\n")
cat(sprintf("%-7s %12s %12s %12s | %14s %14s\n", "design",
            "gh_marg[1]", "taylor[3]", "mismat[4]", "[3]-[1] trunc", "[4]-[1] target"))
for (dn in c("rich", "sparse")) {
  tt  <- if (dn == "rich") RICH else SPARSE
  obs <- lapply(POPS, function(pp) mom_marg(TCL, TCOV, OM, ADD, tt, pp))
  p0  <- c(TCL, TCOV, log(OM), log(ADD))
  for (n in c(10L, 1000L)) {
    m  <- obj_gh_marginal(p0, obs, n, tt)/n
    t3 <- obj_taylor_moments(p0, obs, n, tt)/n
    m4 <- obj_mismatched(p0, obs, n, tt)/n
    cat(sprintf("%-7s %12.5f %12.5f %12.5f | %14.5f %14.5f   (n=%d)\n",
                dn, m, t3, m4, t3 - m, m4 - m, n))
  }
}

## Moment accuracy of construction 3 against construction 1, per population.
cat("\nConstruction 3 vs construction 1, moments at the true parameters (rich):\n")
cat(sprintf("%-12s %8s %8s %14s %14s\n", "population", "mu_a", "sd_a", "rel err E", "rel err V"))
relerr <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-12))
for (k in seq_along(POPS)) {
  ex <- mom_marg(TCL, TCOV, OM, ADD, RICH, POPS[[k]])
  ta <- mom_taylor(TCL, TCOV, OM, ADD, RICH, POPS[[k]])
  cat(sprintf("%-12s %8.2f %8.2f %14.3e %14.3e\n", paste0("pop", k),
              POPS[[k]][["mu"]], POPS[[k]][["sd"]], relerr(ta$E, ex$E), relerr(ta$V, ex$V)))
}
