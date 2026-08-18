## Is the node objective BIASED, or merely a different objective?
##
## The earlier three-way script tested exactly one configuration: node objective
## scored against ONE POOLED (E, V), with V_pred taken CONDITIONAL on the node.
## That confounds two separate things, so it cannot support a blanket claim.
## Here each is isolated. Data are exact marginal moments at the true parameters
## -- no sampling noise -- so any displaced argmin is a property of the OBJECTIVE,
## not a finite-sample effect.
##
##   A  gh-pooled       sum_j w_j NLL(pooled obs; mu(a_j), V_cond(a_j))   [original]
##   B  gh-pooled-Vm    as A but V_pred = the MARGINAL V   -> isolates the mean term
##   C  gh-pooled-Vc    as A but obs V replaced by V_cond  -> isolates the var term
##   D  gh-node-data    sum_j n_j NLL(obs_j; mu(a_j), V_cond(a_j))  -- node DATA
##   M  marginal        the shipped objective                          [reference]
##
## D is the configuration in which the node route is a genuine likelihood.

TCL <- log(1); TV <- log(10); TCOV <- 0.75; OM <- 0.30; ADD <- 0.30; DOSE <- 100
POPS <- list(c(mu = -0.45, sd = 0.30), c(mu = 0.00, sd = 0.55), c(mu = 0.50, sd = 0.35))
RICH <- c(0.5, 1, 1.5, 2, 3, 4, 5, 6, 8)

gh_nodes <- function(k) {
  i <- seq_len(k - 1L)
  J <- matrix(0, k, k); J[cbind(i, i + 1L)] <- sqrt(i); J[cbind(i + 1L, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE); o <- order(e$values)
  list(x = e$values[o], w = e$vectors[1L, o]^2)
}
QE <- gh_nodes(21L); QA <- gh_nodes(21L)
conc <- function(cl, tt) DOSE / exp(TV) * exp(outer(-cl / exp(TV), tt))

mom_cond <- function(a, tcl, tcov, om, add, tt) {
  Y  <- conc(exp(tcl + tcov * a + om * QE$x), tt)
  mu <- as.numeric(crossprod(QE$w, Y)); Yc <- sweep(Y, 2L, mu)
  V  <- t(Yc) %*% (Yc * QE$w); diag(V) <- diag(V) + add^2
  list(E = mu, V = V)
}
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
nll2 <- function(obsE, obsV, pred, n) {
  ch <- tryCatch(chol(pred$V), error = function(e) NULL)
  if (is.null(ch)) return(1e12)
  iv <- chol2inv(ch); r <- obsE - pred$E
  n * (2*sum(log(diag(ch))) + sum(iv * obsV) + as.numeric(t(r) %*% iv %*% r))
}
nodes_of <- function(pop) pop[["mu"]] + pop[["sd"]] * QA$x

## ---- objectives ------------------------------------------------------------
obj_marginal <- function(p, obs, n, tt)
  sum(vapply(seq_along(POPS), function(k) {
    pr <- mom_marg(p[1], p[2], exp(p[3]), exp(p[4]), tt, POPS[[k]])
    nll2(obs[[k]]$E, obs[[k]]$V, pr, n) }, numeric(1)))

obj_gh_pooled <- function(p, obs, n, tt) {          # A
  tot <- 0
  for (k in seq_along(POPS))
    tot <- tot + sum(QA$w * vapply(nodes_of(POPS[[k]]), function(a)
      nll2(obs[[k]]$E, obs[[k]]$V,
           mom_cond(a, p[1],p[2],exp(p[3]),exp(p[4]), tt), n), numeric(1)))
  tot
}
obj_gh_pooled_Vm <- function(p, obs, n, tt) {       # B: node MEAN, marginal V
  tot <- 0
  for (k in seq_along(POPS)) {
    Vm <- mom_marg(p[1], p[2], exp(p[3]), exp(p[4]), tt, POPS[[k]])$V
    tot <- tot + sum(QA$w * vapply(nodes_of(POPS[[k]]), function(a) {
      pr <- mom_cond(a, p[1],p[2],exp(p[3]),exp(p[4]), tt); pr$V <- Vm
      nll2(obs[[k]]$E, obs[[k]]$V, pr, n) }, numeric(1)))
  }
  tot
}
obj_gh_pooled_Vc <- function(p, obs, n, tt) {       # C: node V, obs V made conditional
  tot <- 0
  for (k in seq_along(POPS)) {
    Vc <- mom_cond(POPS[[k]][["mu"]], TCL, TCOV, OM, ADD, tt)$V   # truth-conditional
    tot <- tot + sum(QA$w * vapply(nodes_of(POPS[[k]]), function(a)
      nll2(obs[[k]]$E, Vc,
           mom_cond(a, p[1],p[2],exp(p[3]),exp(p[4]), tt), n), numeric(1)))
  }
  tot
}
obj_node_data <- function(p, obsN, n, tt) {         # D: node objective on node DATA
  tot <- 0
  for (k in seq_along(POPS)) {
    aa <- nodes_of(POPS[[k]])
    for (j in seq_along(aa))
      tot <- tot + nll2(obsN[[k]][[j]]$E, obsN[[k]][[j]]$V,
                        mom_cond(aa[j], p[1],p[2],exp(p[3]),exp(p[4]), tt),
                        n * QA$w[j])
  }
  tot
}

LO <- 0.02; HI <- 2.40
argmin_tcov <- function(objfn, obs, n, tt) {
  prof <- function(tc) optim(c(TCL, log(OM), log(ADD)),
      function(v) objfn(c(v[1], tc, v[2], v[3]), obs, n, tt),
      method = "Nelder-Mead", control = list(maxit = 600, reltol = 1e-11))$value
  o <- optimize(prof, c(LO, HI), tol = 5e-4)
  c(at = o$minimum, bnd = (o$minimum < LO + 0.02 || o$minimum > HI - 0.02))
}

tt   <- RICH
obs  <- lapply(POPS, function(pp) mom_marg(TCL, TCOV, OM, ADD, tt, pp))
obsN <- lapply(POPS, function(pp) lapply(nodes_of(pp), function(a)
          mom_cond(a, TCL, TCOV, OM, ADD, tt)))

cat(sprintf("Truth tcov = %.4f      exact (noiseless) data, 3 populations\n\n", TCOV))
cat(sprintf("%-34s %10s %10s %8s\n", "objective", "argmin", "error", "rel"))
row <- function(lab, fn, dat, n = 100L) {
  a <- argmin_tcov(fn, dat, n, tt)
  cat(sprintf("%-34s %10.4f %+10.4f %7.1f%%%s\n", lab, a[["at"]],
              a[["at"]] - TCOV, 100*(a[["at"]] - TCOV)/TCOV,
              if (a[["bnd"]] > 0) "  <- AT BOUNDARY" else ""))
  invisible(a[["at"]])
}
row("M  marginal (shipped)",             obj_marginal,     obs)
row("A  gh, pooled data, V_cond",        obj_gh_pooled,    obs)
row("B  gh, pooled data, V_marginal",    obj_gh_pooled_Vm, obs)
row("C  gh, conditional data + V_cond",  obj_gh_pooled_Vc, obs)
row("D  gh, NODE data + V_cond",         obj_node_data,    obsN)

## ---- the algebraic decomposition, evaluated numerically --------------------
## sum_j w_j r_j' Vi r_j  =  rbar' Vi rbar  +  tr(Vi * Cov_a(mu))
cat("\nMean-term decomposition at the true parameters (population 2, n=1):\n")
pop <- POPS[[2]]; aa <- nodes_of(pop)
MU  <- vapply(aa, function(a) mom_cond(a, TCL,TCOV,OM,ADD, tt)$E, numeric(length(tt)))
Vi  <- chol2inv(chol(obs[[2]]$V))
rj  <- obs[[2]]$E - MU
pool  <- sum(QA$w * colSums(rj * (Vi %*% rj)))
mubar <- as.numeric(MU %*% QA$w); rb <- obs[[2]]$E - mubar
Cova  <- (sweep(MU, 1L, mubar) * rep(QA$w, each = length(tt))) %*% t(sweep(MU, 1L, mubar))
cat(sprintf("  sum_j w_j r_j' Vi r_j        = %12.6f\n", pool))
cat(sprintf("  rbar' Vi rbar                = %12.6f\n", as.numeric(t(rb) %*% Vi %*% rb)))
cat(sprintf("  tr(Vi %%*%% Cov_a(mu))          = %12.6f\n", sum(Vi * Cova)))
cat(sprintf("  sum of the two               = %12.6f   (matches: %s)\n",
            as.numeric(t(rb) %*% Vi %*% rb) + sum(Vi * Cova),
            isTRUE(all.equal(pool, as.numeric(t(rb) %*% Vi %*% rb) + sum(Vi * Cova),
                             tolerance = 1e-8))))
