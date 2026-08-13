## Why did the development Rmd's GH/Laplace show no bias, while the same
## constructions are biased here?
##
## In that Rmd (lines 815 / 1082 / 1348) the OBSERVED data is regenerated at each
## node:  obsEV <- EV_given_wt(wt) . So obs and pred are BOTH conditional on the
## same wt. Four objectives, to separate that from everything else:
##
##   marginal    obs = marginal (E,V);  pred = marginal moments      <- my method
##   gh          obs = marginal (E,V);  pred = sum_k w_k NLL(a_k)    <- GH on REAL data
##   taylor      obs = marginal (E,V);  pred = 2nd-order approx of gh
##   gh_matched  obs = CONDITIONAL (E,V) at each a_k                 <- the Rmd's setup
##
## gh_matched needs one aggregate dataset PER node -- 17 of them for GH-17 -- which
## a published study does not provide. Run on one population and on three, because
## the single-population likelihood is exactly flat along
## tcl' = tcl + (b-b')*mu_a, omega'^2 = omega^2 + (b^2-b'^2)*sd_a^2.

TCL <- log(1); TV <- log(10); TCOV <- 0.75; OM <- 0.30; ADD <- 0.30; DOSE <- 100
ONE   <- list(c(mu = 0.00, sd = 0.55))
THREE <- list(c(mu = -0.45, sd = 0.30), c(mu = 0.00, sd = 0.55), c(mu = 0.50, sd = 0.35))
TT <- c(0.5, 1, 1.5, 2, 3, 4, 5, 6, 8)
NN <- 100L

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
nll2 <- function(obs, pred, n) {
  ch <- tryCatch(chol(pred$V), error = function(e) NULL)
  if (is.null(ch)) return(1e12)
  iv <- chol2inv(ch); r <- obs$E - pred$E
  n * (2*sum(log(diag(ch))) + sum(iv * obs$V) + as.numeric(t(r) %*% iv %*% r))
}

## obs sets. "marg" = one (E,V) per population. "cond" = one (E,V) per NODE.
obs_marg <- function(pops) lapply(pops, function(p) mom_marg(TCL, TCOV, OM, ADD, TT, p))
obs_cond <- function(pops) lapply(pops, function(p)
  lapply(p[["mu"]] + p[["sd"]] * QA$x, function(a) mom_cond(a, TCL, TCOV, OM, ADD, TT)))

obj_marginal <- function(p, obs, pops)
  sum(vapply(seq_along(pops), function(k)
    nll2(obs[[k]], mom_marg(p[1], p[2], exp(p[3]), exp(p[4]), TT, pops[[k]]), NN), numeric(1)))

obj_gh <- function(p, obs, pops)
  sum(vapply(seq_along(pops), function(k) {
    pop <- pops[[k]]
    sum(QA$w * vapply(pop[["mu"]] + pop[["sd"]] * QA$x, function(a)
      nll2(obs[[k]], mom_cond(a, p[1], p[2], exp(p[3]), exp(p[4]), TT), NN), numeric(1)))
  }, numeric(1)))

obj_taylor <- function(p, obs, pops, hfrac = 0.5)
  sum(vapply(seq_along(pops), function(k) {
    pop <- pops[[k]]; h <- hfrac * pop[["sd"]]
    f <- function(a) nll2(obs[[k]], mom_cond(a, p[1],p[2],exp(p[3]),exp(p[4]), TT), NN)
    f0 <- f(pop[["mu"]]); fp <- f(pop[["mu"]] + h); fm <- f(pop[["mu"]] - h)
    if (!all(is.finite(c(f0, fp, fm)))) return(1e12)
    f0 + 0.5 * pop[["sd"]]^2 * (fp - 2*f0 + fm) / h^2
  }, numeric(1)))

## the Rmd's objective: obs is the node's OWN conditional dataset
obj_gh_matched <- function(p, obs, pops)
  sum(vapply(seq_along(pops), function(k) {
    pop <- pops[[k]]; aa <- pop[["mu"]] + pop[["sd"]] * QA$x
    sum(QA$w * vapply(seq_along(aa), function(j)
      nll2(obs[[k]][[j]], mom_cond(aa[j], p[1],p[2],exp(p[3]),exp(p[4]), TT), NN),
      numeric(1)))
  }, numeric(1)))

LO <- 0.05; HI <- 2.20
argmin <- function(objfn, obs, pops) {
  prof <- function(tc) optim(c(TCL, log(OM), log(ADD)),
      function(v) objfn(c(v[1], tc, v[2], v[3]), obs, pops),
      method = "Nelder-Mead", control = list(maxit = 500, reltol = 1e-11))$value
  o <- optimize(prof, c(LO, HI), tol = 5e-4)
  list(at = o$minimum, bd = (o$minimum < LO + 0.02 || o$minimum > HI - 0.02))
}
## curvature: profile spread over a tcov window. ~0 means NOT identified.
spread <- function(objfn, obs, pops) {
  v <- vapply(c(0.45, 0.60, 0.75, 0.90, 1.05), function(tc)
    optim(c(TCL, log(OM), log(ADD)),
          function(vv) objfn(c(vv[1], tc, vv[2], vv[3]), obs, pops),
          method = "Nelder-Mead", control = list(maxit = 500, reltol = 1e-11))$value,
    numeric(1))
  v - min(v)
}

cat(sprintf("Truth tcov = %.4f, n = %d, exact data (no sampling noise)\n\n", TCOV, NN))
for (nm in c("ONE population", "THREE populations")) {
  pops <- if (nm == "ONE population") ONE else THREE
  om <- obs_marg(pops); oc <- obs_cond(pops)
  cat(sprintf("== %s ==\n", nm))
  cat(sprintf("%-34s %10s %10s   %s\n", "objective", "argmin", "bias",
              "profile spread over tcov 0.45..1.05"))
  for (r in list(
      list("marginal      (obs = marginal)", obj_marginal,   om),
      list("gh            (obs = marginal)", obj_gh,         om),
      list("taylor        (obs = marginal)", obj_taylor,     om),
      list("gh_matched  (obs = CONDITIONAL)", obj_gh_matched, oc))) {
    a <- argmin(r[[2]], r[[3]], pops); s <- spread(r[[2]], r[[3]], pops)
    cat(sprintf("%-34s %10.4f %+10.4f   %s%s\n", r[[1]], a$at, a$at - TCOV,
                paste(sprintf("%8.3f", s), collapse = " "), if (a$bd) "  [boundary]" else ""))
  }
  cat("\n")
}
