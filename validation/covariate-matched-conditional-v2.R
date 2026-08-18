## ============================================================================
## covariate-matched-conditional-v2.R -- CORRECTED re-run
## ============================================================================
##
## WHAT WAS WRONG IN THE ORIGINAL
## ------------------------------
## Original arms: marginal / gh / taylor / gh_matched.
##
##   `gh`     = sum_k w_k nll2(obs_MARGINAL, pred_CONDITIONAL(a_k))   <- BUG
##   `taylor` = a 2nd-order expansion of that same quantity           <- BUG
##
## Both score marginal observations against conditional predictions: a mismatch
## of targets, not a quadrature method. Their argmins (0.0502 at the boundary in
## one population, 0.3045 in three) measure that category error and nothing else.
##
## `marginal` was already the CORRECT Gauss-Hermite implementation (GH over the
## moments, then one likelihood) and is renamed `gh_marginal`.
## `gh_matched` was already correct: per-node observed blocks against per-node
## predicted blocks -- the STRATIFIED construction -- and is renamed
## `stratified` to say so. It is legitimate exactly when the source's data
## genuinely resolve the covariate per node.
##
## ARMS HERE
##   [1] gh_marginal   GH over moments, one likelihood            (correct)
##   [2] stratified    per-node obs vs per-node pred              (correct, needs
##                     one aggregate dataset PER node)
##   [3] taylor_moments 2nd-order expansion of the MOMENTS        (approx of [1])
##   [4] mismatched    marginal obs vs conditional pred           (the bug, shown
##                     ONLY as the error, explicitly labelled)
##
## Run on one population and on three, because the single-population likelihood
## is exactly flat along tcl' = tcl + (b-b')*mu_a,
## omega'^2 = omega^2 + (b^2-b'^2)*sd_a^2.
## ============================================================================

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

## obs sets. "marg" = one (E,V) per population. "cond" = one (E,V) per NODE.
obs_marg <- function(pops) lapply(pops, function(p) mom_marg(TCL, TCOV, OM, ADD, TT, p))
obs_cond <- function(pops) lapply(pops, function(p)
  lapply(p[["mu"]] + p[["sd"]] * QA$x, function(a) mom_cond(a, TCL, TCOV, OM, ADD, TT)))

## [1] correct GH: integrate the moments, evaluate the likelihood once
obj_gh_marginal <- function(p, obs, pops)
  sum(vapply(seq_along(pops), function(k)
    nll2(obs[[k]], mom_marg(p[1], p[2], exp(p[3]), exp(p[4]), TT, pops[[k]]), NN), numeric(1)))

## [3] corrected Taylor: expand the MOMENTS, evaluate the likelihood once
obj_taylor_moments <- function(p, obs, pops, hfrac = 0.5)
  sum(vapply(seq_along(pops), function(k)
    nll2(obs[[k]], mom_taylor(p[1],p[2],exp(p[3]),exp(p[4]), TT, pops[[k]], hfrac), NN),
    numeric(1)))

## [4] THE BUG: marginal obs scored against conditional predictions
obj_mismatched <- function(p, obs, pops)
  sum(vapply(seq_along(pops), function(k) {
    pop <- pops[[k]]
    sum(QA$w * vapply(pop[["mu"]] + pop[["sd"]] * QA$x, function(a)
      nll2(obs[[k]], mom_cond(a, p[1], p[2], exp(p[3]), exp(p[4]), TT), NN), numeric(1)))
  }, numeric(1)))

## [4t] the old `taylor` arm: a 2nd-order expansion OF the bug
obj_mismatched_taylor <- function(p, obs, pops, hfrac = 0.5)
  sum(vapply(seq_along(pops), function(k) {
    pop <- pops[[k]]; h <- hfrac * pop[["sd"]]
    f <- function(a) nll2(obs[[k]], mom_cond(a, p[1],p[2],exp(p[3]),exp(p[4]), TT), NN)
    f0 <- f(pop[["mu"]]); fp <- f(pop[["mu"]] + h); fm <- f(pop[["mu"]] - h)
    if (!all(is.finite(c(f0, fp, fm)))) return(1e12)
    f0 + 0.5 * pop[["sd"]]^2 * (fp - 2*f0 + fm) / h^2
  }, numeric(1)))

## [2] STRATIFIED: obs is the node's OWN conditional dataset, pred is that node's
obj_stratified <- function(p, obs, pops)
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
  cat(sprintf("%-40s %10s %10s   %s\n", "objective", "argmin", "bias",
              "profile spread over tcov 0.45..1.05"))
  for (r in list(
      list("[1] gh_marginal   (obs = marginal)",  obj_gh_marginal,       om),
      list("[3] taylor_moments(obs = marginal)",  obj_taylor_moments,    om),
      list("[2] stratified    (obs = per-NODE)",  obj_stratified,        oc),
      list("[4] MISMATCHED-TARGETS (the bug)  ",  obj_mismatched,        om),
      list("[4t] mismatched_taylor (old arm)  ",  obj_mismatched_taylor, om))) {
    a <- argmin(r[[2]], r[[3]], pops); s <- spread(r[[2]], r[[3]], pops)
    cat(sprintf("%-40s %10.4f %+10.4f   %s%s\n", r[[1]], a$at, a$at - TCOV,
                paste(sprintf("%8.3f", s), collapse = " "), if (a$bd) "  [boundary]" else ""))
  }
  cat("\n")
}
cat("Reading: [1], [2], [3] all target the same thing and should agree when the\n")
cat("design identifies tcov. [4] and [4t] are the category error, kept only to\n")
cat("show its size. In ONE population tcov is not identified from a single\n")
cat("marginal (E,V) -- that is a design statement, not a defect of [1].\n")
