## Second-order marginalisation over CORRELATED covariates: implementation and cost.
##
## The claim under test. Marginalising a_-S out of the aggregate moments, given a
## stratified node a_S = alpha, needs the CONDITIONAL mu_c and Sigma_c. To second
## order the moments are
##
##   E ~= g(mu_c)   + 1/2 tr(Sigma_c Hess g)
##   V ~= Vc(mu_c)  + 1/2 tr(Sigma_c Hess Vc) + J Sigma_c J'      J = dg/da
##
## Both Sigma_c terms are DIRECTIONAL. With Sigma_c = sum_k lam_k v_k v_k',
##
##   tr(Sigma_c Hess g) = sum_k lam_k d2g/dv_k^2
##   J Sigma_c J'       = sum_k lam_k (J v_k)(J v_k)'
##
## so differencing along the EIGENVECTORS of Sigma_c, with step h_k propto
## sqrt(lam_k), gives both terms from 1 + 2p model evaluations -- the same count
## as the independent case. Correlation should therefore be FREE. That is the
## thing to check, against a naive axial expansion that ignores the off-diagonal
## and against an explicit cross-difference stencil that does not rotate.
##
## Run:  Rscript validation/copula-taylor-check.R
say <- function(...) { cat(..., "\n"); utils::flush.console() }
kv  <- function(k, ...) cat(sprintf("%-36s %s\n", k, paste(..., collapse = "  ")))
set.seed(11)

gh <- function(k) { i <- seq_len(k - 1L)
  J <- matrix(0, k, k); J[cbind(i, i+1L)] <- sqrt(i); J[cbind(i+1L, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE); o <- order(e$values)
  list(x = e$values[o], w = e$vectors[1L, o]^2) }

## ============================================================== section 1 ====
## The algebra, on a quadratic where the central second difference is EXACT.
say("======================================================================")
say(" 1. the rotation identity, exact arithmetic (quadratic test function)")
say("======================================================================")
set.seed(3)
for (p in 2:5) {
  A <- matrix(rnorm(p*p), p, p); Sg <- crossprod(A)/p + diag(p)*0.4
  H <- crossprod(matrix(rnorm(p*p), p, p))/p          # Hessian of a quadratic
  m <- rnorm(p)
  q  <- function(a) as.numeric(m %*% a + 0.5 * t(a) %*% H %*% a)
  qv <- function(a) as.numeric(H %*% a)               # exact gradient
  e <- eigen(Sg, symmetric = TRUE); lam <- e$values; Vv <- e$vectors
  ## trace term: FD along eigenvectors vs the closed form
  a0 <- rnorm(p)
  fd <- sum(vapply(seq_len(p), function(k) { h <- 0.5*sqrt(lam[k]); d <- Vv[,k]
    lam[k]*(q(a0+h*d) - 2*q(a0) + q(a0-h*d))/h^2 }, 0))
  ## rank term: sum lam_k (J v_k)(J v_k)' vs J Sigma J'
  Jm <- H                                             # dq/da at a0 is m + H a0; J of the
  Jr <- Reduce(`+`, lapply(seq_len(p), function(k) {  # vector map a -> H a is H
    lam[k]*outer(as.numeric(Jm %*% Vv[,k]), as.numeric(Jm %*% Vv[,k])) }))
  say(sprintf("  p=%d   tr(Sigma H): rot %+.10f  closed %+.10f  |d| %.2e",
              p, fd, sum(Sg*H), abs(fd - sum(Sg*H))))
  say(sprintf("         J Sigma J' : max|rot - closed| %.3e",
              max(abs(Jr - Jm %*% Sg %*% t(Jm)))))
}
say("\n(both are pure linear algebra plus an exact stencil -- machine precision")
say(" here is a necessary condition, not evidence about the real model)")

## ============================================================== section 2 ====
## The real thing: a 1-cmt model whose covariates enter TWO parameters and
## interact, so the Hessian has genuine off-diagonal structure.
say("\n======================================================================")
say(" 2. accuracy vs an exact product-quadrature reference")
say("======================================================================")
TT <- c(0.5, 1, 2, 4, 6, 8, 12); DOSE <- 100
TCL <- log(4); TVV <- log(30); OM <- 0.30; ADD <- 0.20
GE <- gh(21L)
BCL <- c(0.75, -0.40, 0.25, 0.15); BV <- c(0.30, 0.20, 0, 0); BINT <- 0.35

## conditional moments given the covariate vector a
mom_cond <- function(a) {
  p <- length(a)
  s_cl <- sum(BCL[seq_len(p)] * a) + if (p >= 2) BINT * a[1] * a[2] else 0
  s_v  <- sum(BV[seq_len(p)] * a)
  cl <- exp(TCL + s_cl + OM * GE$x); vv <- exp(TVV + s_v)
  Y  <- DOSE/vv * exp(outer(-cl/vv, TT))
  mu <- as.numeric(crossprod(GE$w, Y)); Yc <- sweep(Y, 2L, mu)
  V  <- t(Yc) %*% (Yc * GE$w); diag(V) <- diag(V) + ADD^2
  list(E = mu, V = V) }

## EXACT marginal over a ~ N(mu_c, Sigma_c), product Gauss-Hermite
mom_exact <- function(muc, Sgc, nc) {
  p <- length(muc); Q <- gh(nc); L <- t(chol(Sgc))
  g <- as.matrix(expand.grid(rep(list(seq_len(nc)), p), KEEP.OUT.ATTRS = FALSE))
  w <- apply(matrix(Q$w[g], nrow(g), p), 1L, prod)
  E <- 0; M2 <- 0; EV <- 0
  for (r in seq_len(nrow(g))) {
    a <- as.numeric(muc + L %*% Q$x[g[r, ]]); mm <- mom_cond(a)
    E <- E + w[r]*mm$E; M2 <- M2 + w[r]*outer(mm$E, mm$E); EV <- EV + w[r]*mm$V }
  list(E = E, V = EV + M2 - outer(E, E), nev = nrow(g)) }

## (a) eigen-rotated: 1 + 2p evaluations
mom_taylor_rot <- function(muc, Sgc, hf = 0.5, tol = 1e-8) {
  e <- eigen(Sgc, symmetric = TRUE)
  keep <- e$values > tol * max(e$values)
  lam <- e$values[keep]; Vv <- e$vectors[, keep, drop = FALSE]
  m0 <- mom_cond(muc); E <- m0$E; V <- m0$V; nev <- 1L
  for (k in seq_along(lam)) {
    h <- hf*sqrt(lam[k]); d <- Vv[, k]
    mp <- mom_cond(muc + h*d); mm <- mom_cond(muc - h*d); nev <- nev + 2L
    dE <- (mp$E - mm$E)/(2*h)
    E <- E + 0.5*lam[k]*(mp$E - 2*m0$E + mm$E)/h^2
    V <- V + 0.5*lam[k]*(mp$V - 2*m0$V + mm$V)/h^2 + lam[k]*outer(dE, dE) }
  list(E = E, V = V, nev = nev, ndir = length(lam)) }

## (b) naive axial: uses only diag(Sigma_c), i.e. pretends no correlation
mom_taylor_diag <- function(muc, Sgc, hf = 0.5) {
  p <- length(muc); m0 <- mom_cond(muc); E <- m0$E; V <- m0$V; nev <- 1L
  for (k in seq_len(p)) {
    s2 <- Sgc[k, k]; h <- hf*sqrt(s2); d <- rep(0, p); d[k] <- 1
    mp <- mom_cond(muc + h*d); mm <- mom_cond(muc - h*d); nev <- nev + 2L
    dE <- (mp$E - mm$E)/(2*h)
    E <- E + 0.5*s2*(mp$E - 2*m0$E + mm$E)/h^2
    V <- V + 0.5*s2*(mp$V - 2*m0$V + mm$V)/h^2 + s2*outer(dE, dE) }
  list(E = E, V = V, nev = nev) }

## (c) axial + explicit cross stencil: no rotation, 1 + 2p + p(p-1) evaluations
mom_taylor_cross <- function(muc, Sgc, hf = 0.5) {
  p <- length(muc); h <- hf*sqrt(diag(Sgc)); m0 <- mom_cond(muc)
  ax <- lapply(seq_len(p), function(k) { d <- rep(0,p); d[k] <- 1
    list(p = mom_cond(muc + h[k]*d), m = mom_cond(muc - h[k]*d)) })
  nev <- 1L + 2L*p
  E <- m0$E; V <- m0$V
  for (k in seq_len(p)) {
    dE <- (ax[[k]]$p$E - ax[[k]]$m$E)/(2*h[k])
    E <- E + 0.5*Sgc[k,k]*(ax[[k]]$p$E - 2*m0$E + ax[[k]]$m$E)/h[k]^2
    V <- V + 0.5*Sgc[k,k]*(ax[[k]]$p$V - 2*m0$V + ax[[k]]$m$V)/h[k]^2 +
             Sgc[k,k]*outer(dE, dE) }
  if (p >= 2) for (k in 1:(p-1)) for (l in (k+1):p) {
    dk <- rep(0,p); dk[k] <- h[k]; dl <- rep(0,p); dl[l] <- h[l]
    pp <- mom_cond(muc + dk + dl); mm2 <- mom_cond(muc - dk - dl); nev <- nev + 2L
    den <- 2*h[k]*h[l]
    cE <- (pp$E + mm2$E - ax[[k]]$p$E - ax[[k]]$m$E - ax[[l]]$p$E - ax[[l]]$m$E + 2*m0$E)/den
    cV <- (pp$V + mm2$V - ax[[k]]$p$V - ax[[k]]$m$V - ax[[l]]$p$V - ax[[l]]$m$V + 2*m0$V)/den
    gk <- (ax[[k]]$p$E - ax[[k]]$m$E)/(2*h[k]); gl <- (ax[[l]]$p$E - ax[[l]]$m$E)/(2*h[l])
    E <- E + Sgc[k,l]*cE
    V <- V + Sgc[k,l]*cV + Sgc[k,l]*(outer(gk,gl) + outer(gl,gk)) }
  list(E = E, V = V, nev = nev) }

## (d) the ecological plug-in, for scale
mom_plugin <- function(muc, Sgc) { m <- mom_cond(muc); c(m, list(nev = 1L)) }

relE <- function(a, b) max(abs(a - b)/pmax(abs(b), 1e-12))
relV <- function(a, b) max(abs(a - b))/max(abs(b))

mkSig <- function(p, rho) { R <- matrix(rho, p, p); diag(R) <- 1
  s <- diag(seq(0.45, 0.30, length.out = p), nrow = p); s %*% R %*% s }

say(sprintf("%3s %6s %6s %11s %11s %11s %11s %11s %11s %6s",
            "p", "rho", "arm", "relE", "relV", "", "", "", "", "nev"))
for (p in 2:3) for (rho in c(0.0, 0.5, 0.85)) {
  muc <- rep(0.1, p); Sgc <- mkSig(p, rho)
  ex <- mom_exact(muc, Sgc, if (p == 2) 21L else 13L)
  for (nm in c("rot", "cross", "diag", "plugin")) {
    r <- switch(nm, rot = mom_taylor_rot(muc, Sgc), cross = mom_taylor_cross(muc, Sgc),
                diag = mom_taylor_diag(muc, Sgc), plugin = mom_plugin(muc, Sgc))
    say(sprintf("%3d %6.2f %6s %11.3e %11.3e %47s %6d",
                p, rho, nm, relE(r$E, ex$E), relV(r$V, ex$V), "", r$nev)) }
  say(sprintf("%3s %6s %6s reference: product GH, %d nodes", "", "", "exact", ex$nev)) }
say("\n(rot and cross both carry the full Sigma_c; diag drops the off-diagonal;")
say(" plugin is evaluation at the conditional mean)")

## ============================================================== section 3 ====
say("\n======================================================================")
say(" 3. cost: evaluations and wall time")
say("======================================================================")
say(sprintf("%3s %10s %10s %10s %10s %12s %12s", "p", "quad NC=4", "quad NC=8",
            "taylor rot", "taylor cross", "rot ms", "quad4 ms"))
for (p in 1:4) {
  muc <- rep(0.1, p); Sgc <- mkSig(p, 0.6)
  t_rot <- system.time(for (i in 1:20) mom_taylor_rot(muc, Sgc))[["elapsed"]]/20*1000
  t_q4  <- system.time(mom_exact(muc, Sgc, 4L))[["elapsed"]]*1000
  say(sprintf("%3d %10d %10d %10d %12d %12.1f %12.1f", p, 4^p, 8^p,
              1L+2L*p, 1L+2L*p+p*(p-1L), t_rot, t_q4)) }
say("\nthe rotated count is 1+2p at EVERY rho -- correlation adds no evaluations")

## ============================================================== section 4 ====
say("\n======================================================================")
say(" 4. near-collinear covariates: directions should drop out")
say("======================================================================")
for (rho in c(0.0, 0.9, 0.99, 0.999, 0.99999)) {
  Sgc <- mkSig(3L, rho); r <- mom_taylor_rot(rep(0.1,3), Sgc)
  ev <- eigen(Sgc, symmetric = TRUE)$values
  kv(sprintf("rho=%.5f", rho),
     sprintf("eigenvalues %s  dirs kept %d  nev %d",
             paste(sprintf("%.2e", ev), collapse=" "), r$ndir, r$nev)) }
say("\nStrongly correlated covariates REDUCE the effective dimension, so they")
say("cost less, not more. The tolerance is relative to the largest eigenvalue.")

## ============================================================== section 5 ====
say("\n======================================================================")
say(" 5. the failure mode: marginalising UNCONDITIONALLY when rho != 0")
say("======================================================================")
## Two covariates: stratify on a1 (the source fitted it), marginalise a2.
## Correct: a2 ~ N(mu_c(a1), s2^2(1-rho^2)). Current grid: a2 ~ N(mu2, s2^2).
S1 <- 0.45; S2 <- 0.40; MU2 <- 0.0
say(sprintf("%6s %10s %12s %12s %12s", "rho", "node a1", "relE cond",
            "relE uncond", "ratio"))
for (rho in c(0.0, 0.3, 0.6, 0.85)) {
  for (a1 in c(-0.9, 0.0, 0.9)) {
    muc <- c(a1, MU2 + rho*(S2/S1)*(a1 - 0.0))
    Sgc <- matrix(c(1e-12, 0, 0, S2^2*(1 - rho^2)), 2, 2)
    ex  <- mom_exact(muc, Sgc + diag(c(1e-12, 0)), 21L)
    rc  <- mom_taylor_rot(muc, Sgc)
    muu <- c(a1, MU2); Sgu <- matrix(c(1e-12, 0, 0, S2^2), 2, 2)
    ru  <- mom_taylor_rot(muu, Sgu)
    say(sprintf("%6.2f %10.2f %12.3e %12.3e %12.1f", rho, a1,
                relE(rc$E, ex$E), relE(ru$E, ex$E),
                relE(ru$E, ex$E)/max(relE(rc$E, ex$E), 1e-15))) } }
say("\nAt rho=0 the two coincide. As rho grows the unconditional route predicts")
say("the average-a2 response at every a1 node, while the source's high-a1")
say("subjects actually had high a2 -- a mean error that VARIES ACROSS NODES,")
say("which is the shape that biases the coefficient on the stratified covariate.")
