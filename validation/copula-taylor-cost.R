## Cost and direction-dropping for the eigen-rotated marginalisation.
## Split out of copula-taylor-check.R because both needed repetition to measure:
## a single system.time() call sits under the Windows clock resolution, and the
## drop tolerance has to be set against a target ACCURACY, not picked.
##
## Run:  Rscript validation/copula-taylor-cost.R
say <- function(...) { cat(..., "\n"); utils::flush.console() }
kv  <- function(k, ...) cat(sprintf("%-30s %s\n", k, paste(..., collapse = "  ")))

gh <- function(k) { i <- seq_len(k - 1L)
  J <- matrix(0, k, k); J[cbind(i, i+1L)] <- sqrt(i); J[cbind(i+1L, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE); o <- order(e$values)
  list(x = e$values[o], w = e$vectors[1L, o]^2) }
TT <- c(0.5, 1, 2, 4, 6, 8, 12); DOSE <- 100
TCL <- log(4); TVV <- log(30); OM <- 0.30; ADD <- 0.20
GE <- gh(21L)
BCL <- c(0.75, -0.40, 0.25, 0.15); BV <- c(0.30, 0.20, 0, 0); BINT <- 0.35
NEV <- 0L
mom_cond <- function(a) {
  NEV <<- NEV + 1L; p <- length(a)
  s_cl <- sum(BCL[seq_len(p)]*a) + if (p >= 2) BINT*a[1]*a[2] else 0
  cl <- exp(TCL + s_cl + OM*GE$x); vv <- exp(TVV + sum(BV[seq_len(p)]*a))
  Y  <- DOSE/vv * exp(outer(-cl/vv, TT))
  mu <- as.numeric(crossprod(GE$w, Y)); Yc <- sweep(Y, 2L, mu)
  V  <- t(Yc) %*% (Yc*GE$w); diag(V) <- diag(V) + ADD^2
  list(E = mu, V = V) }
mom_exact <- function(muc, Sgc, nc) {
  p <- length(muc); Q <- gh(nc); L <- t(chol(Sgc))
  g <- as.matrix(expand.grid(rep(list(seq_len(nc)), p), KEEP.OUT.ATTRS = FALSE))
  w <- apply(matrix(Q$w[g], nrow(g), p), 1L, prod)
  E <- 0; M2 <- 0; EV <- 0
  for (r in seq_len(nrow(g))) {
    a <- as.numeric(muc + L %*% Q$x[g[r,]]); mm <- mom_cond(a)
    E <- E + w[r]*mm$E; M2 <- M2 + w[r]*outer(mm$E,mm$E); EV <- EV + w[r]*mm$V }
  list(E = E, V = EV + M2 - outer(E,E)) }
mom_rot <- function(muc, Sgc, hf = 0.5, tol = 0) {
  e <- eigen(Sgc, symmetric = TRUE)
  keep <- e$values > tol*max(e$values)
  lam <- e$values[keep]; Vv <- e$vectors[, keep, drop = FALSE]
  m0 <- mom_cond(muc); E <- m0$E; V <- m0$V
  for (k in seq_along(lam)) { h <- hf*sqrt(lam[k]); d <- Vv[,k]
    mp <- mom_cond(muc + h*d); mm <- mom_cond(muc - h*d)
    dE <- (mp$E - mm$E)/(2*h)
    E <- E + 0.5*lam[k]*(mp$E - 2*m0$E + mm$E)/h^2
    V <- V + 0.5*lam[k]*(mp$V - 2*m0$V + mm$V)/h^2 + lam[k]*outer(dE,dE) }
  list(E = E, V = V, ndir = length(lam)) }
mkSig <- function(p, rho) { R <- matrix(rho, p, p); diag(R) <- 1
  s <- diag(seq(0.45, 0.30, length.out = p), nrow = p); s %*% R %*% s }
relE <- function(a,b) max(abs(a-b)/pmax(abs(b),1e-12))
relV <- function(a,b) max(abs(a-b))/max(abs(b))

say("======================================================================")
say(" cost: model evaluations and wall time (repeated to clear clock noise)")
say("======================================================================")
say(sprintf("%3s %11s %11s %11s | %10s %11s %11s",
            "p", "quad NC=4", "quad NC=8", "taylor rot", "rot ms", "quad4 ms", "quad8 ms"))
for (p in 1:4) {
  muc <- rep(0.1, p); Sgc <- mkSig(p, 0.6)
  R1 <- 200L; NEV <<- 0L
  t1 <- system.time(for (i in seq_len(R1)) mom_rot(muc, Sgc))[["elapsed"]]/R1*1000
  n1 <- NEV/R1
  R2 <- max(1L, as.integer(200/4^p)); NEV <<- 0L
  t2 <- system.time(for (i in seq_len(R2)) mom_exact(muc, Sgc, 4L))[["elapsed"]]/R2*1000
  n2 <- NEV/R2
  R3 <- max(1L, as.integer(400/8^p)); NEV <<- 0L
  t3 <- system.time(for (i in seq_len(R3)) mom_exact(muc, Sgc, 8L))[["elapsed"]]/R3*1000
  n3 <- NEV/R3
  say(sprintf("%3d %11.0f %11.0f %11.0f | %10.3f %11.3f %11.3f", p, n2, n3, n1, t1, t2, t3)) }
say("\ncounts are MEASURED (mom_cond calls), not predicted")

say("\n== the count is invariant to rho, which is the claim ==")
say(sprintf("%3s %8s %10s %12s", "p", "rho", "evals", "wall ms"))
for (p in c(2L, 3L)) for (rho in c(0.0, 0.5, 0.9)) {
  muc <- rep(0.1, p); Sgc <- mkSig(p, rho); NEV <<- 0L
  tt <- system.time(for (i in 1:200) mom_rot(muc, Sgc))[["elapsed"]]/200*1000
  say(sprintf("%3d %8.2f %10.0f %12.3f", p, rho, NEV/200, tt)) }

say("\n======================================================================")
say(" dropping near-null directions: what tolerance, and what does it cost?")
say("======================================================================")
p <- 3L; muc <- rep(0.1, p)
say(sprintf("%9s %9s %28s %6s %11s %11s", "rho", "tol", "eigenvalues", "dirs",
            "relE vs ex", "relV vs ex"))
for (rho in c(0.0, 0.9, 0.99, 0.999)) {
  Sgc <- mkSig(p, rho); ex <- mom_exact(muc, Sgc, 11L)
  ev  <- eigen(Sgc, symmetric = TRUE)$values
  for (tol in c(0, 1e-4, 1e-2)) {
    r <- mom_rot(muc, Sgc, tol = tol)
    say(sprintf("%9.3f %9.0e %28s %6d %11.3e %11.3e", rho, tol,
                paste(sprintf("%.1e", ev), collapse=" "), r$ndir,
                relE(r$E, ex$E), relV(r$V, ex$V))) } }
say("\nA direction contributes in proportion to lam_k, so the sensible tolerance")
say("is set against the target accuracy. At tol=0 nothing is ever dropped and")
say("the cost stays 1+2p; the table shows what each tolerance actually buys.")
