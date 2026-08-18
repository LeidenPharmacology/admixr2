## Can a vine deliver the STRATIFY / MARGINALISE split?
##
## The per-covariate rule needs, for each source and each stratified node alpha:
##     mu_c    = E[a_-S | a_S = alpha]
##     Sigma_c = Cov(a_-S | a_S = alpha)
## Under independence these collapse to the unconditional mean and diag(sd^2) --
## which is what the current covariate grid assumes. This asks whether the vine
## route really produces the CONDITIONAL versions, what the calling convention
## is, and what it costs.
##
## The convention is the whole game: fixing the wrong columns silently conditions
## on the wrong covariate and still returns a plausible answer. So it is probed
## rather than assumed.
##
## Run:  Rscript validation/copula-conditioning-check.R
suppressMessages({library(rvinecopulib)})
kv  <- function(k, ...) cat(sprintf("%-34s %s\n", k, paste(..., collapse = "  ")))
say <- function(...) { cat(..., "\n"); utils::flush.console() }
set.seed(7)
P <- 4L

## ---------------------------------------------------------------- section 1 --
## WHICH column of inverse_rosenblatt() conditions on WHICH variable?
say("======================================================================")
say(" 1. the calling convention (probed, not assumed)")
say("======================================================================")
RH <- matrix(c( 1.00,  0.65,  0.30, -0.20,
                0.65,  1.00,  0.15, -0.10,
                0.30,  0.15,  1.00,  0.40,
               -0.20, -0.10,  0.40,  1.00), P, P, byrow = TRUE)
MU <- rep(0, P); SG <- RH                       # unit margins: Sigma == R
N  <- 200000L
Z  <- matrix(rnorm(N * P), N, P) %*% chol(SG)
vg <- vinecop(pnorm(Z), family_set = "gaussian")
ORD <- vg$structure$order
kv("structure$order", paste(ORD, collapse = "-"))

probe <- function(vc, col, val = 0.99, n = 40000L) {
  set.seed(2); Uq <- matrix(runif(n * P), n, P); Uq[, col] <- val
  colMeans(qnorm(inverse_rosenblatt(Uq, vc)))
}
say(sprintf("\n%6s  %s   <- pin column to qnorm(.99)=2.326", "col",
            paste(sprintf("%8s", paste0("mean_a", 1:P)), collapse = "")))
Pr <- t(vapply(1:P, function(cc) probe(vg, cc), numeric(P)))
for (cc in 1:P) say(sprintf("%6d  %s", cc, paste(sprintf("%8.3f", Pr[cc, ]), collapse = "")))
## the root of the cascade is the column that pins its OWN variable exactly and
## moves every other by exactly rho * 2.326
root <- which.min(abs(diag(Pr) - qnorm(.99)))
kv("\ncascade root (unconditional)",
   sprintf("variable %d (pins itself at %.3f)", root, Pr[root, root]))
kv("closed form  rho_i,root * z", paste(sprintf("%+.3f", RH[, root] * qnorm(.99)), collapse = " "))
kv("observed", paste(sprintf("%+.3f", Pr[root, ]), collapse = " "))
## a column that leaves every OTHER variable at 0 is the LAST in the cascade
leafscore <- vapply(1:P, function(cc) max(abs(Pr[cc, -cc])), numeric(1))
kv("most-conditioned variable",
   sprintf("variable %d (moves others by %.3f)", which.min(leafscore), min(leafscore)))
CASC <- rev(ORD)
kv("=> cascade order", paste(CASC, collapse = " -> "))
kv("=> equals rev(structure order)?", identical(as.integer(CASC), as.integer(rev(ORD))))
kv("=> U columns indexed by", "VARIABLE (natural), not by order")
say("\nCONSEQUENCE: a stratified set S can be conditioned on natively only if it")
say("is a PREFIX of rev(structure order). Fixing any other set of columns")
say("returns a different distribution, with no error raised.")

## ---------------------------------------------------------------- section 2 --
say("\n======================================================================")
say(" 2. conditional moments: three routes to the same quantity")
say("======================================================================")
S <- CASC[1:2]; JJ <- setdiff(seq_len(P), S)
kv("stratify on (cascade prefix)", paste0("a", S, collapse = ","))
kv("marginalise",                  paste0("a", JJ, collapse = ","))

cond_closed <- function(Sig, mu, S, alpha) {
  j <- setdiff(seq_along(mu), S)
  A <- Sig[j, S, drop = FALSE] %*% solve(Sig[S, S, drop = FALSE])
  list(mu  = as.numeric(mu[j] + A %*% (alpha - mu[S])),
       Sig = Sig[j, j, drop = FALSE] - A %*% Sig[S, j, drop = FALSE]) }

## To condition on a_S = alpha (FIXED VALUES, not fixed conditional quantiles)
## the leading u must be the FORWARD Rosenblatt transform of alpha. Fixing
## u = F(alpha) componentwise is only correct for the root.
ustar <- function(vc, alpha, S) {
  a <- rep(0, P); a[S] <- alpha            # trailing entries are irrelevant:
  rosenblatt(matrix(pnorm(a), 1L), vc)[1L, ] }   # the cascade is sequential

invros_cond <- function(vc, alpha, S, n = 200000L, seed = 3L) {
  us <- ustar(vc, alpha, S)
  set.seed(seed); Uq <- matrix(runif(n * P), n, P)
  for (v in S) Uq[, v] <- us[v]
  X <- qnorm(inverse_rosenblatt(Uq, vc))
  j <- setdiff(seq_len(P), S)
  list(mu = colMeans(X)[j], Sig = cov(X[, j, drop = FALSE]),
       pin = colMeans(X)[S]) }

band_cond <- function(vc, alpha, S, n = 4000000L, eps = 0.06, seed = 5L) {
  set.seed(seed); X <- qnorm(inverse_rosenblatt(matrix(runif(n * P), n, P), vc))
  ok <- rep(TRUE, n)
  for (k in seq_along(S)) ok <- ok & abs(X[, S[k]] - alpha[k]) < eps
  Y <- X[ok, , drop = FALSE]; j <- setdiff(seq_len(P), S)
  list(mu = colMeans(Y)[j], Sig = cov(Y[, j, drop = FALSE]), n = sum(ok)) }

say(sprintf("\n%-20s %-8s %11s %11s %10s", "node a_S", "route",
            "max|d mu|", "max|d Sig|", "n_eff"))
for (al in list(c(0, 0), c(1.28, -0.52), c(-1.04, 1.04))) {
  cl <- cond_closed(SG, MU, S, al)
  ir <- invros_cond(vg, al, S); bd <- band_cond(vg, al, S)
  say(sprintf("%-20s %-8s %11.4f %11.4f %10s",
      sprintf("(%+.2f,%+.2f)", al[1], al[2]), "invros",
      max(abs(ir$mu - cl$mu)), max(abs(ir$Sig - cl$Sig)), "2e5"))
  say(sprintf("%-20s %-8s %11.4f %11.4f %10d", "", "band",
      max(abs(bd$mu - cl$mu)), max(abs(bd$Sig - cl$Sig)), bd$n))
  say(sprintf("%-20s %-8s   pinned a_S = %s", "", "",
      paste(sprintf("%+.3f", ir$pin), collapse = " ")))
}
say("(deviation from the closed-form Gaussian conditional -- MC error only)")

## the control: fix the WRONG prefix and see that nothing complains
Sbad <- CASC[3:4]
say(sprintf("\ncontrol -- stratify on a%d,a%d instead (NOT a cascade prefix),", Sbad[1], Sbad[2]))
say("fixing u = F(alpha) componentwise, the natural-looking thing to do:")
badcond <- function(vc, alpha, S, n = 200000L) {
  set.seed(3); Uq <- matrix(runif(n * P), n, P)
  for (k in seq_along(S)) Uq[, S[k]] <- pnorm(alpha[k])
  X <- qnorm(inverse_rosenblatt(Uq, vc)); j <- setdiff(seq_len(P), S)
  list(mu = colMeans(X)[j], Sig = cov(X[, j, drop = FALSE]), pin = colMeans(X)[S]) }
for (al in list(c(1.28, -0.52), c(-1.04, 1.04))) {
  cl <- cond_closed(SG, MU, Sbad, al); bb <- badcond(vg, al, Sbad)
  say(sprintf("  alpha=(%+.2f,%+.2f)  max|d mu| %.4f  max|d Sig| %.4f  (a_S landed at %s)",
      al[1], al[2], max(abs(bb$mu - cl$mu)), max(abs(bb$Sig - cl$Sig)),
      paste(sprintf("%+.2f", bb$pin), collapse = ",")))
}

## ---------------------------------------------------------------- section 3 --
say("\n======================================================================")
say(" 3. does Sigma_c depend on the stratified node?")
say("======================================================================")
NODES <- list(c(0,0), c(1.28,-0.52), c(-1.04,1.04), c(-1.64,-1.64), c(1.64,1.64))
spread <- function(vc, Sset) {
  L <- lapply(NODES, function(a) invros_cond(vc, a, Sset, n = 200000L)$Sig)
  list(mx = max(abs(Reduce(pmax, lapply(L, function(x) abs(x - L[[1]]))))),
       d11 = vapply(L, function(x) x[1,1], 0)) }
gs <- spread(vg, S)
kv("gaussian: max spread of Sigma_c", sprintf("%.4f", gs$mx))
kv("gaussian: Sigma_c[1,1] by node",  paste(sprintf("%.3f", gs$d11), collapse = " "))
kv("gaussian: closed-form Sigma_c",
   paste(sprintf("%.4f", cond_closed(SG, MU, S, c(0,0))$Sig), collapse = " "))

## same Kendall taus, different tree-1 families
TAU1 <- c(0.42, 0.20, 0.10); TAU2 <- c(0.12, 0.06); TAU3 <- 0.03
mkvine <- function(fam1) vinecop_dist(
  pair_copulas = list(
    lapply(TAU1, function(t) bicop_dist(fam1, 0, ktau_to_par(fam1, t))),
    lapply(TAU2, function(t) bicop_dist("gaussian", 0, ktau_to_par("gaussian", t))),
    lapply(TAU3, function(t) bicop_dist("gaussian", 0, ktau_to_par("gaussian", t)))),
  structure = dvine_structure(1:P))
DV <- rev(dvine_structure(1:P)$order); Sd <- DV[1:2]
kv("D-vine cascade", paste(DV, collapse = " -> "))
for (fam in c("gaussian", "clayton", "gumbel")) {
  v <- mkvine(fam); sp <- spread(v, Sd)
  kv(sprintf("%s (same taus): spread", fam), sprintf("%.4f", sp$mx))
  kv(sprintf("  %s Sigma_c[1,1] by node", fam), paste(sprintf("%.3f", sp$d11), collapse = " "))
}
say("\nGaussian: Sigma_c is CONSTANT across nodes -> one eigendecomposition for the")
say("whole grid. A general vine: Sigma_c varies -> one per node (still one-off).")

## ---------------------------------------------------------------- section 4 --
say("\n======================================================================")
say(" 4. MC scaling and cost")
say("======================================================================")
cl0 <- cond_closed(SG, MU, S, c(1.28, -0.52))
say(sprintf("%10s %13s %13s", "n", "max|d mu|", "max|d Sig|"))
for (n in c(2e4, 1e5, 5e5, 2e6)) {
  ir <- invros_cond(vg, c(1.28,-0.52), S, n = as.integer(n), seed = 17L)
  say(sprintf("%10.0f %13.5f %13.5f", n,
              max(abs(ir$mu - cl0$mu)), max(abs(ir$Sig - cl0$Sig)))) }
say("(expect halving per 4x n if the route is unbiased)")

## QMC instead of MC, since the uniforms are ours to choose
sob <- randtoolbox::sobol(200000L, dim = P); sob <- pmin(pmax(sob, 1e-9), 1 - 1e-9)
us <- ustar(vg, c(1.28, -0.52), S); for (v in S) sob[, v] <- us[v]
Xq <- qnorm(inverse_rosenblatt(sob, vg)); jq <- setdiff(seq_len(P), S)
kv("sobol n=2e5: max|d mu|",  sprintf("%.5f", max(abs(colMeans(Xq)[jq] - cl0$mu))))
kv("sobol n=2e5: max|d Sig|", sprintf("%.5f", max(abs(cov(Xq[, jq, drop = FALSE]) - cl0$Sig))))

t0 <- system.time(for (i in 1:5) invros_cond(vg, c(1.28,-0.52), S, n = 100000L))[["elapsed"]]/5
kv("seconds per node (n=1e5)", sprintf("%.3f", t0))
for (nc in c(4L, 8L))
  kv(sprintf("  grid NC=%d, |S|=2 -> %d nodes", nc, nc^2),
     sprintf("%.2f s total, ONCE", nc^2 * t0))
say("\nmu_c and Sigma_c do not depend on the model parameters, so this is paid")
say("once, before the optimiser starts -- never inside the objective.")
