## Implementation check for DEPENDENT covariates in the deterministic paths.
##
##   1. the independent case must be BIT-IDENTICAL to before (no silent drift)
##   2. the u-space product grid must reproduce the dependent distribution
##   3. the eigen-rotated Taylor design must match the closed-form Sigma
##   4. moments from grid and Taylor must agree with a large-sample reference
##
## Run:  Rscript validation/covariate-dependence-impl-check.R
suppressMessages(devtools::load_all(".", quiet = TRUE))
say <- function(...) { cat(..., "\n"); utils::flush.console() }
kv  <- function(k, ...) cat(sprintf("%-44s %s\n", k, paste(..., collapse = "  ")))

say("======================================================================")
say(" 1. independent margins: grid + Taylor design unchanged")
say("======================================================================")
ind <- list(WT = list(mu = 70, sd = 10), AGE = list(mu = 50, sd = 12))
g0  <- admixr2:::.admCovGrid(ind, 5L)
kv("grid rows / cols", sprintf("%d x %d  [%s]", nrow(g0$X), ncol(g0$X),
                               paste(colnames(g0$X), collapse = ",")))
kv("weights sum to 1", sprintf("%.17g", sum(g0$W)))
## the design points a diagonal Sigma must produce: purely axial, in order
td0 <- admixr2:::.admCovTaylorDesign(ind, 0.5)
kv("taylor design points", sprintf("%d (= 1 + 2p)", td0$n_pt))
kv("var == per-covariate variances",
   isTRUE(all.equal(unname(td0$var), c(100, 144), tolerance = 0)))
kv("directions are the coordinate axes",
   isTRUE(all.equal(td0$dir, diag(1, 2), tolerance = 0)))
kv("c sums to 1", sprintf("%.17g", sum(td0$c)))
kv("X rows", paste(apply(round(td0$X, 4), 1L, paste, collapse = "/"), collapse = "  "))

say("\n======================================================================")
say(" 2. u-space product grid reproduces a DEPENDENT distribution")
say("======================================================================")
## Gaussian copula on normal margins: the closed form is exactly known, so the
## grid can be checked against it rather than against another approximation.
RHO <- 0.7
dep <- admixr2:::.admCovDistCanon(
  list(WT  = list(mu = 70, sd = 10, dist = "normal"),
       AGE = list(mu = 50, sd = 12, dist = "normal"),
       cor = RHO))
kv("names after canon", paste(names(dep), collapse = ", "))
for (nn in c(5L, 9L, 15L)) {
  g <- admixr2:::.admCovGrid(dep, nn)
  m <- as.numeric(crossprod(g$W, g$X))
  Xc <- sweep(g$X, 2L, m); S <- crossprod(Xc, g$W * Xc)
  kv(sprintf("  nodes=%2d  rows=%3d", nn, nrow(g$X)),
     sprintf("mean %.4f/%.4f  sd %.4f/%.4f  cor %+.5f",
             m[1], m[2], sqrt(S[1,1]), sqrt(S[2,2]), S[1,2]/sqrt(S[1,1]*S[2,2])))
}
kv("truth", sprintf("mean %.4f/%.4f  sd %.4f/%.4f  cor %+.5f", 70, 50, 10, 12, RHO))

say("\n  lognormal margins (the package default), vs a 400k Sobol reference:")
depl <- admixr2:::.admCovDistCanon(
  list(WT = list(mean = 72, sd = 16), CRCL = list(mean = 90, sd = 25), cor = 0.6))
ref <- admixr2:::.admCovRowsFor(depl, 400000L, 0L)
mr  <- colMeans(ref); Sr <- crossprod(sweep(ref, 2L, mr))/nrow(ref)
for (nn in c(7L, 11L, 21L)) {
  g <- admixr2:::.admCovGrid(depl, nn)
  m <- as.numeric(crossprod(g$W, g$X))
  Xc <- sweep(g$X, 2L, m); S <- crossprod(Xc, g$W * Xc)
  kv(sprintf("  nodes=%2d", nn),
     sprintf("mean %.3f/%.3f  sd %.3f/%.3f  cor %+.4f",
             m[1], m[2], sqrt(S[1,1]), sqrt(S[2,2]), S[1,2]/sqrt(S[1,1]*S[2,2])))
}
kv("  sobol reference (4e5)",
   sprintf("mean %.3f/%.3f  sd %.3f/%.3f  cor %+.4f", mr[1], mr[2],
           sqrt(Sr[1,1]), sqrt(Sr[2,2]), Sr[1,2]/sqrt(Sr[1,1]*Sr[2,2])))
kv("  requested margins", sprintf("mean %.3f/%.3f  sd %.3f/%.3f  cor %+.4f",
                                  72, 90, 16, 25, 0.6))

say("\n======================================================================")
say(" 3. eigen-rotated Taylor design")
say("======================================================================")
td <- admixr2:::.admCovTaylorDesign(dep, 0.5)
mm <- admixr2:::.admCovDistMoments(dep)
kv("measured Sigma (normal margins, rho=0.7)",
   paste(sprintf("%.3f", as.numeric(mm$Sigma)), collapse = " "))
kv("closed form", paste(sprintf("%.3f", c(100, RHO*120, RHO*120, 144)), collapse = " "))
kv("design points", sprintf("%d (= 1 + 2p, unchanged by rho)", td$n_pt))
kv("eigenvalues", paste(sprintf("%.4f", td$var), collapse = "  "))
kv("directions orthonormal",
   isTRUE(all.equal(crossprod(td$dir), diag(1, 2), tolerance = 1e-12)))
## the design must REPRODUCE Sigma: sum_k lam_k v_k v_k' == Sigma
rec <- Reduce(`+`, lapply(seq_along(td$var),
  function(k) td$var[k] * tcrossprod(td$dir[, k])))
kv("sum lam_k v_k v_k' == Sigma", sprintf("max|d| %.2e", max(abs(rec - mm$Sigma))))
kv("c sums to 1", sprintf("%.17g", sum(td$c)))
say("\n  design points (WT, AGE) -- note they are NOT axis-aligned:")
for (i in seq_len(td$n_pt))
  say(sprintf("    %-8s %9.4f %9.4f   c = %+9.4f",
              if (i == 1L) "centre" else sprintf("dir%d%s", (i %/% 2), if (i %% 2 == 0) "+" else "-"),
              td$X[i, 1], td$X[i, 2], td$c[i]))

say("\n======================================================================")
say(" 4. do the two paths agree on the MOMENTS of a nonlinear function?")
say("======================================================================")
## g(a) = exp(0.75*log(WT/70) - 0.4*log(AGE/50)) -- an allometric-style effect,
## nonlinear in both covariates, so the constructions can actually differ.
gfun <- function(X) exp(0.75*log(X[,1]/70) - 0.40*log(X[,2]/50) + 0.15*log(X[,1]/70)*log(X[,2]/50))
ref2 <- admixr2:::.admCovRowsFor(depl, 400000L, 0L)
Eref <- mean(gfun(ref2)); Vref <- mean((gfun(ref2) - Eref)^2)
say(sprintf("%-26s %14s %14s %12s %12s", "arm", "E", "Var", "relE", "relVar"))
say(sprintf("%-26s %14.6f %14.6f %12s %12s", "sobol 4e5 (reference)", Eref, Vref, "-", "-"))
for (nn in c(5L, 7L, 11L, 21L)) {
  g <- admixr2:::.admCovGrid(depl, nn); y <- gfun(g$X)
  E <- sum(g$W*y); V <- sum(g$W*(y - E)^2)
  say(sprintf("%-26s %14.6f %14.6f %12.2e %12.2e",
              sprintf("grid, %d nodes (%d rows)", nn, nrow(g$X)), E, V,
              abs(E-Eref)/abs(Eref), abs(V-Vref)/abs(Vref))) }
## Taylor: E = sum_k c_k g_k ; Cov_a(g) = sum_k lam_k d_k d_k'
tdl <- admixr2:::.admCovTaylorDesign(depl, 0.5)
yk  <- gfun(tdl$X)
Et  <- sum(tdl$c * yk)
dk  <- (yk[tdl$ip] - yk[tdl$im]) / (2*tdl$h)
Vt  <- sum(tdl$var * dk^2)
say(sprintf("%-26s %14.6f %14.6f %12.2e %12.2e",
            sprintf("taylor (%d rows)", tdl$n_pt), Et, Vt,
            abs(Et-Eref)/abs(Eref), abs(Vt-Vref)/abs(Vref)))
say("\n(Var here is Cov_a(g) alone -- the between-subject spread the covariate")
say(" induces, which is the rank-p term. The eta-conditional part is separate.)")

say("\n======================================================================")
say(" 5. refusals that must stay refusals")
say("======================================================================")
tryerr <- function(expr) tryCatch({ expr; "ACCEPTED" },
  error = function(e) paste("refused:", substr(conditionMessage(e), 1, 78)))
kv("discrete covariate + taylor",
   tryerr(admixr2:::.admCovTaylorDesign(list(SEX = list(values = c(0,1))), 0.5)))
kv("zero-variance direction",
   tryerr(admixr2:::.admCovTaylorDesign(list(WT = list(mu = 70, sd = 0)), 0.5)))
kv("perfectly collinear covariates", tryerr({
  cd <- list(A = list(mu = 0, sd = 1, dist = "normal"),
             B = list(mu = 0, sd = 1, dist = "normal"), cor = 1 - 1e-14)
  admixr2:::.admCovTaylorDesign(cd, 0.5) }))
kv("non-PD cor", tryerr(admixr2:::.admCovDistCanon(
  list(A = list(mu=0,sd=1), B = list(mu=0,sd=1), cor = 1.4))))
