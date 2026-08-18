## Two questions the new interface has not been asked yet.
##
##   1. MIXED margins. A real baseline table has a lognormal weight, a normal
##      score, a categorical sex. Does the Gaussian copula still do the right
##      thing, and does the QUADRATURE GRID still weight a discrete margin by
##      its own probabilities rather than by Gauss-Hermite weights?
##   2. A real fitted R-VINE through the whole chain: covDist -> grid ->
##      strata -> conditional marginalisation.
##
## Run:  Rscript validation/covariate-mixed-and-vine.R
suppressMessages(devtools::load_all(".", quiet = TRUE))
say <- function(...) { cat(..., "\n"); utils::flush.console() }
kv  <- function(k, ...) cat(sprintf("%-44s %s\n", k, paste(..., collapse = "  ")))

say("======================================================================")
say(" 1. MIXED margins, independent")
say("======================================================================")
cd <- covDist(WT   = c(mean = 72, sd = 16),                 # lognormal
              SCORE = c(mu = 0, sd = 1),                    # normal
              SEX  = c(female = 0.55, male = 0.45),         # categorical
              AGE  = list(quantile = function(u) qgamma(u, 9, 0.16)))
print(cd)
X <- covDraw(cd, n = 60000L)
say(sprintf("\n%-8s %12s %12s %12s %12s", "cov", "mean", "want", "sd", "want"))
want <- list(WT = c(72, 16), SCORE = c(0, 1), SEX = c(0.45, sqrt(.45*.55)),
             AGE = c(9/0.16, sqrt(9)/0.16))
for (n in colnames(X))
  say(sprintf("%-8s %12.4f %12.4f %12.4f %12.4f", n, mean(X[, n]), want[[n]][1],
              sd(X[, n]), want[[n]][2]))

say("\n-- the QUADRATURE GRID on the same spec --")
g <- admixr2:::.admCovGrid(cd, 5L)
m <- as.numeric(crossprod(g$W, g$X))
Xc <- sweep(g$X, 2L, m); S <- crossprod(Xc, g$W * Xc)
kv("grid rows", sprintf("%d   (5^3 x 2 levels = %d if SEX enumerates)",
                        nrow(g$X), 5L^3 * 2L))
say(sprintf("%-8s %12s %12s", "cov", "grid mean", "grid sd"))
for (j in seq_along(colnames(g$X)))
  say(sprintf("%-8s %12.4f %12.4f", colnames(g$X)[j], m[j], sqrt(S[j, j])))
say("SEX must come out at 0.45/0.4975 -- if the grid gave it Gauss-Hermite")
say("weights instead of its own probabilities it would not.")

say("\n======================================================================")
say(" 2. MIXED margins through a Gaussian COPULA")
say("======================================================================")
cd2 <- covDist(WT = c(mean = 72, sd = 16), SEX = c(female = 0.55, male = 0.45),
               cor = 0.6)
X2 <- covDraw(cd2, n = 60000L)
kv("realised cor(WT, SEX)", sprintf("%.4f  (declared 0.6 on the LATENT scale)",
                                    cor(X2[, "WT"], X2[, "SEX"])))
kv("mean WT | SEX = 0", sprintf("%.3f", mean(X2[X2[, "SEX"] == 0, "WT"])))
kv("mean WT | SEX = 1", sprintf("%.3f", mean(X2[X2[, "SEX"] == 1, "WT"])))
kv("SEX proportions preserved",
   sprintf("%.4f / %.4f  (want 0.55 / 0.45)", mean(X2[, "SEX"] == 0),
           mean(X2[, "SEX"] == 1)))
say("A discrete margin ATTENUATES Pearson correlation by construction --")
say("that is discretisation, not a bug. What must survive is the MARGIN.")
g2 <- admixr2:::.admCovGrid(cd2, 7L)
m2 <- as.numeric(crossprod(g2$W, g2$X))
kv("grid rows", nrow(g2$X))
kv("grid mean WT / SEX", sprintf("%.3f / %.4f", m2[1], m2[2]))

say("\n======================================================================")
say(" 3. a REAL fitted R-vine, end to end")
say("======================================================================")
if (!requireNamespace("rvinecopulib", quietly = TRUE)) {
  say("rvinecopulib not installed -- skipping")
} else {
  set.seed(11)
  n <- 4000L
  R <- matrix(c(1,.65,.30, .65,1,.15, .30,.15,1), 3, 3)
  Z <- matrix(rnorm(n * 3), n, 3) %*% chol(R)
  U <- pnorm(Z)
  vc <- rvinecopulib::vinecop(U, family_set = "all")
  kv("vine structure order", paste(vc$structure$order, collapse = "-"))
  kv("tree-1 families",
     paste(vapply(vc$pair_copulas[[1]], function(b) b$family, ""), collapse = ","))

  ml <- c(log(72), log(90), log(54)); sl <- c(0.22, 0.26, 0.20)
  cl <- function(x) pmin(pmax(x, 1e-12), 1 - 1e-12)
  vine_sampler <- function(u) {
    v <- rvinecopulib::inverse_rosenblatt(cl(u), vc)
    out <- vapply(1:3, function(j) qlnorm(cl(v[, j]), ml[j], sl[j]),
                  numeric(nrow(u)))
    colnames(out) <- c("WT", "CRCL", "AGE"); out
  }
  cdv <- covDist(WT = list(quantile = function(u) qlnorm(u, ml[1], sl[1])),
                 CRCL = list(quantile = function(u) qlnorm(u, ml[2], sl[2])),
                 AGE  = list(quantile = function(u) qlnorm(u, ml[3], sl[3])),
                 joint = vine_sampler)
  print(cdv)

  Xv <- covDraw(cdv, n = 40000L)
  say(sprintf("\n%-14s %10s %10s %10s", "pair", "draws", "grid 7", "vine tau"))
  prs <- list(c(1,2), c(1,3), c(2,3))
  gv <- admixr2:::.admCovGrid(cdv, 7L)
  mv <- as.numeric(crossprod(gv$W, gv$X))
  Xcv <- sweep(gv$X, 2L, mv); Sv <- crossprod(Xcv, gv$W * Xcv)
  for (pp in prs) {
    i <- pp[1]; j <- pp[2]
    say(sprintf("%-14s %10.4f %10.4f",
        paste(colnames(Xv)[c(i, j)], collapse = "-"),
        cor(Xv[, i], Xv[, j]), Sv[i, j] / sqrt(Sv[i, i] * Sv[j, j])))
  }
  kv("grid rows (7^3)", nrow(gv$X))

  say("\n-- strata from the vine: stratify WT, marginalise CRCL and AGE --")
  stv <- covStrata(cdv, stratify = "WT", n_nodes = 4L, n = 200)
  say(sprintf("%8s %10s %8s %12s %12s", "stratum", "WT", "n_k",
              "E[CRCL|k]", "E[AGE|k]"))
  for (k in seq_along(stv)) {
    Xk <- covDraw(stv[[k]]$cov_dist, n = 20000L)
    say(sprintf("%8d %10.2f %8.2f %12.3f %12.3f", k, stv[[k]]$cov$WT,
                stv[[k]]$n, mean(Xk[, "CRCL"]), mean(Xk[, "AGE"])))
  }
  kv("sum n_k", sprintf("%.4f  (want 200)", sum(vapply(stv, function(s) s$n, 0))))
  kv("unconditional means",
     sprintf("CRCL %.3f  AGE %.3f", mean(Xv[, "CRCL"]), mean(Xv[, "AGE"])))
  say("\nIf the conditional means tracked WT, the vine's dependence survived")
  say("the whole chain; if they sat at the unconditional values, it did not.")

  say("\n-- law of total expectation over the strata --")
  w  <- vapply(stv, function(s) s$n / 200, 0)
  ec <- t(vapply(stv, function(s) colMeans(
    covDraw(s$cov_dist, n = 20000L)), numeric(2)))
  kv("sum_k w_k E[CRCL|k]", sprintf("%.3f  vs marginal %.3f",
                                    sum(w * ec[, "CRCL"]), mean(Xv[, "CRCL"])))
  kv("sum_k w_k E[AGE|k]",  sprintf("%.3f  vs marginal %.3f",
                                    sum(w * ec[, "AGE"]), mean(Xv[, "AGE"])))
}
