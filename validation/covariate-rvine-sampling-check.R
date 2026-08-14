## Is inverse_rosenblatt() a valid R-vine sampler, and is it CRN-safe?
##
## admixr2 needs vine draws that are a DETERMINISTIC function of the uniforms it
## supplies: the objective must be a deterministic function of the parameters, or
## a finite difference of it is noise. rvinecop() draws its own randomness and
## takes no `U`, so it cannot be used. inverse_rosenblatt() maps independent
## uniforms to vine-distributed ones, which is what rvinecop() does internally --
## but "should be equivalent" is not evidence, so this checks it.
##
## The vine is fitted the way the reference pipeline does it: kde1d marginals,
## probability-integral transform, then a FREE R-vine (Dissmann) with
## family_set = "all" and no structural constraint.
##
## Run:  Rscript validation/covariate-rvine-sampling-check.R
suppressMessages({library(rvinecopulib); library(kde1d)})
kv <- function(k, v) cat(sprintf("%s\t%s\n", k, paste(v, collapse = ",")))

set.seed(42)
n <- 4000
## five dependent covariates, so the vine has a non-trivial structure (on three
## variables every vine is simultaneously a C- and a D-vine)
mu    <- c(70, 165, 45, 95, 26)
sigma <- c(14,  10, 16, 25,  4)
R <- matrix(c(1.00, 0.62, 0.18, 0.44, 0.71,
              0.62, 1.00, 0.10, 0.30, 0.35,
              0.18, 0.10, 1.00, -0.55, 0.12,
              0.44, 0.30, -0.55, 1.00, 0.28,
              0.71, 0.35, 0.12, 0.28, 1.00), 5, 5, byrow = TRUE)
Sig <- diag(sigma) %*% R %*% diag(sigma)
Z <- sweep(matrix(rnorm(n * 5), n, 5) %*% chol(Sig), 2, mu, "+")
colnames(Z) <- c("WT", "HEIGHT", "AGE", "CRCL", "BMI")

margins <- lapply(seq_len(5), function(j) kde1d::kde1d(Z[, j]))
U <- vapply(seq_len(5), function(j) kde1d::pkde1d(Z[, j], margins[[j]]),
            numeric(n))
colnames(U) <- colnames(Z)

vc <- vinecop(U, family_set = "all")
kv("vine_order", paste(vc$structure$order, collapse = "-"))
kv("vine_trees", length(vc$structure$struct_array))
sm <- summary(vc)
kv("families", paste(unique(as.character(sm$family)), collapse = ","))
## The fit is UNCONSTRAINED -- vinecop() with no structure argument runs
## Dissmann selection over all regular vines, which is what makes it an R-vine
## fit. Which particular structure it lands on is a property of the data, not of
## the method, so the structure is reported rather than classified.
kv("tree1_edges", paste(as.character(sm$edge[sm$tree == 1]), collapse = " "))

## ---- 1. does inverse_rosenblatt reproduce the vine's own distribution? ------
N <- 200000L
set.seed(11); ref  <- rvinecop(N, vc)                        # the package sampler
set.seed(12); Uind <- matrix(runif(N * 5), N, 5)
got <- inverse_rosenblatt(Uind, vc)                          # the CRN-safe route

tau <- function(M) {
  p <- combn(5, 2)
  vapply(seq_len(ncol(p)), function(k)
    cor(M[, p[1, k]], M[, p[2, k]], method = "kendall"), numeric(1))
}
## Kendall's tau is the copula's own dependence measure: margins are uniform by
## construction, so any distributional difference has to show up here.
set.seed(1); i <- sample.int(N, 20000)   # bigger, to shrink MC error on tau
t_ref <- tau(ref[i, ]); t_got <- tau(got[i, ])
kv("kendall_tau_ref", sprintf("%+.4f", t_ref))
kv("kendall_tau_inv", sprintf("%+.4f", t_got))
kv("max_abs_tau_diff", sprintf("%.4f", max(abs(t_ref - t_got))))

ks <- vapply(seq_len(5), function(j)
  suppressWarnings(stats::ks.test(got[, j], "punif")$statistic), numeric(1))
kv("marginal_KS_vs_uniform", sprintf("%.4f", ks))
ks2 <- vapply(seq_len(5), function(j)
  suppressWarnings(stats::ks.test(sample(ref[, j], 20000),
                                  sample(got[, j], 20000))$p.value), numeric(1))
kv("KS_pvalue_ref_vs_inv", sprintf("%.3f", ks2))

## the sharpest check: the vine's own density should agree in expectation
## Mean log density on the FULL samples. This is the sharpest single number: if
## the two samplers targeted different distributions, E[log c(U)] would differ by
## the KL divergence between them. Its Monte Carlo standard error is reported so
## the comparison can be read.
lr <- log(dvinecop(ref, vc)); lg <- log(dvinecop(got, vc))
kv("mean_log_density_ref", sprintf("%.4f (se %.4f)", mean(lr), sd(lr)/sqrt(N)))
kv("mean_log_density_inv", sprintf("%.4f (se %.4f)", mean(lg), sd(lg)/sqrt(N)))
kv("log_density_diff_in_se",
   sprintf("%.2f", abs(mean(lr) - mean(lg)) / sqrt(var(lr)/N + var(lg)/N)))

## ---- 2. is it deterministic in the uniforms? (the CRN requirement) ----------
a <- inverse_rosenblatt(Uind[1:1000, ], vc)
b <- inverse_rosenblatt(Uind[1:1000, ], vc)
kv("deterministic_in_U", isTRUE(all.equal(a, b)))
kv("rvinecop_accepts_U",
   isTRUE(tryCatch({ rvinecop(10, vc, U = Uind[1:10, ]); TRUE },
                   error = function(e) FALSE)))

## ---- 3. does a Sobol stream work through it? -------------------------------
## This is how admixr2 would drive it: low-discrepancy uniforms in, vine draws out.
sob <- randtoolbox::sobol(20000L, dim = 5)
sob <- pmin(pmax(sob, 1e-9), 1 - 1e-9)
qmc <- inverse_rosenblatt(sob, vc)
kv("sobol_tau", sprintf("%+.4f", tau(qmc[sample.int(nrow(qmc), 20000), ])))
kv("sobol_vs_ref_max_tau_diff",
   sprintf("%.4f", max(abs(tau(qmc[sample.int(nrow(qmc), 20000), ]) - t_ref))))

## ---- 4. back on the covariate scale, via the same margins ------------------
to_scale <- function(V) vapply(seq_len(5), function(j)
  kde1d::qkde1d(pmin(pmax(V[, j], 1e-6), 1 - 1e-6), margins[[j]]), numeric(nrow(V)))
Xr <- to_scale(ref[i, ]); Xg <- to_scale(got[i, ])
kv("mean_ref", sprintf("%.2f", colMeans(Xr)))
kv("mean_inv", sprintf("%.2f", colMeans(Xg)))
kv("sd_ref",   sprintf("%.2f", apply(Xr, 2, sd)))
kv("sd_inv",   sprintf("%.2f", apply(Xg, 2, sd)))
kv("cor_WT_HEIGHT", sprintf("ref %.3f  inv %.3f", cor(Xr[,1], Xr[,2]), cor(Xg[,1], Xg[,2])))
kv("cor_AGE_CRCL",  sprintf("ref %.3f  inv %.3f", cor(Xr[,3], Xr[,4]), cor(Xg[,3], Xg[,4])))
