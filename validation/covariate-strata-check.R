## Does .admCovStrata() really condition, or only look like it?
##
## The claim: within stratum k the marginalised covariates follow their
## distribution CONDITIONAL on that stratum. For a Gaussian copula the answer is
## known in closed form, so this is a real check rather than a comparison
## against another approximation.
##
## Run:  Rscript validation/covariate-strata-check.R
suppressMessages(devtools::load_all(".", quiet = TRUE))
say <- function(...) { cat(..., "\n"); utils::flush.console() }
kv  <- function(k, ...) cat(sprintf("%-46s %s\n", k, paste(..., collapse = "  ")))

say("======================================================================")
say(" 1. weights, nodes and shape")
say("======================================================================")
cd <- admixr2:::.admCovDistCanon(
  list(WT   = list(mu = 70, sd = 10, dist = "normal"),
       CRCL = list(mu = 90, sd = 20, dist = "normal"), cor = 0.7))
st <- admixr2:::.admCovStrata(cd, stratify = "WT", n_nodes = 4L)
kv("strata", length(st))
kv("weights sum to 1", sprintf("%.17g", sum(vapply(st, `[[`, 0, "weight"))))
kv("each stratum pins WT", all(vapply(st, function(s) "WT" %in% names(s$cov), TRUE)))
kv("each stratum marginalises CRCL",
   all(vapply(st, function(s) identical(names(s$cov_dist)[1L], "CRCL"), TRUE)))
say(sprintf("\n%8s %10s %10s", "stratum", "WT node", "weight"))
for (k in seq_along(st))
  say(sprintf("%8d %10.3f %10.4f", k, st[[k]]$cov$WT, st[[k]]$weight))

say("\n======================================================================")
say(" 2. is the within-stratum CRCL distribution CONDITIONAL?")
say("======================================================================")
## closed form for a bivariate normal: E[CRCL | WT = w] = mu2 + rho*s2/s1*(w-mu1)
## sd(CRCL | WT) = s2*sqrt(1-rho^2), CONSTANT across strata for a Gaussian copula
RHO <- 0.7; MU1 <- 70; S1 <- 10; MU2 <- 90; S2 <- 20
say(sprintf("%8s %10s %14s %14s %12s %12s", "stratum", "WT node",
            "E[CRCL|WT] got", "closed form", "sd got", "closed form"))
for (k in seq_along(st)) {
  X <- covDraw(st[[k]]$cov_dist, n = 40000L)
  w <- st[[k]]$cov$WT
  say(sprintf("%8d %10.3f %14.4f %14.4f %12.4f %12.4f", k, w,
              mean(X[, "CRCL"]), MU2 + RHO * S2 / S1 * (w - MU1),
              sd(X[, "CRCL"]), S2 * sqrt(1 - RHO^2)))
}
say("\nIf these matched the UNCONDITIONAL mean (90.0) at every stratum, the")
say("conditioning would be the no-op that biases the stratified coefficient.")

say("\n======================================================================")
say(" 3. the law of total expectation, as an end-to-end check")
say("======================================================================")
## sum_k w_k E[a | stratum k] must return the MARGINAL mean, and
## sum_k w_k (Var_k + (E_k - E)^2) must return the marginal variance.
em <- vapply(st, function(s) mean(covDraw(s$cov_dist, n = 40000L)[, "CRCL"]), 0)
vm <- vapply(st, function(s) { X <- covDraw(s$cov_dist, n = 40000L)
  mean((X[, "CRCL"] - mean(X[, "CRCL"]))^2) }, 0)
ww <- vapply(st, `[[`, 0, "weight")
kv("sum_k w_k E_k          (want 90)", sprintf("%.4f", sum(ww * em)))
kv("total variance          (want 400)",
   sprintf("%.4f", sum(ww * vm) + sum(ww * (em - sum(ww * em))^2)))
kv("  within-stratum share", sprintf("%.4f", sum(ww * vm)))
kv("  between-stratum share", sprintf("%.4f", sum(ww * (em - sum(ww * em))^2)))
kv("  = rho^2 * 400 between (want 196)", sprintf("%.1f", RHO^2 * S2^2))

say("\n======================================================================")
say(" 4. independent covariates: conditioning must be a NO-OP")
say("======================================================================")
cdi <- list(WT = list(mu = 70, sd = 10), CRCL = list(mu = 90, sd = 20))
sti <- admixr2:::.admCovStrata(cdi, stratify = "WT", n_nodes = 4L)
emi <- vapply(sti, function(s) mean(covDraw(s$cov_dist, n = 40000L)[, "CRCL"]), 0)
kv("E[CRCL | stratum] across strata", paste(sprintf("%.3f", emi), collapse = "  "))
kv("all equal to the marginal 90?",
   sprintf("max deviation %.4f", max(abs(emi - 90))))

say("\n======================================================================")
say(" 5. a NON-Gaussian sampler (Clayton), conditioned the same way")
say("======================================================================")
## Clayton's conditional inverse is closed form, and its conditional SPREAD
## depends on where you are -- lower-tail dependence. A Gaussian conditional
## would have a constant spread, so this distinguishes them.
clayton <- function(th) function(u) {
  u1 <- u[, 1]
  u2 <- (pmax(u[, 2]^(-th / (1 + th)) - 1, 0) * u1^(-th) + 1)^(-1 / th)
  cbind(A = qnorm(u1), B = qnorm(u2))
}
cdc <- list(A = list(quantile = function(u) qnorm(u)),
            B = list(quantile = function(u) qnorm(u)), joint = clayton(2))
stc <- admixr2:::.admCovStrata(cdc, stratify = "A", n_nodes = 5L)
say(sprintf("%8s %10s %14s %12s", "stratum", "A node", "E[B|A]", "sd(B|A)"))
for (k in seq_along(stc)) {
  X <- covDraw(stc[[k]]$cov_dist, n = 40000L)
  say(sprintf("%8d %10.3f %14.4f %12.4f", k, stc[[k]]$cov$A,
              mean(X[, "B"]), sd(X[, "B"])))
}
say("\nClayton concentrates in the LOWER tail, so sd(B|A) should be SMALLER at")
say("low A and larger at high A -- a Gaussian copula gives a constant spread.")

say("\n======================================================================")
say(" 6. a DISCRETE stratified covariate enumerates its levels")
say("======================================================================")
cdd <- list(SEX = list(values = c(0, 1), probs = c(0.6, 0.4)),
            WT  = list(mu = 70, sd = 10))
std <- admixr2:::.admCovStrata(cdd, stratify = "SEX", n_nodes = 7L)
kv("strata (want 2, not 7)", length(std))
kv("levels", paste(vapply(std, function(s) s$cov$SEX, 0), collapse = "  "))
kv("weights", paste(sprintf("%.3f", vapply(std, `[[`, 0, "weight")), collapse = "  "))
