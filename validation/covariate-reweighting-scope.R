## Which mu-references and functional forms does the importance-reweighting
## construction actually support?
##
## The whole scheme rests on ONE property: the prediction must depend on the
## covariate and the random effect only through their SUM,
##
##     f(a, eta) = F(u),      u = Delta(a) + eta
##
## because then a fixed ensemble of solves indexed by u can be reweighted from
## the eta density to the u density. Everything else -- which link, which
## functional form, how many covariates -- is downstream of that one question.
##
## So this script tests that property DIRECTLY, per model form, rather than
## inferring it from the model text (a syntactic check on `exp(` and the
## covariate name accepts `cl <- exp(tcl + eta) + tcov*WT`, where it is false).
## The gate is: pick (a, eta) pairs with equal u and check the predictions match.
##
## Run:  Rscript validation/covariate-reweighting-scope.R
kv <- function(k, v) cat(sprintf("%s\t%s\n", k, paste(v, collapse = ",")))
TV <- log(10); DOSE <- 100; TT <- c(0.5, 1, 2, 4, 8)
conc <- function(cl, v = exp(TV)) DOSE / v * exp(outer(-cl / v, TT))

gh <- function(k) {
  i <- seq_len(k - 1L)
  J <- matrix(0, k, k); J[cbind(i, i + 1L)] <- sqrt(i); J[cbind(i + 1L, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE); o <- order(e$values)
  list(x = e$values[o], w = e$vectors[1L, o]^2)
}
QA <- gh(40L); QE <- gh(40L)
OM <- 0.30
## covariate: lognormal weight, the realistic case
A_ML <- log(70); A_SL <- 0.25
a_nodes <- exp(A_ML + A_SL * QA$x); a_w <- QA$w
A_REF   <- exp(A_ML + A_SL^2 / 2)

## ---- the model forms ------------------------------------------------------
## each: pred(a, eta) -> vector over TT, and the Delta the scheme would use
FORMS <- list(
  "A log-linear   exp(t + b*a + eta)" = list(
    pred  = function(a, e) conc(exp(0 + 0.010 * a + e)),
    delta = function(a) 0.010 * a),
  "B allometric   exp(t+eta)*(a/70)^b" = list(
    pred  = function(a, e) conc(exp(0 + e) * (a / 70)^0.75),
    delta = function(a) 0.75 * log(a / 70)),
  "C arbitrary h  exp(t+eta)*h(a)" = list(
    pred  = function(a, e) conc(exp(0 + e) * (1 + 0.4 * sin(a / 10))),
    delta = function(a) log(1 + 0.4 * sin(a / 10))),
  "D identity     t + b*a + eta" = list(
    pred  = function(a, e) conc(1.0 + 0.004 * a + e),
    delta = function(a) 0.004 * a),
  "E expit        4*expit(t+b*a+eta)" = list(
    pred  = function(a, e) conc(4 / (1 + exp(-(-1 + 0.02 * a + e)))),
    delta = function(a) 0.02 * a),
  "F ADDITIVE     exp(t+eta) + b*a" = list(
    pred  = function(a, e) conc(exp(0 + e) + 0.006 * a),
    delta = function(a) log((exp(0) + 0.006 * a) / exp(0))),
  "G cov ALSO on v (no eta)" = list(
    pred  = function(a, e) conc(exp(0 + e) * (a / 70)^0.75,
                                v = exp(TV) * (a / 70)),
    delta = function(a) 0.75 * log(a / 70)),
  "H two covariates, one eta" = list(
    pred  = function(a, e) conc(exp(0 + e) * (a / 70)^0.75),
    delta = function(a) 0.75 * log(a / 70)),
  ## the mg/kg dose -- the covariate enters the EVENT, not a parameter. This is
  ## the PAGE case study's shape, and the reason that model is only partly
  ## absorbable.
  "I mg/kg dose   f(centr) <- a" = list(
    pred  = function(a, e) a * conc(exp(0 + e) * (a / 70)^0.75),
    delta = function(a) 0.75 * log(a / 70)),
  ## a covariate whose eta is SHARED with another parameter: shifting that eta
  ## moves v as well as cl, which the covariate does not.
  "J eta shared by cl and v" = list(
    pred  = function(a, e) conc(exp(0 + e) * (a / 70)^0.75, v = exp(TV + e)),
    delta = function(a) 0.75 * log(a / 70)))

## ---- THE GATE: does f depend on (a, eta) only through u = Delta(a)+eta? -----
## Pick several (a, eta) with a common u and compare against the reference point.
gate <- function(fm) {
  worst <- 0
  for (u in c(-0.4, -0.1, 0.15, 0.5)) {
    ref <- fm$pred(A_REF, u - fm$delta(A_REF))
    for (a in quantile(a_nodes, c(0.05, 0.3, 0.7, 0.95))) {
      got <- fm$pred(a, u - fm$delta(a))
      worst <- max(worst, max(abs(got - ref) / abs(ref)))
    }
  }
  worst
}

## ---- exact nested quadrature ------------------------------------------------
exact <- function(fm) {
  m1 <- 0; M2 <- 0
  for (i in seq_along(a_nodes)) for (j in seq_along(QE$x)) {
    w <- a_w[i] * QE$w[j]
    y <- as.numeric(fm$pred(a_nodes[i], OM * QE$x[j]))
    m1 <- m1 + w * y; M2 <- M2 + w * outer(y, y)
  }
  list(E = m1, V = M2 - outer(m1, m1))
}

## ---- reweighting, with the collapse's inflated covariance as the proposal ---
reweight <- function(fm, n = 4000L) {
  dl   <- fm$delta(a_nodes)
  mu_d <- sum(a_w * dl); sd_d <- sqrt(sum(a_w * (dl - mu_d)^2))
  s_p  <- sqrt(OM^2 + sd_d^2)
  z    <- qnorm((seq_len(n) - 0.5) / n)
  u    <- mu_d + s_p * z
  ## ONE ensemble, indexed by u, evaluated at the reference covariate
  Y    <- t(vapply(u, function(uu) as.numeric(fm$pred(A_REF, uu - fm$delta(A_REF))),
                   numeric(length(TT))))
  pu   <- rowSums(vapply(seq_along(dl),
            function(i) a_w[i] * dnorm(u, dl[i], OM), numeric(n)))
  w    <- pu / dnorm(u, mu_d, s_p); w <- w / sum(w)
  E    <- as.numeric(crossprod(w, Y)); Yc <- sweep(Y, 2L, E)
  list(E = E, V = t(Yc) %*% (Yc * w), ess = 1 / sum(w^2) / n)
}

cat(sprintf("%-36s %10s %10s %10s %7s\n",
            "model form", "GATE", "E rel", "V rel", "ESS"))
for (nm in names(FORMS)) {
  fm <- FORMS[[nm]]
  g  <- gate(fm)
  r  <- reweight(fm); x <- exact(fm)
  cat(sprintf("%-36s %10.2e %10.2e %10.2e %6.1f%%%s\n", nm, g,
              max(abs(r$E - x$E) / abs(x$E)),
              max(abs(r$V - x$V)) / max(abs(x$V)),
              100 * r$ess,
              if (g > 1e-8) "   <- GATE FAILS" else ""))
}
cat("\nGATE is max relative disagreement between predictions that share a u.\n")
cat("It is the numerical check that decides support; a form that fails it is\n")
cat("not reweightable no matter how the model text reads.\n")
