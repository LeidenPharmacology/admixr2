## Is the importance-weight degeneracy a per-evaluation problem, or a systematic
## one?
##
## The weights are a DETERMINISTIC function of the parameters -- omega and the
## covariate coefficient both enter them -- and the draws are common random
## numbers across iterations. So the error is not noise that averages out over an
## optimisation: it is a reproducible distortion of the objective SURFACE, and a
## distorted surface has a displaced minimum.
##
## Worse, there is a feedback to worry about. Degeneracy is governed by
## sd(Delta)/omega, and BOTH of those are being estimated. If inflating omega or
## shrinking the covariate coefficient makes the weights better behaved, the
## objective improves for a reason that has nothing to do with fitting the data,
## and the optimiser will take it.
##
## This measures all three things: the per-evaluation error, whether the argmin
## moves, and which way.
##
## Run:  Rscript validation/covariate-reweighting-degeneracy.R
kv <- function(k, v) cat(sprintf("%s\t%s\n", k, paste(v, collapse = ",")))
TV <- log(10); DOSE <- 100; TT <- c(0.5, 1, 2, 4, 8)
conc <- function(cl) DOSE / exp(TV) * exp(outer(-cl / exp(TV), TT))
gh <- function(k) {
  i <- seq_len(k - 1L)
  J <- matrix(0, k, k); J[cbind(i, i + 1L)] <- sqrt(i); J[cbind(i + 1L, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE); o <- order(e$values)
  list(x = e$values[o], w = e$vectors[1L, o]^2)
}
QA <- gh(40L); QE <- gh(40L)

TCL <- 0; ADD <- 0.30; N <- 100L
## THREE populations with differing covariate means AND spreads. One population
## does not identify the coefficient at all -- the likelihood is exactly flat
## along tcl' = tcl + (b-b')*mu_a, omega'^2 = omega^2 + (b^2-b'^2)*sd_a^2 -- so on
## a single population every objective's "optimum" is a point on a ridge and a
## displacement measures nothing. (A first version of this script did that: all
## three objectives returned -209.172, the same ridge value.)
POPS <- list(c(mu = -0.45, sd = 0.30), c(mu = 0.00, sd = 0.55), c(mu = 0.50, sd = 0.35))
TCOV_TRUE <- 0.75; OM_TRUE <- 0.30

delta_of <- function(tcov, pop) tcov * (pop[["mu"]] + pop[["sd"]] * QA$x)

## ---- moments -------------------------------------------------------------
mom_exact <- function(tcov, om, pop) {
  dl <- delta_of(tcov, pop); m1 <- 0; M2 <- 0
  for (i in seq_along(dl)) for (j in seq_along(QE$x)) {
    w <- QA$w[i] * QE$w[j]
    y <- as.numeric(conc(exp(TCL + dl[i] + om * QE$x[j])))
    m1 <- m1 + w * y; M2 <- M2 + w * outer(y, y)
  }
  V <- M2 - outer(m1, m1); diag(V) <- diag(V) + ADD^2
  list(E = m1, V = V, ess = 1)
}
## reweighting from a proposal N(mu_p, s_p^2); `inflate` picks which proposal
mom_is <- function(tcov, om, inflate, pop, n = 4000L) {
  dl   <- delta_of(tcov, pop); pw <- QA$w
  mu_d <- sum(pw * dl); sd_d <- sqrt(sum(pw * (dl - mu_d)^2))
  mu_p <- if (inflate) mu_d else 0
  s_p  <- if (inflate) sqrt(om^2 + sd_d^2) else om
  z    <- qnorm((seq_len(n) - 0.5) / n)
  u    <- mu_p + s_p * z
  Y    <- conc(exp(TCL + u))
  pu   <- rowSums(vapply(seq_along(dl),
            function(i) pw[i] * dnorm(u, dl[i], om), numeric(n)))
  w    <- pu / dnorm(u, mu_p, s_p); sw <- sum(w); w <- w / sw
  E    <- as.numeric(crossprod(w, Y)); Yc <- sweep(Y, 2L, E)
  V    <- t(Yc) %*% (Yc * w); diag(V) <- diag(V) + ADD^2
  list(E = E, V = V, ess = 1 / sum(w^2) / n)
}

nll <- function(obs, pred) {
  ch <- tryCatch(chol(pred$V), error = function(e) NULL)
  if (is.null(ch)) return(1e12)
  iv <- chol2inv(ch); r <- obs$E - pred$E
  N * (2 * sum(log(diag(ch))) + sum(iv * obs$V) + as.numeric(t(r) %*% iv %*% r))
}
OBS <- lapply(POPS, function(p) mom_exact(TCOV_TRUE, OM_TRUE, p))  # exact, at truth

## total objective over the three populations
obj <- function(momfn) sum(vapply(seq_along(POPS),
  function(k) nll(OBS[[k]], momfn(POPS[[k]])), numeric(1)))
ess_of <- function(tcov, om, inflate) mean(vapply(POPS,
  function(p) mom_is(tcov, om, inflate, p)$ess, numeric(1)))
OBJ_EX <- function(tc, om) obj(function(p) mom_exact(tc, om, p))
OBJ_PL <- function(tc, om) obj(function(p) mom_is(tc, om, FALSE, p))
OBJ_IN <- function(tc, om) obj(function(p) mom_is(tc, om, TRUE,  p))

## ---- 1. per-evaluation error, and how it varies over the surface -----------
cat("Per-evaluation objective error, over the region an optimiser walks:\n")
cat(sprintf("%-7s %-7s %12s %12s %10s %8s\n",
            "tcov", "omega", "exact", "IS(plain)", "diff", "ESS"))
for (tc in c(0.30, 0.75, 1.20)) for (om in c(0.20, 0.30, 0.50)) {
  e <- OBJ_EX(tc, om); pv <- OBJ_PL(tc, om)
  cat(sprintf("%-7.2f %-7.2f %12.3f %12.3f %+10.3f %7.1f%%\n",
              tc, om, e, pv, pv - e, 100 * ess_of(tc, om, FALSE)))
}

## ---- 2. does the ARGMIN move? ----------------------------------------------
## Profile tcov with omega re-optimised at each point -- what a fit actually does.
prof <- function(OBJ) {
  f <- function(tc) stats::optimize(
    function(lo) OBJ(tc, exp(lo)), c(log(0.05), log(2)))$objective
  o <- stats::optimize(f, c(0.05, 1.60), tol = 1e-4)
  om <- stats::optimize(function(lo) OBJ(o$minimum, exp(lo)),
                        c(log(0.05), log(2)))$minimum
  c(tcov = o$minimum, omega = exp(om), obj = o$objective)
}
cat("\nProfiled optimum (omega re-optimised at each tcov):\n")
cat(sprintf("%-22s %10s %10s %12s\n", "objective", "tcov", "omega", "value"))
r_ex <- prof(OBJ_EX)
r_pl <- prof(OBJ_PL)
r_in <- prof(OBJ_IN)
cat(sprintf("%-22s %10.4f %10.4f %12.3f\n", "exact quadrature",  r_ex[1], r_ex[2], r_ex[3]))
cat(sprintf("%-22s %10.4f %10.4f %12.3f\n", "reweight, plain",   r_pl[1], r_pl[2], r_pl[3]))
cat(sprintf("%-22s %10.4f %10.4f %12.3f\n", "reweight, inflated",r_in[1], r_in[2], r_in[3]))
kv("truth", sprintf("tcov=%.4f omega=%.4f", TCOV_TRUE, OM_TRUE))
kv("displacement_plain",
   sprintf("tcov %+.4f  omega %+.4f", r_pl[1] - TCOV_TRUE, r_pl[2] - OM_TRUE))
kv("displacement_inflated",
   sprintf("tcov %+.4f  omega %+.4f", r_in[1] - TCOV_TRUE, r_in[2] - OM_TRUE))

## ---- 3. is there a FEEDBACK -- does inflating omega buy a better objective
##         through the weights rather than through the fit? --------------------
## Hold tcov at truth and walk omega; report the objective error alongside ESS.
cat("\nFeedback check: tcov held at truth, omega walked.\n")
cat(sprintf("%-7s %12s %12s %10s %8s\n", "omega", "exact", "IS(plain)", "error", "ESS"))
for (om in c(0.20, 0.30, 0.45, 0.65, 0.90, 1.30)) {
  e <- OBJ_EX(TCOV_TRUE, om); pv <- OBJ_PL(TCOV_TRUE, om)
  cat(sprintf("%-7.2f %12.3f %12.3f %+10.3f %7.1f%%\n",
              om, e, pv, pv - e, 100 * ess_of(TCOV_TRUE, om, FALSE)))
}
cat("\nIf the error becomes more NEGATIVE as omega grows, the weights are paying\n")
cat("the optimiser to inflate omega -- a systematic pull, not per-call noise.\n")
