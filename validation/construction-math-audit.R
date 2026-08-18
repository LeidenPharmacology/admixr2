## AUDIT: the mathematics of the three covariate constructions.
##
## Data are EXACT moments at the true parameters -- no sampling noise -- so a
## displaced argmin is a property of the objective, not an estimate.
##
##   S  stratified   K per-node observations vs K conditional predictions
##                     sum_k n_k NLL(obs_k, pred_k),  n_k = w_k n
##   N  node-pooled  ONE pooled observation vs K conditional predictions
##                     sum_k w_k n NLL(obs, pred_k)          <- the removed one
##   M  marginal     ONE pooled observation vs ONE marginal prediction
##                     n NLL(obs, pred_marg)
##
## The claim to test: S and M are both exact, for different DATA; N is neither.
## And S's exactness should hold at ANY node count, because the MVN kernel is
## minimised at mu = E_obs, V = V_obs, so if the model reproduces every block the
## truth minimises every term simultaneously.
##
## Run:  Rscript validation/construction-math-audit.R
say <- function(...) { cat(..., "\n"); utils::flush.console() }
kv  <- function(k, ...) cat(sprintf("%-40s %s\n", k, paste(..., collapse = "  ")))

gh <- function(k) { i <- seq_len(k - 1L)
  J <- matrix(0, k, k); J[cbind(i, i+1L)] <- sqrt(i); J[cbind(i+1L, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE); o <- order(e$values)
  list(x = e$values[o], w = e$vectors[1L, o]^2) }

DOSE <- 100; TV <- log(10); TCL <- log(1); TCOV <- 0.75; OM <- 0.30; ADD <- 0.30
TT   <- c(0.5, 1, 2, 3, 4, 6, 8)
POP  <- c(mu = 0.0, sd = 0.55)
QE   <- gh(21L)                                   # eta, held fine throughout
NSUB <- 200L

conc <- function(cl) DOSE/exp(TV) * exp(outer(-cl/exp(TV), TT))

## moments CONDITIONAL on covariate a (integrate eta only)
mom_cond <- function(a, tcl, tcov, om, add) {
  Y  <- conc(exp(tcl + tcov*a + om*QE$x))
  mu <- as.numeric(crossprod(QE$w, Y)); Yc <- sweep(Y, 2L, mu)
  V  <- t(Yc) %*% (Yc*QE$w); diag(V) <- diag(V) + add^2
  list(E = mu, V = V) }

## moments MARGINAL over (a, eta), to high accuracy
mom_marg <- function(tcl, tcov, om, add, pop, nc = 41L) {
  Q <- gh(nc); E <- 0; M2 <- 0; EV <- 0
  for (k in seq_along(Q$x)) { a <- pop[["mu"]] + pop[["sd"]]*Q$x[k]
    m <- mom_cond(a, tcl, tcov, om, add)
    E <- E + Q$w[k]*m$E; M2 <- M2 + Q$w[k]*outer(m$E, m$E); EV <- EV + Q$w[k]*m$V }
  list(E = E, V = EV + M2 - outer(E, E)) }

nll2 <- function(obs, pred, n) {
  ch <- tryCatch(chol(pred$V), error = function(e) NULL)
  if (is.null(ch)) return(1e12)
  iv <- chol2inv(ch); r <- obs$E - pred$E
  n*(2*sum(log(diag(ch))) + sum(iv*obs$V) + as.numeric(t(r) %*% iv %*% r)) }

nodes_a <- function(nc, pop) { Q <- gh(nc)
  list(a = pop[["mu"]] + pop[["sd"]]*Q$x, w = Q$w) }

## ---------------------------------------------------------------- TEST 1 ----
say("======================================================================")
say(" T1. the MVN kernel is minimised at mu = E_obs, V = V_obs")
say("======================================================================")
say("If true, S is exact at ANY node set: the truth reproduces every block, so")
say("every term attains its own global minimum simultaneously.")
o <- mom_cond(0.3, TCL, TCOV, OM, ADD)
base <- nll2(o, o, NSUB)
set.seed(1); worse <- 0L; N <- 400L
for (i in seq_len(N)) {
  pe <- o$E * (1 + 0.02*rnorm(length(o$E)))
  Pv <- o$V + 0.02*max(abs(o$V))*crossprod(matrix(rnorm(length(TT)^2), length(TT)))
  if (nll2(o, list(E = pe, V = Pv), NSUB) > base) worse <- worse + 1L }
kv("NLL(obs,obs)", sprintf("%.6f", base))
kv("random perturbations that INCREASE it", sprintf("%d / %d", worse, N))
## and the analytic minimum value
kv("analytic min  n*(log|V| + n_t)",
   sprintf("%.6f", NSUB*(determinant(o$V, logarithm = TRUE)$modulus + length(TT))))

## ---------------------------------------------------------------- TEST 2 ----
say("\n======================================================================")
say(" T2. argmin in tcov, exact data, by construction and node count")
say("======================================================================")
obs_pool <- mom_marg(TCL, TCOV, OM, ADD, POP)

objS <- function(p, nc) { nd <- nodes_a(nc, POP); tot <- 0
  for (k in seq_along(nd$a)) {
    ok <- mom_cond(nd$a[k], TCL, TCOV, OM, ADD)          # per-node DATA
    pk <- mom_cond(nd$a[k], p[1], p[2], exp(p[3]), exp(p[4]))
    tot <- tot + nll2(ok, pk, nd$w[k]*NSUB) }
  tot }
objN <- function(p, nc) { nd <- nodes_a(nc, POP); tot <- 0
  for (k in seq_along(nd$a)) {
    pk <- mom_cond(nd$a[k], p[1], p[2], exp(p[3]), exp(p[4]))
    tot <- tot + nll2(obs_pool, pk, nd$w[k]*NSUB) }      # SAME obs every node
  tot }
objM <- function(p, nc)
  nll2(obs_pool, mom_marg(p[1], p[2], exp(p[3]), exp(p[4]), POP, nc), NSUB)

argmin <- function(fn, nc) {
  prof <- function(tc) optim(c(TCL, log(OM), log(ADD)),
    function(v) fn(c(v[1], tc, v[2], v[3]), nc),
    method = "Nelder-Mead", control = list(maxit = 500, reltol = 1e-12))$value
  optimize(prof, c(0.05, 2.0), tol = 1e-4)$minimum }

say(sprintf("%6s %12s %12s %12s", "nodes", "S (strat)", "N (pooled)", "M (marginal)"))
for (nc in c(3L, 5L, 9L, 15L)) {
  say(sprintf("%6d %12.4f %12.4f %12.4f", nc,
              argmin(objS, nc), argmin(objN, nc), argmin(objM, nc))) }
say(sprintf("%6s %12.4f %12.4f %12.4f", "truth", TCOV, TCOV, TCOV))
say("\nS should be flat in the node count; M should converge; N should not go to")
say("truth at any node count -- its target is displaced, not under-resolved.")

## ---------------------------------------------------------------- TEST 3 ----
say("\n======================================================================")
say(" T3. N's mean term decomposes into a data term and a SPREAD PENALTY")
say("======================================================================")
say("   sum_k w_k r_k' Vi r_k  ==  rbar' Vi rbar + tr(Vi Cov_a(mu))")
for (nc in c(5L, 9L)) {
  nd <- nodes_a(nc, POP)
  mus <- lapply(nd$a, function(a) mom_cond(a, TCL, TCOV, OM, ADD)$E)
  Vi  <- chol2inv(chol(mom_cond(0, TCL, TCOV, OM, ADD)$V))
  lhs <- sum(vapply(seq_along(nd$a), function(k) {
    r <- obs_pool$E - mus[[k]]; nd$w[k]*as.numeric(t(r) %*% Vi %*% r) }, 0))
  mbar <- Reduce(`+`, Map(function(w, m) w*m, nd$w, mus))
  rbar <- obs_pool$E - mbar
  Cv <- Reduce(`+`, Map(function(w, m) w*outer(m - mbar, m - mbar), nd$w, mus))
  rhs <- as.numeric(t(rbar) %*% Vi %*% rbar) + sum(Vi*Cv)
  say(sprintf("  nodes=%2d  lhs %.10f  rhs %.10f  |d| %.2e   penalty share %.3f",
              nc, lhs, rhs, abs(lhs-rhs), sum(Vi*Cv)/rhs)) }
say("\nAt the true parameters the DATA term is ~0 and the penalty is the whole")
say("mean contribution -- it carries no observation, and shrinking tcov shrinks it.")

## ---------------------------------------------------------------- TEST 4 ----
say("\n======================================================================")
say(" T4. S is literally ordinary multi-study fitting")
say("======================================================================")
nc <- 5L; nd <- nodes_a(nc, POP)
studies <- lapply(seq_along(nd$a), function(k)
  list(obs = mom_cond(nd$a[k], TCL, TCOV, OM, ADD), n = nd$w[k]*NSUB, a = nd$a[k]))
p0 <- c(TCL, 0.62, log(OM), log(ADD))
viaS <- objS(p0, nc)
viaStudies <- sum(vapply(studies, function(s)
  nll2(s$obs, mom_cond(s$a, p0[1], p0[2], exp(p0[3]), exp(p0[4])), s$n), 0))
kv("objS(p0)", sprintf("%.10f", viaS))
kv("sum over K independent studies", sprintf("%.10f", viaStudies))
kv("identical?", isTRUE(all.equal(viaS, viaStudies, tolerance = 0)))
kv("sum of n_k (should be n)", sprintf("%.10f  vs %d", sum(nd$w)*NSUB, NSUB))

## ---------------------------------------------------------------- TEST 2b ---
say("\n======================================================================")
say(" T2b. M is exact too -- once the DATA and the PREDICTION use the same rule")
say("======================================================================")
say("T2's M arm scored 41-node data against nc-node predictions. That is a")
say("quadrature mismatch of its own, not a property of the marginal construction.")
objMm <- function(p, nc) nll2(mom_marg(TCL, TCOV, OM, ADD, POP, nc),
                             mom_marg(p[1], p[2], exp(p[3]), exp(p[4]), POP, nc), NSUB)
say(sprintf("%6s %14s %14s", "nodes", "M matched", "M vs 41-node data"))
for (nc in c(3L, 5L, 9L, 15L, 25L))
  say(sprintf("%6d %14.4f %14.4f", nc, argmin(objMm, nc), argmin(objM, nc)))
say(sprintf("%6s %14.4f %14.4f", "truth", TCOV, TCOV))
say("\nSame rule on both sides -> exact at every node count, exactly as for S.")
say("This is ONE principle, not three: the construction that built the data must")
say("be the construction that scores it. N is the case where it cannot be.")

## ---------------------------------------------------------------- TEST 5 ----
say("\n======================================================================")
say(" T5. the omega-inflation identity, and the hybrid built on it")
say("======================================================================")
B1 <- 0.60; B2 <- 0.45; OMU <- 0.26
S1 <- 0.50; S2 <- 0.45                       # covariate SDs (independent)
OMA <- sqrt(OMU^2 + B2^2*S2^2)               # what a source omitting a2 publishes
kv("true omega_U", sprintf("%.4f", OMU))
kv("omega_A = sqrt(om^2 + b2^2 Var(a2))", sprintf("%.4f", OMA))

mom2 <- function(a1, a2, tcl, b1, b2, om, add) {
  Y  <- conc(exp(tcl + b1*a1 + b2*a2 + om*QE$x))
  mu <- as.numeric(crossprod(QE$w, Y)); Yc <- sweep(Y, 2L, mu)
  V  <- t(Yc) %*% (Yc*QE$w); diag(V) <- diag(V) + add^2
  list(E = mu, V = V) }
## marginalise a2 out of a stratum at a1
mom2_marg <- function(a1, tcl, b1, b2, om, add, nc = 9L) {
  Q <- gh(nc); E <- 0; M2 <- 0; EV <- 0
  for (k in seq_along(Q$x)) { m <- mom2(a1, S2*Q$x[k], tcl, b1, b2, om, add)
    E <- E + Q$w[k]*m$E; M2 <- M2 + Q$w[k]*outer(m$E, m$E); EV <- EV + Q$w[k]*m$V }
  list(E = E, V = EV + M2 - outer(E, E)) }

## the identity, checked directly: marginalising a2 out == inflating omega
mq <- mom2_marg(0.4, TCL, B1, B2, OMU, ADD, 41L)
mi <- mom2(0.4, 0, TCL, B1, 0, OMA, ADD)
kv("marginal-over-a2 vs inflated-omega",
   sprintf("relE %.3e   relV %.3e", max(abs(mq$E-mi$E)/abs(mi$E)),
           max(abs(mq$V-mi$V))/max(abs(mi$V))))
say("Exact, because log(cl) is linear in a2: om*z + b2*a2 ~ N(0, om^2 + b2^2 s2^2).")
say("So source A's published blocks ARE the true model marginalised over a2 --")
say("which is why the prediction must marginalise too, and why b2 is recoverable.")

ndA <- nodes_a(7L, c(mu = 0, sd = S1))               # A stratifies a1 only
ndB <- nodes_a(5L, c(mu = 0, sd = S1))               # B stratifies both
ndB2 <- nodes_a(5L, c(mu = 0, sd = S2))
obsA <- lapply(ndA$a, function(a1) mom2(a1, 0, TCL, B1, 0, OMA, ADD))
obsB <- lapply(seq_along(ndB$a), function(i) lapply(seq_along(ndB2$a), function(j)
  mom2(ndB$a[i], ndB2$a[j], TCL, B1, B2, OMU, ADD)))

objHy <- function(v, b2, marg = TRUE) {
  tcl <- v[1]; b1 <- v[2]; om <- exp(v[3]); add <- exp(v[4]); tot <- 0
  for (k in seq_along(ndA$a)) {
    pk <- if (marg) mom2_marg(ndA$a[k], tcl, b1, b2, om, add)
          else      mom2(ndA$a[k], 0, tcl, b1, b2, om, add)   # plug-in a2 = mean
    tot <- tot + nll2(obsA[[k]], pk, ndA$w[k]*NSUB) }
  for (i in seq_along(ndB$a)) for (j in seq_along(ndB2$a))
    tot <- tot + nll2(obsB[[i]][[j]], mom2(ndB$a[i], ndB2$a[j], tcl, b1, b2, om, add),
                      ndB$w[i]*ndB2$w[j]*NSUB)
  if (!is.finite(tot)) 1e12 else tot }

prof_b2 <- function(marg) { LO <- 0.05; HI <- 1.10
  f <- function(b2) optim(c(TCL, B1, log(OMU), log(ADD)),
    function(v) objHy(v, b2, marg), method = "Nelder-Mead",
    control = list(maxit = 500, reltol = 1e-12))$value
  o <- optimize(f, c(LO, HI), tol = 5e-4)
  list(at = o$minimum, bnd = o$minimum < LO+0.02 || o$minimum > HI-0.02) }

h <- prof_b2(TRUE)
kv("hybrid argmin in b2 (truth 0.450)",
   sprintf("%.4f  bias %+.4f%s", h$at, h$at - B2, if (h$bnd) "  [AT BOUNDARY]" else ""))

## ---------------------------------------------------------------- TEST 6 ----
say("\n======================================================================")
say(" T6. the same fit with a PLUG-IN prediction for A (a2 at its mean)")
say("======================================================================")
m <- prof_b2(FALSE)
kv("plug-in argmin in b2 (truth 0.450)",
   sprintf("%.4f  bias %+.4f (%+.1f%%)%s", m$at, m$at - B2,
           100*(m$at - B2)/B2, if (m$bnd) "  [AT BOUNDARY]" else ""))
say("\nA's blocks carry a2's effect inside their inflated omega. A prediction that")
say("plugs a2 in at its mean has no term to match it, so the variance channel is")
say("broken and b2 is left identified only by source B.")

## ---------------------------------------------------------------- TEST 2c ---
say("\n======================================================================")
say(" T2c. why M's argmin wandered: it is NOT identified from one source")
say("======================================================================")
say("log(cl) = tcl + tcov*a + om*z with a ~ N(0,s^2) and z ~ N(0,1), so")
say("   tcov*a + om*z  ~  N(0, tcov^2 s^2 + om^2).")
say("The MARGINAL moments are a functional of cl's marginal law alone, so they")
say("depend on (tcov, om) ONLY through tcov^2 s^2 + om^2. Exactly flat ridge.")
s_a <- POP[["sd"]]; tot0 <- TCOV^2*s_a^2 + OM^2
ref <- mom_marg(TCL, TCOV, OM, ADD, POP, 41L)
say(sprintf("\n%10s %10s %14s %14s", "tcov", "om", "relE vs truth", "relV vs truth"))
for (tc in c(0.20, 0.40, 0.60, 0.75, 0.90)) {
  omk <- sqrt(tot0 - tc^2*s_a^2)
  m <- mom_marg(TCL, tc, omk, ADD, POP, 41L)
  say(sprintf("%10.3f %10.4f %14.3e %14.3e", tc, omk,
              max(abs(m$E-ref$E)/abs(ref$E)), max(abs(m$V-ref$V))/max(abs(ref$V)))) }
kv("objective spread along the ridge",
   sprintf("%.3e", diff(range(vapply(c(0.20,0.40,0.60,0.75,0.90), function(tc)
     nll2(ref, mom_marg(TCL, tc, sqrt(tot0 - tc^2*s_a^2), ADD, POP, 41L), NSUB), 0)))))
say("\nSo T2's and T2b's M columns measure a FLAT DIRECTION, not bias. Withdraw")
say("them as bias figures. S is identified because its per-node data carry a")
say("mean contrast across a; a single marginal source carries none.")

## ---------------------------------------------------------------- TEST 7 ----
say("\n======================================================================")
say(" T7. the paper's actual failure: stratifying a source on a covariate")
say("     its own model never fitted (a FABRICATED null contrast)")
say("======================================================================")
ndA2 <- nodes_a(5L, c(mu = 0, sd = S2))
## A's blocks, stratified on BOTH a1 and a2 -- but A's model has no b2, so every
## a2 node returns the SAME profile. That flat answer is then scored as evidence.
obsAfab <- lapply(seq_along(ndA$a), function(i) lapply(seq_along(ndA2$a), function(j)
  mom2(ndA$a[i], 0, TCL, B1, 0, OMA, ADD)))
objFab <- function(v, b2) {
  tcl <- v[1]; b1 <- v[2]; om <- exp(v[3]); add <- exp(v[4]); tot <- 0
  for (i in seq_along(ndA$a)) for (j in seq_along(ndA2$a))
    tot <- tot + nll2(obsAfab[[i]][[j]],
                      mom2(ndA$a[i], ndA2$a[j], tcl, b1, b2, om, add),
                      ndA$w[i]*ndA2$w[j]*NSUB)
  for (i in seq_along(ndB$a)) for (j in seq_along(ndB2$a))
    tot <- tot + nll2(obsB[[i]][[j]], mom2(ndB$a[i], ndB2$a[j], tcl, b1, b2, om, add),
                      ndB$w[i]*ndB2$w[j]*NSUB)
  if (!is.finite(tot)) 1e12 else tot }
pf <- function(b2) optim(c(TCL, B1, log(OMU), log(ADD)),
  function(v) objFab(v, b2), method = "Nelder-Mead",
  control = list(maxit = 500, reltol = 1e-12))$value
of <- optimize(pf, c(0.02, 1.10), tol = 5e-4)
kv("fabricated-null argmin in b2 (truth 0.450)",
   sprintf("%.4f  bias %+.4f (%+.1f%%)", of$minimum, of$minimum - B2,
           100*(of$minimum - B2)/B2))
kv("  vs hybrid (marginalise a2 for A)", sprintf("%.4f", h$at))
say("\nA's flat answer across a2 is not data -- A never fitted a2. Scoring it as")
say("evidence fights source B's genuine contrast and drags the pooled estimate")
say("down. This is the -68% mechanism, in a setting with no sampling noise.")
