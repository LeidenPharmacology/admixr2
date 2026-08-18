## A CORRECT second-order Taylor approximation of the covariate integral.
##
## The Taylor arm used in `covariate-threeway.R` and `covariate-integration-
## comparison.R` expands the WRONG FUNCTIONAL. It computes
##
##     NLL(mu_a) + 0.5 sd^2 d2/da2 NLL(a)          ~=  E_a[ NLL(a) ]
##
## i.e. a second-order approximation to the `gh` target -- the EXPECTATION OF THE
## OBJECTIVE over the covariate. That target is the stratified one, and it is
## displaced by about -0.45 in tcov on its own. So the old Taylor arm's error is
## overwhelmingly its target, not its truncation, and it cannot be used to say
## anything about Taylor expansion as a technique.
##
## What should be expanded is the MOMENTS, which is what the aggregate likelihood
## actually consumes. With g(a) = E_eta[f(a,eta)] and Vc(a) = Cov_eta(f | a),
##
##     E_marg = E_a[g(a)]           ~=  g(mu) + (sd^2/2) g''(mu)
##     V_marg = E_a[Vc(a)] + Cov_a(g(a))
##                                  ~=  Vc(mu) + (sd^2/2) Vc''(mu) + sd^2 g'(mu) g'(mu)^T
##
## then the NLL is evaluated ONCE at those approximate marginal moments. Three
## solves instead of K nodes.
##
## Note what the last term is. `sd^2 g' g'^T` is the covariate-induced between-
## subject covariance -- to second order it IS the omega'^2 = omega^2 + beta^2 Var(a)
## inflation, appearing here as a rank-one addition to V. The identity the paper
## inverts is the variance term of this expansion.
setwd("C:/package/admixr2/.claude/worktrees/feature-covariate-quadrature")
say <- function(...) { cat(..., "\n"); utils::flush.console() }

TCL <- log(1); TV <- log(10); TCOV <- 0.75; OM <- 0.30; ADD <- 0.30; DOSE <- 100
POPS <- list(c(mu = -0.45, sd = 0.30), c(mu = 0.00, sd = 0.55), c(mu = 0.50, sd = 0.35))
RICH   <- c(0.5, 1, 1.5, 2, 3, 4, 5, 6, 8)
SPARSE <- c(6, 8)
gh_nodes <- function(k) { i <- seq_len(k - 1L)
  J <- matrix(0, k, k); J[cbind(i, i + 1L)] <- sqrt(i); J[cbind(i + 1L, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE); o <- order(e$values)
  list(x = e$values[o], w = e$vectors[1L, o]^2) }
QE <- gh_nodes(21L); QA <- gh_nodes(21L)
conc <- function(cl, tt) DOSE / exp(TV) * exp(outer(-cl / exp(TV), tt))

mom_cond <- function(a, tcl, tcov, om, add, tt) {
  Y  <- conc(exp(tcl + tcov * a + om * QE$x), tt)
  mu <- as.numeric(crossprod(QE$w, Y)); Yc <- sweep(Y, 2L, mu)
  V  <- t(Yc) %*% (Yc * QE$w); diag(V) <- diag(V) + add^2
  list(E = mu, V = V) }
mom_marg <- function(tcl, tcov, om, add, tt, pop) {
  m1 <- 0; M2 <- 0
  for (ia in seq_along(QA$x)) { a <- pop[["mu"]] + pop[["sd"]] * QA$x[ia]
    Y <- conc(exp(tcl + tcov * a + om * QE$x), tt)
    for (ib in seq_along(QE$x)) { w <- QA$w[ia] * QE$w[ib]; y <- Y[ib, ]
      m1 <- m1 + w * y; M2 <- M2 + w * outer(y, y) } }
  V <- M2 - outer(m1, m1); diag(V) <- diag(V) + add^2
  list(E = m1, V = V) }

## ---- the corrected expansion: second order in the MOMENTS --------------------
mom_taylor <- function(tcl, tcov, om, add, tt, pop, hfrac = 0.5) {
  mu <- pop[["mu"]]; s <- pop[["sd"]]; h <- hfrac * s
  m0 <- mom_cond(mu,     tcl, tcov, om, add, tt)
  mp <- mom_cond(mu + h, tcl, tcov, om, add, tt)
  mm <- mom_cond(mu - h, tcl, tcov, om, add, tt)
  dE  <- (mp$E - mm$E) / (2 * h)              # g'(mu)
  d2E <- (mp$E - 2 * m0$E + mm$E) / h^2       # g''(mu)
  d2V <- (mp$V - 2 * m0$V + mm$V) / h^2       # Vc''(mu); add^2 is constant and cancels
  list(E = m0$E + 0.5 * s^2 * d2E,
       V = m0$V + 0.5 * s^2 * d2V + s^2 * outer(dE, dE)) }

nll2 <- function(obs, pred, n) { ch <- tryCatch(chol(pred$V), error = function(e) NULL)
  if (is.null(ch)) return(1e12); iv <- chol2inv(ch); r <- obs$E - pred$E
  n * (2*sum(log(diag(ch))) + sum(iv * obs$V) + as.numeric(t(r) %*% iv %*% r)) }

obj_marginal <- function(p, obs, n, tt) sum(vapply(seq_along(POPS), function(k)
  nll2(obs[[k]], mom_marg(p[1],p[2],exp(p[3]),exp(p[4]), tt, POPS[[k]]), n), numeric(1)))
obj_gh <- function(p, obs, n, tt) { tot <- 0
  for (k in seq_along(POPS)) { pop <- POPS[[k]]
    tot <- tot + sum(QA$w * vapply(pop[["mu"]] + pop[["sd"]] * QA$x, function(a)
      nll2(obs[[k]], mom_cond(a, p[1],p[2],exp(p[3]),exp(p[4]), tt), n), numeric(1))) }
  tot }
## the OLD arm, kept verbatim for the comparison
obj_taylor_old <- function(p, obs, n, tt, hfrac = 0.5) { tot <- 0
  for (k in seq_along(POPS)) { pop <- POPS[[k]]; h <- hfrac * pop[["sd"]]
    f0 <- nll2(obs[[k]], mom_cond(pop[["mu"]],     p[1],p[2],exp(p[3]),exp(p[4]), tt), n)
    fp <- nll2(obs[[k]], mom_cond(pop[["mu"]] + h, p[1],p[2],exp(p[3]),exp(p[4]), tt), n)
    fm <- nll2(obs[[k]], mom_cond(pop[["mu"]] - h, p[1],p[2],exp(p[3]),exp(p[4]), tt), n)
    if (!all(is.finite(c(f0, fp, fm)))) return(1e12)
    tot <- tot + f0 + 0.5 * pop[["sd"]]^2 * (fp - 2*f0 + fm) / h^2 }
  tot }
obj_taylor_new <- function(p, obs, n, tt, hfrac = 0.5) sum(vapply(seq_along(POPS),
  function(k) nll2(obs[[k]], mom_taylor(p[1],p[2],exp(p[3]),exp(p[4]), tt, POPS[[k]], hfrac), n),
  numeric(1)))

LO <- 0.05; HI <- 2.20
argmin_tcov <- function(objfn, obs, n, tt) {
  prof <- function(tc) optim(c(TCL, log(OM), log(ADD)),
      function(v) objfn(c(v[1], tc, v[2], v[3]), obs, n, tt),
      method = "Nelder-Mead", control = list(maxit = 400, reltol = 1e-11))$value
  o <- optimize(prof, c(LO, HI), tol = 5e-4)
  list(at = o$minimum, boundary = (o$minimum < LO + 0.02 || o$minimum > HI - 0.02)) }

say(sprintf("truth tcov = %.4f, exact marginal data, 3 identifying populations\n", TCOV))
say(sprintf("%-7s %5s %10s %10s %10s %10s", "design","n","marginal","gh","taylor-OLD","taylor-NEW"))
for (dn in c("rich","sparse")) { tt <- if (dn=="rich") RICH else SPARSE
  obs <- lapply(POPS, function(pp) mom_marg(TCL,TCOV,OM,ADD,tt,pp))
  for (n in c(10L, 1000L)) {
    a1<-argmin_tcov(obj_marginal,obs,n,tt); a2<-argmin_tcov(obj_gh,obs,n,tt)
    a3<-argmin_tcov(obj_taylor_old,obs,n,tt); a4<-argmin_tcov(obj_taylor_new,obs,n,tt)
    say(sprintf("%-7s %5d %10.4f %10.4f %10.4f %10.4f", dn,n,a1$at,a2$at,a3$at,a4$at)) } }
say("\nbias against truth 0.7500 -- OLD expands E_a[NLL] (the stratified target),")
say("NEW expands the marginal MOMENTS (the target the likelihood actually consumes)")

## ---- how good is the corrected expansion, as a function of regime? ----------
## ratio = tcov*sd_a / om : covariate-induced spread relative to IIV.
say("\n=========== moment accuracy vs the exact marginal, by regime ===========")
say(sprintf("%8s %6s %14s %14s %14s %14s", "ratio","hfrac","relE taylor","relV taylor","relE cond","relV cond"))
relerr <- function(a,b) max(abs(a-b)/pmax(abs(b),1e-12))
POP1 <- c(mu = 0.0, sd = 0.40)
for (ratio in c(0.1,0.25,0.5,1.0,1.5,2.0,3.0)) {
  tcov <- ratio * OM / POP1[["sd"]]
  ex <- mom_marg(TCL,tcov,OM,ADD,RICH,POP1)
  cd <- mom_cond(POP1[["mu"]],TCL,tcov,OM,ADD,RICH)        # the plug-in, for scale
  for (hf in c(0.5)) {
    ta <- mom_taylor(TCL,tcov,OM,ADD,RICH,POP1,hf)
    say(sprintf("%8.2f %6.2f %14.3e %14.3e %14.3e %14.3e", ratio, hf,
      relerr(ta$E,ex$E), relerr(ta$V,ex$V), relerr(cd$E,ex$E), relerr(cd$V,ex$V))) } }
say("\n('cond' = evaluating at the mean covariate, i.e. the ecological plug-in)")

## ---- sensitivity to the differencing step -----------------------------------
say("\n=========== corrected expansion: sensitivity to hfrac (h = hfrac*sd) ===========")
say(sprintf("%8s %14s %14s", "hfrac","relE","relV"))
tcov1 <- 1.0 * OM / POP1[["sd"]]
ex1 <- mom_marg(TCL,tcov1,OM,ADD,RICH,POP1)
for (hf in c(0.05,0.1,0.25,0.5,0.75,1.0)) { ta <- mom_taylor(TCL,tcov1,OM,ADD,RICH,POP1,hf)
  say(sprintf("%8.2f %14.3e %14.3e", hf, relerr(ta$E,ex1$E), relerr(ta$V,ex1$V))) }
say("\n(at ratio = 1.0; the old arm's h was in ABSOLUTE covariate units, this one scales with sd)")
