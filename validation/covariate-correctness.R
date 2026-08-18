## Does admixr2 get the three pieces right?
##
##   L  the LIKELIHOOD itself, against the closed-form MVN log-likelihood
##   M  MARGINALISATION, against a large per-subject reference
##   C  CONDITIONING, against the subjects actually in each stratum
##   R  RECOVERY, from aggregate data generated at known parameters
##
## Everything is checked against something computed independently of admixr2 --
## a formula, or a brute-force sample. Nothing is pinned against admixr2's own
## previous output.
##
## Run:  Rscript validation/covariate-correctness.R
suppressMessages({library(rxode2); devtools::load_all(".", quiet = TRUE)})
say <- function(...) { cat(..., "\n"); utils::flush.console() }
ok  <- function(lbl, pass, detail = "")
  cat(sprintf("  %-52s %s   %s\n", lbl, if (pass) "PASS" else "**FAIL**", detail))

TT <- c(0.5, 1, 2, 4, 8, 12); DOSE <- 100
mod <- function() {
  ini({lcl <- log(3.0); lv <- log(20); bwt <- 0.70; bcr <- 0.55
       eta.cl ~ 0.060; prop.err <- 0.12})
  model({cl <- exp(lcl + eta.cl) * (WT/70)^bwt * (CRCL/90)^bcr
         v  <- exp(lv); d/dt(centr) <- -cl/v*centr; cp <- centr/v
         cp ~ prop(prop.err)})
}
ui   <- suppressMessages(rxode2::rxode2(mod))
ovar <- admixr2:::.admOutputVar(ui)
rxM  <- admixr2:::.admLoadModel(ui)
pop  <- covDist(WT = c(mean = 78, sd = 18), CRCL = c(mean = 95, sd = 26),
                cor = 0.45, dist = "lnorm")

## admixr2's predicted moments for one study spec
adm_moments <- function(study, nodes = 9L) {
  ctl   <- adghControl(studies = list(s = study), grad = "none", print = 0L,
                       covMethod = "none", n_nodes = nodes, cov_nodes = nodes)
  pinfo <- admixr2:::.admDriverPinfo(ui, ctl)
  stu   <- admixr2:::.admCheckCovariates(
             ui, pinfo, admixr2:::.admDriverUnits(list(s = study), ui, ovar)$studies,
             "none", "adgh")
  grid  <- admixr2:::.adghNodeGrid(nodes, pinfo$n_eta)
  pars  <- admixr2:::.admUnpack(admixr2:::.admBuildOptVec(pinfo)$p0, pinfo)
  admixr2:::.adghMoments(pars, pinfo, stu[[1L]], rxM, ovar, grid, 1L)
}
## brute force: solve for n real subjects drawn from cd, add the residual
ref_moments <- function(cd, n = 200000L) {
  cr <- admixr2:::.admCovRowsFor(cd, n, 1L)
  z  <- admixr2:::.admMakeZ(n, list(n_eta = 1L), 1L, "sobol")[[1L]]
  cl <- exp(log(3.0) + sqrt(0.060)*z[, 1L]) * (cr[,"WT"]/70)^0.70 * (cr[,"CRCL"]/90)^0.55
  f  <- DOSE/20 * exp(outer(-cl/20, TT))
  mu <- colMeans(f); V <- crossprod(sweep(f, 2L, mu))/n
  diag(V) <- diag(V) + 0.12^2 * (mu^2 + diag(V))     # law of total variance
  list(E = mu, V = V)
}
rel <- function(a, b) max(abs(a - b)/pmax(abs(b), 1e-12))
relV <- function(a, b) max(abs(a - b))/max(abs(b))

say("=== L. the likelihood itself ===")
set.seed(1)
nt <- 4L; n <- 137
Eo <- c(4.1, 3.2, 2.0, 0.9); Vo <- crossprod(matrix(rnorm(nt*nt), nt))/nt + diag(nt)*0.05
mu <- Eo * 1.04; Vp <- Vo * 1.11 + diag(nt) * 0.01
got  <- admixr2:::nll_cov_cpp(Eo, Vo, mu, Vp, n)
ch   <- chol(Vp); iv <- chol2inv(ch); r <- Eo - mu
want <- n * (2*sum(log(diag(ch))) + sum(iv*Vo) + as.numeric(t(r) %*% iv %*% r))
ok("-2LL == n(log|V| + tr(V^-1 S) + r'V^-1 r)", isTRUE(all.equal(got, want)),
   sprintf("%.10f vs %.10f", got, want))
## and it must be minimised AT the observation
base <- admixr2:::nll_cov_cpp(Eo, Vo, Eo, Vo, n); worse <- 0L
for (i in 1:300) { pe <- Eo*(1+0.02*rnorm(nt))
  Pv <- Vo + 0.02*max(abs(Vo))*crossprod(matrix(rnorm(nt*nt), nt))
  if (admixr2:::nll_cov_cpp(Eo, Vo, pe, Pv, n) > base) worse <- worse + 1L }
ok("minimised at mu = E_obs, V = V_obs", worse == 300L, sprintf("%d/300", worse))

say("\n=== M. marginalisation, vs 200k per-subject solves ===")
st <- list(E = rep(1, length(TT)), V = diag(length(TT)) + 0.01, n = 100L,
           times = TT, ev = rxode2::et(amt = DOSE), cov_dist = pop)
rf <- ref_moments(pop)
for (nd in c(5L, 7L, 11L)) { m <- adm_moments(st, nd)
  ok(sprintf("marginal moments, %2d nodes", nd),
     rel(m$E, rf$E) < 5e-3 && relV(m$V, rf$V) < 5e-3,
     sprintf("relE %.2e  relV %.2e", rel(m$E, rf$E), relV(m$V, rf$V))) }

say("\n=== C. conditioning, vs the subjects actually in each stratum ===")
str <- covStrata(pop, stratify = "WT", n_nodes = 4L, n = 400)
for (k in seq_along(str)) {
  sk <- list(E = rep(1, length(TT)), V = diag(length(TT)) + 0.01,
             n = 100L, times = TT, ev = rxode2::et(amt = DOSE),
             cov_dist = str[[k]]$cov_dist)
  mk <- adm_moments(sk, 9L)
  rk <- ref_moments(str[[k]]$cov_dist, 100000L)
  ok(sprintf("stratum %d (WT %6.1f)", k, str[[k]]$cov$WT),
     rel(mk$E, rk$E) < 8e-3 && relV(mk$V, rk$V) < 2e-2,
     sprintf("relE %.2e  relV %.2e", rel(mk$E, rk$E), relV(mk$V, rk$V)))
}
## law of total expectation across the strata
w  <- vapply(str, function(s) s$weight, 0)
Ek <- vapply(str, function(s) adm_moments(list(E=rep(1,length(TT)),
        V=diag(length(TT))+0.01, n=100L, times=TT, ev=rxode2::et(amt=DOSE),
        cov_dist=s$cov_dist), 9L)$E, numeric(length(TT)))
ok("sum_k w_k E_k == the marginal E", rel(as.numeric(Ek %*% w), rf$E) < 8e-3,
   sprintf("relE %.2e", rel(as.numeric(Ek %*% w), rf$E)))

say("\n=== R. recovery from aggregate data at known parameters ===")
fitm <- function() {
  ini({lcl <- log(2.0); lv <- log(15); bwt <- 0.20; bcr <- 0.20
       eta.cl ~ 0.03; prop.err <- 0.30})
  model({cl <- exp(lcl + eta.cl) * (WT/70)^bwt * (CRCL/90)^bcr
         v  <- exp(lv); d/dt(centr) <- -cl/v*centr; cp <- centr/v
         cp ~ prop(prop.err)})
}
rec <- function(studies, lbl, tol = 0.02) {
  f <- suppressMessages(nlmixr2est::nlmixr2(fitm, admData(), est = "adgh",
        control = adghControl(studies = studies, print = 0L, covMethod = "none")))
  e <- stats::setNames(f$parFixedDf$Estimate, rownames(f$parFixedDf))
  d <- max(abs(c(e[["bwt"]] - 0.70, e[["bcr"]] - 0.55)))
  ok(lbl, d < tol, sprintf("bwt %.4f  bcr %.4f  omega %.4f",
     e[["bwt"]], e[["bcr"]], sqrt(f$omega[1,1])))
}
popB <- covDist(WT = c(mean = 62, sd = 12), CRCL = c(mean = 55, sd = 15),
                cor = 0.45, dist = "lnorm")
g <- function(stratify, pops, ns)
  datagen(stats::setNames(lapply(seq_along(pops), function(i)
    c(list(model = mod, n = ns[i], times = TT, ev = rxode2::et(amt = DOSE),
           cov_dist = pops[[i]]), if (stratify) list(stratify = TRUE))),
    paste0("s", seq_along(pops))),
    control = datagenControl(method = "mc", n_sim = 20000L))
rec(g(TRUE,  list(pop, popB), c(200L, 200L)), "two studies, both stratified")
rec(g(FALSE, list(pop, popB), c(200L, 200L)), "two studies, both marginal")
rec(g(TRUE,  list(pop), 200L),                "one study, stratified on both")
say("\n(a source model with BOTH covariates stratifies on both; the marginal")
say(" arm gives the estimator only two pooled summaries and no contrast)")
