## Dependent aggregate patient statistics via a copula sampler.
##
## Real patient covariates are dependent -- weight with age, weight with height,
## creatinine with age -- and an R-vine copula is the general way to describe
## that: sample the vine on the uniform scale, then push each margin through its
## own quantile function. `cov_dist$joint` is the hook for it.
##
## Checked three ways, because a dependence that is silently dropped still gives
## a plausible-looking fit:
##   1. the drawn covariates reproduce the intended dependence and margins;
##   2. admixr2's predicted moments match an INDEPENDENT reference computed from
##      the same subjects by plain rxode2;
##   3. the dependence actually MATTERS -- the same margins with rho = 0 give
##      materially different aggregate moments, so a check that passed either way
##      would prove nothing.
##
## Run:  Rscript validation/covariate-dependent-copula.R
suppressMessages({library(rxode2)})
suppressMessages(devtools::load_all(".", quiet = TRUE))
kv <- function(k, v) cat(sprintf("%s\t%s\n", k, paste(v, collapse = ",")))

TT <- c(1, 2, 4, 8, 12); EV <- rxode2::et(amt = 20, dur = 0.5)
## weight and age, dependent. A Gaussian copula stands in for the vine here so
## the script has no extra dependency; the interface is identical -- any sampler
## that maps the supplied uniforms to covariate values works, RVineSim included.
MARG <- list(WT  = function(u) stats::qlnorm(u, log(70), 0.25),
             AGE = function(u) stats::qgamma(u, 9, 0.3))
gauss_copula <- function(rho) function(u) {
  z <- stats::qnorm(u)
  z2 <- cbind(z[, 1], rho * z[, 1] + sqrt(1 - rho^2) * z[, 2])
  v  <- stats::pnorm(z2)
  cbind(WT = MARG$WT(v[, 1]), AGE = MARG$AGE(v[, 2]))
}
cov_dist_for <- function(rho) list(
  WT  = list(quantile = MARG$WT),
  AGE = list(quantile = MARG$AGE),
  joint = gauss_copula(rho))

mod <- function() {
  ini({lcl <- log(3.0); lv <- log(20); clwt <- 0.75; agee <- 0.20
       eta.cl ~ 0.09; prop.err <- 0.15})
  model({cl <- exp(lcl + eta.cl) * (WT/70)^clwt * (AGE/30)^agee
         v  <- exp(lv) * (WT/70)
         f(centr) <- WT
         d/dt(centr) <- -cl / v * centr; cp <- centr / v
         cp ~ prop(prop.err)})
}

## ---- 1. do the draws carry the dependence and the margins? -----------------
for (rho in c(0.0, 0.8)) {
  cr <- admixr2:::.admCovRowsFor(cov_dist_for(rho), 20000L, 1L)
  kv(sprintf("draws_rho%.1f", rho),
     sprintf("cor(WT,AGE)=%+.3f  mean WT=%.2f (true %.2f)  mean AGE=%.2f (true %.2f)",
             stats::cor(cr[, "WT"], cr[, "AGE"]), mean(cr[, "WT"]),
             exp(log(70) + 0.25^2 / 2), mean(cr[, "AGE"]), 9 / 0.3))
}

## ---- 2. predicted moments vs an independent per-subject reference ----------
ui   <- suppressMessages(rxode2::rxode2(mod))
ovar <- admixr2:::.admOutputVar(ui)
rx   <- suppressMessages(rxode2::rxode2(ui$simulationModel))
ref_moments <- function(rho, n = 60000L) {
  set.seed(9L)
  u  <- cbind(stats::runif(n), stats::runif(n))
  cr <- gauss_copula(rho)(u)
  p  <- data.frame(lcl = log(3.0), lv = log(20), clwt = 0.75, agee = 0.20,
                   prop.err = 0.15, rxerr.cp = 0,
                   eta.cl = 0.30 * stats::rnorm(n),
                   WT = cr[, "WT"], AGE = cr[, "AGE"])
  o  <- as.data.frame(rxode2::rxSolve(rx, params = p,
         events = EV |> rxode2::et(TT), cores = 1L,
         nDisplayProgress = .Machine$integer.max))
  o  <- o[o$time %in% TT, ]
  f  <- matrix(o$cp, nrow = n, byrow = TRUE)
  mu <- colMeans(f); V <- crossprod(sweep(f, 2L, mu)) / n
  diag(V) <- diag(V) + 0.15^2 * (mu^2 + diag(V))   # law of total variance
  list(E = mu, V = V)
}
adm_moments <- function(rho, n_sim = 60000L) {
  st <- list(s = list(E = rep(1, length(TT)), V = diag(length(TT)) + 0.01,
                      n = 100L, times = TT, ev = EV, cov_dist = cov_dist_for(rho)))
  ctl   <- admControl(studies = st, grad = "sens", n_sim = n_sim, print = 0L,
                      covMethod = "none")
  pinfo <- admixr2:::.admDriverPinfo(ui, ctl)
  u     <- admixr2:::.admCheckCovariates(
             ui, pinfo, admixr2:::.admDriverUnits(st, ui, ovar)$studies, "sens")
  s     <- admixr2:::.admStudyCovRows(u[[1L]], pinfo, n_sim)
  z     <- admixr2:::.admMakeZ(n_sim, pinfo, 1L, "sobol")[[1L]]
  pars  <- admixr2:::.admUnpack(admixr2:::.admBuildOptVec(pinfo)$p0, pinfo)
  eta   <- z %*% t(pars$L); colnames(eta) <- pinfo$eta_col_names
  pm    <- admixr2:::.admMakeParamsList(n_sim, pinfo, 1L)[[1L]]
  cp    <- admixr2:::.admSimulate(admixr2:::.admLoadModel(ui), pars$struct,
                                  pinfo$sigma_names, eta, s, ovar, pm, 1L,
                                  .Machine$integer.max, NULL)
  mu <- colMeans(cp); V <- crossprod(sweep(cp, 2L, mu)) / n_sim
  arr <- admixr2:::.admUnitResidRows(pinfo, ovar, pars$sigma_var, length(TT))
  ap  <- admixr2:::.admResidApply(mu, diag(V), arr, s$times, V)
  list(E = ap$mu, V = admixr2:::.admApplyResidTail(V, ap), path = s$.adm_cov_path)
}
for (rho in c(0.0, 0.8)) {
  a <- adm_moments(rho); r <- ref_moments(rho)
  kv(sprintf("moments_rho%.1f", rho),
     sprintf("path=%s  E rel %.3e  V rel %.3e", a$path,
             max(abs(a$E - r$E) / abs(r$E)),
             max(abs(a$V - r$V)) / max(abs(r$V))))
}

## ---- 3. does the dependence matter? ----------------------------------------
r0 <- ref_moments(0.0); r8 <- ref_moments(0.8)
kv("dependence_effect",
   sprintf("E max rel %.3f  V max rel %.3f  (rho 0 vs 0.8, same margins)",
           max(abs(r8$E - r0$E) / r0$E),
           max(abs(r8$V - r0$V)) / max(abs(r0$V))))

## ---- 4. adgh must refuse rather than integrate the wrong distribution ------
kv("adgh_refuses_joint",
   tryCatch({
     st <- list(s = list(E = rep(1, length(TT)), V = diag(length(TT)), n = 100L,
                         times = TT, ev = EV, cov_dist = cov_dist_for(0.8)))
     ctl <- adghControl(studies = st, grad = "analytical", print = 0L,
                        covMethod = "none")
     pinfo <- admixr2:::.admDriverPinfo(ui, ctl)
     admixr2:::.admCheckCovariates(
       ui, pinfo, admixr2:::.admDriverUnits(st, ui, ovar)$studies,
       "analytical", "adgh")
     "ACCEPTED (wrong -- product grid assumes independence)"
   }, error = function(e) paste("refused:", substr(conditionMessage(e), 1, 60))))
