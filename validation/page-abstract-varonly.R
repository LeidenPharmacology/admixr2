## PAGE case study with VARIANCES ONLY -- no between-timepoint covariances.
##
## This is what a publication actually gives you: mean +- SD at each sampling
## time, never the covariance between the 2 h and the 4 h concentration. So
## `method = "var"` is the realistic aggregate-data case and the full-covariance
## fit is the optimistic one. Same generated studies both ways -- the var fit
## gets diag(V) of the very same data -- so the only difference is how much of it
## the likelihood reads.
##
## THE RESULT: variances alone cannot separate residual error from IIV. IIV
## correlates concentrations ACROSS times; residual error does not, so the
## off-diagonal is the only thing that distinguishes them. Seeing only the
## diagonal, sigma^2 and omega^2 both simply inflate it, and the optimiser drives
## sigma to zero and absorbs it into omega. Measured below as a 337-fold
## difference in how hard each objective resists moving prop.err.
##
## Two traps this script exists to avoid, both of which produced wrong answers
## first time round:
##   * the DRIVER's var fit stops early (nloptr xtol_rel = 1e-4, issue #120) at
##     1008.35 when 1006.58 is reachable -- so it refits with a tolerance that
##     can fire, and reports both.
##   * LBFGS started from the cov solution FAILS on the var objective
##     (NLOPT_FAILURE, max|g| ~ 782, no movement). A run that reports that point
##     is reporting an optimiser failure, not an estimate.
##
## Run:  Rscript validation/page-abstract-varonly.R
suppressMessages({library(rxode2); library(nloptr)})
suppressMessages(devtools::load_all(".", quiet = TRUE))
kv <- function(k, v) cat(sprintf("%s\t%s\n", k, paste(v, collapse = ",")))

ln_from <- function(mu, sd)
  list(meanlog = log(mu^2 / sqrt(sd^2 + mu^2)), sdlog = sqrt(log(1 + sd^2 / mu^2)))
WT_ISS <- ln_from(20.2, 9.2); WT_ALS <- ln_from(18.1, 8.5)
EV_ISS <- rxode2::et(amt = 17.5, dur = 1)
EV_ALS <- rxode2::et(amt = 15, ii = 6, addl = 4, dur = 2)
T_ISS  <- c(0.5, 1, 1.25, 1.5, 2, 3, 4, 5, 6); T_ALS <- c(21, 23.5)
mod_iss <- function() {
  ini({lcl <- log(0.16); lvc <- log(3.86); lvp <- log(0.19); lq <- log(0.13)
       clwt <- 0.97; vpwt <- 1.07; qwt <- 1.19
       eta.cl ~ 0.062; eta.vc ~ 0.048; prop.err <- 0.12})
  model({cl <- exp(lcl + eta.cl) * WT^clwt; vc <- exp(lvc + eta.vc)
         vp <- exp(lvp) * WT^vpwt; q <- exp(lq) * WT^qwt; f(centr) <- WT
         d/dt(centr) <- -(cl + q)/vc*centr + q/vp*peri
         d/dt(peri)  <-        q /vc*centr - q/vp*peri
         cp <- centr/vc; cp ~ prop(prop.err)})
}
mod_als <- function() {
  ini({lcl <- log(2.99); lv <- log(9.55); eta.cl ~ 0.0223; eta.v ~ 0.0134
       prop.err <- 0.119})
  model({cl <- exp(lcl + eta.cl)*(WT/20); v <- exp(lv + eta.v)*(WT/20)
         f(centr) <- WT
         d/dt(centr) <- -cl/v*centr; cp <- centr/v; cp ~ prop(prop.err)})
}
gen <- datagen(list(
  Ayuthaya = list(model = mod_iss, times = T_ISS, ev = EV_ISS, n = 14L,
                  cov_dist = list(WT = WT_ISS)),
  Alsultan = list(model = mod_als, times = T_ALS, ev = EV_ALS, n = 72L,
                  cov_dist = list(WT = WT_ALS))),
  control = datagenControl(n_sim = 40000L, seed = 12345L))
gen_var <- lapply(gen, function(s) { s$V <- diag(as.matrix(s$V)); s })

ui   <- suppressMessages(rxode2::rxode2(mod_iss))
ovar <- admixr2:::.admOutputVar(ui)
setup <- function(studies) {
  ctl   <- adghControl(studies = studies, grad = "analytical", n_nodes = 7L,
                       print = 0L, covMethod = "none")
  pinfo <- admixr2:::.admDriverPinfo(ui, ctl)
  u     <- admixr2:::.admDriverUnits(studies, ui, ovar)$studies
  u     <- admixr2:::.admCheckCovariates(ui, pinfo, u, "analytical", "adgh")
  list(pinfo = pinfo, u = u, ov = admixr2:::.admBuildOptVec(pinfo),
       rx = admixr2:::.admLoadModel(ui),
       sm = tryCatch(admixr2:::.admLoadSensModel(ui), error = function(e) NULL),
       grid = admixr2:::.adghNodeGrid(7L, admixr2:::.admDriverPinfo(ui, ctl)$n_eta))
}
sc <- setup(gen); sv <- setup(gen_var)
kv("method_cov_fit", paste(vapply(sc$u, function(s) s$method, character(1)), collapse = "/"))
kv("method_var_fit", paste(vapply(sv$u, function(s) s$method, character(1)), collapse = "/"))
kv("v_diag_present_var", paste(vapply(sv$u, function(s) !is.null(s$v_diag), logical(1)), collapse = "/"))

mk <- function(th, prop, om1, om2) {
  p <- sc$ov$p0; ns <- length(sc$pinfo$struct_names)
  for (i in seq_len(ns)) p[i] <- th[[sc$pinfo$struct_names[i]]]
  p[ns + 1] <- log(prop^2); p[ns + 2] <- log(om1^2); p[ns + 3] <- log(om2^2); p
}
COV <- mk(list(lcl=-1.6642, lvc=1.0218, lvp=-1.3915, lq=-3.8028, clwt=0.8742,
               vpwt=0.9514, qwt=1.8714), 0.1084, 0.1514, 0.4639)
VAR <- mk(list(lcl=-2.5255, lvc=1.5436, lvp=-5.9163, lq=1.9272, clwt=1.1951,
               vpwt=2.4077, qwt=-0.1432), 0.0164, 0.0260, 0.3862)

fv <- function(p) admixr2:::.adghNLL(p, sv$pinfo, sv$u, sv$rx, ovar, sv$grid, 1L)
gv <- function(p) admixr2:::.adghGrad(p, sv$pinfo, sv$u, sv$sm, sv$rx, ovar, sv$grid, 1L, 1e-4)
kv("var_obj_at_VARsol", sprintf("%.4f", fv(VAR)))
kv("var_obj_at_COVsol", sprintf("%.4f", fv(COV)))
kv("var_grad_maxabs_at_VARsol", sprintf("%.4f", max(abs(gv(VAR)))))

kv("grad_finite_at_COV", all(is.finite(gv(COV))))

## Optimise the VAR objective properly from BOTH starts, with a tolerance that
## can actually fire. nloptr's default xtol_rel = 1e-4 stops early (issue #120),
## which is the leading explanation for a gradient of 3.6 at the reported
## solution.
opt <- function(x0, tag) {
  o <- nloptr(x0 = x0, eval_f = fv, eval_grad_f = gv,
              lb = rep(-Inf, length(x0)), ub = rep(Inf, length(x0)),
              opts = list(algorithm = "NLOPT_LD_LBFGS", maxeval = 6000L,
                          ftol_rel = 1e-14, xtol_rel = 1e-12))
  pr <- admixr2:::.admUnpack(o$solution, sv$pinfo)
  kv(paste0("opt_", tag),
     sprintf("obj=%.4f  evals=%d  max|g|=%.4f  %s", o$objective, o$iterations,
             max(abs(gv(o$solution))), o$message))
  kv(paste0("opt_", tag, "_struct"),
     sprintf("%s=%.4f", sv$pinfo$struct_names, pr$struct))
  kv(paste0("opt_", tag, "_om_sig"),
     sprintf("om=%.4f/%.4f prop=%.4f", sqrt(diag(pr$omega))[1],
             sqrt(diag(pr$omega))[2], sqrt(pr$sigma_var)))
  o
}
o1 <- opt(COV, "from_COV")
o2 <- opt(VAR, "from_VAR")
kv("best_var_objective", sprintf("%.4f", min(o1$objective, o2$objective)))
kv("reported_var_objective", sprintf("%.4f", fv(VAR)))
kv("reported_was_stuck",
   as.integer(min(o1$objective, o2$objective) < fv(VAR) - 1e-3))

## ---- WHY is the var fit degenerate? -----------------------------------------
## Hypothesis: the off-diagonal is what SEPARATES residual error from IIV. IIV
## correlates concentrations across times; residual error does not. Seeing only
## the variances, sigma^2 and omega^2 both just inflate the diagonal, so they are
## confounded -- and the optimiser drives sigma to zero and absorbs it into omega.
##
## Test: fix prop.err at the cov-fit value, re-optimise everything else, and see
## how much each objective degrades. If the var objective barely moves while the
## cov objective moves a lot, sigma is unidentified under var only.
fc <- function(p) admixr2:::.adghNLL(p, sc$pinfo, sc$u, sc$rx, ovar, sc$grid, 1L)
gc_ <- function(p) admixr2:::.adghGrad(p, sc$pinfo, sc$u, sc$sm, sc$rx, ovar, sc$grid, 1L, 1e-4)
ns  <- length(sv$pinfo$struct_names); isig <- ns + 1L

prof_sigma <- function(f, g, x0, sig, tag) {
  free <- setdiff(seq_along(x0), isig)
  fx <- function(q) { p <- x0; p[free] <- q; p[isig] <- log(sig^2); f(p) }
  gx <- function(q) { p <- x0; p[free] <- q; p[isig] <- log(sig^2); g(p)[free] }
  o <- nloptr(x0 = x0[free], eval_f = fx, eval_grad_f = gx,
              lb = rep(-Inf, length(free)), ub = rep(Inf, length(free)),
              opts = list(algorithm = "NLOPT_LD_LBFGS", maxeval = 4000L,
                          ftol_rel = 1e-14, xtol_rel = 1e-12))
  kv(sprintf("profile_%s_sigma%.4f", tag, sig), sprintf("obj=%.4f", o$objective))
  o$objective
}
best_var <- o2$solution
v_free <- prof_sigma(fv, gv, best_var, 0.1084, "var")
v_best <- o2$objective
kv("var_penalty_for_sigma_0.108", sprintf("%.4f  (best %.4f)", v_free - v_best, v_best))

c_best <- fc(COV)
c_free <- prof_sigma(fc, gc_, COV, 0.0082, "cov")
kv("cov_penalty_for_sigma_0.0082", sprintf("%.4f  (best %.4f)", c_free - c_best, c_best))
