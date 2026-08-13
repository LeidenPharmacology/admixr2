## ============================================================================
## PAGE abstract fitted with the 2nd-order TAYLOR combine, for comparison
## ============================================================================
##
##   Rscript validation/page-abstract-taylor.R
##
## Same data, same unifying model as validation/page-abstract-mbma.R, but the
## objective is the node-combination the original analysis used:
##
##   NLL_study = NLL(mu_a) + 0.5 * sd_a^2 * [NLL(mu_a+h) - 2 NLL(mu_a) + NLL(mu_a-h)] / h^2
##
## i.e. a second-order expansion of E_a[NLL(a)], evaluated at three covariate
## values per study with the covariate held FIXED at each. This targets a
## different quantity from the marginal-moment objective: each node's V_pred
## contains IIV only, so the between-subject covariate spread never enters the
## predicted covariance even though it is present in V_obs.
##
## Optimised with xtol_rel passed explicitly. nloptr's default of 1e-4 makes
## BOBYQA report NLOPT_SUCCESS while still on a clear descent direction (issue
## #120), which would otherwise be indistinguishable from the method converging
## somewhere different.

suppressMessages({library(rxode2); library(nloptr)})
suppressMessages(devtools::load_all(".", quiet = TRUE))

H_FRAC <- as.numeric(Sys.getenv("ADM_TAYLOR_H", "0.5"))   # node offset, in SDs
ln_from <- function(mu, sd)
  list(meanlog = log(mu^2 / sqrt(sd^2 + mu^2)), sdlog = sqrt(log(1 + sd^2 / mu^2)))
WT_ISS <- ln_from(20.2, 9.2); WT_ALS <- ln_from(18.1, 8.5)
EV_ISS <- rxode2::et(amt = 17.5, dur = 1)
EV_ALS <- rxode2::et(amt = 15, ii = 6, addl = 4, dur = 2)
T_ISS  <- c(0.5, 1, 1.25, 1.5, 2, 3, 4, 5, 6)
T_ALS  <- c(21, 23.5)

mod_iss <- function() {
  ini({lcl <- log(0.16); lvc <- log(3.86); lvp <- log(0.19); lq <- log(0.13)
       clwt <- 0.97; vpwt <- 1.07; qwt <- 1.19
       eta.cl ~ 0.062; eta.vc ~ 0.048; prop.err <- 0.12})
  model({cl <- exp(lcl + eta.cl) * WT^clwt; vc <- exp(lvc + eta.vc)
         vp <- exp(lvp) * WT^vpwt; q <- exp(lq) * WT^qwt
         f(centr) <- WT
         d/dt(centr) <- -(cl + q)/vc*centr + q/vp*peri
         d/dt(peri)  <-        q /vc*centr - q/vp*peri
         cp <- centr/vc; cp ~ prop(prop.err)})
}
mod_als <- function() {
  ini({lcl <- log(2.99); lv <- log(9.55)
       eta.cl ~ 0.0223; eta.v ~ 0.0134; prop.err <- 0.119})
  model({cl <- exp(lcl + eta.cl) * (WT/20); v <- exp(lv + eta.v) * (WT/20)
         f(centr) <- WT
         d/dt(centr) <- -cl/v*centr; cp <- centr/v; cp ~ prop(prop.err)})
}
mod_unify <- function() {
  ini({lcl <- log(0.30); lvc <- log(2.50); lvp <- log(0.40); lq <- log(0.25)
       clwt <- 0.60; vpwt <- 0.70; qwt <- 0.80
       eta.cl ~ 0.10; eta.vc ~ 0.10; prop.err <- 0.20})
  model({cl <- exp(lcl + eta.cl) * WT^clwt; vc <- exp(lvc + eta.vc)
         vp <- exp(lvp) * WT^vpwt; q <- exp(lq) * WT^qwt
         f(centr) <- WT
         d/dt(centr) <- -(cl + q)/vc*centr + q/vp*peri
         d/dt(peri)  <-        q /vc*centr - q/vp*peri
         cp <- centr/vc; cp ~ prop(prop.err)})
}

## ---- the SAME generated data the marginal fit uses -------------------------
cat("generating study data from the published models ...\n")
gen <- datagen(list(
  Ayuthaya = list(model = mod_iss, times = T_ISS, ev = EV_ISS, n = 14L,
                  cov_dist = list(WT = WT_ISS)),
  Alsultan = list(model = mod_als, times = T_ALS, ev = EV_ALS, n = 72L,
                  cov_dist = list(WT = WT_ALS))),
  control = datagenControl(n_sim = 40000L, seed = 12345L))

## ---- three covariate nodes per study, covariate FIXED at each --------------
## Nodes on the natural weight scale at mu and mu +/- h*sd, matching the
## original. cov_dist is REMOVED: each node is an ordinary fixed-covariate
## study, which is what makes its V_pred contain IIV only.
POP <- list(Ayuthaya = list(mu = 20.2, sd = 9.2),
            Alsultan = list(mu = 18.1, sd = 8.5))
node_studies <- list(); node_meta <- list()
for (nm in names(gen)) {
  p <- POP[[nm]]; h <- H_FRAC * p$sd
  for (k in seq_along(c(0, 1, -1))) {
    off <- c(0, 1, -1)[k]
    s <- gen[[nm]]; s$cov_dist <- NULL
    s$cov <- list(WT = p$mu + off * h)
    node_studies[[paste0(nm, "_", c("0", "p", "m")[k])]] <- s
  }
  node_meta[[nm]] <- list(h = h, sd = p$sd)
}
cat("node weights:\n")
for (nm in names(node_studies))
  cat(sprintf("  %-14s WT = %6.2f\n", nm, node_studies[[nm]]$cov$WT))

ui   <- suppressMessages(rxode2::rxode2(mod_unify))
ovar <- admixr2:::.admOutputVar(ui)
ctl   <- adghControl(studies = node_studies, grad = "none", n_nodes = 7L,
                     print = 0L, covMethod = "none")
pinfo <- admixr2:::.admDriverPinfo(ui, ctl)
u     <- admixr2:::.admDriverUnits(node_studies, ui, ovar)
stu   <- u$studies
rx    <- admixr2:::.admLoadModel(ui)
ov    <- admixr2:::.admBuildOptVec(pinfo)
grid  <- admixr2:::.adghNodeGrid(7L, pinfo$n_eta)

## per-node NLL, one node-study at a time so the three can be combined
nll_node <- function(p, key)
  admixr2:::.adghNLL(p, pinfo, stu[key], rx, ovar, grid, 1L)

n_eval <- 0L
obj <- function(p) {
  n_eval <<- n_eval + 1L
  tot <- 0
  for (nm in names(node_meta)) {
    h  <- node_meta[[nm]]$h; sd <- node_meta[[nm]]$sd
    f0 <- nll_node(p, paste0(nm, "_0"))
    fp <- nll_node(p, paste0(nm, "_p"))
    fm <- nll_node(p, paste0(nm, "_m"))
    if (!all(is.finite(c(f0, fp, fm)))) return(1e12)
    tot <- tot + f0 + 0.5 * sd^2 * (fp - 2*f0 + fm) / h^2
  }
  if (!is.finite(tot)) 1e12 else tot
}

cat(sprintf("\nobjective at start: %.4f\n", obj(ov$p0)))
cat("optimising (BOBYQA, xtol_rel passed explicitly -- see issue #120) ...\n")
t0 <- Sys.time()
o <- nloptr(x0 = ov$p0, eval_f = obj,
            lb = rep(-Inf, length(ov$p0)), ub = rep(Inf, length(ov$p0)),
            opts = list(algorithm = "NLOPT_LN_BOBYQA", maxeval = 4000L,
                        ftol_rel = 1e-12, xtol_rel = 1e-8))
el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("  %.0f s, %d objective calls, %s\n", el, n_eval, o$message))
cat(sprintf("  objective %.4f\n", o$objective))

pars <- admixr2:::.admUnpack(o$solution, pinfo)
est  <- c(pars$struct, prop.err = sqrt(unname(pars$sigma_var[[1]])),
          om.cl = sqrt(pars$omega[1,1]), om.vc = sqrt(pars$omega[2,2]))
cat("\nTaylor-combine estimates:\n"); print(round(est, 4))

## for direct comparison, the marginal-objective fit on the same data
cat("\nmarginal-moment fit on the SAME data (adgh) ...\n")
t1 <- Sys.time()
fitM <- suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
  mod_unify, admData(), est = "adgh",
  control = adghControl(studies = gen, n_nodes = 7L, print = 0L,
                        covMethod = "none", maxeval = 3000L))))
cat(sprintf("  %.0f s, objective %.4f\n",
            as.numeric(difftime(Sys.time(), t1, units = "secs")), fitM$objective))
pm <- setNames(fitM$parFixedDf$Estimate, rownames(fitM$parFixedDf))
estM <- c(pm[c("lcl","lvc","lvp","lq","clwt","vpwt","qwt")],
          prop.err = unname(pm[["prop.err"]]),
          om.cl = sqrt(fitM$omega[1,1]), om.vc = sqrt(fitM$omega[2,2]))

cat("\n", strrep("=", 62), "\n", sep = "")
cat(sprintf("%-10s %12s %12s %12s\n", "param", "Taylor", "marginal", "difference"))
cat(strrep("-", 62), "\n")
for (nm in names(estM))
  cat(sprintf("%-10s %12.4f %12.4f %+12.4f\n", nm, est[[nm]], estM[[nm]],
              est[[nm]] - estM[[nm]]))
cat(sprintf("\nomega(cl): Taylor %.4f vs marginal %.4f  (%+.0f%%)\n",
            est[["om.cl"]], estM[["om.cl"]],
            100*(est[["om.cl"]]/estM[["om.cl"]] - 1)))
cat(sprintf("omega(vc): Taylor %.4f vs marginal %.4f  (%+.0f%%)\n",
            est[["om.vc"]], estM[["om.vc"]],
            100*(est[["om.vc"]]/estM[["om.vc"]] - 1)))
