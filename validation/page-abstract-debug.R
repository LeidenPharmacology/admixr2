## ============================================================================
## PAGE abstract reproduction -- DEBUGGING PASS
## ============================================================================
##
##   Rscript validation/page-abstract-debug.R
##
## Checks every assumption the reproduction rests on, BEFORE any joint fit.
## Output is KEY<TAB>value so each claim can be read off rather than eyeballed.
##
## The frame this is testing, and its one important correction:
##
##   admixr2 generates each study's aggregate data from ITS OWN published model
##   and design, then fits ONE unifying model across all of them. Because the
##   studies come from DIFFERENT published models, the joint fit is a consensus
##   and its optimum is NOT at any single study's parameters. So "recovers the
##   truth" is only a meaningful criterion for a SINGLE-study, same-model
##   simulation-estimation check (block E below). The joint fit has to be judged
##   by whether it describes both studies' data (block F), not by parameter
##   recovery against numbers that no longer define a truth.

suppressMessages(devtools::load_all(".", quiet = TRUE))
kv <- function(k, v) cat(sprintf("%s\t%s\n", k, paste(v, collapse = ",")))

ln_from <- function(mu, sd)
  list(meanlog = log(mu^2 / sqrt(sd^2 + mu^2)), sdlog = sqrt(log(1 + sd^2 / mu^2)))
WT_ISS <- ln_from(20.2, 9.2); WT_ALS <- ln_from(18.1, 8.5)
kv("wt_iss_meanlog_sdlog", sprintf("%.4f", unlist(WT_ISS)))
kv("wt_iss_implied_mean_sd",
   sprintf("%.3f", c(exp(WT_ISS$meanlog + WT_ISS$sdlog^2/2),
                     sqrt((exp(WT_ISS$sdlog^2) - 1) * exp(2*WT_ISS$meanlog + WT_ISS$sdlog^2)))))
kv("wt_als_implied_mean_sd",
   sprintf("%.3f", c(exp(WT_ALS$meanlog + WT_ALS$sdlog^2/2),
                     sqrt((exp(WT_ALS$sdlog^2) - 1) * exp(2*WT_ALS$meanlog + WT_ALS$sdlog^2)))))

mod_iss <- function() {
  ini({lcl <- log(0.16); lvc <- log(3.86); lvp <- log(0.19); lq <- log(0.13)
       clwt <- 0.97; vpwt <- 1.07; qwt <- 1.19
       eta.cl ~ 0.062; eta.vc ~ 0.048
       prop.err <- 0.12})
  model({cl <- exp(lcl + eta.cl) * WT^clwt
         vc <- exp(lvc + eta.vc)
         vp <- exp(lvp) * WT^vpwt
         q  <- exp(lq)  * WT^qwt
         f(centr)    <- WT
         d/dt(centr) <- -(cl + q) / vc * centr + q / vp * peri
         d/dt(peri)  <-        q  / vc * centr - q / vp * peri
         cp <- centr / vc
         cp ~ prop(prop.err)})
}
mod_als <- function() {
  ini({lcl <- log(2.99); lv <- log(9.55)
       eta.cl ~ 0.0223; eta.v ~ 0.0134
       prop.err <- 0.119})
  model({cl <- exp(lcl + eta.cl) * (WT / 20)
         v  <- exp(lv  + eta.v)  * (WT / 20)
         f(centr)    <- WT
         d/dt(centr) <- -cl / v * centr
         cp <- centr / v
         cp ~ prop(prop.err)})
}
EV_ISS <- rxode2::et(amt = 17.5, dur = 1)
EV_ALS <- rxode2::et(amt = 15, ii = 6, addl = 4, dur = 2)
T_ISS  <- c(0.5, 1, 1.25, 1.5, 2, 3, 4, 5, 6)
T_ALS  <- c(21, 23.5)

## ---- A. does f(centr) <- WT actually deliver a mg/kg dose? ------------------
## With dose proportional to WT and Cl proportional to WT^0.97, steady exposure
## (AUC = dose/Cl) should scale as WT^0.03 -- i.e. be nearly flat in weight.
## If f() were ignored, AUC would fall as WT^-0.97.
ui  <- suppressMessages(rxode2::rxode2(mod_iss))
rxM <- admixr2:::.admLoadModel(ui)
pin <- admixr2:::.admParseIniDf(ui$iniDf, ui); pin$cov_map <- admixr2:::.admCovMap(ui)
pin$nDisplayProgress <- .Machine$integer.max
ov  <- admixr2:::.admBuildOptVec(pin); prs <- admixr2:::.admUnpack(ov$p0, pin)
# Zero the allometric exponents: PK then does not depend on weight at all, so
# the ONLY route from WT to the prediction is the dose. cp must scale exactly
# with WT if f() works, and not at all if f() is ignored. Comparing AUC ratios
# against WT^0.03 (an earlier version of this check) could not distinguish
# "correct" from "subtly wrong", because summed concentrations over 6 h are not
# AUC-infinity and Vc is not weight-scaled.
for (.e in c("clwt", "vpwt", "qwt")) prs$struct[[.e]] <- 0
wts <- c(10, 20, 40)
pm  <- admixr2:::.admMakeParamsList(length(wts), pin, 1L)[[1L]]
for (nm in pin$struct_names) pm[, nm] <- prs$struct[[nm]]
pm  <- admixr2:::.admCovCols(pm, rxM$params, NULL,
        matrix(wts, ncol = 1L, dimnames = list(NULL, "WT")))
st_probe <- list(ev_full = admixr2:::.admBuildEvFull(
  list(u = list(ev = EV_ISS, times = T_ISS)))[[1L]]$ev_full, times = T_ISS)
out <- rxode2::rxSolve(rxM, params = as.data.frame(pm),
                       events = st_probe$ev_full, cores = 1L,
                       nDisplayProgress = .Machine$integer.max)
cp <- matrix(out[["cp"]][out[["time"]] %in% T_ISS], nrow = length(wts), byrow = TRUE)
auc <- rowSums(cp)
kv("A_weights", wts)
kv("A_auc", sprintf("%.4f", auc))
kv("A_cp_ratio_40_over_10", sprintf("%.6f", cp[3, ] / cp[1, ]))
kv("A_expected_if_f_works", "4.000000 (exactly, at every time)")
kv("A_expected_if_f_IGNORED", "1.000000")
kv("A_f_scales_dose", as.integer(max(abs(cp[3, ] / cp[1, ] - 4)) < 1e-6))

## ---- B. does the covariate vary PER ROW in a solve? ------------------------
kv("B_cp_varies_across_weights", as.integer(length(unique(round(cp[, 1], 8))) == 3L))

## ---- C. does marginalising over weight change E and V at all? --------------
## If it does not, the covariate is not reaching the aggregation.
g_marg <- datagen(list(s = list(model = mod_iss, times = T_ISS, ev = EV_ISS, n = 14L,
                                cov_dist = list(WT = WT_ISS))),
                  control = datagenControl(n_sim = 20000L, seed = 1L))$s
## the control is the covariate held at the population MEAN, not absent: a model
## that reads WT with nothing supplied cannot be solved at all.
.wt_mean <- exp(WT_ISS$meanlog + WT_ISS$sdlog^2 / 2)
g_fix  <- datagen(list(s = list(model = mod_iss, times = T_ISS, ev = EV_ISS, n = 14L,
                                cov = list(WT = .wt_mean))),
                  control = datagenControl(n_sim = 20000L, seed = 1L))$s
kv("C_marg_E", sprintf("%.4f", g_marg$E))
kv("C_fixedWT_E", sprintf("%.4f", g_fix$E))
kv("C_E_differs", as.integer(max(abs(g_marg$E / g_fix$E - 1)) > 0.01))
kv("C_relsd_marg", sprintf("%.4f", sqrt(diag(g_marg$V)) / g_marg$E))
kv("C_relsd_fixed", sprintf("%.4f", sqrt(diag(g_fix$V)) / g_fix$E))
kv("C_marg_V_larger", as.integer(all(diag(g_marg$V) > diag(g_fix$V))))
kv("C_carries_cov_dist", as.integer(!is.null(g_marg$cov_dist)))
kv("C_carries_cov", as.integer(!is.null(g_marg$cov)))

## ---- D. which covariate path does the unifying model take? ----------------
mod_unify <- mod_iss     # same structure; started elsewhere in the real script
uiU <- suppressMessages(rxode2::rxode2(mod_unify))
kv("D_allCovs", uiU$allCovs)
kv("D_WT_occurrences", admixr2:::.admNameOccurrence(uiU, "WT")[["WT"]])
ctlU  <- admControl(studies = list(s = g_marg), grad = "sens", n_sim = 500L,
                    print = 0L, covMethod = "none")
pinU  <- admixr2:::.admDriverPinfo(uiU, ctlU)
uU    <- admixr2:::.admDriverUnits(list(s = g_marg), uiU, admixr2:::.admOutputVar(uiU))
stU   <- admixr2:::.admCheckCovariates(uiU, pinU, uU$studies, "sens")
kv("D_path", stU[[1L]]$.adm_cov_path)
kv("D_unit_has_cov_dist", as.integer(!is.null(stU[[1L]][["cov_dist"]])))

## ---- E. same model, same population: does the fit see its own data? -------
## The ONLY place "recovers the truth" is a meaningful criterion.
zl <- admixr2:::.admMakeZ(4000L, pinU, 1L, "sobol")
pl <- admixr2:::.admMakeParamsList(4000L, pinU, 1L)
rxU <- admixr2:::.admLoadModel(uiU)
ovU <- admixr2:::.admBuildOptVec(pinU)
f   <- function(pp) admixr2:::.admNLL(pp, pinU, stU, zl, rxU,
                                      admixr2:::.admOutputVar(uiU), pl, 1L)
kv("E_nll_at_truth", sprintf("%.4f", f(ovU$p0)))
k_clwt <- match("clwt", pinU$struct_names)
prof <- vapply(c(0.77, 0.87, 0.97, 1.07, 1.17), function(x) {
  q <- ovU$p0; q[k_clwt] <- x; f(q) }, numeric(1))
kv("E_clwt_grid", c(0.77, 0.87, 0.97, 1.07, 1.17))
kv("E_clwt_profile_minus_min", sprintf("%.3f", prof - min(prof)))
kv("E_clwt_argmin_is_truth", as.integer(which.min(prof) == 3L))

sm <- tryCatch(admixr2:::.admLoadSensModel(uiU), error = function(e) NULL)
kv("E_sens_model_available", as.integer(!is.null(sm)))
ga <- admixr2:::.admGrad(ovU$p0, pinU, stU, zl, rxU,
                         admixr2:::.admOutputVar(uiU), pl, 1L, 1e-4, sensModel = sm)
kv("E_par_names", c(pinU$struct_names, pinU$sigma_names, "omega"))
kv("E_grad_at_truth", sprintf("%.3f", ga))
h <- 1e-5
gf <- vapply(seq_along(ovU$p0), function(k) {
  a <- b <- ovU$p0; a[k] <- a[k] + h; b[k] <- b[k] - h; (f(a) - f(b)) / (2 * h)
}, numeric(1))
kv("E_grad_centralFD", sprintf("%.3f", gf))
kv("E_grad_max_relerr", sprintf("%.2e", max(abs(ga - gf) / pmax(abs(gf), 1e-6))))

## ---- F. the joint frame: is the consensus optimum at either study's truth? --
## Alsultan's data come from a 1-cmt model with its OWN parameters, so the
## unifying 2-cmt fit cannot be at Issaranggoon's values. Quantify the pull.
g_als <- datagen(list(s = list(model = mod_als, times = T_ALS, ev = EV_ALS, n = 72L,
                               cov_dist = list(WT = WT_ALS))),
                 control = datagenControl(n_sim = 20000L, seed = 2L))$s
both  <- list(Issaranggoon = g_marg, Alsultan = g_als)
uB    <- admixr2:::.admDriverUnits(both, uiU, admixr2:::.admOutputVar(uiU))
stB   <- admixr2:::.admCheckCovariates(uiU, pinU, uB$studies, "sens")
zlB   <- admixr2:::.admMakeZ(4000L, pinU, 2L, "sobol")
plB   <- admixr2:::.admMakeParamsList(4000L, pinU, 2L)
fB    <- function(pp) admixr2:::.admNLL(pp, pinU, stB, zlB, rxU,
                                        admixr2:::.admOutputVar(uiU), plB, 1L)
kv("F_joint_nll_at_iss_truth", sprintf("%.4f", fB(ovU$p0)))
profB <- vapply(c(0.77, 0.87, 0.97, 1.07, 1.17), function(x) {
  q <- ovU$p0; q[k_clwt] <- x; fB(q) }, numeric(1))
kv("F_joint_clwt_profile_minus_min", sprintf("%.3f", profB - min(profB)))
kv("F_joint_clwt_argmin_is_iss_truth", as.integer(which.min(profB) == 3L))
kv("F_alsultan_alone_nll_at_iss_truth",
   sprintf("%.4f", fB(ovU$p0) - f(ovU$p0)))
