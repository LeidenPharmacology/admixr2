## Same three-source, covariate-banded MBMA as vignette_covariates_oldvsnew.R,
## but the TRUE kinetics are 2-compartment and ONE of the three published
## analysts (moderate/Ito) fit their trial with a 1-compartment model --
## genuinely misspecified, not just a different covariate choice. normal/Sato
## and mild/Khan fit (correctly) 2-compartment models, matching the pooled
## model's own structure. Compares srcWeight = "n" vs "cov" the same way,
## now in the regime the standalone 1cmt/2cmt test already showed is where
## OLD's mechanism does the most damage.
suppressMessages({ library(nlmixr2est); devtools::load_all(".", quiet = TRUE) })
rxode2::setRxThreads(2L)

set.seed(11)
TIMES <- c(0.5, 1, 2, 4, 8, 12, 24); DOSE <- 200
CL70 <- 5; V1_70 <- 30; QQ <- 8; V2 <- 60; BCRCL <- 0.6; BSEX <- 0.18
OM_CL <- 0.05; SD_ADD <- 0.08
R <- matrix(c(1, .45, .45, 1), 2, 2,
            dimnames = list(c("WT", "CRCL"), c("WT", "CRCL")))

draw_cohort <- function(n, wt_median, wt_cv, crcl_median, crcl_cv, p_male = 0.55) {
  z <- matrix(rnorm(2 * n), nrow = n, ncol = 2) %*% chol(R)
  data.frame(WT = qlnorm(pnorm(z[, 1]), meanlog = log(wt_median), sdlog = wt_cv),
             CRCL = qlnorm(pnorm(z[, 2]), meanlog = log(crcl_median), sdlog = crcl_cv),
             SEX = rbinom(n, size = 1L, prob = p_male))
}
## TRUE 2-cmt IV-bolus biexponential concentration
conc2 <- function(t, dose, cl, v1, q, v2) {
  k10 <- cl/v1; k12 <- q/v1; k21 <- q/v2
  s <- k10 + k12 + k21; dis <- sqrt(max(s^2 - 4*k10*k21, 0))
  a <- (s + dis)/2; b <- (s - dis)/2
  dose/v1 * ((a - k21)/(a - b) * exp(-a*t) + (k21 - b)/(a - b) * exp(-b*t))
}
patient_cl <- function(cohort, cl_70kg, bcrcl, bsex, omega)
  cl_70kg * (cohort$WT/70)^0.75 * (cohort$CRCL/90)^bcrcl * exp(bsex*cohort$SEX) *
    exp(rnorm(nrow(cohort), 0, sqrt(omega)))
run_trial <- function(cohort, first_id) {
  cl <- patient_cl(cohort, CL70, BCRCL, BSEX, OM_CL)
  v1 <- V1_70 * (cohort$WT/70)
  records <- list()
  for (i in seq_len(nrow(cohort))) {
    y <- sapply(TIMES, conc2, dose = DOSE, cl = cl[i], v1 = v1[i], q = QQ, v2 = V2) +
      rnorm(length(TIMES), 0, SD_ADD)
    records[[i]] <- data.frame(
      ID = first_id + i, TIME = c(0, TIMES), DV = c(NA_real_, y),
      AMT = c(DOSE, rep(0, length(TIMES))), EVID = c(1L, rep(0L, length(TIMES))),
      CMT = 1L, WT = cohort$WT[i], CRCL = cohort$CRCL[i], SEX = cohort$SEX[i])
  }
  do.call(rbind, records)
}
cohorts <- list(normal = draw_cohort(260L, 78, .198, 95, .217),
                mild = draw_cohort(210L, 76, .198, 62, .198),
                moderate = draw_cohort(180L, 74, .198, 38, .198))
trial_data <- list(normal = run_trial(cohorts$normal, 10000),
                   mild = run_trial(cohorts$mild, 20000),
                   moderate = run_trial(cohorts$moderate, 30000))
cat("trial data generated (true 2-cmt kinetics for all three cohorts)\n")

## normal / Sato: fits (correctly) 2-cmt, allometry fixed, sex on CL
sato_model <- function() {
  ini({ tcl <- log(4.5); tv1 <- log(28); tq <- log(7); tv2 <- log(55)
        bsex <- 0.05; add.err <- 0.2; eta.cl ~ 0.1 })
  model({ cl <- exp(tcl + eta.cl) * (WT/70)^0.75 * exp(bsex*SEX)
          v1 <- exp(tv1) * (WT/70); q <- exp(tq); v2 <- exp(tv2)
          d/dt(centr) <- -cl/v1*centr - q/v1*centr + q/v2*periph
          d/dt(periph) <- q/v1*centr - q/v2*periph
          cp <- centr/v1; cp ~ add(add.err) })
}
## mild / Khan: fits (correctly) 2-cmt, sex on CL and V1
khan_model <- function() {
  ini({ tcl <- log(4.5); tv1 <- log(28); tq <- log(7); tv2 <- log(55)
        bsex <- 0.05; bsexv <- 0; add.err <- 0.2; eta.cl ~ 0.1 })
  model({ cl <- exp(tcl + eta.cl) * (WT/70)^0.75 * exp(bsex*SEX)
          v1 <- exp(tv1) * (WT/70) * exp(bsexv*SEX); q <- exp(tq); v2 <- exp(tv2)
          d/dt(centr) <- -cl/v1*centr - q/v1*centr + q/v2*periph
          d/dt(periph) <- q/v1*centr - q/v2*periph
          cp <- centr/v1; cp ~ add(add.err) })
}
## moderate / Ito: fits (MISSPECIFIED) 1-cmt, both exponents estimated
ito_model <- function() {
  ini({ tcl <- log(4.5); tv <- log(45); bwtcl <- 0.75; bwtv <- 1
        bsex <- 0.05; add.err <- 0.2; eta.cl ~ 0.1 })
  model({ cl <- exp(tcl + eta.cl) * (WT/70)^bwtcl * exp(bsex*SEX)
          v  <- exp(tv) * (WT/70)^bwtv
          cp <- linCmt(); cp ~ add(add.err) })
}
.ctl <- foceiControl(print = 0L, covMethod = "r", covFull = TRUE)
fit_normal   <- nlmixr2est::nlmixr2(sato_model, trial_data$normal,   est = "focei", control = .ctl)
fit_mild     <- nlmixr2est::nlmixr2(khan_model, trial_data$mild,     est = "focei", control = .ctl)
fit_moderate <- nlmixr2est::nlmixr2(ito_model,  trial_data$moderate, est = "focei", control = .ctl)
cat("three per-trial fits done: normal/mild 2-cmt (correct), moderate 1-cmt (misspecified)\n")

with_om_row <- function(f, rse_om = 0.30, eta = "om.eta.cl") {
  C <- f$cov
  if (eta %in% rownames(C)) return(C)
  nm  <- c(rownames(C), eta)
  out <- matrix(0, length(nm), length(nm), dimnames = list(nm, nm))
  out[rownames(C), rownames(C)] <- C
  out[eta, eta] <- (rse_om * as.numeric(f$omega[1, 1]))^2
  out
}
normal_study <- admStudy(model = fit_normal$ui, cov = with_om_row(fit_normal),
                         population = cohorts$normal, dose = DOSE, times = TIMES,
                         stratify = "SEX")
mild_study <- admStudy(model = fit_mild$ui, cov = with_om_row(fit_mild),
                       population = cohorts$mild, dose = DOSE, times = TIMES,
                       stratify = "SEX")
moderate_study <- admStudy(model = fit_moderate$ui, cov = with_om_row(fit_moderate),
                           population = cohorts$moderate, dose = DOSE, times = TIMES,
                           stratify = "SEX")
studies <- admStudies(normal = normal_study, mild = mild_study, moderate = moderate_study)
cat("studies built\n")

## pooled model: 2-cmt, matching normal/mild's structure, renal term on CL
adm_model <- function() {
  ini({ tcl <- log(4); tv1 <- log(28); tq <- log(7); tv2 <- log(55)
        bcrcl <- 0.3; bsex <- 0.05; add.err <- 0.1; eta.cl ~ 0.1 })
  model({ cl <- exp(tcl + eta.cl) * (WT/70)^0.75 * (CRCL/90)^bcrcl * exp(bsex*SEX)
          v1 <- exp(tv1) * (WT/70); q <- exp(tq); v2 <- exp(tv2)
          d/dt(centr) <- -cl/v1*centr - q/v1*centr + q/v2*periph
          d/dt(periph) <- q/v1*centr - q/v2*periph
          cp <- centr/v1; cp ~ add(add.err) })
}

cat("\n==================== srcWeight = \"n\" =========================\n")
fit_old <- nlmixr2est::nlmixr2(adm_model, admData(), est = "adgh",
  control = adghControl(studies = studies, print = 0L, covMethod = "none", srcWeight = "n"))
po <- fit_old$parFixedDf
cat(sprintf("tcl=%.4f tv1=%.4f tq=%.4f tv2=%.4f bcrcl=%.4f bsex=%.4f\n",
            po["tcl","Estimate"], po["tv1","Estimate"], po["tq","Estimate"],
            po["tv2","Estimate"], po["bcrcl","Estimate"], po["bsex","Estimate"]))
cat(sprintf("natural: CL=%.2f V1=%.2f Q=%.2f V2=%.2f\n",
            exp(po["tcl","Estimate"]), exp(po["tv1","Estimate"]),
            exp(po["tq","Estimate"]), exp(po["tv2","Estimate"])))

cat("\n==================== srcWeight = \"cov\" ========================\n")
fit_new <- nlmixr2est::nlmixr2(adm_model, admData(), est = "adgh",
  control = adghControl(studies = studies, print = 0L, covMethod = "none", srcWeight = "cov"))
pn <- fit_new$parFixedDf
cat(sprintf("tcl=%.4f tv1=%.4f tq=%.4f tv2=%.4f bcrcl=%.4f bsex=%.4f\n",
            pn["tcl","Estimate"], pn["tv1","Estimate"], pn["tq","Estimate"],
            pn["tv2","Estimate"], pn["bcrcl","Estimate"], pn["bsex","Estimate"]))
cat(sprintf("natural: CL=%.2f V1=%.2f Q=%.2f V2=%.2f\n",
            exp(pn["tcl","Estimate"]), exp(pn["tv1","Estimate"]),
            exp(pn["tq","Estimate"]), exp(pn["tv2","Estimate"])))

cat(sprintf("\ntruth: CL=%.1f V1=%.1f Q=%.1f V2=%.1f bcrcl=%.2f bsex=%.2f\n",
            CL70, V1_70, QQ, V2, BCRCL, BSEX))
