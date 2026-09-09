## Faithful reproduction of vignettes/covariates.Rmd -- three published models
## (Sato/Ito/Khan), each a MODEL source with its own reported covariance,
## combined via admStudy()/admStudies(), sex banded, renal function
## marginalised (no trial fit it). Fit twice: srcWeight = "n" (today's
## default, what the vignette's own published numbers reflect) vs
## srcWeight = "cov" (this branch's fix). Point estimate comparison only --
## the SE side is not fixed (see HANDOFF.md), so covMethod = "none" throughout
## to keep the comparison honest about what's actually being tested.
suppressMessages({ library(nlmixr2est); devtools::load_all(".", quiet = TRUE) })
rxode2::setRxThreads(2L)

N_MULT <- as.numeric(Sys.getenv("N_MULT", "1"))
set.seed(11)
TIMES <- c(0.5, 1, 2, 4, 8, 12, 24); DOSE <- 200
CL70 <- 5; V70 <- 50; BCRCL <- 0.6; BSEX <- 0.18
OM_CL <- 0.05; SD_ADD <- 0.08
R <- matrix(c(1, .45, .45, 1), 2, 2,
            dimnames = list(c("WT", "CRCL"), c("WT", "CRCL")))

draw_cohort <- function(n, wt_median, wt_cv, crcl_median, crcl_cv, p_male = 0.55) {
  z      <- matrix(rnorm(2 * n), nrow = n, ncol = 2) %*% chol(R)
  u_wt   <- pnorm(z[, 1]); u_crcl <- pnorm(z[, 2])
  data.frame(WT = qlnorm(u_wt, meanlog = log(wt_median), sdlog = wt_cv),
             CRCL = qlnorm(u_crcl, meanlog = log(crcl_median), sdlog = crcl_cv),
             SEX = rbinom(n, size = 1L, prob = p_male))
}
conc <- function(t, dose, cl, v) dose/v * exp(-cl/v * t)
patient_cl <- function(cohort, cl_70kg, bcrcl, bsex, omega) {
  cl_70kg * (cohort$WT/70)^0.75 * (cohort$CRCL/90)^bcrcl *
    exp(bsex * cohort$SEX) * exp(rnorm(nrow(cohort), 0, sqrt(omega)))
}
run_trial <- function(cohort, first_id) {
  cl <- patient_cl(cohort, CL70, BCRCL, BSEX, OM_CL)
  v  <- V70 * (cohort$WT/70)
  records <- list()
  for (i in seq_len(nrow(cohort))) {
    y <- conc(TIMES, DOSE, cl[i], v[i]) + rnorm(length(TIMES), 0, SD_ADD)
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

sato_model <- function() {
  ini({ tcl <- log(4.5); tv <- log(45); bsex <- 0.05; add.err <- 0.2; eta.cl ~ 0.1 })
  model({ cl <- exp(tcl + eta.cl) * (WT/70)^0.75 * exp(bsex * SEX)
          v  <- exp(tv) * (WT/70); cp <- linCmt(); cp ~ add(add.err) })
}
ito_model <- function() {
  ini({ tcl <- log(4.5); tv <- log(45); bwtcl <- 0.75; bwtv <- 1
        bsex <- 0.05; add.err <- 0.2; eta.cl ~ 0.1 })
  model({ cl <- exp(tcl + eta.cl) * (WT/70)^bwtcl * exp(bsex * SEX)
          v  <- exp(tv) * (WT/70)^bwtv; cp <- linCmt(); cp ~ add(add.err) })
}
khan_model <- function() {
  ini({ tcl <- log(4.5); tv <- log(45); bsex <- 0.05; bsexv <- 0
        add.err <- 0.2; eta.cl ~ 0.1 })
  model({ cl <- exp(tcl + eta.cl) * (WT/70)^0.75 * exp(bsex * SEX)
          v  <- exp(tv) * (WT/70) * exp(bsexv * SEX); cp <- linCmt(); cp ~ add(add.err) })
}
.ctl <- foceiControl(print = 0L, covMethod = "r", covFull = TRUE)
fit_normal   <- nlmixr2est::nlmixr2(sato_model, trial_data$normal,   est = "focei", control = .ctl)
fit_moderate <- nlmixr2est::nlmixr2(ito_model,  trial_data$moderate, est = "focei", control = .ctl)
fit_mild     <- nlmixr2est::nlmixr2(khan_model, trial_data$mild,     est = "focei", control = .ctl)
cat("three per-trial FOCEI fits done\n")

with_om_row <- function(f, rse_om = 0.30, eta = "om.eta.cl") {
  C <- f$cov
  if (eta %in% rownames(C)) return(C)
  nm  <- c(rownames(C), eta)
  out <- matrix(0, length(nm), length(nm), dimnames = list(nm, nm))
  out[rownames(C), rownames(C)] <- C
  out[eta, eta] <- (rse_om * as.numeric(f$omega[1, 1]))^2
  out
}
normal_paper   <- list(COV = with_om_row(fit_normal))
mild_paper     <- list(COV = with_om_row(fit_mild))
moderate_paper <- list(COV = with_om_row(fit_moderate))

sato_published <- fit_normal$ui; ito_published <- fit_moderate$ui; khan_published <- fit_mild$ui
normal_study <- admStudy(model = sato_published, cov = normal_paper$COV,
                         population = cohorts$normal, dose = DOSE, times = TIMES,
                         stratify = "SEX", n = NROW(cohorts$normal) * N_MULT)
moderate_study <- admStudy(model = ito_published, cov = moderate_paper$COV,
                           population = cohorts$moderate, dose = DOSE, times = TIMES,
                           stratify = "SEX")
mild_study <- admStudy(model = khan_published, cov = mild_paper$COV,
                       population = cohorts$mild, dose = DOSE, times = TIMES,
                       stratify = "SEX")
studies <- admStudies(normal = normal_study, moderate = moderate_study, mild = mild_study)
cat("studies built\n")

adm_model <- function() {
  ini({ tcl <- log(4); tv <- log(45); bcrcl <- 0.3; bsex <- 0.05
        add.err <- 0.1; eta.cl ~ 0.1 })
  model({ cl <- exp(tcl + eta.cl) * (WT/70)^0.75 * (CRCL/90)^bcrcl * exp(bsex * SEX)
          v  <- exp(tv) * (WT/70); cp <- linCmt(); cp ~ add(add.err) })
}

cat("\n==================== srcWeight = \"n\" (matches the vignette's own numbers) ====\n")
fit_old <- nlmixr2est::nlmixr2(adm_model, admData(), est = "adgh",
  control = adghControl(studies = studies, print = 0L, covMethod = "none", srcWeight = "n"))
po <- fit_old$parFixedDf
cat(sprintf("tcl=%.10f tv=%.10f bcrcl=%.10f bsex=%.10f\n",
            po["tcl","Estimate"], po["tv","Estimate"], po["bcrcl","Estimate"], po["bsex","Estimate"]))
cat(sprintf("back-transformed CL=%.3f V=%.3f\n", exp(po["tcl","Estimate"]), exp(po["tv","Estimate"])))

cat("\n==================== srcWeight = \"cov\" (this branch's fix) ===================\n")
fit_new <- nlmixr2est::nlmixr2(adm_model, admData(), est = "adgh",
  control = adghControl(studies = studies, print = 0L, covMethod = "none", srcWeight = "cov"))
pn <- fit_new$parFixedDf
cat(sprintf("tcl=%.10f tv=%.10f bcrcl=%.10f bsex=%.10f\n",
            pn["tcl","Estimate"], pn["tv","Estimate"], pn["bcrcl","Estimate"], pn["bsex","Estimate"]))
cat(sprintf("back-transformed CL=%.3f V=%.3f\n", exp(pn["tcl","Estimate"]), exp(pn["tv","Estimate"])))

cat(sprintf("\ntruth: CL=%.1f V=%.1f bcrcl=%.2f bsex=%.2f\n", CL70, V70, BCRCL, BSEX))
