Sys.setenv(NOT_CRAN = "true"); options(warn = 1)
suppressMessages({library(nlmixr2); library(rxode2)
                  pkgload::load_all(".", quiet = TRUE)})
set.seed(31)

## ---- TRUTH: 2 compartments, 200 mg q12h, three covariates ---------------
CL<-5; V<-50; Q<-8; VP<-70; BCRCL<-0.6; BSEX<-0.18; OM<-0.05; SD<-0.10
mod <- rxode2::rxode2({ cl <- z1; v <- z2; q <- z3; vp <- z4; cp <- linCmt() })
Rc  <- matrix(c(1,.45,.45,1), 2, 2)

cohort <- function(n, wt_med, crcl_med, p_male = 0.55) {
  z <- matrix(rnorm(2*n), n, 2) %*% chol(Rc)
  data.frame(WT   = qlnorm(pnorm(z[,1]), log(wt_med),   .198),
             CRCL = qlnorm(pnorm(z[,2]), log(crcl_med), .217),
             SEX  = rbinom(n, 1L, p_male))
}
sim_trial <- function(coh, addl, obs, id0) {
  n   <- nrow(coh)
  cli <- CL * (coh$WT/70)^0.75 * (coh$CRCL/90)^BCRCL * exp(BSEX*coh$SEX) *
         exp(rnorm(n, 0, sqrt(OM)))
  dosing <- if (addl > 0L) rxode2::et(amt=200, ii=12, addl=addl) else rxode2::et(amt=200)
  dose_times <- if (addl > 0L) seq(0, by=12, length.out=addl+1L) else 0
  out <- rxode2::rxSolve(mod, rxode2::add.sampling(dosing, obs),
           params = data.frame(z1=cli, z2=V*(coh$WT/70), z3=Q, z4=VP*(coh$WT/70)),
           returnType="data.frame")
  cps <- split(out$cp, out$sim.id)
  do.call(rbind, lapply(seq_len(n), function(i)
    rbind(data.frame(ID=id0+i, TIME=dose_times, DV=NA_real_, AMT=200, EVID=1L,
                     WT=coh$WT[i], CRCL=coh$CRCL[i], SEX=coh$SEX[i]),
          data.frame(ID=id0+i, TIME=obs, DV=cps[[i]]+rnorm(length(obs),0,SD),
                     AMT=0, EVID=0L,
                     WT=coh$WT[i], CRCL=coh$CRCL[i], SEX=coh$SEX[i]))))
}

## renal function varies BETWEEN trials; each trial also has its own design
OBS_A <- c(239.9, 241)                        # SS trough + peak
OBS_B <- c(0.083, 0.25, 0.5, 1, 2, 4, 8, 12)  # dense, first interval
OBS_C <- c(240.5, 241, 244, 248, 251.9)       # intermediate, SS interval
cohA <- cohort(60, 74, 38)      # moderate impairment, sparse design
cohB <- cohort(60, 78, 95)      # normal function,     dense design
cohC <- cohort(60, 76, 62)      # mild impairment,     intermediate design
dA <- sim_trial(cohA, 20L, OBS_A, 10000)
dB <- sim_trial(cohB,  0L, OBS_B, 20000)
dC <- sim_trial(cohC, 20L, OBS_C, 30000)

## ---- what each analyst publishes: sex yes, renal no ---------------------
one_cmt <- function() {
  ini({ tcl <- log(4); tv <- log(60); bsex <- 0.05; add.err <- 0.2; eta.cl ~ 0.1 })
  model({ cl <- exp(tcl + eta.cl)*(WT/70)^0.75*exp(bsex*SEX); v <- exp(tv)*(WT/70)
          cp <- linCmt(); cp ~ add(add.err) })
}
two_cmt <- function() {
  ini({ tcl <- log(4); tv <- log(40); tq <- log(6); tvp <- log(50)
        bsex <- 0.05; add.err <- 0.2; eta.cl ~ 0.1 })
  model({ cl <- exp(tcl + eta.cl)*(WT/70)^0.75*exp(bsex*SEX); v <- exp(tv)*(WT/70)
          q <- exp(tq); vp <- exp(tvp)*(WT/70); cp <- linCmt(); cp ~ add(add.err) })
}
ct <- foceiControl(print = 0L, covMethod = "r", covFull = TRUE)
fA <- nlmixr2(one_cmt, dA, est="focei", control=ct)   # peak/trough  -> 1 cmt
fB <- nlmixr2(two_cmt, dB, est="focei", control=ct)   # dense early  -> 2 cmt
fC <- nlmixr2(two_cmt, dC, est="focei", control=ct)   # intermediate -> 2 cmt
for (nm in c("A","B","C")) {
  f <- get(paste0("f", nm)); p <- f$parFixedDf
  cat(sprintf("paper %s: CL %.2f  V %.1f  bsex %+.3f%s\n", nm,
      exp(p["tcl","Estimate"]), exp(p["tv","Estimate"]), p["bsex","Estimate"],
      if ("tq" %in% rownames(p)) sprintf("  Q %.2f  Vp %.1f",
          exp(p["tq","Estimate"]), exp(p["tvp","Estimate"])) else ""))
}

## ---- meta-analysis: renal term that NO paper contains -------------------
ev_ss <- function() rxode2::et(amt = 200, ii = 12, addl = 20)
st <- admStudies(
  moderate_sparse = admStudy(model=fA$ui, cov=fA$cov, population=cohA,
                             ev=ev_ss(), times=OBS_A, stratify="SEX"),
  normal_dense    = admStudy(model=fB$ui, cov=fB$cov, population=cohB,
                             ev=rxode2::et(amt=200), times=OBS_B, stratify="SEX"),
  mild_interm     = admStudy(model=fC$ui, cov=fC$cov, population=cohC,
                             ev=ev_ss(), times=OBS_C, stratify="SEX"))
print(st)

pooled <- function() {
  ini({ tcl <- log(4.5); tv <- log(55); tq <- log(6); tvp <- log(60)
        bcrcl <- 0.3; bsex <- 0.05; add.err <- 0.15; eta.cl ~ 0.08 })
  model({ cl <- exp(tcl + eta.cl)*(WT/70)^0.75*(CRCL/90)^bcrcl*exp(bsex*SEX)
          v <- exp(tv)*(WT/70); q <- exp(tq); vp <- exp(tvp)*(WT/70)
          cp <- linCmt(); cp ~ add(add.err) })
}
fit <- nlmixr2(pooled, admData(), est="adgh",
               control=adghControl(studies=st, print=0L, cores=2L))
p <- fit$parFixedDf
cat("\n===== META-ANALYSIS vs TRUTH =====\n")
sl <- function(l,t,e,s) cat(sprintf("  %-6s truth %6.2f   pooled %6.2f  [%6.2f, %6.2f]\n",
                                    l,t,e,e*exp(-1.96*s),e*exp(1.96*s)))
id <- function(l,t,e,s) cat(sprintf("  %-6s truth %6.3f   pooled %6.3f  [%6.3f, %6.3f]\n",
                                    l,t,e,e-1.96*s,e+1.96*s))
sl("CL", CL, exp(p["tcl","Estimate"]),  p["tcl","SE"])
sl("V",  V,  exp(p["tv","Estimate"]),   p["tv","SE"])
sl("Q",  Q,  exp(p["tq","Estimate"]),   p["tq","SE"])
sl("Vp", VP, exp(p["tvp","Estimate"]),  p["tvp","SE"])
id("bcrcl", BCRCL, p["bcrcl","Estimate"], p["bcrcl","SE"])
id("bsex",  BSEX,  p["bsex","Estimate"],  p["bsex","SE"])
cat("  covMethod:", fit$covMethod, "\n")
