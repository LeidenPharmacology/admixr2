Sys.setenv(NOT_CRAN = "true"); options(warn = 1)
suppressMessages({library(nlmixr2); library(rxode2)
                  pkgload::load_all(".", quiet = TRUE)})

NTRIAL <- 150L
REPS   <- 8L
CSV    <- "C:/Users/hidde/.claude/jobs/319faff2/tmp/arm_pinned.csv"
res <- list()
for (.rep in seq_len(REPS)) {
  set.seed(500 + .rep)
  ok <- try({


## ---- TRUTH: 2 compartments, 200 mg q12h, three covariates ---------------
CL<-5; V<-50; Q<-8; VP<-70; BCRCL<-0.6; BSEX<-0.18; OM<-0.05; SD<-0.10
mod <- rxode2::rxode2({ cl <- z1; v <- z2; q <- z3; vp <- z4; cp <- linCmt() })
Rc  <- matrix(c(1,.45,.45,1), 2, 2)

cohort <- function(n, wt_med, crcl_med, p_male = 0.55) {
  z <- matrix(rnorm(2*n), n, 2) %*% chol(Rc)
  data.frame(WT   = qlnorm(pnorm(z[,1]), log(wt_med),   .198),
             CRCL = rep(crcl_med, n),   # NO within-cohort spread
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
OBS_A <- c(240.083, 240.25, 240.5, 241, 242, 244, 248, 251.9)  # rich, SS
OBS_B <- c(0.083, 0.25, 0.5, 1, 2, 4, 8, 12)  # dense, first interval
OBS_C <- c(240.083, 240.25, 240.5, 241, 242, 244, 248, 251.9)  # rich, SS
cohA <- cohort(NTRIAL, 74, 38)      # moderate impairment, sparse design
cohB <- cohort(NTRIAL, 78, 95)      # normal function,     dense design
cohC <- cohort(NTRIAL, 76, 62)      # mild impairment,     intermediate design
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
fA <- nlmixr2(two_cmt, dA, est="focei", control=ct)   # rich -> 2 cmt
fB <- nlmixr2(two_cmt, dB, est="focei", control=ct)   # dense early  -> 2 cmt
fC <- nlmixr2(two_cmt, dC, est="focei", control=ct)   # rich -> 2 cmt
if (FALSE) for (nm in c("A","B","C")) {
  f <- get(paste0("f", nm)); p <- f$parFixedDf
  cat(sprintf("paper %s: CL %.2f  V %.1f  bsex %+.3f%s\n", nm,
      exp(p["tcl","Estimate"]), exp(p["tv","Estimate"]), p["bsex","Estimate"],
      if ("tq" %in% rownames(p)) sprintf("  Q %.2f  Vp %.1f",
          exp(p["tq","Estimate"]), exp(p["tvp","Estimate"])) else ""))
}

## ---- meta-analysis: renal term that NO paper contains -------------------
ev_ss <- function() rxode2::et(amt = 200, ii = 12, addl = 20)
st <- admStudies(
  moderate_sparse = admStudy(model=fA$ui, cov=fA$cov, population=cohA[, c('WT','SEX')], at=list(CRCL=38),
                             ev=ev_ss(), times=OBS_A, stratify="SEX"),
  normal_dense    = admStudy(model=fB$ui, cov=fB$cov, population=cohB[, c('WT','SEX')], at=list(CRCL=95),
                             ev=rxode2::et(amt=200), times=OBS_B, stratify="SEX"),
  mild_interm     = admStudy(model=fC$ui, cov=fC$cov, population=cohC[, c('WT','SEX')], at=list(CRCL=62),
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
    # RECORD WHAT THE SANDWICH ACTUALLY DID. A silent degradation to "r" would
    # otherwise hide inside an averaged coverage number -- the exact failure
    # that cost the earlier render its standard errors.
    row <- data.frame(rep = .rep,
      CL = exp(p["tcl","Estimate"]),  CL_se = p["tcl","SE"],
      V  = exp(p["tv","Estimate"]),   V_se  = p["tv","SE"],
      Q  = exp(p["tq","Estimate"]),   Q_se  = p["tq","SE"],
      Vp = exp(p["tvp","Estimate"]),  Vp_se = p["tvp","SE"],
      bcrcl = p["bcrcl","Estimate"],  bcrcl_se = p["bcrcl","SE"],
      bsex  = p["bsex","Estimate"],   bsex_se  = p["bsex","SE"],
      om = as.numeric(fit$omega[1,1]),
      omA = as.numeric(fA$omega[1,1]), omB = as.numeric(fB$omega[1,1]), omC = as.numeric(fC$omega[1,1]),
      bsexA = fA$parFixedDf["bsex","Estimate"],
      bsexB = fB$parFixedDf["bsex","Estimate"],
      bsexC = fC$parFixedDf["bsex","Estimate"],
      covdimA = nrow(fA$cov), covdimB = nrow(fB$cov), covdimC = nrow(fC$cov),
      om_all = all(vapply(list(fA,fB,fC),
                   function(f) "om.eta.cl" %in% rownames(f$cov), logical(1))),
      sandwich = !is.null(admixr2:::.admSandwichParts(fit)),
      cov_method = fit$covMethod, stringsAsFactors = FALSE)
    res[[length(res)+1L]] <- row
    write.csv(do.call(rbind, res), CSV, row.names = FALSE)
    cat(sprintf("rep %2d ok: bcrcl %.3f  om %.4f  (sources %.3f/%.3f/%.3f)  %s
",
        .rep, row$bcrcl, row$om, row$omA, row$omB, row$omC, row$cov_method))
  }, silent = TRUE)
  if (inherits(ok, "try-error"))
    cat("rep", .rep, "FAILED:", conditionMessage(attr(ok,"condition")), "
")
}
cat("
== done,", length(res), "of", REPS, "
")
