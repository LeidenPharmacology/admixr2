# Parametrised version of the CONTROL arm: every source richly sampled and
# two-compartment, so no structural disagreement. What varies here is the
# NUMBER of sources, the subjects per source, and the WITHIN-source spread of
# the covariate -- the three things that separate "finite-source artefact" from
# "the variance channel" from "inherent".
Sys.setenv(NOT_CRAN = "true"); options(warn = 1)
suppressMessages({library(nlmixr2); library(rxode2)
                  pkgload::load_all(".", quiet = TRUE)})
NSRC   <- as.integer(Sys.getenv("NSRC",   "3"))
NPER   <- as.integer(Sys.getenv("NPER",   "150"))
CRCLSD <- as.numeric(Sys.getenv("CRCLSD", "0.217"))
REPS   <- as.integer(Sys.getenv("REPS",   "8"))
CSV    <- Sys.getenv("CSV")
cat(sprintf("NSRC=%d  NPER=%d  CRCL sdlog=%.3f  reps=%d\n", NSRC, NPER, CRCLSD, REPS))

CL<-5; V<-50; Q<-8; VP<-70; BCRCL<-0.6; BSEX<-0.18; OM<-0.05; SD<-0.10
mod <- rxode2::rxode2({ cl <- z1; v <- z2; q <- z3; vp <- z4; cp <- linCmt() })
Rc  <- matrix(c(1,.45,.45,1), 2, 2)
OBS <- c(240.083, 240.25, 240.5, 241, 242, 244, 248, 251.9)
CRCL_MED <- seq(38, 95, length.out = NSRC)
WT_MED   <- seq(74, 78, length.out = NSRC)

cohort <- function(n, wt_med, crcl_med) {
  z <- matrix(rnorm(2*n), n, 2) %*% chol(Rc)
  data.frame(WT   = qlnorm(pnorm(z[,1]), log(wt_med),   .198),
             CRCL = qlnorm(pnorm(z[,2]), log(crcl_med), CRCLSD),
             SEX  = rbinom(n, 1L, .55))
}
sim_trial <- function(coh, id0) {
  n <- nrow(coh)
  cli <- CL*(coh$WT/70)^0.75*(coh$CRCL/90)^BCRCL*exp(BSEX*coh$SEX)*
         exp(rnorm(n,0,sqrt(OM)))
  dosing <- rxode2::et(amt=200, ii=12, addl=20)
  dt <- seq(0, by=12, length.out=21)
  out <- rxode2::rxSolve(mod, rxode2::add.sampling(dosing, OBS),
           params=data.frame(z1=cli, z2=V*(coh$WT/70), z3=Q, z4=VP*(coh$WT/70)),
           returnType="data.frame")
  cps <- split(out$cp, out$sim.id)
  do.call(rbind, lapply(seq_len(n), function(i)
    rbind(data.frame(ID=id0+i, TIME=dt, DV=NA_real_, AMT=200, EVID=1L,
                     WT=coh$WT[i], CRCL=coh$CRCL[i], SEX=coh$SEX[i]),
          data.frame(ID=id0+i, TIME=OBS, DV=cps[[i]]+rnorm(length(OBS),0,SD),
                     AMT=0, EVID=0L,
                     WT=coh$WT[i], CRCL=coh$CRCL[i], SEX=coh$SEX[i]))))
}
two_cmt <- function() {
  ini({ tcl <- log(4); tv <- log(40); tq <- log(6); tvp <- log(50)
        bsex <- 0.05; add.err <- 0.2; eta.cl ~ 0.1 })
  model({ cl <- exp(tcl+eta.cl)*(WT/70)^0.75*exp(bsex*SEX); v <- exp(tv)*(WT/70)
          q <- exp(tq); vp <- exp(tvp)*(WT/70); cp <- linCmt(); cp ~ add(add.err) })
}
pooled <- function() {
  ini({ tcl <- log(4.5); tv <- log(55); tq <- log(6); tvp <- log(60)
        bcrcl <- 0.3; bsex <- 0.05; add.err <- 0.15; eta.cl ~ 0.08 })
  model({ cl <- exp(tcl+eta.cl)*(WT/70)^0.75*(CRCL/90)^bcrcl*exp(bsex*SEX)
          v <- exp(tv)*(WT/70); q <- exp(tq); vp <- exp(tvp)*(WT/70)
          cp <- linCmt(); cp ~ add(add.err) })
}
ct <- foceiControl(print=0L, covMethod="r", covFull=TRUE)
res <- list()
for (.rep in seq_len(REPS)) {
  set.seed(900 + .rep)
  ok <- try({
    cohs <- lapply(seq_len(NSRC), function(j) cohort(NPER, WT_MED[j], CRCL_MED[j]))
    fits <- lapply(seq_len(NSRC), function(j)
      nlmixr2(two_cmt, sim_trial(cohs[[j]], j*10000), est="focei", control=ct))
    sts <- lapply(seq_len(NSRC), function(j)
      admStudy(model=fits[[j]]$ui, cov=fits[[j]]$cov, population=cohs[[j]],
               ev=rxode2::et(amt=200, ii=12, addl=20), times=OBS, stratify="SEX"))
    names(sts) <- paste0("src", seq_len(NSRC))
    fit <- nlmixr2(pooled, admData(), est="adgh",
                   control=adghControl(studies=do.call(admStudies, sts),
                                       print=0L, cores=2L))
    p <- fit$parFixedDf
    row <- data.frame(rep=.rep, nsrc=NSRC, nper=NPER, crclsd=CRCLSD,
      bcrcl=p["bcrcl","Estimate"], bcrcl_se=p["bcrcl","SE"],
      bsex=p["bsex","Estimate"],   bsex_se=p["bsex","SE"],
      om=as.numeric(fit$omega[1,1]), cov_method=fit$covMethod,
      stringsAsFactors=FALSE)
    res[[length(res)+1L]] <- row
    write.csv(do.call(rbind,res), CSV, row.names=FALSE)
    cat(sprintf("rep %2d ok: bcrcl %.3f  bsex %.3f  om %.4f  %s\n",
                .rep, row$bcrcl, row$bsex, row$om, row$cov_method))
  }, silent = TRUE)
  if (inherits(ok,"try-error"))
    cat("rep", .rep, "FAILED:", conditionMessage(attr(ok,"condition")), "\n")
}
cat("\n== done", length(res), "of", REPS, "\n")
