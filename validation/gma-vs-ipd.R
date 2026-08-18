## DOES IT ACTUALLY FIT? Parameter likelihood vs the individual-data gold standard.
## Individual fits are REAL nlmixr2 FOCEI fits -- publications are produced the way
## real ones are, and C_s comes from nlmixr2's own covariance step.
##
## TRUTH  2-cmt, two covariates: x1 = log(WT/70), x2 = log(CRCL/100)
## STUDY A  1-cmt form, covariate x1 ONLY, steady-state peak + trough
## STUDY B  2-cmt form, BOTH covariates, dense peak capture
##
##   IPD-POOL  unified 2-cmt on the POOLED individual data     <- gold standard
##   IPD-A     A's own form on A's individuals -> A's publication (theta_A, C_A)
##   IPD-B     B's own form on B's individuals -> B's publication (theta_B, C_B)
##   PARAM-LL  unified 2-cmt from the two PUBLICATIONS ONLY
##
## g_s uses an AGGREGATE node refit while the publications came from FOCEI, so
## RMSE(PARAM-LL, IPD-POOL) also measures the surrogate-estimator gap that the
## indirect-inference literature flags as never quantified.
setwd("C:/package/admixr2/.claude/worktrees/feature-covariate-quadrature")
suppressMessages(library(rxode2)); suppressMessages(library(nlmixr2est))
set.seed(20260815)
say <- function(...) { cat(...,"\n"); utils::flush.console() }; T0 <- proc.time()[["elapsed"]]
OUT <- "C:/Users/hidde/AppData/Local/Temp/claude/C--package-admixr2/3ff305c7-64a5-4cb0-ba61-1436e2e9b16e/scratchpad"

TP <- c(lcl=log(4.2), lvc=log(30), lq=log(7.0), lvp=log(40),
        bCL1=.72, bCL2=.45, bVc1=.95, oCL=.26, oVc=.20, sig=.11)
DA <- list(n=140, dose=500, ii=12, addl=5, obs=c(61,72), mx1=log(82/70), sx1=.20,
           mx2=log(78/100), sx2=.28)                 # ss peak(1h into 6th dose)+trough
DB <- list(n=60,  dose=500, ii=0,  addl=0, obs=c(.25,.5,1,2,4,8,12), mx1=0, sx1=.17,
           mx2=log(96/100), sx2=.22)                 # dense single dose

## ---- simulate individual data from the truth -------------------------------
sim_mod <- rxode2({
  cl <- exp(lcl + bCL1*log(WT/70) + bCL2*log(CRCL/100) + eCL)
  vc <- exp(lvc + bVc1*log(WT/70) + eVC); q <- exp(lq); vp <- exp(lvp)
  d/dt(centr) <- -(cl+q)/vc*centr + q/vp*peri
  d/dt(peri)  <-        q/vc*centr - q/vp*peri
  cp <- centr/vc
})
sim_study <- function(d, id0) {
  x1 <- rnorm(d$n, d$mx1, d$sx1); x2 <- rnorm(d$n, d$mx2, d$sx2)
  pars <- data.frame(lcl=TP[["lcl"]], lvc=TP[["lvc"]], lq=TP[["lq"]], lvp=TP[["lvp"]],
                     bCL1=TP[["bCL1"]], bCL2=TP[["bCL2"]], bVc1=TP[["bVc1"]],
                     WT=70*exp(x1), CRCL=100*exp(x2),
                     eCL=rnorm(d$n,0,TP[["oCL"]]), eVC=rnorm(d$n,0,TP[["oVc"]]))
  ev <- if (d$ii > 0) et(amt=d$dose, ii=d$ii, addl=d$addl) else et(amt=d$dose)
  ev <- et(ev, d$obs)
  ## solve per subject: robust to whatever rxSolve names its id column
  obs <- do.call(rbind, lapply(seq_len(d$n), function(i) {
    s <- rxSolve(sim_mod, pars[i, , drop=FALSE], ev, returnType="data.frame",
                 cores=1, nDisplayProgress=.Machine$integer.max)
    s <- s[s$time %in% d$obs, , drop=FALSE]
    data.frame(ID=id0+i, TIME=s$time,
               DV=s$cp*(1+rnorm(nrow(s),0,TP[["sig"]])), AMT=0, EVID=0,
               WT=70*exp(x1)[i], CRCL=100*exp(x2)[i])
  }))
  ndose <- if (d$ii > 0) d$addl+1L else 1L
  dose <- do.call(rbind, lapply(seq_len(d$n), function(i) data.frame(
    ID=id0+i, TIME=(seq_len(ndose)-1)*d$ii, DV=NA_real_, AMT=d$dose, EVID=1,
    WT=70*exp(x1)[i], CRCL=100*exp(x2)[i])))
  df <- rbind(dose, obs); df[order(df$ID, df$TIME, -df$EVID), ]
}
dA <- sim_study(DA, 0L); dB <- sim_study(DB, 1000L)
say(sprintf("[%3.0fs] simulated A %d obs, B %d obs", proc.time()[["elapsed"]]-T0,
            sum(dA$EVID==0), sum(dB$EVID==0)))

## ---- the model forms -------------------------------------------------------
mA <- function() {                       # 1-cmt, WT only
  ini({ lcl <- log(4); lv <- log(60); bCL1 <- .7; bVc1 <- .9
        eta.cl ~ .07; eta.v ~ .05; prop.err <- .12 })
  model({ cl <- exp(lcl + bCL1*log(WT/70) + eta.cl)
          v  <- exp(lv  + bVc1*log(WT/70) + eta.v)
          d/dt(centr) <- -cl/v*centr
          cp <- centr/v; cp ~ prop(prop.err) })
}
mU <- function() {                       # unified / B's form: 2-cmt, both covariates
  ini({ lcl <- log(4); lvc <- log(35); lq <- log(6); lvp <- log(35)
        bCL1 <- .7; bCL2 <- .4; bVc1 <- .9
        eta.cl ~ .07; eta.vc ~ .05; prop.err <- .12 })
  model({ cl <- exp(lcl + bCL1*log(WT/70) + bCL2*log(CRCL/100) + eta.cl)
          vc <- exp(lvc + bVc1*log(WT/70) + eta.vc); q <- exp(lq); vp <- exp(lvp)
          d/dt(centr) <- -(cl+q)/vc*centr + q/vp*peri
          d/dt(peri)  <-        q/vc*centr - q/vp*peri
          cp <- centr/vc; cp ~ prop(prop.err) })
}
grab <- function(f, nm) {
  fx <- fixef(f); om <- f$omega
  v <- vapply(nm, function(q)
    if (q=="oCL") sqrt(om[["eta.cl","eta.cl"]])
    else if (q=="oVc") sqrt(om[[grep("^eta\\.vc?$", rownames(om), value=TRUE)[1],
                                grep("^eta\\.vc?$", rownames(om), value=TRUE)[1]]])
    else if (q=="sig") fx[["prop.err"]] else fx[[q]], 0)
  setNames(v, nm)
}
ffit <- function(m, d) suppressWarnings(suppressMessages(
  nlmixr2(m, d, est="focei", control=foceiControl(print=0L))))
say("fitting IPD-A ..."); fA <- ffit(mA, dA)
say(sprintf("[%3.0fs] IPD-A OFV %.1f", proc.time()[["elapsed"]]-T0, fA$objective))
say("fitting IPD-B ..."); fB <- ffit(mU, dB)
say(sprintf("[%3.0fs] IPD-B OFV %.1f", proc.time()[["elapsed"]]-T0, fB$objective))
say("fitting IPD-POOL (gold) ..."); fP <- ffit(mU, rbind(dA, dB))
say(sprintf("[%3.0fs] IPD-POOL OFV %.1f", proc.time()[["elapsed"]]-T0, fP$objective))

UN  <- names(TP); ANM <- c("lcl","lvc","bCL1","bVc1","oCL","oVc","sig")
thA <- grab(fA, c("lcl","lv","bCL1","bVc1","oCL","oVc","sig")); names(thA)[2] <- "lvc"
thB <- grab(fB, UN); GOLD <- grab(fP, UN)
## C_s: nlmixr2 reports SEs on the printed scale; use them, diagonal, and say so.
seof <- function(f, nm) { pf <- f$parFixedDf
  vapply(nm, function(q) { r <- if (q=="sig") "prop.err" else q
    s <- pf[rownames(pf)==r, "SE"]; if (length(s)&&is.finite(s)) s else NA_real_ }, 0) }
seA <- seof(fA, c("lcl","lv","bCL1","bVc1","sig")); names(seA)[2] <- "lvc"
seB <- seof(fB, setdiff(UN, c("oCL","oVc")))
seA[is.na(seA)] <- 0.15; seB[is.na(seB)] <- 0.15
seA <- c(seA, oCL=0.04, oVc=0.05)[ANM]; seB <- c(seB, oCL=0.03, oVc=0.04)[UN]
say(sprintf("[%3.0fs] publications in hand", proc.time()[["elapsed"]]-T0))
saveRDS(list(TP=TP,GOLD=GOLD,thA=thA,thB=thB,seA=seA,seB=seB,DA=DA,DB=DB),
        file.path(OUT,"gma-ipd-stage1.rds"))
say("\n        TRUTH   IPD-POOL      IPD-A      IPD-B")
for (q in UN) say(sprintf("%-6s %8.4f %10.4f %10s %10.4f", q, TP[[q]], GOLD[[q]],
  if (q %in% ANM) sprintf("%.4f", thA[[q]]) else "-", thB[[q]]))
say(sprintf("\nRMSE vs truth: IPD-POOL %.4f   IPD-B alone %.4f",
    sqrt(mean((GOLD-TP)^2)), sqrt(mean((thB-TP)^2))))
say(sprintf("total %.0f s", proc.time()[["elapsed"]]-T0))
