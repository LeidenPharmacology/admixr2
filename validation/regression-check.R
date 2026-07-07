## Empirical check: with the CORRECT model (IIV on ka, cl, v), do both the var
## and cov aggregate fits recover the individual FOCEI fit on Theophylline?
WT <- "C:/package/admixr2/.claude/worktrees/indometh-validation"
suppressMessages(pkgload::load_all(WT, quiet = TRUE, export_all = TRUE))
suppressMessages({library(rxode2); library(nlmixr2)})
set.seed(1L); DOSE <- 320

e <- new.env(); utils::data("theo_sd", package="nlmixr2data", envir=e); theo <- get("theo_sd", e)
ids <- unique(theo$ID); obs <- theo[theo$EVID==0,]; dz <- theo[theo$EVID!=0,]
amt <- setNames(dz$AMT, dz$ID); K <- 11L
rt <- dvn <- matrix(NA_real_, length(ids), K)
for (i in seq_along(ids)){ s<-obs[obs$ID==ids[i],]; s<-s[order(s$TIME),]; rt[i,]<-s$TIME
  dvn[i,]<-s$DV*(DOSE/amt[[as.character(ids[i])]]) }
times <- round(apply(rt,2,median),3)
E <- colMeans(dvn); V <- cov.wt(dvn, method="ML")$cov; n <- length(ids)

## correct structure: IIV on ka, cl AND v
theo3 <- function() {
  ini({ tka<-log(1.5); tcl<-log(2.8); tv<-log(32)
        prop.sd<-c(0,0.1); add.sd<-c(0,0.3); eta.ka~0.2; eta.cl~0.09; eta.v~0.09 })
  model({ ka<-exp(tka+eta.ka); cl<-exp(tcl+eta.cl); v<-exp(tv+eta.v)
          d/dt(depot)<- -ka*depot; d/dt(center)<-ka*depot-(cl/v)*center
          cp<-center/v; cp~prop(prop.sd)+add(add.sd) })
}

## reference: individual FOCEI
indiv <- theo[,c("ID","TIME","DV","AMT","EVID","CMT")]
t0 <- Sys.time()
ref <- suppressWarnings(suppressMessages(nlmixr2(theo3, indiv, est="focei", control=foceiControl(print=0L))))
th <- ref$theta; om <- ref$omega
refv <- c(ka=exp(th[["tka"]]), cl=exp(th[["tcl"]]), v=exp(th[["tv"]]))
cat(sprintf("FOCEI time %.0fs | ka=%.3f cl=%.3f v=%.3f | om.ka=%.3f om.cl=%.3f om.v=%.3f\n",
    as.numeric(Sys.time()-t0,units="secs"), refv[1],refv[2],refv[3], om[1,1],om[2,2],om[3,3]))

fit_ag <- function(Vmat, tag) {
  t1 <- Sys.time()
  f <- suppressWarnings(suppressMessages(nlmixr2(theo3, admData(), est="adgh",
       control=adghControl(studies=list(t=list(E=E,V=Vmat,n=n,times=times,ev=rxode2::et(amt=DOSE,cmt=1))),
                           maxeval=400L))))
  ex <- f$env$admExtra
  cat(sprintf("%-10s %.0fs | ka=%.3f cl=%.3f v=%.3f  (rel err ka=%.0f%% cl=%.0f%% v=%.0f%%)\n",
      tag, as.numeric(Sys.time()-t1,units="secs"),
      exp(ex$struct[["tka"]]), exp(ex$struct[["tcl"]]), exp(ex$struct[["tv"]]),
      100*abs(exp(ex$struct[["tka"]])/refv[1]-1), 100*abs(exp(ex$struct[["tcl"]])/refv[2]-1),
      100*abs(exp(ex$struct[["tv"]])/refv[3]-1)))
}
fit_ag(diag(diag(V)), "adgh var")
fit_ag(V,             "adgh cov")
