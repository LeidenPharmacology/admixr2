## Can the AGGREGATE method be de-biased by MARGINALISING over the covariates a
## source does not report, instead of stratifying on them?
##
## The failure measured in the overnight study is a construction artefact. `mkgrid`
## builds a PRODUCT grid over both covariates for every source, so a publication
## that never fitted x2 is still evaluated on a grid that varies x2 -- and answers
## "no change" at every node. The unified model is then fitted to that flat answer
## and reads it as bCL2 = 0. The contrast was manufactured by us, not published.
##
## The fix is the matching principle applied PER COVARIATE rather than per study:
##
##   stratify   on the covariates the source's own model conditions on
##   marginalise over the covariates it does not, using the source's reported
##              covariate distribution -- on BOTH sides of the likelihood
##
## Then the source contributes no mean contrast in x2 (correctly -- it has none to
## give), but it still contributes its TOTAL VARIANCE budget, and that budget is
## informative:  oCL_A^2 = oCL_U^2 + bCL2^2 Var_A(x2).  With oCL_U pinned by the
## source that did fit x2, bCL2 is recoverable from the variance channel alone.
##
## Runs off the STORED thA/thB of the overnight study, so AGGm is paired replicate
## by replicate with the AGG already reported -- no refitting, same publications.
setwd("C:/package/admixr2/.claude/worktrees/feature-covariate-quadrature")
OUT <- "C:/Users/hidde/AppData/Local/Temp/claude/C--package-admixr2/3ff305c7-64a5-4cb0-ba61-1436e2e9b16e/scratchpad"
say <- function(...) { cat(...,"\n"); utils::flush.console() }; T0 <- proc.time()[["elapsed"]]
gh <- function(k){ i<-seq_len(k-1L); J<-matrix(0,k,k)
  J[cbind(i,i+1L)]<-sqrt(i); J[cbind(i+1L,i)]<-sqrt(i)
  e<-eigen(J,symmetric=TRUE); o<-order(e$values); list(x=e$values[o],w=e$vectors[1L,o]^2) }
NE<-5L; NC<-4L; GE<-gh(NE); GX<-gh(NC)
TP <- c(lcl=log(4.2),lvc=log(30),lq=log(7.0),lvp=log(40),
        bCL1=.72,bCL2=.45,bVc1=.95,oCL=.26,oVc=.20,sig=.11)
UN <- names(TP)
DA <- list(n=140,dose=500,ii=12,addl=5,obs=c(61,72),mx1=log(82/70),sx1=.20,
           mx2=log(78/100),sx2=.28,times=c(1,12),tau=12)
DB <- list(n=60,dose=500,ii=0,addl=0,obs=c(.25,.5,1,2,4,8,12),mx1=0,sx1=.17,
           mx2=log(96/100),sx2=.22,times=c(.25,.5,1,2,4,8,12),tau=NA)
CELLS <- list(baseline=list(f="2",x2=TRUE), struct=list(f="1",x2=TRUE),
              omit=list(f="2",x2=FALSE),    both=list(f="1",x2=FALSE))
pnames <- function(form,x2) c("lcl","lvc", if (form=="2") c("lq","lvp"),
  "bCL1", if (x2) "bCL2", "bVc1","oCL","oVc","sig")

## ---- grid with a per-covariate stratify/marginalise split -------------------
## strat[i] TRUE  -> covariate i becomes an OUTER node dimension (a separate block)
## strat[i] FALSE -> covariate i joins the INNER expectation, beside the etas, so
##                   it is integrated out of E and folded into V by total variance.
## Row layout r = (k-1)*M + m is identical to the original grid, so `mall` is
## reused unchanged -- the two constructions differ only in this function.
mkgrid <- function(d, strat=c(TRUE,TRUE)){
  ax <- list(list(m=d$mx1,s=d$sx1), list(m=d$mx2,s=d$sx2))
  so <- which(strat); si <- which(!strat)
  if (length(so)) {
    og <- expand.grid(lapply(so, function(i) seq_len(NC)), KEEP.OUT.ATTRS=FALSE)
    wk <- Reduce(`*`, lapply(seq_along(so), function(j) GX$w[og[[j]]]))
  } else { og <- data.frame(a=1L); wk <- 1 }
  wk <- wk/sum(wk); K <- nrow(og)
  ig <- expand.grid(c(list(seq_len(NE), seq_len(NE)),
                      lapply(si, function(i) seq_len(NC))), KEEP.OUT.ATTRS=FALSE)
  we <- GE$w[ig[[1]]]*GE$w[ig[[2]]]
  for (j in seq_along(si)) we <- we*GX$w[ig[[2+j]]]
  we <- we/sum(we); M <- nrow(ig)
  ix <- rep(seq_len(K), each=M); jx <- rep(seq_len(M), K)
  x <- lapply(1:2, function(i)
    if (i %in% so) { j<-match(i,so); (ax[[i]]$m + ax[[i]]$s*GX$x[og[[j]]])[ix] }
    else           { j<-match(i,si); (ax[[i]]$m + ax[[i]]$s*GX$x[ig[[2+j]]])[jx] })
  list(K=K,M=M,wk=wk,we=we,x1=x[[1]],x2=x[[2]],
       z1=GE$x[ig[[1]]][jx], z2=GE$x[ig[[2]]][jx], n=d$n, d=d) }

ssf <- function(e,t) if (is.na(t)) 1 else 1/(1-exp(-e*t))
mall <- function(p,G,form){ d<-G$d
  b2 <- if ("bCL2" %in% names(p)) p[["bCL2"]] else 0
  cl<-exp(p[["lcl"]]+p[["bCL1"]]*G$x1+b2*G$x2+abs(p[["oCL"]])*G$z1)
  vv<-exp(p[["lvc"]]+p[["bVc1"]]*G$x1+abs(p[["oVc"]])*G$z2)
  f <- if(form=="2"){ q<-exp(p[["lq"]]);vp<-exp(p[["lvp"]])
      k10<-cl/vv;k12<-q/vv;k21<-q/vp;s<-k10+k12+k21
      dd<-sqrt(pmax(s^2-4*k10*k21,0));al<-(s+dd)/2;be<-(s-dd)/2
      A<-(d$dose/vv)*(al-k21)/(al-be);Bc<-(d$dose/vv)*(k21-be)/(al-be)
      (A*ssf(al,d$tau))*exp(-outer(al,d$times))+(Bc*ssf(be,d$tau))*exp(-outer(be,d$times))
    } else { k<-cl/vv; outer((d$dose/vv)*ssf(k,d$tau),d$times,function(a,b)a)*exp(-outer(k,d$times)) }
  s2<-abs(p[["sig"]])^2
  lapply(seq_len(G$K),function(k){ r<-((k-1L)*G$M+1L):(k*G$M); fk<-f[r,,drop=FALSE]
    mu<-as.numeric(crossprod(G$we,fk)); fc<-sweep(fk,2L,mu)
    V<-crossprod(fc,fc*G$we); diag(V)<-diag(V)+s2*as.numeric(crossprod(G$we,fk^2))
    list(E=mu,V=V)}) }
nl1 <- function(o,pr,n){ ch<-tryCatch(chol(pr$V),error=function(e) NULL)
  if(is.null(ch)) return(1e12); iv<-chol2inv(ch); r<-o$E-pr$E
  n*(2*sum(log(diag(ch)))+sum(iv*o$V)+as.numeric(t(r)%*%iv%*%r)) }
fitb <- function(BL,GL,FM,nm,start,mx=2500L){
  ob<-function(v){ p<-setNames(v,nm); tot<-0
    for (s in seq_along(BL)){ pr<-mall(p,GL[[s]],FM[s])
      tot<-tot+sum(vapply(seq_along(BL[[s]]),function(k)
        nl1(BL[[s]][[k]],pr[[k]],GL[[s]]$wk[k]*GL[[s]]$n),0)) }
    if(!is.finite(tot)) 1e12 else tot }
  o<-optim(start[nm],ob,method="Nelder-Mead",control=list(maxit=mx,reltol=1e-11))
  p<-setNames(o$par,nm); po<-intersect(nm,c("oCL","oVc","sig")); p[po]<-abs(p[po]); p }

## ---- replay every stored replicate with the marginal construction -----------
d <- readRDS(file.path(OUT,"overnight.rds")); RES <- d$RES
say(sprintf("replaying %d stored replicates\n", length(RES)))
GB <- mkgrid(DB, c(TRUE,TRUE))                  # B always fits both covariates
OUTL <- list()
for (i in seq_along(RES)) { r <- RES[[i]]; cl <- CELLS[[r$cell]]
  ANM <- pnames(cl$f, cl$x2)
  ## A stratifies only on the covariates its OWN publication conditions on
  strA <- c(TRUE, isTRUE(cl$x2))
  GA <- mkgrid(DA, strA)
  AGGm <- try(fitb(list(mall(r$thA,GA,cl$f), mall(r$thB,GB,"2")),
                   list(GA,GB), c("2","2"), UN, r$thB), silent=TRUE)
  if (inherits(AGGm,"try-error")) { say(sprintf("  %-9s rep %2d FAILED", r$cell, r$rep)); next }
  OUTL[[length(OUTL)+1]] <- c(r, list(AGGm=AGGm, strat_x2=strA[2]))
  if (i %% 10 == 0) say(sprintf("[%4.0fs] %3d/%d", proc.time()[["elapsed"]]-T0, i, length(RES)))
}
saveRDS(list(RES=OUTL,TP=TP), file.path(OUT,"aggregate-marginal.rds"))

## ---- report -----------------------------------------------------------------
cells <- vapply(OUTL, `[[`, "", "cell"); rmse <- function(a,b) sqrt(mean((a-b)^2))
ARMS <- c(AGG="AGG", AGGm="AGGm", PLLs="PLLs")
say("\n=========== deviation from the INDIVIDUAL-DATA fit (primary) ===========")
say(sprintf("%-10s %6s %12s %12s %12s", "cell","reps","AGG","AGG-marginal","PLL"))
for (cn in names(CELLS)) { ix <- which(cells==cn); if (!length(ix)) next
  v <- lapply(ARMS, function(a) vapply(ix, function(i) rmse(OUTL[[i]][[a]],OUTL[[i]]$GOLD), 0))
  say(sprintf("%-10s %6d %12.4f %12.4f %12.4f", cn, length(ix),
      mean(v$AGG), mean(v$AGGm), mean(v$PLLs))) }

say("\n=========== deviation from the simulation truth ===========")
say(sprintf("%-10s %6s %12s %12s %12s %12s", "cell","reps","IPD (gold)","AGG","AGG-marginal","PLL"))
for (cn in names(CELLS)) { ix <- which(cells==cn); if (!length(ix)) next
  g <- vapply(ix, function(i) rmse(OUTL[[i]]$GOLD,TP), 0)
  v <- lapply(ARMS, function(a) vapply(ix, function(i) rmse(OUTL[[i]][[a]],TP), 0))
  say(sprintf("%-10s %6d %12.4f %12.4f %12.4f %12.4f", cn, length(ix),
      mean(g), mean(v$AGG), mean(v$AGGm), mean(v$PLLs))) }

say(sprintf("\n=========== bCL2 (truth %.3f) -- the covariate A omits ===========", TP[["bCL2"]]))
say(sprintf("%-10s %6s %14s %14s %14s %14s", "cell","reps","IPD (gold)","AGG","AGG-marginal","PLL"))
for (cn in names(CELLS)) { ix <- which(cells==cn); if (!length(ix)) next
  f <- function(a){ v<-vapply(ix,function(i) OUTL[[i]][[a]][["bCL2"]],0)
    sprintf("%6.3f%+7.3f", mean(v), mean(v)-TP[["bCL2"]]) }
  say(sprintf("%-10s %6d %14s %14s %14s %14s", cn, length(ix),
      f("GOLD"), f("AGG"), f("AGGm"), f("PLLs"))) }
say("\n(mean, and bias against truth)")

say("\n=========== oCL (truth 0.260) -- the variance channel ===========")
say(sprintf("%-10s %14s %14s %14s %14s","cell","IPD (gold)","AGG","AGG-marginal","A's published"))
for (cn in names(CELLS)) { ix <- which(cells==cn); if (!length(ix)) next
  f <- function(a) sprintf("%14.3f", mean(vapply(ix,function(i) OUTL[[i]][[a]][["oCL"]],0)))
  say(sprintf("%-10s %14s %14s %14s %14.3f", cn, f("GOLD"), f("AGG"), f("AGGm"),
      mean(vapply(ix,function(i) OUTL[[i]]$thA[["oCL"]],0)))) }
say(sprintf("\ntotal %.0f s", proc.time()[["elapsed"]]-T0))
