## THREE-WAY, REPLICATED.  AGG vs PARAM-LL vs the IPD gold standard.
##
## Per replicate: simulate individuals -> real FOCEI fits give the two
## "publications" and the pooled gold standard -> then fit the unified model two
## ways from the publications alone:
##
##   AGG       Method 1. Simulate each PUBLICATION's model at covariate nodes,
##             fit the unified model to those blocks. Uses the sources'
##             PREDICTIONS. Study A has no renal term, so its blocks are FLAT in
##             CRCL -- the "absent covariate read as zero effect" trap, live.
##   PARAM-LL  Method 2. Binding function on the published parameter tables.
##
## Scored against IPD-POOL (the target: what the individual data actually
## support) and against TRUTH.
setwd("C:/package/admixr2/.claude/worktrees/feature-covariate-quadrature")
suppressMessages(library(rxode2)); suppressMessages(library(nlmixr2est))
OUT <- "C:/Users/hidde/AppData/Local/Temp/claude/C--package-admixr2/3ff305c7-64a5-4cb0-ba61-1436e2e9b16e/scratchpad"
say <- function(...) { cat(...,"\n"); utils::flush.console() }; T0 <- proc.time()[["elapsed"]]
gh <- function(k){ i <- seq_len(k-1L); J <- matrix(0,k,k)
  J[cbind(i,i+1L)] <- sqrt(i); J[cbind(i+1L,i)] <- sqrt(i)
  e <- eigen(J,symmetric=TRUE); o <- order(e$values); list(x=e$values[o],w=e$vectors[1L,o]^2) }
NE <- 5L; NC <- 4L; GE <- gh(NE); GX <- gh(NC)
TP <- c(lcl=log(4.2), lvc=log(30), lq=log(7.0), lvp=log(40),
        bCL1=.72, bCL2=.45, bVc1=.95, oCL=.26, oVc=.20, sig=.11)
UN <- names(TP); ANM <- c("lcl","lvc","bCL1","bVc1","oCL","oVc","sig")
DA <- list(n=140, dose=500, ii=12, addl=5, obs=c(61,72), mx1=log(82/70), sx1=.20,
           mx2=log(78/100), sx2=.28, times=c(1,12), tau=12)
DB <- list(n=60, dose=500, ii=0, addl=0, obs=c(.25,.5,1,2,4,8,12), mx1=0, sx1=.17,
           mx2=log(96/100), sx2=.22, times=c(.25,.5,1,2,4,8,12), tau=NA)
sim_mod <- rxode2({
  cl <- exp(lcl + bCL1*log(WT/70) + bCL2*log(CRCL/100) + eCL)
  vc <- exp(lvc + bVc1*log(WT/70) + eVC); q <- exp(lq); vp <- exp(lvp)
  d/dt(centr) <- -(cl+q)/vc*centr + q/vp*peri
  d/dt(peri)  <-        q/vc*centr - q/vp*peri
  cp <- centr/vc })
sim_study <- function(d, id0) {
  x1 <- rnorm(d$n,d$mx1,d$sx1); x2 <- rnorm(d$n,d$mx2,d$sx2)
  ev <- if (d$ii>0) et(amt=d$dose, ii=d$ii, addl=d$addl) else et(amt=d$dose)
  ev <- et(ev, d$obs)
  obs <- do.call(rbind, lapply(seq_len(d$n), function(i){
    pr <- data.frame(lcl=TP[["lcl"]],lvc=TP[["lvc"]],lq=TP[["lq"]],lvp=TP[["lvp"]],
      bCL1=TP[["bCL1"]],bCL2=TP[["bCL2"]],bVc1=TP[["bVc1"]],WT=70*exp(x1[i]),
      CRCL=100*exp(x2[i]),eCL=rnorm(1,0,TP[["oCL"]]),eVC=rnorm(1,0,TP[["oVc"]]))
    s <- rxSolve(sim_mod,pr,ev,returnType="data.frame",cores=1,
                 nDisplayProgress=.Machine$integer.max)
    s <- s[s$time %in% d$obs,,drop=FALSE]
    data.frame(ID=id0+i,TIME=s$time,DV=s$cp*(1+rnorm(nrow(s),0,TP[["sig"]])),
               AMT=0,EVID=0,WT=70*exp(x1[i]),CRCL=100*exp(x2[i])) }))
  nd <- if (d$ii>0) d$addl+1L else 1L
  dose <- do.call(rbind,lapply(seq_len(d$n),function(i) data.frame(ID=id0+i,
    TIME=(seq_len(nd)-1)*d$ii,DV=NA_real_,AMT=d$dose,EVID=1,
    WT=70*exp(x1[i]),CRCL=100*exp(x2[i]))))
  df <- rbind(dose,obs); df[order(df$ID,df$TIME,-df$EVID),] }
mA <- function(){ ini({lcl<-log(4);lv<-log(60);bCL1<-.7;bVc1<-.9
    eta.cl~.07;eta.v~.05;prop.err<-.12})
  model({cl<-exp(lcl+bCL1*log(WT/70)+eta.cl); v<-exp(lv+bVc1*log(WT/70)+eta.v)
    d/dt(centr)<--cl/v*centr; cp<-centr/v; cp~prop(prop.err)}) }
mU <- function(){ ini({lcl<-log(4);lvc<-log(35);lq<-log(6);lvp<-log(35)
    bCL1<-.7;bCL2<-.4;bVc1<-.9; eta.cl~.07;eta.vc~.05;prop.err<-.12})
  model({cl<-exp(lcl+bCL1*log(WT/70)+bCL2*log(CRCL/100)+eta.cl)
    vc<-exp(lvc+bVc1*log(WT/70)+eta.vc); q<-exp(lq); vp<-exp(lvp)
    d/dt(centr)<--(cl+q)/vc*centr+q/vp*peri; d/dt(peri)<-q/vc*centr-q/vp*peri
    cp<-centr/vc; cp~prop(prop.err)}) }
grab <- function(f,nm){ fx<-fixef(f); om<-f$omega
  vn <- grep("^eta\\.vc?$",rownames(om),value=TRUE)[1]
  setNames(vapply(nm,function(q) if(q=="oCL") sqrt(om[["eta.cl","eta.cl"]])
    else if(q=="oVc") sqrt(om[[vn,vn]]) else if(q=="sig") fx[["prop.err"]]
    else fx[[q]],0),nm) }
ff <- function(m,d) suppressWarnings(suppressMessages(
  nlmixr2(m,d,est="focei",control=foceiControl(print=0L))))
seof <- function(f,nm,fb){ pf<-f$parFixedDf
  v <- vapply(nm,function(q){ r<-if(q=="sig")"prop.err" else q
    s<-pf[rownames(pf)==r,"SE"]; if(length(s)&&is.finite(s)) s else NA_real_},0)
  v[is.na(v)] <- fb; v }
## ---- analytic harness (shared by AGG and PARAM-LL) -------------------------
ssf <- function(e,t) if (is.na(t)) 1 else 1/(1-exp(-e*t))
mkgrid <- function(d){ cg<-expand.grid(i=seq_len(NC),j=seq_len(NC),KEEP.OUT.ATTRS=FALSE)
  a1<-d$mx1+d$sx1*GX$x[cg$i]; a2<-d$mx2+d$sx2*GX$x[cg$j]
  wk<-GX$w[cg$i]*GX$w[cg$j]; wk<-wk/sum(wk)
  eg<-expand.grid(p=seq_len(NE),q=seq_len(NE),KEEP.OUT.ATTRS=FALSE)
  we<-GE$w[eg$p]*GE$w[eg$q]; we<-we/sum(we); K<-nrow(cg); M<-nrow(eg)
  ix<-rep(seq_len(K),each=M)
  list(K=K,M=M,wk=wk,we=we,x1=a1[ix],x2=a2[ix],z1=GE$x[eg$p][rep(seq_len(M),K)],
       z2=GE$x[eg$q][rep(seq_len(M),K)],n=d$n,d=d) }
mall <- function(p,G,form){ d<-G$d
  b2 <- if("bCL2"%in%names(p)) p[["bCL2"]] else 0
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
fit_blocks <- function(BL,GL,FM,nm,start,mx=2500L){
  ob<-function(v){ p<-setNames(v,nm); tot<-0
    for (s in seq_along(BL)) { pr<-mall(p,GL[[s]],FM[s])
      tot<-tot+sum(vapply(seq_along(BL[[s]]),function(k)
        nl1(BL[[s]][[k]],pr[[k]],GL[[s]]$wk[k]*GL[[s]]$n),0)) }
    if(!is.finite(tot)) 1e12 else tot }
  o<-optim(start[nm],ob,method="Nelder-Mead",control=list(maxit=mx,reltol=1e-11))
  p<-setNames(o$par,nm); po<-intersect(nm,c("oCL","oVc","sig")); p[po]<-abs(p[po]); p }

R <- 6L; RES <- list()
for (rep in seq_len(R)) {
  set.seed(1000+rep); tr <- proc.time()[["elapsed"]]
  dA<-sim_study(DA,0L); dB<-sim_study(DB,1000L)
  fA<-ff(mA,dA); fB<-ff(mU,dB); fP<-ff(mU,rbind(dA,dB))
  thA<-grab(fA,c("lcl","lv","bCL1","bVc1","oCL","oVc","sig")); names(thA)[2]<-"lvc"
  thB<-grab(fB,UN); GOLD<-grab(fP,UN)
  seA<-c(seof(fA,c("lcl","lv","bCL1","bVc1","sig"),.15),oCL=.04,oVc=.05)
  names(seA)[2]<-"lvc"; seA<-seA[ANM]
  seB<-c(seof(fB,setdiff(UN,c("oCL","oVc")),.15),oCL=.03,oVc=.04)[UN]
  GA<-mkgrid(DA); GB<-mkgrid(DB)
  ## AGG: blocks from each PUBLICATION's own model
  blA<-mall(thA,GA,"1"); blB<-mall(thB,GB,"2")
  AGG<-fit_blocks(list(blA,blB),list(GA,GB),c("2","2"),UN,thB)
  ## PARAM-LL
  WARM<-new.env(); WARM$A<-thA[ANM]; WARM$B<-thB[UN]
  g_of<-function(phi,G,form,nm,key){ obs<-mall(phi,G,"2")
    ob<-function(v){ p<-setNames(v,nm); pr<-mall(p,G,form)
      sum(vapply(seq_len(G$K),function(k) nl1(obs[[k]],pr[[k]],G$wk[k]*G$n),0)) }
    o<-optim(WARM[[key]],ob,method="Nelder-Mead",control=list(maxit=1200,reltol=1e-10))
    WARM[[key]]<-o$par; p<-setNames(o$par,nm)
    po<-intersect(nm,c("oCL","oVc","sig")); p[po]<-abs(p[po]); p }
  res<-function(phi) c((thA[ANM]-g_of(phi,GA,"1",ANM,"A"))/seA,
                       (thB[UN]-g_of(phi,GB,"2",UN,"B"))/seB)
  FD<-setNames(rep(6e-3,length(UN)),UN); FD[c("bCL1","bCL2","bVc1")]<-1.2e-2
  phi<-thB[UN]; bs<-Inf; bp<-phi
  for (it in 1:4) { r0<-res(phi); if(sum(r0^2)<bs){bs<-sum(r0^2);bp<-phi}
    J<-vapply(UN,function(q){ph<-phi;ph[[q]]<-ph[[q]]+FD[[q]];(res(ph)-r0)/FD[[q]]},numeric(length(r0)))
    phi<-phi-0.6*as.numeric(solve(crossprod(J)+diag(1e-7,length(UN)),crossprod(J,r0))) }
  PLL<-bp
  RES[[rep]]<-list(GOLD=GOLD,AGG=AGG,PLL=PLL,thA=thA,thB=thB)
  say(sprintf("[%4.0fs] rep %d  |gold-truth| %.4f  AGG %.4f  PLL %.4f  (vs gold: AGG %.4f PLL %.4f)  %.0fs",
    proc.time()[["elapsed"]]-T0, rep, sqrt(mean((GOLD-TP)^2)), sqrt(mean((AGG-TP)^2)),
    sqrt(mean((PLL-TP)^2)), sqrt(mean((AGG-GOLD)^2)), sqrt(mean((PLL-GOLD)^2)),
    proc.time()[["elapsed"]]-tr))
  saveRDS(list(RES=RES,TP=TP), file.path(OUT,"threeway.rds"))
}
M <- function(f) vapply(RES,f,0)
say(sprintf("\n%d replicates\n", length(RES)))
say(sprintf("%-22s %10s %10s", "", "mean RMSE", "sd"))
for (nm in c("GOLD","AGG","PLL")) say(sprintf("%-22s %10.4f %10.4f",
  paste(nm,"vs TRUTH"), mean(M(function(x) sqrt(mean((x[[nm]]-TP)^2)))),
  sd(M(function(x) sqrt(mean((x[[nm]]-TP)^2))))))
for (nm in c("AGG","PLL")) say(sprintf("%-22s %10.4f %10.4f",
  paste(nm,"vs GOLD"), mean(M(function(x) sqrt(mean((x[[nm]]-x$GOLD)^2)))),
  sd(M(function(x) sqrt(mean((x[[nm]]-x$GOLD)^2))))))
say(sprintf("\nbCL2 (renal; ABSENT from A's model)  truth %.3f", TP[["bCL2"]]))
for (nm in c("GOLD","AGG","PLL","thB")) say(sprintf("  %-6s mean %.3f  bias %+.3f  sd %.3f",
  nm, mean(M(function(x) x[[nm]][["bCL2"]])), mean(M(function(x) x[[nm]][["bCL2"]]))-TP[["bCL2"]],
  sd(M(function(x) x[[nm]][["bCL2"]]))))
say(sprintf("total %.0f s", proc.time()[["elapsed"]]-T0))
