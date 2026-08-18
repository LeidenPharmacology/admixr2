## OVERNIGHT SIMULATION STUDY
## Isolates the TWO mechanisms this session identified, and their interaction.
##
## TRUTH  2-cmt, two covariates: x1 = log(WT/70), x2 = log(CRCL/100)
## STUDY B  always 2-cmt with BOTH covariates, dense sampling  (the informative source)
## STUDY A  varies by cell -- it is the source that can go wrong:
##
##   cell            A's structure     A fits x2?     mechanism isolated
##   ------------------------------------------------------------------------
##   baseline        2-cmt             yes            none: methods should agree
##   struct          1-cmt             yes            structural mismatch only
##   omit            2-cmt             NO             absent covariate only
##   both            1-cmt             NO             interaction
##
## METHODS  IPD-POOL (gold, pooled individual records) | AGG (nodes from each
##          publication's model) | PLL-sens (binding function, computed Fisher
##          weight) | PLL-rse (binding function, reported SEs)
##
## Checkpoints after EVERY replicate to overnight.rds, so a partial run is usable.
setwd("C:/package/admixr2/.claude/worktrees/feature-covariate-quadrature")
suppressMessages(library(rxode2)); suppressMessages(library(nlmixr2est))
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
sim_mod <- rxode2({
  cl <- exp(lcl + bCL1*log(WT/70) + bCL2*log(CRCL/100) + eCL)
  vc <- exp(lvc + bVc1*log(WT/70) + eVC); q <- exp(lq); vp <- exp(lvp)
  d/dt(centr) <- -(cl+q)/vc*centr + q/vp*peri
  d/dt(peri)  <-        q/vc*centr - q/vp*peri
  cp <- centr/vc })
sim_study <- function(d,id0){
  x1<-rnorm(d$n,d$mx1,d$sx1); x2<-rnorm(d$n,d$mx2,d$sx2)
  ev <- if (d$ii>0) et(amt=d$dose,ii=d$ii,addl=d$addl) else et(amt=d$dose)
  ev <- et(ev,d$obs)
  obs <- do.call(rbind,lapply(seq_len(d$n),function(i){
    pr<-data.frame(lcl=TP[["lcl"]],lvc=TP[["lvc"]],lq=TP[["lq"]],lvp=TP[["lvp"]],
      bCL1=TP[["bCL1"]],bCL2=TP[["bCL2"]],bVc1=TP[["bVc1"]],WT=70*exp(x1[i]),
      CRCL=100*exp(x2[i]),eCL=rnorm(1,0,TP[["oCL"]]),eVC=rnorm(1,0,TP[["oVc"]]))
    s<-rxSolve(sim_mod,pr,ev,returnType="data.frame",cores=1,
               nDisplayProgress=.Machine$integer.max)
    s<-s[s$time %in% d$obs,,drop=FALSE]
    data.frame(ID=id0+i,TIME=s$time,DV=s$cp*(1+rnorm(nrow(s),0,TP[["sig"]])),
               AMT=0,EVID=0,WT=70*exp(x1[i]),CRCL=100*exp(x2[i])) }))
  nd <- if (d$ii>0) d$addl+1L else 1L
  dose<-do.call(rbind,lapply(seq_len(d$n),function(i) data.frame(ID=id0+i,
    TIME=(seq_len(nd)-1)*d$ii,DV=NA_real_,AMT=d$dose,EVID=1,
    WT=70*exp(x1[i]),CRCL=100*exp(x2[i]))))
  df<-rbind(dose,obs); df[order(df$ID,df$TIME,-df$EVID),] }
## model builders: form (1/2) x covariate set (with/without x2)
mk <- function(form,x2){
  cl <- if (x2) "cl<-exp(lcl+bCL1*log(WT/70)+bCL2*log(CRCL/100)+eta.cl)" else
                "cl<-exp(lcl+bCL1*log(WT/70)+eta.cl)"
  ini2 <- if (x2) "bCL1<-.7;bCL2<-.4;bVc1<-.9" else "bCL1<-.7;bVc1<-.9"
  body <- if (form=="2")
    sprintf("%s; vc<-exp(lvc+bVc1*log(WT/70)+eta.vc); q<-exp(lq); vp<-exp(lvp)
      d/dt(centr)<--(cl+q)/vc*centr+q/vp*peri; d/dt(peri)<-q/vc*centr-q/vp*peri
      cp<-centr/vc; cp~prop(prop.err)", cl)
  else sprintf("%s; vc<-exp(lvc+bVc1*log(WT/70)+eta.vc)
      d/dt(centr)<--cl/vc*centr; cp<-centr/vc; cp~prop(prop.err)", cl)
  ini1 <- if (form=="2") "lcl<-log(4);lvc<-log(35);lq<-log(6);lvp<-log(35)" else
                         "lcl<-log(4);lvc<-log(45)"
  eval(parse(text=sprintf('function(){ini({%s;%s;eta.cl~.07;eta.vc~.05;prop.err<-.12})
    model({%s})}', ini1, ini2, body))) }
pnames <- function(form,x2) c("lcl","lvc", if (form=="2") c("lq","lvp"),
  "bCL1", if (x2) "bCL2", "bVc1","oCL","oVc","sig")
grab <- function(f,nm){ fx<-fixef(f); om<-f$omega
  setNames(vapply(nm,function(q) if(q=="oCL") sqrt(om[["eta.cl","eta.cl"]])
    else if(q=="oVc") sqrt(om[["eta.vc","eta.vc"]]) else if(q=="sig") fx[["prop.err"]]
    else fx[[q]],0),nm) }
ff <- function(m,d) suppressWarnings(suppressMessages(
  nlmixr2(m,d,est="focei",control=foceiControl(print=0L))))
seof <- function(f,nm){ pf<-f$parFixedDf
  v<-vapply(nm,function(q){ r<-if(q=="sig")"prop.err" else q
    s<-pf[rownames(pf)==r,"SE"]; if(length(s)&&is.finite(s)) s else NA_real_},0)
  v[is.na(v)]<-0.15; v[nm %in% c("oCL","oVc")]<-0.04; v }
## ---- analytic harness (binding function + AGG) -----------------------------
ssf <- function(e,t) if (is.na(t)) 1 else 1/(1-exp(-e*t))
mkgrid <- function(d){ cg<-expand.grid(i=seq_len(NC),j=seq_len(NC),KEEP.OUT.ATTRS=FALSE)
  wk<-GX$w[cg$i]*GX$w[cg$j]; wk<-wk/sum(wk)
  eg<-expand.grid(p=seq_len(NE),q=seq_len(NE),KEEP.OUT.ATTRS=FALSE)
  we<-GE$w[eg$p]*GE$w[eg$q]; we<-we/sum(we); K<-nrow(cg); M<-nrow(eg)
  ix<-rep(seq_len(K),each=M)
  list(K=K,M=M,wk=wk,we=we,x1=(d$mx1+d$sx1*GX$x[cg$i])[ix],
       x2=(d$mx2+d$sx2*GX$x[cg$j])[ix],z1=GE$x[eg$p][rep(seq_len(M),K)],
       z2=GE$x[eg$q][rep(seq_len(M),K)],n=d$n,d=d) }
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
## full Gaussian information for a source, from its own design
infoW <- function(th,G,form){
  b<-mall(th,G,form); E0<-do.call(c,lapply(b,`[[`,"E"))
  Vs<-lapply(b,`[[`,"V"); Vi<-lapply(Vs,function(V) chol2inv(chol(V)))
  nm<-names(th); K<-length(nm); dE<-list(); dV<-list()
  for (q in nm){ p<-th; h<-max(abs(p[[q]]),0.1)*1e-2; p[[q]]<-p[[q]]+h
    b1<-mall(p,G,form)
    dE[[q]]<-(do.call(c,lapply(b1,`[[`,"E"))-E0)/h
    dV[[q]]<-lapply(seq_along(Vs),function(k) (b1[[k]]$V-Vs[[k]])/h) }
  M<-matrix(0,K,K,dimnames=list(nm,nm)); nt<-length(b[[1]]$E)
  for (i in seq_len(K)) for (j in i:K){ tot<-0
    for (k in seq_along(Vs)){ idx<-((k-1)*nt+1):(k*nt); nk<-G$wk[k]*G$n
      tot<-tot+nk*as.numeric(t(dE[[i]][idx])%*%Vi[[k]]%*%dE[[j]][idx]) +
        (nk/2)*sum((Vi[[k]]%*%dV[[i]][[k]])*t(Vi[[k]]%*%dV[[j]][[k]])) }
    M[i,j]<-M[j,i]<-tot }
  M }

R <- 25L; RES <- list()
for (cn in names(CELLS)) { cl <- CELLS[[cn]]
  ANM <- pnames(cl$f, cl$x2); mA <- mk(cl$f, cl$x2); mU <- mk("2", TRUE)
  for (rep in seq_len(R)) {
    set.seed(5000 + 100*match(cn,names(CELLS)) + rep); tr <- proc.time()[["elapsed"]]
    ok <- try({
      dA<-sim_study(DA,0L); dB<-sim_study(DB,1000L)
      fA<-ff(mA,dA); fB<-ff(mU,dB); fP<-ff(mU,rbind(dA,dB))
      thA<-grab(fA,ANM); thB<-grab(fB,UN); GOLD<-grab(fP,UN)
      seA<-seof(fA,ANM); seB<-seof(fB,UN)
      GA<-mkgrid(DA); GB<-mkgrid(DB)
      AGG<-fitb(list(mall(thA,GA,cl$f),mall(thB,GB,"2")),list(GA,GB),c("2","2"),UN,thB)
      WARM<-new.env(); WARM$A<-thA; WARM$B<-thB
      g_of<-function(phi,G,form,nm,key){ obs<-mall(phi,G,"2")
        ob<-function(v){ p<-setNames(v,nm); pr<-mall(p,G,form)
          sum(vapply(seq_len(G$K),function(k) nl1(obs[[k]],pr[[k]],G$wk[k]*G$n),0)) }
        o<-optim(WARM[[key]],ob,method="Nelder-Mead",control=list(maxit=1200,reltol=1e-10))
        WARM[[key]]<-o$par; p<-setNames(o$par,nm)
        po<-intersect(nm,c("oCL","oVc","sig")); p[po]<-abs(p[po]); p }
      WA<-infoW(thA,GA,cl$f); WB<-infoW(thB,GB,"2")
      raw<-function(phi) c(thA-g_of(phi,GA,cl$f,ANM,"A"), thB-g_of(phi,GB,"2",UN,"B"))
      FD<-setNames(rep(6e-3,length(UN)),UN); FD[c("bCL1","bCL2","bVc1")]<-1.2e-2
      r0<-raw(thB); J<-vapply(UN,function(q){ph<-thB;ph[[q]]<-ph[[q]]+FD[[q]]
        (raw(ph)-r0)/FD[[q]]},numeric(length(r0)))
      nA<-length(ANM); W1<-matrix(0,length(r0),length(r0))
      W1[1:nA,1:nA]<-WA; W1[(nA+1):length(r0),(nA+1):length(r0)]<-WB
      W2<-diag(c(1/seA^2,1/seB^2))
      solveW<-function(W){ p<-thB
        for (k in 1:3){ r<-if(k==1) r0 else raw(p)
          d<-tryCatch(solve(t(J)%*%W%*%J+diag(1e-7,ncol(J)),t(J)%*%W%*%r),
                      error=function(e) rep(0,ncol(J)))
          p<-p-0.6*as.numeric(d) }; p }
      list(GOLD=GOLD,AGG=AGG,PLLs=solveW(W1),PLLr=solveW(W2),thA=thA,thB=thB)
    }, silent=TRUE)
    if (!inherits(ok,"try-error")) {
      RES[[length(RES)+1]] <- c(ok, list(cell=cn, rep=rep))
      say(sprintf("[%5.0fs] %-9s rep %2d | gold %.4f AGG %.4f PLLs %.4f PLLr %.4f | bCL2 g %.3f A %.3f Ps %.3f (%.0fs)",
        proc.time()[["elapsed"]]-T0, cn, rep,
        sqrt(mean((ok$GOLD-TP)^2)), sqrt(mean((ok$AGG-TP)^2)),
        sqrt(mean((ok$PLLs-TP)^2)), sqrt(mean((ok$PLLr-TP)^2)),
        ok$GOLD[["bCL2"]], ok$AGG[["bCL2"]], ok$PLLs[["bCL2"]],
        proc.time()[["elapsed"]]-tr))
    } else say(sprintf("[%5.0fs] %-9s rep %2d FAILED", proc.time()[["elapsed"]]-T0, cn, rep))
    saveRDS(list(RES=RES,TP=TP,CELLS=CELLS), file.path(OUT,"overnight.rds"))
  }
}
say(sprintf("DONE %d results in %.0f s", length(RES), proc.time()[["elapsed"]]-T0))
