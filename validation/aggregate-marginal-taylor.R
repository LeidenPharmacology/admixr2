## A THIRD arm for the marginal-aggregate study: replace the QUADRATURE over the
## marginalised covariate directions with a second-order TAYLOR expansion of the
## MOMENTS (`taylor-corrected.R`'s `mom_taylor`, lifted into this design).
##
## The split of covariates into stratified vs marginalised is IDENTICAL to AGGm --
## stratify on what the source's own model conditions on, marginalise the rest.
## Only the numerical treatment of the marginalised directions changes:
##
##   AGGm  the marginalised covariate joins the inner expectation as NC GH nodes,
##         so the block moments are exact quadrature over (eta, a).
##   AGGt  the block moments are built from THREE evaluations of the eta-only
##         conditional moments, at a = mu, mu+h, mu-h with h = 0.5 sd:
##
##            E ~= E(mu) + (sd^2/2) E''(mu)
##            V ~= V(mu) + (sd^2/2) V''(mu) + sd^2 dE(mu) dE(mu)^T
##
##         The rank-one term sd^2 dE dE^T is the covariate-induced between-subject
##         covariance -- it is the whole variance channel through which bCL2 is
##         recovered from a source that never fitted x2. Dropping it is a known bug
##         class here, so an AGGt-norank1 diagnostic arm is run to price it.
##
## V(a) is the FULL conditional block covariance including the proportional
## residual diag(sig^2 E_eta[f^2 | a]), which is covariate-dependent, so it is
## expanded along with the rest -- the law of total variance applies to the whole
## conditional block, not just to the structural part.
##
## Where a covariate is STRATIFIED, AGGt must reduce to AGGm exactly. That is
## structural here: with no marginalised direction `mkgridT` emits ONE grid, which
## `.mkg` builds through the same code path as `mkgrid`, so the objective is
## bit-identical and the paired difference is exactly 0. Asserted, not assumed.
##
## Runs off the STORED thA/thB of the overnight study, so all three arms are paired
## replicate by replicate. AGG and AGGm are REFITTED here (not read back) so that
## the wall-clock comparison is like for like; both are checked against their
## stored values.
setwd("C:/package/admixr2/.claude/worktrees/feature-covariate-quadrature")
OUT <- "C:/Users/hidde/AppData/Local/Temp/claude/C--package-admixr2/3ff305c7-64a5-4cb0-ba61-1436e2e9b16e/scratchpad"
say <- function(...) { cat(...,"\n"); utils::flush.console() }; T0 <- proc.time()[["elapsed"]]
gh <- function(k){ i<-seq_len(k-1L); J<-matrix(0,k,k)
  J[cbind(i,i+1L)]<-sqrt(i); J[cbind(i+1L,i)]<-sqrt(i)
  e<-eigen(J,symmetric=TRUE); o<-order(e$values); list(x=e$values[o],w=e$vectors[1L,o]^2) }
NE<-5L; NC<-4L; GE<-gh(NE); GX<-gh(NC); GREF<-gh(21L)   # GREF: reference quadrature
HFRAC <- 0.5
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

## ---- grid construction ------------------------------------------------------
## so = stratified covariates -> OUTER node dimensions (one block each)
## si = marginalised covariates -> either INNER GH nodes (fixed=NULL, the AGGm
##      construction, optionally on a finer node set GQ) or PINNED at the absolute
##      values in `fixed` (the AGGt construction: one grid per Taylor point).
## With si empty both branches are dead and the result is bit-identical to the
## original `mkgrid` of aggregate-marginal.R -- verified by `identical()`.
.mkg <- function(d, so, si, fixed=NULL, GQ=GX){
  ax <- list(list(m=d$mx1,s=d$sx1), list(m=d$mx2,s=d$sx2))
  if (length(so)) {
    og <- expand.grid(lapply(so, function(i) seq_len(NC)), KEEP.OUT.ATTRS=FALSE)
    wk <- Reduce(`*`, lapply(seq_along(so), function(j) GX$w[og[[j]]]))
  } else { og <- data.frame(a=1L); wk <- 1 }
  wk <- wk/sum(wk); K <- nrow(og)
  sq <- if (is.null(fixed)) si else integer(0)
  nq <- length(GQ$x)
  ig <- expand.grid(c(list(seq_len(NE), seq_len(NE)),
                      lapply(sq, function(i) seq_len(nq))), KEEP.OUT.ATTRS=FALSE)
  we <- GE$w[ig[[1]]]*GE$w[ig[[2]]]
  for (j in seq_along(sq)) we <- we*GQ$w[ig[[2+j]]]
  we <- we/sum(we); M <- nrow(ig)
  ix <- rep(seq_len(K), each=M); jx <- rep(seq_len(M), K)
  x <- lapply(1:2, function(i)
    if (i %in% so) { j<-match(i,so); (ax[[i]]$m + ax[[i]]$s*GX$x[og[[j]]])[ix] }
    else if (is.null(fixed)) { j<-match(i,si); (ax[[i]]$m + ax[[i]]$s*GQ$x[ig[[2+j]]])[jx] }
    else { j<-match(i,si); rep(fixed[j], K*M) })
  list(K=K,M=M,wk=wk,we=we,x1=x[[1]],x2=x[[2]],
       z1=GE$x[ig[[1]]][jx], z2=GE$x[ig[[2]]][jx], n=d$n, d=d) }
mkgrid <- function(d, strat=c(TRUE,TRUE), GQ=GX) .mkg(d, which(strat), which(!strat), NULL, GQ)

## Taylor grid bundle: grid 1 is the centre, then (+h, -h) per marginalised covariate.
mkgridT <- function(d, strat=c(TRUE,TRUE), hfrac=HFRAC){
  so <- which(strat); si <- which(!strat)
  ax <- list(list(m=d$mx1,s=d$sx1), list(m=d$mx2,s=d$sx2))
  mu <- vapply(si, function(i) ax[[i]]$m, 0); sd <- vapply(si, function(i) ax[[i]]$s, 0)
  h  <- hfrac*sd
  if (!length(si)) {
    g <- .mkg(d, so, si, NULL)
    return(list(K=g$K, M=g$M, wk=g$wk, n=d$n, d=d, grids=list(g),
                s=numeric(0), h=numeric(0), nmarg=0L)) }
  pts <- list(mu)
  for (j in seq_along(si)) { a<-mu; a[j]<-a[j]+h[j]; pts[[length(pts)+1L]]<-a
                             b<-mu; b[j]<-b[j]-h[j]; pts[[length(pts)+1L]]<-b }
  gl <- lapply(pts, function(v) .mkg(d, so, si, v))
  list(K=gl[[1]]$K, M=gl[[1]]$M, wk=gl[[1]]$wk, n=d$n, d=d, grids=gl,
       s=sd, h=h, nmarg=length(si)) }

## ---- moments ----------------------------------------------------------------
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

## second-order Taylor combination over the marginalised directions.
## rank1 = FALSE drops sd^2 dE dE^T -- the diagnostic arm, NOT the method.
mallT <- function(p,GT,form,rank1=TRUE){
  bs <- lapply(GT$grids, function(G) mall(p,G,form))
  if (!GT$nmarg) return(bs[[1L]])
  lapply(seq_len(GT$K), function(k){
    b0<-bs[[1L]][[k]]; E<-b0$E; V<-b0$V
    for (j in seq_len(GT$nmarg)) {
      bp<-bs[[2L*j]][[k]]; bm<-bs[[2L*j+1L]][[k]]; hj<-GT$h[j]; sj<-GT$s[j]
      dE <-(bp$E-bm$E)/(2*hj)
      E <- E + 0.5*sj^2*(bp$E-2*b0$E+bm$E)/hj^2
      V <- V + 0.5*sj^2*(bp$V-2*b0$V+bm$V)/hj^2
      if (rank1) V <- V + sj^2*outer(dE,dE) }
    list(E=E,V=V) }) }

## ---- likelihood and fitting -------------------------------------------------
nl1 <- function(o,pr,n){ ch<-tryCatch(chol(pr$V),error=function(e) NULL)
  if(is.null(ch)) return(1e12); iv<-chol2inv(ch); r<-o$E-pr$E
  n*(2*sum(log(diag(ch)))+sum(iv*o$V)+as.numeric(t(r)%*%iv%*%r)) }
## one dispatch point: a plain grid -> mall, a Taylor bundle -> mallT
blocks <- function(p,G,form,rank1=TRUE)
  if (is.null(G$grids)) mall(p,G,form) else mallT(p,G,form,rank1)
## rows touched by ONE objective evaluation for a source (the solve-row proxy)
nrows <- function(G) if (is.null(G$grids)) G$K*G$M else length(G$grids)*G$K*G$M
fitb <- function(BL,GL,FM,nm,start,mx=2500L,rank1=TRUE){
  ne <- 0L
  ob<-function(v){ ne<<-ne+1L; p<-setNames(v,nm); tot<-0
    for (s in seq_along(BL)){ pr<-blocks(p,GL[[s]],FM[s],rank1)
      tot<-tot+sum(vapply(seq_along(BL[[s]]),function(k)
        nl1(BL[[s]][[k]],pr[[k]],GL[[s]]$wk[k]*GL[[s]]$n),0)) }
    if(!is.finite(tot)) 1e12 else tot }
  o<-optim(start[nm],ob,method="Nelder-Mead",control=list(maxit=mx,reltol=1e-11))
  p<-setNames(o$par,nm); po<-intersect(nm,c("oCL","oVc","sig")); p[po]<-abs(p[po])
  attr(p,"nev") <- ne; p }

## ---- structural self-checks before any fitting ------------------------------
say("=========== self-checks ===========")
ok <- TRUE
for (st in list(c(TRUE,TRUE),c(TRUE,FALSE),c(FALSE,TRUE),c(FALSE,FALSE)))
  for (dd in list(DA,DB)) {
    GT <- mkgridT(dd,st); Gm <- mkgrid(dd,st)
    if (!any(st==FALSE)) {
      same <- identical(GT$grids[[1L]], Gm)
      ok <- ok && same
      say(sprintf("strat(%s): AGGt bundle == AGGm grid : %s",
                  paste(substr(as.character(st),1,1),collapse=""), same)) } }
set.seed(11)
for (cn in names(CELLS)) { cl<-CELLS[[cn]]; if (!cl$x2) next
  strA<-c(TRUE,isTRUE(cl$x2)); GA<-mkgrid(DA,strA); GAT<-mkgridT(DA,strA)
  for (t in 1:3) { p<-TP*exp(rnorm(length(TP),0,.15))
    if (!identical(mall(p,GA,cl$f), mallT(p,GAT,cl$f))) {
      ok<-FALSE; say(sprintf("  MISMATCH %s draw %d",cn,t)) } } }
say(sprintf("mall == mallT on every fully-stratified cell (3 random draws each): %s", ok))
if (!ok) stop("AGGt does not reduce to AGGm where nothing is marginalised")

## ---- moment accuracy of the Taylor construction ------------------------------
## reference = the same marginal moments with a 21-node inner quadrature, so the
## comparison separates Taylor truncation from the study's own 4-node grid.
##
## TWO sides, and they are not the same problem:
##   obs  A's OWN published model on A's grid. In the marginalising cells that
##        model carries no bCL2 at all, so its moments are FLAT in x2: every
##        construction agrees to machine precision and the rank-one term is
##        exactly 0. That is the premise of the study, not a defect of the check.
##   pred the UNIFIED 2-cmt model (B's fitted vector) on A's grid -- what the
##        objective actually varies. bCL2 is free here, so this is the only place
##        the marginalised covariate enters, and the rank-one term sd^2 dE dE^T
##        is the entire channel through which A's variance budget sees bCL2.
relerr <- function(a,b) max(abs(a-b)/pmax(abs(b),1e-12))
d0 <- readRDS(file.path(OUT,"overnight.rds")); RES0 <- d0$RES
say("\n=========== marginal moment accuracy vs a 21-node reference ===========")
say(sprintf("%-10s %5s %5s %12s %12s %12s %12s %11s", "cell","marg?","side",
            "relE AGGm","relV AGGm","relE AGGt","relV AGGt","rank1 frac"))
MOMCHK <- list()
for (cn in names(CELLS)) { cl<-CELLS[[cn]]; strA<-c(TRUE,isTRUE(cl$x2))
  ix <- which(vapply(RES0,`[[`,"","cell")==cn)
  Gr<-mkgrid(DA,strA,GREF); Gm<-mkgrid(DA,strA); GT<-mkgridT(DA,strA)
  for (side in c("obs","pred")) {
    th <- if (side=="obs") RES0[[ix[1]]]$thA else RES0[[ix[1]]]$thB
    fm <- if (side=="obs") cl$f else "2"
    br<-mall(th,Gr,fm); bm<-mall(th,Gm,fm); bt<-mallT(th,GT,fm)
    e <- c(relE_AGGm=mean(vapply(seq_along(br),function(k) relerr(bm[[k]]$E,br[[k]]$E),0)),
           relV_AGGm=mean(vapply(seq_along(br),function(k) relerr(bm[[k]]$V,br[[k]]$V),0)),
           relE_AGGt=mean(vapply(seq_along(br),function(k) relerr(bt[[k]]$E,br[[k]]$E),0)),
           relV_AGGt=mean(vapply(seq_along(br),function(k) relerr(bt[[k]]$V,br[[k]]$V),0)))
    ## how big is the rank-one term relative to the block covariance it enters?
    r1 <- if (GT$nmarg) { bs<-lapply(GT$grids,function(G) mall(th,G,fm))
        mean(vapply(seq_len(GT$K), function(k){
          dE<-(bs[[2]][[k]]$E-bs[[3]][[k]]$E)/(2*GT$h[1]); R<-GT$s[1]^2*outer(dE,dE)
          sqrt(sum(R^2))/sqrt(sum(bt[[k]]$V^2)) },0)) } else NA_real_
    MOMCHK[[paste(cn,side)]] <- c(e, rank1=r1)
    say(sprintf("%-10s %5s %5s %12.3e %12.3e %12.3e %12.3e %11s", cn,
        if (cl$x2) "no" else "x2", side, e[[1]],e[[2]],e[[3]],e[[4]],
        if (is.na(r1)) "-" else sprintf("%.4f", r1))) } }
say("(relE/relV = max relative deviation of the block moments, averaged over blocks;")
say(" rank1 frac = ||sd^2 dE dE^T||_F / ||V_AGGt||_F, averaged over blocks)")

## ---- replay every stored replicate, all three arms ---------------------------
GB  <- mkgrid(DB, c(TRUE,TRUE))        # B always fits both covariates
GBT <- mkgridT(DB, c(TRUE,TRUE))
RDS <- file.path(OUT,"aggregate-marginal-taylor.rds")
RESUME <- file.exists(RDS) && !identical(Sys.getenv("ADM_REFIT"), "1")
OUTL <- list()
if (RESUME) { OUTL <- readRDS(RDS)$RES
  say(sprintf("\nresuming from a stored run: %d replicates (ADM_REFIT=1 forces a refit)\n",
              length(OUTL)))
} else {
say(sprintf("\nreplaying %d stored replicates, 4 fits each\n", length(RES0)))
for (i in seq_along(RES0)) { r <- RES0[[i]]; cl <- CELLS[[r$cell]]
  strA <- c(TRUE, isTRUE(cl$x2))
  GAf <- mkgrid(DA, c(TRUE,TRUE))      # AGG: full product grid for A too
  GAm <- mkgrid(DA, strA)              # AGGm
  GAT <- mkgridT(DA, strA)             # AGGt
  tk <- function(f) { t0<-proc.time()[["elapsed"]]; v<-try(f(),silent=TRUE)
                      list(v=v, t=proc.time()[["elapsed"]]-t0) }
  a1 <- tk(function() fitb(list(mall(r$thA,GAf,cl$f), mall(r$thB,GB,"2")),
                           list(GAf,GB), c("2","2"), UN, r$thB))
  a2 <- tk(function() fitb(list(mall(r$thA,GAm,cl$f), mall(r$thB,GB,"2")),
                           list(GAm,GB), c("2","2"), UN, r$thB))
  a3 <- tk(function() fitb(list(mallT(r$thA,GAT,cl$f), mallT(r$thB,GBT,"2")),
                           list(GAT,GBT), c("2","2"), UN, r$thB))
  a4 <- tk(function() fitb(list(mallT(r$thA,GAT,cl$f,FALSE), mallT(r$thB,GBT,"2",FALSE)),
                           list(GAT,GBT), c("2","2"), UN, r$thB, rank1=FALSE))
  if (any(vapply(list(a1,a2,a3,a4), function(z) inherits(z$v,"try-error"), TRUE))) {
    say(sprintf("  %-9s rep %2d FAILED", r$cell, r$rep)); next }
  OUTL[[length(OUTL)+1L]] <- c(r, list(
    AGGf=a1$v, AGGm=a2$v, AGGt=a3$v, AGGt0=a4$v, strat_x2=strA[2],
    tim=c(AGG=a1$t,AGGm=a2$t,AGGt=a3$t,AGGt0=a4$t),
    nev=c(AGG=attr(a1$v,"nev"),AGGm=attr(a2$v,"nev"),
          AGGt=attr(a3$v,"nev"),AGGt0=attr(a4$v,"nev")),
    rows=c(AGG=nrows(GAf)+nrows(GB), AGGm=nrows(GAm)+nrows(GB),
           AGGt=nrows(GAT)+nrows(GBT), AGGt0=nrows(GAT)+nrows(GBT)),
    blk=c(AGG=GAf$K+GB$K, AGGm=GAm$K+GB$K, AGGt=GAT$K+GBT$K, AGGt0=GAT$K+GBT$K)))
  say(sprintf("[%5.0fs] %3d/%d %-9s rep %2d | bCL2 gold %.3f AGG %.3f AGGm %.3f AGGt %.3f | %.1f %.1f %.1f s",
      proc.time()[["elapsed"]]-T0, i, length(RES0), r$cell, r$rep, r$GOLD[["bCL2"]],
      a1$v[["bCL2"]], a2$v[["bCL2"]], a3$v[["bCL2"]], a1$t, a2$t, a3$t))
}
say(sprintf("[%5.0fs] %d/%d replicates completed", proc.time()[["elapsed"]]-T0,
            length(OUTL), length(RES0)))
saveRDS(list(RES=OUTL, TP=TP, MOMCHK=MOMCHK, HFRAC=HFRAC), RDS)
}

## ---- reproduction checks against the stored arms -----------------------------
say("\n=========== reproduction of the stored AGG / AGGm ===========")
dm <- tryCatch(readRDS(file.path(OUT,"aggregate-marginal.rds")), error=function(e) NULL)
key <- function(z) paste(z$cell, z$rep)
mxAGG <- max(vapply(OUTL, function(z) max(abs(z$AGGf-z$AGG)), 0))
say(sprintf("AGG  refit vs stored : max abs diff %.3e over %d reps", mxAGG, length(OUTL)))
if (!is.null(dm)) { st <- setNames(lapply(dm$RES,`[[`,"AGGm"), vapply(dm$RES,key,""))
  d2 <- vapply(OUTL, function(z){ s<-st[[key(z)]]
    if (is.null(s)) NA_real_ else max(abs(z$AGGm - s)) }, 0)
  say(sprintf("AGGm refit vs stored : max abs diff %.3e over %d reps",
              max(d2,na.rm=TRUE), sum(!is.na(d2)))) }

## ---- CONTROL: where nothing is marginalised the three arms must coincide ------
say("\n=========== control: AGG vs AGGm vs AGGt, max abs diff over ALL 10 params ===")
say(sprintf("%-10s %6s %6s %16s %16s %16s", "cell","reps","marg?",
            "|AGGm - AGG|","|AGGt - AGGm|","|AGGt - AGG|"))
cells0 <- vapply(OUTL, `[[`, "", "cell")
for (cn in names(CELLS)) { ix <- which(cells0==cn); if (!length(ix)) next
  f <- function(a,b) max(vapply(ix, function(i)
    max(abs(unclass(OUTL[[i]][[a]]) - unclass(OUTL[[i]][[b]]))), 0))
  say(sprintf("%-10s %6d %6s %16.3e %16.3e %16.3e", cn, length(ix),
      if (isTRUE(CELLS[[cn]]$x2)) "none" else "x2",
      f("AGGm","AGGf"), f("AGGt","AGGm"), f("AGGt","AGGf"))) }
say("(baseline/struct stratify BOTH covariates -- the three constructions are the")
say(" same objective there, so anything but exact 0 is a bug)")

## ---- report -----------------------------------------------------------------
cells <- vapply(OUTL, `[[`, "", "cell"); rmse <- function(a,b) sqrt(mean((a-b)^2))
ARMS <- c(AGG="AGG", AGGm="AGGm", AGGt="AGGt", PLLs="PLLs")
say("\n=========== deviation from the INDIVIDUAL-DATA fit (primary) ===========")
say(sprintf("%-10s %6s %10s %10s %10s %10s", "cell","reps","AGG","AGGm","AGGt","PLL"))
for (cn in names(CELLS)) { ix <- which(cells==cn); if (!length(ix)) next
  v <- lapply(ARMS, function(a) vapply(ix, function(i) rmse(OUTL[[i]][[a]],OUTL[[i]]$GOLD), 0))
  say(sprintf("%-10s %6d %10.4f %10.4f %10.4f %10.4f", cn, length(ix),
      mean(v$AGG), mean(v$AGGm), mean(v$AGGt), mean(v$PLLs))) }

say("\n=========== deviation from the simulation truth ===========")
say(sprintf("%-10s %6s %10s %10s %10s %10s %10s", "cell","reps","IPD (gold)","AGG","AGGm","AGGt","PLL"))
for (cn in names(CELLS)) { ix <- which(cells==cn); if (!length(ix)) next
  g <- vapply(ix, function(i) rmse(OUTL[[i]]$GOLD,TP), 0)
  v <- lapply(ARMS, function(a) vapply(ix, function(i) rmse(OUTL[[i]][[a]],TP), 0))
  say(sprintf("%-10s %6d %10.4f %10.4f %10.4f %10.4f %10.4f", cn, length(ix),
      mean(g), mean(v$AGG), mean(v$AGGm), mean(v$AGGt), mean(v$PLLs))) }

say(sprintf("\n=========== bCL2 (truth %.3f) -- the covariate A omits ===========", TP[["bCL2"]]))
say(sprintf("%-10s %5s %13s %13s %13s %13s %13s", "cell","reps","IPD (gold)","AGG","AGGm","AGGt","PLL"))
for (cn in names(CELLS)) { ix <- which(cells==cn); if (!length(ix)) next
  f <- function(a){ v<-vapply(ix,function(i) OUTL[[i]][[a]][["bCL2"]],0)
    sprintf("%6.3f%+7.3f", mean(v), mean(v)-TP[["bCL2"]]) }
  say(sprintf("%-10s %5d %13s %13s %13s %13s %13s", cn, length(ix),
      f("GOLD"), f("AGG"), f("AGGm"), f("AGGt"), f("PLLs"))) }
say("(mean, and bias against truth)")

say("\n=========== oCL (truth 0.260) -- the variance channel ===========")
say(sprintf("%-10s %12s %12s %12s %12s %12s","cell","IPD (gold)","AGG","AGGm","AGGt","A's published"))
for (cn in names(CELLS)) { ix <- which(cells==cn); if (!length(ix)) next
  f <- function(a) sprintf("%12.3f", mean(vapply(ix,function(i) OUTL[[i]][[a]][["oCL"]],0)))
  say(sprintf("%-10s %12s %12s %12s %12s %12.3f", cn, f("GOLD"), f("AGG"), f("AGGm"),
      f("AGGt"), mean(vapply(ix,function(i) OUTL[[i]]$thA[["oCL"]],0)))) }

## ---- PAIRED AGGt vs AGGm ------------------------------------------------------
say("\n=========== PAIRED: AGGt - AGGm ===========")
say(sprintf("%-10s %5s %11s %11s %11s %9s %9s %9s", "cell","reps",
            "d bCL2","se","p (t)","t wins","m wins","p (sign)"))
PAIR <- list()
for (cn in names(CELLS)) { ix <- which(cells==cn); if (!length(ix)) next
  bt <- vapply(ix,function(i) OUTL[[i]]$AGGt[["bCL2"]],0)
  bm <- vapply(ix,function(i) OUTL[[i]]$AGGm[["bCL2"]],0)
  d  <- bt-bm; se <- sd(d)/sqrt(length(d))
  ae <- abs(bt-TP[["bCL2"]]) - abs(bm-TP[["bCL2"]])   # <0 -> AGGt closer to truth
  wt <- sum(ae < 0); wm <- sum(ae > 0); tie <- sum(ae == 0)
  pt <- if (sd(d)>0) t.test(bt,bm,paired=TRUE)$p.value else NA_real_
  ps <- if (wt+wm>0) binom.test(wt,wt+wm)$p.value else NA_real_
  PAIR[[cn]] <- list(d=d, ae=ae, wt=wt, wm=wm, tie=tie, p_t=pt, p_sign=ps)
  say(sprintf("%-10s %5d %11.4f %11.4f %11s %9d %9d %9s", cn, length(ix),
      mean(d), se, if (is.na(pt)) "  (exact 0)" else sprintf("%11.4f",pt),
      wt, wm, if (is.na(ps)) "-" else sprintf("%.4f",ps))) }
say("(d bCL2 = mean paired AGGt - AGGm; wins = closer to truth 0.450; ties excluded from the sign test)")

say("\npaired |error| in bCL2, and paired RMSE-to-gold")
say(sprintf("%-10s %13s %13s %13s %13s", "cell","d|err| bCL2","se","d rmse-gold","se"))
for (cn in names(CELLS)) { ix <- which(cells==cn); if (!length(ix)) next
  ae <- PAIR[[cn]]$ae
  dg <- vapply(ix,function(i) rmse(OUTL[[i]]$AGGt,OUTL[[i]]$GOLD)-
                              rmse(OUTL[[i]]$AGGm,OUTL[[i]]$GOLD),0)
  say(sprintf("%-10s %13.4f %13.4f %13.4f %13.4f", cn, mean(ae), sd(ae)/sqrt(length(ae)),
      mean(dg), sd(dg)/sqrt(length(dg)))) }
say("(negative favours AGGt)")

## ---- the rank-one term is load-bearing ---------------------------------------
say("\n=========== diagnostic: AGGt with the rank-one term DROPPED ===========")
say(sprintf("%-10s %13s %13s %13s %13s", "cell","bCL2 AGGm","bCL2 AGGt","bCL2 AGGt-r1","rmse AGGt-r1"))
for (cn in names(CELLS)) { ix <- which(cells==cn); if (!length(ix)) next
  f <- function(a) mean(vapply(ix,function(i) OUTL[[i]][[a]][["bCL2"]],0))
  say(sprintf("%-10s %13.3f %13.3f %13.3f %13.4f", cn, f("AGGm"), f("AGGt"), f("AGGt0"),
      mean(vapply(ix,function(i) rmse(OUTL[[i]]$AGGt0,OUTL[[i]]$GOLD),0)))) }
r1v <- vapply(MOMCHK, function(z) z[["rank1"]], 0)
r1v <- r1v[grepl(" pred$", names(r1v)) & !is.na(r1v)]
say(sprintf("rank-one Frobenius share of the AGGt block covariance (predicted side): %s",
    paste(sprintf("%s %.4f", sub(" pred$","",names(r1v)), r1v), collapse="  ")))

## ---- cost --------------------------------------------------------------------
say("\n=========== COST ===========")
say(sprintf("%-10s %8s %10s %10s %10s %10s", "cell","arm","rows/obj","blocks","obj calls","s / rep"))
for (cn in names(CELLS)) { ix <- which(cells==cn); if (!length(ix)) next
  for (a in c("AGG","AGGm","AGGt")) {
    say(sprintf("%-10s %8s %10.0f %10.0f %10.0f %10.2f", cn, a,
        mean(vapply(ix,function(i) OUTL[[i]]$rows[[a]],0)),
        mean(vapply(ix,function(i) OUTL[[i]]$blk[[a]],0)),
        mean(vapply(ix,function(i) OUTL[[i]]$nev[[a]],0)),
        mean(vapply(ix,function(i) OUTL[[i]]$tim[[a]],0)))) } }
say("(rows/obj = grid rows evaluated per objective call, both sources summed;")
say(" blocks = Cholesky factorisations per objective call, both sources summed)")
say("\nwhat the cost argument actually rests on -- the multiplier the marginalised")
say("directions put on the inner expectation, quadrature NC^p vs Taylor 1+2p:")
say(sprintf("%12s %12s %12s %12s %12s", "marg covs","quad NC=4","quad NC=8","quad NC=21","taylor"))
for (pp in 1:4) say(sprintf("%12d %12d %12d %12d %12d", pp, 4^pp, 8^pp, 21^pp, 1+2*pp))
say("This study sits at p = 1, NC = 4: the cheapest cell of that table, so it is the")
say("LEAST favourable setting for Taylor there is. Read the speed numbers as a floor.")

## ---- h sweep -----------------------------------------------------------------
## reported because it was run, NOT used to pick h: the study above is hfrac = 0.5
## throughout, the value taylor-corrected.R validated.
say("\n=========== hfrac sweep (first 8 reps of the marginalising cells) ===========")
SW <- list()
say(sprintf("%-10s %7s %12s %12s %12s %12s", "cell","hfrac","bCL2","mean|err|","rmse-gold","relV vs 21n"))
for (cn in c("omit","both")) { cl<-CELLS[[cn]]; strA<-c(TRUE,isTRUE(cl$x2))
  ix <- head(which(cells==cn), 8L)
  Gr <- mkgrid(DA,strA,GREF)
  for (hf in c(0.25,0.5,0.75,1.0)) {
    GAT<-mkgridT(DA,strA,hf)
    ## relV on the PREDICTED side -- A's own model is flat in x2, so measuring the
    ## truncation there would report 1e-15 for every h and say nothing.
    th<-OUTL[[ix[1]]]$thB
    bt<-mallT(th,GAT,"2"); br<-mall(th,Gr,"2")
    rv <- mean(vapply(seq_len(GAT$K), function(k) relerr(bt[[k]]$V, br[[k]]$V), 0))
    fits <- lapply(ix, function(i){ rr<-OUTL[[i]]
      fitb(list(mallT(rr$thA,GAT,cl$f), mallT(rr$thB,GBT,"2")), list(GAT,GBT),
           c("2","2"), UN, rr$thB) })
    b <- vapply(fits, function(p) p[["bCL2"]], 0)
    g <- vapply(seq_along(ix), function(j) rmse(fits[[j]], OUTL[[ix[j]]]$GOLD), 0)
    SW[[paste(cn,hf)]] <- list(bCL2=b, rmse=g, relV=rv)
    say(sprintf("%-10s %7.2f %12.3f %12.3f %12.4f %12.3e", cn, hf, mean(b),
        mean(abs(b-TP[["bCL2"]])), mean(g), rv)) } }
saveRDS(list(RES=OUTL, TP=TP, MOMCHK=MOMCHK, PAIR=PAIR, SWEEP=SW, HFRAC=HFRAC),
        file.path(OUT,"aggregate-marginal-taylor.rds"))
say(sprintf("\ntotal %.0f s", proc.time()[["elapsed"]]-T0))
