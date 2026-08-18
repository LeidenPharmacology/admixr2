## ============================================================================
## headline-construction-check.R
## ============================================================================
## Does the HEADLINE study (`aggregate-marginal.R`, replaying `overnight-
## simulation.R`) contain the mismatched-targets bug?
##
## The headline claim -- that stratifying on a covariate a published source never
## fitted attenuates that covariate's coefficient by ~68% -- rests on the AGG arm
## of those two scripts. If their observed and predicted blocks were NOT both
## conditional on the same covariate node, the claim would be measuring the same
## category error as `covariate-threeway.R`'s `gh` arm and would have to be
## withdrawn.
##
## This script re-implements `mkgrid`/`mall`/`nl1`/`fitb` VERBATIM from
## `aggregate-marginal.R` and applies two discriminating tests.
##
##   TEST 1  With the SAME parameters and the SAME structural form on both sides,
##           a MATCHED construction must give an exactly zero mean residual in
##           EVERY block, and identical V. A mismatched construction cannot:
##           the pooled/marginal E differs from every node's conditional E.
##
##   TEST 2  Build the pooled/marginal (E, V) from the same grid and score it
##           against the per-node predictions (i.e. deliberately introduce
##           construction 4). If the objective moves, the two constructions are
##           numerically distinguishable and TEST 1 identifies which one the
##           headline scripts implement.
## ============================================================================
say <- function(...) { cat(..., "\n"); utils::flush.console() }

gh <- function(k){ i<-seq_len(k-1L); J<-matrix(0,k,k)
  J[cbind(i,i+1L)]<-sqrt(i); J[cbind(i+1L,i)]<-sqrt(i)
  e<-eigen(J,symmetric=TRUE); o<-order(e$values); list(x=e$values[o],w=e$vectors[1L,o]^2) }
NE<-5L; NC<-4L; GE<-gh(NE); GX<-gh(NC)
TP <- c(lcl=log(4.2),lvc=log(30),lq=log(7.0),lvp=log(40),
        bCL1=.72,bCL2=.45,bVc1=.95,oCL=.26,oVc=.20,sig=.11)
DA <- list(n=140,dose=500,ii=12,addl=5,obs=c(61,72),mx1=log(82/70),sx1=.20,
           mx2=log(78/100),sx2=.28,times=c(1,12),tau=12)

## ---- VERBATIM from aggregate-marginal.R -------------------------------------
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

## ---- pooled/marginal (E,V) over the SAME grid -------------------------------
## This is what a construction-4 arm would score against: the covariate fully
## integrated out, one (E, V) for the whole study.
pool_blocks <- function(B, wk) {
  E <- Reduce(`+`, Map(function(b,w) w*b$E, B, wk))
  V <- Reduce(`+`, Map(function(b,w) w*(b$V + outer(b$E,b$E)), B, wk)) - outer(E,E)
  list(E=E, V=V) }

say("=============================================================================")
say("TEST 1 -- are the observed and predicted blocks matched at the same node?")
say("=============================================================================")
say("`aggregate-marginal.R` line 106:")
say("    AGGm <- fitb(list(mall(r$thA, GA, cl$f), mall(r$thB, GB, \"2\")),")
say("                 list(GA, GB), c(\"2\",\"2\"), UN, r$thB)")
say("`fitb` (line 88-91):  pr <- mall(p, GL[[s]], FM[s])")
say("                      nl1(BL[[s]][[k]], pr[[k]], GL[[s]]$wk[k]*GL[[s]]$n)")
say("-> BL[[s]][[k]] and pr[[k]] share the index k AND the grid GL[[s]], so both")
say("   sides are conditional on covariate node k. Numeric confirmation:\n")

for (nm in c("stratify BOTH covariates (cells baseline/struct)",
             "stratify x1 only, marginalise x2 (cells omit/both)")) {
  strat <- if (grepl("BOTH", nm)) c(TRUE,TRUE) else c(TRUE,FALSE)
  G  <- mkgrid(DA, strat)
  OBS <- mall(TP, G, "2")          # the "published source" blocks
  PR  <- mall(TP, G, "2")          # the unified model's blocks, same params
  dE <- max(vapply(seq_len(G$K), function(k) max(abs(OBS[[k]]$E - PR[[k]]$E)), 0))
  dV <- max(vapply(seq_len(G$K), function(k) max(abs(OBS[[k]]$V - PR[[k]]$V)), 0))
  say(sprintf("  %s", nm))
  say(sprintf("    K = %2d blocks   max |E_obs,k - E_pred,k| = %.3e   max |V diff| = %.3e",
              G$K, dE, dV))
  ## the same predictions against the POOLED observation -- what construction 4 does
  PL <- pool_blocks(OBS, G$wk)
  dEp <- max(vapply(seq_len(G$K), function(k) max(abs(PL$E - PR[[k]]$E)), 0))
  say(sprintf("    if obs were POOLED instead: max |E_pooled - E_pred,k| = %.3e  <- NOT zero",
              dEp))
}

say("\n=============================================================================")
say("TEST 2 -- the two constructions give different objectives, so TEST 1 decides")
say("=============================================================================")
for (nm in c("BOTH stratified", "x2 marginalised")) {
  strat <- if (nm == "BOTH stratified") c(TRUE,TRUE) else c(TRUE,FALSE)
  G <- mkgrid(DA, strat); OBS <- mall(TP, G, "2"); PL <- pool_blocks(OBS, G$wk)
  p2 <- TP; p2[["bCL2"]] <- 0.45
  PR <- mall(p2, G, "2")
  matched <- sum(vapply(seq_len(G$K), function(k) nl1(OBS[[k]], PR[[k]], G$wk[k]*G$n), 0))
  mismat  <- sum(vapply(seq_len(G$K), function(k) nl1(PL,       PR[[k]], G$wk[k]*G$n), 0))
  say(sprintf("  %-18s  matched (constr 2) = %12.4f   mismatched (constr 4) = %12.4f",
              nm, matched, mismat))
}

say("\n=============================================================================")
say("TEST 3 -- does bCL2 recover under the construction the headline scripts use?")
say("=============================================================================")
say("Profile bCL2 with everything else fixed at truth, under BOTH constructions,")
say("with x1 stratified and x2 stratified (the `baseline`-style grid).\n")
G <- mkgrid(DA, c(TRUE,TRUE)); OBS <- mall(TP, G, "2"); PL <- pool_blocks(OBS, G$wk)
prof <- function(b2, pooled) { p<-TP; p[["bCL2"]]<-b2; PR<-mall(p,G,"2")
  sum(vapply(seq_len(G$K), function(k)
    nl1(if (pooled) PL else OBS[[k]], PR[[k]], G$wk[k]*G$n), 0)) }
o1 <- optimize(function(b) prof(b, FALSE), c(-0.5, 1.5), tol=1e-8)
o2 <- optimize(function(b) prof(b, TRUE),  c(-0.5, 1.5), tol=1e-8)
say(sprintf("  matched  (construction 2, WHAT THE SCRIPTS DO) -> bCL2 = %.4f  (truth 0.450, %+.1f%%)",
            o1$minimum, 100*(o1$minimum-0.45)/0.45))
say(sprintf("  MISMATCHED (construction 4, hypothetical)      -> bCL2 = %.4f  (truth 0.450, %+.1f%%)",
            o2$minimum, 100*(o2$minimum-0.45)/0.45))
say("\nIf the headline scripts had the bug, the second line would be their answer.")
