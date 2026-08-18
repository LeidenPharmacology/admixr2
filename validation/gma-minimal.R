## MINIMAL Generalized Model Aggregation (GMA) -- CORRECTED.
## Rahmandad, Jalali & Paynabar (2017) PLoS ONE 12(4):e0175111
## = Composite Indirect Inference (Gourieroux & Monfort 2018).
##
## THE FIX over the first version. The binding function g_s(phi) must reproduce
## source s's ESTIMATOR, so it must show that estimator the covariate the way the
## estimator actually saw it:
##
##   s published a covariate coefficient  -> its authors had INDIVIDUAL data
##      resolving weight per subject. So generate PER-NODE (E_k, V_k) across
##      weight nodes with n_k = w_k * n_s, and refit s's form to all of them.
##   s published NO covariate coefficient -> its estimate IS E_{P_s}[.]. So
##      generate ONE pooled (E, V) marginal over P_s and refit to that.
##
## The first version collapsed to a single pooled (E,V) and then fitted B's
## covariate-CARRYING form to it -- recovering a covariate coefficient from a
## marginal summary, which is the configuration the data-processing inequality
## condemns (measured -74% elsewhere in this validation set).
##
## Notation: C_s = covariance of study s's PUBLISHED estimates. `sig` remains the
## residual error, as everywhere else in pharmacometrics.
say <- function(...) { cat(...,"\n"); utils::flush.console() }; T0 <- proc.time()[["elapsed"]]
gh <- function(k){ i <- seq_len(k-1L); J <- matrix(0,k,k)
  J[cbind(i,i+1L)] <- sqrt(i); J[cbind(i+1L,i)] <- sqrt(i)
  e <- eigen(J,symmetric=TRUE); o <- order(e$values); list(x=e$values[o],w=e$vectors[1L,o]^2) }
GE <- gh(5L); GX <- gh(7L)          # eta nodes; covariate nodes

## eta-only grid at a FIXED covariate value -- the conditional moments at node a
grid_at <- function(a){ g <- expand.grid(j=seq_along(GE$x),k=seq_along(GE$x),KEEP.OUT.ATTRS=FALSE)
  w <- GE$w[g$j]*GE$w[g$k]
  list(x=rep(a,nrow(g)), z1=GE$x[g$j], z2=GE$x[g$k], w=w/sum(w)) }
## pooled grid -- covariate marginalised in, for covariate-free sources
grid_marg <- function(d){ xs <- d$mu_x+d$sd_x*GX$x
  g <- expand.grid(i=seq_along(xs),j=seq_along(GE$x),k=seq_along(GE$x),KEEP.OUT.ATTRS=FALSE)
  w <- GX$w[g$i]*GE$w[g$j]*GE$w[g$k]
  list(x=xs[g$i],z1=GE$x[g$j],z2=GE$x[g$k],w=w/sum(w)) }
nodes_of <- function(d) list(a = d$mu_x + d$sd_x*GX$x, w = GX$w/sum(GX$w))

A_th <- c(lcl=log(4.1), lvc=log(28), lq=log(7.5), lvp=log(44),
          bCL=.78, bVc=.95, oCL=.26, oVc=.21, sig=.11)
A_se <- c(lcl=.050,lvc=.075,lq=.170,lvp=.140,bCL=.085,bVc=.105,oCL=.030,oVc=.030,sig=.012)
A_d <- list(times=c(.083,.25,.5,.75,1,1.5,2,3,4,6,8), dose=500, tau=NA, mu_x=0, sd_x=.17, n=45)
B_d <- list(times=c(1,12), dose=500, tau=12, mu_x=log(84/70), sd_x=.22, n=120)
AN <- names(A_th); BN <- c("lcl","lv","bCL","bV","oCL","oV","sig")

ssf <- function(e,tau) if (is.na(tau)) 1 else 1/(1-exp(-e*tau))
f2 <- function(p,gr,d){ cl<-exp(p[["lcl"]]+p[["bCL"]]*gr$x+p[["oCL"]]*gr$z1)
  vc<-exp(p[["lvc"]]+p[["bVc"]]*gr$x+p[["oVc"]]*gr$z2); q<-exp(p[["lq"]]); vp<-exp(p[["lvp"]])
  k10<-cl/vc;k12<-q/vc;k21<-q/vp; s<-k10+k12+k21; dd<-sqrt(pmax(s^2-4*k10*k21,0))
  al<-(s+dd)/2; be<-(s-dd)/2; A<-(d$dose/vc)*(al-k21)/(al-be); B<-(d$dose/vc)*(k21-be)/(al-be)
  (A*ssf(al,d$tau))*exp(-outer(al,d$times))+(B*ssf(be,d$tau))*exp(-outer(be,d$times)) }
f1 <- function(p,gr,d){ cl<-exp(p[["lcl"]]+p[["bCL"]]*gr$x+p[["oCL"]]*gr$z1)
  v<-exp(p[["lv"]]+p[["bV"]]*gr$x+p[["oV"]]*gr$z2); k<-cl/v
  outer((d$dose/v)*ssf(k,d$tau),d$times,function(a,b)a)*exp(-outer(k,d$times)) }
mom <- function(f,gr,sig){ mu<-as.numeric(crossprod(gr$w,f)); fc<-sweep(f,2L,mu)
  V<-crossprod(fc,fc*gr$w); diag(V)<-diag(V)+sig^2*as.numeric(crossprod(gr$w,f^2)); list(E=mu,V=V) }
nll1 <- function(o,pr,n){ ch<-tryCatch(chol(pr$V),error=function(e) NULL)
  if(is.null(ch)) return(1e12); iv<-chol2inv(ch); r<-o$E-pr$E
  n*(2*sum(log(diag(ch)))+sum(iv*o$V)+as.numeric(t(r)%*%iv%*%r)) }

## ---- CONDITIONAL construction: per-node blocks -----------------------------
make_blocks <- function(p, cf, d){ nd <- nodes_of(d)
  lapply(seq_along(nd$a), function(k){ gr <- grid_at(nd$a[k])
    m <- mom(cf(p,gr,d), gr, p[["sig"]])
    list(gr=gr, E=m$E, V=m$V, n=nd$w[k]*d$n) }) }
refit_blocks <- function(bl, cf, d, start, nm){
  ob <- function(v){ p<-setNames(v,nm); po<-intersect(nm,c("oCL","oV","oVc","sig")); p[po]<-abs(p[po])
    sum(vapply(bl, function(b) nll1(b, mom(cf(p,b$gr,d),b$gr,p[["sig"]]), b$n), 0)) }
  o <- optim(start,ob,method="Nelder-Mead",control=list(maxit=8000,reltol=1e-13))
  o <- optim(o$par,ob,method="Nelder-Mead",control=list(maxit=8000,reltol=1e-14))
  p <- setNames(o$par,nm); po<-intersect(nm,c("oCL","oV","oVc","sig")); p[po]<-abs(p[po]); p }

## B PUBLISHED covariate coefficients -> per-node construction
g_A <- function(phi) phi[AN]
g_B <- function(phi) refit_blocks(make_blocks(phi, f2, B_d), f1, B_d,
                                  c(lcl=log(4.2),lv=log(55),bCL=.75,bV=1,oCL=.28,oV=.25,sig=.12), BN)

say("g_B(A) -- per-node (CORRECT) vs pooled (the earlier error):")
gb_cond <- g_B(A_th)
grM <- grid_marg(B_d)
gb_pool <- { ob <- function(v){ p<-setNames(v,BN); p[c("oCL","oV","sig")]<-abs(p[c("oCL","oV","sig")])
    nll1(mom(f2(A_th,grM,B_d),grM,A_th[["sig"]]), mom(f1(p,grM,B_d),grM,p[["sig"]]), B_d$n) }
  o<-optim(c(lcl=log(4.2),lv=log(55),bCL=.75,bV=1,oCL=.28,oV=.25,sig=.12),ob,
           method="Nelder-Mead",control=list(maxit=8000,reltol=1e-13))
  p<-setNames(o$par,BN); p[c("oCL","oV","sig")]<-abs(p[c("oCL","oV","sig")]); p }
say(sprintf("  per-node : CL %.2f V %.1f bCL %.3f bV %.3f", exp(gb_cond[["lcl"]]),
            exp(gb_cond[["lv"]]), gb_cond[["bCL"]], gb_cond[["bV"]]))
say(sprintf("  pooled   : CL %.2f V %.1f bCL %.3f bV %.3f", exp(gb_pool[["lcl"]]),
            exp(gb_pool[["lv"]]), gb_pool[["bCL"]], gb_pool[["bV"]]))
say(sprintf("  bCL differs by %.1f%%, bV by %.1f%%",
  100*(gb_pool[["bCL"]]/gb_cond[["bCL"]]-1), 100*(gb_pool[["bV"]]/gb_cond[["bV"]]-1)))

B_th <- gb_cond + c(.9,-1.1,.8,-.6,1.0,-.7,.9)*c(.045,.070,.120,.150,.045,.055,.011)
B_se <- c(lcl=.045,lv=.070,bCL=.120,bV=.150,oCL=.045,oV=.055,sig=.011)

## ---- fit -------------------------------------------------------------------
iA <- diag(1/A_se[AN]^2); iB <- diag(1/B_se[BN]^2)
res <- function(phi) c((A_th-g_A(phi))/A_se[AN], (B_th-g_B(phi))/B_se[BN])
FD <- c(lcl=5e-3,lvc=5e-3,lq=1.2e-2,lvp=1.2e-2,bCL=1e-2,bVc=1e-2,oCL=5e-3,oVc=5e-3,sig=4e-3)
phi <- A_th; best <- Inf; bphi <- phi; hist <- numeric(0)
for (it in 1:7) {
  r0 <- res(phi); hist <- c(hist, sum(r0^2))
  if (sum(r0^2) < best) { best <- sum(r0^2); bphi <- phi }
  J <- vapply(AN, function(q){ ph<-phi; ph[[q]]<-ph[[q]]+FD[[q]]; (res(ph)-r0)/FD[[q]] },
              numeric(length(r0)))
  phi <- phi - 0.6*as.numeric(solve(crossprod(J)+diag(1e-7,length(AN)), crossprod(J,r0)))
  say(sprintf("[%3.0fs] GN %d SSR %8.4f | CL %.2f Vc %.1f Q %.2f Vp %.1f",
      proc.time()[["elapsed"]]-T0, it, sum(r0^2), exp(phi[["lcl"]]), exp(phi[["lvc"]]),
      exp(phi[["lq"]]), exp(phi[["lvp"]])))
}
phi <- bphi; rF <- res(phi)
J <- vapply(AN, function(q){ ph<-phi; ph[[q]]<-ph[[q]]+FD[[q]]; (res(ph)-rF)/FD[[q]] }, numeric(length(rF)))
SEp <- sqrt(diag(solve(crossprod(J)))); Q <- sum(rF^2); dfQ <- length(rF)-length(AN)
say(sprintf("\nSSR path: %s", paste(sprintf("%.3f",hist), collapse=" -> ")))
say("\n        A published    unified      SE   A's own SE   ratio")
for (q in AN) say(sprintf("%-5s %11.4f %10.4f %7.4f %10.4f %7.2f",
                          q, A_th[[q]], phi[[q]], SEp[[q]], A_se[[q]], SEp[[q]]/A_se[[q]]))
say(sprintf("\nnatural CL %.2f Vc %.1f Q %.2f Vp %.1f  (A alone %.2f %.1f %.2f %.1f)",
  exp(phi[["lcl"]]),exp(phi[["lvc"]]),exp(phi[["lq"]]),exp(phi[["lvp"]]),
  exp(A_th[["lcl"]]),exp(A_th[["lvc"]]),exp(A_th[["lq"]]),exp(A_th[["lvp"]])))
say(sprintf("Q = %.2f on %d df, p = %.3f  [diagnostic only]", Q, dfQ, pchisq(Q,dfQ,lower.tail=FALSE)))
say(sprintf("total %.0f s", proc.time()[["elapsed"]]-T0))
