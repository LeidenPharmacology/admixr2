## PARAMETER LIKELIHOOD on two real-shaped publications. NO TRUTH.
##
## Publication A: a 2-compartment popPK model.   Publication B: a 1-compartment one.
## Same drug. For each we have only: the parameter table, its RSEs, the sampling
## design, and a baseline table (weight mean/SD). No data of any kind.
##
## We fit ONE unified 2-compartment model phi by asking:
##
##     what phi would have made A publish what A published,
##     AND made B publish what B published?
##
##   g_A(phi) : A's form IS the unified form -> the identity (verified, not assumed)
##   g_B(phi) : fit a 1-cmt model to what the unified 2-cmt at phi would produce
##              over B's OWN sampling window and population.
##
## g_B is where the work happens. A 1-cmt fit to a 2-cmt drug returns a V that
## depends on WHEN you sampled: early sampling lands near Vc, late sampling near
## Vc + Vp. So B's published V is not "wrong" -- it is a known functional of the
## 2-cmt parameters, and inverting that functional is what lets a 1-cmt paper
## inform Q and Vp, which its own model does not contain.
OUT <- "C:/Users/hidde/AppData/Local/Temp/claude/C--package-admixr2/3ff305c7-64a5-4cb0-ba61-1436e2e9b16e/scratchpad"
say <- function(...) { cat(..., "\n"); utils::flush.console() }; T0 <- proc.time()[["elapsed"]]
gh <- function(k){ i <- seq_len(k-1L); J <- matrix(0,k,k)
  J[cbind(i,i+1L)] <- sqrt(i); J[cbind(i+1L,i)] <- sqrt(i)
  e <- eigen(J,symmetric=TRUE); o <- order(e$values)
  list(x=e$values[o], w=e$vectors[1L,o]^2) }
GE <- gh(5L); GX <- gh(5L)

## ============ THE TWO PUBLICATIONS (this is all the input there is) =========
PUB_A <- c(lcl=log(4.1), lvc=log(28), lq=log(7.5), lvp=log(44),
           bCL=0.78, bVc=0.95, oCL=0.26, oVc=0.21, sig=0.11)
SE_A  <- c(lcl=.050, lvc=.080, lq=.180, lvp=.150,
           bCL=.090, bVc=.110, oCL=.030, oVc=.030, sig=.012)
DES_A <- list(times=c(.25,.5,1,2,4,8,12), dose=500, mu_x=log(68/70), sd_x=.17, n=45)

PUB_B <- c(lcl=log(4.6), lv=log(61), bCL=0.71, bV=1.05,
           oCL=0.31, oV=0.28, sig=0.14)
SE_B  <- c(lcl=.040, lv=.060, bCL=.070, bV=.100, oCL=.035, oV=.035, sig=.010)
DES_B <- list(times=c(8,12,24,36), dose=500, mu_x=log(82/70), sd_x=.22, n=120)

AN <- names(PUB_A); BN <- names(PUB_B)
grid_of <- function(d){ xs <- d$mu_x + d$sd_x*GX$x
  g <- expand.grid(i=seq_along(xs), j=seq_along(GE$x), k=seq_along(GE$x), KEEP.OUT.ATTRS=FALSE)
  w <- GX$w[g$i]*GE$w[g$j]*GE$w[g$k]
  list(x=xs[g$i], z1=GE$x[g$j], z2=GE$x[g$k], w=w/sum(w)) }
GA <- grid_of(DES_A); GB <- grid_of(DES_B)

conc2 <- function(p, gr, tt, D){
  cl <- exp(p[["lcl"]]+p[["bCL"]]*gr$x+p[["oCL"]]*gr$z1)
  vc <- exp(p[["lvc"]]+p[["bVc"]]*gr$x+p[["oVc"]]*gr$z2)
  q <- exp(p[["lq"]]); vp <- exp(p[["lvp"]])
  k10 <- cl/vc; k12 <- q/vc; k21 <- q/vp
  s <- k10+k12+k21; d <- sqrt(pmax(s^2-4*k10*k21,0)); al <- (s+d)/2; be <- (s-d)/2
  A <- (D/vc)*(al-k21)/(al-be); B <- (D/vc)*(k21-be)/(al-be)
  A*exp(-outer(al,tt)) + B*exp(-outer(be,tt)) }
conc1 <- function(p, gr, tt, D){
  cl <- exp(p[["lcl"]]+p[["bCL"]]*gr$x+p[["oCL"]]*gr$z1)
  v  <- exp(p[["lv"]] +p[["bV"]] *gr$x+p[["oV"]] *gr$z2)
  outer(D/v, tt, function(a,b) a)*exp(-outer(cl/v,tt)) }
mom <- function(f, gr, sig){ mu <- as.numeric(crossprod(gr$w,f)); fc <- sweep(f,2L,mu)
  V <- crossprod(fc, fc*gr$w)
  diag(V) <- diag(V) + sig^2*as.numeric(crossprod(gr$w,f^2)); list(E=mu,V=V) }
nll <- function(obs,pr,n){ ch <- tryCatch(chol(pr$V),error=function(e) NULL)
  if(is.null(ch)) return(1e12); iv <- chol2inv(ch); r <- obs$E-pr$E
  n*(2*sum(log(diag(ch)))+sum(iv*obs$V)+as.numeric(t(r)%*%iv%*%r)) }
fitform <- function(obs, cf, gr, d, start, nm){
  ob <- function(v){ p <- stats::setNames(v,nm)
    p[intersect(nm,c("oCL","oV","oVc","sig"))] <- abs(p[intersect(nm,c("oCL","oV","oVc","sig"))])
    nll(obs, mom(cf(p,gr,d$times,d$dose), gr, p[["sig"]]), d$n) }
  o <- stats::optim(start, ob, method="Nelder-Mead", control=list(maxit=4000,reltol=1e-11))
  o <- stats::optim(o$par, ob, method="Nelder-Mead", control=list(maxit=4000,reltol=1e-12))
  p <- stats::setNames(o$par,nm)
  p[intersect(nm,c("oCL","oV","oVc","sig"))] <- abs(p[intersect(nm,c("oCL","oV","oVc","sig"))]); p }
ST_B <- c(lcl=log(4.5), lv=log(58), bCL=.7, bV=1.0, oCL=.30, oV=.27, sig=.13)
ST_A <- PUB_A

g_A <- function(phi) fitform(mom(conc2(phi,GA,DES_A$times,DES_A$dose),GA,phi[["sig"]]),
                             conc2, GA, DES_A, phi[AN], AN)
g_B <- function(phi) fitform(mom(conc2(phi,GB,DES_B$times,DES_B$dose),GB,phi[["sig"]]),
                             conc1, GB, DES_B, ST_B, BN)

say("checking g_A is the identity (A's form == the unified form) ...")
chk <- g_A(PUB_A); say(sprintf("  max |g_A(phi) - phi| = %.2e", max(abs(chk-PUB_A))))
say("what a 1-cmt paper would report under A's model, over B's design:")
gb0 <- g_B(PUB_A)
say(sprintf("  induced  CL %.2f  V %.1f  bCL %.2f  bV %.2f", exp(gb0[["lcl"]]),
            exp(gb0[["lv"]]), gb0[["bCL"]], gb0[["bV"]]))
say(sprintf("  B PUBLISHED CL %.2f  V %.1f  bCL %.2f  bV %.2f", exp(PUB_B[["lcl"]]),
            exp(PUB_B[["lv"]]), PUB_B[["bCL"]], PUB_B[["bV"]]))
say(sprintf("  (A's Vc = %.0f,  A's Vss = %.0f)", exp(PUB_A[["lvc"]]),
            exp(PUB_A[["lvc"]])+exp(PUB_A[["lvp"]])))

## ---------------- fit the unified 2-cmt model to BOTH publications ----------
resid <- function(phi) c((PUB_A - g_A(phi))/SE_A, (PUB_B - g_B(phi))/SE_B)
FD <- c(lcl=5e-3,lvc=5e-3,lq=1.2e-2,lvp=1.2e-2,bCL=1e-2,bVc=1e-2,oCL=5e-3,oVc=5e-3,sig=4e-3)
phi <- PUB_A
for (it in 1:6) {
  r0 <- resid(phi)
  J <- vapply(AN, function(q){ ph <- phi; ph[[q]] <- ph[[q]]+FD[[q]]
        (resid(ph)-r0)/FD[[q]] }, numeric(length(r0)))
  phi <- phi - 0.7*as.numeric(tryCatch(solve(crossprod(J)+diag(1e-8,length(AN)),
                                             crossprod(J,r0)), error=function(e) rep(0,length(AN))))
  say(sprintf("[%3.0fs] GN %d  SSR %8.3f | CL %.2f Vc %.1f Q %.2f Vp %.1f",
      proc.time()[["elapsed"]]-T0, it, sum(r0^2), exp(phi[["lcl"]]),
      exp(phi[["lvc"]]), exp(phi[["lq"]]), exp(phi[["lvp"]])))
}
rF <- resid(phi); V <- solve(crossprod(J)); SEphi <- sqrt(diag(V))
Q <- sum(rF^2); dfQ <- length(rF) - length(AN)
say("\n         A published    unified      SE   (A's own SE)")
for (q in AN) say(sprintf("%-5s %11.4f %10.4f %7.4f %10.4f", q, PUB_A[[q]], phi[[q]],
                          SEphi[[q]], SE_A[[q]]))
say(sprintf("\nnatural  CL %.2f  Vc %.1f  Q %.2f  Vp %.1f   (A alone: %.2f %.1f %.2f %.1f)",
  exp(phi[["lcl"]]),exp(phi[["lvc"]]),exp(phi[["lq"]]),exp(phi[["lvp"]]),
  exp(PUB_A[["lcl"]]),exp(PUB_A[["lvc"]]),exp(PUB_A[["lq"]]),exp(PUB_A[["lvp"]])))
say(sprintf("\nAGREEMENT  Q = %.2f on %d df,  p = %.4f", Q, dfQ,
            stats::pchisq(Q,dfQ,lower.tail=FALSE)))
say(sprintf("SE shrink vs A alone: Q %.2f -> %.2f   Vp %.2f -> %.2f",
            SE_A[["lq"]], SEphi[["lq"]], SE_A[["lvp"]], SEphi[["lvp"]]))

## ---- induced map traces, for the figure ------------------------------------
say("\ntracing the induced map ...")
QG <- exp(seq(log(2), log(22), length.out = 13))
MAPQ <- vapply(QG, function(qq){ ph <- phi; ph[["lq"]] <- log(qq)
  gg <- g_B(ph); c(V=exp(gg[["lv"]]), CL=exp(gg[["lcl"]])) }, numeric(2))
VPG <- exp(seq(log(12), log(120), length.out = 13))
MAPV <- vapply(VPG, function(vv){ ph <- phi; ph[["lvp"]] <- log(vv)
  gg <- g_B(ph); c(V=exp(gg[["lv"]]), CL=exp(gg[["lcl"]])) }, numeric(2))
saveRDS(list(PUB_A=PUB_A,SE_A=SE_A,PUB_B=PUB_B,SE_B=SE_B,DES_A=DES_A,DES_B=DES_B,
             phi=phi,SEphi=SEphi,Q=Q,dfQ=dfQ,gb0=gb0,QG=QG,MAPQ=MAPQ,VPG=VPG,MAPV=MAPV),
        file.path(OUT,"param-two-pubs.rds"))
say(sprintf("total %.0f s", proc.time()[["elapsed"]]-T0))
