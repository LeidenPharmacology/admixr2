## TWO PUBLICATIONS, same drug, no data of any kind.
##
##  A : 2-compartment, DENSE single-dose sampling that captures the peak and the
##      distribution phase. Resolves Q and Vp.
##  B : 1-compartment, STEADY-STATE peak and trough only (500 mg q12h). Cannot
##      see the distribution phase at all -- which is exactly why its authors
##      fitted one compartment, and why its weight exponents are poorly pinned.
##
## We fit ONE unified 2-cmt model to both parameter tables.
##   g_A(phi) = identity (A's form IS the unified form)
##   g_B(phi) = fit a 1-cmt model to what the unified 2-cmt at phi would produce
##              at STEADY STATE over B's own two sampling times and population.
##
## HOW B's NUMBERS WERE CHOSEN -- stated plainly, because it decides what the Q
## statistic means. We compute g_B(theta_A), i.e. what a ss peak/trough paper
## WOULD have reported if A's model were right, then perturb it by ~1 published SE
## to represent a genuinely separate study. So the two publications are
## constructed to be COMPATIBLE-BUT-NOT-IDENTICAL. Q then tests whether the
## machinery recognises that, rather than testing an incompatibility we built in.
OUT <- "C:/Users/hidde/AppData/Local/Temp/claude/C--package-admixr2/3ff305c7-64a5-4cb0-ba61-1436e2e9b16e/scratchpad"
say <- function(...) { cat(...,"\n"); utils::flush.console() }; T0 <- proc.time()[["elapsed"]]
gh <- function(k){ i <- seq_len(k-1L); J <- matrix(0,k,k)
  J[cbind(i,i+1L)] <- sqrt(i); J[cbind(i+1L,i)] <- sqrt(i)
  e <- eigen(J,symmetric=TRUE); o <- order(e$values)
  list(x=e$values[o], w=e$vectors[1L,o]^2) }
GE <- gh(5L); GX <- gh(5L)

PUB_A <- c(lcl=log(4.1), lvc=log(28), lq=log(7.5), lvp=log(44),
           bCL=0.78, bVc=0.95, oCL=0.26, oVc=0.21, sig=0.11)
SE_A  <- c(lcl=.050, lvc=.075, lq=.170, lvp=.140,
           bCL=.085, bVc=.105, oCL=.030, oVc=.030, sig=.012)
DES_A <- list(times=c(.083,.25,.5,.75,1,1.5,2,3,4,6,8), dose=500, tau=NA,
              mu_x=log(70/70), sd_x=.17, n=45)                  # single dose
DES_B <- list(times=c(1, 12), dose=500, tau=12,
              mu_x=log(84/70), sd_x=.22, n=120)                 # ss peak + trough
## B's design cannot pin exponents: loose SEs there, tight on CL (AUC-driven).
SE_B  <- c(lcl=.045, lv=.070, bCL=.230, bV=.300, oCL=.045, oV=.055, sig=.011)

AN <- names(PUB_A); BN <- names(SE_B)
grid_of <- function(d){ xs <- d$mu_x + d$sd_x*GX$x
  g <- expand.grid(i=seq_along(xs), j=seq_along(GE$x), k=seq_along(GE$x), KEEP.OUT.ATTRS=FALSE)
  w <- GX$w[g$i]*GE$w[g$j]*GE$w[g$k]
  list(x=xs[g$i], z1=GE$x[g$j], z2=GE$x[g$k], w=w/sum(w)) }
GA <- grid_of(DES_A); GB <- grid_of(DES_B)
ss <- function(e, tau) if (is.na(tau)) 1 else 1/(1-exp(-e*tau))   # superposition
conc2 <- function(p, gr, tt, D, tau=NA){
  cl <- exp(p[["lcl"]]+p[["bCL"]]*gr$x+p[["oCL"]]*gr$z1)
  vc <- exp(p[["lvc"]]+p[["bVc"]]*gr$x+p[["oVc"]]*gr$z2)
  q <- exp(p[["lq"]]); vp <- exp(p[["lvp"]])
  k10 <- cl/vc; k12 <- q/vc; k21 <- q/vp
  s <- k10+k12+k21; d <- sqrt(pmax(s^2-4*k10*k21,0)); al <- (s+d)/2; be <- (s-d)/2
  A <- (D/vc)*(al-k21)/(al-be); B <- (D/vc)*(k21-be)/(al-be)
  (A*ss(al,tau))*exp(-outer(al,tt)) + (B*ss(be,tau))*exp(-outer(be,tt)) }
conc1 <- function(p, gr, tt, D, tau=NA){
  cl <- exp(p[["lcl"]]+p[["bCL"]]*gr$x+p[["oCL"]]*gr$z1)
  v  <- exp(p[["lv"]] +p[["bV"]] *gr$x+p[["oV"]] *gr$z2); k <- cl/v
  outer((D/v)*ss(k,tau), tt, function(a,b) a)*exp(-outer(k,tt)) }
mom <- function(f,gr,sig){ mu <- as.numeric(crossprod(gr$w,f)); fc <- sweep(f,2L,mu)
  V <- crossprod(fc, fc*gr$w)
  diag(V) <- diag(V)+sig^2*as.numeric(crossprod(gr$w,f^2)); list(E=mu,V=V) }
nll <- function(o,pr,n){ ch <- tryCatch(chol(pr$V),error=function(e) NULL)
  if(is.null(ch)) return(1e12); iv <- chol2inv(ch); r <- o$E-pr$E
  n*(2*sum(log(diag(ch)))+sum(iv*o$V)+as.numeric(t(r)%*%iv%*%r)) }
fitform <- function(obs, cf, gr, d, start, nm){
  ob <- function(v){ p <- stats::setNames(v,nm); po <- intersect(nm,c("oCL","oV","oVc","sig"))
    p[po] <- abs(p[po]); nll(obs, mom(cf(p,gr,d$times,d$dose,d$tau),gr,p[["sig"]]), d$n) }
  o <- stats::optim(start,ob,method="Nelder-Mead",control=list(maxit=5000,reltol=1e-12))
  o <- stats::optim(o$par,ob,method="Nelder-Mead",control=list(maxit=5000,reltol=1e-13))
  p <- stats::setNames(o$par,nm); po <- intersect(nm,c("oCL","oV","oVc","sig"))
  p[po] <- abs(p[po]); p }
ST_B <- c(lcl=log(4.2), lv=log(55), bCL=.75, bV=1.0, oCL=.28, oV=.25, sig=.12)
g_A <- function(phi) phi[AN]
g_B <- function(phi) fitform(mom(conc2(phi,GB,DES_B$times,DES_B$dose,DES_B$tau),GB,phi[["sig"]]),
                             conc1, GB, DES_B, ST_B, BN)

say("what a ss peak/trough 1-cmt paper would report if A's model were right:")
imp <- g_B(PUB_A)
say(sprintf("  CL %.2f  V %.1f  bCL %.2f  bV %.2f   (A: Vc %.0f, Vss %.0f)",
  exp(imp[["lcl"]]),exp(imp[["lv"]]),imp[["bCL"]],imp[["bV"]],
  exp(PUB_A[["lvc"]]), exp(PUB_A[["lvc"]])+exp(PUB_A[["lvp"]])))
## B's published table: the implied values, displaced by ~1 SE (a separate study)
set.seed(11)
PUB_B <- imp + c(lcl=.9,lv=-1.1,bCL=.8,bV=-.6,oCL=1.0,oV=-.7,sig=.9)[BN]*SE_B[BN]
say(sprintf("  B PUBLISHES        CL %.2f  V %.1f  bCL %.2f  bV %.2f",
  exp(PUB_B[["lcl"]]),exp(PUB_B[["lv"]]),PUB_B[["bCL"]],PUB_B[["bV"]]))

resid <- function(phi) c((PUB_A-g_A(phi))/SE_A, (PUB_B-g_B(phi))/SE_B)
FD <- c(lcl=5e-3,lvc=5e-3,lq=1.2e-2,lvp=1.2e-2,bCL=1e-2,bVc=1e-2,oCL=5e-3,oVc=5e-3,sig=4e-3)
phi <- PUB_A; best <- Inf; bphi <- phi
for (it in 1:9) {
  r0 <- resid(phi); if (sum(r0^2) < best) { best <- sum(r0^2); bphi <- phi }
  J <- vapply(AN, function(q){ ph <- phi; ph[[q]] <- ph[[q]]+FD[[q]]
        (resid(ph)-r0)/FD[[q]] }, numeric(length(r0)))
  st <- tryCatch(solve(crossprod(J)+diag(1e-6,length(AN)), crossprod(J,r0)),
                 error=function(e) rep(0,length(AN)))
  phi <- phi - 0.55*as.numeric(st)
  say(sprintf("[%3.0fs] GN %d SSR %7.3f | CL %.2f Vc %.1f Q %.2f Vp %.1f",
      proc.time()[["elapsed"]]-T0, it, sum(r0^2), exp(phi[["lcl"]]),
      exp(phi[["lvc"]]), exp(phi[["lq"]]), exp(phi[["lvp"]])))
}
phi <- bphi; rF <- resid(phi)
J <- vapply(AN, function(q){ ph <- phi; ph[[q]] <- ph[[q]]+FD[[q]]
      (resid(ph)-rF)/FD[[q]] }, numeric(length(rF)))
SEp <- sqrt(diag(solve(crossprod(J)))); Q <- sum(rF^2); dfQ <- length(rF)-length(AN)
say("\n        A published    unified      SE   A's own SE")
for (q in AN) say(sprintf("%-5s %11.4f %10.4f %7.4f %10.4f", q, PUB_A[[q]], phi[[q]], SEp[[q]], SE_A[[q]]))
say(sprintf("\nnatural  CL %.2f  Vc %.1f  Q %.2f  Vp %.1f    (A alone %.2f %.1f %.2f %.1f)",
  exp(phi[["lcl"]]),exp(phi[["lvc"]]),exp(phi[["lq"]]),exp(phi[["lvp"]]),
  exp(PUB_A[["lcl"]]),exp(PUB_A[["lvc"]]),exp(PUB_A[["lq"]]),exp(PUB_A[["lvp"]])))
say(sprintf("\nAGREEMENT Q = %.2f on %d df, p = %.3f", Q, dfQ, stats::pchisq(Q,dfQ,lower.tail=FALSE)))
say(sprintf("SE vs A alone: %s", paste(sprintf("%s %.3f->%.3f", AN, SE_A[AN], SEp[AN]), collapse="  ")))
say("\ntracing induced map over Q ...")
QG <- exp(seq(log(2), log(25), length.out=13))
MAPQ <- vapply(QG, function(qq){ ph <- phi; ph[["lq"]] <- log(qq); gg <- g_B(ph)
  c(V=exp(gg[["lv"]]), CL=exp(gg[["lcl"]])) }, numeric(2))
saveRDS(list(PUB_A=PUB_A,SE_A=SE_A,PUB_B=PUB_B,SE_B=SE_B,DES_A=DES_A,DES_B=DES_B,
             phi=phi,SEp=SEp,Q=Q,dfQ=dfQ,imp=imp,QG=QG,MAPQ=MAPQ),
        file.path(OUT,"param-peak-trough.rds"))
say(sprintf("total %.0f s", proc.time()[["elapsed"]]-T0))
