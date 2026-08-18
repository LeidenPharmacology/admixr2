## PAGE: how much does the WEIGHT MATRIX matter?
##
## Four ways of weighting the parameter-likelihood residuals, from most to least
## information required of the publication:
##
##   rse        assigned RSEs (15% log-params, 0.12 exponents, 20% variance)
##              -- the PAGE abstracts report NO standard errors, so these are
##              invented. Included as the baseline to be improved on.
##   sens       C_s^-1 proportional to n_s S' Sigma^-1 S, the Fisher information
##              from ONE sensitivity solve. Needs only model, design, n.
##   natural    scale-only: log-params in log units (unit change = 100%),
##              exponents in 1/sd(log WT) units since they enter as beta*log(WT).
##              No solve at all.
##   none       identity. Consistent but scale-dominated; included to show why.
##
## All four are consistent (any PD weight is, in indirect inference); they differ
## in efficiency. The question here is whether the POINT ESTIMATE moves.
##
## METHOD NOTE: one Gauss-Newton Jacobian is computed at Ayuthaya's published
## values and SHARED across the four schemes, so the comparison isolates the
## weight rather than confounding it with different optimiser paths. Two damped
## steps each. This is a linearisation and is reported as one.
suppressMessages({library(rxode2)}); suppressMessages(devtools::load_all(".",quiet=TRUE))
OUT <- "C:/Users/hidde/AppData/Local/Temp/claude/C--package-admixr2/3ff305c7-64a5-4cb0-ba61-1436e2e9b16e/scratchpad"
say <- function(...) { cat(...,"\n"); utils::flush.console() }; T0 <- proc.time()[["elapsed"]]
ln_from <- function(mu,sd) list(meanlog=log(mu^2/sqrt(sd^2+mu^2)), sdlog=sqrt(log(1+sd^2/mu^2)))
WT_ISS <- ln_from(20.2,9.2); WT_ALS <- ln_from(18.1,8.5)
EV_ALS <- et(amt=15,ii=6,addl=4,dur=2); T_ALS <- c(21,23.5)
EV_ISS <- et(amt=17.5,dur=1); T_ISS <- c(0.5,1,1.25,1.5,2,3,4,5,6)
UN <- c("lcl","lvc","lvp","lq","clwt","vpwt","qwt","oCL","oVc","sig")
AN <- c("lcl","lv","oCL","oV","sig")
TH_ISS <- c(lcl=log(0.16),lvc=log(3.86),lvp=log(0.19),lq=log(0.13),
            clwt=0.97,vpwt=1.07,qwt=1.19,oCL=sqrt(0.062),oVc=sqrt(0.048),sig=0.12)
TH_ALS <- c(lcl=log(2.99),lv=log(9.55),oCL=sqrt(0.0223),oV=sqrt(0.0134),sig=0.119)
mk_uni <- function(p) eval(parse(text=sprintf('function(){
  ini({lcl<-%.8f; lvc<-%.8f; lvp<-%.8f; lq<-%.8f; clwt<-%.8f; vpwt<-%.8f; qwt<-%.8f
       eta.cl~%.8f; eta.vc~%.8f; prop.err<-%.8f})
  model({cl<-exp(lcl+eta.cl)*WT^clwt; vc<-exp(lvc+eta.vc)
         vp<-exp(lvp)*WT^vpwt; q<-exp(lq)*WT^qwt; f(centr)<-WT
         d/dt(centr)<--(cl+q)/vc*centr+q/vp*peri
         d/dt(peri)<-q/vc*centr-q/vp*peri
         cp<-centr/vc; cp~prop(prop.err)})}',
  p[["lcl"]],p[["lvc"]],p[["lvp"]],p[["lq"]],p[["clwt"]],p[["vpwt"]],p[["qwt"]],
  max(p[["oCL"]],1e-4)^2,max(p[["oVc"]],1e-4)^2,max(abs(p[["sig"]]),1e-4))))
mod_als <- function(){ ini({lcl<-log(2.99); lv<-log(9.55)
    eta.cl~0.0223; eta.v~0.0134; prop.err<-0.119})
  model({cl<-exp(lcl+eta.cl)*(WT/20); v<-exp(lv+eta.v)*(WT/20); f(centr)<-WT
    d/dt(centr)<--cl/v*centr; cp<-centr/v; cp~prop(prop.err)}) }
gen1 <- function(mod,tt,ev,n,wt,seed=99L) datagen(list(A=list(model=mod,times=tt,ev=ev,
  n=n, cov_dist=list(WT=wt))), control=datagenControl(n_sim=20000L,seed=seed))
g_als <- function(phi){ gg <- gen1(mk_uni(phi),T_ALS,EV_ALS,72L,WT_ALS)
  f <- suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(mod_als,admData(),est="adgh",
       control=adghControl(studies=gg,n_nodes=7L,print=0L,covMethod="none",maxeval=2000L))))
  fx <- nlmixr2est::fixef(f); om <- f$omega
  c(lcl=fx[["lcl"]],lv=fx[["lv"]],oCL=sqrt(om[["eta.cl","eta.cl"]]),
    oV=sqrt(om[["eta.v","eta.v"]]),sig=fx[["prop.err"]]) }

## ---- SENSITIVITY weights: the FULL Gaussian information --------------------
## An aggregate-data likelihood scores BOTH the mean and the covariance, so the
## information has two terms. A first version used only the mean part and gave
## the residual-error parameter a weight of exactly zero (dE/dsig = 0) and the
## variance components orders of magnitude too little. For sufficient statistics
## from n subjects -- mean ~ N(E, V/n), scatter ~ Wishart(V, n) -- the
## Joreskog mean-and-covariance-structure information is
##
##   I_ij = n * (dE_i)' V^-1 (dE_j)  +  (n/2) * tr( V^-1 dV_i V^-1 dV_j )
##
## Both derivative sets come from the same finite differences, so the second term
## is free. V^-1 also removes the need for the ad hoc residual precision the
## first version had to invent.
sens_W <- function(th, mod_of, tt, ev, n, wt) {
  g0 <- gen1(mod_of(th),tt,ev,n,wt)$A; E0 <- g0$E; V0 <- g0$V
  Vi <- chol2inv(chol(V0)); nm <- names(th)
  dE <- list(); dV <- list()
  for (q in nm) { p <- th; h <- max(abs(p[[q]]),0.1)*1e-2; p[[q]] <- p[[q]]+h
    g1 <- gen1(mod_of(p),tt,ev,n,wt)$A
    dE[[q]] <- g1$E - E0; dV[[q]] <- (g1$V - V0)/h; dE[[q]] <- dE[[q]]/h }
  K <- length(nm); M <- matrix(0,K,K,dimnames=list(nm,nm))
  P <- lapply(nm, function(q) Vi %*% dV[[q]])
  for (i in seq_len(K)) for (j in i:K) {
    mean_t <- n * as.numeric(t(dE[[i]]) %*% Vi %*% dE[[j]])
    cov_t  <- (n/2) * sum(P[[i]] * t(P[[j]]))
    M[i,j] <- M[j,i] <- mean_t + cov_t }
  M }
say(sprintf("[%3.0fs] sensitivity weights (full information) ...", proc.time()[["elapsed"]]-T0))
W_sens_I <- sens_W(TH_ISS, mk_uni, T_ISS, EV_ISS, 14L, WT_ISS)
## Alsultan's model rebuilt AT p -- a builder ignoring p makes every derivative
## zero and the whole information block vanish, which is what happened first.
mk_als <- function(p) eval(parse(text=sprintf('function(){
  ini({lcl<-%.8f; lv<-%.8f; eta.cl~%.8f; eta.v~%.8f; prop.err<-%.8f})
  model({cl<-exp(lcl+eta.cl)*(WT/20); v<-exp(lv+eta.v)*(WT/20); f(centr)<-WT
    d/dt(centr)<--cl/v*centr; cp<-centr/v; cp~prop(prop.err)})}',
  p[["lcl"]],p[["lv"]],max(p[["oCL"]],1e-4)^2,max(p[["oV"]],1e-4)^2,
  max(abs(p[["sig"]]),1e-4))))
W_sens_A <- sens_W(TH_ALS, mk_als, T_ALS, EV_ALS, 72L, WT_ALS)
dimnames(W_sens_I) <- list(UN,UN); dimnames(W_sens_A) <- list(AN,AN)
say(sprintf("   Ayuthaya info diag: %s",
  paste(sprintf("%s %.1e",UN,pmax(diag(W_sens_I),0)), collapse=" ")))
say(sprintf("   Alsultan info diag: %s",
  paste(sprintf("%s %.1e",AN,pmax(diag(W_sens_A),0)), collapse=" ")))
## implied SEs, to sanity-check against the invented ones
ise <- function(M) sqrt(pmax(diag(tryCatch(solve(M+diag(1e-10,ncol(M))),
                                           error=function(e) diag(NA,ncol(M)))),0))
say(sprintf("   implied SE (Ayuthaya): %s",
  paste(sprintf("%s %.3f",UN,ise(W_sens_I)), collapse=" ")))

## ---- shared Jacobian --------------------------------------------------------
raw <- function(phi) c(TH_ISS-phi, TH_ALS[AN]-g_als(phi))
FD <- setNames(rep(8e-3,length(UN)),UN); FD[c("clwt","vpwt","qwt")] <- 1.5e-2
r0 <- raw(TH_ISS)
say(sprintf("[%3.0fs] shared Jacobian (%d evaluations) ...", proc.time()[["elapsed"]]-T0, length(UN)))
J <- vapply(UN, function(q){ ph<-TH_ISS; ph[[q]]<-ph[[q]]+FD[[q]]
  (raw(ph)-r0)/FD[[q]] }, numeric(length(r0)))

SD_I <- WT_ISS$sdlog
SEr <- c(rep(.15,4),rep(.12,3),.20*sqrt(0.062),.20*sqrt(0.048),.20*0.12)
SEa <- c(.15,.15,.20*sqrt(0.0223),.20*sqrt(0.0134),.20*0.119)
nat <- c(rep(1,4), rep(SD_I,3), 1,1,1)                # exponents -> log units
WS <- list(
  rse     = diag(c(1/SEr^2, 1/SEa^2)),
  natural = diag(c(nat^2, rep(1,5))),
  none    = diag(1, length(r0)))
WS$sens <- { M <- matrix(0,length(r0),length(r0))
  M[1:10,1:10] <- W_sens_I; M[11:15,11:15] <- W_sens_A; M }
step <- function(W){ p <- TH_ISS
  for (k in 1:2) { r <- if (k==1) r0 else raw(p)
    Jk <- J
    d <- tryCatch(solve(t(Jk)%*%W%*%Jk + diag(1e-8,ncol(Jk)), t(Jk)%*%W%*%r),
                  error=function(e) rep(0,ncol(Jk)))
    p <- p - 0.7*as.numeric(d) }
  p }
say(sprintf("[%3.0fs] solving under each weight ...", proc.time()[["elapsed"]]-T0))
FITS <- lapply(WS[c("rse","sens","natural","none")], step)
AGGv <- tryCatch(readRDS(file.path(OUT,"page-fw.rds")), error=function(e) NULL)
say(sprintf("\n%-24s %8s %8s %8s","weight scheme","clwt","vpwt","qwt"))
say(sprintf("%-24s %8.3f %8.3f %8.3f","Ayuthaya published",
            TH_ISS[["clwt"]],TH_ISS[["vpwt"]],TH_ISS[["qwt"]]))
if (!is.null(AGGv)) say(sprintf("%-24s %8.3f %8.3f %8.3f","AGGREGATE (joint)",
            AGGv$joint[["clwt"]],AGGv$joint[["vpwt"]],AGGv$joint[["qwt"]]))
for (nm in names(FITS)) say(sprintf("%-24s %8.3f %8.3f %8.3f",
  paste("PARAM-LL /",nm), FITS[[nm]][["clwt"]],FITS[[nm]][["vpwt"]],FITS[[nm]][["qwt"]]))
sp <- vapply(c("clwt","vpwt","qwt"), function(q)
  diff(range(vapply(FITS,function(p) p[[q]],0))), 0)
say(sprintf("\nspread of the point estimate ACROSS weight schemes: %s",
            paste(sprintf("%s %.3f",names(sp),sp), collapse="  ")))
saveRDS(list(FITS=FITS, TH_ISS=TH_ISS, AGG=AGGv, W_sens_I=W_sens_I),
        file.path(OUT,"page-weights.rds"))
say(sprintf("total %.0f s", proc.time()[["elapsed"]]-T0))
