## PAGE vancomycin under the matching-principle framework.
##
## Two PUBLISHED models, real ones:
##   Ayuthaya  2-cmt, single 1h infusion 17.5 mg/kg, 9 times, n=14,  WT 20.2+-9.2
##             ESTIMATES weight exponents: clwt 0.97, vpwt 1.07, qwt 1.19
##   Alsultan  1-cmt, 5 infusions q6h x 2h of 15 mg/kg, 2 times, n=72, WT 18.1+-8.5
##             cl ~ (WT/20), v ~ (WT/20), dose ~ WT  -> k = cl/v and C0 = dose/v
##             are BOTH exactly weight-invariant. No weight parameter is estimated.
##
## THE PREDICTION. Alsultan's generated summaries are flat in weight, not because
## its model omits weight but because its parameterisation and mg/kg dosing cancel
## it. Under the aggregate construction that flatness is evidence of NO weight
## effect -- carrying 72 subjects against Ayuthaya's 14. So the weight exponents
## should attenuate toward the values Alsultan implies (linear on cl and v, i.e.
## clwt -> 1 with everything else unconstrained), NOT toward Ayuthaya's published
## values, and the design diagnostic should have said so before any fit.
suppressMessages({library(rxode2)}); suppressMessages(devtools::load_all(".", quiet=TRUE))
say <- function(...) { cat(...,"\n"); utils::flush.console() }; T0 <- proc.time()[["elapsed"]]
ln_from <- function(mu,sd) list(meanlog=log(mu^2/sqrt(sd^2+mu^2)), sdlog=sqrt(log(1+sd^2/mu^2)))
WT_ISS <- ln_from(20.2,9.2); WT_ALS <- ln_from(18.1,8.5)
EV_ISS <- et(amt=17.5, dur=1); EV_ALS <- et(amt=15, ii=6, addl=4, dur=2)
T_ISS <- c(0.5,1,1.25,1.5,2,3,4,5,6); T_ALS <- c(21,23.5)
mod_iss <- function() {
  ini({lcl<-log(0.16); lvc<-log(3.86); lvp<-log(0.19); lq<-log(0.13)
       clwt<-0.97; vpwt<-1.07; qwt<-1.19
       eta.cl~0.062; eta.vc~0.048; prop.err<-0.12})
  model({cl<-exp(lcl+eta.cl)*WT^clwt; vc<-exp(lvc+eta.vc)
         vp<-exp(lvp)*WT^vpwt; q<-exp(lq)*WT^qwt; f(centr)<-WT
         d/dt(centr)<--(cl+q)/vc*centr+q/vp*peri
         d/dt(peri)<-q/vc*centr-q/vp*peri
         cp<-centr/vc; cp~prop(prop.err)}) }
mod_als <- function() {
  ini({lcl<-log(2.99); lv<-log(9.55); eta.cl~0.0223; eta.v~0.0134; prop.err<-0.119})
  model({cl<-exp(lcl+eta.cl)*(WT/20); v<-exp(lv+eta.v)*(WT/20); f(centr)<-WT
         d/dt(centr)<--cl/v*centr; cp<-centr/v; cp~prop(prop.err)}) }
mod_unify <- function() {
  ini({lcl<-log(0.30); lvc<-log(2.50); lvp<-log(0.40); lq<-log(0.25)
       clwt<-0.60; vpwt<-0.70; qwt<-0.80
       eta.cl~0.10; eta.vc~0.10; prop.err<-0.20})
  model({cl<-exp(lcl+eta.cl)*WT^clwt; vc<-exp(lvc+eta.vc)
         vp<-exp(lvp)*WT^vpwt; q<-exp(lq)*WT^qwt; f(centr)<-WT
         d/dt(centr)<--(cl+q)/vc*centr+q/vp*peri
         d/dt(peri)<-q/vc*centr-q/vp*peri
         cp<-centr/vc; cp~prop(prop.err)}) }

## ---- DIAGNOSTIC 1: can each design see weight at all? ----------------------
## Slice each study's own weight distribution into K bands and measure how much
## the predicted mean profile moves across them. This is the check that should
## precede any fit.
strat <- function(spec, K) setNames(lapply(seq_len(K), function(k){ force(k)
  list(quantile=function(u) qlnorm(((k-1)+u)/K, spec$meanlog, spec$sdlog)) }),
  paste0("b",seq_len(K)))
spread_of <- function(mod, tt, ev, n, spec, K=8L) {
  st <- setNames(lapply(strat(spec,K), function(q)
    list(model=mod, times=tt, ev=ev, n=n/K, cov_dist=list(WT=q))), paste0("s",seq_len(K)))
  g <- datagen(st, control=datagenControl(n_sim=20000L, seed=7L))
  m <- vapply(g, function(u) mean(u$E), 0)
  100*(max(m)-min(m))/mean(m)
}
sp_iss <- spread_of(mod_iss, T_ISS, EV_ISS, 14L, WT_ISS)
sp_als <- spread_of(mod_als, T_ALS, EV_ALS, 72L, WT_ALS)
say(sprintf("[%3.0fs] DESIGN DIAGNOSTIC -- spread of predicted mean across weight bands",
            proc.time()[["elapsed"]]-T0))
say(sprintf("   Ayuthaya (2-cmt, estimates exponents) : %6.2f %%   n = 14", sp_iss))
say(sprintf("   Alsultan (1-cmt, mg/kg, linear CL & V): %6.2f %%   n = 72", sp_als))
say(sprintf("   -> Alsultan carries %.1fx the subjects and %s about weight.",
            72/14, if (sp_als < 1) "says NOTHING" else "says something"))

## ---- the joint aggregate fit ----------------------------------------------
gen <- datagen(list(
  Ayuthaya=list(model=mod_iss, times=T_ISS, ev=EV_ISS, n=14L, cov_dist=list(WT=WT_ISS)),
  Alsultan=list(model=mod_als, times=T_ALS, ev=EV_ALS, n=72L, cov_dist=list(WT=WT_ALS))),
  control=datagenControl(n_sim=40000L, seed=12345L))
ff <- function(st) suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
  mod_unify, admData(), est="adgh",
  control=adghControl(studies=st, n_nodes=7L, print=0L, covMethod="none", maxeval=3000L))))
say(sprintf("[%3.0fs] fitting BOTH studies jointly ...", proc.time()[["elapsed"]]-T0))
fB <- ff(gen)
say(sprintf("[%3.0fs] fitting Ayuthaya ALONE (the only weight-informative source) ...",
            proc.time()[["elapsed"]]-T0))
fA <- ff(gen["Ayuthaya"])
gp <- function(f) unlist(f$env$admExtra$struct[c("clwt","vpwt","qwt")])
PUB <- c(clwt=0.97, vpwt=1.07, qwt=1.19)
say("\n            Ayuthaya published   Ayuthaya alone    JOINT (both)")
for (q in names(PUB)) say(sprintf("%-6s %17.3f %17.3f %15.3f",
  q, PUB[[q]], gp(fA)[[q]], gp(fB)[[q]]))
say(sprintf("\nshift toward zero when Alsultan is added: %s",
  paste(sprintf("%s %+.1f%%", names(PUB),
        100*(gp(fB)[names(PUB)]/gp(fA)[names(PUB)]-1)), collapse="  ")))
saveRDS(list(sp_iss=sp_iss, sp_als=sp_als, joint=gp(fB), alone=gp(fA), pub=PUB),
  "C:/Users/hidde/AppData/Local/Temp/claude/C--package-admixr2/3ff305c7-64a5-4cb0-ba61-1436e2e9b16e/scratchpad/page-fw.rds")
say(sprintf("total %.0f s", proc.time()[["elapsed"]]-T0))
