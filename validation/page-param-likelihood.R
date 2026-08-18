## PAGE vancomycin: the PARAMETER likelihood, and the weight-gradient dosing plot.
##
## The previous script (page-framework-test.R) ran the AGGREGATE construction only.
## This one runs the binding-function construction on the same two publications and
## puts both on the dosing figure.
##
## THE PER-SOURCE RULE, applied. Alsultan publishes NO weight exponent -- its
## (WT/20) scaling is fixed, not estimated -- so under the framework it contributes
## MARGINALLY over its own weight distribution, and its parameter vector contains
## nothing that can vote an exponent toward zero. Ayuthaya's form IS the unified
## form, so g_Ayuthaya is the identity. Prediction: the exponents stay at
## Ayuthaya's values instead of attenuating ~10% as they do under AGG.
##
## CAVEAT, stated because it matters: the PAGE abstracts do not report standard
## errors. The C_s below are ASSIGNED (15% RSE on log-scale structural parameters,
## 0.12 absolute on exponents, 20% on variance components). The point estimates are
## insensitive to this within reason; any uncertainty statement would not be.
suppressMessages({library(rxode2); library(ggplot2)})
suppressMessages(devtools::load_all(".", quiet=TRUE))
OUT <- "C:/Users/hidde/AppData/Local/Temp/claude/C--package-admixr2/3ff305c7-64a5-4cb0-ba61-1436e2e9b16e/scratchpad"
say <- function(...) { cat(...,"\n"); utils::flush.console() }; T0 <- proc.time()[["elapsed"]]
ln_from <- function(mu,sd) list(meanlog=log(mu^2/sqrt(sd^2+mu^2)), sdlog=sqrt(log(1+sd^2/mu^2)))
WT_ISS <- ln_from(20.2,9.2); WT_ALS <- ln_from(18.1,8.5)
EV_ISS <- et(amt=17.5,dur=1); EV_ALS <- et(amt=15,ii=6,addl=4,dur=2)
T_ISS <- c(0.5,1,1.25,1.5,2,3,4,5,6); T_ALS <- c(21,23.5)
UN <- c("lcl","lvc","lvp","lq","clwt","vpwt","qwt","oCL","oVc","sig")
AN <- c("lcl","lv","oCL","oV","sig")
TH_ISS <- c(lcl=log(0.16),lvc=log(3.86),lvp=log(0.19),lq=log(0.13),
            clwt=0.97,vpwt=1.07,qwt=1.19,oCL=sqrt(0.062),oVc=sqrt(0.048),sig=0.12)
TH_ALS <- c(lcl=log(2.99),lv=log(9.55),oCL=sqrt(0.0223),oV=sqrt(0.0134),sig=0.119)
SE_ISS <- c(lcl=.15,lvc=.15,lvp=.15,lq=.15,clwt=.12,vpwt=.12,qwt=.12,
            oCL=.20*sqrt(0.062),oVc=.20*sqrt(0.048),sig=.20*0.12)
SE_ALS <- c(lcl=.15,lv=.15,oCL=.20*sqrt(0.0223),oV=.20*sqrt(0.0134),sig=.20*0.119)

mk_uni <- function(p) eval(parse(text=sprintf('function(){
  ini({lcl<-%.8f; lvc<-%.8f; lvp<-%.8f; lq<-%.8f
       clwt<-%.8f; vpwt<-%.8f; qwt<-%.8f
       eta.cl~%.8f; eta.vc~%.8f; prop.err<-%.8f})
  model({cl<-exp(lcl+eta.cl)*WT^clwt; vc<-exp(lvc+eta.vc)
         vp<-exp(lvp)*WT^vpwt; q<-exp(lq)*WT^qwt; f(centr)<-WT
         d/dt(centr)<--(cl+q)/vc*centr+q/vp*peri
         d/dt(peri)<-q/vc*centr-q/vp*peri
         cp<-centr/vc; cp~prop(prop.err)})}',
  p[["lcl"]],p[["lvc"]],p[["lvp"]],p[["lq"]],p[["clwt"]],p[["vpwt"]],p[["qwt"]],
  max(p[["oCL"]],1e-4)^2, max(p[["oVc"]],1e-4)^2, max(abs(p[["sig"]]),1e-4))))
mod_als <- function(){ ini({lcl<-log(2.99); lv<-log(9.55)
    eta.cl~0.0223; eta.v~0.0134; prop.err<-0.119})
  model({cl<-exp(lcl+eta.cl)*(WT/20); v<-exp(lv+eta.v)*(WT/20); f(centr)<-WT
    d/dt(centr)<--cl/v*centr; cp<-centr/v; cp~prop(prop.err)}) }
grabA <- function(f){ fx<-nlmixr2est::fixef(f); om<-f$omega
  c(lcl=fx[["lcl"]],lv=fx[["lv"]],oCL=sqrt(om[["eta.cl","eta.cl"]]),
    oV=sqrt(om[["eta.v","eta.v"]]),sig=fx[["prop.err"]]) }

## g_Alsultan: Alsultan published NO exponent -> MARGINAL construction (one block)
g_als <- function(phi){
  gg <- datagen(list(A=list(model=mk_uni(phi), times=T_ALS, ev=EV_ALS, n=72L,
                            cov_dist=list(WT=WT_ALS))),
                control=datagenControl(n_sim=20000L, seed=99L))
  f <- suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(mod_als, admData(),
        est="adgh", control=adghControl(studies=gg, n_nodes=7L, print=0L,
                                        covMethod="none", maxeval=2000L))))
  grabA(f) }
say(sprintf("[%3.0fs] g_Alsultan at Ayuthaya's parameters:", proc.time()[["elapsed"]]-T0))
g0 <- g_als(TH_ISS)
say(sprintf("   induced   %s", paste(sprintf("%s %.3f",AN,g0), collapse="  ")))
say(sprintf("   PUBLISHED %s", paste(sprintf("%s %.3f",AN,TH_ALS[AN]), collapse="  ")))
say("   note: Alsultan's vector contains NO exponent -> it cannot vote one to zero")

res <- function(phi) c((TH_ISS-phi)/SE_ISS, (TH_ALS[AN]-g_als(phi))/SE_ALS[AN])
FD <- setNames(rep(8e-3,length(UN)),UN); FD[c("clwt","vpwt","qwt")] <- 1.5e-2
phi <- TH_ISS; bs <- Inf; bp <- phi
for (it in 1:4) {
  r0 <- res(phi); if (sum(r0^2)<bs){bs<-sum(r0^2); bp<-phi}
  J <- vapply(UN, function(q){ph<-phi; ph[[q]]<-ph[[q]]+FD[[q]]; (res(ph)-r0)/FD[[q]]},
              numeric(length(r0)))
  phi <- phi - 0.6*as.numeric(solve(crossprod(J)+diag(1e-7,length(UN)),crossprod(J,r0)))
  say(sprintf("[%4.0fs] GN %d SSR %8.3f | clwt %.3f vpwt %.3f qwt %.3f",
      proc.time()[["elapsed"]]-T0, it, sum(r0^2), phi[["clwt"]],phi[["vpwt"]],phi[["qwt"]])) }
PLL <- bp
AGGv <- readRDS(file.path(OUT,"page-fw.rds"))
say("\n           Ayuthaya pub   AGG (joint)   PARAM-LL")
for (q in c("clwt","vpwt","qwt")) say(sprintf("%-5s %13.3f %13.3f %11.3f",
  q, TH_ISS[[q]], AGGv$joint[[q]], PLL[[q]]))
saveRDS(list(PLL=PLL,TH_ISS=TH_ISS,TH_ALS=TH_ALS,g0=g0,AGG=AGGv), file.path(OUT,"page-pll.rds"))

## ---- the weight-gradient dosing plot ---------------------------------------
## Cmax (21 h) and trough (23.5 h) against weight, on the Alsultan design.
WTG <- seq(6, 40, by = 1)
prof <- function(mod, lab){
  do.call(rbind, lapply(WTG, function(w){
    g <- datagen(list(S=list(model=mod, times=T_ALS, ev=EV_ALS, n=1L,
                             cov_dist=list(WT=list(values=w)))),
                 control=datagenControl(n_sim=4000L, seed=5L))
    data.frame(wt=w, model=lab, Cmax=g$S$E[1], Trough=g$S$E[2]) })) }
say(sprintf("[%4.0fs] building the weight-gradient panel ...", proc.time()[["elapsed"]]-T0))
mod_iss_pub <- mk_uni(TH_ISS); mod_pll <- mk_uni(PLL)
df <- rbind(prof(mod_iss_pub,"Ayuthaya (published, 2C)"),
            prof(mod_als,"Alsultan (published, 1C)"),
            prof(mod_pll,"parameter likelihood (2C)"))
dl <- rbind(data.frame(df[,c("wt","model")], y=df$Cmax, panel="Cmax (21 h)"),
            data.frame(df[,c("wt","model")], y=df$Trough, panel="Trough (23.5 h)"))
p <- ggplot(dl, aes(wt, y, colour=model)) + geom_line(linewidth=1.1) +
  facet_wrap(~panel, scales="free_y") +
  scale_colour_manual(values=c("Ayuthaya (published, 2C)"="#1b6ca8",
    "Alsultan (published, 1C)"="#c1121f","parameter likelihood (2C)"="#2a9d8f")) +
  labs(x="body weight (kg)", y="concentration (mg/L)", colour=NULL,
       title="Vancomycin exposure vs weight, Alsultan design (15 mg/kg q6h, 2 h infusion)",
       subtitle="Alsultan's model is exactly flat: mg/kg dosing with CL and V both linear in weight") +
  theme_bw(base_size=11) + theme(legend.position="bottom")
ggsave(file.path(OUT,"page-dosing.png"), p, width=10, height=4.6, dpi=150)
say(sprintf("total %.0f s", proc.time()[["elapsed"]]-T0))
