## PAGE cross-prediction panels: each study's DESIGN, every model on it.
##
## The point of the panel is which model describes which data. A model applied to
## a design it was not developed on is the honest test, and the synthesis methods
## are judged by whether they describe BOTH.
suppressMessages({library(rxode2); library(ggplot2)})
suppressMessages(devtools::load_all(".",quiet=TRUE))
OUT <- "C:/Users/hidde/AppData/Local/Temp/claude/C--package-admixr2/3ff305c7-64a5-4cb0-ba61-1436e2e9b16e/scratchpad"
say <- function(...) { cat(...,"\n"); utils::flush.console() }; T0 <- proc.time()[["elapsed"]]
ln_from <- function(mu,sd) list(meanlog=log(mu^2/sqrt(sd^2+mu^2)), sdlog=sqrt(log(1+sd^2/mu^2)))
WT_ISS <- ln_from(20.2,9.2); WT_ALS <- ln_from(18.1,8.5)
EV_ISS <- et(amt=17.5,dur=1); EV_ALS <- et(amt=15,ii=6,addl=4,dur=2)
T_ISS <- c(0.5,1,1.25,1.5,2,3,4,5,6); T_ALS <- c(21,23.5)
UN <- c("lcl","lvc","lvp","lq","clwt","vpwt","qwt","oCL","oVc","sig")
TH_ISS <- c(lcl=log(0.16),lvc=log(3.86),lvp=log(0.19),lq=log(0.13),
            clwt=0.97,vpwt=1.07,qwt=1.19,oCL=sqrt(0.062),oVc=sqrt(0.048),sig=0.12)
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
mod_iss <- mk_uni(TH_ISS)

## ---- refit AGG to recover its FULL parameter vector -------------------------
say(sprintf("[%3.0fs] refitting the aggregate joint model ...", proc.time()[["elapsed"]]-T0))
gen <- datagen(list(
  Ayuthaya=list(model=mod_iss,times=T_ISS,ev=EV_ISS,n=14L,cov_dist=list(WT=WT_ISS)),
  Alsultan=list(model=mod_als,times=T_ALS,ev=EV_ALS,n=72L,cov_dist=list(WT=WT_ALS))),
  control=datagenControl(n_sim=40000L,seed=12345L))
mod_start <- mk_uni(c(lcl=log(0.30),lvc=log(2.50),lvp=log(0.40),lq=log(0.25),
                      clwt=0.60,vpwt=0.70,qwt=0.80,oCL=0.32,oVc=0.32,sig=0.20))
fA <- suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(mod_start,admData(),est="adgh",
  control=adghControl(studies=gen,n_nodes=7L,print=0L,covMethod="none",maxeval=3000L))))
s <- unlist(fA$env$admExtra$struct); om <- fA$omega
AGG <- c(lcl=s[["lcl"]],lvc=s[["lvc"]],lvp=s[["lvp"]],lq=s[["lq"]],
         clwt=s[["clwt"]],vpwt=s[["vpwt"]],qwt=s[["qwt"]],
         oCL=sqrt(om[["eta.cl","eta.cl"]]),oVc=sqrt(om[["eta.vc","eta.vc"]]),
         sig=nlmixr2est::fixef(fA)[["prop.err"]])
W <- readRDS(file.path(OUT,"page-weights.rds"))
PLL <- W$FITS$sens                                   # the reference weighting
say(sprintf("   AGG  clwt %.3f vpwt %.3f qwt %.3f", AGG[["clwt"]],AGG[["vpwt"]],AGG[["qwt"]]))
say(sprintf("   PLL  clwt %.3f vpwt %.3f qwt %.3f", PLL[["clwt"]],PLL[["vpwt"]],PLL[["qwt"]]))

## ---- profiles on BOTH designs ----------------------------------------------
TG_ISS <- seq(0.1, 8, by=0.1); TG_ALS <- seq(18.05, 24, by=0.05)
prof <- function(mod, tt, ev, wt) {
  g <- datagen(list(S=list(model=mod,times=tt,ev=ev,n=1L,cov_dist=list(WT=wt))),
               control=datagenControl(n_sim=8000L,seed=3L))
  as.numeric(g$S$E) }
MODS <- list("Ayuthaya published (2C)"=mod_iss, "Alsultan published (1C)"=mod_als,
             "aggregate joint (2C)"=mk_uni(AGG), "parameter likelihood (2C)"=mk_uni(PLL))
say(sprintf("[%3.0fs] building cross-prediction profiles ...", proc.time()[["elapsed"]]-T0))
df <- do.call(rbind, c(
  lapply(names(MODS), function(m) data.frame(t=TG_ISS, y=prof(MODS[[m]],TG_ISS,EV_ISS,WT_ISS),
    model=m, panel="Ayuthaya design: 17.5 mg/kg, 1 h infusion, single dose")),
  lapply(names(MODS), function(m) data.frame(t=TG_ALS, y=prof(MODS[[m]],TG_ALS,EV_ALS,WT_ALS),
    model=m, panel="Alsultan design: 15 mg/kg q6h, 2 h infusion, 5 doses"))))
obs <- rbind(
  data.frame(t=T_ISS, y=prof(mod_iss,T_ISS,EV_ISS,WT_ISS),
             panel="Ayuthaya design: 17.5 mg/kg, 1 h infusion, single dose"),
  data.frame(t=T_ALS, y=prof(mod_als,T_ALS,EV_ALS,WT_ALS),
             panel="Alsultan design: 15 mg/kg q6h, 2 h infusion, 5 doses"))
CO <- c("Ayuthaya published (2C)"="#1b6ca8","Alsultan published (1C)"="#c1121f",
        "aggregate joint (2C)"="#8a8f98","parameter likelihood (2C)"="#2a9d8f")
p <- ggplot(df, aes(t,y,colour=model)) +
  geom_line(linewidth=1.05) +
  geom_point(data=obs, aes(t,y), inherit.aes=FALSE, size=2.4, shape=21,
             fill="black", colour="white", stroke=0.9) +
  facet_wrap(~panel, scales="free") +
  scale_colour_manual(values=CO) +
  labs(x="time (h)", y="population mean concentration (mg/L)", colour=NULL,
       title="Cross-prediction: every model on both study designs",
       subtitle="Points mark each study's own sampling times, generated from its own published model") +
  theme_bw(base_size=11) + theme(legend.position="bottom")
ggsave(file.path(OUT,"page-cross.png"), p, width=11, height=4.8, dpi=150)
## RMS deviation from each study's own published model, on its own design
rms <- function(a,b) sqrt(mean(((a-b)/b)^2))*100
say("\nrelative RMS deviation from the study's OWN published model, on its own design (%)")
say(sprintf("%-28s %14s %14s","model","Ayuthaya design","Alsultan design"))
for (m in names(MODS)) say(sprintf("%-28s %13.1f%% %13.1f%%", m,
  rms(prof(MODS[[m]],TG_ISS,EV_ISS,WT_ISS), prof(mod_iss,TG_ISS,EV_ISS,WT_ISS)),
  rms(prof(MODS[[m]],TG_ALS,EV_ALS,WT_ALS), prof(mod_als,TG_ALS,EV_ALS,WT_ALS))))
say(sprintf("total %.0f s", proc.time()[["elapsed"]]-T0))
