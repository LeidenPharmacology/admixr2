## THE DECISIVE TEST for the original motivation.
##
## The parameter likelihood was built on the idea that a published coefficient is
## WITHIN-study knowledge (estimated from that study's individual data), while
## aggregate MBMA sees only BETWEEN-study contrast. To find out whether that
## leverage is real, kill the between-study channel: give every source the SAME
## covariate distribution, changing nothing else.
##
##   AGG  can then reach bCL only through the within-study SPREAD channel, which
##        trades off against omega -> expected to degrade badly.
##   PAR  still has five published bCL as direct data -> expected to barely notice.
##
## If PAR does NOT hold up here, the original motivation is not what makes it work.
setwd("C:/package/admixr2/.claude/worktrees/feature-covariate-quadrature")
source("validation/mc-averaging-harness.R")
OUT <- "C:/Users/hidde/AppData/Local/Temp/claude/C--package-admixr2/3ff305c7-64a5-4cb0-ba61-1436e2e9b16e/scratchpad"
say <- function(...) { cat(..., "\n"); utils::flush.console() }; T0 <- proc.time()[["elapsed"]]

TRUTH <- c(lcl=log(4.2), lv=log(48), bCL=0.75, bV=1.00, oCL=0.28, oV=0.20, sig=0.13)
K <- 4L
FULL <- PAR_NM; NOBV <- setdiff(PAR_NM,"bV"); NOCOV <- setdiff(PAR_NM,c("bCL","bV"))
BASE <- list(
  A = list(mu_x=log(68/70), sd_x=0.17, n= 40, times=c(.5,1,2,4,8), dose=500, free=FULL,  fix=NULL),
  B = list(mu_x=log(84/70), sd_x=0.22, n=150, times=c(10,12),      dose=750, free=NOBV,  fix=c(bV=0)),
  C = list(mu_x=log(61/70), sd_x=0.14, n= 60, times=c(.5,1,2,4,8), dose=500, free=FULL,  fix=NULL),
  D = list(mu_x=log(92/70), sd_x=0.24, n=200, times=c(2,8),        dose=750, free=NOCOV, fix=c(bCL=0,bV=0)),
  E = list(mu_x=log(105/70),sd_x=0.16, n= 30, times=c(.5,1,2,4,8), dose=500, free=NOBV,  fix=c(bV=0)))
## SAME population for every source: the only change.
SAME <- lapply(BASE, function(s) { s$mu_x <- log(78/70); s$sd_x <- 0.20; s })
SE <- list(A=c(lcl=.04,lv=.05,bCL=.05,bV=.07,oCL=.03,oV=.03,sig=.012),
           B=c(lcl=.05,lv=.06,bCL=.20,        oCL=.03,oV=.04,sig=.010),
           C=c(lcl=.04,lv=.05,bCL=.06,bV=.08,oCL=.03,oV=.03,sig=.013),
           D=c(lcl=.06,lv=.07,                oCL=.04,oV=.04,sig=.009),
           E=c(lcl=.05,lv=.06,bCL=.06,        oCL=.04,oV=.04,sig=.014))
START <- c(lcl=log(3),lv=log(60),bCL=.4,bV=.6,oCL=.4,oV=.3,sig=.2)
FDp <- c(lcl=5e-3,lv=5e-3,bCL=1e-2,bV=1e-2,oCL=5e-3,oV=5e-3,sig=4e-3)

run_cell <- function(SRC, lab) {
  bo <- function(p,s) make_blocks(p, SRC[[s]], K)
  PS <- lapply(names(SRC), function(s) fit_form(bo(TRUTH,s), SRC[[s]]$free, START, SRC[[s]]$fix))
  names(PS) <- names(SRC)
  agg <- function(pub) fit_form(unlist(lapply(names(SRC), function(s) bo(pub[[s]],s)),
                                       recursive=FALSE), FULL, START)
  gm  <- function(phi,s) fit_form(bo(phi,s), SRC[[s]]$free, phi, SRC[[s]]$fix)[SRC[[s]]$free]
  rs  <- function(phi,pub) unlist(lapply(names(SRC), function(s)
           (pub[[s]][SRC[[s]]$free] - gm(phi,s)) / SE[[s]][SRC[[s]]$free]))
  r0 <- rs(TRUTH, PS)
  Jp <- vapply(FULL, function(q){ ph<-TRUTH; ph[[q]]<-ph[[q]]+FDp[[q]]
         (rs(ph,PS)-r0)/FDp[[q]] }, numeric(length(r0)))
  a0 <- agg(PS)
  Ja <- do.call(cbind, lapply(names(SRC), function(s)
        vapply(SRC[[s]]$free, function(q){ pb<-PS; h<-FDp[[q]]
          pb[[s]][[q]] <- pb[[s]][[q]]+h; (agg(pb)-a0)/h }, numeric(7))))
  sev <- unlist(lapply(names(SRC), function(s) SE[[s]][SRC[[s]]$free]))
  Ap  <- solve(crossprod(Jp) + diag(1e-10,7))
  say(sprintf("[%4.0fs] %-12s  AGG at truth-pubs: bCL %.3f bV %.3f | cond(Ja) %.1e",
      proc.time()[["elapsed"]]-T0, lab, a0[["bCL"]], a0[["bV"]], kappa(Ja)))
  set.seed(4242); R <- 800L
  M <- array(NA_real_, c(R,2,2), dimnames=list(NULL,c("AGG","PAR-FE"),c("bCL","bV")))
  for (r in seq_len(R)) {
    d <- stats::rnorm(length(sev), 0, sev)
    pA <- a0 + as.numeric(Ja %*% d)
    pP <- TRUTH - as.numeric(Ap %*% crossprod(Jp, r0 + (-d/sev)))
    M[r,"AGG",]    <- c(pA[["bCL"]], pA[["bV"]]) - TRUTH[c("bCL","bV")]
    M[r,"PAR-FE",] <- c(pP[["bCL"]], pP[["bV"]]) - TRUTH[c("bCL","bV")]
  }
  list(M = M, a0 = a0, PS = PS)
}
res <- list(`differing pops` = run_cell(BASE, "differing"),
            `same pop`       = run_cell(SAME, "SAME POP"))
say("\n=================== between-study channel removed ===================")
say(sprintf("%-16s %-8s %9s %9s %9s %9s", "cell","method","bCL bias","bCL SD","bV bias","bV SD"))
for (cn in names(res)) for (m in c("AGG","PAR-FE")) {
  M <- res[[cn]]$M
  say(sprintf("%-16s %-8s %9.4f %9.4f %9.4f %9.4f", cn, m,
      mean(M[,m,"bCL"]), stats::sd(M[,m,"bCL"]),
      mean(M[,m,"bV"]),  stats::sd(M[,m,"bV"])))
}
say("\nSD inflation when the between-study contrast is removed:")
for (m in c("AGG","PAR-FE")) say(sprintf("  %-8s bCL %5.2fx   bV %5.2fx", m,
  stats::sd(res$`same pop`$M[,m,"bCL"])/stats::sd(res$`differing pops`$M[,m,"bCL"]),
  stats::sd(res$`same pop`$M[,m,"bV"]) /stats::sd(res$`differing pops`$M[,m,"bV"])))
saveRDS(res, file.path(OUT,"mc-nocontrast.rds"))
say(sprintf("total %.0f s", proc.time()[["elapsed"]]-T0))
