## Model averaging / MBMA: aggregate likelihood vs parameter likelihood.
## Harness validated against admixr2 in mc-averaging-validate.R (4.9e-04).
setwd("C:/package/admixr2/.claude/worktrees/feature-covariate-quadrature")
source("validation/mc-averaging-harness.R")
OUT <- "C:/Users/hidde/AppData/Local/Temp/claude/C--package-admixr2/3ff305c7-64a5-4cb0-ba61-1436e2e9b16e/scratchpad"
say <- function(...) { cat(..., "\n"); utils::flush.console() }; T0 <- proc.time()[["elapsed"]]

TRUTH <- c(lcl = log(4.2), lv = log(48), bCL = 0.75, bV = 1.00,
           oCL = 0.28, oV = 0.20, sig = 0.13)
K <- 4L
FULL <- PAR_NM; NOBV <- setdiff(PAR_NM,"bV"); NOCOV <- setdiff(PAR_NM,c("bCL","bV"))
SRC <- list(
  A = list(mu_x=log(68/70), sd_x=0.17, n= 40, times=c(.5,1,2,4,8), dose=500, free=FULL,  fix=NULL),
  B = list(mu_x=log(84/70), sd_x=0.22, n=150, times=c(10,12),      dose=750, free=NOBV,  fix=c(bV=0)),
  C = list(mu_x=log(61/70), sd_x=0.14, n= 60, times=c(.5,1,2,4,8), dose=500, free=FULL,  fix=NULL),
  D = list(mu_x=log(92/70), sd_x=0.24, n=200, times=c(2,8),        dose=750, free=NOCOV, fix=c(bCL=0,bV=0)),
  E = list(mu_x=log(105/70),sd_x=0.16, n= 30, times=c(.5,1,2,4,8), dose=500, free=NOBV,  fix=c(bV=0)))
## published SEs. ANTI: small rich studies precise, big sparse ones loose.
SE_ANTI <- list(
  A=c(lcl=.04,lv=.05,bCL=.05,bV=.07,oCL=.03,oV=.03,sig=.012),
  B=c(lcl=.05,lv=.06,bCL=.20,        oCL=.03,oV=.04,sig=.010),
  C=c(lcl=.04,lv=.05,bCL=.06,bV=.08,oCL=.03,oV=.03,sig=.013),
  D=c(lcl=.06,lv=.07,                oCL=.04,oV=.04,sig=.009),
  E=c(lcl=.05,lv=.06,bCL=.06,        oCL=.04,oV=.04,sig=.014))
## ALIGNED: precision follows n, as a naive analyst would assume.
SE_ALIGN <- lapply(names(SRC), function(s) {
  v <- SE_ANTI[[s]]; v[] <- 0.10 * sqrt(120 / SRC[[s]]$n); v })
names(SE_ALIGN) <- names(SRC)

START <- c(lcl=log(3),lv=log(60),bCL=.4,bV=.6,oCL=.4,oV=.3,sig=.2)
blocks_of <- function(p, s) make_blocks(p, SRC[[s]], K)
## pseudo-true publication: fit source s's OWN form to truth-generated blocks
say("computing pseudo-true publications ...")
PSEUDO <- lapply(names(SRC), function(s)
  fit_form(blocks_of(TRUTH, s), SRC[[s]]$free, START, SRC[[s]]$fix))
names(PSEUDO) <- names(SRC)
for (s in names(SRC)) say(sprintf("  %s  %s", s,
  paste(sprintf("%s %.3f", SRC[[s]]$free, PSEUDO[[s]][SRC[[s]]$free]), collapse="  ")))

## ---------------- estimators -------------------------------------------------
agg_fit <- function(pub) {                      # AGGREGATE: blocks from each pub
  bl <- unlist(lapply(names(SRC), function(s) blocks_of(pub[[s]], s)), recursive = FALSE)
  fit_form(bl, FULL, START)
}
g_map <- function(phi, s)                       # induced map
  fit_form(blocks_of(phi, s), SRC[[s]]$free, phi, SRC[[s]]$fix)[SRC[[s]]$free]
par_resid <- function(phi, pub, SE) unlist(lapply(names(SRC), function(s)
  (pub[[s]][SRC[[s]]$free] - g_map(phi, s)) / SE[[s]][SRC[[s]]$free]))

## ---------------- linearisation (once per cell) ------------------------------
## Both estimators are smooth maps from the published numbers to phi_hat, so both
## are linearised about the truth and treated IDENTICALLY. Validated against full
## nonlinear fits on a subsample below.
FDp <- c(lcl=5e-3,lv=5e-3,bCL=1e-2,bV=1e-2,oCL=5e-3,oV=5e-3,sig=4e-3)
build_jac <- function(SE) {
  r0 <- par_resid(TRUTH, PSEUDO, SE)
  Jp <- vapply(FULL, function(q) { ph <- TRUTH; ph[[q]] <- ph[[q]] + FDp[[q]]
    (par_resid(ph, PSEUDO, SE) - r0)/FDp[[q]] }, numeric(length(r0)))
  a0 <- agg_fit(PSEUDO)
  Ja <- do.call(cbind, lapply(names(SRC), function(s)
    vapply(SRC[[s]]$free, function(q) { pb <- PSEUDO; h <- FDp[[q]]
      pb[[s]][[q]] <- pb[[s]][[q]] + h; (agg_fit(pb) - a0)/h }, numeric(7))))
  list(Jp = Jp, r0 = r0, a0 = a0, Ja = Ja)
}
flat <- function(pub) unlist(lapply(names(SRC), function(s) pub[[s]][SRC[[s]]$free]))
SEvec <- function(SE) unlist(lapply(names(SRC), function(s) SE[[s]][SRC[[s]]$free]))

## ---------------- metrics ----------------------------------------------------
REG <- list(INTERP   = list(mu_x=log(78/70), sd_x=0.20),
            `EXTRAP-LO` = list(mu_x=log(45/70), sd_x=0.15),
            `EXTRAP-HI` = list(mu_x=log(125/70),sd_x=0.18))
NEWD <- list(times = c(1,3,6,12,24), dose = 600)     # a regimen nobody ran
div <- function(p, rg) {                             # -2LL of truth's moments
  gr <- grid_of(rg$mu_x, rg$sd_x); tr <- moments(TRUTH, gr, NEWD$times, NEWD$dose)
  pr <- moments(p, gr, NEWD$times, NEWD$dose)
  ch <- tryCatch(chol(pr$V), error=function(e) NULL); if (is.null(ch)) return(NA_real_)
  iv <- chol2inv(ch); r <- tr$E - pr$E
  2*sum(log(diag(ch))) + sum(iv*tr$V) + as.numeric(t(r)%*%iv%*%r)
}
dose_err <- function(p, rg)                          # % error in dose for target AUC
  100*(exp((p[["lcl"]]-TRUTH[["lcl"]]) + (p[["bCL"]]-TRUTH[["bCL"]])*rg$mu_x) - 1)
score <- function(p) c(
  vapply(names(REG), function(r) div(p, REG[[r]]) - div(TRUTH, REG[[r]]), 0),
  stats::setNames(vapply(names(REG), function(r) dose_err(p, REG[[r]]), 0),
                  paste0("dose.", names(REG))),
  bCL = p[["bCL"]] - TRUTH[["bCL"]], bV = p[["bV"]] - TRUTH[["bV"]])

## ---------------- the study --------------------------------------------------
CELLS <- list(`anti/tau=0`    = list(SE=SE_ANTI,  tau=0.00, disc=FALSE),
              `anti/tau=.10`  = list(SE=SE_ANTI,  tau=0.10, disc=FALSE),
              `aligned/tau=0` = list(SE=SE_ALIGN, tau=0.00, disc=FALSE),
              `discrepant`    = list(SE=SE_ANTI,  tau=0.00, disc=TRUE))
R <- 600L; set.seed(20260815); RES <- list()
for (cn in names(CELLS)) {
  cl <- CELLS[[cn]]; say(sprintf("[%4.0fs] cell %s: jacobians ...", proc.time()[["elapsed"]]-T0, cn))
  JJ <- build_jac(cl$SE)
  Ap <- solve(crossprod(JJ$Jp) + diag(1e-10, 7))     # PAR: GN solve
  sev <- SEvec(cl$SE); mu0 <- flat(PSEUDO)
  ## BEST-SINGLE: the source whose solo projection is closest to truth
  solo <- lapply(names(SRC), function(s) fit_form(blocks_of(PSEUDO[[s]], s), FULL, START))
  names(solo) <- names(SRC)
  best <- names(SRC)[which.min(vapply(solo, function(p) div(p, REG$INTERP), 0))]
  out <- array(NA_real_, c(R, 5, 8),
               dimnames=list(NULL, c("AGG","PAR-FE","PAR-RE","NAIVE","BEST-SINGLE"), names(score(TRUTH))))
  for (r in seq_len(R)) {
    b <- if (cl$tau > 0) unlist(lapply(names(SRC), function(s)
           rep(stats::rnorm(1,0,cl$tau), length(SRC[[s]]$free)))) else 0
    y <- mu0 + b + stats::rnorm(length(mu0), 0, sev)
    if (cl$disc) { i <- grep("^D\\.|^D$", names(y)); y[startsWith(names(y),"lcl")][4] <- y[startsWith(names(y),"lcl")][4] }
    yy <- y; if (cl$disc) yy[which(names(mu0)=="lcl")[4]] <- yy[which(names(mu0)=="lcl")[4]] + 0.45
    d  <- yy - mu0
    phiP <- TRUTH - as.numeric(Ap %*% crossprod(JJ$Jp, JJ$r0 + (-d/sev)))
    W  <- 1/(sev^2); tau2 <- max(0, (sum(W*(d - 0)^2) - length(d))/sum(W))
    w2 <- 1/(sev^2 + tau2)
    phiR <- TRUTH - as.numeric(solve(crossprod(JJ$Jp, (w2*sev^2)*JJ$Jp) + diag(1e-10,7),
                                     crossprod(JJ$Jp, (w2*sev^2)*(JJ$r0 + (-d/sev)))))
    phiA <- JJ$a0 + as.numeric(JJ$Ja %*% d)
    nv <- TRUTH; for (q in c("bCL","bV")) { k <- which(names(mu0)==q)
      nv[[q]] <- if (length(k)) sum(yy[k]/sev[k]^2)/sum(1/sev[k]^2) else TRUTH[[q]] }
    nv[c("lcl","lv","oCL","oV","sig")] <- vapply(c("lcl","lv","oCL","oV","sig"), function(q){
      k <- which(names(mu0)==q); sum(yy[k]/sev[k]^2)/sum(1/sev[k]^2)}, 0)
    out[r,"AGG",] <- score(phiA); out[r,"PAR-FE",] <- score(phiP)
    out[r,"PAR-RE",] <- score(phiR); out[r,"NAIVE",] <- score(nv)
    out[r,"BEST-SINGLE",] <- score(solo[[best]])
  }
  RES[[cn]] <- list(out = out, best = best)
  say(sprintf("[%4.0fs]   done (best single = %s)", proc.time()[["elapsed"]]-T0, best))
  saveRDS(list(RES=RES, TRUTH=TRUTH, PSEUDO=PSEUDO, SRC=SRC, REG=REG), file.path(OUT,"mc-averaging.rds"))
}
say(sprintf("total %.0f s", proc.time()[["elapsed"]]-T0))
