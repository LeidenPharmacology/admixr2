## Can the parameter likelihood combine models with DIFFERENT STRUCTURE?
##
## Truth: a 2-compartment drug. Three publications, ALL of them 1-compartment,
## each fitted to its own sampling window and its own population. No data, only
## their parameter tables + RSEs + covariate summary statistics.
##
## THE MECHANISM. The induced map g_s(phi) -- "what would study s have published
## if the truth were phi" -- depends on the study's DESIGN. Fit a 1-cmt model to
## 2-cmt data sampled EARLY and its V lands near Vc; sample LATE and its V lands
## near Vss = Vc + Vp. So the same 2-cmt truth forces DIFFERENT 1-cmt numbers out
## of differently-designed studies, and that disagreement is exactly what
## identifies Q and Vp -- parameters no source model contains.
OUT <- "C:/Users/hidde/AppData/Local/Temp/claude/C--package-admixr2/3ff305c7-64a5-4cb0-ba61-1436e2e9b16e/scratchpad"
say <- function(...) { cat(..., "\n"); utils::flush.console() }; T0 <- proc.time()[["elapsed"]]

gh <- function(k) { i <- seq_len(k-1L)
  J <- matrix(0,k,k); J[cbind(i,i+1L)] <- sqrt(i); J[cbind(i+1L,i)] <- sqrt(i)
  e <- eigen(J, symmetric=TRUE); o <- order(e$values)
  list(x = e$values[o], w = e$vectors[1L,o]^2) }
GE <- gh(5L); GX <- gh(5L)

TRUE2 <- c(lcl=log(4.0), lvc=log(30), lq=log(8.0), lvp=log(40),
           bCL=0.75, bVc=1.00, oCL=0.25, oVc=0.20, sig=0.10)
P1NM  <- c("lcl","lv","bCL","bV","oCL","oV","sig")     # 1-cmt publication

## three 1-cmt publications: same drug, different sampling WINDOWS + populations
SRC <- list(
  early = list(times=c(.25,.5,1,2,4),  dose=500, mu_x=log(68/70), sd_x=.18, n= 45),
  mid   = list(times=c(1,4,8,12),      dose=500, mu_x=log(80/70), sd_x=.22, n= 90),
  late  = list(times=c(8,12,24,36),    dose=500, mu_x=log(95/70), sd_x=.20, n= 60))
SE <- c(lcl=.05, lv=.06, bCL=.08, bV=.09, oCL=.03, oV=.03, sig=.010)

grid_of <- function(s) {
  xs <- s$mu_x + s$sd_x * GX$x
  g <- expand.grid(i=seq_along(xs), j=seq_along(GE$x), k=seq_along(GE$x),
                   KEEP.OUT.ATTRS=FALSE)
  list(x=xs[g$i], z1=GE$x[g$j], z2=GE$x[g$k],
       w=(GX$w[g$i]*GE$w[g$j]*GE$w[g$k]) / sum(GX$w[g$i]*GE$w[g$j]*GE$w[g$k]))
}
conc2 <- function(p, gr, tt, D) {            # 2-cmt IV bolus, biexponential
  cl <- exp(p[["lcl"]] + p[["bCL"]]*gr$x + p[["oCL"]]*gr$z1)
  vc <- exp(p[["lvc"]] + p[["bVc"]]*gr$x + p[["oVc"]]*gr$z2)
  q  <- exp(p[["lq"]]); vp <- exp(p[["lvp"]])
  k10 <- cl/vc; k12 <- q/vc; k21 <- q/vp
  s <- k10+k12+k21; d <- sqrt(pmax(s^2 - 4*k10*k21, 0))
  al <- (s+d)/2; be <- (s-d)/2
  A <- (D/vc)*(al-k21)/(al-be); B <- (D/vc)*(k21-be)/(al-be)
  A*exp(-outer(al,tt)) + B*exp(-outer(be,tt))
}
conc1 <- function(p, gr, tt, D) {            # 1-cmt IV bolus
  cl <- exp(p[["lcl"]] + p[["bCL"]]*gr$x + p[["oCL"]]*gr$z1)
  v  <- exp(p[["lv"]]  + p[["bV"]] *gr$x + p[["oV"]] *gr$z2)
  outer(D/v, tt, function(a,b) a) * exp(-outer(cl/v, tt))
}
mom <- function(f, gr, sig) {
  mu <- as.numeric(crossprod(gr$w, f)); fc <- sweep(f,2L,mu)
  V  <- crossprod(fc, fc*gr$w)
  diag(V) <- diag(V) + sig^2 * as.numeric(crossprod(gr$w, f^2))
  list(E=mu, V=V)
}
nll <- function(obs, pr, n) {
  ch <- tryCatch(chol(pr$V), error=function(e) NULL); if (is.null(ch)) return(1e12)
  iv <- chol2inv(ch); r <- obs$E - pr$E
  n*(2*sum(log(diag(ch))) + sum(iv*obs$V) + as.numeric(t(r)%*%iv%*%r))
}
## fit a 1-cmt model to a source's aggregate moments
fit1 <- function(obs, s, gr, start) {
  o <- stats::optim(start, function(v){ p <- v
    p[c("oCL","oV","sig")] <- abs(p[c("oCL","oV","sig")])
    nll(obs, mom(conc1(p,gr,s$times,s$dose), gr, p[["sig"]]), s$n) },
    method="Nelder-Mead", control=list(maxit=4000, reltol=1e-11))
  o <- stats::optim(o$par, function(v){ p <- v
    p[c("oCL","oV","sig")] <- abs(p[c("oCL","oV","sig")])
    nll(obs, mom(conc1(p,gr,s$times,s$dose), gr, p[["sig"]]), s$n) },
    method="Nelder-Mead", control=list(maxit=4000, reltol=1e-12))
  p <- o$par; p[c("oCL","oV","sig")] <- abs(p[c("oCL","oV","sig")]); p
}
ST1 <- c(lcl=log(3), lv=log(35), bCL=.6, bV=.9, oCL=.3, oV=.25, sig=.15)
GR  <- lapply(SRC, grid_of)

## g_s(phi): what source s would publish if the truth were the 2-cmt model phi
g_s <- function(phi, s) {
  o <- mom(conc2(phi, GR[[s]], SRC[[s]]$times, SRC[[s]]$dose), GR[[s]], phi[["sig"]])
  fit1(o, SRC[[s]], GR[[s]], ST1)
}
say("computing what each 1-cmt study would publish under the 2-cmt truth ...")
PUB <- lapply(names(SRC), function(s) g_s(TRUE2, s)); names(PUB) <- names(SRC)
for (s in names(SRC)) say(sprintf("  %-6s CL %.3f  V %.1f  (Vc=%.0f Vss=%.0f)", s,
  exp(PUB[[s]][["lcl"]]), exp(PUB[[s]][["lv"]]),
  exp(TRUE2[["lvc"]]), exp(TRUE2[["lvc"]])+exp(TRUE2[["lvp"]])))

## ---- the parameter likelihood ----------------------------------------------
UN <- names(TRUE2)
resid <- function(phi) unlist(lapply(names(SRC), function(s) (PUB[[s]]-g_s(phi,s))/SE))
FD <- c(lcl=5e-3,lvc=5e-3,lq=1e-2,lvp=1e-2,bCL=1e-2,bVc=1e-2,oCL=5e-3,oVc=5e-3,sig=4e-3)
phi <- TRUE2; phi[c("lq","lvp")] <- c(log(3.0), log(80))   # deliberately wrong start
phi[c("lcl","lvc")] <- c(log(3.2), log(38))
for (it in 1:6) {
  r0 <- resid(phi)
  J <- vapply(UN, function(q){ ph <- phi; ph[[q]] <- ph[[q]]+FD[[q]]
        (resid(ph)-r0)/FD[[q]] }, numeric(length(r0)))
  st <- tryCatch(solve(crossprod(J)+diag(1e-8,length(UN)), crossprod(J,r0)),
                 error=function(e) rep(0,length(UN)))
  phi <- phi - 0.8*as.numeric(st)
  say(sprintf("[%3.0fs] GN %d SSR %8.3f | CL %.2f Vc %.1f Q %.2f Vp %.1f bCL %.3f",
      proc.time()[["elapsed"]]-T0, it, sum(r0^2), exp(phi[["lcl"]]), exp(phi[["lvc"]]),
      exp(phi[["lq"]]), exp(phi[["lvp"]]), phi[["bCL"]]))
}
SEphi <- sqrt(diag(solve(crossprod(J))))
say("\n            truth   estimate      SE")
for (q in UN) say(sprintf("%-6s %9.4f %10.4f %7.4f", q, TRUE2[[q]], phi[[q]], SEphi[[q]]))
say(sprintf("\nnatural: CL %.2f (%.2f)  Vc %.1f (%.1f)  Q %.2f (%.2f)  Vp %.1f (%.1f)",
  exp(phi[["lcl"]]), exp(TRUE2[["lcl"]]), exp(phi[["lvc"]]), exp(TRUE2[["lvc"]]),
  exp(phi[["lq"]]),  exp(TRUE2[["lq"]]),  exp(phi[["lvp"]]), exp(TRUE2[["lvp"]])))

## ---- the induced-map curve, for the figure ---------------------------------
say("\ntracing the induced map over Q ...")
QG <- exp(seq(log(2), log(20), length.out = 11))
MAP <- vapply(QG, function(qq) { ph <- TRUE2; ph[["lq"]] <- log(qq)
  vapply(names(SRC), function(s) exp(g_s(ph, s)[["lv"]]), 0) }, numeric(3))
saveRDS(list(TRUE2=TRUE2, PUB=PUB, phi=phi, SE=SEphi, SRC=SRC, QG=QG, MAP=MAP,
             GR=GR, conc2=conc2, conc1=conc1, mom=mom, SEs=SE),
        file.path(OUT,"param-1cmt2cmt.rds"))
say(sprintf("total %.0f s", proc.time()[["elapsed"]]-T0))
