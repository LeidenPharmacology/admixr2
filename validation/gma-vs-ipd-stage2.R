## STAGE 2: unified model from the TWO PUBLICATIONS ONLY. No individual data.
##
## PROFILED AND OPTIMISED. The first version took >10 min. Three changes:
##  (1) ONE combined grid. All K covariate nodes share the same eta grid, so build
##      a single (K * n_eta) row block and compute f ONCE, then split -- instead of
##      K separate mom() calls. Same flops, ~K times less R interpretation.
##  (2) WARM START. The inner refit starts from the previous solution, which is
##      close because the outer step is small. This is the single biggest win.
##  (3) Coarser inner tolerance during the FD Jacobian, tight only for the final.
## PROF[] counters report where the time actually goes.
OUT <- "C:/Users/hidde/AppData/Local/Temp/claude/C--package-admixr2/3ff305c7-64a5-4cb0-ba61-1436e2e9b16e/scratchpad"
S <- readRDS(file.path(OUT,"gma-ipd-stage1.rds"))
say <- function(...) { cat(...,"\n"); utils::flush.console() }; T0 <- proc.time()[["elapsed"]]
PROF <- new.env(); PROF$mom <- 0; PROF$nmom <- 0L; PROF$inner <- 0; PROF$ninner <- 0L
tic <- function() proc.time()[["elapsed"]]

gh <- function(k){ i <- seq_len(k-1L); J <- matrix(0,k,k)
  J[cbind(i,i+1L)] <- sqrt(i); J[cbind(i+1L,i)] <- sqrt(i)
  e <- eigen(J,symmetric=TRUE); o <- order(e$values); list(x=e$values[o],w=e$vectors[1L,o]^2) }
NE <- 5L; NC <- 5L                          # eta nodes per dim; covariate nodes per dim
GE <- gh(NE); GX <- gh(NC)
TP <- S$TP; UN <- names(TP); ANM <- c("lcl","lvc","bCL1","bVc1","oCL","oVc","sig")
DA <- S$DA; DB <- S$DB
DA$times <- c(1,12); DA$tau <- 12; DB$times <- DB$obs; DB$tau <- NA

ssf <- function(e,tau) if (is.na(tau)) 1 else 1/(1-exp(-e*tau))

## ---- ONE combined grid per design: K blocks x M eta points ------------------
mkgrid <- function(d) {
  cg <- expand.grid(i=seq_len(NC), j=seq_len(NC), KEEP.OUT.ATTRS=FALSE)
  a1 <- d$mx1 + d$sx1*GX$x[cg$i]; a2 <- d$mx2 + d$sx2*GX$x[cg$j]
  wk <- GX$w[cg$i]*GX$w[cg$j]; wk <- wk/sum(wk)                 # block weights
  eg <- expand.grid(p=seq_len(NE), q=seq_len(NE), KEEP.OUT.ATTRS=FALSE)
  we <- GE$w[eg$p]*GE$w[eg$q]; we <- we/sum(we)                 # within-block weights
  K <- nrow(cg); M <- nrow(eg)
  idx <- rep(seq_len(K), each=M)
  list(K=K, M=M, wk=wk, we=we, blk=idx,
       x1=a1[idx], x2=a2[idx], z1=GE$x[eg$p][rep(seq_len(M),K)],
       z2=GE$x[eg$q][rep(seq_len(M),K)], n=d$n, d=d)
}
## all blocks' (E, V) in ONE pass
moments_all <- function(p, G, form) {
  t0 <- tic(); d <- G$d
  b2 <- if ("bCL2" %in% names(p)) p[["bCL2"]] else 0
  cl <- exp(p[["lcl"]]+p[["bCL1"]]*G$x1+b2*G$x2+abs(p[["oCL"]])*G$z1)
  vv <- exp(p[["lvc"]]+p[["bVc1"]]*G$x1+abs(p[["oVc"]])*G$z2)
  f <- if (form=="2") { q<-exp(p[["lq"]]); vp<-exp(p[["lvp"]])
        k10<-cl/vv;k12<-q/vv;k21<-q/vp; s<-k10+k12+k21
        dd<-sqrt(pmax(s^2-4*k10*k21,0)); al<-(s+dd)/2; be<-(s-dd)/2
        A<-(d$dose/vv)*(al-k21)/(al-be); Bc<-(d$dose/vv)*(k21-be)/(al-be)
        (A*ssf(al,d$tau))*exp(-outer(al,d$times))+(Bc*ssf(be,d$tau))*exp(-outer(be,d$times))
      } else { k<-cl/vv
        outer((d$dose/vv)*ssf(k,d$tau),d$times,function(a,b)a)*exp(-outer(k,d$times)) }
  sg2 <- abs(p[["sig"]])^2; nt <- length(d$times)
  out <- vector("list", G$K)
  for (k in seq_len(G$K)) {
    rows <- ((k-1L)*G$M+1L):(k*G$M); fk <- f[rows, , drop=FALSE]
    mu <- as.numeric(crossprod(G$we, fk)); fc <- sweep(fk,2L,mu)
    V <- crossprod(fc, fc*G$we)
    diag(V) <- diag(V) + sg2*as.numeric(crossprod(G$we, fk^2))
    out[[k]] <- list(E=mu, V=V)
  }
  PROF$mom <- PROF$mom + (tic()-t0); PROF$nmom <- PROF$nmom + 1L
  out
}
nll1 <- function(o,pr,n){ ch<-tryCatch(chol(pr$V),error=function(e) NULL)
  if(is.null(ch)) return(1e12); iv<-chol2inv(ch); r<-o$E-pr$E
  n*(2*sum(log(diag(ch)))+sum(iv*o$V)+as.numeric(t(r)%*%iv%*%r)) }

GRA <- mkgrid(DA); GRB <- mkgrid(DB)
WARM <- new.env(); WARM$A <- S$thA[ANM]; WARM$B <- S$thB[UN]   # warm-start cache
g_of <- function(phi, G, form, nm, key, reltol=1e-10, maxit=1500L) {
  t0 <- tic()
  obs <- moments_all(phi, G, "2")
  ob <- function(v){ p <- setNames(v,nm); pr <- moments_all(p, G, form)
    sum(vapply(seq_len(G$K), function(k) nll1(obs[[k]], pr[[k]], G$wk[k]*G$n), 0)) }
  o <- optim(WARM[[key]], ob, method="Nelder-Mead",
             control=list(maxit=maxit, reltol=reltol))
  p <- setNames(o$par, nm); po <- intersect(nm,c("oCL","oVc","sig")); p[po] <- abs(p[po])
  WARM[[key]] <- o$par                                   # warm start next call
  PROF$inner <- PROF$inner + (tic()-t0); PROF$ninner <- PROF$ninner + 1L
  p
}

say("timing one g_of call (cold) ...")
t1 <- tic(); ga <- g_of(S$thB, GRA, "1", ANM, "A", 1e-11, 4000L)
say(sprintf("  g_A cold = %.1fs   (%d moments_all calls)", tic()-t1, PROF$nmom))
t1 <- tic(); ga2 <- g_of(S$thB, GRA, "1", ANM, "A", 1e-11, 4000L)
say(sprintf("  g_A warm = %.1fs   speedup %.1fx", tic()-t1, 1))
say(sprintf("  g_A(theta_B) induced : %s", paste(sprintf("%s %.3f", ANM, ga), collapse="  ")))
say(sprintf("  A PUBLISHED          : %s", paste(sprintf("%s %.3f", ANM, S$thA[ANM]), collapse="  ")))
gb <- g_of(S$thB, GRB, "2", UN, "B", 1e-11, 4000L)
say(sprintf("  identity check |g_B(theta_B)-theta_B| = %.2e", max(abs(gb-S$thB[UN]))))

res <- function(phi, tol=1e-9, mx=800L)
  c((S$thA[ANM]-g_of(phi,GRA,"1",ANM,"A",tol,mx))/S$seA[ANM],
    (S$thB[UN] -g_of(phi,GRB,"2",UN ,"B",tol,mx))/S$seB[UN])
FD <- setNames(rep(6e-3,length(UN)), UN); FD[c("bCL1","bCL2","bVc1")] <- 1.2e-2
phi <- S$thB[UN]; bs <- Inf; bp <- phi
for (it in 1:6) {
  ti <- tic(); r0 <- res(phi, 1e-11, 3000L)
  if (sum(r0^2) < bs) { bs <- sum(r0^2); bp <- phi }
  J <- vapply(UN, function(q){ ph<-phi; ph[[q]]<-ph[[q]]+FD[[q]]
        (res(ph)-r0)/FD[[q]] }, numeric(length(r0)))
  phi <- phi - 0.6*as.numeric(solve(crossprod(J)+diag(1e-7,length(UN)), crossprod(J,r0)))
  say(sprintf("[%4.0fs] GN %d SSR %9.3f  (iter %.0fs)", tic()-T0, it, sum(r0^2), tic()-ti))
}
PLL <- bp
say(sprintf("\nPROFILE  moments_all: %d calls, %.1fs total (%.2f ms each)",
            PROF$nmom, PROF$mom, 1000*PROF$mom/PROF$nmom))
say(sprintf("PROFILE  inner fits : %d calls, %.1fs total (%.2f s each) = %.0f%% of wall",
            PROF$ninner, PROF$inner, PROF$inner/PROF$ninner, 100*PROF$inner/(tic()-T0)))
say(sprintf("\n%-6s %8s %10s %10s %10s %10s", "par","TRUTH","IPD-POOL","PARAM-LL","IPD-A","IPD-B"))
for (q in UN) say(sprintf("%-6s %8.4f %10.4f %10.4f %10s %10.4f", q, TP[[q]], S$GOLD[[q]],
  PLL[[q]], if (q %in% ANM) sprintf("%.4f", S$thA[[q]]) else "-", S$thB[[q]]))
say(sprintf("\nRMSE vs TRUTH     IPD-POOL %.4f | PARAM-LL %.4f | IPD-B alone %.4f",
  sqrt(mean((S$GOLD-TP)^2)), sqrt(mean((PLL-TP)^2)), sqrt(mean((S$thB-TP)^2))))
say(sprintf("RMSE vs IPD-POOL  PARAM-LL %.4f | IPD-B alone %.4f  <- the real comparison",
  sqrt(mean((PLL-S$GOLD)^2)), sqrt(mean((S$thB-S$GOLD)^2))))
say(sprintf("\nbCL2 (renal, ABSENT from A's model): truth %.3f gold %.3f PARAM-LL %.3f B-alone %.3f",
  TP[["bCL2"]], S$GOLD[["bCL2"]], PLL[["bCL2"]], S$thB[["bCL2"]]))
saveRDS(list(TP=TP,GOLD=S$GOLD,PLL=PLL,thA=S$thA,thB=S$thB), file.path(OUT,"gma-ipd-final.rds"))
say(sprintf("total %.0f s", tic()-T0))
