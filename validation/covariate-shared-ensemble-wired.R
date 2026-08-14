## Does the wired shared ensemble agree with the per-study path, and is it faster?
suppressMessages({library(rxode2)})
suppressMessages(devtools::load_all("C:/package/admixr2/.claude/worktrees/feature-covariate-quadrature", quiet = TRUE))
kv <- function(k, v) cat(sprintf("%s\t%s\n", k, paste(v, collapse = ",")))
TT <- c(0.5, 1, 2, 4, 8); EV <- rxode2::et(amt = 100)
POPS <- list(A = list(meanlog = log(50), sdlog = 0.22),
             B = list(meanlog = log(70), sdlog = 0.30),
             C = list(meanlog = log(95), sdlog = 0.18),
             D = list(meanlog = log(60), sdlog = 0.35),
             E = list(meanlog = log(80), sdlog = 0.25))
mod <- function() {
  ini({tcl <- log(1); tv <- log(10); tcov <- 0.75; add.err <- 0.3; eta.cl ~ 0.09})
  model({cl <- exp(tcl + eta.cl) * (WT / 70)^tcov; v <- exp(tv)
         d/dt(centr) <- -cl / v * centr; cp <- centr / v; cp ~ add(add.err)})
}
ui <- suppressMessages(rxode2::rxode2(mod)); ovar <- admixr2:::.admOutputVar(ui)
rxm <- admixr2:::.admLoadModel(ui)
## REALISTIC observed data: each study's own exact marginal moments. Dummy data
## (E = 1, V = I) makes the objective enormously sensitive to tiny moment
## differences and an earlier version of this script reported a 1e5 NLL gap that
## was entirely that, not a disagreement between the paths.
gh <- function(k){i<-seq_len(k-1L);J<-matrix(0,k,k);J[cbind(i,i+1L)]<-sqrt(i)
  J[cbind(i+1L,i)]<-sqrt(i);e<-eigen(J,symmetric=TRUE);o<-order(e$values)
  list(x=e$values[o],w=e$vectors[1L,o]^2)}
QA<-gh(50L);QE<-gh(50L)
ref_of <- function(pp){m1<-0;M2<-0
  for(i in seq_along(QA$x)){wt<-exp(pp$meanlog+pp$sdlog*QA$x[i])
    for(j in seq_along(QE$x)){w<-QA$w[i]*QE$w[j]
      cl<-exp(0+0.30*QE$x[j])*(wt/70)^0.75
      y<-100/10*exp(-cl/10*TT); m1<-m1+w*y; M2<-M2+w*outer(y,y)}}
  V<-M2-outer(m1,m1); diag(V)<-diag(V)+0.09; list(E=m1,V=V)}
REFM <- lapply(POPS, ref_of)
mk <- function(k) setNames(lapply(seq_len(k), function(i)
  list(E = REFM[[i]]$E, V = REFM[[i]]$V, n = 100L,
       times = TT, ev = EV, cov_dist = list(WT = POPS[[i]]))), names(POPS)[seq_len(k)])

run <- function(S, rw, n_sim = 4000L) {
  st  <- mk(S)
  ctl <- admControl(studies = st, grad = "sens", n_sim = n_sim, print = 0L,
                    covMethod = "none", cov_reweight = rw)
  pin <- admixr2:::.admDriverPinfo(ui, ctl)
  u   <- admixr2:::.admCheckCovariates(
           ui, pin, admixr2:::.admDriverUnits(st, ui, ovar)$studies, "sens")
  ov  <- admixr2:::.admBuildOptVec(pin)
  if (rw) u <- admixr2:::.admCovPromoteReweight(
    ui, u, admixr2:::.admUnpack(ov$p0, pin), pin, rxm, 1L)
  zl <- admixr2:::.admMakeZ(n_sim, pin, length(u), "sobol")
  pl <- admixr2:::.admMakeParamsList(n_sim, pin, length(u))
  f <- function() admixr2:::.admNLL(ov$p0, pin, u, zl, rxm, ovar, pl, 1L)
  for (i in 1:2) f()
  t0 <- Sys.time(); for (i in 1:5) v <- f()
  list(nll = v, sec = as.numeric(difftime(Sys.time(), t0, units = "secs")) / 5,
       paths = paste(unique(vapply(u, function(s) s$.adm_cov_path %||% "-",
                                   character(1))), collapse = ","))
}
cat(sprintf("%-3s %-10s %14s %10s %10s\n", "S", "path", "NLL", "s/eval", "vs base"))
for (S in c(1L, 3L, 5L)) {
  a <- run(S, FALSE); b <- run(S, TRUE)
  cat(sprintf("%-3d %-10s %14.6f %10.4f %10s\n", S, a$paths, a$nll, a$sec, ""))
  cat(sprintf("%-3s %-10s %14.6f %10.4f %9.2fx   dNLL %.2e\n", "", b$paths,
              b$nll, b$sec, a$sec / b$sec, abs(b$nll - a$nll)))
}
