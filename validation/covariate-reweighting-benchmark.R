## Covariate reweighting in the package: correctness first, then speed.
##
## Correctness is checked against exact nested Gauss-Hermite on a closed-form
## model, so the reference carries no approximation of its own. Speed is measured
## as solves per objective evaluation and as wall time per evaluation, which is
## the number that survives changes in machine and model.
##
## Run:  Rscript validation/covariate-reweighting-benchmark.R
suppressMessages({library(rxode2)})
suppressMessages(devtools::load_all(".", quiet = TRUE))
kv <- function(k, v) cat(sprintf("%s\t%s\n", k, paste(v, collapse = ",")))

TT <- c(0.5, 1, 2, 4, 8); EV <- rxode2::et(amt = 100)
TCL <- log(1); TV <- log(10); TCOV <- 0.75; OM <- 0.30; ADD <- 0.30
WT <- list(meanlog = log(70), sdlog = 0.25)

## allometric: NOT a bare theta*COV product, so the closed-form collapse does not
## apply and this would otherwise take the per-subject path
mod <- function() {
  ini({tcl <- log(1); tv <- log(10); tcov <- 0.75; add.err <- 0.3
       eta.cl ~ 0.09})
  model({cl <- exp(tcl + eta.cl) * (WT / 70)^tcov
         v  <- exp(tv)
         d/dt(centr) <- -cl / v * centr; cp <- centr / v
         cp ~ add(add.err)})
}

## ---- exact reference: nested quadrature on the closed form -----------------
gh <- function(k) {
  i <- seq_len(k - 1L)
  J <- matrix(0, k, k); J[cbind(i, i + 1L)] <- sqrt(i); J[cbind(i + 1L, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE); o <- order(e$values)
  list(x = e$values[o], w = e$vectors[1L, o]^2)
}
QA <- gh(60L); QE <- gh(60L)
m1 <- 0; M2 <- 0
for (i in seq_along(QA$x)) {
  wt <- exp(WT$meanlog + WT$sdlog * QA$x[i])
  for (j in seq_along(QE$x)) {
    w  <- QA$w[i] * QE$w[j]
    cl <- exp(TCL + OM * QE$x[j]) * (wt / 70)^TCOV
    y  <- 100 / exp(TV) * exp(-cl / exp(TV) * TT)
    m1 <- m1 + w * y; M2 <- M2 + w * outer(y, y)
  }
}
REF <- list(E = m1, V = { V <- M2 - outer(m1, m1); diag(V) <- diag(V) + ADD^2; V })

ui   <- suppressMessages(rxode2::rxode2(mod))
ovar <- admixr2:::.admOutputVar(ui)
rxm  <- admixr2:::.admLoadModel(ui)
STUDY <- list(s = list(E = REF$E, V = REF$V, n = 100L, times = TT, ev = EV,
                       cov_dist = list(WT = WT)))

setup <- function(est, reweight, n_nodes = 7L, n_sim = 4000L) {
  ctl <- if (est == "adgh")
    adghControl(studies = STUDY, grad = "analytical", n_nodes = n_nodes,
                print = 0L, covMethod = "none", cov_reweight = reweight)
  else
    admControl(studies = STUDY, grad = "sens", n_sim = n_sim, print = 0L,
               covMethod = "none", cov_reweight = reweight)
  pinfo <- admixr2:::.admDriverPinfo(ui, ctl)
  u  <- admixr2:::.admCheckCovariates(
          ui, pinfo, admixr2:::.admDriverUnits(STUDY, ui, ovar)$studies,
          if (est == "adgh") "analytical" else "sens", est)
  ov <- admixr2:::.admBuildOptVec(pinfo)
  if (reweight)
    u <- admixr2:::.admCovPromoteReweight(
      ui, u, admixr2:::.admUnpack(ov$p0, pinfo), pinfo, rxm, 1L)
  list(ctl = ctl, pinfo = pinfo, u = u, ov = ov,
       grid = admixr2:::.adghNodeGrid(n_nodes, pinfo$n_eta), n_sim = n_sim)
}

## ---- 1. does the path change, and is it still right? -----------------------
cat(sprintf("%-22s %-10s %10s %10s %12s %10s\n",
            "estimator", "path", "E rel", "V rel", "s/eval", "solve rows"))
res <- list()
for (est in c("adgh", "admc")) for (rw in c(FALSE, TRUE)) {
  d <- setup(est, rw)
  s1 <- d$u[[1L]]
  if (est == "adgh") {
    pars <- admixr2:::.admUnpack(d$ov$p0, d$pinfo)
    g <- admixr2:::.adghGrid(pars, d$pinfo, d$grid, s1)
    mm <- admixr2:::.adghMoments(pars, d$pinfo, s1, rxm, ovar, d$grid, 1L)
    rows <- nrow(g$eta)
    f <- function() admixr2:::.adghNLL(d$ov$p0, d$pinfo, d$u, rxm, ovar, d$grid, 1L)
  } else {
    zl <- admixr2:::.admMakeZ(d$n_sim, d$pinfo, 1L, "sobol")
    pl <- admixr2:::.admMakeParamsList(d$n_sim, d$pinfo, 1L)
    rows <- d$n_sim
    mm <- NULL
    f <- function() admixr2:::.admNLL(d$ov$p0, d$pinfo, d$u, zl, rxm, ovar, pl, 1L)
  }
  for (i in 1:3) v <- f()          # warm-up, NOT timed
  t0 <- Sys.time(); for (i in 1:10) v <- f()
  sec <- as.numeric(difftime(Sys.time(), t0, units = "secs")) / 10
  erel <- vrel <- NA_real_
  if (!is.null(mm)) {
    erel <- max(abs(mm$E - REF$E) / abs(REF$E))
    vrel <- max(abs(mm$V - REF$V)) / max(abs(REF$V))
  }
  cat(sprintf("%-22s %-10s %10s %10s %12.4f %10d\n",
              paste0(est, if (rw) " reweight" else " baseline"),
              s1$.adm_cov_path,
              if (is.na(erel)) "-" else sprintf("%.2e", erel),
              if (is.na(vrel)) "-" else sprintf("%.2e", vrel),
              sec, rows))
  res[[paste(est, rw)]] <- list(sec = sec, rows = rows, nll = v)
}

cat("\n")
kv("adgh_speedup",  sprintf("%.2fx  (%d -> %d solve rows)",
   res[["adgh FALSE"]]$sec / res[["adgh TRUE"]]$sec,
   res[["adgh FALSE"]]$rows, res[["adgh TRUE"]]$rows))
kv("admc_speedup",  sprintf("%.2fx  (%d -> %d solve rows)",
   res[["admc FALSE"]]$sec / res[["admc TRUE"]]$sec,
   res[["admc FALSE"]]$rows, res[["admc TRUE"]]$rows))
kv("adgh_nll_agreement", sprintf("baseline %.6f  reweight %.6f  diff %.2e",
   res[["adgh FALSE"]]$nll, res[["adgh TRUE"]]$nll,
   abs(res[["adgh TRUE"]]$nll - res[["adgh FALSE"]]$nll)))
kv("admc_nll_agreement", sprintf("baseline %.6f  reweight %.6f  diff %.2e",
   res[["admc FALSE"]]$nll, res[["admc TRUE"]]$nll,
   abs(res[["admc TRUE"]]$nll - res[["admc FALSE"]]$nll)))

## ---- 2. admc: the gain is VARIANCE, not solves -----------------------------
## Same rows either way, but the covariate integral is exact in the weight, so
## only u is sampled. Compare accuracy at matched n_sim.
cat("\nadmc accuracy at matched n_sim (the covariate dimension leaves the sample):\n")
cat(sprintf("%-8s %14s %14s\n", "n_sim", "baseline |dNLL|", "reweight |dNLL|"))
d_ex <- setup("adgh", FALSE, n_nodes = 25L)
NLL_REF <- admixr2:::.adghNLL(d_ex$ov$p0, d_ex$pinfo, d_ex$u, rxm, ovar,
                              d_ex$grid, 1L)
for (ns in c(250L, 1000L, 4000L, 16000L)) {
  row <- c(ns, NA, NA)
  for (k in seq_along(c(FALSE, TRUE))) {
    rw <- c(FALSE, TRUE)[k]
    d  <- setup("admc", rw, n_sim = ns)
    zl <- admixr2:::.admMakeZ(ns, d$pinfo, 1L, "sobol")
    pl <- admixr2:::.admMakeParamsList(ns, d$pinfo, 1L)
    v  <- admixr2:::.admNLL(d$ov$p0, d$pinfo, d$u, zl, rxm, ovar, pl, 1L)
    row[k + 1L] <- abs(v - NLL_REF)
  }
  cat(sprintf("%-8d %14.4f %14.4f\n", ns, row[2], row[3]))
}

## ---- 3. a bigger grid, where the solve count actually dominates -------------
## rxSolve costs ~11 ms before it integrates anything, so a 7-row solve is not
## 11x cheaper than a 77-row one -- the fixed call cost swamps it. The saving
## only shows when the grid is large enough to matter: with TWO random effects
## the product grid is n_node^2 * cov_nodes.
mod2 <- function() {
  ini({tcl <- log(1); tv <- log(10); tcov <- 0.75; add.err <- 0.3
       eta.cl ~ 0.09; eta.v ~ 0.04})
  model({cl <- exp(tcl + eta.cl) * (WT / 70)^tcov
         v  <- exp(tv + eta.v)
         d/dt(centr) <- -cl / v * centr; cp <- centr / v
         cp ~ add(add.err)})
}
ui2   <- suppressMessages(rxode2::rxode2(mod2))
ovar2 <- admixr2:::.admOutputVar(ui2)
rxm2  <- admixr2:::.admLoadModel(ui2)
cat("\nTwo random effects (grid = n_node^2 x cov_nodes):\n")
cat(sprintf("%-12s %-10s %10s %12s %8s\n",
            "n_nodes", "path", "rows", "s/eval", "speedup"))
for (k in c(5L, 7L, 9L)) {
  tm <- c(NA_real_, NA_real_); rw_rows <- c(NA_integer_, NA_integer_)
  for (q in 1:2) {
    rw  <- c(FALSE, TRUE)[q]
    ctl <- adghControl(studies = STUDY, grad = "analytical", n_nodes = k,
                       print = 0L, covMethod = "none", cov_reweight = rw)
    pin <- admixr2:::.admDriverPinfo(ui2, ctl)
    uu  <- admixr2:::.admCheckCovariates(
             ui2, pin, admixr2:::.admDriverUnits(STUDY, ui2, ovar2)$studies,
             "analytical", "adgh")
    ov2 <- admixr2:::.admBuildOptVec(pin)
    if (rw) uu <- admixr2:::.admCovPromoteReweight(
      ui2, uu, admixr2:::.admUnpack(ov2$p0, pin), pin, rxm2, 1L)
    gr <- admixr2:::.adghNodeGrid(k, pin$n_eta)
    pr <- admixr2:::.admUnpack(ov2$p0, pin)
    rw_rows[q] <- nrow(admixr2:::.adghGrid(pr, pin, gr, uu[[1L]])$eta)
    f <- function() admixr2:::.adghNLL(ov2$p0, pin, uu, rxm2, ovar2, gr, 1L)
    for (i in 1:3) f()
    t0 <- Sys.time(); for (i in 1:8) f()
    tm[q] <- as.numeric(difftime(Sys.time(), t0, units = "secs")) / 8
    cat(sprintf("%-12s %-10s %10d %12.4f %8s\n",
                if (q == 1) sprintf("%d", k) else "",
                if (rw) "reweight" else "rows", rw_rows[q], tm[q],
                if (q == 2) sprintf("%.2fx", tm[1] / tm[2]) else ""))
  }
}
