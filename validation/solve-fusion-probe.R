# ---------------------------------------------------------------------------
# validation/solve-fusion-probe.R
#
# Two "could we halve the solve count?" hypotheses, MEASURED rather than argued.
#
# (1) admc pays TWO solves per study per optimiser iteration: .admNLL solves the
#     simulation model, .admGrad solves the sensitivity model, at the same
#     parameters. Fusing them (scoring the NLL off the sens solve's rx_pred_)
#     would halve admc's solve count -- IF the two models' predictions are
#     bit-identical. The sens model carries extra variational compartments, so
#     lsoda's adaptive stepping sees a different system.
#
# (2) Every estimator loops over studies, one rxSolve each. rxSolve accepts a
#     multi-ID event table, so N studies could in principle become ONE call.
#     CLAUDE.md records a REVERTED attempt that merged same-`ev` studies over the
#     UNION of their times -- that changed the output grid. A per-ID event table
#     does not. Does it change the predictions anyway?
#
# CLAUDE.md's warning that "a single-subject probe misleadingly shows
# 0.000e+00" applies to both, so both use many subjects spread over etas.
# ---------------------------------------------------------------------------

suppressMessages({
  library(rxode2)
  pkgload::load_all(".", quiet = TRUE, export_all = FALSE)
})

mod_fn <- function() {
  ini({
    tcl  <- log(4); tv <- log(30); tq <- log(2); tvp <- log(40)
    add.err <- 0.05; prop.err <- 0.1
    eta.cl ~ 0.09
    eta.v  ~ 0.04
  })
  model({
    cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
    q  <- exp(tq); vp <- exp(tvp)
    d/dt(central) <-  q/vp*periph - q/v*central - cl/v*central
    d/dt(periph)  <- -q/vp*periph + q/v*central
    cp <- central / v
    cp ~ add(add.err) + prop(prop.err)
  })
}

ui    <- suppressMessages(rxode2::rxode2(mod_fn))
pinfo <- admixr2:::.admParseIniDf(ui$iniDf, ui)
sens  <- admixr2:::.admLoadSensModel(ui)
rxMod <- admixr2:::.admLoadModel(ui)
rxode2::rxLoad(rxMod)

times <- c(0.25, 0.5, 1, 2, 4, 8, 12, 24)
mkstudy <- function(times, dose) {
  s <- admixr2:::.admNormaliseStudy(
    list(E = rep(1, length(times)), V = diag(length(times)), n = 200L,
         times = times, ev = rxode2::et(amt = dose)), "s")
  s$ev_full <- s$ev |> rxode2::et(s$times)
  s
}
s1 <- mkstudy(times, 100)
s2 <- mkstudy(c(0.5, 1, 3, 6, 12, 24), 300)

pars   <- admixr2:::.admUnpack(admixr2:::.admBuildOptVec(pinfo)$p0, pinfo)
n_sim  <- 500L
z      <- admixr2:::.admMakeZ(n_sim, pinfo, 1L, "sobol")[[1L]]
pm     <- admixr2:::.admMakeParamsList(n_sim, pinfo, 1L)[[1L]]
eta    <- z %*% t(pars$L)
colnames(eta) <- pinfo$eta_col_names

cat("\n(1) admc: sens-model rx_pred_ vs simulation-model prediction\n")
cat("    500 subjects spread over 2 etas, 8 observation times\n")
sim  <- admixr2:::.admSimulate(rxMod, pars$struct, pinfo$sigma_names, eta, s1,
                               "cp", pm, 1L, .Machine$integer.max, pinfo$sigdig)
sns  <- admixr2:::.admSimulateSens(sens, pars$struct, pinfo$sigma_names, eta, s1,
                                   1L, .Machine$integer.max, pars$sigma_var,
                                   pinfo$sigdig)
a <- as.numeric(sim); b <- as.numeric(sns$cp_mat)
cat(sprintf("    identical      : %s\n", identical(sim, sns$cp_mat)))
cat(sprintf("    max rel diff   : %.3e\n", max(abs(a - b) / pmax(abs(a), 1e-300))))
cat(sprintf("    n rows differing: %d of %d\n",
            sum(rowSums(abs(sim - sns$cp_mat)) > 0), nrow(sim)))
# what the fusion would cost the OBJECTIVE
arr <- admixr2:::.admResidRows(pinfo, "cp", pars$sigma_var, ncol(sim))
nll <- function(M) {
  mu <- colMeans(M); cc <- sweep(M, 2L, mu)
  V  <- crossprod(cc) / nrow(M)
  ap <- admixr2:::.admResidApply(mu, diag(V), arr, s1$times, V)
  Vp <- V * tcrossprod(ap$ms); diag(Vp) <- ap$dv
  admixr2:::nll_cov_cpp(s1$E, s1$V, ap$mu, Vp, s1$n)
}
cat(sprintf("    -2LL from sim  : %.17g\n", nll(sim)))
cat(sprintf("    -2LL from sens : %.17g\n", nll(sns$cp_mat)))
cat(sprintf("    -2LL delta     : %.3e\n", abs(nll(sim) - nll(sns$cp_mat))))

cat("\n(2) two studies in ONE multi-ID solve vs two separate solves\n")
cat("    (different doses AND different observation times)\n")
one <- function(s, id0) {
  ev <- as.data.frame(s$ev_full)
  ev$id <- id0
  ev
}
sep1 <- admixr2:::.admSimulate(rxMod, pars$struct, pinfo$sigma_names, eta, s1,
                               "cp", pm, 1L, .Machine$integer.max, pinfo$sigdig)
sep2 <- admixr2:::.admSimulate(rxMod, pars$struct, pinfo$sigma_names, eta, s2,
                               "cp", pm, 1L, .Machine$integer.max, pinfo$sigdig)

# merged: 2*n_sim subjects, ids 1..n_sim on study 1's events, the rest on study 2's
evm <- rbind(do.call(rbind, lapply(seq_len(n_sim), function(i) one(s1, i))),
             do.call(rbind, lapply(seq_len(n_sim), function(i) one(s2, n_sim + i))))
pmm <- rbind(pm, pm)
for (nm in names(pars$struct)) pmm[, nm] <- pars$struct[[nm]]
pmm[, pinfo$eta_col_names] <- rbind(eta, eta)
for (nm in pinfo$sigma_names) pmm[, nm] <- 0
outm <- suppressWarnings(rxode2::rxSolve(rxMod, params = as.data.frame(pmm),
                                         events = evm, cores = 1L,
                                         nDisplayProgress = .Machine$integer.max,
                                         sigdig = pinfo$sigdig))
outm <- as.data.frame(outm)
g1 <- outm[outm$id %in% seq_len(n_sim) & outm$time %in% s1$times, ]
g2 <- outm[!(outm$id %in% seq_len(n_sim)) & outm$time %in% s2$times, ]
m1 <- matrix(g1$cp, nrow = n_sim, byrow = TRUE)
m2 <- matrix(g2$cp, nrow = n_sim, byrow = TRUE)
for (lbl in c("study 1", "study 2")) {
  A <- if (lbl == "study 1") sep1 else sep2
  B <- if (lbl == "study 1") m1   else m2
  cat(sprintf("    %s: identical=%-5s  max rel diff = %.3e\n", lbl,
              identical(A, B),
              max(abs(as.numeric(A) - as.numeric(B)) /
                    pmax(abs(as.numeric(A)), 1e-300))))
}
