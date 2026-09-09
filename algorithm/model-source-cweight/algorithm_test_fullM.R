## The FULL fix: reuse the REAL .admSrcJac (E and V jointly, via datagen() --
## so stratification/conditioning are inherited automatically, not hand-rolled)
## instead of the mean-only rxSolve shortcut. Then empirically test whether M
## needs to be "live" -- i.e. whether D, built by perturbing the SOURCE's own
## model at its own reported theta_src, would meaningfully differ if evaluated
## somewhere else, which is the only way "live" could matter here: D is a
## property of the source model at theta_src, not of the analysis model's
## candidate theta at all (that's what's different from the SE's own M, which
## needs G to project into analysis-parameter space; the objective's residual
## is already in moment space and needs no such projection).
suppressMessages({ library(nlmixr2est); devtools::load_all(".", quiet = TRUE) })
rxode2::setRxThreads(2L)

TRUE_TCL <- log(5); TRUE_TV1 <- log(30); TRUE_Q <- log(8); TRUE_V2 <- log(60)
TRUE_ADD <- 0.10; TRUE_OM <- 0.06
DOSE <- 200; TAU <- 12

SRC_TIMES <- c(0.5, 1, 2, 4, 8, 11.5)
mk_src_ui <- function(tcl, tv, add.err, om) {
  fn <- sprintf(paste0(
    "function(){\n",
    "  ini({ tcl <- %.6f; tv <- %.6f; add.err <- %.6f\n",
    "        eta.cl ~ %.6f })\n",
    "  model({ cl <- exp(tcl + eta.cl); v <- exp(tv)\n",
    "          cp <- linCmt(); cp ~ add(add.err) }) }"),
    tcl, tv, add.err, om)
  suppressMessages(rxode2::rxode2(eval(parse(text = fn))))
}
SD_SRC <- c(tcl = 0.070, tv = 0.050, add.err = 0.004, eta.cl = 0.010)
CSRC_1C <- diag(SD_SRC^2); dimnames(CSRC_1C) <- list(names(SD_SRC), names(SD_SRC))

gen_src_1c <- function(theta, n_model, seed) {
  ui <- mk_src_ui(theta[["tcl"]], theta[["tv"]], theta[["add.err"]], theta[["eta.cl"]])
  ev <- rxode2::et(amt = DOSE, ii = TAU, ss = 1, addl = 0)
  sp <- list(times = SRC_TIMES, ev = ev, n = n_model, model_cov = CSRC_1C)
  suppressWarnings(suppressMessages(datagen(
    list(src = sp), model = ui, control = datagenControl(method = "gh", seed = seed))))
}

DENSE_TIMES <- c(0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 12, 24)
mk_2c_ui <- function(tcl, tv, q, v2, add.err, om) {
  fn <- sprintf(paste0(
    "function(){\n",
    "  ini({ tcl <- %.6f; tv <- %.6f; q <- %.6f; v2 <- %.6f; add.err <- %.6f\n",
    "        eta.cl ~ %.6f })\n",
    "  model({ cl <- exp(tcl + eta.cl); v1 <- exp(tv); qq <- exp(q); vv2 <- exp(v2)\n",
    "          d/dt(centr) <- -cl/v1*centr - qq/v1*centr + qq/vv2*periph\n",
    "          d/dt(periph) <- qq/v1*centr - qq/vv2*periph\n",
    "          cp <- centr/v1; cp ~ add(add.err) }) }"),
    tcl, tv, q, v2, add.err, om)
  suppressMessages(rxode2::rxode2(eval(parse(text = fn))))
}
.fit_2c <- function() {
  ini({ tcl <- log(4); tv <- log(25); q <- log(5); v2 <- log(45)
        add.err <- 0.15; eta.cl ~ 0.1 })
  model({ cl <- exp(tcl + eta.cl); v1 <- exp(tv); qq <- exp(q); vv2 <- exp(v2)
          d/dt(centr) <- -cl/v1*centr - qq/v1*centr + qq/vv2*periph
          d/dt(periph) <- qq/v1*centr - qq/vv2*periph
          cp <- centr/v1; cp ~ add(add.err) })
}
gen_data_2c <- function(n_data) {
  ui <- mk_2c_ui(TRUE_TCL, TRUE_TV1, TRUE_Q, TRUE_V2, TRUE_ADD, TRUE_OM^2)
  ev <- rxode2::et(amt = DOSE)
  sp <- list(times = DENSE_TIMES, ev = ev, n = n_data)
  suppressWarnings(suppressMessages(datagen(
    list(dat = sp), model = ui, control = datagenControl(method = "gh", seed = 2L))))
}

## ---- build M via the REAL .admSrcJac (E + V jointly) -----------------------
pinv_reg <- function(M, tol = 1e-8) {
  ev <- eigen((M + t(M)) / 2, symmetric = TRUE)
  keep <- ev$values > max(ev$values) * tol
  list(Mi = ev$vectors[, keep, drop = FALSE] %*%
         diag(1 / ev$values[keep], sum(keep)) %*%
         t(ev$vectors[, keep, drop = FALSE]),
       rank = sum(keep), values = ev$values)
}

build_M <- function(studies) {
  grp <- admixr2:::.admSrcGroups(studies)
  idx <- grp[[1L]]
  D   <- admixr2:::.admSrcJac(studies, idx)
  prov <- studies[[idx[1L]]][[".adm_src"]]
  Dm  <- D[[1L]]                       # single (unstratified) source study
  M   <- Dm %*% prov$cov[colnames(Dm), colnames(Dm), drop = FALSE] %*% t(Dm)
  M
}

## ---- patched objective, full tau (E+V), M built once from .admSrcJac -------
make_nll_fullM <- function(Mpinv) {
  force(Mpinv)
  function(p, pinfo, studies, rxMod, out_var, grid, cores) {
    pars <- tryCatch(admixr2:::.admUnpack(p, pinfo), error = function(e) NULL)
    if (is.null(pars)) return(Inf)
    if (!admixr2:::.admParsFinite(pars, pinfo)) return(Inf)
    total <- 0
    grp <- admixr2:::.admSrcGroups(studies)
    skip_idx <- if (length(grp)) unlist(grp) else integer(0)
    for (i in seq_along(studies)) {
      if (i %in% skip_idx) next
      s <- studies[[i]]
      m <- admixr2:::.adghMoments(pars, pinfo, s, rxMod, if (is.null(s$output)) out_var else s$output, grid, cores)
      nll <- nll_cov_cpp(s$E, s$V, m$E, m$V, s$n)
      if (!is.finite(nll)) return(Inf)
      total <- total + nll
    }
    for (nm in names(grp)) {
      idx <- grp[[nm]]
      s   <- studies[[idx[1L]]]
      m   <- admixr2:::.adghMoments(pars, pinfo, s, rxMod, if (is.null(s$output)) out_var else s$output, grid, cores)
      tau_pred <- admixr2:::.admTauVec(m$E, m$V, s)
      tau_src  <- admixr2:::.admTauVec(s$E, s$V, s)
      r <- tau_pred - tau_src
      total <- total + as.numeric(t(r) %*% Mpinv %*% r)
    }
    total
  }
}
.orig_adghNLL <- get(".adghNLL", envir = asNamespace("admixr2"))
fit_ctl <- function(g, fn_new) {
  assignInNamespace(".adghNLL", fn_new, ns = "admixr2")
  on.exit(assignInNamespace(".adghNLL", .orig_adghNLL, ns = "admixr2"))
  suppressMessages(nlmixr2est::nlmixr2(.fit_2c, admData(), est = "adgh",
    control = adghControl(studies = g, print = 0L, cores = 2L, covMethod = "none", grad = "fd")))
}

cat("==================== full E+V fix, via real .admSrcJac =================\n")
g_data <- gen_data_2c(300L)
theta_off <- c(tcl = TRUE_TCL - 0.15, tv = TRUE_TV1 - 0.10, add.err = TRUE_ADD, eta.cl = TRUE_OM^2)
g_src <- gen_src_1c(theta_off, 300L, seed = 55L)
studies0 <- c(g_data, g_src)

M0 <- build_M(studies0)
pv0 <- pinv_reg(M0)
cat("M eigenvalues (built at theta_src):", format(pv0$values, digits = 3), "\n")
cat("rank kept:", pv0$rank, "of", length(pv0$values), "\n")

nll_fullM <- make_nll_fullM(pv0$Mi)
f <- fit_ctl(studies0, nll_fullM)
p <- f$parFixedDf
cat(sprintf("tcl=%.4f  tv=%.4f  q=%.4f  v2=%.4f\n",
            p["tcl","Estimate"], p["tv","Estimate"], p["q","Estimate"], p["v2","Estimate"]))
cat(sprintf("truth:  tcl=%.4f tv=%.4f q=%.4f v2=%.4f\n", TRUE_TCL, TRUE_TV1, TRUE_Q, TRUE_V2))

## ---- "live M": does D actually change if evaluated somewhere other than theta_src? ----
cat("\n==================== does D depend on the reference point? ============\n")
## Build D at theta_src (as above), vs at a data-informed alternative -- a
## point some distance away, to see whether the LOCAL sensitivity D actually
## differs enough to matter.
theta_alt <- c(tcl = TRUE_TCL, tv = TRUE_TV1, add.err = TRUE_ADD, eta.cl = TRUE_OM^2)
g_src_alt <- gen_src_1c(theta_alt, 300L, seed = 999L)   # only used to rebuild D at a different theta
studies_alt <- c(g_data, g_src_alt)
M_alt <- build_M(studies_alt)
pv_alt <- pinv_reg(M_alt)
cat("M eigenvalues (built at truth instead of theta_src):", format(pv_alt$values, digits = 3), "\n")
## compare the two M's directly (on their shared rank-2 subspace) via principal angles
cat(sprintf("Frobenius norm of (M0 - M_alt) relative to |M0|: %.4f\n",
            norm(M0 - M_alt, "F") / norm(M0, "F")))
