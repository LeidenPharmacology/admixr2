## Does using a diagonal-only C_src (all a paper's SE column gives you) instead
## of the FULL reported covariance change the POINT ESTIMATE under the NEW
## method -- not just the reported SE, which is the old (already-established)
## concern. Same theta_off (what the source reported) in both arms; only the
## covariance passed as model_cov differs, so any difference in the fit is
## attributable ENTIRELY to the ignored correlation.
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
RHO <- -0.6   # a realistic, moderate-strong CL/V correlation from a real fit
R   <- diag(4); R[1,2] <- R[2,1] <- RHO
dimnames(R) <- list(names(SD_SRC), names(SD_SRC))
CSRC_FULL <- diag(SD_SRC) %*% R %*% diag(SD_SRC); dimnames(CSRC_FULL) <- dimnames(R)
CSRC_DIAG <- diag(SD_SRC^2); dimnames(CSRC_DIAG) <- dimnames(R)
cat("CSRC_FULL:\n"); print(round(CSRC_FULL, 5))
cat("corr(tcl,tv) in FULL:", cov2cor(CSRC_FULL)[1,2], "  in DIAG:", 0, "\n\n")

gen_src_1c <- function(theta, n_model, seed, csrc) {
  ui <- mk_src_ui(theta[["tcl"]], theta[["tv"]], theta[["add.err"]], theta[["eta.cl"]])
  ev <- rxode2::et(amt = DOSE, ii = TAU, ss = 1, addl = 0)
  sp <- list(times = SRC_TIMES, ev = ev, n = n_model, model_cov = csrc)
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

pinv_reg <- function(M, tol = 1e-8) {
  ev <- eigen((M + t(M)) / 2, symmetric = TRUE)
  keep <- ev$values > max(ev$values) * tol
  ev$vectors[, keep, drop = FALSE] %*% diag(1 / ev$values[keep], sum(keep)) %*%
    t(ev$vectors[, keep, drop = FALSE])
}
build_M <- function(studies) {
  grp <- admixr2:::.admSrcGroups(studies)
  idx <- grp[[1L]]
  D   <- admixr2:::.admSrcJac(studies, idx)[[1L]]
  prov <- studies[[idx[1L]]][[".adm_src"]]
  D %*% prov$cov[colnames(D), colnames(D), drop = FALSE] %*% t(D)
}
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

g_data <- gen_data_2c(300L)
theta_off <- c(tcl = TRUE_TCL - 0.15, tv = TRUE_TV1 - 0.10, add.err = TRUE_ADD, eta.cl = TRUE_OM^2)

cat("==================== FULL C_src (real correlation) =====================\n")
g_src_full <- gen_src_1c(theta_off, 300L, seed = 55L, csrc = CSRC_FULL)
studies_full <- c(g_data, g_src_full)
Mf <- pinv_reg(build_M(studies_full))
ff <- fit_ctl(studies_full, make_nll_fullM(Mf))
pf <- ff$parFixedDf
cat(sprintf("tcl=%.4f  tv=%.4f  q=%.4f  v2=%.4f\n",
            pf["tcl","Estimate"], pf["tv","Estimate"], pf["q","Estimate"], pf["v2","Estimate"]))

cat("\n==================== DIAGONAL C_src (correlation dropped) ==============\n")
g_src_diag <- gen_src_1c(theta_off, 300L, seed = 55L, csrc = CSRC_DIAG)
studies_diag <- c(g_data, g_src_diag)
Md <- pinv_reg(build_M(studies_diag))
fd <- fit_ctl(studies_diag, make_nll_fullM(Md))
pd <- fd$parFixedDf
cat(sprintf("tcl=%.4f  tv=%.4f  q=%.4f  v2=%.4f\n",
            pd["tcl","Estimate"], pd["tv","Estimate"], pd["q","Estimate"], pd["v2","Estimate"]))

cat(sprintf("\ntruth:  tcl=%.4f tv=%.4f q=%.4f v2=%.4f\n", TRUE_TCL, TRUE_TV1, TRUE_Q, TRUE_V2))
cat(sprintf("source reported: tcl=%.4f tv=%.4f\n", theta_off[["tcl"]], theta_off[["tv"]]))
cat(sprintf("\ndifference (full - diag): tcl=%.4f tv=%.4f q=%.4f v2=%.4f\n",
            pf["tcl","Estimate"]-pd["tcl","Estimate"], pf["tv","Estimate"]-pd["tv","Estimate"],
            pf["q","Estimate"]-pd["q","Estimate"], pf["v2","Estimate"]-pd["v2","Estimate"]))
