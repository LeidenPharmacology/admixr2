## Test: replace a model source's n-scaled Wishart contribution to the OBJECTIVE
## with a direct quadratic penalty on shared named parameters, weighted by
## C_src -- no n anywhere. S = identity special case (source and analysis model
## share the same ini() names), matching every fixture used this session.
##
## Method: monkey-patch .adghNLL for the duration of the "NEW" fit only, force
## grad = "fd" so the (unpatched) analytic gradient never runs against a
## mismatched objective, leave covMethod = "r,s" untouched (the SE machinery is
## validated separately -- this tests the WEIGHT only).
suppressMessages({ library(nlmixr2est); devtools::load_all(".", quiet = TRUE) })
rxode2::setRxThreads(2L)

TIMES <- c(0.5, 1, 2, 4, 8, 12, 24); DOSE <- 200
TRUE_TCL <- log(5); TRUE_TV <- log(50); TRUE_BWT <- 0.75
TRUE_ADD <- 0.08; TRUE_OM <- 0.05
WTDIST <- covDist(WT = c(mean = 78, sd = 16), dist = "lnorm")
SD  <- c(tcl = 0.080, tv = 0.060, bwt = 0.055, add.err = 0.004, eta.cl = 0.010)
CSRC <- diag(SD^2); dimnames(CSRC) <- list(names(SD), names(SD))
J_STRATA <- 5L

mk_ui <- function(tcl, tv, bwt, add.err, om) {
  fn <- sprintf(paste0(
    "function(){\n",
    "  ini({ tcl <- %.6f; tv <- %.6f; bwt <- %.6f; add.err <- %.6f\n",
    "        eta.cl ~ %.6f })\n",
    "  model({ cl <- exp(tcl + eta.cl) * (WT/70)^bwt; v <- exp(tv)\n",
    "          cp <- linCmt(); cp ~ add(add.err) }) }"),
    tcl, tv, bwt, add.err, om)
  suppressMessages(rxode2::rxode2(eval(parse(text = fn))))
}
.ms_fit_wt <- function() {
  ini({ tcl <- log(4); tv <- log(45); bwt <- 0.4
        eta.cl ~ 0.1; add.err <- 0.1 })
  model({ cl <- exp(tcl + eta.cl) * (WT/70)^bwt; v <- exp(tv)
          cp <- linCmt(); cp ~ add(add.err) })
}
gen_data <- function(n_data) {
  ui <- mk_ui(TRUE_TCL, TRUE_TV, TRUE_BWT, TRUE_ADD, TRUE_OM^2)
  sp <- list(times = TIMES, ev = rxode2::et(amt = DOSE), n = n_data, cov_dist = WTDIST,
             stratify = "WT", strata_nodes = J_STRATA, cov_range = list(WT = c(50, 115)))
  suppressWarnings(suppressMessages(datagen(
    list(dat = sp), model = ui, control = datagenControl(method = "gh", seed = 1L))))
}
gen_src <- function(theta, n_model, seed) {
  ui <- mk_ui(theta[["tcl"]], theta[["tv"]], theta[["bwt"]], theta[["add.err"]], theta[["eta.cl"]])
  sp <- list(times = TIMES, ev = rxode2::et(amt = DOSE), n = n_model,
             model_cov = CSRC, cov_dist = WTDIST,
             stratify = "WT", strata_nodes = J_STRATA, cov_range = list(WT = c(50, 115)))
  suppressWarnings(suppressMessages(datagen(
    list(src = sp), model = ui, control = datagenControl(method = "gh", seed = seed))))
}

## ---- the patched objective ---------------------------------------------------
.adghNLL_srcfix <- function(p, pinfo, studies, rxMod, out_var, grid, cores) {
  pars <- tryCatch(admixr2:::.admUnpack(p, pinfo), error = function(e) NULL)
  if (is.null(pars)) return(Inf)
  if (!admixr2:::.admParsFinite(pars, pinfo)) return(Inf)
  total <- 0
  grp <- admixr2:::.admSrcGroups(studies)
  skip_idx <- if (length(grp)) unlist(grp) else integer(0)
  for (i in seq_along(studies)) {
    if (i %in% skip_idx) next
    s <- studies[[i]]
    if (isTRUE(s$is_joint)) {
      m   <- admixr2:::.adghMomentsJoint(pars, pinfo, s, rxMod, grid, cores)
      nll <- nll_cov_cpp(s$E, s$V, m$E, m$V, s$n)
    } else {
      m <- admixr2:::.adghMoments(pars, pinfo, s, rxMod, if (is.null(s$output)) out_var else s$output, grid, cores)
      nll <- if (identical(s$method, "var"))
        nll_var_cpp(s$E, s$v_diag, m$E, diag(m$V), s$n)
      else
        nll_cov_cpp(s$E, s$V, m$E, m$V, s$n)
    }
    if (!is.finite(nll)) return(Inf)
    total <- total + nll
  }
  ## model-source groups: ONE quadratic penalty per source, on shared named
  ## parameters, weighted by C_src -- no n anywhere.
  th_now <- admixr2:::.admFullTheta(pars, pinfo)
  for (nm in names(grp)) {
    idx  <- grp[[nm]]
    prov <- studies[[idx[1L]]][[".adm_src"]]
    shared <- intersect(names(th_now), prov$par)
    if (!length(shared)) next
    delta <- th_now[shared] - unlist(prov$theta[shared])
    Ci <- tryCatch(solve(prov$cov[shared, shared, drop = FALSE]), error = function(e) NULL)
    if (is.null(Ci)) return(Inf)
    total <- total + as.numeric(t(delta) %*% Ci %*% delta)
  }
  total
}

.orig_adghNLL <- get(".adghNLL", envir = asNamespace("admixr2"))

fit_ctl <- function(g, use_fix, cm = "r,s") {
  if (use_fix) assignInNamespace(".adghNLL", .adghNLL_srcfix, ns = "admixr2")
  on.exit(if (use_fix) assignInNamespace(".adghNLL", .orig_adghNLL, ns = "admixr2"))
  suppressMessages(nlmixr2est::nlmixr2(.ms_fit_wt, admData(), est = "adgh",
    control = adghControl(studies = g, print = 0L, cores = 2L, covMethod = cm,
                          grad = "fd")))
}

cat("==================== smoke: one fit, OLD vs NEW ====================\n")
g_data <- gen_data(300L)
theta_off <- c(tcl = TRUE_TCL, tv = TRUE_TV, bwt = 0.60, add.err = TRUE_ADD, eta.cl = TRUE_OM^2)
g_src <- gen_src(theta_off, 300L, seed = 99L)

f_old <- fit_ctl(c(g_data, g_src), use_fix = FALSE)
cat("OLD  bwt_hat =", f_old$parFixedDf["bwt","Estimate"], " se =", f_old$parFixedDf["bwt","SE"], "\n")

f_new <- fit_ctl(c(g_data, g_src), use_fix = TRUE)
cat("NEW  bwt_hat =", f_new$parFixedDf["bwt","Estimate"], " se =", f_new$parFixedDf["bwt","SE"], "\n")

cat("\n==================== n-invariance sweep, OLD vs NEW ====================\n")
for (nm in c(50L, 300L, 1500L)) {
  g_src_n <- gen_src(theta_off, nm, seed = 99L)
  fo <- fit_ctl(c(g_data, g_src_n), use_fix = FALSE)
  fn <- fit_ctl(c(g_data, g_src_n), use_fix = TRUE)
  cat(sprintf("n_model=%-5d   OLD bwt_hat=%.4f (se %.4f)   NEW bwt_hat=%.4f (se %.4f)\n",
              nm, fo$parFixedDf["bwt","Estimate"], fo$parFixedDf["bwt","SE"],
              fn$parFixedDf["bwt","Estimate"], fn$parFixedDf["bwt","SE"]))
}
