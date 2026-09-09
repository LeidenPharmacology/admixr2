## The CORRECTED general fix: same functional form nll_cov_cpp already uses
## (a quadratic residual against a covariance), but the covariance is built
## from the SOURCE's own Jacobian + C_src -- M_E = D_E C_src D_E' -- projected
## with a pseudo-inverse (handles the rank deficiency: a 2-parameter source
## informs at most 2 directions of a longer mean-profile vector). No n. No
## parameter-name matching -- the comparison lives entirely in PREDICTED
## MOMENT SPACE (concentration at the source's own design points), which is
## meaningful for any two models regardless of structure.
##
## Deliberately scoped to the MEAN residual only (drop the V-matching term) to
## keep this a small, legible test -- and because M_E needs only the SOURCE's
## own Jacobian, not the analysis model's live sensitivity, so it can be built
## ONCE and does not need to move during optimization.
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

## ---- D_E: source's own mean-profile Jacobian w.r.t. (tcl, tv), at its own SS design
solve_src_mean <- function(tcl, tv) {
  ui <- mk_src_ui(tcl, tv, TRUE_ADD, TRUE_OM^2)
  ev <- rxode2::et(amt = DOSE, ii = TAU, ss = 1, addl = 0)
  ev <- rxode2::add.sampling(ev, SRC_TIMES)
  s  <- suppressMessages(rxode2::rxSolve(ui, ev, params = c(eta.cl = 0)))
  as.numeric(s$cp)
}
h <- 1e-4
D_E <- cbind(
  (solve_src_mean(TRUE_TCL + h, TRUE_TV1 - 0.10) - solve_src_mean(TRUE_TCL - h, TRUE_TV1 - 0.10)) / (2*h),
  (solve_src_mean(TRUE_TCL - 0.15, TRUE_TV1 + h) - solve_src_mean(TRUE_TCL - 0.15, TRUE_TV1 - h)) / (2*h)
)
M_E <- D_E %*% CSRC_1C[c("tcl","tv"), c("tcl","tv")] %*% t(D_E)
ev_ME <- eigen((M_E + t(M_E))/2, symmetric = TRUE)
cat("M_E eigenvalues:", format(ev_ME$values, digits = 3), "\n")
keep <- ev_ME$values > max(ev_ME$values) * 1e-8
cat("rank kept:", sum(keep), "of", length(ev_ME$values), "\n")
Mpinv <- ev_ME$vectors[, keep, drop = FALSE] %*%
  diag(1 / ev_ME$values[keep], sum(keep)) %*%
  t(ev_ME$vectors[, keep, drop = FALSE])

## ---- patched objective: mean-only M_E-weighted term for the source block ---
.adghNLL_momentM <- function(p, pinfo, studies, rxMod, out_var, grid, cores) {
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
    r   <- as.numeric(m$E) - as.numeric(s$E)
    total <- total + as.numeric(t(r) %*% Mpinv %*% r)
  }
  total
}
.orig_adghNLL <- get(".adghNLL", envir = asNamespace("admixr2"))
fit_ctl <- function(g, fn_new) {
  assignInNamespace(".adghNLL", fn_new, ns = "admixr2")
  on.exit(assignInNamespace(".adghNLL", .orig_adghNLL, ns = "admixr2"))
  suppressMessages(nlmixr2est::nlmixr2(.fit_2c, admData(), est = "adgh",
    control = adghControl(studies = g, print = 0L, cores = 2L, covMethod = "none", grad = "fd")))
}

cat("\n==================== moment-space M_E fix, no n at all ================\n")
g_data <- gen_data_2c(300L)
theta_off <- c(tcl = TRUE_TCL - 0.15, tv = TRUE_TV1 - 0.10, add.err = TRUE_ADD, eta.cl = TRUE_OM^2)
g_src <- gen_src_1c(theta_off, 300L, seed = 55L)   # n_model is now IRRELEVANT to this term
f <- fit_ctl(c(g_data, g_src), .adghNLL_momentM)
p <- f$parFixedDf
cat(sprintf("tcl=%.4f  tv=%.4f  q=%.4f  v2=%.4f\n",
            p["tcl","Estimate"], p["tv","Estimate"], p["q","Estimate"], p["v2","Estimate"]))
cat(sprintf("truth:  tcl=%.4f tv=%.4f q=%.4f v2=%.4f\n", TRUE_TCL, TRUE_TV1, TRUE_Q, TRUE_V2))
cat(sprintf("source reported: tcl=%.4f tv=%.4f\n", theta_off[["tcl"]], theta_off[["tv"]]))
cat("\n(compare against OLD's n=1500 collapse: q=-1.264 (Q~0.28), v2=0.632 (V2~1.9))\n")

cat("\n==================== n-invariance check (M_E fix never reads s$n) =====\n")
for (nm in c(50L, 300L, 1500L, 10000L)) {
  g_src_n <- gen_src_1c(theta_off, nm, seed = 55L)
  fN <- fit_ctl(c(g_data, g_src_n), .adghNLL_momentM)
  pN <- fN$parFixedDf
  cat(sprintf("n_model=%-6d  tcl=%.4f  tv=%.4f  q=%.4f  v2=%.4f\n",
              nm, pN["tcl","Estimate"], pN["tv","Estimate"], pN["q","Estimate"], pN["v2","Estimate"]))
}

cat("\n==================== replicate check: random source draws =============\n")
R <- 12L
out <- matrix(NA_real_, R, 4); colnames(out) <- c("tcl","tv","q","v2")
for (r in seq_len(R)) {
  set.seed(7000L + r)
  th_r <- c(tcl = TRUE_TCL + SD_SRC[["tcl"]] * rnorm(1),
            tv  = TRUE_TV1 + SD_SRC[["tv"]]  * rnorm(1),
            add.err = TRUE_ADD, eta.cl = TRUE_OM^2)
  g_src_r <- gen_src_1c(th_r, 300L, seed = 6000L + r)
  fr <- try(fit_ctl(c(g_data, g_src_r), .adghNLL_momentM), silent = TRUE)
  if (inherits(fr, "try-error")) { cat(sprintf("rep %2d FAILED\n", r)); next }
  pr <- fr$parFixedDf
  out[r, ] <- c(pr["tcl","Estimate"], pr["tv","Estimate"], pr["q","Estimate"], pr["v2","Estimate"])
  cat(sprintf("rep %2d  src_tcl=%.3f src_tv=%.3f -> tcl=%.4f tv=%.4f q=%.4f v2=%.4f\n",
              r, th_r[["tcl"]], th_r[["tv"]], out[r,1], out[r,2], out[r,3], out[r,4]))
}
truth <- c(TRUE_TCL, TRUE_TV1, TRUE_Q, TRUE_V2)
cat("\nbias:  ", sprintf("%.4f", colMeans(out, na.rm = TRUE) - truth), "\n")
cat("sd:    ", sprintf("%.4f", apply(out, 2, sd, na.rm = TRUE)), "\n")
