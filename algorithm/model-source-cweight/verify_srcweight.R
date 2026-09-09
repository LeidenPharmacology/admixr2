## Step 1 verification only: does .admSrcWeight() reproduce the M already
## computed and validated in algorithm_test_fullM.R? No wiring into .adghNLL
## yet -- this just checks the new helper function in isolation.
suppressMessages({ library(nlmixr2est); devtools::load_all(".", quiet = TRUE) })
rxode2::setRxThreads(2L)

TRUE_TCL <- log(5); TRUE_TV1 <- log(30)
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
theta_off <- c(tcl = TRUE_TCL - 0.15, tv = TRUE_TV1 - 0.10, add.err = 0.10, eta.cl = 0.06^2)
g_src <- gen_src_1c(theta_off, 300L, seed = 55L)

grp <- admixr2:::.admSrcGroups(g_src)
idx <- grp[[1L]]
w   <- admixr2:::.admSrcWeight(g_src, idx)
cat("rank:", w$rank, "\n")
cat("eigenvalues (top 6):", format(head(sort(w$values, decreasing = TRUE), 6), digits = 3), "\n")
cat("\nexpected from algorithm_test_fullM.R: rank 4, values 4.95e-01 3.92e-01 8.48e-02 2.75e-06\n")
