## First real (non-monkeypatched) test of srcWeight = "cov", through the
## genuine public nlmixr2()/adghControl() API.
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

g_data <- gen_data_2c(300L)
theta_off <- c(tcl = TRUE_TCL - 0.15, tv = TRUE_TV1 - 0.10, add.err = TRUE_ADD, eta.cl = TRUE_OM^2)

cat("==================== srcWeight = \"cov\", n=300, REAL API =============\n")
g_src <- gen_src_1c(theta_off, 300L, seed = 55L)
f <- suppressMessages(nlmixr2est::nlmixr2(.fit_2c, admData(), est = "adgh",
  control = adghControl(studies = c(g_data, g_src), print = 0L, cores = 2L,
                        covMethod = "none", grad = "fd", srcWeight = "cov")))
p <- f$parFixedDf
cat(sprintf("tcl=%.4f  tv=%.4f  q=%.4f  v2=%.4f\n",
            p["tcl","Estimate"], p["tv","Estimate"], p["q","Estimate"], p["v2","Estimate"]))
cat("expected (from the monkey-patched test): tcl=1.6088 tv=3.4012 q=2.0793 v2=4.0950\n\n")

cat("==================== n-invariance under the REAL wiring ================\n")
for (nm in c(50L, 1500L, 10000L)) {
  g_src_n <- gen_src_1c(theta_off, nm, seed = 55L)
  fn <- suppressMessages(nlmixr2est::nlmixr2(.fit_2c, admData(), est = "adgh",
    control = adghControl(studies = c(g_data, g_src_n), print = 0L, cores = 2L,
                          covMethod = "none", grad = "fd", srcWeight = "cov")))
  pn <- fn$parFixedDf
  cat(sprintf("n_model=%-6d  tcl=%.4f  tv=%.4f  q=%.4f  v2=%.4f\n",
              nm, pn["tcl","Estimate"], pn["tv","Estimate"], pn["q","Estimate"], pn["v2","Estimate"]))
}

cat("\n==================== srcWeight = \"n\" (default, unchanged?) ===========\n")
fo <- suppressMessages(nlmixr2est::nlmixr2(.fit_2c, admData(), est = "adgh",
  control = adghControl(studies = c(g_data, g_src), print = 0L, cores = 2L,
                        covMethod = "none", grad = "fd")))   # srcWeight defaults to "n"
po <- fo$parFixedDf
cat(sprintf("tcl=%.4f  tv=%.4f  q=%.4f  v2=%.4f\n",
            po["tcl","Estimate"], po["tv","Estimate"], po["q","Estimate"], po["v2","Estimate"]))
