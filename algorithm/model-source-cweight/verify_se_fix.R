## Verifies the two SE fixes for srcWeight = "cov" (HANDOFF.md item #1):
##  (a) covMethod = "r"'s Hessian no longer goes through the Mpinv-blind
##      .adghGradNLL gradient-FD path for a cov-route group.
##  (b) covMethod = "r,s"'s sandwich meat for a cov-route group is now built
##      from .admSrcMeatCov() (4 A' Mpinv A), not the n-route's .admSrcMeat().
##
## The check: under correct specification J = 2H (documented invariant, see
## R/adgh.R around .admSandwich()). A fit that is well-specified and has ONLY
## a cov-route source contribution should show sandwich_HJ$J ~= 2 * sandwich_HJ$H
## elementwise, and covMethod = "r" vs "r,s" standard errors should closely
## agree (both come from the same, now-consistent, objective).
suppressMessages({ library(nlmixr2est); devtools::load_all(".", quiet = TRUE) })
rxode2::setRxThreads(2L)

TRUE_TCL <- log(5); TRUE_TV1 <- log(30)
TRUE_ADD <- 0.10; TRUE_OM <- 0.06
DOSE <- 200; TAU <- 12
SRC_TIMES <- c(0.5, 1, 2, 4, 8, 11.5)
DENSE_TIMES <- c(0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 12, 24)

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

## Fit model is the SAME 1-cmt structure as the source -- correctly specified,
## so J = 2H is the right check (a misspecified fit need not satisfy it).
.fit_1c <- function() {
  ini({ tcl <- log(4); tv <- log(25); add.err <- 0.15; eta.cl ~ 0.1 })
  model({ cl <- exp(tcl + eta.cl); v <- exp(tv)
          cp <- linCmt(); cp ~ add(add.err) })
}
gen_data_1c <- function(n_data) {
  ui <- mk_src_ui(TRUE_TCL, TRUE_TV1, TRUE_ADD, TRUE_OM^2)
  ev <- rxode2::et(amt = DOSE, ii = TAU, ss = 1, addl = 0)
  sp <- list(times = DENSE_TIMES, ev = ev, n = n_data)
  suppressWarnings(suppressMessages(datagen(
    list(dat = sp), model = ui, control = datagenControl(method = "gh", seed = 2L))))
}

g_data <- gen_data_1c(300L)
theta_off <- c(tcl = TRUE_TCL - 0.10, tv = TRUE_TV1 - 0.08, add.err = TRUE_ADD, eta.cl = TRUE_OM^2)
g_src <- gen_src_1c(theta_off, 150L, seed = 55L)

cat("==================== covMethod = \"r\", srcWeight = \"cov\" ==============\n")
f_r <- suppressMessages(nlmixr2est::nlmixr2(.fit_1c, admData(), est = "adgh",
  control = adghControl(studies = c(g_data, g_src), print = 0L, cores = 2L,
                        covMethod = "r", grad = "fd", srcWeight = "cov")))
se_r <- f_r$parFixedDf[c("tcl", "tv", "add.err"), "SE"]
cat("SE (r):    "); print(se_r)

cat("\n================== covMethod = \"r,s\", srcWeight = \"cov\" =============\n")
f_rs <- suppressMessages(nlmixr2est::nlmixr2(.fit_1c, admData(), est = "adgh",
  control = adghControl(studies = c(g_data, g_src), print = 0L, cores = 2L,
                        covMethod = "r,s", grad = "fd", srcWeight = "cov")))
se_rs <- f_rs$parFixedDf[c("tcl", "tv", "add.err"), "SE"]
cat("SE (r,s):  "); print(se_rs)
cat("ratio r,s / r (expect close to 1 under correct specification):\n")
print(se_rs / se_r)

hj <- f_rs$env$admExtra$sandwich
if (is.null(hj)) {
  cat("\nsandwich did not engage (degraded to \"r\") -- see warnings above.\n")
} else {
  H <- hj$H; J <- hj$J
  cat("\nJ / (2H) elementwise, diagonal (expect ~1 under correct specification):\n")
  print(diag(J) / diag(2 * H))
  cat("max |J - 2H| / max |2H|:",
      max(abs(J - 2 * H)) / max(abs(2 * H)), "\n")
}
