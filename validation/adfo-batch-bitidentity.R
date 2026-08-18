# ---------------------------------------------------------------------------
# validation/adfo-batch-bitidentity.R
#
# Does stacking N structural-theta configurations as ROWS of one rxSolve give
# BIT-IDENTICAL per-row results to solving each configuration on its own?
#
# This is the premise behind every "batch it into one solve" change in adfo.
# CLAUDE.md records a REVERTED experiment that merged same-`ev` studies over the
# UNION of their observation times: that changed the output grid, so rxode2's
# adaptive solver took different steps and predictions moved ~1e-7.  Adding
# SUBJECTS (rows) is a different operation -- each subject keeps its own event
# table and the same output grid -- but "different operation" is a hypothesis,
# not a result, so it is measured here rather than assumed.
#
# CLAUDE.md also warns that a single-subject probe misleadingly shows 0.000e+00,
# so this compares MANY configurations spread over a wide range of thetas, and
# reports the worst relative difference AND an exact identical() verdict.
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
sens2 <- admixr2:::.admLoadSensModel(ui, order = 2L)
rxMod <- admixr2:::.admLoadModel(ui)
rxode2::rxLoad(rxMod)

times <- c(0.25, 0.5, 1, 2, 4, 8, 12, 24)
s <- admixr2:::.admNormaliseStudy(
  list(E = rep(1, length(times)), V = diag(length(times)), n = 200L,
       times = times, ev = rxode2::et(amt = 100)), "s")
s$ev_full <- s$ev |> rxode2::et(s$times)

set.seed(3)
base   <- c(tcl = 4, tv = 30, tq = 2, tvp = 40)
n_cfg  <- 40L
sm     <- matrix(rep(base, each = n_cfg), n_cfg, 4L,
                 dimnames = list(NULL, names(base)))
sm     <- sm * exp(matrix(rnorm(n_cfg * 4L, 0, 0.35), n_cfg, 4L))

pars_sig <- admixr2:::.admUnpack(admixr2:::.admBuildOptVec(pinfo)$p0, pinfo)$sigma_var

report <- function(label, big, one_of) {
  ok_all <- TRUE; worst <- 0
  for (k in seq_len(n_cfg)) {
    a <- big[[k]]; b <- one_of(k)
    for (fld in c("mu", "J")) {
      x <- a[[fld]]; y <- b[[fld]]
      if (!identical(x, y)) ok_all <- FALSE
      d <- max(abs(as.numeric(x) - as.numeric(y)) /
                 pmax(abs(as.numeric(y)), 1e-300))
      worst <- max(worst, d)
    }
    # second-order block too, when present
    if (!is.null(a$dJ)) {
      for (nm in names(a$dJ)) {
        if (!identical(a$dJ[[nm]], b$dJ[[nm]])) ok_all <- FALSE
        worst <- max(worst, max(abs(a$dJ[[nm]] - b$dJ[[nm]]) /
                                  pmax(abs(b$dJ[[nm]]), 1e-300)))
      }
    }
  }
  cat(sprintf("%-28s identical=%-5s  worst rel diff = %.3e\n",
              label, ok_all, worst))
  ok_all
}

cat("\n40 configurations, 2-cmt ODE, 8 observation times\n\n")

# ---- sens-model branch (the one adfo actually uses) ------------------------
big_sens <- admixr2:::.adfoGetMuJBatch(sm, pinfo, s, sens2, rxMod, "cp", 1L, pars_sig)
report("sens batch vs one-at-a-time", big_sens, function(k)
  admixr2:::.adfoGetMuJBatch(sm[k, , drop = FALSE], pinfo, s, sens2, rxMod, "cp",
                             1L, pars_sig)[[1L]])

# ---- FD fallback branch (no sens model) ------------------------------------
big_fd <- admixr2:::.adfoGetMuJBatch(sm, pinfo, s, NULL, rxMod, "cp", 1L, NULL)
report("FD batch vs one-at-a-time", big_fd, function(k)
  admixr2:::.adfoGetMuJBatch(sm[k, , drop = FALSE], pinfo, s, NULL, rxMod, "cp",
                             1L, NULL)[[1L]])

# ---- does batch SIZE change the answer? ------------------------------------
half <- admixr2:::.adfoGetMuJBatch(sm[1:20, , drop = FALSE], pinfo, s, sens2,
                                   rxMod, "cp", 1L, pars_sig)
cat(sprintf("%-28s identical=%s\n", "size 40 vs size 20 (rows 1:20)",
            identical(lapply(big_sens[1:20], function(z) z[c("mu", "J")]),
                      lapply(half,           function(z) z[c("mu", "J")]))))
