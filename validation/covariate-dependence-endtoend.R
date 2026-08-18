## END-TO-END: dependent covariates through a real fit.
##
## The unit checks show the grid and the sampler describe the same covariate
## distribution. This asks the question that actually matters: do adgh
## (quadrature grid), adgh (Taylor) and admc (per-subject draws) agree on the
## PREDICTED MOMENTS of a real rxode2 model when the covariates are correlated,
## and does the dependence change the answer enough that agreeing means
## something?
##
## Run:  Rscript validation/covariate-dependence-endtoend.R
suppressMessages({library(rxode2); devtools::load_all(".", quiet = TRUE)})
say <- function(...) { cat(..., "\n"); utils::flush.console() }
kv  <- function(k, ...) cat(sprintf("%-46s %s\n", k, paste(..., collapse = "  ")))

TT <- c(0.5, 1, 2, 4, 8, 12); EV <- rxode2::et(amt = 100)
mod <- function() {
  ini({lcl <- log(3.0); lv <- log(20); clwt <- 0.75; clcr <- 0.35
       eta.cl ~ 0.09; prop.err <- 0.15})
  model({cl <- exp(lcl + eta.cl) * (WT/70)^clwt * (CRCL/90)^clcr
         v  <- exp(lv) * (WT/70)
         d/dt(centr) <- -cl / v * centr; cp <- centr / v
         cp ~ prop(prop.err)})
}
ui   <- suppressMessages(rxode2::rxode2(mod))
ovar <- admixr2:::.admOutputVar(ui)

cov_dist_for <- function(rho)
  admixr2:::.admCovDistCanon(list(WT   = list(mean = 72, sd = 16),
                                  CRCL = list(mean = 90, sd = 25), cor = rho))

## ---- reference: per-subject solve over a large dependent sample -------------
rx <- suppressMessages(rxode2::rxode2(ui$simulationModel))
ref_moments <- function(rho, n = 120000L) {
  cr <- admixr2:::.admCovRowsFor(cov_dist_for(rho), n, 1L)
  z  <- admixr2:::.admMakeZ(n, list(n_eta = 1L), 1L, "sobol")[[1L]]
  p  <- data.frame(lcl = log(3.0), lv = log(20), clwt = 0.75, clcr = 0.35,
                   prop.err = 0.15, rxerr.cp = 0, eta.cl = 0.30 * z[, 1L],
                   WT = cr[, "WT"], CRCL = cr[, "CRCL"])
  o  <- as.data.frame(rxode2::rxSolve(rx, params = p,
         events = EV |> rxode2::et(TT), cores = 1L,
         nDisplayProgress = .Machine$integer.max))
  o  <- o[o$time %in% TT, ]
  f  <- matrix(o$cp, nrow = n, byrow = TRUE)
  mu <- colMeans(f); V <- crossprod(sweep(f, 2L, mu)) / n
  diag(V) <- diag(V) + 0.15^2 * (mu^2 + diag(V))     # law of total variance
  list(E = mu, V = V)
}

## ---- adgh, through the real driver plumbing --------------------------------
adgh_moments <- function(rho, integration, nodes = 9L) {
  st <- list(s = list(E = rep(1, length(TT)), V = diag(length(TT)) + 0.01,
                      n = 100L, times = TT, ev = EV, cov_dist = cov_dist_for(rho)))
  ctl <- adghControl(studies = st, grad = "analytical", print = 0L,
                     covMethod = "none", n_nodes = 9L, cov_nodes = nodes,
                     cov_integration = integration)
  pinfo <- admixr2:::.admDriverPinfo(ui, ctl)
  u  <- admixr2:::.admCheckCovariates(
          ui, pinfo, admixr2:::.admDriverUnits(st, ui, ovar)$studies,
          "analytical", "adgh")
  s  <- u[[1L]]
  gr <- admixr2:::.adghNodeGrid(9L, pinfo$n_eta)
  pars <- admixr2:::.admUnpack(admixr2:::.admBuildOptVec(pinfo)$p0, pinfo)
  m  <- admixr2:::.adghMoments(pars, pinfo, s, admixr2:::.admLoadModel(ui),
                               ovar, gr, 1L)
  c(m, list(path = s$.adm_cov_path))
}

relE <- function(a, b) max(abs(a - b) / abs(b))
relV <- function(a, b) max(abs(a - b)) / max(abs(b))

say("======================================================================")
say(" 1. does adgh now ACCEPT a dependent covariate distribution?")
say("======================================================================")
r <- tryCatch(adgh_moments(0.6, "quadrature"), error = function(e)
  paste("REFUSED:", substr(conditionMessage(e), 1, 70)))
kv("adgh + cor = 0.6", if (is.character(r)) r else paste("accepted, path =", r$path))

say("\n======================================================================")
say(" 2. predicted moments vs the per-subject reference")
say("======================================================================")
say(sprintf("%6s %14s %8s %13s %13s", "rho", "arm", "nodes", "relE", "relV"))
for (rho in c(0.0, 0.6, 0.85)) {
  rf <- ref_moments(rho)
  for (nn in c(5L, 9L, 15L)) {
    a <- adgh_moments(rho, "quadrature", nn)
    say(sprintf("%6.2f %14s %8d %13.3e %13.3e", rho, "quadrature", nn,
                relE(a$E, rf$E), relV(a$V, rf$V)))
  }
  a <- adgh_moments(rho, "taylor")
  say(sprintf("%6.2f %14s %8s %13.3e %13.3e", rho, "taylor", "1+2p",
              relE(a$E, rf$E), relV(a$V, rf$V)))
}
say("\n(reference = 120k per-subject solves through the SAME joint sampler;")
say(" its own MC error is ~1e-3, so that is the floor for every row)")

say("\n======================================================================")
say(" 3. does the dependence MATTER? (else agreement proves nothing)")
say("======================================================================")
r0 <- ref_moments(0.0); r8 <- ref_moments(0.85)
kv("reference, rho 0 vs 0.85 (same margins)",
   sprintf("E rel %.4f   V rel %.4f", relE(r8$E, r0$E), relV(r8$V, r0$V)))
q0 <- adgh_moments(0.0, "quadrature", 9L); q8 <- adgh_moments(0.85, "quadrature", 9L)
kv("adgh quadrature, rho 0 vs 0.85",
   sprintf("E rel %.4f   V rel %.4f", relE(q8$E, q0$E), relV(q8$V, q0$V)))
t0 <- adgh_moments(0.0, "taylor"); t8 <- adgh_moments(0.85, "taylor")
kv("adgh taylor, rho 0 vs 0.85",
   sprintf("E rel %.4f   V rel %.4f", relE(t8$E, t0$E), relV(t8$V, t0$V)))
say("\nIf a path IGNORED the correlation these would be ~0 while the reference")
say("moved -- which is exactly what the old product-over-margins grid did.")

say("\n======================================================================")
say(" 4. cost")
say("======================================================================")
say(sprintf("%14s %8s %10s %12s", "arm", "nodes", "grid rows", "ms/moments"))
for (nn in c(5L, 9L, 15L)) {
  cg <- admixr2:::.admCovGrid(cov_dist_for(0.6), nn)
  tm <- system.time(for (i in 1:5) adgh_moments(0.6, "quadrature", nn))[["elapsed"]]/5*1000
  say(sprintf("%14s %8d %10d %12.1f", "quadrature", nn, nrow(cg$X), tm)) }
td <- admixr2:::.admCovTaylorDesign(cov_dist_for(0.6), 0.5)
tm <- system.time(for (i in 1:5) adgh_moments(0.6, "taylor"))[["elapsed"]]/5*1000
say(sprintf("%14s %8s %10d %12.1f", "taylor", "1+2p", td$n_pt, tm))
