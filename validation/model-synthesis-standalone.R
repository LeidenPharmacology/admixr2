## =============================================================================
## Combining published PK models of DIFFERENT STRUCTURE into one unified model
## =============================================================================
##
##   Rscript validation/model-synthesis-standalone.R
##
## Self-contained. No arguments, no input files. Writes four figures and a
## results table to the directory given by ADM_FIG_DIR (default: current dir).
##
## THE PROBLEM. Two published vancomycin models in paediatric patients. They do
## not share a structure and they do not share a parameterisation:
##
##   source A (2-compartment)  CL = 0.16*WT^0.97,  Vc fixed, Vp = 0.19*WT^1.07,
##                             Q = 0.13*WT^1.19;  17.5 mg/kg over 1 h, 9 samples
##                             over 0.5-6 h.  n = 14.
##   source B (1-compartment)  CL = 2.99*(WT/20),  V = 9.55*(WT/20);
##                             15 mg/kg q6h x5, samples at 21 and 23.5 h. n = 72.
##
## Vp and Q have NO counterpart in source B. And source B's profile is
## weight-INVARIANT by algebra: mg/kg dosing with CL and V both linear in weight
## cancels weight exactly, so B carries no covariate information at all. Both
## facts are demonstrated by the script rather than asserted.
##
## WHAT WE HAVE, AND WHAT IT IS FOR.
##   the model FORM        -> simulate predictions from it
##   the parameter VALUES  -> generate (E, V) over a population
##   covariate mean and SD -> DEFINE that population (cov_dist)
##   parameter SEs         -> the only honest uncertainty; propagate them
##
## There are no observations anywhere here, so there is no likelihood to be the
## likelihood OF. What we are doing is a PROJECTION: choose the unified model's
## parameters to best reproduce what the sources predict, over the patients and
## designs we care about. Numerically that is a weighted minimum-divergence fit,
## which is what admixr2's objective computes when fed datagen() output.
##
## TWO CHOICES THE ANALYST MUST MAKE DELIBERATELY, both exposed below:
##   * K, the number of covariate strata each source is projected over. K = 1
##     matches a source only at its population average; K > 1 matches it ACROSS
##     the covariate range, which is what you want if the unified model must
##     serve light and heavy patients alike. Strata add FIDELITY, not evidence.
##   * n_k, which are WEIGHTS, not sample sizes. They are kept summing to each
##     source's real enrolment so its influence matches its evidence. Inflating
##     them would shrink every standard error without adding information --
##     which is exactly why the uncertainty here comes from route (3) below and
##     never from the fit's own Hessian.
## =============================================================================

suppressMessages({library(rxode2); library(nlmixr2est)})
## Covariate marginalisation is required. Use the installed admixr2 if it has it,
## otherwise load the package this script ships inside -- so the script runs from
## a checkout with nothing installed.
.has_cov <- function()
  isTRUE(try(requireNamespace("admixr2", quietly = TRUE), silent = TRUE)) &&
  exists(".admCovRowsFor", envir = asNamespace("admixr2"), inherits = FALSE)
if (!.has_cov()) {
  .root <- tryCatch(dirname(dirname(normalizePath(
    sys.frame(1)$ofile %||% "validation/x"))), error = function(e) "..")
  if (!file.exists(file.path(.root, "DESCRIPTION"))) .root <- ".."
  if (!file.exists(file.path(.root, "DESCRIPTION"))) .root <- "."
  message("admixr2 with covariate support not installed; loading from ", .root)
  suppressMessages(devtools::load_all(.root, quiet = TRUE))
} else suppressMessages(library(admixr2))
if (!.has_cov()) stop("admixr2 with covariate marginalisation is required.")
FIG <- Sys.getenv("ADM_FIG_DIR", ".")
say <- function(...) { cat(...,"\n"); utils::flush.console() }
M_UNC  <- as.integer(Sys.getenv("ADM_M", "20"))   # source-uncertainty draws
K_MAIN <- 4L                                       # strata for the main fit
RSE    <- 0.15                                     # assumed RSE on source params

## ---- 1. the two published sources -------------------------------------------
lnp <- function(mu, sd)
  list(meanlog = log(mu^2 / sqrt(sd^2 + mu^2)), sdlog = sqrt(log(1 + sd^2/mu^2)))
SRC <- list(
  A = list(wt = lnp(20.2, 9.2), n = 14L, times = c(0.5,1,1.25,1.5,2,3,4,5,6),
           ev = rxode2::et(amt = 17.5, dur = 1)),
  B = list(wt = lnp(18.1, 8.5), n = 72L, times = c(21, 23.5),
           ev = rxode2::et(amt = 15, ii = 6, addl = 4, dur = 2)))
## source A: 2-compartment, weight on CL, Vp and Q
mk_A <- function(p) eval(parse(text = sprintf('function() {
  ini({lcl <- %.8f; lvc <- %.8f; lvp <- %.8f; lq <- %.8f
       clwt <- %.8f; vpwt <- %.8f; qwt <- %.8f
       eta.cl ~ 0.062; eta.vc ~ 0.048; prop.err <- 0.12})
  model({cl <- exp(lcl + eta.cl) * WT^clwt; vc <- exp(lvc + eta.vc)
         vp <- exp(lvp) * WT^vpwt; q <- exp(lq) * WT^qwt; f(centr) <- WT
         d/dt(centr) <- -(cl + q)/vc*centr + q/vp*peri
         d/dt(peri)  <-        q /vc*centr - q/vp*peri
         cp <- centr/vc; cp ~ prop(prop.err)})}',
  p[1], p[2], p[3], p[4], p[5], p[6], p[7])))
P_A <- c(log(0.16), log(3.86), log(0.19), log(0.13), 0.97, 1.07, 1.19)
## source B: 1-compartment, weight linear on CL and V  -> weight CANCELS
mk_B <- function(p) eval(parse(text = sprintf('function() {
  ini({lcl <- %.8f; lv <- %.8f; eta.cl ~ 0.0223; eta.v ~ 0.0134
       prop.err <- 0.119})
  model({cl <- exp(lcl + eta.cl)*(WT/20); v <- exp(lv + eta.v)*(WT/20)
         f(centr) <- WT
         d/dt(centr) <- -cl/v*centr; cp <- centr/v; cp ~ prop(prop.err)})}',
  p[1], p[2])))
P_B <- c(log(2.99), log(9.55))
## the TARGET: one unified 2-compartment model, started away from both sources
mod_unified <- function() {
  ini({lcl <- log(0.16) + 0.35; lvc <- log(3.86) - 0.25; lvp <- log(0.19) + 0.3
       lq <- log(0.13) - 0.3; clwt <- 0.60; vpwt <- 0.70; qwt <- 0.80
       eta.cl ~ 0.09; eta.vc ~ 0.07; prop.err <- 0.15})
  model({cl <- exp(lcl + eta.cl) * WT^clwt; vc <- exp(lvc + eta.vc)
         vp <- exp(lvp) * WT^vpwt; q <- exp(lq) * WT^qwt; f(centr) <- WT
         d/dt(centr) <- -(cl + q)/vc*centr + q/vp*peri
         d/dt(peri)  <-        q /vc*centr - q/vp*peri
         cp <- centr/vc; cp ~ prop(prop.err)})
}
PARS <- c("clwt", "vpwt", "qwt")           # what only source A can speak to
REF  <- c(clwt = 0.97, vpwt = 1.07, qwt = 1.19)

## ---- 2. project each source over K covariate strata --------------------------
## Each stratum keeps its OWN (truncated) covariate distribution -- never a point
## value. Plugging in the stratum mean is the ecological plug-in and biases the
## coefficient upward. n_k sum to the source's real enrolment.
strata_of <- function(mkfun, p, s, K, tag) {
  nk <- diff(round(seq(0, s$n, length.out = K + 1)))
  keep <- which(nk > 0)
  stats::setNames(lapply(keep, function(k) {
    qf <- local({ k <- k; K <- K; w <- s$wt
      function(u) stats::qlnorm(((k - 1) + u)/K, w$meanlog, w$sdlog) })
    list(model = mkfun(p), times = s$times, ev = s$ev, n = as.integer(nk[k]),
         cov_dist = list(WT = list(quantile = qf)))
  }), paste0(tag, "_s", keep))
}
build <- function(K, pA = P_A, pB = P_B)
  datagen(c(strata_of(mk_A, pA, SRC$A, K, "A"),
            strata_of(mk_B, pB, SRC$B, K, "B")),
          control = datagenControl(n_sim = 30000L, seed = 4321L))
project <- function(studies) {
  f <- try(suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
    mod_unified, admData(), est = "adgh",
    control = adghControl(studies = studies, n_nodes = 7L, grad = "analytical",
                          print = 0L, covMethod = "none", maxeval = 3000L)))),
    silent = TRUE)
  if (inherits(f, "try-error")) return(NULL)
  f$env$admExtra$struct
}

## ---- 3. does each source carry covariate information at all? -----------------
say("\n== does each source's profile depend on weight? ==")
spread <- function(K, tag) {
  g <- build(K); u <- g[grepl(paste0("^", tag), names(g))]
  E <- t(vapply(u, function(s) as.numeric(s$E), numeric(length(u[[1]]$E))))
  100 * max(apply(E, 2, function(z) diff(range(z))/mean(z)))
}
sprA <- spread(4L, "A"); sprB <- spread(4L, "B")
say(sprintf("  source A: %.1f%% spread across weight strata", sprA))
say(sprintf("  source B: %.1f%% spread across weight strata  %s", sprB,
            if (sprB < 1) "<- weight cancels; carries NO covariate information" else ""))

## ---- 4. fidelity: how many strata does the projection need? ------------------
say("\n== projection at increasing covariate resolution ==")
KS <- c(1L, 2L, 4L, 8L); fitK <- list()
for (K in KS) {
  s <- project(build(K)); fitK[[as.character(K)]] <- s
  say(sprintf("  K = %-2d  %s", K, if (is.null(s)) "FAILED" else
    paste(sprintf("%s %5.2f", PARS, vapply(PARS, function(p) s[[p]], 0)),
          collapse = "  ")))
}

## ---- 5. uncertainty: propagate the SOURCES' parameter uncertainty ------------
## The fit's own Hessian is meaningless here -- n_k are weights we chose, so
## doubling them would halve every SE. The only real uncertainty is the sources'.
## Published RSEs are assumed (see RSE) because the abstracts do not report them.
say(sprintf("\n== %d draws from the sources' parameter uncertainty (RSE %.0f%%) ==",
            M_UNC, 100*RSE))
## Draw every perturbation UP FRONT. datagen() sets the seed internally, so a
## draw taken inside the loop is reset to the same value on every iteration --
## which silently produced 20 identical "draws" and a 0.01-wide interval.
##
## And perturb on the right scale. lcl/lvc/lvp/lq are LOG parameters, so an RSE
## of r on the natural parameter is an ADDITIVE sd of r on the log scale
## (d log X = dX / X). The exponents are already on their natural scale, so they
## take a proportional perturbation.
set.seed(20260814)
is_log  <- c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE)   # P_A layout
pert_A  <- t(vapply(seq_len(M_UNC), function(m) {
  z <- stats::rnorm(length(P_A))
  ifelse(is_log, P_A + RSE * z, P_A * (1 + RSE * z))
}, numeric(length(P_A))))
pert_B  <- t(vapply(seq_len(M_UNC), function(m)
  P_B + RSE * stats::rnorm(length(P_B)), numeric(length(P_B))))
draws <- matrix(NA_real_, M_UNC, length(PARS), dimnames = list(NULL, PARS))
for (m in seq_len(M_UNC)) {
  pA <- pert_A[m, ]
  pB <- pert_B[m, ]
  s  <- project(build(K_MAIN, pA, pB))
  if (!is.null(s)) draws[m, ] <- vapply(PARS, function(p) s[[p]], 0)
  say(sprintf("  draw %2d/%d  %s", m, M_UNC, if (is.null(s)) "FAILED" else
    paste(sprintf("%s %5.2f", PARS, draws[m, ]), collapse = "  ")))
  saveRDS(list(fitK = fitK, draws = draws, sprA = sprA, sprB = sprB),
          file.path(FIG, "model-synthesis.rds"))
}
say("\n== unified model, with source-uncertainty intervals ==")
for (p in PARS) {
  v <- draws[, p]; v <- v[is.finite(v)]
  say(sprintf("  %-5s %5.2f   90%% interval [%.2f, %.2f]   (source A: %.2f)",
              p, stats::median(v), stats::quantile(v, .05),
              stats::quantile(v, .95), REF[[p]]))
}

## ---- 6. figures --------------------------------------------------------------
TD <- list(A = seq(0.1, 6, by = 0.1), B = seq(18, 24, by = 0.1))
pred_of <- function(struct, om = c(0.09, 0.07), sg = 0.15) {
  mk <- eval(parse(text = sprintf('function() {
    ini({lcl <- %.8f; lvc <- %.8f; lvp <- %.8f; lq <- %.8f
         clwt <- %.8f; vpwt <- %.8f; qwt <- %.8f
         eta.cl ~ %.8f; eta.vc ~ %.8f; prop.err <- %.8f})
    model({cl <- exp(lcl + eta.cl) * WT^clwt; vc <- exp(lvc + eta.vc)
           vp <- exp(lvp) * WT^vpwt; q <- exp(lq) * WT^qwt; f(centr) <- WT
           d/dt(centr) <- -(cl + q)/vc*centr + q/vp*peri
           d/dt(peri)  <-        q /vc*centr - q/vp*peri
           cp <- centr/vc; cp ~ prop(prop.err)})}',
    struct[["lcl"]], struct[["lvc"]], struct[["lvp"]], struct[["lq"]],
    struct[["clwt"]], struct[["vpwt"]], struct[["qwt"]], om[1], om[2], sg)))
  g <- datagen(list(
    A = list(model = mk, times = TD$A, ev = SRC$A$ev, n = SRC$A$n,
             cov_dist = list(WT = SRC$A$wt)),
    B = list(model = mk, times = TD$B, ev = SRC$B$ev, n = SRC$B$n,
             cov_dist = list(WT = SRC$B$wt))),
    control = datagenControl(n_sim = 20000L, seed = 99L))
  list(A = as.numeric(g$A$E), B = as.numeric(g$B$E))
}
src_pred <- function() {
  gA <- datagen(list(A = list(model = mk_A(P_A), times = TD$A, ev = SRC$A$ev,
                              n = SRC$A$n, cov_dist = list(WT = SRC$A$wt))),
                control = datagenControl(n_sim = 20000L, seed = 99L))
  gB <- datagen(list(B = list(model = mk_B(P_B), times = TD$B, ev = SRC$B$ev,
                              n = SRC$B$n, cov_dist = list(WT = SRC$B$wt))),
                control = datagenControl(n_sim = 20000L, seed = 99L))
  list(A = as.numeric(gA$A$E), B = as.numeric(gB$B$E))
}
SP <- src_pred(); UP <- pred_of(fitK[[as.character(K_MAIN)]])
COLK <- grDevices::colorRampPalette(c("#9ecae1", "#08306b"))(length(KS))

grDevices::png(file.path(FIG, "synthesis-1-sources.png"), 1900, 800, res = 150)
op <- par(mfrow = c(1,2), mar = c(4.6,4.6,3.4,1), mgp = c(2.7,.7,0), family="sans")
for (s in c("A","B")) {
  plot(TD[[s]], SP[[s]], type = "l", lwd = 3.5, col = "#c1121f", bty = "l",
       xlab = "time (h)", ylab = "vancomycin (mg/L)", ylim = c(0, max(SP[[s]])*1.1),
       main = sprintf("source %s design (%s)", s,
                      if (s=="A") "2-cmt, rich, single dose" else "1-cmt, sparse, steady state"))
  lines(TD[[s]], UP[[s]], lwd = 3.5, col = "#2a9d8f", lty = 1)
  points(SRC[[s]]$times, approx(TD[[s]], SP[[s]], SRC[[s]]$times)$y, pch = 21,
         bg = "black", col = "white", cex = 1.5, lwd = 2)
  legend("topright", c("published source model", "unified model (projection)",
                       "sampling times"),
         col = c("#c1121f","#2a9d8f","black"), lwd = c(3.5,3.5,NA),
         pch = c(NA,NA,21), pt.bg = "black", bty = "n", cex = .85)
}
invisible(par(op)); grDevices::dev.off()

grDevices::png(file.path(FIG, "synthesis-2-fidelity.png"), 1500, 800, res = 150)
par(mar = c(4.8,4.8,3.6,7.4), mgp = c(2.9,.7,0), family = "sans")
rel <- t(vapply(KS, function(K) {
  s <- fitK[[as.character(K)]]
  if (is.null(s)) rep(NA, 3) else vapply(PARS, function(p) s[[p]], 0)/REF[PARS]
}, numeric(3)))
PC <- c("#1b6ca8","#c1121f","#2a9d8f")
plot(NA, xlim = c(0.85, 10), ylim = range(rel, 1, na.rm = TRUE) * c(.95,1.05),
     log = "x", xaxt = "n", bty = "l", xlab = "K  (covariate strata per source)",
     ylab = "projected / source A value",
     main = "Strata buy FIDELITY of the projection, not evidence")
axis(1, at = KS, labels = KS); abline(h = 1, lty = 2, lwd = 2.5, col = "grey25")
for (j in seq_along(PARS)) {
  lines(KS, rel[, j], col = PC[j], lwd = 3.2)
  points(KS, rel[, j], pch = 21, bg = PC[j], col = "white", cex = 1.9, lwd = 2.4)
  text(10.4, rel[length(KS), j], PARS[j], adj = 0, col = PC[j], font = 2, xpd = NA)
}
grDevices::dev.off()

grDevices::png(file.path(FIG, "synthesis-3-uncertainty.png"), 1500, 800, res = 150)
par(mar = c(4.8,4.8,3.6,1.2), mgp = c(2.9,.7,0), family = "sans")
bx <- lapply(PARS, function(p) draws[is.finite(draws[,p]), p])
plot(NA, xlim = c(.5, 3.5), ylim = range(unlist(bx), REF, na.rm = TRUE)*c(.9,1.1),
     bty = "l", xaxt = "n", xlab = "", ylab = "unified model value",
     main = "Uncertainty comes from the SOURCES, not from the fit")
for (j in seq_along(PARS)) {
  graphics::boxplot(bx[[j]], at = j, add = TRUE, boxwex = .5, col = "#cfe3f2",
                    border = "grey25", outline = FALSE, yaxt = "n")
  points(jitter(rep(j, length(bx[[j]])), amount = .1), bx[[j]], pch = 16,
         col = grDevices::adjustcolor("black", .4), cex = .8)
  points(j, REF[[PARS[j]]], pch = 23, bg = "#c1121f", col = "white", cex = 2, lwd = 2)
}
axis(1, at = seq_along(PARS), labels = PARS, font = 2)
legend("topleft", c("projection over source-parameter draws", "source A value"),
       pch = c(22,23), pt.bg = c("#cfe3f2","#c1121f"), col = c("grey25","white"),
       pt.cex = 1.8, bty = "n", cex = .9)
grDevices::dev.off()

grDevices::png(file.path(FIG, "synthesis-4-information.png"), 1500, 780, res = 150)
par(mar = c(4.8,5.2,3.6,1.2), mgp = c(3.0,.7,0), family = "sans")
bp <- barplot(c(A = sprA, B = sprB), col = c("#2a9d8f","#8a8f98"), border = NA,
              ylim = c(0, max(sprA, 1) * 1.25),
              ylab = "spread of predicted mean across\nweight strata (%)",
              names.arg = c("source A\n(2-cmt)", "source B\n(1-cmt, mg/kg)"),
              main = "Only source A carries covariate information")
text(bp, c(sprA, sprB) + max(sprA)*0.04, sprintf("%.1f%%", c(sprA, sprB)),
     font = 2, cex = 1.1)
text(bp[2], max(sprA)*0.45,
     "mg/kg dosing with CL and V\nboth linear in weight\n=> weight cancels exactly",
     cex = .9, col = "grey30")
grDevices::dev.off()

say(sprintf("\nwrote 4 figures + model-synthesis.rds to %s", normalizePath(FIG)))
