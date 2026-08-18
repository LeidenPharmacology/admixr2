## Are the covariate gradients correct, and is the implementation efficient?
##
##   G  analytical gradient vs a central finite difference of the SAME objective
##   S  is the sensitivity model actually used, or is it silently degrading?
##   E  cost: objective, gradient, and how they scale with the covariate work
##
## Run:  Rscript validation/covariate-gradient-efficiency.R
suppressMessages({library(rxode2); devtools::load_all(".", quiet = TRUE)})
say <- function(...) { cat(..., "\n"); utils::flush.console() }
ok  <- function(lbl, pass, d = "")
  cat(sprintf("  %-46s %s  %s\n", lbl, if (pass) "PASS" else "**FAIL**", d))

TT <- c(0.5, 1, 2, 4, 8, 12); DOSE <- 100
mod <- function() {
  ini({lcl <- log(3.0); lv <- log(20); bwt <- 0.70; bcr <- 0.55
       eta.cl ~ 0.060; prop.err <- 0.12})
  model({cl <- exp(lcl + eta.cl) * (WT/70)^bwt * (CRCL/90)^bcr
         v  <- exp(lv); d/dt(centr) <- -cl/v*centr; cp <- centr/v
         cp ~ prop(prop.err)})
}
ui   <- suppressMessages(rxode2::rxode2(mod))
ovar <- admixr2:::.admOutputVar(ui)
rxM  <- admixr2:::.admLoadModel(ui)
sM   <- tryCatch(admixr2:::.admLoadSensModel(ui), error = function(e) NULL)
say(sprintf("sensitivity model available: %s\n", !is.null(sM)))

rel  <- function(a, b) max(abs(a - b)/pmax(abs(b), 1e-12))
relV <- function(a, b) max(abs(a - b))/max(abs(b))

pop  <- covDist(WT = c(mean = 78, sd = 18), CRCL = c(mean = 95, sd = 26),
                cor = 0.45, dist = "lnorm")
Eobs <- c(4.42, 3.95, 3.20, 2.15, 1.02, 0.51)

setup <- function(studies, integ = "quadrature", nodes = 7L) {
  ctl <- adghControl(studies = studies, grad = "analytical", print = 0L,
                     covMethod = "none", n_nodes = nodes, cov_nodes = nodes,
                     cov_integration = integ)
  pin <- admixr2:::.admDriverPinfo(ui, ctl)
  stu <- admixr2:::.admCheckCovariates(
           ui, pin, admixr2:::.admDriverUnits(studies, ui, ovar)$studies,
           "analytical", "adgh")
  list(pin = pin, stu = stu,
       grid = admixr2:::.adghNodeGrid(nodes, pin$n_eta),
       p0 = admixr2:::.admBuildOptVec(pin)$p0)
}
mkstudy <- function(cd, strat = FALSE, n = 200L) {
  s <- list(E = Eobs, V = diag(0.05, length(TT)) + 0.01, n = n,
            times = TT, ev = rxode2::et(amt = DOSE), cov_dist = cd)
  s
}
gradcheck <- function(d, lbl, integ = "quadrature", h = 1e-5) {
  S <- setup(d, integ)
  f  <- function(p) admixr2:::.adghNLL(p, S$pin, S$stu, rxM, ovar, S$grid, 1L)
  ga <- admixr2:::.adghGrad(S$p0, S$pin, S$stu, sM, rxM, ovar, S$grid, 1L, 1e-4)
  gf <- vapply(seq_along(S$p0), function(j) {
    a <- b <- S$p0; a[j] <- a[j] + h; b[j] <- b[j] - h; (f(a) - f(b))/(2*h) }, 0)
  rel <- max(abs(ga - gf)/pmax(abs(gf), 1e-6))
  ok(lbl, rel < 2e-2, sprintf("max rel %.2e   |g| %.3f", rel, sqrt(sum(ga^2))))
  invisible(list(ga = ga, gf = gf, S = S, f = f))
}

say("=== G. analytical gradient vs central difference ===")
gradcheck(list(s = mkstudy(pop)),                     "marginal, quadrature")
gradcheck(list(s = c(mkstudy(pop), list(stratify = TRUE))),
          "study carrying `stratify` (expanded by datagen)")
gradcheck(list(s = mkstudy(pop)), "marginal, TAYLOR", integ = "taylor")
gradcheck(list(a = mkstudy(pop), b = mkstudy(
  covDist(WT = c(mean = 62, sd = 12), CRCL = c(mean = 55, sd = 15),
          cor = 0.45, dist = "lnorm"))), "two studies, quadrature")
## and a per-stratum study, as datagen would emit
st1 <- covStrata(pop, stratify = "WT", n_nodes = 3L, n = 200)
gradcheck(stats::setNames(lapply(st1, function(s)
  mkstudy(s$cov_dist, n = s$n)), paste0("k", seq_along(st1))),
  "3 strata, as datagen emits them")

say("\n=== E. cost of one objective and one gradient ===")
tm <- function(fn, r = 10L) system.time(for (i in seq_len(r)) fn())[["elapsed"]]/r*1000
bench <- function(d, lbl, integ = "quadrature", nodes = 7L) {
  S  <- setup(d, integ, nodes)
  to <- tm(function() admixr2:::.adghNLL(S$p0, S$pin, S$stu, rxM, ovar, S$grid, 1L))
  tg <- tm(function() admixr2:::.adghGrad(S$p0, S$pin, S$stu, sM, rxM, ovar,
                                          S$grid, 1L, 1e-4), 5L)
  np <- length(S$p0)
  say(sprintf("  %-38s obj %6.1f ms   grad %7.1f ms   ratio %5.1f  (FD would be %d)",
              lbl, to, tg, tg/to, np + 1L))
}
bench(list(s = mkstudy(pop)),              "1 study, marginal, 7 nodes")
bench(list(s = mkstudy(pop)),              "1 study, marginal, 11 nodes", nodes = 11L)
bench(list(s = mkstudy(pop)),              "1 study, TAYLOR", integ = "taylor")
bench(stats::setNames(lapply(st1, function(s) mkstudy(s$cov_dist, n = s$n)),
                      paste0("k", seq_along(st1))), "3 strata")
st2 <- covStrata(pop, stratify = "WT", n_nodes = 5L, n = 200)
bench(stats::setNames(lapply(st2, function(s) mkstudy(s$cov_dist, n = s$n)),
                      paste0("k", seq_along(st2))), "5 strata (the default)")

say("\n=== E. where the covariate work actually lands ===")
say(sprintf("%-30s %10s %10s", "grid rows per objective", "quadrature", "taylor"))
gq <- admixr2:::.admCovGrid(pop, 7L)
tq <- admixr2:::.admCovTaylorDesign(pop, 0.5)
say(sprintf("%-30s %10d %10d", "covariate points", nrow(gq$X), tq$n_pt))
say(sprintf("%-30s %10d %10d", "x eta nodes (7^1)", nrow(gq$X)*7L, tq$n_pt*7L))
say("\nstrata multiply the SOLVE count directly: each stratum is its own study.")

## ===========================================================================
## The same two questions for admc (Monte Carlo, per-subject covariate draws)
## ===========================================================================
say("\n\n=== admc: gradient and cost ===")
mc_setup <- function(studies, n_sim = 4000L) {
  ctl <- admControl(studies = studies, grad = "sens", n_sim = n_sim,
                    print = 0L, covMethod = "none")
  pin <- admixr2:::.admDriverPinfo(ui, ctl)
  stu <- admixr2:::.admCheckCovariates(
           ui, pin, admixr2:::.admDriverUnits(studies, ui, ovar)$studies, "sens")
  list(pin = pin, stu = stu,
       zl = admixr2:::.admMakeZ(n_sim, pin, length(studies), "sobol"),
       pl = admixr2:::.admMakeParamsList(n_sim, pin, length(studies)),
       p0 = admixr2:::.admBuildOptVec(pin)$p0,
       path = vapply(stu, function(s) s$.adm_cov_path %||% "-", ""))
}
mc_grad <- function(d, lbl, n_sim = 4000L, h = 1e-5) {
  S <- mc_setup(d, n_sim)
  f <- function(p) admixr2:::.admNLL(p, S$pin, S$stu, S$zl, rxM, ovar, S$pl, 1L)
  ga <- admixr2:::.admGrad(S$p0, S$pin, S$stu, S$zl, rxM, ovar, S$pl, 1L, 1e-4,
                           sensModel = sM)
  gf <- vapply(seq_along(S$p0), function(j) {
    a <- b <- S$p0; a[j] <- a[j] + h; b[j] <- b[j] - h; (f(a) - f(b))/(2*h) }, 0)
  rel <- max(abs(ga - gf)/pmax(abs(gf), 1e-6))
  ok(sprintf("%s  [path %s]", lbl, paste(unique(S$path), collapse = ",")),
     rel < 5e-2, sprintf("max rel %.2e", rel))
}
mc_grad(list(s = mkstudy(pop)), "1 study, marginal")
mc_grad(stats::setNames(lapply(st1, function(s) mkstudy(s$cov_dist, n = s$n)),
                        paste0("k", seq_along(st1))), "3 strata")

say("")
mc_bench <- function(d, lbl, n_sim = 4000L) {
  S  <- mc_setup(d, n_sim)
  to <- tm(function() admixr2:::.admNLL(S$p0, S$pin, S$stu, S$zl, rxM, ovar,
                                        S$pl, 1L), 5L)
  tg <- tm(function() admixr2:::.admGrad(S$p0, S$pin, S$stu, S$zl, rxM, ovar,
                                         S$pl, 1L, 1e-4, sensModel = sM), 3L)
  say(sprintf("  %-38s obj %6.1f ms   grad %7.1f ms   ratio %5.1f",
              lbl, to, tg, tg/to))
}
mc_bench(list(s = mkstudy(pop)), "1 study, n_sim 2000",  2000L)
mc_bench(list(s = mkstudy(pop)), "1 study, n_sim 4000",  4000L)
mc_bench(list(s = mkstudy(pop)), "1 study, n_sim 16000", 16000L)
mc_bench(stats::setNames(lapply(st1, function(s) mkstudy(s$cov_dist, n = s$n)),
                         paste0("k", seq_along(st1))), "3 strata, n_sim 4000")

say("\n=== do adgh and admc agree on the moments? ===")
A <- setup(list(s = mkstudy(pop)))
ma <- admixr2:::.adghMoments(admixr2:::.admUnpack(A$p0, A$pin), A$pin,
                             A$stu[[1L]], rxM, ovar, A$grid, 1L)
for (ns in c(4000L, 32000L)) {
  M <- mc_setup(list(s = mkstudy(pop)), ns)
  pr <- admixr2:::.admUnpack(M$p0, M$pin)
  cr <- admixr2:::.admStudyCovRows(M$stu[[1L]], M$pin, ns)
  eta <- M$zl[[1L]] %*% t(pr$L); colnames(eta) <- M$pin$eta_col_names
  cp <- admixr2:::.admSimulate(rxM, pr$struct, M$pin$sigma_names, eta, cr, ovar,
                               M$pl[[1L]], 1L, .Machine$integer.max, NULL)
  mu <- colMeans(cp); V <- crossprod(sweep(cp, 2L, mu))/ns
  arr <- admixr2:::.admUnitResidRows(M$pin, ovar, pr$sigma_var, length(mu))
  ap  <- admixr2:::.admResidApply(mu, diag(V), arr, cr$times, V)
  Em  <- ap$mu; Vm <- admixr2:::.admApplyResidTail(V, ap)
  ok(sprintf("adgh vs admc moments, n_sim %5d", ns),
     rel(ma$E, Em) < 5e-3 && relV(ma$V, Vm) < 1e-2,
     sprintf("relE %.2e  relV %.2e", rel(ma$E, Em), relV(ma$V, Vm)))
}
