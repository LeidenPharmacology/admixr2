## Where the Monte Carlo saving actually is: ONE ensemble across MANY studies.
##
## The proposal was that "we can do only one set of MC simulations and then
## calculate the likelihood conditional on covariate values by simply weighting".
## Benchmarking it on a SINGLE study misses the point entirely -- one study needs
## one ensemble either way, so there is nothing to reuse and no saving, which is
## what an earlier benchmark here measured and wrongly concluded from.
##
## The reuse is ACROSS STUDIES. And that is not a niche case: between-study
## variation in the covariate DISTRIBUTION is the only thing that identifies a
## covariate coefficient from aggregate data at all (a single population leaves
## the likelihood exactly flat along tcl' = tcl + (b-b')*mu_a). So the studies
## that make a covariate estimable are, by construction, studies differing in
## their covariate distribution -- and when they share a design, one ensemble
## solved at a common reference serves all of them, each with its own weights.
##
## Requires: same times and same dosing, so the ensemble is transferable; and a
## proposal wide enough to cover EVERY study's u-distribution, not just one.
##
## Run:  Rscript validation/covariate-reweighting-shared-ensemble.R
suppressMessages({library(rxode2)})
suppressMessages(devtools::load_all(".", quiet = TRUE))
kv <- function(k, v) cat(sprintf("%s\t%s\n", k, paste(v, collapse = ",")))

TT <- c(0.5, 1, 2, 4, 8); EV <- rxode2::et(amt = 100)
TCL <- log(1); TV <- log(10); TCOV <- 0.75; OM <- 0.30; ADD <- 0.30
## S studies, SAME design, DIFFERENT weight distributions -- the identifying
## MBMA shape.
POPS <- list(
  A = list(meanlog = log(50), sdlog = 0.22),
  B = list(meanlog = log(70), sdlog = 0.30),
  C = list(meanlog = log(95), sdlog = 0.18),
  D = list(meanlog = log(60), sdlog = 0.35),
  E = list(meanlog = log(80), sdlog = 0.25))

mod <- function() {
  ini({tcl <- log(1); tv <- log(10); tcov <- 0.75; add.err <- 0.3
       eta.cl ~ 0.09})
  model({cl <- exp(tcl + eta.cl) * (WT / 70)^tcov
         v  <- exp(tv)
         d/dt(centr) <- -cl / v * centr; cp <- centr / v
         cp ~ add(add.err)})
}
ui   <- suppressMessages(rxode2::rxode2(mod))
ovar <- admixr2:::.admOutputVar(ui)
rxm  <- admixr2:::.admLoadModel(ui)

mk_studies <- function(k) setNames(lapply(seq_len(k), function(i)
  list(E = rep(1, length(TT)), V = diag(length(TT)) + 0.01, n = 100L,
       times = TT, ev = EV, cov_dist = list(WT = POPS[[i]]))), names(POPS)[seq_len(k)])

## ---- exact reference per study, nested quadrature on the closed form --------
gh <- function(k) {
  i <- seq_len(k - 1L)
  J <- matrix(0, k, k); J[cbind(i, i + 1L)] <- sqrt(i); J[cbind(i + 1L, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE); o <- order(e$values)
  list(x = e$values[o], w = e$vectors[1L, o]^2)
}
QA <- gh(50L); QE <- gh(50L)
ref_of <- function(pp) {
  m1 <- 0; M2 <- 0
  for (i in seq_along(QA$x)) {
    wt <- exp(pp$meanlog + pp$sdlog * QA$x[i])
    for (j in seq_along(QE$x)) {
      w  <- QA$w[i] * QE$w[j]
      cl <- exp(TCL + OM * QE$x[j]) * (wt / 70)^TCOV
      y  <- 100 / exp(TV) * exp(-cl / exp(TV) * TT)
      m1 <- m1 + w * y; M2 <- M2 + w * outer(y, y)
    }
  }
  V <- M2 - outer(m1, m1); diag(V) <- diag(V) + ADD^2
  list(E = m1, V = V)
}
REF <- lapply(POPS, ref_of)

## ---- per-study ensembles: what admc does today -----------------------------
per_study <- function(studies, n_sim) {
  ctl   <- admControl(studies = studies, grad = "sens", n_sim = n_sim,
                      print = 0L, covMethod = "none")
  pinfo <- admixr2:::.admDriverPinfo(ui, ctl)
  u  <- admixr2:::.admCheckCovariates(
          ui, pinfo, admixr2:::.admDriverUnits(studies, ui, ovar)$studies, "sens")
  pars <- admixr2:::.admUnpack(admixr2:::.admBuildOptVec(pinfo)$p0, pinfo)
  zl <- admixr2:::.admMakeZ(n_sim, pinfo, length(u), "sobol")
  pl <- admixr2:::.admMakeParamsList(n_sim, pinfo, length(u))
  out <- vector("list", length(u)); rows <- 0L
  for (i in seq_along(u)) {
    s  <- admixr2:::.admStudyCovRows(u[[i]], pinfo, n_sim)
    et <- zl[[i]] %*% t(pars$L); colnames(et) <- pinfo$eta_col_names
    cp <- admixr2:::.admSimulate(rxm, pars$struct, pinfo$sigma_names, et, s,
                                 ovar, pl[[i]], 1L, .Machine$integer.max, NULL)
    rows <- rows + nrow(cp)
    mu <- colMeans(cp); V <- crossprod(sweep(cp, 2L, mu)) / nrow(cp)
    ar <- admixr2:::.admUnitResidRows(pinfo, ovar, pars$sigma_var, length(TT))
    ap <- admixr2:::.admResidApply(mu, diag(V), ar, s$times, V)
    out[[i]] <- list(E = ap$mu, V = admixr2:::.admApplyResidTail(V, ap))
  }
  list(m = out, rows = rows)
}

## ---- ONE shared ensemble, reweighted per study ------------------------------
shared <- function(studies, n_sim) {
  ctl   <- admControl(studies = studies, grad = "sens", n_sim = n_sim,
                      print = 0L, covMethod = "none")
  pinfo <- admixr2:::.admDriverPinfo(ui, ctl)
  u  <- admixr2:::.admCheckCovariates(
          ui, pinfo, admixr2:::.admDriverUnits(studies, ui, ovar)$studies, "sens")
  pars <- admixr2:::.admUnpack(admixr2:::.admBuildOptVec(pinfo)$p0, pinfo)
  ## A COMMON reference covariate for every study. Delta is measured relative to
  ## whatever `cov` the study carries, so leaving each study on its own reference
  ## puts the u-coordinates on different origins and the shared ensemble is then
  ## being reweighted onto the wrong axis -- measured as 20-44% error, with the
  ## effective sample size still at 95-99%, which is one more demonstration that
  ## ESS does not detect a coverage or alignment fault.
  a_common <- mean(vapply(u, function(s)
    as.numeric(s[["cov"]][["WT"]] %||% admixr2:::.admCovMeanOf(s$cov_dist$WT)),
    numeric(1)))
  u <- lapply(u, function(s) { s[["cov"]] <- list(WT = a_common); s })
  sh <- lapply(u, function(s)
    admixr2:::.admCovShift(ui, rxm, pars, pinfo, s, 1L))
  if (any(vapply(sh, is.null, logical(1)))) return(NULL)
  ## ONE proposal covering EVERY study: centred on the pooled mean shift, with
  ## variance large enough for the widest u-distribution plus the spread of the
  ## study means. A proposal fitted to one study would under-cover the others,
  ## and under-coverage -- unlike inefficiency -- is what produces real error.
  mus <- vapply(sh, `[[`, numeric(1), "mu"); sds <- vapply(sh, `[[`, numeric(1), "sd")
  mu_p <- mean(mus)
  s_p  <- sqrt(max(sds)^2 + OM^2 + max((mus - mu_p)^2))
  z  <- admixr2:::.admMakeZ(n_sim, pinfo, 1L, "sobol")[[1L]]
  eta <- matrix(mu_p + s_p * z[, 1L], ncol = 1L,
                dimnames = list(NULL, pinfo$eta_col_names))
  s0 <- u[[1L]]; s0$cov_rows <- NULL
  s0[["cov"]] <- list(WT = sh[[1L]]$a_ref)
  pl <- admixr2:::.admMakeParamsList(n_sim, pinfo, 1L)
  cp <- admixr2:::.admSimulate(rxm, pars$struct, pinfo$sigma_names, eta, s0,
                               ovar, pl[[1L]], 1L, .Machine$integer.max, NULL)
  out <- vector("list", length(u))
  for (i in seq_along(u)) {
    pu <- rowSums(vapply(seq_along(sh[[i]]$delta), function(k)
      sh[[i]]$w[k] * stats::dnorm(eta[, 1L], sh[[i]]$delta[k], OM),
      numeric(n_sim)))
    w  <- pu / stats::dnorm(eta[, 1L], mu_p, s_p)
    w  <- w / sum(w)
    wm <- admixr2:::weighted_meancov_cpp(cp, w)
    mu <- as.numeric(wm$mu); V <- wm$V
    ar <- admixr2:::.admUnitResidRows(pinfo, ovar, pars$sigma_var, length(TT))
    ap <- admixr2:::.admResidApply(mu, diag(V), ar, u[[i]]$times, V)
    out[[i]] <- list(E = ap$mu, V = admixr2:::.admApplyResidTail(V, ap),
                     ess = (sum(w)^2 / sum(w^2)) / n_sim)
  }
  list(m = out, rows = nrow(cp))
}

relE <- function(a, b) max(abs(a$E - b$E) / abs(b$E))
relV <- function(a, b) max(abs(a$V - b$V)) / max(abs(b$V))

cat(sprintf("%-8s %-12s %8s %10s %10s %10s %8s\n",
            "studies", "scheme", "rows", "s/eval", "E rel", "V rel", "min ESS"))
for (S in c(2L, 3L, 5L)) {
  st <- mk_studies(S); rf <- REF[seq_len(S)]
  for (scheme in c("per-study", "shared")) {
    f <- if (scheme == "per-study") function() per_study(st, 4000L)
         else function() shared(st, 4000L)
    r <- f(); if (is.null(r)) { cat("  shared: declined\n"); next }
    for (i in 1:2) f()
    t0 <- Sys.time(); for (i in 1:5) r <- f()
    sec <- as.numeric(difftime(Sys.time(), t0, units = "secs")) / 5
    e <- max(vapply(seq_len(S), function(i) relE(r$m[[i]], rf[[i]]), numeric(1)))
    v <- max(vapply(seq_len(S), function(i) relV(r$m[[i]], rf[[i]]), numeric(1)))
    es <- if (scheme == "shared")
      min(vapply(r$m, function(x) x$ess, numeric(1))) else 1
    cat(sprintf("%-8d %-12s %8d %10.4f %10.2e %10.2e %7.1f%%\n",
                S, scheme, r$rows, sec, e, v, 100 * es))
  }
}
