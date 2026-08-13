## ============================================================================
## u-quantile covariate marginalisation, end to end through admc
## ============================================================================
##
##   Rscript validation/covariate-uq-endtoend.R
##
## TWO models with genuinely different covariate effects, each fitted to three
## study populations with DIFFERENT covariate means and spreads:
##
##   model A  allometric   cl <- exp(tcl + eta.cl) * (WT/70)^tcov
##   model B  exponential  cl <- exp(tcl + eta.cl) * exp(tcov * (WT - 70))
##
## Neither is a bare `theta * WT` product, so neither is collapsible and neither
## is recognised by rxode2's muRefCovariateDataFrame. Both go through the
## u-quantile path, which measures Delta(a) from the model itself.
##
## Reference aggregate moments are built by EXACT nested 2-D Gauss-Hermite over
## (WT, eta) on an analytic 1-cmt solution, in plain R. No part of the reference
## touches the code under test.

suppressMessages(devtools::load_all(".", quiet = TRUE))

TCL <- log(1.0); TV <- log(10); OM <- 0.30; ADD <- 0.30
DOSE <- 100; TIMES <- c(0.5, 1, 2, 4, 8)

conc <- function(cl) DOSE / exp(TV) * exp(outer(-cl / exp(TV), TIMES))

gh <- function(k) {
  i <- seq_len(k - 1L)
  J <- matrix(0, k, k); J[cbind(i, i + 1L)] <- sqrt(i); J[cbind(i + 1L, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE); o <- order(e$values)
  list(x = e$values[o], w = e$vectors[1L, o]^2)
}
Q <- gh(40L)

## Exact reference: nested quadrature over (WT, eta), WT lognormal.
ref_study <- function(p, cl_of) {
  m1 <- 0; M2 <- 0
  for (ia in seq_along(Q$x)) {
    wt <- exp(p$meanlog + p$sdlog * Q$x[ia])
    for (ib in seq_along(Q$x)) {
      w <- Q$w[ia] * Q$w[ib]
      Y <- as.numeric(conc(cl_of(wt, OM * Q$x[ib])))
      m1 <- m1 + w * Y; M2 <- M2 + w * outer(Y, Y)
    }
  }
  V <- M2 - outer(m1, m1); diag(V) <- diag(V) + ADD^2
  list(E = m1, V = V, n = p$n, times = TIMES, ev = rxode2::et(amt = DOSE))
}

## Three populations: different covariate MEANS and different SPREADS.
POPS <- list(
  light = list(meanlog = log(55), sdlog = 0.12, n = 220L),
  mid   = list(meanlog = log(72), sdlog = 0.28, n = 320L),
  heavy = list(meanlog = log(95), sdlog = 0.18, n = 260L))

MODELS <- list(
  allometric = list(
    tcov = 0.75,
    cl   = function(tcov) function(wt, eta) exp(TCL + eta) * (wt / 70)^tcov,
    fn   = function() {
      ini({tcl <- log(1.0); tv <- log(10); tcov <- 0.40
           eta.cl ~ 0.09; add.err <- 0.3})
      model({cl <- exp(tcl + eta.cl) * (WT / 70)^tcov
             v  <- exp(tv)
             d/dt(centr) <- -cl / v * centr
             cp <- centr / v
             cp ~ add(add.err)})
    }),
  exponential = list(
    tcov = 0.012,
    cl   = function(tcov) function(wt, eta) exp(TCL + eta) * exp(tcov * (wt - 70)),
    fn   = function() {
      ini({tcl <- log(1.0); tv <- log(10); tcov <- 0.005
           eta.cl ~ 0.09; add.err <- 0.3})
      model({cl <- exp(tcl + eta.cl) * exp(tcov * (WT - 70))
             v  <- exp(tv)
             d/dt(centr) <- -cl / v * centr
             cp <- centr / v
             cp ~ add(add.err)})
    })
)

## Budget note: the covariate paths refuse gradients, so this is BOBYQA. The
## objective is shallow near the optimum (tcov 0.76 -> 0.80 costs only ~1.5 -2LL
## units) and carries MC noise from n_sim, so a derivative-free search needs a
## real evaluation budget or it stops anywhere in the basin. Profiling the
## objective directly showed its minimum at the truth; an earlier run with
## maxeval = 400 and n_sim = 4000 reported tcov = 0.80 purely from stopping early.
fit_one <- function(mod, studies, n_sim = 12000L) {
  ctl <- admControl(studies = studies, grad = "none", n_sim = n_sim,
                    print = 0L, covMethod = "none", maxeval = 3000L)
  t0 <- Sys.time()
  f <- try(suppressWarnings(suppressMessages(
        nlmixr2est::nlmixr2(mod, admData(), est = "admc", control = ctl))),
       silent = TRUE)
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(f, "try-error"))
    return(list(err = gsub("[[:space:]]+", " ",
                           conditionMessage(attr(f, "condition")))))
  p <- setNames(f$parFixedDf$Estimate, rownames(f$parFixedDf))
  list(tcl = unname(p["tcl"]), tv = unname(p["tv"]), tcov = unname(p["tcov"]),
       om = sqrt(as.numeric(f$omega[1, 1])), secs = el)
}

cat("\n")
for (nm in names(MODELS)) {
  M   <- MODELS[[nm]]
  clf <- M$cl(M$tcov)
  ref <- lapply(POPS, ref_study, cl_of = clf)

  ## u-quantile: reference covariate value + the distribution the subjects span
  st_uq <- Map(function(r, p)
    c(r, list(cov      = list(WT = exp(p$meanlog)),
              cov_dist = list(WT = list(meanlog = p$meanlog, sdlog = p$sdlog)))),
    ref, POPS)
  ## baseline: covariate at its reference, spread ignored
  st_mo <- Map(function(r, p) c(r, list(cov = list(WT = exp(p$meanlog)))), ref, POPS)

  cat(strrep("=", 78), "\n", sep = "")
  cat(sprintf("model %s   (covariate effect NOT a bare theta*WT product)\n", nm))
  cat(strrep("=", 78), "\n", sep = "")
  cat(sprintf("populations: %s\n", paste(vapply(names(POPS), function(k)
    sprintf("%s median WT %.0f, sdlog %.2f", k, exp(POPS[[k]]$meanlog),
            POPS[[k]]$sdlog), character(1)), collapse = " | ")))
  cat(sprintf("%-16s %9s %9s %9s %9s %8s\n",
              "", "tcl", "tv", "tcov", "omega", "secs"))
  cat(sprintf("%-16s %9.4f %9.4f %9.4f %9.4f\n",
              "TRUTH", TCL, TV, M$tcov, OM))
  for (lbl in c("u-quantile", "mean-only")) {
    r <- fit_one(M$fn, if (lbl == "u-quantile") st_uq else st_mo)
    if (!is.null(r$err)) cat(sprintf("%-16s ERROR: %s\n", lbl, substr(r$err, 1, 60)))
    else cat(sprintf("%-16s %9.4f %9.4f %9.4f %9.4f %8.1f\n",
                     lbl, r$tcl, r$tv, r$tcov, r$om, r$secs))
  }
  cat("\n")
}
