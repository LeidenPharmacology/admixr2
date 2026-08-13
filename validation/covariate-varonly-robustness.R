## How robust is covariate MBMA when studies report only SDs?
##
## page-abstract-varonly.R showed that on ONE two-study problem the variance-only
## likelihood cannot separate residual error from IIV. That was a single dataset
## with no sampling noise and a misspecified joint model, so it cannot say
## whether the COVARIATE COEFFICIENT itself survives -- which is the question
## that matters for MBMA, since sigma and omega are nuisance parameters there.
##
## This is a simulation-estimation study over replicates WITH sampling noise:
## each replicate generates its (E, V) from a finite sample of n subjects
## (n_sim = n), so the summaries carry the noise a real study's would. Fitted
## twice, cov and var, on the SAME data.
##
## Three populations with differing covariate means AND spreads, because one
## population does not identify the coefficient at all (the likelihood is exactly
## flat along tcl' = tcl + (b-b')*mu_a, omega'^2 = omega^2 + (b^2-b'^2)*sd_a^2).
##
## Run:  Rscript validation/covariate-varonly-robustness.R
suppressMessages({library(rxode2); library(nloptr)})
suppressMessages(devtools::load_all(".", quiet = TRUE))
kv <- function(k, v) cat(sprintf("%s\t%s\n", k, paste(v, collapse = ",")))

TCL <- log(1); TV <- log(10); TCOV <- 0.75; OM <- 0.30; ADD <- 0.30
EV   <- rxode2::et(amt = 100)
POPS <- list(p1 = c(mu = -0.45, sd = 0.30),
             p2 = c(mu =  0.00, sd = 0.55),
             p3 = c(mu =  0.50, sd = 0.35))
RICH   <- c(0.5, 1, 1.5, 2, 3, 4, 5, 6, 8)
SPARSE <- c(2, 6)
NREP   <- as.integer(Sys.getenv("NREP", "20"))

truth <- function() {
  ini({tcl <- log(1); tv <- log(10); tcov <- 0.75; add.err <- 0.3
       eta.cl ~ 0.09})
  model({cl <- exp(tcl + tcov * AAA + eta.cl); v <- exp(tv)
         d/dt(centr) <- -cl / v * centr; cp <- centr / v
         cp ~ add(add.err)})
}
fitmod <- function() {
  ini({tcl <- log(1.3); tv <- log(11); tcov <- 0.45; add.err <- 0.35
       eta.cl ~ 0.16})
  model({cl <- exp(tcl + tcov * AAA + eta.cl); v <- exp(tv)
         d/dt(centr) <- -cl / v * centr; cp <- centr / v
         cp ~ add(add.err)})
}

one_rep <- function(rep, times, n) {
  ## finite-sample summaries: n_sim = n, so (E, V) carry a real study's noise
  gen <- datagen(
    setNames(lapply(POPS, function(p)
      list(model = truth, times = times, ev = EV, n = n,
           cov_dist = list(AAA = list(mu = p[["mu"]], sd = p[["sd"]])))),
      names(POPS)),
    control = datagenControl(n_sim = n, seed = 1000L + rep, sampling = "rnorm"))
  gen_var <- lapply(gen, function(s) { s$V <- diag(as.matrix(s$V)); s })

  grab <- function(studies) {
    f <- tryCatch(suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
      fitmod, admData(), est = "adgh",
      control = adghControl(studies = studies, grad = "analytical",
                            n_nodes = 7L, maxeval = 2000L, print = 0L,
                            covMethod = "none")))), error = function(e) NULL)
    if (is.null(f)) return(rep(NA_real_, 5L))
    ax <- f$env$admExtra
    c(tcl = unname(ax$struct[["tcl"]]), tv = unname(ax$struct[["tv"]]),
      tcov = unname(ax$struct[["tcov"]]),
      add = unname(sqrt(ax$sigma_var[[1]])), om = unname(sqrt(ax$omega[1, 1])))
  }
  rbind(cov = grab(gen), var = grab(gen_var))
}

DESIGNS <- list(
  list(lbl = "rich  n=100",  times = RICH,   n = 100L),
  list(lbl = "rich  n=20",   times = RICH,   n =  20L),
  list(lbl = "sparse n=100", times = SPARSE, n = 100L))

TRUE_V <- c(tcl = TCL, tv = TV, tcov = TCOV, add = ADD, om = OM)
for (d in DESIGNS) {
  res <- list(cov = NULL, var = NULL)
  for (r in seq_len(NREP)) {
    o <- tryCatch(one_rep(r, d$times, d$n), error = function(e) NULL)
    if (is.null(o)) next
    res$cov <- rbind(res$cov, o["cov", ]); res$var <- rbind(res$var, o["var", ])
  }
  cat(sprintf("\n== %s, %d replicates ==\n", d$lbl, NROW(res$cov)))
  cat(sprintf("%-6s %-8s %9s %9s %9s %9s %7s\n",
              "meth", "param", "truth", "mean", "bias", "emp SD", "n_ok"))
  for (m in c("cov", "var")) {
    M <- res[[m]]
    for (p in names(TRUE_V)) {
      x <- M[, p]; x <- x[is.finite(x)]
      cat(sprintf("%-6s %-8s %9.4f %9.4f %+9.4f %9.4f %7d\n",
                  m, p, TRUE_V[[p]], mean(x), mean(x) - TRUE_V[[p]],
                  stats::sd(x), length(x)))
    }
  }
  ## the MBMA question: is the covariate coefficient still usable?
  for (m in c("cov", "var")) {
    x <- res[[m]][, "tcov"]; x <- x[is.finite(x)]
    kv(sprintf("tcov_%s_%s", m, gsub("[^A-Za-z0-9]", "", d$lbl)),
       sprintf("bias %+.4f  SD %.4f  RMSE %.4f",
               mean(x) - TCOV, stats::sd(x),
               sqrt(mean((x - TCOV)^2))))
  }
}
