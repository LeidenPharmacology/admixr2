## The node methods on PER-NODE aggregate data -- the development workflow's
## construction, reproduced inside admixr2.
##
## Each node carries its own (E, V), simulated at that node's covariate value,
## and the study objective is sum_k c_k * NLL_k. This is what
## `covariate workflow.Rmd` does (obsEV <- EV_given_wt(wt) inside the objective),
## what datagen(covariate=) generates, and what a publication reporting summaries
## BY COVARIATE STRATUM gives you.
##
## Run:  Rscript validation/covariate-node-data.R
suppressMessages({library(rxode2); library(nloptr)})
suppressMessages(devtools::load_all(".", quiet = TRUE))
kv <- function(k, v) cat(sprintf("%s\t%s\n", k, paste(v, collapse = ",")))

TCL <- log(1); TV <- log(10); TCOV <- 0.75; OM <- 0.30; ADD <- 0.30
TT   <- c(0.5, 1, 1.5, 2, 3, 4, 5, 6, 8)
EV   <- rxode2::et(amt = 100)
NSUB <- 100L
## ONE population -- the development workflow's setting. Per-node data supplies
## the covariate values directly, so the identifiability ridge that makes a
## single population useless for the marginal method does not apply here.
CMU <- 0.0; CSD <- 0.55

truth <- function() {
  ini({tcl <- log(1); tv <- log(10); tcov <- 0.75; add.err <- 0.3
       eta.cl ~ 0.09})
  model({cl <- exp(tcl + tcov * AAA + eta.cl); v <- exp(tv)
         d/dt(centr) <- -cl / v * centr; cp <- centr / v
         cp ~ add(add.err)})
}
fitmod <- function() {
  ini({tcl <- log(1.4); tv <- log(11); tcov <- 0.40; add.err <- 0.35
       eta.cl ~ 0.16})
  model({cl <- exp(tcl + tcov * AAA + eta.cl); v <- exp(tv)
         d/dt(centr) <- -cl / v * centr; cp <- centr / v
         cp ~ add(add.err)})
}

METHODS <- list(gh     = list(m = "gh",     n_nodes = 9L, order = 2L),
                gl     = list(m = "gl",     n_nodes = 9L, order = 2L),
                taylor = list(m = "taylor", n_nodes = 3L, order = 2L))

res <- list()
for (nm in names(METHODS)) {
  sp  <- METHODS[[nm]]
  ## per-node data: datagen() builds one sub-study per node, simulated at that
  ## node's covariate value, and attaches the node's combination coefficient.
  gen <- datagen(
    list(pop = list(model = truth, times = TT, ev = EV, n = NSUB,
                    covariate = list(AAA = list(mu = CMU, sd = CSD)))),
    control = datagenControl(n_sim = 200000L, seed = 20260813L),
    quad_method = sp$m, n_nodes = sp$n_nodes, h = CSD / 2, order = sp$order)
  kv(sprintf("nodes_%s", nm),
     sprintf("%d studies; coefs %s (sum %.6f)", length(gen),
             paste(sprintf("%+.4f", vapply(gen, function(s) s$weight %||% 1,
                                           numeric(1))), collapse = " "),
             sum(vapply(gen, function(s) s$weight %||% 1, numeric(1)))))

  t0  <- Sys.time()
  fit <- tryCatch(suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
    fitmod, admData(), est = "adgh",
    control = adghControl(studies = gen, grad = "analytical", n_nodes = 7L,
                          maxeval = 400L, print = 0L, covMethod = "none")))),
    error = function(e) e)
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(fit, "error")) {
    kv(sprintf("fit_%s_ERROR", nm), conditionMessage(fit)); next
  }
  ax <- fit$env$admExtra
  res[[nm]] <- c(tcov = unname(ax$struct[["tcov"]]),
                 tcl  = unname(ax$struct[["tcl"]]),
                 om   = unname(sqrt(ax$omega[1, 1])),
                 add  = unname(sqrt(ax$sigma_var[[1]])),
                 secs = dt)
  kv(sprintf("fit_%s", nm),
     sprintf("tcov=%.4f tcl=%.4f om=%.4f add=%.4f %.0fs",
             res[[nm]][["tcov"]], res[[nm]][["tcl"]], res[[nm]][["om"]],
             res[[nm]][["add"]], dt))
}

cat("\n")
cat(sprintf("%-8s %9s %9s %9s %9s %7s\n",
            "method", "tcov", "bias", "tcl", "om", "secs"))
cat(sprintf("%-8s %9.4f %9s %9.4f %9.4f %7s\n", "TRUTH", TCOV, "", TCL, OM, ""))
for (nm in names(res)) {
  r <- res[[nm]]
  cat(sprintf("%-8s %9.4f %+9.4f %9.4f %9.4f %7.0f\n",
              nm, r[["tcov"]], r[["tcov"]] - TCOV, r[["tcl"]], r[["om"]],
              r[["secs"]]))
}
