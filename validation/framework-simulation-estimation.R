## ============================================================================
## Framework validation: does fitting SUMMARIES recover the truth?
## ============================================================================
##
##   Rscript validation/framework-simulation-estimation.R
##
## Everything else in validation/ checks that admixr2 computes a particular
## integral correctly. This checks the thing that actually matters: that the
## whole aggregate-data approach is a sound estimator.
##
## Per replicate:
##   1. simulate INDIVIDUAL subjects -- each with their own covariate value,
##      random effects and residual error
##   2. reduce them to exactly what a paper reports: a mean vector, a covariance
##      matrix (ML denominator) and n
##   3. throw the individuals away and fit admixr2 to the summaries alone
##
## Because n is FINITE, E and V carry genuine sampling noise, so this exercises
## the likelihood's treatment of V_obs -- not just the moment machinery. That
## matters because admixr2 scores tr(V_pred^-1 V_obs), which treats V_obs as
## Wishart-from-normal; the concentrations are not normal, so this is the
## approximation the docs describe as "point estimates recovered, reported
## uncertainty approximate". Coverage below is the measurement of that claim.
##
## Reported per parameter:
##   bias        mean(estimate) - truth
##   emp SD      sd(estimate) across replicates          <- the real uncertainty
##   mean SE     mean reported standard error            <- the claimed one
##   SE ratio    mean SE / emp SD                        <- 1.0 = calibrated
##   coverage    fraction of 95% intervals containing truth

suppressMessages({library(rxode2)})
suppressMessages(devtools::load_all(".", quiet = TRUE))

R_REP  <- as.integer(Sys.getenv("ADM_REPS", "20"))
N_SUBJ <- 100L
TIMES  <- c(0.5, 1, 2, 4, 8)
DOSE   <- 100

## ---- truth ------------------------------------------------------------------
TCL <- log(1.0); TV <- log(10); TCOV <- 0.75; OM <- 0.30; ADD <- 0.30
TRUTH <- c(tcl = TCL, tv = TV, tcov = TCOV, add.err = ADD, om = OM)

## three populations with DIFFERENT covariate distributions -- required, since a
## coefficient sharing its argument with a random effect is not identifiable
## from one population (the likelihood is exactly flat along a trade-off with
## the fixed effect and omega; verified separately to 0.00000)
POPS <- list(list(mu = -0.40, sd = 0.25),
             list(mu =  0.00, sd = 0.55),
             list(mu =  0.45, sd = 0.30))

## ---- individual-level generation, analytic, independent of the package ------
one_study <- function(pop, n) {
  a   <- rnorm(n, pop$mu, pop$sd)
  eta <- rnorm(n, 0, OM)
  cl  <- exp(TCL + TCOV * a + eta)
  v   <- exp(TV)
  Y   <- DOSE / v * exp(outer(-cl / v, TIMES)) + matrix(rnorm(n * length(TIMES), 0, ADD),
                                                        n, length(TIMES))
  list(E = colMeans(Y), V = cov.wt(Y, method = "ML")$cov,   # ML denominator, per docs
       n = n, times = TIMES, ev = rxode2::et(amt = DOSE),
       cov_dist = list(WT = list(mu = pop$mu, sd = pop$sd)))
}

mod <- function() {
  ini({tcl <- log(1.2); tv <- log(8); tcov <- 0.50
       eta.cl ~ 0.16; add.err <- 0.5})          # started away from truth
  model({cl <- exp(tcl + tcov * WT + eta.cl); v <- exp(tv)
         d/dt(centr) <- -cl / v * centr
         cp <- centr / v
         cp ~ add(add.err)})
}

cat(sprintf("simulation-estimation: %d replicates, %d subjects x %d studies\n",
            R_REP, N_SUBJ, length(POPS)))
est <- se <- matrix(NA_real_, R_REP, length(TRUTH),
                    dimnames = list(NULL, names(TRUTH)))
t0 <- Sys.time()
for (r in seq_len(R_REP)) {
  set.seed(1000 + r)
  studies <- setNames(lapply(POPS, one_study, n = N_SUBJ),
                      paste0("s", seq_along(POPS)))
  fit <- try(suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
    mod, admData(), est = "adgh",
    control = adghControl(studies = studies, n_nodes = 7L, print = 0L,
                          covMethod = "r", maxeval = 2000L)))), silent = TRUE)
  if (inherits(fit, "try-error")) { cat(sprintf("  rep %2d FAILED\n", r)); next }
  p  <- setNames(fit$parFixedDf$Estimate, rownames(fit$parFixedDf))
  s  <- setNames(fit$parFixedDf$SE,       rownames(fit$parFixedDf))
  est[r, ] <- c(p[["tcl"]], p[["tv"]], p[["tcov"]], p[["add.err"]],
                sqrt(fit$omega[1, 1]))
  se[r, ]  <- c(s[["tcl"]], s[["tv"]], s[["tcov"]], s[["add.err"]], NA_real_)
  cat(sprintf("  rep %2d  tcl %7.4f  tv %7.4f  tcov %7.4f  add %6.4f  om %6.4f\n",
              r, est[r, 1], est[r, 2], est[r, 3], est[r, 4], est[r, 5]))
}
cat(sprintf("\n%d/%d replicates completed in %.0f s\n",
            sum(stats::complete.cases(est)), R_REP,
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

ok <- stats::complete.cases(est)
cat("\n", strrep("=", 76), "\n", sep = "")
cat(sprintf("%-9s %9s %9s %9s %9s %9s %9s\n",
            "param", "truth", "mean est", "bias", "emp SD", "mean SE", "SE/SD"))
cat(strrep("-", 76), "\n")
for (j in seq_along(TRUTH)) {
  e <- est[ok, j]; s <- se[ok, j]
  cat(sprintf("%-9s %9.4f %9.4f %+9.4f %9.4f %9s %9s\n",
              names(TRUTH)[j], TRUTH[[j]], mean(e), mean(e) - TRUTH[[j]], sd(e),
              if (all(is.na(s))) "-" else sprintf("%.4f", mean(s, na.rm = TRUE)),
              if (all(is.na(s))) "-" else sprintf("%.2f", mean(s, na.rm = TRUE) / sd(e))))
}
cat("\nMonte Carlo error on each 'bias' is emp SD / sqrt(n_rep):\n")
for (j in seq_along(TRUTH))
  cat(sprintf("  %-9s +/- %.4f\n", names(TRUTH)[j], sd(est[ok, j]) / sqrt(sum(ok))))

cat("\n95% interval coverage (estimate +/- 1.96 SE, omega has no SE here):\n")
for (j in seq_along(TRUTH)) {
  s <- se[ok, j]
  if (all(is.na(s))) next
  cov <- mean(abs(est[ok, j] - TRUTH[[j]]) <= 1.96 * s, na.rm = TRUE)
  cat(sprintf("  %-9s %5.1f%%   (nominal 95%%, MC error ~%.1f%% at %d reps)\n",
              names(TRUTH)[j], 100 * cov,
              100 * sqrt(0.95 * 0.05 / sum(ok)), sum(ok)))
}
