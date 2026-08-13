## ============================================================================
## PAGE abstract, reproduced in admixr2: two published models -> one unifying fit
## ============================================================================
##
##   Rscript validation/page-abstract-mbma.R
##
## The ADM idiom, end to end. Each study contributes data generated from ITS OWN
## published model and ITS OWN design -- different structural model, different
## dosing, different sampling times, different population -- and a single
## unifying model is then fitted across both, without individual patient data.
##
##   Issaranggoon   2-compartment, single 1 h infusion of 17.5 mg/kg,
##                  9 sampling times, n = 14, weight 20.2 +/- 9.2 kg
##   Alsultan       1-compartment, 5 infusions q6h x 2 h of 15 mg/kg,
##                  2 sampling times at 21 and 23.5 h, n = 72,
##                  weight 18.1 +/- 8.5 kg
##
##   unifying       2-compartment, allometric weight on Cl, Vp and Q
##
## Weight enters FOUR places and is the reason this exercises the general
## covariate path rather than either shortcut:
##   * Cl  -- carries a random effect  (a shift of a mu-referenced argument)
##   * Vp  -- no random effect at all  (nothing to carry a shift)
##   * Q   -- no random effect at all
##   * the DOSE itself, mg/kg, through bioavailability
##
## NOT reproduced, and both are admixr2 limitations rather than choices:
##   * IOV on the Alsultan arm (CV 9.4%) -- admixr2 has no IOV support
##   * the original's second nloptr restart, which exists to work around the
##     premature-convergence bug now filed as issue #120

suppressMessages(devtools::load_all(".", quiet = TRUE))
set.seed(1)

## ---- populations: lognormal matched to the reported mean and SD -------------
## A normal weight distribution at these CVs (46% and 47%) puts mass at negative
## body weight, where an allometric term is undefined.
ln_from <- function(mu, sd)
  list(meanlog = log(mu^2 / sqrt(sd^2 + mu^2)), sdlog = sqrt(log(1 + sd^2 / mu^2)))
WT_ISS <- ln_from(20.2, 9.2)
WT_ALS <- ln_from(18.1, 8.5)

cv2sd <- function(cv) sqrt(log(cv^2 + 1))

## ---- study 1: Issaranggoon's published model (2-cmt) ------------------------
mod_iss <- function() {
  ini({lcl <- log(0.16); lvc <- log(3.86); lvp <- log(0.19); lq <- log(0.13)
       clwt <- 0.97; vpwt <- 1.07; qwt <- 1.19
       eta.cl ~ 0.062; eta.vc ~ 0.048
       prop.err <- 0.12})
  model({cl <- exp(lcl + eta.cl) * WT^clwt
         vc <- exp(lvc + eta.vc)
         vp <- exp(lvp) * WT^vpwt
         q  <- exp(lq)  * WT^qwt
         f(centr)    <- WT                       # dose is mg/kg
         d/dt(centr) <- -(cl + q) / vc * centr + q / vp * peri
         d/dt(peri)  <-        q  / vc * centr - q / vp * peri
         cp <- centr / vc
         cp ~ prop(prop.err)})
}

## ---- study 2: Alsultan's published model (1-cmt, linear in weight) ----------
mod_als <- function() {
  ini({lcl <- log(2.99); lv <- log(9.55)
       eta.cl ~ 0.0223; eta.v ~ 0.0134        # CV 15% and 11.6% on the log scale
       prop.err <- 0.119})
  model({cl <- exp(lcl + eta.cl) * (WT / 20)
         v  <- exp(lv  + eta.v)  * (WT / 20)
         f(centr)    <- WT
         d/dt(centr) <- -cl / v * centr
         cp <- centr / v
         cp ~ prop(prop.err)})
}

## ---- the unifying analysis model -------------------------------------------
## Same structure as Issaranggoon's, deliberately: it is the richer of the two,
## and the Alsultan arm then contributes information about a model it was not
## itself built on -- which is the point of the exercise.
mod_unify <- function() {
  ini({lcl <- log(0.30); lvc <- log(2.50); lvp <- log(0.40); lq <- log(0.25)
       clwt <- 0.60; vpwt <- 0.70; qwt <- 0.80
       eta.cl ~ 0.10; eta.vc ~ 0.10
       prop.err <- 0.20})
  model({cl <- exp(lcl + eta.cl) * WT^clwt
         vc <- exp(lvc + eta.vc)
         vp <- exp(lvp) * WT^vpwt
         q  <- exp(lq)  * WT^qwt
         f(centr)    <- WT
         d/dt(centr) <- -(cl + q) / vc * centr + q / vp * peri
         d/dt(peri)  <-        q  / vc * centr - q / vp * peri
         cp <- centr / vc
         cp ~ prop(prop.err)})
}

TRUTH <- c(lcl = log(0.16), lvc = log(3.86), lvp = log(0.19), lq = log(0.13),
           clwt = 0.97, vpwt = 1.07, qwt = 1.19, prop.err = 0.12)

## ---- generate each study from its own model and design ---------------------
cat("\ngenerating study data from the published models ...\n")
t0 <- Sys.time()
gen <- datagen(
  list(
    Issaranggoon = list(
      model    = mod_iss,
      times    = c(0.5, 1, 1.25, 1.5, 2, 3, 4, 5, 6),
      ev       = rxode2::et(amt = 17.5, dur = 1),
      n        = 14L,
      cov_dist = list(WT = WT_ISS)),
    Alsultan = list(
      model    = mod_als,
      times    = c(21, 23.5),
      ev       = rxode2::et(amt = 15, ii = 6, addl = 4, dur = 2),
      n        = 72L,
      cov_dist = list(WT = WT_ALS))),
  control = datagenControl(n_sim = 40000L, seed = 12345L))
cat(sprintf("  done in %.1f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

for (nm in names(gen))
  cat(sprintf("  %-13s n=%3d  times=%s\n     E = %s\n", nm, gen[[nm]]$n,
              paste(gen[[nm]]$times, collapse = ", "),
              paste(sprintf("%.3f", gen[[nm]]$E), collapse = "  ")))

## ---- fit ONE unifying model across both --------------------------------------
cat("\nfitting the unifying model jointly (analytic gradients) ...\n")
t0 <- Sys.time()
fit <- suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
  mod_unify, admData(), est = "admc",
  control = admControl(studies = gen, n_sim = 3000L, print = 25L,
                       covMethod = "none", maxeval = 2000L))))
el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

est <- setNames(fit$parFixedDf$Estimate, rownames(fit$parFixedDf))
cat(sprintf("\n  fitted in %.1f s, objective %.3f\n", el, fit$objective))

cat("\n", strrep("=", 68), "\n  Parameter recovery\n", strrep("=", 68), "\n", sep = "")
cat(sprintf("%-10s %12s %12s %10s\n", "", "truth", "estimate", "rel err"))
for (p in names(TRUTH)) {
  e <- est[[p]]
  cat(sprintf("%-10s %12.4f %12.4f %9.1f%%\n", p, TRUTH[[p]], e,
              100 * (e - TRUTH[[p]]) / abs(TRUTH[[p]])))
}
om <- sqrt(diag(fit$omega))
cat(sprintf("%-10s %12.4f %12.4f %9.1f%%\n", "om(cl)", sqrt(0.062), om[[1]],
            100 * (om[[1]] - sqrt(0.062)) / sqrt(0.062)))
cat(sprintf("%-10s %12.4f %12.4f %9.1f%%\n", "om(vc)", sqrt(0.048), om[[2]],
            100 * (om[[2]] - sqrt(0.048)) / sqrt(0.048)))

cat("\nOn the natural scale:\n")
for (p in c("lcl", "lvc", "lvp", "lq"))
  cat(sprintf("  %-4s truth %8.4f   estimate %8.4f\n",
              sub("^l", "", p), exp(TRUTH[[p]]), exp(est[[p]])))
cat("\n")
