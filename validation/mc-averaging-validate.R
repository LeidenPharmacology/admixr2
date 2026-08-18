## VALIDATION GATE. The harness must be the same estimator as admixr2's adgh on
## the `rows` path, not merely close to it. Same blocks in -> same estimate out.
## Nothing downstream runs until this passes.
setwd("C:/package/admixr2/.claude/worktrees/feature-covariate-quadrature")
suppressMessages(devtools::load_all(".", quiet = TRUE)); suppressMessages(library(rxode2))
source("validation/mc-averaging-harness.R")
say <- function(...) { cat(..., "\n"); utils::flush.console() }

TRUTH <- c(lcl = log(4.2), lv = log(48), bCL = 0.75, bV = 1.00,
           oCL = 0.28, oV = 0.20, sig = 0.13)
SRC <- list(mu_x = log(68/70), sd_x = 0.17, n = 120,
            times = c(0.5, 1, 2, 4, 8), dose = 500)
K <- 3L

## ---- harness blocks ---------------------------------------------------------
B <- make_blocks(TRUTH, SRC, K)
say("harness E[block1]:", paste(round(B[[1]]$E, 4), collapse = " "))

## ---- the same blocks handed to admixr2 --------------------------------------
wt_q <- function(k) { force(k)
  function(u) 70 * exp(SRC$mu_x + SRC$sd_x * stats::qnorm(((k - 1) + u)/K)) }
studies <- stats::setNames(lapply(seq_len(K), function(k)
  list(E = B[[k]]$E, V = B[[k]]$V, n = SRC$n / K, times = SRC$times,
       ev = rxode2::et(amt = SRC$dose),
       cov_dist = list(WT = list(quantile = wt_q(k))))), paste0("s", seq_len(K)))

mod <- function() {
  ini({ lcl <- log(4.2); lv <- log(48); bCL <- 0.75; bV <- 1.00
        eta.cl ~ 0.0784; eta.v ~ 0.0400; prop.err <- 0.13 })
  model({ cl <- exp(lcl + bCL * log(WT/70) + eta.cl)
          v  <- exp(lv  + bV  * log(WT/70) + eta.v)
          d/dt(centr) <- -cl/v * centr
          cp <- centr/v
          cp ~ prop(prop.err) })
}
START <- c(lcl = log(3.0), lv = log(60), bCL = 0.40, bV = 0.60,
           oCL = 0.40, oV = 0.30, sig = 0.20)
mod_s <- function() {
  ini({ lcl <- log(3.0); lv <- log(60); bCL <- 0.40; bV <- 0.60
        eta.cl ~ 0.16; eta.v ~ 0.09; prop.err <- 0.20 })
  model({ cl <- exp(lcl + bCL * log(WT/70) + eta.cl)
          v  <- exp(lv  + bV  * log(WT/70) + eta.v)
          d/dt(centr) <- -cl/v * centr
          cp <- centr/v
          cp ~ prop(prop.err) })
}
f <- suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
  mod_s, admData(), est = "adgh",
  control = adghControl(studies = studies, n_nodes = 5L, grad = "analytical",
                        print = 0L, covMethod = "none"))))
fx <- nlmixr2est::fixef(f); om <- f$omega
pkg <- c(lcl = fx[["lcl"]], lv = fx[["lv"]], bCL = fx[["bCL"]], bV = fx[["bV"]],
         oCL = sqrt(om[["eta.cl","eta.cl"]]), oV = sqrt(om[["eta.v","eta.v"]]),
         sig = fx[["prop.err"]])
har <- fit_form(B, PAR_NM, START)

say("\n           ", paste(sprintf("%8s", PAR_NM), collapse = ""))
say("truth      ", paste(sprintf("%8.4f", TRUTH[PAR_NM]), collapse = ""))
say("admixr2    ", paste(sprintf("%8.4f", pkg[PAR_NM]),   collapse = ""))
say("harness    ", paste(sprintf("%8.4f", har[PAR_NM]),   collapse = ""))
say("|diff|     ", paste(sprintf("%8.1e", abs(pkg[PAR_NM] - har[PAR_NM])), collapse = ""))
say(sprintf("\nmax |pkg - harness| = %.2e     max |harness - truth| = %.2e",
            max(abs(pkg[PAR_NM] - har[PAR_NM])), max(abs(har[PAR_NM] - TRUTH[PAR_NM]))))
say(if (max(abs(pkg[PAR_NM] - har[PAR_NM])) < 5e-3) "GATE: PASS" else "GATE: FAIL")
