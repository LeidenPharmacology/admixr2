## PAGE case study, four covariate constructions.
##
## Two published vancomycin models -- Issaranggoon (2-cmt, single infusion,
## Ayuthaya design) and Alsultan (1-cmt, multiple dose) -- are used to GENERATE
## aggregate data for their own designs, and one unifying 2-cmt model is then fit
## to both. Weight enters cl (which carries an eta), vp and q (which do not), and
## the mg/kg dose through f(centr).
##
##   marginal   ONE aggregate (E, V) per study, marginal over its weight
##              distribution -- what a publication reports -- scored once against
##              covariate-marginal moments.
##   gh/gl/taylor
##              PER-NODE aggregate data: one (E, V) per weight node, each scored
##              against a prediction at its own weight, combined as sum_k c_k NLL_k.
##
## The two data shapes are NOT the same dataset, and that is the point: the node
## methods need summaries by weight stratum, the marginal method needs the pooled
## summary. Both are run at their own best, so what the table compares is the
## whole approach, data requirement included.
##
## Run:  Rscript validation/page-abstract-methods.R
suppressMessages({library(rxode2); library(nloptr)})
suppressMessages(devtools::load_all(".", quiet = TRUE))
kv <- function(k, v) cat(sprintf("%s\t%s\n", k, paste(v, collapse = ",")))

ln_from <- function(mu, sd)
  list(meanlog = log(mu^2 / sqrt(sd^2 + mu^2)), sdlog = sqrt(log(1 + sd^2 / mu^2)))
WT_ISS <- ln_from(20.2, 9.2); WT_ALS <- ln_from(18.1, 8.5)
EV_ISS <- rxode2::et(amt = 17.5, dur = 1)
EV_ALS <- rxode2::et(amt = 15, ii = 6, addl = 4, dur = 2)
T_ISS  <- c(0.5, 1, 1.25, 1.5, 2, 3, 4, 5, 6); T_ALS <- c(21, 23.5)
N_ISS  <- 14L; N_ALS <- 72L

mod_iss <- function() {
  ini({lcl <- log(0.16); lvc <- log(3.86); lvp <- log(0.19); lq <- log(0.13)
       clwt <- 0.97; vpwt <- 1.07; qwt <- 1.19
       eta.cl ~ 0.062; eta.vc ~ 0.048; prop.err <- 0.12})
  model({cl <- exp(lcl + eta.cl) * WT^clwt; vc <- exp(lvc + eta.vc)
         vp <- exp(lvp) * WT^vpwt; q <- exp(lq) * WT^qwt; f(centr) <- WT
         d/dt(centr) <- -(cl + q)/vc*centr + q/vp*peri
         d/dt(peri)  <-        q /vc*centr - q/vp*peri
         cp <- centr/vc; cp ~ prop(prop.err)})
}
mod_als <- function() {
  ini({lcl <- log(2.99); lv <- log(9.55); eta.cl ~ 0.0223; eta.v ~ 0.0134
       prop.err <- 0.119})
  model({cl <- exp(lcl + eta.cl)*(WT/20); v <- exp(lv + eta.v)*(WT/20)
         f(centr) <- WT
         d/dt(centr) <- -cl/v*centr; cp <- centr/v; cp ~ prop(prop.err)})
}
## the unifying model: one 2-cmt structure for both designs
mod_fit <- function() {
  ini({lcl <- log(0.16); lvc <- log(3.86); lvp <- log(0.19); lq <- log(0.13)
       clwt <- 0.97; vpwt <- 1.07; qwt <- 1.19
       eta.cl ~ 0.062; eta.vc ~ 0.048; prop.err <- 0.12})
  model({cl <- exp(lcl + eta.cl) * WT^clwt; vc <- exp(lvc + eta.vc)
         vp <- exp(lvp) * WT^vpwt; q <- exp(lq) * WT^qwt; f(centr) <- WT
         d/dt(centr) <- -(cl + q)/vc*centr + q/vp*peri
         d/dt(peri)  <-        q /vc*centr - q/vp*peri
         cp <- centr/vc; cp ~ prop(prop.err)})
}

DG <- datagenControl(n_sim = 40000L, seed = 12345L)

## ---- marginal: one pooled (E, V) per study ---------------------------------
gen_marg <- datagen(list(
  Ayuthaya = list(model = mod_iss, times = T_ISS, ev = EV_ISS, n = N_ISS,
                  cov_dist = list(WT = WT_ISS)),
  Alsultan = list(model = mod_als, times = T_ALS, ev = EV_ALS, n = N_ALS,
                  cov_dist = list(WT = WT_ALS))), control = DG)

## ---- node methods: per-node (E, V) per study -------------------------------
gen_nodes <- function(meth, n_nodes, ord = 2L) {
  a <- datagen(list(Ayuthaya = list(model = mod_iss, times = T_ISS, ev = EV_ISS,
                                    n = N_ISS, covariate = list(WT = WT_ISS))),
               control = DG, quad_method = meth, n_nodes = n_nodes,
               h = WT_ISS$sdlog / 2, order = ord)
  b <- datagen(list(Alsultan = list(model = mod_als, times = T_ALS, ev = EV_ALS,
                                    n = N_ALS, covariate = list(WT = WT_ALS))),
               control = DG, quad_method = meth, n_nodes = n_nodes,
               h = WT_ALS$sdlog / 2, order = ord)
  c(a, b)
}

RUNS <- list(marginal = gen_marg,
             gh       = gen_nodes("gh", 9L),
             gl       = gen_nodes("gl", 9L),
             taylor   = gen_nodes("taylor", 3L))
for (nm in names(RUNS))
  kv(sprintf("data_%s", nm),
     sprintf("%d studies; coef sum %.6f", length(RUNS[[nm]]),
             sum(vapply(RUNS[[nm]], function(s) s$weight %||% 1, numeric(1)))))

## ---- fit each ---------------------------------------------------------------
PARS <- c("lcl", "lvc", "lvp", "lq", "clwt", "vpwt", "qwt")
res  <- list()
for (nm in names(RUNS)) {
  t0  <- Sys.time()
  fit <- tryCatch(suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
    mod_fit, admData(), est = "adgh",
    control = adghControl(studies = RUNS[[nm]], grad = "analytical",
                          n_nodes = 7L, maxeval = 600L, print = 0L,
                          covMethod = "none")))),
    error = function(e) e)
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(fit, "error")) {
    kv(sprintf("fit_%s_ERROR", nm), conditionMessage(fit)); next
  }
  ax <- fit$env$admExtra
  res[[nm]] <- list(
    struct = vapply(PARS, function(p) unname(ax$struct[[p]]), numeric(1)),
    om     = sqrt(diag(ax$omega)), prop = sqrt(ax$sigma_var[[1]]),
    obj    = fit$objective, secs = dt)
  kv(sprintf("fit_%s", nm), sprintf("obj=%.3f %.0fs", fit$objective, dt))
}

## ---- report ----------------------------------------------------------------
## Reference: the two generating models' own weight coefficients. Neither is
## "truth" for a joint 2-cmt fit -- the studies disagree structurally, which is
## what the case study is about -- but clwt is on the same footing in all of
## them, so it is the row to read.
cat("\n")
nms <- names(res)
cat(sprintf("%-8s", "param"));  for (n in nms) cat(sprintf(" %11s", n)); cat("\n")
for (p in PARS) {
  cat(sprintf("%-8s", p))
  for (n in nms) cat(sprintf(" %11.4f", res[[n]]$struct[[p]]))
  cat("\n")
}
for (k in seq_len(2L)) {
  cat(sprintf("%-8s", paste0("om", k)))
  for (n in nms) cat(sprintf(" %11.4f", res[[n]]$om[[k]]))
  cat("\n")
}
cat(sprintf("%-8s", "prop"))
for (n in nms) cat(sprintf(" %11.4f", res[[n]]$prop)); cat("\n")
cat(sprintf("%-8s", "secs"))
for (n in nms) cat(sprintf(" %11.0f", res[[n]]$secs)); cat("\n")
cat("\nGenerating models: Issaranggoon clwt 0.97 vpwt 1.07 qwt 1.19;",
    "\nAlsultan is linear in WT (clwt = 1 by construction, 1-cmt).\n")
cat("Objectives are NOT comparable across columns -- different functionals",
    "of\ndifferent data shapes.\n")
