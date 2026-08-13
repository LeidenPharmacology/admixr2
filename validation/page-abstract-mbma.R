## ============================================================================
## PAGE abstract in admixr2: two published models -> one unifying fit
## ============================================================================
##
##   Rscript validation/page-abstract-mbma.R
##
## The ADM idiom end to end. Each study contributes data generated from ITS OWN
## published model and ITS OWN design; a single unifying model is then fitted
## across both, without individual patient data.
##
##   Issaranggoon   2-compartment, single 1 h infusion of 17.5 mg/kg,
##                  9 sampling times, n = 14, weight 20.2 +/- 9.2 kg
##   Alsultan       1-compartment, 5 infusions q6h x 2 h of 15 mg/kg,
##                  2 sampling times (21, 23.5 h), n = 72, weight 18.1 +/- 8.5 kg
##
## Weight enters FOUR places, which is why this needs the general covariate path:
## Cl (carries an eta), Vp and Q (no eta at all), and the mg/kg DOSE itself.
##
## HOW TO READ THE RESULT. The two studies come from DIFFERENT published models,
## so the joint fit is a CONSENSUS and its optimum is not at either study's
## parameters. Parameter "recovery" against Issaranggoon's values is therefore
## not a correctness criterion -- an earlier version of this script printed one
## and it was meaningless. The fit is judged by how well it describes BOTH
## studies' data (the RMS table and Panel A below).
##
## Identifiability note: a covariate effect is identified only by BETWEEN-STUDY
## variation in the covariate distribution. Within one population the likelihood
## is exactly flat along tcl' = tcl + (tcov-tcov')*mu_a,
## omega'^2 = omega^2 + (tcov^2-tcov'^2)*sd_a^2 -- verified to 0.00000. Here the
## two studies differ in both mean weight (20.2 vs 18.1) and spread (9.2 vs 8.5),
## so the allometric exponents are identified, if weakly.
##
## NOT reproduced, both admixr2 limitations: IOV on the Alsultan arm, and the
## original's steady-state closed form (explicit multiple dosing is used).

suppressMessages({library(rxode2); library(ggplot2)})
suppressMessages(devtools::load_all(".", quiet = TRUE))
set.seed(1)

FIG <- Sys.getenv("ADM_FIG_DIR", ".")
ln_from <- function(mu, sd)
  list(meanlog = log(mu^2 / sqrt(sd^2 + mu^2)), sdlog = sqrt(log(1 + sd^2 / mu^2)))
WT_ISS <- ln_from(20.2, 9.2); WT_ALS <- ln_from(18.1, 8.5)

EV_ISS <- rxode2::et(amt = 17.5, dur = 1)
EV_ALS <- rxode2::et(amt = 15, ii = 6, addl = 4, dur = 2)
T_ISS  <- c(0.5, 1, 1.25, 1.5, 2, 3, 4, 5, 6)
T_ALS  <- c(21, 23.5)

## ---- the two published models ----------------------------------------------
mod_iss <- function() {
  ini({lcl <- log(0.16); lvc <- log(3.86); lvp <- log(0.19); lq <- log(0.13)
       clwt <- 0.97; vpwt <- 1.07; qwt <- 1.19
       eta.cl ~ 0.062; eta.vc ~ 0.048; prop.err <- 0.12})
  model({cl <- exp(lcl + eta.cl) * WT^clwt; vc <- exp(lvc + eta.vc)
         vp <- exp(lvp) * WT^vpwt; q <- exp(lq) * WT^qwt
         f(centr) <- WT
         d/dt(centr) <- -(cl + q)/vc*centr + q/vp*peri
         d/dt(peri)  <-        q /vc*centr - q/vp*peri
         cp <- centr/vc; cp ~ prop(prop.err)})
}
mod_als <- function() {
  ini({lcl <- log(2.99); lv <- log(9.55)
       eta.cl ~ 0.0223; eta.v ~ 0.0134; prop.err <- 0.119})
  model({cl <- exp(lcl + eta.cl) * (WT/20); v <- exp(lv + eta.v) * (WT/20)
         f(centr) <- WT
         d/dt(centr) <- -cl/v*centr; cp <- centr/v; cp ~ prop(prop.err)})
}
## the unifying model: Issaranggoon's structure, started well away from it
mod_unify <- function() {
  ini({lcl <- log(0.30); lvc <- log(2.50); lvp <- log(0.40); lq <- log(0.25)
       clwt <- 0.60; vpwt <- 0.70; qwt <- 0.80
       eta.cl ~ 0.10; eta.vc ~ 0.10; prop.err <- 0.20})
  model({cl <- exp(lcl + eta.cl) * WT^clwt; vc <- exp(lvc + eta.vc)
         vp <- exp(lvp) * WT^vpwt; q <- exp(lq) * WT^qwt
         f(centr) <- WT
         d/dt(centr) <- -(cl + q)/vc*centr + q/vp*peri
         d/dt(peri)  <-        q /vc*centr - q/vp*peri
         cp <- centr/vc; cp ~ prop(prop.err)})
}

## ---- generate each study from its own model and design ----------------------
cat("generating study data from the published models ...\n")
gen <- datagen(list(
  Ayuthaya = list(model = mod_iss, times = T_ISS, ev = EV_ISS, n = 14L,
                  cov_dist = list(WT = WT_ISS)),
  Alsultan = list(model = mod_als, times = T_ALS, ev = EV_ALS, n = 72L,
                  cov_dist = list(WT = WT_ALS))),
  control = datagenControl(n_sim = 40000L, seed = 12345L))

## ---- fit ONE unifying model across both ------------------------------------
## Estimator chosen by benchmark, not habit. On this covariate structure, against
## exact nested quadrature, adgh's product grid at 7 nodes gives E 1.1e-06 /
## V 2.3e-06 per evaluation where admc needs 16000 draws for E 4.0e-05 /
## V 3.1e-03 -- and V is what this objective mostly compares. adgh is also
## deterministic, so the optimiser sees no Monte Carlo noise.
EST <- Sys.getenv("ADM_EST", "adgh")
stopifnot(EST %in% c("adgh", "admc"))
cat(sprintf("fitting the unifying model jointly (est = %s) ...\n", EST))
t0  <- Sys.time()
## braces required: at top level, `if (x) expr` <newline> `else` is a parse error
ctl <- if (EST == "adgh") {
  adghControl(studies = gen, n_nodes = 7L, print = 50L,
              covMethod = "none", maxeval = 3000L)
} else {
  admControl(studies = gen, n_sim = 4000L, print = 50L,
             covMethod = "none", maxeval = 3000L)
}
fit <- suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
  mod_unify, admData(), est = EST, control = ctl)))
cat(sprintf("  %.0f s, objective %.3f\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs")), fit$objective))
est <- setNames(fit$parFixedDf$Estimate, rownames(fit$parFixedDf))
om  <- sqrt(diag(fit$omega))
cat("\nFitted unifying model (a consensus across two different published models,\n")
cat("so these are NOT expected to equal either study's parameters):\n")
print(round(c(est, om.cl = om[[1]], om.vc = om[[2]]), 4))

## ---- three parameter sets, the MBMA one taken FROM the fit ------------------
P_ISS <- list(two = TRUE, prop = 0.12,
              th = c(lcl = log(0.16), lvc = log(3.86), lvp = log(0.19),
                     lq = log(0.13), clwt = 0.97, vpwt = 1.07, qwt = 1.19),
              om = c(sqrt(0.062), sqrt(0.048)))
P_ALS <- list(two = FALSE, prop = 0.119,
              th = c(lcl = log(2.99), lv = log(9.55)),
              om = c(sqrt(0.0223), sqrt(0.0134)))
P_MBMA <- list(two = TRUE, prop = unname(est[["prop.err"]]),
               th = c(lcl = est[["lcl"]], lvc = est[["lvc"]], lvp = est[["lvp"]],
                      lq = est[["lq"]], clwt = est[["clwt"]],
                      vpwt = est[["vpwt"]], qwt = est[["qwt"]]),
               om = c(om[[1]], om[[2]]))

m2c <- rxode2::rxode2({
  cl <- exp(lcl + eta.cl) * WT^clwt; vc <- exp(lvc + eta.vc)
  vp <- exp(lvp) * WT^vpwt; q <- exp(lq) * WT^qwt
  f(centr) <- WT
  d/dt(centr) <- -(cl + q)/vc*centr + q/vp*peri
  d/dt(peri)  <-        q /vc*centr - q/vp*peri
  cp <- centr/vc
})
m1c <- rxode2::rxode2({
  cl <- exp(lcl + eta.cl) * (WT/20); v <- exp(lv + eta.v) * (WT/20)
  f(centr) <- WT
  d/dt(centr) <- -cl/v*centr; cp <- centr/v
})

N_SIM <- 800L
draws <- function(P, wt, n = N_SIM) {
  u  <- (seq_len(n) - 0.5) / n
  z1 <- qnorm(u); z2 <- qnorm((u + 0.5) %% 1)
  if (P$two)
    data.frame(t(P$th), eta.cl = P$om[1]*z1, eta.vc = P$om[2]*z2, WT = wt)
  else
    data.frame(t(P$th), eta.cl = P$om[1]*z1, eta.v = P$om[2]*z2, WT = wt)
}
profile_of <- function(P, ev, times, wtd) {
  wt <- qlnorm(rev((seq_len(N_SIM) - 0.5)/N_SIM), wtd$meanlog, wtd$sdlog)
  out <- rxode2::rxSolve(if (P$two) m2c else m1c, params = draws(P, wt),
                         events = ev |> rxode2::et(times), cores = 1L,
                         nDisplayProgress = .Machine$integer.max)
  cpm <- matrix(out[["cp"]][out[["time"]] %in% times], nrow = N_SIM, byrow = TRUE)
  mu  <- colMeans(cpm)
  data.frame(time = times, mean = mu,
             sd = sqrt(pmax(apply(cpm, 2L, var) + P$prop^2 * mu^2, 0)))
}

## ---- Panel A: obs vs pred per study ----------------------------------------
cat("\nPanel A ...\n")
DES <- list(Ayuthaya = list(ev = EV_ISS, tg = seq(0.05, 6, by = 0.05),
                            wt = WT_ISS, own = P_ISS, cross = P_ALS),
            Alsultan = list(ev = EV_ALS, tg = seq(18, 24, by = 0.05),
                            wt = WT_ALS, own = P_ALS, cross = P_ISS))
df_pred <- do.call(rbind, lapply(names(DES), function(sn) {
  d <- DES[[sn]]
  rbind(cbind(profile_of(d$own,   d$ev, d$tg, d$wt), model = "Original", study = sn),
        cbind(profile_of(d$cross, d$ev, d$tg, d$wt), model = "Cross",    study = sn),
        cbind(profile_of(P_MBMA,  d$ev, d$tg, d$wt), model = "MBMA",     study = sn))
}))
df_pred$model <- factor(df_pred$model, levels = c("Original", "Cross", "MBMA"))
df_obs <- do.call(rbind, lapply(names(gen), function(sn) {
  g <- gen[[sn]]
  data.frame(study = sn, time = g$times, mean = as.numeric(g$E),
             sd = sqrt(pmax(diag(g$V), 0)))
}))
df_obs$lo <- df_obs$mean - df_obs$sd; df_obs$hi <- df_obs$mean + df_obs$sd

cols <- c(Original = "#E69F00", Cross = "#56B4E9", MBMA = "#1B7837")
pA <- ggplot() +
  geom_ribbon(data = df_pred, aes(time, ymin = mean - sd, ymax = mean + sd,
                                  fill = model), alpha = 0.15, colour = NA) +
  geom_line(data = df_pred, aes(time, mean, colour = model), linewidth = 0.9) +
  geom_pointrange(data = df_obs, aes(time, mean, ymin = lo, ymax = hi),
                  colour = "black", size = 0.3, linewidth = 0.6) +
  facet_wrap(~study, scales = "free", nrow = 1) +
  scale_colour_manual(values = cols, name = NULL) +
  scale_fill_manual(values = cols, name = NULL) +
  labs(x = "Time (h)", y = "Vancomycin (mg/L)",
       subtitle = "Marginal over each study's weight distribution; residual included") +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"),
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "bottom",
        plot.subtitle = element_text(size = 8, colour = "grey40"))
ggsave(file.path(FIG, "page-panelA.png"), pA, width = 11, height = 5, dpi = 150)

cat("\nRMS relative deviation from the observed means:\n")
for (sn in names(DES)) {
  ob <- df_obs[df_obs$study == sn, ]
  for (m in levels(df_pred$model)) {
    pr <- df_pred[df_pred$study == sn & df_pred$model == m, ]
    v  <- approx(pr$time, pr$mean, xout = ob$time)$y
    cat(sprintf("  %-9s %-9s %6.1f%%\n", sn, m, 100*sqrt(mean((v/ob$mean - 1)^2))))
  }
}

## ---- weight-gradient panel --------------------------------------------------
cat("\nWeight-gradient panel ...\n")
WT_SEQ <- seq(5, 55, by = 1); TP <- c(21, 23.5)
wt_stats <- function(P) do.call(rbind, lapply(WT_SEQ, function(w) {
  out <- rxode2::rxSolve(if (P$two) m2c else m1c, params = draws(P, rep(w, N_SIM)),
                         events = EV_ALS |> rxode2::et(TP), cores = 1L,
                         nDisplayProgress = .Machine$integer.max)
  cpm <- matrix(out[["cp"]][out[["time"]] %in% TP], nrow = N_SIM, byrow = TRUE)
  data.frame(wt = w, tp = c("Cmax (21 h)", "Trough (23.5 h)"),
             mean = colMeans(cpm),
             p10 = apply(cpm, 2L, quantile, 0.10),
             p90 = apply(cpm, 2L, quantile, 0.90))
}))
df_wt <- rbind(cbind(wt_stats(P_MBMA), model = "admixr2 joint (2C)"),
               cbind(wt_stats(P_ISS),  model = "Original Issaranggoon (2C)"),
               cbind(wt_stats(P_ALS),  model = "Original Alsultan (1C)"))
df_wt$model <- factor(df_wt$model, levels = c("admixr2 joint (2C)",
  "Original Issaranggoon (2C)", "Original Alsultan (1C)"))
wcols <- c("admixr2 joint (2C)" = "#0072B2",
           "Original Issaranggoon (2C)" = "#E69F00",
           "Original Alsultan (1C)" = "#56B4E9")
pW <- ggplot(df_wt, aes(wt, colour = model, fill = model)) +
  geom_ribbon(aes(ymin = p10, ymax = p90), alpha = 0.15, colour = NA) +
  geom_line(aes(y = mean), linewidth = 0.9) +
  geom_vline(data = data.frame(wt = c(20.2, 18.1),
                               label = c("Issaranggoon", "Alsultan")),
             aes(xintercept = wt, linetype = label), colour = "grey40") +
  facet_wrap(~tp, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = wcols, name = NULL) +
  scale_fill_manual(values = wcols, name = NULL) +
  scale_linetype_manual(values = c(Issaranggoon = "dashed", Alsultan = "dotted"),
                        name = "Study mean WT") +
  labs(x = "Body weight (kg)", y = "Vancomycin (mg/L)",
       subtitle = "Alsultan design: 15 mg/kg q6h (2 h infusion); ribbon = P10-P90 of IIV") +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"),
        strip.text = element_text(size = 11, face = "bold"),
        legend.position = "bottom")
ggsave(file.path(FIG, "page-weight-gradient.png"), pW, width = 10, height = 5, dpi = 150)
cat("wrote both figures to", normalizePath(FIG), "\n")
