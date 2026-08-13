## ============================================================================
## Weight-gradient panel, to compare against the PAGE poster figure
## ============================================================================
##
##   Rscript validation/page-abstract-weightplot.R
##
## Cmax (t = 21 h) and trough (t = 23.5 h) against body weight, on the Alsultan
## design (15 mg/kg q6h, 2 h infusion, 5 doses), for three models:
##
##   Original Issaranggoon (2C)   its published parameters
##   Original Alsultan (1C)       its published parameters
##   admixr2 joint (2C)           the unifying model fitted across BOTH studies
##
## Ribbons are the 10th-90th percentile across BETWEEN-SUBJECT variability only
## (the random effects). Residual error is NOT included, so these are the spread
## of true individual exposures rather than of measured concentrations.

suppressMessages({library(rxode2); library(ggplot2)})
set.seed(1)

N_SIM  <- 600L
WT_SEQ <- seq(5, 55, by = 1)
TIMES  <- c(21, 23.5)
EV_ALS <- rxode2::et(amt = 15, ii = 6, addl = 4, dur = 2) |> rxode2::et(TIMES)

## ---- the three parameter sets ----------------------------------------------
## Issaranggoon's published 2-cmt model
P_ISS <- list(lcl = log(0.16), lvc = log(3.86), lvp = log(0.19), lq = log(0.13),
              clwt = 0.97, vpwt = 1.07, qwt = 1.19,
              om = c(cl = sqrt(0.062), vc = sqrt(0.048)))
## admixr2's joint fit across both studies (validation/page-abstract-mbma.R)
P_JNT <- list(lcl = -1.6655, lvc = 1.0203, lvp = -1.3883, lq = -3.8081,
              clwt = 0.8746, vpwt = 0.9505, qwt = 1.8737,
              om = c(cl = 0.1517, vc = 0.4656))
## Alsultan's published 1-cmt model
P_ALS <- list(lcl = log(2.99), lv = log(9.55),
              om = c(cl = sqrt(0.0223), v = sqrt(0.0134)))

m2c <- rxode2::rxode2({
  cl <- exp(lcl + eta.cl) * WT^clwt
  vc <- exp(lvc + eta.vc)
  vp <- exp(lvp) * WT^vpwt
  q  <- exp(lq)  * WT^qwt
  f(centr)    <- WT
  d/dt(centr) <- -(cl + q) / vc * centr + q / vp * peri
  d/dt(peri)  <-        q  / vc * centr - q / vp * peri
  cp <- centr / vc
})
m1c <- rxode2::rxode2({
  cl <- exp(lcl + eta.cl) * (WT / 20)
  v  <- exp(lv  + eta.v)  * (WT / 20)
  f(centr)    <- WT
  d/dt(centr) <- -cl / v * centr
  cp <- centr / v
})

## one solve per model: every (weight, subject) pair is a row
sim_grid <- function(mod, p, two_cmt) {
  z1 <- qnorm((seq_len(N_SIM) - 0.5) / N_SIM)
  z2 <- qnorm((seq(N_SIM) - 0.5 + 0.5) %% N_SIM / N_SIM + 1e-9)   # decorrelated
  g  <- expand.grid(i = seq_len(N_SIM), wt = WT_SEQ)
  pars <- if (two_cmt)
    data.frame(lcl = p$lcl, lvc = p$lvc, lvp = p$lvp, lq = p$lq,
               clwt = p$clwt, vpwt = p$vpwt, qwt = p$qwt,
               eta.cl = p$om[["cl"]] * z1[g$i], eta.vc = p$om[["vc"]] * z2[g$i],
               WT = g$wt)
  else
    data.frame(lcl = p$lcl, lv = p$lv,
               eta.cl = p$om[["cl"]] * z1[g$i], eta.v = p$om[["v"]] * z2[g$i],
               WT = g$wt)
  out <- rxode2::rxSolve(mod, params = pars, events = EV_ALS, cores = 1L,
                         nDisplayProgress = .Machine$integer.max)
  keep <- out[["time"]] %in% TIMES
  cp   <- matrix(out[["cp"]][keep], ncol = length(TIMES), byrow = TRUE)
  do.call(rbind, lapply(seq_along(TIMES), function(k) {
    m <- matrix(cp[, k], nrow = N_SIM)          # N_SIM x length(WT_SEQ)
    data.frame(wt   = WT_SEQ,
               mean = colMeans(m),
               p10  = apply(m, 2L, quantile, 0.10),
               p90  = apply(m, 2L, quantile, 0.90),
               tp   = c("Cmax (21 h)", "Trough (23.5 h)")[k])
  }))
}

cat("simulating ...\n")
df <- rbind(
  cbind(sim_grid(m2c, P_JNT, TRUE),  model = "admixr2 joint (2C)"),
  cbind(sim_grid(m2c, P_ISS, TRUE),  model = "Original Issaranggoon (2C)"),
  cbind(sim_grid(m1c, P_ALS, FALSE), model = "Original Alsultan (1C)"))
df$model <- factor(df$model, levels = c("admixr2 joint (2C)",
                                        "Original Issaranggoon (2C)",
                                        "Original Alsultan (1C)"))
df$tp <- factor(df$tp, levels = c("Cmax (21 h)", "Trough (23.5 h)"))

study_wts <- data.frame(wt = c(20.2, 18.1),
                        label = c("Issaranggoon", "Alsultan"))
cols <- c("admixr2 joint (2C)"         = "#0072B2",
          "Original Issaranggoon (2C)" = "#E69F00",
          "Original Alsultan (1C)"     = "#56B4E9")

p <- ggplot(df, aes(x = wt, colour = model, fill = model)) +
  geom_ribbon(aes(ymin = p10, ymax = p90), alpha = 0.15, colour = NA) +
  geom_line(aes(y = mean), linewidth = 0.9) +
  geom_vline(data = study_wts, aes(xintercept = wt, linetype = label),
             colour = "grey40", linewidth = 0.5) +
  facet_wrap(~tp, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = cols, name = NULL) +
  scale_fill_manual(values = cols, name = NULL) +
  scale_linetype_manual(values = c("Issaranggoon" = "dashed",
                                   "Alsultan" = "dotted"),
                        name = "Study mean WT") +
  labs(x = "Body weight (kg)", y = "Vancomycin (mg/L)",
       subtitle = "Alsultan design: 15 mg/kg q6h (2 h infusion), 5 doses; ribbon = P10-P90 of IIV") +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"),
        strip.text = element_text(size = 11, face = "bold"),
        legend.position = "bottom")

out <- file.path(Sys.getenv("ADM_FIG_DIR", "."), "page-weight-gradient.png")
ggsave(out, p, width = 10, height = 5, dpi = 150)
cat("wrote", normalizePath(out), "\n")

## numbers at the two study mean weights, for a quick read-off
cat("\nMean concentration at the study mean weights:\n")
for (w in c(18.1, 20.2)) {
  d <- df[abs(df$wt - round(w)) < 0.5, ]
  for (tp in levels(df$tp)) {
    dd <- d[d$tp == tp, ]
    cat(sprintf("  WT %4.1f  %-16s %s\n", w, tp,
                paste(sprintf("%-28s %6.2f", dd$model, dd$mean), collapse = " | ")))
  }
}
