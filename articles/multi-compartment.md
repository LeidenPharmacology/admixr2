# Several observed compartments (plasma and brain)

## The problem

You are developing a CNS drug and you want one number: **how much of it
reaches the brain?** You don’t have patient-level data. What you *do*
have is a published paper with two figures — a **plasma** and a
**brain** concentration–time curve, each drawn as a mean with error
bars. You digitise them into a mean and an SD at every sampling time.

That aggregate summary is exactly what admixr2 is built for. It fits a
population PK model *directly* to means-and-covariances, so a digitised
figure becomes a fittable dataset — no individual records required. This
vignette builds the analysis in two steps:

1.  **Plasma only** — the ordinary single-output workflow, to set the
    scene.
2.  **Plasma + brain** — add the brain as a second observed output and
    read the brain-penetration ratio straight off the fit.

``` r

library(admixr2)
library(rxode2)
library(nlmixr2)
library(ggplot2)

ev <- rxode2::et(amt = 100, cmt = "central")   # single 100-unit dose, shared throughout
```

## Step 1 — Plasma only

Start where every PK analysis starts: the plasma curve. In admixr2 a
**study** is the digitised summary bundled with its design — the mean
vector `E`, its variance `V` (here `SD^2`, read as a diagonal
covariance, exactly what error bars give you), the sample size `n`, the
sampling `times`, and the dosing `ev`.

``` r

plasma_times <- c(0.5, 1, 2, 4, 8, 12)
plasma_mean  <- c(8.793, 7.812, 6.850, 5.597, 3.985, 2.910)
plasma_sd    <- c(1.151, 0.911, 0.765, 0.649, 0.770, 0.828)

plasma_study <- list(
  E = plasma_mean, V = plasma_sd^2, n = 60L,
  times = plasma_times, ev = ev
)
```

The model is an ordinary two-compartment PK model with one observed
output, `cp`. Fitting is a single call to `nlmixr2()` with
`est = "adgh"`, admixr2’s Gauss–Hermite estimator:

``` r
pk_plasma <- function() {
  ini({
    tcl <- log(1);  tv1 <- log(10);  tq <- log(3);  tv2 <- log(8)
    prop.cp <- 0.05
    eta.cl ~ 0.09
    eta.v1 ~ 0.04
  })
  model({
    cl <- exp(tcl + eta.cl); v1 <- exp(tv1 + eta.v1)
    q  <- exp(tq);           v2 <- exp(tv2)
    d/dt(central) <- -(cl/v1)*central - (q/v1)*central + (q/v2)*periph
    d/dt(periph)  <-  (q/v1)*central - (q/v2)*periph
    cp <- central / v1
    cp ~ prop(prop.cp)
  })
}

fit_plasma <- nlmixr2(pk_plasma, admData(), est = "adgh",
                      control = adghControl(studies = list(trial = plasma_study)))
fit_plasma
 [1m──  [34mnlmix [39m [31mr² [39m  [33madgh [39m ── [22m

         OBJF      AIC      BIC Log-likelihood
adgh 229.6289 243.6289 270.8316      -114.8144

 [1m── Time (sec  [33mfit_plasma [39m [34m$time [39m): ── [22m

        optimize covariance elapsed other
elapsed    3.705      0.356   4.061  0.08

 [1m── Population Parameters ( [33mfit_plasma [39m [34m$parFixed [39m or  [33mfit_plasma [39m [34m$parFixedDf [39m): ── [22m

             [1mEst. [0m       [1m [1mSE [0m [0m   [1m%RSE [0m  [1mBack-transformed(95%CI) [0m  [1mBSV(CV%) [0m  [1mShrink(SD)% [0m
 [1m [1mtcl [0m [0m      0.03319 0.03098 93.36   1.034 (0.9728, 1.098)     27.9            
 [1m [1mtv1 [0m [0m        2.292 0.04476 1.953     9.893 (9.062, 10.8)     15.1            
 [1m [1mtq [0m [0m        0.8808  0.3692 41.91     2.413 (1.17, 4.975)                     
 [1m [1mtv2 [0m [0m        0.867   0.153 17.64     2.38 (1.763, 3.212)                     
 [1m [1mprop.cp [0m [0m 0.004104                              0.004104                     
 
  Covariance Type ( [33mfit_plasma [39m [1m [34m$covMethod [39m [22m):  [1mr [22m
  No correlations in between subject variability (BSV) matrix
  Full BSV covariance ( [33mfit_plasma [39m [1m [34m$omega [39m [22m) 
    or correlation ( [33mfit_plasma [39m [1m [34m$omegaR [39m [22m; diagonals=SDs)
  Distribution stats (mean/skewness/kurtosis/p-value) available in  [1m [34m$shrink [39m [22m 
  Censoring ( [33mfit_plasma [39m [1m [34m$censInformation [39m [22m): No censoring
  Minimization message ( [33mfit_plasma [39m [1m [34m$message [39m [22m):  
    NLOPT_XTOL_REACHED: Optimization stopped because xtol_rel or xtol_abs (above) was reached. 
```

This is a perfectly good plasma model — but look at what it *cannot*
answer. Its `periph` compartment is a mathematical distribution store:
we never measured it, and nothing connects it to the brain. To quantify
brain exposure we need brain data **and** a model with a real brain
compartment.

## Step 2 — Add the brain

Here are the brain concentrations digitised from the same paper:

``` r

brain_times <- c(1, 2, 4, 8, 12)
brain_mean  <- c(3.004, 3.394, 3.018, 2.157, 1.551)
brain_sd    <- c(0.353, 0.349, 0.309, 0.369, 0.405)
```

Now swap the anonymous peripheral compartment for a **mechanistic brain
compartment**. Drug moves plasma → brain with influx clearance `qin` and
back brain → plasma with efflux clearance `qout`. The steady-state
brain:plasma ratio is the quantity we want:

``` math
K_{p,uu} = \frac{q_{in}}{q_{out}}
```

The model now has **two** observed outputs — plasma `cp` and brain `cb`
— so it carries a residual-error term for each. (`vb`, the brain volume,
is a fixed physiological constant, not an estimated parameter.)

``` r

pk_cns <- function() {
  ini({
    tcl  <- log(1);  tv1  <- log(10)
    tqin <- log(3);  tqout <- log(6)
    prop.cp <- 0.05      # plasma residual (proportional)
    add.cb  <- 0.02      # brain residual (additive)
    eta.cl ~ 0.09
    eta.v1 ~ 0.04
  })
  model({
    cl <- exp(tcl + eta.cl); v1 <- exp(tv1 + eta.v1)
    qin <- exp(tqin);        qout <- exp(tqout)
    vb  <- 5
    d/dt(central) <- -(cl/v1)*central - (qin/v1)*central + (qout/vb)*brain
    d/dt(brain)   <-  (qin/v1)*central - (qout/vb)*brain
    cp <- central / v1        # plasma concentration
    cb <- brain   / vb        # brain concentration
    cp ~ prop(prop.cp)
    cb ~ add(add.cb)
  })
}
```

Two observed outputs means the study needs two summaries. Instead of a
single `E`/`V`, give it an **`observations` list** — one named entry per
observed compartment, each pairing a model output with its own `times`,
`E` and `V`:

``` r

cns_study <- list(
  n = 60L, ev = ev,
  observations = list(
    plasma = list(output = "cp", times = plasma_times, E = plasma_mean, V = plasma_sd^2),
    brain  = list(output = "cb", times = brain_times,  E = brain_mean,  V = brain_sd^2)
  )
)
```

The only other change is telling
[`admData()`](https://leidenpharmacology.github.io/admixr2/reference/admData.md)
which outputs to expect. Then the fit call is identical to Step 1:

``` r
fit_cns <- nlmixr2(pk_cns, admData(c("cp", "cb")), est = "adgh",
                   control = adghControl(studies = list(lit = cns_study)))
fit_cns
 [1m──  [34mnlmix [39m [31mr² [39m  [33madgh [39m ── [22m

          OBJF       AIC       BIC Log-likelihood
adgh -88.58817 -72.58817 -36.65025       44.29408

 [1m── Time (sec  [33mfit_cns [39m [34m$time [39m): ── [22m

  optimize covariance elapsed
1     3.72      0.797   4.517

 [1m── Population Parameters ( [33mfit_cns [39m [34m$parFixed [39m or  [33mfit_cns [39m [34m$parFixedDf [39m): ── [22m

            [1mEst. [0m       [1m [1mSE [0m [0m    [1m%RSE [0m  [1mBack-transformed(95%CI) [0m  [1mBSV(CV%) [0m  [1mShrink(SD)% [0m
 [1m [1mtcl [0m [0m     0.04145 0.01875  45.22    1.042 (1.005, 1.081)     27.0            
 [1m [1mtv1 [0m [0m       2.269 0.01028 0.4529    9.672 (9.479, 9.868)     13.9            
 [1m [1mtqin [0m [0m      1.083 0.03604  3.328    2.953 (2.752, 3.169)                     
 [1m [1mtqout [0m [0m     1.777 0.04141   2.33    5.915 (5.454, 6.415)                     
 [1m [1mprop.cp [0m [0m 0.04853                                0.04853                     
 [1m [1madd.cb [0m [0m  0.01999                                0.01999                     
 
  Covariance Type ( [33mfit_cns [39m [1m [34m$covMethod [39m [22m):  [1mr [22m
  No correlations in between subject variability (BSV) matrix
  Full BSV covariance ( [33mfit_cns [39m [1m [34m$omega [39m [22m) 
    or correlation ( [33mfit_cns [39m [1m [34m$omegaR [39m [22m; diagonals=SDs)
  Distribution stats (mean/skewness/kurtosis/p-value) available in  [1m [34m$shrink [39m [22m 
  Censoring ( [33mfit_cns [39m [1m [34m$censInformation [39m [22m): No censoring
  Minimization message ( [33mfit_cns [39m [1m [34m$message [39m [22m):  
    NLOPT_XTOL_REACHED: Optimization stopped because xtol_rel or xtol_abs (above) was reached. 
```

### Model against data

A single fit now describes both compartments. The figure overlays the
observed summaries (points, mean ± SD) with the fitted population
prediction — the mean curve and the ±SD band implied by the estimated
between-subject variability (1000 simulated subjects, residual error
excluded). Plasma and brain are distinguished by colour.

``` r

# 1. Pull the fitted estimates from the fit (standard nlmixr2 accessors).
theta <- fit_cns$theta   # fixed effects (tcl, tv1, tqin, tqout, residuals)
omega <- fit_cns$omega   # between-subject covariance (rows/cols already named)

# 2. Simulate the fitted population. We solve a plain, residual-free copy of the
#    model so the band shows between-subject variability alone -- and because a
#    multi-endpoint fit can only be re-solved with per-endpoint DVID/CMT tags,
#    whereas this gives cp and cb directly.
sim_model <- rxode2::rxode2({
  cl <- exp(tcl + eta.cl); v1 <- exp(tv1 + eta.v1)
  qin <- exp(tqin); qout <- exp(tqout); vb <- 5
  d/dt(central) <- -(cl/v1)*central - (qin/v1)*central + (qout/vb)*brain
  d/dt(brain)   <-  (qin/v1)*central - (qout/vb)*brain
  cp <- central / v1
  cb <- brain / vb
})
grid <- seq(0.25, 13, by = 0.25)
sim  <- rxode2::rxSolve(sim_model,
  params = theta[c("tcl", "tv1", "tqin", "tqout")],
  omega  = omega, nSub = 1000L,
  events = ev |> rxode2::et(grid), returnType = "data.frame")

# 3. Summarise model and data as mean +/- SD per time, per compartment.
band <- function(value)
  data.frame(time = sort(unique(sim$time)),
             mean = tapply(value, sim$time, mean),
             sd   = tapply(value, sim$time, sd))

pred <- rbind(cbind(compartment = "Plasma", band(sim$cp)),
              cbind(compartment = "Brain",  band(sim$cb)))
obs  <- rbind(
  data.frame(compartment = "Plasma", time = plasma_times, mean = plasma_mean, sd = plasma_sd),
  data.frame(compartment = "Brain",  time = brain_times,  mean = brain_mean,  sd = brain_sd))

# 4. Model (line + band) over data (points + error bars), coloured by compartment.
pal <- c(Plasma = "#0072B2", Brain = "#D55E00")   # Okabe-Ito, colour-blind safe

ggplot() +
  geom_ribbon(data = pred,
              aes(time, ymin = mean - sd, ymax = mean + sd, fill = compartment),
              alpha = 0.15) +
  geom_line(data = pred,
            aes(time, mean, colour = compartment), linewidth = 1) +
  geom_errorbar(data = obs,
                aes(time, ymin = mean - sd, ymax = mean + sd, colour = compartment),
                width = 0.4, linewidth = 0.6, show.legend = FALSE) +
  geom_point(data = obs,
             aes(time, mean, fill = compartment),
             shape = 21, colour = "white", size = 3, stroke = 0.7) +
  scale_colour_manual(values = pal, breaks = c("Plasma", "Brain")) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_x_continuous(breaks = seq(0, 12, 2),
                     expand = expansion(mult = c(0.01, 0.03))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  coord_cartesian(ylim = c(0, NA), clip = "off") +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.4))) +
  labs(x = "Time (h)", y = "Concentration (mg/L)",
       title = "Fitted population prediction vs. observed data",
       subtitle = "Line + band: model (mean ± SD)   ·   Points + bars: data (mean ± SD)",
       colour = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.justification = "left",
    legend.margin = margin(b = 2),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(colour = "grey40", size = 9, margin = margin(b = 9)),
    axis.title = element_text(colour = "grey25"),
    axis.title.x = element_text(margin = margin(t = 6)),
    axis.title.y = element_text(margin = margin(r = 6)),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.4),
    plot.margin = margin(12, 16, 10, 12)
  )
```

![](multi-compartment_files/figure-html/brain-plot-1.png)

The prediction tracks both compartments. Now the payoff — read the brain
penetration directly off the estimates:

``` r

Kp_uu <- exp(theta[["tqin"]]) / exp(theta[["tqout"]])
round(Kp_uu, 2)
#> [1] 0.5
```

`Kp,uu ≈ 0.5`: at steady state the brain sees about half the plasma
concentration. This is the whole reason the brain data was needed — with
plasma alone, `qin` and `qout` are not separately identifiable and
`Kp,uu` cannot be estimated. The brain measurements resolve it.

### Built-in diagnostics

The overlay above was assembled by hand for a custom figure, but you
don’t have to: calling
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) on the fit
draws observed-vs-predicted panels directly — one per observed output.

``` r

plot(fit_cns, which = "mean")
```

![](multi-compartment_files/figure-html/plot-fit-1.png)![](multi-compartment_files/figure-html/plot-fit-2.png)

## Notes

- **Estimators.** `adgh` (used here), `adfo` and `admc` all support
  several observed outputs, each with its own analytical / sensitivity
  gradient. `adirmc` does not — use one of the other three. See the
  [estimator
  comparison](https://leidenpharmacology.github.io/admixr2/articles/estimator-comparison.md)
  vignette to choose.
- **Structural vs. observed compartments.** How many compartments the
  ODEs contain is irrelevant to admixr2; what matters is how many
  outputs you *observe* and fit — one in Step 1, two in Step 2.
- **Hard-coded constants** such as `vb <- 5` keep their value; not every
  physiological constant has to be an estimated parameter.
- **Same-subject data.** Here plasma and brain came from separate
  figures, so they are treated as independent likelihood blocks. If they
  had instead been measured in the *same* subjects, you can supply the
  plasma–brain cross-covariance for a joint fit (a per-output-pair
  `cross` list, or a study-level joint `V`); see
  [`?admControl`](https://leidenpharmacology.github.io/admixr2/reference/admControl.md).
  \`\`\`
