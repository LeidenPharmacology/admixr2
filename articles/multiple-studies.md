# Multiple studies

## Why multiple studies?

Passing several studies to
[`admControl()`](https://leidenpharmacology.github.io/admixr2/reference/admControl.md)
fits them simultaneously, minimising the sum of per-study NLLs under one
shared set of population parameters. This is **meta-analysis**, the core
use case: summary statistics from several trials that may differ in
dose, size or schedule, and one population model consistent with all of
them. Each study’s `E`, `V` and `n` can come from a digitised figure
([`vignette("aggregate-data")`](https://leidenpharmacology.github.io/admixr2/articles/aggregate-data.md))
or from its own published model
([`vignette("datagen")`](https://leidenpharmacology.github.io/admixr2/articles/datagen.md)).

## Splitting examplomycin into two cohorts

Split the 500 examplomycin subjects into two cohorts of 250, with
separate aggregate statistics for each:

``` r

library(admixr2)
library(rxode2)
library(nlmixr2)
library(ggplot2)

dv_mat <- admVignetteDvMatrix()       # 500 subjects x 9 times
times  <- as.numeric(colnames(dv_mat))

# Alternate subjects into two equal cohorts, then take E, V and n for each
cohort1 <- admVignetteStats(dv_mat, seq(1, nrow(dv_mat), by = 2))
cohort2 <- admVignetteStats(dv_mat, seq(2, nrow(dv_mat), by = 2))

E1 <- cohort1$E; V1 <- cohort1$V; n1 <- cohort1$n
E2 <- cohort2$E; V2 <- cohort2$V; n2 <- cohort2$n
```

Both helpers are defined in this vignette’s setup file: the first
reshapes `examplomycin` into one row per subject, the second takes `E`,
`V` and `n` off a set of its rows.

## Comparing observed profiles across cohorts

Check the raw summary statistics first: both cohorts come from the same
population here, so they should be comparable.

``` r

df_obs <- rbind(
  data.frame(cohort = "Cohort 1", time = times,
             mean = E1,
             lo   = E1 - sqrt(diag(V1)),
             hi   = E1 + sqrt(diag(V1))),
  data.frame(cohort = "Cohort 2", time = times,
             mean = E2,
             lo   = E2 - sqrt(diag(V2)),
             hi   = E2 + sqrt(diag(V2)))
)

ggplot(df_obs, aes(x = time, y = mean, colour = cohort, fill = cohort)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_x_log10(breaks = times, labels = times) +
  scale_colour_manual(values = c("Cohort 1" = "#0072B2", "Cohort 2" = "#D55E00")) +
  scale_fill_manual(  values = c("Cohort 1" = "#0072B2", "Cohort 2" = "#D55E00")) +
  labs(title    = "Observed mean ± 1 SD by cohort",
       x        = "Time (h, log scale)",
       y        = "Concentration",
       colour   = NULL, fill = NULL) +
  theme_bw()
```

![Observed mean ± 1 SD for each cohort on a log time
axis.](multiple-studies_files/figure-html/obs-compare-1.png)

Observed mean ± 1 SD for each cohort on a log time axis.

## Model definition

The two-compartment model from [Getting
started](https://leidenpharmacology.github.io/admixr2/articles/admixr2.md),
supplied by this vignette’s setup file and fitted to both cohorts at
once:

``` r

pk_model <- admVignetteModel
```

## Fitting with two studies

Pass both cohorts as a named list. Each entry may independently specify
`times`, `ev`, `V`, `n`, and `method`:

``` r
fit_multi <- nlmixr2(
  pk_model, admData(), est = "admc",
  control = admControl(
    studies = list(
      cohort1 = list(E = E1, V = V1, n = n1,
                     times = times, ev = rxode2::et(amt = 100)),
      cohort2 = list(E = E2, V = V2, n = n2,
                     times = times, ev = rxode2::et(amt = 100))
    ),
    n_sim     = 5000L,
    cov_n_sim = 10000L,
    maxeval   = 300L,
    seed      = 1L
  )
)

print(fit_multi)
── nlmixr² admc ──

          OBJF       AIC       BIC Log-likelihood
admc -3690.262 -3668.262 -3597.732       1845.131

── Time (sec fit_multi$time): ──

  optimize covariance other elapsed
1    45.84     20.181     0  66.021

── Population Parameters (fit_multi$parFixed or fit_multi$parFixedDf): ──

                                  Parameter    Est.       SE  %RSE
tcl                    Log clearance (L/hr)   1.602  0.02000 1.248
tv1                  Log central volume (L)   2.328   0.1320 5.670
tv2               Log peripheral volume (L)   3.397  0.05338 1.571
tq        Log inter-compartmental CL (L/hr)   2.276  0.02729 1.199
tka     Log absorption rate constant (1/hr) 0.02979   0.1230 412.9
prop.sd      Proportional residual error SD  0.1895 0.003282 1.732
        Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
tcl        4.963 (4.772, 5.161)    32.62         NaN
tv1        10.26 (7.920, 13.29)    33.15         NaN
tv2        29.89 (26.92, 33.18)    31.81         NaN
tq         9.738 (9.231, 10.27)    33.68         NaN
tka       1.030 (0.8096, 1.311)    31.87         NaN
prop.sd 0.1895 (0.1831, 0.1960)                     
 
  Covariance Type (fit_multi$covMethod): r,s
  No correlations in between subject variability (BSV) matrix
  Full BSV covariance (fit_multi$omega) 
    or correlation (fit_multi$omegaR; diagonals=SDs)
  Distribution stats (mean/skewness/kurtosis/p-value) available in $shrink 
  Censoring (fit_multi$censInformation): No censoring
  Minimization message (fit_multi$message):  
    NLOPT_XTOL_REACHED: Optimization stopped because xtol_rel or xtol_abs (above) was reached. 
```

## Per-study diagnostic plots

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) produces
separate panels per study, named `mean_<study>` and `cov_<study>`:

``` r

plots <- plot(fit_multi, which = "mean")
```

![Mean diagnostics for both cohorts (one panel per
study).](multiple-studies_files/figure-html/diag-1.png)

Mean diagnostics for both cohorts (one panel per study).

![Mean diagnostics for both cohorts (one panel per
study).](multiple-studies_files/figure-html/diag-2.png)

Mean diagnostics for both cohorts (one panel per study).

``` r

names(plots)
#>  [1] "mean_cohort1"           "mean_cohort1_obs"       "mean_cohort1_pred"     
#>  [4] "mean_cohort1_resid"     "mean_cohort1_std_resid" "mean_cohort2"          
#>  [7] "mean_cohort2_obs"       "mean_cohort2_pred"      "mean_cohort2_resid"    
#> [10] "mean_cohort2_std_resid"
```

Access individual panels to compare studies side by side:

``` r

plots$mean_cohort1
plots$mean_cohort2

# Combine with patchwork if installed
if (requireNamespace("patchwork", quietly = TRUE)) {
  patchwork::wrap_plots(plots, ncol = 1)
}
```

## Different doses and schedules

Studies may differ in any aspect — a typical development programme:

``` r

fit_program <- nlmixr2(
  pk_model, admData(), est = "admc",
  control = admControl(
    studies = list(
      phase1_50mg  = list(E = E_50,  V = V_50,  n = 30L,
                          times = c(1, 2, 4, 8),
                          ev    = rxode2::et(amt = 50)),
      phase2_100mg = list(E = E_100, V = V_100, n = 120L,
                          times = c(0.5, 1, 2, 4, 8, 12),
                          ev    = rxode2::et(amt = 100)),
      phase2_200mg = list(E = E_200, V = V_200, n = 115L,
                          times = c(0.5, 1, 2, 4, 8, 12),
                          ev    = rxode2::et(amt = 200))
    ),
    n_sim   = 5000L,
    maxeval = 1000L,
    seed    = 1L
  )
)
```

A study with a diagonal `V` (or a plain vector of variances) gets
`method = "var"`, skipping the O(n_t³) Cholesky solve there is no
off-diagonal structure to justify.

## See also

- [From a published figure to E, V and
  n](https://leidenpharmacology.github.io/admixr2/articles/aggregate-data.md)
  — prepare each study’s `E`, `V` and `n`
- [Several observed
  compartments](https://leidenpharmacology.github.io/admixr2/articles/multi-compartment.md)
  — several outputs per study
- [Estimator
  comparison](https://leidenpharmacology.github.io/admixr2/articles/estimator-comparison.md)
