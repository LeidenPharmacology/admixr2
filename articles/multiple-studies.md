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

## Two trials that differ

The point of a meta-analysis is that the studies are *not*
interchangeable, so this one combines a trial we hold individual records
for with a trial we have only summary statistics from — at a different
dose, a sparser schedule and a different size.

``` r

library(admixr2)
library(rxode2)
library(nlmixr2)
library(ggplot2)

# Trial A: 500 subjects at 100 mg, richly sampled. Individual records, reduced
# to the E, V and n a publication would print.
dv_mat  <- admVignetteDvMatrix()      # 500 subjects x 9 times
trial_a <- admVignetteStats(dv_mat)
times_a <- trial_a$times

# Trial B: published as summary statistics only -- 120 subjects at 200 mg, four
# sampling times. datagen() turns its model into the aggregate data it implies,
# which is what makes a published study a direct input; see vignette("datagen").
times_b <- c(0.5, 2, 6, 12)
trial_b <- datagen(
  studies = list(b = list(times = times_b,
                          ev    = rxode2::et(amt = 200),
                          n     = 120L)),
  model   = admVignetteModel,
  control = datagenControl(n_sim = 10000L, seed = 1L)
)$b
```

`admVignetteDvMatrix()` and `admVignetteStats()` come from this
vignette’s setup file: the first reshapes `examplomycin` into one row
per subject, the second takes `E`, `V` and `n` off it.

Trial B is simulated here so the vignette has a second study to combine;
in practice its `E`, `V` and `n` would be digitised from the paper
([`vignette("aggregate-data")`](https://leidenpharmacology.github.io/admixr2/articles/aggregate-data.md))
or derived from the model the paper published.

## Comparing the two observed profiles

The two are on different doses, so they should *not* lie on top of each
other — that separation is the information a joint fit uses.

``` r

df_obs <- rbind(
  data.frame(trial = "Trial A (100 mg, n = 500)", time = times_a,
             mean = trial_a$E,
             lo   = trial_a$E - sqrt(diag(trial_a$V)),
             hi   = trial_a$E + sqrt(diag(trial_a$V))),
  data.frame(trial = "Trial B (200 mg, n = 120)", time = times_b,
             mean = trial_b$E,
             lo   = trial_b$E - sqrt(diag(trial_b$V)),
             hi   = trial_b$E + sqrt(diag(trial_b$V)))
)

ggplot(df_obs, aes(x = time, y = mean, colour = trial, fill = trial)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_x_log10(breaks = sort(unique(c(times_a, times_b)))) +
  scale_colour_manual(values = c("Trial A (100 mg, n = 500)" = "#0072B2",
                                 "Trial B (200 mg, n = 120)" = "#D55E00")) +
  scale_fill_manual(  values = c("Trial A (100 mg, n = 500)" = "#0072B2",
                                 "Trial B (200 mg, n = 120)" = "#D55E00")) +
  labs(title  = "Observed mean ± 1 SD by trial",
       x      = "Time (h, log scale)",
       y      = "Concentration",
       colour = NULL, fill = NULL) +
  theme_bw()
```

![Observed mean ± 1 SD for each trial on a log time axis. Trial B is at
twice the dose and four sampling
times.](multiple-studies_files/figure-html/obs-compare-1.png)

Observed mean ± 1 SD for each trial on a log time axis. Trial B is at
twice the dose and four sampling times.

## Model definition

The two-compartment model from [Getting
started](https://leidenpharmacology.github.io/admixr2/articles/admixr2.md),
supplied by this vignette’s setup file and fitted to both trials at
once:

``` r

pk_model <- admVignetteModel
```

## Fitting with two studies

Pass both trials as a named list. Each entry carries its own `times`,
`ev`, `V`, `n` and `method`, which is what lets them differ:

``` r
fit_multi <- nlmixr2(
  pk_model, admData(), est = "admc",
  control = admControl(
    studies = list(
      trial_a = list(E = trial_a$E, V = trial_a$V, n = trial_a$n,
                     times = times_a, ev = rxode2::et(amt = 100)),
      trial_b = list(E = trial_b$E, V = trial_b$V, n = trial_b$n,
                     times = times_b, ev = rxode2::et(amt = 200))
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
admc -3351.844 -3329.844 -3258.199       1675.922

── Time (sec fit_multi$time): ──

  optimize covariance other elapsed
1    43.15     25.036     0  68.186

── Population Parameters (fit_multi$parFixed or fit_multi$parFixedDf): ──

                                  Parameter    Est.       SE  %RSE
tcl                    Log clearance (L/hr)   1.603  0.01771 1.105
tv1                  Log central volume (L)   2.328   0.1293 5.553
tv2               Log peripheral volume (L)   3.398  0.04948 1.456
tq        Log inter-compartmental CL (L/hr)   2.281  0.02542 1.114
tka     Log absorption rate constant (1/hr) 0.02919   0.1192 408.2
prop.sd      Proportional residual error SD  0.1900 0.003277 1.725
        Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
tcl        4.967 (4.797, 5.142)    32.32         NaN
tv1        10.26 (7.963, 13.22)    32.98         NaN
tv2        29.90 (27.14, 32.95)    32.07         NaN
tq         9.789 (9.314, 10.29)    33.66         NaN
tka       1.030 (0.8151, 1.301)    31.79         NaN
prop.sd 0.1900 (0.1835, 0.1964)                     
 
  Covariance Type (fit_multi$covMethod): r,s
  Some strong fixed parameter correlations exist (fit_multi$cor) :
                cor:tv1,tcl             cor:tv2,tcl              cor:tq,tcl 
                 0.342                  -0.508                   0.202  
            cor:tka,tcl         cor:prop.sd,tcl       cor:om.eta.cl,tcl 
                 0.367                  0.0553                 -0.0961  
      cor:om.eta.v1,tcl       cor:om.eta.v2,tcl        cor:om.eta.q,tcl 
                -0.322                  -0.191                   0.249  
      cor:om.eta.ka,tcl             cor:tv2,tv1              cor:tq,tv1 
                 0.305                  -0.856                   0.211  
            cor:tka,tv1         cor:prop.sd,tv1       cor:om.eta.cl,tv1 
                 0.984                 -0.0463                -0.00150  
      cor:om.eta.v1,tv1       cor:om.eta.v2,tv1        cor:om.eta.q,tv1 
                -0.683                  0.0883                   0.500  
      cor:om.eta.ka,tv1              cor:tq,tv2             cor:tka,tv2 
                 0.701                  -0.247                  -0.861  
        cor:prop.sd,tv2       cor:om.eta.cl,tv2       cor:om.eta.v1,tv2 
                0.0142                  0.0456                   0.632  
      cor:om.eta.v2,tv2        cor:om.eta.q,tv2       cor:om.eta.ka,tv2 
                0.0501                  -0.438                  -0.648  
             cor:tka,tq          cor:prop.sd,tq        cor:om.eta.cl,tq 
                 0.253                 -0.0202                 -0.0463  
       cor:om.eta.v1,tq        cor:om.eta.v2,tq         cor:om.eta.q,tq 
                -0.111                  0.0774                   0.337  
       cor:om.eta.ka,tq         cor:prop.sd,tka       cor:om.eta.cl,tka 
                0.0402                 -0.0492                -0.00304  
      cor:om.eta.v1,tka       cor:om.eta.v2,tka        cor:om.eta.q,tka 
                -0.675                  0.0836                   0.528  
      cor:om.eta.ka,tka   cor:om.eta.cl,prop.sd   cor:om.eta.v1,prop.sd 
                 0.689                 -0.0340                 -0.0251  
  cor:om.eta.v2,prop.sd    cor:om.eta.q,prop.sd   cor:om.eta.ka,prop.sd 
                -0.221                  -0.188                 -0.0106  
cor:om.eta.v1,om.eta.cl cor:om.eta.v2,om.eta.cl  cor:om.eta.q,om.eta.cl 
                0.0335                  -0.161                 0.00154  
cor:om.eta.ka,om.eta.cl cor:om.eta.v2,om.eta.v1  cor:om.eta.q,om.eta.v1 
               -0.0291                 -0.0495                  -0.260  
cor:om.eta.ka,om.eta.v1  cor:om.eta.q,om.eta.v2 cor:om.eta.ka,om.eta.v2 
                -0.897                  -0.103                  0.0636  
 cor:om.eta.ka,om.eta.q 
                 0.225  
 

  No correlations in between subject variability (BSV) matrix
  Full BSV covariance (fit_multi$omega) 
    or correlation (fit_multi$omegaR; diagonals=SDs)
  Distribution stats (mean/skewness/kurtosis/p-value) available in $shrink 
  Censoring (fit_multi$censInformation): No censoring
  Minimization message (fit_multi$message):  
    NLOPT_XTOL_REACHED: Optimization stopped because xtol_rel or xtol_abs (above) was reached. 
```

Nothing above assumes the two share a dose or a schedule. The objective
is the sum of per-study negative log-likelihoods under one set of
population parameters, and each study is predicted at its own dosing and
its own times.

## Per-study diagnostic plots

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) produces
separate panels per study, named `mean_<study>` and `cov_<study>`:

``` r

plots <- plot(fit_multi, which = "mean")
```

![Mean diagnostics for both trials (one panel per
study).](multiple-studies_files/figure-html/diag-1.png)

Mean diagnostics for both trials (one panel per study).

![Mean diagnostics for both trials (one panel per
study).](multiple-studies_files/figure-html/diag-2.png)

Mean diagnostics for both trials (one panel per study).

``` r

names(plots)
#>  [1] "mean_trial_a"           "mean_trial_a_obs"       "mean_trial_a_pred"     
#>  [4] "mean_trial_a_resid"     "mean_trial_a_std_resid" "mean_trial_b"          
#>  [7] "mean_trial_b_obs"       "mean_trial_b_pred"      "mean_trial_b_resid"    
#> [10] "mean_trial_b_std_resid"
```

Access individual panels to compare studies side by side:

``` r

plots$mean_trial_a
plots$mean_trial_b

# Combine with patchwork if installed
if (requireNamespace("patchwork", quietly = TRUE)) {
  patchwork::wrap_plots(plots, ncol = 1)
}
```

## Scaling to a programme

A third and fourth study are more entries in the same list, so a
development programme is written the same way the two above are — one
entry per trial, each with its own dose, schedule and size.

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
