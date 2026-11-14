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
stats_a <- admVignetteStats(dv_mat)
times_a <- stats_a$times
trial_a <- admStudy(E = stats_a$E, V = stats_a$V, n = stats_a$n,
                    times = times_a, dose = 100)

# Trial B: 120 subjects at 200 mg, four sampling times. SIMULATED here from
# the model examplomycin came from, so it is scored as 120 sampled patients
# (with standard errors), just as a digitised summary would be.
times_b <- c(0.5, 2, 6, 12)
trial_b <- admStudy(model = admVignetteModel, n = 120, dose = 200,
                    times = times_b, simulated = TRUE)

studies <- admStudies(trial_a, trial_b)
studies
#> admixr2 studies: 2  (0 published models, 1 simulated trial, 1 digitised)
#>   trial_a        data   n = 500    9 times
#>   trial_b        sim    n = 120    4 times
#> 
#> print() a single study to check its transcription.
```

`admVignetteDvMatrix()` and `admVignetteStats()` come from this
vignette’s setup file: the first reshapes `examplomycin` into one row
per subject, the second takes `E`, `V` and `n` off it.

Trial B is simulated here so the vignette has a second study to combine;
in practice its `E`, `V` and `n` would be digitised from the paper
([`vignette("aggregate-data")`](https://leidenpharmacology.github.io/admixr2/articles/aggregate-data.md)),
or the paper’s published model would be entered with
`admStudy(model, est)` — which then contributes as a model, without
standard errors
([`vignette("datagen")`](https://leidenpharmacology.github.io/admixr2/articles/datagen.md)).

## Comparing the two observed profiles

The two are on different doses, so they should *not* lie on top of each
other — that separation is the information a joint fit uses.

[`datagen()`](https://leidenpharmacology.github.io/admixr2/reference/datagen.md)
shows the `E` and `V` each study contributes:

``` r

gen <- datagen(studies)
df_obs <- rbind(
  data.frame(trial = "Trial A (100 mg, n = 500)", time = times_a,
             mean = gen$trial_a$E,
             lo   = gen$trial_a$E - sqrt(diag(gen$trial_a$V)),
             hi   = gen$trial_a$E + sqrt(diag(gen$trial_a$V))),
  data.frame(trial = "Trial B (200 mg, n = 120)", time = times_b,
             mean = gen$trial_b$E,
             lo   = gen$trial_b$E - sqrt(diag(gen$trial_b$V)),
             hi   = gen$trial_b$E + sqrt(diag(gen$trial_b$V)))
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

Pass the collection. Each study carries its own `times`, dose, `V` and
`n`, which is what lets them differ:

``` r
fit_multi <- nlmixr2(
  pk_model, admData(), est = "admc",
  control = admControl(
    studies   = studies,
    n_sim     = 5000L,
    cov_n_sim = 10000L,
    maxeval   = 300L,
    seed      = 1L
  )
)

print(fit_multi)
── nlmixr² admc ──

          OBJF       AIC       BIC Log-likelihood
admc -3351.412 -3329.412 -3257.767       1675.706

── Time (sec fit_multi$time): ──

  optimize covariance other elapsed
1   79.861     33.573     0 113.434

── Population Parameters (fit_multi$parFixed or fit_multi$parFixedDf): ──

                                  Parameter    Est.       SE  %RSE
tcl                    Log clearance (L/hr)   1.603  0.01771 1.105
tv1                  Log central volume (L)   2.328   0.1290 5.542
tv2               Log peripheral volume (L)   3.398  0.04941 1.454
tq        Log inter-compartmental CL (L/hr)   2.281  0.02543 1.115
tka     Log absorption rate constant (1/hr) 0.02909   0.1189 408.8
prop.sd      Proportional residual error SD  0.1900 0.003278 1.725
        Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
tcl        4.967 (4.797, 5.142)    32.33         NaN
tv1        10.26 (7.966, 13.21)    32.99         NaN
tv2        29.90 (27.14, 32.94)    32.08         NaN
tq         9.790 (9.314, 10.29)    33.66         NaN
tka       1.030 (0.8155, 1.300)    31.78         NaN
prop.sd 0.1900 (0.1835, 0.1964)                     
 
  Covariance Type (fit_multi$covMethod): r,s
  Some strong fixed parameter correlations exist (fit_multi$cor) :
                cor:tv1,tcl             cor:tv2,tcl              cor:tq,tcl 
                 0.341                  -0.508                   0.202  
            cor:tka,tcl         cor:prop.sd,tcl       cor:om.eta.cl,tcl 
                 0.366                  0.0553                 -0.0960  
      cor:om.eta.v1,tcl       cor:om.eta.v2,tcl        cor:om.eta.q,tcl 
                -0.321                  -0.192                   0.249  
      cor:om.eta.ka,tcl             cor:tv2,tv1              cor:tq,tv1 
                 0.304                  -0.855                   0.212  
            cor:tka,tv1         cor:prop.sd,tv1       cor:om.eta.cl,tv1 
                 0.984                 -0.0463                -0.00149  
      cor:om.eta.v1,tv1       cor:om.eta.v2,tv1        cor:om.eta.q,tv1 
                -0.682                  0.0879                   0.499  
      cor:om.eta.ka,tv1              cor:tq,tv2             cor:tka,tv2 
                 0.699                  -0.247                  -0.860  
        cor:prop.sd,tv2       cor:om.eta.cl,tv2       cor:om.eta.v1,tv2 
                0.0142                  0.0456                   0.630  
      cor:om.eta.v2,tv2        cor:om.eta.q,tv2       cor:om.eta.ka,tv2 
                0.0508                  -0.437                  -0.646  
             cor:tka,tq          cor:prop.sd,tq        cor:om.eta.cl,tq 
                 0.254                 -0.0203                 -0.0463  
       cor:om.eta.v1,tq        cor:om.eta.v2,tq         cor:om.eta.q,tq 
                -0.110                  0.0773                   0.337  
       cor:om.eta.ka,tq         cor:prop.sd,tka       cor:om.eta.cl,tka 
                0.0399                 -0.0492                -0.00303  
      cor:om.eta.v1,tka       cor:om.eta.v2,tka        cor:om.eta.q,tka 
                -0.673                  0.0831                   0.527  
      cor:om.eta.ka,tka   cor:om.eta.cl,prop.sd   cor:om.eta.v1,prop.sd 
                 0.688                 -0.0340                 -0.0252  
  cor:om.eta.v2,prop.sd    cor:om.eta.q,prop.sd   cor:om.eta.ka,prop.sd 
                -0.221                  -0.188                 -0.0105  
cor:om.eta.v1,om.eta.cl cor:om.eta.v2,om.eta.cl  cor:om.eta.q,om.eta.cl 
                0.0335                  -0.161                 0.00155  
cor:om.eta.ka,om.eta.cl cor:om.eta.v2,om.eta.v1  cor:om.eta.q,om.eta.v1 
               -0.0292                 -0.0491                  -0.259  
cor:om.eta.ka,om.eta.v1  cor:om.eta.q,om.eta.v2 cor:om.eta.ka,om.eta.v2 
                -0.897                  -0.103                  0.0632  
 cor:om.eta.ka,om.eta.q 
                 0.223  
 

  No correlations in between subject variability (BSV) matrix
  Full BSV covariance (fit_multi$omega) 
    or correlation (fit_multi$omegaR; diagonals=SDs)
  Distribution stats (mean/skewness/kurtosis/p-value) available in $shrink 
  Information about run found (fit_multi$runInfo):
   • "r,s" covariance is conditional on the ETAs; setCov() is marginal 
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
