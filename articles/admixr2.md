# Getting started with admixr2

## What is aggregate-data modelling?

`admixr2` fits pharmacometric PK/PD models to **summary-level data**
instead of individual records — a **meta-analysis** framework for
population PK/PD. The input is **digitised aggregate data** (means,
error bars and covariances), **previously published models**, or both,
and the result is one population model with interpretable fixed, random
and covariate effects.

For each study you supply:

- **E** — observed mean vector (one entry per observation time)
- **V** — observed covariance matrix (or variance vector)
- **n** — sample size
- **times** — observation time points
- **ev** — dosing event table

Each study is written down with
[`admStudy()`](https://leidenpharmacology.github.io/admixr2/reference/admStudy.md),
which checks it as it is built and prints it so it can be read back
against the source. The estimators match E and V against their
model-predicted counterparts and return a standard nlmixr2 fit object.
Turning a published figure into `E`, `V` and `n` is covered in
[`vignette("aggregate-data")`](https://leidenpharmacology.github.io/admixr2/articles/aggregate-data.md);
letting a study carry its own published model instead, in
[`vignette("datagen")`](https://leidenpharmacology.github.io/admixr2/articles/datagen.md).

Four estimators are available:

| Estimator | `est =` | Control function | Approach |
|----|----|----|----|
| First-Order | `"adfo"` | [`adfoControl()`](https://leidenpharmacology.github.io/admixr2/reference/adfoControl.md) | First-order Taylor expansion at η = 0; one rxSolve per NLL eval; fastest |
| Monte Carlo | `"admc"` | [`admControl()`](https://leidenpharmacology.github.io/admixr2/reference/admControl.md) | Sample average over η; asymptotically exact |
| Gauss-Hermite | `"adgh"` | [`adghControl()`](https://leidenpharmacology.github.io/admixr2/reference/adghControl.md) | Deterministic quadrature over η; noise-free, unbiased at any IIV |
| Iterative Reweighting MC | `"adirmc"` | [`adirmcControl()`](https://leidenpharmacology.github.io/admixr2/reference/adirmcControl.md) | Proposals fixed per phase; inner loop needs no new rxSolve calls |

Start with `adfo` for screening and initial estimates. `adgh` is
deterministic and most efficient with few random effects, so it is used
for the fit below; `admc` is the workhorse for standard PK models;
`adirmc` suits expensive ODE solves, high-dimensional IIV or poor
starting values. See
[`vignette("estimator-comparison")`](https://leidenpharmacology.github.io/admixr2/articles/estimator-comparison.md).

## The examplomycin dataset

`examplomycin` ships with admixr2: 500 simulated subjects from a
two-compartment PK model with first-order absorption (100 mg oral dose,
sampled at 0.1, 0.25, 0.5, 1, 2, 3, 5, 8, and 12 h). True parameters: CL
= 5 L/h, V1 = 10 L, V2 = 30 L, Q = 10 L/h, ka = 1 h⁻¹; IIV = 0.3 (SD on
log scale) for all parameters; proportional error SD = 0.2.

``` r

library(admixr2)
library(rxode2)
library(nlmixr2)

data("examplomycin")
head(examplomycin[examplomycin$EVID == 0, c("ID", "TIME", "DV")], 9)
#>     ID  TIME    DV
#> 2  460  0.10 0.752
#> 3  460  0.25 1.932
#> 4  460  0.50 3.694
#> 5  460  1.00 3.479
#> 6  460  2.00 4.003
#> 7  460  3.00 3.825
#> 8  460  5.00 1.756
#> 9  460  8.00 1.155
#> 10 460 12.00 0.742
```

## Computing aggregate statistics

Reshape individual records into a subjects × times matrix, then compute
E and V:

``` r

dv_mat <- admVignetteDvMatrix()          # 500 subjects x 9 times

E     <- colMeans(dv_mat)
V     <- cov.wt(dv_mat, method = "ML")$cov
n     <- nrow(dv_mat)                    # 500
times <- as.numeric(colnames(dv_mat))

round(E, 2)
#>  0.1 0.25  0.5    1    2    3    5    8   12 
#> 0.97 1.94 2.79 3.02 2.26 1.65 1.06 0.75 0.51
```

`admVignetteDvMatrix()` is a helper defined in this vignette’s setup
file; it does nothing but reshape `examplomycin`’s individual records
into one row per subject and one column per time.

`V`’s off-diagonal entries capture within-subject correlation across
time. Using the full matrix (`method = "cov"`) typically tightens
estimates against the diagonal-only approximation (`method = "var"`);
admixr2 picks the method from the structure of `V`.

## Model definition

Models use standard nlmixr2 syntax with mu-referenced log-scale
parameters:

``` r

pk_model <- function() {
  ini({
    tcl     <- log(5)  ; label("Log clearance (L/hr)")
    tv1     <- log(10) ; label("Log central volume (L)")
    tv2     <- log(30) ; label("Log peripheral volume (L)")
    tq      <- log(10) ; label("Log inter-compartmental CL (L/hr)")
    tka     <- log(1)  ; label("Log absorption rate constant (1/hr)")
    prop.sd <- c(0, 0.2); label("Proportional residual error SD")
    eta.cl ~ 0.09
    eta.v1 ~ 0.09
    eta.v2 ~ 0.09
    eta.q  ~ 0.09
    eta.ka ~ 0.09
  })
  model({
    cl <- exp(tcl + eta.cl)
    v1 <- exp(tv1 + eta.v1)
    v2 <- exp(tv2 + eta.v2)
    q  <- exp(tq  + eta.q)
    ka <- exp(tka + eta.ka)
    d/dt(depot)      <- -ka * depot
    d/dt(central)    <- ka * depot - (cl/v1 + q/v1) * central + (q/v2) * peripheral
    d/dt(peripheral) <- (q/v1) * central - (q/v2) * peripheral
    cp <- central / v1
    cp ~ prop(prop.sd)
  })
}
```

Writing each parameter as `exp(tcl + eta.cl)` is **mu-referencing**: the
fixed effect and its random effect enter additively on the log scale.
admixr2 uses that pairing to get analytical gradients from sensitivity
equations — see [Advanced
usage](https://leidenpharmacology.github.io/admixr2/articles/advanced.html#mu-referencing-and-sensitivity-equations),
which also covers parameters with no random effect.

## Writing down the study

[`admStudy()`](https://leidenpharmacology.github.io/admixr2/reference/admStudy.md)
takes the statistics and checks them; printing it is the read-back:

``` r
study <- admStudy(
  E     = E,
  V     = V,          # full 9x9 covariance matrix (ML denominator, as computed)
  n     = n,
  times = times,
  dose  = 100         # single 100 mg oral dose; `ev = rxode2::et(...)` for more
)
study
admixr2 study 'study'  -- digitised DATA
  design    n = 500, dose 100, 9 times (0.1 - 12)
  V scale   ml denominator
  reported  mean profile, full covariance
```

## Fitting

Collect the studies with
[`admStudies()`](https://leidenpharmacology.github.io/admixr2/reference/admStudies.md)
and pass them to
[`adghControl()`](https://leidenpharmacology.github.io/admixr2/reference/adghControl.md).
Gauss-Hermite quadrature is deterministic, so the same call always gives
the same fit:

``` r

fit <- nlmixr2(
  pk_model, admData(), est = "adgh",
  control = adghControl(
    studies = admStudies(examplomycin = study),
    maxeval = 300L,
    cores   = 2L
  )
)
```

## Inspecting the fit

``` r
print(fit)
── nlmixr² adgh ──

          OBJF       AIC       BIC Log-likelihood
adgh -3690.048 -3668.048 -3597.518       1845.024

── Time (sec fit$time): ──

  optimize covariance other elapsed
1   10.471      5.618     0  16.089

── Population Parameters (fit$parFixed or fit$parFixedDf): ──

                                  Parameter    Est.       SE  %RSE
tcl                    Log clearance (L/hr)   1.602  0.02017 1.259
tv1                  Log central volume (L)   2.329   0.1364 5.858
tv2               Log peripheral volume (L)   3.397  0.05492 1.617
tq        Log inter-compartmental CL (L/hr)   2.276  0.02737 1.203
tka     Log absorption rate constant (1/hr) 0.03030   0.1272 419.6
prop.sd      Proportional residual error SD  0.1895 0.003282 1.732
        Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
tcl        4.964 (4.771, 5.164)    32.53         NaN
tv1        10.27 (7.858, 13.41)    32.91         NaN
tv2        29.86 (26.81, 33.26)    31.53         NaN
tq         9.739 (9.230, 10.28)    33.66         NaN
tka       1.031 (0.8034, 1.323)    32.10         NaN
prop.sd 0.1895 (0.1831, 0.1960)                     
 
  Covariance Type (fit$covMethod): r,s
  Some strong fixed parameter correlations exist (fit$cor) :
                cor:tv1,tcl             cor:tv2,tcl              cor:tq,tcl 
                 0.382                  -0.540                   0.212  
            cor:tka,tcl         cor:prop.sd,tcl       cor:om.eta.cl,tcl 
                 0.406                  0.0592                 -0.0974  
      cor:om.eta.v1,tcl       cor:om.eta.v2,tcl        cor:om.eta.q,tcl 
                -0.361                  -0.186                   0.275  
      cor:om.eta.ka,tcl             cor:tv2,tv1              cor:tq,tv1 
                 0.348                  -0.859                   0.236  
            cor:tka,tv1         cor:prop.sd,tv1       cor:om.eta.cl,tv1 
                 0.984                 -0.0139                -0.00635  
      cor:om.eta.v1,tv1       cor:om.eta.v2,tv1        cor:om.eta.q,tv1 
                -0.711                  0.0566                   0.507  
      cor:om.eta.ka,tv1              cor:tq,tv2             cor:tka,tv2 
                 0.735                  -0.275                  -0.865  
        cor:prop.sd,tv2       cor:om.eta.cl,tv2       cor:om.eta.v1,tv2 
               -0.0101                  0.0473                   0.663  
      cor:om.eta.v2,tv2        cor:om.eta.q,tv2       cor:om.eta.ka,tv2 
                0.0739                  -0.448                  -0.682  
             cor:tka,tq          cor:prop.sd,tq        cor:om.eta.cl,tq 
                 0.272                -0.00667                 -0.0487  
       cor:om.eta.v1,tq        cor:om.eta.v2,tq         cor:om.eta.q,tq 
                -0.143                  0.0728                   0.324  
       cor:om.eta.ka,tq         cor:prop.sd,tka       cor:om.eta.cl,tka 
                0.0780                 -0.0180                -0.00788  
      cor:om.eta.v1,tka       cor:om.eta.v2,tka        cor:om.eta.q,tka 
                -0.705                  0.0530                   0.532  
      cor:om.eta.ka,tka   cor:om.eta.cl,prop.sd   cor:om.eta.v1,prop.sd 
                 0.724                 -0.0278                 -0.0490  
  cor:om.eta.v2,prop.sd    cor:om.eta.q,prop.sd   cor:om.eta.ka,prop.sd 
                -0.198                  -0.158                  0.0199  
cor:om.eta.v1,om.eta.cl cor:om.eta.v2,om.eta.cl  cor:om.eta.q,om.eta.cl 
                0.0367                  -0.169                -0.00235  
cor:om.eta.ka,om.eta.cl cor:om.eta.v2,om.eta.v1  cor:om.eta.q,om.eta.v1 
               -0.0326                 -0.0335                  -0.282  
cor:om.eta.ka,om.eta.v1  cor:om.eta.q,om.eta.v2 cor:om.eta.ka,om.eta.v2 
                -0.907                  -0.137                  0.0447  
 cor:om.eta.ka,om.eta.q 
                 0.252  
 

  No correlations in between subject variability (BSV) matrix
  Full BSV covariance (fit$omega) or correlation (fit$omegaR; diagonals=SDs) 
  Distribution stats (mean/skewness/kurtosis/p-value) available in fit$shrink 
  Information about run found (fit$runInfo):
   • "r,s" covariance is conditional on the ETAs; setCov() is marginal 
  Censoring (fit$censInformation): No censoring
  Minimization message (fit$message):  
    NLOPT_FTOL_REACHED: Optimization stopped because ftol_rel or ftol_abs (above) was reached. 
```

Key entries in `fit$env$admExtra`:

``` r

fit$objective                    # -2 log-likelihood
#> [1] -3690.048
fit$env$admExtra$struct          # structural parameters (log scale)
#>       tcl       tv1       tv2        tq       tka 
#> 1.6021766 2.3289123 3.3965932 2.2761540 0.0303032
fit$env$admExtra$sigma_var       # residual variance(s)
#>    prop.sd 
#> 0.03592593

logLik(fit)
#> 'log Lik.' 1845.024 (df=11)
AIC(fit)
#> [1] -3668.048
```

The estimated between-subject covariance matrix `Omega`:

``` r

knitr::kable(fit$env$admExtra$omega, digits = 4,
             caption = "Estimated Omega (between-subject covariance).")
```

|        |        |        |        |        |
|-------:|-------:|-------:|-------:|-------:|
| 0.1006 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 0.0000 | 0.1028 | 0.0000 | 0.0000 | 0.0000 |
| 0.0000 | 0.0000 | 0.0948 | 0.0000 | 0.0000 |
| 0.0000 | 0.0000 | 0.0000 | 0.1074 | 0.0000 |
| 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0981 |

Estimated Omega (between-subject covariance). {.table}

## Diagnostic plots

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws up to
four panel types and returns them as a named list of ggplot2 objects:

``` r

plots <- plot(fit, which = c("mean", "nll"))
```

![Mean diagnostics (observed vs predicted, with residuals), then the NLL
convergence trace.](admixr2_files/figure-html/plot-1.png)

Mean diagnostics (observed vs predicted, with residuals), then the NLL
convergence trace.

![Mean diagnostics (observed vs predicted, with residuals), then the NLL
convergence trace.](admixr2_files/figure-html/plot-2.png)

Mean diagnostics (observed vs predicted, with residuals), then the NLL
convergence trace.

All four panel types and their customisation are covered in
[`vignette("diagnostic-plots")`](https://leidenpharmacology.github.io/admixr2/articles/diagnostic-plots.md).

## Where to next

- [From a published figure to E, V and
  n](https://leidenpharmacology.github.io/admixr2/articles/aggregate-data.html)
  — turn published summaries into `E`, `V` and `n`
- [Simulating data & using published
  models](https://leidenpharmacology.github.io/admixr2/articles/datagen.html)
  — the model-as-input path
- [Multiple
  studies](https://leidenpharmacology.github.io/admixr2/articles/multiple-studies.html)
  — meta-analysis across several studies
- [Estimator
  comparison](https://leidenpharmacology.github.io/admixr2/articles/estimator-comparison.html)
  — choosing a backend
- [Choosing a residual error
  model](https://leidenpharmacology.github.io/admixr2/articles/error-models.html)
  — what `cp ~ prop(...)` does to `E` and `V`
