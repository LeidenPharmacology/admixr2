# Getting started with admixr2

## What is aggregate-data modelling?

`admixr2` fits pharmacometric PK/PD models to **summary-level data**
instead of individual records — a **meta-analysis** framework for
population PK/PD. The inputs can be **digitised aggregate data** from
published studies (means, error bars and covariances), **previously
published PK/PD models**, or a mix of both. The result is a single
unified population model with interpretable fixed, random and covariate
effects, recovered without individual patient data.

For each study you supply:

- **E** — observed mean vector (one entry per observation time)
- **V** — observed covariance matrix (or variance vector)
- **n** — sample size
- **times** — observation time points
- **ev** — dosing event table

The estimators match E and V against their model-predicted counterparts
and return a standard nlmixr2 fit object, so established nlmixr2 models
apply to aggregate statistics from publications or internal summaries
where individual records are unavailable. Turning a published figure
into `E`, `V` and `n` is covered in
[`vignette("aggregate-data", package = "admixr2")`](https://leidenpharmacology.github.io/admixr2/articles/aggregate-data.md);
letting each study carry its own published model as the input is covered
in
[`vignette("datagen", package = "admixr2")`](https://leidenpharmacology.github.io/admixr2/articles/datagen.md).

Four estimators are available:

| Estimator | `est =` | Control function | Approach |
|----|----|----|----|
| First-Order | `"adfo"` | [`adfoControl()`](https://leidenpharmacology.github.io/admixr2/reference/adfoControl.md) | First-order Taylor expansion at η = 0; one rxSolve per NLL eval; fastest |
| Monte Carlo | `"admc"` | [`admControl()`](https://leidenpharmacology.github.io/admixr2/reference/admControl.md) | Sample average over η; asymptotically exact |
| Gauss-Hermite | `"adgh"` | [`adghControl()`](https://leidenpharmacology.github.io/admixr2/reference/adghControl.md) | Deterministic quadrature over η; noise-free, unbiased at any IIV |
| Iterative Reweighting MC | `"adirmc"` | [`adirmcControl()`](https://leidenpharmacology.github.io/admixr2/reference/adirmcControl.md) | Proposals fixed per phase; inner loop needs no new rxSolve calls |

`adfo` is the natural starting point for model screening and initial
estimates. `admc` is the workhorse for standard PK models. `adgh` is a
noise-free alternative to `admc`, most efficient when the number of
random effects is small. `adirmc` is preferred for complex ODE systems
with expensive solves, high-dimensional IIV, or poor starting values.
See
[`vignette("estimator-comparison", package = "admixr2")`](https://leidenpharmacology.github.io/admixr2/articles/estimator-comparison.md)
for a detailed comparison.

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

obs   <- examplomycin[examplomycin$EVID == 0, ]
obs   <- obs[order(obs$ID, obs$TIME), ]
times <- sort(unique(obs$TIME))
ids   <- unique(obs$ID)
n     <- length(ids)                     # 500

dv_mat <- matrix(NA_real_, nrow = n, ncol = length(times))
for (i in seq_along(ids)) {
  sub         <- obs[obs$ID == ids[i], ]
  dv_mat[i, ] <- sub$DV[order(sub$TIME)]
}

E <- colMeans(dv_mat)
V <- cov.wt(dv_mat, method = "ML")$cov

round(E, 2)
#> [1] 0.97 1.94 2.79 3.02 2.26 1.65 1.06 0.75 0.51
```

`V` is the 9×9 sample covariance matrix. Its off-diagonal entries
capture within-subject correlation across time; using the full matrix
(`method = "cov"`) typically tightens parameter estimates compared to
the diagonal-only approximation (`method = "var"`). `admixr2`
auto-detects the method from the structure of V.

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

Writing each parameter as `exp(tcl + eta.cl)` is called
**mu-referencing**: the structural fixed effect and its random effect
enter additively on the log scale. `admixr2` exploits this pairing to
compute analytical gradients via sensitivity equations. See the
[Advanced
usage](https://leidenpharmacology.github.io/admixr2/articles/advanced.html#mu-referencing-and-sensitivity-equations)
vignette for details, including how parameters without a random effect
are handled.

## Assembling the study specification

Bundle each study’s statistics into a named list:

``` r

study <- list(
  E     = E,
  V     = V,                       # full 9x9 covariance matrix
  n     = n,
  times = times,
  ev    = rxode2::et(amt = 100)    # single 100 mg oral dose
)
```

## Fitting

Pass one or more named studies to
[`admControl()`](https://leidenpharmacology.github.io/admixr2/reference/admControl.md):

``` r

fit <- nlmixr2(
  pk_model, admData(), est = "admc",
  control = admControl(
    studies   = list(examplomycin = study),
    n_sim     = 5000L,
    cov_n_sim = 10000L,
    maxeval   = 300L,
    seed      = 1L
  )
)
```

## Inspecting the fit

``` r
print(fit)
── nlmixr² admc ──

          OBJF       AIC       BIC Log-likelihood
admc -3690.254 -3668.254 -3597.724       1845.127

── Time (sec fit$time): ──

  optimize covariance elapsed
1    18.19      8.588  26.778

── Population Parameters (fit$parFixed or fit$parFixedDf): ──

                                  Parameter    Est.      SE  %RSE
tcl                    Log clearance (L/hr)   1.602 0.01956 1.221
tv1                  Log central volume (L)   2.328  0.1207 5.183
tv2               Log peripheral volume (L)   3.398 0.05196 1.529
tq        Log inter-compartmental CL (L/hr)   2.276 0.02652 1.165
tka     Log absorption rate constant (1/hr) 0.02992  0.1138 380.5
prop.sd      Proportional residual error SD  0.1895 0.00321 1.694
        Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
tcl        4.963 (4.776, 5.157)     32.6            
tv1           10.26 (8.098, 13)     32.8            
tv2         29.9 (27.01, 33.11)     32.0            
tq          9.74 (9.247, 10.26)     33.7            
tka        1.03 (0.8243, 1.288)     32.2            
prop.sd 0.1895 (0.1832, 0.1958)                     
 
  Covariance Type (fit$covMethod): r
  No correlations in between subject variability (BSV) matrix
  Full BSV covariance (fit$omega) or correlation (fit$omegaR; diagonals=SDs) 
  Distribution stats (mean/skewness/kurtosis/p-value) available in fit$shrink 
  Censoring (fit$censInformation): No censoring
  Minimization message (fit$message):  
    NLOPT_XTOL_REACHED: Optimization stopped because xtol_rel or xtol_abs (above) was reached. 
```

Key entries in `fit$env$admExtra`:

``` r

fit$objective                    # -2 log-likelihood
#> [1] -3690.254
fit$env$admExtra$struct          # structural parameters (log scale)
#>        tcl        tv1        tv2         tq        tka 
#> 1.60193579 2.32808009 3.39796638 2.27623769 0.02991704
fit$env$admExtra$sigma_var       # residual variance(s)
#>    prop.sd 
#> 0.03591429

logLik(fit)
#> 'log Lik.' 1845.127 (df=11)
AIC(fit)
#> [1] -3668.254
```

The estimated between-subject covariance matrix `Omega`:

``` r

knitr::kable(fit$env$admExtra$omega, digits = 4,
             caption = "Estimated Omega (between-subject covariance).")
```

|        |       |        |        |        |
|-------:|------:|-------:|-------:|-------:|
| 0.1011 | 0.000 | 0.0000 | 0.0000 | 0.0000 |
| 0.0000 | 0.102 | 0.0000 | 0.0000 | 0.0000 |
| 0.0000 | 0.000 | 0.0975 | 0.0000 | 0.0000 |
| 0.0000 | 0.000 | 0.0000 | 0.1073 | 0.0000 |
| 0.0000 | 0.000 | 0.0000 | 0.0000 | 0.0986 |

Estimated Omega (between-subject covariance). {.table}

## Diagnostic plots

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) produces up to
four panel types and returns them as a named list of ggplot2 objects:

``` r

plots <- plot(fit, which = c("mean", "nll"))
```

![Left: observed vs predicted mean with residuals. Right: NLL
convergence trace.](admixr2_files/figure-html/plot-1.png)

Left: observed vs predicted mean with residuals. Right: NLL convergence
trace.

![Left: observed vs predicted mean with residuals. Right: NLL
convergence trace.](admixr2_files/figure-html/plot-2.png)

Left: observed vs predicted mean with residuals. Right: NLL convergence
trace.

For a detailed walkthrough of all four panel types and customisation
options, see
[`vignette("diagnostic-plots", package = "admixr2")`](https://leidenpharmacology.github.io/admixr2/articles/diagnostic-plots.md).

## Where to next

- [From a published figure to E, V and
  n](https://leidenpharmacology.github.io/admixr2/articles/aggregate-data.md)
  — turn published summaries into `E`, `V` and `n`
- [Simulating data & using published
  models](https://leidenpharmacology.github.io/admixr2/articles/datagen.md)
  — the model-as-input path
- [Multiple
  studies](https://leidenpharmacology.github.io/admixr2/articles/multiple-studies.md)
  — meta-analysis across several studies
- [Estimator
  comparison](https://leidenpharmacology.github.io/admixr2/articles/estimator-comparison.md)
  — choosing a backend
- [Choosing a residual error
  model](https://leidenpharmacology.github.io/admixr2/articles/error-models.md)
  — what `cp ~ prop(...)` does to `E` and `V`
