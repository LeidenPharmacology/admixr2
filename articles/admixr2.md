# Getting started with admixr2

## What is aggregate-data modelling?

`admixr2` fits pharmacometric PK/PD models to **summary-level data**.
For each clinical study you supply:

- **E** — observed mean vector (one entry per observation time)
- **V** — observed covariance matrix (or variance vector)
- **n** — sample size
- **times** — observation time points
- **ev** — dosing event table

The estimators match E and V against their model-predicted counterparts
and return a standard nlmixr2 fit object. This lets you apply
established nlmixr2 models to aggregate statistics from publications or
internal data summaries where individual records are unavailable.

Four estimators are available:

| Estimator | `est =` | Control function | Approach |
|----|----|----|----|
| First-Order | `"adfo"` | [`adfoControl()`](https://leidenpharmacology.github.io/admixr2/reference/adfoControl.md) | First-order Taylor expansion at η = 0; one rxSolve per NLL eval; fastest |
| Monte Carlo | `"admc"` | [`admControl()`](https://leidenpharmacology.github.io/admixr2/reference/admControl.md) | Sample average over η; asymptotically exact |
| Gauss-Hermite | `"adgh"` | [`adghControl()`](https://leidenpharmacology.github.io/admixr2/reference/adghControl.md) | Deterministic quadrature over η; noise-free, unbiased at any IIV |
| Iterative Reweighting MC | `"adirmc"` | [`adirmcControl()`](https://leidenpharmacology.github.io/admixr2/reference/adirmcControl.md) | Proposals fixed per phase; inner loop needs no new rxSolve calls |

`adfo` is the natural starting point for model screening and initial
estimates. `admc` is the workhorse for standard PK models. `adgh` is a
noise-free alternative to `admc` for models with up to ~4 etas. `adirmc`
is preferred for complex ODE systems with expensive solves,
high-dimensional IIV, or poor starting values. See
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
 [1m──  [34mnlmix [39m [31mr² [39m  [33madmc [39m ── [22m

          OBJF       AIC       BIC Log-likelihood
admc -3690.835 -3668.835 -3598.305       1845.418

 [1m── Time (sec  [33mfit [39m [34m$time [39m): ── [22m

  optimize covariance elapsed
1   23.516      5.081  28.597

 [1m── Population Parameters ( [33mfit [39m [34m$parFixed [39m or  [33mfit [39m [34m$parFixedDf [39m): ── [22m

                                   [1m [1mParameter [0m [0m     [1mEst. [0m       [1m [1mSE [0m [0m    [1m%RSE [0m
 [1m [1mtcl [0m [0m                    Log clearance (L/hr)   1.601 0.01635  1.021
 [1m [1mtv1 [0m [0m                  Log central volume (L)   2.314 0.08719  3.768
 [1m [1mtv2 [0m [0m               Log peripheral volume (L)   3.402 0.04007  1.178
 [1m [1mtq [0m [0m        Log inter-compartmental CL (L/hr)   2.285 0.02132 0.9332
 [1m [1mtka [0m [0m     Log absorption rate constant (1/hr) 0.02423 0.08198  338.4
 [1m [1mprop.sd [0m [0m      Proportional residual error SD  0.1984               
         [1mBack-transformed(95%CI) [0m  [1mBSV(CV%) [0m  [1mShrink(SD)% [0m
 [1m [1mtcl [0m [0m         4.958 (4.802, 5.12)     32.8            
 [1m [1mtv1 [0m [0m           10.12 (8.528, 12)     33.8            
 [1m [1mtv2 [0m [0m        30.03 (27.76, 32.48)     32.0            
 [1m [1mtq [0m [0m          9.822 (9.42, 10.24)     33.4            
 [1m [1mtka [0m [0m       1.025 (0.8725, 1.203)     31.2            
 [1m [1mprop.sd [0m [0m                  0.1984                     
 
  Covariance Type ( [33mfit [39m [1m [34m$covMethod [39m [22m):  [1mr [22m
  No correlations in between subject variability (BSV) matrix
  Full BSV covariance ( [33mfit [39m [1m [34m$omega [39m [22m) or correlation ( [33mfit [39m [1m [34m$omegaR [39m [22m; diagonals=SDs) 
  Distribution stats (mean/skewness/kurtosis/p-value) available in  [33mfit [39m [1m [34m$shrink [39m [22m 
  Censoring ( [33mfit [39m [1m [34m$censInformation [39m [22m): No censoring
  Minimization message ( [33mfit [39m [1m [34m$message [39m [22m):  
    NLOPT_XTOL_REACHED: Optimization stopped because xtol_rel or xtol_abs (above) was reached. 
```

Key entries in `fit$env$admExtra`:

``` r

fit$objective                    # -2 log-likelihood
#> [1] -3690.835
fit$env$admExtra$struct          # structural parameters (log scale)
#>       tcl       tv1       tv2        tq       tka 
#> 1.6010696 2.3142405 3.4022370 2.2845989 0.0242259
fit$env$admExtra$omega           # estimated Omega matrix
#>           [,1]      [,2]       [,3]     [,4]       [,5]
#> [1,] 0.1021589 0.0000000 0.00000000 0.000000 0.00000000
#> [2,] 0.0000000 0.1080018 0.00000000 0.000000 0.00000000
#> [3,] 0.0000000 0.0000000 0.09747445 0.000000 0.00000000
#> [4,] 0.0000000 0.0000000 0.00000000 0.105588 0.00000000
#> [5,] 0.0000000 0.0000000 0.00000000 0.000000 0.09281048
fit$env$admExtra$sigma_var       # residual variance(s)
#>    prop.sd 
#> 0.03937635

logLik(fit)
#> 'log Lik.' 1845.418 (df=11)
AIC(fit)
#> [1] -3668.835
```

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
