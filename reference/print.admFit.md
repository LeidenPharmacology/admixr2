# Print method for admFit objects

Delegates to `print.nlmixr2FitCore` for the standard nlmixr2 coloured
output. `admFit` class is kept on the object during the call so that
`head.admFit` intercepts any `head(fit)` calls that arise in the paged-
output path (R Markdown / notebooks), preventing the
`[.data.frame(.subset2(env, integer))` crash that occurs when an
environment-backed fit is subscripted like a plain list.

## Usage

``` r
# S3 method for class 'admFit'
print(x, ...)
```

## Arguments

- x:

  An `admFit` object.

- ...:

  Passed to `print.nlmixr2FitCore`.

## Value

`x`, invisibly.

## Examples

``` r
# \donttest{
library(rxode2)
library(nlmixr2)

data("examplomycin")
obs    <- examplomycin[examplomycin$EVID == 0, ]
obs    <- obs[order(obs$ID, obs$TIME), ]
times  <- sort(unique(obs$TIME))
ids    <- unique(obs$ID)
dv_mat <- do.call(rbind, lapply(ids, function(i) {
  sub <- obs[obs$ID == i, ]; sub$DV[order(sub$TIME)]
}))
E <- colMeans(dv_mat)
V <- cov.wt(dv_mat, method = "ML")$cov

pk_model <- function() {
  ini({
    tcl <- log(5); tv <- log(30)
    prop.sd <- c(0, 0.2)
    eta.cl ~ 0.09; eta.v ~ 0.04
  })
  model({
    cl <- exp(tcl + eta.cl)
    v  <- exp(tv  + eta.v)
    d/dt(central) <- -(cl/v) * central
    cp <- central / v
    cp ~ prop(prop.sd)
  })
}

fit <- nlmixr2(
  pk_model, admData(), est = "adfo",
  control = adfoControl(
    studies = list(study1 = list(E = E, V = V, n = length(ids),
                                 times = times, ev = et(amt = 100))),
    maxeval = 100L
  )
)
#>  
#>  
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> === admixr2: Aggregate Data Modeling (FO) ===
#>   Obs units: 1 | Params: 5 | Cores: 2 | Grad: Analytical | Restarts: 1
#> +----------+----------+----------+----------+----------+----------+----------+
#> |          |     -2LL |      tcl |       tv |  prop.sd |   eta.cl |    eta.v |
#> +----------+----------+----------+----------+----------+----------+----------+
#> | 0010     |  1768.15 |    4.967 |    29.88 |   0.2587 |   0.0888 |  0.04603 |
#> | 0020     |   862.47 |    6.391 |    37.74 |   0.3864 |  0.08003 |   0.0422 |
#> | 0029 ✓   |   861.90 |    6.384 |    38.03 |     0.39 |  0.08051 |  0.04074 |
#> | 0.7 sec  |          |          |          |          |          |          |
#>   Computing covariance (R method, Analytical-Hessian, 6 gradient evaluations)
#> → compress origData in nlmixr2 object, save 1160
#>  
#>  
print(fit)
#> ── nlmixr² adfo ──
#> 
#>          OBJF      AIC      BIC Log-likelihood
#> adfo 861.8956 871.8956 903.9548      -430.9478
#> 
#> ── Time (sec fit$time): ──
#> 
#>   optimize covariance other elapsed
#> 1    0.652      0.158     0    0.81
#> 
#> ── Population Parameters (fit$parFixed or fit$parFixedDf): ──
#> 
#>           Est.       SE   %RSE Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
#> tcl      1.854  0.01620 0.8742    6.384 (6.184, 6.590)    28.95         NaN
#> tv       3.638  0.01234 0.3391    38.03 (37.13, 38.96)    20.39         NaN
#> prop.sd 0.3900 0.006554  1.681 0.3900 (0.3771, 0.4028)                     
#>  
#>   Covariance Type (fit$covMethod): r
#>   No correlations in between subject variability (BSV) matrix
#>   Full BSV covariance (fit$omega) or correlation (fit$omegaR; diagonals=SDs) 
#>   Distribution stats (mean/skewness/kurtosis/p-value) available in fit$shrink 
#>   Censoring (fit$censInformation): No censoring
#>   Minimization message (fit$message):  
#>     NLOPT_XTOL_REACHED: Optimization stopped because xtol_rel or xtol_abs (above) was reached. 
# }
```
