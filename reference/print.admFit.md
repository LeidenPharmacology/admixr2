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
#>   Obs units: 1 | Params: 5 | Cores: 2 | Grad: none | Restarts: 1
#> +----------+----------+----------+----------+----------+----------+----------+
#> |          |     -2LL |      tcl |       tv |  prop.sd |   eta.cl |    eta.v |
#> +----------+----------+----------+----------+----------+----------+----------+
#> | 0010     |  1.5e+30 |        5 |        1 |      0.2 |     0.09 |     0.04 |
#> | 0020     |  3878.99 |    4.803 |    31.64 |      0.2 |     0.09 |     0.04 |
#> | 0030     |  3303.72 |    4.806 |     30.1 |   0.2068 |  0.09261 |     0.04 |
#> | 0040     |  2828.57 |    4.816 |    28.56 |   0.2164 |  0.09217 |  0.04034 |
#> | 0050     |   931.70 |    6.841 |       37 |    0.381 |   0.1389 |  0.02184 |
#> | 0060     |   828.72 |    6.712 |    39.14 |   0.4072 |   0.1429 |  0.02032 |
#> | 0070     |   822.99 |    6.742 |    39.64 |   0.4167 |   0.1405 |  0.02019 |
#> | 0080     |   821.51 |    6.716 |    39.82 |    0.422 |   0.1389 |  0.02072 |
#> | 0090     |   816.00 |    6.677 |    39.33 |   0.4187 |   0.1361 |  0.02537 |
#> | 0100     |   813.39 |    6.705 |    39.36 |   0.4168 |   0.1346 |  0.02999 |
#> | 0102 ✓   |   813.39 |    6.705 |    39.36 |   0.4168 |   0.1346 |  0.02999 |
#> | 1.0 sec  |          |          |          |          |          |          |
#>   Computing covariance (R method, 19 NLL evaluations)
#>   Note: covMethod='r' computes covariance for structural and sigma parameters only; omega (IIV) SEs are not computed (matching nlmixr2 FOCEI behavior).
#> → compress origData in nlmixr2 object, save 1160
#>  
#>  
print(fit)
#> ── nlmixr² adfo ──
#> 
#>          OBJF      AIC      BIC Log-likelihood
#> adfo 813.3949 823.3949 855.4541      -406.6975
#> 
#> ── Time (sec fit$time): ──
#> 
#>   optimize covariance elapsed
#> 1     0.95      0.088   1.038
#> 
#> ── Population Parameters (fit$parFixed or fit$parFixedDf): ──
#> 
#>           Est.      SE   %RSE Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
#> tcl      1.903 0.01231 0.6467    6.705 (6.546, 6.869)     38.0            
#> tv       3.673 0.01025  0.279    39.36 (38.58, 40.16)     17.4            
#> prop.sd 0.4168                                 0.4168                     
#>  
#>   Covariance Type (fit$covMethod): r
#>   No correlations in between subject variability (BSV) matrix
#>   Full BSV covariance (fit$omega) or correlation (fit$omegaR; diagonals=SDs) 
#>   Distribution stats (mean/skewness/kurtosis/p-value) available in fit$shrink 
#>   Censoring (fit$censInformation): No censoring
#>   Minimization message (fit$message):  
#>     NLOPT_MAXEVAL_REACHED: Optimization stopped because maxeval (above) was reached. 
# }
```
