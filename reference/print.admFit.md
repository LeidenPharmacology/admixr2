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
#>   Studies: 1 | Params: 5 | Cores: 1 | Grad: none | Restarts: 1
#> +----------+----------+----------+----------+----------+----------+----------+
#> |          |     -2LL |      tcl |       tv |  prop.sd |   eta.cl |    eta.v |
#> +----------+----------+----------+----------+----------+----------+----------+
#> | 0010     |  1.5e+30 |        5 |        1 |      0.2 |     0.09 |     0.04 |
#> | 0020     |  3878.99 |    4.803 |    31.64 |      0.2 |     0.09 |     0.04 |
#> | 0030     |  3303.72 |    4.806 |     30.1 |   0.2068 |  0.09261 |     0.04 |
#> | 0040     |  2828.60 |    4.816 |    28.56 |   0.2164 |  0.09217 |  0.04034 |
#> | 0050     |   929.13 |    6.839 |    37.04 |   0.3811 |   0.1391 |  0.02192 |
#> | 0060     |   833.01 |    6.862 |    39.38 |   0.4076 |   0.1484 |  0.02135 |
#> | 0070     |   826.48 |    6.702 |    39.97 |   0.4139 |   0.1411 |  0.02057 |
#> | 0080     |   822.02 |    6.697 |    39.52 |   0.4177 |   0.1395 |  0.02036 |
#> | 0090     |   821.28 |    6.719 |    39.61 |   0.4213 |   0.1374 |  0.02043 |
#> | 0100     |   820.07 |     6.65 |    39.62 |   0.4181 |   0.1294 |  0.02105 |
#> | 0102 ✓   |   819.71 |    6.591 |    39.47 |   0.4158 |   0.1246 |  0.02138 |
#> | 3.5 sec  |          |          |          |          |          |          |
#>   Computing covariance (R method, 19 NLL evaluations)
#>   Note: covMethod='r' computes covariance for structural and sigma parameters only; omega (IIV) SEs are not computed (matching nlmixr2 FOCEI behavior).
#> → compress origData in nlmixr2 object, save 1160
print(fit)
#> ── nlmixr² adfo ──
#> 
#>          OBJF      AIC      BIC Log-likelihood
#> adfo 819.7109 829.7109 861.7701      -409.8555
#> 
#> ── Time (sec fit$time): ──
#> 
#>   optimize covariance elapsed
#> 1    3.503      0.656   4.159
#> 
#> ── Population Parameters (fit$parFixed or fit$parFixedDf): ──
#> 
#>           Est.       SE   %RSE Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
#> tcl      1.886  0.01222 0.6479    6.591 (6.435, 6.751)     36.4            
#> tv       3.676 0.009646 0.2624    39.47 (38.73, 40.22)     14.7            
#> prop.sd 0.4158                                  0.4158                     
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
