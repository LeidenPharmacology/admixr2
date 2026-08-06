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
#> | 0010     |  1.5e+29 |        5 |        1 |      0.2 |     0.09 |     0.04 |
#> | 0020     |  2351.54 |        5 |     23.2 |   0.2232 |     0.09 |    1.588 |
#> | 0030     |  1826.01 |    5.246 |    27.66 |   0.2112 |  0.05607 |   0.9401 |
#> | 0040     |  1693.60 |    5.429 |    29.31 |   0.2266 |   0.0586 |   0.9714 |
#> | 0050     |  1459.82 |    5.752 |     33.9 |   0.2796 |  0.04258 |   0.6821 |
#> | 0060     |  1339.28 |     5.96 |    35.23 |   0.3241 |  0.03856 |   0.5019 |
#> | 0070     |  1323.48 |    5.712 |    40.83 |   0.3624 |   0.0265 |   0.2476 |
#> | 0080     |  1137.41 |    5.855 |    35.91 |   0.3459 |  0.02343 |   0.1307 |
#> | 0090     |  1050.92 |    5.746 |    37.23 |   0.3671 |  0.02109 |  0.06168 |
#> | 0100     |  1018.46 |    5.853 |    37.91 |     0.39 |  0.02121 |  0.04717 |
#> | 0102 ✓   |  1018.46 |    5.853 |    37.91 |     0.39 |  0.02121 |  0.04717 |
#> | 1.1 sec  |          |          |          |          |          |          |
#>   Computing covariance (R method, 51 NLL evaluations)
#> → compress origData in nlmixr2 object, save 1160
#>  
#>  
print(fit)
#> ── nlmixr² adfo ──
#> 
#>          OBJF      AIC      BIC Log-likelihood
#> adfo 1018.459 1028.459 1060.518      -509.2295
#> 
#> ── Time (sec fit$time): ──
#> 
#>   optimize covariance other elapsed
#> 1    1.096      0.106     0   1.202
#> 
#> ── Population Parameters (fit$parFixed or fit$parFixedDf): ──
#> 
#>           Est.       SE   %RSE Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
#> tcl      1.767 0.008944 0.5062    5.853 (5.751, 5.957)    14.64         NaN
#> tv       3.635  0.01258 0.3461    37.91 (36.99, 38.86)    21.98         NaN
#> prop.sd 0.3900 0.006145  1.576 0.3900 (0.3779, 0.4020)                     
#>  
#>   Covariance Type (fit$covMethod): r
#>   No correlations in between subject variability (BSV) matrix
#>   Full BSV covariance (fit$omega) or correlation (fit$omegaR; diagonals=SDs) 
#>   Distribution stats (mean/skewness/kurtosis/p-value) available in fit$shrink 
#>   Information about run found (fit$runInfo):
#>    • adfoCalcCov: the full Hessian including omega was not positive definite; reporting structural and sigma standard errors only. 
#>   Censoring (fit$censInformation): No censoring
#>   Minimization message (fit$message):  
#>     NLOPT_MAXEVAL_REACHED: Optimization stopped because maxeval (above) was reached. 
# }
```
