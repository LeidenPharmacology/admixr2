# Diagnostic plots for an admixr2 fit

Generates up to four diagnostic panels:

## Usage

``` r
# S3 method for class 'admFit'
plot(x, which = c("mean", "cov", "nll", "par"), n_sim = NULL, seed = 1L, ...)
```

## Arguments

- x:

  An `admFit` object returned by `nlmixr2()` with `est = "adfo"`,
  `est = "admc"`, or `est = "adirmc"`.

- which:

  Character vector selecting which panel types to produce. Any subset of
  `c("mean", "cov", "nll", "par")`. Defaults to all four.

- n_sim:

  Number of MC samples for the final prediction. Defaults to the value
  used during fitting. Only used when `"mean"` or `"cov"` is in `which`.

- seed:

  Random seed for reproducibility.

- ...:

  Unused.

## Value

A named list of ggplot2 objects, invisibly. Prints each selected plot.

## Details

1.  `"mean"` – Observed vs predicted mean per study (2x2 grid). Upper
    row: observed and predicted mean lines with +/-1 SD ribbon on a
    shared y scale (black throughout). Lower row: raw residual lollipop
    with +/-2 SE band and standardised residual z-scores with +/-1.96
    reference lines.

2.  `"cov"` – Observed vs predicted (co)variance heatmaps per study (2x2
    grid). Upper row shares a common colour scale (blue-white-red).
    Lower row uses distinct diverging scales: residual (red-white-green)
    and standardised residual (gold-white-purple). Significance stars
    overlaid on the standardised residual panel.

3.  `"nll"` – NLL trace per restart over optimizer evaluations. Restarts
    coloured with the Okabe-Ito palette.

4.  `"par"` – Parameter trace per restart on the natural scale (struct
    thetas back-transformed, sigma as SD, omega diagonal as variance
    labelled `V(eta.x)`). Facets ordered as in the model
    [`ini()`](https://nlmixr2.github.io/rxode2/reference/ini.html)
    block. Restarts coloured with the Okabe-Ito palette.

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
#>  
#>  
#>  
#>  
#> === admixr2: Aggregate Data Modeling (FO) ===
#>   Studies: 1 | Params: 5 | Cores: 1 | Grad: none | Restarts: 1
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
#> | 4.3 sec  |          |          |          |          |          |          |
#>   Computing covariance (R method, 19 NLL evaluations)
#> → compress origData in nlmixr2 object, save 1120
plot(fit)
#>  
#>  




# }
```
