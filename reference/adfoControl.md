# Control settings for the FO (First-Order) estimator

Creates a control object for `nlmixr2(est = "adfo")`. The FO estimator
linearises model predictions at \\\eta = 0\\: it is faster than the MC
estimator but less accurate for models with large IIV or strongly
non-linear individual predictions.

## Usage

``` r
adfoControl(
  studies = list(),
  grad = c("none", "analytical", "fd", "cfd"),
  algorithm = NULL,
  maxeval = 500L,
  ftol_rel = .Machine$double.eps^(1/2),
  print = 10L,
  seed = 12345L,
  cores = 1L,
  nDisplayProgress = .Machine$integer.max,
  grad_h = 1e-04,
  grad_bounds = 5,
  cov_h = 0.001,
  cov_h_outer = .Machine$double.eps^(1/5),
  covMethod = c("r", "none"),
  n_restarts = 1L,
  restart_sd = 0.5,
  workers = 1L,
  rxControl = NULL,
  calcTables = FALSE,
  compress = TRUE,
  ci = 0.95,
  sigdig = 4,
  sigdigTable = NULL,
  addProp = c("combined2", "combined1"),
  optExpression = TRUE,
  sumProd = FALSE,
  literalFix = TRUE,
  returnAdmr = FALSE,
  ...
)
```

## Arguments

- studies:

  Named list of study specifications (same format as
  [`admControl()`](https://leidenpharmacology.github.io/admixr2/reference/admControl.md):
  `E`, `V`, `n`, `times`, `ev`, optional `method`; or an `observations`
  list for multi-compartment fits – see
  [`admControl()`](https://leidenpharmacology.github.io/admixr2/reference/admControl.md)).

- grad:

  Gradient mode. `"none"` (default) uses derivative-free BOBYQA;
  `"analytical"` uses the closed-form FO gradient (requires sensitivity
  equations); `"fd"` uses forward finite differences of the full NLL;
  `"cfd"` uses central finite differences for struct theta gradient
  (more accurate than `"fd"`, roughly twice as many NLL evaluations per
  step).

- algorithm:

  nloptr algorithm, or `NULL` (default) to pick the default that matches
  `grad`: `"NLOPT_LD_LBFGS"` with a gradient, `"NLOPT_LN_BOBYQA"` when
  `grad = "none"`. Any algorithm reported by
  [`nloptr::nloptr.print.options()`](https://astamm.github.io/nloptr/reference/nloptr.print.options.html)
  is accepted. An explicit algorithm is reconciled with `grad`: when
  `grad = "none"` a gradient-based algorithm (`NLOPT_LD_*` /
  `NLOPT_GD_*`) falls back to `"NLOPT_LN_BOBYQA"`; when a gradient is
  requested a derivative-free algorithm (`NLOPT_LN_*` / `NLOPT_GN_*`)
  turns the gradient off. Both emit a message.

- maxeval:

  Maximum function evaluations (default 500).

- ftol_rel:

  Relative tolerance (default `sqrt(.Machine$double.eps)`).

- print:

  Print-frequency for live progress (0 = silent).

- seed:

  Random seed (used for restarts).

- cores:

  OpenMP threads for
  [`rxSolve()`](https://nlmixr2.github.io/rxode2/reference/rxSolve.html)
  (default 1).

- nDisplayProgress:

  Passed to
  [`rxSolve()`](https://nlmixr2.github.io/rxode2/reference/rxSolve.html):
  show the solver's text progress bar only once a single solve exceeds
  this many subjects. The default (`.Machine$integer.max`) keeps it off
  for clean script/vignette output; lower it (e.g. `1000L`) to see
  progress during long fits.

- grad_h:

  Finite-difference step for unpaired struct theta gradient and FD
  Jacobian.

- grad_bounds:

  Box-constraint half-width when using gradients.

- cov_h:

  Inner FD step for the gradient-based Hessian (only used when
  `covMethod = "r"` and `grad != "none"`). Default 1e-3.

- cov_h_outer:

  Outer step scale for NLL-FD Hessian.

- covMethod:

  `"r"` computes covariance via numerical Hessian for structural and
  residual-error parameters only (omega/IIV SEs are not computed,
  consistent with nlmixr2 FOCEI); `"none"` skips it.

- n_restarts:

  Number of optimizer restarts (1 = no multi-start).

- restart_sd:

  Standard deviation for random perturbations of initial struct thetas
  at each restart (\> 1).

- workers:

  Number of parallel PSOCK/fork workers for multi-restart (default 1 =
  sequential).

- rxControl:

  [`rxode2::rxControl()`](https://nlmixr2.github.io/rxode2/reference/rxSolve.html)
  object. Created automatically when `NULL`.

- calcTables, compress, ci, sigdig, sigdigTable, optExpression, sumProd,
  literalFix:

  Passed to
  [`nlmixr2est::foceiControl()`](https://nlmixr2.github.io/nlmixr2est/reference/foceiControl.html)
  for the table/output machinery.

- addProp:

  How combined additive+proportional error is parameterised in the
  nlmixr2 output tables: `"combined2"` (default, variance form) or
  `"combined1"` (SD form). Has no effect on admixr2's own estimation.

- returnAdmr:

  If `TRUE`, return a plain list instead of the full nlmixr2 fit object.

- ...:

  Unused arguments (trigger an error).

## Value

An `adfoControl` object (a named list).

## See also

[`admControl()`](https://leidenpharmacology.github.io/admixr2/reference/admControl.md),
[`adirmcControl()`](https://leidenpharmacology.github.io/admixr2/reference/adirmcControl.md)

## Examples

``` r
# Inspect defaults
ctl <- adfoControl()
ctl$grad
#> [1] "none"
ctl$maxeval
#> [1] 500

# Analytical gradient, more evaluations
ctl2 <- adfoControl(grad = "analytical", maxeval = 1000L)

# \donttest{
library(rxode2)
#> rxode2 5.1.2 using 2 threads (see ?getRxThreads)
#>   no cache: create with `rxCreateCache()`
library(nlmixr2)
#> ── Attaching packages ───────────────────────────────────────── nlmixr2 5.0.0 ──
#> ✔ lotri        1.0.4     ✔ nlmixr2extra 5.1.0
#> ✔ nlmixr2data  2.0.9     ✔ nlmixr2plot  5.0.2
#> ✔ nlmixr2est   6.0.1     
#> ── Optional Packages Loaded/Ignored ─────────────────────────── nlmixr2 5.0.0 ──
#> ✖ babelmixr2     ✖ nonmem2rx
#> ✖ ggPMX     ✖ posologyr
#> ✖ monolix2rx     ✖ shinyMixR
#> ✖ nlmixr2lib     ✖ xpose.nlmixr2
#> ✖ nlmixr2rpt     
#> ── Conflicts ───────────────────────────────────────────── nlmixr2conflicts() ──
#> ✖ nlmixr2est::boxCox()     masks rxode2::boxCox()
#> ✖ nlmixr2est::yeoJohnson() masks rxode2::yeoJohnson()

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
#>   Obs units: 1 | Params: 5 | Cores: 1 | Grad: none | Restarts: 1
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
#> | 3.9 sec  |          |          |          |          |          |          |
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
#> 1    3.889      0.689   4.578
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
