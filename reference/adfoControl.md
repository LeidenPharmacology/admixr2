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
  cores = rxode2::rxCores(),
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
  [`rxSolve()`](https://nlmixr2.github.io/rxode2/reference/rxSolve.html).
  Defaults to
  [`rxode2::rxCores()`](https://nlmixr2.github.io/rxode2/reference/getRxThreads.html).
  When `workers > 1` it is a *total* budget, split across the workers.

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

  Number of parallel workers (mirai daemons) for multi-restart (default
  1 = sequential). Requires the `mirai` package.

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

## Installing memuse

[`rxode2::rxSolve()`](https://nlmixr2.github.io/rxode2/reference/rxSolve.html)
estimates free RAM on every call. When the `memuse` package is not
installed its fallback ends up shelling out to `vm_stat`, a macOS-only
command, so on Windows and Linux every solve spawns a process that can
only fail. Because the FO estimator issues many small solves, this
overhead is measurable (roughly 17% of an FO gradient). Installing
`memuse` makes the fallback unreachable:

    install.packages("memuse")

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
#> ✔ lotri        1.0.4      ✔ nlmixr2extra 5.1.0 
#> ✔ nlmixr2data  2.0.10     ✔ nlmixr2plot  5.0.2 
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
#>   Obs units: 1 | Params: 5 | Cores: 2 | Grad: none | Restarts: 1
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
#> | 1.5 sec  |          |          |          |          |          |          |
#>   Computing covariance (R method, 19 NLL evaluations)
#>   Note: covMethod='r' computes covariance for structural and sigma parameters only; omega (IIV) SEs are not computed (matching nlmixr2 FOCEI behavior).
#> → compress origData in nlmixr2 object, save 1160
#>  
#>  
print(fit)
#> ── nlmixr² adfo ──
#> 
#>          OBJF      AIC      BIC Log-likelihood
#> adfo 819.7109 829.7109 861.7701      -409.8555
#> 
#> ── Time (sec fit$time): ──
#> 
#>         optimize covariance elapsed other
#> elapsed    1.459      0.107   1.566 2.314
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
