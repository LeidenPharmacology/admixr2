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
  resid_nodes = 81L,
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

  `"r"` computes covariance via a numerical Hessian over the structural,
  residual-error and omega parameters; `"none"` skips it. Omega is
  included because excluding it also biases the STRUCTURAL standard
  errors downward – a theta carrying an eta is correlated with that
  eta's variance. If the weakly-identified omega Cholesky makes the
  Hessian non-positive definite, the structural + residual sub-block is
  reported with a warning.

  All three blocks are reported on the scale the ESTIMATES are printed
  on, as `nlmixr2est` does: structural thetas on the log/optimizer
  scale, residual error as an SD, and omega as the variance/covariance
  entries (named `om.<eta>` and `cov.<eta_i>.<eta_j>`). The omega block
  is rotated by the full Jacobian of Omega with respect to the
  log-Cholesky, which is not diagonal once omega is correlated.

  **An adfo standard error describes scatter, not accuracy.** FO
  linearises the model at eta = 0, and on a non-additive residual (or a
  saturating endpoint, or a large omega) the resulting point estimates
  carry a bias of several standard errors – measured 5-20 SE, giving 0%
  coverage for a nominal 95% interval even where the SE itself matches
  the sampling SD. Use `adgh` or `admc` when the uncertainty matters.

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

- resid_nodes:

  Gauss-Hermite nodes used to integrate the RESIDUAL for a
  transform-both-sides endpoint (`boxCox`, `yeoJohnson`, `logitNorm`,
  `probitNorm`), where `y = g(h(f) + sigma*eps)` has no closed-form mean
  and variance. Ignored by every other error model, which has closed
  forms. Default 81. Measured worst-case relative error against an
  independent quadrature, over all four transforms and residual SD of
  0.5, 1, 2 and 3: n = 15 gives 5.7e-2, 31 gives 4.5e-3, 81 gives
  5.0e-5. The error is dominated by large residual SD; at SD \<= 1, n =
  31 already gives 1e-7 or better.

  This is an ACCURACY dial, not a speed one. The quadrature is linear in
  `resid_nodes` in isolation (~50 us at 15, 300 us at 81 for an 8-row
  study) but negligible beside the ODE solve: a full NLL evaluation
  measured 0.750 s per 60 evaluations at BOTH 31 and 81 nodes. Raise it
  if you have a saturating endpoint with a large residual SD; there is
  little to gain by lowering it.

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
#> rxode2 5.1.5 using 2 threads (see ?getRxThreads)
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
#> | 1.2 sec  |          |          |          |          |          |          |
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
#>         optimize covariance elapsed other
#> elapsed    1.245      0.093   1.338 3.337
#> 
#> ── Population Parameters (fit$parFixed or fit$parFixedDf): ──
#> 
#>          Est.       SE   %RSE Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
#> tcl     1.767 0.008944 0.5062    5.853 (5.751, 5.957)     14.6            
#> tv      3.635  0.01258 0.3461    37.91 (36.99, 38.86)     22.0            
#> prop.sd  0.39 0.006145  1.576    0.39 (0.3779, 0.402)                     
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
