# Control settings for the ADM estimator

Constructs a control object for `est = "admc"`, the Monte Carlo
aggregate data modelling estimator.

## Usage

``` r
admControl(
  studies = list(),
  n_sim = 5000L,
  sampling = c("sobol", "halton", "torus", "lhs", "rnorm"),
  algorithm = NULL,
  maxeval = 500L,
  ftol_rel = .Machine$double.eps^2,
  print = 10L,
  seed = 12345L,
  cores = rxode2::rxCores(),
  nDisplayProgress = .Machine$integer.max,
  grad = c("sens", "fd", "cfd", "none"),
  grad_h = 1e-04,
  cov_h = 0.001,
  cov_h_outer = .Machine$double.eps^(1/5),
  grad_bounds = 5,
  covMethod = c("r", "none"),
  cov_n_sim = 10000L,
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

  Named list of study specifications. Each element is a list with:

  - `E` – observed mean vector

  - `V` – observed covariance matrix or variance vector (auto-detected)

  - `n` – sample size

  - `times` – numeric vector of observation times

  - `ev` –
    [`rxode2::et()`](https://nlmixr2.github.io/rxode2/reference/et.html)
    dosing event table

  - `method` – `"cov"` or `"var"` (optional; auto-detected from `V`)

  **Multi-compartment (multiple observed outputs).** To fit several
  observed compartments simultaneously (e.g. plasma and brain/CSF), give
  the study an `observations` list instead of top-level `E`/`V`/`times`.
  Each entry is one observed output with its own `output` (the model
  prediction variable, e.g. `"cp"` or `"cCSF"`), `times`, `E`, `V` and –
  for independent fits – `ev` and `n`. Pass the endpoint names to
  [`admData()`](https://leidenpharmacology.github.io/admixr2/reference/admData.md),
  e.g. `admData(c("cp", "cCSF"))`, so nlmixr2 recognises every endpoint.
  There are two modes:

  - *Independent* – each observed output has its own `n`/`ev` (separate
    experiments / subjects, e.g. a plasma study and a brain study
    combined for meta-analysis). The outputs are independent likelihood
    blocks and the aggregate `-2LL` is their sum.

  - *Joint (same subjects)* – the outputs are measured on the SAME
    subjects. Give the study a shared `n` and `ev`, and a joint
    covariance either as a study-level full matrix `V` (blocks in
    `observations` order) or as per-output marginal `V` plus a `cross`
    list of cross-covariance blocks keyed `"outA:outB"` (each
    `length(times_A)` x `length(times_B)`; omitted pairs are zero). The
    compartments are then scored by a single MVN over the stacked vector
    with shared random effects. `est = "adirmc"` does not support
    multiple observed outputs; use `"admc"`, `"adfo"` or `"adgh"`.

  **Long format (one row per endpoint/time).** As an alternative to the
  `observations` list, a study may carry a `data` frame that keys each
  observed summary by endpoint, the way nlmixr2 keys observations by
  `DVID`/`CMT`. The frame needs an endpoint column (`DVID`, `CMT` or
  `output`), a time column (`TIME`), a mean column (`E`) and – unless a
  joint `V` is given – a variance column (`V`) or an SD column (`SD`).
  It is normalised into exactly the same units as the `observations`
  form, so the two are interchangeable:

      # independent blocks: per-row variances; optional per-endpoint `n` column
      # and per-endpoint `ev` (a list of event tables keyed by endpoint)
      list(n = 60L, ev = ev,
           data = data.frame(DVID = c("cp", "cp", "cCSF"), TIME = c(1, 2, 2),
                             E = c(9.1, 7.4, 2.2), V = c(1.2, 0.9, 0.1)))

      # joint (same subjects): ONE stacked covariance whose rows/cols align with
      # the rows of `data` -- no `cross` blocks to assemble by hand
      list(n = 60L, ev = ev, data = data.frame(DVID = ..., TIME = ..., E = ...),
           V = V_joint)

  A study-level `V` (or an explicit `joint = TRUE`) marks the endpoints
  as same-subject; without one, each endpoint is an independent
  likelihood block. Endpoints are stacked in the order they first appear
  in `data`.

- n_sim:

  Number of Monte Carlo samples per NLL evaluation.

- sampling:

  Sampling method for eta draws: `"sobol"` (Sobol, default), `"halton"`
  (Halton), `"torus"` (Kronecker/torus), `"lhs"` (Latin hypercube), or
  `"rnorm"` (iid normal).

- algorithm:

  nloptr algorithm string, or `NULL` (default) to pick the default that
  matches `grad`: `"NLOPT_LD_LBFGS"` with a gradient,
  `"NLOPT_LN_BOBYQA"` when `grad = "none"`. Any algorithm reported by
  [`nloptr::nloptr.print.options()`](https://astamm.github.io/nloptr/reference/nloptr.print.options.html)
  is accepted (e.g. `"NLOPT_LD_MMA"`, `"NLOPT_LN_NELDERMEAD"`). An
  explicit algorithm is reconciled with `grad`: when `grad = "none"` a
  gradient-based algorithm (`NLOPT_LD_*` / `NLOPT_GD_*`) falls back to
  `"NLOPT_LN_BOBYQA"`; when a gradient is requested a derivative-free
  algorithm (`NLOPT_LN_*` / `NLOPT_GN_*`) turns the gradient off. Both
  emit a message.

- maxeval:

  Maximum number of optimizer function evaluations.

- ftol_rel:

  Relative function-value tolerance for convergence.

- print:

  Print progress every this many evaluations (0 = silent).

- seed:

  Random seed for reproducibility.

- cores:

  Number of OpenMP threads for
  [`rxSolve()`](https://nlmixr2.github.io/rxode2/reference/rxSolve.html).
  Defaults to
  [`rxode2::rxCores()`](https://nlmixr2.github.io/rxode2/reference/getRxThreads.html).
  [`rxSolve()`](https://nlmixr2.github.io/rxode2/reference/rxSolve.html)
  parallelises over subjects, so this is the main speed lever for the MC
  estimators; when `workers > 1` it is a *total* budget, split across
  the workers.

- nDisplayProgress:

  Passed to
  [`rxSolve()`](https://nlmixr2.github.io/rxode2/reference/rxSolve.html):
  the solver shows its text progress bar only once a single solve
  exceeds this many subjects. The default (`.Machine$integer.max`) keeps
  the bar off, which is what you want for scripts, vignettes and logs;
  lower it (e.g. `1000L`) to see solver progress during long interactive
  fits.

- grad:

  Gradient mode: `"sens"` (sensitivity equations, default), `"fd"`
  (forward finite differences), `"cfd"` (central finite differences), or
  `"none"` (derivative-free). A warning is issued when `"sens"` is
  requested but the sensitivity model is unavailable; the estimator then
  falls back to forward finite differences.

- grad_h:

  Step size for finite-difference gradient evaluation during
  optimization (used by `grad = "fd"` or `"cfd"`). The default 1e-4 is
  near the optimal balance between truncation error (grows with `h`) and
  MC noise amplification (grows as `1/h`) for forward FD. Central FD
  (`"cfd"`) has a slightly wider optimum around 1e-3, but 1e-4 works
  well for both.

- cov_h:

  Inner FD step for the gradient-based Hessian (only used when
  `covMethod = "r"` and `grad != "none"`). Each gradient evaluation has
  MC noise of order `sigma / cov_h`; the Hessian divides that noise by
  the outer step, giving total noise
  `sigma / (cov_h * cov_h_outer * |p|)`. `cov_h = 1e-3` balances
  truncation error and noise amplification. Increase to `1e-2` if the
  Hessian is non-positive definite.

- cov_h_outer:

  Outer step scale for the numerical Hessian. The actual step for
  parameter `p` is `max(|p|, 0.1) * cov_h_outer`. Applied to both the
  gradient-FD Hessian (`grad != "none"`) and the NLL-FD Hessian
  (`grad = "none"`). Default `eps^(1/5)` (~2.5e-3) is larger than the
  textbook `eps^(1/4)` to account for MC noise in NLL and gradient
  evaluations; empirically it matches the analytical
  (sensitivity-equation) Hessian ground truth. Increase (e.g. to `5e-3`
  or `1e-2`) if the Hessian is non-positive definite.

- grad_bounds:

  Box-constraint half-width when using gradients.

- covMethod:

  Covariance method: `"r"` (numerical Hessian for structural and
  residual-error parameters only; omega/IIV SEs are not computed,
  consistent with nlmixr2 FOCEI) or `"none"`.

- cov_n_sim:

  Number of MC samples for the covariance (Hessian) step. More samples
  reduce MC noise in NLL evaluations. The NLL-based Hessian
  (`grad = "none"`) uses a central second difference of the NLL with the
  same Sobol sequence (CRN) at every perturbed point, so noise largely
  cancels and `cov_n_sim = 10000` (default) is sufficient for most
  models.

- n_restarts:

  Number of optimization restarts. Runs in parallel when `workers > 1`.

- restart_sd:

  Standard deviation of structural theta perturbations for restart
  initialisation.

- workers:

  Number of parallel workers for multi-restart. `1` (default) runs
  restarts sequentially. Values `> 1` run the restarts on a pool of
  background R processes (mirai daemons), which behaves the same way on
  every platform. Requires the `mirai` package. Workers are stopped
  automatically after the restart phase so all cores are available for
  the Hessian step; if a fit is interrupted,
  [`admStopWorkers()`](https://leidenpharmacology.github.io/admixr2/reference/admStopWorkers.md)
  cleans up any survivors.

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
  `"combined1"` (SD form). Has no effect on admixr2's own estimation;
  passed to
  [`nlmixr2est::foceiControl()`](https://nlmixr2.github.io/nlmixr2est/reference/foceiControl.html)
  for the table/output machinery only.

- returnAdmr:

  If `TRUE`, return a plain list instead of a full nlmixr2 fit object
  (useful for debugging).

- ...:

  Additional arguments (none allowed; triggers an error).

## Value

An object of class `admControl`.

## Examples

``` r
# Minimal control object -- inspect defaults
ctl <- admControl()
ctl$n_sim
#> [1] 5000
ctl$algorithm
#> [1] "NLOPT_LD_LBFGS"

# Override key settings without fitting
ctl2 <- admControl(
  n_sim    = 2000L,
  maxeval  = 300L,
  grad     = "fd",
  seed     = 42L
)

# \donttest{
library(rxode2)
library(nlmixr2)

data("examplomycin")
obs   <- examplomycin[examplomycin$EVID == 0, ]
obs   <- obs[order(obs$ID, obs$TIME), ]
times <- sort(unique(obs$TIME))
ids   <- unique(obs$ID)
dv_mat <- do.call(rbind, lapply(ids, function(i) {
  sub <- obs[obs$ID == i, ]; sub$DV[order(sub$TIME)]
}))
E <- colMeans(dv_mat)
V <- cov.wt(dv_mat, method = "ML")$cov

pk_model <- function() {
  ini({
    tcl <- log(5);  tv1 <- log(12); tv2 <- log(25)
    tq  <- log(12); tka <- log(1.2)
    prop.sd <- c(0, 0.2)
    eta.cl ~ 0.09; eta.v1 ~ 0.09; eta.v2 ~ 0.09
    eta.q  ~ 0.09; eta.ka ~ 0.09
  })
  model({
    cl <- exp(tcl + eta.cl); v1 <- exp(tv1 + eta.v1)
    v2 <- exp(tv2 + eta.v2); q  <- exp(tq  + eta.q)
    ka <- exp(tka + eta.ka)
    d/dt(depot)      <- -ka * depot
    d/dt(central)    <- ka * depot - (cl/v1 + q/v1) * central + (q/v2) * peripheral
    d/dt(peripheral) <- (q/v1) * central - (q/v2) * peripheral
    cp <- central / v1
    cp ~ prop(prop.sd)
  })
}

fit <- nlmixr2(
  pk_model, admData(), est = "admc",
  control = admControl(
    studies  = list(study1 = list(E = E, V = V, n = length(ids),
                                  times = times, ev = et(amt = 100))),
    n_sim    = 1000L,
    maxeval  = 200L
  )
)
#>  
#>  
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> === admixr2: Aggregate Data Modeling (MC) ===
#>   Obs units: 1 | MC samples: 1000 | Params: 11 | Cores: 2 | Grad: Sens | Restarts: 1
#> +----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+
#> |          |     -2LL |      tcl |      tv1 |      tv2 |       tq |      tka |  prop.sd |   eta.cl |   eta.v1 |   eta.v2 |    eta.q |   eta.ka |
#> +----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+
#> | 0010     | -3644.91 |    5.149 |    11.61 |    26.17 |    10.87 |    1.243 |   0.2005 |   0.0903 |  0.09021 |  0.09008 |  0.09075 |   0.0903 |
#> | 0020     | -3689.85 |    4.995 |    10.84 |    29.13 |    9.758 |     1.09 |   0.1994 |   0.1081 |   0.1038 |  0.09577 |   0.1045 |   0.1068 |
#> | 0030     | -3690.56 |    4.967 |     10.3 |    29.77 |    9.843 |    1.045 |   0.1983 |   0.1047 |    0.103 |  0.09964 |   0.1097 |   0.1028 |
#> | 0040     | -3690.69 |    4.951 |     10.2 |    30.01 |    9.823 |    1.034 |   0.1983 |   0.1044 |   0.1128 |   0.1038 |   0.1066 |  0.09513 |
#> | 0050     | -3690.73 |    4.952 |    10.11 |    30.04 |    9.813 |    1.026 |   0.1983 |   0.1044 |   0.1156 |   0.1009 |    0.106 |  0.09254 |
#> | 0052 ✓   | -3690.73 |    4.952 |    10.11 |    30.04 |    9.815 |    1.026 |   0.1983 |   0.1044 |   0.1155 |   0.1009 |    0.106 |  0.09258 |
#> | 6.6 sec  |          |          |          |          |          |          |          |          |          |          |          |          |
#>   Computing covariance (R method, Sens-Hessian, 7 gradient evaluations)
#>   Note: covMethod='r' computes covariance for structural and sigma parameters only; omega (IIV) SEs are not computed (matching nlmixr2 FOCEI behavior).
#> → compress origData in nlmixr2 object, save 1160
#>  
#>  
print(fit)
#> ── nlmixr² admc ──
#> 
#>           OBJF       AIC       BIC Log-likelihood
#> admc -3690.728 -3668.728 -3598.198       1845.364
#> 
#> ── Time (sec fit$time): ──
#> 
#>   optimize covariance elapsed
#> 1      6.6      4.369  10.969
#> 
#> ── Population Parameters (fit$parFixed or fit$parFixedDf): ──
#> 
#>            Est.      SE  %RSE Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
#> tcl         1.6 0.01658 1.037    4.952 (4.793, 5.115)     33.2            
#> tv1       2.314 0.08724 3.771       10.11 (8.522, 12)     35.0            
#> tv2       3.403 0.04058 1.193    30.04 (27.74, 32.53)     32.6            
#> tq        2.284 0.02142 0.938    9.815 (9.411, 10.24)     33.4            
#> tka     0.02529 0.08189 323.7   1.026 (0.8735, 1.204)     31.1            
#> prop.sd  0.1983                                0.1983                     
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
