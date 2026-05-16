# Control settings for the IRMC estimator

Constructs a control object for `est = "adirmc"`, the Importance
Resampling Monte Carlo estimator.

## Usage

``` r
adirmcControl(
  studies = list(),
  n_sim = 2500L,
  outer_iter = 50L,
  sampling = c("sobol", "halton", "torus", "lhs", "rnorm"),
  algorithm = "NLOPT_LN_BOBYQA",
  maxeval = 5000L,
  ftol_rel = .Machine$double.eps,
  print = 1L,
  omega_expansion = 1,
  seed = 12345L,
  cores = 1L,
  grad = c("analytical", "none", "fd"),
  kappa_method = c("first-order", "second-order"),
  grad_h = 1e-04,
  cov_h = 0.001,
  cov_h_outer = .Machine$double.eps^(1/5),
  phases = c(2, 1, 0.5, 0.01),
  convcrit = 0.05,
  max_worse = 3L,
  covMethod = c("r", "none"),
  cov_n_sim = 10000L,
  n_restarts = 1L,
  restart_sd = 0.2,
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

- n_sim:

  Number of Monte Carlo samples per NLL evaluation.

- outer_iter:

  Maximum inner optimiser iterations per phase.

- sampling:

  Sampling method for eta draws: `"sobol"` (Sobol, default), `"halton"`
  (Halton), `"torus"` (Kronecker/torus), `"lhs"` (Latin hypercube), or
  `"rnorm"` (iid normal).

- algorithm:

  nloptr algorithm string. Automatically switched to `"NLOPT_LD_LBFGS"`
  when `grad != "none"`.

- maxeval:

  Maximum number of optimizer function evaluations.

- ftol_rel:

  Relative function-value tolerance for convergence.

- print:

  Print progress every this many evaluations (0 = silent).

- omega_expansion:

  Inflate proposal Omega by this factor (\>= 1).

- seed:

  Random seed for reproducibility.

- cores:

  Number of OpenMP threads for
  [`rxSolve()`](https://nlmixr2.github.io/rxode2/reference/rxSolve.html).

- grad:

  Gradient mode for the inner optimiser: `"analytical"` (default,
  closed-form weight-path gradient), `"none"` (derivative-free BOBYQA),
  or `"fd"` (finite differences). Note: `"sens"` and `"cfd"` are not
  available for the IRMC estimator.

- kappa_method:

  Kappa correction method for non-mu-referenced models: `"first-order"`
  (default, population prediction at eta=0) or `"second-order"`
  (delta-method correction for Jensen's inequality, adds
  `(1/2) sum_k omega_k * d2f/deta_k2` to the reference mean).

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

- phases:

  Numeric vector of box-constraint half-widths, one per phase. Phases
  progressively tighten the search region.

- convcrit:

  Convergence criterion: phase ends when `|approx - exact| < convcrit`.

- max_worse:

  Stop a phase after this many consecutive worsening iterations.

- covMethod:

  Covariance method: `"r"` (numerical Hessian) or `"none"`.

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
  restarts sequentially. Values `> 1` use a PSOCK cluster on Windows and
  fork workers on Unix/macOS. Workers are stopped automatically after
  the restart phase so all cores are available for the Hessian step.

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

An object of class `adirmcControl`.

## Examples

``` r
# Inspect defaults
ctl <- adirmcControl()
ctl$phases
#> [1] 2.00 1.00 0.50 0.01
ctl$omega_expansion
#> [1] 1

# Tighter phases, more restarts
ctl2 <- adirmcControl(
  n_sim           = 1000L,
  omega_expansion = 1.5,
  phases          = c(2, 1, 0.5, 0.01),
  n_restarts      = 3L
)

# \donttest{
library(rxode2)
#> rxode2 5.0.2 using 2 threads (see ?getRxThreads)
#>   no cache: create with `rxCreateCache()`
library(nlmixr2)
#> ── Attaching packages ───────────────────────────────────────── nlmixr2 5.0.0 ──
#> ✔ lotri        1.0.4     ✔ nlmixr2extra 5.0.0
#> ✔ nlmixr2data  2.0.9     ✔ nlmixr2plot  5.0.1
#> ✔ nlmixr2est   5.0.2     
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
obs   <- examplomycin[examplomycin$EVID == 0, ]
obs   <- obs[order(obs$ID, obs$TIME), ]
times <- sort(unique(obs$TIME))
ids   <- unique(obs$ID)
dv_mat <- do.call(rbind, lapply(ids, function(i) {
  sub <- obs[obs$ID == i, ]; sub$DV[order(sub$TIME)]
}))
E <- colMeans(dv_mat)
V <- diag(diag(cov(dv_mat)))

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
  pk_model, admData(), est = "adirmc",
  control = adirmcControl(
    studies = list(study1 = list(E = E, V = V, n = length(ids),
                                 times = times, ev = et(amt = 100))),
    n_sim   = 500L
  )
)
#>  
#>  
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of full model...
#> ✔ done
#> → calculate jacobian
#> → calculate sensitivities
#> → calculate ∂(f)/∂(η)
#> → calculate ∂(R²)/∂(η)
#> → finding duplicate expressions in inner model...
#> → optimizing duplicate expressions in inner model...
#> → finding duplicate expressions in EBE model...
#> → optimizing duplicate expressions in EBE model...
#> → compiling inner model...
#>  
#>  
#> ✔ done
#> → finding duplicate expressions in FD model...
#> → optimizing duplicate expressions in FD model...
#> → compiling EBE model...
#>  
#>  
#> ✔ done
#> → compiling events FD model...
#>  
#>  
#> ✔ done
#>  
#>  
#>  
#>  
#> === admixr2: Aggregate Data Modeling (IR-MC) ===
#>   Studies: 1 | MC samples: 500 | Phases: 4 | Iters/phase: 50 | Expansion: 1.00 | Grad: analytic+Sens-Hessian | Restarts: 1
#> +----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+
#> |          |     -2LL |      tcl |      tv1 |      tv2 |       tq |      tka |  prop.sd |   eta.cl |   eta.v1 |   eta.v2 |    eta.q |   eta.ka |
#> +-- Phase 1: Wide (+/-2.00) -------------------------------------------------------------------------------------------------------------------+
#> | 0001     |  2394.83 |    3.993 |    2.604 |    37.52 |    8.863 |   0.7407 |   0.1931 |   0.2467 |    0.665 |    0.665 |  0.02795 |    0.135 |
#> | 0002     | -1006.70 |    4.629 |    8.384 |     32.1 |     9.94 |   0.8639 |   0.1828 |   0.1074 |     0.09 |    0.115 |   0.2065 |    0.247 |
#> | 0003     | -1254.25 |    4.962 |    8.193 |    31.16 |    8.859 |   0.8156 |   0.1827 |   0.1168 |   0.0762 |  0.07697 |  0.02795 |   0.1441 |
#> | 0004     | -1253.53 |    4.892 |    8.033 |    31.74 |    8.861 |   0.8039 |   0.1754 |   0.1313 |  0.07653 |   0.1217 |  0.02351 |   0.1459 |
#> | 0005     | -1256.51 |    4.947 |     8.06 |    31.52 |    8.898 |   0.8092 |   0.1793 |   0.1202 |  0.07351 |  0.09691 |   0.0236 |   0.1395 |
#> | 0006 ✓   | -1257.40 |    4.948 |    8.176 |    31.23 |     8.93 |   0.8205 |   0.1801 |   0.1234 |  0.07442 |  0.09763 |  0.02351 |   0.1422 |
#> +-- Phase 2: Focused (+/-1.00) ----------------------------------------------------------------------------------------------------------------+
#> | 0007 ✓   | -1257.41 |    4.949 |    8.179 |    31.23 |    8.928 |   0.8201 |     0.18 |   0.1233 |  0.07442 |  0.09764 |  0.02351 |   0.1422 |
#> +-- Phase 3: Fine-tuning (+/-0.50) ------------------------------------------------------------------------------------------------------------+
#> | 0008 ✓   | -1257.41 |    4.947 |    8.174 |    31.23 |     8.93 |     0.82 |   0.1799 |   0.1232 |  0.07441 |  0.09771 |   0.0235 |   0.1419 |
#> +-- Phase 4: Precision (+/-0.01) --------------------------------------------------------------------------------------------------------------+
#> | 0009 ✓   | -1257.41 |    4.948 |    8.172 |    31.23 |    8.932 |   0.8199 |     0.18 |   0.1232 |  0.07444 |  0.09777 |   0.0235 |   0.1419 |
#> | 1.0 sec  |          |          |          |          |          |          |          |          |          |          |          |          |
#>   Computing covariance (R method, MC NLL, Sens-Hessian, 7 gradient evaluations)
#> → compress origData in nlmixr2 object, save 1120
print(fit)
#> ── nlmixr² adirmc ──
#> 
#>             OBJF       AIC       BIC Log-likelihood
#> adirmc -1257.413 -1235.413 -1164.882       628.7063
#> 
#> ── Time (sec fit$time): ──
#> 
#>   optimize covariance elapsed
#> 1    1.007     10.286  11.293
#> 
#> ── Population Parameters (fit$parFixed or fit$parFixedDf): ──
#> 
#>            Est.      SE  %RSE Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
#> tcl       1.599 0.02433 1.521    4.948 (4.717, 5.189)     36.2            
#> tv1       2.101  0.1919 9.137      8.172 (5.61, 11.9)     27.8            
#> tv2       3.441 0.07038 2.045     31.23 (27.2, 35.85)     32.0            
#> tq         2.19 0.05773 2.637       8.932 (7.976, 10)     15.4            
#> tka     -0.1985  0.1772 89.28  0.8199 (0.5793, 1.161)     39.0            
#> prop.sd    0.18                                  0.18                     
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
