# Control settings for the IRMC estimator

Constructs a control object for `est = "adirmc"`, the Iterative
Reweighting Monte Carlo estimator.

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
  kappa_method = c("first-order", "linear"),
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

  Kappa correction method for models with non-mu-referenced struct
  thetas: `"first-order"` (default, population prediction `f(theta, 0)`,
  one rxSolve per inner step) or `"linear"` (first-order Taylor
  expansion — precomputes `J = df/d(theta)` once per outer iteration,
  zero rxSolve per inner step).

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
V <- diag(diag(cov.wt(dv_mat, method = "ML")$cov))

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
#> | 0001     |  1764.39 |    4.278 |    3.685 |    38.18 |    8.767 |   0.8694 |   0.1872 |   0.1653 |    0.665 |    0.665 |  0.05283 |  0.09956 |
#> | 0002     | -1231.42 |    4.817 |    8.687 |     30.5 |    9.044 |   0.8399 |   0.1594 |   0.1168 |   0.1085 |     0.09 |  0.06498 |   0.1675 |
#> | 0003     | -1264.96 |     4.95 |    7.932 |    31.87 |    8.796 |   0.7944 |   0.1752 |   0.1264 |  0.06648 |   0.0982 |  0.02324 |   0.1434 |
#> | 0004     | -1266.20 |    4.915 |    7.996 |    31.74 |    8.866 |   0.8038 |   0.1777 |   0.1222 |  0.07718 |    0.109 |  0.02161 |    0.143 |
#> | 0005 ✓   | -1266.42 |    4.916 |    8.011 |     31.8 |    8.878 |   0.8039 |   0.1773 |   0.1229 |  0.07669 |   0.1089 |   0.0216 |   0.1403 |
#> +-- Phase 2: Focused (+/-1.00) ----------------------------------------------------------------------------------------------------------------+
#> | 0006 ✓   | -1266.42 |    4.913 |    8.008 |    31.79 |    8.879 |   0.8041 |   0.1773 |   0.1229 |   0.0767 |   0.1089 |   0.0216 |   0.1404 |
#> +-- Phase 3: Fine-tuning (+/-0.50) ------------------------------------------------------------------------------------------------------------+
#> | 0007     | -1266.30 |    4.938 |    8.075 |    31.51 |    8.902 |    0.811 |   0.1803 |   0.1215 |  0.07672 |  0.09915 |   0.0213 |   0.1384 |
#> | 0008     | -1266.40 |    4.924 |    8.062 |    31.57 |    8.894 |   0.8097 |     0.18 |   0.1233 |  0.07969 |   0.1039 |  0.02059 |   0.1383 |
#> | 0009 ✓   | -1266.45 |    4.924 |    8.033 |    31.65 |    8.889 |   0.8069 |     0.18 |   0.1228 |  0.07934 |   0.1034 |  0.02057 |   0.1372 |
#> +-- Phase 4: Precision (+/-0.01) --------------------------------------------------------------------------------------------------------------+
#> | 0010 ✓   | -1266.45 |    4.923 |    8.032 |    31.65 |     8.89 |    0.807 |     0.18 |   0.1228 |  0.07934 |   0.1034 |  0.02057 |   0.1372 |
#> | 1.5 sec  |          |          |          |          |          |          |          |          |          |          |          |          |
#>   Computing covariance (R method, MC NLL, Sens-Hessian, 7 gradient evaluations)
#> → compress origData in nlmixr2 object, save 1120
print(fit)
#> ── nlmixr² adirmc ──
#> 
#>             OBJF       AIC       BIC Log-likelihood
#> adirmc -1266.449 -1244.449 -1173.919       633.2246
#> 
#> ── Time (sec fit$time): ──
#> 
#>   optimize covariance elapsed
#> 1    1.506      9.301  10.807
#> 
#> ── Population Parameters (fit$parFixed or fit$parFixedDf): ──
#> 
#>            Est.      SE  %RSE Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
#> tcl       1.594 0.02539 1.593    4.923 (4.684, 5.175)     36.1            
#> tv1       2.083  0.1895 9.093    8.032 (5.541, 11.64)     28.7            
#> tv2       3.455 0.07219  2.09    31.65 (27.47, 36.46)     33.0            
#> tq        2.185 0.05916 2.708     8.89 (7.916, 9.983)     14.4            
#> tka     -0.2145  0.1742  81.2   0.807 (0.5736, 1.135)     38.4            
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
