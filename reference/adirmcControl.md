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
  kappa_method = c("exact", "linearized"),
  grad_h = 1e-04,
  cov_h = 0.001,
  cov_h_outer = .Machine$double.eps^(1/5),
  phases = c(2, 1, 0.5, 0.01),
  convcrit = 1e-05,
  max_worse = 5L,
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
  thetas: `"exact"` (default, re-evaluates population prediction
  `f(theta, 0)` via rxSolve at each inner step) or `"linearized"`
  (precomputes `J = df/d(theta)` once per outer iteration, approximates
  kappa via linear expansion — zero rxSolve per inner step).

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
#> | 0001     |  1761.55 |    4.279 |    3.685 |    38.19 |    8.771 |   0.8692 |   0.1871 |   0.1652 |    0.665 |    0.665 |  0.05266 |  0.09961 |
#> | 0002     | -1243.14 |    4.843 |    8.558 |    30.68 |    8.997 |   0.8286 |   0.1629 |   0.1175 |   0.1081 |     0.09 |  0.05605 |   0.1567 |
#> | 0003     | -1265.36 |    4.947 |    8.729 |     31.2 |    9.204 |   0.8731 |    0.176 |   0.1212 |  0.08796 |   0.1019 |  0.05458 |   0.1281 |
#> | 0004     | -1265.30 |    4.921 |    8.055 |    31.62 |     8.85 |    0.808 |    0.174 |   0.1302 |  0.07557 |   0.1074 |  0.02433 |   0.1461 |
#> | 0005     | -1266.01 |    4.941 |    8.092 |    31.51 |    8.914 |   0.8127 |   0.1793 |   0.1206 |  0.07518 |  0.09887 |  0.02373 |   0.1389 |
#> | 0006     | -1266.43 |    4.942 |    8.164 |    31.31 |    8.933 |   0.8194 |   0.1798 |   0.1229 |  0.07595 |  0.09944 |  0.02365 |   0.1408 |
#> | 0007     | -1266.43 |    4.943 |    8.166 |    31.31 |     8.93 |   0.8191 |   0.1797 |   0.1228 |  0.07594 |  0.09945 |  0.02364 |   0.1407 |
#> | 0008 ✓   | -1266.43 |    4.942 |    8.166 |    31.31 |    8.931 |   0.8191 |   0.1797 |   0.1228 |  0.07594 |  0.09945 |  0.02364 |   0.1407 |
#> +-- Phase 2: Focused (+/-1.00) ----------------------------------------------------------------------------------------------------------------+
#> | 0009     | -1266.43 |    4.941 |    8.163 |    31.32 |    8.932 |   0.8191 |   0.1796 |   0.1228 |  0.07592 |   0.0995 |  0.02364 |   0.1404 |
#> | 0010 ✓   | -1266.43 |    4.941 |    8.162 |    31.32 |    8.932 |   0.8191 |   0.1796 |   0.1228 |  0.07592 |   0.0995 |  0.02364 |   0.1404 |
#> +-- Phase 3: Fine-tuning (+/-0.50) ------------------------------------------------------------------------------------------------------------+
#> | 0011 ✓   | -1266.43 |    4.942 |    8.161 |    31.32 |    8.934 |    0.819 |   0.1797 |   0.1228 |  0.07593 |  0.09952 |  0.02363 |   0.1404 |
#> +-- Phase 4: Precision (+/-0.01) --------------------------------------------------------------------------------------------------------------+
#> | 0012     | -1266.43 |    4.942 |    8.161 |    31.32 |    8.932 |    0.819 |   0.1797 |   0.1228 |  0.07593 |  0.09953 |  0.02363 |   0.1404 |
#> | 0013     | -1266.43 |    4.942 |     8.16 |    31.33 |    8.933 |   0.8189 |   0.1797 |   0.1228 |  0.07594 |  0.09956 |  0.02363 |   0.1404 |
#> | 0014     | -1266.44 |    4.941 |     8.16 |    31.33 |    8.932 |   0.8189 |   0.1797 |   0.1228 |  0.07594 |  0.09956 |  0.02363 |   0.1404 |
#> | 0015     | -1266.44 |    4.941 |    8.159 |    31.33 |    8.933 |   0.8188 |   0.1797 |   0.1228 |  0.07596 |  0.09959 |  0.02362 |   0.1404 |
#> | 0016     | -1266.44 |    4.941 |    8.159 |    31.33 |    8.932 |   0.8187 |   0.1797 |   0.1228 |  0.07596 |  0.09959 |  0.02362 |   0.1404 |
#> | 0017     | -1266.44 |    4.941 |    8.158 |    31.33 |    8.931 |   0.8186 |   0.1797 |   0.1228 |  0.07597 |  0.09962 |  0.02362 |   0.1404 |
#> | 0018     | -1266.44 |    4.941 |    8.157 |    31.33 |    8.931 |   0.8186 |   0.1797 |   0.1228 |  0.07597 |  0.09963 |  0.02362 |   0.1404 |
#> | 0019     | -1266.44 |    4.941 |    8.156 |    31.34 |    8.933 |   0.8186 |   0.1797 |   0.1228 |  0.07598 |  0.09965 |  0.02362 |   0.1404 |
#> | 0020     | -1266.44 |    4.941 |    8.156 |    31.33 |    8.931 |   0.8185 |   0.1797 |   0.1228 |  0.07598 |  0.09965 |  0.02362 |   0.1403 |
#> | 0021     | -1266.44 |    4.941 |    8.155 |    31.34 |    8.932 |   0.8185 |   0.1797 |   0.1228 |  0.07599 |  0.09967 |  0.02361 |   0.1403 |
#> | 0022     | -1266.44 |    4.941 |    8.155 |    31.34 |    8.931 |   0.8184 |   0.1797 |   0.1228 |    0.076 |  0.09968 |  0.02361 |   0.1403 |
#> | 0023     | -1266.44 |    4.941 |    8.154 |    31.34 |    8.932 |   0.8184 |   0.1797 |   0.1228 |  0.07601 |   0.0997 |  0.02361 |   0.1403 |
#> | 0024     | -1266.44 |    4.941 |    8.154 |    31.34 |    8.931 |   0.8184 |   0.1797 |   0.1228 |  0.07601 |  0.09971 |  0.02361 |   0.1403 |
#> | 0025     | -1266.44 |    4.941 |    8.153 |    31.34 |    8.932 |   0.8183 |   0.1797 |   0.1228 |  0.07602 |  0.09974 |   0.0236 |   0.1403 |
#> | 0026     | -1266.44 |     4.94 |    8.153 |    31.34 |    8.931 |   0.8183 |   0.1797 |   0.1228 |  0.07602 |  0.09974 |   0.0236 |   0.1403 |
#> | 0027     | -1266.44 |    4.936 |     8.12 |    31.41 |    8.923 |   0.8153 |   0.1798 |   0.1229 |  0.07679 |   0.1007 |  0.02341 |   0.1397 |
#> | 0028     | -1266.44 |    4.935 |    8.122 |    31.44 |    8.922 |   0.8153 |   0.1797 |   0.1227 |  0.07678 |   0.1007 |   0.0234 |   0.1396 |
#> | 0029     | -1266.44 |    4.935 |    8.121 |    31.44 |    8.922 |   0.8153 |   0.1797 |   0.1227 |  0.07678 |   0.1007 |   0.0234 |   0.1396 |
#> | 0030 ✓   | -1266.44 |    4.935 |    8.121 |    31.44 |    8.922 |   0.8153 |   0.1797 |   0.1227 |  0.07678 |   0.1008 |   0.0234 |   0.1396 |
#> | 2.2 sec  |          |          |          |          |          |          |          |          |          |          |          |          |
#>   Computing covariance (R method, MC NLL, Sens-Hessian, 7 gradient evaluations)
#>   Note: covMethod='r' computes covariance for structural and sigma parameters only; omega (IIV) SEs are not computed (matching nlmixr2 FOCEI behavior).
#> → compress origData in nlmixr2 object, save 1160
print(fit)
#> ── nlmixr² adirmc ──
#> 
#>             OBJF       AIC       BIC Log-likelihood
#> adirmc -1266.442 -1244.442 -1173.911       633.2208
#> 
#> ── Time (sec fit$time): ──
#> 
#>   optimize covariance elapsed
#> 1    2.177     10.364  12.541
#> 
#> ── Population Parameters (fit$parFixed or fit$parFixedDf): ──
#> 
#>            Est.      SE  %RSE Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
#> tcl       1.596 0.02478 1.552    4.935 (4.701, 5.181)     36.1            
#> tv1       2.094  0.1906 9.102     8.121 (5.589, 11.8)     28.3            
#> tv2       3.448 0.07124 2.066    31.44 (27.34, 36.15)     32.6            
#> tq        2.189  0.0581 2.655    8.922 (7.962, 9.999)     15.4            
#> tka     -0.2042  0.1757 86.06  0.8153 (0.5778, 1.151)     38.7            
#> prop.sd  0.1797                                0.1797                     
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
