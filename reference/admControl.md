# Control settings for the ADM estimator

Constructs a control object for `est = "admc"`, the Monte Carlo
aggregate data modelling estimator.

## Usage

``` r
admControl(
  studies = list(),
  n_sim = 5000L,
  sampling = c("sobol", "halton", "torus", "lhs", "rnorm"),
  algorithm = "NLOPT_LN_BOBYQA",
  maxeval = 500L,
  ftol_rel = .Machine$double.eps^2,
  print = 10L,
  seed = 12345L,
  cores = 1L,
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

- n_sim:

  Number of Monte Carlo samples per NLL evaluation.

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

- seed:

  Random seed for reproducibility.

- cores:

  Number of OpenMP threads for
  [`rxSolve()`](https://nlmixr2.github.io/rxode2/reference/rxSolve.html).

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
V <- cov(dv_mat)

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
#>   Studies: 1 | MC samples: 1000 | Params: 11 | Cores: 1 | Grad: Sens | Restarts: 1
#> +----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+
#> |          |     -2LL |      tcl |      tv1 |      tv2 |       tq |      tka |  prop.sd |   eta.cl |   eta.v1 |   eta.v2 |    eta.q |   eta.ka |
#> +----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+
#> | 0010     | -3635.67 |    5.144 |    11.61 |    26.16 |    10.87 |    1.244 |   0.2006 |  0.09031 |  0.09022 |  0.09008 |  0.09076 |  0.09031 |
#> | 0020     | -3681.10 |     4.99 |    10.73 |    29.22 |    9.787 |     1.08 |   0.1992 |   0.1069 |   0.1038 |  0.09625 |   0.1052 |   0.1067 |
#> | 0030     | -3681.63 |    4.965 |    10.27 |    29.85 |    9.855 |    1.043 |   0.1988 |   0.1043 |   0.1069 |   0.1017 |   0.1083 |  0.09879 |
#> | 0040     | -3681.71 |    4.951 |    10.16 |    29.99 |    9.822 |    1.031 |   0.1986 |   0.1046 |    0.114 |   0.1019 |   0.1068 |   0.0942 |
#> | 0046 ✓   | -3681.73 |    4.952 |    10.12 |    30.03 |    9.815 |    1.027 |   0.1985 |   0.1046 |   0.1156 |    0.101 |   0.1063 |  0.09264 |
#> | 9.7 sec  |          |          |          |          |          |          |          |          |          |          |          |          |
#>   Computing covariance (R method, Sens-Hessian, 7 gradient evaluations)
#> → compress origData in nlmixr2 object, save 1120
print(fit)
#> ── nlmixr² admc ──
#> 
#>           OBJF       AIC       BIC Log-likelihood
#> admc -3681.725 -3659.725 -3589.195       1840.863
#> 
#> ── Time (sec fit$time): ──
#> 
#>   optimize covariance elapsed
#> 1    9.676     10.859  20.535
#> 
#> ── Population Parameters (fit$parFixed or fit$parFixedDf): ──
#> 
#>           Est.      SE   %RSE Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
#> tcl        1.6 0.01659  1.037    4.952 (4.793, 5.115)     33.2            
#> tv1      2.315 0.08736  3.774    10.12 (8.528, 12.01)     35.0            
#> tv2      3.402 0.04064  1.194    30.03 (27.73, 32.52)     32.6            
#> tq       2.284 0.02143 0.9384    9.815 (9.411, 10.24)     33.5            
#> tka     0.0262 0.08201  313.1   1.027 (0.8741, 1.206)     31.2            
#> prop.sd 0.1985                                 0.1985                     
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
