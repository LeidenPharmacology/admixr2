# Control settings for the Gauss-Hermite (GH) quadrature estimator

Creates a control object for `nlmixr2(est = "adgh")`. The GH estimator
integrates model predictions against the random-effects prior \\\eta
\sim N(0, \Omega)\\ using a deterministic tensor-product Gauss-Hermite
quadrature grid. It is unbiased at any IIV magnitude (unlike FO),
noise-free (unlike MC), and much faster than MC for models with up to ~4
etas.

## Usage

``` r
adghControl(
  studies = list(),
  n_nodes = 5L,
  grad = c("analytical", "fd", "cfd", "none"),
  algorithm = "NLOPT_LN_BOBYQA",
  maxeval = 500L,
  ftol_rel = .Machine$double.eps^(1/2),
  print = 10L,
  seed = 12345L,
  cores = 1L,
  grad_h = 1e-04,
  grad_bounds = 5,
  cov_h = 0.001,
  cov_h_outer = .Machine$double.eps^(1/4),
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
  `E`, `V`, `n`, `times`, `ev`, optional `method`).

- n_nodes:

  Number of quadrature nodes per eta dimension (default 5). Total nodes
  = `n_nodes^n_eta`. `n_nodes = 5` achieves near-exact covariance
  moments for IIV SD up to ~0.5; `n_nodes = 7` extends coverage to SD
  ~0.7. For models with \>= 5 etas the node count grows steeply;
  consider reducing `n_nodes` or using a different estimator.

- grad:

  Gradient mode. `"analytical"` (default) uses closed-form contractions
  through the sensitivity equations – cheapest and exact. `"fd"` uses
  forward finite differences; `"cfd"` uses central FD. `"none"` uses
  derivative-free BOBYQA.

- algorithm:

  nloptr algorithm. Automatically coerced to `"NLOPT_LD_LBFGS"` when
  `grad != "none"`.

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

- grad_h:

  Finite-difference step for unpaired struct theta gradient and FD
  Jacobian fallback.

- grad_bounds:

  Box-constraint half-width when using gradients.

- cov_h:

  Inner FD step for the gradient-based Hessian (only used when
  `covMethod = "r"` and `grad != "none"`).

- cov_h_outer:

  Outer step scale for numerical Hessian. Default `eps^(1/4)` (tighter
  than admc's `eps^(1/5)` because the GH surface is noise-free).

- covMethod:

  `"r"` computes covariance via numerical Hessian for structural and
  residual-error parameters only; `"none"` skips it.

- n_restarts:

  Number of optimizer restarts (1 = no multi-start).

- restart_sd:

  SD of random perturbations of initial struct thetas at each restart.

- workers:

  Number of parallel PSOCK/fork workers (default 1 = sequential).

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
  nlmixr2 output tables: `"combined2"` (default) or `"combined1"`.

- returnAdmr:

  If `TRUE`, return a plain list instead of the full nlmixr2 fit object.

- ...:

  Unused arguments (trigger an error).

## Value

An `adghControl` object (a named list).

## See also

[`admControl()`](https://leidenpharmacology.github.io/admixr2/reference/admControl.md),
[`adfoControl()`](https://leidenpharmacology.github.io/admixr2/reference/adfoControl.md),
[`adirmcControl()`](https://leidenpharmacology.github.io/admixr2/reference/adirmcControl.md)

## Examples

``` r
ctl <- adghControl()
ctl$n_nodes
#> [1] 5
ctl$grad
#> [1] "analytical"

# More nodes for large IIV, analytical gradient
ctl2 <- adghControl(n_nodes = 7L, grad = "analytical", maxeval = 300L)

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
  pk_model, admData(), est = "adgh",
  control = adghControl(
    studies = list(study1 = list(E = E, V = V, n = length(ids),
                                 times = times, ev = et(amt = 100)))
  )
)
#>  
#>  
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> === admixr2: Aggregate Data Modeling (GH) ===
#>   Studies: 1 | Params: 5 | Nodes: 5^2=25 | Cores: 1 | Grad: Analytical | Restarts: 1
#> +----------+----------+----------+----------+----------+----------+----------+
#> |          |     -2LL |      tcl |       tv |  prop.sd |   eta.cl |    eta.v |
#> +----------+----------+----------+----------+----------+----------+----------+
#> | 0010     |   972.97 |     5.84 |    34.83 |   0.3185 |  0.08983 |  0.05427 |
#> | 0020     |   727.90 |    7.788 |    37.96 |   0.4167 |   0.2292 |  0.04492 |
#> | 0030     |   726.86 |     8.17 |    38.19 |   0.4231 |   0.2708 |  0.04749 |
#> | 0034 ✓   |   726.85 |     8.15 |    38.21 |   0.4229 |    0.269 |  0.04736 |
#> | 0.8 sec  |          |          |          |          |          |          |
#>   Computing covariance (R method, Analytical-Hessian, 4 gradient evaluations)
#>   Note: covMethod='r' computes covariance for structural and sigma parameters only; omega (IIV) SEs are not computed (matching nlmixr2 FOCEI behavior).
#> → compress origData in nlmixr2 object, save 1160
# }
```
