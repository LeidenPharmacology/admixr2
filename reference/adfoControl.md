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
  algorithm = "NLOPT_LN_BOBYQA",
  maxeval = 500L,
  ftol_rel = .Machine$double.eps^(1/2),
  print = 10L,
  seed = 12345L,
  cores = 1L,
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
  `E`, `V`, `n`, `times`, `ev`, optional `method`).

- grad:

  Gradient mode. `"none"` (default) uses derivative-free BOBYQA;
  `"analytical"` uses the closed-form FO gradient (requires sensitivity
  equations); `"fd"` uses forward finite differences of the full NLL;
  `"cfd"` uses central finite differences for struct theta gradient
  (more accurate than `"fd"`, roughly twice as many NLL evaluations per
  step).

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
  Jacobian.

- grad_bounds:

  Box-constraint half-width when using gradients.

- cov_h_outer:

  Outer step scale for NLL-FD Hessian.

- covMethod:

  `"r"` computes covariance via numerical Hessian; `"none"` skips it.

- n_restarts:

  Number of optimizer restarts (1 = no multi-start).

- restart_sd:

  Standard deviation for random perturbations of initial struct thetas
  at each restart (\> 1).

- workers:

  Number of parallel PSOCK/fork workers for multi-restart (default 1 =
  sequential).

- returnAdmr:

  If `TRUE`, return a plain list instead of the full nlmixr2 fit object.

- ...:

  Unused arguments (trigger an error).

## Value

An `adfoControl` object (a named list).

## See also

[`admControl()`](https://leidenpharmacology.github.io/admixr2/reference/admControl.md),
[`adirmcControl()`](https://leidenpharmacology.github.io/admixr2/reference/adirmcControl.md)
