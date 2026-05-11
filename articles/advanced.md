# Advanced usage

## Full covariance vs variance NLL

Each study’s `V` determines which branch of the NLL is used. `admixr2`
auto-detects the method:

| V structure | Auto-detected `method` | NLL cost |
|----|----|----|
| Non-diagonal matrix | `"cov"` | O(n_t³) Cholesky solve |
| Diagonal matrix | `"var"` | O(n_t) element-wise |
| Plain vector (variances only) | `"var"` | O(n_t) element-wise |

Override by setting `method` explicitly in the study list:

``` r

# Force full covariance NLL even for a diagonal V
study_cov <- c(study, list(method = "cov"))

# Diagonal approximation when you only have marginal variances
study_var <- list(
  E      = E,
  V      = diag(diag(V)),   # drop off-diagonal entries
  n      = n,
  times  = times,
  ev     = rxode2::et(amt = 100),
  method = "var"
)
```

Use `"cov"` when you have a full covariance matrix and expect temporal
correlation to inform parameter estimates. Use `"var"` when only
marginal variances are available (e.g. from a published table of means
and SDs) or when runtime is a priority.

## Gradient modes

[`admControl()`](https://leidenpharmacology.github.io/admixr2/reference/admControl.md)
offers four gradient strategies via the `grad` argument:

| `grad =` | Method | Notes |
|----|----|----|
| `"sens"` | Sensitivity equations (default) | Analytical; requires ODE or `linCmt()` model |
| `"fd"` | Forward finite differences | Falls back to this if sens unavailable |
| `"cfd"` | Central finite differences | More accurate FD; ~2× slower than `"fd"` |
| `"none"` | BOBYQA (derivative-free) | No gradient; useful for debugging or simple models |

``` r

# Analytical gradient (default; fastest for ODE models)
fit_sens <- nlmixr2(pk_model, admData(), est = "admc",
  control = admControl(studies = list(s = study), grad = "sens",
                       n_sim = 5000L, seed = 1L))

# Derivative-free — ignores maxeval; use nloptr stopping criteria instead
fit_bobyqa <- nlmixr2(pk_model, admData(), est = "admc",
  control = admControl(studies = list(s = study), grad = "none",
                       n_sim = 5000L, seed = 1L))
```

`"sens"` is recommended for all standard ODE and `linCmt()` models. If
the sensitivity model is unavailable (e.g. rare ODE features), `admixr2`
falls back to `"fd"` automatically with a warning.

## Mu-referencing and sensitivity equations

**Mu-referencing** is the nlmixr2 convention of expressing each
structural parameter as a fixed effect plus a random effect through a
known back-transformation:

``` r

cl <- exp(tcl + eta.cl)   # mu-referenced: tcl paired with eta.cl
```

`admixr2` uses this pairing to classify parameters into three cases,
each handled differently by `grad = "sens"`:

**Paired (mu-referenced) parameters** — rxode2 augments the ODE system
with sensitivity equations `d(pred)/d(eta_i)` for each random effect.
Because `tcl` and `eta.cl` enter additively on the log scale, the
gradient of the NLL with respect to `tcl` equals the gradient with
respect to `eta.cl`. Both are recovered from a single ODE solve, with no
extra function evaluations.

**Non-mu-referenced parameters** — when a structural parameter and its
random effect are written separately
(e.g. `cl <- exp(tcl) * exp(eta.cl)`), `admixr2` cannot recover
`d(NLL)/d(tcl)` from the sensitivity equations alone. The sensitivity
model is still used for all etas, omega, and sigma; only the unpaired
structural parameters require a targeted CRN finite difference solve. A
warning is issued but `grad = "sens"` is kept.

**Parameters without a random effect** — structural parameters with no
corresponding eta (e.g. `v2 <- exp(tv2)` with IIV dropped from V2) are
always handled by CRN finite differences. The perturbation uses the same
quasi-random seed as the nominal draw so that Monte Carlo noise largely
cancels in the difference.

## IRMC kappa correction

The IRMC estimator draws proposals from an inflated omega. For nonlinear
models, the mean of simulated predictions does not equal the prediction
at eta = 0 — a Jensen’s inequality gap. The kappa correction shifts the
predicted mean to account for this bias:

- **`"first-order"`** (default): baseline is the population prediction
  `f(theta, 0)`. Fast; sufficient for most models.
- **`"second-order"`**: adds a curvature term
  `(1/2) * sum_k(omega_k * d²f/deta_k²)`. More accurate for strongly
  nonlinear models.

``` r

adirmcControl(..., kappa_method = "second-order")
```

## Parallel restarts

Multi-restart fitting guards against local optima. `workers > 1` runs
restarts in parallel; `admixr2` manages the worker pool internally:

``` r

fit_par <- nlmixr2(
  pk_model, admData(), est = "admc",
  control = admControl(
    studies    = list(examplomycin = study),
    n_sim      = 5000L,
    n_restarts = 4L,
    workers    = 4L,      # fork on Unix/macOS; PSOCK cluster on Windows
    cores      = 8L,      # total rxSolve OpenMP threads across workers
    restart_sd = 0.3,     # SD of log-scale perturbation from the starting point
    seed       = 1L
  )
)
```

`cores` is distributed automatically: `floor(total_cores / n_workers)`
threads per worker, with remainder allocated cyclically (e.g. 13 cores /
4 workers → 4, 3, 3, 3). Workers are stopped after the restart phase so
all cores are freed for the covariance Hessian step. Call
[`admStopWorkers()`](https://leidenpharmacology.github.io/admixr2/reference/admStopWorkers.md)
manually if a fit is interrupted before cleanup.

## Quasi-random sampling

The `sampling` argument controls how eta samples are drawn:

| `sampling =` | Method                   |
|--------------|--------------------------|
| `"sobol"`    | Sobol sequence (default) |
| `"halton"`   | Halton sequence          |
| `"torus"`    | Kronecker torus          |
| `"lhs"`      | Latin hypercube          |
| `"rnorm"`    | Plain normal             |

Sobol sequences typically require 2–5× fewer samples than plain normal
draws for equivalent NLL variance. For small `n_sim`, `"lhs"` provides
better coverage than Sobol.

## Parameter uncertainty

`covMethod = "r"` (default) computes a numerical Hessian after
optimisation and reports standard errors in `print(fit)`. A larger
`cov_n_sim` reduces MC noise in the Hessian:

``` r

fit_cov <- nlmixr2(
  pk_model, admData(), est = "admc",
  control = admControl(
    studies     = list(examplomycin = study),
    n_sim       = 5000L,
    cov_n_sim   = 10000L,   # more samples → lower Hessian noise
    cov_h_outer = 2.5e-3,   # outer FD step scale (default: eps^(1/5) ≈ 2.5e-3)
    covMethod   = "r",
    seed        = 1L
  )
)
```

Set `covMethod = "none"` to skip uncertainty estimation during early
model development when you only need the point estimates:

``` r

admControl(..., covMethod = "none")
```

If the Hessian is non-positive-definite (SEs printed as `NA`), increase
`cov_h_outer` or `cov_n_sim`.

## Model comparison: AIC and BIC

Standard information criteria work directly on `admFit` objects. Here we
compare the full model (IIV on all five parameters) against a reduced
model with IIV on CL and V1 only:

``` r

ctl <- admControl(
  studies   = list(examplomycin = study),
  n_sim     = 5000L,
  cov_n_sim   = 10000L,
  maxeval   = 300L,
  seed      = 1L
)

fit_full    <- nlmixr2(pk_model,   admData(), est = "admc", control = ctl)
#> [====|====|====|====|====|====|====|====|====|====] 0:00:07
fit_reduced <- nlmixr2(pk_reduced, admData(), est = "admc", control = ctl)
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00 
#> 
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00 
#> 
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> [====|====|====|====|====|====|====|====|====|====] 0:00:04 
#> 
#> [====|====|====|====|====|====|====|====|====|====] 0:00:08

AIC(fit_full, fit_reduced)
#>             df       AIC
#> fit_full    11 -3659.832
#> fit_reduced  8 -3519.403
BIC(fit_full, fit_reduced)
#>             df       BIC
#> fit_full    11 -3589.302
#> fit_reduced  8 -3468.108
```

Lower AIC/BIC favours the more parsimonious model; a difference \> 10 is
generally considered strong evidence.
