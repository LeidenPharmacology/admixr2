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
and
[`adfoControl()`](https://leidenpharmacology.github.io/admixr2/reference/adfoControl.md)
offer gradient strategies via the `grad` argument.
[`adghControl()`](https://leidenpharmacology.github.io/admixr2/reference/adghControl.md)
offers `"analytical"`, `"fd"`, `"cfd"`, and `"none"`.
[`adirmcControl()`](https://leidenpharmacology.github.io/admixr2/reference/adirmcControl.md)
offers `"analytical"`, `"none"`, and `"fd"` only (no `"sens"` or
`"cfd"`).

**[`admControl()`](https://leidenpharmacology.github.io/admixr2/reference/admControl.md)
gradient modes:**

| `grad =` | Method | Notes |
|----|----|----|
| `"sens"` | Sensitivity equations (default) | Analytical; requires ODE or `linCmt()` model |
| `"fd"` | Forward finite differences | Falls back to this if sens unavailable |
| `"cfd"` | Central finite differences | More accurate FD; ~2× slower than `"fd"` |
| `"none"` | BOBYQA (derivative-free) | No gradient; useful for debugging or simple models |

**[`adfoControl()`](https://leidenpharmacology.github.io/admixr2/reference/adfoControl.md)
gradient modes:**

| `grad =` | Method | Notes |
|----|----|----|
| `"none"` (default) | BOBYQA | Derivative-free; robust starting point for FO |
| `"analytical"` | Chain rule through V_pred | Omega/sigma analytical; struct thetas FD only |
| `"fd"` | Forward FD of full NLL | All parameters; `n_p + 1` NLL evals per step |
| `"cfd"` | Central FD of full NLL | More accurate; `2 n_p` NLL evals per step |

**[`adghControl()`](https://leidenpharmacology.github.io/admixr2/reference/adghControl.md)
gradient modes:**

| `grad =` | Method | Notes |
|----|----|----|
| `"analytical"` (default) | Closed-form contractions through sensitivity equations | Exact and cheapest; one batched sensitivity solve per study over the node grid |
| `"fd"` | Forward FD of full NLL | All parameters |
| `"cfd"` | Central FD of full NLL | More accurate; ~2× slower than `"fd"` |
| `"none"` | BOBYQA (derivative-free) | No gradient |

Because the GH objective is **noise-free** (deterministic quadrature, no
MC draws), its analytical gradient is exact and the resulting Hessian is
well-conditioned — `"analytical"` is the recommended default for `adgh`.

For `adfo`, all gradient modes still use the sensitivity model inside
each NLL evaluation to compute J in a single rxSolve — `grad` controls
only how the *optimizer* gradient is formed.

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

The IRMC estimator shifts IS weights by moving the proposal mean when a
structural parameter changes. This works directly for mu-referenced
parameters (e.g. `cl <- exp(tcl + eta.cl)`) because the weight shift is
a closed-form function of the back-transformed ratio.

For **non-mu-referenced parameters** — structural thetas without a
paired random effect (e.g. `ka <- exp(tka)` with no `eta.ka`) — changing
`tka` during inner optimisation alters the population mean prediction
`f(theta, 0)` without a corresponding IS weight path. The kappa
correction accounts for this by adding
`kappa = f(theta_cand, 0) - f(theta_outer, 0)` to the predicted mean,
keeping the inner NLL anchored to the true prediction at the candidate
point. Kappa is zero when all struct thetas are mu-referenced.

Two methods control how `f(theta_cand, 0)` is evaluated during inner
optimisation:

- **`"exact"`** (default): re-evaluates `f(theta_cand, 0)` via a single
  rxSolve at each inner NLL evaluation. Exact; costs one extra rxSolve
  per inner step.
- **`"linearized"`**: precomputes `J = df/d(theta)` once per outer
  iteration via a small FD batch. Approximates
  `kappa_fn(theta_cand) ≈ f0 + J %*% (theta_cand - theta0)` — pure
  arithmetic per inner step, zero extra rxSolve calls.

`"exact"` is accurate but adds an rxSolve to every inner NLL call.
`"linearized"` is faster and the approximation is good when inner steps
stay small relative to the outer box constraint — which is typical for
converged phases. Prefer `"linearized"` for complex ODE models where
each rxSolve is expensive.

``` r

adirmcControl(..., kappa_method = "exact")       # default
adirmcControl(..., kappa_method = "linearized")  # faster for complex ODE models
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
    workers    = 4L,      # background worker processes (mirai daemons)
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
optimisation and reports standard errors in `print(fit)`. **Important:**
these SEs are computed for structural and residual-error parameters
only; omega (IIV) SEs are not computed (matching nlmixr2 FOCEI
behavior). A larger `cov_n_sim` reduces MC noise in the Hessian:

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

Standard information criteria work directly on `admFit` objects. **AIC
and BIC are comparable only across estimators that evaluate the same
likelihood.** `admc`, `adgh`, and `adirmc` all target the exact
aggregate MVN likelihood (MC, quadrature, and importance-weighted
estimates of the same integral), so their information criteria are
mutually comparable. `adfo` evaluates a *linearised* likelihood that
differs from these, so it must not be compared against the other three.

Here we compare the full model (IIV on all five parameters) against a
reduced model with IIV on CL and V1 only:

``` r

ctl <- admControl(
  studies   = list(examplomycin = study),
  n_sim     = 5000L,
  cov_n_sim   = 10000L,
  maxeval   = 300L,
  seed      = 1L
)

fit_full    <- nlmixr2(pk_model,   admData(), est = "admc", control = ctl)
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

AIC(fit_full, fit_reduced)
#>             df       AIC
#> fit_full    11 -3668.835
#> fit_reduced  8 -3528.356
BIC(fit_full, fit_reduced)
#>             df       BIC
#> fit_full    11 -3598.305
#> fit_reduced  8 -3477.061
```

Lower AIC/BIC favours the more parsimonious model; a difference \> 10 is
generally considered strong evidence.
