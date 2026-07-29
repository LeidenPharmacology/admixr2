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
  algorithm = NULL,
  maxeval = 5000L,
  ftol_rel = .Machine$double.eps,
  print = 1L,
  omega_expansion = 1,
  seed = 12345L,
  cores = rxode2::rxCores(),
  nDisplayProgress = .Machine$integer.max,
  grad = c("analytical", "none", "fd"),
  kappa_method = c("exact", "linearized", "linearized_gh"),
  kappa_n_nodes = 5L,
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
  resid_nodes = 81L,
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

- outer_iter:

  Maximum inner optimiser iterations per phase.

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

- omega_expansion:

  Inflate proposal Omega by this factor (\>= 1).

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

  Gradient mode for the inner optimiser: `"analytical"` (default,
  closed-form weight-path gradient), `"none"` (derivative-free BOBYQA),
  or `"fd"` (finite differences). Note: `"sens"` and `"cfd"` are not
  available for the IRMC estimator.

- kappa_method:

  Kappa correction method for models with non-mu-referenced struct
  thetas: `"exact"` (default, re-evaluates population prediction
  `f(theta, 0)` via rxSolve at each inner step), `"linearized"`
  (precomputes `J = df/d(theta)` once per outer iteration using
  `f(theta, 0)` as baseline — zero rxSolve per inner step), or
  `"linearized_gh"` (same linear approximation but baseline and Jacobian
  use Gauss-Hermite quadrature `E_GH[f(theta, eta)]` instead of
  `f(theta, 0)` — more accurate baseline at any IIV magnitude, still
  zero rxSolve per inner step).

- kappa_n_nodes:

  Number of GH nodes per eta dimension for
  `kappa_method = "linearized_gh"` (default 5). Total quadrature points
  = `kappa_n_nodes^n_eta`. Ignored for other kappa methods.

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

  Covariance method: `"r"` (numerical Hessian over the structural,
  residual-error and omega parameters) or `"none"`. Omega is included
  because excluding it also biases the STRUCTURAL standard errors
  downward – a theta carrying an eta is correlated with that eta's
  variance. If the weakly-identified omega Cholesky makes the Hessian
  non-positive definite, the structural + residual sub-block is reported
  with a warning.

  All three blocks are reported on the scale the ESTIMATES are printed
  on, as `nlmixr2est` does: structural thetas on the log/optimizer
  scale, residual error as an SD, and omega as the variance/covariance
  entries (named `om.<eta>` and `cov.<eta_i>.<eta_j>`). The omega block
  is rotated by the full Jacobian of Omega with respect to the
  log-Cholesky, which is not diagonal once omega is correlated.

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

  Additional arguments (none allowed; triggers an error).

## Value

An object of class `adirmcControl`.

## Details

Multi-compartment fits (a study `observations` list with several
observed outputs) are **not** supported by `adirmc`; use `est = "admc"`,
`"adfo"`, or `"adgh"` for those. Single-output studies are fit as usual.

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
#> → calculate sensitivities
#> → finding duplicate expressions in admixr2 sensitivity model...
#> → optimizing duplicate expressions in admixr2 sensitivity model...
#>  
#>  
#>  
#>  
#>  
#>  
#> === admixr2: Aggregate Data Modeling (IR-MC) ===
#>   Studies: 1 | MC samples: 500 | Phases: 4 | Iters/phase: 50 | Expansion: 1.00 | Grad: analytic+Sens-Hessian | Restarts: 1
#> +----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+
#> |          |     -2LL |      tcl |      tv1 |      tv2 |       tq |      tka |  prop.sd |   eta.cl |   eta.v1 |   eta.v2 |    eta.q |   eta.ka |
#> +-- Phase 1: Wide (+/-2.00) -------------------------------------------------------------------------------------------------------------------+
#> | 0001     |   728.22 |    3.552 |     6.49 |    36.82 |    8.796 |   0.5216 |   0.1867 |   0.4043 |   0.1828 |    0.665 |  0.02393 |   0.2747 |
#> | 0002     | -1247.61 |    4.966 |    6.662 |    33.45 |    7.927 |   0.6548 |   0.1217 |   0.1298 |  0.02474 |   0.1283 |  0.02614 |   0.2091 |
#> | 0003     | -1194.49 |        5 |    8.026 |    33.06 |    8.569 |   0.7865 |   0.1911 |  0.09852 |   0.1406 |  0.07392 | 0.006252 |   0.1188 |
#> | 0004     | -1238.75 |    4.872 |    7.888 |    31.64 |    8.864 |   0.7894 |   0.1762 |   0.1499 |  0.07673 |   0.1169 | 0.006838 |   0.1372 |
#> | 0005     | -1264.67 |    4.955 |    7.829 |    31.56 |    8.832 |     0.79 |   0.1859 |   0.1141 |  0.07267 |  0.08511 | 0.006873 |   0.1292 |
#> | 0006     | -1266.43 |    4.959 |    8.001 |    31.16 |    8.888 |   0.8065 |   0.1863 |    0.119 |  0.07429 |  0.08569 | 0.006851 |   0.1323 |
#> | 0007     | -1266.44 |    4.959 |    8.004 |    31.15 |    8.885 |   0.8061 |   0.1863 |    0.119 |  0.07429 |   0.0857 | 0.006851 |   0.1322 |
#> | 0008 ✓   | -1266.44 |    4.959 |    8.002 |    31.15 |    8.887 |   0.8061 |   0.1863 |    0.119 |  0.07429 |  0.08571 |  0.00685 |   0.1322 |
#> +-- Phase 2: Focused (+/-1.00) ----------------------------------------------------------------------------------------------------------------+
#> | 0009     | -1266.44 |    4.959 |    8.002 |    31.15 |    8.887 |   0.8061 |   0.1863 |    0.119 |  0.07429 |  0.08571 |  0.00685 |   0.1322 |
#> | 0010     | -1266.44 |    4.959 |    7.991 |    31.15 |    8.886 |   0.8054 |   0.1862 |    0.119 |  0.07433 |  0.08583 | 0.006848 |   0.1319 |
#> | 0011     | -1266.44 |    4.959 |    7.991 |    31.16 |    8.886 |   0.8054 |   0.1862 |    0.119 |  0.07433 |  0.08583 | 0.006848 |   0.1319 |
#> | 0012     | -1266.44 |    4.956 |    7.984 |     31.2 |    8.882 |   0.8044 |   0.1863 |   0.1191 |  0.07447 |  0.08599 | 0.006845 |    0.132 |
#> | 0013     | -1266.44 |    4.957 |    7.983 |     31.2 |     8.88 |   0.8043 |   0.1862 |   0.1191 |  0.07448 |    0.086 | 0.006845 |   0.1319 |
#> | 0014     | -1266.44 |    4.956 |    7.982 |     31.2 |    8.881 |   0.8043 |   0.1862 |   0.1191 |  0.07448 |    0.086 | 0.006845 |   0.1319 |
#> | 0015     | -1266.45 |    4.957 |    7.974 |    31.19 |    8.881 |   0.8038 |   0.1861 |   0.1191 |  0.07453 |   0.0861 | 0.006843 |   0.1317 |
#> | 0016     | -1266.45 |    4.957 |    7.974 |    31.19 |    8.881 |   0.8038 |   0.1861 |   0.1191 |  0.07453 |   0.0861 | 0.006843 |   0.1317 |
#> | 0017     | -1266.22 |    4.923 |    7.869 |    31.59 |    8.843 |   0.7939 |    0.185 |   0.1201 |  0.08131 |  0.09585 | 0.006669 |   0.1294 |
#> | 0018     | -1266.41 |    4.939 |    7.867 |    31.58 |    8.843 |   0.7933 |   0.1854 |   0.1188 |  0.07741 |  0.09027 | 0.006751 |   0.1287 |
#> | 0019     | -1266.47 |    4.937 |    7.862 |    31.54 |    8.859 |   0.7941 |   0.1856 |    0.119 |  0.07746 |  0.09029 |  0.00675 |   0.1289 |
#> | 0020     | -1266.48 |    4.939 |    7.872 |    31.51 |    8.851 |   0.7944 |   0.1858 |   0.1191 |  0.07757 |  0.09032 |  0.00675 |   0.1293 |
#> | 0021     | -1266.48 |    4.938 |    7.874 |    31.52 |    8.848 |   0.7945 |   0.1858 |   0.1191 |  0.07757 |  0.09032 |  0.00675 |   0.1293 |
#> | 0022 ✓   | -1266.48 |    4.938 |    7.874 |    31.52 |    8.849 |   0.7944 |   0.1858 |   0.1191 |  0.07758 |  0.09032 |  0.00675 |   0.1293 |
#> +-- Phase 3: Fine-tuning (+/-0.50) ------------------------------------------------------------------------------------------------------------+
#> | 0023     | -1266.48 |    4.938 |    7.874 |    31.52 |     8.85 |   0.7945 |   0.1858 |   0.1191 |  0.07758 |  0.09032 |  0.00675 |   0.1293 |
#> | 0024 ✓   | -1266.48 |    4.938 |    7.874 |    31.52 |     8.85 |   0.7945 |   0.1858 |   0.1191 |  0.07758 |  0.09032 |  0.00675 |   0.1293 |
#> +-- Phase 4: Precision (+/-0.01) --------------------------------------------------------------------------------------------------------------+
#> | 0025     | -1266.48 |    4.938 |    7.883 |    31.51 |     8.85 |   0.7952 |   0.1856 |    0.119 |  0.07768 |  0.09035 | 0.006749 |   0.1295 |
#> | 0026     | -1266.48 |    4.938 |    7.883 |     31.5 |    8.851 |   0.7953 |   0.1856 |   0.1191 |  0.07768 |  0.09036 | 0.006749 |   0.1295 |
#> | 0027     | -1266.48 |    4.939 |    7.884 |     31.5 |     8.85 |   0.7953 |   0.1856 |   0.1191 |  0.07769 |  0.09036 | 0.006749 |   0.1295 |
#> | 0028     | -1266.48 |    4.939 |    7.884 |     31.5 |    8.851 |   0.7953 |   0.1856 |   0.1191 |  0.07769 |  0.09036 | 0.006749 |   0.1295 |
#> | 0029 ✓   | -1266.48 |    4.939 |    7.885 |     31.5 |    8.851 |   0.7954 |   0.1856 |   0.1191 |  0.07769 |  0.09037 | 0.006749 |   0.1295 |
#> | 1.2 sec  |          |          |          |          |          |          |          |          |          |          |          |          |
#>   Computing covariance (R method, MC NLL, Sens-Hessian, 12 gradient evaluations)
#> → compress origData in nlmixr2 object, save 1160
#>  
#>  
print(fit)
#> ── nlmixr² adirmc ──
#> 
#>             OBJF       AIC       BIC Log-likelihood
#> adirmc -1266.481 -1244.481 -1173.951       633.2406
#> 
#> ── Time (sec fit$time): ──
#> 
#>   optimize covariance elapsed
#> 1    1.182      9.057  10.239
#> 
#> ── Population Parameters (fit$parFixed or fit$parFixedDf): ──
#> 
#>           Est.      SE  %RSE Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
#> tcl      1.597 0.04205 2.633    4.939 (4.548, 5.363)     35.6            
#> tv1      2.065  0.3248 15.73     7.885 (4.172, 14.9)     28.4            
#> tv2       3.45  0.1365 3.957     31.5 (24.11, 41.16)     30.8            
#> tq        2.18  0.1072 4.916    8.851 (7.174, 10.92)     8.23            
#> tka     -0.229  0.2938 128.3  0.7954 (0.4472, 1.415)     37.2            
#> prop.sd 0.1856 0.01504   8.1 0.1856 (0.1562, 0.2151)                     
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
