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
  grad = c("sens", "fd", "none"),
  grad_h = 1e-04,
  cov_h = 0.001,
  cov_h_outer = .Machine$double.eps^(1/5),
  grad_bounds = 5,
  covMethod = c("r,s", "r", "none"),
  cov_n_sim = 10000L,
  n_restarts = 1L,
  restart_sd = 0.5,
  workers = 1L,
  rxControl = NULL,
  calcTables = FALSE,
  compress = TRUE,
  ci = 0.95,
  sigdig = NULL,
  sigdigTable = NULL,
  addProp = c("combined2", "combined1"),
  optExpression = TRUE,
  sumProd = FALSE,
  literalFix = TRUE,
  returnAdmr = FALSE,
  resid_nodes = 81L,
  xtol_rel = .Machine$double.eps^(1/2),
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

  - `v_denom` – `"ml"` (default) or `"unbiased"`, declaring which
    denominator the supplied `V` uses. The likelihood is the exact one
    for `n` iid draws only under the ML (`n`) covariance, which is what
    `cov.wt(method = "ML")` and
    [`datagen()`](https://leidenpharmacology.github.io/admixr2/reference/datagen.md)
    produce. A **published** SD is the unbiased (`n - 1`) SD, so a
    digitised figure gives `V = SD^2` on the `n - 1` scale: declare
    `v_denom = "unbiased"` and admixr2 converts it. Declared per study,
    since a meta-analysis routinely mixes a digitised source with a
    model-derived one and the two need not share a denominator. At
    `n = 60` the factor is 1.7%; it matters more the smaller `n` is, and
    more again for any method that scores the reported covariance
    against its own sampling law.

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
  (central finite differences; forward was removed in 0.4.1), or
  `"none"` (derivative-free). A warning is issued when `"sens"` is
  requested but the sensitivity model is unavailable; the estimator then
  falls back to central finite differences.

- grad_h:

  Step size for finite-difference gradient evaluation during
  optimization (used by `grad = "fd"`). This is the FALLBACK step: the
  step is normally measured per parameter by the Shi (2021) procedure,
  and `grad_h` is what a parameter falls back to when that measurement
  cannot be made (a direction the objective is flat in, or a failed
  noise estimate).

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

  Box-constraint half-width when using gradients: the fit is confined to
  `p0 +/- grad_bounds` on the optimizer scale, which for a log-scale
  parameter is a factor of `exp(grad_bounds)` (~148 at the default 5).
  This bound is admixr2's, not the model's – an unbounded parameter has
  no other – and nloptr reports normal convergence at a box corner, so a
  warning is emitted if an estimate finishes on it.

- covMethod:

  `"r,s"` (the DEFAULT) computes the sandwich `H^-1 J H^-1`; `"r"` the
  numerical Hessian alone, `2H^-1`; `"none"` skips the covariance. A
  study generated from a published model defaults to `"none"` and
  refuses an explicit covariance method because it has no sampling law.
  All three span the structural, residual-error and omega parameters.
  Omega is included because excluding it also biases the STRUCTURAL
  standard errors downward – a theta carrying an eta is correlated with
  that eta's variance. If the weakly-identified omega Cholesky makes the
  Hessian non-positive definite, the structural + residual sub-block is
  reported with a warning.

  `"r,s"` adds a sandwich correction, `H^-1 J H^-1`, on the same
  Hessian. The aggregate objective scores the reported mean and
  covariance as though the subjects behind them were multivariate
  normal; they are not, because the model is nonlinear in the random
  effects. `"r,s"` scores that law from the model instead, on a
  quadrature ensemble rather than on this fit's own MC draws, so the
  reported uncertainty carries no sampling noise of its own. Point
  estimates are untouched, and under correct specification it reduces to
  `"r"` exactly. **It is the default because it is the conservative
  choice, not the aggressive one.** Under correct specification `J = 2H`
  and the sandwich returns what `"r"` returns, so defaulting to it costs
  nothing when the normal-theory assumption holds and corrects the
  standard errors when it does not. Anything it cannot build degrades to
  `"r"` and reports `"r"`, so no fit loses its covariance by asking.
  Pass `covMethod = "r"` for the pre-0.4.1 behaviour.

  Applies to every residual family whose conditional law is independent
  across timepoints, which is all of them except
  [`ar()`](https://rdrr.io/r/stats/ar.html): the conditionally-normal
  set (`add`, `prop`, `pow`, `combined1`, `combined2`), the closed-form
  distributional ones (`lnorm`, `pois`, `binom`, `nbinomMu`, `beta`, and
  [`t()`](https://rdrr.io/r/base/t.html) with `nu > 4`), and the
  transform-both-sides ones (`boxCox`, `yeoJohnson`, `logitNorm`,
  `probitNorm`), whose third and fourth conditional moments come off the
  same quadrature that already gives their mean and variance. Refused,
  and degraded to `"r"`: [`ar()`](https://rdrr.io/r/stats/ar.html),
  because it correlates the residual ACROSS timepoints and the cross
  terms the expansion drops are then real;
  [`t()`](https://rdrr.io/r/base/t.html) with `nu <= 4`, whose kurtosis
  does not exist; and `ordinal()` and same-subject `joint` studies,
  which stack several outputs into one covariance the per-output node
  ensemble does not describe. These four are refusals by construction
  rather than failures, so the fit reports the reason as a message and
  falls back to `"r"`; a sandwich that was attempted and could not be
  built still warns.

  **`"r,s"` is more sensitive to an ill-conditioned Hessian than `"r"`
  is.** `"r"` reports `2H^-1` and inverts `H` once; the sandwich reports
  `H^-1 J H^-1` and inverts it twice, so in a direction the data barely
  identifies any gap between `J` and `2H` is amplified quadratically. A
  residual SD contributing 0.01 variance against 1.7 from
  between-subject variability is such a direction: measured on one 1-cmt
  fixture at `cond(H) = 3.5e5`, the reported residual SE moved by a
  factor of 0.11 and two omega entries by 0.59 and 1.55, while the same
  model and design on a study the residual IS identified in
  (`cond(H) = 247`) reproduced `"r"` to four decimals on every
  parameter. Neither number is a correction there – both methods are
  reporting an unidentified direction, and `"r,s"` is louder about it.
  admixr2 says so: when the Hessian's reciprocal condition number falls
  below `eps^(1/4)` – the point at which squaring the conditioning
  reaches the bound a single inversion is already called singular at –
  the fit records a note naming the parameter that loads most heavily on
  the offending direction. It arrives on `fit$runInfo` and is listed by
  `print(fit)`, which is where `nlmixr2est` routes an estimator's
  warnings. The sandwich is still reported, because the well-determined
  parameters of the same fit are unaffected; check the named parameter's
  relative standard error before reading its `"r,s"` value as a finding.

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

- calcTables, compress, ci, sigdigTable, optExpression, sumProd,
  literalFix:

  Passed to
  [`nlmixr2est::foceiControl()`](https://nlmixr2.github.io/nlmixr2est/reference/foceiControl.html)
  for the table/output machinery.

- sigdig:

  Significant digits asked of the ODE solver, or `NULL` (the default) to
  leave rxode2's own solver tolerances alone. When set, it is passed to
  [`rxode2::rxSolve()`](https://nlmixr2.github.io/rxode2/reference/rxSolve.html)'s
  own `sigdig` argument for every solve the estimator issues – rxode2
  owns the mapping to `atol`/`rtol` and has changed it between releases,
  which is why the digits, not the tolerances, are what travels – and to
  [`nlmixr2est::foceiControl()`](https://nlmixr2.github.io/nlmixr2est/reference/foceiControl.html)
  for the post-fit tables.

  It is a speed lever, and an opt-in one because it is not free. The
  estimators finite-difference the solve with steps of the same order:
  `grad_h` (1e-4), `cov_h` (1e-3) and `cov_h_outer` (~2.5e-3), while
  `sigdig = 4` maps to a relative tolerance of ~1e-4 on current rxode2.
  Differencing a solution whose own noise is 1e-4 with a 1e-4 step
  returns noise, and it surfaces as a moved objective and an indefinite
  covariance Hessian (every `SE` reported `NA`) rather than as an error.
  Most worthwhile where the gradient is fully analytic and nothing
  differences the solve – `adfoControl(grad = "analytical")` measured
  ~4.8x faster at `sigdig = 4` with standard errors unchanged to 4
  significant figures. Elsewhere, compare the objective and the standard
  errors against `NULL` before relying on it. Table formatting is
  unaffected either way: `sigdigTable` defaults to 4 regardless.

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

- xtol_rel:

  Relative parameter tolerance for convergence (default
  `sqrt(.Machine$double.eps)`).

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
#> | 0010     | -3667.69 |    4.896 |    11.82 |    27.71 |    9.353 |    1.208 |   0.1949 |  0.09176 |  0.09044 |  0.09008 |  0.09218 |  0.09068 |
#> | 0020     | -3689.45 |    4.992 |    10.83 |    29.16 |    9.664 |     1.08 |     0.19 |   0.1069 |    0.104 |   0.0955 |   0.1058 |   0.1078 |
#> | 0030     | -3690.00 |    4.967 |     10.4 |     29.7 |     9.75 |    1.047 |   0.1896 |   0.1035 |   0.1027 |  0.09876 |   0.1106 |   0.1041 |
#> | 0040     | -3690.05 |    4.958 |    10.37 |    29.81 |    9.743 |    1.043 |   0.1894 |   0.1033 |   0.1087 |   0.1018 |   0.1092 |  0.09935 |
#> | 0050     | -3690.08 |    4.956 |    10.25 |     29.9 |    9.734 |    1.031 |   0.1894 |   0.1034 |   0.1118 |  0.09989 |   0.1081 |  0.09633 |
#> | 0060     | -3690.08 |    4.957 |    10.26 |    29.88 |    9.736 |    1.033 |   0.1894 |   0.1034 |    0.112 |  0.09973 |   0.1085 |  0.09614 |
#> | 0070     | -3690.08 |    4.957 |    10.26 |    29.88 |    9.736 |    1.033 |   0.1894 |   0.1034 |    0.112 |  0.09973 |   0.1085 |  0.09614 |
#> | 0080     | -3690.08 |    4.957 |    10.26 |    29.88 |    9.736 |    1.033 |   0.1894 |   0.1034 |    0.112 |  0.09973 |   0.1085 |  0.09614 |
#> | 0083 ✓   | -3690.08 |    4.957 |    10.26 |    29.88 |    9.736 |    1.033 |   0.1894 |   0.1034 |    0.112 |  0.09973 |   0.1085 |  0.09614 |
#> | 10.4 sec |          |          |          |          |          |          |          |          |          |          |          |          |
#>   Computing covariance (R method, Sens-Hessian, sandwich, 12 gradient evaluations)
#> → compress origData in nlmixr2 object, save 1160
#>  
#>  
print(fit)
#> ── nlmixr² admc ──
#> 
#>          OBJF      AIC      BIC Log-likelihood
#> admc -3690.08 -3668.08 -3597.55        1845.04
#> 
#> ── Time (sec fit$time): ──
#> 
#>   optimize covariance other elapsed
#> 1   10.361     13.062     0  23.423
#> 
#> ── Population Parameters (fit$parFixed or fit$parFixedDf): ──
#> 
#>            Est.       SE  %RSE Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
#> tcl       1.601  0.01976 1.234    4.957 (4.769, 5.153)    33.01         NaN
#> tv1       2.329   0.1294 5.556    10.26 (7.964, 13.22)    34.42         NaN
#> tv2       3.397  0.05201 1.531    29.88 (26.99, 33.09)    32.38         NaN
#> tq        2.276  0.02769 1.217    9.736 (9.222, 10.28)    33.86         NaN
#> tka     0.03217   0.1200 372.9   1.033 (0.8163, 1.306)    31.77         NaN
#> prop.sd  0.1894 0.003291 1.737 0.1894 (0.1830, 0.1959)                     
#>  
#>   Covariance Type (fit$covMethod): r,s
#>   Some strong fixed parameter correlations exist (fit$cor) :
#>                 cor:tv1,tcl             cor:tv2,tcl              cor:tq,tcl 
#>                  0.301                  -0.483                   0.206   
#>             cor:tka,tcl         cor:prop.sd,tcl       cor:om.eta.cl,tcl 
#>                  0.330                  0.0481                   -0.125   
#>       cor:om.eta.v1,tcl       cor:om.eta.v2,tcl        cor:om.eta.q,tcl 
#>                 -0.243                   -0.183                    0.253   
#>       cor:om.eta.ka,tcl             cor:tv2,tv1              cor:tq,tv1 
#>                  0.203                   -0.838                   0.217   
#>             cor:tka,tv1         cor:prop.sd,tv1       cor:om.eta.cl,tv1 
#>                  0.981                 -0.0353                 -0.00544   
#>       cor:om.eta.v1,tv1       cor:om.eta.v2,tv1        cor:om.eta.q,tv1 
#>                 -0.617                  0.0743                    0.510  
#>       cor:om.eta.ka,tv1              cor:tq,tv2             cor:tka,tv2 
#>                  0.629                  -0.269                   -0.846  
#>         cor:prop.sd,tv2       cor:om.eta.cl,tv2       cor:om.eta.v1,tv2 
#>                 0.0123                   0.0571                    0.543  
#>       cor:om.eta.v2,tv2        cor:om.eta.q,tv2       cor:om.eta.ka,tv2 
#>                 0.0582                   -0.450                  -0.556  
#>              cor:tka,tq          cor:prop.sd,tq        cor:om.eta.cl,tq 
#>                  0.257                 -0.00147                  -0.0496   
#>        cor:om.eta.v1,tq        cor:om.eta.v2,tq         cor:om.eta.q,tq 
#>                 -0.131                   0.0657                    0.326  
#>        cor:om.eta.ka,tq         cor:prop.sd,tka       cor:om.eta.cl,tka 
#>              0.0000123                  -0.0383                 -0.00740   
#>       cor:om.eta.v1,tka       cor:om.eta.v2,tka        cor:om.eta.q,tka 
#>                 -0.605                  0.0691                    0.537  
#>       cor:om.eta.ka,tka   cor:om.eta.cl,prop.sd   cor:om.eta.v1,prop.sd 
#>                  0.618                 -0.0349                  -0.0447   
#>   cor:om.eta.v2,prop.sd    cor:om.eta.q,prop.sd   cor:om.eta.ka,prop.sd 
#>                 -0.183                   -0.159                  -0.0172   
#> cor:om.eta.v1,om.eta.cl cor:om.eta.v2,om.eta.cl  cor:om.eta.q,om.eta.cl 
#>                 0.0270                   -0.136                  -0.0149   
#> cor:om.eta.ka,om.eta.cl cor:om.eta.v2,om.eta.v1  cor:om.eta.q,om.eta.v1 
#>                -0.0203                  -0.0181                   -0.266   
#> cor:om.eta.ka,om.eta.v1  cor:om.eta.q,om.eta.v2 cor:om.eta.ka,om.eta.v2 
#>                 -0.642                  -0.120                   0.0408   
#>  cor:om.eta.ka,om.eta.q 
#>                  0.191   
#>  
#> 
#>   No correlations in between subject variability (BSV) matrix
#>   Full BSV covariance (fit$omega) or correlation (fit$omegaR; diagonals=SDs) 
#>   Distribution stats (mean/skewness/kurtosis/p-value) available in fit$shrink 
#>   Censoring (fit$censInformation): No censoring
#>   Minimization message (fit$message):  
#>     NLOPT_FAILURE: Generic failure code. 
# }
```
