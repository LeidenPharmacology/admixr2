# Control settings for the FO (First-Order) estimator

Creates a control object for `nlmixr2(est = "adfo")`. The FO estimator
linearises model predictions at \\\eta = 0\\: it is faster than the MC
estimator but less accurate for models with large IIV or strongly
non-linear individual predictions.

## Usage

``` r
adfoControl(
  studies = list(),
  grad = c("analytical", "none", "fd"),
  algorithm = NULL,
  maxeval = 500L,
  ftol_rel = .Machine$double.eps^(1/2),
  print = 10L,
  seed = 12345L,
  cores = rxode2::rxCores(),
  nDisplayProgress = .Machine$integer.max,
  grad_h = 1e-04,
  grad_bounds = 5,
  cov_h = 0.001,
  cov_h_outer = .Machine$double.eps^(1/5),
  covMethod = c("r,s", "r", "none"),
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
  ...
)
```

## Arguments

- studies:

  Named list of study specifications (same format as
  [`admControl()`](https://leidenpharmacology.github.io/admixr2/reference/admControl.md):
  `E`, `V`, `n`, `times`, `ev`, optional `method`; or an `observations`
  list for multi-compartment fits – see
  [`admControl()`](https://leidenpharmacology.github.io/admixr2/reference/admControl.md)).

- grad:

  Gradient mode. `"analytical"` (default) uses the closed-form FO
  gradient with LBFGS; `"none"` uses derivative-free BOBYQA; `"fd"` uses
  central finite differences of the full NLL. Forward differencing was
  removed in 0.4.1 – it was 10^2 to 10^4 times less accurate than a
  central difference at every site measured, and the one solve per
  parameter it saved did not pay for a gradient the optimizer struggles
  to descend.

  The default was `"none"` up to 0.4.0, because the structural thetas
  were finite-differenced through the whole NLL and the resulting
  gradient was too noisy for a quasi-Newton step to pay off. They are
  now differentiated analytically from a second-order sensitivity model
  (relative error ~1e-7 against a central difference, where the
  finite-difference pass reached 1e-2), so LBFGS on the exact gradient
  is the better default. A model that cannot build that sensitivity
  model falls back to the finite-difference gradient automatically, and
  `grad = "none"` remains available.

- algorithm:

  nloptr algorithm, or `NULL` (default) to pick the default that matches
  `grad`: `"NLOPT_LD_LBFGS"` with a gradient, `"NLOPT_LN_BOBYQA"` when
  `grad = "none"`. Any algorithm reported by
  [`nloptr::nloptr.print.options()`](https://astamm.github.io/nloptr/reference/nloptr.print.options.html)
  is accepted. An explicit algorithm is reconciled with `grad`: when
  `grad = "none"` a gradient-based algorithm (`NLOPT_LD_*` /
  `NLOPT_GD_*`) falls back to `"NLOPT_LN_BOBYQA"`; when a gradient is
  requested a derivative-free algorithm (`NLOPT_LN_*` / `NLOPT_GN_*`)
  turns the gradient off. Both emit a message.

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
  [`rxSolve()`](https://nlmixr2.github.io/rxode2/reference/rxSolve.html).
  Defaults to
  [`rxode2::rxCores()`](https://nlmixr2.github.io/rxode2/reference/getRxThreads.html).
  When `workers > 1` it is a *total* budget, split across the workers.

- nDisplayProgress:

  Passed to
  [`rxSolve()`](https://nlmixr2.github.io/rxode2/reference/rxSolve.html):
  show the solver's text progress bar only once a single solve exceeds
  this many subjects. The default (`.Machine$integer.max`) keeps it off
  for clean script/vignette output; lower it (e.g. `1000L`) to see
  progress during long fits.

- grad_h:

  Finite-difference step for unpaired struct theta gradient and FD
  Jacobian.

- grad_bounds:

  Box-constraint half-width when using gradients: the fit is confined to
  `p0 +/- grad_bounds` on the optimizer scale, which for a log-scale
  parameter is a factor of `exp(grad_bounds)` (~148 at the default 5).
  This bound is admixr2's, not the model's – an unbounded parameter has
  no other – and nloptr reports normal convergence at a box corner, so a
  warning is emitted if an estimate finishes on it.

- cov_h:

  Inner FD step for the gradient-based Hessian (only used when
  `covMethod = "r"` and `grad != "none"`). Default 1e-3.

- cov_h_outer:

  Outer step scale for NLL-FD Hessian.

- covMethod:

  `"r,s"` (the DEFAULT) computes the sandwich `H^-1 J H^-1`; `"r"` the
  numerical Hessian alone, `2H^-1`; `"none"` skips the covariance. All
  three span the structural, residual-error and omega parameters. Omega
  is included because excluding it also biases the STRUCTURAL standard
  errors downward – a theta carrying an eta is correlated with that
  eta's variance. If the weakly-identified omega Cholesky makes the
  Hessian non-positive definite, the structural + residual sub-block is
  reported with a warning.

  `"r,s"` adds a sandwich correction, `H^-1 J H^-1`, on the same
  Hessian. For FO this does more than correct kurtosis:
  `V = J Omega J' + Sigma` is the covariance of an exactly normal
  individual law, so the reported standard errors otherwise answer to
  the linearisation rather than to the model. The correction scores the
  FO fit against the model's true nonlinear law, built post-fit on a
  quadrature ensemble, and so absorbs part of the linearisation error as
  well. Point estimates are untouched. **Transform-both-sides
  endpoints.** `adfo` composes the residual by a second-order expansion
  about the linearised moments, because FO carries no node ensemble to
  compose over – that is what the method is. `adgh` and `admc` compose
  exactly at their nodes/draws, so an `adfo` fit of a `boxCox`,
  `yeoJohnson`, `logitNorm` or `probitNorm` endpoint differs from theirs
  by the expansion's truncation: roughly 0.3% in `V` at moderate
  between-subject variability, rising to ~3% for a tightly-bounded
  `logit`/`probit` at high variability. That is a property of the
  estimator, not a discrepancy.

  **It is the default because it is the conservative choice, not the
  aggressive one.** Under correct specification `J = 2H` and the
  sandwich returns what `"r"` returns, so defaulting to it costs nothing
  when the normal-theory assumption holds and corrects the standard
  errors when it does not. Anything it cannot build degrades to `"r"`
  and reports `"r"`, so no fit loses its covariance by asking. Pass
  `covMethod = "r"` for the pre-0.4.1 behaviour.

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

  **An adfo standard error describes scatter, not accuracy.** FO
  linearises the model at eta = 0, and on a non-additive residual (or a
  saturating endpoint, or a large omega) the resulting point estimates
  carry a bias of several standard errors – measured 5-20 SE, giving 0%
  coverage for a nominal 95% interval even where the SE itself matches
  the sampling SD. Use `adgh` or `admc` when the uncertainty matters.

- n_restarts:

  Number of optimizer restarts (1 = no multi-start).

- restart_sd:

  Standard deviation for random perturbations of initial struct thetas
  at each restart (\> 1).

- workers:

  Number of parallel workers (mirai daemons) for multi-restart (default
  1 = sequential). Requires the `mirai` package.

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
  `"combined1"` (SD form). Has no effect on admixr2's own estimation.

- returnAdmr:

  If `TRUE`, return a plain list instead of the full nlmixr2 fit object.

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

  Unused arguments (trigger an error).

## Value

An `adfoControl` object (a named list).

## Installing memuse

[`rxode2::rxSolve()`](https://nlmixr2.github.io/rxode2/reference/rxSolve.html)
estimates free RAM on every call. When the `memuse` package is not
installed its fallback ends up shelling out to `vm_stat`, a macOS-only
command, so on Windows and Linux every solve spawns a process that can
only fail. Because the FO estimator issues many small solves, this
overhead is measurable (roughly 17% of an FO gradient). Installing
`memuse` makes the fallback unreachable:

    install.packages("memuse")

## See also

[`admControl()`](https://leidenpharmacology.github.io/admixr2/reference/admControl.md),
[`adirmcControl()`](https://leidenpharmacology.github.io/admixr2/reference/adirmcControl.md)

## Examples

``` r
# Inspect defaults
ctl <- adfoControl()
ctl$grad
#> [1] "analytical"
ctl$maxeval
#> [1] 500

# Analytical gradient, more evaluations
ctl2 <- adfoControl(grad = "analytical", maxeval = 1000L)

# \donttest{
library(rxode2)
#> rxode2 5.1.6 using 2 threads (see ?getRxThreads)
#>   no cache: create with `rxCreateCache()`
library(nlmixr2)
#> ── Attaching packages ───────────────────────────────────────── nlmixr2 7.0.1 ──
#> ★ lotri        1.0.4      ★ nlmixr2est   7.0.2 
#> ★ nlmixr2data  2.0.10     ★ nlmixr2extra 5.2.0 
#> ★ nlmixr2save  0.2.0      ★ nlmixr2plot  5.1.0 
#> ── Optional Packages Not Installed ──────────────────────────── nlmixr2 7.0.1 ──
#> ✖ babelmixr2     ✖ nlmixr2targets
#> ✖ FME     ✖ nonmem2rx
#> ✖ ggPMX     ✖ pmxNODE
#> ✖ monolix2rx     ✖ PopED
#> ✖ nlmixr2auto     ✖ posologyr
#> ✖ nlmixr2autoinit     ✖ shinyMixR
#> ✖ nlmixr2lib     ✖ xpose.nlmixr2
#> ✖ nlmixr2rpt     
#> ── Conflicts ───────────────────────────────────────────── nlmixr2conflicts() ──
#> ✖ nlmixr2est::boxCox()     masks rxode2::boxCox()
#> ✖ nlmixr2est::yeoJohnson() masks rxode2::yeoJohnson()

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
  pk_model, admData(), est = "adfo",
  control = adfoControl(
    studies = list(study1 = list(E = E, V = V, n = length(ids),
                                 times = times, ev = et(amt = 100))),
    maxeval = 100L
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
#> → calculate sensitivities
#> → calculate sensitivities
#> → finding duplicate expressions in admixr2 sensitivity model...
#> → optimizing duplicate expressions in admixr2 sensitivity model...
#>  
#>  
#>  
#>  
#>  
#>  
#> === admixr2: Aggregate Data Modeling (FO) ===
#>   Obs units: 1 | Params: 5 | Cores: 2 | Grad: Analytical | Restarts: 1
#> +----------+----------+----------+----------+----------+----------+----------+
#> |          |     -2LL |      tcl |       tv |  prop.sd |   eta.cl |    eta.v |
#> +----------+----------+----------+----------+----------+----------+----------+
#> | 0010     |  1768.15 |    4.967 |    29.88 |   0.2587 |   0.0888 |  0.04603 |
#> | 0020     |   862.47 |    6.391 |    37.74 |   0.3864 |  0.08003 |   0.0422 |
#> | 0029 ✓   |   861.90 |    6.384 |    38.03 |     0.39 |  0.08051 |  0.04074 |
#> | 0.6 sec  |          |          |          |          |          |          |
#>   Computing covariance (R method, Analytical-Hessian, sandwich, 6 gradient evaluations)
#> → compress origData in nlmixr2 object, save 1160
#>  
#>  
print(fit)
#> ── nlmixr² adfo ──
#> 
#>          OBJF      AIC      BIC Log-likelihood
#> adfo 861.8956 871.8956 903.9548      -430.9478
#> 
#> ── Time (sec fit$time): ──
#> 
#>         optimize covariance other elapsed other
#> elapsed    0.551      0.331     0   0.882 6.172
#> 
#> ── Population Parameters (fit$parFixed or fit$parFixedDf): ──
#> 
#>           Est.       SE   %RSE Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
#> tcl      1.854  0.01961  1.058    6.384 (6.143, 6.634)    28.95         NaN
#> tv       3.638  0.01689 0.4641    38.03 (36.80, 39.31)    20.39         NaN
#> prop.sd 0.3900 0.009106  2.335 0.3900 (0.3721, 0.4078)                     
#>  
#>   Covariance Type (fit$covMethod): r,s
#>   No correlations in between subject variability (BSV) matrix
#>   Full BSV covariance (fit$omega) or correlation (fit$omegaR; diagonals=SDs) 
#>   Distribution stats (mean/skewness/kurtosis/p-value) available in fit$shrink 
#>   Censoring (fit$censInformation): No censoring
#>   Minimization message (fit$message):  
#>     NLOPT_XTOL_REACHED: Optimization stopped because xtol_rel or xtol_abs (above) was reached. 
# }
```
