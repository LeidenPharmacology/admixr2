# Banding a source that contributed a MODEL, not data

Handoff for admixr2. Concerns `stratify` / `covStrata()` in `R/covariate.R`.
Evidence: a 4,000-replicate simulation (`kleb_jres_mc.R`, base R, outside the
package) plus an independent adversarial review that derived the failure modes
and verified them numerically. Nothing here refutes the construction — point
estimates are sound. What is wrong is the *inference*, and the fixes are local.

## What is already right — do not change it

- **`n_k = w_k * n`, deliberately unrounded.** Correct, and the comment saying so
  should stay. The weights must sum to `n` at every resolution.
- **`stratify = TRUE` derived from the source model's `allCovs`.** Correct
  principle, and the error when the model reads none of the declared covariates
  is exactly right.
- **Every stratum gets its own generated observation.** This is the property that
  makes the fabricated null contrast unreachable through this path. Keep it.
- **Point estimates.** Confirmed unbiased at every resolution and every weighting:
  J=1 0.75186, J=3 0.75102, J=9 0.75101 against a truth of 0.7500, over 4,000
  replicates. The Gaussian discrepancy has a strict unique minimum at
  `(ybar_j, V_j)`, so a model source's block sits at its floor to `0.000e+00` at
  any J. No estimation change is needed.

## Fix 1 — bands must be bounded by the source's REPORTED covariate range

**The bug.** `.admCovStrata(cov_dist, stratify, n_nodes)` cuts strata from the
analyst's `cov_dist` over its full support. A published model is then evaluated
in covariate bands where that study enrolled nobody, and its coefficient is
credited as evidence there. The inflation is exact:

```
information_claimed / information_earned  =  var_assumed / var_enrolled
```

Verified to 3 decimals: a source that enrolled +/- 1 SD is credited **3.43x**
what it earned; +/- 0.5 SD, **12.38x**. This is not bias — estimates stay
correct — it is **false confidence**, and it only bites once a second source
disagrees, which is when it matters.

**The fix.** Accept a reported range per stratified covariate, truncate
`cov_dist` to it, renormalise `w_k` over the truncated support, keep
`sum(n_k) = n`. Suggested surface, matching the existing style:

```r
list(..., stratify = TRUE, cov_range = list(WT = c(52, impossible_max)))
```

Publications routinely report min-max or median/IQR, so this is available. When
`cov_range` is absent, warn rather than silently using full support — the
default is the leaky case.

## Fix 2 — refuse to band a covariate whose coefficient the source FIXED

**The bug.** `allCovs` reports which covariates a model *reads*, not which
coefficients it *estimated*. A model containing `clwt <- fix(0.75)` — the
allometric convention, so this is common — reads `WT`, passes the current
check, and gets banded. But it carries **zero** evidence about the covariate:
`information_earned = 0` and the ratio above is **unbounded** at every J > 1.

**The fix.** For each covariate in `keep`, locate the theta(s) that reach it and
refuse if every one is fixed. `ui$iniDf$fix` carries this. Error text should say
the source asserted the coefficient rather than estimating it, so there is no
contrast to extract, and the study is generated marginal instead — parallel to
the existing "reads none of the covariates" error.

A partially-fixed set (one covariate estimated, another fixed) should band on the
estimated one only.

## Fix 3 — `strata_nodes` is a CONVERGENCE parameter, not a modelling choice

**The bug.** It defaults to `5L` and reads as a user preference. It is not. The
`J -> infinity` limit is the well-defined object (`N * E_a[l]`); intermediate J
is an approximation to it, and under misspecification the answer can jump
between basins — a verified counterexample flips from `beta = 1.57` at J=4 to
`beta = 1.6e-05` at J=5.

**The fix.** Document it as a quadrature resolution to be driven to convergence,
not chosen; raise the default; and ideally offer a convergence check that refits
at `2 * strata_nodes` and reports the movement, in the manner of the existing
node-count guidance. The joint collapse (see `HANDOFF-joint-collapse.md`) is what
makes a high count affordable.

**Consequence to document:** the objective is J-dependent — measured at **441
units** across J = 1 to 500 with the model correct. So OFV, AIC, BIC and any LRT
are comparable **only at fixed J**. Two fits at different `strata_nodes` cannot
be compared. This deserves a hard warning, and arguably a stored `strata_nodes`
on the fit so a comparison across differing values can be refused outright.

## Fix 4 — standard errors for a model source

**The bug.** All J blocks derived from one published model are deterministic
functions of that model's parameter vector, so the sandwich meat has rank at most
`dim(theta_src)` however many blocks there are. Treating them as independent
contributions overstates the information. The point estimate is unaffected;
`covMethod = "r"` is invalid for such a source at every J.

**The fix.** A published model comes *with* its parameter uncertainty, and that
is the right meat: use the source's reported `Var(theta_hat_src)` for its blocks
rather than a block-independence sandwich. This also answers the standing
question of what a notional-N standard error is an SE *of*. Interacts with the
`covMethod = "r,s"` work — see `algorithm/adf/HANDOFF-INFERENCE.md`.

If the source's covariance is not supplied, the honest outcome is to refuse a
standard error for those parameters rather than report an over-confident one.

## Priority

1 and 2 are correctness-of-inference and cheap. 3 is documentation plus a default
change plus a comparison guard. 4 is the largest and can follow.

## What the simulation could NOT test, so do not assume it

- The leak of Fix 1 was invisible to `kleb_jres_mc.R` because the model source's
  covariate law *was* the analyst's law by construction. It comes from the review's
  derivation, verified numerically there but not reproduced here.
- A single source cannot test a weighting claim at all: scaling every `n_k` of one
  source by a constant leaves the argmin exactly invariant. Any future test needs
  two disagreeing sources.
- No PK-scale demonstration of any of this exists yet.
