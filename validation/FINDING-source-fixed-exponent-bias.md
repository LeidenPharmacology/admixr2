# A source that fixes an exponent while omitting a CORRELATED covariate biases the meta-analysis

Measured 2026-08-27/31 on `feature/covariate-quadrature`.

> **This document replaces `FINDING-marginal-covariate-attenuation.md`, whose
> central claim was WRONG.** That draft concluded that a covariate marginalised
> in every source is attenuated ~10% with coverage 0.5-0.7, and that this was
> inherent to identifying a coefficient from the between-source contrast. It is
> not. The attenuation was an artefact of the simulated sources, and it
> disappears when the covariates are uncorrelated. Two intermediate claims were
> also withdrawn along the way -- see "Retractions" at the end. The mistake in
> each case was concluding from a comparison I went looking for, without the
> control that would have refuted it.

## The result

Truth: `CL = 5 (WT/70)^0.75 (CRCL/90)^0.6 exp(0.18 SEX) exp(eta)`, omega 0.05.
Three cohorts at CRCL 38 / 62 / 95, n = 150. Every published source model holds
weight at the conventional `^0.75` and has **no renal term** -- within one
cohort CRCL cannot be told from IIV. The pooled model has the renal term. The
only thing varied below is the WT-CRCL correlation and whether the sources are
FITTED (focei on simulated patients) or written down EXACTLY at the values an
analyst would converge to.

| sources | rho(WT,CRCL) | reps | mean | bias +/- SE | SE/sd | coverage |
|---|---|---|---|---|---|---|
| fitted | 0.45 | 10 | 0.546 | **-0.054** +/- 0.009 | 1.16 | **0.50** |
| exact  | 0.45 |  8 | 0.571 | -0.029 +/- 0.006 | 0.64 | 0.25 |
| **fitted** | **0** | 8 | **0.622** | **+0.022** +/- 0.012 | **1.06** | **1.00** |
| **exact**  | **0** | 8 | **0.601** | **+0.001** +/- 0.007 | 0.58 | 0.88 |

**With uncorrelated covariates the estimator is unbiased and the interval is
correctly sized.** `SE/sd` = 1.06 and coverage 1.00 in the realistic arm. The
covariance machinery was never implicated: `covMethod` resolved to `"r,s"` on
every replicate of every arm.

## The mechanism

`CRCL` correlates with `WT` (rho = 0.45 on the latent scale). Regressing the
truth's log-CL on log-WT alone therefore gives an *induced* exponent

    0.75 + 0.6 * rho * (sd_logCRCL / sd_logWT)
      = 0.75 + 0.6 * 0.45 * (0.217/0.198)  ~=  1.05

but every analyst pinned weight at 0.75. Their published models consequently
misrepresent their own cohorts' CL-WT relationship, and admixr2 -- correctly --
generates `(E, V)` from what they published. The pooled model, which
marginalises over the joint (WT, CRCL) distribution and so DOES capture the
induced relationship, cannot reconcile that with the sources, and the renal
coefficient absorbs the discrepancy.

Note the sign flips: at rho = 0 the fitted arm sits at +0.022, at rho = 0.45 at
-0.054. Two source-side effects of opposite sign (the induced-exponent
mismatch, and how focei apportions the omitted variance between omega and
residual error) partially cancel, which is why intermediate diagnostics looked
inconsistent.

## What this means in practice

**The hazard is real, not a simulation artefact.** Weight and creatinine
clearance are correlated in essentially every PK dataset -- Cockcroft-Gault
contains weight -- and fixing allometry at 0.75 is standard practice. So the
configuration that fails here is the ordinary one.

- **Prefer a source that ESTIMATED its allometric exponent** over one that
  fixed it, when you intend to add a covariate correlated with weight. In the
  main experiment one analyst (Ito) did estimate both exponents; that is the
  model to trust for this purpose.
- **Banding remains the strongest remedy** where available: coverage 0.93 from
  the 720-replicate campaign, one banded source as good as three.
- **A marginalised covariate is NOT intrinsically problematic.** That was the
  wrong lesson. Where the sources are internally consistent, it recovers
  unbiased with correct coverage.

## Also established, and unaffected

- **Structural disagreement between sources is essentially free.** A
  one-compartment source fitted to steady-state peak and trough only, pooled
  with two-compartment sources, recovers every PK parameter: CL 4.97/5.0,
  V 50.08/50.0, Q 7.98/8.0, Vp 71.4/70, coverage 0.89-1.00. A source publishing
  `V = 57` (its one-compartment volume absorbing the peripheral compartment)
  still yields a pooled `V` of 50.1.
- **A covariate banded in its sources is textbook**: `bsex` unbiased at
  coverage 0.95 in every arm.
- **Not quadrature**: `cov_nodes` 7 -> 15 (121 -> 529 and 196 -> 900 design
  points) moves the estimate by 1.8e-11.

## Retractions, in order

1. *"~10% attenuation is inherent to a between-source contrast."* Refuted by
   the rho = 0 arms.
2. *"Use a source that reports AT a covariate value."* The pinned arm changed
   the TRUTH (no within-cohort spread), so it diagnosed a mechanism rather than
   licensing a recipe. The spread is real; declaring it away misrepresents the
   source.
3. *"More sources fixes it, O(1/k)."* Read off k = 3 and k = 6 (p = 0.053);
   k = 9 did not replicate (p = 0.43 vs k = 3).

## A caveat on the exact-source arms

`SE/sd` is 0.58-0.64 there, but that comparison is not a coverage assessment:
the exact generator computes `(E, V)` deterministically and adds **no
patient-level sampling noise**, which is precisely what admixr2's `V/n` term
describes, while its only randomness -- the cohort draw -- is something admixr2
conditions on as data. The FITTED arms are the valid coverage assessments.

## Drivers

`study-heterogeneous-designs.R`, `study-source-count-arm.R` (env vars
`NSRC`/`NPER`/`CRCLSD`/`RHO`), `study-exact-source-arm.R` (env vars
`RHO`/`CSV`), `study-pinned-covariate-arm.R`. Per-replicate results in
`sscov_reps2.csv`, `sscov_control.csv`, `sscov_nodes.csv`, `arm_pinned.csv`,
`exp_6src.csv`, `exp_9src.csv`, `exp_narrow.csv`, `exp_exact.csv`,
`exp_exact_rho0.csv`, `exp_fit_rho0.csv`.
