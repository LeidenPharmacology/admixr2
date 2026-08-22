# Between-source disagreement lands on the variance parameters

Measured 2026-08-22. Answers "are all the SEs correct under MBMA" with: the SE
machinery is correct, and the model is missing a term.

## What is exact

Every quantity with an analytically known answer checks out. These are
deterministic -- a departure would be a bug, not sampling noise -- and each
costs a handful of fits rather than hours of replication.

| check | result |
|---|---|
| lone model source | `SE = sqrt(C_src)` exactly |
| k identical sources | x1.000 / 0.707 / 0.577 / 0.500 at k = 1..4 |
| invariance in `n` | exact (0.03536 at n = 75, 150, 600) |
| subgroup-reported source (`by =`) | `= C_src`, ratio 1.0000 |
| meat additive across different covariate centres | 2.9e-13, bread bit-identical |
| **mixed model + data** | `Var(mixed) = Var(data) + Var(model)` to 8.0e-15 |
| Jacobian `dV11/d(add.err)` | 0.160000 across FIVE decades of step, `= 2*sigma` |
| pooling lands on the average | deviations 9e-5 .. 2e-3 |
| coverage `tcl` / `tv` / `om`, sources agreeing | 0.950 / 0.973 / 0.970 |

## What breaks, and why it is not an SE problem

With TWO sources whose drawn `theta_src` differ, `add.err` covers at 0.533 with
a +11.5% bias and SE/SD 0.29. Varying one family of parameters at a time
isolates the cause completely:

| what varies across sources | `add.err` bias | SE/SD | coverage |
|---|---|---|---|
| the residual only | +0.00027 | 0.964 | 0.930 |
| **the means only** | **+0.00849** | 0.436 | 0.570 |
| everything | +0.00867 | 0.380 | 0.540 |

**Varying only the MEANS reproduces 98% of the bias.** Two sources with
identical designs report different `E`; the model predicts one profile, so
nothing in it can explain a between-study difference in the mean, and the only
free knob that widens `V` enough to reconcile them is the residual error. Omega
is hit too (bias +0.0171, SE/SD 0.40), and the mean parameters' own scatter
inflates beyond first-order prediction (`tcl` SE/SD 0.66).

A first-order `G C_src G'` propagation cannot see any of this: absorption is
driven by `|theta_1 - theta_2|`, which is SECOND order in the source draws.

### Four mechanisms eliminated by measurement, not argument

- **the variance formula** -- every exact check above passes;
- **the Jacobian `D`** -- stable to six figures across five decades of FD step,
  and equal to its analytic value `2*sigma`;
- **optimizer convergence noise** -- five different starting values give
  bit-identical estimates and the same objective;
- **finite-sample / non-normal `theta_hat`** -- `tcl` covers 0.950 with SE/SD
  1.03 in the very same fits.

### One separate, smaller effect: SAMPLING DESIGN

Separating a constant additive residual from a time-varying omega needs an
early sample. On `c(1, 3, 8, 24)` the residual-only arm gives SE/SD 0.780 and
coverage 0.820; on log-spaced times from 0.5 h it gives 0.964 and 0.930, and is
insensitive to whether there are 4 or 10 of them. Design, not defect -- but
worth saying in the docs, because it is invisible.

## What this means

The SEs are right for what the model says. The model has no BETWEEN-STUDY
variance component, so when published sources disagree by more than their
reported `C_src` -- the normal case, and the stated premise of the covariates
vignette ("published models disagree, and reconciling them is not free") --
that disagreement has nowhere to go and lands on the residual and omega.

This is the `tau^2` item already on record as unimplemented. No refinement of
the sandwich substitutes for it: the term is absent from the model, not
mis-estimated.

Three possible responses, increasing in cost:

1. **Detect it.** A Cochran's-Q-style check on whether sources disagree by more
   than `C_src` explains, and warn. Cheap, and it is what classical
   meta-analysis has always done. admixr2 is currently a meta-analysis package
   with no heterogeneity diagnostic.
2. **Report the residual and omega SEs as CONDITIONAL** on homogeneity.
3. **Model it** -- a between-study `tau^2`. Correct, and much larger.

Not chosen here; it is a scope decision.

## Reproducing

`scratchpad/{pool_exact,exact_suite,absorb,djac,optnoise,cov_exactSE}.R`.
`absorb.R` takes `MODE` in {mean, resid, all}, `NT` (timepoint count) and
`NREP`; it is the one that isolates the cause.
