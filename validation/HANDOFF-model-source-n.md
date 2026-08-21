# What is `n` for a source that contributed a MODEL?

Handoff back, for theory. Concerns fix 4 of `HANDOFF-model-source-bands.md`.
Fixes 1–3 are implemented (`b67eb8d`, `e3e582d`); this is the part that needs
deciding before it can be.

## Settled empirically — do not re-derive

- **Fix 3 was not the default.** The objective moved 975 units in J because
  banding ran the *pooled* branch, which cuts equiprobable bins and evaluates at
  a representative point — a midpoint rule, `O(1/J)`. The Gauss–Hermite branch
  was gated on `cor` having been supplied, so any study declaring no dependence
  fell to the slow rule. Independence is a known latent structure, so it now
  takes the fast branch: **749 units of J-dependence → 0.2**, flat from J = 9.
  Both routes head to the same limit.
- **Point estimates are J-invariant** at every resolution and both branches
  (0.750 exactly, J = 3 to 25), confirming the original finding.
- Banding preserves total weight: `n_k = w_k·n`, `Σ n_k = n`. So the Hessian
  does not grow with J and `covMethod = "r"` should not shrink with J. **The
  J-inflation fix 4 describes appears to be prevented by the weighting already**
  — which, if right, means fix 4 is not about J at all.

## The open question

For a **data** source, `−2LL = n(log|V| + tr(V⁻¹V_obs) + r'V⁻¹r)` is the
normal/Wishart likelihood of sufficient statistics from `n` subjects. `n` is a
sample size and precision follows from it.

For a **model** source, `(E, V)` are exact functions of `θ_src`. There is no
sampling distribution — the block sits at its floor (`0.000e+00`), as the
original handoff observed. So the criterion is a **weighted discrepancy, not a
likelihood**, and `n` is the weight on it. Asking what standard error it yields
is asking what a weighted least-squares criterion yields, which is whatever the
weight says.

That makes fix 4 a **competition**, not a gap. `n` and `Var(θ̂_src)` are two
statements of the same thing, and using both double-counts. Which should own it?

1. **`n` owns it.** Read `n` as the sample size of the data the source model was
   fitted to. Then fix 4 is redundant for a single source, and what is actually
   needed is guidance on choosing `n`, plus a check.
2. **`Var(θ̂_src)` owns it.** Then `n` is a pure weight and must not also carry
   precision — which changes what `n` means depending on the source type, inside
   one objective.

## Why this matters more than the rank algebra

`n` is user-supplied and SE falls like `1/√n`. For digitised data `n` is
checkable against the paper. **For a model input there is nothing anchoring it**,
so anyone can buy arbitrary confidence by typing a larger number. That hazard is
not in the original handoff and looks more consequential than the rank point.

Two things `n` structurally cannot express, which is where fix 4 would earn its
place regardless:

- **Per-parameter precision.** A source may pin CL tightly and its covariate
  exponent loosely. One scalar cannot say that; a covariance can.
- **Correlation between the source's parameters**, which propagates into ours.

## Worth considering

The cleanest formulation may be that a published model is a **posterior summary,
not pseudo-data**: it contributes a term like `(θ − θ_src)' Var(θ̂_src)⁻¹
(θ − θ_src)` rather than a discrepancy against generated `(E, V)` at a notional
`n`. Then `n` disappears for model sources, the double-count cannot arise, and
per-parameter precision comes for free.

The costs to think through: it only works when the source's model is nested in
(or mappable onto) ours — the generated-data route works for *any* structural
mismatch, which is much of the point of ADM. And mixing the two source types in
one objective needs the relative weighting to be defensible.

## What would settle the empirical half

Compare, for one source and one parameter:

- the SE the **source itself** reported from individual data (the yardstick — a
  summary of a model cannot make us more certain than the analyst who had every
  patient);
- the ADM SE at `n` = the source's sample size, across J and across
  `covMethod`, and for `admc` as well as `adgh`.

A ratio near 1 says `n` is doing the job; well below says it is not.

### Result

Source: 400 patients, individual data, fitted with FOCEI. Its own SE for the
covariate exponent is **0.0558**. That model was then banded and re-fitted
through ADM at `n` = 400.

| method | SE | vs the source's own |
|---|---|---|
| `adgh`, `covMethod = "r"` | 0.0509 | **0.91×** |
| `adgh`, `covMethod = "r,s"` | 0.0746 | **1.34×** |
| `admc`, `covMethod = "r"` | 0.0509–0.0511 | 0.91× |

**Invariant in J** — identical to three decimals at J = 5, 9, 25 and 50 — which
confirms that `Σ n_k = n` prevents the J-inflation fix 4 describes. The rank
concern does not materialise through the weighting.

**So `n` is doing the job**, to within 9%. With `n` set to the source's actual
sample size, `covMethod = "r"` lands just inside the source's own SE — slightly
over-confident, as expected from ignoring the source's parameter uncertainty
entirely, but not by an order of magnitude. The sandwich is 34% *conservative*
by the same yardstick.

That points at option 1: **`n` owns the precision, and fix 4 is not needed to
rescue it for a single source.** What remains genuinely unaddressed is the
per-parameter shape, the correlation between a source's parameters, and the
unanchored-`n` hazard above.

**admc needs no special setting.** n_sim = 2000 and 8000 give 0.0511 and 0.0509,
so the banded studies are not sampling-limited; the cost is J × n_sim solves,
which makes `adgh` the better fit for heavily banded sources.

### One more reason the branch mattered

On the pooled branch the *estimate* is not J-invariant either: 0.6810 at J = 5,
9, 15 and 100, 0.6803 at 25, and **0.6968 at J = 50**. That is the basin-jumping
the original handoff warned about, reproduced. On the Gauss–Hermite branch the
estimate is 0.681000 at every J tested. The pooled branch also had not converged
at J = 100 — still 51 units from the GH value, having come down from 452 — and
approaches it non-monotonically.
