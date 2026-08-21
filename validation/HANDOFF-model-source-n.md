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
- **On the Gauss–Hermite branch the OBJECTIVE is J-invariant too** — 0.03 units
  across J = 5 to 100 — so OFV, AIC, BIC and likelihood ratios are comparable
  across resolutions there. The J-dependence is a property of the banding RULE,
  not of the estimator or of ADM.
- **Point estimates are J-invariant on that branch** (0.681000 at every J).
  They are *not* on the pooled branch: 0.6810 at J = 5, 9, 15, 100 and 0.6968 at
  J = 50 — the basin-jumping the original handoff cited, reproduced and
  localised.
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


## Open, and the largest remaining gap in the banding path

**When does the slow rule actually run?** Two triggers, and the second is the
common one:

1. an opaque user `joint` sampler — nothing there can be conditioned;
2. **any declared discrete covariate**, even when the covariate being banded is
   continuous and independent of it.

Trigger 2 is most real covariate models — sex, genotype, formulation — and it
costs **515 units of J-dependence across J = 5 to 50, still moving**, against
0.03 for the same study without the discrete covariate.

Relaxing it was attempted and **reverted**, because it produced a result that is
not understood rather than one that is wrong in an obvious way:

- the machinery admits discrete margins cleanly — `.admCovQuantile` maps a
  latent normal onto declared levels, and the generated strata carried the right
  level probability (0.550 against a declared 0.55) in *every* band, with the
  banded covariate an exact point per stratum as that branch intends;
- J-dependence fell from 515 units to **0.0**;
- but the recovered coefficient for the discrete covariate moved from **0.150**
  (its truth, which the pooled rule returns exactly at every J) to **−0.052**,
  stably, at every J.

**That reading was wrong, and the follow-up settles it.** Profiling the
objective in that coefficient, everything else refitted at each value:

```
bsex fixed at:   -0.10  -0.06     0    0.05   0.10   0.15   0.20
exact (GH):      0.004  0.065  0.003  0.003  0.002  0.000  0.010   spread  0.07
banded (pooled):18.264 12.911  6.677  3.033  0.797  0.000  0.866   spread 18.26
```

Under exact conditioning the surface is **flat** — it *is* non-identifiability.
A deterministic optimiser on a flat ridge stops in the same place every time,
which is why it looked stable. Under banding the same coefficient has 18 units
of curvature and is cleanly identified, with its minimum at the truth.

Two consequences:

- **The collapse is not involved.** Disabling both `.admCovCollapse` and
  `.admJointCollapse` moves the estimate by 1e-4 (−0.0595 → −0.0594). Worth
  recording because it was the natural suspicion after the recent collapse work.
- **Routing discrete covariates through exact conditioning would silently
  destroy their identifiability.** Exact conditioning removes the within-band
  spread of the *banded* covariate, and that spread is what separates the
  discrete coefficient's contribution to `V` from `omega`; with the level
  proportions identical in every stratum there is no between-stratum contrast to
  fall back on either.

### RESOLVED. The pool was the bug, in two separate ways

Both halves are now implemented and measured.

**1. Discreteness never was the obstacle — the SAMPLER was.** A discrete margin
cannot ride a quadrature grid when the sampler mixes the uniforms before mapping
them to levels, because fixing the input uniform does not fix the output level.
But a margin latently INDEPENDENT of the others has `chol(R)[, j] = e_j`, so the
copula leaves `z_j = qnorm(u_j)` and the level is monotone in its own uniform
after all. `.admCovDiscExact()` records those margins; `.admCovGrid()` crosses
them with the Gauss-Hermite grid over the rest, and `.admCovStrata()` enumerates
them as strata when stratified and as exact level nodes when not.

That is the case the handoff called "most real covariate models" — one declared
sex, genotype or formulation:

| | J-dependence, J = 5 to 50 |
|---|---|
| before | **515 units, still rising** |
| after | **0.009 units** |

`sum(n_k) = 400.0` at every J, and the weight-carrying estimate (`bwt`) is
J-invariant to four decimals.

**2. The pool was not just slower — it was SHARED, and that manufactured
evidence.** An equal-weight pool is a deterministic function of `cov_dist`, which
is the property common random numbers need. It also means `datagen()` and the
fit draw the SAME rows, so whatever that one finite ensemble happens to contain
is generated into `(E, V)` and read straight back out as if it were population
structure.

That is where the −0.052 came from, and the reading in the previous section was
still not right. Replacing SEX by its *exact declared marginal* while leaving
the band's weight spread untouched — same `(E, V)`, different design — collapses
the profile:

```
bsex fixed at:        -0.10      0    0.05    0.10   0.15   0.20   spread
pooled, as shipped:  18.264  6.677   3.033   0.797  0.000  0.866   18.264
pooled + exact SEX:   0.001  0.002   0.001   0.019  0.000  0.008    0.019
```

So the pooled route's 18.3 units of curvature were **the fit reading its own
sampling noise**. It was not carrying information the exact route lost; it was
hiding the absence of information. Two candidate mechanisms were measured and
ruled out first: the per-band level proportion is flat to 4e-4 and does not move
with pool size, and the pool has no within-band covariate association (±2e-4,
against ±4e-2 for an ordinary pseudo-random pool of the same size).

**3. So a marginalised discrete covariate is genuinely not identified**, and the
remedy is a design change, not an integration one. Its effect shows only through
the mixture it induces — level probabilities shift `E`, between-level spread adds
to `V` — which is what a random effect on the same parameter does, and with the
same level distribution in every band there is no between-band contrast either.

Stratifying on it fixes it completely, and costs one study per level with no
quadrature nodes:

| design | `bsex` (truth 0.150) | across J = 5 to 50 |
|---|---|---|
| sex marginalised | −0.0595 | invariant to 1e-4 |
| sex **stratified** | **0.1516** | invariant to 1e-4 |

`stratify = TRUE` already produces the second, since it bands every covariate
whose coefficient the source's own model estimated. For an explicit `stratify`,
`.admCovDiscContrast()` warns and names the remedy — the point being that the
failure is silent otherwise: a deterministic optimizer on a flat ridge stops in
the same place every run, so −0.059 reads as a converged answer.

**Refuted, do not resurrect:** that exact point conditioning is what loses the
coefficient. It is not — the clean product design *with* full within-band weight
spread is equally flat (0.019). Conditioning point-vs-interval is a question
about which estimand a banded source reports; it has nothing to do with this.
