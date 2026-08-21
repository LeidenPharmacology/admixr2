# The SE path, for every kind of ADM / MBMA setup

Written 2026-08-21, from `HANDOFF-model-source-n-REPLY.md` plus the measurements
in `HANDOFF-model-source-n.md`. Estimator-independent: adfo, admc, adgh and
adirmc differ in how they compute `H`, not in any of the logic below.

## The one principle

    Var(theta_hat) = H^-1 . MEAT . H^-1

`H` is the Hessian of the objective actually minimised — the same one
`covMethod = "r"` inverts, unchanged, shared by every case here. What differs is
the MEAT, which is the variance of the score, and therefore depends entirely on
**what is random in each block**.

Three kinds of block, and that is the whole taxonomy:

| block came from | what is random in it | meat contribution |
|---|---|---|
| digitised aggregate DATA | the sufficient statistics `(ybar, vech V)`, from `n_s` real patients | `G_s Omega_s(n_s) G_s'` — the ADF/Browne weight, what `"r,s"` builds today |
| a published MODEL, with reported uncertainty | `theta_src_hat` only | `(sum_j G_j D_j) C_src (sum_j G_j D_j)'`, `D_j = dt_j/dtheta_src` |
| a published model treated as KNOWN, or a simulation truth | **nothing** | **zero** |

The third row is not a degenerate case, it is a real and useful one — see
"Simulation and design work" below.

`n` appears in exactly one place: it sets each source's RELATIVE WEIGHT in the
objective, which is what determines `H` and the `G_s`. It is not a precision
statement. For a lone model source it divides straight out and the estimate does
not depend on it at all (measured: `bwt` = 0.75000 at n = 100 / 400 / 1600 /
6400).

## Setup by setup

### 1. All sources are digitised aggregate data — classic MBMA

Nothing changes. `n_s` is a real sample size, the blocks really are sample
statistics from `n_s` patients, and the existing machinery is correct:

- `covMethod = "r"` — `2H^-1`, right under correct specification.
- `covMethod = "r,s"` — `H^-1 J H^-1` with the ADF/Browne weight, right under
  misspecification too. Prefer it wherever the marginal is non-normal, which is
  most of the time.

This is the case everything in `R/adfweight.R` was built for and it needs no
work.

### 2. All sources are published MODELS — "models as input"

    Var(theta_hat) = sum_m G_m C_src,m G_m' ,      G_m = d theta_hat / d theta_src,m

Requires `C_src` per source. `n` is relative weight only.

Where our model can reproduce a lone source's, `G = I` and `Var = C_src`
exactly. That is the sanity check to ship: **with one model source, the ADM SE
must not come out smaller than that source's own reported SE.** A summary of a
model cannot make us more certain than the analyst who had every patient.

### 3. Mixed — some data sources, some model sources

The realistic meta-analysis, and the reason the meat is defined per block rather
than per fit:

    MEAT = sum_{data s} G_s Omega_s(n_s) G_s'  +  sum_{model m} G_m C_src,m G_m'

The two scale DIFFERENTLY in their own `n`, which is why they cannot share a
formula: a data source's score is an average over `n_s` patients, so its term is
linear in `n_s`; a model source's `Var(s_m)` does not involve `n_m` at all, so
its term is quadratic. Nothing about this is optional — using the data weight on
a model block is the current defect.

### 4. One model source BANDED into J strata

The J strata are **one source**, not J. `C_src` is applied ONCE across the
stacked `D`:

    MEAT_m = (sum_{j=1..J} G_j D_j) C_src (sum_j G_j D_j)'

so the contribution has rank at most `dim(theta_src)` however fine the banding.
Applying `C_src` per stratum would give J independent copies and inflate the
apparent evidence in the anti-conservative direction — raising the resolution
would silently buy confidence. `Sum n_k = n` already prevents the analogous
error in the objective; this is the same rule for the covariance.

### 5. A model source that reports SEs but no correlations

Diagonal `C_src`. Usable, and **flag it**. The cost is exactly zero at the
source's own reference and grows with extrapolation distance — measured
diagonal/full, SEs 0.080 and 0.060 referenced at CRCL 90:

| rho | CRCL 90 | CRCL 60 | CRCL 40 | CRCL 22 |
|---|---|---|---|---|
| -0.9 | 1.00 | 0.82 | 0.75 | 0.73 |
| -0.5 | 1.00 | 0.88 | 0.83 | 0.82 |
| +0.5 | 1.00 | 1.18 | 1.34 | 1.41 |
| +0.9 | 1.00 | 1.42 | 2.23 | 3.14 |

Below 1 means the diagonal is OVER-confident. At the reference
`log(x/xref) = 0` kills the cross term, so no correlation can matter there
whatever it is.

**Mitigation worth taking:** the intercept-slope correlation is largely an
artefact of WHERE the source referenced its model. Re-centre the source's model
at its own covariate median before generating blocks and the two are
near-orthogonal, which makes the diagonal nearly exact.

### 6. A model source that reports NO uncertainty at all

There is no honest SE available from that source. Report `NA` for the affected
directions rather than the naive number, which is wrong by a factor the analyst
chooses: measured 1.69 / 0.84 / 0.42 / 0.21 times the source's own SE at
`n` = 100 / 400 / 1600 / 6400, and the direction is not predictable.

The point estimate is unaffected and should still be reported.

### 7. Parameters the source ASSERTED rather than estimated

A fixed allometric 0.75, a `fix()`ed theta. Those rows and columns of `C_src`
are **zero** — the source claims no uncertainty about them, so they contribute
none. The contribution's rank drops accordingly.

This is the same rule `.admCovCoefThetas()` already applies on the design side,
where banding on an asserted coefficient is refused because it credits evidence
that was never earned. Same principle, applied to the covariance.

### 8. Simulation and design work — the third row of the table

Generating blocks from a KNOWN true model, to evaluate a method or size a study.
`theta_src` is not an estimate, so `C_src = 0`, so the correct
`Var(theta_hat)` is **zero** — and that is right: the blocks are deterministic,
and ADM recovers the truth exactly (measured, 0.75000 at every `n` and every
design).

The SE people want here is a different quantity: *"if this were a real study of
`n` patients with this design, how precise would I be?"* That is exactly what
`2(nh)^-1` answers, and it answers it correctly. So the naive covariance is not
useless — it is a **design** statistic, and should be labelled as one rather
than as evidence.

    naive SE  ->  a DESIGN question: how much would a study this size tell me?
    C_src SE  ->  an EVIDENCE question: how much does this source actually know?

Both are legitimate. Reporting the first when the user asked the second is the
defect.

### 9. The double-counting hazard, which no formula catches

If a trial's published MODEL and its published DATA are both entered, the same
patients contribute twice — once through `C_src` and once through
`Omega_s(n_s)`. The SE comes out too small and nothing in the algebra notices,
because the two blocks look independent. This is a data-entry discipline
problem, not a covariance problem, and belongs in the documentation next to
`n`.

## What the naive covariance IS, so it is not misread

`2(n h)^-1 = (1/n) . 2h^-1`. `n` multiplies the WHOLE matrix by a scalar, so its
correlation structure is invariant in `n` — verified exactly: correlations
identical to 0.000e+00 between `n` = 100 and 1600, elementwise ratio exactly 16
in every cell.

So the set of answers the naive route can produce over every `n` is a
**one-parameter ray through a single shape**. `C_src` is a general PD matrix and
generally sits off that ray, which is why no choice of `n` reaches it, and why
the `n` that would match one parameter (3743) differs from the one that would
match another (448) within the same source.

The shape it carries is real — it describes how the GENERATED DESIGN identifies
the parameters. It simply is not a description of how the source's own data
identified them.

## Implementation consequences

- The bread is shared, so `H` is needed in every case and nothing about it
  changes. `covMethod = "r,s"` becomes **source-aware** rather than gaining a
  sibling: same sandwich, meat chosen per block.
- `D_j = dt_j/dtheta_src` costs `dim(theta_src)` re-aggregations per source.
  Note `datagen()` DOES solve — the reply's "no ODE solve, only re-aggregation"
  is not true of this codebase — but there is no refit, which is the expensive
  and noisy part of differentiating an argmin.
- `C_src` must be on the SOURCE MODEL'S OWN `ini()` scale, since that is what
  gets perturbed: log for a mu-referenced theta, natural for a covariate
  exponent. This needs a checked contract, not a documented convention.
- **`cond(H)` must be checked before differentiating.** A near-singular `H` lets
  the refit wander along the flat direction and returns a large arbitrary number
  instead of a clean infinity — the same failure mode as a marginalised discrete
  covariate on a flat ridge (`.admCovDiscContrast`), caught at inference time
  rather than design time.
- Where the sandwich cannot be formed, the existing convention applies: degrade,
  and report what the covariance IS rather than what was asked for. But a model
  source degrading to `"r"` must WARN loudly, because unlike the data case the
  fallback number is not merely approximate — it is unrelated to the question.
