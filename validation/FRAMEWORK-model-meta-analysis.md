# Meta-analysis of published models with covariate effects

A framework, consolidating what this session established. Each claim is marked
**[proved]**, **[measured]**, **[literature]** or **[open]**. Nothing is asserted
that is only plausible.

---

## 1. The input, and the two currencies

A publication gives you one of two things, and they are not interchangeable:

| currency | what you receive | what it supports |
|---|---|---|
| **DATA** | aggregate `(E, V, n)` per arm/time — digitised figures, summary tables | an exact likelihood of sufficient statistics |
| **MODEL** | a parameter table + RSEs + structure + design + baseline covariate table | a binding-function likelihood |

Most literature reviews yield a mixture. The framework scores each study in the
currency it was actually reported in, rather than forcing one shape onto the other
(digitising a model back into pseudo-data, or reducing data to a fitted vector).

**Everything downstream depends on which currency a study is in.** In particular,
the covariate question has *opposite* answers in the two cases (§3).

---

## 2. The target

### 2.1 Conditional, not marginal **[proved]**

For a source with covariate law `P_s`, two candidate targets:

```
(M)  D( marginal_s(y) || marginal_phi(y) )
(C)  E_{a~P_s} [ D( p_s(y|a) || p_phi(y|a) ) ]
```

By the data-processing inequality for f-divergences (Csiszar 1963; Ali & Silvey
1966; Liese & Vajda 2006), marginalisation is a deterministic kernel, so

```
(M)  <=  (C)      pointwise in phi
```

with equality **iff** marginalisation is a *sufficient statistic* for the pair —
which fails precisely when the covariate discriminates the two conditionals.
So (M) is a **lower bound**: minimising it minimises a bound, and its minimum is
attained on a flat manifold. Measured symptom: a covariate coefficient of 0.196
against a truth of 0.75 (**-74%**), with omega inflated ~40% **[measured]**.

Two corollaries worth stating in any write-up:

- `(C)` **equals the joint KL** over `(a, y)` when `P_s` is used on both sides —
  the chain-rule term `D(P(a)||Q(a))` is identically zero. So it is not an ad hoc
  average of divergences; it *is* a divergence. **[proved]**
- `(C)` is the expected **log score of the conditional forecaster**, hence a
  strictly proper scoring rule. **[proved]**

### 2.2 The framing that will land with a pharmacometrics audience

> Both objectives are strictly proper — **for different forecasting tasks**.
> - marginal is proper for: *predict for a randomly drawn patient whose weight you do not know.*
> - conditional is proper for: *predict for a patient whose weight you do know.*
>
> A unified model used for **dosing** is always deployed in the second task.

Present the DPI as the derivation; present this as the argument.

### 2.3 Direction: forward KL, `D(source || unified)` **[literature]**

Three converging reasons:

1. Forward KL is the **population limit of maximum likelihood on data simulated
   from the source** — so the synthesis agrees with the fit you would have got
   had the individual data been available. Reverse KL has no such reading.
2. Over an unrestricted family, `argmin sum_s w_s D(p_s||q)` is the **linear
   opinion pool**; `argmin sum_s w_s D(q||p_s)` is the **logarithmic** pool,
   which is zero wherever *any* source is zero — one source vetoes all others.
   That is the opposite of synthesis.
3. Forward KL over-disperses, reverse under-disperses. Under-stated predictive
   variance means over-confident dosing intervals — the dangerous direction.

---

## 3. THE PER-SOURCE RULE (the core of the framework)

**Score each source at the granularity it actually resolves.** This is the single
most important design decision and it is where the naive approach fails.

| source | contributes via |
|---|---|
| fitted the covariate | **conditional** — node quadrature over its `P_s` |
| **did not** fit it | **marginal**, over its own reported `P_s` |
| reported aggregate DATA | the ordinary aggregate likelihood |

**Why the middle row is not optional.** A model that omits covariate `a` is **not**
the model `beta_a = 0`. It is that source's conditional *already marginalised* over
`a`'s distribution in its own population. Forcing a conditional on it by
substituting zero injects a false conditional. Measured cost of getting this wrong:
a covariate coefficient biased **-29%**, and a second coefficient **-43%**, *at zero
noise* — because the largest study in the set had no covariate term and was read as
asserting no effect **[measured]**.

Absence of evidence is not evidence of absence, and no pharmacometric paper states
this for model omission (the ISoP covariate guidance says it only about
non-significant hypothesis tests) **[literature]**.

---

## 4. The objective

```
-2LL(phi)  =   sum_{s in DATA}   sum_k  n_sk * [ log|V| + tr(V^-1 V_obs) + r' V^-1 r ]
             + sum_{s in MODEL}  ( theta_s - g_s(phi) )' Sigma_s^-1 ( theta_s - g_s(phi) )
```

`g_s(phi)` is the **binding function**: what study `s` would have published, under
*its own* model form, design and population, if the truth were `phi`. Computed by
simulating the unified model at `phi` over `s`'s design and refitting `s`'s form.

**Names, so this is positioned rather than reinvented [literature]:**

- the second term alone = **indirect inference**, Wald/parameter-matching form
  (Gourieroux, Monfort & Renault 1993), with a fixed exogenous design per source
- multiple heterogeneous auxiliary models = **Composite Indirect Inference**
  (Gourieroux & Monfort 2018, *Econometrics and Statistics* 7:30-45)
- the whole construction, published outside pharmacometrics =
  **Generalized Model Aggregation**, Rahmandad, Jalali & Paynabar (2017),
  *PLoS ONE* 12(4):e0175111 — cite it, do not re-derive it
- for a biostatistics audience, the citable name is **"the indirect method"**,
  Jiang & Turnbull (2004), *Statistical Science* 19(2):239-263 — framed exactly
  for the case where the auxiliary estimate is *not* consistent for the target
- the sum of the two terms = **composite likelihood** (Varin, Reid & Firth 2011);
  closest precedent Diao et al. (2026) for hybrid IPD + aggregate
- the multi-source structure = **multiparameter evidence synthesis**
  (Ades & Sutton 2006)

**Contrast to state explicitly:** two-stage multivariate meta-analysis
(`mvmeta`) requires a *common parametrisation* across studies. The binding
function is exactly what removes that requirement. That is the cleanest one-line
statement of what is new.

### 4.1 Structural facts that simplify our case

- **`g_s` is the identity** whenever a source's form equals the unified form.
  Verified to `0.00e+00` **[measured]**. Only structurally reduced sources need an
  inner solve, which is what makes the nesting affordable.
- **Our `Sigma` is block-diagonal, and that is the easy case.** In CII all
  auxiliary models are fitted to the *same* data, so the optimal weight carries
  non-zero cross-blocks. Our sources are independent samples, so the criterion
  above is **exactly** CII-optimal, not an approximation **[literature]**.
- **One currency per study.** If a study offers both data and a model, use the
  data and discard its model. This sidesteps double counting entirely, which the
  sandwich alone would not fix (it corrects the variance, not the pull on the
  estimate).

---

## 5. The covariate integral

### 5.1 Node quadrature, and when it is valid

The node route is a genuine likelihood **iff each node carries its own `(E_k, V_k)`**.

| data | prediction | verdict |
|---|---|---|
| one pooled `(E,V)`, reused per node | per node | **invalid**, -59% **[measured]** |
| per node | per node | **exact**, 0.7500 vs truth **[measured]** |
| per stratum | marginalised within stratum | valid |
| per stratum | stratum mean plugged in | +17% at K=2 **[measured]** |

In the MODEL currency you *generate* `(E_k, V_k)` at each node, so you are always
in the valid regime. The invalid regime requires a pooled observation, which
models-as-input never supplies. **[proved]**

The distinguishing question is therefore *where the data come from*, never *how
many nodes there are*.

### 5.2 Mechanics

- **Point nodes**, not interval strata. A point node is the direct quadrature of
  `E_a[...]`; an interval averages within itself and reintroduces a slice of the
  marginal collapse.
- **Two nested integrals with different jobs**: the inner (over eta) builds the
  conditional-on-`a` aggregate moments — the "data" at that node, carrying IIV and
  residual only; the outer (over `a`) is the expectation the target names.
- **Weights**: `n_k = w_k * n_s`, so `sum_k n_k = n_s`. Total information is the
  study's real size, not `K * n_s`. A per-stratum list with `n_k = n_s` at every
  node overstates information by `K`x.
- **Node count is a convergence diagnostic**, not a tuning choice — same target
  throughout, so it must settle. Failure to settle means a non-identified
  parameter, and should be reported as such.
- **Nodes on the scale the covariate is normal on.** Moment-matching a normal to
  a lognormal weight puts outer nodes at negative weight, where `WT^clwt` does not
  exist — the PAGE failure. Build on `log(WT)` and exponentiate.
- **Truncate to the source's enrolled range**, then rebuild the rule on the
  truncated law (Golub-Welsch; use modified moments, raw moments are
  ill-conditioned — Gautschi 1970). A node outside the studied range asks the
  source's power law to extrapolate and then treats the answer as evidence. The
  covariate exponent is the parameter most sensitive to lever-arm at the extremes,
  so this is where it bites hardest. **[open]** — no literature on this.

### 5.3 Dependent covariates: use the uniform scale

A product grid assumes independence **on the covariate scale**. It does not on the
**uniform** scale: the inverse Rosenblatt transform maps any vine to independent
uniforms *by construction*.

```
tensor GH grid on z   ->   u = Phi(z)   ->   a = R^-1(u)  via the copula
```

Weights stay a clean product because the `u` genuinely are independent. Exact, not
an approximation. `.admCovGrid()` already builds nodes as
`spec$quantile(pnorm(gh$x))` per margin — i.e. it is already a product grid on the
uniform scale — so the change is to push the whole uniform vector through
`joint(u)` rather than applying margins separately. That makes `adgh`'s blanket
refusal of dependent covariates unnecessary.

**Caveat to measure, not assume:** `f(R^-1(u))` can be far less smooth in `u` than
`f(a)` is in `a`, especially under tail dependence, so Gauss rules may need more
nodes than the independent case. **[open]**

Beyond two covariates a tensor grid is `n^d`; sparse grids carry **negative
weights**, which for a nonnegative divergence integrand can return a negative
objective — guard explicitly, or use QMC.

**And the binding constraint is not numerical.** Publications report marginal
mean +/- SD per covariate and essentially never the joint. `P_s` is *not identified*
from the publication. No quadrature repairs that: treat the dependence as a
sensitivity analysis with an assumed copula and say so.

---

## 6. Weighting and uncertainty

- **Do not read SEs off the fit.** With manufactured nodes the objective is an
  M-estimator, not a likelihood; the Hessian reports whatever precision you asked
  for. Measured coverage of a nominal 95% interval: **0.58** **[measured]**.
- **The composite objective needs the Godambe sandwich `H^-1 J H^-1`**, with the
  **study** as the clustering unit: `J = sum_s U_s U_s'` from *total* per-study
  scores, so an (A)+(B) overlap is absorbed automatically **[literature]**.
  admixr2's `2*H^-1` is not the right estimator for a combined objective.
- **`(1 + 1/S)` inflation** on the indirect-inference variance, where `S` is the
  simulated population size relative to `n_s` — not the number of studies. With
  deterministic quadrature `S -> Inf` and the factor vanishes **[literature]**.
- **Bootstrap `Sigma_s` rather than trusting published RSEs.** We already simulate
  and refit each source, so bootstrapping that refit gives the *full* sampling
  covariance under our own model — no publication needs to report correlations.
  Published RSEs are inverse-Hessian SEs from a fit we know to be misspecified, so
  they are probably **too small**, which systematically over-weights the MODEL
  term **[open]**.
  - **Prerequisite, learned the hard way:** the binding function's dimension must
    be supported by the source's design. Matching 7 parameters from a 2-timepoint
    design gave `Sigma_B` with a 0.96 off-diagonal correlation and SEs of 14.7
    against a true 0.05 **[measured]**. If a source publishes more parameters than
    its design identifies, drop them from the match — do not regularise past it.
  - Fallback when bootstrapping is not possible: diag + a single shared
    correlation (Riley, Thompson & Abrams 2008).
- **AIC/BIC and LRT are invalid** on a composite objective without the
  Chandler-Bate adjustment (the composite LR statistic is a weighted sum of
  chi-sq_1) **[literature]**.
- **Source weighting is genuinely open.** No data => no AIC weights => no
  performance weights. Every scheme in climate, epidemiology and pharmacometric
  model averaging needs observations. `n_s`? study quality? equal? **[open]**

---

## 7. Diagnostics — all of them, every time

1. **Can the design see the covariate?** Spread of the predicted mean across
   nodes. mg/kg dosing with CL and V both linear in weight cancels algebraically —
   measured 0.0% spread, so that source says nothing about weight at any `n`
   **[measured]**. Run this *before* fitting.
2. **Node convergence.** Sweep the node count; a parameter still moving is not
   identified.
3. **`g_s` identity check** for any source whose form matches the unified one. It
   must return `phi` unchanged; anything else is an implementation fault.
4. **Each source projected alone**, compared with the joint fit. Solo answers that
   disagree badly mean one structure does not hold — misspecification, not
   uncertainty.
5. **Over-identification `Q`**, `df = sum_s p_s - dim(phi)`. Report it as a
   **diagnostic, never a gate**: the asymptotic chi-sq is known to be poor, more
   matched coefficients *lowers* power (4 -> 69.6%, 27 -> 21.6%), and with large
   `n_s` it rejects almost any approximate model **[literature]**. Simulate the
   null instead — the machinery is already there.
6. **Conflict diagnostic**: fit DATA-only, MODEL-only, test agreement at `phi`
   (Presanis et al. 2013). Cheap, and it is what a reviewer will ask for.

---

## 8. Heterogeneity

`theta_s ~ N(g_s(phi) + u_s, Sigma_s)`, `u_s ~ N(0, tau^2)`.

**Not established anywhere for a binding function** — five independent literature
sweeps found no formulation. That makes it a contribution, and also unsupported.

Two honest caveats:

- **`u_s` *is* the misspecification of `g_s`, re-described as random.** `tau^2`
  buys variance inflation, not diagnosis; whether the discrepancy is exchangeable
  across studies is unfalsifiable from the data.
- **Do not add `tau^2` on top of a diagonalised `Sigma_s`.** `tau^2` would then
  absorb both genuine heterogeneity and the corruption introduced by diagonalising,
  with no way to separate them. Bootstrap `Sigma_s` first.

Estimator: Paule-Mandel or REML with **modified** HKSJ intervals (Rover, Knapp &
Friede 2015) — plain HKSJ fails when source precisions vary widely, which is our
regime. Restrict the between-source covariance (diagonal or single shared
correlation); a full matrix from 5-10 sources is not estimable.

Serious alternative worth considering instead: Armstrong & Kolesar (2021) — a
bounded *deterministic* discrepancy set giving honest intervals, rather than an
estimated random-effect distribution.

---

## 9. The limitation that must be stated

**The estimand is not invariant to the study set.** Under misspecification the
pseudo-true value depends on which moments are matched (Hall, Inoue, Nason & Rossi
2012), so adding or dropping a source changes what `phi_hat` converges to. This is
not a robustness weakness to be minimised in discussion — it is a property of the
estimator and belongs in the methods section **[literature]**.

Related: no bounded influence. One badly wrong source enters linearly with weight
`(J'Sigma^-1 J)^-1 J' Sigma^-1`. The literature's answers are *select*
(information criteria for which moments to match) or *downweight*
(bounded-influence II, Genton & Ronchetti 2003). There is no free robustness.

---

## 10. Implementation status in admixr2

**Available today, unchanged:** point-node quadrature. A study with a point `cov`
and `n = w_k * n_s` reproduces the node objective bit-for-bit in all four
estimators; `n` is already a `double` linear multiplier in every C++ kernel.
`.admCovNodesFor()` / `.admCovGrid()` are already better node generators than the
removed `admBuildQuadrature` (no `statmod`, log-scale lognormals, discrete and
arbitrary-quantile specs).

**Needs fixing:**

1. `R/utils.R:49` — `as.integer(s$n)` **truncates** fractional `n_k`, silently
   zeroing tail nodes in AIC/BIC. The only genuine bug blocking this workflow.
2. `.admRefuseNodeStudies()`'s message advises "Do NOT give a stratum a point
   `cov`" — right for a *published* stratum spanning a range, wrong for a *chosen
   evaluation node* whose within-node distribution is degenerate. It also
   contradicts `datagen()`'s own message, which is correct. Keep the refusal
   (`weight`/`cov_method` must stay rejected); fix the advice.
3. The `R/covariate.R` header enumerates two data shapes and never considers
   models-as-input.
4. `.admCovGrid()` does not filter `rho`/`Sigma`/`joint` from `names(cov_dist)` —
   latent bug: `adgh` with a bare `rho` builds a zero-row grid.
5. `man/plot.admFit.Rd` documents a `"covariate"` panel `plot.admFit()` no longer
   accepts. Stale today, independent of any of this.

**Wanted:** an exported `covNodes(cov_dist, n_nodes)` beside `covDraw()` — ~10
lines over `.admCovDistCanon()` + `.admCovGrid()`; the copula route of §5.3; and
analytic `dg_s/dphi` by the implicit function theorem on the inner first-order
condition, since finite-differencing through a nested Nelder-Mead leaves the outer
gradient at the inner solver's tolerance and makes Gauss-Newton oscillate
**[measured]**.

**Do not restore:** `gl`, `taylor` (taylor at `h = sqrt(3)*sigma` *is* GH-3 to
2e-16, every other `h` a worse rule), or `.admCovExpandNodes()` — that is the
pooled-data path the -59% condemns.

---

## 11. Novelty, honestly

- **Ancestor**: Valitalo (2021) — the aggregate likelihood with `n` propagated.
  Contains the word "covariate" **zero** times (verified in full text). The
  covariate extension is unclaimed in this line.
- **Competitor**: Suzuki et al. (Nov 2025) "M-cubed" — 19 published vancomycin
  models -> virtual patients -> pooled refit. Claims the same novelty, does not
  cite Valitalo, **discards published RSEs entirely** and pools pseudo-individuals
  as if real. Runs the opposite direction to the binding function.
- **Published elsewhere**: GMA (Rahmandad 2017) is our estimator, in systems
  science. Cite it.
- **Genuinely unclaimed**: the covariate treatment — the per-source granularity
  rule of §3, node truncation to the studied range, and `tau^2` on a binding
  function.
- **Unchecked**: PAGE/ACoP/ISoP abstract archives. Pharmacometrics methodology
  surfaces there first. Check manually before asserting absence in print.
