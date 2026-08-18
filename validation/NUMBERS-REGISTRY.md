# NUMBERS REGISTRY — validation/

Authoritative source of computed results for the **aggregate/ADM likelihood** paper.
Compiled 2026-08-16 by reading every script in `validation/` and loading every `.rds`
in the scratchpad with R 4.5.3.

**Revised 2026-08-18** against `RERUN-REPORT.md` (the mismatched-targets bug) and the
corrected-Taylor work. Four constructions are now distinguished throughout — see the
box in §0. Corrected or relabelled: `NODE.*`, `MATCHED.*`, `THREEWAY.*`, `INTEG.*`,
§4.1. Withdrawn: §4.7. Added: `HEADLINE.construction` (§2.11), `TAYLOR2.*` (§2.12),
`TAYLOR3.*` (§2.13). Retired conflict: §7.1. **`AGGMARG.*`, `OVN.*`, `THEORY.T2`/`T3`
and the headline are unaffected and were re-verified.**

**Scope rule applied throughout.** The paper now covers ONLY the aggregate/ADM
construction: re-simulate each published source model over a covariate grid to make
`(E, V, n)` blocks, fit one unified model to them. The **parameter likelihood /
binding function** method (PLL, PARAM-LL, GMA, indirect inference) is being removed.
Its numbers are catalogued and labelled **OUT OF SCOPE**; they are not wrong, they
are simply not in this paper.

Conventions used below:

- **`[VERIFIED]`** — value read from a `.rds` written by the script on disk, or
  re-run by me today and reproduced.
- **`[UNVERIFIED]`** — the script exists but no artefact exists and I did not run it,
  or the number appears only in prose with no producing code on disk.
- Every citable value carries a **label** in `CAPS.dotted` form. Cite the label, not
  a loose recollection of the number.

Paths:
`V/` = `C:\package\admixr2\.claude\worktrees\feature-covariate-quadrature\validation\`
`S/` = `C:\Users\hidde\AppData\Local\Temp\claude\C--package-admixr2\3ff305c7-64a5-4cb0-ba61-1436e2e9b16e\scratchpad\`

---

## 0. THE HEADLINE, STATED ONCE

Four findings dominate everything else in this directory, and three of them
invalidate numbers that were treated as settled a day earlier.

> **FOUR CONSTRUCTIONS — read this before any number in §2.4–§2.7 or §2.12.**
> `RERUN-REPORT.md` (2026-08-16) found that several arms labelled `gh` or `taylor`
> in the older covariate scripts are not methods at all. Let `f(a, η)` be the
> individual model, `g(a) = E_η[f|a]`, `Vc(a) = Cov_η(f|a)`:
>
> | # | construction | definition | verdict |
> |---|---|---|---|
> | **1** | **marginal / correct GH** | `E = E_a[g(a)]`, `V = E_a[Vc(a)] + Cov_a(g(a))`, then **one** `nll2` | correct; this is the shipped objective |
> | **2** | **stratified / matched** | `Σ_k w_k·n·nll2(obs_k, pred_k)`, both conditional on node *k* | correct, a genuine likelihood; exact |
> | **3** | **Taylor (corrected)** | 2nd-order expansion of the **moments**, then **one** `nll2` | usable approximation to 1, within a regime (§2.12) |
> | **4** | **mismatched targets** | `Σ_k w_k·nll2(obs_MARGINAL, pred_COND(a_k), n)` | **a category error, not a method** |
>
> Construction 4's displacement is a difference of *targets*, not a truncation
> error: the observed `(E, V)` has the covariate integrated out while each
> predicted block pretends the whole study sat at one covariate value, so the
> objective picks up a data-free `tr(V⁻¹ Cov_a(μ))` term that is minimised by
> driving the covariate coefficient to zero. **Wherever a number below measures
> construction 4, it is labelled as the size of a category error and must never
> be quoted as the bias of the aggregate method, of Gauss-Hermite, or of Taylor
> expansion.**

1. **The "aggregate method cannot recover a covariate a source omits" result is a
   CONSTRUCTION ARTEFACT, and it is fixed.** Every AGG arm written before
   2026-08-16 builds a *product* covariate grid for every source, including sources
   whose published model never contained that covariate. Such a source answers "no
   change" at every node, and the unified fit reads that manufactured flatness as
   `beta = 0`. `V/aggregate-marginal.R` replays the identical 100 replicates with a
   per-covariate stratify/marginalise split and the bias largely disappears
   (`bCL2` bias `-0.306 → +0.087` in `omit`, `-0.310 → -0.041` in `both`).
   See §2.1. **This supersedes the AGG arm of `overnight-*`, `three-way-replicates`,
   `mc-averaging-*` and `model-synthesis-two-likelihoods`.**

2. **Theory test T1 is invalid twice over** — not merely "known invalid". Its design
   is a tautology (both arms are correctly specified, so both must return truth), and
   the stored numbers are literally the optimiser's *start values*: `fit_marg` and
   `fit_cond` are started at the truth, Nelder-Mead never improves on the initial
   vertex, and `b_m = b_c = 0.75` bit-exactly at all seven `sd(x)`. Verified by
   restarting away from truth (§3.1). **Do not cite any T1 number, including the
   correlation and the "max rel err" lines.**

3. **The marginal/aggregate objective is EXACT when there is between-population
   contrast, and so is genuine stratification. What is badly displaced is the
   *mismatched-targets* construction (4), which is not a method.** Node-level
   predictions scored against one pooled/marginal `(E, V)` land at `−59.4%`,
   `−56.7%` or `−73.3%` depending on setup — that is the size of a category
   error. **Construction 2 (per-node observed blocks against per-node
   predictions) is exact, `0.0%`** (`NODE.argmin` row D). The corrected Taylor
   (construction 3) errs by `−1.9%` in its regime. With one population the
   marginal objective's profile is **flat** (`0 0 0 0 5.912` over
   `tcov ∈ 0.45…1.05`), so the argmin at `0.8956` is non-identification, **not a
   +19.4% bias**. See §2.4–§2.7 and §2.12.

4. **The headline attenuation is NOT the construction-4 bug, and that is now
   verified numerically rather than asserted.** `aggregate-marginal.R` and
   `overnight-simulation.R` are construction 2: observed block *k* is scored
   against predicted block *k* on the same grid, the two sides agree to
   **0.000e+00**, and profiling `bCL2` under that construction returns
   **0.4500 (+0.0%)** against a hypothetical mismatched **0.3515 (−21.9%)** —
   neither of which is the reported −68%. The −68% comes from stratifying a
   source on a covariate its own published model never fitted, which is finding
   1 above. See §2.11 (`HEADLINE.construction`). **Do not re-open this.**

---

## 1. SCRIPT INVENTORY

| # | script | one-line purpose | status | in scope? |
|---|---|---|---|---|
| 1 | `aggregate-marginal.R` | Per-covariate stratify/marginalise split; replays the 100 stored overnight replicates paired | **CURRENT — definitive** | ✅ yes |
| 2 | `overnight-simulation.R` | 2×2 (`baseline/struct/omit/both`) × 25 reps: IPD gold vs AGG vs PLL-sens vs PLL-rse | **CURRENT for IPD/PLL arms; AGG arm SUPERSEDED by #1** | ✅ partly (AGG only, corrected) |
| 3 | `overnight-summary.R` | Reporting layer over `overnight.rds` | **CURRENT** (but reports the superseded AGG) | ✅ partly |
| 4 | `theory-tests.R` | T1 DPI scaling / T2 attenuation identity / T3 omitted-covariate binding function | **T1 INVALID; T2, T3 CURRENT** | ✅ T2,T3 |
| 5 | `covariate-node-retest.R` | Isolates *which* node configuration is biased (A/B/C/D vs marginal) | **ROWS A/B/C ARE CONSTRUCTION 4 — SUPERSEDED by #42** (rows M, D unchanged) | ✅ yes, corrected |
| 6 | `covariate-threeway.R` | marginal vs gh vs taylor, 3 populations, 2 designs, 3 sizes | **`gh`/`taylor` arms ARE CONSTRUCTION 4 — SUPERSEDED by #43** (`marginal` unchanged) | ✅ yes, corrected |
| 7 | `covariate-matched-conditional.R` | Same, but separates ONE population from THREE; adds `gh_matched` | **`gh`/`taylor` arms ARE CONSTRUCTION 4 — SUPERSEDED by #44** (`marginal`, `gh_matched` unchanged) | ✅ yes, corrected |
| 8 | `covariate-integration-comparison.R` | Six-part study of the covariate integral: targets, taylor error, IS/ESS, s-collapse, `h`-scale, regime sweep | **target C and the `h`-sweep ARE CONSTRUCTION 4 — SUPERSEDED by #45** (targets A, s-collapse, ESS unchanged) | ✅ yes, corrected |
| 9 | `model-synthesis-omitted-covariate.R` | Recovering a covariate NO source model contains, via the covariate joint | **CURRENT** | ✅ yes |
| 10 | `mc-averaging-harness.R` | Closed-form re-implementation of admixr2 `adgh` on the `rows` path | **CURRENT** | ✅ yes |
| 11 | `mc-averaging-validate.R` | Gate: harness must equal admixr2 | **CURRENT — PASSES** | ✅ yes |
| 12 | `mc-averaging-study.R` | 5-publication model-averaging MC, AGG vs PAR-FE/RE vs NAIVE vs BEST-SINGLE | **PARTIAL RUN + AGG arm SUPERSEDED by #1** | ⚠️ AGG only, and contested |
| 13 | `mc-averaging-nocontrast.R` | Kills the between-study covariate contrast; measures SD inflation | **AGG arm SUPERSEDED by #1** | ⚠️ contested |
| 14 | `model-synthesis-two-likelihoods.R` | AGGREGATE (real admixr2) vs PARAMETER on the same two publications | **AGG arm SUPERSEDED by #1; PARAM arm out of scope** | ⚠️ partly |
| 15 | `page-framework-test.R` | PAGE vancomycin: design diagnostic + joint aggregate fit (Ayuthaya + Alsultan) | **CURRENT — the in-scope PAGE case study** | ✅ yes |
| 16 | `page-cross-dosing.R` | Each study's design, every model on it (cross-prediction panels) | **CURRENT** (figure script) | ✅ yes |
| 17 | `page-abstract-mbma.R` | Earlier PAGE ADM end-to-end fit | SUPERSEDED by #15 | ✅ but superseded |
| 18 | `page-abstract-panelA.R` | Obs-vs-pred per study, own/cross/MBMA | SUPERSEDED by #16 | ✅ but superseded |
| 19 | `page-abstract-weightplot.R` | Cmax/trough vs weight panel | SUPERSEDED by #16 / `page-param-likelihood.R` figure | ✅ but superseded |
| 20 | `page-abstract-taylor.R` | Same fit under the 2nd-order taylor combine | **CONTAINS A CONSTRUCTION-4 ARM, NOT RE-RUN — see §5.5.** The arm is deliberate (it prices the package's own combine); its finding must be stated as "that code path implements construction 4", never as "Taylor is inaccurate" | ⚠️ needs a construction-3 arm |
| 21 | `page-abstract-debug.R` | Pre-fit assumption checks for the PAGE reproduction | SUPERSEDED by #15's diagnostic | ✅ but superseded |
| 22 | `page-abstract-varonly.R` | PAGE with variances only (no off-diagonals) | **CURRENT** for the single-dataset claim; generalised by #23 | ✅ yes |
| 23 | `covariate-varonly-robustness.R` | Does the covariate coefficient survive var-only, over replicates with noise | **CURRENT — supersedes #22's generality** | ✅ yes |
| 24 | `covariate-collapse-endtoend.R` | Covariate marginalisation through the real estimators, 3 studies | **CURRENT** | ✅ yes |
| 25 | `covariate-uq-endtoend.R` | u-quantile marginalisation end-to-end through `admc`, 2 covariate forms | **CURRENT** | ✅ yes |
| 26 | `covariate-u-distribution.R` | The `u = Delta(a) + eta` framing across 5 functional forms at equal solve budget | **CURRENT** | ✅ yes |
| 27 | `covariate-dependent-copula.R` | Dependent covariates via an R-vine copula, checked 3 ways | **CURRENT** | ✅ yes |
| 28 | `covariate-rvine-sampling-check.R` | Is `inverse_rosenblatt()` a valid, CRN-safe vine sampler | **CURRENT** | ✅ yes |
| 29 | `framework-simulation-estimation.R` | Does fitting SUMMARIES recover the truth (the estimator-soundness check) | **CURRENT** | ✅ yes |
| 30 | `regression-check.R` | Theophylline: var and cov aggregate fits vs individual FOCEI | **CURRENT** | ✅ yes |
| 31 | `model-synthesis-standalone.R` | Self-contained "two published vancomycin models of different structure" narrative + figures | **CURRENT** (narrative/figure artefact) | ✅ yes |
| 32 | `gma-minimal.R` | Minimal GMA / composite indirect inference, per-node binding function | **CURRENT (bug fixed — see §5.1)** | ❌ OUT OF SCOPE |
| 33 | `gma-vs-ipd.R` | Stage 1: real FOCEI publications + IPD gold standard | **CURRENT**, n = 1 replicate | ❌ OUT OF SCOPE |
| 34 | `gma-vs-ipd-stage2.R` | Stage 2: unified model from the two publications only | **CURRENT**, n = 1 replicate; SUPERSEDED by #2 (25 reps/cell) | ❌ OUT OF SCOPE |
| 35 | `three-way-replicates.R` | AGG vs PARAM-LL vs IPD, 6 replicates, one configuration | **SUPERSEDED by #2** (25 reps/cell, 2×2) **and #1** | ❌ mostly |
| 36 | `model-synthesis-formal-likelihood.R` | Is there a formal likelihood for synthesising published models | **CURRENT** (conceptual) | ❌ OUT OF SCOPE |
| 37 | `page-param-likelihood.R` | PAGE binding-function fit + dosing figure | **CURRENT** | ❌ OUT OF SCOPE |
| 38 | `page-weights.R` | Four weight matrices for the parameter likelihood | **CURRENT (both known bugs fixed — §5.1, §5.2)** | ❌ OUT OF SCOPE |
| 39 | `param-likelihood-two-publications.R` | Parameter likelihood, 2 real-shaped publications, no truth | **CURRENT** | ❌ OUT OF SCOPE |
| 40 | `param-likelihood-peak-trough.R` | 2-cmt dense + 1-cmt steady-state peak/trough | **CURRENT** | ❌ OUT OF SCOPE |
| 41 | `param-likelihood-1cmt-2cmt.R` | Three 1-cmt publications → one 2-cmt unified model | **UNVERIFIED — no output artefact exists** | ❌ OUT OF SCOPE |
| 42 | `covariate-node-retest-v2.R` → `out-node-v2.txt` | #5 with every arm classified by construction; adds the corrected-Taylor row T | **CURRENT — supersedes #5** | ✅ yes |
| 43 | `covariate-threeway-v2.R` → `out-threeway-v2.txt` | #6 with `gh_marginal` / `taylor_moments` / `mismatched` / `mismatched_taylor` separated | **CURRENT — supersedes #6** | ✅ yes |
| 44 | `covariate-matched-conditional-v2.R` → `out-matched-v2.txt` | #7 with the same separation, plus the profile spread that settles the one-population case | **CURRENT — supersedes #7** | ✅ yes |
| 45 | `covariate-integration-comparison-v2.R` → `out-integ-v2.txt` | #8 with every error measured against construction 1; corrected `h`-sweep, IS arm and regime sweep | **CURRENT — supersedes #8** | ✅ yes |
| 46 | `headline-construction-check.R` | Re-implements `mkgrid`/`mall`/`nl1`/`fitb` verbatim and proves the headline scripts are construction 2 | **CURRENT — the verification of §0.4** | ✅ yes |
| 47 | `taylor-corrected.R` | The validated corrected 2nd-order Taylor on the MOMENTS: coefficient recovery, regime table, `hfrac` sensitivity | **CURRENT — defines construction 3** | ✅ yes |
| 48 | `aggregate-marginal-taylor.R` → `S/aggregate-marginal-taylor.rds`, `S/taylor-report.log` | #1's 100 replicates with a THIRD arm (`AGGt`, marginalised directions by construction 3) + a rank-one-dropped guard arm | **CURRENT — extends #1, does not supersede it** | ✅ yes |
| 49 | `taylor-arm-summary.R` | Reporting layer over `aggregate-marginal-taylor.rds` | **CURRENT** | ✅ yes |

### Artefacts in `S/` that are NOT produced by anything in `validation/`

`page_method_fits.rds`, `page_method_fits_admc.rds`, `page_method_pred.rds`,
`page_method_pred_admc.rds`, `page_ksweep.rds`, `sim_strata.rds` are written by
ad-hoc scratchpad scripts (`S/page_ksweep.R`, `S/page_method_pred.R`,
`S/sim_strata.R`, `S/fig_*.R`), not by any registered validation script.
**Treat as UNREGISTERED: do not cite without first re-deriving from a
`validation/` script.**

---

## 2. CURRENT SCRIPTS — HEADLINE NUMBERS

### 2.1 `aggregate-marginal.R` → `S/aggregate-marginal.rds` — **THE DEFINITIVE RESULT**

100 replicates (4 cells × 25), paired replicate-by-replicate with the overnight AGG
because it replays the *stored* `thA`/`thB` rather than refitting. `strat_x2` is
`TRUE` in `baseline`/`struct` (A did fit x2) and `FALSE` in `omit`/`both`.

**Internal consistency gate `[VERIFIED]`:** in `baseline` and `struct` the marginal
construction reduces to the original one, and the refit reproduced the stored AGG
with a paired difference of **exactly 0.0000 in all 50 replicates**. The replay is
therefore faithful; every difference below is the construction, not the optimiser.

**`AGGMARG.rmse.gold`** — mean RMSE against the individual-data (IPD) fit, MCSE in
parentheses `[VERIFIED]`:

| cell | reps | AGG (product grid) | **AGG-marginal** | PLL (out of scope) |
|---|---|---|---|---|
| baseline | 25 | 0.0470 (0.0068) | **0.0470 (0.0068)** | 0.0739 (0.0051) |
| struct | 25 | 0.0325 (0.0044) | **0.0325 (0.0044)** | 0.0475 (0.0050) |
| omit | 25 | 0.1305 (0.0061) | **0.0631 (0.0066)** | 0.0788 (0.0083) |
| both | 25 | 0.1024 (0.0048) | **0.0468 (0.0042)** | 0.0549 (0.0054) |

**`AGGMARG.rmse.truth`** — mean RMSE against the simulation truth `[VERIFIED]`:

| cell | IPD (gold) | AGG | **AGG-marginal** | PLL |
|---|---|---|---|---|
| baseline | 0.0575 | 0.0763 | **0.0763** | 0.0944 |
| struct | 0.0547 | 0.0642 | **0.0642** | 0.0717 |
| omit | 0.0490 | 0.1407 | **0.0854** | 0.0968 |
| both | 0.0529 | 0.1222 | **0.0742** | 0.0738 |

**`AGGMARG.bCL2`** — the covariate source A omits. Truth `0.450`. Mean, bias, MCSE
`[VERIFIED]`:

| cell | IPD (gold) | AGG | **AGG-marginal** | PLL |
|---|---|---|---|---|
| baseline | 0.452 (+0.002, ±0.013) | 0.479 (+0.029) | **0.479 (+0.029, ±0.016)** | 0.471 (+0.021) |
| struct | 0.452 (+0.002, ±0.014) | 0.420 (−0.030) | **0.420 (−0.030, ±0.015)** | 0.418 (−0.032) |
| omit | 0.464 (+0.014, ±0.009) | **0.144 (−0.306)** | **0.537 (+0.087, ±0.032)** | 0.560 (+0.110) |
| both | 0.429 (−0.021, ±0.016) | **0.140 (−0.310)** | **0.409 (−0.041, ±0.029)** | 0.443 (−0.007) |

**`AGGMARG.oCL`** — the variance channel. Truth `0.260` `[VERIFIED]`:

| cell | IPD (gold) | AGG | **AGG-marginal** | A's published |
|---|---|---|---|---|
| baseline | 0.260 | 0.276 | 0.276 | 0.275 |
| struct | 0.255 | 0.241 | 0.241 | 0.211 |
| omit | 0.257 | 0.315 | **0.274** | 0.310 |
| both | 0.261 | 0.275 | **0.241** | 0.235 |

**`AGGMARG.paired`** — paired test of AGG-marginal vs AGG on RMSE-to-gold
`[VERIFIED, computed today from the stored replicates]`:

| cell | mean improvement | SE | AGGm better in | paired t p |
|---|---|---|---|---|
| baseline | +0.0000 | 0.0000 | 0/25 (identical) | — |
| struct | +0.0000 | 0.0000 | 0/25 (identical) | — |
| omit | **+0.0674** | 0.0053 | **25/25** | 3.45e−12 |
| both | **+0.0556** | 0.0066 | **23/25** | 1.19e−08 |

> **Caveat to state in the paper.** `AGG-marginal` is still slightly biased low on
> `bCL2` in `both` (−0.041) and slightly high in `omit` (+0.087). The correction is
> large, not perfect, and both remaining biases are within ~1.5–3 MCSE.

> **This script is construction 2 and is NOT affected by the mismatched-targets bug
> — verified, not assumed.** See §2.11. A third arm, `AGGt` (the marginalised
> directions by corrected Taylor instead of quadrature), was added to the identical
> 100 replicates and changes nothing detectable — see §2.13 (`TAYLOR3.*`).

### 2.2 `overnight-simulation.R` / `overnight-summary.R` → `S/overnight.rds`

100 replicates, 25 per cell, all four cells complete `[VERIFIED]`. Truth vector
`TP = (lcl=log 4.2, lvc=log 30, lq=log 7, lvp=log 40, bCL1=.72, bCL2=.45,
bVc1=.95, oCL=.26, oVc=.20, sig=.11)`.

**`OVN.rmse.gold`** — deviation from the IPD fit `[VERIFIED]`:

| cell | reps | AGG | PLL (computed W) | PLL (reported SE) |
|---|---|---|---|---|
| baseline | 25 | 0.0470 (0.0068) | 0.0739 (0.0051) | 0.0732 (0.0058) |
| struct | 25 | 0.0325 (0.0044) | 0.0475 (0.0050) | 0.0498 (0.0052) |
| omit | 25 | 0.1305 (0.0061) | 0.0788 (0.0083) | 0.0811 (0.0057) |
| both | 25 | 0.1024 (0.0048) | 0.0549 (0.0054) | 0.0683 (0.0066) |

**`OVN.ipd.rmse`** — the IPD gold standard's own RMSE against truth, the ceiling for
any summary method `[VERIFIED]`: baseline **0.0575**, struct **0.0547**, omit
**0.0490**, both **0.0529**. *This is the most quotable in-scope number from this
script and it is unaffected by the AGG construction issue.*

**`OVN.beats.ipd`** — how often a summary method lands nearer truth than IPD
`[VERIFIED]`: AGG 6/25, 5/25, 0/25, 0/25 across baseline/struct/omit/both;
PLLs 3/25, 3/25, 2/25, 3/25. Frame as luck, not skill.

> **The AGG column of this table is the product-grid construction.** Use
> `AGGMARG.*` instead for any AGG claim. The IPD and PLL columns are unaffected.

### 2.3 `theory-tests.R` → `S/theory-tests.rds` — T2 and T3 only

**`THEORY.T3`** `[VERIFIED]` — a covariate-free model fitted to data from a
covariate-carrying one returns `omega_s^2 = omega^2 + b^2 Var(x)`.
Truth `omega = 0.30`, `b = 0.75`.

| sd(x) | `om_free` (fitted) | `sqrt(om² + b²·sd²)` (predicted) |
|---|---|---|
| 0.02 | 0.3003747 | 0.3003748 |
| 0.05 | 0.3023346 | 0.3023347 |
| 0.10 | 0.3092328 | 0.3092329 |
| 0.15 | 0.3204000 | 0.3204001 |
| 0.20 | 0.3354101 | 0.3354102 |
| 0.30 | 0.3749998 | 0.3750000 |
| 0.40 | 0.4242641 | 0.4242641 |

**`THEORY.T3.maxdev` = 1.95e−07** `[VERIFIED]`. This is the channel by which a
publication with no covariate term still carries covariate information, and it is
the mechanism `aggregate-marginal.R` exploits.

**`THEORY.T2`** `[VERIFIED, re-run today]` at `sd(x) = 0.30`, true parameters:

- `E_a[r' V⁻¹ r]` = **0.33275866**
- `rbar' V⁻¹ rbar` = **0.00000000** ← the only data-bearing term
- `tr(V⁻¹ Cov_a(mu))` = **0.33275866** ← pure covariate-spread penalty
- identity holds to **1.67e−16**

> Note for the write-up: `rbar' V⁻¹ rbar` is exactly zero *by construction* at the
> true parameters (the marginal mean *is* the weighted average of the conditional
> means), so T2's content is "at truth the whole of `E_a[r'V⁻¹r]` is a penalty
> carrying no data", not "the identity is a surprising decomposition". State it that
> way; it is still the reason the pooled/node-mismatch route drives `b → 0`.

### 2.4 `covariate-node-retest.R` / `covariate-node-retest-v2.R` — which construction is which

`[VERIFIED]`, re-run 2026-08-16 as `-v2` (exit 0; `out-node-v2.txt` on disk). Exact
noiseless data, 3 populations, truth `tcov = 0.7500`.

**`NODE.argmin`** — the three constructions that actually target the aggregate
likelihood. *Previously this label also carried rows A/B/C below; it no longer
does.*

| objective | construction | argmin | error | rel |
|---|---|---|---|---|
| **M** marginal (the shipped objective) | **1** | **0.7500** | +0.0000 | **0.0%** |
| **D** gh, **NODE data** + `V_cond` (genuine stratification) | **2** | **0.7500** | +0.0000 | **0.0%** |
| **T** taylor on the **MOMENTS** (corrected) | **3** | **0.7354** | −0.0146 | **−1.9%** |

**`NODE.mismatched`** — the same script's rows A/B/C, which are **construction 4**.
The numbers reproduce bit-for-bit; only their name changes. [*These now measure the
size of a mismatched-targets category error. They previously appeared inside
`NODE.argmin` and were read as "how biased the node route is" — a route that, done
properly, is exact.*]

| objective | construction | argmin | error | rel |
|---|---|---|---|---|
| A gh, pooled data, `V_cond` | **4** | 0.3045 | −0.4455 | −59.4% |
| B gh, pooled data, `V_marginal` | **4** | 0.4567 | −0.2933 | −39.1% |
| C gh, conditional obs `V` + `V_cond` | **4** | 0.3248 | −0.4252 | −56.7% |

**`NODE.decomp`** (population 2, n = 1, true parameters): `sum_j w_j r_j' V⁻¹ r_j =
0.852572`; `rbar' V⁻¹ rbar = 0.000000`; `tr(V⁻¹ Cov_a(mu)) = 0.852572`; identity
matches. **Reproduces bit-for-bit, but it is the decomposition of *construction 4's*
mean term** — i.e. the mechanism by which the category error drives `b → 0`. [*It is
a diagnosis of the bug, not a property of the shipped objective, which it was
previously read as.*]

**Reading**: constructions 1 and 2 are both exact, and construction 3 costs 1.9%.
Rows A–C are **not candidate objectives**: they need a pooled/marginal observation
scored against conditional predictions, which is a category error rather than a
different likelihood. In the MODEL-as-input currency you always generate
`(E_k, V_k)` per node, so construction 4 cannot arise from the data; it can only be
written into the code. The old framing "is the node objective biased, or merely a
different objective?" is **withdrawn** — it presupposed A was a candidate.

### 2.5 `covariate-matched-conditional.R` / `-v2.R` — one population vs three

`[VERIFIED]`, re-run 2026-08-16 as `-v2` (`out-matched-v2.txt` on disk). Truth
`tcov = 0.7500`, n = 100, exact data. Arm renames: `marginal` → **`gh_marginal`**
(it always *was* the 21-node GH implementation), `gh_matched` → **`stratified`**.

**`MATCHED.threepop`** — profile spread is over `tcov ∈ 0.45…1.05`:

| arm | construction | argmin | bias | profile spread |
|---|---|---|---|---|
| `gh_marginal` (was `marginal`) | **1** | **0.7500** | +0.0000 | 43.903 12.416 0.000 14.797 60.511 |
| `taylor_moments` **(new)** | **3** | **0.7353** | −0.0147 | 38.681 9.427 0.000 15.130 54.903 |
| `stratified` (was `gh_matched`) | **2** | **0.7501** | +0.0001 | 110.949 30.244 0.000 30.369 111.505 |
| `mismatched` (was `gh`) | **4** | 0.3045 | −0.4455 | 0.000 48.552 125.104 223.187 335.985 |
| `mismatched_taylor` (was `taylor`) | **4** | 0.3128 | −0.4372 | 0.000 40.317 101.772 177.100 260.251 |

[*The bottom two rows now measure the size of a mismatched-targets category error
and a 2nd-order expansion of it. They previously stood as "the gh result" and "the
taylor result" — i.e. as verdicts on quadrature and on Taylor expansion, which they
never were.*]

**`MATCHED.onepop`** — truth 0.7500, and **the profile is flat**:

| arm | construction | argmin | bias | profile spread |
|---|---|---|---|---|
| `gh_marginal` | **1** | 0.8956 | +0.1456 | **0.000 0.000 0.000 0.000 5.912** |
| `taylor_moments` **(new)** | **3** | 0.0637 | −0.6863 | 0.000 0.450 2.517 8.575 20.640 [boundary] |
| `stratified` | **2** | **0.7500** | −0.0000 | 34.794 9.422 0.000 9.573 35.655 |
| `mismatched` | **4** | 0.0502 | −0.6998 | [boundary] |
| `mismatched_taylor` | **4** | 0.0502 | −0.6998 | [boundary] |

> **This is the sharpest in-scope statement available, but state it as
> identification, not as bias.** The aggregate/marginal objective is exact *given
> between-population covariate contrast*. Without it the objective is **exactly flat
> from `tcov = 0.45` to at least 0.90**, so `0.8956` is simply where an optimiser
> stops inside a flat region. **Quote the profile `0 0 0 0 5.912`; do NOT quote
> "+19.4%" as a bias magnitude** — see §4.7. Construction 3 landing at the boundary
> in the same cell is that same flatness, not a separate defect. With three
> populations the profiles are strongly curved and constructions 1, 2, 3 all recover
> the truth.

### 2.6 `covariate-threeway.R` / `-v2.R`

`[VERIFIED]`, re-run 2026-08-16 as `-v2` (`out-threeway-v2.txt` on disk).

**`THREEWAY.bias`** — invariant to n (10 / 100 / 1000; the argmin does not move):

| design | `gh_marginal` [1] | `taylor_moments` [3] | `mismatched` [4] | `mismatched_taylor` [4t] |
|---|---|---|---|---|
| rich | **0.7500** (+0.0000) | **0.7353** (−0.0147) | 0.3045 (−0.4455) | 0.3128 (−0.4372) |
| sparse | **0.7501** (+0.0001) | **0.7401** (−0.0099) | 0.3477 (−0.4023) | 0.3624 (−0.3876) |

[*Column 2 is new and REPLACES the old `taylor` column: the corrected Taylor
expands the marginal MOMENTS and errs by −0.0147 / −0.0099, a **30× smaller error**
than the old `taylor` figures −0.4372 / −0.3876, which expanded the category error
instead and are withdrawn as Taylor results (§4.7). Column 3 is the old `gh`
column, unchanged in value, now labelled as construction 4 rather than as "what
Gauss-Hermite gives".*]

**`THREEWAY.perN`** — per-subject objective at the true parameters:

| design | `gh_marginal`/n [1] | `taylor_moments`/n [3] | `mismatched`/n [4] | `[3]−[1]` truncation | `[4]−[1]` target mismatch |
|---|---|---|---|---|---|
| rich | −22.86914 | **−22.76895** | −14.39631 | **+0.10019** | +8.47283 |
| sparse | 2.86045 | **2.86445** | 8.07093 | **+0.00400** | +5.21048 |

[*The truncation column now measures the truncation error of an expansion of the
right functional. The old figures `−1.15938` / `−0.71076` (old `taylor/n`
`−15.55568` / `7.36017`) were the truncation error of an expansion of the WRONG
functional and are withdrawn — see §4.7. Any prose saying "Taylor's truncation
error is large" must be reversed: it is one to two orders of magnitude smaller than
the target mismatch beside it.*]

**`THREEWAY.momacc`** `[VERIFIED, new]` — moment accuracy of construction 3 against
construction 1, at the true parameters, rich design:

| population | `mu_a` | `sd_a` | rel err E | rel err V |
|---|---|---|---|---|
| pop1 | −0.45 | 0.30 | 1.853e−04 | 2.167e−02 |
| pop2 | 0.00 | 0.55 | 2.252e−03 | 1.096e−01 |
| pop3 | 0.50 | 0.35 | 4.057e−04 | 2.737e−02 |

The expansion's weak point is `V`, and it degrades with `sd_a` — quantified as a
regime boundary in `TAYLOR2.regime` (§2.12).

### 2.7 `covariate-integration-comparison.R` / `-v2.R`

`[VERIFIED]`, re-run 2026-08-16 as `-v2` (`out-integ-v2.txt` on disk).

One population, `mu_a = 0`, `sd_a = 1`, `omega = 0.30` — so at truth
`tcov·sd_a/omega = 2.5`. That ratio matters for the Taylor row.

**`INTEG.targets`** — restructured: the old (B) and (C) are not targets a published
aggregate dataset can have. Truth `tcov` 0.7500, omega 0.3000:

| row | construction | minimiser | bias | omega at truth-fixed `tcov` |
|---|---|---|---|---|
| `gh_marginal` (was target A) | **1** | **0.7500** | +0.0000 | **0.3000** (−0.0000) |
| `stratified` **(new)** | **2** | **0.7500** | +0.0000 | **0.3000** (+0.0000) |
| `taylor_moments` **(new)** | **3** | **1.0864** | **+0.3364 (+44.9%)** | **0.5970** (+0.2970) |
| `mismatched` (was target C, `E[NLL]`) | **4** | 0.2000 | −0.5500 (−73.3%) | 1.2000 (+0.9000) |
| mixture (old target B) | 4-family | 0.8390 | +0.0890 (+11.9%) | — |

[*The `mismatched` row is unchanged in value and is the origin of the "−74%"
figure (§4.1), but it now measures the size of a mismatched-targets category error,
not "how biased the aggregate method is". The mixture row is a valid likelihood for
a **different** data model (all subjects share one unknown `a`) — keep it, but never
present it as an option for summary data.*]

> **Construction 3 genuinely fails here, and that is a real finding, not an
> artefact.** At `ratio = 2.5` the covariate spread is 2.5× the IIV and the
> second-order moment expansion is out of its regime: `+44.9%` on `tcov`, `+99%` on
> omega. Contrast `THREEWAY.bias`, where `sd_a ≤ 0.55` and construction 3 errs by
> −1.9%. **The corrected Taylor is not unconditionally safe** — see
> `TAYLOR2.regime` (§2.12) and `INTEG.regime` below for the boundary.

**`INTEG.accuracy`** `[VERIFIED, new]` — absolute error against construction 1
(reference: `gh_marginal` at 80 nodes via the s-collapse), at `tcov = 0.75`:
`gh_mom-5` **3.74e+00**, `gh_mom-9` **1.25e−02**, `taylor_moments` **1.54e+03**,
`mismatched` **5.42e+04**.

**`INTEG.taylor.h` — WITHDRAWN** (see §4.7); replaced by **`INTEG.taylor.hfrac`**
`[VERIFIED]`. The surviving, still-valid point is that the old `h` was in **ABSOLUTE
covariate units**, so a default `h = 2.0` lands arbitrarily for a covariate on any
other scale. Cite it from the corrected sweep, with `h` a fraction of `sd_a`, errors
against construction 1:

| `hfrac` | rel err E | rel err V | NLL err vs [1] |
|---|---|---|---|
| 0.05 | 1.072e−01 | 3.821e−01 | 2142.33 |
| 0.10 | 1.067e−01 | 3.805e−01 | 2119.89 |
| 0.25 | 1.039e−01 | 3.696e−01 | 1970.86 |
| 0.50 | 9.404e−02 | 3.316e−01 | 1534.79 |
| 0.75 | 7.845e−02 | 2.721e−01 | 1043.16 |
| 1.00 | 5.830e−02 | 1.961e−01 | 650.20 |
| 2.00 | 4.056e−02 | 1.758e−01 | 812.96 |

[*Large because `sd_a = 1` puts the expansion out of regime; the point of the table
is that the error is now monotone and interpretable in `sd_a` units instead of
sign-flipping across the sweep as the old one did. For a well-conditioned `hfrac`
sensitivity see `TAYLOR2.hfrac` (§2.12).*]

**`INTEG.ess`** — **UNCHANGED and explicitly re-verified.** IS reweighting of one
20000-draw eta pool to 9 covariate nodes; `ESS/N` decays as `exp(−(D/omega)²)`:
`D/omega = 0.5 → 7.79e−01` (15576 effective); `1.0 → 3.68e−01` (7358);
`2.0 → 1.83e−02` (366); `3.0 → 1.23e−04` (2.5); `4.0 → 1.13e−07` (0.0). The ESS
mechanics do not depend on what is combined afterwards.

**`INTEG.is.err`** — at `tcov = 0.75` the IS absolute error is **979.67** with min
`ESS/N = 5.74e−05`. [*Replaces **4308.06**, which was the error against construction
4. The corrected arm reweights the pool to per-node MOMENTS, combines them into
marginal moments and evaluates one likelihood. **Conclusion unchanged**: IS is
unusable at realistic covariate effect sizes.*] Full sweep: shift/omega
`0.33 → 9.70`, `0.67 → 20.52`, `1.33 → 120.80`, `2.33 → 646.11`, `3.33 → 1102.02`,
`5.00 → 979.67`.

**`INTEG.scollapse`** — **UNCHANGED**, reproduced verbatim; this section never
involved construction 4. For a LINEAR covariate effect and a normal covariate,
`s = theta_cov·a + b ~ N(theta_cov·mu_a, theta_cov²·sd_a² + omega²)` exactly.
Max relative difference vs nested 2-D quadrature: **E 9.99e−16, V 5.28e−14**.
One 1-D integral, no nodes, no extra solves.

**`INTEG.regime`** — **REPLACED**: relative errors are now against **construction
1**, not against target C. `ratio = theta_cov·sd_a / omega`:

| ratio | `tcov` | `gh_mom-9` [1] | `taylor_moments` [3] | IS-moments | `stratified` [2] | mismatch size [4] |
|---|---|---|---|---|---|---|
| 0.10 | 0.030 | 6.50e−16 | 2.72e−10 | 4.30e−08 | 2.56e−03 | 2.63e−03 |
| 0.25 | 0.075 | 3.73e−15 | 3.70e−07 | 1.40e−07 | 1.59e−02 | 1.87e−02 |
| 0.50 | 0.150 | 1.95e−15 | 6.73e−05 | 3.05e−06 | 6.27e−02 | 1.09e−01 |
| 1.00 | 0.300 | 2.81e−15 | 7.78e−03 | 8.99e−04 | 2.46e−01 | 1.12e+00 |
| 2.00 | 0.600 | 1.42e−07 | 5.93e−01 | 2.73e−01 | 1.14e+00 | 2.41e+01 |
| 3.00 | 0.900 | 6.26e−04 | 6.85e+00 | 5.70e+00 | 4.20e+00 | 2.35e+02 |

[*The old row read `gl-9 ≤ 8.24e−04`, `taylor2` reaching `7.52e−01` at ratio 1, and
a "target gap" — all measured against target C, i.e. against construction 4.*]
Readings: **9-node GH over the moments is essentially exact to ratio ≈ 2**
(`1.42e−07`) and still `6.26e−04` at ratio 3. The **corrected Taylor is good to
~1e−02 up to ratio 1 and breaks beyond ratio 2**. The `stratified` column is a
DIFFERENT likelihood on DIFFERENT data — it is not an approximation to [1], so read
it as a magnitude, not an error. The last column (formerly "target gap") is **the
size of the category error**, `2.63e−03 → 2.35e+02` as the covariate effect grows;
that is the number to quote when explaining why construction 4 must not be used.

### 2.8 `model-synthesis-omitted-covariate.R` `[VERIFIED, re-run today]`

Two published models, **NEITHER containing CRCL**; recovered through the covariate
JOINT. Truth `b_WT = 0.7500`, `b_CRCL = 0.4000`, IIV sd `0.3000`.

Published (derived): source A (`rho = 0.20`) `b_WT 0.8433`, IIV `0.3299`;
source B (`rho = 0.70`) `b_WT 1.0767`, IIV `0.3162`.

**`OMITCOV.recovery`**:

| projection | `b_WT` | `b_CRCL` | IIV sd |
|---|---|---|---|
| joint covariate law, K = 2 | **0.7500** | **0.4000** | **0.3000** |
| joint covariate law, K = 4 | **0.7500** | **0.4000** | **0.3000** |
| joint covariate law, K = 8 | **0.7500** | **0.4000** | **0.3000** |
| CRCL drawn UNCONDITIONALLY, K = 4 (control) | 0.9615 | **−0.2609** | 0.3118 |

> Exact recovery at every K, and the control gets the *sign* of `b_CRCL` wrong. This
> is arguably the strongest positive in-scope result in the directory and it is
> currently under-used.

### 2.9 `mc-averaging-validate.R` — **the harness gate** `[VERIFIED, re-run today]`

**`GATE.pkg_vs_harness`**:

| | lcl | lv | bCL | bV | oCL | oV | sig |
|---|---|---|---|---|---|---|---|
| truth | 1.4351 | 3.8712 | 0.7500 | 1.0000 | 0.2800 | 0.2000 | 0.1300 |
| admixr2 | 1.4351 | 3.8712 | 0.7495 | 1.0002 | 0.2800 | 0.2000 | 0.1300 |
| harness | 1.4351 | 3.8712 | 0.7500 | 1.0000 | 0.2800 | 0.2000 | 0.1300 |

`max |pkg − harness| = **4.91e−04**`; `max |harness − truth| = **2.80e−06**`;
**GATE: PASS**. This is the number quoted as "4.9e-04" in `mc-averaging-study.R`'s
header — confirmed.

### 2.10 `page-framework-test.R` → `S/page-fw.rds` — the PAGE vancomycin case study

`[VERIFIED from the rds]`.

**`PAGE.diagnostic`** — spread of the predicted mean across 8 weight bands:

- Ayuthaya (2-cmt, estimates exponents, n = 14): **29.5506 %**
- Alsultan (1-cmt, mg/kg dosing, CL and V both linear in WT, n = 72):
  **6.72e−08 %** — i.e. **exactly flat**; the design says *nothing* about weight at
  any n, and it carries 5.1× the subjects.

**`PAGE.exponents`**:

| | `clwt` | `vpwt` | `qwt` |
|---|---|---|---|
| Ayuthaya published | 0.97 | 1.07 | 1.19 |
| Ayuthaya ALONE (aggregate fit) | 0.9698 | 1.0614 | 1.1636 |
| **JOINT aggregate (both studies)** | **0.8742** | **0.9515** | **1.8714** |

Shift on adding Alsultan: `clwt −9.9%`, `vpwt −10.4%`, **`qwt +60.8%`**.

> ⚠️ **Conflict with the prose.** `page-framework-test.R`'s own header and
> `page-param-likelihood.R`'s header both describe this as exponents "attenuating
> ~10%" / "shifting toward zero". That is true of `clwt` and `vpwt` **only**.
> `qwt` moves the *opposite* way, by +60.8%, which is six times the size of the
> effect being narrated. Any sentence claiming uniform attenuation is contradicted
> by this artefact. Either report `qwt` explicitly, or restrict the claim to the
> CL and Vp exponents.

### 2.11 `headline-construction-check.R` — **the headline is construction 2, VERIFIED**

`[VERIFIED, re-run 2026-08-18, exit 0, reproduced line for line]`. This exists so
the question is never re-opened. It re-implements `mkgrid`/`mall`/`nl1`/`fitb`
**verbatim** from `aggregate-marginal.R` (whose call at line 106 is the identical
construction to `overnight-simulation.R:152`) and asks whether the observed and
predicted blocks are conditional on the same covariate node.

**`HEADLINE.construction`**:

| test | result |
|---|---|
| obs block *k* vs pred block *k*, same parameters, both grids, K = 16 (both covariates stratified) | `max abs diff` E **0.000e+00**, V **0.000e+00** |
| same, K = 4 (x1 stratified, x2 marginalised) | `max abs diff` E **0.000e+00**, V **0.000e+00** |
| the same predictions against a POOLED `(E, V)` from the same grid | `max abs diff` E **1.177e+01** (K = 16) / **6.833e+00** (K = 4) — plainly not zero |
| objective, matched vs mismatched, both stratified | **803.4109** vs 1083.5678 |
| objective, matched vs mismatched, x2 marginalised | **837.2395** vs 1013.2661 |
| profile `bCL2`, everything else at truth, **matched (what the scripts do)** | **0.4500** (truth 0.450, **+0.0%**) |
| profile `bCL2`, same grid, **mismatched (hypothetical)** | 0.3515 (**−21.9%**) |

> **Reading, and it is final.** `aggregate-marginal.R` and `overnight-simulation.R`
> are **construction 2**: `BL[[s]][[k]]` and `pr[[k]]` share both the block index
> `k` and the grid object `GL[[s]]`, and there is no pooled `(E, V)` anywhere in
> either objective. The matched construction recovers `bCL2` exactly; a hypothetical
> mismatched version of the same code would give −21.9%, which is **neither the
> reported −68% nor anything these scripts print**. The −68% attenuation therefore
> comes from what §0.1 and §4.3 already document — **stratifying a source on a
> covariate its own published model never fitted, so every node returns the same
> flat answer** — and not from the mismatched-targets bug. Two different defects;
> only the first is the headline's subject. **`AGGMARG.*`, `OVN.*` and the §7
> conflict-3 resolution stand as written.**

### 2.12 `taylor-corrected.R` — **construction 3, defined and priced**

`[VERIFIED, re-run 2026-08-18, exit 0, every table reproduced]`. Base R, no rxode2.
Expands the **moments** (which is what the aggregate likelihood consumes), not the
objective:

```
E_marg ≈ g(mu) + (sd²/2)·g''(mu)
V_marg ≈ Vc(mu) + (sd²/2)·Vc''(mu) + sd²·g'(mu)g'(mu)ᵀ
```

then evaluates the NLL **once** at those approximate marginal moments — three solves
instead of K nodes. The rank-one term `sd²·g'g'ᵀ` is the covariate-induced
between-subject covariance, i.e. to second order it *is* the
`omega'² = omega² + beta²·Var(a)` inflation of `THEORY.T3`, appearing as a rank-one
addition to `V`.

**`TAYLOR2.recovery`** — truth `tcov = 0.7500`, exact marginal data, 3 identifying
populations; argmin invariant to n (10 and 1000 give the same four decimals):

| design | `marginal` [1] | `gh` = mismatched [4] | `taylor-OLD` [4t] | **`taylor-NEW`** [3] |
|---|---|---|---|---|
| rich | 0.7500 | 0.3045 | 0.3128 | **0.7353** |
| sparse | 0.7501 | 0.3477 | 0.3624 | **0.7401** |

[*`taylor-NEW` is the corrected expansion and is the only Taylor number that should
be cited. `taylor-OLD` is kept verbatim in the script solely to show what the
retired arm was: a 2nd-order expansion of `E_a[NLL]`, i.e. of the category error,
whose displacement is overwhelmingly its target and not its truncation.*]

**`TAYLOR2.regime`** — moment accuracy of construction 3 against the exact marginal,
`hfrac = 0.5`, one population `mu_a = 0`, `sd_a = 0.40`;
`ratio = tcov·sd_a/omega`. `cond` = evaluating at the mean covariate, i.e. the
ecological plug-in, shown for scale:

| ratio | rel err E (taylor) | rel err V (taylor) | rel err E (cond) | rel err V (cond) |
|---|---|---|---|---|
| 0.10 | 7.447e−08 | 1.008e−05 | 1.024e−04 | 1.094e−02 |
| 0.25 | 2.895e−06 | 3.729e−04 | 6.386e−04 | 6.495e−02 |
| 0.50 | 4.556e−05 | 5.001e−03 | 2.534e−03 | 2.205e−01 |
| 1.00 | 6.819e−04 | 4.730e−02 | 9.887e−03 | 5.486e−01 |
| 1.50 | 3.089e−03 | 1.344e−01 | 2.145e−02 | 7.538e−01 |
| 2.00 | 8.413e−03 | 2.435e−01 | 3.617e−02 | 8.628e−01 |
| 3.00 | 3.352e−02 | 4.413e−01 | 7.012e−02 | 9.511e−01 |

**`TAYLOR2.boundary`** — **the citable regime statement.** Construction 3 is
accurate to **~1e−02 up to `ratio ≈ 1`** (`7.78e−03` on the NLL, `INTEG.regime`),
degrades through ratio 2 (`5.93e−01`), and **breaks at `ratio ≈ 2.5`, where it
returns `tcov = 1.0864`, a `+44.9%` error, and omega `0.5970`, `+99%`**
(`INTEG.targets`). Below the boundary it costs −1.9% (`NODE.argmin` row T) to
−0.0147 in `tcov` (`THREEWAY.bias`). Use it where the covariate-induced spread is at
most comparable to the IIV; use 9-node GH over the moments otherwise, which is
exact to `1.42e−07` at ratio 2.

**`TAYLOR2.hfrac`** — sensitivity to the differencing step at `ratio = 1.0`
(`h = hfrac·sd`, so scale-free, unlike the retired `INTEG.taylor.h`):
`0.05 → E 7.490e−04 / V 4.997e−02`; `0.10 → 7.469e−04 / 4.989e−02`;
`0.25 → 7.327e−04 / 4.933e−02`; `0.50 → 6.819e−04 / 4.730e−02`;
`0.75 → 5.980e−04 / 4.393e−02`; `1.00 → 4.820e−04 / 3.918e−02`. **Flat across a
20× range of steps** — the expansion's error is its truncation, not its
differencing.

### 2.13 `aggregate-marginal-taylor.R` / `taylor-arm-summary.R` — construction 3 inside the 100-replicate study

`[VERIFIED]` from `S/aggregate-marginal-taylor.rds` and `S/taylor-report.log`
(100 replicates, 357 s). Same 4 cells × 25 replicates and the same stored `thA`/`thB` as
§2.1, so all arms are paired replicate by replicate. Three arms:

- **AGG** — stratified on ALL covariates (the product-grid construction of §0.1);
- **AGGm** — per-covariate split, marginalised directions by **quadrature** (§2.1);
- **AGGt** — same split, marginalised directions by **construction 3** (`hfrac = 0.5`);
- **AGGt-r1** — AGGt with the rank-one `sd²·dE·dEᵀ` term dropped, a deliberate guard.

**`TAYLOR3.controls`** — the exact-zero controls, all `[VERIFIED]`. In `baseline`
and `struct` **both** covariates are stratified, so the three constructions are the
same objective and anything but an exact zero is a bug:

| check | result |
|---|---|
| `strat(TT)`: AGGt bundle == AGGm grid | TRUE (both sources) |
| `mall == mallT` on every fully-stratified cell (3 random draws each) | TRUE |
| AGG refit vs the stored overnight AGG, 100 reps | `max abs diff` **0.000e+00** |
| AGGm refit vs the stored AGGm, 100 reps | `max abs diff` **0.000e+00** |
| max abs diff AGGm−AGG, AGGt−AGGm, AGGt−AGG, cells `baseline`/`struct` | **0.000e+00** on all 10 parameters |
| the same three, cell `omit` | 6.382e−01 / 2.399e−01 / 6.379e−01 |
| the same three, cell `both` | 4.759e−01 / 1.107e−01 / 4.916e−01 |

**`TAYLOR3.rmse`** — deviation from the individual-data (IPD) fit, and from truth:

| cell | reps | AGG | AGGm | **AGGt** | PLL | IPD (gold) vs truth | AGG | AGGm | **AGGt** |
|---|---|---|---|---|---|---|---|---|---|
| baseline | 25 | 0.0470 | 0.0470 | **0.0470** | 0.0739 | 0.0575 | 0.0763 | 0.0763 | **0.0763** |
| struct | 25 | 0.0325 | 0.0325 | **0.0325** | 0.0475 | 0.0547 | 0.0642 | 0.0642 | **0.0642** |
| omit | 25 | 0.1305 | 0.0631 | **0.0679** | 0.0788 | 0.0490 | 0.1407 | 0.0854 | **0.0888** |
| both | 25 | 0.1024 | 0.0468 | **0.0474** | 0.0549 | 0.0529 | 0.1222 | 0.0742 | **0.0747** |

(columns 3–6 are RMSE-to-gold; columns 7–10 are RMSE-to-truth.)

**`TAYLOR3.bCL2`** — truth 0.450, mean and bias:

| cell | IPD (gold) | AGG | AGGm | **AGGt** | PLL |
|---|---|---|---|---|---|
| baseline | 0.452 (+0.002) | 0.479 (+0.029) | 0.479 (+0.029) | **0.479 (+0.029)** | 0.471 (+0.021) |
| struct | 0.452 (+0.002) | 0.420 (−0.030) | 0.420 (−0.030) | **0.420 (−0.030)** | 0.418 (−0.032) |
| omit | 0.464 (+0.014) | 0.144 (−0.306) | 0.537 (+0.087) | **0.539 (+0.089)** | 0.560 (+0.110) |
| both | 0.429 (−0.021) | 0.140 (−0.310) | 0.409 (−0.041) | **0.414 (−0.036)** | 0.443 (−0.007) |

`oCL` (truth 0.260): AGGt `0.276 / 0.241 / 0.275 / 0.242` across the four cells,
against AGGm `0.276 / 0.241 / 0.274 / 0.241`.

**`TAYLOR3.paired`** — the load-bearing comparison, AGGt − AGGm on `bCL2`, paired:

| cell | reps | mean d(bCL2) | se | p (t) | AGGt wins | m wins | p (sign) |
|---|---|---|---|---|---|---|---|
| baseline | 25 | 0.0000 | 0.0000 | (exact 0) | 0 | 0 | — |
| struct | 25 | 0.0000 | 0.0000 | (exact 0) | 0 | 0 | — |
| omit | 25 | **+0.0013** | 0.0165 | **0.94** | 10 | 15 | 0.4244 |
| both | 25 | **+0.0046** | 0.0099 | **0.65** | 13 | 12 | 1.0000 |

Paired `|error|` in `bCL2` and paired RMSE-to-gold (negative favours AGGt): `omit`
`−0.0055` (se 0.0142) and `+0.0048` (se 0.0040); `both` `−0.0001` (se 0.0099) and
`+0.0006` (se 0.0033). **Reading: replacing the quadrature over the marginalised
directions by construction 3 changes nothing detectable at 25 replicates.**

**`TAYLOR3.momacc`** — block-moment accuracy against a 21-node reference, averaged
over blocks. **The observed side is exactly zero for both arms** (`0.000e+00` in
`baseline`/`struct`; ~1e−15 in `omit`/`both`) — the approximation enters only through
the prediction:

| cell | side | relE AGGm | relV AGGm | relE AGGt | relV AGGt | rank-one Frobenius share |
|---|---|---|---|---|---|---|
| omit | obs | 9.863e−16 | 1.427e−15 | 7.688e−16 | 1.298e−15 | **0** |
| omit | pred | 1.370e−09 | 9.673e−07 | 1.030e−04 | 1.512e−02 | **0.2678** |
| both | obs | 1.108e−15 | 1.120e−15 | 1.050e−15 | 1.212e−15 | **0** |
| both | pred | 2.847e−09 | 2.127e−06 | 1.782e−04 | 1.853e−02 | **0.2819** |

**`TAYLOR3.rank1`** — the rank-one term `sd²·dE·dEᵀ` carries **26.78%** (`omit`) and
**28.19%** (`both`) of the Frobenius norm of the AGGt block covariance on the
**predicted** side, and **exactly 0** on the observed side. Dropping it moves `bCL2`
from 0.539 → **0.560** (`omit`) and 0.414 → **0.440** (`both`), and RMSE-to-gold from
0.0679 → 0.0730 and 0.0474 → 0.0500. It is the variance channel through which the
coefficient is recovered from a source that never fitted the covariate; it is not
optional.

**`TAYLOR3.cost`** — per objective call, both sources summed
(`rows/obj` = grid rows evaluated, `blocks` = Cholesky factorisations):

| cell | arm | rows/obj | blocks | obj calls | s / rep |
|---|---|---|---|---|---|
| omit | AGG / AGGm / **AGGt** | 800 / 800 / **700** | 32 / 20 / **20** | 2501 / 2498 / 2501 | 9.14 / 6.75 / **7.66** |
| both | AGG / AGGm / **AGGt** | 800 / 800 / **700** | 32 / 20 / **20** | 2501 / 2423 / 2480 | 7.70 / 5.47 / **6.37** |

(`baseline`/`struct` are identical across arms by construction: 800 / 32 and
~10 s/rep.)

**`TAYLOR3.scaling`** — the multiplier that `p` marginalised covariate directions put
on the inner expectation, quadrature `NC^p` against Taylor `1 + 2p`:

| marg. covariates `p` | quad `NC = 4` | quad `NC = 8` | quad `NC = 21` | taylor |
|---|---|---|---|---|
| 1 | 4 | 8 | 21 | **3** |
| 2 | 16 | 64 | 441 | **5** |
| 3 | 64 | 512 | 9261 | **7** |
| 4 | 256 | 4096 | 194481 | **9** |

> **This study sits at `p = 1`, `NC = 4`** — the cheapest cell of that table and
> therefore **the LEAST favourable setting for Taylor there is**. Read the measured
> speed numbers as a floor, not as the achievable gain.

**`TAYLOR3.hfrac`** — first 8 replicates of the two marginalising cells, `bCL2` and
`relV` against a 21-node reference: `omit` `0.25 → 0.523 / 1.579e−02`,
`0.50 → 0.563 / 1.512e−02`, `0.75 → 0.549 / 1.399e−02`, `1.00 → 0.567 / 1.243e−02`;
`both` `0.25 → 0.486 / 1.934e−02`, `0.50 → 0.481 / 1.853e−02`,
`0.75 → 0.478 / 1.718e−02`, `1.00 → 0.493 / 1.528e−02`. Consistent with
`TAYLOR2.hfrac`: the step is not the binding constraint.

---

## 3. OUT-OF-SCOPE NUMBERS (parameter likelihood / binding function)

Catalogued for completeness. **Do not put these in the aggregate-only paper.**

### 3.1 `page-param-likelihood.R` → `S/page-pll.rds` `[VERIFIED]`

**`PLL.page.exponents`** (label kept for traceability): `clwt 0.9230`,
`vpwt 1.0049`, `qwt 1.2171`. Induced map at Ayuthaya's parameters
`g_Alsultan(TH_ISS) = (lcl 1.2408, lv 2.5087, oCL 0.2774, oV 0.0095, sig 0.1302)`
against Alsultan's published `(1.0953, 2.2565, 0.1493, 0.1158, 0.1190)`.

### 3.2 `page-weights.R` → `S/page-weights.rds` `[VERIFIED]`

**`PLL.page.weights`** — four weight schemes, one shared Gauss-Newton Jacobian:

| scheme | `clwt` | `vpwt` | `qwt` |
|---|---|---|---|
| Ayuthaya published | 0.9700 | 1.0700 | 1.1900 |
| AGGREGATE (joint) | 0.8742 | 0.9515 | 1.8714 |
| PARAM-LL / rse | 0.9285 | 1.0025 | 1.2182 |
| PARAM-LL / sens | 0.8414 | 1.4460 | 1.6232 |
| PARAM-LL / natural | 0.9408 | 0.9545 | 1.2149 |
| PARAM-LL / none | 0.9389 | 0.9877 | 1.2023 |

**`PLL.page.spread`** — range of the point estimate across weight schemes:
`clwt 0.0994`, `vpwt 0.4915`, `qwt 0.4209`. The `sens` scheme is the outlier on
`vpwt` and `qwt`; the header's implied conclusion ("the point estimate barely
moves") holds for `clwt` only.

Consistency check: `PLL.page.exponents` (4 damped GN steps, damping 0.6) and the
`rse` row here (2 shared-Jacobian steps, damping 0.7) agree to ≤0.006 — as expected
for the same objective under different linearisations.

### 3.3 `model-synthesis-two-likelihoods.R` → `S/two-lik.rds` `[VERIFIED]`

| | lcl | lv | clwt | vwt | om.cl | prop.err |
|---|---|---|---|---|---|---|
| A published | 1.4351 | 3.8712 | 0.7500 | 1.0000 | 0.0800 | 0.1300 |
| **AGGREGATE** | 1.5720 | 3.8694 | **0.7832** | **0.6306** | 0.0602 | 0.1548 |
| PARAM-LL (fe) | 1.5907 | 3.8480 | 0.7229 | 0.8951 | 0.0630 | 0.1137 |
| PARAM-LL (re) | 1.5853 | 3.8269 | 0.6960 | 0.9737 | 0.0321 | 0.1097 |

`Q = 23.427` on 5 df, `p = 0.0003`; `tau = 0.0870`.

> **In-scope reading, and it is the artefact again.** Source B's model has no `vwt`
> term, yet the AGG arm stratifies B on weight anyway → B's blocks are flat in
> weight for V → AGG `vwt = 0.6306`, i.e. **−36.9%**. `clwt` (which B *does*
> estimate) is fine at **+4.4%**. This is `aggregate-marginal.R`'s mechanism showing
> up in a real admixr2 fit, at a per-parameter granularity. Worth re-running with
> the per-covariate split before citing anything from the AGG row.

### 3.4 `mc-averaging-study.R` → `S/mc-averaging.rds` — **PARTIAL RUN**

⚠️ **Only 3 of the intended cells are stored**: `anti/tau=0`, `anti/tau=.10`,
`aligned/tau=0`. `discrepant` never completed, and `aligned/tau=.10` — one of the
four cells in `DESIGN-synthesis-mc-study.md` §4 — **is not in the script at all**.
`R = 600` (design says 1000) and `K = 4` (design says 8). No `ORACLE` arm; a
`BEST-SINGLE` arm was added. **The design document does not describe the executed
study.**

`[VERIFIED]` `anti/tau=0` cell, means over 600 replicates:

| metric | AGG | PAR-FE | PAR-RE | NAIVE | BEST-SINGLE |
|---|---|---|---|---|---|
| bias `bCL` | **−0.2174** | 0.0007 | 0.0006 | 0.0730 | 0.0000 |
| bias `bV` | **−0.4220** | −0.0016 | −0.0016 | 0.0002 | 0.0000 |
| SD `bCL` | 0.0842 | 0.0339 | 0.0340 | 0.0331 | 0.0000 |
| SD `bV` | 0.1015 | 0.0469 | 0.0471 | 0.0542 | 0.0000 |
| dose err INTERP (%) | −3.24 | 0.04 | 0.04 | −1.87 | 0.00 |
| dose err EXTRAP-LO (%) | +9.07 | 0.02 | 0.02 | −5.72 | −0.00 |
| dose err EXTRAP-HI (%) | −12.53 | 0.09 | 0.09 | +1.58 | 0.00 |

`aligned/tau=0`: AGG bias `bCL −0.2215`, `bV −0.4394`; SD `bCL 0.1456`, `bV 0.3942`.
`anti/tau=.10`: AGG bias `bCL −0.2135`, `bV −0.4335`; PAR-RE recovers where PAR-FE
blows up (RMS INTERP `0.4088` vs `363.96`), which is the one clean, in-design result
this cell delivers.

> **The AGG column is the product-grid artefact.** Sources B, D and E do not
> estimate `bV` (D estimates neither coefficient), yet `make_blocks()` in
> `mc-averaging-harness.R` builds K covariate strata for every source regardless.
> D carries n = 200. The `bV` bias of −0.42 is therefore the same mechanism as
> `AGGMARG.bCL2`. **Do not report "AGG loses 29% of the covariate coefficient" from
> this study.**

### 3.5 `mc-averaging-nocontrast.R` → `S/mc-nocontrast.rds` `[VERIFIED]`

| cell | method | `bCL` bias | `bCL` SD | `bV` bias | `bV` SD |
|---|---|---|---|---|---|
| differing pops | AGG | −0.2147 | 0.0838 | −0.4278 | 0.1007 |
| differing pops | PAR-FE | 0.0009 | 0.0344 | −0.0015 | 0.0450 |
| same pop | AGG | −0.2218 | 0.1205 | −0.6382 | 0.1392 |
| same pop | PAR-FE | 0.0005 | 0.0259 | −0.0020 | 0.0457 |

SD inflation when the between-study contrast is removed: AGG `bCL 1.44×`,
`bV 1.38×`; PAR-FE `bCL 0.75×`, `bV 1.02×`. AGG at truth-pubs: differing
`bCL 0.5327 / bV 0.5784`; same pop `bCL 0.5264 / bV 0.3651`.

> Same caveat as §3.4 for the AGG arm. The *SD inflation* result (1.44× / 1.38×)
> is about identification and would likely survive the correction — but that has
> **not been measured**, so mark it UNVERIFIED-UNDER-CORRECTION if reused.

### 3.6 `three-way-replicates.R` → `S/threeway.rds` — **SUPERSEDED**

6 replicates, one configuration (A has no renal term ≈ the `omit` cell)
`[VERIFIED]`: RMSE vs truth IPD **0.0442**, AGG **0.1076**, PLL **0.0641**;
RMSE vs gold AGG **0.1219**, PLL **0.0502**; `bCL2` truth 0.450, gold **0.492**,
AGG **0.167**, PLL **0.489**.

Superseded by `overnight-*` (25 reps/cell across a 2×2) on power, and by
`aggregate-marginal.R` on the AGG construction. Its AGG `bCL2 = 0.167` is the same
artefact as `overnight`'s `0.144`.

### 3.7 `gma-vs-ipd.R` / `gma-vs-ipd-stage2.R` → `S/gma-ipd-stage1.rds`,
`S/gma-ipd-final.rds` — **n = 1 replicate** `[VERIFIED]`

`RMSE(GOLD, truth) = 0.0825`; `RMSE(PLL, truth) = 0.0988`;
`RMSE(PLL, GOLD) = 0.0275`. `bCL2`: truth 0.450, gold 0.496, PLL 0.471.
A single replicate cannot support a comparison; superseded by `overnight-*`.

### 3.8 `param-likelihood-*` `[VERIFIED for two of three]`

`S/param-two-pubs.rds`: `phi = (lcl 1.3633, lvc 3.2785, lq 2.1336, lvp 3.3870,
bCL 0.7969, bVc 0.8165, oCL 0.2710, oVc 0.2299, sig 0.1229)`; `Q = 175.23` on 7 df.
`S/param-peak-trough.rds`: `phi = (lcl 1.4179, lvc 3.3251, lq 2.0124, lvp 3.7723,
bCL 0.8027, bVc 0.9365, oCL 0.2654, oVc 0.2128, sig 0.1102)`; `Q = 3.332` on 7 df.
`param-likelihood-1cmt-2cmt.R` writes `param-1cmt2cmt.rds`, **which does not
exist** → **UNVERIFIED**.

### 3.9 `gma-minimal.R` — no artefact, **UNVERIFIED**

The corrected (per-node) version IS on disk — confirmed by reading it (§5.1). But it
writes no `.rds` and prints only. **Every number it reports, including the
per-node-vs-pooled `bCL`/`bV` percentage difference (the "33% coefficient change"
that motivated the correction), is UNVERIFIED.** Re-run it if any of it is needed;
it is out of scope for this paper regardless.

---

## 4. DO NOT USE

### 4.1 `FRAMEWORK-model-meta-analysis.md` §2.1 — "0.196 against a truth of 0.75 (−74%), with omega inflated ~40%"

- The coefficient half traces to `covariate-integration-comparison.R` target C.
  **Corrected value: `INTEG.targets` gives the minimiser at `0.2000`, bias
  `−0.5500`, i.e. `−73.3%`**, not `0.196 / −74%`. No script on disk produces
  `0.196`; grep across `validation/*.R` returns nothing. The `0.196` is from a
  superseded revision.
- **The omega half is flatly wrong.** Measured omega under target C is **1.2000**
  against a truth of `0.3000` — an inflation of **+300% (4×)**, not "~40%".
  **Corrected label: `INTEG.targets` omega row.**
- Also note the `0.2000` looks like a search-grid value; treat `−73.3%` as
  "large and negative", not as a precise figure. `NODE.mismatched` row A
  (**−59.4%**) is the better-conditioned measurement of the same phenomenon.
- **And the phenomenon must be RENAMED.** Both numbers are **construction 4**:
  they measure the size of a *mismatched-targets category error*, not the bias of
  the aggregate method or of quadrature. Suggested wording: *"scoring a marginal
  observation against conditional predictions displaces the coefficient by −59% to
  −73% and inflates omega 4×; it is a mismatch of targets, not a quadrature
  error."* [*Previously this bullet, and every sentence built on it, read as a
  verdict on the node/quadrature route; the genuine node route (construction 2) is
  exact — `NODE.argmin` row D.*]

### 4.2 Any number from theory test T1

`theory-tests.rds` columns `b_m` and `b_c` are **bit-exactly 0.75 at all seven
`sd(x)`** because `fit_marg`/`fit_cond` start Nelder-Mead *at the truth* and the
initial simplex never improves. Probing from a displaced start
(`b_start = 0.45`) returns `0.6166` at `sd(x) = 0.10` and `0.4677` at
`sd(x) = 0.40` — non-convergence, not bias, since both arms are correctly
specified and their argmin is at truth by construction.
**Do not cite**: the T1 correlation (`NA` — zero variance), "marginal err at
smallest/largest sd(x) = 0.00%", "node-objective max rel err = 0.000%", or any
DPI-scaling claim resting on T1. **Use `NODE.argmin` and `MATCHED.*` instead**
(the `-v2` versions in §2.4/§2.5) — they test the real question with the
constructions named, so a displacement cannot be misread as a method's bias.

### 4.3 Any AGG number computed on a product covariate grid

Superseded by `AGGMARG.*`. Specifically retire:

| do not use | source | corrected value |
|---|---|---|
| AGG `bCL2 = 0.144` (bias −0.306), cell `omit` | `overnight.rds` | **0.537 (bias +0.087)** — `AGGMARG.bCL2` |
| AGG `bCL2 = 0.140` (bias −0.310), cell `both` | `overnight.rds` | **0.409 (bias −0.041)** — `AGGMARG.bCL2` |
| AGG RMSE-to-gold `0.1305` (`omit`) | `overnight.rds` | **0.0631** — `AGGMARG.rmse.gold` |
| AGG RMSE-to-gold `0.1024` (`both`) | `overnight.rds` | **0.0468** — `AGGMARG.rmse.gold` |
| AGG RMSE-to-truth `0.1407` / `0.1222` | `overnight.rds` | **0.0854 / 0.0742** — `AGGMARG.rmse.truth` |
| AGG `bCL2 = 0.167`, AGG RMSE-to-gold `0.1219` | `threeway.rds` | no direct correction; use `AGGMARG.*` |
| AGG `bCL` bias `−0.2174`, `bV` bias `−0.4220` | `mc-averaging.rds` | **not corrected — do not report as a property of the aggregate method** |
| AGG `vwt = 0.6306` (−36.9%) | `two-lik.rds` | **not corrected — re-run with the per-covariate split first** |
| `FRAMEWORK` §3 "**−29%**" and "**−43%**" | `mc-nocontrast.rds` (`−0.2147/0.75 = −28.6%`, `−0.4278 = −42.8%`) | correct as a measurement of the **wrong** construction; **must be relabelled as such**, not as the aggregate method's performance |

### 4.4 `FRAMEWORK-model-meta-analysis.md` §4.1 — "`g_s` is the identity. Verified to `0.00e+00`"

`[VERIFIED contradiction]`. `S/two-lik.rds$gA_check` vs `PHI` gives
`max |g_A(phi) − phi| = **5.284e−04**` (worst element `clwt`, −5.28e−04; then `vwt`
−3.14e−04, `om.cl` −1.72e−04). **Corrected value: 5.3e−04, not 0.00e+00.** The claim
"the identity check passes to numerical tolerance" is fine; the specific figure is
not. (`param-likelihood-two-publications.R` also prints an identity check but writes
no artefact for it — that one is UNVERIFIED.)

### 4.5 `DESIGN-synthesis-mc-study.md` as a description of what was run

The pre-registration says 4 cells × 1000 replicates, `K = 8`, arms
`ORACLE/AGG/PAR-FE/PAR-RE/NAIVE`. What exists is 3 cells × 600 replicates,
`K = 4`, arms `AGG/PAR-FE/PAR-RE/NAIVE/BEST-SINGLE`, with `aligned/tau=.10`
never coded and `discrepant` never completed. **Do not cite the design document's
numbers as the study's parameters.** Tier 2 ("40 replicates fitting source forms to
simulated individual data") was never run at all in this harness.

### 4.6 Uniform-attenuation language for the PAGE exponents

See the ⚠️ in §2.10. `qwt` moves **+60.8%** on adding Alsultan, against `clwt`
−9.9% and `vpwt` −10.4%.

### 4.7 Anything measuring **construction 4** and named as a method

`RERUN-REPORT.md` (2026-08-16) re-derived every conclusion that rested on an arm
scoring **marginal observations against conditional predictions**. The values below
mostly reproduce bit-for-bit — what is withdrawn is the claim each was carrying.

| withdrawn | source | replacement |
|---|---|---|
| **`INTEG.taylor.h`, the whole six-value sweep** (`h=0.05 → −18758.5` … `h=3.00 → +14249.0`) | `covariate-integration-comparison.R` §5 | **No replacement as a Taylor result — it measured the differencing step of an expansion of construction 4, i.e. an implementation artefact.** The one surviving point (the step was in ABSOLUTE covariate units, so `h = 2.0` lands arbitrarily) is now carried by **`INTEG.taylor.hfrac`** (§2.7), and the well-conditioned step sensitivity by **`TAYLOR2.hfrac`** (§2.12) |
| **`THREEWAY.perN` "taylor truncation error" `−1.15938` (rich) / `−0.71076` (sparse)**, and the `taylor/n` values `−15.55568` / `7.36017` they came from | `covariate-threeway.R` | **`THREEWAY.perN` `[3]−[1]` = +0.10019 / +0.00400**, with `taylor_moments/n` = −22.76895 / 2.86445 (§2.6). Any prose saying "Taylor's truncation error is large" must be **reversed** |
| **`THREEWAY.bias` `taylor` column, −0.4372 / −0.3876**, as a Taylor result | `covariate-threeway.R` | **−0.0147 / −0.0099** (`THREEWAY.bias`, construction 3) — a 30× smaller error. The old numbers survive only as `mismatched_taylor`, a 2nd-order expansion of the category error |
| **`MATCHED.*` `gh` and `taylor` rows** as "gh"/"taylor" results | `covariate-matched-conditional.R` | the values stand as `mismatched` / `mismatched_taylor` (construction 4) in §2.5; as method results, replaced by `gh_marginal` and `taylor_moments` |
| **`NODE.argmin`'s `−59.4%`** (and rows B, C) as a property of the node route | `covariate-node-retest.R` | **No replacement as a node-route bias — it measures construction 4, and genuine stratification is EXACT at 0.0%** (`NODE.argmin` row D). The number itself survives, correctly named, as `NODE.mismatched` row A: the size of a category error |
| **`MATCHED.onepop`'s `+19.4%` as a bias magnitude** | `covariate-matched-conditional.R` | **No replacement as a bias — it is a NON-IDENTIFICATION result.** The `gh_marginal` profile over `tcov ∈ 0.45…1.05` is `0 0 0 0 5.912`, i.e. **exactly flat from 0.45 to at least 0.90**, so `0.8956` is where an optimiser stops inside a flat region. **Quote the flat profile.** The registry's framing ("an identification statement about the design, not a defect of the method") was and remains correct |
| **`INTEG.targets` target C `−73.3%` / omega `+300%`** as "how biased the aggregate method is" | `covariate-integration-comparison.R` | values stand; **relabel as the size of a category error** — see §4.1 |
| §7 conflict 1's inclusion of **`−43.7%`** in the list of mismatch measurements | `covariate-threeway.R` `taylor` | **removed** — it is an expansion of construction 4, superseded by the corrected Taylor's −1.9%, and leaving it in implies Taylor expansion is a 44%-biased technique, which the re-run refutes. See §7.1, now retired |

> **What is NOT withdrawn.** `THEORY.T2` computes the same mismatched quadratic form
> **on purpose**, as an algebraic identity, with no fit and no optimiser: it is the
> *proof* that construction 4 decomposes into `rbar'V⁻¹rbar + tr(V⁻¹Cov_a(μ))` and
> therefore attenuates `b`. It is a diagnosis of the bug, not an instance of it, and
> it is correctly labelled as written. Likewise `gma-minimal.R`'s `gb_pool` is
> pooled-marginal on **both** sides — a construction-1 control, deliberately printed
> as "the earlier error".

---

## 5. KNOWN BUGS — WHERE THEY WERE FIXED, AND WHERE THEY WERE NOT

### 5.1 The `function(p) mod_als` closure that ignored its argument

**Status: FIXED, and confined to one script.** `page-weights.R:89–94` now builds
Alsultan's model with a `sprintf`-parsed `mk_als(p)`, with an in-source comment
naming the failure ("a builder ignoring `p` makes every derivative zero and the
whole information block vanish"). Grep across `page-*.R` finds no surviving
argument-ignoring builder. `page-framework-test.R` and `page-cross-dosing.R` use a
bare `mod_als` correctly — as a *fixed target model* being refitted, not as a
`p`-dependent builder, so they were never affected. `page-param-likelihood.R`
correctly uses `mk_uni(phi)` inside `g_als`.

**Also fixed:** `gma-minimal.R` on disk is the **corrected, per-node** version — it
builds `make_blocks()` per weight node with `n_k = w_k · n_s` and refits the source's
form to all of them, and prints the pooled variant side by side as the documented
error. The buggy pooled-only version is not on disk.

### 5.2 Mean-only Fisher information (residual-error parameter weighted exactly zero)

**Status: FIXED in both places that compute a Fisher weight.**

- `page-weights.R:66` — `I_ij = n·(dE_i)' V⁻¹ (dE_j) + (n/2)·tr(V⁻¹ dV_i V⁻¹ dV_j)`,
  implemented at lines 79–83 with `P <- Vi %*% dV[[q]]` and
  `cov_t <- (n/2) * sum(P[[i]] * t(P[[j]]))`.
- `overnight-simulation.R`'s `infoW()` — same two-term form, per node, weighted
  `n_k = wk·n` (`nk*... + (nk/2)*sum(...)`). So the `PLLs` arm of the overnight study
  used the corrected information.

No other script computes a Fisher weight. **No surviving instance of the bug.**

### 5.3 Fitting a covariate-carrying model to pooled rather than per-node `(E, V)`

**Status: fixed in `gma-minimal.R` (§5.1); NOT the same bug as §0.1.** The §0.1
artefact is the mirror image — building per-node blocks for a source that has no
covariate coefficient to give. `aggregate-marginal.R` is the fix for that direction
and it is not yet applied to `mc-averaging-harness.R`, `three-way-replicates.R`, or
`model-synthesis-two-likelihoods.R`.

### 5.4 Dead code found today (harmless, but noted)

`mc-averaging-study.R:109` — `y[startsWith(names(y),"lcl")][4] <-
y[startsWith(names(y),"lcl")][4]` is a self-assignment on a temporary and has no
effect; `i <- grep(...)` on the same line is unused. The real discrepancy injection
is line 110 (into `yy`). Only reachable in the `discrepant` cell, which never
completed, so **no stored result is affected.**


### 5.5 The mismatched-targets bug (construction 4) — where it still lives

**Status: corrected in four scripts via `-v2` siblings (§1 rows 42–45); SURVIVES in
exactly one, deliberately.** All 36 remaining `.R` files in `validation/` were read
and classified for `RERUN-REPORT.md` §4. Originals are untouched; each `-v2` file is
headed by a comment naming the wrong arm and why.

**Still contains a construction-4 arm, and was NOT re-run:**

| script | arm | lines | why not re-run |
|---|---|---|---|
| `page-abstract-taylor.R` | `obj()`, fed by `nll_node()` | objective **111–123**, node construction **78–90** | needs `rxode2` + `nlmixr2est` + a full `datagen()` + two BOBYQA fits (minutes to tens of minutes) — out of the pure-base-R scope of the re-run |

Both construction-4 tells are present at once: lines 83–87 take the **marginal**
generated study (built at 67–72 *with* `cov_dist`, i.e. the covariate already
integrated out of `E`/`V`), strip `cov_dist`, and pin only the **prediction** to a
covariate value; lines 116–120 then finite-difference the **NLL scalars**. That is
`Σ_k w_k·nll2(obs_marginal, pred_cond(a_k))` expanded to second order — the same
object as the old `covariate-threeway.R` `taylor` arm, applied to the PAGE
vancomycin case and to a real `admixr2` `.adghNLL`.

> **This one is deliberate and the script is not confused.** Its header (lines
> 10–17) says it is pricing its own local implementation of a second-order taylor
> combine against the marginal fit at lines 144–149. **But the finding it produces
> must now be stated as "that combine implements construction 4", never as "the
> Taylor approximation is inaccurate"** — the corrected Taylor errs by **−1.9%**, not
> −44%, in the regime where it applies (§2.12).

> **There is no live PACKAGE bug.** `.adm_combine_nll` **does not exist** in
> `admixr2` — verified 2026-08-18, zero hits across `R/` and `src/` in both the
> main tree and this worktree; the only occurrences are comments in
> `page-abstract-taylor.R` and `covariate-integration-comparison.R` (lines 12–13,
> 91), which describe what those scripts' own local combines do. The construction-4
> code is confined to `validation/` scripts. **Nothing shipped needs fixing on
> account of this**, and no claim about a package code path may be sourced from
> those comments without first re-checking the package. The worktree's own
> `R/covariate.R` header (lines 18–40) already records that the node methods
> (`gl`/`gh`/`taylor`) were **removed deliberately** and states the same
> decomposition and the same `0.3045` / `0.4567` figures as `NODE.mismatched` —
> so the package and this registry now agree, in the same words.

**Registry exposure: none.** §6 lists every `page-abstract-*.R` as "NOT RUN TODAY,
no `.rds`", so no registry label depends on this script's output. It is on the §8
re-run list.

Scripts checked and found **clean** (construction 1 and/or 2, obs side a per-node
list indexed with the same `k` and the same grid object as the prediction):
`aggregate-marginal.R`, `overnight-simulation.R`, `aggregate-marginal-taylor.R`
(1 + 2 + **3**), `mc-averaging-{harness,study,nocontrast,validate}.R`,
`three-way-replicates.R`, `gma-vs-ipd-stage2.R`, `gma-minimal.R`,
`model-synthesis-{omitted-covariate,two-likelihoods,standalone}.R`, `theory-tests.R`,
`covariate-{collapse-endtoend,uq-endtoend,varonly-robustness,u-distribution,dependent-copula}.R`,
`framework-simulation-estimation.R`, `page-abstract-{mbma,debug,varonly}.R`,
`page-{framework-test,weights,param-likelihood}.R`, `param-likelihood-*.R`. Scripts
with no covariate-node integration at all: `model-synthesis-formal-likelihood.R`,
`covariate-rvine-sampling-check.R`, `gma-vs-ipd.R`, `regression-check.R`,
`overnight-summary.R`, `page-abstract-{panelA,weightplot}.R`, `page-cross-dosing.R`.

---

## 6. UNVERIFIED — LISTED RATHER THAN GUESSED

These appear in `FRAMEWORK-model-meta-analysis.md` or elsewhere marked
**[measured]**, but I could not tie them to a script + artefact on disk.

| claim | where it appears | status |
|---|---|---|
| "coefficient of **0.196**" | FRAMEWORK §2.1 | **UNVERIFIED** — no script produces it. Nearest current value `0.2000` (`INTEG.targets`). |
| "omega inflated **~40%**" | FRAMEWORK §2.1 | **CONTRADICTED** — measured `+300%` (§4.1) |
| "stratum mean plugged in → **+17% at K=2**" | FRAMEWORK §5.1 table | **UNVERIFIED** — no script on disk computes a stratum-mean-plug-in arm |
| "coverage of a nominal 95% interval: **0.58**" | FRAMEWORK §6 | **UNVERIFIED** — no coverage study exists in `validation/` |
| "`Sigma_B` with a **0.96** off-diagonal correlation and SEs of **14.7** against a true 0.05" | FRAMEWORK §6 | **UNVERIFIED** — no artefact stores a `Sigma_B`; out of scope anyway |
| "`g_s` identity verified to **0.00e+00**" | FRAMEWORK §4.1 | **CONTRADICTED** — 5.28e−04 (§4.4) |
| "over-identification power **4 → 69.6%**, **27 → 21.6%**" | FRAMEWORK §7.5 | marked [literature], not measured here — fine, but do not relabel as measured |
| "`R/utils.R:49` `as.integer(s$n)` truncates fractional `n_k`" | FRAMEWORK §10 | **UNVERIFIED against the current package source** — re-check before repeating |
| `gma-minimal.R` per-node vs pooled `bCL`/`bV` % difference (the "33% change") | task brief + `gma-minimal.R` prose | **UNVERIFIED** — script prints only, no artefact |
| `param-likelihood-1cmt-2cmt.R` results | script §41 | **UNVERIFIED** — `param-1cmt2cmt.rds` absent |
| `covariate-collapse-endtoend.R`, `covariate-uq-endtoend.R`, `covariate-u-distribution.R`, `covariate-dependent-copula.R`, `covariate-rvine-sampling-check.R`, `covariate-varonly-robustness.R`, `framework-simulation-estimation.R`, `regression-check.R`, `model-synthesis-standalone.R`, `model-synthesis-formal-likelihood.R`, `page-abstract-*.R` | — | **NOT RUN TODAY, no `.rds`.** Scripts look current; all are console/figure output. Any number from them must be regenerated before citing. |

---

## 7. OPEN CONFLICTS — STATED, NOT SMOOTHED

1. ~~**How biased is the mismatched node route?**~~ **RESOLVED 2026-08-16 —
   CONFLICT RETIRED. There was never a quantity here to reconcile.** The five
   figures once listed — `−59.4%` (`NODE.argmin` A), `−56.7%` (row C), `−43.7%`
   (`THREEWAY.bias` taylor), `−73.3%` (`INTEG.targets` target C), `−93.3%`
   (`MATCHED.onepop` gh, at the boundary) — are **all construction 4, or a
   second-order expansion of it**: marginal observations scored against conditional
   predictions. They measure the size of a **category error**, which naturally
   varies with the number of populations, with whether `V_pred` is conditional or
   marginal, and with where the optimiser's boundary sits — so their disagreement
   was never evidence about a method. **The genuine node/stratified route is
   construction 2, and it is EXACT: `NODE.argmin` row D, `0.7500`, `0.0%`**;
   `MATCHED.*` `stratified` agrees at `0.7500` / `0.7501`, as does
   `INTEG.targets` `stratified` (`tcov` 0.7500, omega 0.3000). If the size of the
   category error is needed — e.g. to justify why the construction must not be
   used — quote `NODE.mismatched` row A (**−59.4%**, the best-conditioned of the
   set, noiseless and agreeing with `THREEWAY.bias` `mismatched` to four decimals)
   or the `INTEG.regime` mismatch-size column (`2.63e−03 → 2.35e+02`), and **name
   it as a mismatch of targets**. `−43.7%` is **removed** from the set entirely: it
   was an expansion of the category error and is superseded by the corrected
   Taylor's −1.9%, so leaving it in implied that Taylor expansion is a 44%-biased
   technique, which the re-run refutes (§4.7).

2. **Is the marginal/aggregate objective biased at all?** `theory-tests.R` T1 says
   no (invalid, §4.2). `covariate-node-retest.R`, `covariate-threeway.R` and
   `covariate-matched-conditional.R` (3 pops) all say **no, exactly 0.0%** — under
   constructions 1 and 2, with construction 3 costing −1.9%. With ONE population
   the objective is **exactly flat** over `tcov ∈ 0.45…0.90` (profile
   `0 0 0 0 5.912`), so `0.8956` is non-identification, not bias. **Believe that as
   the caveat**: the marginal objective is exact given between-population covariate
   contrast, and *unidentified* without it. That is an identification statement and
   it belongs in the paper — but **stated as the flat profile, not as "+19.4%"**
   (§4.7).

3. **Does the aggregate method fail on a covariate a source omits?**
   `overnight.rds` / `threeway.rds` / `mc-averaging.rds` / `two-lik.rds` all say yes,
   badly. `aggregate-marginal.rds` says the failure was manufactured by the grid
   construction and mostly disappears once each source is stratified only on the
   covariates its own model conditions on. **Believe `aggregate-marginal.rds`**: it
   is the newest, it is paired replicate-by-replicate on the *identical*
   publications, it reproduces the old AGG bit-for-bit in the two cells where the
   two constructions coincide (so the replay is not a different estimator), and the
   improvement is 25/25 and 23/25 with p < 1e−7. The older AGG arms have **not**
   been re-run under the correction and their AGG columns should be treated as
   measuring the artefact, not the method.
   **Additionally settled 2026-08-18**: this conflict is not contaminated by the
   mismatched-targets bug. Both `aggregate-marginal.R` and `overnight-simulation.R`
   are construction 2, verified numerically — the two sides differ by `0.000e+00`
   and profiling gives matched `0.4500 (+0.0%)` against a hypothetical mismatched
   `0.3515 (−21.9%)`, neither of which is the reported −68%. See §2.11
   (`HEADLINE.construction`). **`AGGMARG.*` and `OVN.*` stand as written.**

4. **`overnight` vs `aggregate-marginal` on PLL.** Both report the same PLLs column
   (0.0739 / 0.0475 / 0.0788 / 0.0549) — no conflict; `aggregate-marginal.R` reads
   PLLs straight from the stored replicates. But note that under the correction
   **AGG-marginal now beats PLL on RMSE-to-gold in `omit` (0.0631 vs 0.0788) and
   `both` (0.0468 vs 0.0549)**, reversing the ordering the overnight study
   reported. Since PLL is being removed from the paper this matters mainly as a
   reason not to repeat the old ordering in prose.

---

## 8. WHAT SHOULD BE RE-RUN BEFORE THE PAPER IS FINALISED

Ordered by how much a paper claim depends on it.

1. `mc-averaging-harness.R::make_blocks()` under the per-covariate stratify/
   marginalise split, then `mc-averaging-study.R` and `mc-averaging-nocontrast.R`.
   Until then §3.4 and §3.5's AGG arms are uncitable.
2. `model-synthesis-two-likelihoods.R` with the same correction — it is the only
   place the artefact is visible in a **real admixr2 fit** (`vwt` −36.9%), so the
   corrected version would be the strongest possible demonstration.
3. `page-framework-test.R` with Alsultan marginalised rather than stratified on
   weight. Alsultan is *exactly* flat in weight (6.72e−08 %), which is the extreme
   case of the artefact, and the corrected joint fit is the natural PAGE headline.
4. `page-abstract-taylor.R` with a **construction-3 arm beside its current
   construction-4 one** (§5.5). It is the only remaining script carrying the
   mismatched-targets construction, and the only one that runs it against a real
   `admixr2` `.adghNLL` on the PAGE case, so the corrected version is what decides
   whether a second-order taylor combine is worth offering at all. (The package
   itself no longer has the node methods — `R/covariate.R` removed them
   deliberately — so this decides a future option, not a live defect.) No registry
   label depends on it today, so it blocks only a PAGE Taylor claim.
5. The `discrepant` and `aligned/tau=.10` cells, if the MC study is kept at all.
6. Anything in §6 that the paper wants to assert.
