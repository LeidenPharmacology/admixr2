# Covariate support in admixr2 — handoff

Branch `feature/covariate-quadrature`, worktree
`C:\package\admixr2\.claude\worktrees\feature-covariate-quadrature`.
11 commits, **nothing pushed**. Tier 1 2108 pass / 0 fail; integration covariate 46;
`-grad` 34, `-nll` 5, `-cov` 53, `-pipeline` 154 — all unchanged from main.

Run R as `"C:/Program Files/R/R-4.5.3/bin/Rscript.exe"`. Use
`devtools::load_all(".")`, **never** `library(admixr2)` — the installed package has
none of this and fails with "the following parameter(s) are required for solving".
`test_dir()` from a plain Rscript needs `TESTTHAT_PARALLEL=false`.

---

## 1. What the feature is

A study may carry `cov_dist`, the distribution its subjects' covariates span:

```r
list(E = ..., V = ..., n = 300L, times = ..., ev = ...,
     cov_dist = list(WT = list(meanlog = log(72), sdlog = 0.28)))
```

Supported specs: `mu`+`sd` (normal), `meanlog`+`sdlog` (lognormal),
`values`(+`probs`) (discrete/categorical), `quantile` (a function of a uniform).
`cov` (the value the model is solved at) is derived from `cov_dist` when absent, via
`.admCovMeanOf()`.

The likelihood needs the **covariate-marginal moments**: `mu = E_{a,eta}[f]`,
`V = Cov_{a,eta}(f) + residual`. Covariate-induced between-subject variability
belongs *inside* `V_pred`, because that is where the observed `V_obs` carries it.

## 2. Three paths, chosen per study by `.admCheckCovariates()`

Field `.adm_cov_path` on each study unit.

| path | precondition | cost |
|---|---|---|
| `collapse` | bare `theta*COV`, normal, covariate appears once | closed form, free |
| `uq` | appears once, shares a parameter assignment with an eta, has a `cov` reference | `n_sim` rows + one ~11 ms probe |
| `rows` | **none** | `n_sim` rows (admc) / product grid (adgh) |

`collapse` inflates `Omega -> Omega + J Sigma_a J'` (`.admCovInflateL`). `uq` measures
`Delta(a)` from the model itself (`.admCovDelta`, one rxSolve at eta=0 across a
covariate grid) and replaces the affected eta column with `u = F_u^-1(Phi(z_j))`.
`rows` gives every subject its own covariate value; rxode2 evaluates whatever the
model contains, so any functional form / several parameters / no-eta parameters /
dose scaling all work.

**Routing is gradient-aware**: with any gradient mode, only `rows` is valid, so it
is selected. That is deliberate — every estimator defaults to a gradient, and
erroring would make covariate models fail out of the box.

## 3. Estimator support

| | admc | adgh | adfo | adirmc |
|---|---|---|---|---|
| collapse / uq / rows | yes / yes / yes | yes / →grid / yes | refused | refused |
| analytic gradient | yes | yes | — | — |
| Hessian (SEs) | via NLL-FD | via NLL-FD | — | — |
| joint units | refused | refused | — | — |

`adfo`/`adirmc`/joint units refuse via `.admRefuseCovariates()`. **Refusing matters**:
every study also carries a covariate *value*, so an unwired estimator does not fail —
it silently solves at the covariate mean and inflates omega (measured 0.30 → 0.44).

`.admCalcCov`/`.adghCalcCov` downgrade `use_grad` for covariate studies, because
`.admGradBatch` (five hand-built params frames) and `.adghGradNLL` are not covariate
aware; the NLL-FD Hessian route is. Same precedent as the beta-endpoint downgrade.

## 4. Established results

**Accuracy vs exact nested Gauss-Hermite**, model with WT on both `cl` and `v`, eta
on only one:

```
adgh product grid   E 1.06e-06   V 2.19e-06     0.024 s/eval
admc rows (n 16000) E 3.98e-05   V 3.11e-03     1.03  s/eval
```

adgh is ~1400x more accurate on V at ~40x less cost. **Use adgh for covariate work.**

**Analytic gradients** vs central FD of the same NLL, away from the optimum:
admc max rel 1e-4 (sens) / 1e-5 (fd); adgh max rel 1.05e-04 on gradients of
magnitude ~1900.

**Framework simulation-estimation** (`validation/framework-simulation-estimation.R`)
— individual subjects simulated, reduced to (E, V, n), individuals discarded:

```
param      truth   mean est      bias   emp SD   mean SE   SE/SD   coverage
tcl       0.0000    -0.0001   -0.0001   0.0266    0.0214    0.80      90.0%
tv        2.3026     2.3023   -0.0003   0.0016    0.0016    1.02     100.0%
tcov      0.7500     0.7339   -0.0161   0.0493    0.0446    0.90      85.0%
add.err   0.3000     0.2986   -0.0014   0.0053    0.0066    1.23     100.0%
om        0.3000     0.3048   +0.0048   0.0272         -       -          -
```

Point estimates unbiased (all within ~1.5 MC SE). **Reported SEs are optimistic** —
coverage 85–90% against nominal 95%. Predicted by `tr(V_pred^-1 V_obs)` treating
`V_obs` as Wishart-from-normal. 20 reps → ~4.9% MC error per coverage figure;
**rerun at ~100 reps, and across admc/adgh/adfo, before quoting.**

**Value of the covariate** — same data, with and without:
`delta -2LL 227.4` on 1 df; omega 0.3134 → 0.4942 (+58%); residual unchanged.

## 5. THE definitive result (`scratchpad/threeway.R`)

Three objectives, exact marginal data at truth, 3 identifying populations,
`tcl`/`omega`/`add.err` profiled out:

```
design       n   marginal         gh     taylor |  marg bias    gh bias  tayl bias
rich        10     0.7500     0.3045     0.3128 |    +0.0000    -0.4455    -0.4372
rich      1000     0.7500     0.3045     0.3128 |    +0.0000    -0.4455    -0.4372
sparse      10     0.7501     0.3477     0.3624 |    +0.0001    -0.4023    -0.3876
sparse    1000     0.7501     0.3477     0.3624 |    +0.0001    -0.4023    -0.3876
```

1. **Marginal moments recover truth exactly.**
2. **Node constructions are biased ~58% low**, and GH (exact quadrature) is *as
   biased as* Taylor — truncation is <2% of the error. More nodes do not help; it is
   the **target**, not the approximation.
3. **The bias is completely independent of n** (identical to 4 dp over 100x). More
   data converges more precisely on the wrong value.

Per subject at truth: `gh - marginal` = 8.47 (rich) / 5.21 (sparse);
`taylor - gh` = -1.16 / -0.71.

## 6. Identifiability (do not lose this)

With the covariate sharing an argument with an eta, the likelihood is **exactly**
flat along

```
theta' = theta + (b - b')*mu_a ,   omega'^2 = omega^2 + (b^2 - b'^2)*sd_a^2
```

Verified `spread 0.00000`. **One population identifies nothing**, and design
variation does not help — two studies with different doses/times/regimens but the
same covariate distribution are still exactly flat. Only differing covariate
**means** or **spreads** between studies identify it.

```
design                                0.45     0.60     0.75     0.90     1.05
A one population                     0.000    0.000    0.000    0.000    0.000
B 3 pops, MEANS differ, sd same     52.564   13.372    0.000   13.386   52.897
C 3 pops, SDs differ, mean same     28.078    8.811    0.000   10.595   46.121
D 3 pops, means AND sds differ     109.499   31.647    0.000   34.143  133.254
```

A covariate on a parameter with **no** random effect is not on this ridge.

An identifiability **warning** was drafted (`.admWarnCovIdentifiability`,
`.admCovSdOf`) but **NOT committed** — HvdB pushed back that multi-study identifies
it, and although the ridge test says design variation alone does not break it, the
disagreement is unresolved. Do not add it without settling that.

## 7. PAGE case study

`validation/page-abstract-mbma.R` (set `ADM_EST=adgh|admc`), plus
`page-abstract-taylor.R`, `page-abstract-debug.R`. Two published models
(Issaranggoon 2-cmt single infusion; Alsultan 1-cmt multiple dose) → `datagen()` per
design → one unifying 2-cmt fit. Weight enters `cl` (has eta), `vp`, `q` (no eta) and
the mg/kg dose via `f(centr) <- WT`.

```
              time    objective      RMS Ayuthaya / Alsultan
admc (MC)     470 s     868.199          15.2% / 2.3%
adgh (grid)    41 s     868.202          15.2% / 2.3%
```

Cross-method agreement (0.003 objective units, ~0.4% worst parameter) ⇒ the 15.2%
Ayuthaya deviation is **genuine structural conflict** between a 2-cmt and a 1-cmt
generating model, not numerical error. Confirmed converged: LBFGS restarted from the
Taylor solution returns to 868.2020.

Taylor on the same data converges (297 s) to a materially different point —
`Vp` 5.5x different, omega +5%/+20% — but only 5.16 -2LL units worse under the
marginal objective, because the peripheral parameters are weakly identified.
**This case study demonstrates the workflow, not the method's advantage**: most of
the weight effect is on parameters with no random effect, so the absorption
mechanism barely operates. Do not use it as evidence for the construction; use §4–5.

Objective decomposition (`scratchpad/investigate.R`) shows Taylor trades `logdet`
against `trace`: MAR Ayuthaya 282.3/127.6/10.4 vs TAY 299.2/111.8/9.6.

## 8. Manuscript

`C:\Users\hidde\Nextcloud\1. Aggregate data\covariates\covariates-adm.Rmd`
(+ `references.bib`), formatted like `C:\package\ferx\letter-jpkpd\preprint.Rmd`.
Presents **only** the marginal construction — no comparison to node methods, per
HvdB (unpublished, so readers do not know they existed). A scope note explaining
this sits in an HTML comment, invisible in the PDF.

**It does not knit yet**: four chunks expect CSVs from `analysis/accuracy_sweep.R`,
`simulation_estimation.R`, `covariate_value.R`, `casestudy_vancomycin.R`. Working
versions of all four exist in `validation/`; they need adapting to write CSVs.

The sketch (`sketch_aggdata_covariates.docx`) claims simulating individuals and
pseudo-group averaging are "equivalent". §5 shows they are not. That claim must not
carry over.

## 9. Open items

1. **Coverage at ~100 replicates**, across admc/adgh/adfo. If the under-coverage
   appears in all three it is a property of the aggregate likelihood, which is a
   much stronger statement. Highest value.
2. **`uq` gradient** — the Newton solve is done (`.admCovUQuantile`, residual 9e-16,
   `d(u)/d(s)` matches FD to 7e-11, so the implicit-function derivative is ready);
   the chain rule into `.admGrad` is not. Until then `uq` never runs with gradients.
3. **adfo** covariate support (needs `J` per covariate value).
4. **Documentation** — `cov_dist` is in no Rd and no vignette.
5. **Issue #120** (nloptr `xtol_rel` default 1e-4) — evidence posted, fix not made.
   Do **not** set `xtol_rel = 0` (BOBYQA hits ROUNDOFF_LIMITED); pick a finite value
   against the measured noise floor. Affects every derivative-free fit.
6. **HvdB says the admr "naive method" did not work** — unresolved. `genopts.R`
   lines 196–225 draw covariates per simulated subject from Sobol dims *after* the
   Omega dims (the same construction as `.admCovRowsFor`), so on its face it should
   be the marginal method. Worth establishing what "did not work" meant.
7. Run HvdB's original PAGE script for a genuine external cross-check.

## 10. Process notes — these cost real time

**Silent no-ops, four times.** A discarded return value (`.admCheckCovariates`
annotations); `library(admixr2)` instead of `load_all`; a python `replace` that did
not match; an `assert` that did fire but was chained with a backgrounded run so its
failure scrolled past. All produced plausible output from code that was not doing
what was claimed. **Make the code announce what it did** (the PAGE script now prints
`est = ...`) and verify the edit landed before running anything.

**Experiments whose parameter is not identified, twice.** The original design table
and the first `n` sweep both used a single population and produced spurious
"biases" that were just the optimiser wandering a flat ridge. **Check
identifiability before designing the experiment.**

**Test at a point where the quantity is non-zero.** The first adgh gradient check
evaluated at the true parameters, where every component is ~0; it passed while
proving nothing. Regression tests now assert `max|grad| > 100` first.

**Independent references, not self-comparison.** Every real defect this session was
silent and produced plausible numbers. What caught them: exact nested quadrature as
a reference, or a contrast that *must* differ (`mean-only` beside the real path).

**`$` partial-matches on lists.** `study$cov` returned `cov_dist` whenever `cov` was
absent, and `NULL` when both `cov_dist` and `cov_rows` existed. All covariate reads
are `[["cov"]]`. Grep for `$cov` before adding new fields with that prefix.
