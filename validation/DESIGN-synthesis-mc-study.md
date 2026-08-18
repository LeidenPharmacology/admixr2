# MC study: aggregate likelihood vs parameter likelihood

Pre-registered design. Written **before** the harness was run; hypotheses in §6 are
committed to in advance so they can be wrong.

## 1. Truth `M*`

One-compartment IV bolus, individual level:

```
log CL_i = lcl + b_CL * x_i + eta_CL,i        x_i = log(WT_i / 70)
log V_i  = lv  + b_V  * x_i + eta_V,i
eta ~ N(0, diag(om_CL^2, om_V^2))             y_ij = f(t_ij) * (1 + eps), eps ~ N(0, sig^2)
```

| lcl | lv | b_CL | b_V | om_CL | om_V | sig |
|---|---|---|---|---|---|---|
| log(4.2) | log(48) | 0.75 | 1.00 | 0.28 | 0.20 | 0.13 |

`M*` is the oracle: every metric is reported **relative to `M*` itself**, so the
floor is known and a method can be scored on how much it loses, not on an
arbitrary scale.

## 2. The publications (S = 5)

Each source has its own design, covariate distribution, size, and **model form**.
Forms are deliberately heterogeneous — this is the "models that do not share
structure" case, and it is what makes the induced map `g_s` non-trivial.

| s | n | WT mean (SD) | sampling | form |
|---|---|---|---|---|
| 1 | 40 | 68 (12) | rich, single dose | full (`b_CL`, `b_V`) |
| 2 | 150 | 84 (19) | trough only, steady state | reduced: no `b_V` |
| 3 | 60 | 61 (9) | rich | full |
| 4 | 200 | 92 (22) | sparse (2 pts) | reduced: no covariate at all |
| 5 | 30 | 105 (17) | rich | reduced: no `b_V` |

`n` and precision are deliberately **anti-aligned** in the base cell: the small
studies sampled richly and publish tight exponents; the large ones took troughs and
publish loose ones. This is the configuration in which the two objectives separate;
the aligned configuration is a factor level (§4).

## 3. Methods compared

| arm | what it is |
|---|---|
| `ORACLE` | `M*` itself. The floor. |
| `AGG` | aggregate marginal likelihood on stratified re-simulated blocks, `K = 8` |
| `PAR-FE` | parameter likelihood, fixed effect: `theta_s ~ N(g_s(phi), Sigma_s)` |
| `PAR-RE` | as `PAR-FE` with between-source `tau^2` (Paule-Mandel) |
| `NAIVE` | precision-weighted average of the published coefficients, **ignoring `g_s`** |

`NAIVE` is the strawman that isolates what the induced map contributes: it is what
someone would do by hand, and it treats a reduced-form source's coefficient as if it
estimated the same quantity as a full-form one.

## 4. Factors

| factor | levels |
|---|---|
| A. `n` vs precision | anti-aligned (base) / aligned |
| B. between-source heterogeneity | `tau = 0` / `tau = 0.10` on the log-scale structural parameters |

2 x 2 = 4 cells. Heterogeneity is drawn per source per replicate as
`b_s ~ N(0, tau^2)` added to the pseudo-true value — i.e. studies genuinely differ
for reasons no covariate explains.

## 5. Evaluation

Three target populations, reported **separately** (never pooled):

| region | WT mean (SD) | relation to sources |
|---|---|---|
| `INTERP` | 75 (14) | inside the source range |
| `EXTRAP-LO` | 45 (8) | below every source |
| `EXTRAP-HI` | 120 (20) | above every source |

Four metrics, all four requested:

1. **Predictive `-2LL`** (primary). Draw new subjects from `M*` at a standard
   design; score their observations under each synthesised model, marginal over
   `eta`. Reported as **excess over `ORACLE`, per subject**. A proper scoring rule,
   so it grades the variance model too — which matters, since the two methods
   disagreed most about `om` and `sig`.
2. **Dosing error.** % error in the dose achieving a target AUC at the target
   population's median weight. `AUC = D / CL`, so this is a direct read on the
   clearance-covariate relationship. The interpretable headline.
3. **Aggregate moment error.** Relative error in predicted `(E, V)` for an unseen
   study design + population. Note this is close to what `AGG` optimises, so it is
   expected to favour `AGG`; reported for that reason, not despite it.
4. **Parameter recovery.** Bias and RMSE on `(lcl, lv, b_CL, b_V, om, sig)` vs `M*`.

Also recorded: `Q` and its p-value from `PAR-FE`, and the estimated `tau`, to check
the goodness-of-fit test detects factor B.

## 6. Hypotheses (committed in advance)

- **H1** `PAR-FE` beats `AGG` on parameter recovery when `n` and precision are
  anti-aligned, and they tie when aligned. *(Partly established: ~2x efficiency.)*
- **H2** `AGG` wins or ties on `INTERP` predictive `-2LL` — it minimises predictive
  discrepancy weighted by `n`, which is nearly the evaluation metric.
- **H3** `PAR-FE` wins on `EXTRAP-LO` and `EXTRAP-HI`, because extrapolation is
  governed by the covariate coefficients it estimates more efficiently.
- **H4** Under `tau > 0`, `PAR-FE` degrades (its likelihood is then misspecified)
  and `PAR-RE` recovers most of the loss; `AGG` degrades more gracefully but emits
  no warning.
- **H5** `NAIVE` is badly biased in every cell, because sources 2, 4 and 5 do not
  estimate the same quantity as sources 1 and 3.

H2 and H3 disagree with each other by design. If they both hold, the recommendation
is: interpolate -> aggregate, extrapolate -> parameter. If neither holds, the study
has told us something.

## 7. Two tiers

| tier | publications generated by | replicates | purpose |
|---|---|---|---|
| 1 | `theta_s ~ N(g_s(M*) + b_s, Sigma_s)`, asymptotic | 1000 / cell | power |
| 2 | fitting source `s`'s form to simulated **individual** data by marginal ML | 40 / cell | shows tier 1 is not rigged |

Tier 1's draw makes `PAR`'s link function correct **by construction**, which is a
thumb on the scale in its favour. Tier 2 exists solely to check that the tier-1
conclusions survive when publications are produced the way real ones are. **Any
conclusion that holds in tier 1 but not tier 2 is a tier-1 artefact and must be
reported as such.**

## 8. Harness

Tier 1 runs on a closed-form 1-compartment model with Gauss-Hermite moments, not on
rxode2 — an `rxSolve` costs ~11 ms before it integrates anything, and 1000
replicates x 4 cells x 5 arms is not reachable through it. A separate
cross-check fits a handful of replicates through admixr2 itself and confirms the
analytic harness and the package agree on `AGG` to within MC error.
