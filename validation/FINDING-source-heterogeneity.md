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

## What this means -- FIRST ANSWER, SUPERSEDED

> The reading below was that the model lacks a between-study variance component
> and needs `tau^2`. **It is wrong, and the arithmetic of the simulation says so
> plainly: `absorb.R` DECLARES exactly the `C_src` it draws each source from, so
> `tau^2` is zero by construction and there was never any heterogeneity to
> absorb.** A multiplicative heterogeneity factor was built on this reading and
> did roughly half the job for the wrong reason -- it widened an interval to
> cover a bias. The actual cause is in the next section.
>
> Kept because the DIAGNOSIS above it -- means reproduce 98% of the bias, four
> mechanisms eliminated by measurement -- is what led to the right answer, and
> because "sources disagree by more than they report" remains a real and
> separate thing that `tau^2` would be for.

The SEs are right for what the model says. The model has no BETWEEN-STUDY
variance component, so when published sources disagree by more than their
reported `C_src`, that disagreement has nowhere to go and lands on the residual
and omega. Three possible responses, increasing in cost: detect it with a
Cochran's-Q-style check; report the residual and omega SEs as CONDITIONAL on
homogeneity; or model `tau^2`.

## What it actually was: the OBJECTIVE's weight, not a missing variance

The aggregate -2LL splits exactly along its two sufficient statistics,

    n log|Vt| + n tr(Vt^-1 V) + n r'Vt^-1 r
      = [(n-1) log|Vt| + n tr(Vt^-1 V)]     Wishart:  n V ~ W(n-1, Vt)
      + [log|Vt| + n r'Vt^-1 r]             mean:     ybar ~ N(mu, Vt/n)

so a study's MEAN is scored against `Vt/n` -- the sampling error of `n`
patients, and nothing else. True for a digitised study. False for a model
source: its moments were generated at the published `theta_src_hat`, so with
`D = dE/d(theta_src)` the mean's real covariance is

    Vt/n  +  D C_src D'

and the second term does not shrink in `n` while the first does. At `n = 150`
the optimizer was being told the residuals were two orders of magnitude more
surprising than they are. Nothing in a model with one profile per study explains
a between-source difference in the mean, so the only knobs that widen `Vt`
absorbed it -- `add.err` and `omega`, exactly as the decomposition above found.

**It is a WEIGHT error, which is why no post-hoc covariance could undo it.** A
sandwich rescales an interval; it cannot move a biased point estimate. That is
also the general lesson: `project_adf_proper_ofv`'s "a post-fit sandwich fixes
it with bit-identical estimates, do NOT rewrite the objective" was about a case
where the ESTIMATES were right. Here they were not. Do not conflate the two.

### The fix, and why it is source-level

`.admSrcMeanCorr()` adds the term. `D C_src D'` is not estimated: it is data --
the source's own reported covariance pushed through its own model at its own
published parameters -- so it is fixed across the fit and costs nothing per
iteration.

A source BANDED into J strata is ONE source; every stratum's mean is shifted by
the same `delta`, so the weight is dense across the strata, and applying it per
stratum would count one draw J times -- the same error already fixed on the meat
side. Woodbury collapses the joint form to one `dim(theta_src)`-square term:

    corr = log|I + C M| - q' S q,   M = sum_j n_j D_j' Vt_j^-1 D_j
                                    q = sum_j n_j D_j' Vt_j^-1 r_j
                                    S = (I + C M)^-1 C

zero at `C = 0`, so a fit with no model source is bit-for-bit unchanged.

It has a consequence worth stating on its own. Because `D` is STACKED, a
contrast between two strata of the SAME source has variance
`c'Ac + |(D_j - D_k) delta|^2`, and `D_j ~ D_k` because the source's parameters
move all its strata together -- so **the shared `delta` cancels and a
within-source covariate contrast keeps its full weight**, while a between-source
contrast is correctly down-weighted. A per-study `Lambda_j` would have
attenuated both, which is the "explaining the covariate away" failure. It is
structurally impossible in this form.

`tau^2` is the one place that worry does bite: with a covariate varying only
BETWEEN sources it competes with the coefficient for the same `k-1` degrees of
freedom, and at `k = 2` they are exactly confounded. That is a first-class
reason to leave `tau^2` out rather than a detail, and if it is ever added it
must be structured so it cannot reach within-source contrasts.

### Measured, same seeds as the table above

| | before | after |
|---|---|---|
| `add.err` bias | +0.00849 | **+0.00023** |
| `tcl` bias | +0.02100 | **+0.00045** |
| `om` bias | +0.01710 | **+0.00059** |

Gradients re-verified against central FD after the change: adgh 1.4e-09, admc
4.6e-09, adfo 1.3e-08, in both the `cov` and `var` branches, with `C_src = 0`
reducing to the old objective exactly.

### The SEs that follow, measured

100 replicates, two model sources drawn from exactly the `C_src` they declare,
`covMethod = "r,s"`. The target for each is the source's own SE pooled over two
sources, `C_src/sqrt(2)`:

| | median SE | target | |
|---|---|---|---|
| `tcl` | 0.03529 | 0.03536 | -0.2% |
| `add.err` | 0.00354 | 0.00354 | exact |
| `om.eta.cl` | 0.00857 | 0.00849 | +0.9% |

`tcl`: bias +0.00045, SE/SD 1.093, **coverage 0.980** at MC SE 0.022 -- about
1.4 MC SE above nominal, i.e. mildly conservative and never anti-conservative.
Zero of 100 replicates exceeded 5x the median SE.

And in the HONEST regime -- `MODE=all`, every source parameter drawn from the
`C_src` it declares, which is the case a real meta-analysis is in:

| | bias before -> after | SE/SD before -> after | coverage before -> after |
|---|---|---|---|
| `add.err` | +0.00867 -> **+0.00055** | 0.380 -> **0.966** | 0.540 -> **0.930** |
| `tcl` | +0.02159 -> **+0.00083** | 0.655 -> **1.073** | 0.930 -> **0.980** |
| `om.eta.cl` | +0.01806 -> **+0.00095** | 0.325 -> **0.997** | 0.860 -> **0.970** |

Nominal 0.950 at MC SE 0.022, so all three are within about 1.4 MC SE. The
lowest, `add.err` at 0.930, is exactly what the residual-only control gives on
this sampling design -- `c(1, 3, 8, 24)` has no early sample, so a constant
residual and a time-varying omega are poorly separated. Design, not weight; see
"One separate, smaller effect: SAMPLING DESIGN" above.

The residual-only arm is the control that the fix does not disturb what it
should not: SE/SD 0.964 -> 0.965 and coverage 0.930 -> 0.930, unchanged. That is
mechanical rather than lucky -- only `add.err` is drawn there, and for additive
error `add.err` does not enter `E`, so its column of `D_E` is ~0 and the mean
correction has nothing to carry.

> **A false alarm worth recording, because it names the invariant.** An earlier
> run of the same simulation reported a mean SE of 3.63 against an empirical SD
> of 0.033. It was not a heavy tail and not `.admHeteroFactor` (phi comes out at
> 1 here, and switching it off changes nothing to five figures). It was a run
> launched from a tree where the OBJECTIVE carried the source-mean term and
> `.admScoreCross` did not yet. **The bread and the meat must come from the same
> objective**; a sandwich assembled from two of them is not approximately
> anything. It also shows up as a lone source returning 34.0 where its own
> `C_src` says 0.012, which is the cheap way to detect it.

### What it costs

Exact reproduction of a lone source becomes O(1/n) rather than exact -- 4.8e-03
at `n = 100`, 1.9e-04 at `n = 6400`. `log|I + C M|` depends on `Vt`, so it pulls
on the variance parameters with an O(1) force against the objective's own O(n)
restoring force. The log-determinant is NOT optional: without it
`E[dF/dPsi] != 0` at the truth, because differentiating a quadratic form whose
weight moves with `Psi` leaves `tr(W^-1 dW/dPsi)` behind. Exactness at every `n`
was an artefact of the mean and the covariance sharing one weight -- the thing
that was wrong.

adgh, admc and adfo carry the term. **adirmc REFUSES** a model source with
`C_src`: its `mu` is the importance-weighted mean computed inside
`irmc_inner_nll_cpp` and never exists in R, so the term cannot be formed.
Refusing beats scoring without it.

Cost in solves: one extra moment pass per gradient for adgh and admc (adfo needs
none -- Pass 1 already caches every study's moments), and admc's batch paths
fall back per configuration. Paid only when a model source declared a `C_src`.

## Reproducing

`scratchpad/{pool_exact,exact_suite,absorb,djac,optnoise,cov_exactSE}.R`.
`absorb.R` takes `MODE` in {mean, resid, all}, `NT` (timepoint count) and
`NREP`; it is the one that isolates the cause.

---

# The corrected chi-square and TIC are BOTH essential under MBMA

Measured 2026-08-22, 150 replicates under the null (beta = 0), three model
sources with `C_src`, MC SE 0.018 on a rate.

| statistic | measured | target |
|---|---|---|
| `dOFV` mean | 12.215 | 1.000 (chi2_1) |
| `dOFV` median | 3.632 | 0.455 |
| KS of `dOFV` vs chi2_1 | D = 0.512, p = 0.0000 | -- |
| **naive LRT at 0.05** | **0.493** | 0.050 |
| **corrected LRT at 0.05** | **0.068** | 0.050 |
| **AIC selects the full model** | **0.667** | 0.157 |
| **TIC selects the full model** | **0.180** | 0.157 |
| `bcrcl` estimate | -0.0020 | 0 |

Under MBMA the uncorrected statistics are not slightly off, they are unusable: a
49% false-positive rate on the likelihood ratio test, and AIC choosing a
covariate that is not there two times in three. Both corrections land within
about one Monte Carlo standard error of nominal.

**Why TIC works, exactly.** `Delta p_eff = 155.79 - 143.57 = 12.22`, against a
mean `dOFV` of 12.215. The penalty matches the inflation to three significant
figures, which is what Takeuchi's correction is for. TIC selects the full model
when `dOFV > 2 * Delta p_eff = 24.4`; with `dOFV ~ 12.2 * chi2_1` that is
`P(chi2_1 > 2.0) = 0.157` -- the same threshold AIC would use if the objective
were a proper likelihood.

**So `p_eff` needs no repair.** It is not a parameter count under MBMA -- it
scales with `n`, measured at exactly 2.000x for a doubling -- and it does not
need to be one. It is the optimism correction, and it carries the right
magnitude. An earlier draft of this work concluded that "not a parameter count"
implied "not usable" and suppressed TIC; that was wrong on the evidence here,
and the suppression was reverted. The only thing misleading is the NAME:
"effective parameters" invites comparison with `p`, and under MBMA the two are
unrelated.

The `n`-scaling is likewise not a defect. The objective scales with `n` and so
does the penalty, so `Delta TIC` between two models is scale-consistent and
selection is unaffected.

---

# Calibration after the weight fix: the LRT holds, TIC's PENALTY does not

Re-measured 2026-08-22 under the corrected objective, because none of the
earlier calibration transfers -- changing the weight changes `H`, changes `J`,
and therefore changes both the eigenvalues the corrected test uses and the
`tr(H^-1 J)` TIC penalises with. 140 null replicates (beta = 0), three model
sources at different renal centres, `n = 120`, 80 of them carrying the
per-replicate eigenvalue.

**The target is not "dOFV ~ chi2_1".** Under misspecification the statistic is a
weighted sum of chi-squares and cannot be chi2_1. Three separable claims:

## 1. The corrected LRT is CALIBRATED -- holds, and it matches standard ADM

Reported in the format `algorithm/adf/DERIVATION-DOFV.md` section 8 uses, so the
MBMA row is directly comparable with the standard-ADM row it must reproduce.
360 replicates carrying the eigenvalue:

| arm | mean | var | P(>3.841) | 95% CI | KS p |
|---|---|---|---|---|---|
| target chi2_1 | 1.000 | 2.000 | 0.050 | | |
| unscaled `dOFV` | 4.823 | 54.127 | 0.392 | [0.341, 0.442] | 0.0000 |
| ORACLE-scaled | 1.000 | 2.327 | 0.042 | [0.021, 0.062] | 0.9558 |
| **SANDWICH-scaled** | **0.962** | **2.123** | **0.039** | **[0.019, 0.059]** | **0.6950** |

Standard ADM, pooled to 700 replicates, for comparison: mean 0.986, var 2.087,
P 0.057 [0.041, 0.077], KS p 0.552. **The MBMA case reproduces it.**

    c_sandwich 5.1411 +- 0.0442     c_oracle 4.8234 +- 0.3878
    gap 0.318 = 0.81 SE             not significant

Rejection rates at 400 replicates: 0.0075 / 0.0350 / 0.0800 / 0.1675 against
0.01 / 0.05 / 0.10 / 0.20, KS of the corrected p against U(0,1) D = 0.038,
p = 0.610. Naive: 0.235 / 0.390 / 0.468 / 0.558, KS p = 0.

> **A PHANTOM WAS CHASED HERE, AND IT IS THE ONE DERIVATION-DOFV.md ALREADY
> RECORDS.** At 240-400 replicates `c_sandwich/c_oracle` read 1.296, then 1.108,
> and a second-order plug-in-bias story was built on it -- that the sandwich is
> biased upward by `sum_m tr(C Cov(g_m))`. **The inequality is real algebra but
> its magnitude here is zero to within noise.** `c_oracle` is variance-like:
> resolving a 5% gap at 3 SE needs ~7200 replicates. The standard-ADM derivation
> hit exactly this (18% high at 400, 1% agreement at 700) and its section 8 says
> in terms: *anything computed from `var(gam_hat)` at 300-400 replicates should
> be treated as having a 3-SE window, not read as a point estimate.* Read that
> before running another oracle comparison.
>
> What the detour did leave behind is worth keeping, because both CONFIRM the
> machinery rather than question it: the sandwich reproduces a refit-based
> `sum_m g_m' C g_m` to **six significant figures** (0.011104 vs 0.0111037), and
> the Omega check below.

## 1b. Omega = D C_src D' -- the MBMA analogue of probe 30

`DERIVATION-DOFV.md` verifies `J = sum_s G_s Omega_s G_s'` by checking `Omega_s`
against the empirical `Cov(t_s)` with no fitting at all -- 20000 draws, relative
Frobenius error 0.86-1.65%, contraction 1.012. A MODEL source replaces `Omega_s`
by `D C_src D'`, a FIRST-ORDER delta approximation, and nothing analogous had
ever been run. It is measurable to ~1% because both sides are deterministic
functions of `theta_src`:

| | |
|---|---|
| relative Frobenius error | 0.0487 (MC floor at 400 draws: 0.0707) |
| trace ratio model/empirical | 1.0015 |
| MEAN block | Frobenius 0.046, trace ratio 0.9930 |
| **contraction `g'Omega g / g'Cov_emp g`** | **median 1.0111**, IQR [0.975, 1.048] |

**1.0111 against probe 30's 1.012.** The delta approximation for a source's
moments is as good as the normal-theory `Omega` is for a data source.
`scratchpad/omegachk.R`.

## 1c. Every link, side by side

| link | standard ADM | MBMA |
|---|---|---|
| derivation steps (1)-(4) | oracle-scaled chi2_1, KS 0.97 | KS 0.956 |
| quadratic expansion (2) | rel.err 1.8e-3 | `Var(g)/E[naive]` 3.939 vs `E[dOFV]` 3.939 |
| sandwich formula | 0.966 vs `var(gam_hat)` | 6 s.f. vs refit |
| `Omega = Cov(t)` | 1.012 | 1.0111 |
| `c_sandwich = c_oracle` | 1% at 700 reps | 0.81 SE at 360 |
| `dOFV / c ~ chi2_1` | KS p 0.552 | KS p 0.695 |

## 2. Superseded reading of the uniformity check

## 3. E[delta p_eff] = lambda -- REFUTED, and it is what TIC rests on

    E[delta p_eff] 3.514 +- 0.262      E[lambda] 5.038 +- 0.100     gap 5.4 SE
    correlation(lambda, delta p_eff) = 0.103     delta p_eff < 0 in 9.3%

**Correlation 0.10.** `delta p_eff` is not a noisy estimate of `lambda`, it is a
different quantity. The reason is structural: `tr(H^-1 J)` is the sum of the
eigenvalues, so adding a parameter contributes its own `lambda` AND moves every
other one, since both `H` and `J` change. The difference of traces is
`lambda_new + sum(lambda_i^full - lambda_i^red)` and that second sum has no sign
constraint -- hence the 9.3% of fits where the penalty is negative and the
larger model wins by construction.

| null selection rate | |
|---|---|
| AIC | 0.507 |
| TIC (delta p_eff) | 0.236 |
| **TIC (lambda)** | **0.125** (target 0.157, MC SE 0.041) |

> **And the repair turns out to be something already in the table.** Selecting
> when `dOFV > 2 lambda` is `dOFV/lambda > 2`, i.e. `chi2_1 > 2`, i.e. the
> corrected `p` below 0.157. **A lambda-penalised criterion IS the corrected
> LRT.** So no new column is warranted: for nested models the `p` column already
> carries the right answer, and TIC keeps its role only where there is no
> lambda to have -- comparisons that are not nested.
>
> `anova.admFit()` warns when `delta p_eff < 0`, which is the case where reading
> the TIC column would silently mislead.

Reproduce: `scratchpad/ticbig.R` (chunked, checkpointed; `SEED0`/`TAG` per
chunk) and `scratchpad/ticpool.R` (pools every chunk and runs all three tests).
