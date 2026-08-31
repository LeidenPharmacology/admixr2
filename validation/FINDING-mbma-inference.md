# Inference for a meta-analysis of MODEL SOURCES: SE, LRT, AIC/TIC

Measured 2026-08-26 on `feature/covariate-quadrature`. 520 + 200 replicates,
adgh, klebsiella. Scripts: `mbma_sim.R`, `mbma_tic.R` (job scratch).

The prior campaign (`algorithm/adf/`, probes 01–42) established all of this for
**digitised** aggregate data. None of it had been checked where the studies are
**model sources** carrying their own `model_cov`, which is this branch's feature
and the case where `J` is a completely different object.

---

## Design

Three published models, each fitted on its own cohort and supplied through
`datagen()` with its own `C_src`, stratified on `SEX` (which all three
estimated), marginalising `WT` and a covariate `Z`.

The sources are **structurally correct** for the truth — they contain every term
except `Z`, whose true coefficient is 0 — so this isolates the model-source
question from the structural-misfit question. The randomness is the sources'
published estimates, drawn `theta_hat ~ N(theta_0, C_src)` per replicate, which
is the only thing that is random about a model source.

Nesting is by `bz <- fix(0)`. Verified equivalent to a model that never mentions
`Z` (objectives agree to 9.8e-11) and `bz` is correctly excluded from the reduced
fit's parameter set. It is also *required*: admixr2 refuses a `cov_dist` naming a
covariate the model does not read.

---

## 1. Standard errors — `covMethod = "r,s"` is right

400 replicates, n = 300.

| param | truth | empSD | mean SE | SE/SD | coverage |
|---|---|---|---|---|---|
| `bz` (null) | 0 | 0.1052 | 0.1020 | 0.969 | 0.910 [0.882, 0.938] |
| `bsex` | 0.18 | 0.0321 | 0.0313 | 0.974 | 0.948 [0.926, 0.969] |
| `tcl` | 1.609 | 0.1168 | 0.1139 | 0.975 | 0.908 [0.879, 0.936] |

SE/SD ≈ 0.97 throughout, so the scale is right. Coverage for `tcl` and `bz` is
0.91, about 4 MC SE below nominal — with the scale correct that is mildly heavy
tails, not a mis-sized SE. Not chased further.

Separately, on a lone reproducible source the sandwich returns the declared
`C_src` **exactly and independently of n** (0.08000/0.06000/0.00400/0.01000 at
both n = 100 and n = 6400), where `covMethod = "r"` falls as 1/sqrt(n).

## 2. Likelihood-ratio test — corrected works, naive is unusable

400 replicates, null covariate.

```
alpha = 0.05   corrected 0.0600 [0.0367, 0.0833]     naive 0.7950 [0.7554, 0.8346]
alpha = 0.01   corrected 0.0200 [0.0063, 0.0337]     naive 0.7100 [0.6655, 0.7545]

KS vs Uniform(0,1)   corrected D = 0.031, p = 0.851
                     naive     D = 0.745, p = 0.000
lambda  median 44.5  IQR 41.2 - 47.0  range 23.5 - 54.5
```

The corrected p-values are **uniform under the null**, which is stronger than
correct size at a chosen alpha. `anova()`'s `p` matched lambda computed directly
from the stored `H` and `J` to 2.2e-16.

`lambda ~ 45` here against 1.08–3.65 in the digitised-data work, because a model
source's `J` block is `G C_src G'` and bears no relation to `2H`.

## 3. AIC fails; TIC does not

200 replicates, both fits' `p_eff` recorded.

```
penalty for ONE extra parameter
  AIC   2.00 (fixed)
  TIC   92.08  [IQR 88.2 - 97.1]
  d(p_eff) / lambda = 1.019       <- TIC and the LRT agree on scale

picking the WRONG (larger) model
  AIC          0.8400 [0.789, 0.891]
  TIC          0.1450 [0.096, 0.194]    (nominal ~0.157 for an AIC-type rule)
  LRT p<0.05   0.0750 [0.039, 0.112]
  naive dOFV   0.7550 [0.695, 0.815]

absolute p_eff: full 125.2, reduced 78.9, for 6 and 5 estimated parameters
```

**AIC's `2p` is a derived quantity, not a convention.** The optimism of a
maximised objective is `tr(H^-1 J)/2`; the information equality `J = 2H` makes
that equal `p`, and that is the only reason AIC has the form it does. Here the
optimism is 125 where `p` is 6, so AIC under-penalises by a factor of ~21 and
takes the larger model 84% of the time.

**TIC is the corrected version and it works.** `d(p_eff)/lambda = 1.019` says
TIC's penalty and the LRT's rescaling are the same quantity, which is what the
theory predicts; TIC then selects at 14.5% against its own nominal 15.7%.

`p_eff` scales **exactly linearly in n** (slope of `log p_eff` on `log n` =
1.002; medians 41.9 / 168.4 / 669.6 at n = 100 / 400 / 1600). So `p_eff` is not
an effective parameter count and must not be read as one — but the *difference*
between nested models is the right penalty, and that is all TIC uses.

---

## 4. What lambda actually is

The objective never reads `model_cov` -- only the sandwich does -- so `H` is
independent of `C_src` while `J = sum_s G_s C_src G_s'` is linear in it. Both
consequences were measured on the vignette design, scaling `C_src` over a 64x
range:

```
 C scale     lambda  lam/scale    bcrcl est     bcrcl SE        p_eff
    0.25     0.8497     3.3986     0.553290     0.021592        17.23
    1.00     3.3986     3.3986     0.553290     0.043184        68.91
    4.00    13.5945     3.3986     0.553290     0.086369       275.66
   16.00    54.3782     3.3986     0.553290     0.172738      1102.63
```

`lambda` is **exactly linear in `C_src`** (ratio constant to 4 decimals),
the estimates are **bit-identical** across all four, SE scales as `sqrt(C_src)`
and `p_eff` linearly. With `p_eff ∝ n` from §3, the rule is

    lambda  ∝  n * C_src

so lambda is not a nuisance constant: it is the share of the objective's
apparent information that is really the sources' own uncertainty.

This explains the two values seen in practice. The simulation used hand-chosen
source SEs about 3.6x larger than the vignette's, which are derived from focei
fits on 180-260 subjects -- 3.6^2 ≈ 13, and 45/3.4 ≈ 13. **It is not structural
misfit**, which was the initial (wrong) guess.

**`lambda < 1` is a real regime.** At scale 0.25 it is 0.85: where the sources
are very precise the naive test is CONSERVATIVE rather than liberal, and the
sandwich SE is smaller than the naive one. The correction goes both ways.

## Settings for an MBMA

1. `covMethod = "r,s"` on every fit. It is what makes the SE valid and it is
   what the LRT and TIC are built from.
2. Compare models with `anova(full, reduced)` or on **TIC**. Never on AIC, never
   on raw dOFV.
3. Nest by fixing the coefficient, not by deleting the term.
4. `stratify` the covariates each source's own model estimated.
5. adgh or admc — adfo and adirmc refuse covariate marginalisation.

## Wired in

`anova()` now reports **`p_eff` beside `Npar`**, and adds a footer when
`p_eff/Npar > 1.5`:

> NOTE: p_eff is 20.9x Npar, so the information equality H = J fails here and
> AIC's 2p penalty is too small. Compare on TIC or on the corrected p, not on
> AIC or on raw dOFV.

AIC is *not* suppressed: users compare AIC across papers, and dropping the
column silently would be its own surprise. The gate is wide (1.5x) because the
ratio was ~21 where it mattered.

## Not covered

- Boundary restrictions (an omega to zero) — `anova()` refuses; use TIC or
  simulate under the null.
- Non-nested comparison — needs Vuong, not provided.
- BIC under model sources — untested here.
- Structural misfit (sources simpler than the fitted model) was deliberately
  excluded from this design. See `FINDING-source-heterogeneity.md`.
