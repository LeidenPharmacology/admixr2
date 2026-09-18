# Cut a covariate distribution into strata

Shows the strata
[`datagen()`](https://leidenpharmacology.github.io/admixr2/reference/datagen.md)
builds for a study declaring `stratify` — the covariate value each
stratum is pinned at, its effective sample size, and the distribution
the remaining covariates follow *within* it. It is the way to see what a
stratification actually describes before a fit depends on it.

## Usage

``` r
covStrata(
  cov_dist,
  stratify,
  n_nodes = .ADM_STRATA_NODES,
  n = 1,
  n_pool = 32768L,
  cov_range = NULL
)
```

## Arguments

- cov_dist:

  A covariate specification, as given to a study — see
  [`covDraw()`](https://leidenpharmacology.github.io/admixr2/reference/covDraw.md)
  for the full grammar.

- stratify:

  Character vector naming the covariates to stratify on. The rest are
  left in each stratum's own `cov_dist`, to be marginalised over their
  distribution **conditional** on that stratum.

- n_nodes:

  Strata per stratified covariate (default 9). A discrete covariate
  ignores it and is cut at its levels, exactly, with the level
  probabilities as weights.

  For a continuous one it sets how finely the covariate is resolved, and
  the rule depends on what is known about the joint distribution. Where
  the latent structure is known — an explicit `cor`, or no declared
  dependence, which IS independence — each stratum is a Gauss-Hermite
  node, held at a POINT, and the remaining covariates follow their exact
  conditional law within it. An opaque `joint` sampler cannot be
  conditioned, and there the covariate is cut into that many
  equiprobable BINS instead, each described by the pool members that
  landed in it.

  **This is a convergence parameter, not a modelling choice.** The
  well-defined object is the limit as the count grows; any finite value
  approximates it, and under misspecification the answer can jump
  between basins rather than drift. Raise it until the estimates stop
  moving rather than picking a value you like. Cost is `n_nodes^p`
  studies for `p` stratified covariates, so this is cheap in one
  covariate and expensive in three.

  How much the objective depends on it turns on which rule applies.
  Where the latent structure is known — an explicit `cor`, or no
  declared dependence, which is independence — conditioning uses a
  Gauss-Hermite grid and the objective is stable: 0.03 units across
  counts of 5 to 100, so objective, AIC, BIC and likelihood ratios are
  comparable across resolutions.

  A declared discrete covariate no longer forces the slow rule: as long
  as it is latently independent of the others it enumerates exactly at
  its levels, which took a study declaring one sex covariate from 515
  units of movement, still rising at a count of 50, to 0.009.

  An opaque `joint` sampler cannot be conditioned and falls back to
  equiprobable bins, where it is not stable: 51 units at a count of 100
  and still moving, with the estimate itself wandering (0.681 at 5, 9,
  15 and 100; 0.697 at 50). There, raise the count until the estimates
  settle, and note that [`anova()`](https://rdrr.io/r/stats/anova.html)
  refuses to compare two fits built at different ones.

- n:

  Total sample size to divide among the strata. The default of `1`
  returns the raw stratum weights.

- n_pool:

  Minimum size of the deterministic pool the strata are cut from. Each
  stratum's covariate distribution is that pool restricted to its own
  bin, so this sets how finely the within-stratum distribution is
  resolved — and, because no stratum can produce a covariate value
  beyond the pool's extremes, how far into the tails it reaches. The
  pool is grown automatically with the number of strata, and a stratum
  left with too few draws is an error rather than a silently coarse
  answer; `n_pool_cell` in the result reports what each one got.

  It applies only where a pool is used at all: an opaque `joint`
  sampler, or a stratified covariate latently correlated with a
  marginalised one. Where the conditional law is closed form the stratum
  carries the DECLARED specs themselves, no sample stands in for them,
  and `n_pool_cell` is `NA`.

- cov_range:

  Optional named list giving the range each stratified covariate was
  actually ENROLLED over, e.g. `list(WT = c(52, 118))` — publications
  routinely report a min-max or a median with an IQR.

  Without it the strata are cut over the whole declared distribution,
  and a source is credited with evidence at covariate values it never
  sampled. The overstatement is exactly `var(declared) / var(enrolled)`:
  3.4x for a source spanning ±1 SD, 12.4x at ±0.5 SD. Estimates stay
  correct — what inflates is confidence, and it only shows once a second
  source disagrees. Omitting it warns for that reason.

  More strata resolve the covariate range more finely but do not buy
  accuracy, and each one costs a solve: on a matched one-covariate fit
  the coefficient came back at 0.7000 / 0.7002 / 0.7005 / 0.7005 for 3 /
  4 / 10 / 16 strata against a true 0.700, while the fit took 7.2 / 4.1
  / 8.0 / 11.5 seconds. Raise it when the covariate effect is strongly
  nonlinear over the range; the only hard limit is that a stratum must
  keep enough pooled draws to describe itself, which admixr2 checks.

## Value

A list with one element per stratum, each containing `cov` (the value
each stratified covariate is held at — a quadrature node, a level, or a
bin mean), `cov_dist` (the distribution of ALL covariates *within* that
stratum, carrying the full joint dependence), `n` and `weight`.

## Details

Stratify on the covariates a source's own published model conditions on,
and leave the rest to be marginalised. A source that never fitted a
covariate has no contrast in it to report, and evaluating it at nodes
that vary that covariate manufactures a null one — measured, with no
sampling noise, as an attenuation from 0.450 to 0.211 in the fitted
coefficient.

## See also

[`covDraw()`](https://leidenpharmacology.github.io/admixr2/reference/covDraw.md)
to inspect a covariate specification,
[`datagen()`](https://leidenpharmacology.github.io/admixr2/reference/datagen.md)
to generate the strata as studies.

## Examples

``` r
# a source that fitted weight but not renal function, the two correlated
st <- covStrata(list(WT   = list(mean = 72, sd = 15),
                     CRCL = list(mean = 90, sd = 22), cor = 0.7),
                stratify = "WT", n_nodes = 4L, n = 300)
#> Warning: admixr2: conditioning on ‘WT’ over the FULL declared distribution, because no range was given for it. The source is then credited with evidence at covariate values it may never have enrolled -- the overstatement is var(declared)/var(enrolled), which is 3.4x for a source spanning +/-1 SD. Supply the reported span: `range = list(WT = c(min, max))` in admStudy(), or `cov_range` here.
vapply(st, function(s) s$cov$WT, numeric(1))     # where each stratum sits
#> [1] 107.01621  83.12946  60.87054  36.98379
vapply(st, function(s) s$n, numeric(1))          # and how big it is
#> [1]  13.76276 136.23724 136.23724  13.76276

# within a stratum, CRCL follows its CONDITIONAL distribution: a heavier
# stratum has a higher creatinine clearance, which is the whole point
vapply(st, function(s) mean(covDraw(s$cov_dist, n = 2000L)[, "CRCL"]),
       numeric(1))
#> [1] 125.96663 101.44289  78.59041  54.06667
```
