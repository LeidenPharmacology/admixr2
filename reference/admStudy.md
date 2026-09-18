# Write down a published study

One study, transcribed from one paper. It holds what the paper reported
and nothing else, is checked as you build it, and prints so the
transcription can be read back against the source before anything is
fitted.

## Usage

``` r
admStudy(
  model = NULL,
  est = NULL,
  E = NULL,
  V = NULL,
  sd = NULL,
  sem = NULL,
  n = NULL,
  times = NULL,
  dose = NULL,
  ev = NULL,
  population = NULL,
  at = NULL,
  by = NULL,
  strata_nodes = NULL,
  range = NULL,
  label = NULL,
  v_denom = NULL
)
```

## Arguments

- model:

  The published model, as an nlmixr2-style function (or a parsed
  `rxUi`). Omit for a digitised study.

- est:

  Named vector of the paper's parameter estimates, on the scale the
  model's [`ini()`](https://nlmixr2.github.io/rxode2/reference/ini.html)
  is written on — so `tcl = log(5.2)` for a clearance reported as 5.2.
  Omitted parameters keep the model's own
  [`ini()`](https://nlmixr2.github.io/rxode2/reference/ini.html) value.

- E, V, sd, sem:

  Digitised aggregate data: the reported mean profile `E`, with its
  spread as `sd` (per timepoint), `sem` (converted using `n`), or a full
  covariance `V`.

- n:

  Number of subjects the study reports on. Taken from `population` when
  that is the cohort's data frame.

- times:

  Observation times.

- dose:

  Dose amount, as shorthand for a single-dose `ev`.

- ev:

  A dosing event table from
  [`rxode2::et()`](https://nlmixr2.github.io/rxode2/reference/et.html),
  for anything `dose` cannot express.

- population:

  The population the study enrolled — see
  [`admPopulation()`](https://leidenpharmacology.github.io/admixr2/reference/admPopulation.md).
  Covariates the model reads are integrated over it. A **data frame** of
  enrolled subjects is accepted directly and read into a table, taking
  `n` from its row count. EVERY column that is not a reserved data
  column (`ID`, `TIME`, `DV`, `AMT`, ...) becomes a covariate, including
  ones this model does not read — a covariate another source fitted is
  still evidence here, and each one adds a quadrature axis, so drop the
  columns you do not want integrated over.

- at:

  Named list pinning a covariate at a value, for a study reported in one
  subgroup, e.g. `at = list(SEX = 1)`. A pinned covariate must be
  omitted from `population`.

- by:

  Covariate the paper reports results SEPARATELY by, e.g. `by = "SEX"`.
  Expands a published model into one study per level; digitised subgroup
  profiles must be supplied as separate studies.

- strata_nodes:

  Quadrature resolution for a conditional CONTINUOUS covariate, i.e. how
  many nodes it is cut into. Default 9. Raise it for precision, at a
  multiplicative cost: a source conditional on two continuous covariates
  becomes `strata_nodes^2` studies, and each is solved. Discrete
  covariates are unaffected — their levels are exact and cost one
  stratum each. Note that the objective VALUE depends on this, so two
  fits are only comparable at the same resolution; admixr2 records it
  per study and [`anova()`](https://rdrr.io/r/stats/anova.html) checks
  it.

- range:

  Optional named list giving the covariate span the source ENROLLED,
  e.g. `range = list(WT = c(52, 118))`. The declared distribution is
  truncated to it, whether that covariate ends up conditional (strata
  are then cut from the truncated margin) or marginal (the truncated
  margin is integrated over). Needed when that distribution is wider
  than the enrolment — a `mean +/- SD` transcribed from a baseline table
  has tails where nobody was — and *not* wanted when `population` is the
  patients themselves, since the margins are then fitted to exactly who
  was enrolled and truncating would score the paper against a
  sub-population of its own. An unnamed value is allowed only when the
  source itself is conditional on exactly one covariate — which
  covariate an unnamed range belongs to is a question about the source,
  so the answer does not change with the model being fitted to it.

- label:

  Optional display name; otherwise taken from the argument name in
  [`admStudies()`](https://leidenpharmacology.github.io/admixr2/reference/admStudies.md).

- v_denom:

  Which denominator the supplied spread uses: `"unbiased"` (`n - 1`) or
  `"ml"` (`n`). Usually leave it unset — the currency you wrote the
  study in already says which it is, and `admStudy()` resolves it and
  shows the answer when the study is printed.

  `sd`/`sem` default to `"unbiased"`, because a **published** spread is
  the `n - 1` one; `V` defaults to `"ml"`, because handing over a
  covariance matrix is the deliberate act of someone who computed it,
  and the likelihood wants the ML denominator. Set it explicitly when
  your source breaks that pattern — a `sd` you computed yourself with
  the ML denominator, say. A published `model` is always `"ml"`: its
  moments are generated at that denominator, so `"unbiased"` is refused
  rather than applied to a `V` that is already right.

## Value

An `admStudy` object.

## Details

A study contributes in one of two currencies, and this takes either:

- **a published model** — `model`, with the paper's parameter table as
  `est`;

- **digitised aggregate data** — `E` with `sd` (or `sem`, or `V`).

A study generated from a model is not a sample: its mean and covariance
are exact functions of that model's parameters, and `n` sets its
RELATIVE WEIGHT against the other studies rather than its precision. No
standard error is reported for a fit that includes one — see
[`admStudies()`](https://leidenpharmacology.github.io/admixr2/reference/admStudies.md).

Nothing is solved here. The study is generated when it reaches the fit,
so building one is cheap and a mistake surfaces on
[`print()`](https://rdrr.io/r/base/print.html) rather than after a long
run.

## Conditional and marginal covariates are derived, not declared

Whether a covariate is **conditional** or **marginal** for a source is a
property of that source's own model, so admixr2 works it out and there
is nothing to set.

A covariate is **conditional** when the source's model ESTIMATED its
coefficient. That paper reports a contrast along it, so the source is
cut into nodes and the contrast is carried into the fit. A covariate the
model merely READS is not conditional — weight at a fixed allometric
exponent carries no fitted effect to recover, and conditioning on it
would buy nodes and no evidence.

A covariate is **marginal** when the model did not estimate it. There is
no contrast in that paper to condition on, and splitting it on a
covariate its model never saw would manufacture evidence, so admixr2
integrates over the population instead. That covers both the covariate
the model never mentions and the one it reads at an asserted
coefficient; [`print()`](https://rdrr.io/r/base/print.html) names them
separately, because only the second is easy to mistake for conditional.

Which covariates a source is conditional on comes from that source.
Which of them admixr2 has to cut into nodes is narrowed to the ones the
**analysis** model reads, and that narrowing is exact: if the model's
prediction does not move across a source's nodes, the mixture those
nodes collapse to is a sufficient statistic for it, so the objective is
unchanged to eight decimal places while the study count falls. Measured
on a source conditional on two covariates with an analysis model reading
one, collapsing the unread one moved the objective by 0.00004 and
collapsing the read one by 73.6.

This matters and is why it is not left to the caller: a covariate every
source marginalises is not identified against a random effect on the
same parameter. Measured over 720 replicates, coverage was 0.933 with
one source conditional and 0.925 with all of them conditional, against
**0.817 with none**.

A source given as digitised `E`/`V` has no model to derive from and is
marginal over whatever it declares. Use one `admStudy(..., at = ...)`
per reported subgroup to enter a paper that published by subgroup.

`stratify` was removed;
[`covStrata()`](https://leidenpharmacology.github.io/admixr2/reference/covStrata.md)
still takes it for working with a covariate distribution directly.
`strata_nodes` and `range` remain: the first is a precision setting and
the second is a fact about the source, neither of which is a statement
about which covariates are conditional.

## See also

[`admStudies()`](https://leidenpharmacology.github.io/admixr2/reference/admStudies.md)
to collect several,
[`admPopulation()`](https://leidenpharmacology.github.io/admixr2/reference/admPopulation.md)
for the baseline table.
