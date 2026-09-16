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
  stratify = NULL,
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

- stratify:

  Band the source into strata, so a covariate it fitted contributes a
  CONTRAST rather than one pooled number. `TRUE` bands over every
  covariate the source's own model ESTIMATED a coefficient for, and
  leaves the rest marginalised — derived from the model, so nothing has
  to be restated. A covariate the model merely READS is not banded:
  weight at a fixed allometric exponent carries no fitted effect to
  recover, and banding on it would buy strata and no evidence. Name a
  covariate explicitly to override that judgement. Available for
  published model sources only. Measured over 720 replicates: coverage
  0.933 with one source banded and 0.925 with all of them, against 0.817
  with none — one banded source is as good as three. See
  [`covStrata()`](https://leidenpharmacology.github.io/admixr2/reference/covStrata.md).

- strata_nodes, range:

  Resolution and the enrolled range for `stratify`.

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

## See also

[`admStudies()`](https://leidenpharmacology.github.io/admixr2/reference/admStudies.md)
to collect several,
[`admPopulation()`](https://leidenpharmacology.github.io/admixr2/reference/admPopulation.md)
for the baseline table.
