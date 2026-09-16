# Describe the population a study enrolled

Written the way a baseline demographics table reads. Each covariate
takes whichever summary the paper printed — `mean`/`sd`, `median`/`iqr`,
a `cv` as a percent, or a proportion for a binary one — and correlations
are given for the PAIRS that were reported, everything else being
independent.

## Usage

``` r
admPopulation(..., cor = NULL, dist = c("lnorm", "normal"), data = NULL)
```

## Arguments

- ...:

  Named covariates. A continuous one takes a named vector, e.g.
  `WT = c(mean = 75, sd = 16)` or
  `CRCL = c(median = 92, iqr = c(62, 118))`. A binary one takes a single
  named proportion, e.g. `SEX = c(male = 0.55)`, which becomes levels
  `0`/`1` with that probability on `1`.

- cor:

  Correlations between covariate PAIRS, named `A.B`, e.g.
  `cor = c(WT.CRCL = 0.45)`. Pairs not named are independent, so a
  partial table needs no identity padding. A full matrix is accepted
  too.

- dist:

  `"lnorm"` (default) or `"normal"`, for the continuous margins.
  Lognormal is the usual choice for a positive covariate — a normal
  margin wide enough to matter puts mass at or below zero, which is
  `NaN` inside any power or log term.

- data:

  A data frame of individual covariates to derive the table FROM,
  instead of typing it out — a digitised baseline listing, or the cohort
  itself in a simulation study. Each numeric column becomes a margin (a
  0/1 column becomes a proportion, everything else a continuous margin
  matching that column's mean and SD), and every continuous PAIR gets
  its correlation — taken on the LATENT scale, so on the logs for a
  lognormal margin, which is the conversion easiest to get wrong by
  hand. Anything named in `...` or `cor` overrides what the data would
  have given, so a column you would rather state yourself simply gets
  stated.

## Value

A covariate specification, as
[`covDist()`](https://leidenpharmacology.github.io/admixr2/reference/covDist.md)
returns.

## See also

[`admStudy()`](https://leidenpharmacology.github.io/admixr2/reference/admStudy.md),
which takes one;
[`covDraw()`](https://leidenpharmacology.github.io/admixr2/reference/covDraw.md)
to inspect it.
