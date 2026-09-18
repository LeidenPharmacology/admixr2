# Observed and predicted aggregate moments, tidied

The numbers every diagnostic panel is drawn from, as one row per study
and observation time.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws a fixed
set of panels; this returns the moments behind them so you can draw your
own.

## Usage

``` r
admMoments(fit, n_sim = NULL, seed = 1L, by = c("source", "stratum"))
```

## Arguments

- fit:

  An `admFit` object.

- n_sim, seed:

  Simulation size and seed for the predicted moments. Defaults to the
  fit's own.

- by:

  `"source"` to collapse a conditional source back together, `"stratum"`
  to keep its strata apart.

## Value

A data frame with one row per study and observation time.

## What the columns mean

`obs_sd` and `pred_sd` are the SD of **one observation across subjects**
– not between-subject variability. They carry BSV, residual error, and,
wherever a study marginalises a covariate, the spread that covariate
induces. `struct_sd` is the same quantity **before** residual error is
composed on, so `pred_sd - struct_sd` is what sigma contributes: a
predicted spread that misses the observed one can then be attributed.
`struct_sd` is `NA` for a transforming error model, which moves the mean
as well and leaves the pre-sigma variance on a different scale.

`z` is the mean standardised residual,
`(obs - pred) / sqrt(pred_var / n)`.

## Sources and strata

`stratify` cuts a source into one study per covariate level, because a
covariate a study marginalises over is not identified against a random
effect on the same parameter. Those strata are the unit of the
likelihood, not a unit a reader recognises, so `by = "source"` (the
default) puts them back together: the mean is the n-weighted mean, and
the variance is the law of total variance – within plus **between**, the
between term being the covariate effect conditioning created.

`by = "stratum"` returns them separately, which is what you want when
the between-level contrast is the thing you are looking at.

## Examples

``` r
if (FALSE) { # \dontrun{
m <- admMoments(fit)
library(ggplot2)
ggplot(m, aes(pred_sd, obs_sd, colour = study)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point()
} # }
```
