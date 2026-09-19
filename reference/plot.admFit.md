# Diagnostic plots for an admixr2 fit

Generates up to five diagnostic panels:

## Usage

``` r
# S3 method for class 'admFit'
plot(
  x,
  which = c("mean", "cov", "covariate", "nll", "par"),
  n_sim = NULL,
  seed = 1L,
  ...
)
```

## Arguments

- x:

  An `admFit` object returned by `nlmixr2()` with `est = "adfo"`,
  `est = "admc"`, `est = "adgh"`, or `est = "adirmc"`.

- which:

  Character vector selecting which panel types to produce. Any subset of
  `c("mean", "cov", "covariate", "nll", "par")`. Defaults to all five.
  Note `"cov"` (predicted vs observed covariance) and `"covariate"`
  (covariate effect) are different panels.

- n_sim:

  Number of MC samples for the final prediction. Defaults to the value
  used during fitting. Only used when `"mean"`, `"cov"` or `"covariate"`
  is in `which`.

- seed:

  Random seed for reproducibility.

- ...:

  Unused.

## Value

A named list of ggplot2 objects, invisibly. Prints each selected
top-level panel. For the `"mean"` and `"cov"` panels the returned list
also contains each sub-panel individually so a single panel (or a few)
can be extracted in code without reprinting the whole grid. Elements can
be pulled out by name – `plot(fit, which = "mean")$mean_study1_pred` or
`plot(fit, which = "cov")$cov_study1_std_resid` – or by position, with
the combined 2x2 grid stored first per study
(`plot(fit, which = "mean")[[1]]` is the full grid, `[1]` the length-1
named sub-list). The sub-panel keys are `<type>_<source>_obs`, `_pred`,
`_resid`, and `_std_resid`; the combined grid stays under
`<type>_<source>`. `<source>` is the study name, and for a source cut
into nodes it is the name of the SOURCE rather than of a node – these
panels are about the paper, so the nodes are put back together first and
there is no `mean_study1_s1`. The extra sub-panel keys are not printed
on their own.

## Details

1.  `"mean"` – Observed vs predicted mean per study (2x2 grid). Upper
    row: observed and predicted mean lines with +/-1 SD ribbon on a
    shared y scale (black throughout). Lower row: raw residual lollipop
    with +/-2 SE band and standardised residual z-scores with +/-1.96
    reference lines.

2.  `"cov"` – Observed vs predicted (co)variance heatmaps per study (2x2
    grid). Upper row shares a common colour scale (blue-white-red).
    Lower row uses distinct diverging scales: residual (red-white-green)
    and standardised residual (gold-white-purple). Significance stars
    overlaid on the standardised residual panel.

3.  `"covariate"` – Two BETWEEN-study covariate panels, one facet per
    covariate any source describes – including one the model does not
    read, which is exactly the case `covariate_resid` exists for.
    `covariate_effect` sweeps each covariate across the range the
    sources between them cover and draws every model parameter that
    reads it, at the fitted thetas, as ONE line per parameter with the
    other covariates held at the pooled centre; the region no source
    sampled is shaded as extrapolation, and a discrete covariate is
    drawn on its levels rather than swept through the values between
    them. `covariate_resid` plots each study's mean standardised
    residual against the covariate value it sits at, with an `lm` trend
    on a continuous axis: a slope there is a mis-specified covariate
    form. Produced only for a fit whose studies declare covariates, so
    it is silently absent otherwise.

4.  `"nll"` – NLL trace per restart over optimizer evaluations. Restarts
    coloured with the Okabe-Ito palette.

5.  `"par"` – Parameter trace per restart on the natural scale (struct
    thetas back-transformed, sigma as SD, omega diagonal as variance
    labelled `V(eta.x)`). Facets ordered as in the model
    [`ini()`](https://nlmixr2.github.io/rxode2/reference/ini.html)
    block. Restarts coloured with the Okabe-Ito palette.

## Why the covariate panels are separate

With aggregate data a covariate effect is identified BETWEEN studies:
the renal exponent in
[`vignette("covariates")`](https://leidenpharmacology.github.io/admixr2/articles/covariates.md)
is recovered from three cohorts sitting at three different
creatinine-clearance medians, none of which fitted a renal term at all.
The `"mean"` and `"cov"` panels are per study and in isolation, so every
source can sit beautifully on its own panel while the relation tying
them together is wrong. That contrast is what `"covariate"` draws.

## Marginal and conditioned sources

Both panels distinguish how each study enters a covariate, because the
two are different statements about the source rather than degrees of the
same one. A **marginalised** covariate – one the estimator integrates
over a declared distribution – is drawn as a round point at the study's
median with a 10th-90th percentile bar and a 2.5th-97.5th whisker: the
source constrains the effect through the spread it induces, and the
whole span is what it speaks for. A **conditional** covariate – cut into
nodes, or `at` or `by`, so the model is solved at one value – is drawn
as a diamond at that value with no spread, because it has none; what
that source reports along it is a relationship, drawn as its own
regression over the range it covers.

A SOURCE IS ONE MARK, not one per stratum. Conditioning cuts a source
into one study per quadrature node, and a node is an internal
discretisation of the very distribution the source already stands for:
drawn straight it became nine small dots for one paper, each carrying a
ninth of its patients, with the outermost node setting the axis. On a
LEVEL axis the positions are values the paper actually reported, so a
source appears once per level.

## Reading `covariate_resid`

A **continuous** covariate is read as a trend: the dashed `lm` across
studies, where a slope is a mis-specified covariate form.

A **discrete** one is read as a contrast instead, and is drawn that way.
Its axis is ticked at its levels and nowhere else, no regression is
fitted through it – a line across the levels of a factor reports as a
slope what is a difference between groups – and the strata `stratify`
cut from one source are joined, because that pairing is the evidence
conditioning creates: the same study, one covariate moved. Several
sources tilting the same way is the mis-specification. A regression over
the pooled cloud of strata cannot show it, since it averages the pairs
away.

With few sources, covariates whose study medians happen to move together
cannot be told apart here – three cohorts whose weights and renal
function both decline will show a slope in both facets whichever one is
mis-specified. The panel localises the problem to a set of covariates,
not to one; separating them needs a source that breaks the pattern.

## Aggregate data slot

Every admixr2 fit also carries the observed and predicted aggregate data
in `fit$env$aggData`, a named list with one entry per study. Each entry
holds the observation `times`, the study `n`, and two moment sets –
`obs` (from the data) and `pred` (predicted at the fitted parameters) –
each a list with the mean vector `E` and the (co)variance matrix `V`:


      fit$env$aggData$study1$obs$E    # observed mean vector
      fit$env$aggData$study1$obs$V    # observed covariance matrix
      fit$env$aggData$study1$pred$E   # predicted mean vector
      fit$env$aggData$study1$pred$V   # predicted covariance matrix

The predicted moments are computed by one MC simulation at the fitted
parameters using the fit's own `n_sim` and a fixed seed. The slot is
absent only when the fit cannot be simulated (no simulation model
available).

It is **per study**, which for a source cut into nodes means per node –
where the `"mean"` and `"cov"` panels are per SOURCE, the nodes put back
together by the mixture law. So the two no longer line up entry for
entry on such a fit, and the panel keys are `mean_<source>` rather than
`mean_<source>_s1`.
[`admMoments()`](https://leidenpharmacology.github.io/admixr2/reference/admMoments.md)
returns whichever of the two you want, tidied: `by = "source"` matches
the panels, `by = "stratum"` matches this slot.

## nlmixr2 `traceplot()`

admixr2 fits also plug into the nlmixr2 `traceplot()` generic. During
fitting the parameter iteration history of the best restart is stored on
the fit in the standard `parHistData` slot (natural scale), so
`traceplot(fit)` produces the familiar per-parameter, free-y facetted
trace used elsewhere in the nlmixr2 ecosystem. There is no burn-in
marker (admixr2 records optimizer evaluations, not SAEM iterations), and
only the best restart is shown – the per-restart overlay and the NLL
trace remain available via `plot(fit, which = c("par", "nll"))`. The
trace stores only improving evaluations (steps that lowered the best
NLL), so the `iter` axis indexes those improvement steps rather than raw
optimizer iterations.

## Examples

``` r
# \donttest{
library(rxode2)
library(nlmixr2)

data("examplomycin")
obs    <- examplomycin[examplomycin$EVID == 0, ]
obs    <- obs[order(obs$ID, obs$TIME), ]
times  <- sort(unique(obs$TIME))
ids    <- unique(obs$ID)
dv_mat <- do.call(rbind, lapply(ids, function(i) {
  sub <- obs[obs$ID == i, ]; sub$DV[order(sub$TIME)]
}))
E <- colMeans(dv_mat)
V <- cov.wt(dv_mat, method = "ML")$cov

pk_model <- function() {
  ini({
    tcl <- log(5); tv <- log(30)
    prop.sd <- c(0, 0.2)
    eta.cl ~ 0.09; eta.v ~ 0.04
  })
  model({
    cl <- exp(tcl + eta.cl)
    v  <- exp(tv  + eta.v)
    d/dt(central) <- -(cl/v) * central
    cp <- central / v
    cp ~ prop(prop.sd)
  })
}

fit <- nlmixr2(
  pk_model, admData(), est = "adfo",
  control = adfoControl(
    studies = list(study1 = list(E = E, V = V, n = length(ids),
                                 times = times, ev = et(amt = 100))),
    maxeval = 100L
  )
)
#>  
#>  
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> === admixr2: Aggregate Data Modeling (FO) ===
#>   Obs units: 1 | Params: 5 | Cores: 2 | Grad: Analytical | Restarts: 1
#> +----------+----------+----------+----------+----------+----------+----------+
#> |          |     -2LL |      tcl |       tv |  prop.sd |   eta.cl |    eta.v |
#> +----------+----------+----------+----------+----------+----------+----------+
#> | 0010     |  1768.15 |    4.967 |    29.88 |   0.2587 |   0.0888 |  0.04603 |
#> | 0020     |   862.47 |    6.391 |    37.74 |   0.3864 |  0.08003 |   0.0422 |
#> | 0029 ✓   |   861.90 |    6.384 |    38.03 |     0.39 |  0.08051 |  0.04074 |
#> | 0.5 sec  |          |          |          |          |          |          |
#>   Computing covariance (R method, Analytical-Hessian, sandwich, 6 gradient evaluations)
#> → compress origData in nlmixr2 object, save 1160
#>  
#>  
plot(fit)




# }
```
