# Describe the covariate distribution a study's subjects span

Builds the `cov_dist` a study carries, from what a publication actually
reports: a mean and a standard deviation per covariate, and optionally a
correlation between them. Validation happens here, where an error can
name the covariate, rather than at the first objective evaluation of a
fit.

## Usage

``` r
covDist(..., cor = NULL, joint = NULL, dist = c("normal", "lnorm"))
```

## Arguments

- ...:

  One argument per covariate, named as the model reads it. Each is one
  of:

  `c(mean = , sd = )`

  :   Mean and SD on the covariate's own scale — what a baseline table
      reports. Normal by default; `dist = "lnorm"` moment-matches a
      lognormal margin instead, which is what you want for a covariate
      that must stay positive (see `dist`).

  `c(mu = , sd = )` / `c(meanlog = , sdlog = )`

  :   A normal or lognormal margin given directly on its own scale.

  `c(label = prob, ...)`

  :   A categorical covariate: names are the level labels, values their
      proportions. Levels are coded `0, 1, ...` in the order given.

  a `list(...)`

  :   The canonical form, for anything else — including
      `list(quantile = f)` for an arbitrary margin.

  Alternatively a single data.frame with `covariate`, `mean` and `sd`
  columns (and optionally `dist`): a baseline-characteristics table
  transcribed as-is.

- cor:

  Correlation between the covariates: a scalar for two of them, or a
  correlation matrix (named, in any order). Realised through a Gaussian
  copula on the declared margins.

- joint:

  Dependence a single correlation cannot express: a function receiving
  the matrix of uniforms admixr2 supplies and returning one named column
  per covariate, on each covariate's OWN scale. Overrides `cor`.

  The sampler owns its margins. A copula returns dependence on the
  *copula scale* only, so a sampler that hands those uniforms back
  returns covariate values on `(0, 1)` — every covariate mean 0.50 and
  SD 0.29, whatever margins were declared, with no error anywhere. Push
  them through the quantile function you intend.

  `stratify` needs nothing further: an opaque sampler is banded by
  binning its output rather than by fixing its input uniforms.

- dist:

  Default margin for the `c(mean = , sd = )` form: `"normal"` (default)
  or `"lnorm"`. A per-covariate `dist` wins over it.

  Choose `"lnorm"` when the model needs the covariate to stay positive.
  A normal margin is unbounded below and the quadrature reaches 3.75
  standard deviations, so any covariate with a coefficient of variation
  above about 0.27 (the guard is `mu - 3.75 * sd <= 0`) gets a node at
  or below zero, where an allometric or log term is `NaN`. `covDist()`
  warns when that would happen.

## Value

A validated `cov_dist`, ready to pass to a study. Printing it shows each
covariate's realised mean, SD and type.

## See also

[`covDraw()`](https://leidenpharmacology.github.io/admixr2/reference/covDraw.md)
to see what a specification describes,
[`covStrata()`](https://leidenpharmacology.github.io/admixr2/reference/covStrata.md)
to cut it into strata.

## Examples

``` r
# what a baseline-characteristics table reports
cd <- covDist(WT = c(mean = 72, sd = 16), CRCL = c(mean = 90, sd = 25),
              cor = 0.6)
#> Warning: admixr2: covariate ‘CRCL’ has a NORMAL margin with mean 90 and sd 25, so the quadrature reaches -3.75 -- at or below zero. If the model uses it in a power, allometric or log term that is NaN, and you want dist = "lnorm", which is positive by construction. Ignore this if the covariate is genuinely centred.
cd
#> <covDist> 2 covariate(s)
#>  covariate   type   mean     sd
#>         WT normal 72.002 15.968
#>       CRCL normal 89.993 24.938
#>   dependence: cor = 0.6 
#>     realised cor(WT, CRCL) = +0.599

# a categorical covariate: labels are the levels, values their proportions
covDist(SEX = c(female = 0.55, male = 0.45))
#> <covDist> 1 covariate(s)
#>  covariate        type mean    sd
#>        SEX categorical 0.45 0.498
#>   SEX levels: female=0 (55%), male=1 (45%)

# or transcribe the table itself
covDist(data.frame(covariate = c("WT", "CRCL"),
                   mean = c(72, 90), sd = c(16, 25)))
#> Warning: admixr2: covariate ‘CRCL’ has a NORMAL margin with mean 90 and sd 25, so the quadrature reaches -3.75 -- at or below zero. If the model uses it in a power, allometric or log term that is NaN, and you want dist = "lnorm", which is positive by construction. Ignore this if the covariate is genuinely centred.
#> <covDist> 2 covariate(s)
#>  covariate   type   mean     sd
#>         WT normal 72.002 15.968
#>       CRCL normal 89.989 24.959
#>   dependence: independent 
#>     realised cor(WT, CRCL) = -0.002
```
