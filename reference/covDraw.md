# Draw covariate values from a `cov_dist` specification

Simulates the covariate values admixr2 itself would use for `n`
subjects, from the same specification a study carries in `cov_dist`. It
is the way to see what a specification actually describes — to check
that a reported mean and standard deviation were transcribed correctly,
that a correlation points the way round you meant, or that a copula
sampler returns what you think it does — before a fit depends on it.

## Usage

``` r
covDraw(cov_dist, n = 1000L, n_eta = 0L)
```

## Arguments

- cov_dist:

  A covariate specification, as given to a study. Each element names a
  covariate and describes its distribution, either as `mean` and `sd` on
  the covariate's own scale — NORMAL by default, so pass
  `dist = "lnorm"` for the positive margin an allometric term needs
  (only
  [`covDist()`](https://leidenpharmacology.github.io/admixr2/reference/covDist.md)
  defaults to lognormal) — or as a `quantile` function, or as `values`
  (with optional `probs`) for a discrete covariate. A `cor` entry — a
  scalar for two covariates, or a correlation matrix — links them
  through a Gaussian copula. A `joint` function takes the matrix of
  uniforms and returns one column per covariate, which is how an
  arbitrary copula, an R-vine included, is supplied.

- n:

  Number of subjects to draw.

- n_eta:

  Number of random effects in the model the specification belongs to.
  The draws are deterministic and come from Sobol dimensions after the
  random effects', so passing the model's own `n_eta` reproduces exactly
  the values a fit would use. The default of `0` is right for inspecting
  a specification on its own.

## Value

A numeric matrix with `n` rows and one named column per covariate.

## Examples

``` r
# a baseline-characteristics table, transcribed directly
X <- covDraw(list(WT   = list(mean = 72, sd = 16),
                  CRCL = list(mean = 90, sd = 25),
                  cor  = 0.6), n = 500)
colMeans(X)
#>       WT     CRCL 
#> 72.05359 90.14087 
stats::cor(log(X[, "WT"]), log(X[, "CRCL"]))
#> [1] 0.5752932
```
