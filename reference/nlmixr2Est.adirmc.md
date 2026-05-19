# Fit an aggregate data model via Iterative Reweighting MC (adirmc estimator)

Called automatically by
`nlmixr2(model, admData(), est = "adirmc", control = adirmcControl(...))`.
Not typically called directly.

## Usage

``` r
# S3 method for class 'adirmc'
nlmixr2Est(env, ...)
```

## Arguments

- env:

  nlmixr2 environment containing `ui` and `control`.

- ...:

  Unused.

## Value

An `admFit` nlmixr2 fit object.
