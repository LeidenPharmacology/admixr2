# Fit an aggregate data model via Monte Carlo (admc estimator)

Called automatically by
`nlmixr2(model, admData(), est = "admc", control = admControl(...))`.
Not typically called directly.

## Usage

``` r
# S3 method for class 'admc'
nlmixr2Est(env, ...)
```

## Arguments

- env:

  nlmixr2 environment containing `ui` and `control`.

- ...:

  Unused.

## Value

An `admFit` nlmixr2 fit object.
