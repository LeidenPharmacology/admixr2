# Fit an aggregate data model via Gauss-Hermite quadrature

Called automatically by
`nlmixr2(model, admData(), est = "adgh", control = adghControl(...))`.
Not typically called directly.

## Usage

``` r
# S3 method for class 'adgh'
nlmixr2Est(env, ...)
```

## Arguments

- env:

  nlmixr2 environment containing `ui` and `control`.

- ...:

  Unused.

## Value

An `admFit` nlmixr2 fit object.
