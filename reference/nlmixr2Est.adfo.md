# Fit an aggregate data model via First-Order (FO) approximation

Called automatically by
`nlmixr2(model, admData(), est = "adfo", control = adfoControl(...))`.
Not typically called directly.

## Usage

``` r
# S3 method for class 'adfo'
nlmixr2Est(env, ...)
```

## Arguments

- env:

  nlmixr2 environment containing `ui` and `control`.

- ...:

  Unused.

## Value

An `admFit` nlmixr2 fit object.
