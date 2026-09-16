# Collect the studies a meta-analysis draws on

Takes
[`admStudy()`](https://leidenpharmacology.github.io/admixr2/reference/admStudy.md)
objects and names them. Names come from the argument names where given,
so `admStudies(smith2019, jones2021)` labels them by the objects they
were built as.

## Usage

``` r
admStudies(...)
```

## Arguments

- ...:

  [`admStudy()`](https://leidenpharmacology.github.io/admixr2/reference/admStudy.md)
  objects, optionally named.

## Value

An `admStudies` object, to pass as `studies` to any of
[`admControl()`](https://leidenpharmacology.github.io/admixr2/reference/admControl.md),
[`adghControl()`](https://leidenpharmacology.github.io/admixr2/reference/adghControl.md),
[`adfoControl()`](https://leidenpharmacology.github.io/admixr2/reference/adfoControl.md)
or
[`adirmcControl()`](https://leidenpharmacology.github.io/admixr2/reference/adirmcControl.md).

## Details

Nothing is generated here either. The studies are materialised once,
inside the fit, so the whole specification stays cheap to build and
inspect.

## See also

[`admStudy()`](https://leidenpharmacology.github.io/admixr2/reference/admStudy.md),
[`admPopulation()`](https://leidenpharmacology.github.io/admixr2/reference/admPopulation.md).
