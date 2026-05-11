# Dummy data frame for nlmixr2 dispatch

Returns a minimal NONMEM-style data frame that satisfies nlmixr2's data
argument requirement. All `DV` values are `NA` so nlmixr2 adds zero
log(2pi) constants to OBJF, keeping `fit$objective == our -2LL` exactly.

## Usage

``` r
admData()
```

## Value

A data frame with columns `ID`, `TIME`, `DV`, `AMT`, `EVID`, `CMT`.

## Examples

``` r
admData()
#>   ID TIME DV AMT EVID CMT
#> 1  1    0 NA 100  101   1
#> 2  1    1 NA   0    0   2

if (FALSE) { # \dontrun{
# Pass to nlmixr2() as the data argument when using admixr2 estimators:
# fit <- nlmixr2(model, admData(), est = "admc", control = admControl(...))
} # }
```
