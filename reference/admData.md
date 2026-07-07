# Dummy data frame for nlmixr2 dispatch

Returns a minimal NONMEM-style data frame that satisfies nlmixr2's data
argument requirement. All `DV` values are `NA` so nlmixr2 adds zero
log(2pi) constants to OBJF, keeping `fit$objective == our -2LL` exactly.

## Usage

``` r
admData(outputs = NULL)
```

## Arguments

- outputs:

  Optional character vector of observed output (endpoint) names for a
  multi-compartment model with several prediction lines (e.g.
  `c("cp", "cCSF")`). One `NA` observation row is emitted per endpoint,
  keyed by name in the `CMT` column, so nlmixr2's data translation
  recognises every endpoint. When `NULL` (default) the single-endpoint
  dummy frame is returned unchanged.

## Value

A data frame with columns `ID`, `TIME`, `DV`, `AMT`, `EVID`, `CMT`.

## Examples

``` r
admData()
#>   ID TIME DV AMT EVID CMT
#> 1  1    0 NA 100  101   1
#> 2  1    1 NA   0    0   2
admData(c("cp", "cCSF"))
#>   ID TIME DV AMT EVID CMT DVID
#> 1  1    0 NA 100    1   1 <NA>
#> 2  1    1  1   0    0  NA   cp
#> 3  1    2  1   0    0  NA cCSF
```
