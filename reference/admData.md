# Dummy data frame for nlmixr2 dispatch

Returns a minimal NONMEM-style data frame that satisfies nlmixr2's data
argument requirement. The single observation row carries a non-`NA`
placeholder `DV` (`1`); the dose row keeps `DV = NA`. The placeholder is
purely for dispatch and output construction and never enters the
reported objective – each estimator overwrites `fit$env$objective` (and
OBJF/logLik/AIC/BIC) with its own aggregate -2LL. A non-`NA` observation
is required because nlmixr2's post-fit output construction
(`nlmixr2CreateOutputFromUi`) solves the model over this frame, and
rxode2's event-table translation rejects a dataset with no non-`NA`
observation rows ("no rows in event table or input data"), the same
reason the multi-endpoint frame below uses a placeholder `DV`.

## Usage

``` r
admData(outputs = NULL)
```

## Arguments

- outputs:

  Optional character vector of observed output (endpoint) names for a
  multi-compartment model with several prediction lines (e.g.
  `c("cp", "cCSF")`). One observation row is emitted per endpoint, keyed
  by name in the `DVID` column (nlmixr2's endpoint identifier; `CMT` is
  `NA` on those rows), so nlmixr2's data translation recognises every
  endpoint. These rows carry a non-`NA` placeholder `DV` (`1`) because
  nlmixr2's multi-endpoint translator rejects an all-`NA`-DV dataset;
  the placeholder is purely for dispatch and never enters the reported
  objective (each estimator overwrites it with its own aggregate -2LL).
  When `NULL` (default) the single-endpoint dummy frame is returned
  unchanged.

## Value

A data frame with columns `ID`, `TIME`, `DV`, `AMT`, `EVID`, `CMT`
(single-endpoint), plus a `DVID` endpoint column when `outputs` is
given.

## Examples

``` r
admData()
#>   ID TIME DV AMT EVID CMT
#> 1  1    0 NA 100  101   1
#> 2  1    1  1   0    0   2
admData(c("cp", "cCSF"))
#>   ID TIME DV AMT EVID CMT DVID
#> 1  1    0 NA 100    1   1 <NA>
#> 2  1    1  1   0    0  NA   cp
#> 3  1    2  1   0    0  NA cCSF
```
