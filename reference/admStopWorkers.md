# Stop parallel workers

Stops any worker processes (mirai daemons) started by a parallel-restart
fit (`admControl(workers = N)`). Workers are stopped automatically after
the restart phase completes, so this function is only needed if a fit
was interrupted before cleanup could run.

## Usage

``` r
admStopWorkers()
```

## Value

`NULL`, invisibly.

## Examples

``` r
# Safe to call at any time; no-op if no workers are running
admStopWorkers()
#> No admixr2 workers running.
```
