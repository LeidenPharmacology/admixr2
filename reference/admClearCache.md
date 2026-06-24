# Clear the admixr2 model cache

Removes all cached simulation and sensitivity models from both the
session-level in-memory cache and the qs2 disk files written to
[`rxode2::rxTempDir()`](https://nlmixr2.github.io/rxode2/reference/rxTempDir.html).
Call this in long-running sessions to free memory and disk space after
fitting many distinct models.

## Usage

``` r
admClearCache()
```

## Value

Invisibly returns the number of in-memory objects removed.
