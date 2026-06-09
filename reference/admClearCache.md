# Clear the admixr2 model cache

Removes all cached simulation models and pinned foceiModel objects from
the session-level cache. Call this in long-running sessions to free
memory after fitting many distinct models.

## Usage

``` r
admClearCache()
```

## Value

Invisibly returns the number of objects removed.
