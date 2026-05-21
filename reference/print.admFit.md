# Print method for admFit objects

Delegates to `print.nlmixr2FitCore` for the standard nlmixr2 coloured
output. `admFit` class is kept on the object during the call so that
`head.admFit` intercepts any `head(fit)` calls that arise in the paged-
output path (R Markdown / notebooks), preventing the
`[.data.frame(.subset2(env, integer))` crash that occurs when an
environment-backed fit is subscripted like a plain list.

## Usage

``` r
# S3 method for class 'admFit'
print(x, ...)
```

## Arguments

- x:

  An `admFit` object.

- ...:

  Passed to `print.nlmixr2FitCore`.

## Value

`x`, invisibly.

## Examples

``` r
# \donttest{
fit <- nlmixr2(model, admData(), est = "admc", control = admControl(...))
#>  
#>  
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> Error: rxode2 model function requires one 'model({})' block
print(fit)   # or just: fit
#> Error: object 'fit' not found
# }
```
