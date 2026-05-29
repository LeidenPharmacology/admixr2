# Control parameters for [`datagen()`](https://leidenpharmacology.github.io/admixr2/reference/datagen.md)

Control parameters for
[`datagen()`](https://leidenpharmacology.github.io/admixr2/reference/datagen.md)

## Usage

``` r
datagenControl(
  n_sim = 5000L,
  sampling = c("sobol", "halton", "torus", "lhs", "rnorm"),
  seed = 12345L,
  cores = 1L,
  add_residual_error = TRUE,
  return_samples = FALSE
)
```

## Arguments

- n_sim:

  Number of Monte Carlo samples used to approximate population moments.

- sampling:

  Quasi-random sampling method: `"sobol"` (default), `"halton"`,
  `"torus"`, `"lhs"`, or `"rnorm"`.

- seed:

  Integer seed. Applied before stochastic methods (`"rnorm"`, `"lhs"`).

- cores:

  Number of `rxSolve` threads.

- add_residual_error:

  Add residual-error variance to the diagonal of `V` (`TRUE` by
  default), matching the admixr2 NLL convention.

- return_samples:

  Include the raw `n_sim x length(times)` prediction matrix as
  `$samples` in each study's output.

## Value

A list of class `"datagenControl"`.

## See also

[`datagen()`](https://leidenpharmacology.github.io/admixr2/reference/datagen.md)

## Examples

``` r
ctrl <- datagenControl(n_sim = 2000L)
ctrl$sampling  # "sobol"
#> [1] "sobol"
```
