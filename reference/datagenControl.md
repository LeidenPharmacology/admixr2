# Control parameters for [`datagen()`](https://leidenpharmacology.github.io/admixr2/reference/datagen.md)

Control parameters for
[`datagen()`](https://leidenpharmacology.github.io/admixr2/reference/datagen.md)

## Usage

``` r
datagenControl(
  method = c("mc", "fo", "gh"),
  n_sim = 5000L,
  n_nodes = 5L,
  sampling = c("sobol", "halton", "torus", "lhs", "rnorm"),
  seed = 12345L,
  cores = 1L,
  return_samples = FALSE
)
```

## Arguments

- method:

  Moment approximation used to generate `E` and `V`: `"mc"` (default)
  draws Monte Carlo samples over the IIV distribution, as in
  `est = "admc"`; `"fo"` uses the deterministic First-Order expansion
  (`mu = f(theta, 0)`, `V = J Omega J' + Sigma`), matching
  `est = "adfo"`; `"gh"` uses deterministic Gauss-Hermite quadrature
  over the random-effects prior, matching `est = "adgh"` – unbiased at
  any IIV magnitude and noise-free. Use `"fo"` or `"gh"` for design
  evaluation where the data-generating and data-analytic models must
  coincide.

- n_sim:

  Number of Monte Carlo samples used to approximate population moments.
  Ignored when `method = "fo"` or `"gh"`.

- n_nodes:

  Number of Gauss-Hermite nodes per eta dimension for `method = "gh"`
  (default 5). Total nodes = `n_nodes^n_eta`. Ignored for `"mc"` and
  `"fo"`.

- sampling:

  Quasi-random sampling method: `"sobol"` (default), `"halton"`,
  `"torus"`, `"lhs"`, or `"rnorm"`. Ignored when `method = "fo"` or
  `"gh"`.

- seed:

  Integer seed. Applied before stochastic methods (`"rnorm"`, `"lhs"`).
  Ignored when `method = "fo"` or `"gh"`.

- cores:

  Number of `rxSolve` threads.

- return_samples:

  Include the raw `n_sim x length(times)` prediction matrix as
  `$samples` in each study's output. No effect when `method = "fo"` or
  `"gh"` (those methods draw no samples).

## Value

A list of class `"datagenControl"`.

## See also

[`datagen()`](https://leidenpharmacology.github.io/admixr2/reference/datagen.md)

## Examples

``` r
ctrl <- datagenControl(n_sim = 2000L)
ctrl$sampling  # "sobol"
#> [1] "sobol"

# Deterministic FO moments for design evaluation:
datagenControl(method = "fo")$method  # "fo"
#> [1] "fo"

# GH quadrature moments (unbiased, noise-free):
datagenControl(method = "gh", n_nodes = 5L)$n_nodes
#> [1] 5
```
