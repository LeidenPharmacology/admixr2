# Generate aggregate study data from (possibly different) pharmacometric models

Simulates population mean vectors (`E`) and covariance matrices (`V`)
for each study using Monte Carlo integration over the IIV distribution.
Each study may specify its own PK/PD model (as would be the case when
digitising data from several published studies, each fit with a
different structural model). True parameter values are taken from the
[`ini()`](https://nlmixr2.github.io/rxode2/reference/ini.html) block of
each study's model. Each element of the returned list is ready to supply
directly to `admControl(studies = ...)`.

## Usage

``` r
datagen(studies, model = NULL, control = datagenControl())
```

## Arguments

- studies:

  A named list of study specifications. Each element is a list with:

  `model`

  :   An nlmixr2-style model function with
      [`ini()`](https://nlmixr2.github.io/rxode2/reference/ini.html) and
      [`model()`](https://nlmixr2.github.io/rxode2/reference/model.html)
      blocks. Serves as the data-generating model for this study. May
      differ between studies. Can be omitted if a top-level default is
      supplied via the `model` argument.

  `times`

  :   Numeric vector of observation times.

  `ev`

  :   A dosing event table created with
      [`rxode2::et()`](https://nlmixr2.github.io/rxode2/reference/et.html).

  `n`

  :   (Optional) integer sample size; stored as metadata and used when
      supplying the result to
      [`admControl()`](https://leidenpharmacology.github.io/admixr2/reference/admControl.md).

- model:

  Optional default model function used for any study that does not
  supply its own `model` element. At least one of `model` or each
  study's `model` must be non-`NULL`.

- control:

  A
  [`datagenControl`](https://leidenpharmacology.github.io/admixr2/reference/datagenControl.md)
  object.

## Value

A named list with one element per study. Each element contains:

- `E`:

  Population mean vector at `times`.

- `V`:

  Population covariance matrix (`length(times)` x `length(times)`, ML
  denominator `n_sim`). Residual error is added to the diagonal when
  `control$add_residual_error = TRUE`.

- `n`:

  Sample size (`NA_integer_` if not supplied).

- `times`:

  Observation times.

- `ev`:

  Dosing event table.

- `samples`:

  Raw `n_sim x length(times)` prediction matrix (only when
  `control$return_samples = TRUE`).

## Details

Population moments are computed via the same Monte Carlo engine as
`est = "admc"`: \$\$E_t = \bar{f}\_s(\hat\theta_s, \eta_i, t)\$\$
\$\$V\_{ts} = \widehat{\mathrm{Cov}}\_\eta\[f\_{s,t}, f\_{s,s'}\] +
\Sigma_s\$\$ where \\f_s\\ and \\\hat\theta_s\\ are the model and
initial estimates from the
[`ini()`](https://nlmixr2.github.io/rxode2/reference/ini.html) block of
study \\s\\, the sample covariance uses the ML denominator `n_sim`, and
\\\Sigma_s\\ is diagonal with entries determined by that study model's
residual error type (additive, proportional, or log-normal).

Models are compiled and cached on first use (keyed by model expression
digest), so repeated calls or multiple studies sharing the same model
incur only a single compilation.

## See also

[`datagenControl`](https://leidenpharmacology.github.io/admixr2/reference/datagenControl.md),
[`admControl`](https://leidenpharmacology.github.io/admixr2/reference/admControl.md)
