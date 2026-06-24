# Generate aggregate study data from (possibly different) pharmacometric models

Generates population mean vectors (`E`) and covariance matrices (`V`)
for each study by integrating over the IIV distribution – either by
Monte Carlo (the default) or by a deterministic First-Order expansion
(`method = "fo"`, see
[`datagenControl()`](https://leidenpharmacology.github.io/admixr2/reference/datagenControl.md)).
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
  [`datagenControl()`](https://leidenpharmacology.github.io/admixr2/reference/datagenControl.md)
  object.

## Value

A named list with one element per study. Each element contains:

- `E`:

  Population mean vector at `times`.

- `V`:

  Population covariance matrix (`length(times)` x `length(times)`; ML
  denominator `n_sim` for `method = "mc"`, the analytical FO covariance
  for `method = "fo"`, or the GH weighted covariance for
  `method = "gh"`). The diagonal carries the model's residual-error
  variance; to generate residual-free (IIV-only) moments, omit the error
  term from the model.

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

With `control = datagenControl(method = "mc")` (the default) population
moments are computed via the same Monte Carlo engine as `est = "admc"`:
\$\$E_t = \bar{f}\_s(\hat\theta_s, \eta_i, t)\$\$ \$\$V\_{ts} =
\widehat{\mathrm{Cov}}\_\eta\[f\_{s,t}, f\_{s,s'}\] + \Sigma_s\$\$ where
\\f_s\\ and \\\hat\theta_s\\ are the model and initial estimates from
the [`ini()`](https://nlmixr2.github.io/rxode2/reference/ini.html) block
of study \\s\\, the sample covariance uses the ML denominator `n_sim`,
and \\\Sigma_s\\ is diagonal with entries determined by that study
model's residual error type (additive, proportional, or log-normal).

With `method = "fo"` the moments are instead the deterministic
First-Order expansion used by `est = "adfo"`: \$\$E = f_s(\hat\theta_s,
0)\$\$ \$\$V = J \Omega_s J^\top + \Sigma_s, \quad J\_{tj} = \partial
f\_{s,t}/\partial \eta_j \|\_{\eta = 0}\$\$ with the Jacobian \\J\\
obtained from the sensitivity model (or finite differences if that is
unavailable). This is the natural choice for design evaluation and
optimal design: the moments are fast and reproducible, and because the
data-generating and data-analytic models coincide, the FO Hessian of the
log-likelihood (the expected information matrix) is evaluated at the
true maximum rather than at a point that is not an MLE of the generated
data. Note `est = "adfo"` always adds \\\Sigma\\ to its predicted
covariance, so for a consistent FIM keep the residual error in the
generating model; omit it only when residual-free (IIV-only) moments are
genuinely what you want.

With `method = "gh"` the moments are computed by deterministic
Gauss-Hermite quadrature over the random-effects prior \\\eta \sim N(0,
\Omega)\\: \$\$E = \sum_q w_q f(\hat\theta, \eta_q), \quad V = \sum_q
w_q (f_q - E)(f_q - E)^\top + \Sigma\$\$ where \\(\eta_q, w_q)\\ are the
Cholesky-scaled tensor-product GH nodes and weights. Unlike FO this is
unbiased at any IIV magnitude; unlike MC the result is noise-free and
exactly reproducible. Matching the moments of `est = "adgh"` makes
`method = "gh"` the natural choice for optimal design with that
estimator.

Models are compiled and cached on first use (keyed by model expression
digest), so repeated calls or multiple studies sharing the same model
incur only a single compilation.

## See also

[`datagenControl()`](https://leidenpharmacology.github.io/admixr2/reference/datagenControl.md),
[`admControl()`](https://leidenpharmacology.github.io/admixr2/reference/admControl.md)

## Examples

``` r
# \donttest{
library(rxode2)

pk_model <- function() {
  ini({
    tcl <- log(5); tv <- log(30)
    prop.sd <- c(0, 0.2)
    eta.cl ~ 0.09; eta.v ~ 0.04
  })
  model({
    cl <- exp(tcl + eta.cl)
    v  <- exp(tv  + eta.v)
    d/dt(central) <- -(cl/v) * central
    cp <- central / v
    cp ~ prop(prop.sd)
  })
}

study_data <- datagen(
  studies = list(
    study1 = list(times = c(1, 2, 4, 8, 12, 24),
                  ev = rxode2::et(amt = 100), n = 200L)
  ),
  model   = pk_model,
  control = datagenControl(n_sim = 2000L)
)
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments

# E and V plug directly into admControl(studies = ...)
round(study_data$study1$E, 2)
#>    1    2    4    8   12   24 
#> 2.83 2.37 1.68 0.88 0.48 0.10 
# }
```
