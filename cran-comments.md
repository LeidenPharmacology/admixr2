# CRAN submission comments — admixr2 0.1.0

## Second resubmission

This is a second resubmission addressing the following point raised by the
reviewer after the first resubmission:

- **`\dontrun{}`**: `plot.admFit()` and `print.admFit()` examples still used
  `\dontrun{}`; replaced with `\donttest{}`.

Additional improvements made proactively:

- **Missing examples**: added `@examples` to `adfoControl()`, `datagen()`, and
  `datagenControl()`. Fast constructor calls run unconditionally; fitting
  examples are wrapped in `\donttest{}`.

## First resubmission (addressed points)

- **Title**: removed redundant "in R".
- **Description**: expanded the PK/PD acronym to
  "pharmacokinetic/pharmacodynamic (PK/PD)"; added DOI references for
  Välitalo (2021) <doi:10.1007/s10928-021-09760-1> and
  van de Beek et al. (2025) <doi:10.1007/s10928-025-10011-w>.
- **`\dontrun{}`**: replaced with `\donttest{}` in `admControl()` and
  `adirmcControl()` examples. Blocks containing only commented-out code were
  removed entirely (`admData()`, `nlmixr2Est.admc()`, `nlmixr2Est.adirmc()`).
- **Commented-out code in examples**: removed from all affected files
  (`admData.Rd`, `admControl.Rd`, `adirmcControl.Rd`).

## Test environments

- Local: Windows 11, R 4.4.x
- GitHub Actions: ubuntu-latest (R release, R devel), windows-latest (R release), macOS-latest (R release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Notes on dependencies

- `nlmixr2est` and `rxode2` are on CRAN and provide the model specification
  and ODE-solving infrastructure this package integrates with.
- `qs2` is used for caching compiled rxode2 model objects between runs.
- `furrr` and `future` (Suggests) are used only when `workers > 1` is set in
  the control object; they are loaded conditionally via `requireNamespace()`.
- `patchwork` (Suggests) is loaded conditionally via `requireNamespace()` for
  optional 2x2 diagnostic plot layouts.
- `knitr` (Suggests) is used for a `knit_print.admFit` S3 method registered
  at package load time.
- `expm` (Suggests) is a fallback for Hessian inversion when Cholesky and
  `solve()` both fail; loaded conditionally via `requireNamespace()`.

## Vignettes

One vignette ("Getting started with admixr2") is included and fully executed
during package building. Five additional vignettes (diagnostic plots, multiple
studies, estimator comparison, advanced usage, data generation) are excluded
from the tarball via `.Rbuildignore` because they each require long-running
fits. They are available on the package website at
<https://leidenpharmacology.github.io/admixr2>.

## Integration tests

Integration tests (those requiring `rxode2` model compilation and ODE solving)
are skipped on CRAN via `skip_on_cran()`. All remaining tests pass with 0
failures.
