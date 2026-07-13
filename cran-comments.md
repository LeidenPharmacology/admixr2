# CRAN submission comments — admixr2 0.2.0

## Summary

This is a minor update from the current CRAN version (0.1.0). Notable changes:

- New deterministic Gauss-Hermite quadrature estimator (`est = "adgh"`), a
  noise-free exact estimator configured via `adghControl()`.
- Fix for an infinite recursion ("evaluation nested too deeply" / "node stack
  overflow") that could abort the first model fit of a session when a covariance
  matrix was requested (`covMethod = "r"`).
- Minimum versions declared for the imported 'rxode2' (>= 5.1.2) and
  'nlmixr2est' (>= 6.0.1), and for suggested 'nlmixr2' (>= 5.0.0).

See NEWS.md for the full list of changes.

## Test environments

- Local: Windows 11, R 4.5.3
- GitHub Actions: ubuntu-latest (R release, R devel), windows-latest (R release),
  macOS-latest (R release)
- win-builder: R-devel and R-release

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no reverse dependencies on CRAN.

## Notes on dependencies

- `nlmixr2est` and `rxode2` are on CRAN and provide the model specification
  and ODE-solving infrastructure this package integrates with. Minimum versions
  are declared in `Imports`.
- `qs2` is used for caching compiled rxode2 model objects between runs.
- `nlmixr2` (Suggests) is used in examples and tests; a minimum version is
  declared.
- `mirai` (Suggests) is used only when `workers > 1` is set in the control
  object; it is loaded conditionally via `requireNamespace()`.
- `patchwork` (Suggests) is loaded conditionally via `requireNamespace()` for
  optional 2x2 diagnostic plot layouts.
- `knitr` (Suggests) is used for a `knit_print.admFit` S3 method registered
  at package load time.
- `expm` (Suggests) is a fallback for Hessian inversion when Cholesky and
  `solve()` both fail; loaded conditionally via `requireNamespace()`.

## Vignettes

One vignette ("Getting started with admixr2") is included and fully executed
during package building. The remaining vignettes (diagnostic plots, multiple
studies, estimator comparison, advanced usage, data generation) are excluded
from the tarball via `.Rbuildignore` because they each require long-running
fits. They are available on the package website at
<https://leidenpharmacology.github.io/admixr2>.

## Integration tests

Integration tests (those requiring `rxode2` model compilation and ODE solving)
are skipped on CRAN via `skip_on_cran()`. All remaining tests pass with 0
failures.
