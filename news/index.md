# Changelog

## admixr2 0.1.0

CRAN release: 2026-06-02

- Initial release.
- Monte Carlo estimator (`est = "admc"`) via
  [`admControl()`](https://leidenpharmacology.github.io/admixr2/reference/admControl.md).
- Iterative Reweighting Monte Carlo estimator (`est = "adirmc"`) via
  [`adirmcControl()`](https://leidenpharmacology.github.io/admixr2/reference/adirmcControl.md).
- Analytical CRN gradient with sensitivity equations (`grad = "sens"`).
- Multi-restart parallelism via `furrr`/`future`.
- Diagnostic plots: observed vs predicted mean/covariance, NLL trace,
  parameter trace.
- Integrates with the nlmixr2/rxode2 ecosystem.
