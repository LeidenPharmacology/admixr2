# admixr2 0.1.0

* Initial release.
* Monte Carlo estimator (`est = "admc"`) via `admControl()`.
* Iterative Reweighting Monte Carlo estimator (`est = "adirmc"`) via `adirmcControl()`.
* Analytical CRN gradient with sensitivity equations (`grad = "sens"`).
* Multi-restart parallelism via `furrr`/`future`.
* Diagnostic plots: observed vs predicted mean/covariance, NLL trace, parameter trace.
* `traceplot()` support: admixr2 fits populate the standard `parHistData` slot,
  so the nlmixr2 `traceplot()` generic works natively (best restart, natural
  scale, no burn-in marker).
* Integrates with the nlmixr2/rxode2 ecosystem.
