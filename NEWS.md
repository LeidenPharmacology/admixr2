# admixr2 0.2.0

## New features

* New estimator `est = "adgh"`: deterministic Gauss-Hermite quadrature over the
  random-effects prior, configured via `adghControl()`. The objective is
  noise-free (no Monte Carlo draws), the analytical gradient is exact, and it is
  unbiased at any IIV magnitude. For models with up to ~4 random effects it is
  the fastest exact estimator (#65).
* `datagen()` gains FO-approximated population moments (`method = "fo"`, matching
  `est = "adfo"`) for design evaluation and optimal-design work (#56).
* `adirmcControl(kappa_method = "linearized_gh")`: GH-averaged kappa baseline for
  the IRMC inner loop.
* `admClearCache()` prunes the session-level compiled-model cache (#10).
* Control objects now accept any `nloptr` algorithm; the default is chosen from
  the gradient mode, and `grad`/`algorithm` are reconciled automatically (#70).

## Bug fixes

* Use the ML denominator (`1/n_sim`) consistently in the MC gradient kernels,
  matching the NLL (#48).
* Fix parallel multi-restart dispatch for fork/PSOCK, and fix `adirmc`
  multi-restart (#45).
* Guard non-positive predicted variance in the diagonal-NLL paths (#57).
* Correct the FO diagonal omega gradient scaling, plus assorted plot,
  output-variable detection, caching, and worker-serialization fixes.

## Documentation

* Add Gauss-Hermite sections across the vignettes and fix the pkgdown reference
  index so the documentation site builds (#79).

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
