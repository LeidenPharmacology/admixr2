# Changelog

## admixr2 (development version)

### New features

- **Parallel restarts now run on `mirai` daemons.** `workers > 1` starts
  a pool of background R processes instead of dispatching through
  `future`/`furrr`. This replaces the previous fork (Unix/macOS) vs
  PSOCK (Windows/RStudio) split with a single code path that behaves
  identically on every platform, and the pool lives on its own mirai
  compute profile so it never disturbs daemons the user has set up for
  their own code. `furrr` and `future` are no longer used; `mirai` moves
  into `Suggests`. Workers are still stopped automatically after the
  restart phase (and now also on error/interrupt, via
  [`on.exit()`](https://rdrr.io/r/base/on.exit.html)), so all cores are
  free for the covariance step;
  [`admStopWorkers()`](https://leidenpharmacology.github.io/admixr2/reference/admStopWorkers.md)
  remains available.

- **Multi-compartment fitting (multiple observed outputs).** A study may
  now observe several model outputs at once (e.g. plasma and brain/CSF)
  via an `observations` list – one entry per observed output with its
  own `output`, `times`, `E` and `V`. Two modes
  ([\#85](https://github.com/LeidenPharmacology/admixr2/issues/85)):

  - *Independent* – each output has its own `n`/`ev` (separate
    experiments, e.g. literature meta-analysis); the aggregate `-2LL` is
    the sum of the per-output likelihood blocks. Fit with full
    analytical / sensitivity gradients.
  - *Joint (same subjects)* – outputs measured on the same subjects,
    with a shared `n`/`ev` and a joint covariance given either as a
    study-level full `V` or as per-output marginal `V` plus a `cross`
    list of cross-covariance blocks. Scored by a single MVN over the
    stacked vector with shared random effects and the full
    **analytical** gradient in all three estimators (any number of
    compartments; the assembled joint covariance is checked for
    positive-definiteness).

  Supported by `est = "admc"`, `"adfo"` and `"adgh"`;
  [`datagen()`](https://leidenpharmacology.github.io/admixr2/reference/datagen.md)
  generates multi-output aggregate data and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) renders one
  panel set per compartment. Pass the endpoint names to
  [`admData()`](https://leidenpharmacology.github.io/admixr2/reference/admData.md),
  e.g. `admData(c("cp", "cCSF"))`. `est = "adirmc"` does not support
  multiple observed outputs.

- **`nDisplayProgress` control argument** for every estimator
  ([`admControl()`](https://leidenpharmacology.github.io/admixr2/reference/admControl.md),
  [`adfoControl()`](https://leidenpharmacology.github.io/admixr2/reference/adfoControl.md),
  [`adghControl()`](https://leidenpharmacology.github.io/admixr2/reference/adghControl.md),
  [`adirmcControl()`](https://leidenpharmacology.github.io/admixr2/reference/adirmcControl.md)),
  passed through to the
  [`rxSolve()`](https://nlmixr2.github.io/rxode2/reference/rxSolve.html)
  calls that drive fitting. It sets how many subjects a single solve
  must exceed before the solver shows its text progress bar. The default
  (`.Machine$integer.max`) keeps the bar off, so it no longer leaks into
  scripts, logs or rendered vignettes; lower it (e.g. `1000L`) to watch
  progress during long interactive fits.

### Bug fixes

- The solver progress bar no longer appears during covariance/gradient
  batches. Most internal
  [`rxSolve()`](https://nlmixr2.github.io/rxode2/reference/rxSolve.html)
  calls already suppressed it, but the covariance and batched-gradient
  solves in `admc` hard-coded a low `nDisplayProgress` (1000), so the
  bar printed once a chunk exceeded 1000 solves. All solves now honour
  the new `nDisplayProgress` control argument (default off).

- Hard-coded numeric constants in a model’s `model({})` block (e.g. a
  fixed brain volume `vb <- 5`, common in PBPK/CNS models) are no longer
  zeroed. admixr2 used to hand-fill every model parameter it did not set
  with `0`, clobbering such a constant’s default and producing an
  `NA`/non-finite objective (e.g. a `qout / vb` divide-by-zero). It now
  supplies only the parameters it varies and lets
  [`rxSolve()`](https://nlmixr2.github.io/rxode2/reference/rxSolve.html)
  fill the rest from the model’s own defaults, so constants and
  covariate defaults keep their value.

- `adgh` now computes gradients for non-mu-referenced (unpaired)
  structural thetas. The unpaired-parameter set was derived from the
  eta-indexed `struct_eta_idx`, so it was always empty and those thetas
  silently received a zero gradient; it now uses the struct-indexed
  `struct_has_eta`.

## admixr2 0.2.0

CRAN release: 2026-07-02

### New features

- New estimator `est = "adgh"`: deterministic Gauss-Hermite quadrature
  over the random-effects prior, configured via
  [`adghControl()`](https://leidenpharmacology.github.io/admixr2/reference/adghControl.md).
  The objective is noise-free (no Monte Carlo draws), the analytical
  gradient is exact, and it is unbiased at any IIV magnitude. For models
  with up to ~4 random effects it is the fastest exact estimator
  ([\#65](https://github.com/LeidenPharmacology/admixr2/issues/65)).
- [`datagen()`](https://leidenpharmacology.github.io/admixr2/reference/datagen.md)
  gains FO-approximated population moments (`method = "fo"`, matching
  `est = "adfo"`) for design evaluation and optimal-design work
  ([\#56](https://github.com/LeidenPharmacology/admixr2/issues/56)).
- `adirmcControl(kappa_method = "linearized_gh")`: GH-averaged kappa
  baseline for the IRMC inner loop.
- [`admClearCache()`](https://leidenpharmacology.github.io/admixr2/reference/admClearCache.md)
  prunes the session-level compiled-model cache
  ([\#10](https://github.com/LeidenPharmacology/admixr2/issues/10)).
- Control objects now accept any `nloptr` algorithm; the default is
  chosen from the gradient mode, and `grad`/`algorithm` are reconciled
  automatically
  ([\#70](https://github.com/LeidenPharmacology/admixr2/issues/70)).

### Bug fixes

- Fix an infinite recursion (“evaluation nested too deeply” / “node
  stack overflow”) that aborted the first fit of an R session when a
  covariance matrix was requested (`covMethod = "r"`). Accessing
  `ui$simulationModel` left a self-referential compiled-model object in
  `ui$meta`, which nlmixr2’s ui-cloning during fit assembly could not
  traverse. admixr2 now clears that transient artifact in
  `.admLoadModel()`, keeping the ui in the canonical state nlmixr2
  expects. Affected all four estimators (`adfo`/`admc`/`adgh`/`adirmc`)
  ([\#81](https://github.com/LeidenPharmacology/admixr2/issues/81)).
- Use the ML denominator (`1/n_sim`) consistently in the MC gradient
  kernels, matching the NLL
  ([\#48](https://github.com/LeidenPharmacology/admixr2/issues/48)).
- Fix parallel multi-restart dispatch for fork/PSOCK, and fix `adirmc`
  multi-restart
  ([\#45](https://github.com/LeidenPharmacology/admixr2/issues/45)).
- Guard non-positive predicted variance in the diagonal-NLL paths
  ([\#57](https://github.com/LeidenPharmacology/admixr2/issues/57)).
- Correct the FO diagonal omega gradient scaling, plus assorted plot,
  output-variable detection, caching, and worker-serialization fixes.

### Documentation

- Add Gauss-Hermite sections across the vignettes and fix the pkgdown
  reference index so the documentation site builds
  ([\#79](https://github.com/LeidenPharmacology/admixr2/issues/79)).

### Dependencies

- Declare minimum versions for the imported `rxode2 (>= 5.1.2)` and
  `nlmixr2est (>= 6.0.1)`, and for the suggested `nlmixr2 (>= 5.0.0)`
  (used in examples and tests).

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
- `traceplot()` support: admixr2 fits populate the standard
  `parHistData` slot, so the nlmixr2 `traceplot()` generic works
  natively (best restart, natural scale, no burn-in marker).
- Integrates with the nlmixr2/rxode2 ecosystem.
