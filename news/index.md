# Changelog

## admixr2 0.3.0

### New features

- **Analytical gradients for non-mu-referenced (“unpaired”) structural
  thetas.** A structural theta with no mu-referencing eta (`tka` with no
  `eta.ka`, or the `exp(tcl) * exp(eta.cl)` writing style rxode2 does
  not mu-reference) used to cost an extra finite-difference `rxSolve`
  per gradient call. admixr2 now emits its own first-order sensitivity
  model over an explicit direction set (one direction per random effect
  plus one per unpaired theta), compiled with `eventSens = "jump"` so
  dosing-modifier (`f`/`lag`/`rate`/`dur`) sensitivities are no longer
  silently zero. This mirrors the scheme nlmixr2est’s fast-focei uses
  (`.foceiAnalyticDirections`) but first-order only, and is
  cross-validated against nlmixr2est’s inner model to ~1e-13 across ODE,
  linCmt, dosing modifiers, initial conditions, covariates, if/else and
  multi-endpoint models. Consumed by `admc`, `adgh` (including joint
  multi-output studies); `adfo` keeps finite differences (its
  `V_pred = J Omega J' + Sigma` needs a second derivative). Measured
  2.5-3.8x faster and ~100x more accurate than the previous
  finite-difference path on a 2-compartment model. This adds `symengine`
  (already a hard dependency of `nlmixr2est`, so always installed
  alongside admixr2) to `Imports`, used to emit the linCmt direction
  derivatives. The feature degrades gracefully on rxode2 without
  `eventSens = "jump"` support (it falls back to the finite-difference
  path), so no minimum-version bump is required.

- **Residual error models: `pow()`, `addPow()` and `combined1()` are now
  supported, with analytical gradients**
  ([\#84](https://github.com/LeidenPharmacology/admixr2/issues/84)).
  admixr2 previously supported only `add`, `prop` and `lnorm`. The
  residual error model is now read from `ui$predDf`
  (`errType`/`errTypeF`/`transform`/`addProp`) rather than from
  `iniDf$err` alone, and every estimator evaluates it through one shared
  specification:

  | form                                   | variance                 |
  |----------------------------------------|--------------------------|
  | `combined2` (default for `add + prop`) | `a^2 + b^2 * f^(2c)`     |
  | `combined1`                            | `(a + b * f^c)^2`        |
  | `lnorm`                                | moment-matched lognormal |

  with `c = 1` recovering `prop` and `b = 0` recovering `add`.
  Analytical `d(var)/d(sigma)`, `d(mu)/d(sigma)` and `d(var)/d(f)` are
  supplied for all of them, so residual parameters keep an exact
  gradient under `grad = "sens"`/`"analytical"`.

  Existing `add`/`prop`/`lnorm` fits are unaffected: the aggregate
  `-2LL` is bit-for-bit identical, and their gradients change only by
  floating-point reassociation (~1 ulp).

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

- The aggregate-data estimators (`adfo`, `adgh`, `adirmc`, `admc`) now
  carry `type` and `description` attributes classifying them as “Model
  Based Meta Analysis” methods, so they appear in the category-grouped
  estimation-method list nlmixr2est prints for an unsupported `est=` (or
  a bare `nlmixr2()` call)
  ([\#107](https://github.com/LeidenPharmacology/admixr2/issues/107)).

### Bug fixes

- **`pow()` models no longer fit the wrong residual model, silently.**
  `pow(b, c)` produces two `iniDf` rows – the coefficient
  (`err = "pow"`) and the *exponent* (`err = "pow2"`). admixr2
  recognised neither, warned once, and then treated **both as additive
  variances**: the exponent was stored as `2*log(c)` and optimized as a
  variance contributing `exp(2*log(c))` to `diag(V)`. A `pow` model
  therefore ran to completion and reported plausible estimates for a
  model it was not fitting. Residual parameters now carry a role, and a
  `pow` exponent is estimated on its own (unconstrained, identity)
  scale.

- **`combined1()` is honoured.** `predDf$addProp` selects SD-additive
  (`combined1`) versus variance-additive (`combined2`) residual error.
  admixr2 ignored it and always computed `combined2`, dropping the
  `2*a*b*f` cross term. (`combined2` is nlmixr2’s default, so only
  models that explicitly asked for `combined1()` were affected.)

- **An unrepresentable residual model is now refused rather than
  approximated.** Error types admixr2 cannot express as a Gaussian
  aggregate MVN (`logitNorm`, `probitNorm`, Box-Cox/Yeo-Johnson
  transforms, `t`/`cauchy`, `propF`/`powF`) previously emitted a
  one-time warning and were then **treated as additive**, so the fit
  proceeded with the wrong residual model. They now
  [`stop()`](https://rdrr.io/r/base/stop.html). This is a behaviour
  change: a model that “worked” before may now error.

- **`propT`/`propF`, `norm`/`dnorm` and `dlnorm`/`logn`/`dlogn` no
  longer emit spurious “modelled as …” approximation warnings.** These
  are aliases, not approximations: `norm` *is* `add`, `logn` *is*
  `lnorm`, and on an untransformed model `propT` (which scales by the
  transformed prediction) *is* exactly `prop`, because there the
  transformed and untransformed predictions are the same quantity. The
  warnings claimed an inaccuracy that did not exist.

- **Lognormal residual error is now applied to the plotted predicted
  mean.**
  [`plot.admFit()`](https://leidenpharmacology.github.io/admixr2/reference/plot.admFit.md)’s
  aggregate-data helper added the lnorm variance to the predicted
  covariance but never applied the `exp(s/2)` mean scaling to the
  predicted `E`, so lnorm fits plotted a mean the NLL does not use.

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

- **Parallel restarts under `devtools::load_all()` warn once about the
  installed package.** In dev mode the admixr2 namespace is locked, so
  worker daemons run the *installed* package rather than the loaded
  source; if it is stale the parallel objective silently diverges from
  the sequential one. `.admRunRestarts` now emits a one-time warning in
  this case telling you to `devtools::install()`. It never fires in
  production (installed package == source).

### Internal changes

- **Model loading and per-fit memory now follow nlmixr2est’s own
  conventions.** admixr2 previously pinned each fit’s `foceiModel`
  companion objects in a package-level environment (a Windows
  GC-finalizer heap-corruption guard) and reclaimed rxode2’s global
  model registry with a bespoke snapshot/teardown after every fit. Both
  are gone: the companion objects are no longer pinned (the guard proved
  unnecessary – verified by running the `covMethod = "r"` fit path
  repeatedly under aggressive GC with no crash), and each estimator now
  frees memory the way nlmixr2est does, with
  [`gc(); rxode2::rxUnloadAll()`](https://rdrr.io/r/base/gc.html). The
  disk model cache continues to use `qs2` + `digest`, exactly like
  rxode2/nlmixr2est; the in-memory pin cache was removed (same-model
  reloads come from the `qs2` files). Net: ~290 fewer lines, no
  admixr2-specific memory machinery, and fit results are unchanged.

- **`admClearCache()` is removed; use
  [`rxode2::rxClean()`](https://nlmixr2.github.io/rxode2/reference/rxClean.html).**
  admixr2’s `qs2` caches live in
  [`rxode2::rxTempDir()`](https://nlmixr2.github.io/rxode2/reference/rxTempDir.html)
  alongside rxode2’s and nlmixr2est’s, so
  [`rxode2::rxClean()`](https://nlmixr2.github.io/rxode2/reference/rxClean.html)
  – rxode2’s standard cache wipe (unload all models + clear the temp
  dir), which nlmixr2est itself calls to reset – already clears
  admixr2’s cache too. The package-specific `admClearCache()` is
  therefore redundant.

- **[`print()`](https://rdrr.io/r/base/print.html) on a fit no longer
  writes into rmarkdown’s namespace.** `print.admFit` temporarily
  overwrote `rmarkdown:::print.paged_df` via
  [`assignInNamespace()`](https://rdrr.io/r/utils/getFromNamespace.html)
  (restoring it `on.exit`) to steer nlmixr2est away from its paged-table
  branch. That branch is in fact unreachable: nlmixr2est decides between
  paged and console output by *probing behaviour* – it prints a
  `paged_df`-classed frame into
  [`capture.output()`](https://rdrr.io/r/utils/capture.output.html) and
  infers “a paged renderer consumed my output” from zero captured lines
  – but `rmarkdown:::print.paged_df` returns its `knit_asis` object
  visibly and no `print.knit_asis` method exists, so the probe always
  collects output, always returns `FALSE`, and the console branch is
  always taken. The stub therefore changed nothing except skipping the
  discarded probe render (~20 ms per `print(fit)`), at the cost of
  mutating a foreign namespace – fragile, unsafe under concurrent
  rendering, and a CRAN-policy grey area. Printed output is unchanged,
  byte for byte.
  ([\#58](https://github.com/LeidenPharmacology/admixr2/issues/58))

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
- `admClearCache()` prunes the session-level compiled-model cache
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
