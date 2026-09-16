# admixr2 0.4.1

## New features

* **A covariate the model does not read is dropped, not refused**, so a nested
  pair shares one `studies` object and needs no `fix(0)`.

* **`covMethod = "r,s"`: standard errors scored against the model's own
  sampling law** -- an ADF sandwich `H^-1 J H^-1`, on all four estimators.

* **`v_denom` declares which denominator a study's `V` uses**, per study,
  rather than a hand-applied `(n-1)/n` correction. Default `"ml"`.

* **`sigdig` now controls the fit, not just the output tables**, and is
  opt-in: the default `NULL` leaves rxode2's own tolerances alone.

* **adfo differentiates its structural thetas analytically**, so `grad =
  "analytical"` (LBFGS) is now its default.

* **`linCmt()` models are supported at second order**, by promotion to the
  explicit ODE form.

* **Finite-difference steps are measured per parameter** (Shi 2021), replacing
  the fixed `pmax(abs(p), 0.1) * h` scale.

* **`anova()` on nested fits**: the likelihood-ratio test on the objective
  difference against a chi-squared reference.

* **A study can contribute as a published MODEL**, not only as digitised
  aggregate data. No standard error is reported for a fit containing one.

* **Covariate marginalisation over a declared distribution** (`covDist()`),
  for `admc` and `adgh`, instead of solving at the covariate mean.

* **A sparse-grid route for several covariates**: `cov_integration = "sparse"`
  uses a Smolyak grid -- 49 points against the product grid's 81.

* **A paper-shaped study API**: `admStudy()`/`admStudies()` describe a source
  as a publication does, and `admPopulation()` reads a baseline table.

* **The covariate integral is collapsed onto the directions it actually has**
  -- same answers, fewer design points.

* **`plot(fit, which = "covariate")`**: the fitted covariate effect, and each
  study's mean residual against it -- a slope is a bad covariate form.

## Changes that can move an existing fit

Several changes in this release alter results for scripts that do not name a new
argument. None is a bug fix, so all are listed here rather than below.

* **The sandwich's `G` is evaluated at `tau`, not at the observed summary**,
  so every `covMethod = "r,s"` standard error moves slightly.

* **`covMethod` now defaults to `"r,s"`**, so reported standard errors change
  for every script that does not name it.

* **Transform-both-sides endpoints are composed exactly**, by quadrature, and
  their estimates move.

* **The per-parameter Shi (2021) step is not optional**, so any fit relying on
  a finite difference moves slightly.

* **`gill` is removed from all four controls**, superseded by the Shi (2021)
  step search.

* **Forward finite differences are removed**: `grad = "fd"` is now a central
  difference, and `grad = "cfd"` is gone.

* **`adirmcControl(grad = "fd")` differences with `grad_h`**, not a hard-coded
  `1e-6`.

* **`adfoControl()`'s new `grad = "analytical"` default brings the
  `grad_bounds` box with it.**

## Bug fixes

* **A discrete covariate latently correlated with another margin is refused**,
  rather than integrated as if it were independent.

* **A joint collapse probes Omega as well as the structural parameters**
  before it freezes the design's rank.

* **Derivative-free fits no longer inherit nloptr's loose `xtol_rel = 1e-4`.**

* **Parallel restarts could fail to read the compiled-model cache** whenever a
  second R session was using admixr2 at the same time.

* **adfo could report `NA` for every standard error** on a fit that converged
  normally.

* **A joint study normalised before the model was known kept `NULL` block
  outputs.**

* **`.admCacheWrite()` could delete another session's valid cache entry.**

* **A sensitivity model that failed to build was reported as quietly as one
  refused by design.**

* **A fit that stops on the gradient box constraint now says so audibly.**

* **The order-2 `linCmt()` promotion did not run for a linCmt assigned to a
  variable**, so adfo kept finite-differencing its structural thetas there.

* **`adghControl()` accepted an invalid nloptr algorithm**, and would hand a
  derivative-free one a gradient.

* **`adirmcControl()` validated neither `ci` nor `returnAdmr`.**

* **A cache write that fails no longer discards the model it just compiled**,
  or kills the fit.

* **The session-ownership guard rejected nlmixr2est's own sensitivity model**
  unconditionally.

* **The order-2 `linCmt()` promotion could write a theta's value into the
  wrong `THETA[k]` slot.**

* **The gradient-box warning judged the fit against the wrong point**, and
  stayed silent for the parameters most likely to need it.

* **An explicit `adfoControl(grad = "analytical")` that cannot build a
  sensitivity model warns again.**

* **Normalising a study twice no longer leaves its endpoint unset.**

* **Dev-mode parallel restarts could not see any function this release
  introduced.**

* **Generated models are built under role-tagged names**, in their own
  directory, and a cached one is checked before it is trusted.

* **Normalising a study twice turned it into a joint (same-subject) study.**

* **Non-finite parameters no longer reach the ODE solver.**

* **A cache-key collision solved fits at another model's fixed value.**

* **A stale sensitivity cache entry could survive a change to what it
  caches.**

* **`linCmt()` second-order promotion built its direction set from the
  pre-promotion model.**

* **A struct theta missing from the cached direction map crashed the fit.**

* **A transformed endpoint no longer pays for second-order compartments it
  cannot use.**

* **A parallel worker no longer walks on from a model it could not load.**

* **`.admNLL()` gained the non-finite screen the other estimators got.**

* **Compiled models are held in a session cache**, instead of being reloaded
  from the disk cache on every call.

* **The diagnostic panels solved a covariate study at the covariate mean**, so
  its predicted moments lost the spread the observed ones carry.

* **A covariate the model no longer reads keeps its declared distribution**,
  so the residual panel can still plot a deleted term against it.

* **Stratified studies are titled by the covariate value they condition at**,
  not a bare `_s1`/`_s2` index.

* **A covariate facet is kept when the model is flat at one conditioned level**
  and varies at another, instead of being dropped on the first level alone.

* **A marginalised level mix sits at its mean, not its median**, so an even
  binary split no longer greys a studied level or loses its residual facet.

* **The `lm` trend is guarded per covariate**, not across the figure, so a
  facet with two studies gets no line through its two points.

* **A source keeps one colour across both covariate panels**, and the effect
  panel's curves split only on covariates some study conditions at a point.

* **Merging a conditioned mark into a marginal one no longer invents a spread**
  for it; the merged mark is drawn as the conditioned one it contains.

* **A shaded discrete level no longer costs the axis its level-only ticks**,
  which had a three-level factor ticked at 0.5 and 1.5.

* **The covariate panels' legends keep a fixed order.** ggplot2 sorts
  equal-priority guides by a hash, which is not stable across sessions.

## Internal changes

* **The covariate panels are built by `.admCovEffectPanel()` and
  `.admCovResidPanel()`**, not inline in `plot.admFit()`, which loses 170 lines.

* **Visual regression tests for the covariate panels**, on synthetic studies
  and no fit. Requires the new `vdiffr` suggested dependency.

* **`print.admFit()` reaches nlmixr2est's printer through `getS3method()`**,
  not an `asNamespace()` lookup.

* **The sensitivity-model builder takes an `order` argument**: `1L` is the
  existing direction set, `2L` adds the cross block adfo needs.

* **CI: `R-CMD-check` gained a `workflow_dispatch` trigger** and a dependency
  cache-version bump.

# admixr2 0.4.0

## New features

* **Student-t residual error (`cp ~ add(a) + t(nu)`) is supported**, as
  nlmixr2's scale family: the residual is `scale * T_nu`.

* **The transform-both-sides transforms call rxode2's own kernel**, instead of
  an inline reimplementation.

* **`.admBackTransform()` uses `rxode2::probitInv()`** instead of an inline
  `low + (high - low) * pnorm(p)`. Numerically identical.

* **New `resid_nodes` control argument** on all four estimators: the
  Gauss-Hermite node count for a transform-both-sides residual integral.

* **New vignette: "Choosing a residual error model"**
  (`vignette("error-models", package = "admixr2")`).

## Bug fixes

* **Dropped the `qs2` dependency**: the compiled-model and sensitivity disk
  caches use `saveRDS()`/`readRDS()`.

* **IRMC importance-sampling shift was wrong for every non-`exp` mu-referenced
  theta.**

* **A `fix()`ed prediction-dependent residual lost its gradient.**

* **`binom(20L, p)` was refused as a non-constant size.**

* **A non-positive `nbinomMu` size now gives a clear domain error.**

* **`beta` precision denominator is guarded against a zero draw.**

* **Standard errors: sigma SEs were uninitialised memory, and omega was
  excluded.**

* **Omega and sigma standard errors are reported**, on the scale the estimates
  are printed on.

* **A printed standard error now belongs to the parameter it is printed
  beside.**

* **Count endpoints could not be fitted with the default gradient.**

* **The covariance Hessian used the starting lambda for a transformed
  endpoint.**

* **`beta()` endpoints were only ever right on the plain NLL path.**

* **The `ar()` and `ordinal` guards judged every study, not the affected
  one.**

* **Ordinal categories were grouped by exact floating-point time equality.**

* **The moment expansion and its derivative capped the same pole
  differently.**

* **A parallel worker could invert a transform with another model's lambda.**

* **`plot()` back-transformed three residual roles on the wrong scale.**

* **A `binom` size written as a model constant was refused as non-constant.**

* **A count or beta endpoint alongside another endpoint is now refused.**

* **`datagen()` refuses an ordinal endpoint.**

* **`resid_nodes` no longer changes what a positional call means.**

* **The ordinal same-time grouping is defined once.**

* **Endpoints transformed differently from one another refused the sensitivity
  model.**

* **A joint (same-subject) study had no aggregate diagnostics.**

* **Documented: an `adfo` standard error describes scatter, not accuracy.**

* **A failed covariance is no longer silent.**

* **A study `ev` containing observation records now warns.**

* **Residual parameters fixed with `fix()` were silently dropped.**

* **A `prop()`/`pow()` term on a transform-both-sides endpoint contributed
  nothing.**

* **The post-fit covariance was a Hessian of the wrong objective for several
  error models.**

* **`adfo` dropped `ar()` from its objective while keeping it in the
  gradient.**

* **An out-of-support transform aborted the whole fit.**

* **The sensitivity-model cache could serve a stale transform spec.**

* **`0^negative` in the moment expansion.**

* **`ordinal` endpoints are supported.**

* **`dv()` is now refused.**

* **`ar()` combined with `prop()`/`pow()`/combined is now refused.**

* **Known upstream issue: simulating an `ar()` fit will not reproduce its
  covariance.**

* **Prediction-dependent residual error is composed correctly** (`prop()`,
  `pow()`, `lnorm()`, combined).

* **`lnorm()` analytic gradients were computed against the wrong quantity.**

* **`delay()` (DDE) models get an accurate sensitivity solve.**

## Internal changes

* **The post-fit covariance's reported-scale rotation and its non-PD omega
  fallback are now shared helpers.**

* **The residual variance's dependence on `(mu, var_f)` is computed once per
  study/unit**, not three times.

* **The residual V-composition tail is one helper, `.admApplyResidTail()`.**

# admixr2 0.3.0

## New features

* **Analytical gradients for non-mu-referenced ("unpaired") structural
  thetas.**

* **Residual error models `pow()`, `addPow()` and `combined1()` are
  supported**, with analytical gradients.

* **Multi-compartment fitting (multiple observed outputs).**

* **Parallel restarts run on `mirai` daemons.**

* **New `nDisplayProgress` control argument.**

* **The aggregate-data estimators carry `type` and `description` attributes**,
  classifying them as Model Based Meta Analysis.

## Bug fixes

* **`pow()` models no longer fit the wrong residual model, silently.**

* **`combined1()` is honoured.**

* **An unrepresentable residual model is refused rather than approximated.**

* **`propT`/`propF`, `norm`/`dnorm` and `dlnorm`/`logn`/`dlogn` no longer emit
  spurious approximation warnings.**

* **Lognormal residual error is applied to the plotted predicted mean.**

* **The solver progress bar no longer appears during covariance/gradient
  batches.**

* **Hard-coded numeric constants in a `model({})` block are no longer
  zeroed.**

* **`adgh` computes gradients for non-mu-referenced (unpaired) structural
  thetas.**

* **Parallel restarts under `devtools::load_all()` warn once about the
  installed package.**

## Internal changes

* **`adgh` gradient-mode fits are about twice as fast**: the objective and the
  gradient share one solve.

* **Model loading and per-fit memory follow nlmixr2est's own conventions.**

* **`admClearCache()` is removed; use `rxode2::rxClean()`.**

* **`print()` on a fit no longer writes into rmarkdown's namespace.**

# admixr2 0.2.0

## New features

* New estimator `est = "adgh"`: deterministic Gauss-Hermite quadrature over
  the random-effects prior, via `adghControl()`. Noise-free and exact (#65).
* `datagen()` gains FO-approximated population moments (`method = "fo"`) for
  design evaluation and optimal-design work (#56).
* `adirmcControl(kappa_method = "linearized_gh")`: GH-averaged kappa baseline
  for the IRMC inner loop.
* `admClearCache()` prunes the session-level compiled-model cache (#10).
* Control objects accept any `nloptr` algorithm; the default is chosen from
  the gradient mode, and `grad`/`algorithm` are reconciled (#70).

## Bug fixes

* Fix an infinite recursion that aborted the first fit of an R session when a
  covariance matrix was requested. All four estimators (#81).
* Use the ML denominator (`1/n_sim`) consistently in the MC gradient kernels,
  matching the NLL (#48).
* Fix parallel multi-restart dispatch for fork/PSOCK, and fix `adirmc`
  multi-restart (#45).
* Guard non-positive predicted variance in the diagonal-NLL paths (#57).
* Correct the FO diagonal omega gradient scaling, plus assorted plot,
  output-variable detection, caching and worker-serialization fixes.

## Documentation

* Add Gauss-Hermite sections across the vignettes and fix the pkgdown
  reference index so the documentation site builds (#79).

## Dependencies

* Declare minimum versions for `rxode2 (>= 5.1.2)`, `nlmixr2est (>= 6.0.1)`
  and the suggested `nlmixr2 (>= 5.0.0)`.

# admixr2 0.1.0

* Initial release.
* Monte Carlo estimator (`est = "admc"`) via `admControl()`.
* Iterative Reweighting Monte Carlo estimator (`est = "adirmc"`) via
  `adirmcControl()`.
* Analytical CRN gradient with sensitivity equations (`grad = "sens"`).
* Multi-restart parallelism via `furrr`/`future`.
* Diagnostic plots: observed vs predicted mean/covariance, NLL trace,
  parameter trace.
* `traceplot()` support: admixr2 fits populate the standard `parHistData`
  slot, so the nlmixr2 generic works natively.
* Integrates with the nlmixr2/rxode2 ecosystem.
