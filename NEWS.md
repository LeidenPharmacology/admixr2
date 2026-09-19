# admixr2 0.4.1

## New features

* **A conditional source is cut on the span of the two models, not on the
  product grid.** The objective reaches the covariates only through the source
  direction and the analysis direction, so a source conditional on three
  covariates needs nodes on that plane rather than the full `strata_nodes^3`
  cube: 729 studies become 81, exactly rather than approximately. Truncated
  margins keep the reduction through Caratheodory recombination, which draws
  its atoms from a fine cloud and reproduces that cloud's moments in every
  projected direction. The reduction is admitted only where the analysis
  model's covariate loading is parameter-invariant -- fixed allometric
  exponents, in practice -- because the nodes are cut once at materialisation
  and an estimated exponent would carry the loading off the plane they sit on.
  Anything else keeps the product grid.

* **One `rxSolve` per group of studies rather than one per study.** A call
  costs 0.0251 s to enter and 1.6e-06 s per subject, so at the sizes a node
  expansion produces it is almost all overhead; the objective, the gradient and
  the post-fit diagnostics now stack their parameter frames and solve together,
  across different doses and schedules alike. A four-study fit went from 126.4 s
  to 38.9 s, and the sampled profile from 83.97 s to 41.78 s, with the objective
  unchanged to every digit.

* **The mean and covariance panels are per SOURCE, not per stratum.** A conditional
  source is collapsed by the mixture law, so `_s1`/`_s2` never reaches a figure.

* **The panels name what `V` contains** -- `BSV + covariate spread + sigma`,
  read off the fit -- and the predicted ribbon shows the pre-sigma part inside.

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

* **The effect panel compares the estimated effect against its sources**: a
  dotted line, against each source's own published model.

* **A conditional source draws its OWN regression over the range it covers**;
  a marginal one draws a whisker, because what it reported is a distribution.

* **One mark per source per facet**, not one per stratum, so a source conditional on
  sex no longer draws twice on every other covariate's panel.

* **Both covariate panels key colour on the source**, so a source keeps one
  colour across the figure; its strata are joined in grey instead.

* **A residual facet is dropped when the sources' contrast is sampling noise**,
  measured against a typical within-study 10th-90th.

* **Point area is the study's sample size on both covariate panels**, with a
  legend. The residual panel encoded it already and said so nowhere.

* **Whether a covariate is conditional or marginal is derived, not declared**:
  conditional when the source's own model ESTIMATED its coefficient.

* **`admMoments(fit)` gives the observed and predicted first two moments per
  source**, with the structural variance share and the standardised residual.

## Changes that can move an existing fit

Several changes in this release alter results for scripts that do not name a new
argument. None is a bug fix, so all are listed here rather than below.

* **`covDist(joint = )` and a `population` carrying its own sampler are refused**,
  pending the vine-copula work; `cor` is the supported route to dependence.

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

* **`admStudy(stratify = )` is removed**, conditioning being derived: a script
  naming fewer covariates than its model uses moves. `covStrata()` still takes it.

* **`range` truncates a MARGINAL covariate's declared distribution too**, not
  only a conditional one's: the enrolled span holds however the covariate is used.

* **Nodes need an ESTIMATED coefficient, not just a covariate the model
  reads**, so a fixed allometric exponent leaves weight marginal.

* **`fit$env$strataNodes` is now one entry per covariate the model reads**,
  rather than one number per fit, and `anova()` compares where the two overlap.

## Bug fixes

* **`adgh` on a no-IIV (`n_eta = 0`) model failed under covariate
  marginalisation**, with a dimnames-length error from a phantom `"eta."` name.

* **A DISCRETE covariate latently correlated with any other margin is refused**,
  rather than integrated as if it were independent; a level is not a point.

* **The estimated-effect curve and the source marks were different quantities**
  for a staged model, so every source drew a constant factor off the fitted line.

* **A `by =` source got no mark and no regression line on `covariate_effect`**,
  its own model being reachable under no name the panels use.

* **`print()` told a `by =` source its estimated coefficient was asserted**, and
  that it conditions on nothing, contradicting its own `reported by` line.

* **The no-enrolled-range warning fired for DISCRETE conditional covariates**,
  where truncation is a no-op; it asked for a span for a two-level factor.

* **A `range` given as `c(lo =, hi =)` was applied by the fit and ignored by the
  plot**, and one keyed by covariate names without being a list is now refused.

* **An unnamed `range` is resolved against the SOURCE's conditional covariates**,
  not the analysis-narrowed set: one `studies` object now serves a nested pair.

* **The source marks and the source regression line took the FIRST assignment
  to a parameter**, where the solve uses the last; a staged model drew a zig-zag.

* **An enrolled `range` is refused against a `cov_dist` with its own `joint`
  sampler**, rather than accepted and silently ignored.

* **`anova()` skipped the resolution check when only one fit's stamp was
  named**, which is the pre-rename fit the unnamed comparison exists to catch.

* **`fit$env$strataNodes` tells a marginalised covariate from one pinned at a
  single node**, and reports the whole set of node counts rather than its maximum.

* **A `cov` entry longer than one took the `covariate_resid` panel out**, through
  a label the panel computed and never read.

* **A free-scaled covariate facet could be ticked at another covariate's
  levels**; the breaks are read per covariate now.

* **The `diagnostic-plots` article called an undefined function** and could not
  be built.

* **`datagen()` cut every study that declared a covariate distribution into
  nodes**, having no `model` of its own to derive from. It needs `stratify` now.

* **`stratify = FALSE` stopped reaching the spec**, so a study that refused
  conditioning had one derived for it -- the opt-out became its opposite.

* **`strata_nodes` was recorded only when something was cut into nodes**, so it
  could vary between studies of one fit -- which `anova()` refuses to compare.

* **`anova()` refused the nested pair a covariate test is made of**, the null
  having dropped the term and so having no nodes; measured, they agree to 5e-05.

* **`range` was silently dropped by a source conditional on nothing**, the
  transcribed `mean +/- SD` it exists for.

* **The rebuilt stratum sampler was discarded one line later**, leaving
  correlated conditional margins to be drawn independently of each other.

* **A continuous covariate conditional at `strata_nodes <= 8` was drawn on a
  discrete axis**, a dot per quadrature node and the ticks on that grid.

* **An unnamed `range` crashed the `covariate_effect` panel out of existence**,
  the error being caught and the panel reported as absent.

* **A stratum's source is recorded rather than recovered by regex**, so two
  studies a user named `a_s1` and `a_s2` are no longer merged into one.

* **`admMoments()` returned `NULL` where it documents an empty data frame.**

* **A source conditional on a covariate the analysis model does not read could
  not be reduced**: the sampler's inputs are kept, so it rebuilds on the subset.

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

* **No effect panel for a covariate the model does not estimate** -- a fixed
  allometric exponent is not a finding to check agreement on.

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

* **Coincident marks from different sources are no longer merged**, which
  invented labels like `3 studies` for studies that share only a position.

* **An eighth source gets its own colour** instead of the first one's, via an
  HCL ramp past the seven Okabe-Ito hues.

* **Black is reserved for the fit** in the covariate panels, so no source is
  drawn in the colour of the thing it is compared against.

## Internal changes

* **The covariate panels are built by `.admCovEffectPanel()` and
  `.admCovResidPanel()`**, not inline in `plot.admFit()`, which loses 170 lines.

* **Visual regression tests for the covariate panels**, on synthetic studies
  and no fit. Requires the new `vdiffr` suggested dependency.

* **`vignette("diagnostic-plots")` renders the covariate panels** from its own
  three-source fit, instead of describing them in prose.

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
