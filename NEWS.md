# admixr2 0.4.1

## New features

* **`sigdig` now controls the fit, not just the output tables -- and it is
  opt-in.** The `sigdig` and `rxControl` arguments were documented as solver
  controls, but the object they built only ever reached nlmixr2's *post-fit*
  table solves: every optimizer solve ran at rxode2's own default tolerances, so
  setting `sigdig` changed nothing about the fit. It is now passed to
  `rxode2::rxSolve()`'s own `sigdig` argument at every solve the estimators
  issue.

  **The default is `sigdig = NULL`, so default fit results are unchanged.**
  `NULL` means "leave rxode2's own tolerances alone", which is exactly the
  numerics every admixr2 fit had before. It is the default because a looser
  solve is not free: the estimators finite-difference these solves with steps of
  the same order -- `grad_h` 1e-4, `cov_h` 1e-3, `cov_h_outer` ~2.5e-3 -- while
  rxode2 5.1.5 maps `sigdig = 4` to `rtol = 1e-4`. Differencing a solution whose
  own relative noise is 1e-4 with a 1e-4 step returns noise, and it surfaces as a
  moved objective and an indefinite covariance Hessian (every `SE` reported `NA`)
  rather than as an error. Turning that on by default would have changed the
  numerics of every existing script silently, for a knob that looked like table
  formatting before this release.

  Set explicitly, it is the lever for trading solver accuracy against speed.
  Measured on a 1-cmt oral ODE model with two studies, at a fixed iteration
  count, `sigdig = 4` makes **adfo 4.8x faster**, `admc` 1.3x, and `adgh`
  unchanged (its batched quadrature solve is not integration-bound); the
  objective moves by 5e-09 relative and the `covMethod = "r"` standard errors are
  unchanged to four significant figures. It is most worthwhile where the gradient
  is fully analytic and nothing differences the solve, which after this release
  is `adfoControl(grad = "analytical")`. Elsewhere, compare the objective and the
  standard errors against `NULL` before relying on it.

  Passing the digits rather than re-deriving tolerances keeps the mapping
  rxode2's business, which matters because rxode2 has changed it between
  releases: `sigdig = 4` is `atol = rtol = 5e-07` on rxode2 5.1.4 but
  `rtol = 1e-04` on 5.1.5. `sigdig = NULL` is the one setting whose meaning does
  not move under an upgrade -- and no single `sigdig` value reproduces rxode2's
  defaults anyway, since they are asymmetric (`atol` 1e-8 against `rtol` 1e-6)
  while the `sigdig` map is one-dimensional. Table formatting is unaffected
  either way: `sigdigTable` falls back to 4 when `sigdig` is `NULL`.

  `plot()` now solves at the tolerance the fit used, so the diagnostic panels
  describe the same integration the objective was minimised on. `datagen()`
  deliberately does not -- it generates the reference, and integrates at
  rxode2's own tolerances regardless.

* **adfo differentiates its structural thetas analytically, and
  `grad = "analytical"` (LBFGS) is now the default.** adfo was the last estimator
  with a finite-difference component: its `V_pred = J Omega J' + resid` depends on
  a structural theta through `J`, so the gradient needs `dJ/dtheta`, a *second*
  derivative that the first-order sensitivity model does not carry. It therefore
  finite-differenced the whole objective -- differencing a log-determinant and a
  quadratic form, the noisiest construction in the package.

  `.admBuildThetaSens()` now emits a second-order **cross block**
  `d2f/(d eta d dir)` on request, and `.adfoGrad()` contracts it against the same
  `dNLL/dV` matrix the omega path already uses, so the two cannot drift apart.
  Against a central difference of the objective the structural gradient is
  accurate to 2e-07..2e-06, where the finite-difference pass it replaces reached
  8e-04..1e-02 (worst case `tv`, 1.4%).

  Because the gradient is now exact, LBFGS on it beats the derivative-free
  BOBYQA that `grad = "none"` used, so the default changed. `grad = "none"`
  remains available, and any model whose second-order model cannot be built
  falls back to the previous finite-difference gradient automatically.

  **That default flip changes three more things than the gradient**, because
  `grad != "none"` is the switch for all four. Spelled out, since only the first
  is obvious:

  1. The gradient itself, as above.
  2. **A box constraint.** A gradient fit is confined to `p0 +/- grad_bounds`
     (default 5) on the optimizer scale -- a factor of ~148 on the log scale.
     adfo fits were unbounded before. nloptr reports normal convergence at a box
     corner, so admixr2 now *warns* when an estimate finishes on that bound and
     the model itself declared none; `admc`/`adgh` have always run this way and
     gained the same warning.
  3. **The covariance method.** `covMethod = "r"` builds its Hessian by
     forward-differencing the gradient rather than the objective when a gradient
     is available. That is now gated on the struct-theta gradient being genuinely
     analytic, not merely on `grad != "none"`: with an order-1 fallback or a
     transformed endpoint the gradient is itself a finite difference, and
     differencing it again produced a singular Hessian and "standard errors are
     unavailable for this fit". Those cases keep the Gill objective-FD Hessian
     that 0.4.0 used.
  4. **Whether a sensitivity model is asked for at all.** `.admLoadSensModel()`
     returns `NULL` by design for a fixed-effects-only model, an ordinal
     endpoint, and mixed transformed/untransformed endpoints. Each used to run
     BOBYQA silently and briefly warned on every fit; they now emit a single
     plain message saying the gradient is finite-differenced.

  The startup line distinguishes the two analytic levels: `Grad: Analytical`
  means the struct thetas come from the second-order block, `Analytical (struct
  FD)` means the omega/sigma blocks are analytic and the struct thetas are not.

* **`linCmt()` models are supported at second order, by promotion.** `linCmt()`
  has no second derivative -- `rxFromSE()` cannot emit the nested `linCmtB`
  derivative, which is why nlmixr2est refuses `linCmt()` outright for its own
  analytic gradient and covariance. admixr2 instead promotes the model to its
  explicit ODE form with the exported `rxode2::linToOde()` and builds the
  second-order block from that. The promoted solve reproduces the analytic
  `linCmt()` prediction to 1.8e-08 relative.

  Only the `order = 2` request promotes: `admc`/`adgh` continue to use the fast
  solved form, which is all their first-order moments need.

* **Gill (1983) finite-difference steps (`gill = TRUE`), for the covariance
  Hessian and for the optimizer's gradient.** Every finite difference in admixr2
  took its step from the same heuristic -- `pmax(abs(p), 0.1) * h`, with `h` a
  fixed constant. That is a single guess about how much noise the objective
  carries, applied identically to every parameter, and it is the guess behind
  the "Hessian not positive definite ... try increasing `cov_h_outer`" warning:
  too small a step and the difference is swamped by solver noise, too large and
  it is swamped by curvature the second-order term does not model. A parameter
  the objective is flat in and one it is sharp in want different steps, and the
  right step moves with the ODE tolerance.

  `gill = TRUE` measures instead of guessing. It calls
  `nlmixr2est::nlmixr2Gill83()` -- exported, so no `:::`, and the same algorithm
  FOCEI uses for its own steps -- which probes the objective and returns, per
  parameter, the step where condition error and truncation error balance.
  Available on all four estimator controls. Try it first when the covariance
  comes back indefinite.

  Applied wherever admixr2 finite-differences the OBJECTIVE: the post-fit
  covariance Hessian, and the optimizer's gradient. The gradient steps are
  measured ONCE, at the starting values, and reused for the rest of the fit --
  the mechanism FOCEI uses, whose `numericGrad` runs Gill83 at its first
  evaluation (`nF == 1`) and finite-differences with the resulting `aEps`/`rEps`
  thereafter. Only the parameters actually differenced are probed: under
  `grad = "fd"`/`"cfd"` that is all of them, otherwise only those the sensitivity
  model could not supply a column for, and when it supplies all of them the probe
  is skipped entirely. A theta with no mu-referencing eta does **not** qualify by
  itself -- `.admBuildThetaSens()` emits a direction for each, so its gradient is
  exact and no step is chosen for it.

  Two places keep the fixed scale on purpose. admc's ETA perturbations
  (`eta_hi[, j] + h`, over `n_sim` rows) have no parameter to key on and
  difference the *prediction*, not the objective; under common random numbers the
  two solves share their draws, so most of the noise a step choice would guard
  against cancels anyway. And the gradient-FD covariance branch differences the
  analytic gradient, whose error is not what Gill83 measures.

  What it buys, measured on a 1-cmt ODE fit with `grad = "fd"`, against the
  analytical gradient: it **flattens** the error rather than uniformly shrinking
  it. Worst-parameter relative error `2.2e-03 -> 7.9e-04` (2.8x better); median
  `2.1e-04 -> 3.9e-04` (1.8x worse). The chosen steps span 20x across parameters
  (5.4e-05 to 1.7e-03) where the heuristic spans 3x, so it fixes the outliers --
  the structural thetas, which one shared constant suited badly -- and gives a
  little back on the parameters that constant already suited. For a quasi-Newton
  step the worst component is the binding one, but "more accurate everywhere"
  would be an overclaim. Cost: 20 objective evaluations for 5 parameters (~14
  NLL-equivalents), once; the end-to-end fit measured 4.4 s -> 4.1 s, i.e. the
  probe pays for itself rather than being a tax.

  Default `FALSE`, so nothing moves unless asked. Because it now reaches the
  gradient, `gill = TRUE` can change the fit and not only the standard errors --
  but only when the gradient has a finite-differenced parameter to begin with.

  One thing worth knowing if you read a step back: the acceptance rule keeps
  Gill83's `hf` for every return code except `"Not Assessed"`, `"Constant Grad"`
  and `"Odd/Linear Grad"`. `"High Grad Error"` is *kept*, because with the
  upstream default `fTol = 0` it is the ordinary return, not a failure -- a
  noiseless quadratic with an exact second derivative reports it. Rejecting
  everything but `"Good"` makes the whole option a silent no-op.

## Changes that can move an existing fit

Two changes in this release alter results for scripts that do not name a new
argument. Neither is a bug fix, so both are listed here rather than below.

* **`adirmcControl(grad = "fd")` now differences with `grad_h`, not a hard-coded
  `1e-6`.** The IRMC *inner* gradient ignored `grad_h` entirely -- it was the one
  finite difference in the package that could not be tuned, which is why the new
  `gill = TRUE` option could not reach it either. It now honours the argument, and
  under `gill = TRUE` takes Gill83's measured steps.

  **The default `grad_h` is `1e-4`, so an existing `grad = "fd"` adirmc fit
  differences with a step 100x larger than before** and can return a different
  objective and different estimates. Pass `grad_h = 1e-6` to reproduce 0.4.0
  exactly. Every other estimator already used `grad_h` here, so this also removes
  a discrepancy: the same control meant something different for adirmc than for
  the other three.

* **`adfoControl()`'s new `grad = "analytical"` default brings the `grad_bounds`
  box with it.** The box constraint (`p0 +/- grad_bounds`, default 5 on the
  optimizer scale, a factor of ~148 for a log-scale theta) applies only to
  gradient-based fits, so under the previous `grad = "none"` default an adfo fit
  was unconstrained. A default `adfoControl(studies = ...)` call is now confined
  to that box.

  This is rarely reachable -- it takes a starting value off by more than ~148x --
  and a fit that stops on the box now says so. Set `grad_bounds = Inf` for the
  old behaviour with the new gradient, or `grad = "none"` for the old behaviour
  entirely.

  The box itself is not an adfo peculiarity: `admControl()` (`grad = "sens"`) and
  `adghControl()` (`grad = "analytical"`) have always defaulted to a gradient
  *and* `grad_bounds = 5`, so this aligned adfo with them rather than singling it
  out. That is why the default stays and the REPORTING is what changed: the
  bounds notice is now emitted as a `message()` as well as a `warning()`.
  nlmixr2est muffles conditions inside `nlmixr2Est.*`, so the warning reaches
  `fit$warnings` -- where `print(fit)` surfaces it -- but never `warnings()`. A
  batch script that writes coefficients to disk without printing the fit would
  have seen nothing at all; the message goes to the same channel as the live
  progress table, which such a script does see.

## Bug fixes

* **adfo could report `NA` for every standard error on a fit that converged
  normally.** The driver decided whether to build the covariance Hessian by
  forward-differencing the *gradient* from the sensitivity model's shape alone,
  while `.adfoGrad()` re-derives that decision at run time with stricter
  requirements -- every study's cached `dJ` present, every theta's direction
  resolvable. When they disagreed, the gradient was itself finite-differenced and
  the Hessian then finite-differenced *that*, which is exactly the nested FD the
  gate exists to prevent. `.adfoGrad()` now reports what it actually did and the
  driver believes that, falling back to the Gill NLL-FD Hessian otherwise.

* **A joint (same-subject) study normalised before the model was known kept
  `NULL` block outputs.** Only the driver's pass carries the endpoint name, and
  the short-circuit for an already-normalised study skipped joint units
  altogether, so each block's `cmt` tag stayed empty: the joint sensitivity solve
  either dropped the fit to finite differences or read an untagged compartment,
  giving a finite but wrong joint objective with no warning.

* **`.admCacheWrite()` could delete another session's valid cache entry.** The
  cleanup that removes a half-written file ran on any `saveRDS` failure,
  including one that fails at open time and leaves a complete pre-existing entry
  untouched. A concurrent fit of the same model could therefore have its compiled
  model removed underneath its parallel workers, which then fail with
  "parallel restart N failed". It now only removes a file that call created.

* **A sensitivity model that failed to build was reported as quietly as one
  refused by design.** `.admLoadSensModel()` returns `NULL` both for models that
  cannot have one (no random effects, ordinal, mixed or unlike endpoint
  transforms) and for genuine failures such as an unwritable `rxode2::rxTempDir()`.
  The second case now warns rather than messages, and says what to check --
  previously the same script silently produced a coarser gradient, and different
  estimates and standard errors, on a machine with a read-only cache directory.

* **`adirmcControl()`'s `grad_h` now defaults to `1e-6`, not `1e-4`.** Making the
  IRMC inner finite difference honour `grad_h` was right, but inheriting the
  common default made every `grad = "fd"` adirmc fit converge on a step 100x
  coarser than the inner loop was tuned for. The IRMC inner NLL is deterministic
  given fixed proposals, so it wants a finer step than the sampling estimators,
  whose coarser default exists to step over Monte Carlo noise.

* **A fit that stops on the gradient box constraint now says so audibly.**
  nlmixr2est muffles conditions raised inside `nlmixr2Est.*`, so the warning
  reached `fit$warnings` -- where `print(fit)` shows it -- but never `warnings()`.
  A script that writes coefficients to disk without printing the fit saw nothing
  at all. The notice is now also a `message()`, on the same channel as the live
  progress table.

* **The order-2 `linCmt()` promotion did not run for a linCmt assigned to a
  variable, so adfo kept finite-differencing its structural thetas there.**
  `linCmt()` carries no second derivative, so an order-2 request promotes the
  model to explicit ODE form and builds from that. The gate detecting a
  solved-form model read `ui$predDf$linCmt`, and on rxode2 5.1.4 that column
  depends on how the model is *written*:

  | model line | `predDf$linCmt` | promoted before | promoted now |
  |---|---|---|---|
  | `linCmt() ~ add(a)` | `TRUE` | yes | yes |
  | `cp <- linCmt(); cp ~ add(a)` | `FALSE` | **no** | **yes** |
  | `cp <- 2 * linCmt(); cp ~ add(a)` | `FALSE` | no | no -- see below |

  The assigned form is the common way to write it, and there the promotion was
  never reached: `.admLoadSensModel(order = 2L)` served an order-1 model and adfo
  silently kept the forward-FD struct-theta pass (8e-04..1e-02 relative, against
  ~1e-09 for the analytic block). A correct fit, just the slow noisy one.

  Detection now uses the exported `rxode2::testRxLinCmt()`, which checks
  `ui$.linCmtM` as well and is `TRUE` for all three forms. The third still yields
  no cross block, because `rxode2::linToOde()` hands a derived `linCmt` back
  unchanged; the `linCmtB` text backstop then correctly refuses and the caller
  falls back to order 1. So this widens the fix rather than completing it.

* **`adghControl()` accepted an invalid nloptr algorithm, and would hand a
  derivative-free one a gradient.** It had its own two-line algorithm rule instead
  of the shared `.admResolveAlgorithm()` the other three controls use, and that
  rule was one-directional and unvalidated: `adghControl(algorithm =
  "NOT_AN_ALGO")` was accepted and surfaced as a cryptic nloptr error mid-fit,
  and `adghControl(grad = "analytical", algorithm = "NLOPT_LN_NELDERMEAD")` kept
  both -- paying for a gradient the algorithm discards on every iteration. It now
  goes through the shared reconciliation, so `grad == "none"` if and only if the
  algorithm is derivative-free, as documented. `algorithm` now defaults to `NULL`
  ("match `grad`"). `adgh`'s `cov_h_outer` default stays `eps^(1/4)` rather than
  the other three's `eps^(1/5)` -- that difference is deliberate, since the
  quadrature surface is noise-free.

  **One combination changes, and it is the one worth knowing about.** The old
  two-line rule existed to special-case exactly `"NLOPT_LN_BOBYQA"`, upgrading it
  to LBFGS whenever a gradient was requested -- because BOBYQA was `adghControl`'s
  own *default*, so naming it could not be distinguished from leaving it alone.
  With the default now `NULL`, naming a derivative-free algorithm is unambiguous
  and is honoured:

  | `adghControl(...)` | 0.4.0 | 0.4.1 |
  |---|---|---|
  | `grad = "analytical"`, `algorithm = "NLOPT_LN_BOBYQA"` | `analytical` + LBFGS | **`none` + BOBYQA** |
  | `grad = "fd"`, `algorithm = "NLOPT_LN_BOBYQA"` | `fd` + LBFGS | **`none` + BOBYQA** |

  So a script that explicitly restated the old default now gets a
  derivative-free fit where it had a quasi-Newton one. It says so (the
  reconciliation emits a message), but if you wrote `algorithm =
  "NLOPT_LN_BOBYQA"` meaning "the default", **delete the argument** -- `NULL`
  now picks LBFGS for you. Every other combination is unchanged, including the
  four pinned in `test-adgh-nodes.R`.

* **`adirmcControl()` validated neither `ci` nor `returnAdmr`.** `ci = 99` reached
  the interval columns as a nonsense level and `returnAdmr = "x"` made the
  driver's `isTRUE()` quietly `FALSE`, returning a full fit where a plain list was
  requested. Both are now checked, as in the other three controls.

* **A cache write that fails no longer discards the model it just compiled, or
  kills the fit.** Both disk caches wrote with a bare `saveRDS()`. The cache is an
  optimisation -- by the time it is written the model is compiled and loaded -- so
  an unwritable or full `rxTempDir()` (a locked-down HPC home, a cache directory
  owned by another user) should cost speed, not correctness. Instead
  `.admLoadModel()` propagated the error and failed the whole `nlmixr2()` call
  with `cannot open the connection`, while `.admLoadSensModel()`'s callers wrap
  the build in `tryCatch(error = function(e) NULL)`, so a write failure threw away
  a *successfully compiled* sensitivity model and dropped adfo from its order-2
  analytic structural gradient to forward FD, silently. Both now warn once per
  file and carry on with the model in hand. (The previous release swallowed the
  error entirely, which was also wrong -- the parallel restart workers find these
  models by reading exactly these files, so a silent failure resurfaced much later
  as every restart failing to read the cache.)

* **The session-ownership guard on a cached model rejected nlmixr2est's own
  sensitivity model unconditionally.** `.admRxLoadAll()` requires an artifact
  under a session-local `*Sens` build directory to belong to the running session,
  and tested that by comparing against `.admModDir()` -- which is
  `<tempdir>/admixr2Sens`, and so can never equal nlmixr2est 7.x's
  `<tempdir>/nlmixr2estSens`. A cached `.admSensFromInner()` result was therefore
  reported stale on every call *in the session that built it*, and the model
  recompiled (~3 s) for every fit -- the endless-recompile failure the guard
  exists to prevent, caused by the guard. It now tests membership of the current
  session's `tempdir()`, which is what identifies the session and covers both
  build directories.

* **The order-2 `linCmt()` promotion could write a theta's value into the wrong
  `THETA[k]` slot.** An order-2 request on a solved-form `linCmt()` model promotes
  it to explicit ODE form, and `.admBuildThetaSens()` numbers its emitted
  derivative directions from the *promoted* `iniDf`; `.admLoadSensModel()` built
  the `rename_map` that fills those columns at solve time from the *original* one.
  Any difference across the promotion -- a renumbered `ntheta`, an inserted or
  dropped row, a reordered eta -- meant each theta was differentiated in one slot
  and filled in another. The solve still succeeds and `use_d2` skips adfo's FD
  cross-check, so the fit would converge to wrong estimates and wrong standard
  errors with no error and no warning. Both are now derived from the same frame by
  construction (`.admSensNameMaps()`). Latent on every model measured here --
  `linToOde()` does preserve the `iniDf` on those -- but not guaranteed.

* **The gradient-box warning judged the fit against the wrong point, and stayed
  silent for the parameters most likely to need it.** It differenced the solution
  against the fit's `p0`, but each restart's box is centred on its own perturbed
  starting value, so a restart pinned to its box was not reported (its distance
  from `p0` never reaches `grad_bounds`) while an interior one could be reported
  spuriously. It also suppressed any hit on a parameter whose model declares a
  bound on that side -- a residual-error parameter always does -- even when that
  bound was nowhere near and admixr2's box was what stopped the fit. It now
  reconstructs the box actually given to nloptr, centred on the winning restart's
  own init, and reports a hit only when that box, rather than a model-declared
  bound, is the binding edge.

* **An explicit `adfoControl(grad = "analytical")` that cannot build a
  sensitivity model warns again.** 0.4.1 demoted this to a `message()`, which is
  right when `grad` was left at its (new) default -- an unavailable sensitivity
  model is routine and unactionable for a fixed-effects-only model or an ordinal
  endpoint -- but wrong when the user named the argument: a message is swallowed
  by `suppressMessages()`, by a knitr chunk with `message = FALSE`, and by any
  stderr-capturing wrapper, leaving no record that the fit used the gradient the
  control asked it not to use.

  Where it survives is worth stating precisely, because the obvious answer is
  wrong: nlmixr2est intercepts and muffles conditions raised inside
  `nlmixr2Est.*`, so this warning does **not** reach `warnings()` and
  `options(warn = 2)` does **not** turn it into an error. It is recorded on
  `fit$warnings`, which `print(fit)` displays. That is a durable record where a
  `message()` left none, which is the point -- but do not rely on
  `options(warn = 2)` to catch it.

* **Normalising a study twice no longer leaves its endpoint unset.** The
  idempotence guard added in this release returned before the point where a unit's
  `output` is filled from the caller's default, so a study first normalised
  without one (which is what the test fixtures do) kept `output = NULL`
  permanently -- and for a multi-endpoint model `.admBuildEvFull(tag_cmt = TRUE)`
  then has nothing to tag `cmt` with, so the unit reads the wrong compartment's
  trajectory. A second pass now fills what is still missing before returning.

* **Dev-mode parallel restarts could not see any function this release
  introduced.** `utils::assignInNamespace()` can replace a binding in a daemon's
  locked installed namespace but cannot *add* one, and the failure was swallowed,
  so a newly introduced helper was simply absent in the worker while everything
  looked healthy -- `.admGH()`/`.admGH0()`, called from every finite-difference
  site in `.admGrad()`/`.admGradBatch()`, would have broken every dev-mode
  `workers > 1` restart. Dev functions are now injected via a patch *environment*
  that patched closures are re-parented onto, so new and existing names resolve
  alike and a future helper needs no special handling. The same dispatch stopped
  shipping the parent's model-cache environments to each daemon: an environment
  serialises by value, so every dev-mode restart was copying compiled rxode2
  models that the worker has to rebuild anyway.

* **Generated models are built under role-tagged names, in their own directory,
  and a cached one is checked before it is trusted.** rxode2 names an anonymous
  model's `.c`/`.so` from the parsed model text alone, but the emitted C also
  depends on inputs that text cannot see -- above all the event-sensitivity code,
  which is injected afterwards. Two builds of one text that differ there land on
  a single artifact, and because entry points resolve by NAME (`R_GetCCallable`)
  a model bound to the earlier one silently starts executing the replacement
  (nlmixr2/rxode2#1171). admixr2 is exposed in the worst way of any package:
  `.admSensFromInner()` recompiles *nlmixr2est's own* inner model text with
  `eventSens = "jump"`, where nlmixr2est built the same text with a different
  one. Generated models now carry a name folding in the role and `eventSens`
  alongside the parsed md5, and are built in a session-local directory rather
  than the persistent `rxTempDir()`.

  The `.rds` caches stay in `rxTempDir()`, which persists, so a cache entry
  written by an earlier session necessarily references an artifact this session
  does not have. `.admRxLoadAll()` therefore checks that a cached model's DLL
  exists **and** belongs to this session's build directory, and reports the entry
  as stale otherwise so it is rebuilt. Both halves are load-bearing:
  `rxode2::rxLoad()` on a vanished DLL does not reliably error -- it returns
  quietly and the model then solves to garbage (a prediction frozen at its `t = 0`
  value, and `NA` structural gradients) -- and R removes its temp directory only
  on a clean exit, so a killed session leaves one behind that satisfies
  `file.exists()` indefinitely.

* **Normalising a study twice turned it into a joint (same-subject) study.**
  `.admNormaliseStudy()` was not idempotent, and the second pass changed the
  likelihood. Normalising a legacy single-output study ADDS an `observations`
  list while KEEPING its top-level `V` -- which is exactly the signature the
  joint branch tests for (`!is.null(s$observations) && !is.null(s$V)`), so a
  second pass collapsed it into one joint unit:

  ```
  pass 1  is_joint = FALSE      pass 2  is_joint = TRUE      pass 3  TRUE
  ```

  No error, no warning, and a perfectly plausible fit -- down a different
  likelihood path, and with adfo's `have_d2` forced `FALSE` (it requires
  `!any_joint`), so the order-2 analytical structural gradient this release adds
  quietly turned itself off. Each estimator normalises exactly once, so a normal
  fit never reached it; the test fixtures hand out pre-normalised studies that
  the driver then normalises again, which is how it was found -- meaning a number
  of end-to-end tests had been exercising the joint path while appearing to test
  the ordinary one. `.admNormaliseStudy()` now marks what it has normalised and
  returns such a study untouched. Genuine joint studies are detected exactly as
  before.

* **Non-finite parameters no longer reach the ODE solver.** The screen that
  rejects an unusable parameter vector before a solve tested that the omega
  diagonal was positive -- and `Inf > 0` is `TRUE`. A covariance probe that
  perturbs a residual parameter to `exp(1e5/2)` therefore handed `Inf` to
  `rxSolve()`, which integrated garbage and emitted on the order of 190,000
  `intdy -- t = <denormal> illegal` and `lsoda -- h too small` warnings before
  the caller discarded the result anyway. The parameter vector is now also
  checked for finiteness, at the objective *and* gradient entry points of all
  three affected estimators (the gradients unpack the optimizer vector
  themselves, so the objective's guard did not cover them). Aside from the
  console noise, this removes a guaranteed-useless ODE solve every time the
  optimizer or the covariance step overflows a parameter.

* **A cache-key collision solved fits at another model's fixed value.** A
  `fix()`ed parameter never reaches the optimizer, so it travels to the solve as
  data -- either baked into `$simulationModel` as that parameter's default, or
  carried on the cached sensitivity object. Both caches were keyed on the
  `model({})` block, which does not distinguish `theta <- fix(0.5)` from
  `fix(0.9)`. Two such fits therefore shared one compiled model, and the second
  silently solved at the first's fixed value: a plausible objective, plausible
  estimates and plausible standard errors for a model the user never wrote, with
  no error and no warning, persisting across sessions because the cache directory
  does. Fixed values now key both caches; ordinary *starting* values deliberately
  do not, since a starting value is optimizer state and keying it would force a
  recompile for nothing.

  The parallel workers are why this needed more than a longer key: a worker has
  no `ui` to re-derive from, and it recomputed the cache path itself. The parent
  now sends the path on `pinfo` (which travels by value, so no worker signature
  changes). One consequence for developers only: a `devtools::load_all()` parent
  and an older INSTALLED admixr2 derive that path differently, so a daemon
  started from a stale install cannot find the file -- run `devtools::install()`
  before testing parallel restarts, as the contributor notes already say. The
  worker's error message names that as the first thing to check.

* **A stale sensitivity cache entry could survive a change to what it caches.**
  The order-1 fallback is stored under the order-2 key, so editing what the
  order-2 build emits has to invalidate the entry -- otherwise the previously
  compiled model is served and `.adfoGrad()` contracts its second-order columns
  against the new code's direction map: a finite, plausible, silently wrong
  structural gradient with a normal-looking objective. That used to rest on
  editing a schema-tag string by hand, and was forgotten once during this work.
  The key now carries the package version *and* a digest of the emitter's own
  source, so an edit between releases invalidates it too, with nothing to
  remember.

* **`linCmt()` second-order promotion built its direction set from the
  pre-promotion model.** `.admBuildThetaSens()` swapped in the `linToOde()`
  model but kept the parameter rows derived from the original, despite a comment
  claiming otherwise. Any `iniDf` difference across promotion -- a renumbered
  `ntheta`, an added or dropped row, a different eta ordering -- would have made
  the emitted `rx_f1_THETA_k_`/`rx_f2_ETA_i_THETA_k_` differentiate a different
  parameter, and with `grad = "analytical"` now the default and the FD pass
  skipped, adfo would descend a gradient computed for the wrong theta. Latent
  rather than firing (promotion preserves the `iniDf` on the models measured
  here), but nothing enforced it.

* **A struct theta missing from the cached direction map crashed the fit.** The
  intended fallback -- turn the analytic pass off and finite-difference -- was
  unreachable, because `[[` with an unmatched name on an atomic vector throws
  rather than returning `NULL`. `.adfoGrad()` is not wrapped in a `tryCatch`
  there, so the whole nloptr run died with a bare `subscript out of bounds`.

* **A transformed endpoint no longer pays for second-order compartments it
  cannot use.** The solve paths deliberately discard the second-order block for
  an `lnorm`/`boxCox`/`yeoJohnson`/`logit`/`probit` endpoint (chaining a second
  derivative through the transform needs terms the first-order chain does not
  carry), but nothing stopped it being *built*: a 2-state, 2-eta, 1-unpaired-theta
  `lnorm` model integrated 20 states instead of 8 on every solve and then threw
  the extra away. Those endpoints now build the order-1 model directly.

* **A parallel worker no longer walks on from a model it could not load.** The
  worker's re-load step discarded its own failure return, so a model whose shared
  library had been unloaded was handed to the estimator as if live and the first
  `rxSolve()` dereferenced a dead pointer -- an opaque error deep in the restart,
  or a heap-corruption crash on Windows. It now stops with the cache path and the
  likely cause. A cache file whose contents are not a compiled model at all is
  also detected again and rebuilt, rather than being reported as loaded.

* **`.admNLL()` gained the non-finite screen the other estimators got.** admc's
  objective -- the function nloptr calls as `eval_f` -- still carried only the
  omega-diagonal test described below, which `Inf` passes, so admc users kept
  seeing the console flood that adfo and adgh users no longer do.

* **Compiled models are held in a session cache.** `.admLoadModel()` and
  `.admLoadSensModel()` reloaded their model from the disk cache on EVERY call --
  a `readRDS()` of a compiled model plus a `dyn.load()`, for every fit and every
  test -- and handed back a fresh wrapper object each time. A repeat load now
  costs a hash lookup instead: measured at roughly 50 ms to 0.5 ms.

  The mechanism is `nlmixr2est`'s, deliberately rather than invented: an
  `emptyenv()`-parented environment per purpose, a composite key covering
  everything that changes the emitted model, and a wholesale wipe at 64 entries
  to bound retained compiled models -- the same shape as its
  `.foceiAnalyticAugCache`. The load step matches `rxUiGet.foceiModel()` too,
  which re-loads EVERY `rxode2` element of a cached object rather than one by
  name.

  A cached model is only served while its disk cache file still exists, so
  `rxode2::rxClean()` still forces a genuine recompile, and the metadata the key
  cannot capture (a `boxCox`/`yeoJohnson` lambda's VALUE) is re-derived on every
  hit exactly as the disk path already did.

  Note what this does NOT change: memory. A session's footprint is set by how
  many DISTINCT models it compiles and loads -- measured at roughly 4-10 MB and
  two shared libraries each, and nothing unloads them -- so a session fitting
  many different models grows regardless of caching. Caching changes how often
  the same model is re-read, not how many are resident.

## Internal changes

* **`print.admFit()` reaches nlmixr2est's printer through `getS3method()`.** It
  used `get("print.nlmixr2FitCore", envir = asNamespace("nlmixr2est"))`, which is
  semantically a `:::` call that merely evades `R CMD check`'s syntactic scan --
  and carries exactly the upstream-refactor fragility the package's no-`:::`
  policy exists to avoid. The function is a registered S3 method, so method
  lookup is the supported public route to it.

* **The sensitivity-model builder takes an `order` argument.**
  `.admBuildThetaSens()`/`.admLoadSensModel()` default to `order = 1L` -- the
  existing first-order direction set, unchanged, which is what `admc`/`adgh`
  read. `order = 2L` additionally emits the eta x direction cross block that
  `adfo` needs. The block is deliberately asymmetric: `rxode2::rxExpandSens2_()`
  accepts two different direction sets, so no theta x theta compartment is
  generated, and admixr2 needs none of the residual-variance chains that
  dominate nlmixr2est's own second-order build (`errmodel.R` derives the
  residual analytically). Second-order initial conditions are emitted too,
  without which a parameter-dependent IC leaves the cross compartment at zero.

  The eta x eta half of that block is SYMMETRIC and `rxExpandSens2_()` does not
  know it: asked for the full rectangle it emits `d2/(d eta_1 d eta_2)` and
  `d2/(d eta_2 d eta_1)` as two variational compartments carrying the same
  equation. The block is therefore requested one eta row at a time, against only
  the directions at or after it, and the name matrix mirrors the duplicate cell
  onto the canonical one -- exact, since mixed partials of a smooth prediction
  commute, and it means the redundant chain expression is not emitted either.
  Saves `n_states * n_eta(n_eta - 1)/2` integrated states: 20 -> 18 on a
  2-state/2-eta/1-theta model, 63 -> 54 on a 3-state/3-eta/2-theta one. A joint
  fit now also asks for order 1, since `have_d2` excludes joint units and the
  cross block it used to compile was integrated on every solve for a result
  nothing read.

  The sensitivity cache key includes the order; see the cache-invalidation fix
  above for how a change to what the order-2 build emits invalidates an existing
  entry.

  A `fix()`ed parameter's VALUE now keys the cache too. A fixed parameter never
  reaches the optimizer, so it travels to the solve as data carried on the cached
  object -- and a parallel worker, which reads that file and has no `ui` to
  re-derive from, used the value it found. Two fits of the same model differing
  only in `theta <- fix(0.5)` versus `fix(0.9)` therefore shared one cache entry,
  and every parallel restart solved at the other fit's fixed value: silently, and
  across sessions, since the cache directory persists. Starting values
  deliberately do not key the cache -- they are optimizer state, and invalidating
  on them would force a recompile for nothing.

* **CI: `R-CMD-check` gained a `workflow_dispatch` trigger** and a dependency
  cache-version bump. The RcppParallel/TBB -> stringfish -> qs2 -> rxode2 stack
  has broken twice from CRAN-side rebuilds alone, with no commit of this
  package's involved, so being able to ask "does the current CRAN state still
  build?" without pushing a dummy commit is worth two lines. Pair it with a
  cache-version bump for a genuinely cold resolve -- a warm dependency cache is
  what made macOS look healthy right through the RcppParallel 6.0.0 break.

# admixr2 0.4.0

## New features

* **Student-t residual error (`cp ~ add(a) + t(nu)`) is now supported.** nlmixr2
  writes Student-t residuals as a *scale family* -- residual = scale * T_nu, with
  the scale being whatever `add()`/`prop()`/`pow()`/combined structure the
  endpoint already has -- so on the aggregate scale it is exactly the normal
  variance times `nu/(nu-2)`. admixr2 moment-matches it: the mean is unchanged and
  the variance is multiplied, which is exact for every existing residual form and
  works with all four estimators. Previously refused outright.

  **`nu` cannot be estimated from aggregate data and should be fixed**
  (`nu <- fix(5)`). This is structural, not a small-sample issue: `nu` reaches the
  aggregate moments only through the multiplier `nu/(nu-2)`, so it is aliased with
  the scale and only the product `a^2*nu/(nu-2)` is identified -- an estimated
  `nu` reflects its starting value, not the data. admixr2 warns when `nu` is left
  free. `t()` is intended for carrying a *known* `nu` through an aggregate
  analysis (for instance a published model supplied via `datagen()`), not for
  estimating tail weight; a fitted t model is observationally equivalent to a
  normal one with the same total residual variance.

  Note that moment matching makes the *mean* term of the objective correct, but
  admixr2 also scores `tr(V_pred^-1 V_obs)`, which treats the observed covariance
  as arising from a normal; that term stays approximate for heavy-tailed
  residuals, increasingly so as `nu` approaches 2.

* **The transform-both-sides transforms now call rxode2's own kernel.** `boxCox`,
  `yeoJohnson`, `logitNorm` and `probitNorm` used to be evaluated by a
  line-by-line R port of rxode2's C `_powerD`/`_powerDi` -- about ninety lines of
  branch order, clamps and short-circuits. They now call
  `rxode2::.rxTransform()`, which is what rxode2's own `boxCox()`/`yeoJohnson()`/
  `logit()`/`probit()` call and which bottoms out in the very C routine the solve
  transforms with, so admixr2 and rxode2 cannot drift apart. (The port had agreed
  with it exactly -- 0 mismatches over every transform x lambda x bounds
  combination -- but only by re-deriving it.)

  The residual quadrature now evaluates the transform for the whole node grid in
  one call instead of once per node, which made the switch a **speed-up rather
  than a cost**: 3.00 -> 0.60 ms per moment evaluation for `boxCox` at 81 nodes,
  1.50 -> 0.30 for `yeoJohnson`, 0.95 -> 0.30 for `logitNorm` (8 observations).
  The per-node accumulation is still a loop, deliberately, so the summation order
  and hence the objective are unchanged.

  Two pieces stay admixr2's own, for stated reasons: the derivative of the
  INVERSE transform (rxode2 exposes no equivalent -- and its `_powerDD` has a sign
  error on the Yeo-Johnson negative branch that admixr2 does not reproduce), and
  `dim()` restoration, which `.rxTransform()` drops.

* `.admBackTransform()`/`.admLogBackTransform()` now use `rxode2::probitInv()`
  instead of an inline `low + (high - low) * pnorm(p)`. Numerically identical
  across the whole range, but it is the kernel rxode2 itself transforms with, and
  it matches the neighbouring `expit` branch, which had always called rxode2.

* **New `resid_nodes` control argument** on all four estimators. A
  transform-both-sides endpoint (`boxCox`, `yeoJohnson`, `logitNorm`,
  `probitNorm`) has no closed-form mean and variance -- `y = g(h(f) + sigma*eps)`
  -- so admixr2 integrates the residual by Gauss-Hermite quadrature.
  `resid_nodes` sets that node count (default 81); every other error model has
  closed forms and ignores it.

  It is an **accuracy** dial, not a speed one. Worst-case relative error against
  an independent quadrature, over all four transforms and residual SD in
  {0.5, 1, 2, 3}: 5.7e-2 at 15 nodes, 4.5e-3 at 31, 5.0e-5 at 81. The error is
  dominated entirely by the largest SD; at SD <= 1 (a realistic residual on a
  transformed scale) 31 nodes already gives 1e-7 or better. Cost is linear in the
  node count in isolation (~50 us at 15, ~300 us at 81 for an eight-row study) but
  negligible beside the ODE solve -- a full NLL evaluation measured 0.750 s per 60
  evaluations at *both* 31 and 81 nodes. Raise it for a saturating endpoint with a
  large residual SD; there is little to gain by lowering it.

  `datagenControl()` takes it too, with the same default, so a study generated by
  `datagen()` and the fit that consumes it integrate the residual identically
  unless you deliberately change one of them.

* **New vignette: "Choosing a residual error model"**
  (`vignette("error-models", package = "admixr2")`). What the `cp ~ prop(...)`
  line actually does once your data are a mean and a covariance; a side-by-side
  fit of the same study with the right and the wrong residual model (the
  structural parameters survive, the IIV does not); the full menu of supported
  models; reading the covariance diagnostic panel, which is the aggregate-data
  substitute for a residual-vs-predicted plot; `resid_nodes`; the parameters
  aggregate data cannot identify (`t()`'s `nu`, an estimated `binom` size); and
  the combinations that are refused, with the reason for each.

## Bug fixes

* **Dropped the `qs2` dependency.** The compiled-model and sensitivity disk
  caches under `rxode2::rxTempDir()` are written with `saveRDS()`/`readRDS()`
  instead of `qs2`, and the files are named `adm-sim-*.rds` / `adm-sens-*.rds`.
  Base R serialization does this job, so the dependency bought nothing; it came
  up while tracking down rxode2 reverse-dependency failures in a check library
  that did not contain `qs2`. The caches are keyed by a model digest and live in
  the session temporary directory, so nothing needs migrating -- a leftover
  `adm-*.qs2` is simply a cache miss and the model is recompiled.

  Note this does not change which packages get *loaded*: `rxode2` itself imports
  `qs2`, and R loads a package's `Imports` with its namespace, so `qs2` (and
  `stringfish`) still enter the session behind `library(admixr2)`.

* **IRMC importance-sampling shift was wrong for every non-`exp` mu-referenced
  theta.** For a paired parameter `param <- h(theta + eta)`, `eta` and `theta`
  enter the transform through the *same* argument, so shifting `theta` by `Delta`
  shifts the importance-sampling target mean of `eta` by exactly `Delta` --
  `theta_new - theta_orig` -- for *any* `h`. The code computed that shift as
  `log(back(theta))`, which equals `theta` only for `exp` (`log(exp(theta))`);
  for a bounded (`expit`/`probit`) paired theta it used a natural-scale-log form
  and for an additive one (`emax <- temax + eta.emax`, the standard Emax writing
  style) it used `log(theta)`. Both biased the estimate and its analytical
  gradient, and the additive case went `-Inf`/`NaN` once the parameter passed
  through zero. Measured against a direct `adgh` evaluation, the `expit` shift
  drove the IRMC objective ~140 `-2LL` units off within a few tenths of the
  proposal point; the shift is now the identity `theta_new - theta_orig` for all
  transforms and matches the direct objective to importance-sampling noise
  (~0.06). Only `adirmc` fits used this path; the other three estimators integrate
  the random effects directly and were unaffected.

* **A `fix()`ed prediction-dependent residual lost its gradient.** A single
  endpoint whose only residual parameter is `fix()`ed -- `cp ~ prop(b)` with
  `b <- fix(0.2)`, or a fixed `lnorm`/`boxCox` coefficient -- is still
  prediction-dependent: `Var(y|eta)` moves with the prediction. `.admResidDeriv`
  early-returned `d(var)/df = 0` and `d(V_pred)/d(V_struct) = 1` whenever there
  were no *estimated* residual parameters, dropping that dependence from both the
  structural-theta and omega gradients, so under the default analytic gradient the
  optimizer descended a direction the objective did not follow. The early return
  now fires only for a genuinely additive residual (where those defaults are
  correct); every prediction-dependent form runs the full derivative even with no
  estimated residual parameter. FD-verified across `adfo`/`adgh`/`admc` for fixed
  `prop`, `lnorm` and combined residuals.

* **`binom(20L, p)` was refused as a non-constant size.** An integer literal
  deparses to `"20L"`, and `as.numeric("20L")` is `NA`, so a binomial size written
  with the integer suffix was misclassified as non-constant and refused with advice
  to `fix()` a parameter that does not exist -- the only difference from
  `binom(20, p)` being the suffix. A bare numeric literal (double or integer) is
  now read straight from the model AST.

* **A non-positive `nbinomMu` size now gives a clear domain error.** The size is
  estimated on the log scale, so a start `<= 0` made `log(size)` `-Inf`/`NaN` and
  the first NLL evaluation `NaN` with no explanation. It now refuses at parse with
  a domain message, matching the sibling `ar()` correlation and `t()` degrees-of-
  freedom guards.

* **`beta` precision denominator is guarded against a zero draw.** `.admSimulate`
  computed a `beta` endpoint's derived mean `b1/(b1+b2)` without the zero-
  denominator guard its three sibling solve paths already carry, so a draw with
  `b1 + b2 = 0` would have produced a `NaN` objective; it now floors the
  denominator the same way.

* **Standard errors: sigma SEs were uninitialised memory, and omega was excluded.**
  All three `CalcCov` functions built the Hessian over structural *and* residual
  parameters but returned only the structural corner. nlmixr2est's C++ `popDf`
  builder then read past the end of the matrix, so every sigma row of
  `parFixedDf$SE` printed a denormal (6.953178e-310, `%RSE` ~1e+307) instead of
  `NA` -- while the discarded sigma SEs were in fact good (reported SE / empirical
  sampling SD 0.90-0.94). The Hessian now also spans **omega**: excluding it made
  the *structural* SEs too small, because a theta carrying an eta is correlated
  with that eta's variance. Reported SE / empirical SD for the eta-carrying theta
  went from 0.67 to 1.17 (`prop`) and 0.67 to 1.06 (`lnorm`); a purely additive
  model barely moved. Under-stated SEs give over-confident intervals, so that was
  the dangerous direction of error. If the weakly-identified omega Cholesky makes
  the full Hessian indefinite, the struct+sigma sub-block is reported with a
  warning rather than nothing at all.

* **Omega and sigma standard errors are now reported, on the scale the estimates
  are printed on.** `fit$cov` previously covered the structural thetas alone.
  It now spans structural thetas, residual error *and* omega, delta-transformed
  out of the optimizer's parameterisation the way `nlmixr2est` does it -- so
  `Estimate +- 1.96*SE` is meaningful for every row. Residual error is reported as
  an SD (from `log(sigma^2)`), with the other `sigma_role`s mapped through their
  own derivatives (`t` degrees of freedom from `log(nu - 2)`, an `ar()` correlation
  from its logit, a negative-binomial size from its log). Omega is reported as the
  variance/covariance entries, named exactly as nlmixr2est names them
  (`om.<eta>`, `cov.<eta_i>.<eta_j>`).

  Omega is the one block that is not a per-row rescaling: the optimizer holds the
  log-Cholesky, `Omega = L L'`, and `d(Omega_ij)/d(L_ab)` is dense once omega is
  correlated. The new `.admOmegaJacobian()` builds that Jacobian in full and
  rotates both the omega block and its cross-covariance with struct/sigma; it
  agrees with a finite difference to 3e-10 on a correlated two-eta model.
  Calibration against the empirical sampling SD over 40 simulated studies gives
  reported SE / empirical SD = 1.13 for an IIV variance.

  Note for anyone touching this: nlmixr2est's C++ `foceiFitCpp_` re-dimnames the
  covariance from its own theta-name vector and blanks the omega rows (upstream
  ships `.impmapNameCov()` to repair the same thing for its importance-sampling
  estimator). It does so **in place**, which also blanks the driver's own copy, so
  the names are snapshot before the matrix is handed over and restored afterwards
  by `.admRestoreCovNames()`.

* **A printed standard error now belongs to the parameter it is printed beside.**
  nlmixr2est fills `parFixedDf$SE` *positionally*: it walks the thetas in `iniDf`
  order and takes the next entry of `sqrt(diag(fit$cov))` for each one it is not
  skipping. admixr2 builds its covariance in optimizer order -- structural thetas
  first, then residual error -- and hands over a matrix that also carries the
  residual parameters, so two things had to be said explicitly:

  - `.admCovThetaOrder()` puts the theta rows back in `iniDf` order. A model that
    declares its residual parameter first, `ini({ a <- 0.1; tcl <- log(3); tv <-
    log(30) })`, previously printed `a` with `tcl`'s SE, `tcl` with `tv`'s and
    `tv` with `a`'s -- a silent rotation, every number finite and plausible.
  - `.admCovSkip()` tells nlmixr2est which thetas the matrix actually carries,
    derived from the matrix itself rather than from a convention. nlmixr2est's own
    default is version-dependent: 6.2.0 skips only *fixed* thetas, while earlier
    versions (including 6.0.1, current on CRAN) also skip every residual-error
    theta, because FOCEI's covariance genuinely does not include them. Without
    this, those versions printed `NA` for every residual SD and read the structural
    SEs off the wrong rows.

  Verified on nlmixr2est 6.0.1 and 6.2.0, with the residual declared first and
  last, and with a `fix()`ed structural theta (which correctly stays `NA`).

* **Count endpoints could not be fitted with the default gradient.** `y ~ pois(cp)`
  and `y ~ nbinomMu(k, cp)` emit `rx_pred_ = llikPois(DV, ...)` -- the
  log-likelihood, not the mean -- and sensitivity columns that differentiate it,
  both of which need `DV`, which an aggregate fit does not have. `.adghGrad()` and
  `.admGradBatch()` returned all-`NA`, so `adgh` died at iteration 0 with
  "gradient of objective in x0 returns NA" and `admc` silently produced a **zero**
  Hessian and therefore no standard errors. admixr2 now emits sensitivities of the
  count MEAN (the distribution's argument), exactly as it already did for `beta`;
  gradients agree with a finite difference to 1.7e-05.

* **The covariance Hessian used the starting lambda for a transformed endpoint.**
  `.admGradBatch()` -- the evaluator behind `covMethod = "r"` when a gradient is
  available -- inherited the transform back-transform but not the
  estimated-lambda fix: an estimated `boxCox`/`yeoJohnson` lambda is a *sigma*
  name, so the zero-fill of the solve frame handed rxode2 lambda = 0 (a plain log
  transform) while the inverse used the model's STARTING lambda, held constant
  across every configuration. That is the same mismatch documented elsewhere here
  as making the sensitivity gradient ~60x wrong, driving the Hessian: every
  reported SE came from the gradient of a different function, and lambda's own row
  was insensitive to lambda. Each configuration now writes and inverts with its
  own lambda; measured against `.admGrad()` at a lambda well away from its start,
  the batch gradient went from 67% wrong to exact.

* **`beta()` endpoints were only ever right on the plain NLL path.** The
  prediction of `y ~ beta(b1, b2)` is the derived mean `b1/(b1+b2)` and its
  variance needs the SOLVED precision `phi = b1 + b2`, and every other path read
  the raw first shape parameter, or dropped `phi`, or both: the `covMethod = "r"`
  objective evaluator scored a different model from the fit, the finite-difference
  gradient returned all-`NA`, `datagen()` emitted an `E` that was a shape
  parameter and a `V` of `NA`s, and `plot()` gave an all-`NA` predicted covariance
  after a perfectly ordinary fit. None of it raised anything. Every path that
  turns a solve into a prediction now combines the pair and carries `phi`.

  A beta fit is also now driven **derivative-free**, with a message: a structural
  theta reaches the objective through `phi` as well as through the mean, and every
  gradient path chains through the mean alone. `datagen(method = "fo")` refuses a
  beta endpoint for the related reason that FO has no path to `phi` at all.

* **The `ar()` and `ordinal` guards judged every study, not the affected one.**
  Both decided from a model-level scan and then rejected every flattened unit, so
  `cp ~ add(a) + ar(rho); ct ~ add(a2)` refused a `ct` study whose `V` happened to
  be diagonal, and a PK + ordinal model could never be fitted at all -- the
  ordinary `cp` study is neither joint nor supplies one block per category. Each
  guard now looks only at units that observe the endpoint it is about.

* **Ordinal categories were grouped by exact floating-point time equality.** The
  row times come from the per-category blocks, i.e. from independent user inputs:
  `seq(0.1, 0.7, by = 0.2)` and `c(0.1, 0.3, 0.5, 0.7)` are the same grid to a
  reader and differ in the last bit to `match()`, which put the two categories in
  different groups and silently dropped the `-p_j*p_k` cross-covariance for those
  rows -- the term a joint ordinal fit exists to capture. Grouped by tolerance now.

* **The moment expansion and its derivative capped the same pole differently.**
  `.admMomF()` (what the NLL scores) caps the divergent `mu^(k-2)` correction
  against the leading term; `.admMomFd()` (what the gradient chains through) zeroed
  it past a magnitude threshold instead. For `pow(b, c)` with `c < 1` near a zero
  prediction the two differed by orders of magnitude, so the optimizer was handed a
  direction that does not descend the function it is minimising. The derivatives
  are now the derivatives of the capped expression, piecewise, and agree with a
  finite difference of `.admMomF()` across the capped and uncapped regimes alike.

* **A parallel worker could invert a transform with another model's lambda.** The
  sensitivity cache key covers the `model({})` block, the `iniDf` names, the
  `fix()` flags and the `err` column -- but not the estimates, so two models
  differing only in the VALUE of a `fix()`ed lambda share one file. The parent
  re-derives `pred_tbs` on a cache hit; the worker could not, and used the file's.
  The parallel restarts then minimised a different objective from the sequential
  ones, invisibly, because the NLL itself is bit-identical. The worker now
  re-derives it from `pinfo`, which it already holds.

* **`plot()` back-transformed three residual roles on the wrong scale.** The trace
  panel special-cased `pow_exp` and `t_df` and let `ar_cor`, `nb_size` and
  `tbs_lam` fall through to the generic `exp(v/2)` variance rule: a converged
  `ar()` correlation of 0.6 plotted as **1.22**, outside its own support and
  disagreeing with what `print(fit)` reports. The display map now comes from
  `.admSigmaNat()` itself, so a new `sigma_role` cannot be added in one place and
  forgotten in the other.

* **A `binom` size written as a model constant was refused as non-constant.**
  `nt <- 20; y ~ binom(nt, p)` -- a genuinely constant number of trials, and how
  one is usually written -- hard-errored with advice to `fix()` a parameter that
  does not exist. A bare numeric assignment in `model({})` now resolves; an
  estimated size is still refused, since it has no gradient path.

* **A count or beta endpoint alongside another endpoint is now refused.**
  Multi-endpoint solves route observations by compartment, and a count endpoint is
  read through its distribution's ARGUMENT -- a model variable, not a compartment
  -- so the tagged records matched nothing and the objective came back `Inf` with
  no explanation. Relatedly, the dummy frame handed to nlmixr2est now carries the
  ENDPOINT names rather than the solve columns, which is what its `dvid`->`cmt`
  translation expects.

* **`datagen()` refuses an ordinal endpoint** instead of emitting a study without
  the cross-category covariance: its categories are one joint observation, and
  `datagen()` derives each observed output separately.

* **`resid_nodes` no longer changes what a positional call means.** It was added
  as the SECOND argument of every estimator control, so `adghControl(studies, 7L)`
  -- which had always meant `n_nodes = 7` -- set `resid_nodes = 7` instead, passed
  its own validation, and left `n_nodes` at its default, changing the eta
  quadrature grid the whole fit is built on with no message. `admControl(studies,
  20000L)` and `adirmcControl(studies, 2000L)` were the same for `n_sim`, and
  `datagenControl` shifted `sampling`/`seed`/`cores`. `resid_nodes` is now the
  last argument of all five, where it cannot capture a positional one.

* **The ordinal same-time grouping is now defined once.** `.admResidApply()`,
  `.admResidVChain()` and `.admResidMuCoupling()` each grouped the rows
  themselves. When the tolerance-based grouping above was first added, it went
  into one of the three -- which is worse than the exact-match bug it replaced:
  wrong-but-consistent became objective-and-gradient-disagree, so the optimizer
  descended a direction the objective does not follow. All three now call
  `.admOrdTimeGroup()`.

* **Endpoints transformed differently from one another refused the sensitivity
  model.** The guard caught transformed-vs-*un*transformed mixtures only, while
  the back-transform spec is a single one taken from the first endpoint. So
  `cp ~ lnorm(a); ct ~ boxCox(b, lam)` applied `exp()` to `ct`'s Box-Cox rows, two
  `logitNorm` endpoints with different bounds shared the first one's bounds, and
  two `boxCox` endpoints shared the first one's lambda -- the residual path being
  per-endpoint already, the gradient then described a different function from the
  one the objective scored. Any non-identical set now falls back to finite
  differences.

* **A joint (same-subject) study had no aggregate diagnostics.** `.admAggData()`
  solved one output for the whole stacked unit and applied the first endpoint's
  residual spec to every row; it then died on the dimnames (the row count is the
  stacked total, the labels were one block's times) and, being guarded, left
  `fit$env$aggData` unset -- so `plot(fit)`'s mean/cov panels had nothing to show
  and said nothing about it. It now uses the estimators' own shared-eta solve and
  per-row-output residual, and labels rows `<endpoint>@<time>`.

* **Documented: an `adfo` standard error describes scatter, not accuracy.** FO
  linearises at eta = 0, so on a non-additive residual or a large omega the point
  estimate carries a bias of several standard errors -- measured 5-20 SE, giving
  0% coverage for a nominal 95% interval even where the SE itself matches the
  sampling SD. Prefer `adgh`/`admc` when the uncertainty matters.

* **A failed covariance is no longer silent.** When the Hessian was singular the
  covariance came back `NULL`, `covMethod` was set to `""` and every SE was `NA`
  with no warning reaching the user. The drivers now say so.

* **A study `ev` containing observation records now warns.** `ev` is dosing-only;
  observation rows in it were appended a second time by the study's own `times`,
  silently duplicating every time point.

* **Residual parameters fixed with `fix()` were silently dropped.** `add(a)` with
  `a <- fix(0.7)` fitted with **no residual variance at all**; `add(a) + prop(b)`
  with a fixed `b` lost the proportional term; and `pow(b, c)` with a fixed `c`
  reverted to `prop()`. `.admParseIniDf()` removes fixed rows from the optimizer,
  and the residual spec indexed only the estimated ones -- `tdf_fixed`,
  `ar_fixed` and `lam_fixed` existed for exactly this reason but `add`/`prop`/`pow`
  had no equivalent. They now carry `add_fixed`/`prop_fixed`/`pow_fixed`. Fixing a
  residual parameter is routine (it is what this package's own `t()` advice tells
  you to do for `nu`), so this was reachable in ordinary use.

* **A `prop()`/`pow()` term on a transform-both-sides endpoint contributed
  nothing.** For `boxCox`/`yeoJohnson`/`logitNorm`/`probitNorm` the quadrature used
  only the additive parameter, so `cp ~ add(a) + prop(b) + boxCox(lam)` scored
  identically with and without `b`: the parameter entered the optimizer, had an
  exactly-zero gradient, and was reported back at its starting value. rxode2 emits
  `rx_r_ ~ (a)^2 + (rx_pred_f_)^2*(b)^2` for that model, and admixr2 now builds the
  transformed-scale residual SD from the same expression, including `propT()`/
  `powT()` (which scale by the transformed prediction) and `combined1()`.

* **The post-fit covariance was a Hessian of the wrong objective for several error
  models.** `.admNLLBatch()` -- the evaluator `covMethod = "r"` differentiates --
  called the fused C++ kernels unconditionally. Those implement additive,
  proportional, combined and lnorm only, so transform-both-sides, count, beta,
  ordinal and `ar()` models were scored as `combined2`: standard errors and RSEs
  came from a different model than the one fitted (measured on a boxCox model,
  190.28 against 49.46). It now applies the same `.admResidCppOK()` gate `.admNLL()`
  uses. `adirmc` cannot take that route (its kernel forms the importance-weighted
  mean internally) and now refuses those models with a message.

* **`adfo` dropped `ar()` from its objective while keeping it in the gradient.**
  `.adfoVpred()` never received the observation times and never added the residual
  correlation, so the FO objective was exactly invariant in `rho` while
  `.adfoGrad()` returned a non-zero `rho` gradient -- the optimizer walked a
  direction the objective could not move along, and `adfo` reported a different
  objective from `adgh`/`admc` on identical data. Relatedly,
  `adfoControl(grad = "analytical")` warned that it was falling back to finite
  differences when no sensitivity model was available but did not actually do so.

* **An out-of-support transform aborted the whole fit.** `any(ap$ms != 1)` was not
  NaN-guarded in nine places. `.admTBSi()` legitimately returns `NaN` outside a
  transform's support, and the default `grad_bounds = 5` lets a line search reach
  it, so `any(NaN != 1)` -- which is `NA` -- raised "missing value where
  TRUE/FALSE needed" instead of the optimizer simply rejecting the point. This
  killed every `yeoJohnson` fit.

* **The sensitivity-model cache could serve a stale transform spec.** The cache key
  digests `ui$lstExpr`, the `model({})` block only, but a Box-Cox lambda's starting
  value and its `fix()` status live in `ini({})` -- so `lam <- fix(0.5)` and
  `lam <- 0.5` collided. `pred_tbs` is what tells the solve which lambda to use and
  how to back-transform, and it was not re-derived on a cache hit (unlike
  `rename_map`/`fixed_theta`, which are, for the same reason). Gradients came back
  wrong by 10^2-10^4x with one component of the wrong sign, while the objective
  stayed bit-identical, so nothing warned and the fit simply stalled.

* **`0^negative` in the moment expansion.** `pow(b, c)` with `c < 1` at a structural
  prediction of exactly zero -- routine for a depot model observed at `t = 0` --
  produced a **negative variance** (measured -3.4e+20 at `c = 0.25`) or a
  plausible-looking 2.3e+05 at `c = 0.75`. The second-order term has a genuine pole
  there and is now dropped rather than evaluated at machine epsilon. The C++ twin
  `adm_mom_f()` had no guard at all and returned `NaN` where the R path returned a
  finite value, so the same model fitted or did not depending on the estimator.

* **`ordinal` endpoints are now supported** (`y ~ c(p1, p2)`), as a joint
  same-subject unit with one observation block per category. The spec is registered
  under every category probability (only the first was, leaving the others with no
  residual variance), and the same-time cross-category covariance correctly
  *replaces* the structural covariance rather than adding to it -- by the law of
  total covariance `Cov(1_j, 1_k) = -E[p_j]E[p_k]` exactly, the structural term
  cancelling. Verified against a multinomial simulation with between-subject
  variability.

* **`dv()` is now refused.** It scales the residual by the observed DV, an
  individual-level quantity an aggregate mean and covariance cannot recover.
  rxode2's simulation ignores `dv()`, so admixr2 had been silently fitting the
  prediction-scaled model instead.

* **`ar()` combined with `prop()`/`pow()`/combined is now refused**, as is `ar()`
  inside a joint multi-output study. rxode2's innovation scaling leaves the marginal
  variance equal to `rx_r_` only when `rx_r_` is constant; with a prediction-
  dependent variance the process is non-stationary and admixr2's covariance was
  measured 2.4-12x too high.

* **Known upstream issue -- simulating an `ar()` fit will not reproduce its
  covariance.** rxode2 has two `ar()` emitters and they do not agree with each
  other. Its *estimation* lines are the prediction-error decomposition
  (`rx_pred_ + phi*prev_resid`, `rx_r_ * (1 - phi^2)`), whose implied marginal
  variance is the stationary AR(1) admixr2 scores. Its *simulation* is not
  stationary when a dose record precedes the first observation: the first
  observation carries up to 2x the nominal residual variance. A **zero-amount**
  dose reproduces it and a plain `add()` model does not, so it is record-driven
  and specific to `ar()`; nlmixr2's own focei cannot recover `rho` from rxode2's
  own simulation either (0.4617 against a truth of 0.60 on individual-level data,
  with no admixr2 involved). admixr2 keeps the stationary form -- matching the
  simulator would put it at odds with nlmixr2's estimator and would break when
  this is fixed upstream. Every other error model round-trips (simulate from the
  fitted model, aggregate, and recover the fitted mean and covariance) to within
  Monte-Carlo noise.

* **Prediction-dependent residual error is now composed correctly (`prop()`,
  `pow()`, `lnorm()`, combined).** admixr2 built the predicted covariance as
  `Var_eta(f) + Sigma(mu_pred)` -- evaluating the residual variance at the
  *population mean* prediction rather than averaging it over individual
  predictions. That is exact only for additive error. The predicted covariance is
  now the **law of total variance**,
  `Var_eta(E[y|eta]) + E_eta[Var(y|eta)]`, which for a proportional model adds the
  previously missing `b^2 * Var_eta(f)` to the diagonal, and for `lnorm()` also
  scales the **off-diagonals** by `exp(s)` (its conditional mean is `f*exp(s/2)`,
  so the whole covariance is scaled, not just its diagonal). Validated against
  individual-level simulation: the old formulas carried fixed biases of ~15-20%
  that did not shrink with sample size, while the new ones converge to the
  empirical moments.

  **This changes results for every `prop()`, `pow()` and `lnorm()` model.**
  Objective values, residual-error and IIV estimates and all standard errors move
  -- for a proportional model with 30-50% IIV, the residual SD and omega were both
  biased upward by roughly 2-4%; for `lnorm()` the effect is larger. Purely
  additive (`add()`) models are unchanged, bit for bit. Refits are expected to
  differ from results produced by earlier versions.

* **`lnorm()` analytic gradients were computed against the wrong quantity.** For a
  log-transformed endpoint the sensitivity model returns `rx_pred_ = log(f)` while
  the NLL path reads the natural-scale prediction, so `grad = "sens"`/`"analytical"`
  differentiated `log(f)` while the objective scored `f`. The sensitivity paths now
  back-transform with the chain rule. This affected every `lnorm()` fit using an
  analytic gradient and went unnoticed because `lnorm()` appeared in no gradient
  test; a finite-difference gradient check across all estimators and error models
  has been added.

* **`delay()` (DDE) models get an accurate sensitivity solve.** A delay model's
  sensitivity system -- the base ODEs plus one variational compartment per state
  per direction, all delayed -- is stiff enough to trip rxode2's `hasDelay`
  AutoSwitch composite (`dop853`+`ros4`) into its `ros4` leg, whose dense
  delay-history is inaccurate for this system. The failure is silent: the
  sensitivity model's predictions match the ordinary solve for the first
  observations and then drift once `delay()` begins reading the recorded (solved)
  history, so `grad = "sens"` gradients on a DDE model could be wrong without any
  error or warning. Sensitivity solves for a delay model are now forced onto pure
  `dop853` (dense, no `ros4` secondary), whose 8th-order dense output reproduces
  the ordinary solve. Non-delay models are untouched, and their solves are
  unchanged byte for byte. Found by porting the equivalent fix from nlmixr2est's
  own augmented-sensitivity solve.

## Internal changes

* **The post-fit covariance's reported-scale rotation and its non-PD omega
  fallback are now single shared helpers.** The ~46-line block that rotates the
  optimizer-scale covariance onto the printed scale (residual delta factors plus
  the omega Jacobian) was byte-identical in all three `CalcCov` functions, and the
  "drop to the struct+sigma sub-block when omega makes the Hessian indefinite"
  fallback was duplicated in `adfo`/`admc` with an already-divergent invert-first
  variant in `adgh`. Both are now `.admScaleReportedCov()` and
  `.admReduceNpdOmega()` in `utils.R`, so a change to how residual/omega SEs reach
  the printed scale, or to the fallback threshold, is made in one place rather than
  three. The `adgh` fallback converges onto the same eigenvalue threshold the other
  two use; results are unchanged (the full pipeline and covariance suites pass
  identically).

* **The residual variance's dependence on `(mu, var_f)` is computed once per
  study/unit instead of three times.** `.admResidVChain()`, `.admSigmaGrad()` and
  `.admResidMuCoupling()` each recomputed `.admResidDeriv()` internally, in every
  estimator's hot gradient loop -- three `resid_nodes` (default 81) quadratures per
  observation row for a transform-both-sides endpoint. They now accept the
  precomputed derivative as an optional last argument, which the estimators (which
  call all three on the same inputs) pass, cutting that to one. Gradients are
  bit-identical (the same computation, reused); the optional argument defaults to
  recomputing, so every other caller is unchanged.

* **The residual V-composition tail is one helper, `.admApplyResidTail()`.** The
  three-line `V <- V * tcrossprod(ms); diag(V) <- dv; V <- V + rmat` that composes a
  structural covariance with the residual (lnorm/TBS off-diagonal scale, the
  composed diagonal, an `ar()` correlation matrix) was hand-copied at eleven sites
  across every estimator's moment/objective path, `plot.R`, `datagen.R` and
  `.admJointResidual`. Adding an off-diagonal residual channel meant editing all of
  them, and missing one silently dropped that endpoint's off-diagonal predicted
  covariance on that path. It is now written once, including the load-bearing
  `na.rm` guard that keeps a NaN from a transform's out-of-support tail from
  aborting the fit. Objective and gradients are bit-identical.

# admixr2 0.3.0

## New features

* **Analytical gradients for non-mu-referenced ("unpaired") structural thetas.**
  A structural theta with no mu-referencing eta (`tka` with no `eta.ka`, or the
  `exp(tcl) * exp(eta.cl)` writing style rxode2 does not mu-reference) used to
  cost an extra finite-difference `rxSolve` per gradient call. admixr2 now emits
  its own first-order sensitivity model over an explicit direction set (one
  direction per random effect plus one per unpaired theta), compiled with
  `eventSens = "jump"` so dosing-modifier (`f`/`lag`/`rate`/`dur`) sensitivities
  are no longer silently zero. This mirrors the scheme nlmixr2est's fast-focei
  uses (`.foceiAnalyticDirections`) but first-order only, and is cross-validated
  against nlmixr2est's inner model to ~1e-13 across ODE, linCmt, dosing
  modifiers, initial conditions, covariates, if/else and multi-endpoint models.
  Consumed by `admc`, `adgh` (including joint multi-output studies); `adfo` keeps
  finite differences (its `V_pred = J Omega J' + Sigma` needs a second
  derivative). Measured 2.5-3.8x faster and ~100x more accurate than the previous
  finite-difference path on a 2-compartment model. This adds `symengine` (already
  a hard dependency of `nlmixr2est`, so always installed alongside admixr2) to
  `Imports`, used to emit the linCmt direction derivatives. The feature degrades
  gracefully on rxode2 without `eventSens = "jump"` support (it falls back to the
  finite-difference path), so no minimum-version bump is required.

* **Residual error models: `pow()`, `addPow()` and `combined1()` are now
  supported, with analytical gradients** (#84). admixr2 previously supported
  only `add`, `prop` and `lnorm`. The residual error model is now read from
  `ui$predDf` (`errType`/`errTypeF`/`transform`/`addProp`) rather than from
  `iniDf$err` alone, and every estimator evaluates it through one shared
  specification:

  | form | variance |
  |------|----------|
  | `combined2` (default for `add + prop`) | `a^2 + b^2 * f^(2c)` |
  | `combined1` | `(a + b * f^c)^2` |
  | `lnorm` | moment-matched lognormal |

  with `c = 1` recovering `prop` and `b = 0` recovering `add`. Analytical
  `d(var)/d(sigma)`, `d(mu)/d(sigma)` and `d(var)/d(f)` are supplied for all of
  them, so residual parameters keep an exact gradient under
  `grad = "sens"`/`"analytical"`.

  Existing `add`/`prop`/`lnorm` fits are unaffected: the aggregate `-2LL` is
  bit-for-bit identical, and their gradients change only by floating-point
  reassociation (~1 ulp).

* **Multi-compartment fitting (multiple observed outputs).** A study may now
  observe several model outputs at once (e.g. plasma and brain/CSF) via an
  `observations` list -- one entry per observed output with its own `output`,
  `times`, `E` and `V`. Two modes (#85):
    * *Independent* -- each output has its own `n`/`ev` (separate experiments,
      e.g. literature meta-analysis); the aggregate `-2LL` is the sum of the
      per-output likelihood blocks. Fit with full analytical / sensitivity
      gradients.
    * *Joint (same subjects)* -- outputs measured on the same subjects, with a
      shared `n`/`ev` and a joint covariance given either as a study-level full
      `V` or as per-output marginal `V` plus a `cross` list of cross-covariance
      blocks. Scored by a single MVN over the stacked vector with shared random
      effects and the full **analytical** gradient in all three estimators (any
      number of compartments; the assembled joint covariance is checked for
      positive-definiteness).

  Supported by `est = "admc"`, `"adfo"` and `"adgh"`; `datagen()` generates
  multi-output aggregate data and `plot()` renders one panel set per compartment.
  Pass the endpoint names to `admData()`, e.g. `admData(c("cp", "cCSF"))`.
  `est = "adirmc"` does not support multiple observed outputs.

* **Parallel restarts now run on `mirai` daemons.** `workers > 1` starts a pool
  of background R processes instead of dispatching through `future`/`furrr`.
  This replaces the previous fork (Unix/macOS) vs PSOCK (Windows/RStudio) split
  with a single code path that behaves identically on every platform, and the
  pool lives on its own mirai compute profile so it never disturbs daemons the
  user has set up for their own code. `furrr` and `future` are no longer used;
  `mirai` moves into `Suggests`. Workers are still stopped automatically after
  the restart phase (and now also on error/interrupt, via `on.exit()`), so all
  cores are free for the covariance step; `admStopWorkers()` remains available.

* **`nDisplayProgress` control argument** for every estimator (`admControl()`,
  `adfoControl()`, `adghControl()`, `adirmcControl()`), passed through to the
  `rxSolve()` calls that drive fitting. It sets how many subjects a single solve
  must exceed before the solver shows its text progress bar. The default
  (`.Machine$integer.max`) keeps the bar off, so it no longer leaks into scripts,
  logs or rendered vignettes; lower it (e.g. `1000L`) to watch progress during
  long interactive fits.

* The aggregate-data estimators (`adfo`, `adgh`, `adirmc`, `admc`) now carry
  `type` and `description` attributes classifying them as "Model Based Meta
  Analysis" methods, so they appear in the category-grouped estimation-method
  list nlmixr2est prints for an unsupported `est=` (or a bare `nlmixr2()` call)
  (#107).

## Bug fixes

* **`pow()` models no longer fit the wrong residual model, silently.** `pow(b, c)`
  produces two `iniDf` rows -- the coefficient (`err = "pow"`) and the *exponent*
  (`err = "pow2"`). admixr2 recognised neither, warned once, and then treated
  **both as additive variances**: the exponent was stored as `2*log(c)` and
  optimized as a variance contributing `exp(2*log(c))` to `diag(V)`. A `pow`
  model therefore ran to completion and reported plausible estimates for a model
  it was not fitting. Residual parameters now carry a role, and a `pow` exponent
  is estimated on its own (unconstrained, identity) scale.

* **`combined1()` is honoured.** `predDf$addProp` selects SD-additive
  (`combined1`) versus variance-additive (`combined2`) residual error. admixr2
  ignored it and always computed `combined2`, dropping the `2*a*b*f` cross term.
  (`combined2` is nlmixr2's default, so only models that explicitly asked for
  `combined1()` were affected.)

* **An unrepresentable residual model is now refused rather than approximated.**
  Error types admixr2 cannot express as a Gaussian aggregate MVN (`logitNorm`,
  `probitNorm`, Box-Cox/Yeo-Johnson transforms, `t`/`cauchy`, `propF`/`powF`)
  previously emitted a one-time warning and were then **treated as additive**,
  so the fit proceeded with the wrong residual model. They now `stop()`. This is
  a behaviour change: a model that "worked" before may now error.

* **`propT`/`propF`, `norm`/`dnorm` and `dlnorm`/`logn`/`dlogn` no longer emit
  spurious "modelled as ..." approximation warnings.** These are aliases, not
  approximations: `norm` *is* `add`, `logn` *is* `lnorm`, and on an untransformed
  model `propT` (which scales by the transformed prediction) *is* exactly `prop`,
  because there the transformed and untransformed predictions are the same
  quantity. The warnings claimed an inaccuracy that did not exist.

* **Lognormal residual error is now applied to the plotted predicted mean.**
  `plot.admFit()`'s aggregate-data helper added the lnorm variance to the
  predicted covariance but never applied the `exp(s/2)` mean scaling to the
  predicted `E`, so lnorm fits plotted a mean the NLL does not use.

* The solver progress bar no longer appears during covariance/gradient batches.
  Most internal `rxSolve()` calls already suppressed it, but the covariance and
  batched-gradient solves in `admc` hard-coded a low `nDisplayProgress` (1000),
  so the bar printed once a chunk exceeded 1000 solves. All solves now honour the
  new `nDisplayProgress` control argument (default off).

* Hard-coded numeric constants in a model's `model({})` block (e.g. a fixed brain
  volume `vb <- 5`, common in PBPK/CNS models) are no longer zeroed. admixr2 used
  to hand-fill every model parameter it did not set with `0`, clobbering such a
  constant's default and producing an `NA`/non-finite objective (e.g. a
  `qout / vb` divide-by-zero). It now supplies only the parameters it varies and
  lets `rxSolve()` fill the rest from the model's own defaults, so constants and
  covariate defaults keep their value.

* `adgh` now computes gradients for non-mu-referenced (unpaired) structural
  thetas. The unpaired-parameter set was derived from the eta-indexed
  `struct_eta_idx`, so it was always empty and those thetas silently received a
  zero gradient; it now uses the struct-indexed `struct_has_eta`.

* **Parallel restarts under `devtools::load_all()` warn once about the installed
  package.** In dev mode the admixr2 namespace is locked, so worker daemons run
  the *installed* package rather than the loaded source; if it is stale the
  parallel objective silently diverges from the sequential one. `.admRunRestarts`
  now emits a one-time warning in this case telling you to `devtools::install()`.
  It never fires in production (installed package == source).

## Internal changes

* **`adgh` gradient-mode fits are about twice as fast: the objective and the
  gradient now share one solve** (#76). `nloptr` asks for the objective and the
  gradient as two separate calls, but LBFGS always asks at the same parameter
  vector, and `.adghGrad` already builds exactly the moments the negative
  log-likelihood needs -- so the objective's solve was duplicate work. It is now
  memoised onto the gradient's solve. Measured on a 3-compartment, 5-eta,
  40-timepoint fit with a full covariance `V`: `rxSolve` calls per fit drop from
  58 to 23 (`n_nodes = 3`) and 63 to 25 (`n_nodes = 5`), roughly halving wall
  time. Applies to `grad = "analytical"` only (including multi-restart fits);
  `grad = "fd"`/`"cfd"`/`"none"` are unchanged, as are all gradient values.

  Note for anyone comparing objectives across versions: the reported objective
  now comes from the sensitivity solve rather than the plain one. Both integrate
  the same underlying model, but the augmented system makes rxode2's adaptive
  stepper land a little differently -- about 5e-11 relative on the objective,
  well inside the solver's own tolerance, and parameter estimates are unchanged
  (identical to six decimal places in testing). As a side effect the objective
  and its gradient are now computed from a single trajectory, where previously
  they came from two slightly different ones.

* **Model loading and per-fit memory now follow nlmixr2est's own conventions.**
  admixr2 previously pinned each fit's `foceiModel` companion objects in a
  package-level environment (a Windows GC-finalizer heap-corruption guard) and
  reclaimed rxode2's global model registry with a bespoke snapshot/teardown after
  every fit. Both are gone: the companion objects are no longer pinned (the guard
  proved unnecessary -- verified by running the `covMethod = "r"` fit path
  repeatedly under aggressive GC with no crash), and each estimator now frees
  memory the way nlmixr2est does, with `gc(); rxode2::rxUnloadAll()`. The disk
  model cache continues to use `qs2` + `digest`, exactly like rxode2/nlmixr2est;
  the in-memory pin cache was removed (same-model reloads come from the `qs2`
  files). Net: ~290 fewer lines, no admixr2-specific memory machinery, and fit
  results are unchanged.

* **`admClearCache()` is removed; use `rxode2::rxClean()`.** admixr2's `qs2`
  caches live in `rxode2::rxTempDir()` alongside rxode2's and nlmixr2est's, so
  `rxode2::rxClean()` -- rxode2's standard cache wipe (unload all models + clear
  the temp dir), which nlmixr2est itself calls to reset -- already clears
  admixr2's cache too. The package-specific `admClearCache()` is therefore
  redundant.

* **`print()` on a fit no longer writes into rmarkdown's namespace.**
  `print.admFit` temporarily overwrote `rmarkdown:::print.paged_df` via
  `assignInNamespace()` (restoring it `on.exit`) to steer nlmixr2est away from
  its paged-table branch. That branch is in fact unreachable: nlmixr2est decides
  between paged and console output by *probing behaviour* -- it prints a
  `paged_df`-classed frame into `capture.output()` and infers "a paged renderer
  consumed my output" from zero captured lines -- but `rmarkdown:::print.paged_df`
  returns its `knit_asis` object visibly and no `print.knit_asis` method exists,
  so the probe always collects output, always returns `FALSE`, and the console
  branch is always taken. The stub therefore changed nothing except skipping the
  discarded probe render (~20 ms per `print(fit)`), at the cost of mutating a
  foreign namespace -- fragile, unsafe under concurrent rendering, and a
  CRAN-policy grey area. Printed output is unchanged, byte for byte. (#58)

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

* Fix an infinite recursion ("evaluation nested too deeply" / "node stack
  overflow") that aborted the first fit of an R session when a covariance matrix
  was requested (`covMethod = "r"`). Accessing `ui$simulationModel` left a
  self-referential compiled-model object in `ui$meta`, which nlmixr2's ui-cloning
  during fit assembly could not traverse. admixr2 now clears that transient
  artifact in `.admLoadModel()`, keeping the ui in the canonical state nlmixr2
  expects. Affected all four estimators (`adfo`/`admc`/`adgh`/`adirmc`) (#81).
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

## Dependencies

* Declare minimum versions for the imported `rxode2 (>= 5.1.2)` and
  `nlmixr2est (>= 6.0.1)`, and for the suggested `nlmixr2 (>= 5.0.0)` (used in
  examples and tests).

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
