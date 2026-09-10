# Handoff: weighting a model source by its own reported covariance, not by `n`

Branch `experiment/model-source-cweight`, forked from `feature/covariate-quadrature`
at `bc2c721` (the model-source sandwich SE machinery -- `.admSrcMeat`,
`.admSrcJac`, `.admSandwich`'s skip/extra split -- is already committed there;
this branch does not depend on anything uncommitted in that PR).

## The defect

A model source's **point estimate** influence is set by `nll_cov_cpp(..., n)`
-- the same n-scaled Wishart term real aggregate data gets, where `n` is a
number the analyst types when generating the study. For real data `n` and
precision are the same thing by construction; for a model source they are not
-- `C_src` (what the source actually reported) never enters the objective at
all, only the post-fit sandwich SE. Measured three ways in
`algorithm/covariate-shift/`: `implied_n.R` (a source's own real, correctly-
paired `n` and `C_src` disagree by up to 564x across its own parameters),
`ms_coverage.R`'s Test W (the point estimate moves with `n` while every real
input -- `E`, `V`, `C_src` -- is held fixed), and `toy_weighting_fix.R`
(36%/99% coverage from the same source's two parameters, closed form).

## The fix

For a model-source study, replace `nll_cov_cpp(E_src, V_src, E_pred(theta),
V_pred(theta), n)` with `r' Mpinv r`, where:

* `r = tau_pred(theta) - tau_src`, `tau = (E, vech V)` stacked -- same
  comparison OLD makes (predicted profile vs the source's reported one), no
  new predicted quantity.
* `M = D C_src D'`, `D` = the SOURCE's own Jacobian (`.admSrcJac`, already
  shipped for the SE) -- how the source's OWN predicted profile would move if
  its OWN parameter estimates had come out slightly differently. Built ONCE,
  from the source model alone; does not depend on the analysis model's
  candidate theta (unlike `.admSrcMeat`'s `S = sum G D`, which needs `G` to
  project into the analysis model's parameter space -- this term never leaves
  moment space, so it needs no such projection).
* `Mpinv` = a pseudo-inverse of `M`, restricted to the directions `C_src`
  actually informs. `M` is rank <= dim(theta_src) regardless of how many
  moments are stacked (the rank-2-over-90-timepoints result); a naive inverse
  does not exist. This is the moment-space analogue of `Vt/n` in
  `.admSandwich`, which exists for the identical reason.

No `n`. No parameter-name matching -- the comparison lives in predicted-
concentration space, which is meaningful for two structurally different
models the way no shared parameter name is (see "what was wrong first",
below).

## What's implemented (real code, not a monkey-patch)

* `.admSrcWeight(studies, idx)` (R/adfweight.R) -- builds `Mpinv` for one
  source group. Row-stacks `D` across every member study (a banded source
  stacks J blocks; this is what keeps a banded source's cross-band
  correlations, the same "one contribution, not J" property `.admSandwich`
  already enforces for the SE).
* `.admAttachSrcWeight(studies)` (R/adfweight.R) -- calls the above once per
  source group, attaches `Mpinv` to the group's first member's `.adm_src`.
  Refuses (stops) rather than silently falling back if any group's weight
  can't be built.
* `adghControl(srcWeight = c("n", "cov"))` -- new, LAST-positioned control
  argument. `"n"` (default) is today's shipped behaviour, byte-identical
  (verified below). `"cov"` opts into the fix.
* `nlmixr2Est.adgh` -- calls `.admAttachSrcWeight` once, after studies are
  finalised, when `srcWeight = "cov"`.
* `.adghNLL` -- a source group carrying `Mpinv` leaves the per-study
  n-weighted loop entirely (ALL its member studies, found via the cheap
  `.admSrcGroups`, not just the one `Mpinv` is attached to) and is scored
  once, after the loop, via `r' Mpinv r`.

## What's verified

* `verify_srcweight.R` -- `.admSrcWeight()` reproduces the eigenvalues
  (`0.495, 0.392, 0.0848, 2.75e-6`, rank 4) already validated by the
  monkey-patched prototype in `algorithm_test_fullM.R`.
* `test_real_wiring.R` -- through the real `nlmixr2()`/`adghControl()` API,
  no monkey-patching:
  * `srcWeight = "cov"` reproduces the monkey-patched prototype's point
    estimate exactly (`tcl=1.6088 tv=3.4012 q=2.0793 v2=4.0950`).
  * n-invariance holds under the real wiring: identical across
    n = 50 / 1500 / 10000.
  * `srcWeight = "n"` (default) reproduces OLD's result exactly
    (`tcl=1.9458 tv=3.3872 q=0.1275 v2=1.9975`) -- **zero regression** to
    existing behaviour.
* `test-adfweight.R` and `test-integration-model-source.R` (the existing,
  shipped test suite) both pass unchanged against these edits.

## What was wrong first, and why (read before "fixing" this a third way)

1. `algorithm_test_srcfix.R` -- named-parameter penalty (`S = identity`,
   directly compare shared `ini()` names). Works ONLY when source and
   analysis model share the exact same structural form (the wt-covariate
   case). Validated there: n-invariance, closed-form GLS match
   (`toy_weighting_fix.R`).
2. `algorithm_test_srcfix_1c2c.R` -- the SAME named-parameter approach,
   applied to a 1-cmt source informing a 2-cmt model. Looked clean (Q/V2
   protected, n-invariant) but was WRONG: pooling `tcl`/`tv` by name assumes
   the misspecified 1-cmt source's own "CL"/"V" are the same estimand as the
   2-cmt model's -- they are not, because model misspecification contaminates
   ALL of a source's own parameter estimates, not just the ones without a
   structural counterpart. Q/V2 being protected was trivially true of any
   version that excludes them; it was not evidence the mechanism was right.
3. `algorithm_test_momentM.R` / `algorithm_test_fullM.R` -- the corrected,
   moment-space version above. Validated in the same 1cmt/2cmt scenario:
   n-invariant, Q/V2 correctly bounded (not by exclusion -- by the source's
   actual rank), full E+V (via the REAL `.admSrcJac`) gives rank 4 matching
   the source's real parameter count and lands closer to truth than
   mean-only.
4. `algorithm_test_diagvsfull.R` -- a diagonal-only `C_src` (all a paper's SE
   column gives you) changes the POINT ESTIMATE now, not just the reported
   SE, since `C_src` feeds the objective. Measured, real, modest in this
   scenario (data arm dominated); direction/magnitude not established --
   would need a replicate study, not a single draw.

## Update: srcWeight = "cov" was a silent no-op under the DEFAULT gradient mode

Found reproducing `vignettes/covariates.Rmd` faithfully (three published
models -- Sato/Ito/Khan -- each `admStudy()`-built with its own `cov`, sex
banded, renal effect marginalised, `est = "adgh"`): `srcWeight = "n"` and
`"cov"` gave BIT-IDENTICAL results to 10 decimal places, and -- the tell --
BOTH moved identically when one source's declared `n` was deliberately
mis-set 20x, which is impossible if `"cov"` genuinely never reads `n`.

Root cause: `.adghGradNLL` (see below) drives the optimizer under the
DEFAULT `grad = "analytical"`, not `.adghNLL` -- so `srcWeight = "cov"`
never reached the function it patches unless the caller also set
`grad = "fd"` explicitly, which nothing enforced or even documented. Every
earlier test in this file happened to set `grad = "fd"`, which is why none
of them caught this.

**Fixed**: `adghControl()` now forces `grad = "fd"` (with a message) whenever
`srcWeight = "cov"` is requested under the default analytical gradient --
same pattern already used for a beta() endpoint's own grad override. See
`vignette_covariates_oldvsnew.R`, the script that found this, kept as the
regression case.

## The 1cmt/2cmt mismatch, inside the real 3-source vignette

`vignette_2cmt_mismatch.R`: same three-source setup, but the true kinetics
are 2-compartment and one analyst (moderate/Ito) fit their trial with a
1-compartment model -- genuinely misspecified, not just a different
covariate choice; normal/Sato and mild/Khan fit (correctly) 2-compartment
models matching the pooled model. Natural scale:

           CL      V1      Q      V2    bcrcl
  truth   5.00   30.00   8.00   60.00   0.60
  old     6.54   31.68   5.17   40.45   0.38
  new     5.11   30.09   8.34   65.64   0.62

`new` recovers every parameter within a few percent. `old`'s damage is not
contained to Q/V2 (35%/33% understated, matching the standalone test) -- it
leaks into `bcrcl`, the renal effect the whole vignette exists to recover
(36% understated), because the fit is joint: a bad pull on Q/V2 distorts
what `bcrcl` has to be to explain the rest of the data.

## What's NOT done -- in priority order

1. **FIXED.** `.adghGradNLL` -- a structurally SEPARATE, hand-coded copy of
   the objective -- drove BOTH the analytical-gradient point estimate (worked
   around earlier by forcing `grad = "fd"`) AND, independently, the Hessian
   for `covMethod = "r"`/`"r,s"` (the `grad = "fd"` forcing did not reach
   this: `use_grad_cov` at the covariance call site was computed separately
   and stayed on the gradient-FD path). Two fixes, both in `.adghCalcCov`/
   `.admSandwichCov`:
   * `.adghCalcCov`'s `use_grad` guard now also forces NLL-FD (Mpinv-aware)
     whenever any study carries an `Mpinv` -- the same test `.adghNLL`
     itself uses to route a source group, so this cannot drift out of sync
     with the objective the way the `grad = "fd"` forcing alone did.
   * `.admSandwichCov` gained `.admSrcMeatCov()`/`.admTauVecDeriv()`, an
     Mpinv-aware sandwich-meat term (`4 A' Mpinv A`, from the score of
     `r' Mpinv r`) for cov-route groups, alongside the unchanged
     `.admSrcMeat()` n-route path -- mixed n/cov fits work with no new
     plumbing.
   Verified on a correctly-specified 1-cmt fit (`verify_se_fix.R`):
   `J ~= 2H` (diagonal ratios 0.999-1.055) and `covMethod = "r"` vs `"r,s"`
   SEs agree to ~1%. `adghControl(covMethod = "none")` is no longer required
   for `srcWeight = "cov"`.
   * Remaining cost, not yet measured or surfaced: `srcWeight = "cov"` still
     forces `grad = "fd"` for the OPTIMIZER (unchanged from before this fix)
     -- so a fit using it loses the analytical gradient entirely. If
     `srcWeight = "cov"` is to fully replace the n route, `.adghGradNLL`
     itself needs the group term (`2 A' Mpinv r`, the same `A` this fix
     already derives) so the analytical path stops being FD-only. Not
     started.
2. **Stratification/covariates + the CORRECTED moment-space fix, together:**
   never tested. The only stratified scenario tested (wt-covariate) used the
   now-invalidated named-parameter version. `.admSrcWeight`'s row-stacking
   across a group's member studies should handle it -- theory says so,
   nobody has run it.
3. **`M` is frozen, not live.** Built once from `theta_src` (fixed, since `D`
   is a property of the source model alone). Measured sensitivity: ~36%
   relative Frobenius-norm change in `M` for a realistic parameter shift
   between `theta_src` and the final combined estimate. Didn't hurt the
   tested fit (data arm dominated); would matter more where the model source
   dominates. A proper fix would iterate: fit, rebuild `M` at the resulting
   theta-hat, refit, check convergence (IRLS). Not implemented.
4. **`adfo`/`admc` untouched entirely** -- only `adgh`.
5. **Multiple data + multiple model sources together, and interaction with
   the existing profiled-lambda_s multi-study weighting**: untested.
6. **Replicate/bias check for the full E+V fix**: only single-draw tests so
   far (`algorithm_test_fullM.R`, `algorithm_test_diagvsfull.R`); the
   mean-only version (`algorithm_test_momentM.R`) got a 12-replicate
   bias/SD check, full E+V hasn't had the equivalent.
7. **Diagonal-vs-full `C_src` direction/magnitude**: one measured example,
   not a replicate study -- can't yet say whether ignoring correlation is a
   systematic bias or just noise in that one draw.

## Literature grounding

* Välitalo, *J Pharmacokinet Pharmacodyn* 48:623-638 (2021), "Pharmacometric
  estimation methods for aggregate data, including data simulated from other
  pharmacometric models" -- almost certainly the direct methodological basis
  for admixr2's own aggregate-data / model-source design (`.adm*` naming
  matches "aggregate data modelling" closely); has a 2025 follow-up
  describing an R package.
* Wang, Kim, Quinney, Zhou & Li, *BMC Systems Biology* 4(S1):S8 (2010),
  "Non-compartment model to compartment model pharmacokinetics transformation
  meta-analysis" -- closest REAL precedent for the structurally-mismatched
  case tested here (10 published midazolam studies, NCA-to-compartmental).
* Gisleskog, Karlsson & Beal, *J Pharmacokinet Pharmacodyn* 29:473-505
  (2002) -- NONMEM's `$PRIOR`, the closest production-software mechanical
  precedent; assumes matching parameterisation, which is exactly the gap the
  moment-space construction here covers.
* Gourieroux, Monfort & Renault, *J Applied Econometrics* 8:S85-S118 (1993),
  "Indirect Inference" -- the binding-function identification condition
  (auxiliary model needs >= as many parameters as the structural model)
  explains the rank-2-cannot-identify-4-parameters result exactly.
* Neuenschwander, Capkun-Niggli, Branson & Spiegelhalter, *Clinical Trials*
  7(1):5-18 (2010) -- meta-analytic-predictive priors, the Bayesian-borrowing
  analogue of the same weight-by-reported-precision principle.
