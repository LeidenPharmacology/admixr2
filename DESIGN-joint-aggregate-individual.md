# Joint estimation from aggregate summaries and individual records

A design note, not an implementation. It says what the combined objective is, what
in admixr2 already supports it, what breaks, and the order I would build it in.

The motivating case is ordinary: you have your own trial's individual records and
several published summaries for the same drug. Today you must choose — fit the
individual data in nlmixr2 and ignore the literature, or digitise the literature
and ignore your own subject-level detail.

## 1. The objective

Both likelihoods are already log-densities of the same population model
`(theta, Omega, sigma)`, so they simply add:

```
-2LL_joint(theta, Omega, sigma) = sum_s -2LL_agg(study s) + sum_i -2LL_ind(subject i)
```

That is not an approximation or a weighting scheme. Independent data sets
contribute additively to a log-likelihood, and the two terms describe disjoint
subjects. **No tuning constant belongs here**, and one should be refused if
proposed: `lambda * LL_agg + LL_ind` is not a likelihood for any `lambda != 1`, and
its maximiser has no interpretation. If the aggregate studies appear to dominate,
that is information content, not a bug — see §4.

The two terms differ in what they integrate over:

- **Individual**: `-2 log integral p(y_i | eta_i) p(eta_i) d(eta_i)` per subject,
  which is FOCEI/Laplace/SAEM territory — nlmixr2est's business, and admixr2 has
  no reason to reimplement it.
- **Aggregate**: `n_s (log|V_pred| + tr(V_pred^-1 V_obs) + r' V_pred^-1 r)` per
  study, computed from the marginal moments — admixr2's business.

## 2. What already fits

More than I expected when I started writing this.

- **The parameterisation is shared.** `.admParseIniDf()` builds the optimizer
  vector from `ui$iniDf`, the same `iniDf` nlmixr2est parameterises from. Both
  sides already speak `(struct thetas on the log scale, Omega as log-Cholesky,
  sigma by role)`.
- **The residual layer is shared and is the harder half.** `R/errmodel.R` derives
  its specs from `ui$predDf`, which is nlmixr2's own object, and it already knows
  every error model nlmixr2 supports. Crucially it computes `E[Var(y|eta)]` —
  the aggregate quantity — from the same spec that gives the conditional density.
- **Covariates now work on both sides.** The per-subject path treats the covariate
  as data, exactly as an individual-level fit does. A joint fit where the
  individual records carry real weights and the published studies carry weight
  DISTRIBUTIONS is the natural case, and both are already expressible.
- **The gradient composes.** `d(-2LL_joint) = d(-2LL_agg) + d(-2LL_ind)`, and
  admixr2's analytic gradients are already in the optimizer's coordinates.

## 3. What breaks, in the order it will bite

**(a) Two optimizers, two conventions.** nlmixr2est drives its own optimisation
(`foceiFit`), and admixr2 drives nloptr. A joint fit needs ONE optimizer calling
both objectives. The cheapest honest route is to keep admixr2's nloptr loop and
call nlmixr2est for the individual term as a *function evaluation*, which means
finding an entry point that returns an objective and gradient at a given parameter
vector without running a fit. If none exists, the fallback is to evaluate the
individual likelihood in admixr2 directly for the subset of error models already
implemented — duplication, but bounded and testable.

**(b) Eta handling is genuinely different.** The individual term needs per-subject
EBEs (an inner optimisation, or a Laplace expansion, per subject per outer step).
The aggregate term needs no EBEs at all — it integrates eta out to form moments.
These are not reconcilable into one machinery, and they should not be: the joint
objective is a SUM, so each term computes its own way and only the parameter vector
is shared. Attempting to share eta machinery is the design error to avoid.

**(c) Residual sharing is a modelling decision that must be surfaced.** If your
assay and a published study's assay differ, a shared `sigma` is wrong. The natural
handling is nlmixr2's existing per-endpoint residuals: give the aggregate studies
their own endpoint. This must be a user choice, not a default, and the default
should probably be SEPARATE, since assuming a shared assay across a literature
meta-analysis is the stronger claim.

**(d) The covariance report spans both.** `.admScaleReportedCov()` and
`.admOmegaJacobian()` already rotate an optimizer-scale Hessian to the reported
scale; a joint Hessian is the sum of the two blocks and rotates identically. This
part should be nearly free.

**(e) Objective-value reporting.** `fit$objective` currently holds the aggregate
-2LL. A joint fit must report the sum, and the two components separately, or no
one can tell which data are driving the answer.

## 4. The thing to measure first, before building anything

**How much does an aggregate study weigh against individual records?**

An aggregate study contributes `n_s` times a per-subject term; a set of individual
records contributes one term per subject. Nominally comparable — but the aggregate
term scores `V_obs` as Wishart-from-normal, and the simulation-estimation study
here found its standard errors already run 0.84–0.95 of the empirical spread
(90–94% coverage against nominal 95%). If the aggregate term is mildly
over-confident on its own, then in a joint fit it will pull harder than its true
information content warrants, and the individual data — the data you trust most —
will be under-weighted.

That is measurable before a line of joint code is written: fit the same simulated
population three ways (individual only, aggregate only, joint) and compare each
parameter's empirical spread against what each fit reports. If the joint fit's
spread sits closer to the aggregate-only fit than subject counts predict, the
over-confidence is real and needs addressing — most likely by correcting the
aggregate term's effective sample size rather than by weighting it.

**Do this experiment first.** It decides whether joint estimation is worth building
in this form, and it needs no new machinery: simulate individuals, fit them with
nlmixr2est, reduce the same individuals to summaries, fit those with admixr2, and
compare.

## 5. Build order

1. The §4 experiment. It can invalidate the rest.
2. A joint objective behind one nloptr call, aggregate term from admixr2,
   individual term from nlmixr2est, no gradient (derivative-free) — enough to
   prove the objective is right on a simulated case where both data sets come from
   one known truth.
3. Analytic gradient by summation, checked against a finite difference of the
   joint objective.
4. Joint covariance, reusing `.admScaleReportedCov()`.
5. Per-endpoint residual separation for aggregate versus individual data.

## 6. What not to do

- Do not introduce a weight on either term. See §1.
- Do not convert individual data to summaries and fit everything through the
  aggregate likelihood. That discards exactly the information the individual
  records were included for, and the aggregate likelihood's Wishart assumption is
  a worse description of a small hand-made summary than of a published one.
- Do not convert aggregate data to pseudo-individuals. Simulating subjects
  consistent with a reported mean and covariance and then fitting them as though
  they were observed manufactures information that was never measured; the
  resulting standard errors describe a study that does not exist.
- Do not share eta machinery between the two terms (§3b).
