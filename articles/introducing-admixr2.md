# The mathematics of aggregate data modelling

## Why aggregate data?

Meta-analysis in pharmacometrics means combining evidence across studies
for estimates no single study could support alone. The classical
approach pools individual patient data (IPD) into one dataset —
scientifically ideal, and rare: data sharing agreements, proprietary
restrictions and the fact that most published work reports only summary
statistics make IPD pooling the exception.

Model-based meta-analysis (MBMA) fits pharmacometric models directly to
aggregate summaries instead. Those summaries come from two sources.
*Directly published*: many studies report a mean concentration-time
profile with the variance at each time point, sometimes the full
covariance matrix. *Model-derived*: a published model and its estimates
(fixed effects, random-effects variance, residual error) generate the
expected mean vector and covariance by simulation, without the original
data. The covariance carries information about the shape of individual
profiles that a variance-only summary discards — two studies can share
every variance and differ entirely in correlation structure, reflecting
different sources of variability. Both types analyse within the same
framework, which makes the full literature — competitor compounds,
historical datasets, regulatory submissions — available for dose
selection, trial design and disease-progression modelling.

The aggregate data likelihood, first systematically described by
[Välitalo (2021)](https://doi.org/10.1007/s10928-021-09760-1), provides
the statistical foundation for MBMA: given a nonlinear mixed-effects
model, what is the probability of observing the reported mean and
covariance? The first R implementation was admr, described in [van de
Beek et al. (2025)](https://doi.org/10.1007/s10928-025-10011-w). admixr2
re-implements the same approach within the nlmixr2/rxode2 ecosystem,
extending it with new algorithmic ideas that are the focus of this post.

------------------------------------------------------------------------

## The aggregate likelihood

### Observed data

For a single study with $`n`$ subjects, each measured at observation
times $`t_1, \ldots, t_T`$, the observed data are two sufficient
statistics:

``` math
\bar{y} \in \mathbb{R}^T, \qquad S \in \mathbb{R}^{T \times T}
```

where $`\bar{y}`$ is the observed mean vector and $`S`$ is the observed
covariance matrix (ML denominator $`n`$, not $`n - 1`$).

### The model

A nonlinear mixed-effects model specifies that subject $`i`$’s
observations arise from

``` math
y_i = f(\theta, \eta_i) + \varepsilon_i, \qquad
\eta_i \sim \mathcal{N}(0, \Omega), \quad
\varepsilon_i \sim \mathcal{N}(0, \Sigma)
```

where $`f(\theta, \cdot)`$ is an ODE system (or closed-form prediction
function) with structural parameters $`\theta`$ and individual random
effects $`\eta_i`$. $`\Omega \in \mathbb{R}^{d \times d}`$ is the
between-subject covariance, $`\Sigma \in \mathbb{R}^{T \times T}`$ is
the residual error covariance (typically a diagonal function of the
structural prediction, e.g. proportional error).

### Population moments

Under this model the population-level distribution of $`y_i`$ has mean
and covariance

``` math
\mu_\text{pred} = \mathbb{E}_\eta\bigl[f(\theta, \eta)\bigr], \qquad
V_\text{pred}   = \operatorname{Var}_\eta\bigl[f(\theta, \eta)\bigr]
                + \mathbb{E}_\eta\bigl[\operatorname{Var}(y \mid \eta)\bigr].
```

This is the law of total variance, and the second term is written out
rather than abbreviated to $`\Sigma`$ for a reason. Only for
**additive** error is the residual variance free of the prediction,
giving $`\mathbb{E}_\eta[\operatorname{Var}(y\mid\eta)] = \Sigma`$. With
`prop()`, `pow()` or `lnorm()` it is a function of $`f`$, and averaging
it over subjects is not evaluating it at the population mean — for
proportional error the two differ by exactly
$`b^2\operatorname{Var}_\eta(f)`$. admixr2 averages. Writing
$`+\,\Sigma`$ throughout is the “no $`\eta`$–$`\varepsilon`$
interaction” convention of individual-level fitting, and it understates
$`\operatorname{diag}(V_\text{pred})`$ here; see [Choosing a residual
error
model](https://leidenpharmacology.github.io/admixr2/articles/error-models.md).
The sections below write $`\Sigma`$ for brevity, with that averaging
understood.

### The aggregate log-likelihood

If the $`n`$ subjects are independent and drawn from the same
population, the sample mean $`\bar{y}`$ follows approximately
$`\mathcal{N}(\mu_\text{pred},\, V_\text{pred}/n)`$, and the sample
covariance $`S`$ concentrates around $`V_\text{pred}`$. Under a
multivariate normal approximation to the joint density of
$`(\bar{y}, S)`$, the log-likelihood reduces to

``` math
-2\ell = n \Bigl[
  \log |V_\text{pred}| +
  \operatorname{tr}(V_\text{pred}^{-1} S) +
  r^\top V_\text{pred}^{-1} r
\Bigr], \qquad
r = \bar{y} - \mu_\text{pred}.
```

This is the central formula in admixr2. Every estimator minimises it;
they differ only in how $`\mu_\text{pred}`$ and $`V_\text{pred}`$ are
computed.

The trace and quadratic terms go through a Cholesky factorisation of
$`V_\text{pred}`$, avoiding explicit inversion and staying stable even
for near-singular covariances. The log-determinant falls out as twice
the sum of the factor’s log-diagonal.

------------------------------------------------------------------------

## The integration problem

For a nonlinear ODE model, $`\mu_\text{pred}`$ and $`V_\text{pred}`$ are
integrals over $`\eta \sim \mathcal{N}(0, \Omega)`$. There is no closed
form. The four estimators differ in how they handle this intractability:
`adfo` linearises the integrand, `admc` samples it randomly, `adgh`
integrates it on a deterministic quadrature grid, and `adirmc` reuses a
fixed sample, iteratively reweighting it as the parameters move.

------------------------------------------------------------------------

## First-Order linearisation (adfo)

The simplest approach approximates $`f(\theta, \eta)`$ by its
first-order Taylor expansion around $`\eta = 0`$:

``` math
f(\theta, \eta) \approx f(\theta, 0) + J\,\eta,
\qquad J_{t,j} = \left.\frac{\partial f_t}{\partial \eta_j}\right|_{\eta=0}.
```

Substituting into the population moment definitions gives closed-form
expressions:

``` math
\mu_\text{pred} = f(\theta, 0), \qquad
V_\text{pred}   = J\,\Omega\,J^\top + \Sigma.
```

The Jacobian $`J`$ is the sensitivity matrix of the ODE output with
respect to the random effects, evaluated at $`\eta = 0`$. In admixr2 it
is obtained by **augmenting the ODE system with sensitivity equations**

``` math
\frac{d}{dt}\frac{\partial x}{\partial \eta_j} =
\frac{\partial g}{\partial x}\,\frac{\partial x}{\partial \eta_j} +
\frac{\partial g}{\partial \eta_j},
```

where $`g`$ is the ODE right-hand side. One solver call delivers $`f`$
and every column of $`J`$ together, with no extra solve per
random-effect dimension. Its predecessor admr computed this Jacobian by
**finite differences** on $`f`$, one extra solve per dimension.

**Gradient of the FO NLL.** Because
$`V_\text{pred} = J\Omega J^\top + \Sigma`$ is explicit in $`\Omega`$
and $`\Sigma`$, the gradient with respect to these parameters is
available analytically via the matrix derivative identities
$`d \log|A| = \operatorname{tr}(A^{-1} dA)`$ and
$`d \operatorname{tr}(A^{-1}B) = -\operatorname{tr}(A^{-1} dA\, A^{-1} B)`$.
Structural parameter gradients need the second derivative
$`dJ/d\theta`$, which admixr2 takes from an order-2 sensitivity model,
so they too are analytical; central finite differences are the fallback
for when that model cannot be built.

**Accuracy.** Exact when $`f`$ is linear in $`\eta`$. For nonlinear
models at large IIV it underestimates
$`\operatorname{Var}_\eta[f(\theta, \eta)]`$, so $`V_\text{pred}`$ comes
out too small and $`\Omega`$ negatively biased. In practice the bias is
negligible below a CV of 20–30 % on each PK parameter, on a model weakly
nonlinear in $`\eta`$.

------------------------------------------------------------------------

## Monte Carlo simulation (admc)

The MC estimator makes no approximation to $`f`$. It draws $`N`$
random-effect samples and computes sample moments:

``` math
\hat{\mu} = \frac{1}{N}\sum_{i=1}^N f(\theta, \eta_i), \qquad
\hat{V}   = \frac{1}{N}\sum_{i=1}^N
  \bigl(f(\theta,\eta_i) - \hat{\mu}\bigr)
  \bigl(f(\theta,\eta_i) - \hat{\mu}\bigr)^\top + \Sigma,
\qquad \eta_i \sim \mathcal{N}(0, \Omega).
```

It converges to the true population moments as $`N \to \infty`$, and its
AIC is directly interpretable.

**Quasi-random sampling.** admixr2 draws from Sobol sequences (also in
admr) rather than pseudo-random normal deviates. Sobol points are
*low-discrepancy*, placed to fill gaps in the sample space rather than
clustering at random. On the smooth integrands of pharmacometric MC
integration that uniform coverage cuts the integration error at a given
$`N`$, for no extra model evaluations.

**Numerical stability.** The sample covariance accumulates centred
products in a single pass through fused C++ kernels, avoiding an
intermediate $`N \times T`$ allocation, and the resulting $`\hat{V}`$ is
positive semi-definite by construction.

------------------------------------------------------------------------

## Gauss-Hermite quadrature (adgh)

MC replaces the population integral with a *random* average; GH replaces
it with a *deterministic* one, a tensor-product Gauss-Hermite rule
evaluating the same expectation over
$`\eta \sim \mathcal{N}(0, \Omega)`$ at a small set of fixed nodes. It
is new in admixr2, with no counterpart in admr.

For a single random effect, the probabilists’ Gauss-Hermite rule
approximates an expectation under the standard normal by a weighted sum
over $`m`$ nodes,

``` math
\mathbb{E}_{z\sim\mathcal{N}(0,1)}[g(z)] \approx \sum_{q=1}^m w_q\, g(x_q),
```

with nodes $`x_q`$ and weights $`w_q`$ (normalised so $`\sum_q w_q = 1`$
and $`\sum_q w_q x_q^2 = 1`$) obtained from the Golub-Welsch algorithm —
the eigendecomposition of the symmetric tridiagonal Jacobi matrix of the
Hermite recurrence, requiring no external tables. For $`n_\eta`$ random
effects the rule is applied as a tensor product, giving a grid of
$`Q = m^{n_\eta}`$ nodes $`x_q \in \mathbb{R}^{n_\eta}`$ with product
weights. The Cholesky factor $`L`$ ($`\Omega = LL^\top`$) maps the
standard-normal nodes into the correlated random-effect space,
$`\eta_q = L x_q`$, so the population moments become

``` math
\mu_\text{pred} = \sum_{q=1}^Q w_q\, f(\theta, L x_q), \qquad
V_\text{pred} = \sum_{q=1}^Q w_q\,
  \bigl(f(\theta, L x_q) - \mu_\text{pred}\bigr)
  \bigl(\cdots\bigr)^\top + \Sigma.
```

Structurally this is the MC estimator with a fixed deterministic node
grid in place of $`N`$ random draws and quadrature weights in place of
uniform ones, so it plugs into exactly the same aggregate MVN
$`-2\ell`$.

**Noise-free objective.** Fixed nodes and weights make the objective a
smooth deterministic function of $`(\theta, \Omega, \Sigma)`$ — no Monte
Carlo noise to average away, no `n_sim` to tune. The analytical gradient
(closed-form contractions through the same ODE sensitivity outputs the
MC estimator uses) is therefore exact, and the numerical Hessian behind
the standard errors is well-conditioned. Like `admc` it approximates
nothing about $`f`$, so it is unbiased at any IIV, on the same
likelihood scale, and its AIC compares directly to `admc` and `adirmc`.

**Cost.** $`Q = m^{n_\eta}`$ grows exponentially in the number of random
effects. Up to ~4 etas it stays small ($`5^4 = 625`$) and `adgh` is
substantially faster than MC at equal accuracy; beyond that the node
count becomes prohibitive and `admc` or `adirmc` are preferable.
`n_nodes = 5` gives near-exact covariance moments to IIV SD ~0.5, and 7
extends that to ~0.7.

------------------------------------------------------------------------

## Iterative Reweighting MC (adirmc)

MC’s central bottleneck is that *every* NLL evaluation takes $`N`$ ODE
solves, one per sample — expensive on a complex model, over the hundreds
of evaluations an L-BFGS optimisation may need. IRMC decouples sample
generation from optimisation.

### Proposal distribution

At each outer phase, $`N`$ proposals are drawn once from an inflated
prior:

``` math
\tilde{\eta}_i \sim \mathcal{N}\!\left(0,\, \alpha\,\Omega_0\right),
\qquad \alpha \geq 1,
```

where $`\Omega_0`$ is the current Omega estimate and $`\alpha`$ the
expansion factor. The predictions $`f(\theta_0, \tilde{\eta}_i)`$ are
computed once, stored, and reused for the rest of the inner optimisation
without further ODE solves.

### Importance weights

A move to a candidate $`(\theta, \Omega)`$ leaves the stored predictions
inexact. For the random-effects covariance, the proposals are reweighted
by the ratio of target to proposal density:

``` math
w_i \propto \frac{p(\tilde{\eta}_i \mid \Omega)}{q(\tilde{\eta}_i \mid \alpha\,\Omega_0)}
= \exp\!\left[
  -\tfrac{1}{2}\tilde{\eta}_i^\top
  (\Omega^{-1} - (\alpha\,\Omega_0)^{-1})
  \tilde{\eta}_i
\right].
```

Weights are normalised by softmax, and the weighted mean and covariance
of the predictions are the IRMC estimates of $`\mu_\text{pred}`$ and
$`V_\text{pred}`$:

``` math
\mu_\text{pred} = \sum_i w_i\, f(\theta_0, \tilde\eta_i), \qquad
V_\text{pred} = \sum_i w_i\,
  (f(\theta_0, \tilde\eta_i) - \mu_\text{pred})
  (\cdots)^\top + \Sigma.
```

With fixed proposals the inner objective is a smooth deterministic
function of $`(\theta, \Omega)`$, optimisable to high precision.
Proposals refresh between phases, whose box constraints tighten
progressively to guide global convergence, so the number of ODE solves
scales with phases rather than inner steps.

------------------------------------------------------------------------

## Gradient computation: from finite differences to sensitivity equations

The gradient is the single most important ingredient in efficient
nonlinear optimisation, and it is where admr and admixr2 differ most
consequentially.

### Gradient in admr

admr stores a fixed Sobol base matrix `biseq` for the lifetime of an
optimisation run. Random-effect samples are formed as
$`\eta_i = \texttt{biseq}
\cdot L(\Omega)`$, so the NLL is a smooth *deterministic* function of
the parameter vector for fixed `biseq`. The gradient for the MC
estimator is computed by forward finite differences of this
deterministic NLL:

``` math
\frac{\partial \ell}{\partial \theta_k} \approx
\frac{\ell(\theta + h\,e_k) - \ell(\theta)}{h}, \qquad h = 10^{-6}.
```

On a model with $`p`$ parameters each gradient call therefore costs
$`p`$ extra NLL evaluations, $`p \times N`$ additional ODE solves. That
being expensive, admr defaults to gradient-free BOBYQA, with L-BFGS and
an FD gradient available but off (`use_grad = FALSE`).

Under IRMC, proposals are fixed inside a closure before the inner
optimisation, so each inner NLL evaluation is matrix operations only.
admr’s optional FD inner gradient costs $`p`$ extra importance-weighted
recomputations per step and no ODE solves — but the inner optimizer
still defaults to BOBYQA, and the FD approximation carries truncation
error of order $`O(h)`$.

### CRN analytical gradient in admixr2

admixr2 fixes the sample set $`\{\eta_i\}`$ and computes the gradient of
the frozen deterministic MC NLL analytically — an application of the
Common Random Numbers (CRN) identity.

For a mu-referenced parameterisation
$`\psi_k = \exp(\theta_k + \eta_k)`$, the chain rule gives

``` math
\frac{\partial \hat{\mu}}{\partial \theta_k}
= \frac{1}{N}\sum_{i=1}^N
  \frac{\partial f(\theta, \eta_i)}{\partial \eta_k} \cdot
  \underbrace{\frac{\partial \psi_k}{\partial \theta_k}}_{\psi_k}
  \cdot \underbrace{\frac{\partial \eta_k}{\partial \psi_k}}_{1/\psi_k}
= \frac{1}{N}\sum_{i=1}^N
  \frac{\partial f(\theta, \eta_i)}{\partial \eta_k},
```

i.e. the sample average of the ODE sensitivity output
$`\partial f / \partial
\eta_k`$ — which is already available from the NLL evaluation via the
augmented ODE system. No additional ODE solves are needed.

The gradient of $`\hat{V}`$ with respect to $`\Omega`$ follows
similarly: the Cholesky factor $`L`$ ($`\Omega = LL^\top`$) enters
through the sample covariance, and differentiating through the MVN chain
rule gives closed-form expressions in the stored sensitivity outputs.

So the **full MC gradient costs the same as a single NLL evaluation** —
no extra ODE solves, and exact rather than a finite-difference
approximation. On a model with $`p = 10`$ parameters and $`N = 5000`$
samples, admr’s FD gradient needs $`50\,000`$ additional solves per
optimizer step; admixr2’s CRN gradient needs none.

### Analytical IRMC inner gradient

The importance-weighted inner objective has a tractable analytical
gradient with respect to both $`\Omega`$ and $`\theta`$. The gradient
passes through the softmax normalisation via the identity

``` math
\frac{\partial}{\partial p} \sum_i w_i(\,p\,)\, g_i
= \sum_i w_i \left(\frac{\partial g_i}{\partial p} +
  \frac{\partial \log w_i}{\partial p}\,(g_i - \bar{g})\right),
```

where $`\bar{g} = \sum_i w_i g_i`$. Differentiating the MVN
log-likelihood of the weighted mean and covariance through this identity
gives closed-form expressions in the pre-computed ODE solutions and
sensitivity outputs. Fixed proposals mean no inner gradient evaluation
needs an ODE solve under any implementation, so the gain over admr’s FD
approach is exactness and a $`p\times`$ reduction in inner-level matrix
operations.

The inner optimizer therefore runs L-BFGS on an exact gradient,
converging in far fewer iterations than BOBYQA.

### Kappa correction for non-mu-referenced parameters

Both packages handle structural parameters $`\beta_k`$ entering $`f`$
with no paired random effect
($`\mathrm{CL} = \exp(\beta_1)\cdot\exp(\eta_1)`$, but $`V_\text{max}`$
has no $`\eta`$). Changing one during the inner IRMC optimisation leaves
the stored predictions $`f(\theta_0, \tilde\eta_i)`$ inexact even in the
mean, so both apply a kappa correction:

``` math
\kappa(\beta) = f(\beta_\text{new}, 0) - f(\beta_\text{orig}, 0),
```

shifting $`\mu_\text{pred}`$ by the predicted change at $`\eta = 0`$.

admixr2 adds a **linearized kappa** (`kappa_method = "linearized"`) that
computes $`J_\kappa = \partial f / \partial \beta_\text{single}`$ once
per outer iteration from one batched FD solve, then takes the correction
as $`J_\kappa (\beta_\text{new} - \beta_\text{orig})`$ through the inner
loop — matrix work only, no solver calls. That matters most on complex
ODE systems: exact kappa costs one extra rxSolve per inner NLL step for
the single-beta parameters, while linearized kappa leaves the inner loop
entirely free of them.

------------------------------------------------------------------------

## Parameter space geometry

### Cholesky parameterisation of $`\Omega`$

$`\Omega`$ must be positive definite, and unconstrained optimisation
over symmetric matrices can step outside that.

Both packages parameterise it through a Cholesky factor $`L`$
($`\Omega = LL^\top`$), and differ in how the entries are encoded. admr
uses a log transform on the diagonal and a **covlogit** off it, mapping
each $`L_{jk}`$’s corresponding correlation through a logit after
normalisation. admixr2 uses:

``` math
p_k = \log(\Omega_{kk}) = 2\log(L_{kk})
\quad\text{(diagonal)}, \qquad
p_{jk} = L_{jk}
\quad\text{(off-diagonal, raw)}.
```

The raw form simplifies the gradient substantially —
$`\partial \Omega / \partial
p_{jk}`$ is a rank-1 matrix with no logit derivative to chain through —
and avoids the instability of logit transforms as correlations approach
±1.

Encoding $`p = \log(\Omega_{kk})`$ rather than $`\log(L_{kk})`$
equalises gradient sensitivity across parameter types: a unit optimizer
step changes $`\Omega_{kk}`$ by $`e^2`$, matching structural parameters
on the log scale.

### Parameter preconditioning

The unconstrained parameter vector is pre-scaled by a diagonal $`C`$
before reaching L-BFGS, each entry estimating the curvature in that
direction: unity for exponential transforms, a derivative-based
magnitude for logit/probit, and $`\max(|L_{jk}|, 0.1)`$ off the Cholesky
diagonal. This lowers the condition number and speeds convergence
wherever structural parameters and variance components differ by orders
of magnitude.

------------------------------------------------------------------------

## Summary of improvements over admr

admr (van de Beek et al., 2025) established the core aggregate data
workflow in R. admixr2 keeps its three estimators (FO, MC, IRMC) and the
same aggregate MVN likelihood, adds a fourth deterministic Gauss-Hermite
estimator, and replaces the surrounding architecture and several
algorithmic components.

| Aspect | admr | admixr2 |
|----|----|----|
| Ecosystem | Standalone R package | nlmixr2/rxode2 integration |
| Model syntax | Custom `genopts()` + prediction function | nlmixr2 [`ini()`](https://nlmixr2.github.io/rxode2/reference/ini.html) / [`model()`](https://nlmixr2.github.io/rxode2/reference/model.html) blocks |
| Fit object | Plain list | nlmixr2 fit object (AIC, logLik, plot, print) |
| Estimators | FO, MC, IRMC | FO, MC, GH, IRMC |
| Gauss-Hermite quadrature estimator | Not available | New deterministic, noise-free estimator |
| FO Jacobian ($`\partial f/\partial\eta`$) | C++ finite differences | ODE sensitivity equations |
| MC gradient | FD of frozen MC NLL ($`p \times N`$ extra solves) | CRN analytical (0 extra ODE solves) |
| IRMC inner optimizer | BOBYQA (default) or FD gradient | L-BFGS with analytical gradient |
| Linearized kappa | Not available | Optional; eliminates ODE calls in inner loop |
| $`\Omega`$ off-diagonal encoding | covlogit (correlation logit) | Raw Cholesky entry $`L_{jk}`$ |
| Parameter preconditioning | No | Yes (diagonal scaling) |
| Parallel restarts | Sequential chains | `mirai` workers |

The gradient improvements are the most consequential. Going from
$`p \times N`$ extra ODE solves per step (admr FD) to none (admixr2 CRN)
cuts wall-clock fitting time by one to two orders of magnitude for
gradient-based algorithms on complex models. The analytical IRMC inner
gradient compounds it: the inner optimisation converges in fewer
iterations, each iteration is faster, and linearized kappa lets the
whole inner loop run without touching the ODE solver.

------------------------------------------------------------------------

## Further reading

The papers linked at the top give the mathematical foundations in full.
[Estimator
comparison](https://leidenpharmacology.github.io/admixr2/articles/estimator-comparison.md)
works all four estimators through the included `examplomycin` dataset,
and [Advanced
usage](https://leidenpharmacology.github.io/admixr2/articles/advanced.md)
covers gradient modes, multi-restart fitting and model comparison via
AIC.
