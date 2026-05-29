# The mathematics of aggregate data modelling

## Why aggregate data?

Pharmacometric analyses traditionally begin with individual patient
records. Each observation row carries a subject identifier, a time
stamp, and a measured drug concentration. The statistical model is
fitted to this raw dataset using methods such as FOCE or SAEM.

In many practically important situations, individual-level records are
simply not available. Clinical pharmacology databases, published
meta-analyses, and regulatory submissions frequently report only
*summary statistics*: a mean concentration profile and a between-subject
variance structure. Reconstructing individual data from these summaries
is at best approximate, and often impossible.

The aggregate data likelihood, first systematically described by
[Välitalo et al. (2021)](https://doi.org/10.1007/s10928-021-09760-1),
bypasses individual records entirely. It works directly with the mean
vector and covariance matrix that a study reports and asks: how likely
are these summaries under a given nonlinear mixed-effects model?

The first R implementation was admr, described in [van de Beek et
al. (2025)](https://doi.org/10.1007/s10928-025-10011-w). admixr2
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
V_\text{pred}   = \operatorname{Var}_\eta\bigl[f(\theta, \eta)\bigr] + \Sigma.
```

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

This is the central formula in admixr2. Each estimator minimises the
same objective; they differ only in how $`\mu_\text{pred}`$ and
$`V_\text{pred}`$ are computed.

The trace and quadratic terms are evaluated via Cholesky factorisation
of $`V_\text{pred}`$, which avoids explicit matrix inversion and is
numerically stable even for near-singular covariances. The
log-determinant falls out as twice the sum of log-diagonal elements of
the Cholesky factor.

------------------------------------------------------------------------

## The integration problem

For a nonlinear ODE model, $`\mu_\text{pred}`$ and $`V_\text{pred}`$ are
integrals over $`\eta \sim \mathcal{N}(0, \Omega)`$. There is no closed
form. The three estimators differ in how they handle this
intractability.

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

where $`g`$ is the ODE right-hand side. A single solver call delivers
both $`f`$ and all columns of $`J`$ simultaneously — no additional ODE
solves per random-effect dimension. In the predecessor admr, this
Jacobian was computed by **finite differences** on $`f`$, requiring one
extra ODE solve per random-effect dimension.

**Gradient of the FO NLL.** Because
$`V_\text{pred} = J\Omega J^\top + \Sigma`$ is explicit in $`\Omega`$
and $`\Sigma`$, the gradient with respect to these parameters is
available analytically via the matrix derivative identities
$`d \log|A| = \operatorname{tr}(A^{-1} dA)`$ and
$`d \operatorname{tr}(A^{-1}B) = -\operatorname{tr}(A^{-1} dA\, A^{-1} B)`$.
Structural parameter gradients use forward finite differences through
the full FO NLL, reusing the sensitivity solve.

**Accuracy.** The FO approximation is exact when $`f`$ is linear in
$`\eta`$. For nonlinear models with large IIV the approximation
underestimates $`\operatorname{Var}_\eta[f(\theta, \eta)]`$:
$`V_\text{pred}`$ is too small and the estimated $`\Omega`$ is
negatively biased. Practically, bias is negligible when the coefficient
of variation for each PK parameter is below 20–30 % and the model is
weakly nonlinear in $`\eta`$.

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

The estimator converges to the true population moments as
$`N \to \infty`$; its AIC is directly interpretable.

**Quasi-random sampling.** admixr2 draws samples using Sobol sequences
(also available in admr) rather than pseudo-random normal deviates.
Sobol sequences are *low-discrepancy*: successive points are placed to
fill gaps in the sample space, rather than clustering randomly. For the
smooth integrands that arise in pharmacometric MC integration, this
uniform coverage reduces the integration error for a given $`N`$
compared to i.i.d. draws, at the cost of no additional model
evaluations.

**Numerical stability.** The sample covariance is computed in fused C++
kernels that accumulate centred products in a single pass, avoiding an
intermediate $`N \times T`$ allocation. The resulting $`\hat{V}`$ is
guaranteed positive semi-definite by construction.

------------------------------------------------------------------------

## Iterative Reweighting MC (adirmc)

The central bottleneck of the MC estimator is that *every* NLL
evaluation requires $`N`$ ODE solves — one per sample. For a complex
model each solve is expensive, and an L-BFGS optimisation may require
hundreds of NLL evaluations.

IRMC decouples sample generation from optimisation.

### Proposal distribution

At each outer phase, $`N`$ proposals are drawn once from an inflated
prior:

``` math
\tilde{\eta}_i \sim \mathcal{N}\!\left(0,\, \alpha\,\Omega_0\right),
\qquad \alpha \geq 1,
```

where $`\Omega_0`$ is the current Omega estimate and $`\alpha`$ is the
expansion factor. The predictions $`f(\theta_0, \tilde{\eta}_i)`$ are
computed once and stored. For the remainder of the inner optimisation
they are reused without further ODE solves.

### Importance weights

When the parameter vector moves to a candidate $`(\theta, \Omega)`$, the
stored predictions are no longer exact. For the random-effects
covariance $`\Omega`$, the proposals are reweighted by the ratio of
target density to proposal density:

``` math
w_i \propto \frac{p(\tilde{\eta}_i \mid \Omega)}{q(\tilde{\eta}_i \mid \alpha\,\Omega_0)}
= \exp\!\left[
  -\tfrac{1}{2}\tilde{\eta}_i^\top
  (\Omega^{-1} - (\alpha\,\Omega_0)^{-1})
  \tilde{\eta}_i
\right].
```

Weights are normalised via softmax. The weighted mean and covariance of
the predictions then provide the IRMC estimates of $`\mu_\text{pred}`$
and $`V_\text{pred}`$:

``` math
\mu_\text{pred} = \sum_i w_i\, f(\theta_0, \tilde\eta_i), \qquad
V_\text{pred} = \sum_i w_i\,
  (f(\theta_0, \tilde\eta_i) - \mu_\text{pred})
  (\cdots)^\top + \Sigma.
```

Given fixed proposals this inner objective is a smooth, deterministic
function of $`(\theta, \Omega)`$ and can be optimised to high precision.
Between phases, proposals are refreshed and box constraints are
progressively tightened to guide global convergence. The number of ODE
solves scales with the number of phases, not the number of inner
optimizer steps.

------------------------------------------------------------------------

## Gradient computation: from finite differences to sensitivity equations

The gradient is the single most important ingredient for efficient
nonlinear optimisation. The improvements in gradient computation
represent the most consequential algorithmic differences between admr
and admixr2.

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

For a model with $`p`$ parameters, each gradient call therefore costs
$`p`$ extra NLL evaluations, i.e. $`p \times N`$ additional ODE solves.
Because this is expensive, admr defaults to gradient-free BOBYQA; L-BFGS
with FD gradient is available but off by default (`use_grad = FALSE`).

For the IRMC estimator, proposals are fixed inside a closure before the
inner optimisation begins, so each inner NLL evaluation requires only
matrix operations — no ODE solves. admr’s optional FD inner gradient
therefore incurs $`p`$ extra importance-weighted recomputations per
step, but still no extra ODE solves. Nonetheless, the inner optimizer
defaults to BOBYQA (no gradient) and the FD approximation introduces
truncation error of order $`O(h)`$.

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
rule yields closed-form expressions in the stored sensitivity outputs.

The net result is that the **full MC gradient costs the same as a single
NLL evaluation** — it requires no additional ODE solves and is exact
rather than a finite-difference approximation. For a model with
$`p = 10`$ parameters and $`N = 5000`$ samples, admr’s FD gradient
requires $`50\,000`$ additional ODE solves per optimizer step; admixr2’s
CRN gradient requires none.

### Analytical IRMC inner gradient

The importance-weighted inner objective has a tractable analytical
gradient with respect to both $`\Omega`$ and $`\theta`$. The gradient
passes through the softmax normalisation via the identity

``` math
\frac{\partial}{\partial p} \sum_i w_i(\,p\,)\, g_i
= \sum_i w_i \left(\frac{\partial g_i}{\partial p} +
  \frac{\partial \log w_i}{\partial p}\,(g_i - \bar{g})\right),
```

where $`\bar{g} = \sum_i w_i g_i`$. Differentiating the multivariate
normal log-likelihood of the weighted mean and covariance through this
identity yields closed-form expressions in the pre-computed ODE
solutions and sensitivity outputs. Since proposals are fixed, the inner
gradient evaluation requires no ODE solves regardless of the
implementation — the admixr2 advantage over admr’s FD approach is
exactness and a $`p\times`$ reduction in inner-level matrix operations.

The inner optimizer therefore runs L-BFGS with an exact gradient and
converges in far fewer iterations than BOBYQA.

### Kappa correction for non-mu-referenced parameters

Both admr and admixr2 handle models where some structural parameters
$`\beta_k`$ enter $`f`$ without a paired random effect
(e.g. $`\mathrm{CL} = \exp(\beta_1)
\cdot \exp(\eta_1)`$, but $`V_\text{max}`$ has no $`\eta`$). When such
parameters change during the inner IRMC optimisation, the stored
predictions $`f(\theta_0,
\tilde\eta_i)`$ are no longer exact even for the mean. Both packages
apply a kappa correction:

``` math
\kappa(\beta) = f(\beta_\text{new}, 0) - f(\beta_\text{orig}, 0),
```

shifting $`\mu_\text{pred}`$ by the predicted change at $`\eta = 0`$.

admixr2 introduces a **linearized kappa** option
(`kappa_method = "linearized"`) that computes the Jacobian
$`J_\kappa = \partial f / \partial \beta_\text{single}`$ once per outer
iteration (via a single batched FD solve) and approximates the
correction as $`J_\kappa (\beta_\text{new} - \beta_\text{orig})`$ during
the inner loop — a pure matrix operation requiring no ODE solver calls.
This is particularly valuable for complex ODE systems where each kappa
evaluation is otherwise expensive: with the exact kappa, every inner NLL
step requires one extra rxSolve for the single-beta parameters; with
linearized kappa, the inner loop is completely free of ODE solver calls.

------------------------------------------------------------------------

## Parameter space geometry

### Cholesky parameterisation of $`\Omega`$

The between-subject covariance $`\Omega`$ must be positive definite.
Unconstrained optimisation in the space of symmetric matrices can
produce non-positive-definite steps.

Both admr and admixr2 parameterise $`\Omega`$ through a Cholesky factor
$`L`$ ($`\Omega = LL^\top`$). The packages differ in how the Cholesky
entries are encoded. admr uses a log transform for diagonal entries and
a **covlogit** transform for off-diagonal entries (the correlation
corresponding to each $`L_{jk}`$ is mapped through a logit after
normalisation). admixr2 uses:

``` math
p_k = \log(\Omega_{kk}) = 2\log(L_{kk})
\quad\text{(diagonal)}, \qquad
p_{jk} = L_{jk}
\quad\text{(off-diagonal, raw)}.
```

The raw parameterisation simplifies the gradient computation
substantially: the Jacobian $`\partial \Omega / \partial p_{jk}`$ is a
simple rank-1 matrix with no logit derivative to chain through. It also
avoids the numerical instability of logit-based transforms when
correlations approach ±1.

The specific encoding $`p = \log(\Omega_{kk})`$ rather than
$`p = \log(L_{kk})`$ equalises gradient sensitivity across parameter
types: a unit optimizer step changes $`\Omega_{kk}`$ by a factor of
$`e^2`$, matching the sensitivity of structural parameters on the log
scale.

### Parameter preconditioning

The unconstrained parameter vector is pre-scaled by a diagonal matrix
$`C`$ before being passed to L-BFGS. Each diagonal entry estimates the
curvature in that direction: unity for exponential transforms, a
derivative-based magnitude for logit/probit transforms, and
$`\max(|L_{jk}|, 0.1)`$ for off-diagonal Cholesky entries.
Preconditioning reduces condition number and accelerates convergence
when structural parameters and variance components differ by orders of
magnitude.

------------------------------------------------------------------------

## Summary of improvements over admr

admr (van de Beek et al., 2025) established the core aggregate data
workflow in R. admixr2 preserves the same three estimators and the same
aggregate MVN likelihood while replacing the surrounding architecture
and several algorithmic components.

| Aspect | admr | admixr2 |
|----|----|----|
| Ecosystem | Standalone R package | nlmixr2/rxode2 integration |
| Model syntax | Custom `genopts()` + prediction function | nlmixr2 [`ini()`](https://nlmixr2.github.io/rxode2/reference/ini.html) / [`model()`](https://nlmixr2.github.io/rxode2/reference/model.html) blocks |
| Fit object | Plain list | nlmixr2 fit object (AIC, logLik, plot, print) |
| Estimators | FO, MC, IRMC | FO, MC, IRMC |
| FO Jacobian ($`\partial f/\partial\eta`$) | C++ finite differences | ODE sensitivity equations |
| MC gradient | FD of frozen MC NLL ($`p \times N`$ extra solves) | CRN analytical (0 extra ODE solves) |
| IRMC inner optimizer | BOBYQA (default) or FD gradient | L-BFGS with analytical gradient |
| Linearized kappa | Not available | Optional; eliminates ODE calls in inner loop |
| $`\Omega`$ off-diagonal encoding | covlogit (correlation logit) | Raw Cholesky entry $`L_{jk}`$ |
| Parameter preconditioning | No | Yes (diagonal scaling) |
| Parallel restarts | Sequential chains | `furrr` workers |

The gradient improvements are the most consequential. Moving from a
gradient that requires $`p \times N`$ extra ODE solves per step (admr
FD) to one that requires none (admixr2 CRN) reduces wall-clock fitting
time by one to two orders of magnitude for gradient-based algorithms on
complex models. The analytical IRMC inner gradient compounds this
saving: the inner optimisation converges in fewer iterations, each
iteration is faster, and the linearized kappa option allows the entire
inner loop to run without any ODE solver calls.

------------------------------------------------------------------------

## Further reading

The mathematical foundations are described in detail in the papers
linked at the top of this post. The [Estimator
comparison](https://leidenpharmacology.github.io/admixr2/estimator-comparison.md)
vignette shows worked examples on the included `examplomycin` dataset
with all three estimators. The [Advanced
usage](https://leidenpharmacology.github.io/admixr2/advanced.md)
vignette covers gradient modes, multi-restart fitting, and model
comparison via AIC.
