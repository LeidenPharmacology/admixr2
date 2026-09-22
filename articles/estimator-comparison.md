# Estimator comparison: adfo, admc, adgh and adirmc

## Four estimators, one interface

All four take the same model and study specification and return the same
`admFit` object. They share one likelihood and differ only in how they
approximate the predicted population mean and covariance.

## The shared integration problem

For each study the observed data are a mean vector $`\bar{y}`$ and a
covariance matrix $`S`$ computed from $`n`$ subjects. Under a
multivariate normal approximation the aggregate log-likelihood is

``` math
-2\ell = n \left[
  \log |V_\text{pred}| +
  \operatorname{tr}(V_\text{pred}^{-1} S) +
  r^\top V_\text{pred}^{-1} r
\right], \qquad r = \bar{y} - \mu_\text{pred}
```

where $`\mu_\text{pred}`$ and $`V_\text{pred}`$ are the model-predicted
population mean and covariance of the observations. For a nonlinear
model $`f(\theta,\eta)`$ these are integrals over the random-effects
distribution $`\eta \sim \mathcal{N}(0, \Omega)`$:

``` math
\mu_\text{pred} = \mathbb{E}_\eta\bigl[\mathbb{E}[y \mid \eta]\bigr], \qquad
V_\text{pred} = \operatorname{Var}_\eta\bigl(\mathbb{E}[y \mid \eta]\bigr)
             + \mathbb{E}_\eta\bigl[\operatorname{Var}(y \mid \eta)\bigr]
```

that is, the **law of total variance**. Writing $`V_\text{pred} =
\operatorname{Var}_\eta[f] + \Sigma`$ is correct only when the residual
variance does not depend on the prediction — i.e. for purely
**additive** error. With `prop()`, `pow()` or `lnorm()` the residual
variance is a function of $`f`$, so
$`\mathbb{E}_\eta[\operatorname{Var}(y\mid\eta)] \neq
\operatorname{Var}(y \mid \eta = \bar\eta)`$, and `lnorm()` additionally
scales the conditional *mean* by $`e^{s/2}`$, which scales the whole
covariance rather than just its diagonal. Evaluating the residual at the
population mean instead — the NONMEM “no $`\eta`$–$`\varepsilon`$
interaction” convention — understates
$`\operatorname{diag}(V_\text{pred})`$ by
$`b^2\operatorname{Var}_\eta(f)`$ for proportional error, and admixr2
does **not** do this.

For nonlinear $`f`$ these integrals have no closed form, and evaluating
them is the only thing the four estimators disagree about.

## First-Order (FO): Taylor expansion at $`\eta = 0`$

`adfo` sidesteps the integral entirely by approximating
$`f(\theta, \eta)`$ with its first-order Taylor expansion around
$`\eta = 0`$:

``` math
f(\theta, \eta) \approx f(\theta, 0) + J\,\eta, \qquad
J_{t,j} = \left.\frac{\partial f_t}{\partial \eta_j}\right|_{\eta=0}
```

Substituting into the population moments gives closed-form expressions:

``` math
\mu_\text{pred} = f(\theta, 0), \qquad
V_\text{pred} = J\,\Omega\,J^\top + \Sigma
```

$`J`$ comes from one rxSolve, via sensitivity equations that augment the
ODE system with $`\partial f / \partial \eta_j`$ — no finite
differences, no extra solve per eta. That makes `adfo` the fastest and
fully deterministic: one rxSolve per NLL evaluation.

**When it holds:** exactly, when $`f`$ is linear in $`\eta`$. Otherwise
the true population covariance exceeds $`J\Omega J^\top`$, so FO
underestimates $`\Omega`$, and the bias grows with the nonlinearity.

**AIC:** `adfo` evaluates a *linearised* likelihood, on a different
scale from the other three. Never compare its objective across
estimators.

## Monte Carlo (MC): sample average over $`\eta`$

`admc` estimates the population integrals directly by drawing $`N`$
samples $`\eta_i \sim \mathcal{N}(0, \Omega)`$ and computing sample
moments:

``` math
\mu_\text{pred} = \frac{1}{N} \sum_{i=1}^N f(\theta, \eta_i), \qquad
V_\text{pred} = \frac{1}{N} \sum_{i=1}^N
  \bigl(f(\theta,\eta_i) - \mu_\text{pred}\bigr)
  \bigl(f(\theta,\eta_i) - \mu_\text{pred}\bigr)^\top + \Sigma
```

Nothing is approximated about $`f`$, so the estimator is asymptotically
exact as $`N \to \infty`$. Sobol sequences (the default) typically need
2–5× fewer samples than plain normal draws for the same precision.

The cost is that **every NLL evaluation takes $`N`$ rxSolve calls**, one
per sample. The analytical gradient (sensitivity equations for
mu-referenced parameters, an added first-order direction set for the
rest) adds none, but the base cost per optimizer step stays $`N \times`$
one rxSolve.

Asymptotically exact, with AIC comparable across `admc` fits, and
reliable on standard 1–2 compartment models at moderate IIV. Against
that: slower than `adfo`, and at low `n_sim` the MC noise in the
gradient can set the optimiser oscillating.

## Iterative Reweighting MC (IRMC): proposals fixed, inner loop free

`adirmc` attacks `admc`’s bottleneck — $`N`$ fresh rxSolve calls at
every optimizer step — by **decoupling proposal generation from
optimisation**.

Each outer phase draws $`N`$ proposals $`\eta_i`$ once and holds them
fixed for the whole inner optimisation, so evaluating the NLL is
**matrix operations only**. Reweighting the proposals by their
likelihood under the current $`\Omega`$ keeps the inner objective valid
as the parameters move:

``` math
w_i \propto \frac{p(\eta_i \mid \Omega)}{q(\eta_i)}, \qquad
\mu_\text{pred} = \sum_i w_i\, f(\theta, \eta_i), \qquad
V_\text{pred} = \sum_i w_i\,
  \bigl(f(\theta,\eta_i) - \mu_\text{pred}\bigr)
  \bigl(\cdots\bigr)^\top + \Sigma
```

The inner optimizer therefore runs to convergence for almost nothing;
proposals refresh only between phases, whose box constraints tighten
progressively to guide global convergence.

That suits **complex ODE models** where each solve is expensive: the
total number of solves scales with phases rather than optimizer steps.
It also tolerates poor starting values, the inner loop being
deterministic rather than re-sampled, and its AIC is comparable to
`admc`’s. The cost is per-phase overhead, which does not pay off on a
simple, well-initialised problem with cheap solves.

## Gauss-Hermite (GH): deterministic quadrature over $`\eta`$

`adgh` evaluates the integrals exactly, up to the quadrature rule, on a
tensor-product Gauss-Hermite grid over
$`\eta \sim \mathcal{N}(0, \Omega)`$. Golub-Welsch gives the nodes
$`\eta_q`$ and weights $`w_q`$; the Cholesky factor maps standard-normal
nodes into the correlated $`\eta`$ space:

``` math
\mu_\text{pred} = \sum_{q=1}^Q w_q\, f(\theta, \eta_q), \qquad
V_\text{pred} = \sum_{q=1}^Q w_q\,
  \bigl(f(\theta,\eta_q) - \mu_\text{pred}\bigr)
  \bigl(\cdots\bigr)^\top + \Sigma
```

With $`m`$ nodes per dimension and $`n_\eta`$ random effects the grid
holds $`Q = m^{n_\eta}`$ points, and the objective is **fully
deterministic** — no MC noise, so a clean gradient and a
well-conditioned Hessian.

It is reproducible without fixing `n_sim`, makes no approximation to
$`f`$ (so unlike FO it is unbiased at any IIV), and its AIC is
comparable to `admc` and `adirmc`. Below about four etas $`Q`$ stays
small — $`5^4 = 625`$ — and `adgh` beats `admc` at equal accuracy. Above
that the grid grows exponentially ($`5^5 = 3125`$, $`5^6 = 15625`$), and
there is no inner-loop saving as in `adirmc`: every NLL evaluation runs
all $`Q`$ solves. For high-dimensional IIV, use `admc` or `adirmc`.

## Common setup

``` r

library(admixr2)
library(rxode2)
library(nlmixr2)
library(ggplot2)

# The examplomycin study and model from `vignette("admixr2")`, supplied by the
# shared helpers in this vignette's setup file. All four fits below differ only
# in the estimator, so the model and the data are held fixed here.
pk_model <- admVignetteModel
study    <- c(admVignetteStats(), list(ev = rxode2::et(amt = 100)))
```

It is mu-referenced (`cl <- exp(tcl + eta.cl)`), which is what gives all
four estimators an analytical gradient from sensitivity equations — see
[Advanced
usage](https://leidenpharmacology.github.io/admixr2/articles/advanced.html#mu-referencing-and-sensitivity-equations).

## Fitting with adfo

The fastest, and the natural starting point for screening or for initial
estimates to hand to MC.

``` r

fit_fo <- nlmixr2(
  pk_model, admData(), est = "adfo",
  control = adfoControl(
    studies  = list(examplomycin = study),
    maxeval  = 500L,
    seed     = 1L
  )
)
```

## Fitting with admc

``` r

fit_mc <- nlmixr2(
  pk_model, admData(), est = "admc",
  control = admControl(
    studies   = list(examplomycin = study),
    n_sim     = 5000L,
    cov_n_sim = 10000L,
    maxeval   = 300L,
    seed      = 1L
  )
)
```

## Fitting with adirmc

Slower per evaluation, but more robust to poor starting values and a
high-dimensional $`\Omega`$. These settings keep the vignette build
fast; raise `n_sim` for real work.

``` r

fit_irmc <- nlmixr2(
  pk_model, admData(), est = "adirmc",
  control = adirmcControl(
    studies         = list(examplomycin = study),
    n_sim           = 2000L,
    phases          = c(2, 1, 0.5, 0.01),
    cov_n_sim       = 10000L,
    omega_expansion = 1.5,
    seed            = 1L
  )
)
```

## Fitting with adgh

A drop-in alternative to `admc` when there are few etas. The five-eta
model here gives $`5^5 = 3125`$ nodes, about the practical limit. Start
at `n_nodes = 3` (243 nodes) and go to 5 or 7 only if IIV is large (SD
\> 0.4).

``` r

fit_gh <- nlmixr2(
  pk_model, admData(), est = "adgh",
  control = adghControl(
    studies  = list(examplomycin = study),
    n_nodes  = 5L,
    maxeval  = 300L,
    seed     = 1L
  )
)
```

## Comparing parameter estimates

Both should land near the truth (CL = 5, V1 = 10, V2 = 30, Q = 10, ka =
1; IIV variance 0.09 throughout; prop.sd = 0.2). examplomycin has
moderate IIV and is near-linear in $`\eta`$, so FO bias is small here.

``` r

get_pars <- function(fit) {
  s  <- fit$env$admExtra$struct
  om <- diag(fit$env$admExtra$omega)
  sg <- sqrt(unlist(fit$env$admExtra$sigma_var))
  c(exp(s), om, sg)
}

tbl <- data.frame(
  Parameter = c(
    paste0("exp(", names(fit_fo$env$admExtra$struct), ")"),
    paste0("var(", fit_fo$env$admExtra$eta_col_names, ")"),
    "prop.sd"
  ),
  True = c(5, 10, 30, 10, 1, rep(0.09, 5), 0.2),
  adfo = round(get_pars(fit_fo), 4),
  admc = round(get_pars(fit_mc), 4)
)

knitr::kable(tbl, caption = "Parameter estimates vs true values")
```

| Parameter   |  True |    adfo |    admc |
|:------------|------:|--------:|--------:|
| exp(tcl)    |  5.00 |  4.9100 |  4.9628 |
| exp(tv1)    | 10.00 |  7.3161 | 10.2582 |
| exp(tv2)    | 30.00 | 32.3637 | 29.8857 |
| exp(tq)     | 10.00 | 10.0853 |  9.7381 |
| exp(tka)    |  1.00 |  0.8043 |  1.0302 |
| var(eta.cl) |  0.09 |  0.0968 |  0.1011 |
| var(eta.v1) |  0.09 |  0.1427 |  0.1042 |
| var(eta.v2) |  0.09 |  0.0905 |  0.0964 |
| var(eta.q)  |  0.09 |  0.0938 |  0.1075 |
| var(eta.ka) |  0.09 |  0.0777 |  0.0967 |
| prop.sd     |  0.20 |  0.1908 |  0.1895 |

Parameter estimates vs true values {.table}

## Comparing objectives

`admc` and `adirmc` evaluate the same likelihood, so their -2LL values
compare directly. `adfo`’s is linearised, on a different scale, and
**must not** be set against them or used for cross-estimator AIC:

``` r

cat(sprintf("adfo  -2LL = %.2f   AIC = %.2f\n", fit_fo$objective, AIC(fit_fo)))
#> adfo  -2LL = -3676.41   AIC = -3654.41
cat(sprintf("admc  -2LL = %.2f   AIC = %.2f\n", fit_mc$objective, AIC(fit_mc)))
#> admc  -2LL = -3690.26   AIC = -3668.26
```

Use AIC only within the same estimator for model selection.

## NLL convergence traces

``` r

plots_fo <- plot(fit_fo, which = "nll")
```

![adfo and admc NLL convergence
traces.](estimator-comparison_files/figure-html/nll-trace-1.png)

adfo and admc NLL convergence traces.

``` r

plots_mc <- plot(fit_mc, which = "nll")
```

![adfo and admc NLL convergence
traces.](estimator-comparison_files/figure-html/nll-trace-2.png)

adfo and admc NLL convergence traces.

An `adirmc` trace shows outer phases rather than individual L-BFGS
steps, so the jumps where box constraints tighten are expected.

## When to use each

| Situation | Recommendation |
|----|----|
| Rapid model screening, many candidate models | `adfo` — fastest per evaluation |
| Initial estimates before MC/GH refinement | `adfo` then hand off to `admc` or `adgh` |
| Weak IIV (CV \< 20 %) and near-linear model | `adfo` estimates reliable for inference |
| Standard 1–2 compartment PK, ≤ 4 etas | `adgh` — noise-free, faster than MC at equivalent accuracy |
| Standard 1–2 compartment PK, ≥ 5 etas or large IIV | `admc` with `grad = "sens"` |
| Initial exploration or poor starting values | `admc` with `n_restarts >= 3` |
| Complex ODE system with expensive solves | `adirmc` — inner loop needs no new rxSolve calls |
| High-dimensional Omega (≥ 5 etas) | `adirmc` — inner loop scales with phases not steps |
| Non-Gaussian or bounded IIV | `adirmc` with `omega_expansion > 1` |
| Optimal design / design evaluation (matched moments) | `adgh` with `datagenControl(method = "gh")` |
| Need exact likelihood for AIC comparison across models | `admc`, `adgh`, or `adirmc` (not `adfo`) |
| Maximum reproducibility, no MC noise | `adfo` or `adgh` (both deterministic) |

Every argument, with its default, is in
[`?adfoControl`](https://leidenpharmacology.github.io/admixr2/reference/adfoControl.md),
[`?admControl`](https://leidenpharmacology.github.io/admixr2/reference/admControl.md),
[`?adghControl`](https://leidenpharmacology.github.io/admixr2/reference/adghControl.md)
and
[`?adirmcControl`](https://leidenpharmacology.github.io/admixr2/reference/adirmcControl.md).

## See also

- [Choosing a residual error
  model](https://leidenpharmacology.github.io/admixr2/articles/error-models.md)
  — which error models each estimator supports
- [Advanced
  usage](https://leidenpharmacology.github.io/admixr2/articles/advanced.md)
  — gradient modes, restarts, and AIC/BIC model comparison
- [Getting
  started](https://leidenpharmacology.github.io/admixr2/articles/admixr2.md)
