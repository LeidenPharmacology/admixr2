# Dimension collapse for the aggregate covariate integral

Handoff for the covariates paper. Implemented on `feature/covariate-quadrature`
(PR #121), commits `b944feb` → `b23b847`.

## The problem

The aggregate likelihood needs population moments, so every prediction is
marginalised over the random effects **and** over the covariate distribution
in the population. Done as a product quadrature that costs `n^(n_eta + p)`
nodes for `n_eta` random effects and `p` covariates — and every node is a row
in the ODE solve. Three covariates and two random effects at 7 nodes is 16,807
solve rows per study per objective evaluation.

## The observation

That product integrates in `n_eta + p` dimensions, but the integrand rarely
lives there. Write the latent vector as `ξ = (η_std, z)`, standard normal in
both blocks, with `z` the Gaussian-copula scores of the covariates. For the
multiplicative and allometric forms that dominate PK/PD,

```
log θ = a + B'ξ        B's η rows scaled by L = chol(Ω)
```

so the parameters — and therefore the prediction — depend on `ξ` **only
through** `w = U'ξ`, where `U` spans the column space of `B`. Since `ξ` is
Gaussian, `w` is Gaussian, and the integral is exactly `r = rank(B)`
dimensional. The design becomes an `r`-dimensional Gauss–Hermite rule with
preimage `ξ = Uw`.

**This is not an approximation.** It is a change of variables; the reduced
design integrates the same function against the same law.

The bound that matters for the paper: `rank(B) ≤` the number of **parameters**
the latents reach. A two-parameter model never needs more than a
two-dimensional design, however many random effects and covariates it carries.
Cost stops scaling with the number of covariates and starts scaling with the
number of parameters.

## Evidence

Two random effects, three covariates on the random-effect-carrying parameter,
three studies. Reference is a `15² × 15³ = 759,375`-row product grid.

| design | rows/study | decomposition | −2LL error vs reference |
|---|---|---|---|
| full product quadrature | 8,575 | `5² × 7³` | not run |
| shipping route | 525 | `5² × 21` | 1.1e-02 … 5.0e-02 |
| joint collapse | **324** | `18²` | **2.1e-09 … 4.2e-08** |

Every count decomposes. `n_nodes = 5` (random effects), `cov_nodes = 7`
(covariates), both package defaults. `21 = ⌈7·3/1⌉` is the covariate-only
collapse's per-direction count at rank 1; `18 = ⌈7·5/2⌉` is the joint one's at
rank 2 over five latent dimensions.

**The middle row is not a product quadrature** and must not be labelled one —
the shipping route already collapses the covariate block, so it is the random-
effect grid crossed with a rank-1 covariate design. The genuine product is the
top row and was not run; the reference for the error column is a `15² × 15³ =
759,375`-row grid.

Fewer rows *and* six orders of magnitude more accurate, at four parameter
points. Structural moments agree with the reference to 1e-12…3e-09 across six
points. The analytic gradient matches central differences to 9e-10…2e-08.

For the covariate block alone (no random effects), a rank-1 collapse of three
covariates uses 21 design points against a 343-point product grid, and is more
accurate at every point tested.

## What the construction requires

1. **A certificate.** The map must be affine in `ξ`, or a *single index* — some
   function of one linear combination, which is enough because a Gauss–Hermite
   rule integrates any composition with a normal index exactly to degree
   `2n−1`. Affine is the identity-link special case. Both are detected
   numerically from a probe; anything else is refused.
2. **Verification of the design, not the certificate.** The reduced design must
   reproduce the law of every covariate-reading assignment — checked against a
   large probe on the first two moments and the reciprocal (a PK parameter
   enters the prediction as both `f` and `1/f`). Proxies for this were tried
   and are not sharp enough.
3. **Node counts must scale up, not across.** A collapsed direction carries the
   *combined* spread of the axes it absorbed, so it needs `~p/r` times the
   nodes of a single axis. Giving it the per-axis count makes the collapse
   **worse than the grid it replaces** — measured at 2.3e-02 against the grid's
   4.5e-05.

## The subtlety worth a paragraph in the paper

**The rotation is not a property of the data.** `B` depends on the covariate
coefficients, which are estimated, and on `Ω`, which is estimated. So `U` moves
as the optimiser moves, and the design must be re-aimed at every objective
evaluation. A design fixed at the starting values pins the covariates to the
wrong line in latent space and misses the orthogonal variation entirely —
measured at 53 to 163 −2LL units for a 0.1 move in one coefficient, while every
*non*-covariate parameter stayed exact to 1e-5.

This is worth stating because it is invisible to the obvious checks: moment and
gradient tests evaluate at a single parameter point, and at that point a stale
design is correct by construction. It took a full simulation–estimation run to
surface.

What is fixed at admission is structural only — the rank, the node count, and
which route each column's loading is read through. None may be re-derived per
call: rank or node count moving would change the *number* of design points
mid-fit and step the objective discontinuously; a route is a threshold decision
a borderline column could flip.

Re-aiming costs no ODE solves (it is arithmetic on the model's parameter
assignments) — 0.8 ms per study per evaluation, 6% of the total.

## Limitations to state honestly

- Discrete covariates are strata, not directions: they have no latent normal to
  rotate into and are enumerated exactly instead.
- A covariate-by-random-effect interaction makes the loading depend on `η`,
  which the construction detects and refuses.
- **Quadrature construction, and it does not transfer to the Monte-Carlo
  estimator at all.** The same change of variables lets it sample in the
  reduced subspace, and that is exactly right, but it was tried and removed for
  two reasons. The accuracy benefit is not measurable. With randomised QMC (Cranley–Patterson, 24 replicates, paired):
  the covariate collapse buys 1.12×–1.47× on the covariance and the joint one
  adds 1.01×–1.17×, with the interquartile range spanning 1.0 in every cell.
  Do not claim a speed or accuracy figure for the sampler.

  *An earlier draft of this handoff claimed ~2×. That came from medians taken
  across randtoolbox `seed` values, and `seed` is inert — the sequences are
  bit-identical — so it was one realisation with no error bar. `scrambling` in
  the same package is disabled outright. Both are worth knowing before quoting
  any QMC comparison.*

  And it breaks **common random numbers**. The rotation depends on the
  covariate coefficients, which are estimated, so re-aiming it — which
  correctness requires — makes the *draws* move with the parameters: a step of
  1e-6 in one coefficient shifted the sampled covariate values by 1.1e-4, so a
  CRN finite difference in that direction carries a design change on top of the
  parameter change. The sampler is deterministic in the covariate distribution
  alone, which is data, and the gradient depends on that. Worth a sentence in
  the paper: the reduction is a property of *quadrature*, not of the change of
  variables, because a design can be re-aimed between evaluations and a common
  random number cannot.

  Dropping a weak direction from the sampler is separately **not** viable: it
  biases and floors (8.3e-03 where keeping it converges to 3.0e-04).

## Open

- *(closed)* The Monte-Carlo estimator's batched Hessian paths were listed here
  as needing per-block rotations. They do not: the sampler no longer collapses
  at all, so every path takes the same deterministic product draw and there is
  nothing to serve per block. The deterministic estimator's Hessian goes through
  the ordinary objective and is re-aimed there.


## For Appendix B: differentiating through a parameter-dependent design

A reviewer will ask this, because the quadrature nodes depend on the parameters
being differentiated. The argument is short.

The computed objective is `F(ψ) = Q(U(ψ), ψ)`, where `Q` is the `r`-dimensional
rule and `U(ψ)` an orthonormal basis of `col(B(ψ))`. The implementation
differentiates only the explicit `ψ`, holding the design fixed, so the question
is whether the omitted `(∂Q/∂U)(∂U/∂ψ)` matters.

It does not, and the reason is an invariance. `θ` depends on `ξ` only through
`u = B'ξ`, and with `ξ ~ N(0, I)` we have `u ~ N(0, B'B)`. The design maps the
product Gauss–Hermite grid through `M = (U'B)'`, and

```
M M' = B'U U'B = B'B         since col(U) = col(B), so U U'B = B
```

for **every** orthonormal basis `U` of `col(B)`. So each admissible `U` yields a
different matrix square root of the same covariance, and the rules they generate
agree on everything inside the exactness class of the product rule. They differ
only outside it — that is, by terms of the order of the quadrature error itself.
Hence `∂Q/∂U = O(ε_quad)`, and the omitted gradient term is `O(ε_quad)` too:
the same order as the objective's own error, not a separate approximation. Here
`ε_quad ≈ 2e-09`, and the analytic gradient matches central differences of the
re-aimed objective to 9e-10…2e-08 — consistent.

**State the scope precisely**, because the invariance is about *re-choosing a
basis within a fixed span*, not about using a wrong span. A design whose `U`
spans `col(B(ψ₀))` while the parameters have moved to `ψ` is not covered by it:
the loading acquires a component orthogonal to the design's subspace, and the
error is `O(1)`, not `O(ε_quad)`. That is the failure mode that made re-aiming
necessary in the first place, and the two statements should appear together or
the first will read as licence for the second.

One caveat worth not overstating: a natural guess is that the stale-design error
should be second order in `‖ψ − ψ₀‖`, since a loading component of size `ε`
contributes `O(ε²)` in variance. The measured errors (5e-3, 7.8e-3, 2.3e-2
relative at coefficient displacements of 0.1, 0.3, 0.6) are **sub**-quadratic,
so do not assert an order for it — report the invariance argument and the
finite-difference check, which are what actually hold.
