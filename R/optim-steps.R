# Optimizer plumbing: nloptr algorithm selection, Shi (2021) finite-difference
# step selection, and the box-constraint warning.
#
# Split out of utils.R (see R/covreport.R for why that is
# safe). The two step accessors .admGH()/.admGH0() look collapsible and are not:
# indexing a per-parameter measured-step vector by an ETA number would recycle
# across the n_sim rows of an eta perturbation and return a plausible wrong
# gradient, which is exactly what .admGH0() exists to prevent.


# -- nloptr algorithm selection ------------------------------------------------

# Valid nloptr algorithm names, queried from the installed nloptr so the set
# always matches the user's version (no hardcoded list to go stale). Returns
# character(0) if the query fails (unexpected nloptr internals) -- callers then
# defer validation to nloptr itself at fit time.
.admNloptrAlgorithms <- function() {
  algs <- tryCatch({
    o  <- nloptr::nloptr.get.default.options()
    pv <- o[o$name == "algorithm", "possible_values"]
    trimws(strsplit(as.character(pv), ",")[[1]])
  }, error = function(e) character(0))
  algs[nzchar(algs)]
}

# TRUE if the algorithm consumes a user-supplied gradient (the _LD_ / _GD_
# NLopt families); FALSE for the derivative-free _LN_ / _GN_ families.
.admAlgoNeedsGrad <- function(algorithm) grepl("_(LD|GD)_", algorithm)

# Default nloptr algorithm for a gradient mode: BOBYQA when gradless, LBFGS
# otherwise. The default pairing is LBFGS (gradient) <-> BOBYQA (gradless).
.admDefaultAlgorithm <- function(grad)
  if (grad == "none") "NLOPT_LN_BOBYQA" else "NLOPT_LD_LBFGS"

# Index a finite-difference step that may be a SCALAR -- the fixed fallback, one
# constant shared by every parameter -- or a per-parameter VECTOR of Shi21's
# measured steps.
#
# Every FD site in the package reads its step through this or .admGH0(), which is
# what lets measured steps reach the optimizer's gradient without a single
# signature change: the driver hands the gradient a vector where it used to hand
# it a number. A scalar takes the identity branch, which is the path taken
# whenever the measurement could not be made.
.admGH <- function(h, idx) if (length(h) == 1L) h else h[idx]

# The step for a difference that has NO parameter index to key on.
#
# Not every finite difference in admixr2 perturbs a parameter of the objective.
# .admGrad also differences in ETA space -- `eta_hi[, j] <- eta_hi[, j] + h` over
# n_sim rows -- to build d(pred)/d(eta), which the analytic kernels then contract
# into the omega and paired-theta gradient. There is no parameter index there, so
# a per-parameter vector must NOT be indexed with the eta number: it would pick
# an unrelated parameter's step, and used bare it would RECYCLE across the n_sim
# rows without a warning (n_sim %% n_par is virtually always 0), handing every
# draw a different perturbation and returning a plausible, wrong gradient.
#
# So a measured-step vector carries the scalar it was built from as attr(, "h0"),
# and index-less sites read that. The eta step therefore stays the fixed one,
# which is also the honest answer: the probe measured the OBJECTIVE, and
# d(pred)/d(eta) is a different function.
.admGH0 <- function(h)
  if (length(h) == 1L) h else {
    h0 <- attr(h, "h0", exact = TRUE)
    if (is.null(h0)) h[[1L]] else h0
  }

# Per-parameter FD steps for the covariance HESSIAN.
#
# `pmax(abs(p[idx]), 0.1) * cov_h_outer` -- one constant applied to every
# parameter -- is the FALLBACK now, not the behaviour. It is a guess about how
# much noise the objective carries, and it is the guess behind the "Hessian not
# positive definite ... try increasing cov_h_outer" advice: too fine a step and
# the difference is noise, too coarse and it is curvature the second-order term
# does not capture. Shi21 measures the objective instead of guessing, per
# parameter (see .admShi21Central).
#
# Single-sourced because this expression used to be written out in
# .adfoCalcCov(), .adghCalcCov() and .admCalcCov(), byte-identical but for the
# label -- and the fixed-step form appeared TWICE inside each copy, so the
# heuristic was restated six times across three files. The fallback exists so a
# failed probe lands exactly where the fixed step would have.
#
# NOT .admShi21Steps(). That returns the optimum for a FIRST central derivative,
# `h ~ (3 eps_f/|f'''|)^(1/3)`, and the Hessian takes a SECOND difference, whose
# error is `(h^2/12)|f''''| + 4 eps_f/h^2` and whose optimum therefore scales as
# `eps_f^(1/4)` -- around ten times larger at a machine-precision objective. Too
# fine a step in a second difference amplifies noise as `4 eps_f/h^2`, which is
# precisely what tips a marginal Hessian out of positive-definiteness. Measured
# against an exact second derivative, the first-derivative step was 6x to 385x
# worse across noise levels from 1e-15 to 1e-7.
#
# So what is reused here is the MEASUREMENT -- `eps_f` from ECnoise, the thing
# `cov_h_outer` could only guess -- fed to the second-difference rule rather than
# the first-difference one.
# `cov_h_outer` SCALES the measured step; it does not merely back it up. Making
# it a pure fallback was tried and is wrong: the measurement almost always
# succeeds, so the argument would be inert in practice -- and it is the escape
# hatch the docs point users at ("Hessian not positive definite ... try
# increasing cov_h_outer"). An argument that looks like it does something and
# does not is the exact failure mode `nlmixr2Gill83`'s dropped tuning arguments
# are an example of. So the step is the measured one at the DEFAULT
# `cov_h_outer`, and moves proportionally when the user changes it.
.ADM_COV_H_REF <- .Machine$double.eps^(1/5)

.admHessSteps <- function(fn, p, idx, cov_h_outer, .var.name = "CalcCov") {
  fixed <- pmax(abs(p[idx]), 0.1) * cov_h_outer
  if (length(idx) == 0L) return(fixed)
  eps_f <- tryCatch(.admEcNoise(fn, p, idx[[1L]]), error = function(e) NA_real_)
  f0    <- tryCatch(abs(fn(p)), error = function(e) NA_real_)
  if (!is.finite(eps_f) || eps_f <= 0 || !is.finite(f0)) {
    warning(sprintf(
      "%s: could not estimate the objective's noise level -- falling back to the fixed step scale.",
      .var.name), call. = FALSE)
    return(fixed)
  }
  # adgh deliberately ships a finer cov_h_outer than the sampling estimators
  # (its quadrature objective is noise-free); expressing the user's value as a
  # ratio to one reference keeps that relative intent instead of erasing it.
  pmax(abs(p[idx]), 0.1) * (eps_f / max(f0, 1))^(1/4) *
    (cov_h_outer / .ADM_COV_H_REF)
}

# Per-parameter FD steps for the OPTIMIZER's gradient, measured once at the point
# the fit starts from.
#
# The covariance Hessian can afford to call .admShi21Steps() directly: it runs
# once, post-fit, at a known parameter vector. A gradient cannot -- the optimizer
# calls it thousands of times. So measure once and REUSE, which is what FOCEI
# does too (its gill83 runs at the first gradient evaluation, `nF == 1`, and
# every later difference uses the stored per-parameter step). Here the
# measurement happens in the driver, before the optimizer is handed anything, and
# travels as the `grad_h` argument that was already there.
#
# The measured step is ABSOLUTE at `p`. The estimators express theirs
# differently -- adfo/adgh scale by `pmax(abs(p), 0.1)`, admc uses the raw
# number -- so it is divided by whatever that site will multiply it back by
# (`scaled`). The step is therefore exactly the measured one at `p`, and tracks
# the parameter afterwards under the convention that site already had, rather
# than freezing an absolute number that stops making sense once the optimizer
# has moved a decade.
#
# `idx` is the set of parameters that will actually be finite-differenced;
# everything else keeps the constant and costs nothing. Passing `integer(0)` --
# a fit whose gradient is fully analytic -- skips the probe entirely.
.admShi21GradH <- function(fn, p, idx, grad_h, scaled = TRUE, .var.name = "grad") {
  if (length(idx) == 0L) return(grad_h)
  out <- rep_len(as.numeric(grad_h), length(p))
  base <- if (scaled) pmax(abs(p[idx]), 0.1) else rep(1, length(idx))
  hf <- .admShi21Steps(fn, p, idx, fallback = base * grad_h, .var.name = .var.name)
  out[idx] <- hf / base
  # The scalar every index-less difference keeps using -- see .admGH0().
  attr(out, "h0") <- as.numeric(grad_h)[[1L]]
  out
}

# Warn when the fit finished ON the gradient-mode box constraint.
#
# A gradient fit is run inside `p0 +/- grad_bounds` on the optimizer scale, a
# constraint the user did not write: .admBuildOptVec() returns -Inf/Inf for
# struct thetas and omega unless the model declares explicit bounds. nloptr
# reports normal convergence at a box corner, and a finite estimate and a finite
# SE are printed, so a parameter pinned 5 optimizer units from its starting value
# is indistinguishable from a converged one. On the log scale that is a factor of
# exp(5) ~ 148: fit `tv <- log(20)` to data whose true V is 5000 and V is clamped
# at ~2968, silently.
#
# admc/adgh/adirmc have always run with a gradient by default and so have always
# had this; adfo acquired it in 0.4.1 when its default gradient mode changed.
# Reporting it is the cheap half -- the fix is the user's (widen grad_bounds, or
# start closer), and it is only actionable if they are told.
#
# `p` is the final optimizer-scale solution, `p0` the start. Only entries whose
# model-declared bound is infinite are reported: a user-written bound reached is
# the user's own constraint, not this one.
.admWarnOnBounds <- function(p, p0, ov, grad_bounds, pinfo) {
  if (is.null(p) || is.null(p0) || !is.finite(grad_bounds) || grad_bounds <= 0)
    return(invisible(character(0)))
  n <- min(length(p), length(p0))
  if (n == 0L) return(invisible(character(0)))
  p  <- p[seq_len(n)]; p0 <- p0[seq_len(n)]
  lo <- if (is.null(ov$lower)) rep(-Inf, n) else ov$lower[seq_len(n)]
  hi <- if (is.null(ov$upper)) rep(Inf,  n) else ov$upper[seq_len(n)]
  # Reconstruct the box nloptr was actually given, and ask whether the solution
  # sits on it. Within 0.1% of the half-width counts as "on" -- nloptr stops just
  # inside.
  #
  # `lb > lo` / `ub < hi` is the whole point: it says the binding edge is
  # ADMIXR2'S box and not a bound the model itself declared, which is the only
  # case worth warning about. Testing `!is.finite(lo)` instead -- i.e. "warn only
  # if the model declared no bound on that side at all" -- silently drops every
  # hit on a parameter that has one, even when that bound is nowhere near and the
  # box is what actually stopped the fit. A residual-error parameter carries a
  # lower bound, so exactly the parameters most likely to run away were the ones
  # that could never report it.
  tol <- grad_bounds * 1e-3
  lb  <- pmax(lo, p0 - grad_bounds)
  ub  <- pmin(hi, p0 + grad_bounds)
  hit <- ((p - lb) <= tol & lb > lo) | ((ub - p) <= tol & ub < hi)
  hit[is.na(hit)] <- FALSE
  if (!any(hit)) return(invisible(character(0)))
  # Same order .admBuildOptVec() packs the vector in (parse.R).
  nms <- tryCatch(c(pinfo$struct_names, pinfo$sigma_names, pinfo$omega_par_names),
                  error = function(e) NULL)
  lab <- if (!is.null(nms) && length(nms) >= n) nms[seq_len(n)][hit]
         else paste0("p", which(hit))
  .txt <- sprintf(
    paste0("admixr2: %s finished on the gradient box constraint (grad_bounds = %g ",
           "from the starting value), not at an interior optimum. The reported ",
           "estimate and SE are those of a constrained fit. Widen grad_bounds, or ",
           "start closer to the expected value."),
    paste(lab, collapse = ", "), grad_bounds)
  # BOTH channels, deliberately.
  #
  # warning() is what lands in fit$warnings, so print(fit) surfaces it later --
  # but nlmixr2est muffles conditions inside nlmixr2Est.*, so it never reaches
  # warnings() and a batch script that only writes coefficients to disk sees
  # nothing. That matters more since 0.4.1: adfo's grad default moved
  # "none" -> "analytical", which also flips want_grad and therefore imposes this
  # box where the derivative-free path used -Inf/Inf. A start value far from the
  # optimum now converges to a corner with nloptr reporting success and a finite
  # estimate and SE printed. (admc and adgh have always defaulted to a gradient
  # AND grad_bounds = 5, so this aligned adfo with them rather than singling it
  # out -- which is why the box stays and the reporting is what changes.)
  #
  # message() is the channel the live progress table already uses, so it is
  # visible as the fit runs. Redundant in an interactive session; the only way to
  # see it in a non-interactive one.
  message(.txt)
  warning(.txt, call. = FALSE)
  invisible(lab)
}

# Reconcile a user-chosen nloptr algorithm with the gradient mode.
#   * algorithm = NULL (unset) -> pick the default matching `grad` (no message).
#   * grad == "none" but a gradient-based algorithm was chosen -> there is no
#     gradient to give nloptr, so fall back to BOBYQA (with a message).
#   * grad != "none" but a derivative-free algorithm was chosen -> the gradient
#     cannot be used, so turn it off (with a message).
# Validates explicit algorithm names against the installed nloptr.
# Returns list(algorithm = <chr>, grad = <chr>).
.admResolveAlgorithm <- function(algorithm, grad, .var.name = "algorithm") {
  # Unset -> the default that matches the gradient mode; always consistent.
  if (is.null(algorithm)) return(list(algorithm = .admDefaultAlgorithm(grad),
                                       grad = grad))

  checkmate::assertString(algorithm, .var.name = .var.name)
  # Validate early against the installed nloptr when we can; if the query failed
  # (empty), defer to nloptr -- it rejects bad names and lists the valid ones.
  valid <- .admNloptrAlgorithms()
  if (length(valid) && !algorithm %in% valid)
    stop(sprintf(
      "%s: '%s' is not a valid nloptr algorithm. See nloptr::nloptr.print.options() for the full list.",
      .var.name, algorithm), call. = FALSE)

  # AUGLAG / MLSL are meta-algorithms requiring a subsidiary local optimiser
  # (local_opts) that the control objects do not expose -- warn up front rather
  # than surface a cryptic nloptr error at fit time.
  if (grepl("AUGLAG|MLSL", algorithm))
    warning(sprintf(
      "%s: '%s' needs a subsidiary local optimiser (local_opts) that admixr2 does not configure; it may fail.",
      .var.name, algorithm), call. = FALSE)

  algo_grad <- .admAlgoNeedsGrad(algorithm)

  # grad == "none" -> derivative-free optimisation. A gradient-based algorithm
  # has no gradient to consume, so fall back to BOBYQA.
  if (grad == "none" && algo_grad) {
    message(sprintf(
      "%s: '%s' is gradient-based but grad = 'none'; using 'NLOPT_LN_BOBYQA'.",
      .var.name, algorithm))
    algorithm <- "NLOPT_LN_BOBYQA"

  # grad != "none" -> a gradient is computed. A derivative-free algorithm cannot
  # use it, so turn the gradient off.
  } else if (grad != "none" && !algo_grad) {
    message(sprintf(
      "%s: '%s' is derivative-free; gradient ('%s') is unused (grad set to 'none').",
      .var.name, algorithm, grad))
    grad <- "none"
  }

  list(algorithm = algorithm, grad = grad)
}

# -- Shi (2021) adaptive CENTRAL-difference intervals ---------------------------
#
# Gill83, above, answers "what forward step balances condition error against
# curvature". Shi/Xie/Xu/Nocedal (2021) answer the same question for a CENTRAL
# difference, and measurement says they answer it far better on this package's
# objectives. Scored against the analytic gradient (exact; cross-checked to 3e-9
# against a high-accuracy central difference), max relative error over the five
# parameters of the integration model:
#
#                       adirmc inner NLL   adfo NLL
#   forward fixed 1e-4       7.3e-04        8.6e-04
#   forward fixed 1e-6       7.3e-06        9.8e-06
#   gill83 forward           8.1e-04        7.9e-04
#   shi21 central            7.0e-08        9.5e-08
#
# Gill83 is not merely beaten, it is beaten by the FIXED step it exists to
# improve on: its measured steps come out 5e-5..1.7e-3 where the objective wants
# ~1e-8..1e-6. It is not mis-implemented here -- FOCEI's own defaults are what
# the exported wrapper reaches, and those are
# tuned for a per-subject objective carrying real solver noise, not for an
# aggregate objective evaluated to near machine precision.
#
# WHY REIMPLEMENTED rather than called. nlmixr2est HAS this algorithm, as
# `shi21CentralWrap`, but it is not exported -- and admixr2 makes zero `:::`
# calls into nlmixr2est, a rule that covers reaching in via asNamespace() just as
# much as the `:::` token. Reimplementing also buys the thing the gill83 path
# cannot have: `eps_f` is a real argument here. It is the single input that
# matters (h* scales as eps_f^(1/3)), and nlmixr2Gill83's wrapper drops every
# tuning argument it accepts, so its noise assumption is unreachable. Measured
# against the upstream routine as an ORACLE in test-optim-steps-shi.R.
#
# The maths. For a central difference the two error terms are
#     truncation  (h^2/6)|f'''|      noise  eps_f/h
# minimised together at
#     h* = (3 * eps_f / |f'''|)^(1/3)
# and the whole job is estimating |f'''| without knowing it. The symmetric third
# difference does that:
#     D3(h) = f(p+2h) - 2f(p+h) + 2f(p-h) - f(p-2h)  ~  2 h^3 f'''
# It is used only when it stands clear of the noise floor: its four evaluations
# carry coefficients (1,2,2,1), so noise in D3 has scale sqrt(1+4+4+1) = 3.16
# eps_f, and a D3 below a multiple of that is measuring nothing but noise -- the
# signal that h is too SMALL, which is the one failure a fixed step cannot detect
# and the reason a too-fine step degrades so sharply (a fixed 1e-6 on an
# objective with 1e-8 relative noise measured 2.55 relative error, i.e. no
# correct digits at all).
#
# Returns list(h, gr): the chosen step and the central-difference derivative at
# it, per requested index. `gr` is free -- the last iterate already evaluated it.
.admShi21Central <- function(fn, p, k, eps_f, h0 = NULL, maxiter = 10L) {
  scale <- max(abs(p[k]), 0.1)
  # Start where a unit third derivative would put the optimum, but never finer
  # than the scale-relative cube root of machine precision.
  h <- if (!is.null(h0)) h0 else
    max((3 * eps_f)^(1/3), scale * .Machine$double.eps^(1/3))
  at <- function(d) { q <- p; q[k] <- q[k] + d; fn(q) }
  # Do NOT iterate h to a fixed point of h* = (3 eps_f/|f'''|)^(1/3). That was
  # tried and it CANNOT converge, for a reason worth recording: at the optimum
  # h*^3 = 3 eps_f/|f'''|, so the third difference there is d3 ~ 2 h*^3 |f'''| =
  # 6 eps_f -- only about twice the sqrt(1+4+4+1) = 3.16 eps_f of noise it
  # carries. Any floor high enough to trust d3 therefore REJECTS h* itself, and
  # the iteration limit-cycles: grow because d3 is under the floor, refine back
  # down to h*, get rejected, grow again. It ran to maxiter every time, at 4
  # evaluations an iteration, and returned the starting guess unrefined.
  #
  # |f'''| is a property of the FUNCTION, not of h, so estimate it ONCE at a
  # probe step deliberately coarse enough for d3 to stand far clear of the noise
  # (100x), then evaluate the formula. Monotone, terminating, and it uses the
  # regime where the third difference is actually informative.
  d3_noise <- 3.1623 * eps_f
  f3 <- NA_real_
  for (i in seq_len(maxiter)) {
    fp2 <- at(2 * h); fp1 <- at(h); fm1 <- at(-h); fm2 <- at(-2 * h)
    if (!all(is.finite(c(fp2, fp1, fm1, fm2)))) { h <- h / 4; next }
    d3 <- fp2 - 2 * fp1 + 2 * fm1 - fm2
    if (abs(d3) >= 100 * d3_noise) { f3 <- abs(d3) / (2 * h^3); break }
    h <- h * 4
    if (h > 1e3 * scale) break     # flat in this parameter; keep the last h
  }
  # `measured` is load-bearing, not decoration. When the third difference never
  # clears the noise floor there is no |f'''| and hence no h*, and the loop exits
  # holding whatever h it had grown to -- which is the CAP, i.e. a step the size
  # of the parameter. Returning that silently would be the worst outcome of the
  # three: a converged NLL is near-quadratic, so |f'''| is smallest exactly where
  # fits end up, and a step that size perturbs omega out of positive-definiteness
  # or sigma to Inf. The caller substitutes its fixed fallback instead.
  measured <- is.finite(f3) && f3 > 0
  hs <- if (measured) (3 * eps_f / f3)^(1/3) else h
  # Never so fine it is pure cancellation, nor larger than the parameter itself.
  hs <- min(max(hs, scale * 1e-12), scale)
  q1 <- p; q1[k] <- q1[k] + hs; q2 <- p; q2[k] <- q2[k] - hs
  list(h = hs, gr = (fn(q1) - fn(q2)) / (2 * hs), measured = measured)
}

# Per-parameter central steps: a numeric vector
# the length of `idx`, never NA, never non-positive, falling back per parameter
# rather than all-or-nothing.
#
# `eps_f` is the ABSOLUTE noise level of `fn`. Supplied by the caller when it
# knows one; otherwise estimated by .admEcNoise(). Do NOT default it to a
# relative quantity -- h* scales as eps_f^(1/3), so a value wrong by 10^6 moves
# the step by 10^2, which is the whole difference between the gill83 and shi21
# rows in the table above.
.admShi21Steps <- function(fn, p, idx = seq_along(p), fallback = NULL,
                           eps_f = NULL, .var.name = "cov") {
  if (is.null(fallback))
    fallback <- pmax(abs(p[idx]), 0.1) * .Machine$double.eps^(1/3)
  fallback <- rep_len(as.numeric(fallback), length(idx))
  if (is.null(eps_f))
    eps_f <- tryCatch(.admEcNoise(fn, p, idx[[1L]]), error = function(e) NA_real_)
  if (!is.finite(eps_f) || eps_f <= 0) {
    warning(sprintf(
      "%s: could not estimate the objective's noise level -- falling back to the fixed step scale.",
      .var.name), call. = FALSE)
    return(fallback)
  }
  h <- vapply(seq_along(idx), function(j) {
    r <- tryCatch(.admShi21Central(fn, p, idx[[j]], eps_f), error = function(e) NULL)
    if (is.null(r) || !isTRUE(r$measured) || !is.finite(r$h) || r$h <= 0)
      fallback[[j]] else r$h
  }, numeric(1))
  h
}

# Moré & Wild (2011) ECnoise: the noise level of `fn` from a short line of
# equally spaced evaluations.
#
# Needed because eps_f is the one input .admShi21Central() cannot guess. The
# construction: differences of order j annihilate any polynomial of degree < j,
# so once j exceeds the local smoothness of `fn` the j-th differences are pure
# noise, amplified by a known factor gamma_j = (j!)^2/(2j)!. The estimate is the
# level at which successive estimates stop moving AND the differences alternate
# in sign -- both, because either alone is fooled (a still-smooth level looks
# stable; a genuinely noisy one need not alternate every single time).
.admEcNoise <- function(fn, p, k, h = NULL, m = 6L) {
  scale <- max(abs(p[k]), 0.1)
  if (is.null(h)) h <- scale * .Machine$double.eps^(1/3)
  d <- seq(-m/2, m/2) * h
  f <- vapply(d, function(dd) { q <- p; q[k] <- q[k] + dd; fn(q) }, numeric(1))
  if (!all(is.finite(f))) return(NA_real_)
  # A perfectly constant slice carries no noise information at all.
  if (all(f == f[[1L]])) return(NA_real_)
  tab <- f
  gamma <- 1
  est <- rep(NA_real_, m)
  alt <- rep(NA_real_, m)
  for (j in seq_len(m)) {
    tab <- diff(tab)
    gamma <- gamma * (j^2) / (2 * j * (2 * j - 1))
    est[j] <- sqrt(gamma * mean(tab^2))
    s <- sign(tab); s <- s[s != 0]
    alt[j] <- if (length(s) < 2L) 0 else mean(s[-1] != s[-length(s)])
  }
  for (j in seq_len(m - 2L)) {
    trio <- est[j:(j + 2L)]
    if (any(!is.finite(trio)) || any(trio <= 0)) next
    if (max(trio) / min(trio) <= 100 && alt[j] >= 0.5) return(est[[j]])
  }
  # Nothing stabilised: the smallest positive estimate is the least-bad answer,
  # and is at worst an OVER-estimate, which errs toward a coarser step.
  ok <- is.finite(est) & est > 0
  if (!any(ok)) NA_real_ else min(est[ok])
}

# The central-difference Hessian, given steps and a scalar objective.
#
# Diagonal by the second difference and off-diagonal by the four-point cross:
#
#   H_kk = [ f(p + h_k) - 2 f(p) + f(p - h_k) ] / h_k^2
#   H_ij = [ f(++) - f(+-) - f(-+) + f(--) ] / (4 h_i h_j)
#
# Extracted because .adfoCalcCov and .adghCalcCov carried it BYTE-IDENTICALLY
# over nineteen lines. The stencil and its normalisation are a decision -- the
# `-2 f(p)` on the diagonal, the 4 in the cross denominator, and the fact that
# both use the SAME per-parameter step vector .admHessSteps measured -- so two
# copies are two chances for one estimator's reported covariance to be computed
# differently from another's.
#
# .admCalcCov does NOT call this: it builds the same point set but hands it to
# .admNLLBatch, because an MC objective is far cheaper evaluated in one batched
# solve than point by point. Same stencil, different evaluation strategy; if
# the stencil changes here it has to change there too.
#
# `nll_fn` must already be memoised or cheap -- .adfoCalcCov primes a memo over
# these exact points first, which is why it is passed in rather than built.
.admFDHessian <- function(p_hat, cov_idx, h_fd, nll0, nll_fn) {
  np <- length(cov_idx)
  H  <- matrix(0, np, np)
  for (k in seq_len(np)) {
    ki  <- cov_idx[k]; hk <- h_fd[k]
    p_p <- p_hat; p_p[ki] <- p_p[ki] + hk
    p_m <- p_hat; p_m[ki] <- p_m[ki] - hk
    H[k, k] <- (nll_fn(p_p) - 2 * nll0 + nll_fn(p_m)) / hk^2
  }
  if (np > 1L)
    for (i in seq_len(np - 1L)) {
      for (j in seq(i + 1L, np)) {
        ii <- cov_idx[i]; ji <- cov_idx[j]
        hi <- h_fd[i];    hj <- h_fd[j]
        p_pp <- p_hat; p_pp[ii] <- p_pp[ii] + hi; p_pp[ji] <- p_pp[ji] + hj
        p_pm <- p_hat; p_pm[ii] <- p_pm[ii] + hi; p_pm[ji] <- p_pm[ji] - hj
        p_mp <- p_hat; p_mp[ii] <- p_mp[ii] - hi; p_mp[ji] <- p_mp[ji] + hj
        p_mm <- p_hat; p_mm[ii] <- p_mm[ii] - hi; p_mm[ji] <- p_mm[ji] - hj
        H[i, j] <- H[j, i] <-
          (nll_fn(p_pp) - nll_fn(p_pm) - nll_fn(p_mp) + nll_fn(p_mm)) /
          (4 * hi * hj)
      }
    }
  H
}
