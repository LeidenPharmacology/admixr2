# Optimizer plumbing: nloptr algorithm selection, Gill (1983) finite-difference
# step selection, and the box-constraint warning.
#
# Split out of utils.R; contents unchanged (see R/covreport.R for why that is
# safe). The two step accessors .admGH()/.admGH0() look collapsible and are not:
# indexing a per-parameter Gill vector by an ETA number would silently recycle
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

# Per-parameter finite-difference steps by Gill, Murray, Saunders & Wright (1983),
# via nlmixr2est's own implementation.
#
# Everything in admixr2 that finite-differences picks its step from the same
# heuristic: `pmax(abs(p), 0.1) * h`, with h a fixed constant (`cov_h_outer`
# defaults to eps^(1/5)). That is a guess about one thing -- how much noise the
# objective carries -- applied identically to every parameter, and it is the
# guess behind the "Hessian not positive definite ... try increasing
# cov_h_outer" warning: too small a step and the difference is dominated by
# solver noise; too large and it is dominated by curvature the second-order term
# does not capture. A parameter the objective is flat in and one it is sharp in
# want different steps, and the right step also moves with the ODE tolerance.
#
# Gill83 measures instead of guessing: it probes the objective, estimates the
# condition error and the second derivative, and returns the step where the two
# error sources balance. This is the algorithm FOCEI uses to choose its own
# steps, so admixr2 and nlmixr2est now agree about what a finite difference is.
#
# nlmixr2est::nlmixr2Gill83 is EXPORTED (no ::: call), and returns a data.frame
# with ONE ROW PER PARAMETER of `p` -- `hf` (the forward step), `df` (the
# derivative at it), `err`, and `info`, a factor whose non-"Good" levels flag a
# parameter whose step could not be assessed. Those fall back to the caller's
# heuristic rather than being trusted.
#
# Two things about the upstream function worth knowing:
#  * `which` MASKS rather than subsets -- an excluded parameter still occupies a
#    row, as "Not Assessed" with hf = NA -- so the result is indexed by `idx`,
#    not consumed in order.
#  * its wrapper drops the tuning arguments it accepts (gillRtol / gillK /
#    gillStep / gillFtol are re-supplied as defaults in the inner call), so only
#    the defaults are reachable through the exported entry point. That is fine
#    here -- they are FOCEI's own defaults -- but do not add a control for them
#    expecting it to have any effect.
#
# Cost: up to gillK (10) extra objective evaluations per parameter, ONCE. So it
# must never run INSIDE a gradient the optimizer calls thousands of times -- but
# it can perfectly well choose the step that gradient then reuses, which is
# exactly what FOCEI does (`numericGrad` runs gill83 on the first evaluation,
# `nF == 1`, stores aEps/rEps per parameter and finite-differences with them for
# the rest of the fit). See .admGillGradH().
#
# Returns a numeric vector the length of `idx`, never NA, never non-positive.
.admGillSteps <- function(fn, p, idx = seq_along(p), fallback = NULL,
                          .var.name = "cov") {
  if (is.null(fallback))
    fallback <- pmax(abs(p[idx]), 0.1) * .Machine$double.eps^(1/5)
  fallback <- rep_len(as.numeric(fallback), length(idx))
  if (!requireNamespace("nlmixr2est", quietly = TRUE)) return(fallback)
  .which <- rep(FALSE, length(p))
  .which[idx] <- TRUE
  g <- tryCatch(
    nlmixr2est::nlmixr2Gill83(fn, p, which = .which),
    error = function(e) NULL)
  # One row per parameter of `p` -- NOT one per requested index. A parameter
  # `which` excludes still gets a row, with info = "Not Assessed" and hf = NA.
  # (Checked against nlmixr2est rather than assumed: reading it the other way
  # silently returned the fallback for every subset request.)
  if (is.null(g) || is.null(g$hf) || length(g$hf) != length(p)) {
    warning(sprintf(
      "%s: Gill83 step selection failed -- falling back to the fixed step scale.",
      .var.name), call. = FALSE)
    return(fallback)
  }
  h <- as.numeric(g$hf)[idx]
  # `info` is a factor. Which levels to accept is NOT "only Good" -- that reading
  # made the whole option a no-op and it took a probe to notice:
  #
  #   nlmixr2Gill83(function(p) sum((p - 1)^2), 0.5)
  #   #>   info            hf         df df2
  #   #>   High Grad Error 0.000149505 -0.9998505 2
  #
  # A noiseless quadratic, an exact second derivative, a textbook-perfect step --
  # reported as "High Grad Error", because that level means "the derivative
  # estimate's error is `fTol` or more of the derivative" and `fTol` defaults to
  # ZERO. Any error at all trips it, so it is the ORDINARY return, not a failure.
  # FOCEI agrees in the only way that matters: `numericGrad` takes gill83's `hf`
  # into aEps/rEps unconditionally, whatever the return code.
  #
  # The levels that genuinely leave no step to use are the ones describing a
  # function the algorithm cannot fit a step to at all:
  #   "Not Assessed"    -- masked out by `which`; hf is NA anyway
  #   "Constant Grad"   -- flat in this parameter, so no h is better than another
  #   "Odd/Linear Grad" -- df2 ~ 0, and h balances condition error against df2
  # "Grad changes quickly" IS kept: FOCEI's response to it is to use the step now
  # and re-measure later (repeatGill), and using it beats the blanket constant --
  # but we do not re-measure, so it is the least confident step we accept.
  ok <- is.finite(h) & h > 0
  if (!is.null(g$info))
    ok <- ok & !(as.character(g$info)[idx] %in%
                   c("Not Assessed", "Constant Grad", "Odd/Linear Grad"))
  h[!ok] <- fallback[!ok]
  h
}

# Index a finite-difference step that may be a SCALAR -- the fixed heuristic, one
# constant shared by every parameter -- or a per-parameter VECTOR of Gill83's
# measured steps.
#
# Every FD site in the package reads its step through this or .admGH0(), which is
# what lets `gill = TRUE` reach the optimizer's gradient without a single
# signature change: the driver hands the gradient a vector where it used to hand
# it a number. A scalar takes the identity branch, so `gill = FALSE` (the
# default) is bit-for-bit what it was.
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
# So a Gill vector carries the scalar it was built from as attr(, "h0"), and
# index-less sites read that. The eta step therefore stays exactly what it was
# before `gill = TRUE`, which is also the honest answer: Gill83 measured the
# OBJECTIVE, and d(pred)/d(eta) is a different function.
.admGH0 <- function(h)
  if (length(h) == 1L) h else {
    h0 <- attr(h, "h0", exact = TRUE)
    if (is.null(h0)) h[[1L]] else h0
  }

# Per-parameter FD steps for the OPTIMIZER's gradient, measured once at the point
# the fit starts from.
#
# The covariance Hessian can afford to call .admGillSteps() directly: it runs
# once, post-fit, at a known parameter vector. A gradient cannot -- the optimizer
# calls it thousands of times. FOCEI's answer, which this mirrors, is to measure
# once and REUSE: gill83 runs at the first gradient evaluation, its `hf` is
# folded into per-parameter aEps/rEps, and every later difference is taken with
# those. Here the measurement happens in the driver, before the optimizer is
# handed anything, and the result travels as the `grad_h` argument that was
# already there.
#
# `hf` is an ABSOLUTE step at `p`. The estimators express theirs differently --
# adfo/adgh scale by `pmax(abs(p), 0.1)`, admc uses the raw number -- so the
# measured step is divided by whatever that site will multiply it back by
# (`scaled`). The step is therefore exactly Gill's at `p`, and tracks the
# parameter afterwards under the convention that site already had, rather than
# freezing an absolute number that stops making sense once the optimizer has
# moved a decade.
#
# One deliberate difference from FOCEI, which tracks with `h = |p|*rEps + aEps`
# (aEps = rEps = hf/(|p0| + 1), so `hf * (|p| + 1)/(|p0| + 1)`) against
# admixr2's `hf * max(|p|, 0.1)/max(|p0|, 0.1)`. Both reproduce `hf` exactly at
# the point it was measured and both scale mildly and monotonically away from it;
# what Gill83 actually contributes -- the MAGNITUDE, which moves by orders of
# magnitude between parameters -- is identical either way. Keeping admixr2's own
# law means `gill = TRUE` changes one thing rather than two, and every FD site
# keeps a single scaling convention.
#
# `idx` is the set of parameters that will actually be finite-differenced;
# everything else keeps the constant and costs nothing. Passing `integer(0)` --
# a fit whose gradient is fully analytic -- skips the probe entirely.
.admGillGradH <- function(fn, p, idx, grad_h, scaled = TRUE, .var.name = "grad") {
  if (length(idx) == 0L) return(grad_h)
  out <- rep_len(as.numeric(grad_h), length(p))
  base <- if (scaled) pmax(abs(p[idx]), 0.1) else rep(1, length(idx))
  hf <- .admGillSteps(fn, p, idx, fallback = base * grad_h, .var.name = .var.name)
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
