# Optimizer plumbing: nloptr algorithm selection, Shi (2021) finite-difference
# step selection, and box-constraint warnings.
# .admGH()/.admGH0() separate per-parameter steps from index-less eta perturbations
# (which must not recycle across n_sim draws).


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

# Return scalar fallback or per-parameter step for parameter `idx`.
.admGH <- function(h, idx) if (length(h) == 1L) h else h[idx]

# Step for differences without a parameter index (e.g., eta perturbation in .admGrad).
# Uses attr(, "h0") from measured-step vectors to avoid recycling across draws.
.admGH0 <- function(h)
  if (length(h) == 1L) h else {
    h0 <- attr(h, "h0", exact = TRUE)
    if (is.null(h0)) h[[1L]] else h0
  }

# Per-parameter FD steps for the covariance Hessian.
# Uses ECnoise-measured eps_f with second-difference scaling h ~ eps_f^(1/4)
# (scaled proportionally by cov_h_outer / .ADM_COV_H_REF).
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
  # adgh uses finer default cov_h_outer; scale relative to reference keeps user intent
  pmax(abs(p[idx]), 0.1) * (eps_f / max(f0, 1))^(1/4) *
    (cov_h_outer / .ADM_COV_H_REF)
}

# Measure per-parameter FD steps once at initial `p` and reuse across optimizer iterations.
# `idx` specifies the parameters to finite-difference; empty vector skips the probe.
.admShi21GradH <- function(fn, p, idx, grad_h, scaled = TRUE, .var.name = "grad") {
  if (length(idx) == 0L) return(grad_h)
  out <- rep_len(as.numeric(grad_h), length(p))
  base <- if (scaled) pmax(abs(p[idx]), 0.1) else rep(1, length(idx))
  hf <- .admShi21Steps(fn, p, idx, fallback = base * grad_h, .var.name = .var.name)
  out[idx] <- hf / base
  # Scalar step for index-less differences (.admGH0)
  attr(out, "h0") <- as.numeric(grad_h)[[1L]]
  out
}

# Warn when a gradient-mode fit finishes on the synthetic box constraint [p0 +/- grad_bounds]
# rather than an interior optimum or user-declared model bound.
.admWarnOnBounds <- function(p, p0, ov, grad_bounds, pinfo) {
  if (is.null(p) || is.null(p0) || !is.finite(grad_bounds) || grad_bounds <= 0)
    return(invisible(character(0)))
  n <- min(length(p), length(p0))
  if (n == 0L) return(invisible(character(0)))
  p  <- p[seq_len(n)]; p0 <- p0[seq_len(n)]
  lo <- if (is.null(ov$lower)) rep(-Inf, n) else ov$lower[seq_len(n)]
  hi <- if (is.null(ov$upper)) rep(Inf,  n) else ov$upper[seq_len(n)]
  # Check if solution sits within 0.1% of admixr2's synthetic box (lb > lo / ub < hi)
  tol <- grad_bounds * 1e-3
  lb  <- pmax(lo, p0 - grad_bounds)
  ub  <- pmin(hi, p0 + grad_bounds)
  hit <- ((p - lb) <= tol & lb > lo) | ((ub - p) <= tol & ub < hi)
  hit[is.na(hit)] <- FALSE
  if (!any(hit)) return(invisible(character(0)))
  # Parameter names matching .admBuildOptVec ordering
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
  # Emit to both message (live display) and warning (fit$warnings)
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
# Balances truncation error (h^2/6)|f'''| against noise eps_f/h via
#   h* = (3 * eps_f / |f'''|)^(1/3)
# where |f'''| is estimated from the symmetric third difference D3(h).
# Beats gill83's forward step by orders of magnitude on this package's
# objectives (max relative error vs analytic gradient: gill83 ~8e-4,
# shi21 central ~1e-7). Reimplemented rather than calling nlmixr2est's
# unexported shi21CentralWrap -- admixr2 makes zero `:::`-equivalent calls,
# and the wrapper drops the `eps_f` tuning argument this code needs.
# Returns list(h, gr, measured).
.admShi21Central <- function(fn, p, k, eps_f, h0 = NULL, maxiter = 10L) {
  scale <- max(abs(p[k]), 0.1)
  # Start where a unit third derivative would put the optimum, bounded below
  # by scale * eps^(1/3)
  h <- if (!is.null(h0)) h0 else
    max((3 * eps_f)^(1/3), scale * .Machine$double.eps^(1/3))
  at <- function(d) { q <- p; q[k] <- q[k] + d; fn(q) }
  # Estimate |f'''| at a probe step coarse enough for D3 to clear the noise floor
  # (3.16 * eps_f) by 100x; fixed-point iteration on h* does not converge (at h*
  # itself D3 ~ 6 eps_f, only ~2x the noise floor, so it gets rejected and the
  # loop regrows/shrinks to maxiter without refining).
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
  # When D3 never clears noise, fall back to caller's default rather than the loop cap
  measured <- is.finite(f3) && f3 > 0
  hs <- if (measured) (3 * eps_f / f3)^(1/3) else h
  # Guard against excessive cancellation or oversized perturbation
  hs <- min(max(hs, scale * 1e-12), scale)
  q1 <- p; q1[k] <- q1[k] + hs; q2 <- p; q2[k] <- q2[k] - hs
  list(h = hs, gr = (fn(q1) - fn(q2)) / (2 * hs), measured = measured)
}

# Per-parameter central steps with per-parameter fallback.
# eps_f is the absolute noise level of fn (estimated via .admEcNoise if NULL).
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

# Moré & Wild (2011) ECnoise: estimates noise level eps_f of fn from equally
# spaced evaluations where differences stop moving and alternate in sign.
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
  # Smallest positive estimate if noise did not stabilise (errs toward coarser step)
  ok <- is.finite(est) & est > 0
  if (!any(ok)) NA_real_ else min(est[ok])
}
