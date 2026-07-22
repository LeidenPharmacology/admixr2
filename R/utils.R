#' @importFrom stats cov dnorm pnorm qnorm rnorm runif setNames
#' @importFrom utils assignInNamespace head
#' @importFrom Rcpp sourceCpp
#' @useDynLib admixr2, .registration = TRUE
NULL

# Suppress R CMD check notes for ggplot2 NSE column names used in plot.admFit.
utils::globalVariables(c(
  "time", "pred_lo", "pred_hi", "pred_mean",
  "obs_lo", "obs_hi", "obs_mean",
  "t_col", "t_row", "value",
  "nll", "restart", "iter",
  "z", "z_label", "z_vjust",
  "resid", "lo", "hi"
))

`%||%` <- function(x, y) if (is.null(x)) y else x

# Jacobian of the REPORTED omega entries w.r.t. the optimizer's log-Cholesky.
#
# The optimizer holds L (lower-triangular) with the diagonal stored as
# p = log(Omega_ii) and the off-diagonal as the raw L_ij; .admFullTheta() reports
# the Omega entries themselves. Omega = L L', so
#
#   d(Omega)/d(L_ab)  = E_ab L' + L E_ab'   =>   d(Omega_ij)/d(L_ab)
#                                              = delta_ia L_jb + delta_ja L_ib
#   d(L_ab)/d(p)      = L_aa/2 on the diagonal (L_aa = exp(p/2)), 1 off it.
#
# This is NOT diagonal once omega is correlated, which is why a per-row delta
# factor -- fine for sigma -- cannot be used here. Returns an (n_report x n_par)
# matrix so the covariance transforms as J %*% cov %*% t(J).
#
# WHY NOT REUSE nlmixr2est/rxode2 -- checked, and the answer is structural:
#
#   * nlmixr2est never needs this Jacobian. `.foceiCalcRanalytic()` builds its R
#     matrix on NATURAL-scale directions from the start (via `.omegaVarCovDeriv()`
#     -> `rxode2::rxOmegaVarCovDeriv()`, the derivatives of Omega^-1 and
#     log|Omega| w.r.t. the Omega ELEMENTS), so its result is already on the
#     reported scale -- hence the variable is literally called `.covNat`. There is
#     no transform of ours to borrow. We cannot borrow the assembler either: it
#     needs FOCEI's EBEs, `etaObf` and individual-level `dataSav`, none of which
#     exist in an aggregate fit. admixr2's covariance is a NUMERICAL Hessian of the
#     OPTIMIZER-scale objective, so a change of variables is unavoidable.
#   * `rxode2::rxOmegaVarCovDeriv()` differentiates the wrong way round (w.r.t.
#     Omega's own entries, not w.r.t. a Cholesky parameter).
#   * `rxode2::rxSymInvCholCreate()` is a Cholesky of **Omega^-1** with
#     `diag.xform` in {sqrt, log, identity} -- a different parameterisation from
#     admixr2's (Cholesky of Omega, diagonal held as log(Omega_ii)), so its thetas
#     are not ours and its derivatives are not the ones we need.
#
# What IS shared with upstream, deliberately: the (row, col) ENUMERATION of the
# free Omega entries -- pinned to `rxOmegaVarCovDeriv()$elements` in
# test-cov-reporting.R -- and the resulting names (`.foceiOmegaCovNames`'s
# `om.<eta>` / `cov.<eta_i>.<eta_j>`). The delta-transform PATTERN is theirs too:
# `.postEstimationBoundedTransformJacobian()` does `env$cov <- J cov J'` for
# bounded structural thetas, with a diagonal J; the sigma factors in
# .admSigmaReportJac() are that same pattern, and omega is the case where J is
# not diagonal.
.admOmegaJacobian <- function(pinfo, L) {
  n_o <- length(pinfo$omega_par)
  if (n_o == 0L) return(NULL)
  ci <- pinfo$chol_i; cj <- pinfo$chol_j; cd <- pinfo$chol_diag
  J  <- matrix(0, n_o, n_o)
  for (r in seq_len(n_o)) {           # reported entry Omega[i, j]
    i <- ci[r]; j <- cj[r]
    for (k in seq_len(n_o)) {         # parameter p_k -> L[a, b]
      a <- ci[k]; b <- cj[k]
      dOm <- (if (i == a) L[j, b] else 0) + (if (j == a) L[i, b] else 0)
      J[r, k] <- dOm * (if (cd[k]) L[a, a] / 2 else 1)
    }
  }
  J
}

# Names for the reported omega entries, matching nlmixr2est's own convention
# exactly (.foceiOmegaCovNames): `om.<eta>` on the diagonal, `cov.<eta_i>.<eta_j>`
# off it, with i the LATER eta -- which is how admixr2 stores chol_i/chol_j too
# (they come from iniDf's neta1/neta2, i.e. the lower triangle).
.admOmegaReportNames <- function(pinfo) {
  et <- pinfo$eta_names %||% pinfo$eta_col_names
  if (is.null(et)) return(pinfo$omega_par_names)
  vapply(seq_along(pinfo$omega_par), function(r) {
    i <- pinfo$chol_i[r]; j <- pinfo$chol_j[r]
    if (i == j) paste0("om.", et[i]) else paste0("cov.", et[i], ".", et[j])
  }, character(1))
}

# Put the dimnames back on fit$env$cov after nlmixr2est has been through it.
#
# nlmixr2est's C++ `foceiFitCpp_` re-dimnames whatever covariance it finds in the
# fit environment using its OWN parameter-name vector, which knows about the
# thetas only -- so the omega rows we append come back named "". This is a known
# shape upstream, not a bug of ours: nlmixr2est ships `.impmapNameCov()` to
# repair exactly these blanks for its importance-sampling estimator, reading the
# omega names off the model. We do the same thing, except we already hold the
# authoritative names (we built the block), so we restore them verbatim.
#
# `nms` must be SNAPSHOT with .admCovNames() before the matrix is handed to
# nlmixr2est: foceiFitCpp_ sets the dimnames attribute IN PLACE on the same SEXP
# (no R-level copy happens when a matrix is merely assigned into an environment),
# so by the time we get here the driver's own `.cov` has been blanked as well.
# Reading the names back off it would restore nothing.
#
# Guarded on the length matching: if a future nlmixr2est returns a covariance of
# a different shape, leaving it untouched is the safe outcome -- a wrongly
# labelled SE is far worse than an unlabelled one.
.admCovNames <- function(cov) if (is.matrix(cov)) rownames(cov) else NULL

# Put the theta rows of the covariance in iniDf's OWN order.
#
# nlmixr2est fills its `SE` column POSITIONALLY: it walks the thetas in iniDf
# order and takes the next entry of `sqrt(diag(cov))` for each one it did not
# skip. admixr2 builds the covariance in OPTIMIZER order -- every structural
# theta, then every residual parameter -- and those two orders agree only when
# the model happens to declare its residual parameters last.
#
#   ini({ a <- 0.1; tcl <- log(3); tv <- log(30) })   # residual declared FIRST
#
# printed `a` with tcl's SE, tcl with tv's and tv with a's: a silent rotation,
# every value finite and plausible. So this is not cosmetic ordering -- it is
# what makes the SE belong to the parameter it is printed beside.
#
# Rows that are not thetas (the appended omega block) keep their position at the
# end. Anything unrecognised is left alone: a covariance we cannot map is better
# reported in the order we built it than permuted on a guess.
#
# The other half of the contract is .admCovSkip(), which tells nlmixr2est WHICH
# thetas this matrix carries -- without it, nlmixr2est < 6.2.0 skips every
# residual-error theta (FOCEI computes its covariance without them) and so reads
# the residual's row as the first structural theta's standard error.
.admCovThetaOrder <- function(cov, ui) {
  if (!is.matrix(cov) || is.null(rownames(cov))) return(cov)
  .th <- .admThetaIniDf(ui)
  if (is.null(.th)) return(cov)
  nms  <- rownames(cov)
  # the thetas this matrix carries, in iniDf's order; everything else (the omega
  # block) keeps its position after them, which is where nlmixr2est stops looking
  want <- .th$name[.th$name %in% nms]
  if (!length(want)) return(cov)
  ord <- c(match(want, nms), which(!(nms %in% .th$name)))
  cov[ord, ord, drop = FALSE]
}

# iniDf's theta rows, in ntheta order. NULL when there is nothing to order by.
.admThetaIniDf <- function(ui) {
  iniDf <- tryCatch(ui$iniDf, error = function(e) NULL)
  if (is.null(iniDf) || is.null(iniDf$ntheta)) return(NULL)
  .th <- iniDf[!is.na(iniDf$ntheta), , drop = FALSE]
  if (!nrow(.th)) return(NULL)
  .th[order(.th$ntheta), , drop = FALSE]
}

# The `skipCov` vector to hand nlmixr2est: TRUE for every theta this covariance
# does NOT carry a standard error for.
#
# Derived from the matrix itself rather than from a convention, so the two cannot
# drift apart. nlmixr2est's own default is version-dependent -- 6.2.0 skips only
# fixed thetas, earlier versions skip every residual-error theta as well, because
# FOCEI's covariance genuinely does not include them. admixr2's does, so saying so
# is what makes the residual SE print at all on the older host, and what stops the
# structural SEs being read off the wrong rows there.
#
# NULL (leave nlmixr2est's default alone) when the model has no thetas or the
# covariance is missing, and when we do not hold every theta nlmixr2est expects to
# find -- the same refusal as .admCovThetaOrder(), for the same reason.
.admCovSkip <- function(cov, ui) {
  if (!is.matrix(cov) || is.null(rownames(cov))) return(NULL)
  .th <- .admThetaIniDf(ui)
  if (is.null(.th)) return(NULL)
  # Indexed BY ntheta, not by row position: nlmixr2est checks the length against
  # max(ntheta) and, if it disagrees, silently substitutes a vector of its own.
  skip <- rep(TRUE, max(.th$ntheta))
  skip[.th$ntheta] <- !(.th$name %in% rownames(cov))
  if (all(skip)) return(NULL)
  skip
}

.admRestoreCovNames <- function(fit, nms) {
  if (is.null(nms)) return(invisible(NULL))
  e <- tryCatch(fit$env, error = function(e) NULL)
  if (is.null(e) || is.null(e$cov) || !is.matrix(e$cov)) return(invisible(NULL))
  if (nrow(e$cov) != length(nms)) return(invisible(NULL))
  dimnames(e$cov) <- list(nms, nms)
  invisible(NULL)
}

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

# Map an rxSolve output column name to the one rxSolve actually returns.
# Internal nlmixr2 linCmt names (rxLinCmt, linCmtB, ...) don't appear in the
# simulation model rxSolve output -- use ipredSim which is always present.
.admOutputColName <- function(var)
  if (startsWith(var, "rx") || startsWith(var, "linCmt")) "ipredSim" else var

# The model variable an endpoint's predictions actually live in.
#
# For a residual-error endpoint this is predDf$var itself (`cp ~ add(a)` -> "cp").
# For a COUNT endpoint it is NOT: `y ~ pois(cp)` has predDf$var == "y", the DV
# name, while the quantity that is solved and that admixr2 must read is the
# distribution's ARGUMENT, `cp`. Following predDf$var there sent every solve
# looking for a column that does not exist, which is what made count endpoints
# unreachable. .admCountSpec() recovers the argument from the model line.
.admEndpointVar <- function(ui, i = 1L) {
  pd <- tryCatch(ui$predDf, error = function(e) NULL)
  if (is.null(pd) || !"var" %in% names(pd) || nrow(pd) < i) return("cp")
  v <- as.character(pd$var[i])
  d <- if ("distribution" %in% names(pd)) as.character(pd$distribution[i]) else "norm"
  if (d %in% c("pois", "dpois", "binom", "dbinom", "nbinomMu", "dnbinomMu")) {
    cs <- tryCatch(.admCountSpec(ui, v, d), error = function(e) NULL)
    if (!is.null(cs) && nzchar(cs$mean_var %||% "")) return(cs$mean_var)
  }
  if (d %in% c("beta", "dbeta")) {
    bs <- tryCatch(.admBetaSpec(ui, v), error = function(e) NULL)
    if (!is.null(bs)) return(bs$b1)      # the pair travels on the study/unit
  }
  v
}

# The two shape columns a beta endpoint's prediction is derived from, or NULL.
# Attached to each study as `out_pair` so the solve paths can combine them.
.admBetaPair <- function(ui) {
  pd <- tryCatch(ui$predDf, error = function(e) NULL)
  if (is.null(pd) || !"distribution" %in% names(pd)) return(NULL)
  w <- which(as.character(pd$distribution) %in% c("beta", "dbeta"))
  if (length(w) == 0L) return(NULL)
  bs <- tryCatch(.admBetaSpec(ui, as.character(pd$var[w[1L]])), error = function(e) NULL)
  if (is.null(bs)) return(NULL)
  c(bs$b1, bs$b2)
}

# Detect the primary/default output variable name from ui$predDf (default "cp").
# Used as the fallback output for studies/observations that don't name one.
.admOutputVar <- function(ui) {
  var <- tryCatch(.admEndpointVar(ui, 1L), error = function(e) "cp")
  .admOutputColName(var)
}

# All observable output variable names from ui$predDf (one per model endpoint).
# A multi-endpoint model (e.g. `cp ~ ...; cCSF ~ ...`) has several predDf rows.
.admOutputVars <- function(ui) {
  vars <- tryCatch({
    pd <- ui$predDf
    if (!is.null(pd) && "var" %in% names(pd))
      vapply(seq_len(nrow(pd)), function(i) .admEndpointVar(ui, i), character(1))
    else "cp"
  }, error = function(e) "cp")
  unique(vapply(vars, .admOutputColName, character(1), USE.NAMES = FALSE))
}

# Logical selector over pinfo$sigma_names: which residual-error parameters
# belong to `output`. When the sigma->output mapping is unknown (single-output
# model, or Tier-1 mock iniDf with no `condition` column) every sigma is treated
# as belonging to the one output -- preserving legacy single-output behaviour.
.admSigmaSel <- function(pinfo, output) {
  so <- pinfo$sigma_output
  n  <- length(pinfo$sigma_names)
  if (n == 0L) return(logical(0))
  if (is.null(so) || all(is.na(so)) || is.null(output) || is.na(output))
    return(rep(TRUE, n))
  # `so` holds endpoint names from iniDf$condition (e.g. "cp", "rxLinCmt");
  # `output` is the rxSolve column name. Map endpoints to their column name so
  # linCmt endpoints ("rxLinCmt"/"linCmt*" -> "ipredSim") match. See
  # .admOutputColName().
  so_col <- vapply(so, function(x)
    if (is.na(x)) NA_character_ else .admOutputColName(x), character(1),
    USE.NAMES = FALSE)
  sel <- so_col == output
  sel[is.na(sel)] <- FALSE
  sel
}

# Add each output's residual error to the correct rows of a joint (same-subject)
# predicted covariance. `mu_struct`/`V_pred` are the structural stacked mean and
# covariance; each block's own sigma(s) act only on that block's rows. Returns
# the residual-adjusted mean (`mu`, lnorm-corrected) and covariance (`V`).
.admJointResidual <- function(mu_struct, V_pred, unit, pinfo, sigma_var) {
  n_t <- length(mu_struct)
  arr <- .admResidRows(pinfo, .admRowOutput(unit, n_t), sigma_var, n_t)
  # Row times and the STRUCTURAL covariance are both needed by residual forms that
  # reach the off-diagonal: ar() correlates by time gap, and ordinal's same-time
  # cross-category entry replaces V_struct outright (see .admResidApply). This used
  # to pass neither, so ap$rmat was silently discarded for joint units.
  rt  <- .admRowTimes(unit, n_t)
  ap  <- .admResidApply(mu_struct, diag(V_pred), arr, rt, V_pred)
  if (any(ap$ms != 1, na.rm = TRUE)) V_pred <- V_pred * tcrossprod(ap$ms)   # lnorm off-diagonals
  diag(V_pred) <- ap$dv
  if (!is.null(ap$rmat)) V_pred <- V_pred + ap$rmat
  list(mu = ap$mu, V = V_pred)
}

# Observation time governing each row of a unit's stacked mean vector -- the
# companion of .admRowOutput(). Joint units carry per-output blocks each with
# their own times; a plain unit has one time vector.
.admRowTimes <- function(unit, n_t) {
  if (!is.null(unit$blocks) && length(unit$blocks) > 0L) {
    rt <- rep(NA_real_, n_t)
    # A hand-built unit (Tier-1 tests) may carry blocks with `rows` but no `times`;
    # leave those rows NA rather than erroring. Every consumer guards on
    # length(times) == length(mu), so NA times simply disable the off-diagonal
    # residual forms -- which is right: without times there is no time structure.
    for (blk in unit$blocks)
      if (length(blk$times) == length(blk$rows)) rt[blk$rows] <- blk$times
    return(rt)
  }
  if (!is.null(unit$times) && length(unit$times) == n_t) return(as.numeric(unit$times))
  rep(NA_real_, n_t)
}

# Output column name governing each row of a unit's stacked mean vector.
# Joint units carry per-output blocks; a plain unit is a single output.
.admRowOutput <- function(unit, n_t) {
  if (!is.null(unit$blocks) && length(unit$blocks) > 0L) {
    ro <- rep(NA_character_, n_t)
    for (blk in unit$blocks) ro[blk$rows] <- blk$output
    return(ro)
  }
  rep(unit$output %||% NA_character_, n_t)
}

# Normalise one observed-compartment unit: coerce E, coerce V to matrix,
# auto-detect diagonal, set method + v_diag, inherit n/ev/output from study-level
# `defaults`, and validate dimensions. Returns a self-contained unit.
# V as vector -> treated as variances -> expand to diag matrix, force "var".
# V as matrix with all off-diagonal zeros -> force "var" unless user said "cov".
# V must use ML denominator n (not n-1); use cov.wt(dv_mat, method="ML")$cov.
.admNormaliseObs <- function(ob, label, defaults = list()) {
  ob$n      <- ob$n      %||% defaults$n
  ob$ev     <- ob$ev     %||% defaults$ev
  ob$output <- ob$output %||% defaults$output
  for (f in c("n", "E", "V", "times"))
    if (is.null(ob[[f]])) stop(sprintf("Study '%s' missing '%s'", label, f), call. = FALSE)
  ob$E <- as.numeric(ob$E)
  if (is.vector(ob$V) && !is.list(ob$V)) {
    if (identical(ob$method, "cov"))
      warning(sprintf("Study '%s': V is a vector (variances only) but method='cov' requested -- using method='var'", label), call. = FALSE)
    vv        <- as.numeric(ob$V)
    # diag(x) with a length-1 x treats x as a DIMENSION, not a diagonal value --
    # build the 1x1 matrix explicitly for single-timepoint observations.
    ob$V      <- if (length(vv) == 1L) matrix(vv, 1L, 1L) else diag(vv)
    ob$method <- "var"
  } else {
    ob$V     <- unname(as.matrix(ob$V))
    is_diag <- all(ob$V[lower.tri(ob$V)] == 0) && all(ob$V[upper.tri(ob$V)] == 0)
    ob$method <- if (is_diag && is.null(ob$method)) "var" else
      match.arg(ob$method %||% "cov", c("cov", "var"))
    if (!is_diag && ob$method == "var")
      warning(sprintf("Study '%s': V has non-zero off-diagonal entries but method='var' -- off-diagonal entries will be ignored", label), call. = FALSE)
  }
  n_t <- length(ob$times)
  if (length(ob$E) != n_t)
    stop(sprintf("Study '%s': length(E) (%d) must equal length(times) (%d)",
                 label, length(ob$E), n_t), call. = FALSE)
  if (nrow(ob$V) != n_t || ncol(ob$V) != n_t)
    stop(sprintf("Study '%s': V must be %d x %d to match times", label, n_t, n_t),
         call. = FALSE)
  if (identical(ob$method, "var")) ob$v_diag <- diag(ob$V)
  ob$label <- label
  ob
}

# Coerce a per-observation V spec to a full matrix (matrix as-is; length-1 vector
# -> 1x1; longer vector -> diagonal). Used when assembling a joint covariance.
.admObsVmat <- function(v) {
  if (is.null(v)) return(NULL)
  if (is.vector(v) && !is.list(v)) {
    vv <- as.numeric(v)
    if (length(vv) == 1L) matrix(vv, 1L, 1L) else diag(vv)
  } else unname(as.matrix(v))
}

# Build a single JOINT (same-subject) unit from a study whose observed
# compartments are measured on the SAME subjects: one shared n and ev, a stacked
# observation vector [E_1, ..., E_K] and a joint covariance across all
# compartments. The joint V is supplied either as a study-level full matrix
# (`s$V`, blocks in `observations` order) or assembled from per-observation
# marginal V on the diagonal plus optional cross-covariance blocks (`s$cross`, a
# named list keyed "outA:outB"). Missing cross pairs are zero (block-diagonal).
# Each output is simulated with the SAME random effects and scored by one MVN.
.admBuildJointUnit <- function(s, nm, default_output) {
  onames <- names(s$observations)
  if (is.null(onames) || any(!nzchar(onames)))
    onames <- paste0("obs", seq_along(s$observations))

  blocks <- vector("list", length(s$observations))
  E_list <- vector("list", length(s$observations))
  Vmarg  <- vector("list", length(s$observations))
  row_output <- integer(0); offset <- 0L
  for (k in seq_along(s$observations)) {
    ob     <- s$observations[[k]]
    output <- ob$output %||% default_output
    if (is.null(ob$E) || is.null(ob$times))
      stop(sprintf("Study '%s' observation '%s': joint fits need `E` and `times`.",
                   nm, onames[k]), call. = FALSE)
    tk  <- as.numeric(ob$times); ord <- order(tk); tk <- tk[ord]
    Ek  <- as.numeric(ob$E)
    if (length(Ek) != length(tk))
      stop(sprintf("Study '%s.%s': length(E) (%d) != length(times) (%d)",
                   nm, onames[k], length(Ek), length(tk)), call. = FALSE)
    Ek  <- Ek[ord]
    Vk  <- .admObsVmat(ob$V)
    if (!is.null(Vk)) {
      if (nrow(Vk) != length(tk) || ncol(Vk) != length(tk))
        stop(sprintf("Study '%s.%s': V must be %d x %d to match times",
                     nm, onames[k], length(tk), length(tk)), call. = FALSE)
      Vk <- Vk[ord, ord, drop = FALSE]
    }
    nk <- length(tk)
    blocks[[k]] <- list(name = onames[k], output = output, times = tk,
                        rows = offset + seq_len(nk))
    E_list[[k]] <- Ek; Vmarg[[k]] <- Vk
    row_output  <- c(row_output, rep.int(k, nk))
    offset      <- offset + nk
  }
  T_total   <- offset
  E_stacked <- unlist(E_list, use.names = FALSE)

  if (!is.null(s$V)) {
    V <- unname(as.matrix(s$V))
    if (nrow(V) != T_total || ncol(V) != T_total)
      stop(sprintf("Study '%s': joint `V` must be %d x %d (sum of per-output times)",
                   nm, T_total, T_total), call. = FALSE)
  } else {
    V <- matrix(0, T_total, T_total)
    for (k in seq_along(blocks)) {
      if (is.null(Vmarg[[k]]))
        stop(sprintf("Study '%s.%s': observation needs its own `V` when no study-level joint `V` is given.",
                     nm, onames[k]), call. = FALSE)
      r <- blocks[[k]]$rows; V[r, r] <- Vmarg[[k]]
    }
    for (cn in names(s$cross)) {
      parts <- strsplit(cn, ":", fixed = TRUE)[[1]]
      if (length(parts) != 2L)
        stop(sprintf("Study '%s': cross name '%s' must be 'outA:outB'.", nm, cn),
             call. = FALSE)
      ia <- match(parts[1], onames); ib <- match(parts[2], onames)
      if (is.na(ia) || is.na(ib))
        stop(sprintf("Study '%s': cross pair '%s' must name two observations (%s).",
                     nm, cn, paste(onames, collapse = ", ")), call. = FALSE)
      Cab <- unname(as.matrix(s$cross[[cn]]))
      ra  <- blocks[[ia]]$rows; rb <- blocks[[ib]]$rows
      if (nrow(Cab) != length(ra) || ncol(Cab) != length(rb))
        stop(sprintf("Study '%s': cross block '%s' must be %d x %d.",
                     nm, cn, length(ra), length(rb)), call. = FALSE)
      # cross blocks are given rows=first output's times, cols=second's; reorder
      # to the sorted-time order used for the stacked vector.
      Cab <- Cab[order(as.numeric(s$observations[[ia]]$times)),
                 order(as.numeric(s$observations[[ib]]$times)), drop = FALSE]
      V[ra, rb] <- Cab; V[rb, ra] <- t(Cab)
    }
  }
  V <- (V + t(V)) / 2
  if (is.null(tryCatch(chol(V), error = function(e) NULL)))
    stop(sprintf("Study '%s': the assembled joint covariance is not positive-definite; check the cross / V blocks.",
                 nm), call. = FALSE)

  list(is_joint = TRUE, label = nm, n = s$n, ev = s$ev,
       output = blocks[[1L]]$output,   # any valid endpoint, for cmt-tagging
       times  = sort(unique(unlist(lapply(blocks, `[[`, "times")))),
       method = "cov", E = E_stacked, V = V, blocks = blocks,
       row_output = row_output, n_total = T_total)
}

# -- long-format study input ---------------------------------------------------

# Accepted column-name synonyms in a long-format study `data` frame. The output
# column names the model endpoint each row belongs to (nlmixr2 keys observations
# the same way, by DVID / CMT).
.adm_long_cols <- list(
  output = c("DVID", "dvid", "CMT", "cmt", "output"),
  time   = c("TIME", "time", "t"),
  E      = c("E", "mean", "DV", "dv"),
  var    = c("V", "var", "variance"),
  sd     = c("SD", "sd"),
  n      = c("n", "N")
)

.admLongCol <- function(df, key) {
  hit <- intersect(.adm_long_cols[[key]], names(df))
  if (length(hit) == 0L) NULL else hit[[1L]]
}

# Rewrite a long-format study into the canonical `observations` form, so the rest
# of the pipeline is untouched. The study gives one row per observed
# (endpoint, time) pair:
#
#   list(n = 60L, ev = ev,
#        data = data.frame(DVID = c("cp","cp","cb"), TIME = c(1,2,1),
#                          E = c(...), V = c(...)))     # V column = variances
#
# Same-subject (joint) studies instead carry ONE study-level covariance matrix
# whose rows/cols align with the rows of `data`; the per-row variance column is
# then unnecessary:
#
#   list(n = 60L, ev = ev, data = data.frame(DVID = ..., TIME = ..., E = ...),
#        V = Vjoint)
#
# A study-level `V` (or an explicit `joint = TRUE`) means the endpoints were
# measured on the SAME subjects and the whole stacked vector is scored by one
# MVN. Without it each endpoint is an independent likelihood block -- a separate
# experiment, so it may carry its own `n` (an `n` column) and its own dosing (a
# named list of event tables in `ev`, keyed by endpoint).
.admExpandLongStudy <- function(s, nm) {
  df <- s$data
  if (!is.data.frame(df) || nrow(df) == 0L)
    stop(sprintf("Study '%s': `data` must be a non-empty data frame", nm),
         call. = FALSE)

  c_out <- .admLongCol(df, "output"); c_t <- .admLongCol(df, "time")
  c_E   <- .admLongCol(df, "E")
  for (cc in list(list(c_out, "an endpoint column (DVID / CMT / output)"),
                  list(c_t, "a time column (TIME)"),
                  list(c_E, "a mean column (E)")))
    if (is.null(cc[[1L]]))
      stop(sprintf("Study '%s': long-format `data` needs %s", nm, cc[[2L]]),
           call. = FALSE)

  c_var <- .admLongCol(df, "var"); c_sd <- .admLongCol(df, "sd")
  c_n   <- .admLongCol(df, "n")

  out  <- as.character(df[[c_out]])
  tvec <- as.numeric(df[[c_t]])
  # Endpoints keep the order they first appear in (factor levels win if given).
  onames <- if (is.factor(df[[c_out]])) levels(droplevels(df[[c_out]])) else unique(out)
  joint  <- s$joint %||% (!is.null(s$V) || !is.null(s$cross))

  key <- paste(out, tvec, sep = "@")
  if (anyDuplicated(key))
    stop(sprintf("Study '%s': duplicate endpoint/time rows in `data` (%s)", nm,
                 paste(unique(key[duplicated(key)]), collapse = ", ")),
         call. = FALSE)

  # Row order used for the stacked E / joint V: endpoint block order, then time.
  perm <- order(match(out, onames), tvec)

  # Per-endpoint variances: from the study-level joint V's diagonal blocks when
  # given, otherwise from the variance (or SD) column.
  Vfull <- NULL
  if (!is.null(s$V)) {
    Vfull <- unname(as.matrix(s$V))
    if (nrow(Vfull) != nrow(df) || ncol(Vfull) != nrow(df))
      stop(sprintf("Study '%s': `V` must be %d x %d to match the rows of `data`",
                   nm, nrow(df), nrow(df)), call. = FALSE)
    Vfull <- Vfull[perm, perm, drop = FALSE]
  } else if (is.null(c_var) && is.null(c_sd)) {
    stop(sprintf(paste("Study '%s': long-format `data` needs a variance column",
                       "(V) or an SD column (SD), or a study-level joint `V`."),
                 nm), call. = FALSE)
  }
  vrow <- if (!is.null(c_var)) as.numeric(df[[c_var]]) else
    if (!is.null(c_sd)) as.numeric(df[[c_sd]])^2 else NULL

  # Per-endpoint dosing: one shared `ev`, or a list of event tables keyed by
  # endpoint. rxode2 event tables inherit from data.frame -- test for that, not
  # for is.list().
  ev_per <- !is.null(s$ev) && is.list(s$ev) && !is.data.frame(s$ev)
  if (ev_per && joint)
    stop(sprintf(paste("Study '%s': a joint (same-subject) study shares one `ev`",
                       "-- per-endpoint event tables describe separate experiments."),
                 nm), call. = FALSE)

  offset <- 0L
  obs <- lapply(onames, function(o) {
    idx <- which(out[perm] == o)                       # rows of this endpoint,
    j   <- perm[idx]                                   # in study/original order
    ob  <- list(output = o, times = tvec[j], E = as.numeric(df[[c_E]])[j])
    # Independent blocks each keep their own covariance: the endpoint's diagonal
    # block of a supplied V, else its per-row variances. A joint study scores the
    # whole stack with one V, so the per-endpoint copies would be redundant.
    if (!is.null(Vfull)) {
      if (!joint) ob$V <- Vfull[idx, idx, drop = FALSE]
    } else if (!is.null(vrow)) ob$V <- vrow[j]
    if (!is.null(c_n)) {
      nk <- unique(as.numeric(df[[c_n]])[j])
      if (length(nk) != 1L)
        stop(sprintf("Study '%s': endpoint '%s' has more than one `n` (%s)", nm, o,
                     paste(nk, collapse = ", ")), call. = FALSE)
      ob$n <- nk
    }
    if (ev_per) {
      if (is.null(s$ev[[o]]))
        stop(sprintf("Study '%s': `ev` is a per-endpoint list but has no entry for '%s'",
                     nm, o), call. = FALSE)
      ob$ev <- s$ev[[o]]
    }
    ob
  })
  names(obs) <- onames

  if (joint) {
    ns <- unique(vapply(obs, function(o) o$n %||% NA_real_, numeric(1)))
    ns <- ns[!is.na(ns)]
    if (length(ns) > 1L)
      stop(sprintf(paste("Study '%s': a joint (same-subject) study has one shared",
                         "`n`, but `data` gives several (%s)."),
                   nm, paste(ns, collapse = ", ")), call. = FALSE)
    if (length(ns) == 1L) s$n <- s$n %||% ns
    for (k in seq_along(obs)) obs[[k]]$n <- NULL
    if (!is.null(Vfull)) s$V <- Vfull
    # No study-level V: the joint covariance is assembled from the per-endpoint
    # marginal variances plus any `cross` blocks, as in the observations form.
  } else {
    s$V <- NULL
  }
  if (ev_per) s$ev <- NULL
  s$data <- NULL
  s$joint <- joint
  s$observations <- obs
  s
}

# Normalise one study spec into a list of observed-compartment units.
#
# Long-format input (study carries a `data` frame with one row per observed
# endpoint/time, plus an optional study-level joint `V`) is rewritten into the
# `observations` form first -- see .admExpandLongStudy().
#
# Multi-compartment forms (study carries an `observations` list):
#   * Independent blocks -- each observed output has its own n/ev/times/E/V and is
#     summed as a separate likelihood block (separate experiments / subjects; no
#     cross-compartment covariance).
#   * Joint same-subject -- outputs measured on the SAME subjects; the study gives
#     a joint covariance (study-level `V`, or per-output marginal `V` + a `cross`
#     list) and shared n/ev. Collapsed to ONE joint unit scored by a single MVN
#     with shared random effects. See .admBuildJointUnit().
#
# Legacy single-output form: the study's E/V/n/times fields describe one implicit
# observation. Top-level normalised fields (V, method, v_diag) are preserved for
# backward compatibility; `$observations` holds the single unit either way.
.admNormaliseStudy <- function(s, nm, default_output = NULL) {
  if (!is.null(s$data)) s <- .admExpandLongStudy(s, nm)
  if (!is.null(s$observations) &&
      (isTRUE(s$joint) || !is.null(s$cross) || !is.null(s$V))) {
    if (!is.list(s$observations) || length(s$observations) == 0L)
      stop(sprintf("Study '%s': `observations` must be a non-empty list", nm),
           call. = FALSE)
    if (is.null(s$n) || is.null(s$ev))
      stop(sprintf(paste("Study '%s': a joint (same-subject) study needs a shared",
                         "`n` and `ev` at the study level (measured on the same subjects)."),
                   nm), call. = FALSE)
    s$observations <- setNames(list(.admBuildJointUnit(s, nm, default_output)), nm)
    s$multi <- TRUE; s$joint <- TRUE
  } else if (!is.null(s$observations)) {
    if (!is.list(s$observations) || length(s$observations) == 0L)
      stop(sprintf("Study '%s': `observations` must be a non-empty list", nm),
           call. = FALSE)
    onames   <- names(s$observations)
    if (is.null(onames) || any(!nzchar(onames)))
      onames <- paste0("obs", seq_along(s$observations))
    defaults <- list(n = s$n, ev = s$ev, output = s$output %||% default_output)
    s$observations <- setNames(lapply(seq_along(s$observations), function(k)
      .admNormaliseObs(s$observations[[k]], paste0(nm, ".", onames[k]), defaults)),
      onames)
    s$multi <- TRUE
  } else {
    unit <- .admNormaliseObs(
      list(E = s$E, V = s$V, n = s$n, times = s$times, ev = s$ev,
           method = s$method, output = s$output %||% default_output), nm)
    # Preserve top-level normalised fields (legacy callers / tests read these).
    s$E <- unit$E; s$V <- unit$V; s$method <- unit$method
    s$v_diag <- unit$v_diag; s$output <- unit$output
    s$observations <- setNames(list(unit), nm)
    s$multi <- FALSE
  }
  s
}

# Flatten normalised studies into a single list of independent observation units.
# Each unit is self-contained (output/ev/times/E/V/n/method/v_diag/label); the
# aggregate -2LL is the sum over units, reusing the multi-study summation loop.
.admFlattenStudies <- function(studies) {
  units <- unlist(lapply(studies, `[[`, "observations"), recursive = FALSE,
                  use.names = FALSE)
  setNames(units, vapply(units, function(u) u$label, character(1)))
}

# Attach `ev_full` (dosing merged with observation times) to each unit. Defaults
# to a 100-unit bolus into compartment 1 when a unit gives no `ev`.
#
# tag_cmt: when TRUE (multi-compartment fits) each unit's observation records are
# tagged with its output compartment. nlmixr2's simulation model for a
# multi-endpoint model routes observations by compartment; untagged observations
# are ambiguous across endpoints and the solve errors. Single-output fits keep
# untagged observations (unchanged behaviour; also handles linCmt where the
# output resolves to "ipredSim", which is not a valid dosing/observation cmt).
.admBuildEvFull <- function(units, tag_cmt = FALSE) {
  lapply(units, function(u) {
    ev <- if (!is.null(u$ev)) u$ev else rxode2::et(amt = 100)
    # `ev` is documented as DOSING-only. If a user also puts observation rows in it,
    # the et() calls below append the study times a SECOND time and every point is
    # silently duplicated -- a badly wrong fit with no indication. Warn rather than
    # silently rewriting the event table: reconstructing `ev` from a filtered
    # data.frame loses event attributes rxode2 needs (it broke the sensitivity
    # solve outright), so telling the user is both safer and clearer.
    if (isTRUE(getOption("admixr2.warn.ev.obs", TRUE))) {
      .nobs <- tryCatch({
        .d <- as.data.frame(ev)
        if ("evid" %in% names(.d)) sum(.d$evid == 0) else 0L
      }, error = function(e) 0L)
      if (.nobs > 0L)
        warning("study event table `ev` contains ", .nobs, " observation record(s). ",
                "`ev` should carry DOSING only -- the study's `times` are added ",
                "separately, so those rows will be duplicated.", call. = FALSE)
    }
    # Joint units are always multi-endpoint -> always tag. A single tag (the
    # first output) is enough: the multi-endpoint solve returns every output
    # column at the observation times, and each block is extracted by name.
    u$ev_full <- if ((tag_cmt || isTRUE(u$is_joint)) && !is.null(u$output))
      ev |> rxode2::et(u$times, cmt = u$output)
    else
      ev |> rxode2::et(u$times)
    # Joint units also need a per-block event table (obs tagged with that block's
    # output cmt at its own times) so the sensitivity model can return each
    # output's prediction + sensitivities for the analytical joint gradient.
    if (isTRUE(u$is_joint))
      u$blocks <- lapply(u$blocks, function(blk) {
        blk$ev_full <- ev |> rxode2::et(blk$times, cmt = blk$output)
        blk
      })
    u
  })
}

# TRUE when the flattened units observe more than one distinct output variable
# (i.e. a genuine multi-compartment fit needing per-output residual error and
# the FD gradient path rather than the single-output analytical path).
.admIsMultiOutput <- function(units, default_output) {
  outs <- vapply(units, function(u) u$output %||% default_output %||% NA_character_,
                 character(1))
  length(unique(outs[!is.na(outs)])) > 1L
}

# Assemble fullTheta vector in iniDf row order (thetas + sigma SDs + omega entries).
.admFullTheta <- function(pars, pinfo) {
  .ini <- pinfo$iniDf
  setNames(vapply(seq_len(nrow(.ini)), function(i) {
    nm <- .ini$name[i]
    if (.ini$fix[i])                      return(.ini$est[i])
    if (nm %in% names(pars$struct))       return(unname(pars$struct[nm]))
    # Residual params report on their iniDf scale: variances as an SD; a pow()
    # exponent and a t() degrees-of-freedom as themselves (.admSigmaNat has
    # already mapped log(nu - 2) back to nu).
    if (nm %in% names(pars$sigma_var)) {
      .k <- match(nm, pinfo$sigma_names)
      .v <- unname(pars$sigma_var[nm])
      return(if (.admSigmaRole(pinfo)[.k] == "var") sqrt(.v) else .v)
    }
    if (!is.na(.ini$neta1[i]))
      return(pars$omega[.ini$neta1[i], .ini$neta2[i]])
    .ini$est[i]
  }, double(1)), .ini$name)
}

# Compute objective function statistics for a completed fit.
# nobs = sum(n_subjects * n_times) across studies, matching nlmixr2's individual-level convention.
.admCalcObjStats <- function(objective, npar, studies) {
  nobs <- sum(vapply(studies, function(s)
    as.integer(s$n) * (s$n_total %||% length(s$times)), integer(1)))
  ll   <- -objective / 2
  attr(ll, "df")   <- npar
  attr(ll, "nobs") <- nobs
  class(ll) <- "logLik"
  objDf <- data.frame(
    OBJF             = objective,
    AIC              = objective + 2 * npar,
    BIC              = objective + log(nobs) * npar,
    "Log-likelihood" = as.numeric(ll),
    check.names      = FALSE
  )
  list(ll = ll, nobs = nobs, npar = npar, objDf = objDf)
}

# Bridge admControl/adirmcControl fields into foceiControl for nlmixr2 table machinery.
#
# `skip_cov` says which thetas nlmixr2est should NOT expect a standard error for.
# It fills its SE column by walking the thetas in iniDf order and taking the next
# entry of sqrt(diag(cov)) for each one it is not skipping, so this vector and the
# row order of the covariance are two halves of the same contract -- see
# .admCovSkip()/.admCovThetaOrder(). Left to nlmixr2est's own default when NULL.
.admToFoceiControl <- function(ctl, skip_cov = NULL) {
  .args <- list(
    rxControl          = ctl$rxControl,
    maxOuterIterations = 0L,
    maxInnerIterations = 0L,
    covMethod          = 0L,
    sumProd            = ctl$sumProd,
    optExpression      = ctl$optExpression,
    literalFix         = ctl$literalFix,
    scaleTo            = 0,
    calcTables         = ctl$calcTables,
    addProp            = ctl$addProp,
    interaction        = 0L,
    compress           = ctl$compress,
    ci                 = ctl$ci,
    sigdigTable        = ctl$sigdigTable)
  if (!is.null(skip_cov)) .args$skipCov <- skip_cov
  do.call(nlmixr2est::foceiControl, .args)
}

# LHS: one sample per stratum per dimension, independently permuted.
.lhsSample <- function(n, d) {
  m <- matrix(0, nrow = n, ncol = d)
  for (j in seq_len(d))
    m[, j] <- (sample.int(n) - 1L + runif(n)) / n
  m
}

# Generate z matrices -- one per study (seed must be set by caller).
# sampling: "sobol" (quasi-random), "lhs" (Latin hypercube), "rnorm" (iid normal).
.admMakeZ <- function(n_sim, pinfo, n_studies, sampling = "sobol") {
  replicate(n_studies, {
    if (pinfo$n_eta == 0L)
      return(matrix(0, nrow = n_sim, ncol = 1L))
    switch(sampling,
      sobol  = {
        z <- qnorm(randtoolbox::sobol(n = n_sim, dim = pinfo$n_eta))
        if (!is.matrix(z)) z <- matrix(z, ncol = 1L)
        z
      },
      halton = {
        z <- qnorm(randtoolbox::halton(n = n_sim, dim = pinfo$n_eta))
        if (!is.matrix(z)) z <- matrix(z, ncol = 1L)
        z
      },
      torus  = {
        z <- qnorm(randtoolbox::torus(n = n_sim, dim = pinfo$n_eta))
        if (!is.matrix(z)) z <- matrix(z, ncol = 1L)
        z
      },
      lhs    = qnorm(.lhsSample(n_sim, pinfo$n_eta)),
      rnorm  = matrix(rnorm(n_sim * pinfo$n_eta), nrow = n_sim),
      stop("admMakeZ: unknown sampling method '", sampling, "'", call. = FALSE)
    )
  }, simplify = FALSE)
}

# Pre-allocate params matrix list -- one per study.
# Matrix avoids data.frame list COW overhead: first col-write copies once,
# subsequent writes modify in-place. as.data.frame() wraps at rxSolve call site.
#
# Includes one residual-error placeholder column per observed output
# (rxerr.<output>). rxSolve requires every endpoint's rxerr parameter to be
# present in the params frame (its value is immaterial -- admixr2 adds residual
# error analytically and reads the structural prediction). There is exactly one
# rxerr per prediction line regardless of the error model (add/prop/lnorm or a
# combined add+prop), so the endpoint set is `unique(sigma_output)`; the mapping
# holds for named endpoints and linCmt (`rxerr.rxLinCmt`). All other model
# parameters (CMT, hard-coded constants) are left for rxSolve to fill from the
# model's own defaults -- do not add them here.
.admMakeParamsList <- function(n_sim, pinfo, n_studies = 1L) {
  so    <- unique(pinfo$sigma_output[!is.na(pinfo$sigma_output)])
  rxerr <- if (length(so)) paste0("rxerr.", so) else "rxerr.cp"
  col_names <- c(pinfo$struct_names, pinfo$eta_col_names,
                 pinfo$sigma_names, rxerr)
  replicate(n_studies, {
    m <- matrix(0, nrow = n_sim, ncol = length(col_names),
                dimnames = list(NULL, col_names))
    m[, rxerr] <- 1L
    m
  }, simplify = FALSE)
}

# -- FOCEI-style aligned progress output ---------------------------------------

# Column names: -2LL, back-transformed struct thetas, sigma SDs, omega diagonal variances.
.admProgressNames <- function(pinfo) {
  omega_diag_nms <- if (pinfo$n_eta > 0L)
    pinfo$eta_names[pinfo$chol_i[pinfo$chol_diag]]
  else character(0)
  c("-2LL", pinfo$struct_names, pinfo$sigma_names, omega_diag_nms)
}

# Bordered header block (separator + header row + separator). iter_w sets the
# label column width; data columns are max(name_width, 10). bottom=FALSE omits
# the closing separator so a phase divider can follow immediately.
.admProgressHeader <- function(pinfo, iter_w = 8L, bottom = TRUE) {
  nms <- .admProgressNames(pinfo)
  cws <- pmax(nchar(nms), 8L)
  sep <- paste0("+", strrep("-", iter_w + 2L), "+",
                paste(vapply(cws, function(w) strrep("-", w + 2L), character(1)),
                      collapse = "+"), "+")
  hdr <- paste0("| ", formatC("", width = iter_w), " | ",
                paste(mapply(formatC, nms, width = cws), collapse = " | "), " |")
  if (bottom) paste0(sep, "\n", hdr, "\n", sep) else paste0(sep, "\n", hdr)
}

# Full-width section divider matching .admProgressHeader outer width; pads label right with dashes.
.admProgressDivider <- function(label, pinfo, iter_w = 8L) {
  nms     <- .admProgressNames(pinfo)
  cws     <- pmax(nchar(nms), 8L)
  inner_w <- iter_w + 2L + sum(cws) + 3L * length(cws)
  n_dash  <- max(inner_w - 2L - nchar(label), 0L)
  paste0("+--", label, strrep("-", n_dash), "+")
}

.admProgressPhase <- function(phase_idx, phase_name, ph_step, pinfo, iter_w = 8L)
  .admProgressDivider(sprintf(" Phase %d: %s (+/-%.2f) ", phase_idx, phase_name, ph_step), pinfo, iter_w)

.admProgressRestart <- function(r, n_r, pinfo, iter_w = 8L)
  .admProgressDivider(sprintf(" Restart %d / %d ", r, n_r), pinfo, iter_w)

# One bordered data row aligned to the same column widths as .admProgressHeader.
.admProgressRow <- function(iter_label, nll, p, pinfo, iter_w = 8L) {
  pars <- tryCatch(.admUnpack(p, pinfo), error = function(e) NULL)
  if (is.null(pars)) return(NULL)
  nms <- .admProgressNames(pinfo)
  cws <- pmax(nchar(nms), 8L)
  struct_vals <- vapply(pinfo$struct_names, function(nm)
    .admBackTransform(pars$struct[[nm]], pinfo$struct_transforms[[nm]]), double(1))
  sigma_vals  <- sqrt(pars$sigma_var)
  omega_vals  <- if (pinfo$n_eta > 0L)
    diag(pars$omega)[pinfo$chol_i[pinfo$chol_diag]]
  else numeric(0)
  nll_str  <- local({
    s <- formatC(nll, format = "f", digits = 2)
    if (nchar(s) <= cws[1L]) formatC(nll, format = "f", digits = 2, width = cws[1L])
    else                     formatC(nll, format = "e", digits = 1, width = cws[1L])
  })
  par_strs <- mapply(function(x, w) formatC(x, format = "g", digits = 4, width = w),
                     c(struct_vals, sigma_vals, omega_vals), cws[-1L])
  val_strs <- c(nll_str, par_strs)
  paste0("| ", formatC(iter_label, width = iter_w, flag = "-"), " | ",
         paste(val_strs, collapse = " | "), " |")
}

# Timing row: label column shows elapsed time, all data columns blank.
.admProgressTimingRow <- function(sec, pinfo, iter_w = 8L) {
  nms    <- .admProgressNames(pinfo)
  cws    <- pmax(nchar(nms), 8L)
  blanks <- vapply(cws, function(w) formatC("", width = w), character(1))
  time_label <- if (sec >= 60) sprintf("%.1f min", sec / 60) else sprintf("%.1f sec", sec)
  paste0("| ", formatC(time_label, width = iter_w, flag = "-"), " | ",
         paste(blanks, collapse = " | "), " |")
}
