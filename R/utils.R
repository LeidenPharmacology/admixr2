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

# Detect output variable name from ui$predDf (default "cp").
.admOutputVar <- function(ui) {
  var <- tryCatch(
    { pd <- ui$predDf; if (!is.null(pd) && "var" %in% names(pd)) pd$var[1] else "cp" },
    error = function(e) "cp")
  # Internal nlmixr2 linCmt names (rxLinCmt, linCmtB, ...) don't appear in
  # simulation model rxSolve output — use ipredSim which is always present.
  if (startsWith(var, "rx") || startsWith(var, "linCmt")) "ipredSim" else var
}

# Normalise one study spec: coerce V to matrix, auto-detect diagonal, set method.
# V as vector -> treated as variances -> expand to diag matrix, force "var".
# V as matrix with all off-diagonal zeros -> force "var" unless user said "cov".
# V must use ML denominator n (not n-1); use cov.wt(dv_mat, method="ML")$cov.
.admNormaliseStudy <- function(s, nm) {
  for (f in c("n", "E", "V", "times"))
    if (is.null(s[[f]])) stop(sprintf("Study '%s' missing '%s'", nm, f), call. = FALSE)
  s$E <- as.numeric(s$E)
  if (is.vector(s$V) && !is.list(s$V)) {
    if (identical(s$method, "cov"))
      warning(sprintf("Study '%s': V is a vector (variances only) but method='cov' requested -- using method='var'", nm), call. = FALSE)
    s$V      <- diag(as.numeric(s$V))
    s$method <- "var"
  } else {
    s$V     <- unname(as.matrix(s$V))
    is_diag <- all(s$V[lower.tri(s$V)] == 0) && all(s$V[upper.tri(s$V)] == 0)
    s$method <- if (is_diag && is.null(s$method)) "var" else
      match.arg(s$method %||% "cov", c("cov", "var"))
    if (!is_diag && s$method == "var")
      warning(sprintf("Study '%s': V has non-zero off-diagonal entries but method='var' -- off-diagonal entries will be ignored", nm), call. = FALSE)
  }
  if (identical(s$method, "var")) s$v_diag <- diag(s$V)
  s
}

# Assemble fullTheta vector in iniDf row order (thetas + sigma SDs + omega entries).
.admFullTheta <- function(pars, pinfo) {
  .ini <- pinfo$iniDf
  setNames(vapply(seq_len(nrow(.ini)), function(i) {
    nm <- .ini$name[i]
    if (.ini$fix[i])                      return(.ini$est[i])
    if (nm %in% names(pars$struct))       return(unname(pars$struct[nm]))
    if (nm %in% names(pars$sigma_var))    return(sqrt(unname(pars$sigma_var[nm])))
    if (!is.na(.ini$neta1[i]))
      return(pars$omega[.ini$neta1[i], .ini$neta2[i]])
    .ini$est[i]
  }, double(1)), .ini$name)
}

# Compute objective function statistics for a completed fit.
# nobs = sum(n_subjects * n_times) across studies, matching nlmixr2's individual-level convention.
.admCalcObjStats <- function(objective, npar, studies) {
  nobs <- sum(vapply(studies, function(s) as.integer(s$n) * length(s$times), integer(1)))
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
.admToFoceiControl <- function(ctl) {
  nlmixr2est::foceiControl(
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
.admMakeParamsList <- function(n_sim, pinfo, n_studies = 1L) {
  col_names <- c(pinfo$struct_names, pinfo$eta_col_names,
                 pinfo$sigma_names, "rxerr.cp")
  replicate(n_studies, {
    m <- matrix(0, nrow = n_sim, ncol = length(col_names),
                dimnames = list(NULL, col_names))
    m[, "rxerr.cp"] <- 1L
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
