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

# Map an rxSolve output column name to the one rxSolve actually returns.
# Internal nlmixr2 linCmt names (rxLinCmt, linCmtB, ...) don't appear in the
# simulation model rxSolve output -- use ipredSim which is always present.
.admOutputColName <- function(var)
  if (startsWith(var, "rx") || startsWith(var, "linCmt")) "ipredSim" else var

# Detect the primary/default output variable name from ui$predDf (default "cp").
# Used as the fallback output for studies/observations that don't name one.
.admOutputVar <- function(ui) {
  var <- tryCatch(
    { pd <- ui$predDf; if (!is.null(pd) && "var" %in% names(pd)) pd$var[1] else "cp" },
    error = function(e) "cp")
  .admOutputColName(var)
}

# All observable output variable names from ui$predDf (one per model endpoint).
# A multi-endpoint model (e.g. `cp ~ ...; cCSF ~ ...`) has several predDf rows.
.admOutputVars <- function(ui) {
  vars <- tryCatch(
    { pd <- ui$predDf; if (!is.null(pd) && "var" %in% names(pd)) as.character(pd$var) else "cp" },
    error = function(e) "cp")
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
  ap  <- .admResidApply(mu_struct, diag(V_pred), arr)
  diag(V_pred) <- ap$dv
  list(mu = ap$mu, V = V_pred)
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
    # Residual params report on their iniDf scale: variances as an SD, a pow()
    # exponent as itself.
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
