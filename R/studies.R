# Study and observation-unit handling: resolving a model's endpoints, normalising
# a study specification, flattening it to independent observation units (or one
# joint same-subject unit), and the per-row maps the estimators read off a unit.
#
# Split out of utils.R; contents unchanged (see R/covreport.R for why that is
# safe). This is the largest single concern that file held, and the one with the
# most invariants: `multi_out` is MODEL-level while `is_joint` is study-level, a
# joint unit routes per ROW (row_output) and so has no endpoint of its own --
# though it still carries one copied off blocks[[1]] for cmt-tagging, and its
# BLOCKS each have a real one -- and normalising twice must be idempotent WITHOUT
# being inert: a second pass still has to fill an output the first pass had no
# default for, in the blocks as well as in the plain units.


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

# The ENDPOINT names, as nlmixr2 knows them -- predDf$var verbatim.
#
# NOT the same thing as .admOutputVars(), and the difference matters for exactly
# the endpoints that made .admEndpointVar() necessary: `y ~ pois(lam)` is SOLVED
# through `lam` but nlmixr2 knows the endpoint as `y`. These names go in the DVID
# column of the dummy frame handed to nlmixr2CreateOutputFromUi(), and its
# dvid->cmt translation rejects a name that is not an endpoint -- so passing the
# solve variable there made a converged multi-endpoint count fit die at the
# output-building step with "'dvid'->'cmt' ... on a undefined compartment".
.admEndpointNames <- function(ui) {
  nms <- tryCatch(as.character(ui$predDf$var), error = function(e) NULL)
  if (is.null(nms) || !length(nms)) return(.admOutputVars(ui))
  unique(nms)
}

# A count or beta endpoint cannot share a model with other endpoints.
#
# Multi-endpoint solves route observations by COMPARTMENT: .admBuildEvFull() tags
# each unit's records with `cmt = unit$output`, and rxode2 resolves that against
# the model's endpoints. A count endpoint's output is its distribution's ARGUMENT
# (`y ~ pois(lam)` is read through `lam`), which is an ordinary model variable and
# not an endpoint at all, so the tagged records match nothing: the solve returns no
# rows for that unit and the objective silently comes back Inf -- there is no
# wrong-but-plausible number, but there is also nothing telling the user why.
#
# Single-endpoint count/beta models are unaffected (no tagging happens) and are
# what the count/beta support was built for. An ordinal endpoint is ONE predDf row
# whose categories are separate outputs, so it is not "mixed" either.
.admCheckMixedEndpoints <- function(ui) {
  pd <- tryCatch(ui$predDf, error = function(e) NULL)
  if (is.null(pd) || nrow(pd) < 2L || !"distribution" %in% names(pd))
    return(invisible(NULL))
  d <- as.character(pd$distribution)
  w <- which(d %in% c("pois", "dpois", "binom", "dbinom", "nbinomMu", "dnbinomMu",
                      "beta", "dbeta"))
  if (!length(w)) return(invisible(NULL))
  stop("A ", d[w[1L]], "() endpoint cannot be combined with other endpoints in ",
       "one model.\n",
       "  Endpoint ", sQuote(as.character(pd$var[w[1L]])), " is read through its ",
       "distribution's argument,\n",
       "  which is a model variable rather than a compartment -- so the ",
       "multi-endpoint solve,\n",
       "  which routes observations by compartment, cannot deliver its ",
       "observations.\n",
       "  Fit that endpoint in a model of its own.", call. = FALSE)
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
  m <- .admResidMoments(mu_struct, diag(V_pred), arr, V_pred, rt)
  list(mu = m$mu, V = m$V)
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
  # Covariates describe the study's SUBJECTS, so every observed output of that
  # study inherits them: `cov` is the value written into the solve, `cov_dist`
  # the distribution those subjects span (which drives the Omega collapse).
  ob[["cov"]]      <- ob[["cov"]]      %||% defaults[["cov"]]
  ob[["cov_dist"]] <- ob[["cov_dist"]] %||% defaults[["cov_dist"]]
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
       cov = s[["cov"]], cov_dist = s[["cov_dist"]],
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
# Convert a study's reported covariance to the ML (denominator n) convention the
# likelihood requires.
#
# The two input types admixr2 serves disagree about what `V` IS, and until now
# the difference was a footnote the user had to act on:
#
#   a digitised figure  ->  SD is the UNBIASED (n-1) sample SD, so V = SD^2 is
#                           an (n-1) covariance
#   datagen / own data  ->  cov.wt(method = "ML"), an n covariance
#
# Eq. (1) is the exact log-likelihood of n iid draws only for the ML form, so a
# published SD is strictly V = SD^2 * (n-1)/n. At n = 60 that is 1.7% and was
# reasonably ignored. It stops being ignorable the moment the summary is scored
# against its own sampling law: the same factor reappears there as the alignment
# of tau with E[t], where getting it wrong is measurably WORSE than not
# correcting at all.
#
# So it becomes a declaration rather than a convention, PER STUDY -- a
# meta-analysis routinely mixes a digitised figure with a model-derived source,
# and the two do not share a denominator.
#
# Idempotent: `v_denom` is set to "ml" once applied, so normalising twice cannot
# apply it twice.
.admVDenom <- function(s, nm) {
  vd <- s[["v_denom"]] %||% "ml"
  if (!is.character(vd) || length(vd) != 1L || !vd %in% c("ml", "unbiased"))
    stop(sprintf("Study '%s': `v_denom` must be \"ml\" or \"unbiased\"", nm),
         call. = FALSE)
  if (identical(vd, "ml")) return(s)
  conv <- function(V, n, what) {
    if (is.null(V)) return(NULL)
    if (is.null(n)) stop(sprintf(
      "Study '%s': `v_denom = \"unbiased\"` needs `n` to convert %s to the ML denominator",
      nm, what), call. = FALSE)
    n <- as.numeric(n)[[1L]]
    if (!is.finite(n) || n <= 1)
      stop(sprintf("Study '%s': `v_denom = \"unbiased\"` needs n > 1 (got %s)",
                   nm, format(n)), call. = FALSE)
    V * (n - 1) / n
  }
  s$V <- conv(s$V, s$n, "V")
  if (!is.null(s$observations))
    s$observations <- lapply(s$observations, function(ob) {
      ob$V <- conv(ob$V, ob$n %||% s$n, "an observation's V"); ob
    })
  if (!is.null(s$cross))
    s$cross <- lapply(s$cross, function(x) conv(x, s$n, "a cross block"))
  s$v_denom <- "ml"
  s
}

.admNormaliseStudy <- function(s, nm, default_output = NULL) {
  # IDEMPOTENT, and it has to be stated rather than assumed.
  #
  # Normalising a legacy single-output study ADDS an `observations` list while
  # KEEPING its top-level `V` -- which is precisely the signature the joint
  # (same-subject) branch below tests for. So a second pass over an
  # already-normalised study silently collapsed it into ONE JOINT unit:
  #
  #   u <- .admFlattenStudies(list(.admNormaliseStudy(raw,  "s")))  # is_joint FALSE
  #   u <- .admFlattenStudies(list(.admNormaliseStudy(once, "s")))  # is_joint TRUE
  #
  # No error, no warning, a plausible fit -- down a different likelihood path,
  # and for adfo with `have_d2` forced FALSE (it requires `!any_joint`), so the
  # order-2 analytical struct-theta gradient quietly turns itself off.
  #
  # Each driver normalises exactly once, so this was not reachable from a normal
  # fit. It WAS reachable from the test fixtures, which hand out pre-normalised
  # studies that then get normalised again by the driver -- which is how it was
  # found. Guarding here rather than in the fixtures because "normalise a study"
  # should not be an operation you can only safely perform once.
  # ... but idempotent is not the same as INERT. The first pass may have run
  # without a `default_output` (nothing but the driver knows the model's endpoint,
  # and the fixtures normalise before there is a model), which leaves every unit
  # with output = NULL. Short-circuiting outright made the driver's later pass --
  # the one that DOES carry output_var -- a no-op, so the NULL was permanent: for a
  # multi-endpoint model .admBuildEvFull(tag_cmt = TRUE) then has nothing to tag
  # `cmt` with and the unit silently reads the wrong compartment's trajectory.
  # Repeating the normalisation was masking that; so fill what is still missing,
  # and only then return.
  if (isTRUE(s$.adm_normalised)) {
    if (!is.null(default_output)) {
      if (is.null(s$output)) s$output <- default_output
      s$observations <- lapply(s$observations, function(u) {
        if (!isTRUE(u$is_joint)) {
          if (is.null(u$output)) u$output <- default_output
          return(u)
        }
        # A JOINT unit routes per ROW, so it carries no endpoint of its own --
        # except the one .admBuildJointUnit() copies off blocks[[1]] purely for
        # cmt-tagging. Skipping joint units ENTIRELY here was too strong: their
        # BLOCKS each do have an output, taken from `ob$output %||% default_output`
        # at construction, so a study normalised before the model was known (the
        # fixtures do exactly this) leaves every blk$output NULL -- and nothing
        # later fills it, because this short-circuit is the only second pass.
        #
        # .admBuildEvFull() then runs `et(blk$times, cmt = blk$output)` per block
        # with cmt = NULL, so the joint sens solve either errors out of
        # .admSimulateJointSens() -- dropping the fit to FD -- or reads an
        # untagged compartment, giving a finite but wrong joint objective with no
        # warning.
        #
        # `row_output` needs nothing: it holds block INDICES, not names.
        if (!is.null(u$blocks))
          u$blocks <- lapply(u$blocks, function(blk) {
            if (is.null(blk$output)) blk$output <- default_output
            blk
          })
        # Restore the tagging endpoint the constructor would have set, rather
        # than stamping default_output over a stack that names its own endpoints.
        if (is.null(u$output) && length(u$blocks))
          u$output <- u$blocks[[1L]]$output
        u
      })
    }
    return(.admStampStudy(s, nm))
  }
  # BEFORE any branch reads V: the joint constructor assembles its own matrix
  # from the raw per-observation blocks and never passes through
  # .admNormaliseObs, so converting there would miss it.
  s <- .admVDenom(s, nm)
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
    defaults <- list(n = s$n, ev = s$ev, output = s$output %||% default_output,
                     cov = s[["cov"]], cov_dist = s[["cov_dist"]])
    s$observations <- setNames(lapply(seq_along(s$observations), function(k)
      .admNormaliseObs(s$observations[[k]], paste0(nm, ".", onames[k]), defaults)),
      onames)
    s$multi <- TRUE
  } else {
    unit <- .admNormaliseObs(
      list(E = s$E, V = s$V, n = s$n, times = s$times, ev = s$ev,
           method = s$method, output = s$output %||% default_output,
           cov = s[["cov"]], cov_dist = s[["cov_dist"]]), nm)
    # Preserve top-level normalised fields (legacy callers / tests read these).
    s$E <- unit$E; s$V <- unit$V; s$method <- unit$method
    s$v_diag <- unit$v_diag; s$output <- unit$output
    s$observations <- setNames(list(unit), nm)
    s$multi <- FALSE
  }
  s <- .admStampStudy(s, nm)
  s$.adm_normalised <- TRUE
  s
}

# Stamp every observation unit with the STUDY it came from.
#
# `label` identifies a unit; `study` identifies the trial. For a single-output
# study they coincide, and for a multi-output or joint study several units share
# one `study`. Nothing today reads it, and it is here because the things that
# will read it need it to have been recorded from the start:
#
#   * WITHIN- vs BETWEEN-study covariate effects. A covariate coefficient
#     estimated across studies is confounded with everything else that differs
#     between them -- measured on a confounded three-trial design, a true 0.75
#     comes back at 1.14. Separating a within-study from a between-study
#     coefficient, or giving each study its own baseline so the between-study
#     differences stop being billed to the covariate, both need to know which
#     units belong to the same trial. Units alone cannot say: covariate STRATA of
#     one trial are separate units that must share a baseline, and reporting them
#     as separate studies is exactly the mistake (measured: strata alone recover
#     almost nothing, 1.11 against the same 0.75, because the model has no
#     parameter for "trial" and the slope absorbs it regardless).
#   * A random study effect (tau^2) on the baseline, which is the same grouping.
#   * Joint individual + aggregate fits, where an aggregate study and an
#     individual dataset from the same trial share study-level parameters.
#
# Recording it costs nothing and cannot be reconstructed later: once units are
# flattened, `label` is all that survives, and a label is deliberately allowed to
# be anything. Callers that build units by hand (fixtures) may not set it, so
# .admUnitStudy() falls back to the label rather than failing -- one unit per
# study is the right reading of a unit that never declared otherwise.
.admStampStudy <- function(s, nm) {
  if (is.null(s$observations)) return(s)
  s$observations <- lapply(s$observations, function(u) {
    if (is.null(u$study)) u$study <- nm
    u
  })
  s
}

# The study a unit belongs to. See .admStampStudy() for why this is not `label`.
.admUnitStudy <- function(u) u$study %||% u$label

# Units grouped by study, in first-appearance order. The grouping every
# study-level parameter will index off.
.admStudyGroups <- function(units) {
  k <- vapply(units, .admUnitStudy, character(1))
  split(seq_along(units), factor(k, levels = unique(k)))
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
