# =============================================================================
# Covariate marginalisation for aggregate data modelling
# =============================================================================
#
# A study's subjects span a distribution of covariate values (`cov_dist`), and
# the aggregate (E, V) it reports is MARGINAL over it, so the prediction must be
# too:
#
#   mu = E_{a,eta}[f]        V = Cov_{a,eta}(f) + residual
#
# Covariate-induced between-subject variability belongs INSIDE V_pred, because
# that is where the observed V carries it. Two paths compute those moments --
# the general per-row/product-grid path ("rows") and the shift path ("shift") --
# chosen per study by .admCheckCovariates(). adgh's `cov_integration` selects
# between the product grid, its second-order Taylor surrogate and the shift.
#
# THE RULE, of which everything below is a consequence: a block's PREDICTION
# must integrate over the same covariate distribution its DATA were aggregated
# over. Pooled data need a pooled prediction; per-stratum data need the
# prediction marginalised within that stratum. Both mismatches cost, in opposite
# directions, and neither announces itself.
#
# THREE EARLIER PATHS WERE REMOVED. Do not resurrect them:
#
#  * "collapse" folded a normal covariate's variance into Omega in closed form.
#    It is the GAUSSIAN SPECIAL CASE OF THE SHIFT -- for a linear Delta and a
#    normal covariate, u = Delta(a) + eta is normal with variance
#    Omega + beta^2 Var(a) -- verified at 6.7e-16 on the mean and 5.0e-13 on the
#    covariance with the ODE solver taken out of the comparison. It also needed
#    grad = "none", and every estimator defaults to a gradient, so it was
#    unreachable in a real fit.
#  * "uq" replaced the affected eta column with quantiles of u = Delta(a) + eta.
#    Right idea, and it is what "shift" does -- but it INFERRED the property
#    from the model TEXT, with four measured silent-wrong-answer modes: a
#    discrete or multi-modal covariate (13-20% on the mean, 3-5x on the
#    covariance), a non-log link or additive effect (16-21% on the covariance),
#    two covariates on one eta (result bit-identical to the second alone), and a
#    user `quantile` integrated by 32 equal-weight midpoints (coefficient ~7%
#    high). .admShiftVerify() now checks numerically instead.
#  * NODE METHODS (gl / gh / taylor) scored the study at fixed covariate values
#    and combined the per-node -2LL linearly. Given per-node DATA that is exact
#    -- ordinary multi-study fitting with quadrature weights standing in for
#    stratum sizes, recovering a true 0.75 as 0.7500. Given the ONE pooled
#    (E, V) a publication reports it is not a likelihood of anything: the sum is
#    -2 log of an unnormalised geometric mean of densities, and its mean term
#    decomposes as
#
#      sum_k w_k r_k' Vi r_k  =  rbar' Vi rbar + tr(Vi %*% Cov_a(mu))
#
#    whose second term is a covariate-spread penalty carrying no data -- at the
#    true parameters the first is 0 and the second is the whole objective. It
#    pulls 0.75 to 0.3045 (-59%), or 0.4567 (-39%) even with V_pred supplied
#    marginally, so it is the log-density pooling and not the variance handling.
#    A displaced TARGET, not estimator bias: the construction is exact the
#    moment it is given the data it describes.
#
# So fit stratum summaries as ordinary studies, each with its own `n` and its
# own `cov_dist` -- NOT `cov`. A publication reporting strata also reports their
# real n_k, which beats quadrature weights.
#
# `cov` is a POINT value, and plugging in the stratum mean is the ecological
# plug-in: the aggregate relation equals the individual one only if the model is
# affine over that stratum's support, or the stratum is degenerate. Against a
# true 0.75, marginalising within the stratum recovers 0.7500 at every K, while
# plugging in the mean gives 0.875 at K = 2 (+17%), 0.782 at K = 4 and 0.757 at
# K = 8 -- dying only as K grows, and biasing UPWARD, the mirror of the node
# route's downward pull. A subgroup table reports the stratum's own mean and SD,
# which is exactly the truncated `cov_dist` this needs.
#
# (`%||%` is defined in utils.R.)

# Append this study's covariate columns to an rxSolve params frame (matrix or
# data.frame). Returns `mat` unchanged when there is nothing to add.
#
# ONLY names present in `cov_s` are added -- deliberately NOT
# setdiff(mod$params, colnames(mat)) with a 0 default. That blanket zero-fill is
# the pattern this package has been bitten by twice: it clobbered genuine
# hard-coded model constants (qout/vb -> /0 -> an NA objective), and it handed
# the solve lambda = 0 for an ESTIMATED boxCox/yeoJohnson, making the sens
# gradient ~60x wrong for boxCox and NaN for yeoJohnson. Everything that is not
# a covariate is rxSolve's to default from the model itself.
#
# Named .adm* on purpose: the dev-mode daemon payload is collected by a
# /^\.(adm|adfo|adgh|adirmc)/i regex, so a helper named .cov_fills would be
# missing inside every parallel-restart worker.
.admCovCols <- function(mat, mod_params, cov_s, cov_rows = NULL) {
  # PER-ROW values (general path): each simulated subject carries its own
  # covariate, so rxode2 evaluates whatever functional form the model contains.
  if (!is.null(cov_rows)) {
    nms <- setdiff(intersect(colnames(cov_rows), mod_params), colnames(mat))
    if (nrow(cov_rows) != nrow(mat))
      stop(".admCovCols: cov_rows has ", nrow(cov_rows), " rows but the params ",
           "frame has ", nrow(mat), ". Recycling here would hand subjects the ",
           "wrong covariate values silently.", call. = FALSE)
    if (length(nms)) {
      add <- as.matrix(cov_rows[, nms, drop = FALSE])
      mat <- if (is.data.frame(mat))
        cbind(mat, as.data.frame(add, check.names = FALSE)) else cbind(mat, add)
    }
    # ... and then FALL THROUGH to the fixed values. A study may declare a
    # distribution for one covariate and a constant for another; returning here
    # dropped every constant one, so a model reading both could not solve at all
    # ("parameter(s) are required for solving: SEX"), swallowed by .admNLL's
    # tryCatch into an Inf objective everywhere.
    cov_s <- cov_s[setdiff(names(cov_s), colnames(cov_rows))]
  }
  if (is.null(cov_s) || length(cov_s) == 0L) return(mat)
  nms <- setdiff(intersect(names(cov_s), mod_params), colnames(mat))
  if (length(nms) == 0L) return(mat)
  add <- matrix(rep(as.numeric(unlist(cov_s[nms], use.names = FALSE)),
                    each = nrow(mat)),
                nrow(mat), length(nms), dimnames = list(NULL, nms))
  if (is.data.frame(mat)) cbind(mat, as.data.frame(add, check.names = FALSE))
  else                    cbind(mat, add)
}

# Covariate columns for a params frame that stacks `n_blk` blocks of `n_sim`
# rows, each block a PERTURBATION OF THE SAME SUBJECTS (the finite-difference
# frames in .admGrad). The covariate rows are tiled per block so every block sees
# the same subjects' covariates -- which is exactly what makes the difference a
# common-random-numbers one. Tiling with the wrong stride would give each
# perturbation different subjects and quietly turn the gradient into noise.
.admCovColsTiled <- function(mat, mod_params, s, n_sim, n_blk) {
  cr <- s[["cov_rows"]]
  if (!is.null(cr)) cr <- cr[rep(seq_len(n_sim), times = n_blk), , drop = FALSE]
  .admCovCols(mat, mod_params, s[["cov"]], cr)
}

# Attach per-row covariate values to a study for the GENERAL path. Returns the
# study unchanged on the shift path (where the covariate is held at its
# reference and its whole effect rides in the shifted eta column) and when no
# distribution is declared.
.admStudyCovRows <- function(s, pinfo, n_row) {
  if (!identical(s$.adm_cov_path, "rows")) return(s)
  # THE FULL PRODUCT DRAW, ALWAYS -- deliberately, and not for want of a
  # cheaper one.
  #
  # Sampling in the collapsed subspace is an exact change of variables and was
  # tried. Two things killed it. It buys nothing measurable: with randomised
  # QMC (Cranley-Patterson, 24 replicates, paired within replicate) the
  # covariance error ratio is 1.12x [0.89,1.84] at n = 1000, 1.24x [0.90,1.53]
  # at 4000 and 1.47x [0.78,2.17] at 16000 -- interquartile range spanning 1.0
  # in every cell. (An earlier "~2x" came from medians across randtoolbox
  # `seed` values, and `seed` is inert -- the sequences are bit-identical.)
  #
  # And it breaks COMMON RANDOM NUMBERS. The rotation depends on the covariate
  # coefficients, which are estimated, so re-aiming it -- which correctness
  # requires -- makes the DRAWS move with the parameters. Measured, a step of
  # 1e-6 in one coefficient shifted the sampled covariate values by 1.1e-4, so
  # a CRN-FD difference in that direction carries a design change on top of the
  # parameter change. .admCovRowsFor is deterministic in `cov_dist` alone, which
  # is data, and that is the property the gradient depends on.
  #
  # The quadrature estimators are unaffected: their gain is in DESIGN POINTS,
  # they re-aim per objective call, and they have no random numbers to hold
  # common.
  # ...AND THE SAME DRAW EVERY TIME, WHICH IS WORTH CACHING. `.admCovRowsFor`
  # is deterministic in `cov_dist` alone -- that is the CRN property above --
  # so every objective evaluation was re-running an identical Sobol draw.
  # Measured on 3 correlated covariates and 2 etas: 9.5 ms at n_row = 1000 and
  # 5000, 79 ms at 50000, PER STUDY PER EVALUATION. An rxSolve costs ~11 ms
  # before it integrates anything, so at cov_n_sim this was another whole solve
  # per study spent redrawing numbers that cannot change. Caching also cannot
  # weaken CRN: it hands back the identical rows by construction.
  #
  # Keyed on the covariate distribution, the row count and the eta count --
  # every input the draw has. An opaque user `joint` is NOT cached: digesting a
  # closure means digesting its environment, which can hold an arbitrarily
  # large object (a fitted vine, say), and paying that per call could cost more
  # than the draw it saves. That case behaves exactly as it did before.
  if (is.function(s$cov_dist[["joint"]])) {
    s$cov_rows <- .admCovRowsFor(s$cov_dist, n_row, pinfo$n_eta)
    return(s)
  }
  key <- digest::digest(list(s$cov_dist, n_row, pinfo$n_eta))
  hit <- .adm_covrows_env[[key]]
  if (is.null(hit)) {
    # a plain bound, so a long session fitting many models cannot grow it
    # without limit; the entries are per (study x n_row), so a handful per fit
    if (length(ls(.adm_covrows_env, all.names = TRUE)) > 64L)
      rm(list = ls(.adm_covrows_env, all.names = TRUE), envir = .adm_covrows_env)
    hit <- .admCovRowsFor(s$cov_dist, n_row, pinfo$n_eta)
    .adm_covrows_env[[key]] <- hit
  }
  s$cov_rows <- hit
  s
}

# Memo for the covariate row draws. An ENVIRONMENT, so it is caught by the
# generic is.environment() filter that keeps the daemon payload from shipping
# the parent's caches -- each worker keeps its own, which is what you want
# since the draw is cheap to rebuild and heavy to send.
.adm_covrows_env <- new.env(parent = emptyenv())

# Refuse `cov_dist` for an estimator that has no covariate path.
#
# Silence here is the dangerous outcome, not an error: every study carries a
# covariate VALUE as well (derived from the distribution when not given), so an
# unwired estimator does not fail -- it solves at the covariate mean and reports
# a fit whose omega has quietly absorbed the between-subject covariate spread.
# Measured on the general path's own test model, that is omega 0.30 -> 0.44.
.admRefuseCovariates <- function(studies, est) {
  # A fully stratified study still carries a cov_dist, but only of point specs
  # -- there is nothing to marginalise, so nothing for adfo or adirmc to refuse.
  has <- vapply(studies, function(s)
    !is.null(s[["cov_dist"]]) && !.admCovDistDegenerate(s[["cov_dist"]]),
    logical(1))
  if (!any(has)) return(invisible(NULL))
  stop("admixr2: `", est, "` does not support covariate marginalisation. ",
       "Stud", if (sum(has) > 1L) "ies " else "y ",
       paste(sQuote(names(studies)[has]), collapse = ", "),
       " declare(s) `cov_dist`. Use `admc` or `adgh`, both of which integrate ",
       "the covariate distribution -- dependent ones included; running here ",
       "would silently solve at the covariate mean and inflate omega instead.",
       call. = FALSE)
}

# A MARGINALISED DISCRETE COVARIATE WITH NO CONTRAST IS NOT IDENTIFIED, and it
# fails silently -- an ordinary-looking coefficient, finite SE and all.
#
# Marginalising a discrete covariate leaves its effect visible only through the
# MIXTURE it induces: the level probabilities shift E, and the spread between
# levels adds to V. Both are exactly what a random effect on the same parameter
# does, so where the level distribution is the SAME in every study there is no
# between-study contrast to break the confounding either, and the two are
# separated only by the shape difference between a two-point mixture and a
# lognormal -- fourth order, and worth almost nothing.
#
# Measured on a 400-subject source, one weight-banded study, sex declared
# 0.45/0.55 in every band and its coefficient estimated: profiling the
# objective in that coefficient moves it 0.019 units across its whole plausible
# range, against 18.3 units of curvature for the same coefficient once sex is
# STRATIFIED. The optimizer settled at -0.059 against a truth of +0.150, at
# every resolution, because a deterministic optimizer on a flat ridge stops in
# the same place every time -- which reads as a stable, converged answer.
#
# The remedy is cheap and is what a source that FITTED a sex effect actually
# supports: stratify on it. A discrete covariate needs no quadrature nodes, so
# stratifying multiplies the study count by its number of levels and nothing
# else, and the coefficient then comes back at 0.1500 at every resolution.
# `stratify = TRUE` already does this, since it bands every covariate whose
# coefficient the source's own model estimated.
# Does a random effect reach the same parameter this covariate modulates?
#
# Followed TRANSITIVELY, because a model routinely splits the two apart:
#   cl0 <- exp(tcl + eta.cl)
#   cl  <- cl0 * exp(bsex * SEX)
# is the same confounding as writing them on one line. Anything that cannot be
# parsed answers TRUE -- the caller only warns, and the failure this guards is
# silent.
.admCovMeetsEta <- function(ui, cov) {
  lst <- tryCatch(ui$lstExpr, error = function(e) NULL)
  ini <- tryCatch(ui$iniDf,   error = function(e) NULL)
  if (is.null(lst) || is.null(ini)) return(TRUE)
  etas <- unique(stats::na.omit(ini$name[!is.na(ini$neta1)]))
  if (!length(etas)) return(FALSE)
  asg <- Filter(function(e) is.call(e) && length(e) == 3L &&
                  as.character(e[[1L]]) %in% c("<-", "=") && is.name(e[[2L]]),
                lst)
  carry <- etas
  repeat {
    add <- character(0)
    for (e in asg) {
      v <- all.vars(e[[3L]])
      if (any(v %in% carry) && !as.character(e[[2L]]) %in% carry)
        add <- c(add, as.character(e[[2L]]))
    }
    if (!length(add)) break
    carry <- c(carry, add)
  }
  for (e in asg) {
    v <- all.vars(e[[3L]])
    if (cov %in% v && (any(v %in% carry) ||
                       as.character(e[[2L]]) %in% carry)) return(TRUE)
  }
  FALSE
}

.admCovDiscContrast <- function(.ui, studies, nms_has) {
  lvl <- function(sp) {
    p <- sp[["probs"]]
    if (is.null(p)) p <- rep(1, length(sp[["values"]]))
    paste(sprintf("%.10g", c(as.numeric(sp[["values"]]),
                             as.numeric(p) / sum(p))), collapse = "|")
  }
  seen <- list()
  for (nm in nms_has) {
    cd <- studies[[nm]]$cov_dist
    for (cv in .admCovSpecNames(cd)) {
      if (is.null(cd[[cv]][["values"]])) next
      seen[[cv]] <- c(seen[[cv]], lvl(cd[[cv]]))
    }
  }
  if (!length(seen)) return(invisible(NULL))
  # A study that PINS the covariate at a value carries a contrast in it, and so
  # does one that never declares it -- both are a different design point.
  for (cv in names(seen)) {
    pinned <- unique(unlist(lapply(studies, function(s) {
      v <- s[["cov"]][[cv]]
      if (is.null(v) || is.null(s$cov_dist) ||
          !is.null(s$cov_dist[[cv]][["values"]])) NULL else
        sprintf("%.10g", as.numeric(v)) })))
    if (length(unique(seen[[cv]])) > 1L || length(pinned) ||
        length(seen[[cv]]) < length(studies)) next
    est <- tryCatch(length(.admCovCoefThetas(.ui, cv, NULL)) > 0L,
                    error = function(e) TRUE)
    if (!isTRUE(est)) next
    # The confounding needs a PARTNER. Where no random effect reaches the
    # parameter the covariate modulates, the mixture is the only thing putting
    # spread on it and the coefficient is identified from V after all -- so
    # warning there would be a false positive, and a warning users learn to
    # ignore is worse than none.
    if (!isTRUE(.admCovMeetsEta(.ui, cv))) next
    warning("admixr2: '", cv, "' is a DISCRETE covariate whose coefficient ",
            "this model estimates, but every study marginalises over it with ",
            "the same level distribution, so nothing in the data separates ",
            "its effect from a random effect on the same parameter. The ",
            "coefficient will still be reported, and it will look converged: ",
            "measured on one source, the objective moves 0.019 units across ",
            "its whole range and the optimizer settled at -0.059 against a ",
            "truth of +0.150.\n",
            "  STRATIFY on it instead -- `stratify = ",
            if (length(names(seen)) > 1L)
              paste0("c(", paste(sQuote(names(seen)), collapse = ", "), ")") else
              sQuote(cv),
            "` in the study, or `stratify = TRUE` to derive it from the ",
            "source's own model. A discrete covariate needs no quadrature ",
            "nodes, so that costs one study per level and nothing else.\n",
            "  Or fix() the coefficient, if the source asserted it rather than ",
            "estimating it.", call. = FALSE)
  }
  invisible(NULL)
}

# Is there anything left to MARGINALISE, or is every margin a point?
#
# .admExpandStrata() gives each stratum a cov_dist carrying the stratified
# covariates as degenerate point specs (deliberately -- see pt_spec). Asking
# `!is.null(cov_dist)` therefore reports "this study marginalises" for a study
# that marginalises nothing, which made adfo and adirmc refuse a fully
# stratified source outright and made adgh and admc drop the gradient-based
# Hessian for one.
.admCovDistDegenerate <- function(cd) {
  nm <- .admCovSpecNames(cd)
  length(nm) > 0L &&
    all(vapply(nm, function(k) isTRUE(cd[[k]][[".point"]]), logical(1)))
}

# `cov` as a LIST, on every study.
#
# `cov` is documented as name -> value, and a named NUMERIC VECTOR is the
# natural way to write that -- but for an atomic vector `x[["absent"]]` is an
# ERROR ("subscript out of bounds") where a list returns NULL. Every consumer
# that asks "does this study pin covariate cv?" does so with `[[ ]]`, so the
# coercion has to happen before ANY of them run, not part-way down one function.
# It used to sit ten lines below .admCovDiscContrast() and was absent from
# .admWarnCovIdentifiability() entirely -- which every driver runs -- so a study
# written the documented way (`cov = c(WT = 70)` beside a cov_dist for SEX)
# aborted the fit with a message naming nothing.
# Does this covariate reach the model anywhere the probes cannot see?
#
# .admCovCollapse, .admJointCollapse and .admShiftSpec all decide which lines
# read a covariate by scanning TOP-LEVEL ASSIGNMENTS. `quote(if (a) b else c)`
# has length 4 and `quote(if (a) b)` has `if` as its head, so an if() is not an
# assignment either way and its covariate is invisible to all three.
#
# The consequence is not a refusal, it is a wrong design: with
# `if (CRCL < 30) fr <- 0.5 else fr <- 1` beside an allometric CL, the collapse
# is ADMITTED at r = 1 and CRCL takes exactly one value across the whole
# design -- its median -- so P(CRCL < 30) is 0 under the design against 0.0228
# declared and the renal branch is never exercised. Neither ver() nor the const
# re-probe can see it: both score only the `hit` columns, which the design
# reproduces exactly. On the shift path it is worse, since .admShiftVerify
# probes a single covariate row and never crosses a threshold at all.
#
# A branch is not affine in the covariate, so there is nothing to rotate onto
# even in principle. Refuse, and the study takes the ordinary product grid,
# which solves the whole model per row and handles if() correctly.
.admCovInBranch <- function(lst, covs) {
  if (is.null(lst) || !length(covs)) return(FALSE)
  for (e in lst) {
    if (!length(e)) next
    .asgn <- is.call(e) && length(e) == 3L &&
      as.character(e[[1L]])[1L] %in% c("<-", "=")
    if (.asgn) next
    if (any(covs %in% all.vars(e))) return(TRUE)
  }
  FALSE
}

.admCovAsList <- function(studies) {
  lapply(studies, function(s) {
    if (is.list(s) && !is.null(s[["cov"]]) && !is.list(s[["cov"]]))
      s[["cov"]] <- as.list(s[["cov"]])
    s
  })
}

# `est` goes LAST, and defaults to NULL meaning "build everything" -- the
# historical behaviour, which every Tier-1 mock and direct call relies on.
.admCheckCovariates <- function(.ui, pinfo, studies, est = NULL) {
  has <- vapply(studies, function(s) !is.null(s$cov_dist), logical(1))
  if (!any(has)) return(studies)
  bad <- function(...) stop("admixr2: ", ..., call. = FALSE)

  # Expand the friendly grammar ONCE, before anything routes on it: `cor`
  # becomes a `joint` sampler, and every downstream consumer (the per-row
  # draws, the product grid, the Taylor design, the shift grid) reads `joint`
  # rather than the margins. Expanding later would let a path integrate margins
  # whose correlation it had silently dropped.
  for (nm in names(studies)[has])
    studies[[nm]]$cov_dist <- .admCovDistCanon(studies[[nm]]$cov_dist)

  covs <- tryCatch(.ui$allCovs, error = function(e) character(0))
  # A JOINT unit has no per-row covariate path -- the shared-eta solve produces
  # one prediction set for every output at once -- so marginalising over a
  # covariate would silently solve at its mean. Refused HERE, at admission,
  # rather than only inside .adghMomentsJoint: that refusal sits on the NLL
  # path, so a fit could converge and then fail at the covariance, or take a
  # different route (batched moments, the plot's aggregate data) that never
  # asks. One admission gate covers all of them.
  .jt <- vapply(studies[has], function(s) isTRUE(s$is_joint), logical(1))
  if (any(.jt))
    stop("admixr2: covariate marginalisation is not supported for a JOINT ",
         "(same-subject, multi-output) unit -- study ",
         paste(sQuote(names(studies)[has][.jt]), collapse = ", "),
         ". The shared-eta joint solve has no per-row covariate path, so this ",
         "would silently solve at the covariate mean.", call. = FALSE)
  # A DESIGN IS A PROPERTY OF QUADRATURE, and admc has none: it samples through
  # .admCovRowsFor, which is deterministic in `cov_dist` alone -- the property
  # common random numbers depend on -- and re-aiming a design between
  # evaluations would make the DRAWS move with the parameters. So admc reads
  # none of .adm_cov_sparse / .adm_cov_collapse / .adm_cov_joint. Building them
  # anyway cost a probe and an SVD per study and, worse, printed "using 196
  # design points rather than ..." at a user whose fit uses neither number.
  .no_design <- identical(est, "admc")
  studies <- .admCovAsList(studies)      # BEFORE anything reads cov[[...]]
  .admCovDiscContrast(.ui, studies, names(studies)[has])

  for (nm in names(studies)[has]) {
    cd <- studies[[nm]]$cov_dist
    # `rho` and `Sigma` are metadata siblings of the covariate specs, not
    # covariates. Looping over them made a correlated spec fail with
    # "declares cov_dist for 'rho', which the model never reads" -- a message
    # pointing at the wrong thing.
    for (cv in .admCovSpecNames(cd)) {
      if (!cv %in% covs)
        bad("study '", nm, "' declares `cov_dist` for '", cv,
            "', which the model never reads. Model covariates: ",
            if (length(covs)) paste(covs, collapse = ", ") else "(none)", ".")
      sp <- cd[[cv]]
      if (!is.list(sp)) bad("`cov_dist` for '", cv, "' must be a list.")
      normal <- !is.null(sp$mu) && !is.null(sp$sd) &&
                is.finite(sp$mu) && is.finite(sp$sd) && sp$sd > 0
      lnorm  <- !is.null(sp$meanlog) && !is.null(sp$sdlog) &&
                is.finite(sp$meanlog) && is.finite(sp$sdlog) && sp$sdlog > 0
      disc   <- !is.null(sp$values) && length(sp$values) > 0L &&
                all(is.finite(sp$values))
      # A DROPPED LEVEL IS A TYPO THAT TWO ESTIMATORS READ DIFFERENTLY.
      # values = c(0, 1, 2) with probs = c(0.5, 0.5): admc's .admCovQuantile
      # renormalises over 2 entries, so cumsum is (0.5, 1) and level 2 is
      # UNREACHABLE; adgh's .admCovNodesFor returns 3 nodes against 2 weights,
      # which .admCovGrid recycles to (0.5, 0.5, 0.5) and normalises to a
      # UNIFORM 3-level covariate. Both finite, both plausible, neither the
      # declared distribution -- and the raw `values`(+`probs`) form is what
      # the error message just below invites.
      if (disc && !is.null(sp$probs)) {
        if (length(sp$probs) != length(sp$values))
          bad("`cov_dist` for '", cv, "' in study '", nm, "' gives ",
              length(sp$values), " `values` but ", length(sp$probs),
              " `probs`. Give one probability per level.")
        if (!all(is.finite(sp$probs)) || any(sp$probs < 0) || !sum(sp$probs) > 0)
          bad("`cov_dist` for '", cv, "' in study '", nm, "' has `probs` that ",
              "are not non-negative and summing to something positive.")
      }
      userq  <- is.function(sp$quantile)
      if (!(normal || lnorm || disc || userq))
        bad("`cov_dist` for '", cv, "' in study '", nm, "' is not a supported ",
            "distribution. Give one of: `mu`+`sd` (normal), `meanlog`+`sdlog` ",
            "(lognormal), `values`(+`probs`) (discrete/categorical), or ",
            "`quantile` (a function of a uniform).")

      # Fill the solve value from the distribution when the study omits it.
      if (is.null(studies[[nm]][["cov"]][[cv]])) {
        m <- .admCovMeanOf(sp)
        if (!is.null(m) && is.finite(m)) studies[[nm]][["cov"]][[cv]] <- m
      }
    }

    # THE DEFAULT PATH IS "rows", and it assumes nothing at all: every simulated
    # subject (admc) or grid point (adgh) carries its own covariate value, so
    # rxode2 evaluates the whole model -- covariate on several parameters, on a
    # parameter with no random effect, or interacting with one. It also carries
    # a gradient with no new chain rule, because on this path the covariate is
    # DATA (a per-row params column) and the existing sensitivity directions
    # differentiate exactly the function the NLL scores.
    #
    # It is the fallback for every refusal below, so a refusal costs solve rows
    # and never accuracy: both paths are differentiable.
    studies[[nm]]$.adm_cov_path <- "rows"

    # The covariate SHIFT was removed here. It pinned the covariate at its
    # reference and folded its whole contribution into one eta column, so the
    # integral collapsed onto u = Delta(a) + eta. .admJointCollapse below
    # finds the same structure -- rank 1 on the certified single-eta case --
    # without needing a certificate, and measured against a high-resolution
    # reference it is both cheaper and more accurate than the shift wherever
    # both applied: 14 rows at 2.05e-08 against 10 rows at 1.03e-06 at one
    # eta, 121 at 2.35e-07 against 50 at 2.58e-06 at two. Where the joint
    # refuses, the full product grid runs, which is the correct integral.

    # cov_integration = "sparse" replaces the product grid on the "rows" path
    # (and only there -- a study that took the shift has no product grid left to
    # expand). Its refusals are properties of `cov_dist` alone, so build the
    # grid once HERE, where the message can name the study, rather than letting
    # the first objective evaluation error out of the middle of a fit.
    if (identical(studies[[nm]]$.adm_cov_path, "rows") &&
        identical(pinfo$cov_integration %||% "on", "sparse") &&
        !.no_design) {
      .sg <- tryCatch(.admCovSparseGrid(cd, pinfo$cov_sparse_level %||% 3L,
                                        pinfo$cov_nodes %||% 7L),
                      error = function(e) conditionMessage(e))
      if (is.character(.sg))
        bad("study '", nm, "': ", sub("^admixr2: ", "", .sg))
      # KEEP it. The grid is a pure function of `cov_dist` and the level, both
      # DATA, but .adghGrid runs inside the objective -- so rebuilding it there
      # would repeat the combination technique's point merge on every
      # evaluation for a design that cannot change. Numeric only (no closures),
      # so it serialises to a parallel-restart worker by value like the rest of
      # the study.
      studies[[nm]]$.adm_cov_sparse <- .sg
    }
  }

  # EVERY covariate the ANALYSIS model reads must be described by a study that
  # has opted into covariate handling -- either a distribution to integrate
  # over, or a `cov` value if it genuinely does not vary in that study.
  #
  # Without this a covariate the model reads but the study never mentions is
  # silently held at whatever rxSolve defaults it to, which is the ecological
  # plug-in wearing a fit's clothes: finite, plausible, and biased.
  #
  # Checked LAST, and only for studies declaring some `cov_dist`. Last because
  # a mistyped covariate name fails BOTH this and the "declares a covariate the
  # model never reads" check above, and that one names the typo directly. Only
  # for opted-in studies so a fit handling covariates entirely through fixed
  # `cov` values is untouched.
  for (nm in names(studies)[has]) {
    dcl <- c(.admCovSpecNames(studies[[nm]]$cov_dist),
             names(studies[[nm]][["cov"]] %||% list()))
    miss <- setdiff(covs, dcl)
    if (length(miss))
      bad("study '", nm, "' does not describe covariate(s) ",
          paste(sQuote(miss), collapse = ", "),
          ", which the model reads. Give each one a `cov_dist` entry (the ",
          "distribution its subjects span) or, if it does not vary in that ",
          "study, a fixed `cov` value. Leaving it out solves at whatever ",
          "value rxSolve defaults to, which is the ecological plug-in.")
  }

  studies
}

# =============================================================================
# General covariate marginalisation (ANY functional form)
# =============================================================================
#
# The DEFAULT path, and the fallback for every refusal the shift path makes.
# (WT/70)^theta, exp(theta*WT), allometric-in-log, Emax, if/else, categorical, a
# covariate on several parameters or on one with no random effect -- all of it
# goes through here, with no structural precondition of any kind.
#
# The trick is that we never have to KNOW the functional form. A covariate
# reaches rxode2 only as a column of the params frame, so if each simulated
# subject carries its OWN covariate value, rxode2 evaluates whatever the model
# contains and the aggregate moments come out of the pooled ensemble. That is
# also why it costs nothing extra for admc: n_sim subjects still means n_sim
# rows, each now carrying its own (covariate, eta) pair rather than a shared
# covariate.

# Quantile function for one covariate spec. Supported:
#   list(mu, sd)              normal
#   list(meanlog, sdlog)      lognormal
#   list(values, probs)       discrete / categorical (probs default to uniform)
#   list(quantile = f(u))     anything else, supplied by the caller
# Canonical form for a user-supplied `cov_dist`.
#
# What a paper reports is a MEAN and an SD on the natural scale, and a
# CORRELATION. What the internals consume is `meanlog`/`sdlog` (or `mu`/`sd`)
# and a `joint` sampler. Every caller used to bridge that gap by hand -- each
# one re-deriving the lognormal moment match and writing its own Gaussian
# copula closure -- which is both tedious and a place to get the algebra
# quietly wrong.
#
# This is the ONE place the friendly grammar is expanded, and it is idempotent,
# so it can be applied defensively wherever a spec enters. It runs BEFORE
# routing on the estimator path, which matters: `cor` becomes a `joint`
# sampler, and `joint` is what adgh refuses. Expanding later would let adgh
# integrate a product grid over margins whose correlation it had silently
# dropped.
#
#   cov_dist = list(WT   = list(mean = 72, sd = 16),
#                   CRCL = list(mean = 90, sd = 25),
#                   cor  = 0.6)
#
# `dist` selects the margin: "lnorm" (default -- covariates in this setting are
# positive, and a normal margin puts mass at and below zero, where a power
# covariate model returns NaN) or "normal".
#
# `cor` is the GAUSSIAN COPULA correlation -- the correlation of the latent
# normals, not the Pearson correlation of the covariates themselves. The two
# coincide for `dist = "normal"` margins and differ slightly otherwise: with
# the lognormal defaults above, cor = 0.6 realises a Pearson correlation of
# 0.592. That is a property of the copula construction, not an approximation,
# and the quadrature grid and the per-subject sampler agree on it to 4 decimal
# places, so the two estimator families see the SAME distribution. Supply a
# `joint` sampler directly if a specific Pearson correlation must be matched.
.admCovMomentMatch <- function(m, sd, nm, dist) {
  # sd == 0 IS refused, not just a negative one: a constant covariate has no
  # spread to integrate, every design point coincides, and the failure surfaces
  # much later and much less legibly. The guard was added to .admPopFromData
  # only, so both typed-out routes still reached it.
  if (is.null(sd) || !is.finite(sd) || sd <= 0)
    stop("admixr2: covariate ", sQuote(nm), " gives `mean` without a finite ",
         "POSITIVE `sd`.", call. = FALSE)
  if (identical(dist, "normal")) return(list(mu = m, sd = sd))
  if (!is.finite(m) || m <= 0)
    stop("admixr2: covariate ", sQuote(nm), " has mean ", m, ", which a ",
         "lognormal margin cannot represent. Give a positive mean, or ",
         'dist = "normal" if the covariate really is unbounded below.',
         call. = FALSE)
  # match the natural-scale mean and SD exactly
  list(meanlog = log(m^2 / sqrt(sd^2 + m^2)),
       sdlog   = sqrt(log(1 + sd^2 / m^2)))
}

.admCovDistCanon <- function(cov_dist) {
  if (!is.list(cov_dist) || !length(cov_dist)) return(cov_dist)
  nms <- .admCovSpecNames(cov_dist)
  if (!length(nms)) return(cov_dist)

  for (nm in nms) {
    sp <- cov_dist[[nm]]
    if (!is.list(sp) || is.null(sp[["mean"]])) next
    # an explicit canonical form always wins; `mean` is only a shorthand
    if (!is.null(sp$quantile) || !is.null(sp$values) ||
        !is.null(sp$meanlog)  || !is.null(sp$mu)) next
    d <- tolower(sp$dist %||% "normal")
    if (!d %in% c("lnorm", "lognormal", "normal", "norm"))
      stop("admixr2: covariate ", sQuote(nm), ' has dist = "', sp$dist,
           '"; only "lnorm" and "normal" are understood as shorthands. ',
           "Supply a `quantile` function for anything else.", call. = FALSE)
    mm <- .admCovMomentMatch(sp[["mean"]], sp[["sd"]], nm,
                             if (d %in% c("normal", "norm")) "normal" else "lnorm")
    sp[["mean"]] <- NULL; sp[["dist"]] <- NULL
    if (!is.null(mm$meanlog)) sp[["sd"]] <- NULL
    cov_dist[[nm]] <- utils::modifyList(sp, mm)
  }

  # THREE spellings of the same statement, and they must all reach every path.
  # `cor` is the user-facing one; `rho`/`Sigma` are the retired collapse's older
  # normal-scale form, still accepted because published specs are written that
  # way. Until this was unified, a spec written with `rho` built its Gaussian
  # copula for the collapse and a DIAGONAL grid for everything else -- the
  # correlation silently present in one path and absent in the other.
  d  <- length(nms)
  cr <- cov_dist[["cor"]]
  cov_dist[["cor"]] <- NULL
  R  <- NULL
  if (!is.null(cr)) {
    if (d < 2L)
      stop("admixr2: `cor` needs at least two covariates; ", d, " declared.",
           call. = FALSE)
    R <- if (length(cr) == 1L) {
      if (d != 2L)
        stop("admixr2: a scalar `cor` is only unambiguous for two covariates; ",
             d, " declared, so give a ", d, " x ", d, " correlation matrix.",
             call. = FALSE)
      matrix(c(1, cr, cr, 1), 2L, 2L)
    } else as.matrix(cr)
    if (!identical(dim(R), c(d, d)))
      stop("admixr2: `cor` is ", nrow(R), " x ", ncol(R), " but ", d,
           " covariates are declared.", call. = FALSE)
    # Order the matrix to the declared covariates when it says which is which;
    # a named matrix in a different order would otherwise correlate the wrong
    # pair. This has to happen HERE, in the branch that can receive a named
    # matrix from the user -- doing it after the branches ran once against a
    # matrix whose dimnames had already been stripped, i.e. never at all.
    if (!is.null(rownames(R)) && all(nms %in% rownames(R)))
      R <- R[nms, nms, drop = FALSE]
  } else if (!is.null(cov_dist[["Sigma"]])) {
    S <- as.matrix(cov_dist[["Sigma"]])
    if (!identical(dim(S), c(d, d)))
      stop("admixr2: `Sigma` is ", nrow(S), " x ", ncol(S), " but ", d,
           " covariates are declared.", call. = FALSE)
    if (!is.null(rownames(S)) && all(nms %in% rownames(S)))
      S <- S[nms, nms, drop = FALSE]
    if (any(diag(S) <= 0))
      stop("admixr2: `Sigma` has a non-positive diagonal entry, so it ",
           "describes no distribution.", call. = FALSE)
    R <- stats::cov2cor(S)
  } else if (!is.null(cov_dist[["rho"]]) &&
             any(abs(as.numeric(cov_dist[["rho"]])) > 1e-12)) {
    rr <- as.numeric(cov_dist[["rho"]])[1L]
    if (d < 2L)
      stop("admixr2: `rho` needs at least two covariates; ", d, " declared.",
           call. = FALSE)
    # scalar rho applies to EVERY pair
    R <- matrix(rr, d, d); diag(R) <- 1
  }
  if (is.null(R)) return(cov_dist)
  # An explicit sampler is the more specific statement; do not override it.
  if (is.function(cov_dist[["joint"]])) return(cov_dist)
  dimnames(R) <- NULL                    # ordered above, in the branch that can
  Lc <- tryCatch(chol(R), error = function(e) NULL)   # receive a named matrix
  if (is.null(Lc))
    stop("admixr2: `cor` is not positive definite, so it describes no ",
         "distribution.", call. = FALSE)
  # admixr2 built this copula, so it also knows its DENSITY -- and supplying it
  # lets covStrata() condition on an exact covariate value by
  # sampling-importance-resampling rather than by binning, which removes the
  # attenuation a band carries. Gaussian copula density on the copula scale:
  #   c(u) = |R|^-1/2 exp(-0.5 z' (R^-1 - I) z),   z = qnorm(u)
  # NOT named `cor_matrix`: `$` PARTIAL-MATCHES, so `cov_dist$cor` would have
  # silently returned this matrix everywhere `cor` was read.
  cov_dist[["latentR"]] <- R
  cov_dist[["discExact"]] <- .admCovDiscExact(cov_dist, nms, R)
  margins <- lapply(nms, function(nm) cov_dist[[nm]])
  # Marked as OURS. covStrata() truncates the margins and re-canonicalises, and
  # the early return at "an explicit sampler is the more specific statement"
  # would otherwise hand back a closure still holding the UNTRUNCATED margins --
  # so bands were cut over the full declared support with no warning (the
  # cov_range gate had been satisfied), which is the var(declared)/var(enrolled)
  # overstatement that gate exists to prevent. A user's own `joint` is still
  # never rebuilt: only a closure carrying this flag is discarded.
  cov_dist[["jointOwn"]] <- TRUE
  cov_dist[["joint"]] <- local({
    nms <- nms; margins <- margins; Lc <- Lc; d <- d
    function(u) {
      # ONE tolerance, shared with every other site -- .admCovU(). This closure
      # used to clamp at 1e-12, and .admCovGrid hands it uniforms it has ALREADY
      # clamped at .Machine$double.eps, so the looser bound silently won on the
      # dependent path: the same design point sat at qnorm -7.03 through a
      # correlated spec and -8.13 through an independent one.
      cl  <- function(x) pmin(pmax(x, .Machine$double.eps),
                              1 - .Machine$double.eps)
      z   <- stats::qnorm(cl(u)) %*% Lc
      out <- vapply(seq_len(d), function(j)
        .admCovQuantile(margins[[j]], .admCovU(z[, j])),
        numeric(nrow(u)))
      out <- matrix(out, nrow = nrow(u), ncol = d)
      colnames(out) <- nms
      out
    }
  })
  cov_dist
}

.admCovQuantile <- function(spec, u) {
  # `values` FIRST, matching .admCovNodesFor and .admCovMeanOf,
  # .admCovTruncSpec and .admCovDiscExact. Testing `quantile` first made a spec
  # carrying both integrate as a CONTINUOUS margin under admc while adgh
  # enumerated its discrete levels -- two estimators, two distributions, both
  # finite and plausible. (A degenerate point spec carries only `quantile`, so
  # it is unaffected.)
  if (!is.null(spec$values)) {
    pr <- spec$probs %||% rep(1 / length(spec$values), length(spec$values))
    pr <- pr / sum(pr)
    return(as.numeric(spec$values)[findInterval(u, cumsum(pr), rightmost.closed = TRUE) + 1L])
  }
  if (is.function(spec$quantile)) return(as.numeric(spec$quantile(u)))
  if (!is.null(spec$meanlog)) return(stats::qlnorm(u, spec$meanlog, spec$sdlog))
  stats::qnorm(u, spec$mu, spec$sd)
}

# Deterministic per-row covariate values for `n` simulated subjects.
#
# Deterministic on purpose: the covariate distribution is DATA, not a parameter,
# so the same rows must come back on every objective evaluation or the optimizer
# sees noise. Common random numbers therefore hold with no seed plumbing.
#
# The uniforms are taken from Sobol dimensions AFTER the eta dimensions.
# sobol(n, dim = k)[, 1:j] is exactly sobol(n, dim = j) (verified), so this
# yields dimensions genuinely different from the ones .admMakeZ used for eta.
# Drawing a separate halton/sobol sequence instead would NOT: every low-
# discrepancy family starts from the same base-2 van der Corput sequence, so
# covariate column 1 would have been a copy of eta column 1.
# Per-subject covariate draws for the general path.
#
# DEPENDENT COVARIATES. `cov_dist$joint` is a function taking the n x d matrix of
# uniforms admixr2 draws and returning the n x d matrix of covariate values. That
# is exactly the shape a copula produces -- an R-vine included: sample the vine
# on the uniform scale, then push each column through its own marginal quantile
# function. Aggregate patient statistics are dependent (weight with age, weight
# with height, creatinine with age), and the per-covariate branch below draws
# each margin independently, so it cannot represent that.
#
# The uniforms come from ADMIXR2's Sobol stream, deliberately, and a user sampler
# must consume them rather than draw its own. The stream is what makes the
# objective a deterministic function of the parameters: common random numbers
# across optimizer iterations is what lets a finite difference of it mean
# anything, and a sampler calling RVineSim() internally would reseed every
# evaluation and turn the objective into noise. The dimensions sit AFTER the
# random-effect dimensions so the eta draws are unchanged by adding a covariate.
#
#   cov_dist = list(
#     WT  = list(quantile = function(u) qlnorm(u, log(70), 0.25)),
#     AGE = list(quantile = function(u) qgamma(u, 9, 0.3)),
#     joint = function(u) {                       # u is n x 2, columns WT, AGE
#       v <- VineCopula::RVineSim(nrow(u), RVM, U = u)   # dependence
#       cbind(WT  = qlnorm(v[, 1], log(70), 0.25),
#             AGE = qgamma(v[, 2], 9, 0.3))
#     })
#
# Passing `U` keeps the vine driven by admixr2's stream. A sampler that ignores
# its argument still runs, and still gives the right MARGINAL answer in
# expectation, but the objective becomes stochastic and gradients degrade -- so
# the contract is stated here rather than enforced, since only the caller knows
# whether its sampler is deterministic in `u`.
.admCovRowsFor <- function(cov_dist, n, n_eta) {
  # `rho`/`Sigma` are dependence metadata; `joint` is the sampler itself.
  cov_dist <- .admCovDistCanon(cov_dist)
  nms <- .admCovSpecNames(cov_dist)
  d   <- length(nms)
  u   <- randtoolbox::sobol(n, dim = n_eta + d)
  if (!is.matrix(u)) u <- matrix(u, nrow = n)
  u   <- u[, n_eta + seq_len(d), drop = FALSE]
  # sobol emits an exact 0 in its first row; qnorm(0) is -Inf.
  u   <- pmin(pmax(u, .Machine$double.eps), 1 - .Machine$double.eps)
  colnames(u) <- nms

  jf <- cov_dist[["joint"]]
  if (is.function(jf)) {
    out <- tryCatch(as.matrix(jf(u)), error = function(e)
      stop("cov_dist$joint failed on the ", n, " x ", d, " uniform matrix: ",
           conditionMessage(e), call. = FALSE))
    if (!is.matrix(out) || nrow(out) != n)
      stop("cov_dist$joint must return a matrix with one ROW per subject (",
           n, "); got ", nrow(out), ".", call. = FALSE)
    if (is.null(colnames(out)) || !setequal(colnames(out), nms))
      stop("cov_dist$joint must return columns named ",
           paste(sQuote(nms), collapse = ", "), "; got ",
           if (is.null(colnames(out))) "none" else
             paste(sQuote(colnames(out)), collapse = ", "), ".", call. = FALSE)
    out <- out[, nms, drop = FALSE]
    if (!all(is.finite(out)))
      stop("cov_dist$joint returned non-finite covariate values.", call. = FALSE)
    return(out)
  }

  .admCovXFromU(cov_dist, nms, u)
}

# =============================================================================
# Shared covariate-design primitives
# =============================================================================
#
# Five designs -- the product grid, the sparse grid, the strata, the covariate
# collapse and the joint collapse -- do the same four things to a `cov_dist`:
# enumerate its discrete margins, refuse a discrete margin that is latently
# correlated with a continuous one, factorise the continuous correlation, and
# push latent normal nodes through each margin's quantile function. Each used to
# spell all four out.
#
# That is the shape CLAUDE.md names as the source of the covariate bugs: "the
# residual row arrays are rebuilt by SIX consumers, and that is where the bugs
# are". It has already cost one here -- .admCovSparseGrid was written with the
# enumeration and the factorisation copied in but the REFUSAL left out, so a
# discrete covariate correlated with a continuous one was integrated as if it
# were independent until a converted test caught it.

# The uniform a latent normal node maps to.
#
# ONE tolerance. pnorm() saturates to exactly 0 or 1 in the tails -- the grid
# reaches |z| ~ 8 at 21 nodes and a copula's mixing step pushes that further --
# after which a margin's quantile function returns +/-Inf. Nine sites clamped at
# .Machine$double.eps and one at 1e-12, which put the extreme node at qnorm
# -8.13 against -7.03, i.e. the same design point in a different place
# depending on which function built it.
.admCovU <- function(z)
  pmin(pmax(stats::pnorm(z), .Machine$double.eps), 1 - .Machine$double.eps)

# Latent normal nodes -> covariate values, one column per margin.
#
# Z is n x length(cn) on the LATENT scale (already rotated and scaled by the
# caller); the result is n x length(cn) on each covariate's own scale, named.
# The matrix() is not decoration: vapply drops to a vector at n == 1 or at one
# covariate, and four of the six callers carried their own `if (!is.matrix(x))`
# repair for exactly that.
.admCovXFromZ <- function(cd, cn, Z)
  .admCovXFromU(cd, cn, .admCovU(as.matrix(Z)))

# The same push, from UNIFORMS rather than latent normals -- what a copula's
# output needs, since inverse_rosenblatt() and admixr2's own sampler both hand
# back the copula scale and neither knows the margins.
.admCovXFromU <- function(cd, cn, U) {
  U <- as.matrix(U)
  matrix(vapply(seq_along(cn), function(k)
    .admCovQuantile(cd[[cn[k]]], U[, k]), numeric(nrow(U))),
    nrow(U), length(cn), dimnames = list(NULL, cn))
}

# Exact enumeration of the discrete margins: their levels and probabilities ARE
# the integration rule, so they are crossed whole rather than put on any rule.
#
# Returns the pieces every caller wanted between them -- the per-margin levels
# and normalised probabilities, the crossed cells with their probabilities, and
# the cells as named lists for the probe. Reads each margin through
# .admCovNodesFor, which is where "probs, defaulting to uniform, normalised"
# already lived; five sites re-implemented that line in three spellings.
.admCovDiscCells <- function(cd, dn) {
  if (!length(dn))
    return(list(lv = list(), pr = list(), cells = matrix(numeric(0), 1L, 0L),
                pcell = 1, cell_list = list(list())))
  nodes <- lapply(dn, function(n) .admCovNodesFor(cd[[n]], 1L))
  lv <- lapply(nodes, `[[`, "x")
  pr <- lapply(nodes, `[[`, "w")
  cells <- as.matrix(expand.grid(lv, KEEP.OUT.ATTRS = FALSE))
  colnames(cells) <- dn
  pcell <- apply(as.matrix(expand.grid(pr, KEEP.OUT.ATTRS = FALSE)), 1L, prod)
  list(lv = lv, pr = pr, cells = cells, pcell = as.numeric(pcell),
       cell_list = lapply(seq_len(nrow(cells)), function(i)
         as.list(stats::setNames(cells[i, ], dn))))
}

# The latent structure a design may use, or NULL if it may not.
#
# Refuses a DISCRETE margin latently correlated with a CONTINUOUS one: a level
# is then a truncation of the latent normal rather than a point, so the
# continuous conditional differs cell to cell and one shared design is the wrong
# design in every cell. Then subsets the correlation to the continuous block and
# factorises it. `R` is indexed POSITIONALLY -- latentR carries no dimnames --
# so `nms` must be the full declared order.
.admCovLatentBlock <- function(cd, nms, cn, dn, R) {
  pc <- length(cn)
  ic <- match(cn, nms); id <- match(dn, nms)
  if (length(dn) && !is.null(R) && any(abs(R[id, ic, drop = FALSE]) > 0))
    return(NULL)
  Rc <- if (is.null(R)) diag(1, pc) else R[ic, ic, drop = FALSE]
  Lc <- tryCatch(chol(Rc), error = function(e) NULL)
  if (is.null(Lc)) return(NULL)
  c(list(ic = ic, id = id, Rc = Rc, Lc = Lc), .admCovDiscCells(cd, dn))
}

# Deterministic quadrature nodes + weights for one covariate (adgh path).
# Discrete specs enumerate exactly; normal/lognormal use Gauss-Hermite.
.admCovNodesFor <- function(spec, n_nodes) {
  if (!is.null(spec$values)) {
    pr <- spec$probs %||% rep(1 / length(spec$values), length(spec$values))
    return(list(x = as.numeric(spec$values), w = pr / sum(pr)))
  }
  # A DEGENERATE POINT NEEDS ONE NODE, NOT n_nodes COPIES OF ITSELF. A
  # stratified covariate is held at a point inside its stratum (pt_spec), and a
  # point spec carries only `quantile`, so it fell through to the branch below
  # and was handed the whole standard-normal axis -- n_nodes rows with
  # IDENTICAL covariate values and Gauss-Hermite weights summing to one. The
  # answer was right and the grid was n_nodes times too big PER STRATIFIED
  # COVARIATE: 49x on a band cut in two of them, every extra row an exact
  # duplicate solve.
  if (isTRUE(spec[[".point"]]))
    return(list(x = as.numeric(spec$quantile(0.5)), w = 1))
  g <- .adghNodes1(n_nodes)                       # standard-normal nodes/weights
  # A user-supplied quantile function. E_a[h(a)] = E_z[h(F^-1(Phi(z)))] for
  # z ~ N(0,1), so pushing the standard-normal nodes through Phi and then F^-1
  # is an exact quadrature for ANY margin -- which is also what makes this the
  # hook a copula-based joint sampler plugs into.
  #
  # Without this branch the function fell through to `spec$mu + spec$sd * g$x`,
  # and a quantile spec has neither: x came back numeric(0) against 7 weights,
  # .admCovGrid built a zero-row grid, and adgh -- which VALIDATES the spec as
  # supported in .admCheckCovariates -- died on the first objective evaluation.
  # admc was unaffected, because it goes through .admCovQuantile instead.
  # CLAMP before a user quantile function sees it. pnorm() returns exactly 1
  # from |z| >= 8.30, which this grid reaches at 31 nodes, and an unbounded
  # quantile is then infinite -- measured: qweibull(1, 1.5, 55) = Inf at 41
  # nodes, which propagates silently into the moments. The clamp costs nothing
  # below 31 nodes, where nothing is clamped at all. (The joint branch of
  # .admCovGrid has the same guard for the same reason.)
  if (is.function(spec$quantile))
    return(list(x = as.numeric(spec$quantile(
                     .admCovU(g$x))), w = g$w))
  # the closed forms are safe: they are built from g$x directly, never from a
  # probability, so no saturation can occur
  if (!is.null(spec$meanlog)) list(x = exp(spec$meanlog + spec$sdlog * g$x), w = g$w)
  else                        list(x = spec$mu + spec$sd * g$x,              w = g$w)
}


# Metadata siblings of the per-covariate specs. Named once, because every
# consumer that enumerated `names(cov_dist)` and forgot them produced a
# different wrong answer: .admCovGrid tried to read `$values` off a numeric
# `rho` ("$ operator is invalid for atomic vectors") and off the `joint`
# CLOSURE ("object of type 'closure' is not subsettable"), and
# .admCheckCovariates reported "declares cov_dist for 'rho', which the model
# never reads".
.ADM_COV_META <- c("rho", "Sigma", "cor", "joint", "jointOwn", "latentR",
                   "discExact")


# Discrete margins a `joint` sampler maps STRAIGHT FROM THEIR OWN UNIFORM.
#
# A discrete margin normally cannot ride a quadrature grid once the covariates
# are dependent: the sampler mixes the uniforms before mapping them to levels,
# so fixing the INPUT uniform does not fix the OUTPUT level and the cell
# weights are not the level probabilities (measured: a covariate declared
# 0.55/0.45 came off the grid at 0.477, and the error does not shrink with
# `cov_nodes`, because it is a property of the mixing rather than of the
# resolution).
#
# That is a property of THE SAMPLER, not of discreteness. Whenever column j is
# latently independent of the rest, chol(R)[, j] is e_j, so the copula's
# `z <- qnorm(u) %*% Lc` leaves `z[, j] = qnorm(u[, j])` and the level is a
# monotone function of that margin's own uniform after all. Such a margin
# enumerates EXACTLY at its levels, and it is common: one declared sex,
# genotype or formulation alongside a correlated (WT, CRCL) pair.
#
# Recorded by whoever BUILDS the sampler, because only they can know it. A user
# `joint` is opaque and never gets the flag.
.admCovDiscExact <- function(cov_dist, nms, R, tol = 1e-12) {
  keep <- vapply(seq_along(nms), function(j) {
    if (is.null(cov_dist[[nms[j]]][["values"]])) return(FALSE)
    all(abs(R[j, -j, drop = TRUE]) < tol)
  }, logical(1))
  nms[keep]
}

.admCovSpecNames <- function(cov_dist) setdiff(names(cov_dist), .ADM_COV_META)

# =============================================================================
# Covariate STRATA -- the per-covariate stratify/marginalise split
# =============================================================================
#
# A published source conditions its model on SOME covariates and not others,
# and the matching rule is per covariate, not per study:
#
#   STRATIFY   on the covariates the source's own model fitted
#   MARGINALISE over the ones it did not, using the source's covariate
#              distribution -- on BOTH the observed and the predicted side
#
# A full product grid over ALL covariates for EVERY source is the failure this
# exists to prevent: a model with no term in x is evaluated at nodes that vary
# x, answers "no change" at every one, and that FABRICATED null contrast is
# scored as evidence against sources that have a real one. Measured with no
# sampling noise, against a true 0.450: 0.2107 (-53.2%) fabricated versus
# 0.4500 matched.
#
# GENERATION-SIDE, NOT AN ESTIMATOR SETTING, because stratifying needs
# per-stratum OBSERVATIONS. Scoring one pooled (E, V) against K conditional
# predictions is -2 log of an unnormalised geometric mean, which is not the
# likelihood of anything. So strata arrive as ordinary studies, each with its
# own n and data, and the estimator needs no knowledge of where they came from
# (verified bit-identical to fitting them as unrelated studies).
#
# CONDITIONING. Within stratum k the covariates the source did NOT fit are
# marginalised over their distribution CONDITIONAL on that stratum. The
# unconditional shortcut predicts the average-x2 response at every x1 node while
# the source's high-x1 subjects had high x2 -- a mean error that VARIES ACROSS
# NODES, the shape that biases the stratified covariate's own coefficient.
# Relative error of the block mean: 1e-3 conditional against 0.20 at rho = 0.3
# and 0.96 at rho = 0.85. They coincide exactly at rho = 0, which is why an
# independent-covariate test never showed it.
#
# Done in U-SPACE, which is what makes it work for an arbitrary sampler: a
# copula maps INDEPENDENT uniforms to dependent values, so holding the leading
# uniforms fixed and varying the rest samples the conditional distribution --
# for a vine, its Rosenblatt coordinates. Nothing inverts the sampler and no
# forward Rosenblatt is needed, because a stratum is defined BY ITS NODE in
# u-space and the covariate value there is read off the sampler.
#
# THAT ROUTE IS GONE FOR AN OPAQUE SAMPLER, and this paragraph describes what
# replaced it only for admixr2's own `cor` copula, which is unwrapped to
# margins plus `latentR` before it gets here. Fixing the stratified
# covariates' INPUT uniforms requires them to lead the sampler's own
# conditioning cascade, and fails silently when they do not -- so a user
# `joint` bins the sampler's OUTPUT instead, which has no such requirement.
# See the binning branch below for the measurement that forced it.
#
# Do not reintroduce a cascade contract on the strength of the paragraph
# above: it was read that way once and cost a column reversal that bound
# covariate k to the sampler's variable d-k+1.
# Strata per stratified covariate.
#
# A CONVERGENCE PARAMETER, NOT A MODELLING CHOICE -- the old default of 5 read
# like a preference. The well-defined object is the J -> infinity limit
# (N * E_a[l]); any finite J approximates it, and under misspecification the
# answer can jump between basins rather than drift: a verified counterexample
# flips from beta = 1.57 at J = 4 to 1.6e-05 at J = 5. Drive it up until the
# answer stops moving. 9 is a starting point; cost is J^p studies for p
# stratified covariates, so affordable in one and expensive in three.
#
# How J-dependent the objective is depends on the rule, and on the default it
# is barely at all -- measured across J = 5 to 100 on one banded source:
# Gauss-Hermite 0.03 units, pooled bins (an opaque `joint` only) 51 and still
# moving. So OFV, AIC, BIC and a likelihood ratio ARE comparable across
# resolutions on the default path and are not on the pooled one, where the
# estimate also wanders (0.6810 at J = 5, 9, 15 and 100 but 0.6968 at J = 50).
# The value is stamped onto every generated study and carried onto the fit so
# anova() can refuse the comparison.
.ADM_STRATA_NODES <- 9L

# Truncate one covariate's margin to the range a source actually enrolled.
#
# Strata are cut from the analyst's `cov_dist` over its FULL support, so a
# published model gets evaluated -- and credited as evidence -- in covariate
# bands where that study enrolled nobody. The inflation is exact:
#
#     information_claimed / information_earned  =  var_assumed / var_enrolled
#
# measured at 3.43x for a source that enrolled +/- 1 SD, and 12.38x at +/- 0.5
# SD. It is not bias, estimates stay correct; it is FALSE CONFIDENCE, and it
# only bites once a second source disagrees -- which is exactly when it matters.
#
# Truncating the SPEC rather than either branch of .admCovStrata means both
# inherit it: the exact-conditioning route and the pooled route both reach the
# margin through .admCovQuantile.
.admCovTruncSpec <- function(spec, rng, nm) {
  bad <- function(...) stop("admixr2: ", ..., call. = FALSE)
  if (is.null(rng)) return(spec)
  rng <- sort(as.numeric(rng))
  if (length(rng) != 2L || !all(is.finite(rng)) || rng[1L] >= rng[2L])
    bad("`cov_range` for ", sQuote(nm), " must be two finite increasing ",
        "values, e.g. c(52, 118).")
  # DISCRETE: a range keeps the levels inside it and renormalises
  if (!is.null(spec[["values"]])) {
    v  <- as.numeric(spec[["values"]])
    pr <- spec[["probs"]] %||% rep(1 / length(v), length(v))
    k  <- v >= rng[1L] & v <= rng[2L]
    if (!any(k)) bad("`cov_range` for ", sQuote(nm), " excludes every level.")
    return(list(values = v[k], probs = pr[k] / sum(pr[k])))
  }
  # CONTINUOUS: re-map u onto [F(a), F(b)], which is the truncated quantile
  # function. Weights renormalise automatically over the truncated support, so
  # sum(n_k) = n still holds.
  cdf <- if (!is.null(spec[["meanlog"]]))
    function(x) stats::plnorm(x, spec[["meanlog"]], spec[["sdlog"]])
  else if (!is.null(spec[["mu"]]))
    function(x) stats::pnorm(x, spec[["mu"]], spec[["sd"]])
  else NULL
  qf <- if (!is.null(spec[["meanlog"]]))
    function(u) stats::qlnorm(u, spec[["meanlog"]], spec[["sdlog"]])
  else if (!is.null(spec[["mu"]]))
    function(u) stats::qnorm(u, spec[["mu"]], spec[["sd"]])
  else if (is.function(spec[["quantile"]])) spec[["quantile"]]
  else bad("`cov_range` was given for ", sQuote(nm), " but its margin is not ",
           "one admixr2 can invert (declare it as normal, lognormal, or with ",
           "a `quantile` function).")
  if (is.null(cdf)) {
    # a user quantile function: invert on a fine grid, which is enough since
    # this only sets the two endpoints of the truncation
    gu <- seq(1e-6, 1 - 1e-6, length.out = 20001L)
    gx <- as.numeric(qf(gu))
    cdf <- function(x) stats::approx(gx, gu, xout = x, rule = 2L)$y
  }
  pa <- as.numeric(cdf(rng[1L])); pb <- as.numeric(cdf(rng[2L]))
  if (!is.finite(pa) || !is.finite(pb) || pb - pa < 1e-8)
    bad("`cov_range` for ", sQuote(nm), " covers essentially none of its ",
        "declared distribution -- check the range and the margin agree on ",
        "units.")
  list(quantile = function(u) qf(pa + u * (pb - pa)))
}

.admCovStrata <- function(cov_dist, stratify, n_nodes = 5L,
                          n_pool = 32768L, cov_range = NULL) {
  cov_dist <- .admCovDistCanon(cov_dist)
  nms <- .admCovSpecNames(cov_dist)
  bad <- function(...) stop("admixr2: ", ..., call. = FALSE)
  if (!length(stratify)) return(NULL)
  if (!is.character(stratify))
    bad("`stratify` must be a character vector of covariate names.")
  if (!all(stratify %in% nms))
    bad("`stratify` names covariate(s) ",
        paste(sQuote(setdiff(stratify, nms)), collapse = ", "),
        " that this study's `cov_dist` does not declare. Declared: ",
        paste(sQuote(nms), collapse = ", "), ".")
  checkmate::assertCount(n_nodes, positive = TRUE)
  # BAND ONLY WHERE THE SOURCE ENROLLED. Without a reported range the strata are
  # cut over the whole declared distribution, and the source is credited with
  # evidence from bands it never sampled -- see .admCovTruncSpec(). Warn rather
  # than do it silently, because the silent case is the leaky one.
  if (!is.null(cov_range)) {
    if (!is.list(cov_range) || is.null(names(cov_range)))
      bad("`cov_range` must be a NAMED list, e.g. list(WT = c(52, 118)).")
    if (!all(names(cov_range) %in% nms))
      bad("`cov_range` names covariate(s) ",
          paste(sQuote(setdiff(names(cov_range), nms)), collapse = ", "),
          " that this study's `cov_dist` does not declare.")
    for (nm in names(cov_range))
      cov_dist[[nm]] <- .admCovTruncSpec(cov_dist[[nm]], cov_range[[nm]], nm)
    # Discard the derived fields FIRST: the canon short-circuits on an existing
    # `joint`, so re-running it over truncated margins was a no-op and the
    # sampler kept the full declared support. The correlation has to be handed
    # BACK as `cor` -- the first canon consumed it (cov_dist[["cor"]] <- NULL)
    # and left only `latentR` -- or dropping the derived fields would drop the
    # dependence with them.
    if (isTRUE(cov_dist[["jointOwn"]])) {
      .R <- cov_dist[["latentR"]]
      cov_dist[c("joint", "jointOwn", "latentR", "discExact")] <- NULL
      if (!is.null(.R) && length(nms) > 1L &&
          identical(dim(.R), c(length(nms), length(nms)))) {
        dimnames(.R) <- list(nms, nms)
        cov_dist[["cor"]] <- .R
      }
    }
    cov_dist <- .admCovDistCanon(cov_dist)
  }
  .no_rng <- setdiff(stratify, names(cov_range))
  if (length(.no_rng))
    warning("admixr2: stratifying on ", paste(sQuote(.no_rng), collapse = ", "),
            " over the FULL declared distribution, because no `cov_range` was ",
            "given for ", if (length(.no_rng) == 1L) "it" else "them",
            ". The source is then credited with evidence in covariate bands it ",
            "may never have enrolled -- the overstatement is ",
            "var(declared)/var(enrolled), which is 3.4x for a source spanning ",
            "+/-1 SD. Supply the reported range, e.g. `cov_range = list(",
            .no_rng[1L], " = c(min, max))`.", call. = FALSE)

  d  <- length(nms)
  iS <- match(stratify, nms)

  # === EXACT point conditioning, for a copula admixr2 BUILT ==================
  # Fixing the stratified covariates' uniforms and varying the rest only
  # conditions correctly when they lead the sampler's own cascade -- and when
  # they do not, it does not even hold them fixed, because a variable earlier
  # in the cascade feeds them. No choice of u repairs that for an OPAQUE
  # sampler; the pool route below handles that case instead.
  #
  # But when the dependence came from `cor`/`rho`/`Sigma`, admixr2 built the
  # copula and the conditional is closed form on the latent scale:
  #     z_-S | z_S ~ N(A z_S, Sig_-S - A Sig_S,-S),  A = Sig_-S,S Sig_SS^-1
  # Every draw lands on it exactly -- nothing to resample, no importance
  # weights to collapse -- and the stratified covariate is a genuine POINT.
  # Verified against the closed-form conditional to 8e-4 on the mean and
  # 14.278 vs 14.283 on the SD.
  # INDEPENDENCE IS A KNOWN LATENT STRUCTURE, not a missing one. `latentR` is
  # recorded only when `cor` was supplied, so a study that declared no
  # dependence fell to the pooled branch below -- and that branch cuts
  # EQUIPROBABLE BINS evaluated at a representative point, which is a midpoint
  # rule. The objective then converges in J at O(1/J), which is not good enough
  # to be a likelihood: measured on an independent pair, the pooled route moves
  # the objective 749 units between J = 3 and J = 25 and is still moving, while
  # the Gauss-Hermite route moves 0.2 units in total and is flat to three
  # decimals from J = 9. Same distribution, same estimates (both give the
  # exponent as 0.750); only the rule differs.
  #
  # OFV, AIC, BIC and any likelihood ratio are comparable only once this has
  # converged, so the fast rule has to be the default wherever it applies. It
  # applies whenever the latent structure is KNOWN -- an explicit `cor`, or no
  # declared dependence at all, which IS independence. An opaque user `joint`
  # is the one case that stays on the pooled branch, since nothing can be
  # conditioned there.
  Rm <- cov_dist[["latentR"]]
  if (is.null(Rm) && is.null(cov_dist[["joint"]])) Rm <- diag(length(nms))
  disc <- vapply(nms, function(n) !is.null(cov_dist[[n]][["values"]]),
                 logical(1))
  # A DISCRETE COVARIATE RIDES THIS BRANCH ONLY IF IT IS LATENTLY INDEPENDENT
  # of the rest. Then chol(R)[, j] is e_j, its level is a monotone function of
  # its own uniform, and it enumerates EXACTLY -- as a stratum when it is
  # stratified, and as exact level nodes inside every stratum when it is not.
  # Correlated with a continuous covariate it is not a point but a TRUNCATION
  # of the latent, and conditioning on it is a different construction; that
  # case keeps the pooled route below.
  #
  # Getting here matters more than it looks. Discreteness used to send the
  # whole study to the pooled route, and that is most real covariate models --
  # one declared sex, genotype or formulation. It cost 515 units of
  # J-dependence across J = 5 to 50, still moving, against 0.03 without it.
  #
  # It is not only a rate. The pooled route represents a band by an equal-weight
  # SAMPLE, and that sample is a deterministic function of cov_dist, so
  # datagen() and the fit draw the same rows. Whatever that one finite ensemble
  # happens to contain is generated into (E, V) and then read back out as if it
  # were population structure. Measured: it manufactured 18.3 units of
  # curvature in a MARGINALISED sex coefficient and returned it at its truth,
  # where the same study scored against SEX's exact declared marginal is flat
  # to 0.019 -- i.e. that coefficient is not identified at all and the pool was
  # hiding it behind its own sampling noise.
  disc_sep <- !any(disc) || (!is.null(Rm) && all(vapply(which(disc), function(j)
    all(abs(Rm[j, -j, drop = TRUE]) < 1e-12), logical(1))))
  if (!is.null(Rm) && disc_sep) {
    iSc <- iS[!disc[iS]]; iSd <- iS[disc[iS]]
    iM  <- setdiff(seq_len(d), iS)
    iMc <- iM[!disc[iM]]
    # --- the stratified block ------------------------------------------------
    # Gauss-Hermite over the continuous stratified covariates, CROSSED with the
    # levels of the discrete ones. Weights multiply: the discrete block is
    # latently independent, so the cell probability factorises exactly.
    # ROTATED BY chol(R_SS), or the bands are laid out as if the stratified
    # covariates were INDEPENDENT. .adghNodeGrid returns a product grid over
    # standard normals, and the declared correlation was read only for the
    # marginalised block's conditional mean just below -- so the code held the
    # block and did not use it here. Measured on WT/AGE at rho = 0.9 with
    # strata_nodes = 5: the weight-weighted correlation across the 25 strata was
    # 0.0000, and the (max WT, min AGE) corner -- essentially impossible under
    # N(0, R_SS) at that rho -- carried weight 1.27e-04, which .admExpandStrata
    # turns into a real study with a real n at a covariate combination the
    # declared distribution says does not occur.
    #
    # chol(R) is upper triangular with U'U = R, so Z %*% U has covariance R and
    # every COLUMN is still standard normal -- which is what keeps pnorm() the
    # right uniform for each margin's own quantile, and what makes A %*% zC the
    # conditional mean it already claims to be. The GH weights are unchanged;
    # this is the same construction the eta block uses with L.
    if (length(iSc)) {
      .ng <- .adghNodeGrid(n_nodes, length(iSc))
      zC  <- .ng$X; wC <- as.numeric(.ng$W / sum(.ng$W))
      .Rss <- Rm[iSc, iSc, drop = FALSE]
      if (length(iSc) > 1L &&
          !isTRUE(all.equal(unname(.Rss), diag(length(iSc)), tolerance = 1e-12))) {
        .U <- tryCatch(chol(.Rss), error = function(e) NULL)
        if (is.null(.U))
          stop("admixr2: the declared correlation among the stratified ",
               "covariates is not positive definite, so it describes no ",
               "distribution to band.", call. = FALSE)
        zC <- zC %*% .U
      }
    } else { zC <- matrix(0, 1L, 0L); wC <- 1 }
    nC <- nrow(zC)
    xC <- matrix(if (length(iSc)) vapply(seq_along(iSc), function(k)
      .admCovQuantile(cov_dist[[nms[iSc[k]]]],
        .admCovU(zC[, k])), numeric(nC)) else 0, nC, length(iSc))
    if (length(iSd)) {
      .dz <- .admCovDiscCells(cov_dist, nms[iSd])   # levels + probs, shared
      lvD <- .dz$lv; prD <- .dz$pr
      lgD <- as.matrix(expand.grid(lapply(lvD, seq_along), KEEP.OUT.ATTRS = FALSE))
      wD  <- Reduce(`*`, lapply(seq_along(iSd), function(k) prD[[k]][lgD[, k]]))
    } else { lvD <- list(); lgD <- matrix(0L, 1L, 0L); wD <- 1 }
    nD <- nrow(lgD)
    # --- the marginalised block ----------------------------------------------
    # Conditioning is a NO-OP wherever the latent blocks do not couple, and the
    # stratum's remaining covariates then keep their DECLARED specs -- exact,
    # with no sample standing in for them. That is the common case, and it is
    # what lets .admCovGrid enumerate a discrete margin at its levels and put
    # Gauss-Hermite nodes on a continuous one.
    A <- if (length(iMc) && length(iSc))
      Rm[iMc, iSc, drop = FALSE] %*% solve(Rm[iSc, iSc, drop = FALSE]) else
      matrix(0, length(iMc), length(iSc))
    exact_marg <- !length(A) || all(abs(A) < 1e-12)
    Rmm <- Rm[iM, iM, drop = FALSE]
    dep_marg <- length(iM) > 1L &&
      !isTRUE(all.equal(unname(Rmm), diag(length(iM)), tolerance = 1e-12))
    # a stratified covariate is held at a POINT, carried as a degenerate spec so
    # covStrata() still shows it and every consumer sees the same name set
    # `.point` is what .admCovDistDegenerate() reads. A stratum keeps the
    # stratified covariate in its own cov_dist ON PURPOSE -- a stratum is a
    # range, and downstream code needs every covariate in scope -- so presence
    # of cov_dist cannot mean "there is something to marginalise".
    pt_spec <- function(v) list(.point = TRUE, quantile = local({ v <- as.numeric(v)
      function(u) rep(v, length.out = length(u)) }))
    if (!exact_marg) {
      Sc <- Rm[iMc, iMc, drop = FALSE] - A %*% Rm[iSc, iMc, drop = FALSE]
      Lc <- t(chol(Sc + diag(1e-12, nrow(Sc))))
      nd <- max(as.integer(n_pool), 8192L)
      q  <- randtoolbox::sobol(nd, dim = length(iMc))
      if (!is.matrix(q)) q <- matrix(q, ncol = length(iMc))
      Z0 <- stats::qnorm(pmin(pmax(q, 1e-12), 1 - 1e-12))
    }
    return(unlist(lapply(seq_len(nD), function(kd) lapply(seq_len(nC), function(kc) {
      cd_k <- list()
      for (k in seq_along(iSc)) cd_k[[nms[iSc[k]]]] <- pt_spec(xC[kc, k])
      for (k in seq_along(iSd)) cd_k[[nms[iSd[k]]]] <- pt_spec(lvD[[k]][lgD[kd, k]])
      if (exact_marg) {
        for (j in iM) cd_k[[nms[j]]] <- cov_dist[[nms[j]]]
        if (dep_marg) {
          Rk <- diag(length(nms)); Rk[iM, iM] <- Rmm
          dimnames(Rk) <- list(nms, nms)
          cd_k[["cor"]] <- Rk[names(cd_k), names(cd_k), drop = FALSE]
        }
      } else {
        # The closed-form Gaussian conditional of the CONTINUOUS block, drawn
        # once as a low-discrepancy sample. Discrete margins do not go through
        # it: they are independent of the conditioning, so they keep their
        # exact specs and are mapped from their own uniform, which is what
        # `discExact` tells .admCovGrid it may enumerate.
        mu <- as.numeric(A %*% zC[kc, ])
        Xc <- sweep(Z0 %*% t(Lc), 2L, mu, "+")
        Xc <- matrix(vapply(seq_along(iMc), function(k)
          .admCovQuantile(cov_dist[[nms[iMc[k]]]],
            .admCovU(Xc[, k])), numeric(nd)), nd, length(iMc),
          dimnames = list(NULL, nms[iMc]))
        # sorted before it is handed out: the rows arrive in Sobol order and the
        # sampler that draws from them runs another Sobol stream, so indexing
        # rows by that stream correlates the two and returns a structured
        # subsample (measured: sd 13.094 against a true 14.283)
        Xc <- Xc[order(Xc[, 1L]), , drop = FALSE]
        for (k in seq_along(iMc)) cd_k[[nms[iMc[k]]]] <-
          list(quantile = local({ v <- sort(Xc[, k])
            function(u) v[1L + pmin(floor(u * length(v)), length(v) - 1L)] }))
        for (j in setdiff(iM, iMc)) cd_k[[nms[j]]] <- cov_dist[[nms[j]]]
        nk  <- names(cd_k)
        cn  <- nms[iMc]
        jc  <- match(cn[1L], nk)
        cd_k[["joint"]] <- local({ Xc <- Xc; nk <- nk; jc <- jc; cn <- cn
          cds <- cd_k[nk]; other <- setdiff(seq_along(nk), match(cn, nk))
          function(u) {
            i <- 1L + pmin(floor(u[, jc] * nrow(Xc)), nrow(Xc) - 1L)
            out <- matrix(0, nrow(u), length(nk), dimnames = list(NULL, nk))
            out[, cn] <- Xc[i, , drop = FALSE]
            for (k in other) out[, k] <- .admCovQuantile(cds[[k]], u[, k])
            out
          }})
        cd_k[["discExact"]] <- setdiff(nk, cn)
      }
      cv <- c(stats::setNames(as.list(xC[kc, ]), nms[iSc]),
              stats::setNames(lapply(seq_along(iSd), function(k)
                lvD[[k]][lgD[kd, k]]), nms[iSd]))
      list(cov = cv[nms[iS]], cov_dist = cd_k,
           weight = wC[kc] * wD[kd],
           n_pool_cell = if (exact_marg) NA_integer_ else nd)
    })), recursive = FALSE))
  }

  # --- the pool -------------------------------------------------------------
  # ONE deterministic draw from the study's own covariate distribution, through
  # exactly the sampler the fit will use. Every stratum is then a SUBSET of
  # this pool, so the strata cannot describe a different population from the
  # one that was declared, and their weights are counts rather than an
  # assumption.
  #
  # Sobol, so this is a fixed function of `cov_dist` -- the covariate
  # distribution is DATA and must not move between objective evaluations.
  #
  # THE COST OF THIS ROUTE IS THE TAIL. A stratum resamples pool members, so no
  # stratum can produce a covariate value beyond the pool's own extremes, and
  # the outermost strata are described by the fewest draws. A Sobol pool of
  # 32768 reaches about the 99.997th percentile of a lognormal margin, so the
  # truncation is far out -- but it is real, and it is why the pool is sized by
  # the CELL count below rather than by the population, and why a thin cell is
  # an error rather than a shrug.
  # THE POOL MUST BE SIZED FOR THE CELLS, NOT THE POPULATION. Every stratum is
  # a subset, so a pool that is ample overall can still be thin in a corner
  # cell -- and it is the corner cells that carry the covariate extremes the
  # coefficient is estimated from. Two covariates at 4 bins each is 16 cells;
  # three is 64. Scale with the cell count and keep a floor per cell.
  n_cell_max <- prod(vapply(seq_along(iS), function(k) {
    v <- cov_dist[[nms[iS[k]]]][["values"]]
    if (is.null(v)) as.integer(n_nodes) else length(v) }, integer(1)))
  n_pool <- max(as.integer(n_pool), 4096L * n_cell_max)
  pool <- .admCovRowsFor(cov_dist, n_pool, 0L)
  pool <- pool[, nms, drop = FALSE]

  # --- the bins -------------------------------------------------------------
  # A discrete covariate is cut at its levels; a continuous one into
  # EQUIPROBABLE bins, so every stratum carries the same number of subjects and
  # the weights need no quadrature theory to justify. Cutting on the pool's own
  # quantiles also means a skewed covariate gets narrow bins where it is dense.
  cutone <- function(j) {
    x <- pool[, j]
    if (!is.null(cov_dist[[nms[j]]][["values"]])) {
      lv <- sort(unique(x))
      return(list(id = match(x, lv), k = length(lv)))
    }
    br <- stats::quantile(x, probs = seq(0, 1, length.out = n_nodes + 1L),
                          names = FALSE, type = 7)
    br[1L] <- -Inf; br[length(br)] <- Inf
    br <- unique(br)
    if (length(br) < 3L)
      bad("covariate ", sQuote(nms[j]), " takes too few distinct values to ",
          "cut into ", n_nodes, " strata.")
    list(id = as.integer(cut(x, breaks = br, labels = FALSE,
                             include.lowest = TRUE)), k = length(br) - 1L)
  }
  cuts <- lapply(iS, cutone)
  key  <- Reduce(function(a, b) (a - 1L) * b$k + b$id,
                 cuts[-1L], cuts[[1L]]$id)
  present <- sort(unique(key))
  if (length(present) < 2L)
    bad("stratifying on ", paste(sQuote(stratify), collapse = ", "),
        " produced a single stratum, so there is no contrast to gain from it.")

  # A cell resamples from the pool members that landed in it, so a thin cell
  # resolves its own conditional distribution coarsely -- and the thin cells are
  # the EXTREME ones, which is where a covariate coefficient gets its leverage.
  # This is the same failure the SIR samplers report as a collapsed effective
  # sample size. Reported rather than silently tolerated.
  cnt <- tabulate(match(key, present), length(present))
  if (min(cnt) < 200L)
    bad("stratifying on ", paste(sQuote(stratify), collapse = ", "), " into ",
        length(present), " strata leaves only ", min(cnt), " of ", n_pool,
        " pooled draws in the smallest one, which is too few to resolve its ",
        "covariate distribution -- and the thin strata are the extreme ones, ",
        "where the coefficient gets its leverage.
",
        "  Use fewer strata (`strata_nodes`), stratify on fewer covariates, ",
        "or raise `n_pool`.")

  # --- one stratum per non-empty cell --------------------------------------
  # The stratum's covariate distribution is the pool RESTRICTED to its cell --
  # the empirical conditional, carrying every covariate jointly, so all of the
  # dependence survives with no assumption about the sampler's structure.
  #
  # This is the rejection/SIR idea used to sample vines conditionally, applied
  # to the sampler's OUTPUT rather than to a density: admixr2 is handed a
  # sampler, not a `dvinecop`, so it cannot importance-weight, but it can bin
  # what the sampler produced. The earlier route -- fixing the stratified
  # covariates' INPUT uniforms -- only works when they lead the sampler's own
  # conditioning cascade, and failed silently otherwise: on a vine ordered
  # AGE -> WT -> CRCL, stratifying on WT left AGE at its unconditional 55.09 in
  # every stratum against a true 50.1 to 60.4. Binning has no such requirement.
  lapply(present, function(kk) {
    idx <- which(key == kk)
    sub <- pool[idx, , drop = FALSE]
    # a deterministic resampler over the cell: whole ROWS, so the joint
    # dependence within the stratum is carried exactly. One uniform column
    # indexes; a Sobol stream therefore walks the cell evenly.
    cd_k <- stats::setNames(lapply(nms, function(n) {
      list(quantile = local({ v <- sort(sub[, n])
        function(u) v[1L + pmin(floor(u * length(v)), length(v) - 1L)] })
      )}), nms)
    # sorted for the same reason as the exact branch above: the pool arrives in
    # Sobol order and is indexed by another Sobol stream
    jsort <- setdiff(seq_len(d), iS); if (!length(jsort)) jsort <- 1L
    sub <- sub[order(sub[, jsort[1L]]), , drop = FALSE]
    cd_k[["joint"]] <- local({ sub <- sub; nms <- nms
      function(u) {
        i <- 1L + pmin(floor(u[, 1L] * nrow(sub)), nrow(sub) - 1L)
        out <- sub[i, , drop = FALSE]; colnames(out) <- nms; out
      }})
    list(cov = stats::setNames(as.list(colMeans(sub[, iS, drop = FALSE])),
                               nms[iS]),
         cov_dist = cd_k,
         weight = length(idx) / nrow(pool),
         n_pool_cell = length(idx))
  })
}

#' Cut a covariate distribution into strata
#'
#' Shows the strata `datagen()` builds for a study declaring `stratify` --- the
#' covariate value each stratum is pinned at, its effective sample size, and the
#' distribution the remaining covariates follow *within* it. It is the way to
#' see what a stratification actually describes before a fit depends on it.
#'
#' Stratify on the covariates a source's own published model conditions on, and
#' leave the rest to be marginalised. A source that never fitted a covariate has
#' no contrast in it to report, and evaluating it at nodes that vary that
#' covariate manufactures a null one --- measured, with no sampling noise, as an
#' attenuation from 0.450 to 0.211 in the fitted coefficient.
#'
#' @param cov_dist A covariate specification, as given to a study --- see
#'   [covDraw()] for the full grammar.
#' @param stratify Character vector naming the covariates to stratify on. The
#'   rest are left in each stratum's own `cov_dist`, to be marginalised over
#'   their distribution **conditional** on that stratum.
#' @param n_nodes Strata per stratified covariate (default 9). A discrete
#'   covariate ignores it and is cut at its levels, exactly, with the level
#'   probabilities as weights.
#'
#'   For a continuous one it sets how finely the covariate is resolved, and the
#'   rule depends on what is known about the joint distribution. Where the
#'   latent structure is known --- an explicit `cor`, or no declared dependence,
#'   which IS independence --- each stratum is a Gauss-Hermite node, held at a
#'   POINT, and the remaining covariates follow their exact conditional law
#'   within it. An opaque `joint` sampler cannot be conditioned, and there the
#'   covariate is cut into that many equiprobable BINS instead, each described
#'   by the pool members that landed in it.
#'
#'   **This is a convergence parameter, not a modelling choice.** The
#'   well-defined object is the limit as the count grows; any finite value
#'   approximates it, and under misspecification the answer can jump between
#'   basins rather than drift. Raise it until the estimates stop moving rather
#'   than picking a value you like. Cost is `n_nodes^p` studies for `p`
#'   stratified covariates, so this is cheap in one covariate and expensive in
#'   three.
#'
#'   How much the objective depends on it turns on which rule applies. Where the
#'   latent structure is known — an explicit `cor`, or no declared dependence,
#'   which is independence — banding conditions on a Gauss-Hermite grid and the
#'   objective is stable: 0.03 units across counts of 5 to 100, so objective,
#'   AIC, BIC and likelihood ratios are comparable across resolutions.
#'
#'   A declared discrete covariate no longer forces the slow rule: as long as it
#'   is latently independent of the others it enumerates exactly at its levels,
#'   which took a study declaring one sex covariate from 515 units of movement,
#'   still rising at a count of 50, to 0.009.
#'
#'   An opaque `joint` sampler cannot be conditioned and falls back to
#'   equiprobable bins, where it is not stable: 51 units at a count of 100 and
#'   still moving, with the estimate itself wandering (0.681 at 5, 9, 15 and
#'   100; 0.697 at 50). There, raise the count until the estimates settle, and
#'   note that `anova()` refuses to compare two fits built at different ones.
#' @param cov_range Optional named list giving the range each stratified
#'   covariate was actually ENROLLED over, e.g. `list(WT = c(52, 118))` —
#'   publications routinely report a min-max or a median with an IQR.
#'
#'   Without it the strata are cut over the whole declared distribution, and a
#'   source is credited with evidence in covariate bands it never sampled. The
#'   overstatement is exactly `var(declared) / var(enrolled)`: 3.4x for a source
#'   spanning ±1 SD, 12.4x at ±0.5 SD. Estimates stay correct — what inflates is
#'   confidence, and it only shows once a second source disagrees. Omitting it
#'   warns for that reason.
#'
#'   More strata resolve the covariate range more finely but do not buy
#'   accuracy, and each one costs a solve: on a matched one-covariate fit the
#'   coefficient came back at 0.7000 / 0.7002 / 0.7005 / 0.7005 for 3 / 4 / 10 /
#'   16 strata against a true 0.700, while the fit took 7.2 / 4.1 / 8.0 / 11.5
#'   seconds. Raise it when the covariate effect is strongly nonlinear over the
#'   range; the only hard limit is that a stratum must keep enough pooled draws
#'   to describe itself, which admixr2 checks.
#' @param n Total sample size to divide among the strata. The default of `1`
#'   returns the raw stratum weights.
#' @param n_pool Minimum size of the deterministic pool the strata are cut
#'   from. Each stratum's covariate distribution is that pool restricted to its
#'   own bin, so this sets how finely the within-stratum distribution is
#'   resolved --- and, because no stratum can produce a covariate value beyond
#'   the pool's extremes, how far into the tails it reaches. The pool is grown
#'   automatically with the number of strata, and a stratum left with too few
#'   draws is an error rather than a silently coarse answer; `n_pool_cell` in
#'   the result reports what each one got.
#'
#'   It applies only where a pool is used at all: an opaque `joint` sampler, or
#'   a stratified covariate latently correlated with a marginalised one. Where
#'   the conditional law is closed form the stratum carries the DECLARED specs
#'   themselves, no sample stands in for them, and `n_pool_cell` is `NA`.
#'
#' @return A list with one element per stratum, each containing `cov` (the
#'   value each stratified covariate is held at --- a quadrature node, a level,
#'   or a bin mean), `cov_dist` (the distribution of ALL covariates *within*
#'   that stratum, carrying the full joint dependence), `n` and `weight`.
#'
#' @examples
#' # a source that fitted weight but not renal function, the two correlated
#' st <- covStrata(list(WT   = list(mean = 72, sd = 15),
#'                      CRCL = list(mean = 90, sd = 22), cor = 0.7),
#'                 stratify = "WT", n_nodes = 4L, n = 300)
#' vapply(st, function(s) s$cov$WT, numeric(1))     # where each stratum sits
#' vapply(st, function(s) s$n, numeric(1))          # and how big it is
#'
#' # within a stratum, CRCL follows its CONDITIONAL distribution: a heavier
#' # stratum has a higher creatinine clearance, which is the whole point
#' vapply(st, function(s) mean(covDraw(s$cov_dist, n = 2000L)[, "CRCL"]),
#'        numeric(1))
#'
#' @seealso [covDraw()] to inspect a covariate specification, [datagen()] to
#'   generate the strata as studies.
#' @export
covStrata <- function(cov_dist, stratify, n_nodes = .ADM_STRATA_NODES, n = 1,
                      n_pool = 32768L, cov_range = NULL) {
  checkmate::assertNumber(n, lower = 0, finite = TRUE)
  st <- .admCovStrata(cov_dist, stratify, n_nodes, n_pool, cov_range)
  lapply(st, function(s)
    list(cov = s$cov, cov_dist = s$cov_dist, n = n * s$weight,
         weight = s$weight, n_pool_cell = s$n_pool_cell))
}

# Which ESTIMATED thetas parameterise a covariate's effect?
#
# `allCovs` reports which covariates a model READS, not which coefficients it
# ESTIMATED, and the difference decides whether that source carries any evidence
# about the covariate at all. A model containing `(WT/70)^0.75`, or
# `clwt <- fix(0.75)`, reads WT while ASSERTING its coefficient -- the
# allometric convention, so this is the common case, not a corner one. Banding
# such a source credits it with evidence it never earned: the information it
# actually contributes about the covariate is ZERO, so the ratio of claimed to
# earned information is unbounded at every stratum count.
#
# Structural parsing cannot answer this reliably -- the coefficient may be an
# exponent, a multiplier, a slope inside a link, a spline knot -- and the
# tempting shortcut of "which thetas appear in the same assignment" is wrong for
# the commonest form of all: in `cl <- exp(tcl + eta.cl) * (WT/70)^0.75` the
# theta `tcl` shares the expression with WT and has nothing to do with it.
#
# So it is answered NUMERICALLY, which is exact for any expression: theta
# parameterises the covariate's effect iff the MIXED second difference
#
#     [f(cov+h, th+d) - f(cov, th+d)] - [f(cov+h, th) - f(cov, th)]
#
# is non-zero -- that is, iff changing theta changes what the covariate DOES.
# A literal exponent gives exactly zero; an estimated one does not.
.admCovCoefThetas <- function(ui, cov, cov_dist = NULL, tol = 1e-8) {
  lst <- tryCatch(ui$lstExpr, error = function(e) NULL)
  ini <- tryCatch(ui$iniDf,   error = function(e) NULL)
  if (is.null(lst) || is.null(ini)) return(NULL)
  th <- ini[is.na(ini$neta1) & is.na(ini$err) & !ini$fix, , drop = FALSE]
  if (!nrow(th)) return(character(0))
  etas <- unique(stats::na.omit(ini$name[!is.na(ini$neta1)]))
  is_asgn <- vapply(lst, function(e) is.call(e) && length(e) == 3L &&
                      (identical(e[[1L]], quote(`<-`)) ||
                       identical(e[[1L]], quote(`=`))), logical(1))
  hit <- which(is_asgn & vapply(lst, function(e)
    is.call(e) && length(e) == 3L && cov %in% all.vars(e[[3L]]), logical(1)))
  if (!length(hit)) return(character(0))
  # a nominal value for every covariate: the declared median where there is one
  covs <- tryCatch(ui$allCovs, error = function(e) character(0))
  base_cov <- stats::setNames(lapply(covs, function(nm) {
    sp <- if (!is.null(cov_dist)) cov_dist[[nm]] else NULL
    v  <- if (!is.null(sp)) tryCatch(.admCovQuantile(sp, 0.5),
                                     error = function(e) NULL) else NULL
    if (is.null(v) || !is.finite(v) || v == 0) 1 else as.numeric(v)
  }), covs)
  ev_at <- function(cov_val, th_over) {
    ev <- new.env(parent = asNamespace("rxode2"))
    for (i in seq_len(nrow(ini)))
      if (is.na(ini$neta1[i])) assign(ini$name[i], ini$est[i], ev)
    for (nm in names(th_over)) assign(nm, th_over[[nm]], ev)
    for (e in etas) assign(e, 0, ev)
    for (nm in covs) assign(nm, base_cov[[nm]], ev)
    assign(cov, cov_val, ev)
    out <- numeric(0)
    for (ii in seq_along(lst)) {
      if (!isTRUE(is_asgn[ii])) next
      v <- tryCatch(eval(lst[[ii]][[3L]], ev), error = function(e) NULL)
      if (is.null(v)) next
      assign(as.character(lst[[ii]][[2L]]), v, ev)
      if (ii %in% hit && length(v) == 1L && is.finite(v)) out <- c(out, v)
    }
    out
  }
  x0 <- base_cov[[cov]] %||% 1
  x1 <- x0 * 1.05
  f00 <- ev_at(x0, list()); f10 <- ev_at(x1, list())
  if (!length(f00) || length(f00) != length(f10)) return(character(0))
  # ON THE LOG SCALE, and that is the whole test rather than a detail. These
  # models are multiplicative, so a SCALE parameter does change the covariate's
  # ABSOLUTE effect -- in cl = exp(tcl + eta) * (WT/70)^0.75, raising tcl
  # raises the number of L/h a 10% weight change buys. Differencing f would
  # therefore flag tcl as WT's coefficient, which it is not. What tcl leaves
  # alone is the RELATIVE effect: d(log f)/d(log WT) is 0.75 whatever tcl is,
  # and 1 * bwt when the exponent is estimated. So the mixed difference is
  # taken on log f, and only a theta that moves the covariate's PROPORTIONAL
  # effect counts as parameterising it.
  lg <- function(a, b, c_, d_) {           # log where every arm is positive
    if (all(c(a, b, c_, d_) > 0)) list(log(a), log(b), log(c_), log(d_))
    else list(a, b, c_, d_)
  }
  keep <- vapply(seq_len(nrow(th)), function(i) {
    nm <- th$name[i]
    ov <- stats::setNames(list(th$est[i] + 1e-3 * max(abs(th$est[i]), 1)), nm)
    f01 <- ev_at(x0, ov); f11 <- ev_at(x1, ov)
    if (length(f01) != length(f00) || length(f11) != length(f00)) return(FALSE)
    z <- lg(f00, f10, f01, f11)
    d0 <- z[[2L]] - z[[1L]]; d1 <- z[[4L]] - z[[3L]]
    any(abs(d1 - d0) > tol * pmax(abs(d0), 1e-12))
  }, logical(1))
  th$name[keep]
}

# Expand every study carrying `stratify` into one ordinary study per stratum.
#
# The output is plain studies -- own `n`, own `cov`, own `cov_dist` -- so the
# generator and the estimator both treat them as unrelated sources, which is
# exactly what they are once built. That is also the safety property: there is
# no way to reach the invalid construction (one pooled observation scored
# against several conditional predictions) through this path, because every
# stratum gets its OWN observation generated for it.
#
# `n_k = w_k * n` is the stratum's effective size. It is deliberately not
# rounded: it is a quadrature weight standing in for a stratum size, and
# rounding it would break sum(n_k) = n. A source that publishes REAL subgroups
# reports their real sizes, which beat quadrature weights -- pass those as
# ordinary studies instead of using `stratify` at all.
.admExpandStrata <- function(studies, study_names, model = NULL) {
  has <- vapply(studies, function(s) !is.null(s[["stratify"]]), logical(1))
  if (!any(has)) return(list(studies = studies, names = study_names))
  out <- list(); nms <- character(0)
  for (i in seq_along(studies)) {
    s <- studies[[i]]; nm <- study_names[[i]]
    if (!isTRUE(has[[i]])) { out[[length(out) + 1L]] <- s; nms <- c(nms, nm); next }
    # `stratify = TRUE` derives the split FROM THE SOURCE MODEL, which is the
    # only thing that actually knows it: a covariate is stratified when this
    # study's own data-generating model conditions on it, and marginalised when
    # it does not. Asking the user to restate that is asking them to duplicate
    # information already in the model, and the failure mode of getting it
    # wrong -- stratifying on a covariate the source never fitted -- is the
    # fabricated null contrast, measured at 0.450 -> 0.211.
    if (isTRUE(s[["stratify"]])) {
      m <- s[["model"]] %||% model
      # AN rxUi COUNTS, and admStudy() only ever supplies one: it parses the
      # model at construction so the transcription can be checked, and hands
      # the ui down. rxode2::rxode2() is idempotent on a ui -- the next line
      # calls it either way -- and datagen() accepts both for the same reason,
      # so demanding a function here rejected the whole admStudy() route.
      if (!is.function(m) && !inherits(m, "rxUi"))
        stop("admixr2: study '", nm, "' asks for `stratify = TRUE`, which is ",
             "derived from that study's own data-generating model, but no ",
             "`model` was supplied for it.", call. = FALSE)
      cvs <- tryCatch(suppressMessages(rxode2::rxode2(m))$allCovs,
                      error = function(e) NULL)
      if (is.null(cvs))
        stop("admixr2: study '", nm, "' asks for `stratify = TRUE` but its ",
             "model could not be parsed to find which covariates it reads.",
             call. = FALSE)
      keep <- intersect(.admCovSpecNames(s[["cov_dist"]]), cvs)
      if (!length(keep))
        stop("admixr2: study '", nm, "' asks for `stratify = TRUE`, but its ",
             "data-generating model reads none of the covariates its ",
             "`cov_dist` declares", if (length(cvs))
               paste0(" (model reads: ", paste(cvs, collapse = ", "), ")"),
             ". There is no contrast to stratify on; drop `stratify`, and the ",
             "study is generated marginal over the distribution instead.",
             call. = FALSE)
      # READING a covariate is not the same as having ESTIMATED its
      # coefficient. A source that asserted it -- `(WT/70)^0.75`, or a fix()ed
      # theta -- carries no evidence about that covariate, and banding on it
      # credits information that was never earned. Band only where a free theta
      # actually modulates the covariate; a partially-fixed set keeps the
      # estimated members and drops the rest.
      .ui_s <- tryCatch(suppressMessages(rxode2::rxode2(m)),
                        error = function(e) NULL)
      if (!is.null(.ui_s)) {
        .est <- vapply(keep, function(cv) length(.admCovCoefThetas(
          .ui_s, cv, s[["cov_dist"]])) > 0L, logical(1))
        if (!any(.est))
          stop("admixr2: study '", nm, "' asks for `stratify = TRUE`, but its ",
               "model ASSERTS the coefficient of ",
               paste(sQuote(keep), collapse = ", "),
               " rather than estimating it (a fixed exponent or a fix()ed ",
               "theta). A source that fixed a covariate's coefficient carries ",
               "no evidence about that covariate, so there is no contrast to ",
               "extract; drop `stratify` and the study is generated marginal ",
               "over the distribution instead.", call. = FALSE)
        if (!all(.est))
          message("admixr2: study '", nm, "': not stratifying on ",
                  paste(sQuote(keep[!.est]), collapse = ", "),
                  " -- the model asserts ",
                  if (sum(!.est) == 1L) "its coefficient" else "their coefficients",
                  " rather than estimating ",
                  if (sum(!.est) == 1L) "it." else "them.")
        keep <- keep[.est]
      }
      s[["stratify"]] <- keep
    }
    if (is.null(s[["cov_dist"]]))
      stop("admixr2: study '", nm, "' declares `stratify` but no `cov_dist`. ",
           "Stratifying needs the covariate distribution the strata are cut ",
           "from.", call. = FALSE)
    if (!is.null(s[["observations"]]))
      stop("admixr2: study '", nm, "' declares both `stratify` and ",
           "`observations`. A multi-output study is already a list of blocks ",
           "and nothing downstream defines the product of the two.",
           call. = FALSE)
    n_tot <- s[["n"]]
    if (is.null(n_tot) || !is.finite(n_tot) || n_tot <= 0)
      stop("admixr2: study '", nm, "' declares `stratify` but has no positive ",
           "`n` to divide among the strata.", call. = FALSE)
    stl <- .admCovStrata(s[["cov_dist"]], s[["stratify"]],
                         s[["strata_nodes"]] %||% .ADM_STRATA_NODES,
                         cov_range = s[["cov_range"]])
    # `stratify = character(0)` reached .admCovStrata's `if (!length(stratify))
    # return(NULL)` and the loop below then ran zero times, DROPPING the study
    # -- the one malformed `stratify` that did not error. Every other spelling
    # is refused loudly, so this was refused silently.
    if (!length(stl))
      stop("admixr2: study '", nm, "' declares `stratify` but it names no ",
           "covariate, so the study would be dropped rather than banded. Give ",
           "a covariate name, or remove `stratify`.", call. = FALSE)
    .Jk <- s[["strata_nodes"]] %||% .ADM_STRATA_NODES
    for (k in seq_along(stl)) {
      sk <- s
      sk[["stratify"]] <- NULL; sk[["strata_nodes"]] <- NULL
      sk[["cov_range"]] <- NULL
      # carried so the fit can refuse to be compared against one built at a
      # different resolution -- the objective is J-dependent
      sk[[".adm_strata_nodes"]] <- .Jk
      # a stratum's own point value for the stratified covariates, merged over
      # any `cov` the study already set for covariates it does not stratify on
      sk[["cov"]] <- utils::modifyList(as.list(s[["cov"]] %||% list()),
                                       stl[[k]]$cov)
      sk[["cov_dist"]] <- stl[[k]]$cov_dist          # NULL if all stratified
      sk[["n"]] <- n_tot * stl[[k]]$weight
      out[[length(out) + 1L]] <- sk
      nms <- c(nms, sprintf("%s_s%d", nm, k))
    }
  }
  names(out) <- nms
  list(studies = out, names = nms)
}

# Product grid over several covariates: every combination, weights multiplied.
#
# DEPENDENT COVARIATES. A `joint` sampler maps INDEPENDENT uniforms to
# dependent covariate values -- that is what a copula is, and for a vine the
# uniforms are its Rosenblatt coordinates, which are independent by
# construction WHATEVER the dependence. So a product grid in u-space is exact
# for any joint distribution: the weights genuinely factorise there, and the
# dependence is carried entirely by the map.
#
# That is the SAME construction the independent branch already uses -- it
# pushes pnorm(z) through each margin's own quantile function -- with the per
# margin quantiles replaced by the joint map. Nothing about the grid, the
# weights or the row stride changes, so every consumer downstream is unaffected.
# A SMOLYAK SPARSE GRID over the continuous margins, crossed with the discrete
# ones enumerated exactly.
#
# WHY THIS REPLACED THE "taylor" DESIGN. That design was derived as a
# second-order moment expansion, but it is a cubature rule and its own comment
# recorded which one: at the moment-matched radius h = sqrt(3 lambda) it
# "coincides with 3-point Gauss-Hermite". Measured, at p = 1 it was
# BIT-IDENTICAL to .admCovGrid at n_nodes = 3, and above that it was the axial
# subset of the 3^p product grid. In cubature terms it was exactly Smolyak
# LEVEL 2 -- so the way to make it more accurate is the next LEVEL, not a
# bigger radius and not more nodes along each axis.
#
# Three measurements settle that (lognormal margins, allometric + saturable
# integrand, reference = product-15/21):
#
#   * MORE NODES PER AXIS SATURATES. axial-3 -> 5 -> 7 -> 9 is unchanged from 5
#     onward, and from 3 onward once the covariates are correlated: what is
#     left is the mixed terms an axial rule structurally cannot see.
#   * LEVEL 3 IS A LARGE AND CHEAP GAIN. At p = 4, rho = 0.85 it is 49 points
#     against the product grid's 81 and is 44x more accurate on the mean and
#     33x on the covariance -- cheaper AND better. At p = 3 it matches
#     product-3's cost at ~30x the accuracy.
#   * CORRELATION FLIPS SIGN BETWEEN THE LEVELS. Level 2 gets WORSE with it
#     (relE 9.3e-05 at rho = 0 against 1.9e-03 at 0.5); level 3 gets BETTER
#     (6.6e-07 against 4.8e-08). The cheap rule that handles correlation well
#     is real -- it is level 3, and it was never level 2.
#
# The price is SIGNED weights. sum(W) is exactly 1 at every level, but sum|W|
# is 1.0 at level 2, 2.5 at level 3 for p = 3 and 4.1 for p = 4, so the answer
# is a difference of terms several times its own size and any solver noise in a
# design point is amplified by that factor. .admSandwichCov checks Om for
# indefiniteness because of it. Level 2 is not innocent either: its centre
# weight is already negative at p = 4 (min -1/3).
#
# Returns list(X, W, z) exactly as .admCovGrid does, so every consumer -- the
# moments, the omega chain, the ADF weight, the shift certificate -- takes it
# unchanged. That is the whole reason this is a GRID rather than a design
# carrying its own derivative machinery, which is what the Taylor path was and
# where two of its defects lived.
.ADM_SPARSE_GROWTH <- c(1L, 3L, 5L, 7L, 9L)

.admCovSparseGrid <- function(cov_dist, level = 3L, n_nodes = NULL) {
  cov_dist <- .admCovDistCanon(cov_dist)
  nms <- .admCovSpecNames(cov_dist)
  d   <- length(nms)
  if (!d) return(NULL)
  level <- as.integer(level)
  if (!is.finite(level) || level < 2L || level > length(.ADM_SPARSE_GROWTH))
    stop("admixr2: `cov_sparse_level` must be an integer between 2 and ",
         length(.ADM_SPARSE_GROWTH), "; got ", format(level), ".",
         call. = FALSE)
  # Discrete margins are enumerated at their levels and crossed in, exactly as
  # the product grid does. A sparse rule is a statement about the CONTINUOUS
  # dimensions; a level probability is not an approximation of anything.
  disc <- vapply(nms, function(n) !is.null(cov_dist[[n]][["values"]]),
                 logical(1))
  cn <- nms[!disc]; dn <- nms[disc]
  dc <- length(cn)
  if (!dc) return(.admCovGrid(cov_dist, n_nodes %||% 7L))
  # An opaque `joint` sampler is refused for the reason the Taylor design was:
  # the canoniser early-returns before recording `latentR`, so the rotation
  # below would read a dependent distribution as independent.
  Rz <- cov_dist[["latentR"]]
  if (is.null(Rz) && is.function(cov_dist[["joint"]]) && dc > 1L)
    stop("admixr2: cov_integration = \"sparse\" cannot integrate a covariate ",
         "distribution supplied as a `joint` sampler -- the sparse rule needs ",
         "the latent correlation and an opaque closure does not report one. ",
         "Declare it with `cor` (or `rho`/`Sigma`), which admixr2 builds its ",
         "own sampler from, or use cov_integration = \"quadrature\".",
         call. = FALSE)
  # The refusal, the correlation block and the discrete enumeration all come
  # from .admCovLatentBlock, which the two collapses use. Writing them out here
  # is what let the REFUSAL be omitted from the first version of this function:
  # a discrete margin latently correlated with a continuous one was integrated
  # as if independent until a test caught it.
  .lb <- .admCovLatentBlock(cov_dist, nms, cn, dn, Rz)
  if (is.null(.lb))
    stop("admixr2: this covariate distribution cannot be integrated on a ",
         "sparse grid -- either a DISCRETE covariate is latently correlated ",
         "with a continuous one (a level is then a truncation of the latent ",
         "normal rather than a point, so the continuous conditional differs ",
         "from cell to cell and one shared design is the wrong design in ",
         "every cell), or the declared correlation is not positive definite. ",
         "Use cov_integration = \"quadrature\", or declare the discrete ",
         "covariate independent of the continuous ones.", call. = FALSE)
  Rc <- .lb$Rc
  # Rotate onto the eigenvectors of the latent correlation, so the rule runs
  # along the directions the distribution actually varies in. At Rc = I this is
  # the identity and the grid is the ordinary axis-aligned one.
  Arot <- diag(1, dc)
  if (dc > 1L && max(abs(Rc[lower.tri(Rc)])) > 0) {
    ez <- eigen(Rc, symmetric = TRUE)
    if (min(ez$values) <= -1e-8)
      stop("admixr2: the covariate correlation implied by `cov_dist` is not ",
           "positive semi-definite (smallest eigenvalue ",
           format(min(ez$values)), "), so it describes no distribution.",
           call. = FALSE)
    Arot <- ez$vectors %*% diag(sqrt(pmax(ez$values, 0)), dc)
  }
  sm <- .admSparseNodes(dc, level)
  Z  <- sm$X %*% t(Arot)
  Xc <- .admCovXFromZ(cov_dist, cn, Z)
  if (!all(is.finite(Xc)))
    stop("admixr2: the sparse grid produced non-finite covariate values -- a ",
         "margin's quantile function saturated at a tail node. Lower ",
         "`cov_sparse_level`, or give a margin that is finite there.",
         call. = FALSE)
  Wc <- sm$W
  if (length(dn)) {
    nq <- nrow(Xc); nl <- nrow(.lb$cells)
    ix <- rep(seq_len(nq), times = nl); il <- rep(seq_len(nl), each = nq)
    X  <- cbind(Xc[ix, , drop = FALSE],
                .lb$cells[il, , drop = FALSE])[, nms, drop = FALSE]
    W  <- Wc[ix] * .lb$pcell[il]
    Zo <- Z[ix, , drop = FALSE]
  } else {
    X <- Xc[, nms, drop = FALSE]; W <- Wc; Zo <- Z
  }
  list(X = X, W = W / sum(W), z = Zo)
}

# The Smolyak combination rule over `d` standard-normal dimensions.
#
#   A(d, L) = sum over multi-indices i with L-d+1 <= |i| <= L of
#             (-1)^(L-|i|) * choose(d-1, L-|i|) * (U^{i_1} x ... x U^{i_d})
#
# with L = d + level - 1 and U^k the .ADM_SPARSE_GROWTH[k]-point Gauss-Hermite
# rule. Level 2 reproduces the axial rule exactly, which is what makes the
# retired Taylor design a special case rather than a separate method.
#
# Points are merged BY VALUE, and that is why the centre node has to be snapped
# to exactly zero: .adghNodes1 returns it as ~1e-16 out of the eigen
# decomposition, so without the snap the same point arrives under several keys
# and the grid comes back with spurious rows carrying split weights (measured:
# 7 rows where the rule has 5, at p = 2 level 2).
.admSparseNodes <- function(d, level) {
  L   <- d + level - 1L
  idx <- as.matrix(expand.grid(rep(list(seq_len(level)), d),
                               KEEP.OUT.ATTRS = FALSE))
  idx <- idx[rowSums(idx) >= L - d + 1L & rowSums(idx) <= L, , drop = FALSE]
  gs  <- lapply(.ADM_SPARSE_GROWTH[seq_len(level)], function(n) {
    g <- .adghNodes1(n); g$x[abs(g$x) < 1e-12] <- 0; g })
  env <- new.env(hash = TRUE, parent = emptyenv())
  for (r in seq_len(nrow(idx))) {
    i  <- idx[r, ]
    cf <- (-1)^(L - sum(i)) * choose(d - 1L, L - sum(i))
    if (cf == 0) next
    gk <- gs[i]
    gg <- as.matrix(expand.grid(lapply(gk, function(g) seq_along(g$x)),
                                KEEP.OUT.ATTRS = FALSE))
    for (t in seq_len(nrow(gg))) {
      x  <- vapply(seq_len(d), function(j) gk[[j]]$x[gg[t, j]], numeric(1))
      w  <- prod(vapply(seq_len(d), function(j) gk[[j]]$w[gg[t, j]],
                        numeric(1)))
      ky <- paste(sprintf("%.12g", x), collapse = "|")
      env[[ky]] <- (env[[ky]] %||% 0) + cf * w
    }
  }
  ks <- ls(env)
  W  <- vapply(ks, function(k) env[[k]], numeric(1))
  ok <- abs(W) > 1e-14
  X  <- do.call(rbind, lapply(ks[ok], function(k)
    as.numeric(strsplit(k, "|", fixed = TRUE)[[1L]])))
  list(X = matrix(X, sum(ok), d), W = unname(W[ok]))
}

# Call a user-supplied `joint` sampler and hold it to its contract.
#
# Two sites in .admCovGrid did this identically: the discExact branch, which
# crosses an enumerated discrete margin with the sampler, and the plain joint
# branch. The contract is a DECISION -- one row per supplied uniform, named
# columns, exactly the declared name set -- so a copy that validated less would
# let a malformed sampler through on one path and not the other.
.admCovJointEval <- function(jf, u, nms) {
  d <- length(nms)
  X <- tryCatch(as.matrix(jf(u)), error = function(e)
    stop("admixr2: cov_dist$joint failed on the ", nrow(u), " x ", d,
         " grid of uniforms: ", conditionMessage(e), call. = FALSE))
  if (!is.matrix(X) || nrow(X) != nrow(u) || is.null(colnames(X)) ||
      !setequal(colnames(X), nms))
    stop("admixr2: cov_dist$joint must return one ROW per supplied uniform ",
         "row, with columns ", paste(sQuote(nms), collapse = ", "), ".",
         call. = FALSE)
  X[, nms, drop = FALSE]
}

.admCovGrid <- function(cov_dist, n_nodes) {
  cov_dist <- .admCovDistCanon(cov_dist)
  nms <- .admCovSpecNames(cov_dist)
  jf  <- cov_dist[["joint"]]
  if (is.function(jf)) {
    d  <- length(nms)
    # A DISCRETE margin cannot ride the Gauss-Hermite grid once the covariates
    # are DEPENDENT. The sampler mixes the uniforms before mapping them to
    # levels -- admixr2's own Gaussian copula computes qnorm(u) %*% L and then
    # pnorm() per margin -- so fixing a discrete covariate's INPUT uniform does
    # not fix its OUTPUT level, and the cell weights are not the level
    # probabilities. Enumerating the slabs anyway was tried and measured: a
    # covariate declared 0.55/0.45 came off the grid at 0.477, and the error
    # does NOT shrink with cov_nodes because it is a property of the mixing.
    #
    # So that case switches integration rule rather than losing the covariate:
    # an equal-weight low-discrepancy pool drawn through the sampler itself,
    # which needs no smoothness and reproduces levels and dependence together.
    # It converges more slowly than Gauss-Hermite on a smooth integrand, which
    # is why it is used only where GH is inapplicable.
    disc <- vapply(nms, function(n) !is.null(cov_dist[[n]][["values"]]),
                   logical(1))
    # ... unless the sampler maps that margin straight from its own uniform,
    # which admixr2 records when it builds the copula. See .admCovDiscExact().
    # intersected with `disc` deliberately: a stratum carries its stratified
    # covariates as DEGENERATE specs and lists them in `discExact` too, and
    # those have no `values` to enumerate.
    ex   <- disc & nms %in% (cov_dist[["discExact"]] %||% character(0))
    if (any(disc & !ex)) {
      # THE POOL IS A LAST RESORT, AND IT IS NOT JUST SLOWER -- IT IS SHARED.
      # An equal-weight pool is a deterministic function of `cov_dist`, so
      # datagen() and the fit draw the SAME rows: any idiosyncrasy of that one
      # finite ensemble is generated into (E, V) and then read back out as if
      # it were population structure. Measured on a study banded on WT with an
      # independent SEX declared, that manufactured 18.3 units of curvature in
      # the sex coefficient and "recovered" it at its truth, while the same
      # study fitted against SEX's exact declared marginal is flat to 0.019 --
      # i.e. the coefficient is not identified at all and the pool was hiding
      # it. Every margin that can be enumerated must be.
      npool <- max(as.integer(n_nodes)^d, 4096L)
      X <- .admCovRowsFor(cov_dist, npool, 0L)[, nms, drop = FALSE]
      # no latent normal score exists for a discrete pool, and none is wanted:
      # a discrete covariate can never make Delta normal.
      return(list(X = X, W = rep(1 / nrow(X), nrow(X)), z = NULL))
    }
    if (any(ex)) {
      # Cross the joint grid over the remaining margins with an EXACT
      # enumeration of the separable discrete ones. A pass-through column is
      # monotone in its own uniform, so feeding the MIDPOINT of a level's
      # probability interval selects exactly that level -- verified below
      # rather than assumed, since `discExact` is a claim about a closure.
      iE <- which(ex); iJ <- which(!ex)
      .dz <- .admCovDiscCells(cov_dist, nms[iE])    # levels + probs, shared
      lv <- .dz$lv; pr <- .dz$pr
      lg <- as.matrix(expand.grid(lapply(lv, seq_along), KEEP.OUT.ATTRS = FALSE))
      wE <- Reduce(`*`, lapply(seq_along(iE), function(k) pr[[k]][lg[, k]]))
      uE <- vapply(seq_along(iE), function(k) {
        cp <- cumsum(pr[[k]]); (cp - pr[[k]] / 2)[lg[, k]] }, numeric(nrow(lg)))
      uE <- matrix(uE, nrow(lg), length(iE))
      if (length(iJ)) {
        .ng <- .adghNodeGrid(n_nodes, length(iJ))
        zJ  <- .ng$X; wJ <- .ng$W / sum(.ng$W)
      } else { zJ <- matrix(0, 1L, 0L); wJ <- 1 }
      nJ <- nrow(zJ); nE <- nrow(lg)
      u  <- matrix(0, nJ * nE, d, dimnames = list(NULL, nms))
      rj <- rep(seq_len(nJ), times = nE)      # joint block fastest
      re <- rep(seq_len(nE), each  = nJ)
      if (length(iJ))
        u[, iJ] <- .admCovU(zJ[rj, , drop = FALSE])
      u[, iE] <- uE[re, , drop = FALSE]
      X <- .admCovJointEval(jf, u, nms)
      want <- vapply(seq_along(iE), function(k) lv[[k]][lg[re, k]],
                     numeric(nrow(u)))
      if (!isTRUE(all.equal(unname(X[, iE, drop = FALSE]),
                            unname(matrix(want, nrow(u), length(iE))))))
        stop("admixr2: cov_dist$discExact names ",
             paste(sQuote(nms[iE]), collapse = ", "), ", but the sampler does ",
             "not map ", if (length(iE) == 1L) "it" else "them",
             " from the matching input uniform, so the enumerated levels are ",
             "not the ones it returns. Drop `discExact`.", call. = FALSE)
      if (!all(is.finite(X)))
        stop("admixr2: cov_dist$joint returned non-finite covariate values on ",
             "the quadrature grid.", call. = FALSE)
      # z is the latent score behind the CONTINUOUS block only; a discrete
      # covariate has none, and can never make Delta normal.
      return(list(X = X, W = as.numeric(wJ[rj] * wE[re]), z = NULL))
    }
    .ng <- .adghNodeGrid(n_nodes, d)
    z   <- .ng$X
    W   <- .ng$W
    u  <- stats::pnorm(z)
    u  <- pmin(pmax(u, .Machine$double.eps), 1 - .Machine$double.eps)
    colnames(u) <- nms
    X <- .admCovJointEval(jf, u, nms)
    if (!all(is.finite(X))) {
      # Almost always saturation, not a broken sampler. The grid's extreme
      # nodes reach |z| ~ 6.4 at 15 nodes, a copula's mixing step scales that
      # by up to (rho + sqrt(1-rho^2)) ~ 1.4, and pnorm() of the result rounds
      # to exactly 1 -- after which qlnorm(1) is Inf. The cure is to clamp
      # AFTER pnorm, not only on the way in, which is what admixr2's own
      # `cor` sampler does.
      nb <- which(!is.finite(rowSums(X)))
      stop("admixr2: cov_dist$joint returned non-finite covariate values at ",
           length(nb), " of ", nrow(X), " quadrature nodes (first at u = ",
           paste(sprintf("%.3e", u[nb[1L], ]), collapse = ", "), ").\n",
           "  The grid's tail nodes are more extreme than a per-subject ",
           "sample, and a copula amplifies them further, so a sampler that ",
           "does `qnorm(u)` then `pnorm(...)` can saturate at exactly 0 or 1 ",
           "and return +/-Inf from its quantile function.\n",
           "  Clamp inside the sampler AFTER the pnorm step, e.g.\n",
           "    v <- pmin(pmax(pnorm(z), 1e-12), 1 - 1e-12)\n",
           "  or declare the dependence with `cor` and let admixr2 build the ",
           "sampler, or lower `cov_nodes`.", call. = FALSE)
    }
    return(list(X = X, W = W / sum(W), z = z))
  }
  one <- lapply(nms, function(n) .admCovNodesFor(cov_dist[[n]], n_nodes))
  X   <- as.matrix(expand.grid(lapply(one, `[[`, "x"), KEEP.OUT.ATTRS = FALSE))
  W   <- Reduce(`*`, lapply(seq_along(one), function(k)
           rep(rep(one[[k]]$w, each = prod(vapply(one[seq_len(k - 1L)],
                 function(o) length(o$x), integer(1)))),
               length.out = nrow(X))))
  colnames(X) <- nms
  # The latent normal scores behind X, assembled with the SAME stride. Each
  # continuous margin is a standard-normal node pushed through Phi and then the
  # margin's quantile function, so `z` is what X is an image of -- which is what
  # certifies Delta's law (see .admShiftAffineResid).
  #
  # A margin given by `values` is discrete and has no score, but it does not
  # void the others: the certificate regresses Delta on the columns that DO
  # exist, so a Delta reaching a discrete covariate simply fails to be affine in
  # them and is refused, while one that does not reach it still certifies. That
  # matters for the ordinary case of a continuous covariate on the shifted
  # parameter beside a discrete one somewhere else in the model -- returning
  # NULL for the whole grid sent it down the mixture route for no reason.
  # Expanded over EVERY margin, so the stride matches X row for row, then
  # narrowed to the continuous columns. Expanding only the kept margins gives a
  # z with fewer rows than X, and .admShiftAffineResid then declines on the
  # length check -- which looks like a refusal to certify rather than a bug.
  keep <- !vapply(nms, function(n) !is.null(cov_dist[[n]][["values"]]),
                  logical(1))
  z <- if (!any(keep)) NULL else {
    ixz <- as.matrix(expand.grid(lapply(one, function(o) seq_along(o$x)),
                                 KEEP.OUT.ATTRS = FALSE))
    do.call(cbind, lapply(which(keep), function(k)
      .adghNodes1(length(one[[k]]$x))$x[ixz[, k]]))
  }
  list(X = X, W = W / sum(W), z = z)
}


# =============================================================================
# Shift path -- the covariate dimension leaves the solver entirely
# =============================================================================
#
# When the covariates influence the model ONLY through a mu-referenced argument
# they act as a pure shift of that argument's random effect:
#
#     f(a, eta) == f(a_ref, eta + Delta(a))
#
# so the (a, eta) integral collapses onto an integral over u = Delta(a) + eta.
# The covariate never reaches the solver: it is held at its reference and the
# affected eta column carries the covariate's whole contribution.  The solve
# therefore costs n_u * (nodes for the OTHER etas) rows instead of
# n_node^n_eta * n_cov^p, which is CONSTANT in the number of covariates.
#
# Two properties make this cheap as well as exact:
#
#  * Delta is a function of the structural thetas and the covariates only -- no
#    ODE, no eta, no time -- so it is evaluated by running the model's own
#    parameter assignment in R over the covariate nodes as VECTORS.  Zero solves.
#  * u's law is the mixture sum_j w_j N(Delta_j, omega^2), whose mean and
#    variance are known exactly (E[Delta], Var(Delta) + omega^2), so its
#    quantiles are found by Newton from a moment-matched normal start rather
#    than by a grid or a bracketing solve.
#
# THE PRECONDITION IS CHECKED NUMERICALLY, NEVER READ OFF THE MODEL TEXT.
# .admShiftVerify() evaluates the identity above against the compiled model at
# several etas and covariate values; the path is granted only if it holds to
# tolerance.  This is the positive check the earlier `uq` route lacked -- its
# four silent-wrong-answer modes (additive effect, non-log link, two covariates
# per eta, discrete covariate) were all consequences of inferring the property
# from syntax.  Measured separation on the probe set: valid forms 1e-12..1e-14,
# invalid forms 3e-02..7e-01.

# Which assignment carries the covariates, and which eta shares it?
# Returns NULL unless EXACTLY one eta appears with them -- m > 1 (a covariate on
# two mu-referenced parameters) is representable but needs a vector u, which
# this path does not build.
.admShiftSpec <- function(ui, cov_names, eta_names) {
  lst <- tryCatch(ui$lstExpr, error = function(e) NULL)
  if (is.null(lst) || !length(cov_names)) return(NULL)
  # a covariate inside an if() is invisible to the scan below -- see
  # .admCovInBranch()
  if (.admCovInBranch(lst, cov_names)) return(NULL)
  hit <- list()
  for (e in lst) {
    if (!is.call(e) || length(e) < 3L) next
    if (!as.character(e[[1L]])[1L] %in% c("<-", "=")) next
    v <- all.vars(e[[3L]])
    if (!any(cov_names %in% v)) next
    et <- intersect(eta_names, v)
    # an assignment carrying a covariate but NO random effect has nothing to
    # shift; two random effects in one assignment cannot be separated either
    if (length(et) != 1L) return(NULL)
    hit[[length(hit) + 1L]] <- list(param = as.character(e[[2L]])[1L],
                                    eta = et, rhs = e[[3L]])
  }
  if (!length(hit)) return(NULL)
  # More than one assignment may carry covariates -- a weight term on both CL
  # and V, say. Each contributes its own shifted column, so m = the number of
  # DISTINCT affected random effects; two assignments sharing one eta cannot.
  ets <- vapply(hit, `[[`, "", "eta")
  if (anyDuplicated(ets)) return(NULL)
  # The link of the mu-referenced argument, so Delta can be measured on the
  # scale eta lives on. rxode2 reports this per THETA; blank means it could not
  # classify the expression, and `exp` is the pharmacometric default that the
  # numerical verification will reject if wrong.
  ce <- tryCatch(ui$muRefCurEval, error = function(e) NULL)
  md <- tryCatch(ui$muRefDataFrame, error = function(e) NULL)
  lks <- vapply(ets, function(et) {
    lk <- ""
    if (!is.null(ce) && !is.null(md) && NROW(md)) {
      th <- md$theta[md$eta == et]
      if (length(th) && th[1L] %in% ce$parameter)
        lk <- as.character(ce$curEval[ce$parameter == th[1L]])[1L]
    }
    if (!is.null(lk) && nzchar(lk)) lk else "exp" }, "")
  list(eta = ets, link = unname(lks),
       param = vapply(hit, `[[`, "", "param"),
       rhs = lapply(hit, `[[`, "rhs"))
}

# Delta at a matrix of covariate values, relative to `aref`. ONE vectorised
# evaluation of the model's own assignment -- the node set goes through as
# columns, not row by row.
.admShiftDelta <- function(spec, struct, X, aref) {
  m  <- length(spec$eta)
  # The rxode2 namespace is the EVAL PARENT here, not a source of named
  # internals -- `spec$rhs` is the user's own model expression and rxode2 is the
  # scope it was written against, so expit/logit/probit and anything else the
  # assignment names resolve the way they do in the model. That is the reason
  # CLAUDE.md's "grep asNamespace( as well as :::" rule does not bite: the rule
  # is about REACHING INTO a namespace for a specific unexported function, which
  # is what plot.R used to do. Narrowing this to an explicit list of exports
  # would silently move any model naming something outside that list onto the
  # product grid.
  ev <- new.env(parent = asNamespace("rxode2"))
  for (k in names(struct)) assign(k, struct[[k]], ev)
  for (e in spec$eta) assign(e, 0, ev)
  out <- matrix(0, nrow(X), m, dimnames = list(NULL, spec$eta))
  for (k in seq_len(m)) {
    inv <- switch(spec$link[k],
                  exp = , log = log,
                  expit = , logit = rxode2::logit,
                  probitInv = , probit = rxode2::probit,
                  identity)
    for (j in names(aref)) assign(j, aref[[j]], ev)
    p0 <- eval(spec$rhs[[k]], ev)
    # as.list() on a matrix row drops dimnames, so assign COLUMNS by name --
    # anything else silently leaves every node at the reference and Delta is
    # identically zero, which is finite, plausible and wrong.
    for (j in colnames(X)) assign(j, as.numeric(X[, j]), ev)
    d <- inv(eval(spec$rhs[[k]], ev)) - inv(p0)
    if (length(d) != nrow(X) || !all(is.finite(d))) return(NULL)
    out[, k] <- as.numeric(d)
  }
  out
}

# Quantiles of u = Delta + eta at Gauss-Hermite probabilities.
#
# u's law is the mixture sum_j W_j N(Delta_j, om^2). Its mean and variance are
# exact, so Newton starts from the moment-matched normal and converges in a few
# steps; each step is one n_u x n_node matrix. A grid was tried and is slower
# (512 x n_node evaluations of pnorm against n_u x n_node per Newton step).
# Beyond ~8 standard deviations pnorm() saturates at exactly 0 or 1, so the
# target probability carries no information and the quantile there is not
# determined. Those nodes have Gauss-Hermite weight ~1e-23 and cannot move a
# moment, so they are left as Newton leaves them rather than special-cased.
# Is the weighted law of Delta itself NORMAL?
#
# When it is, u = Delta(a) + eta is normal too -- the convolution of two
# normals -- with variance Var(Delta) + omega^2, and its quantiles are the
# closed form that the Newton inversion below merely converges to. That is the
# retired "collapse" path, recovered here as a special case of the shift rather
# than as a separate route with its own preconditions.
#
# The test is on (Delta, W) ALONE: no model text, no muRefCovariateDataFrame.
# That matters twice over. rxode2 splits covariate forms across THREE frames --
# `tcov*WT` lands in muRefCovariateDataFrame, `tcov*WT/70` in muRefExtra and
# `tcov*log(WT/70)` in mu2RefCovariateReplaceDataFrame (measured on 5.1.4) --
# so any single frame sees a third of the cases and the union is still a
# syntactic guess. And the property that actually matters is not "the covariate
# is normal" but "Delta(a) is normal", which admits the allometric case a
# distributional test catches for free: for lognormal WT, `tcov*log(WT/70)` is
# affine in the latent normal score and so exactly qualifies.
#
# Standardised central moments 3..6 against N(0,1). Delta arrives on a
# Gauss-Hermite rule, which reproduces a normal's moments to machine precision
# up to degree 2n-1, so a qualifying Delta scores ~1e-13 and a non-qualifying
# one ~1e0: measured separation over the probe set is 9.9e-14 (normal+linear),
# 7.3e-14 (lognormal+log) against 2.6e+00 (normal+square), 2.8e+00
# (normal+log), 5.9e+00 (lognormal+linear) -- thirteen orders, so the threshold
# is not a tuned quantity.
.ADM_SHIFT_GAUSS_TOL <- 1e-8

.admShiftGaussResid <- function(D, W) {
  W  <- W / sum(W)
  m  <- sum(W * D)
  s2 <- sum(W * (D - m)^2)
  # A CONSTANT Delta (no covariate spread) is degenerate-normal: u is then just
  # eta and the closed form is exact, so admit it rather than dividing by zero.
  # A NON-FINITE one is not: scoring it 0 would route NaN through the closed
  # form and return a node set full of NaN wearing a valid shape.
  if (!is.finite(s2)) return(Inf)
  if (s2 <= 0) return(0)
  zz  <- (D - m) / sqrt(s2)
  mom <- vapply(3:6, function(k) sum(W * zz^k), numeric(1))
  max(abs(mom - c(0, 3, 0, 15)) / c(1, 3, 1, 15))
}




# Is Delta an AFFINE image of the latent normal scores?
#
# This is the certificate, and it replaces asking whether Delta's own moments
# look normal. admixr2 builds every continuous covariate as X = F^-1(Phi(z))
# from a jointly normal z, so an affine Delta = c + B z is exactly normal, and
# jointly so across coordinates -- an affine image of a Gaussian is Gaussian in
# any dimension. That covers both cases worth having by construction rather
# than by inspection: a normal covariate entering linearly, and a LOGNORMAL one
# entering through log(), where log(F^-1(Phi(z))) is affine in z.
#
# A moment test on Delta cannot do this job above one dimension. Normality of
# every margin plus finitely many fixed projections does not certify joint
# normality -- Cramer-Wold needs ALL projections -- so it can only ever fail to
# find a counterexample. Affinity is checkable, and sufficient.
#
# It is one-sided: a jointly normal Delta that is not affine in z is refused and
# takes the mixture route, which is correct but slower. That is the safe
# direction, and the two forms this exists for are both affine.
.admShiftAffineResid <- function(D, W, z) {
  if (is.null(z)) return(Inf)
  D <- as.matrix(D); z <- as.matrix(z)
  if (nrow(z) != nrow(D) || !all(is.finite(z)) || !all(is.finite(D)))
    return(Inf)
  Wn <- W / sum(W); sw <- sqrt(Wn)
  A  <- cbind(1, z) * sw
  q  <- tryCatch(qr(A), error = function(e) NULL)
  if (is.null(q) || q$rank < ncol(A)) return(Inf)
  r <- 0
  for (k in seq_len(ncol(D))) {
    m   <- sum(Wn * D[, k])
    sd0 <- sqrt(max(sum(Wn * (D[, k] - m)^2), 0))
    if (sd0 <= 0) next                  # constant column is affine trivially
    y   <- D[, k] * sw
    r   <- max(r, sqrt(sum((y - qr.fitted(q, y))^2)) / sd0)
  }
  r
}

# One predicate, so the objective and the gradient cannot key off different
# tests. `m` is the shift dimension: at m = 1 a moment test on Delta is a valid
# fallback when no latent score is to hand; above that only the affine
# certificate will do.
.admShiftGaussOK <- function(D, W, z, m) {
  if (!is.null(z) &&
      .admShiftAffineResid(D, W, z) < .ADM_SHIFT_GAUSS_TOL) return(TRUE)
  if (m > 1L) return(FALSE)
  D <- as.matrix(D)
  .admShiftGaussResid(D[, 1L], W) < .ADM_SHIFT_GAUSS_TOL
}










# Every structural theta by name -- ESTIMATED ones at their current value,
# fix()ed ones at the value in iniDf. Delta is evaluated in a plain R
# environment, so a fix()ed covariate coefficient that is absent from
# `struct` would be looked up in the search path instead of the model.
.admShiftStruct <- function(pinfo, struct = NULL) {
  out <- list()
  d <- pinfo$iniDf
  if (!is.null(d) && NROW(d)) {
    keep <- is.na(d$neta1) & (is.na(d$err) | d$err == "")
    out <- as.list(stats::setNames(as.numeric(d$est[keep]),
                                   as.character(d$name[keep])))
  }
  cur <- if (is.null(struct)) pinfo$struct_init else struct
  for (k in names(cur)) out[[k]] <- unname(cur[[k]])
  out
}

# The reference covariate values the shift is measured against: each
# covariate's own mean, which is where the study is solved.
.admShiftRef <- function(cov_dist, cov_names) {
  r <- lapply(cov_names, function(n) .admCovMeanOf(cov_dist[[n]]))
  names(r) <- cov_names
  if (any(vapply(r, function(z) is.null(z) || !is.finite(z), logical(1))))
    return(NULL)
  r
}

# =============================================================================
# Covariate spec moments, and the identifiability warning
# =============================================================================

# Which parameter assignment does a covariate enter, and which eta shares it?
#
# The ONLY consumer is .admWarnCovIdentifiability(), which needs to know whether
# a covariate sits on the flat (theta, omega, beta) ridge -- it does exactly when
# it shares a mu-referenced argument with a random effect. That is a question
# about whether a warning applies, so a syntactic answer is adequate: the cost of
# being wrong is a missing or spurious warning.
#
# It is NOT adequate for ROUTING, and must never be used for it again. The
# retired `uq` path made this decision from the same syntax and was silently
# wrong in four measured ways (see the file header). .admShiftSpec() +
# .admShiftVerify() are the routing pair: same question, answered against the
# compiled model.
.admCovParamEta <- function(ui, cov, eta_names) {
  lst <- tryCatch(ui$lstExpr, error = function(e) NULL)
  if (is.null(lst)) return(NULL)
  asg <- Filter(function(e)
    is.call(e) && length(e) >= 3L &&
      as.character(e[[1L]])[1L] %in% c("<-", "="), lst)
  reads <- vapply(asg, function(e) cov %in% all.vars(e[[3L]]), logical(1))
  # The ridge argument requires the model to see the covariate ONLY through the
  # mu-referenced sum: it is what makes the prediction a function of
  # `gamma*a + omega*b` rather than of `(a, b)` separately. A covariate read by
  # a SECOND assignment -- a second parameter, or a scaling the design applies
  # -- restores the separate dependence and is identified from one population,
  # so there is no ridge and no warning to give. Reading only the first
  # assignment that paired the covariate with an eta declared such a model
  # unidentified when it is not.
  if (sum(reads) != 1L) return(NULL)
  e  <- asg[[which(reads)]]
  et <- intersect(eta_names, all.vars(e[[3L]]))
  if (length(et) != 1L) return(NULL)
  list(param = as.character(e[[2L]])[1L], eta = et)
}

# The covariate value the model is SOLVED at when the study does not name one.
#
# The shift path measures Delta(a) relative to it (.admShiftRef), and any study
# mixing a fixed `cov` with a distributed one needs a value for the fixed part.
# Deriving it from `cov_dist` rather than making the user restate it removes a
# way for the two to disagree silently.
.admCovMeanOf <- function(spec) {
  if (!is.null(spec$values)) {
    pr <- spec$probs %||% rep(1, length(spec$values))
    return(sum(as.numeric(spec$values) * pr) / sum(pr))
  }
  if (is.function(spec$quantile))
    return(mean(spec$quantile((seq_len(1024L) - 0.5) / 1024L)))
  if (!is.null(spec$meanlog)) return(exp(spec$meanlog + spec$sdlog^2 / 2))
  spec$mu
}

# Var(a) on the covariate's OWN scale -- which is NOT .admCovSdOf()^2.
#
# .admCovSdOf reports the spread "on the scale the identifiability ridge lives
# on", so for a lognormal covariate it returns `sdlog`, a LOG-scale number,
# while .admCovMeanOf returns the natural-scale mean. Squaring it and pairing it
# with that mean would expand about the right point with the wrong second
# moment. The Taylor expansion needs both moments on the scale the model reads
# the covariate on, so it gets its own function rather than a caller remembering
# which of the two conventions it is holding.
.admCovVarOf <- function(spec) {
  if (!is.null(spec$values)) {
    pr <- spec$probs %||% rep(1, length(spec$values)); pr <- pr / sum(pr)
    v  <- as.numeric(spec$values); m <- sum(v * pr)
    return(sum(pr * (v - m)^2))
  }
  # The POPULATION variance of the same 1024 midpoints .admCovMeanOf averages,
  # not stats::sd()^2: these are quadrature nodes for E_a[.], not a sample.
  if (is.function(spec$quantile)) {
    q <- as.numeric(spec$quantile((seq_len(1024L) - 0.5) / 1024L))
    return(mean((q - mean(q))^2))
  }
  if (!is.null(spec$meanlog))
    return((exp(spec$sdlog^2) - 1) * exp(2 * spec$meanlog + spec$sdlog^2))
  as.numeric(spec$sd)^2
}

# Spread of a covariate spec, on the scale the identifiability ridge lives on.
.admCovSdOf <- function(spec) {
  if (!is.null(spec$values)) {
    pr <- spec$probs %||% rep(1, length(spec$values)); pr <- pr / sum(pr)
    m <- sum(as.numeric(spec$values) * pr)
    return(sqrt(sum(pr * (as.numeric(spec$values) - m)^2)))
  }
  if (is.function(spec$quantile)) {
    # POPULATION sd of the 1024 midpoints, not stats::sd() -- these are
    # quadrature nodes for E_a[.], not a sample, and the `values` branch above
    # already uses the population form.
    q <- spec$quantile((seq_len(1024L) - 0.5) / 1024L)
    return(sqrt(mean((q - mean(q))^2)))
  }
  if (!is.null(spec$meanlog)) return(spec$sdlog)
  spec$sd
}


# Warn when a covariate coefficient cannot be identified from the data supplied.
#
# When a covariate shares a mu-referenced argument with a random effect, the
# model sees only u = Delta(a) + eta, so ONE population determines just two
# quantities -- u's mean and variance -- against three parameters. The likelihood
# is then EXACTLY flat along
#
#     theta'   = theta + (b - b')*mu_a
#     omega'^2 = omega^2 + (b^2 - b'^2)*sd_a^2
#
# (verified: the objective is bit-identical across b from 0.40 to 1.10). Only
# BETWEEN-STUDY variation in the covariate distribution breaks it -- differing
# means break the first equation, differing spreads the second. Note this is
# variation in the DATA, not a between-study random effect: admixr2 has none, and
# a tau^2 would in fact compete with the covariate for the same signal.
#
# A covariate on a parameter with NO random effect is not affected: there is no
# omega for its variance to be absorbed into, so it is identified by shape.
.admWarnCovIdentifiability <- function(.ui, pinfo, studies) {
  # CANONICALISE FIRST. This is the one entry point that reads the RAW study
  # list -- every driver calls it before normalising, deliberately -- so it is
  # also the one that sees the user's shorthand un-expanded. `mean`/`sd` has no
  # branch in .admCovMeanOf, which returned NULL and made the mean NA, so
  # between-study variation was invisible and the "not identifiable" warning
  # fired on exactly the contrast that identifies the coefficient.
  studies <- .admCovAsList(lapply(studies, function(s) {
    if (is.list(s) && !is.null(s[["cov_dist"]]))
      s[["cov_dist"]] <- tryCatch(.admCovDistCanon(s[["cov_dist"]]),
                                  error = function(e) s[["cov_dist"]])
    s
  }))
  # A study declares a covariate one of two ways: as the distribution it is
  # marginalised over (`cov_dist`), or as the value it is conditioned at
  # (`cov`). These are the same object at two resolutions -- a conditioned
  # stratum is a degenerate distribution -- and BOTH break the ridge, because
  # each supplies its own equation in (theta, gamma). Reading `cov_dist` alone
  # warned that a source reporting by stratum was unidentified, which is the
  # one case that identifies the coefficient WITHIN a single population.
  .decl <- function(s, cv) {
    d <- s[["cov_dist"]][[cv]]
    if (!is.null(d))
      return(c(.admCovMeanOf(d) %||% NA_real_, .admCovSdOf(d) %||% NA_real_))
    v <- s[["cov"]][[cv]]
    if (!is.null(v) && is.numeric(v) && length(v) == 1L && is.finite(v))
      return(c(as.numeric(v), 0))
    NULL
  }
  cvs <- unique(unlist(lapply(studies, function(s)
    c(.admCovSpecNames(s[["cov_dist"]]), names(s[["cov"]])))))
  for (cv in cvs) {
    # only covariates that share an argument with an eta sit on the ridge
    if (is.null(.admCovParamEta(.ui, cv, pinfo$eta_col_names))) next
    sp <- Filter(Negate(is.null), lapply(studies, .decl, cv = cv))
    if (length(sp) == 0L) next
    mu <- vapply(sp, `[`, numeric(1), 1L)
    sd <- vapply(sp, `[`, numeric(1), 2L)
    varies <- (length(unique(signif(mu, 10))) > 1L) ||
              (length(unique(signif(sd, 10))) > 1L)
    # THE RIDGE IS FLAT ONLY WHERE u = Delta(a) + eta IS NORMAL, which by
    # Cramer is exactly where Delta is -- the same certificate the shift path
    # routes on. Then u's law is pinned by (beta + gamma mu_a,
    # gamma^2 sigma_a^2 + omega^2) and any (gamma, omega) holding that pair
    # fixed predicts identically. For a covariate that does NOT certify, u is a
    # MIXTURE, whose law is not determined by its first two moments, and f is
    # nonlinear -- so the aggregate V separates the pair and the coefficient is
    # identified from ONE pooled source. Measured along the ridge, -2LL moving
    # from its centre:
    #
    #   normal              0.000000  exactly flat, machine precision
    #   binary {0,1}        4.9  3.5  80.1
    #   4-level discrete    2.0  1.4  19.8
    #   lognormal, raw      167  63   115
    #
    # Warning on those told the user a design could not identify something it
    # identifies strongly -- which sends them hunting for a second source they
    # do not need, or to discard the analysis. Sex, formulation, food status
    # and genotype are all in that group and all routinely reported as one
    # pooled summary, so this was not a corner case.
    if (!varies && !isTRUE(.admCovRidgeFlat(.ui, pinfo, studies, cv)))
      varies <- TRUE
    if (!varies)
      warning("admixr2: the coefficient on covariate '", cv, "' is not ",
              "identifiable from these data. It enters the same argument as a ",
              "random effect, and every study declaring it has the SAME ",
              "covariate distribution, so the likelihood is exactly flat along ",
              "a trade-off between that coefficient, the corresponding fixed ",
              "effect and omega. Identification needs sources whose covariate ",
              "MEANS or SPREADS differ, or one source reporting its summary ",
              "metrics BY covariate stratum.", call. = FALSE)
  }
  invisible(NULL)
}


#' Draw covariate values from a `cov_dist` specification
#'
#' Simulates the covariate values admixr2 itself would use for `n` subjects,
#' from the same specification a study carries in `cov_dist`. It is the way to
#' see what a specification actually describes --- to check that a reported mean
#' and standard deviation were transcribed correctly, that a correlation points
#' the way round you meant, or that a copula sampler returns what you think it
#' does --- before a fit depends on it.
#'
#' @param cov_dist A covariate specification, as given to a study. Each element
#'   names a covariate and describes its distribution, either as `mean` and `sd`
#'   on the covariate's own scale --- NORMAL by default, so pass
#'   `dist = "lnorm"` for the positive margin an allometric term needs (only
#'   [admPopulation()] defaults to lognormal) --- or as a `quantile` function,
#'   or as `values` (with optional
#'   `probs`) for a discrete covariate. A `cor` entry --- a scalar for two
#'   covariates, or a correlation matrix --- links them through a Gaussian
#'   copula. A `joint` function takes the matrix of uniforms and returns one
#'   column per covariate, which is how an arbitrary copula, an R-vine included,
#'   is supplied.
#' @param n Number of subjects to draw.
#' @param n_eta Number of random effects in the model the specification belongs
#'   to. The draws are deterministic and come from Sobol dimensions after the
#'   random effects', so passing the model's own `n_eta` reproduces exactly the
#'   values a fit would use. The default of `0` is right for inspecting a
#'   specification on its own.
#'
#' @return A numeric matrix with `n` rows and one named column per covariate.
#'
#' @examples
#' # a baseline-characteristics table, transcribed directly
#' X <- covDraw(list(WT   = list(mean = 72, sd = 16),
#'                   CRCL = list(mean = 90, sd = 25),
#'                   cor  = 0.6), n = 500)
#' colMeans(X)
#' stats::cor(log(X[, "WT"]), log(X[, "CRCL"]))
#' @export
covDraw <- function(cov_dist, n = 1000L, n_eta = 0L) {
  if (!is.list(cov_dist) || !length(cov_dist))
    stop("admixr2: `cov_dist` must be a non-empty list naming each covariate.",
         call. = FALSE)
  .admCovRowsFor(cov_dist, as.integer(n), as.integer(n_eta))
}

# =============================================================================
# covDist() -- the user-facing constructor
# =============================================================================
#
# What a user HAS is a baseline-characteristics table: a covariate name, a mean
# and an SD per row, sometimes a correlation. What the internals consume is a
# nested list of canonical specs. covDist() is the bridge, and it exists so
# that the bridge is crossed ONCE, at construction, where an error can name the
# covariate -- rather than at the first objective evaluation, where it cannot.
#
# A NAMED numeric vector is required for the two-number forms, deliberately.
# `WT = c(72, 16)` would be terser, but `SEX = c(0, 1)` is then ambiguous
# between "mean 0, sd 1" and "the levels 0 and 1" -- and the two differ by a
# factor no fit can detect. Names make the intent explicit and cost four
# characters.

# LIFTED OUT OF R/study-api.R, where it also serves admPopulation().
# .admCovSpecFromVec() below calls it, so `covDist(WT = c(mean = 72,
# cv = 22))` -- the documented spelling -- does not load without it. The
# study API lands AFTER this PR and must not re-add it.
# Convert a baseline-table entry into a covDist margin.
#
# Papers report a covariate in whichever way suited the journal, so accept the
# forms that actually appear rather than one canonical pair.
.admPopSpec <- function(v, nm, dist) {
  bad <- function(...) stop("admixr2: covariate '", nm, "': ", ..., call. = FALSE)
  if (is.list(v) && (!is.null(v$quantile) || !is.null(v$values))) return(v)
  v <- unlist(v)
  g <- function(k) if (k %in% names(v)) unname(v[[k]]) else NULL
  # `c(median = 92, iqr = c(62, 118))` flattens to names median, iqr1, iqr2 --
  # c() renames the elements of a vector argument -- so a two-value entry has
  # to be collected by PREFIX, not by exact name.
  gg <- function(k) { i <- grepl(paste0("^", k, "[0-9]*$"), names(v))
                      if (any(i)) unname(v[i]) else NULL }
  # a proportion: `SEX = c(male = 0.55)` is how a baseline table prints it
  .known <- c("mean", "median", "sd", "cv", "meanlog", "sdlog")
  if (length(v) == 1L && !is.null(names(v)) && !names(v) %in% .known &&
      !grepl("^(iqr|range)[0-9]*$", names(v))) {
    p <- unname(v[[1L]])
    if (!is.finite(p) || p < 0 || p > 1)
      bad("'", names(v), " = ", p, "' must be a proportion between 0 and 1.")
    return(list(values = c(0, 1), probs = c(1 - p, p), .level = names(v)))
  }
  if (!is.null(g("meanlog")) && !is.null(g("sdlog")))
    return(list(meanlog = g("meanlog"), sdlog = g("sdlog")))
  ctr <- g("mean") %||% g("median")
  if (is.null(ctr)) bad("needs a `mean` or a `median` (or `meanlog`/`sdlog`).")
  # spread, in whichever currency the paper used
  sd <- g("sd")
  if (is.null(sd) && !is.null(g("cv"))) sd <- ctr * g("cv") / 100
  if (is.null(sd) && !is.null(gg("iqr"))) {
    q <- gg("iqr")
    if (length(q) != 2L) bad("`iqr` needs two values, e.g. iqr = c(62, 88).")
    # exact for the assumed shape rather than a rule of thumb: on the scale the
    # margin is modelled on, IQR = 1.349 sd
    q <- sort(q)
    # A REPORTED (median, IQR) NEED NOT BE CONSISTENT with the assumed shape,
    # and usually is not: real quartiles are rarely symmetric about the median
    # on any scale. Only two of the three numbers can be honoured, so honour
    # the median (the headline) and the IQR WIDTH (the spread), and say when
    # the third disagrees rather than quietly reproducing neither quartile.
    .chk <- function(mid, tol, what) {
      if (abs(mid - ctr) / ctr > tol)
        message("admixr2: covariate '", nm, "': median ", signif(ctr, 4),
                " and iqr c(", signif(q[1L], 4), ", ", signif(q[2L], 4),
                ") are not consistent with a ", what, " margin -- the ",
                "quartiles centre on ", signif(mid, 4), " instead. Matching ",
                "the median and the IQR WIDTH; the quartiles will come out ",
                "shifted. Give `sd =` if the paper reports one.")
    }
    if (identical(dist, "lnorm")) {
      .chk(sqrt(q[1L] * q[2L]), 0.02, "lognormal")
      return(list(meanlog = log(ctr), sdlog = diff(log(q)) / 1.34898))
    }
    .chk(mean(q), 0.02, "normal")
    sd <- diff(q) / 1.34898
  }
  if (is.null(sd) && !is.null(gg("range"))) {
    r <- gg("range")
    if (length(r) != 2L) bad("`range` needs two values, e.g. range = c(35, 160).")
    # A MIN-MAX IS NOT A SPREAD, and the conversion is a rule of thumb whose
    # error grows with n -- said out loud rather than folded in silently,
    # because the covariate spread is what a covariate coefficient is estimated
    # from and getting it wrong is not visible in any diagnostic.
    message("admixr2: covariate '", nm, "': deriving the spread from a min-max ",
            "via range/4, which is a rule of thumb, not an identity -- it is ",
            "reasonable near n = 30 and too WIDE for a large study. Prefer ",
            "`iqr = c(lo, hi)` or `sd =` when the paper gives one.")
    sd <- diff(sort(r)) / 4
  }
  if (is.null(sd))
    bad("needs a spread: `sd`, `cv` (as a percent), `iqr = c(lo, hi)` or ",
        "`range = c(min, max)`.")
  if (identical(dist, "lnorm")) {
    # match the reported mean and sd on the natural scale exactly
    s2 <- log1p((sd / ctr)^2)
    list(meanlog = log(ctr) - s2 / 2, sdlog = sqrt(s2))
  } else list(mu = ctr, sd = sd)
}

# `dist` goes LAST, and it is the covDist() argument -- NOT hard-coded "norm".
# The delegation below converts median/iqr, mean/cv and mean/range straight to
# mu/sd, which removes the `mean` field the caller's dist default keys on, so
# covDist(dist = "lnorm") was silently discarded for every vocabulary except
# c(mean=, sd=) and an allometric model got an unbounded-below margin.
.admCovSpecFromVec <- function(v, nm, dist = "normal") {
  bad <- function(...) stop("admixr2: covariate ", sQuote(nm), " ", ...,
                            call. = FALSE)
  if (is.list(v) || is.function(v)) return(v)
  if (!is.numeric(v)) bad("must be numeric.")
  k <- names(v)
  if (is.null(k) || any(!nzchar(k)))
    bad("needs NAMED values, e.g. c(mean = 72, sd = 16), c(mu = , sd = ) or ",
        "c(meanlog = , sdlog = ). An unnamed pair is ambiguous: c(0, 1) could ",
        "be a mean and an SD, or the two levels of a binary covariate.")
  k <- tolower(k)
  has <- function(...) all(c(...) %in% k)
  g <- function(x) as.numeric(v[[which(k == x)[1L]]])
  if (has("mean", "sd"))       return(list(mean = g("mean"), sd = g("sd")))
  if (has("mu", "sd"))         return(list(mu = g("mu"), sd = g("sd")))
  if (has("meanlog", "sdlog")) return(list(meanlog = g("meanlog"),
                                           sdlog = g("sdlog")))
  # CONTINUOUS-SUMMARY VOCABULARY IS NOT A SET OF LEVELS. `c(median = 92,
  # iqr = c(62, 118))` and `c(mean = 72, cv = 22)` are the forms admPopulation()
  # documents beside this function, and both fell through to the categorical
  # branch below: the first became a 3-level covariate with probabilities
  # 92/272, 62/272, 118/272, the second a 2-level one. The whole quadrature then
  # integrates over a covariate that does not exist, with no error. Delegate to
  # the same parser admPopulation() uses rather than teach two vocabularies.
  if (any(k %in% c("mean", "median", "sd", "cv", "meanlog", "sdlog")) ||
      any(grepl("^(iqr|range)[0-9]*$", k))) {
    sp <- .admPopSpec(stats::setNames(as.numeric(v), k), nm,
                      if (identical(dist, "lnorm")) "lnorm" else "norm")
    # [[ ]], NOT $. .admPopSpec returns meanlog/sdlog on the lnorm branch, and
    # BOTH `$mean` and `$sd` PARTIAL-MATCH those -- so this line silently
    # reinterpreted a lognormal's log-scale parameters as a normal's natural
    # ones: covDist(WT = c(mean = 72, cv = 22), dist = "lnorm") came back as
    # list(mu = 4.25, sd = 0.217), a normal margin centred at 4.25 kg, and every
    # quadrature node sat there instead of near 72. .admPopSpec never returns
    # `mean`, so with [[ ]] this is the no-op it was always meant to be.
    if (!is.null(sp[["mean"]])) sp <- list(mu = sp[["mean"]], sd = sp[["sd"]])
    return(sp)
  }
  # anything else named is read as a CATEGORICAL covariate: the names are the
  # levels' labels and the values their proportions, which is how a baseline
  # table reports sex, a genotype or a dosing band
  if (any(v < 0)) bad("has negative probabilities, reading it as categorical.")
  if (!sum(v) > 0) bad("has probabilities summing to zero.")
  # A SINGLE named proportion is a BINARY covariate, not a one-level one.
  # `SEX = c(male = 0.55)` is exactly how a baseline table prints it, and it is
  # what admPopulation() documents -- but normalising it like a level vector
  # divides 0.55 by itself, giving one level at probability 1, i.e. a CONSTANT.
  # The covariate then contributes no variation, its coefficient is not
  # identified, and covDraw() returns the same value for every subject. That is
  # a silent wrong answer on the documented spelling, and .admPopSpec (which
  # admPopulation routes through) already handled it correctly -- so the two
  # parsers disagreed on the one form the docs teach most.
  # A SINGLE value is only a proportion if it IS one. `c(kg = 70)` is not a
  # probability, and reading it as one gives P(absent) = -69: the multi-level
  # branch below is normalised by its sum, so it absorbed any scale silently,
  # but a lone value has no second number to be normalised against. Refuse it
  # and name both spellings -- the alternative is a negative probability
  # reaching the quadrature weights, where it is a wrong answer rather than an
  # error.
  if (length(v) == 1L) {
    p1 <- as.numeric(v)
    if (!is.finite(p1) || p1 < 0 || p1 > 1)
      bad("is a single value ", format(p1), ", which is read as the ",
          "proportion in level ", sQuote(names(v)), " and so must lie in ",
          "[0, 1]. For a continuous covariate give two numbers, e.g. ",
          "c(mean = ", format(p1), ", sd = ...); for a fixed value use the ",
          "study's `cov`, not `cov_dist`.")
    return(list(values = c(0, 1), probs = c(1 - p1, p1), labels = names(v)))
  }
  list(values = seq_along(v) - 1, probs = as.numeric(v) / sum(v),
       labels = names(v))
}

.admCovDistFromDf <- function(df) {
  low <- tolower(names(df))
  pick <- function(...) { i <- which(low %in% c(...))
    if (length(i)) i[1L] else NA_integer_ }
  inm <- pick("covariate", "cov", "name", "parameter")
  if (is.na(inm))
    stop("admixr2: a data.frame passed to covDist() needs a column naming the ",
         "covariate -- one of `covariate`, `cov`, `name` or `parameter`. Got: ",
         paste(names(df), collapse = ", "), ".", call. = FALSE)
  im <- pick("mean", "m"); isd <- pick("sd", "s", "stdev", "std")
  if (is.na(im) || is.na(isd))
    stop("admixr2: a data.frame passed to covDist() needs `mean` and `sd` ",
         "columns. Got: ", paste(names(df), collapse = ", "), ".", call. = FALSE)
  idst <- pick("dist", "distribution")
  stats::setNames(lapply(seq_len(nrow(df)), function(i) {
    sp <- list(mean = as.numeric(df[[im]][i]), sd = as.numeric(df[[isd]][i]))
    if (!is.na(idst)) sp$dist <- as.character(df[[idst]][i])
    sp
  }), as.character(df[[inm]]))
}

#' Describe the covariate distribution a study's subjects span
#'
#' Builds the `cov_dist` a study carries, from what a publication actually
#' reports: a mean and a standard deviation per covariate, and optionally a
#' correlation between them. Validation happens here, where an error can name
#' the covariate, rather than at the first objective evaluation of a fit.
#'
#' @param ... One argument per covariate, named as the model reads it. Each is
#'   one of:
#'   \describe{
#'     \item{`c(mean = , sd = )`}{Mean and SD on the covariate's own scale ---
#'       what a baseline table reports. Normal by default; `dist = "lnorm"`
#'       moment-matches a lognormal margin instead, which is what you want for
#'       a covariate that must stay positive (see `dist`).}
#'     \item{`c(mu = , sd = )` / `c(meanlog = , sdlog = )`}{A normal or
#'       lognormal margin given directly on its own scale.}
#'     \item{`c(label = prob, ...)`}{A categorical covariate: names are the
#'       level labels, values their proportions. Levels are coded `0, 1, ...`
#'       in the order given.}
#'     \item{a `list(...)`}{The canonical form, for anything else --- including
#'       `list(quantile = f)` for an arbitrary margin.}
#'   }
#'   Alternatively a single data.frame with `covariate`, `mean` and `sd`
#'   columns (and optionally `dist`): a baseline-characteristics table
#'   transcribed as-is.
#' @param cor Correlation between the covariates: a scalar for two of them, or
#'   a correlation matrix (named, in any order). Realised through a Gaussian
#'   copula on the declared margins.
#' @param joint Dependence a single correlation cannot express: a function
#'   receiving the matrix of uniforms admixr2 supplies and returning one named
#'   column per covariate, on each covariate's OWN scale. Overrides `cor`.
#'
#'   The sampler owns its margins. A copula returns dependence on the *copula
#'   scale* only, so a sampler that hands those uniforms back returns covariate
#'   values on `(0, 1)` --- every covariate mean 0.50 and SD 0.29, whatever
#'   margins were declared, with no error anywhere. Push them through the
#'   quantile function you intend.
#'
#'   `stratify` needs nothing further: an opaque sampler is banded by binning
#'   its output rather than by fixing its input uniforms.
#' @param dist Default margin for the `c(mean = , sd = )` form: `"normal"`
#'   (default) or `"lnorm"`. A per-covariate `dist` wins over it.
#'
#'   Choose `"lnorm"` when the model needs the covariate to stay positive. A
#'   normal margin is unbounded below and the quadrature reaches 3.75 standard
#'   deviations, so any covariate with a coefficient of variation above about
#'   0.27 (the guard is `mu - 3.75 * sd <= 0`) gets a node at or below zero,
#'   where an allometric or log term is `NaN`.
#'   `covDist()` warns when that would happen.
#'
#' @return A validated `cov_dist`, ready to pass to a study. Printing it shows
#'   each covariate's realised mean, SD and type.
#'
#' @examples
#' # what a baseline-characteristics table reports
#' cd <- covDist(WT = c(mean = 72, sd = 16), CRCL = c(mean = 90, sd = 25),
#'               cor = 0.6)
#' cd
#'
#' # a categorical covariate: labels are the levels, values their proportions
#' covDist(SEX = c(female = 0.55, male = 0.45))
#'
#' # or transcribe the table itself
#' covDist(data.frame(covariate = c("WT", "CRCL"),
#'                    mean = c(72, 90), sd = c(16, 25)))
#'
#' @seealso [covDraw()] to see what a specification describes,
#'   [covStrata()] to cut it into strata.
#' @export
covDist <- function(..., cor = NULL, joint = NULL,
                    dist = c("normal", "lnorm")) {
  dist <- match.arg(dist)
  a <- list(...)
  if (length(a) == 1L && is.null(names(a)) && is.data.frame(a[[1L]]))
    out <- .admCovDistFromDf(a[[1L]])
  else {
    if (!length(a) || is.null(names(a)) || any(!nzchar(names(a))))
      stop("admixr2: covDist() needs one NAMED argument per covariate, named ",
           "as the model reads it -- e.g. covDist(WT = c(mean = 72, sd = 16)).",
           call. = FALSE)
    out <- stats::setNames(
      lapply(seq_along(a), function(i)
        .admCovSpecFromVec(a[[i]], names(a)[i], dist)),
      names(a))
  }
  # the default margin applies only where the covariate did not choose one, and
  # only to the mean/sd shorthand -- every other form already names its scale
  for (nm in names(out)) {
    sp <- out[[nm]]
    if (is.list(sp) && !is.null(sp[["mean"]]) && is.null(sp[["dist"]]))
      out[[nm]][["dist"]] <- dist
  }
  if (!is.null(cor))   out[["cor"]]   <- cor
  if (!is.null(joint)) out[["joint"]] <- joint
  # Canonicalise NOW, so an impossible moment match or a non-PD correlation is
  # reported here, naming the covariate, rather than surfacing from inside the
  # first objective evaluation of a fit.
  out <- .admCovDistCanon(out)
  # A NORMAL margin is unbounded below, and the quadrature reaches |z| = 5.19
  # at the default 7 nodes -- so any covariate with a CV above about 0.27 gets
  # a node at or below zero, and an allometric or log term evaluated there is
  # NaN. That covers most PK covariates (weight 72+/-16 is CV 0.22), so it is
  # worth saying out loud rather than leaving to a solver failure. It is only a
  # WARNING: a centred covariate -- a log-ratio, a z-score, a change from
  # baseline -- is legitimately negative and needs no fixing.
  # .admCovSpecNames, NOT names(out): the metadata siblings are a function
  # (`joint`) and a matrix (`latentR`), and sp[["mu"]] on either is an error,
  # not a NULL. This is the third distinct bug from enumerating names(cov_dist)
  # directly -- it is what the accessor exists for.
  for (nm in .admCovSpecNames(out)) {
    sp <- out[[nm]]
    if (!is.list(sp) || is.null(sp[["mu"]]) || is.null(sp[["sd"]])) next
    if (sp[["mu"]] - 3.75 * sp[["sd"]] <= 0 && sp[["mu"]] > 0)
      warning("admixr2: covariate ", sQuote(nm), " has a NORMAL margin with ",
              "mean ", format(sp[["mu"]]), " and sd ", format(sp[["sd"]]),
              ", so the quadrature reaches ",
              format(round(sp[["mu"]] - 3.75 * sp[["sd"]], 2)),
              " -- at or below zero. If the model uses it in a power, ",
              "allometric or log term that is NaN, and you want ",
              'dist = "lnorm", which is positive by construction. Ignore this ',
              "if the covariate is genuinely centred.", call. = FALSE)
  }
  # Remember HOW dependence was declared. Canonicalisation turns `cor` into a
  # joint sampler, so without this the print method reports every correlated
  # spec as an opaque sampler and hides the number the user actually typed.
  structure(out, class = c("covDist", "list"),
            declared = if (!is.null(joint)) "joint sampler"
                       else if (is.null(cor)) "independent"
                       else if (length(cor) == 1L) sprintf("cor = %s", format(cor))
                       else "correlation matrix")
}

#' @export
print.covDist <- function(x, ...) {
  nms <- .admCovSpecNames(x)
  cat("<covDist>", length(nms), "covariate(s)\n")
  X <- try(covDraw(unclass(x), n = 4000L), silent = TRUE)
  if (inherits(X, "try-error")) {
    cat("  could not be drawn from:",
        conditionMessage(attr(X, "condition")), "\n")
    return(invisible(x))
  }
  d <- data.frame(
    covariate = nms,
    type = vapply(nms, function(n) {
      sp <- x[[n]]
      if (!is.null(sp$values)) "categorical"
      else if (!is.null(sp$meanlog)) "lognormal"
      else if (is.function(sp$quantile)) "quantile fn"
      else "normal" }, character(1)),
    mean = round(colMeans(X[, nms, drop = FALSE]), 3),
    sd = round(apply(X[, nms, drop = FALSE], 2L, stats::sd), 3),
    row.names = NULL)
  print(d, row.names = FALSE)
  for (n in nms) if (!is.null(x[[n]]$labels))
    cat(sprintf("  %s levels: %s\n", n,
        paste(sprintf("%s=%.0f (%.0f%%)", x[[n]]$labels, x[[n]]$values,
                      100 * x[[n]]$probs), collapse = ", ")))
  if (length(nms) > 1L) {
    dep <- attr(x, "declared") %||%
           (if (is.function(x[["joint"]])) "joint sampler" else "independent")
    cat("  dependence:", dep, "\n")
    if (length(nms) == 2L)
      cat(sprintf("    realised cor(%s, %s) = %+.3f\n", nms[1L], nms[2L],
                  stats::cor(X[, nms[1L]], X[, nms[2L]])))
  }
  invisible(x)
}




# Is the (coefficient, fixed effect, omega) ridge EXACTLY flat for this
# covariate? True only where Delta certifies as Gaussian -- see the note in
# .admWarnCovIdentifiability().
#
# Two independent things can break flatness, and they are decidable at different
# cost. The covariate's own LAW may make u non-normal whatever the model does
# (a discrete margin has no latent normal score at all); or Delta may be a
# non-affine function of a normal latent. Only the second needs Delta, so the
# first is settled first and for free.
#
# When Delta cannot be evaluated -- a hand-built pinfo with no structural
# values, an unusual model -- fall back to the MARGIN. A normal margin under an
# identity link is the textbook flat case and the one worth warning about;
# anything else is left unwarned, because a nearly-flat likelihood surfaces as
# an enormous standard error that covMethod = "r,s" now reports honestly,
# whereas a false "not identifiable" gets acted on and cannot be recovered from.
.admCovRidgeFlat <- function(ui, pinfo, studies, cv) {
  cd <- NULL
  for (s in studies) if (!is.null(s[["cov_dist"]][[cv]])) {
    cd <- s[["cov_dist"]]; break
  }
  # Declared only as a conditioning VALUE, and the caller has already
  # established every study uses the SAME one. Then the covariate enters as a
  # CONSTANT, so its coefficient is confounded with the fixed effect sharing its
  # argument -- exactly flat along (theta + d*a, gamma - d), whatever the
  # covariate's law would have been. A different flat direction from the
  # (gamma, omega) ridge below, and it needs no certificate.
  if (is.null(cd)) return(TRUE)
  spec <- cd[[cv]]
  # a discrete margin has no latent normal score and can never make Delta
  # normal -- decided without building anything
  if (!is.null(spec[["values"]])) return(FALSE)
  cert <- tryCatch({
    cn <- .admCovSpecNames(cd)
    sp <- .admShiftSpec(ui, cn, pinfo$eta_col_names)
    if (is.null(sp)) NA else {
      ar <- .admShiftRef(cd, cn)
      if (is.null(ar)) NA else {
        g <- .admCovGrid(cd, pinfo$cov_nodes %||% 7L)
        D <- .admShiftDelta(sp, .admShiftStruct(pinfo), g$X, ar)
        if (is.null(D)) NA else {
          D <- as.matrix(D)
          .admShiftGaussOK(D, g$W, g$z, ncol(D))
        }
      }
    }
  }, error = function(e) NA)
  if (!is.na(cert)) return(isTRUE(cert))
  # Delta unavailable: the margin is all there is to go on.
  !is.null(spec[["mu"]]) && !is.null(spec[["sd"]])
}
