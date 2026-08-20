# =============================================================================
# Covariate marginalisation for aggregate data modelling
# =============================================================================
#
# A study's subjects span a distribution of covariate values (`cov_dist`), and
# the aggregate (E, V) it reports is MARGINAL over that distribution. So the
# prediction has to be marginal too:
#
#   mu = E_{a,eta}[f]        V = Cov_{a,eta}(f) + residual
#
# Covariate-induced between-subject variability belongs INSIDE V_pred, because
# that is where the observed V carries it. Two paths compute those moments --
# the general per-row/product-grid path ("rows") and the shift path ("shift") --
# chosen per study by .admCheckCovariates(); see each section below. adgh's
# `cov_integration` selects between the product grid, its second-order Taylor
# surrogate and the shift.
#
# TWO EARLIER PATHS WERE REMOVED, and the reasons generalise:
#
#  * "collapse" folded a normal covariate's variance into Omega in closed form
#    (Omega + J Sigma_a J'). Exact, but it needed a bare `theta * COV` product,
#    a NORMAL covariate, the study solved at the covariate mean, AND grad =
#    "none" -- and every estimator defaults to a gradient, so it was unreachable
#    in a real fit. It was gated on rxode2's muRefCovariateDataFrame, which
#    recognises 1 of 8 realistic covariate parameterisations (the table that
#    used to sit here).
#
#    It is the GAUSSIAN SPECIAL CASE OF THE SHIFT, not a separate method: for a
#    linear Delta and a normal covariate, u = Delta(a) + eta is itself normal
#    with variance Omega + beta^2 Var(a), which is precisely the inflation.
#    Verified with the ODE solver taken out of the comparison (analytic 1-cmt f
#    over both node sets): 6.7e-16 on the mean and 5.0e-13 on the covariance at
#    31 covariate nodes, 2.5e-10 / 2.7e-8 at the default 7 -- an order of
#    magnitude under the solver's own ~1e-6 floor, which is why the same
#    comparison run THROUGH the solver bottoms out at 1.5e-06 / 1.7e-05 and
#    says nothing. "shift" covers the case with no distributional restriction
#    and with a gradient.
#  * "uq" replaced the affected eta column with quantiles of u = Delta(a) + eta.
#    That is the right idea, and it is what "shift" does -- but it INFERRED the
#    shift property from the model TEXT, which has four measured silent-wrong-
#    answer modes (discrete/multi-modal covariate 13-20% on the mean and 3-5x on
#    the covariance; a non-log link or additive effect 16-21% on the covariance;
#    two covariates on one eta leaving the result bit-for-bit the second
#    covariate alone; a user `quantile` integrated by 32 equal-weight midpoints
#    biasing the coefficient ~7% high). It had been disabled (`ok_uq <- FALSE`,
#    unconditional) rather than deleted, pending exactly the positive numerical
#    check that .admShiftVerify() now performs.
#
# THE RULE, of which everything below is a consequence: a block's PREDICTION
# must integrate over the same covariate distribution its DATA were aggregated
# over. Pooled data need a pooled (marginal) prediction; per-stratum data need
# the prediction marginalised within that stratum. Both mismatches cost, in
# opposite directions, and neither announces itself.
#
# NODE METHODS (gl / gh / taylor) WERE REMOVED, deliberately -- but NOT because
# they estimate wrongly. They score the study at fixed covariate values and
# combine the per-node -2LL linearly. Given per-node DATA that is a genuine
# likelihood and is EXACT: sum_k c_k n NLL_k is just sum_k n_k NLL_k with
# n_k = c_k n, i.e. ordinary multi-study fitting with quadrature weights standing
# in for stratum sizes, and it recovers a true coefficient of 0.75 as 0.7500.
#
# What it cannot take is the ONE pooled (E, V) a publication reports. There the
# same sum is sum_k w_k * NLL(same obs) = E_a[-2LL] -- an average of
# log-densities, so -2 log of an unnormalised geometric mean of densities, which
# is not a likelihood of anything. Its mean term decomposes exactly as
#
#   sum_k w_k r_k' Vi r_k  =  rbar' Vi rbar  +  tr(Vi %*% Cov_a(mu))
#
# and the second term is a covariate-spread penalty carrying no data at all: at
# the true parameters the first term is 0 and the second is the whole objective.
# It pulls the coefficient to 0.3045 (-59%), or 0.4567 (-39%) even when V_pred is
# supplied marginally, so it is the log-density pooling and not the variance
# handling that does it. That is a displaced TARGET, not estimator bias -- the
# same construction is exact the moment it is given the data it describes.
# Measured in validation/covariate-threeway.R (configuration A) and separated
# case by case in validation/covariate-node-retest.R.
#
# So the machinery was redundant where it was valid -- a publication reporting
# strata also reports their real n_k, which beats quadrature weights -- and
# inapplicable where it was not. Fit stratum summaries as ordinary studies, each
# with its own `n` and its own `cov_dist` -- NOT `cov`.
#
# `cov` is a POINT value, and plugging the stratum's mean covariate in is the
# ecological plug-in, which reintroduces the same bias at stratum scale: the
# aggregate relation equals the individual one only if the model is affine over
# that stratum's covariate support, or the stratum is degenerate. Measured
# against a true coefficient of 0.75, marginalising within the stratum recovers
# 0.7500 at every K, while plugging in the stratum mean gives 0.875 at K = 2
# (+17%), 0.782 at K = 4 (+4.3%) and 0.757 at K = 8 -- dying only as K grows.
# Note the sign: the plug-in biases UPWARD, the mirror of the node route's
# downward pull. A subgroup table reports the stratum's own mean and SD, which is
# exactly the truncated `cov_dist` this needs.
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

# Per-subject covariate draws from the COLLAPSED subspace.
#
# admc draws sobol(n, dim = n_eta + p) and keeps the covariate columns. QMC
# error grows with dimension, so when the covariates reach the model through r
# scalars the sequence should be n_eta + r wide, not n_eta + p. Drawing w in r
# dimensions and mapping z = U Lr w reproduces the parameter law exactly -- the
# same construction the quadrature design uses, sampled rather than gridded.
#
# The saving is in ACCURACY, not rows: admc pays the same number of solves
# either way.
#
# AND IT IS SMALL. An earlier note here claimed ~2x on the covariance. That was
# measured across randtoolbox `seed` values, and `seed` is INERT -- sobol(n, d,
# seed = 1) and seed = 999 are bit-identical -- so it was a median over
# identical numbers: one realisation with no error bar. (`scrambling` is
# disabled outright in that package too; check before trusting either.)
#
# Re-measured with Cranley-Patterson shifts, 24 replicates, paired within
# replicate, the covariance error ratio is 1.12x [0.89, 1.84] at n = 1000,
# 1.24x [0.90, 1.53] at 4000 and 1.47x [0.78, 2.17] at 16000 -- interquartile
# range spanning 1.0 in every cell. So: a change of variables that is exactly
# right and costs nothing, with no benefit this measurement can distinguish
# from zero. Kept because it is exact and free, NOT because it is faster.
#
# A DISCRETE covariate keeps a dimension of its own: its level is genuinely
# random per subject and there is nothing to project. Still narrower than the
# full product whenever r < pc.
.admCovRowsCollapsed <- function(co, n, n_eta) {
  pc <- length(co$cn); nd <- length(co$dn); r <- co$r
  u  <- randtoolbox::sobol(n, dim = n_eta + r + nd)
  if (!is.matrix(u)) u <- matrix(u, nrow = n)
  u  <- u[, n_eta + seq_len(r + nd), drop = FALSE]
  u  <- pmin(pmax(u, .Machine$double.eps), 1 - .Machine$double.eps)
  w  <- matrix(stats::qnorm(u[, seq_len(r), drop = FALSE]), nrow = n)
  Z  <- w %*% co$Lr %*% t(co$U)                       # z = U Lr w
  # the same clamp .admCovCollapse applies, for the same reason: the uniforms
  # are clamped above but qnorm() of the clamp is already +-8.1, and the
  # rotation carries that further still, so pnorm() saturates at exactly 1 and
  # an unbounded margin quantile comes back infinite
  X  <- vapply(seq_len(pc), function(k)
    .admCovQuantile(co$cd[[co$cn[k]]],
                    pmin(pmax(stats::pnorm(Z[, k]), .Machine$double.eps),
                         1 - .Machine$double.eps)), numeric(n))
  if (!is.matrix(X)) X <- matrix(X, n, pc)
  colnames(X) <- co$cn
  if (nd) {
    D <- vapply(seq_len(nd), function(j) {
      sp <- co$cd[[co$dn[j]]]
      lv <- as.numeric(sp[["values"]])
      pr <- sp[["probs"]] %||% rep(1 / length(lv), length(lv))
      lv[findInterval(u[, r + j], cumsum(pr / sum(pr))) + 1L]
    }, numeric(n))
    if (!is.matrix(D)) D <- matrix(D, n, nd)
    colnames(D) <- co$dn
    X <- cbind(X, D)
  }
  if (!all(is.finite(X))) return(NULL)
  X[, co$nms, drop = FALSE]
}

# Attach per-row covariate values to a study for the GENERAL path. Returns the
# study unchanged on the shift path (where the covariate is held at its
# reference and its whole effect rides in the shifted eta column) and when no
# distribution is declared.
.admStudyCovRows <- function(s, pinfo, n_row, struct = NULL) {
  if (!identical(s$.adm_cov_path, "rows")) return(s)
  co <- s[[".adm_cov_collapse"]]
  # The collapsed subspace is aimed by the CURRENT structural thetas, so a
  # caller that cannot supply them gets the full product draw instead. That is
  # not a fallback for convenience: without them the rotation is the one the
  # STARTING values implied, and using it would be silently wrong rather than
  # merely slower. The batched Hessian paths are the case -- they hold many
  # parameter vectors at once, so no single rotation serves them.
  if (!is.null(co) && !is.null(struct))
    co <- .admCovRefresh(co, .admShiftStruct(pinfo, struct))
  else if (!is.null(co)) co <- NULL
  # The collapsed subspace when one was certified, the full product otherwise.
  # r == pc is admitted by the collapse for the sake of adgh's node search, but
  # it is no reduction in DIMENSION, and dimension is the whole of admc's
  # interest here -- taking it would rotate every draw for no gain. Measured:
  # a sampler that drops a direction outright is worse at every n and FLOORS
  # (8.3e-3 on V, against 3.0e-4 still falling), so admc keeps every direction
  # a nonzero singular value earns.
  cr <- if (!is.null(co) && isTRUE(co$r < co$pc))
    tryCatch(.admCovRowsCollapsed(co, n_row, pinfo$n_eta),
             error = function(e) NULL) else NULL
  s$cov_rows <- cr %||% .admCovRowsFor(s$cov_dist, n_row, pinfo$n_eta)
  s
}

# Refuse `cov_dist` for an estimator that has no covariate path.
#
# Silence here is the dangerous outcome, not an error: every study carries a
# covariate VALUE as well (derived from the distribution when not given), so an
# unwired estimator does not fail -- it solves at the covariate mean and reports
# a fit whose omega has quietly absorbed the between-subject covariate spread.
# Measured on the general path's own test model, that is omega 0.30 -> 0.44.
.admRefuseCovariates <- function(studies, est) {
  has <- vapply(studies, function(s) !is.null(s[["cov_dist"]]), logical(1))
  if (!any(has)) return(invisible(NULL))
  stop("admixr2: `", est, "` does not support covariate marginalisation. ",
       "Stud", if (sum(has) > 1L) "ies " else "y ",
       paste(sQuote(names(studies)[has]), collapse = ", "),
       " declare(s) `cov_dist`. Use `admc` or `adgh`, both of which integrate ",
       "the covariate distribution -- dependent ones included; running here ",
       "would silently solve at the covariate mean and inflate omega instead.",
       call. = FALSE)
}

# `grad` and `est` are both VESTIGIAL as of the collapse/uq removal -- neither
# enters the routing any more (see the note on the default path below) -- and
# they are kept only because all four drivers call this positionally. Drop them
# when those call sites are next touched.
.admCheckCovariates <- function(.ui, pinfo, studies) {
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

  for (nm in names(studies)[has]) {
    cd <- studies[[nm]]$cov_dist
    # `cov` is documented as name -> value, and a named NUMERIC VECTOR is the
    # natural way to write that. It must be a LIST here: for an atomic vector
    # `x[["absent"]]` is an error ("subscript out of bounds") where a list
    # returns NULL, so the fill-from-distribution test below died on any study
    # that pinned one covariate with `cov` and described another with
    # `cov_dist`. Coerced once, at the only place that writes into it.
    if (!is.null(studies[[nm]][["cov"]]) && !is.list(studies[[nm]][["cov"]]))
      studies[[nm]][["cov"]] <- as.list(studies[[nm]][["cov"]])
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

    # SHIFT. If the covariates act only as a shift of a mu-referenced argument
    # -- f(a, eta) == f(a_ref, eta + Delta(a)) -- the whole covariate dimension
    # leaves the solve: the covariates are held at their reference and the
    # affected eta column carries quantiles of u = Delta(a) + eta instead. Cost
    # becomes CONSTANT in the number of covariates, against the product grid's
    # cov_nodes^p.
    #
    # ADMITTED BY NUMERICAL VERIFICATION, never by reading the model text. The
    # retired `uq` route inferred exactly this property from syntax and had four
    # measured silent-wrong-answer modes as a result (see the file header);
    # .admShiftVerify() evaluates the identity against the COMPILED model and
    # separates valid forms (1e-12..1e-14) from invalid ones (3e-02..7e-01) by
    # ten orders of magnitude.
    #
    #   "shift"  demand it; a refusal is an ERROR naming the reason, because the
    #            user asked for this path specifically.
    #   "auto"   try it; a refusal falls back to the product grid with a
    #            message. Never an error -- "auto" is a speed lever, and the
    #            fallback is the more accurate path, so a refusal costs solve
    #            rows and nothing else.
    #
    # `"auto"` is NOT the default. Turning it on changes every existing
    # covariate fit's numbers at the ~1e-5 level (that is the tolerance the
    # shift-vs-grid agreement test pins), and a default that silently moves
    # results is the class of change this package keeps paying for. Opt in.
    .ci <- pinfo$cov_integration %||% "quadrature"
    if (identical(studies[[nm]]$.adm_cov_path, "rows") &&
        .ci %in% c("shift", "auto")) {
      .cn  <- .admCovSpecNames(cd)
      .sp  <- .admShiftSpec(.ui, .cn, pinfo$eta_col_names)
      .ar  <- if (is.null(.sp)) NULL else .admShiftRef(cd, .cn)
      # A DISCRETE covariate makes u = Delta(a) + eta a multi-modal mixture.
      # The shift IDENTITY still holds -- Delta is well defined at each level --
      # but Gauss-Hermite nodes on u cannot resolve separated modes: measured
      # 13% on V and 21.4 -2LL units against the product grid at identical cost,
      # and it does NOT converge away with more nodes (1891 rows still 1.2% out).
      #
      # Scoped to the covariates that actually REACH the shifted argument. One
      # that is discrete but enters the model elsewhere -- a SEX term on a
      # parameter the covariate never shifts -- makes no mode in u, and
      # .admShiftVerify is what decides whether holding it at its reference is
      # legitimate. Testing every declared covariate instead refused the shift
      # for models where nothing about it was multi-modal.
      .used <- if (is.null(.sp)) .cn else
        intersect(.cn, unlist(lapply(.sp$rhs, all.vars)))
      .disc <- vapply(.used, function(n) !is.null(cd[[n]][["values"]]),
                      logical(1))
      .why <- NULL
      if (is.null(.sp)) .why <- paste0(
        "the covariates do not all enter one model assignment together with ",
        "exactly one random effect, so no single shifted column can carry them")
      # ... and reported through `.why` rather than as an error, so that `auto`
      # falls back to the product grid the way it does for every other
      # disqualification. It used to stop() unconditionally, which made a
      # discrete covariate the one property that could not be auto-detected
      # around.
      # A DISCRETE covariate reaching the shifted argument is no longer a
      # disqualification, it is a STRATIFICATION instruction: .admCovGrid
      # enumerates its levels exactly, so conditioning on them leaves each cell
      # with only the continuous covariates varying -- the mild sub-mixture the
      # quadrature resolves well -- and the cells recombine with their own
      # probabilities. See .admShiftNodesStrat(). What that buys is keeping the
      # CONTINUOUS covariates off the product grid, which is where the n_cov^p
      # cost lives; the discrete ones cost K cells either way.
      #
      # It is still refused where the cells cannot be formed: a JOINT sampler
      # mixes the uniforms before mapping them to levels, so a grid row's
      # discrete value is not a cell label and the weights are not the level
      # probabilities (measured: a 0.55/0.45 covariate came off the grid at
      # 0.477, and the error does not shrink with cov_nodes).
      else if (length(.used) && any(.disc) && is.function(cd[["joint"]]))
        .why <- paste0(
          "covariate(s) ", paste(sQuote(.used[.disc]), collapse = ", "),
          " are discrete, reach the shifted argument, and are drawn through a ",
          "JOINT sampler, so the grid carries no exact level cells to ",
          "condition on")
      # ALL discrete: the shift buys nothing. Its saving is eliminating the
      # CONTINUOUS covariate dimension; the discrete levels cost K cells on
      # either route, and n_u is inflated above n_nodes to resolve the widened
      # u, so the stratified shift comes out with MORE rows than the product
      # grid (measured 26 vs 18 for one binary covariate, 48 vs 36 for a
      # four-level one) -- and as an approximation where the grid enumerates
      # exactly. With a continuous covariate alongside it is the other way
      # round by a mile: 28 vs 162.
      else if (length(.used) && all(.disc)) .why <- paste0(
        "every covariate reaching the shifted argument is discrete, and the ",
        "product grid enumerates those exactly at no greater cost")
      else if (is.null(.ar)) .why <- "a covariate has no finite mean to shift against"
      # The grid and Delta cost no solve and are needed by the off-diagonal
      # test below as well as by the sizing further down, so build them once.
      .g <- .D0 <- NULL
      if (is.null(.why)) {
        .g  <- .admCovGrid(cd, pinfo$cov_nodes %||% 7L)
        .D0 <- tryCatch(.admShiftDelta(.sp, .admShiftStruct(pinfo), .g$X, .ar),
                        error = function(e) NULL)
      }
      # Delta = c + B z is absorbed into Omega as Omega + P, P_SS = B B', which
      # keeps every correlation a substituted column would have dropped.
      .abs <- !is.null(.D0) && .admShiftGaussOK(as.matrix(.D0), .g$W, .g$z,
                                                ncol(as.matrix(.D0)))
      # HOW A CORRELATED OMEGA IS CARRIED. The plain shift replaces ONE eta
      # column and rebuilds the others as sqrt(Omega_kk) * z, which silently
      # discards every Cholesky off-diagonal -- including between two etas the
      # covariate never touches. Measured on a 3-eta model with
      # cov(eta.v, eta.ka) estimated: 2.5e-2 on E, 78x on V and 951 -2LL units;
      # the same model with that covariance zeroed is 6.5e-10 / 3.2e-07.
      #
      # There are two ways not to drop it, and which one applies is decided by
      # the certificate, not by the correlation:
      #
      #   Delta certifies  -> ABSORPTION. No column is substituted at all: the
      #                       whole eta vector is drawn from Omega + P, so
      #                       Cov(u_S, eta_O) stays Omega_SO exactly.
      #   Delta does not   -> CONDITIONING. eta_O comes off its own grid on
      #                       chol(Omega_OO) and u_S from the conditional law
      #                       given it. See .admCondShiftParts().
      #
      # A DIAGONAL Omega keeps the substitution, bit-identically -- conditioning
      # reduces to it there (Omega_SO = 0 gives K = 0), but running the general
      # code would move existing fits by rounding alone.
      .cond <- FALSE
      if (is.null(.why))
        .cond <- !.abs && length(which(!pinfo$chol_diag)) > 0L
      if (is.null(.why)) {
        .w <- .admShiftVerify(.sp, .ui, NULL, pinfo, studies[[nm]], NULL,
                              .cn, .ar)
        if (is.na(.w) || .w > 1e-8) .why <- paste0(
          "the shift identity f(a, eta) == f(a_ref, eta + Delta(a)) does not ",
          "hold for this model (worst relative discrepancy ",
          if (is.na(.w)) "could not be evaluated" else format(.w, digits = 3),
          "); an additive covariate effect, a link the covariate sits outside ",
          "of, or a covariate on a parameter with no random effect all do this")
      }
      if (!is.null(.why) && identical(.ci, "shift"))
        bad("study '", nm, "': cov_integration = \"shift\" is not applicable -- ",
            .why, ". Use cov_integration = \"quadrature\".")
      if (!is.null(.why)) {
        # Recorded on the study as well as messaged: a message is easy to miss
        # in a fit's output, and "why did auto not take the fast path" is the
        # question this answers.
        studies[[nm]]$.adm_cov_shift_why <- .why
        message("admixr2: study '", nm, "': cov_integration = \"auto\" is ",
                "using the covariate product grid -- ", .why, ".")
      } else {
        # n_u is fixed HERE, from DATA, so the objective cannot step when
        # the optimizer moves omega. It used to scale with the current
        # sqrt((Var(Delta) + omega^2)/omega^2), so the node count changed
        # mid-fit and the objective jumped 0.078 -2LL units across each
        # switch -- enough to send an FD Hessian entry from -3169 to
        # -256961. Sized once at the initial omega, then held.
        .nn0 <- as.integer(pinfo$n_nodes %||% 7L)
        .D0  <- tryCatch(.admShiftDelta(.sp, .admShiftStruct(pinfo), .g$X,
                                        .ar), error = function(e) NULL)
        .od  <- tryCatch(diag(as.matrix(pinfo$omega_init))[[
                  match(.sp$eta[1L], pinfo$eta_col_names)]],
                  error = function(e) NA_real_)
        .om0 <- sqrt(if (is.finite(.od) && .od > 0) .od else 0.09)
        .vD0 <- if (is.null(.D0)) 0 else
          max(colSums(.g$W * as.matrix(.D0)^2) -
              colSums(.g$W * as.matrix(.D0))^2, 0)
        .nu  <- min(101L, max(.nn0, as.integer(ceiling(
                  .nn0 * sqrt((.vD0 + .om0^2) / .om0^2)))))
        studies[[nm]]$.adm_cov_shift <- list(
          spec = .sp, aref = .ar, cov_names = .cn, cond = .cond,
          # which grid rows share a discrete cell; NULL when none is discrete,
          # and then every node path below is exactly what it was
          strata = .admShiftStrata(cd, .g$X),
          eta_idx = match(.sp$eta, pinfo$eta_col_names),
          m = length(.sp$eta), X = .g$X, W = .g$W, z = .g$z, n_u = .nu,
          # Fixed at admission for the same reason n_u is: a path that could
          # flip mid-fit would step the objective.
          #
          # Taken wherever the column substitution would cost a FINITE-DIFFERENCE
          # gradient: a correlated Omega (which it cannot represent at all) and a
          # VECTOR shift (whose later coordinates move through the Rosenblatt
          # posterior weights, a chain .admShiftDu does not carry). The
          # absorption has no recursion and its Cholesky differential is exact,
          # so both become analytic. A scalar shift on a diagonal Omega keeps its
          # own closed form, which is already analytic and already verified.
          absorb = .abs && (length(which(!pinfo$chol_diag)) > 0L ||
                            length(.sp$eta) > 1L))
        studies[[nm]]$.adm_cov_path <- "shift"
      }
    }

    # cov_integration = "taylor" replaces the product grid on the "rows" path
    # (and only there -- a study that took the shift has no product grid left to
    # expand). Its refusals are properties of `cov_dist` alone, so build the
    # design once HERE, where the message can name the study, rather than
    # letting the first objective evaluation error out of the middle of a fit.
    if (identical(studies[[nm]]$.adm_cov_path, "rows") &&
        identical(.ci, "taylor")) {
      .td <- tryCatch(.admCovTaylorDesign(cd, pinfo$cov_taylor_h %||% 1),
                      error = function(e) conditionMessage(e))
      if (is.character(.td))
        bad("study '", nm, "': ", sub("^admixr2: ", "", .td))
      # KEEP it. The design is a pure function of `cov_dist` and `cov_taylor_h`,
      # both of which are DATA, but .adghGrid runs inside the objective -- so
      # rebuilding it there re-ran .admCovDistMoments' 8192-row Sobol pass
      # through the joint sampler on every evaluation: 8.6 ms a call, ~21 s of
      # pure overhead across a fit, for a number that cannot change. The
      # independent path was 0.6 ms and never noticed.
      #
      # Numeric only (no closures), so it serialises to a parallel-restart
      # worker by value like the rest of the study.
      studies[[nm]]$.adm_cov_taylor <- .td
    }
  }

  # DIMENSION COLLAPSE. Where the covariates reach the model through a single
  # scalar -- p covariates on one parameter, the allometric case -- the integral
  # is ONE-dimensional however many covariates there are, and the product grid
  # was integrating it in p. This is the shift's argument with the random effect
  # removed, so it applies exactly where the shift refuses for want of an eta.
  #
  # Cached for the same reason the Taylor design is: it is a pure function of
  # `cov_dist` and the model, both fixed for the fit, and .adghGrid runs inside
  # the objective. Numeric only, so it serialises to a restart worker by value.
  #
  # Only where the shift did NOT take the study: the shift removes the
  # dimension rather than reducing it, so it is cheaper still.
  for (nm in names(studies)) {
    s_nm <- studies[[nm]]
    if (is.null(s_nm[["cov_dist"]])) next
    if (identical(s_nm$.adm_cov_path, "shift")) next
    if (identical(pinfo$cov_integration %||% "quadrature", "taylor")) next
    .co <- tryCatch(.admCovCollapse(.ui, pinfo, s_nm[["cov_dist"]],
                                    pinfo$cov_nodes %||% 7L,
                                    cov_fixed = s_nm[["cov"]]),
                    error = function(e) NULL)
    if (!is.null(.co)) {
      studies[[nm]]$.adm_cov_collapse <- .co
      message("admixr2: study '", nm, "': the ", .co$p, " covariates reach the ",
              "model through ", .co$r,
              if (.co$r == 1L) " scalar" else " independent scalars",
              ", so the integral is ", .co$r, "-dimensional -- using ",
              nrow(.co$X), " design points rather than a ", .co$p,
              "-way product grid.")
    }
    # JOINT: the etas are latent normal directions too, and the design crosses
    # them with the covariate design as if the two were independent. Where an
    # eta and a covariate index reach the model through the same sum they are
    # ONE direction, and the rank is bounded by how many PARAMETERS the latents
    # reach -- not by how many etas and covariates there are.
    #
    # Tried after the covariate collapse and preferred over it where it holds:
    # it subsumes that reduction and adds the eta block. adgh only -- this is a
    # quadrature construction, and admc reaches its etas through .admMakeZ.
    if (pinfo$n_eta > 0L && !isTRUE(s_nm$is_joint)) {
      .jc <- tryCatch({
        .j0 <- .admJointCollapse(.ui, pinfo, s_nm[["cov_dist"]],
                                 pinfo$cov_nodes %||% 7L, s_nm, NULL,
                                 cov_fixed = s_nm[["cov"]])
        .p0 <- .admBuildOptVec(pinfo)$p0
        .pr <- .admUnpack(.p0, pinfo)
        .admJointAdmit(.j0, .admShiftStruct(pinfo, .pr$struct), .pr$L)
      }, error = function(e) NULL)
      # AND IT MUST BE CHEAPER THAN WHAT IT REPLACES. Subsuming the covariate
      # collapse on RANK does not make it cheaper: where the latents share no
      # directions the joint rank is the whole latent dimension, and the
      # per-direction cap then applies to every one of them. Measured on 3 etas
      # with 2 covariates on a parameter carrying none -- rank 4 of 5 -- the
      # joint design is 6561 rows against 1750, 3.75x WORSE, and the absolute
      # max_rows cap is far too loose to catch it.
      #
      # The alternative is the eta grid crossed with whatever covariate design
      # would otherwise be used. No discrete covariates can be present here --
      # .admJointCollapse refuses those -- so the fallback grid is cov_nodes^pc.
      if (!is.null(.jc)) {
        .alt <- (pinfo$n_nodes %||% 5L)^pinfo$n_eta *
                (if (!is.null(.co)) nrow(.co$X)
                 else (pinfo$cov_nodes %||% 7L)^.jc$pc *
                      max(.jc$n_cell %||% 1L, 1L))
        if (.jc$m^.jc$r * max(.jc$n_cell %||% 1L, 1L) >= .alt) .jc <- NULL
      }
      if (!is.null(.jc)) {
        studies[[nm]]$.adm_cov_joint <- .jc
        message("admixr2: study '", nm, "': the ", pinfo$n_eta,
                " random effects and ", .jc$pc, " covariates span ", .jc$r,
                if (.jc$r == 1L) " direction" else " directions",
                " jointly -- using ", .jc$m^.jc$r,
                " design points rather than crossing an ", pinfo$n_eta,
                "-way eta grid with the covariate design.")
      }
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
  if (is.null(sd) || !is.finite(sd) || sd < 0)
    stop("admixr2: covariate ", sQuote(nm), " gives `mean` without a finite ",
         "non-negative `sd`.", call. = FALSE)
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
  margins <- lapply(nms, function(nm) cov_dist[[nm]])
  cov_dist[["joint"]] <- local({
    nms <- nms; margins <- margins; Lc <- Lc; d <- d
    function(u) {
      cl  <- function(x) pmin(pmax(x, 1e-12), 1 - 1e-12)
      z   <- stats::qnorm(cl(u)) %*% Lc
      out <- vapply(seq_len(d), function(j)
        .admCovQuantile(margins[[j]], cl(stats::pnorm(z[, j]))),
        numeric(nrow(u)))
      out <- matrix(out, nrow = nrow(u), ncol = d)
      colnames(out) <- nms
      out
    }
  })
  cov_dist
}

.admCovQuantile <- function(spec, u) {
  if (is.function(spec$quantile)) return(as.numeric(spec$quantile(u)))
  if (!is.null(spec$values)) {
    pr <- spec$probs %||% rep(1 / length(spec$values), length(spec$values))
    pr <- pr / sum(pr)
    return(as.numeric(spec$values)[findInterval(u, cumsum(pr), rightmost.closed = TRUE) + 1L])
  }
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

  out <- vapply(seq_len(d), function(k) .admCovQuantile(cov_dist[[nms[k]]], u[, k]),
                numeric(n))
  if (!is.matrix(out)) out <- matrix(out, nrow = n)
  colnames(out) <- nms
  out
}

# Deterministic quadrature nodes + weights for one covariate (adgh path).
# Discrete specs enumerate exactly; normal/lognormal use Gauss-Hermite.
.admCovNodesFor <- function(spec, n_nodes) {
  if (!is.null(spec$values)) {
    pr <- spec$probs %||% rep(1 / length(spec$values), length(spec$values))
    return(list(x = as.numeric(spec$values), w = pr / sum(pr)))
  }
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
                     pmin(pmax(stats::pnorm(g$x), .Machine$double.eps),
                          1 - .Machine$double.eps))), w = g$w))
  # the closed forms are safe: they are built from g$x directly, never from a
  # probability, so no saturation can occur
  if (!is.null(spec$meanlog)) list(x = exp(spec$meanlog + spec$sdlog * g$x), w = g$w)
  else                        list(x = spec$mu + spec$sd * g$x,              w = g$w)
}

# =============================================================================
# Second-order Taylor expansion of the MARGINAL MOMENTS
# =============================================================================
#
# The product grid above costs n_nodes^p solves-worth of rows. The alternative
# is to stop integrating the covariate and expand the two things the aggregate
# likelihood actually consumes -- the marginal mean and covariance -- about the
# covariate mean. With g(a) = E_eta[f(a, eta)] and Vc(a) = Cov_eta(f | a),
#
#   E_marg = E_a[g(a)]                ~= g(m) + (1/2) sum_j v_j g''_j(m)
#   V_marg = E_a[Vc(a)] + Cov_a(g(a)) ~= Vc(m) + (1/2) sum_j v_j Vc''_j(m)
#                                          + sum_j v_j g'_j(m) g'_j(m)'
#
# for v_j = Var(a_j), differences taken per covariate at h_j = hfrac * sd_j.
# That is 1 + 2p covariate points instead of n_nodes^p, and the likelihood is
# evaluated ONCE at the approximate marginal moments rather than per node.
#
# THE THIRD TERM IS THE POINT. `sum_j v_j g' g''` is Cov_a(g(a)) -- the
# between-subject covariance the covariate itself induces, and to second order
# it IS the omega'^2 = omega^2 + beta^2 Var(a) inflation a linear covariate
# effect produces in closed form, appearing here as a rank-p addition to V.
# Dropping it leaves
# E_a[Vc(a)], i.e. the average WITHIN-covariate covariance, which is the same
# class of mistake as writing V_struct + Sigma(mu) for the residual: a number
# that is finite, plausible and biased low, with the covariate's whole
# contribution to V missing. `test-covariate.R` fails without it.
#
# Everything here is a moment expansion, so it needs only that `a` has a finite
# mean and variance and that g is smooth over the covariate -- not that `a` is
# normal. What it does NOT have is cross terms: d^2 g/(da_j da_k) and
# Cov(a_j, a_k) both enter at the same order, so a dependent covariate
# distribution is REFUSED rather than integrated as if it were independent.
# Refusals are listed in .admCovTaylorDesign() and every one of them is an
# error, never a silent approximation.
#
# The design points, the combination coefficients and the row map, for one
# study's covariate distribution. Depends only on `cov_dist` -- which is DATA --
# so it is a pure function of the study and can be built up front (that is what
# .admCheckCovariates() does with it) as well as inside the objective.
.admCovTaylorDesign <- function(cov_dist, hfrac = 1) {
  cov_dist <- .admCovDistCanon(cov_dist)
  bad <- function(...) stop("admixr2: ", ..., call. = FALSE)
  if (!is.numeric(hfrac) || length(hfrac) != 1L || !is.finite(hfrac) || hfrac <= 0)
    bad("`cov_taylor_h` must be one positive finite number; got ",
        paste(format(hfrac), collapse = ", "), ".")
  nms <- .admCovSpecNames(cov_dist)
  d   <- length(nms)
  if (d == 0L) return(NULL)

  # A DISCRETE covariate is not expanded, it is ENUMERATED. A step function has
  # no curvature -- design points would land ON levels and the differences would
  # measure the gap between two levels rather than a derivative -- but it does
  # not need one: its levels and their probabilities ARE the integration rule,
  # exactly. So the design is a product of an exact discrete enumeration with
  # the continuous cubature, and the two are combined by the law of total
  # variance across the cells:
  #
  #   E   = sum_l p_l E_l
  #   V   = sum_l p_l [ within-cell expansion ]          (within)
  #       + sum_l p_l (E_l - E)(E_l - E)'                (between, EXACT)
  #
  # The between-cell term is a quadratic form in a LINEAR functional of the
  # solved rows, exactly like the rank-p derivative term, so it rides the same
  # (Dw, var) machinery and neither the moment assembly nor the analytical
  # gradient needs to know it exists. Cost is L * (1 + 2 p_continuous) points.
  disc <- vapply(nms, function(n) !is.null(cov_dist[[n]][["values"]]), logical(1))
  cn   <- nms[!disc]                     # continuous
  dn   <- nms[disc]                      # discrete
  dc   <- length(cn)

  Rz <- cov_dist[["latentR"]]
  if (length(dn) && !is.null(Rz) && any(abs(Rz[lower.tri(Rz)]) > 0)) {
    Rd <- Rz[disc, !disc, drop = FALSE]
    if (length(Rd) && any(abs(Rd) > 0))
      bad("covariate(s) ", paste(sQuote(dn), collapse = ", "), " are discrete ",
          "and DEPENDENT on continuous covariate(s) ",
          paste(sQuote(cn[apply(abs(Rd) > 0, 2L, any)]), collapse = ", "),
          '. cov_integration = "taylor" enumerates a discrete covariate at its ',
          "levels, but under dependence a level is a TRUNCATION of the latent ",
          "normal rather than a point, so the continuous covariates' ",
          "conditional mean and spread differ from cell to cell and expanding ",
          "them about the marginal mean would be the wrong expansion in every ",
          'cell. Use cov_integration = "quadrature", which resolves this with ',
          "a pooled sampler.")
    if (length(dn) > 1L) {
      Rdd <- Rz[disc, disc, drop = FALSE]
      if (any(abs(Rdd[lower.tri(Rdd)]) > 0))
        bad("discrete covariates ", paste(sQuote(dn), collapse = ", "),
            " are declared DEPENDENT on each other, so their joint cell ",
            "probabilities are not the product of their level probabilities ",
            'and the enumeration would be wrong. Use cov_integration = ',
            '"quadrature".')
    }
  }

  # ---- discrete cells: exact product enumeration -----------------------------
  if (length(dn)) {
    lev <- lapply(dn, function(n) as.numeric(cov_dist[[n]][["values"]]))
    prb <- lapply(dn, function(n) { s <- cov_dist[[n]]
      p <- s[["probs"]] %||% rep(1 / length(s[["values"]]), length(s[["values"]]))
      p / sum(p) })
    cells <- as.matrix(expand.grid(lev, KEEP.OUT.ATTRS = FALSE))
    colnames(cells) <- dn
    pcell <- apply(as.matrix(expand.grid(prb, KEEP.OUT.ATTRS = FALSE)), 1L, prod)
  } else {
    cells <- matrix(numeric(0), 1L, 0L)
    pcell <- 1
  }
  L <- length(pcell)

  # ---- continuous cubature ---------------------------------------------------
  # DEPENDENCE IS FREE, in the eigenbasis of the latent correlation. Both
  # covariance-weighted terms of the expansion are DIRECTIONAL: with
  # Sigma = sum_k lam_k v_k v_k',
  #
  #   tr(Sigma Hess g) = sum_k lam_k * d2g/dv_k^2
  #   J Sigma J'       = sum_k lam_k * (J v_k)(J v_k)'
  #
  # so differencing along the EIGENVECTORS gives both from 1 + 2p points -- the
  # same count as the independent case, at ANY correlation. This is the
  # third-degree spherical-radial cubature rule, i.e. the unscented transform's
  # sigma points, and the step h_k propto sqrt(lam_k) is its scaling.
  #
  # THE EXPANSION RUNS IN THE LATENT NORMAL, NOT ON THE COVARIATE SCALE. Every
  # margin is x = F^-1(Phi(z)), so differencing in z and mapping each design
  # point back through F^-1 keeps it inside the covariate's own support by
  # construction; differencing on the covariate scale does not (a lognormal
  # weight with mean 2 and sd 8 put a design point at -2, and the model returned
  # NaN from (WT/70)^theta). For a NORMAL margin the two are the same expansion,
  # which is why this changed no existing number.
  #
  # INDEPENDENT margins keep the AXIAL design bit-for-bit: eigen() of a diagonal
  # matrix would reorder the directions by eigenvalue, which changes nothing
  # mathematically but would move every stored number in the tests.
  lam <- numeric(0); dir <- matrix(0, 0L, 0L); h <- numeric(0)
  Xc  <- matrix(numeric(0), 1L, 0L)
  n_cpt <- 1L
  Rc  <- diag(1, nrow = max(dc, 1L))[seq_len(dc), seq_len(dc), drop = FALSE]
  if (dc > 0L) {
    mm <- .admCovDistMoments(cov_dist)
    muc <- mm$mu[cn]
    if (anyNA(muc) || !all(is.finite(muc)))
      bad("covariate(s) ", paste(sQuote(cn[!is.finite(muc)]), collapse = ", "),
          " have no finite mean, which the Taylor expansion expands about.")
    ic <- match(cn, nms)
    vr <- mm$Sigma[cbind(ic, ic)]
    if (anyNA(vr) || !all(is.finite(vr)) || any(vr <= 0))
      bad("covariate(s) ", paste(sQuote(cn[!(is.finite(vr) & vr > 0)]),
                                 collapse = ", "),
          " have no spread, so every design point would coincide. A covariate ",
          "that does not vary should be supplied as a fixed `cov` value instead.")
    Rc  <- if (is.null(Rz)) diag(1, dc) else Rz[!disc, !disc, drop = FALSE]
    lam <- stats::setNames(rep(1, dc), cn)
    dir <- diag(1, nrow = dc)
    if (dc > 1L && max(abs(Rc[lower.tri(Rc)])) > 0) {
      ez <- eigen(Rc, symmetric = TRUE)
      if (min(ez$values) <= -1e-8)
        bad("the covariate correlation implied by `cov_dist` is not positive ",
            "semi-definite (smallest eigenvalue ", format(min(ez$values)),
            "), so it describes no distribution.")
      if (min(ez$values) <= 1e-8)
        bad("covariate(s) ", paste(sQuote(cn), collapse = ", "), " are ",
            "(near-)perfectly collinear: the smallest eigenvalue of their ",
            "correlation matrix is ", format(min(ez$values)), ". The expansion ",
            "would difference along that direction with a step near zero, which ",
            "is catastrophic cancellation, not a derivative. Drop one of the ",
            "covariates, or express the pair as one covariate plus its ",
            "(independent) residual.")
      lam <- pmax(ez$values, 0)
      dir <- ez$vectors
    }
    # THE RADIUS IS A CUBATURE NODE, NOT A FINITE-DIFFERENCE STEP, so it is
    # chosen to match moments rather than to be small. With weights
    # c(+-) = lam/(2h^2) and c(0) = 1 - sum(lam/h^2) the design reproduces
    # E[1], E[z], E[z^2] and E[z^3] at ANY h -- but E[z^4] comes out as h^2/lam
    # against a true 3, so only h = sqrt(3 lam) matches it. There the design
    # coincides with 3-point Gauss-Hermite (nodes 0, +-sqrt(3), weights 2/3,
    # 1/6, 1/6) and is exact through degree 5 instead of degree 3, at the same
    # three points. It is also the classical unscented transform's scaling,
    # kappa = 3 - d.
    #
    # Measured against the exact marginal, at rho = beta*sd_a/omega = 1:
    # relE 9.6e-04 -> 1.4e-05, relV 2.6e-02 -> 5.4e-03.
    #
    # It fixes the conditioning too. The centre weight is 1 - d/h^2, which at
    # the old hfrac = 0.5 was 1 - 4d: -3, -7, -11, -15 for d = 1..4, so the
    # answer was a difference of terms carrying coefficients up to 11x its own
    # size and any solver noise in a design point was amplified by that factor.
    # At the matched radius it is 1 - d/3.
    #
    # `hfrac` (cov_taylor_h) is a MULTIPLIER on that radius, defaulting to 1.
    # Below 1 the points are pulled back toward the covariate mean -- the scaled
    # unscented transform -- which costs the fourth-moment match but keeps the
    # first three, and matters for a strongly skewed margin, where the matched
    # radius sits 1.73 latent SDs out and maps a long way into the tail.
    h     <- hfrac * sqrt(3 * lam)
    n_cpt <- 1L + 2L * dc
    Z     <- matrix(0, n_cpt, dc)
    for (k in seq_len(dc)) {
      Z[2L * k,        ] <-  h[k] * dir[, k]
      Z[2L * k + 1L,   ] <- -h[k] * dir[, k]
    }
    U  <- pmin(pmax(stats::pnorm(Z), 1e-12), 1 - 1e-12)
    Xc <- matrix(vapply(seq_len(dc), function(j)
      .admCovQuantile(cov_dist[[cn[j]]], U[, j]), numeric(n_cpt)),
      n_cpt, dc, dimnames = list(NULL, cn))
  }

  # E_marg = sum_k c_k E_k with sum(c) = 1: the second difference written as a
  # SIGNED weighted sum over the design points. At the default hfrac = 1 and one
  # continuous covariate that is c = (2/3, 1/6, 1/6) -- 3-point Gauss-Hermite.
  ccont       <- numeric(n_cpt)
  ccont[1L]   <- 1 - sum(lam / h^2)
  if (dc > 0L) { ccont[2L * seq_len(dc)] <- lam / (2 * h^2)
                 ccont[2L * seq_len(dc) + 1L] <- lam / (2 * h^2) }

  # ---- assemble: cell SLOWEST, cubature point FASTEST ------------------------
  n_pt <- L * n_cpt
  X <- matrix(0, n_pt, d, dimnames = list(NULL, nms))
  if (dc > 0L) X[, cn] <- Xc[rep(seq_len(n_cpt), times = L), , drop = FALSE]
  if (length(dn)) X[, dn] <- cells[rep(seq_len(L), each = n_cpt), , drop = FALSE]
  cvec <- rep(pcell, each = n_cpt) * rep(ccont, times = L)
  ip   <- if (dc > 0L) as.integer(outer(2L * seq_len(dc), (seq_len(L) - 1L) * n_cpt, "+")) else integer(0)
  im   <- ip + 1L
  ic   <- (seq_len(L) - 1L) * n_cpt + 1L      # the centre row of each cell
  # var carried by each FIRST-difference direction, cell by cell: p_l * lam_j.
  dvar <- if (dc > 0L) rep(pcell, each = dc) * rep(lam, times = L) else numeric(0)
  # ... and by each SECOND difference: p_l * lam_j^2 / 2.
  #
  # Cov_a[y~(a)] was the rank-p term sum_j lam_j g'_j g'_j', i.e. y~(a)
  # LINEARISED about the covariate mean. Keeping the quadratic term of the same
  # expansion, y~ = y~(mu) + g'(a-mu) + (1/2) g''(a-mu)^2, gives for Gaussian
  # latent a
  #
  #   Var(y~) = lam g'^2 + (1/4) g''^2 Var[(a-mu)^2] = lam g'^2 + (1/2) lam^2 g''^2
  #
  # and g'' is ALREADY computed -- it is what produces the mean correction
  # (1/2) lam g''. So the term costs no extra model evaluation, and it has the
  # same shape as everything else here: var times a quadratic form in a linear
  # functional of the solved rows. Measured relV at rho = 1: 5.4e-03 -> 8.3e-04.
  #
  # In more than one direction this is the PURE-direction part of the quadratic
  # contribution. The mixed term sum_{j<k} lam_j lam_k g''_jk g''_jk' needs the
  # cross partials, which an axial design does not carry (that is the
  # 1 + 2p + p(p-1) coordinate-basis rule). Omitting it is the same
  # approximation the design already makes, not a new one.
  dvar2 <- if (dc > 0L) 0.5 * rep(pcell, each = dc) * rep(lam^2, times = L)
           else numeric(0)

  list(X = X, c = cvec, ip = ip, im = im, ic = ic,
       var = dvar, var2 = dvar2, h = rep(h, times = L), mu = X[1L, ],
       n_pt = n_pt, names = nms, dir = dir, Sigma = Rc,
       cont = cn, n_cell = L, p_cell = pcell, n_cpt = n_cpt,
       c_cont = ccont,
       lam = lam)
}

# Expand a design into the row map the moment assembly works from, for a params
# frame that stacks `n_pt` blocks of `nq` eta-quadrature rows (block index
# slowest, eta fastest -- the same stride .admCovGrid's branch uses).
#
#   W    signed combined weights c_k * w_q, so crossprod(W, cp) IS E_marg
#   rows the row indices of each design point, for BLOCK-WISE centring: Vc_k is
#        centred at its own E_k, not at the pooled mean
#   Dw   n_row x m, column m a LINEAR functional of the solved rows, paired with
#        var[m] so that sum_m var_m (Dw_m' cp)(Dw_m' cp)' is the whole rank
#        correction to V. Three kinds of column:
#          - one per (discrete cell, continuous direction): the central
#            difference for g'_j inside that cell, carrying var = p_l * lam_j
#          - one per (discrete cell, continuous direction): the SECOND
#            difference for g''_j, carrying var = p_l * lam_j^2 / 2, which is
#            the quadratic term of Cov_a[y~(a)] and costs no extra solve
#          - one per discrete cell, present only when there IS more than one:
#            E_l - E, carrying var = p_l, which is the BETWEEN-cell term of the
#            law of total variance and is exact
#        Both are quadratic forms in a linear functional of cp, which is why the
#        moment assembly and the analytical gradient handle them with the same
#        two lines and neither needs to know a discrete covariate exists.
.admCovTaylorRows <- function(td, w_eta, nq) {
  n_pt <- td$n_pt
  L    <- td$n_cell %||% 1L
  ncpt <- td$n_cpt  %||% n_pt
  nd   <- length(td$ip)                       # L * (continuous directions)
  dc   <- if (L > 0L) nd %/% L else 0L
  rows <- lapply(seq_len(n_pt), function(k) (k - 1L) * nq + seq_len(nq))
  nb   <- if (L > 1L) L else 0L
  Dw   <- matrix(0, n_pt * nq, 2L * nd + nb)
  for (m in seq_len(nd)) {
    # first difference -> g'_j, paired with var = p_l lam_j
    Dw[rows[[td$ip[m]]], m] <-  w_eta / (2 * td$h[m])
    Dw[rows[[td$im[m]]], m] <- -w_eta / (2 * td$h[m])
    # second difference -> g''_j, paired with var = p_l lam_j^2 / 2. Uses the
    # centre row of THIS cell, not the global centre.
    l <- ((m - 1L) %/% dc) + 1L
    Dw[rows[[td$ip[m]]],   nd + m] <-      w_eta / td$h[m]^2
    Dw[rows[[td$im[m]]],   nd + m] <-      w_eta / td$h[m]^2
    Dw[rows[[td$ic[l]]],   nd + m] <- -2 * w_eta / td$h[m]^2
  }
  if (nb) {
    cc <- td$c_cont
    for (l in seq_len(L)) for (k in seq_len(ncpt)) {
      pt <- (l - 1L) * ncpt + k
      Dw[rows[[pt]], 2L * nd + l] <- w_eta * cc[k]
    }
    # subtract the global weights once, so column l is exactly E_l - E
    for (pt in seq_len(n_pt))
      Dw[rows[[pt]], 2L * nd + seq_len(L)] <-
        Dw[rows[[pt]], 2L * nd + seq_len(L)] - w_eta * td$c[pt]
  }
  list(n_pt = n_pt, nq = nq, rows = rows, w = w_eta, c = td$c,
       var = c(td$var, td$var2, if (nb) td$p_cell else numeric(0)),
       h = td$h, Dw = Dw, names = td$names)
}

# Metadata siblings of the per-covariate specs. Named once, because every
# consumer that enumerated `names(cov_dist)` and forgot them produced a
# different wrong answer: .admCovGrid tried to read `$values` off a numeric
# `rho` ("$ operator is invalid for atomic vectors") and off the `joint`
# CLOSURE ("object of type 'closure' is not subsettable"), and
# .admCheckCovariates reported "declares cov_dist for 'rho', which the model
# never reads".
.ADM_COV_META <- c("rho", "Sigma", "cor", "joint", "latentR")

.admCovSpecNames <- function(cov_dist) setdiff(names(cov_dist), .ADM_COV_META)

# Mean vector and covariance MATRIX of a covariate distribution, on the scale
# the model reads the covariates on.
#
# Independent margins have a diagonal Sigma and closed-form entries. A `joint`
# sampler has neither, so its moments are measured by pushing a deterministic
# Sobol block through the sampler itself -- the same map the fit will use, so
# the moments cannot disagree with the draws.
#
# This is DATA, not a parameter: `cov_dist` does not move during a fit, so the
# whole thing is computed once, up front, and never inside an objective call.
.admCovDistMoments <- function(cov_dist, n_sobol = 8192L) {
  cov_dist <- .admCovDistCanon(cov_dist)
  nms <- .admCovSpecNames(cov_dist)
  d   <- length(nms)
  if (!d) return(NULL)
  mu <- vapply(nms, function(n)
    as.numeric(.admCovMeanOf(cov_dist[[n]]) %||% NA_real_), numeric(1))
  jf <- cov_dist[["joint"]]
  if (!is.function(jf)) {
    vr <- vapply(nms, function(n) .admCovVarOf(cov_dist[[n]]), numeric(1))
    S  <- diag(vr, nrow = d)                 # nrow= : diag(scalar) is an IDENTITY
    dimnames(S) <- list(nms, nms)
    return(list(mu = mu, Sigma = S, names = nms, diagonal = TRUE))
  }
  u <- randtoolbox::sobol(n_sobol, dim = d)
  if (!is.matrix(u)) u <- matrix(u, nrow = n_sobol)
  u <- pmin(pmax(u, .Machine$double.eps), 1 - .Machine$double.eps)
  colnames(u) <- nms
  X <- tryCatch(as.matrix(jf(u)), error = function(e)
    stop("admixr2: cov_dist$joint failed while measuring the covariate ",
         "moments: ", conditionMessage(e), call. = FALSE))
  if (!is.matrix(X) || nrow(X) != n_sobol || is.null(colnames(X)) ||
      !setequal(colnames(X), nms))
    stop("admixr2: cov_dist$joint must return a ", n_sobol, " x ", d,
         " matrix with columns ", paste(sQuote(nms), collapse = ", "), ".",
         call. = FALSE)
  X <- X[, nms, drop = FALSE]
  if (!all(is.finite(X)))
    stop("admixr2: cov_dist$joint returned non-finite covariate values.",
         call. = FALSE)
  m <- colMeans(X)
  Xc <- sweep(X, 2L, m)
  S  <- crossprod(Xc) / nrow(X)              # ML denominator: these are nodes
  dimnames(S) <- list(nms, nms)
  list(mu = m, Sigma = S, names = nms,
       diagonal = max(abs(S[lower.tri(S)])) <=
                  1e-10 * max(diag(S), .Machine$double.eps))
}

# =============================================================================
# Covariate STRATA -- the per-covariate stratify/marginalise split
# =============================================================================
#
# A published source conditions its model on SOME covariates and not others.
# The matching rule is per covariate, not per study:
#
#   STRATIFY   on the covariates the source's own model fitted
#   MARGINALISE over the ones it did not, using the source's covariate
#              distribution -- on BOTH the observed and the predicted side
#
# Building the grid as a full product over ALL covariates for EVERY source is
# the failure this exists to prevent: a model with no term in x is evaluated at
# nodes that vary x, answers "no change" at every one, and that FABRICATED null
# contrast is then scored as evidence against sources that have a real one.
# Measured, with no sampling noise at all, against a true coefficient of 0.450:
# 0.2107 (-53.2%) fabricated versus 0.4500 (-0.0%) matched.
#
# WHY THIS IS A GENERATION-SIDE FUNCTION, NOT AN ESTIMATOR SETTING. Stratifying
# needs per-stratum OBSERVATIONS. Scoring ONE pooled (E, V) against K
# conditional predictions is a different object entirely -- a weighted average
# of log-densities of the same data under different models, i.e. -2 log of an
# unnormalised geometric mean, which is not the likelihood of anything. So
# there is deliberately no "stratify this study" flag on the estimator side:
# strata arrive as ordinary studies, each with its own n and its own data, and
# the estimator needs no knowledge of where they came from (verified
# bit-identical to fitting them as unrelated studies).
#
# CONDITIONING. Within stratum k the covariates the source did NOT fit must be
# marginalised over their distribution CONDITIONAL on that stratum, not their
# marginal one. With correlated covariates the unconditional shortcut predicts
# the average-x2 response at every x1 node while the source's high-x1 subjects
# actually had high x2 -- a mean error that VARIES ACROSS NODES, which is
# precisely the shape that biases the stratified covariate's own coefficient.
# Measured relative error of the block mean: 1e-3 conditional against 0.20 at
# rho = 0.3 and 0.96 at rho = 0.85 unconditional. At rho = 0 they coincide
# exactly, which is why this never showed up in an independent-covariate test.
#
# The conditioning is done in U-SPACE, and that is what makes it work for an
# arbitrary sampler. A copula maps INDEPENDENT uniforms to dependent values, so
# holding the leading uniforms fixed and varying the rest samples the
# conditional distribution -- for a vine those uniforms are its Rosenblatt
# coordinates. Those are the INDEPENDENT uniforms the sampler's
# inverse-Rosenblatt step consumes, not the copula's own margins -- a copula's
# margins are uniform but DEPENDENT, since the copula is the dependence.
# Nothing here has to invert the sampler, and no forward Rosenblatt transform
# is needed, because a stratum is defined BY ITS NODE IN U-SPACE and the
# covariate value at that node is read off the sampler rather than solved for.
#
# The contract on a user `joint`, stated because it cannot be checked: its
# uniform columns must act as a conditioning cascade in DECLARED COVARIATE
# ORDER, so that fixing the columns of the stratified covariates conditions on
# them. inverse_rosenblatt() satisfies this when the vine is ordered with those
# covariates leading; admixr2's own `cor` sampler satisfies it by construction
# (it multiplies by a Cholesky factor in declared order).
.admCovStrata <- function(cov_dist, stratify, n_nodes = 5L,
                          n_pool = 32768L) {
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
  Rm <- cov_dist[["latentR"]]
  if (!is.null(Rm) && all(vapply(nms, function(n)
        is.null(cov_dist[[n]][["values"]]), logical(1)))) {
    # the ordinary product Gauss-Hermite grid -- bit-identical to building it
    # from indices, checked across node counts and dimensions
    .ng <- .adghNodeGrid(n_nodes, length(iS))
    K   <- nrow(.ng$X)
    wk  <- .ng$W / sum(.ng$W)
    zS  <- .ng$X
    iM  <- setdiff(seq_len(d), iS)
    if (length(iM)) {
      A  <- Rm[iM, iS, drop = FALSE] %*% solve(Rm[iS, iS, drop = FALSE])
      Sc <- Rm[iM, iM, drop = FALSE] - A %*% Rm[iS, iM, drop = FALSE]
      Lc <- t(chol(Sc + diag(1e-12, nrow(Sc))))
    }
    nd <- max(as.integer(n_pool), 8192L)
    Z0 <- if (length(iM)) {
      q <- randtoolbox::sobol(nd, dim = length(iM))
      if (!is.matrix(q)) q <- matrix(q, ncol = length(iM))
      stats::qnorm(pmin(pmax(q, 1e-12), 1 - 1e-12))
    } else NULL
    return(lapply(seq_len(K), function(k) {
      Z <- matrix(0, nd, d)
      Z[, iS] <- rep(zS[k, ], each = nd)
      if (length(iM))
        Z[, iM] <- sweep(Z0 %*% t(Lc), 2L, as.numeric(A %*% zS[k, ]), "+")
      U  <- pmin(pmax(stats::pnorm(Z), 1e-12), 1 - 1e-12)
      Xk <- matrix(vapply(seq_len(d), function(j)
        .admCovQuantile(cov_dist[[nms[j]]], U[, j]), numeric(nd)),
        nd, d, dimnames = list(NULL, nms))
      # sorted before it is handed out: the rows arrive in Sobol order and the
      # sampler that draws from them runs another Sobol stream, so indexing
      # rows by that stream correlates the two and returns a structured
      # subsample (measured: sd 13.094 against a true 14.283)
      srt <- if (length(iM)) order(Xk[, iM[1L]]) else seq_len(nd)
      Xk  <- Xk[srt, , drop = FALSE]
      cd_k <- stats::setNames(lapply(nms, function(n)
        list(quantile = local({ v <- sort(Xk[, n])
          function(u) v[1L + pmin(floor(u * length(v)), length(v) - 1L)] }))),
        nms)
      cd_k[["joint"]] <- local({ sub <- Xk; nms <- nms
        function(u) { i <- 1L + pmin(floor(u[, 1L] * nrow(sub)), nrow(sub) - 1L)
          out <- sub[i, , drop = FALSE]; colnames(out) <- nms; out } })
      list(cov = stats::setNames(as.list(Xk[1L, iS]), nms[iS]),
           cov_dist = cd_k, weight = wk[k], n_pool_cell = nd)
    }))
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
#' @param n_nodes Strata per stratified covariate (default 5). The covariate
#'   is cut into that many EQUIPROBABLE bins, so every stratum carries the same
#'   number of subjects; a discrete covariate ignores this and is cut at its
#'   levels.
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
#' @return A list with one element per stratum, each containing `cov` (the
#'   stratum's mean for every stratified covariate), `cov_dist` (the
#'   distribution of ALL covariates *within* that stratum --- the empirical
#'   conditional, carrying the full joint dependence), `n` and `weight`.
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
covStrata <- function(cov_dist, stratify, n_nodes = 5L, n = 1,
                      n_pool = 32768L) {
  checkmate::assertNumber(n, lower = 0, finite = TRUE)
  st <- .admCovStrata(cov_dist, stratify, n_nodes, n_pool)
  lapply(st, function(s)
    list(cov = s$cov, cov_dist = s$cov_dist, n = n * s$weight,
         weight = s$weight, n_pool_cell = s$n_pool_cell))
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
      if (!is.function(m))
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
                         s[["strata_nodes"]] %||% 5L)
    for (k in seq_along(stl)) {
      sk <- s
      sk[["stratify"]] <- NULL; sk[["strata_nodes"]] <- NULL
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
    if (any(disc)) {
      npool <- max(as.integer(n_nodes)^d, 4096L)
      X <- .admCovRowsFor(cov_dist, npool, 0L)[, nms, drop = FALSE]
      # no latent normal score exists for a discrete pool, and none is wanted:
      # a discrete covariate can never make Delta normal.
      return(list(X = X, W = rep(1 / nrow(X), nrow(X)), z = NULL))
    }
    .ng <- .adghNodeGrid(n_nodes, d)
    z   <- .ng$X
    W   <- .ng$W
    u  <- stats::pnorm(z)
    u  <- pmin(pmax(u, .Machine$double.eps), 1 - .Machine$double.eps)
    colnames(u) <- nms
    X <- tryCatch(as.matrix(jf(u)), error = function(e)
      stop("admixr2: cov_dist$joint failed on the ", nrow(u), " x ", d,
           " grid of uniforms: ", conditionMessage(e), call. = FALSE))
    if (!is.matrix(X) || nrow(X) != nrow(u) || is.null(colnames(X)) ||
        !setequal(colnames(X), nms))
      stop("admixr2: cov_dist$joint must return one ROW per supplied uniform ",
           "row, with columns ", paste(sQuote(nms), collapse = ", "), ".",
           call. = FALSE)
    X <- X[, nms, drop = FALSE]
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

.admShiftNodes <- function(Delta, W, om, n_u, tol = 1e-12, maxit = 30L,
                           ftol = 1e-14, z = NULL, gauss_ok = TRUE) {
  # NON-FINITE INPUT is refused rather than iterated on, and BEFORE the Gaussian
  # test: .admShiftDelta screens its own output, but this is also called
  # directly and per-coordinate from the Rosenblatt recursion, and a NaN reaches
  # `if (max(abs(st)) < tol)` as "missing value where TRUE/FALSE needed" rather
  # than as the NULL every caller is written to handle.
  if (!all(is.finite(Delta)) || !all(is.finite(W)) || !is.finite(om))
    return(NULL)
  g  <- .adghNodes1(n_u)
  tg <- stats::pnorm(g$x)
  mD <- sum(W * Delta); vD <- max(sum(W * Delta^2) - mD^2, 0)
  u  <- mD + sqrt(vD + om^2) * g$x
  # THE GAUSSIAN CASE IS EXACT AND COSTS NOTHING. `u` above is already the
  # answer; the loop would spend up to `maxit` mixture evaluations rediscovering
  # it. Measured on the default grid: 0.02 ms against 1.95 ms, and the Newton
  # route is not merely slower but LESS accurate as sd(Delta)/omega grows --
  # 8.5e-04 at ratio 4 and 8.8e-02 at ratio 16, against <= 1.3e-12 here at every
  # ratio, because a widely-separated mixture is what n_u nodes resolve worst.
  # .admShiftDu applies the matching closed-form derivatives, keyed off the same
  # test, so the objective and the gradient cannot disagree about which node set
  # is in play.
  # The certificate when the latent scores are available, the moment test when
  # they are not. In ONE dimension the moment test is answering exactly the
  # right question -- is this univariate law normal -- so it stays as the
  # fallback for hand-built specs and direct calls. It is not sound above one
  # dimension, which is why .admShiftNodesMulti requires `z` instead.
  # `gauss_ok = FALSE` forces the mixture inversion even where the closed form
  # would apply. The Rosenblatt recursion needs it: its derivatives come from
  # the DISCRETE mixture identity F(u) = Phi(z), and a discretised affine Delta
  # is a quadrature of a normal rather than a normal, so the closed-form u does
  # not satisfy that identity exactly. Nodes from one construction and
  # derivatives from the other is the objective-and-gradient disagreement this
  # file exists to avoid.
  if (gauss_ok && .admShiftGaussOK(Delta, W, z, 1L))
    return(list(u = u, w = g$w, gauss = TRUE))
  # CONVERGENCE IS TESTED ON THE RESIDUAL, NOT ON THE STEP. `max(abs(st)) < tol`
  # never fired: at the outermost node the CDF is saturated (target within
  # 1.3e-14 of 1), so its Newton step stays finite no matter how exactly the
  # root is found, and one such node held the whole vector short of the test.
  # The loop therefore ran its full `maxit` on EVERY call -- measured 30 of 30
  # iterations while the residual it actually cares about had reached 1.1e-16
  # within a few. That is what made node construction cost 1.95 ms at the
  # default grid and 23 ms at cov_nodes = 127.
  for (it in seq_len(maxit)) {
    Z  <- outer(u, Delta, "-") / om
    Fu <- as.numeric(stats::pnorm(Z) %*% W)
    if (max(abs(Fu - tg)) < ftol) break
    fu <- as.numeric(stats::dnorm(Z) %*% W) / om
    st <- (Fu - tg) / pmax(fu, 1e-300)
    # safeguarded: an unbounded Newton step on a mixture CDF can leap past the
    # support of a neighbouring component and stall in a flat tail
    u  <- u - sign(st) * pmin(abs(st), 2 * om)
    if (max(abs(st)) < tol) break
  }
  # A clamped Newton on a multi-modal mixture can exhaust maxit while still far
  # from the root, and returning those nodes silently is a wrong answer wearing
  # a valid shape (measured |u - u_exact| = 3.7 on a well-separated mixture).
  # Report the failure instead; the caller refuses the path.
  Zc <- outer(u, Delta, "-") / om
  if (max(abs(as.numeric(stats::pnorm(Zc) %*% W) - tg)) > 1e-6) return(NULL)
  list(u = u, w = g$w)
}

# d(Delta)/d(theta) for every structural theta, by central difference.
#
# Delta costs no solve -- two vectorised evaluations of the parameter assignment
# -- so a fixed step is right here and a measured one (.admShi21Steps) would be
# pure overhead. The step lives in ONE place: it was written out at three sites,
# and the objective and its gradient must not be able to drift apart on it.
#
# Returns a named list of (n_cov x m) matrices, or NULL if any evaluation fails,
# which every caller treats as "no analytic derivative available".
.admShiftDDelta <- function(spec, struct, X, aref, h = 1e-6) {
  out <- lapply(names(struct), function(k) {
    s1 <- struct; s1[[k]] <- s1[[k]] + h
    s2 <- struct; s2[[k]] <- s2[[k]] - h
    d1 <- .admShiftDelta(spec, s1, X, aref)
    d2 <- .admShiftDelta(spec, s2, X, aref)
    if (is.null(d1) || is.null(d2)) return(NULL)
    (as.matrix(d1) - as.matrix(d2)) / (2 * h)
  })
  stats::setNames(out, names(struct))
}

# Derivatives of the u nodes with respect to the parameters they move with.
#
# u_i is defined implicitly by F_u(u_i) = Phi(z_i) with
# F_u(t) = sum_j W_j Phi((t - Delta_j)/om), so differentiating that identity
# gives both chain factors as CONDITIONAL EXPECTATIONS under the same covariate
# quadrature that built F_u -- no solve, no extra nodes:
#
#   du_i/dtheta = E[ dDelta/dtheta | u = u_i ]
#   du_i/dom    = E[ (u_i - Delta)/om | u = u_i ]
#
# dDelta/dtheta is finite-differenced rather than derived symbolically because
# Delta costs no solve: two vectorised evaluations of the parameter assignment
# per theta, against one ODE solve for anything else.
#
# The typical value of a mu-referenced parameter CANCELS out of Delta (it is a
# difference of the same expression at two covariate values), so its column
# comes back zero and the chain correctly contributes nothing through u.
.admShiftDu <- function(spec, struct, X, aref, Delta, W, om, u, h = 1e-6,
                        z = NULL) {
  # Scalar shift only. With m > 1 the later coordinates' nodes also move through
  # the POSTERIOR weights of the Rosenblatt recursion, which is a second chain
  # this does not carry -- .adghGrad finite-differences that case instead.
  Delta <- as.numeric(as.matrix(Delta)[, 1L])
  nm <- names(struct)
  # dDelta/dtheta by central difference, shared by both branches below.
  dDs <- lapply(.admShiftDDelta(spec, struct, X, aref, h),
                function(m) if (is.null(m)) NULL else m[, 1L])

  # THE GAUSSIAN BRANCH, matched to .admShiftNodes' closed form. Keyed off the
  # SAME predicate, evaluated on the same (Delta, W): if the two ever disagreed
  # the gradient would describe a node set the objective is not using, which is
  # the silent class of failure this file exists to avoid. Here
  # u_i = mD + s*x_i with s = sqrt(vD + om^2), so
  #   du_i/dom    = (om/s) * x_i
  #   du_i/dtheta = dmD + (x_i/(2s)) * dvD
  # with dmD = sum(W dDelta) and dvD = 2*sum(W Delta dDelta) - 2*mD*dmD.
  # THE SAME PREDICATE .admShiftNodes USED, on the same inputs -- which requires
  # the same `z`. Testing moments here while the nodes were chosen by the affine
  # certificate is not equivalent: Gauss-Hermite on n nodes is exact to degree
  # 2n-1 and this checks degrees 3..6, so at cov_nodes = 3 an exactly affine
  # Delta certifies (9.5e-15) and FAILS the moment test (4.0e-01). The nodes were
  # then the closed form while the derivatives differentiated the mixture
  # identity at them, and the analytic gradient came back 4.7e-03 off a central
  # difference against 2.1e-09 at 4 nodes or more.
  if (.admShiftGaussOK(matrix(Delta, ncol = 1L), W, z, 1L)) {
    Wn <- W / sum(W)
    mD <- sum(Wn * Delta); vD <- max(sum(Wn * Delta^2) - mD^2, 0)
    s  <- sqrt(vD + om^2)
    x  <- (u - mD) / s
    du_dom <- (om / s) * x
    dth <- vapply(nm, function(k) {
      dD <- dDs[[k]]
      if (is.null(dD)) return(rep(NA_real_, length(u)))
      dmD <- sum(Wn * dD)
      dvD <- 2 * sum(Wn * Delta * dD) - 2 * mD * dmD
      dmD + x * dvD / (2 * s)
    }, numeric(length(u)))
    if (!is.matrix(dth)) dth <- matrix(dth, length(u), length(nm))
    dimnames(dth) <- list(NULL, nm)
    return(list(du_dtheta = dth, du_domega = du_dom))
  }

  # Both derivatives are conditional expectations over the mixture components.
  #
  # A central difference of the NODES is not a valid check on them at the
  # outermost node: its target probability sits 1.3e-14 from 1, where the
  # mixture CDF is a numerical plateau -- moving u by 1e-5 there changes F by
  # exactly 0 -- so u is undetermined across a wide interval and its difference
  # quotient is noise. Nodes 3..20 match a central difference to 1.0000; node 1
  # does not, and the analytic value is the trustworthy one. Do not "fix" that
  # gap, and do not use a node FD as the reference for it.
  Z   <- outer(u, Delta, "-") / om
  Ph  <- stats::dnorm(Z)
  den <- pmax(as.numeric(Ph %*% W), 1e-300)
  du_dom <- as.numeric((Ph * Z) %*% W) / den
  dth <- vapply(nm, function(k) {
    dD <- dDs[[k]]
    if (is.null(dD)) return(rep(NA_real_, length(u)))
    as.numeric((Ph %*% (W * dD)) / den)
  }, numeric(length(u)))
  if (!is.matrix(dth)) dth <- matrix(dth, length(u), length(nm))
  dimnames(dth) <- list(NULL, nm)
  list(du_dtheta = dth, du_domega = du_dom)
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
# Quadrature over a VECTOR shift, u = Delta(a) + eta in R^m, WITH derivatives.
#
# u's law is the mixture sum_j W_j N(Delta_j, diag(om^2)). The affected etas are
# mutually uncorrelated (a correlated Omega absorbs instead, or is refused), so
# the components factor given the mixture index -- but not marginally, because
# they share it. Rosenblatt handles exactly that:
#
#   u_1        ~ sum_j W_j N(Delta_j1, om_1^2)                    invert its CDF
#   u_2 | u_1  ~ sum_j W_j^(u_1) N(Delta_j2, om_2^2),
#                W_j^(v) proportional to W_j phi((v - Delta_j1)/om_1)
#
# so each conditional is again a ONE-dimensional mixture and the same Newton
# inversion serves every level. Cost is n_u^m nodes and n_u^(m-1) inversions,
# against n_cov^p * n_node^m for the product grid -- independent of the number
# of covariates, which is the point of the shift.
#
# THE DERIVATIVES ARE CARRIED ALONGSIDE, which is what this used to lack. Every
# u_k is defined implicitly by F_k(u_k; W^(k)) = Phi(z), so differentiating that
# identity gives, for a direction moving Delta by dD and om by dom,
#
#   du_k = sum_j p_j (dDelta_jk + z_j dom_k)  -  om_k sum_j dW_j Phi(z_j) / SA
#
# with p_j the posterior W_j phi(z_j)/SA and SA its normaliser. The SECOND term
# is the chain the old implementation did not have: the later coordinates move
# because the posterior weights move, not only because Delta does. Propagating
# it needs the weight derivative itself,
#
#   dA_j = phi_j (dW_j - W_j z_j dz_j),   dW'_j = (dA_j - W'_j sum_l dA_l) / SA
#
# formed on the UNNORMALISED A_j = W_j phi_j so that a zero weight never divides.
#
# Without this the whole objective had to be finite-differenced, which is a cost
# reserved for a failed sensitivity model rather than a routine path.
.admShiftNodesMultiD <- function(D, W, om, n_u, dirs = NULL) {
  D <- as.matrix(D); m <- ncol(D)
  nd <- length(dirs)
  if (m == 1L) {
    un <- .admShiftNodes(D[, 1L], W, om[1L], n_u, gauss_ok = FALSE)
    if (is.null(un)) return(NULL)
    out <- list(u = matrix(un$u, ncol = 1L), w = un$w, gauss = un$gauss)
    if (nd) out$du <- .admShiftDuLevel(D, W, om, un$u, 1L, dirs,
                                       matrix(0, nrow(D), nd))$du
    return(out)
  }
  rec <- function(k, Wc, dWc) {
    un <- .admShiftNodes(D[, k], Wc, om[k], n_u, gauss_ok = FALSE)
    if (is.null(un)) return(NULL)
    lv <- .admShiftDuLevel(D, Wc, om, un$u, k, dirs, dWc)
    if (k == m)
      return(list(u = matrix(un$u, ncol = 1L), w = un$w, du = lv$du))
    ou <- list(); ow <- numeric(0); od <- list()
    for (i in seq_along(un$u)) {
      sub <- rec(k + 1L, lv$Wnext[[i]], lv$dWnext[[i]])
      if (is.null(sub)) return(NULL)
      nr <- nrow(sub$u)
      ou[[i]] <- cbind(un$u[i], sub$u)
      ow <- c(ow, un$w[i] * sub$w)
      if (nd) {
        a <- array(0, c(nr, m - k + 1L, nd))
        a[, 1L, ] <- rep(lv$du[i, 1L, ], each = nr)
        a[, -1L, ] <- sub$du
        od[[i]] <- a
      }
    }
    list(u = do.call(rbind, ou), w = ow,
         du = if (nd) .admBindDu(od) else NULL)
  }
  r <- rec(1L, W, if (nd) matrix(0, nrow(D), nd) else NULL)
  if (is.null(r)) return(NULL)
  list(u = r$u, w = r$w / sum(r$w), du = r$du)
}

# One Rosenblatt level: the nodes' derivatives, and the posterior weights (with
# THEIR derivatives) that the next level conditions on.
.admShiftDuLevel <- function(D, Wc, om, u, k, dirs, dWc) {
  nd <- length(dirs); nu <- length(u)
  du <- if (nd) array(0, c(nu, 1L, nd)) else NULL
  Wn <- vector("list", nu); dWn <- vector("list", nu)
  Z  <- outer(u, D[, k], "-") / om[k]
  Ph <- stats::dnorm(Z); CP <- stats::pnorm(Z)
  for (i in seq_len(nu)) {
    A  <- Wc * Ph[i, ]; SA <- sum(A)
    if (!is.finite(SA) || SA <= 0) { SA <- .Machine$double.xmin; }
    p  <- A / SA
    Wn[[i]] <- p
    if (!nd) next
    z <- Z[i, ]
    for (d in seq_len(nd)) {
      dD <- dirs[[d]]$dD; dom <- dirs[[d]]$dom
      duv <- sum(p * dD[, k]) + sum(p * z) * dom[k] -
             om[k] * sum(dWc[, d] * CP[i, ]) / SA
      du[i, 1L, d] <- duv
      dz  <- (duv - dD[, k]) / om[k] - z * dom[k] / om[k]
      dA  <- Ph[i, ] * (dWc[, d] - Wc * z * dz)
      dWn[[i]] <- cbind(dWn[[i]], (dA - p * sum(dA)) / SA)
    }
  }
  list(du = du, Wnext = Wn, dWnext = dWn)
}

# rbind a list of (rows x cols x dir) arrays along rows.
.admBindDu <- function(lst) {
  d <- dim(lst[[1L]]); nr <- sum(vapply(lst, function(a) dim(a)[1L], 0L))
  out <- array(0, c(nr, d[2L], d[3L])); off <- 0L
  for (a in lst) {
    n <- dim(a)[1L]
    out[off + seq_len(n), , ] <- a
    off <- off + n
  }
  out
}

# The covariate absorbed into Omega.
#
# When Delta = c + B z, the joint law of (u_S, eta_O) is normal with covariance
# Omega + P -- P zero except P_SS = B B' -- and mean c on the affected block.
# The covariate then needs no node set of its own: it is a rank-<= p additive
# term on Omega, and the ORDINARY eta grid carries it, correlations with the
# etas the covariate never touches included. That is what the retired collapse
# path computed as Omega + J Sigma_a J', now with a certificate rather than a
# muRefCovariateDataFrame lookup, and without its normality-of-the-covariate,
# solve-at-the-mean or grad = "none" preconditions.
#
# Used only where the column substitution cannot go -- a correlated Omega --
# because it gives up the analytic shift gradient (d(eta)/d(L_ab) is a Cholesky
# differential of chol(Omega + P), not the column the existing chain reads), and
# .adghGrad finite-differences instead.
# Omega + B B' and the mean shift c, with no node grid. This is the whole of the
# absorption; .admShiftAbsorb adds a grid on top, and adfo/adirmc need only this.
.admAbsorbFit <- function(D, W, z, omega, j) {
  D <- as.matrix(D)
  if (!.admShiftGaussOK(D, W, z, ncol(D))) return(NULL)
  Wn <- W / sum(W); sw <- sqrt(Wn)
  q  <- tryCatch(qr(cbind(1, z) * sw), error = function(e) NULL)
  if (is.null(q)) return(NULL)
  co <- qr.coef(q, D * sw)
  if (!all(is.finite(co))) return(NULL)
  B  <- t(co[-1L, , drop = FALSE])
  Om <- as.matrix(omega)
  Om[j, j] <- Om[j, j] + tcrossprod(B)
  Lt <- tryCatch(t(chol(Om)), error = function(e) NULL)
  if (is.null(Lt)) return(NULL)
  list(Lt = Lt, omega = Om, cc = as.numeric(co[1L, ]), B = B,
       q = q, sw = sw, j = j)
}

.admShiftAbsorb <- function(D, W, z, omega, j, n_u, nn0) {
  ab <- .admAbsorbFit(D, W, z, omega, j)
  if (is.null(ab)) return(NULL)
  Om <- ab$omega; L <- ab$Lt
  # The affected axes carry the covariate's spread as well as the eta's, so
  # they keep the wider node count the shift sized for them; the rest keep
  # theirs. Node COUNTS, not the covariance, so this cannot move with omega.
  nn <- rep(as.integer(nn0), nrow(Om)); nn[j] <- as.integer(n_u)
  gs <- lapply(nn, .adghNodes1)
  ix <- as.matrix(expand.grid(lapply(nn, seq_len), KEEP.OUT.ATTRS = FALSE))
  X  <- vapply(seq_along(nn), function(k) gs[[k]]$x[ix[, k]],
               numeric(nrow(ix)))
  W2 <- Reduce(`*`, lapply(seq_along(nn), function(k) gs[[k]]$w[ix[, k]]))
  mu <- numeric(nrow(Om)); mu[j] <- ab$cc
  c(ab, list(eta = sweep(X %*% t(L), 2L, mu, "+"), W = W2 / sum(W2), X = X))
}

# d(chol(S))/d(param), for S = Lt Lt'.
#
#   dLt = Lt * Phi(Lt^-1 dS Lt^-T),   Phi(M) = lower triangle, diagonal halved
#
# Verified against a central difference of chol() at 5e-9 relative, both for a
# general symmetric direction and for the dS = E_ij L' + L E_ij' that an omega
# Cholesky entry produces.
.admCholDiff <- function(Lt, dS) {
  Li <- backsolve(Lt, diag(nrow(Lt)), upper.tri = FALSE)
  M  <- Li %*% dS %*% t(Li)
  M[upper.tri(M)] <- 0
  diag(M) <- diag(M) / 2
  Lt %*% M
}

# The absorption's eta-path derivatives wrt the structural thetas.
#
# theta moves eta twice over: through mu = c and through Omega + B B'. Both come
# from d(Delta)/d(theta), which costs no solve -- two vectorised evaluations of
# the parameter assignment, exactly as .admShiftDu pays for the scalar shift.
.admShiftAbsorbDeriv <- function(spec, struct, Xcov, aref, ab, n_eta, h = 1e-6) {
  nm <- names(struct)
  dmu <- matrix(0, n_eta, length(nm), dimnames = list(NULL, nm))
  dP  <- stats::setNames(vector("list", length(nm)), nm)
  dDs <- .admShiftDDelta(spec, struct, Xcov, aref, h)
  for (k in nm) {
    dD <- dDs[[k]]
    if (is.null(dD)) return(NULL)
    dco <- qr.coef(ab$q, dD * ab$sw)
    if (!all(is.finite(dco))) return(NULL)
    dB  <- t(dco[-1L, , drop = FALSE])
    dmu[ab$j, k] <- dco[1L, ]
    M <- matrix(0, n_eta, n_eta)
    M[ab$j, ab$j] <- tcrossprod(dB, ab$B) + tcrossprod(ab$B, dB)
    dP[[k]] <- M
  }
  list(dmu = dmu, dP = dP)
}

# d(f)/d(param) for the absorption: eta = mu + X Lt', so every eta dimension
# moves, not just one column. The existing chain is the special case
# dLt = E_ij (and dmu = 0), which collapses this sum to Jl[[i]] * X[, j].
# d(f)/d(param) for a vector shift: sum the per-coordinate node derivatives
# against that coordinate's own sensitivity column.
.admShiftBase <- function(Jl, idx, du) {
  out <- NULL
  for (a in seq_along(idx)) {
    cl <- du[, a, 1L]
    if (all(cl == 0)) next
    out <- if (is.null(out)) Jl[[idx[a]]] * cl else out + Jl[[idx[a]]] * cl
  }
  out
}

.admAbsorbBase <- function(Jl, X, dLt, dmu) {
  Xd <- X %*% t(dLt)
  out <- NULL
  for (a in seq_along(Jl)) {
    cl <- Xd[, a] + dmu[a]
    if (all(cl == 0)) next
    out <- if (is.null(out)) Jl[[a]] * cl else out + Jl[[a]] * cl
  }
  out
}

# Does the shift identity actually hold for THIS model? Compiled model, a few
# rows, once per fit. Returns the worst relative discrepancy.
.admShiftVerify <- function(spec, ui, rxMod, pinfo, s, out_var, cov_names,
                            aref, cores = 1L) {
  if (is.null(out_var)) out_var <- tryCatch(.admOutputVar(ui), error = function(e) NULL)
  if (is.null(rxMod))   rxMod   <- tryCatch(.admLoadModel(ui), error = function(e) NULL)
  if (is.null(rxMod) || is.null(out_var)) return(NA_real_)
  cd <- s[["cov_dist"]]
  # .admCovVarOf, not .admCovSdOf: the latter is deliberately the LOG-scale
  # spread for a lognormal margin, so pairing it with a natural-scale mean
  # moved the probe by 0.23 kg instead of 18 and a real violation registered
  # ~80x smaller than it is.
  hi <- vapply(cov_names, function(n) {
    m <- .admCovMeanOf(cd[[n]]); v <- .admCovVarOf(cd[[n]])
    if (is.null(m) || !is.finite(m)) return(NA_real_)
    if (is.null(v) || !is.finite(v) || v <= 0) m else m + sqrt(v)
  }, 0)
  if (anyNA(hi)) return(NA_real_)
  X <- matrix(hi, 1L, length(cov_names), dimnames = list(NULL, cov_names))
  struct <- .admShiftStruct(pinfo)
  # .admSimulate indexes struct_theta with `[nm]` and writes it into a params
  # MATRIX column, so it needs the named numeric vector of ESTIMATED thetas --
  # the full list (which also carries fix()ed ones, for evaluating Delta in R)
  # is not addressable that way and silently fails the solve.
  sv <- pinfo$struct_init
  D <- .admShiftDelta(spec, struct, X, aref)
  if (is.null(D)) return(NA_real_)
  # THE PROBE MUST MAKE Delta MATERIAL, or the identity holds VACUOUSLY. A
  # covariate coefficient initialised at 0 gives Delta == 0, at which
  # f(a, eta) == f(a_ref, eta) for ANY model -- including ones the shift does
  # not describe. The optimizer then walks the coefficient off zero and the
  # objective is wrong in the direction that looks better (measured on an
  # ADDITIVE effect granted at b1 = 0: -211 units at 0.002, -70622 at 0.05).
  # So if Delta is negligible at the initial values, re-probe at a nominal
  # coefficient; the identity is a property of the model's STRUCTURE, not of
  # the parameter value, so any material Delta tests the same thing.
  if (max(abs(D)) < 1e-6) {
    st2 <- struct
    for (k in names(st2)) {
      if (isTRUE(abs(st2[[k]]) < 1e-8)) st2[[k]] <- 0.5
    }
    D2 <- .admShiftDelta(spec, st2, X, aref)
    if (is.null(D2) || max(abs(D2)) < 1e-6) return(NA_real_)
    struct <- st2; D <- D2
    sv <- utils::modifyList(as.list(sv), st2[intersect(names(st2), names(sv))])
    sv <- unlist(sv[names(pinfo$struct_init)])
  }
  j  <- match(spec$eta, pinfo$eta_col_names)
  if (anyNA(j)) return(NA_real_)
  es <- c(0, 0.4, -0.6)
  worst <- 0
  for (e in es) {
    e1 <- matrix(0, 1L, pinfo$n_eta); e1[1L, j] <- e
    e2 <- e1; e2[1L, j] <- e + as.numeric(D[1L, ])
    colnames(e1) <- colnames(e2) <- pinfo$eta_col_names
    sA <- s; sA$cov_rows <- X
    sR <- s; sR$cov_rows <- matrix(unlist(aref[cov_names]), 1L, length(cov_names),
                                   dimnames = list(NULL, cov_names))
    pm <- .admMakeParamsList(1L, pinfo, 1L)[[1L]]
    p1 <- tryCatch(.admSimulate(rxMod, sv, pinfo$sigma_names, e1, sA, out_var,
                                pm, cores, .Machine$integer.max, pinfo$sigdig),
                   error = function(e) NULL)
    p2 <- tryCatch(.admSimulate(rxMod, sv, pinfo$sigma_names, e2, sR, out_var,
                                pm, cores, .Machine$integer.max, pinfo$sigdig),
                   error = function(e) NULL)
    if (is.null(p1) || is.null(p2)) return(NA_real_)
    worst <- max(worst, max(abs(p1 - p2) / pmax(abs(p1), 1e-12)))
  }
  worst
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
    # already uses the population form. .admCovVarOf documents the same point.
    q <- spec$quantile((seq_len(1024L) - 0.5) / 1024L)
    return(sqrt(mean((q - mean(q))^2)))
  }
  if (!is.null(spec$meanlog)) return(spec$sdlog)
  spec$sd
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
  studies <- lapply(studies, function(s) {
    if (!is.null(s[["cov_dist"]]))
      s[["cov_dist"]] <- tryCatch(.admCovDistCanon(s[["cov_dist"]]),
                                  error = function(e) s[["cov_dist"]])
    s
  })
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

# Refuse the removed node-quadrature inputs, loudly.
#
# gl/gh/taylor covariate marginalisation was removed (see the file header). Old
# code passing `weight`, `cov_method` or a `quadrature` attribute would
# otherwise be accepted and silently fitted as an ordinary unweighted study
# list, which is a DIFFERENT objective -- a plausible number from a model the
# user did not ask for, which is the failure mode this package keeps meeting.
.admRefuseNodeStudies <- function(studies) {
  # A multi-output study can carry these per observation, so look there too.
  .fields <- function(s, f)
    c(list(s[[f]]), lapply(s$observations %||% list(), `[[`, f))
  bad_w <- vapply(studies, function(s)
    any(vapply(.fields(s, "weight"), function(x)
      !is.null(x) && is.numeric(x) && length(x) &&
        !isTRUE(all.equal(as.numeric(x)[[1L]], 1)), logical(1))), logical(1))
  bad_m <- vapply(studies, function(s)
    any(vapply(.fields(s, "cov_method"), function(x)
      !is.null(x) && !identical(as.character(x)[[1L]], "marginal"),
      logical(1))), logical(1))
  if (!any(bad_w) && !any(bad_m)) return(invisible(NULL))
  which_f <- if (any(bad_m)) "cov_method" else "weight"
  stop("admixr2: node-quadrature covariate marginalisation (gl / gh / taylor) ",
       "was removed.\n",
       "  Stud", if (sum(bad_w | bad_m) > 1L) "ies " else "y ",
       paste(sQuote(names(studies)[bad_w | bad_m]), collapse = ", "),
       " set `", which_f, "`.\n",
       "  For ONE pooled (E, V) per study -- what a publication reports -- give ",
       "the study its\n  `cov_dist` and let admixr2 marginalise over it.\n",
       "  For summaries BY COVARIATE STRATUM, pass the strata as ordinary ",
       "studies, each with\n  its own `n` and its own `cov_dist` -- that ",
       "stratum's own mean and SD, which is\n  what a subgroup table reports. ",
       "That is the same likelihood the node route\n  computed, with the real ",
       "stratum sizes instead of quadrature weights standing in\n  for them. Do ",
       "NOT give a stratum a point `cov`: plugging in the stratum mean is\n  the ",
       "ecological plug-in, and it biases the coefficient upward (+17% at 2 ",
       "strata,\n  +4.3% at 4, for a lognormal covariate).",
       call. = FALSE)
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
#'   on the covariate's own scale (lognormal by default, `dist = "normal"` for a
#'   normal margin), or as a `quantile` function, or as `values` (with optional
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

.admCovSpecFromVec <- function(v, nm) {
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
  # anything else named is read as a CATEGORICAL covariate: the names are the
  # levels' labels and the values their proportions, which is how a baseline
  # table reports sex, a genotype or a dosing band
  if (any(v < 0)) bad("has negative probabilities, reading it as categorical.")
  if (!sum(v) > 0) bad("has probabilities summing to zero.")
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
#' @param joint An arbitrary sampler, for dependence a correlation cannot
#'   express. It receives the matrix of uniforms admixr2 supplies and returns
#'   one named column per covariate --- the shape a copula or an R-vine
#'   produces. Overrides `cor`.
#' @param dist Default margin for the `c(mean = , sd = )` form: `"normal"`
#'   (default) or `"lnorm"`. A per-covariate `dist` wins over it.
#'
#'   Choose `"lnorm"` when the model needs the covariate to stay positive. A
#'   normal margin is unbounded below and the quadrature reaches 3.75 standard
#'   deviations, so any covariate with a coefficient of variation above about
#'   0.19 --- which includes most of them, weight at 72 +/- 16 being 0.22 ---
#'   gets a node at or below zero, where an allometric or log term is `NaN`.
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
      lapply(seq_along(a), function(i) .admCovSpecFromVec(a[[i]], names(a)[i])),
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

# -- Correlated Omega through a NON-certified shift ----------------------------
#
# The column substitution rebuilds the unaffected etas from the DIAGONAL, which
# drops every Omega off-diagonal -- including between two etas the covariate
# never touches. Absorption avoids that, but only when Delta certifies as
# Gaussian. Everything else used to fall back to the product grid, at
# n_cov^p * n_node^m: exponential in the number of covariates.
#
# It need not. Conditional on the covariate node j, shifting by a constant does
# not change the covariance, so
#
#     (u_S, eta_O) | j  ~  N( (Delta_j, 0), Omega )
#
# -- a mixture of normals ALL SHARING ONE Omega, differing only in the S-block
# mean. Condition the other way round and the mixture collapses into one block:
#
#     eta_O          ~ N(0, Omega_OO)                    ordinary grid
#     u_S | eta_O, j ~ N(Delta_j + K eta_O, Sigma_c)     K = Omega_SO Omega_OO^-1
#                                                        Sigma_c = Omega_SS - K Omega_OS
#
# Rotating w = Ls^-1 (u_S - K eta_O) with Ls = chol(Sigma_c) gives
# w | j ~ N(Ls^-1 Delta_j, I): a UNIT-covariance mixture, which is exactly what
# .admShiftNodes / .admShiftNodesMultiD already invert. Two consequences make
# this cheap rather than a new path -- w's law does not involve eta_O, so the
# inversion still runs ONCE, and the node count is unchanged.
#
# Reduces to the existing construction exactly when Omega is diagonal:
# Omega_SO = 0 gives K = 0 and Sigma_c = diag(omega_S^2).
#
# Verified against 4e6 Monte Carlo draws on the hardest cell -- m = 2 shifted
# etas, one untouched, Delta quadratic in a discrete covariate, correlated
# Omega -- at the MC floor (1.2e-3 against a 1.0e-3 reference error), where the
# diagonal shortcut is 0.217 off and zeroes Cov(u_S, eta_O) outright.
.admCondShiftParts <- function(omega, j) {
  omega <- as.matrix(omega)
  n     <- nrow(omega)
  O     <- setdiff(seq_len(n), j)
  Oss   <- omega[j, j, drop = FALSE]
  if (!length(O))
    return(list(O = O, K = NULL, Ls = tryCatch(t(chol(Oss)),
                                               error = function(e) NULL),
                Lo = NULL, Sc = Oss))
  Ooo <- omega[O, O, drop = FALSE]
  Oso <- omega[j, O, drop = FALSE]
  Lo  <- tryCatch(t(chol(Ooo)), error = function(e) NULL)
  if (is.null(Lo)) return(NULL)
  Ki  <- tryCatch(chol2inv(chol(Ooo)), error = function(e) NULL)
  if (is.null(Ki)) return(NULL)
  K   <- Oso %*% Ki
  Sc  <- Oss - K %*% t(Oso)
  Sc  <- (Sc + t(Sc)) / 2
  Ls  <- tryCatch(t(chol(Sc)), error = function(e) NULL)
  if (is.null(Ls)) return(NULL)
  list(O = O, K = K, Ls = Ls, Lo = Lo, Sc = Sc, Ooo_inv = Ki)
}

# Derivatives of the conditional pieces wrt ONE perturbation dOmega of Omega.
#
#   dK   = (dOmega_SO - K dOmega_OO) Omega_OO^-1
#   dSc  = dOmega_SS - dK Omega_OS - K dOmega_OS
#   dLo  = cholDiff(Lo, dOmega_OO)      dLs = cholDiff(Ls, dSc)
#
# and the rotated means move because Ls does:  dDw = -Ls^-1 dLs Dw.
.admCondShiftDeriv <- function(cp, omega, j, dOm, Dw) {
  omega <- as.matrix(omega); dOm <- as.matrix(dOm)
  dSS <- dOm[j, j, drop = FALSE]
  if (!length(cp$O)) {
    dLs <- .admCholDiff(cp$Ls, dSS)
    return(list(dLo = NULL, dK = NULL, dLs = dLs,
                dDw = -t(solve(cp$Ls, dLs %*% t(Dw)))))
  }
  O   <- cp$O
  dOO <- dOm[O, O, drop = FALSE]
  dSO <- dOm[j, O, drop = FALSE]
  dLo <- .admCholDiff(cp$Lo, dOO)
  dK  <- (dSO - cp$K %*% dOO) %*% cp$Ooo_inv
  dSc <- dSS - dK %*% t(omega[j, O, drop = FALSE]) - cp$K %*% t(dSO)
  dSc <- (dSc + t(dSc)) / 2
  dLs <- .admCholDiff(cp$Ls, dSc)
  list(dLo = dLo, dK = dK, dLs = dLs,
       dDw = -t(solve(cp$Ls, dLs %*% t(Dw))))
}

# d(f)/d(direction) for a grid whose eta is NOT X L'.
#
# The omega chain in .adghGradNLL forms d(f)/d(L_ab) as Jl[[a]] * X[, b], which
# assumes eta = X L'. Under conditioning it is
#   eta_O = X_O Lo',  eta_S = w Ls' + eta_O K'
# so every eta column responds to every direction and the contraction has to be
# taken in full -- the same shape .admAbsorbBase uses for the absorption.
.admShiftCondBase <- function(Jl, dEta) {
  out <- NULL
  for (a in seq_along(Jl)) {
    if (all(dEta[, a] == 0)) next
    out <- if (is.null(out)) Jl[[a]] * dEta[, a] else out + Jl[[a]] * dEta[, a]
  }
  out
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

# -- Discrete covariates on the shift path -------------------------------------
#
# A discrete covariate reaching the shifted argument used to disqualify the
# shift outright, on the grounds that u = Delta(a) + eta becomes a multi-modal
# mixture that Gauss-Hermite nodes cannot resolve. The diagnosis is right and
# the remedy was too broad: it dropped the WHOLE study onto the product grid,
# including any continuous covariates riding along, which then cost n_cov^p
# again -- the very thing the shift exists to avoid.
#
# .admShiftNodes places nodes by inverting the mixture CDF at Gauss-Hermite
# probability points but KEEPS the Gaussian weights, so it is exact for a normal
# target and degrades as the target departs -- 8.5e-04 at sd(Delta)/omega = 4,
# 8.8e-02 at 16. A well-separated two-component mixture is the worst case.
#
# Stratifying fixes it at the source. Condition on the discrete levels, which
# .admCovGrid already enumerates EXACTLY with their true probabilities: within a
# stratum only the continuous covariates vary, so the sub-mixture is the mild
# one the quadrature handles well, and the strata recombine with weights pi_l.
# Cost is K * n_u nodes -- linear in the number of discrete CELLS and still
# constant in the number of continuous covariates.
#
# Which grid rows share a discrete cell. NULL when no covariate is discrete, so
# the caller keeps the single-mixture path unchanged.
.admShiftStrata <- function(cov_dist, X) {
  cd  <- .admCovDistCanon(cov_dist)
  nms <- .admCovSpecNames(cd)
  dsc <- nms[vapply(nms, function(n) !is.null(cd[[n]][["values"]]), logical(1))]
  dsc <- intersect(dsc, colnames(X))
  if (!length(dsc)) return(NULL)
  key <- do.call(paste, c(lapply(dsc, function(n) X[, n]), sep = "\r"))
  match(key, unique(key))
}

# The mixture inversion, run PER discrete cell and recombined.
#
# Each cell is its own quadrature problem with its own n_u nodes; the weights
# carry pi_l so the whole set still integrates to one. Derivatives come back
# stratum by stratum and stack the same way -- d(u)/d(psi) is local to the cell
# a node came from, because the cell's own sub-mixture is what produced it.
.admShiftNodesStrat <- function(D, W, om, n_u, strata, dirs = NULL) {
  if (is.null(strata)) return(.admShiftNodesMultiD(D, W, om, n_u, dirs))
  D  <- as.matrix(D)
  ids <- unique(strata)
  us <- ws <- vector("list", length(ids)); ds <- vector("list", length(ids))
  for (i in seq_along(ids)) {
    r  <- which(strata == ids[i])
    pi <- sum(W[r])
    if (!is.finite(pi) || pi <= 0) next
    dl <- if (is.null(dirs)) NULL else lapply(dirs, function(d) list(
      dD = as.matrix(d$dD)[r, , drop = FALSE], dom = d$dom))
    un <- .admShiftNodesMultiD(D[r, , drop = FALSE], W[r] / pi, om, n_u, dl)
    if (is.null(un)) return(NULL)
    us[[i]] <- un$u; ws[[i]] <- un$w * pi; ds[[i]] <- un$du
  }
  keep <- !vapply(us, is.null, logical(1))
  if (!any(keep)) return(NULL)
  us <- us[keep]; ws <- ws[keep]; ds <- ds[keep]
  du <- if (is.null(dirs) || any(vapply(ds, is.null, logical(1)))) NULL else
    .admBindDu(ds)
  w  <- unlist(ws)
  list(u = do.call(rbind, us), w = w / sum(w), du = du)
}


# The index DIRECTION for a parameter that is not affine in the latent scores.
#
# An ordinary least-squares fit of v on z estimates the AVERAGE DERIVATIVE,
# which for a single-index model is proportional to the index direction whatever
# the link does afterwards. It is only an estimate, and it is not verified here:
# .admCovCollapse verifies the DESIGN instead, by checking that it reproduces
# the parameter law -- which is the property actually needed, rather than a
# proxy for it.
#
# Returns NULL for a link with no first-order signal. A SYMMETRIC one (v even in
# u) has zero average derivative, so this declines it and the product grid
# stands -- correct, if conservative.
.admIndexDir <- function(v, Z) {
  cf <- tryCatch(stats::lm.fit(cbind(1, Z), v)$coefficients,
                 error = function(e) NULL)
  if (is.null(cf) || !all(is.finite(cf))) return(NULL)
  b <- cf[-1L]
  nb <- sqrt(sum(b^2))
  if (!is.finite(nb) || nb <= 0) return(NULL)
  # A CONSTANT has no direction. Guarding on the ratio alone gets this
  # backwards: sd(v) = 0 sends it to infinity and a floating-point 1e-16
  # coefficient is returned as though it meant something.
  sv <- stats::sd(v)
  if (!is.finite(sv) || sv <= 0) return(NULL)
  # and a direction whose signal is swamped is noise, not a direction
  if (nb / sv < 1e-8) return(NULL)
  b
}

# -- Dimension collapse: cost scales with the RANK, not the covariate count ----
#
# p covariates reaching the model through r < p independent scalars make an
# r-DIMENSIONAL integral, whatever p is: the model cannot tell two covariate
# vectors apart when they give every parameter the same value. The product grid
# integrates it in p dimensions at n^p points.
#
# This is the shift's argument with the random effect removed. The shift
# integrates u = Delta(a) + eta and needs an eta to substitute into; here there
# is none, so it integrates Delta(a) alone. A covariate on a parameter with NO
# random effect was the largest case the shift refused, and it is the allometric
# one -- CL and V on weight and creatinine clearance.
#
# THE CERTIFICATE IS THE SAME ONE THE SHIFT ROUTES ON. Every continuous
# covariate is X = F^-1(Phi(z)) from a standard normal latent z, so if each
# covariate-reading assignment is affine in z -- log-affine for the
# multiplicative forms that dominate -- the model depends on z only through
# B'z, with B the p x m matrix of loadings.
#
# THE BASIS IS WHAT MAKES IT GENERIC. Take the SVD of B and keep the
# ORTHONORMAL basis U_r of its column space. Then
#
#     w = U_r' z  ~  N(0, I_r)      exactly, since z ~ N(0, I_p)
#
# so the design is the ORDINARY r-dimensional Gauss-Hermite grid -- no
# covariance to factor, no rescaling -- and z = U_r w is an exact preimage,
# minimum-norm and as good as any other because the model sees only U_r' z.
#
# rank(B) is the whole story: three covariates on one parameter give r = 1,
# three on two parameters r = 2, and three on three separate parameters r = 3,
# where there is nothing to gain and this correctly declines.
#
# Measured against a genuine 21^3 product grid, three covariates on ONE
# parameter at CV = 0.5:
#
#     grid n=7        2401 rows   E 5.48e-10   V 1.26e-06
#     COLLAPSED n=11    77 rows   E 1.18e-11   V 3.27e-08
#
# Returns a design in .admCovGrid's shape, so nothing downstream changes, or
# NULL when it does not apply and the product grid stands.
# Evaluate the covariate-reading assignments at a given set of structural
# thetas, covariate values, eta value and discrete cell.
#
# Standalone rather than a closure inside .admCovCollapse, because the SAME
# evaluation has to be redone on every objective call at the CURRENT thetas --
# see .admCovRefresh() for why.
.admCovProbeAt <- function(pr, st, eta_at, cell, AA) {
  nrw <- nrow(AA)
  ev  <- new.env(parent = asNamespace("rxode2"))
  for (k in names(st)) assign(k, st[[k]], ev)
  # eta_at is a SCALAR for the covariate collapse (etas held, only their
  # invariance is being checked) and a MATRIX for the joint one, where the etas
  # are part of the latent vector being probed
  if (is.matrix(eta_at)) {
    for (j in seq_along(pr$eta_names)) assign(pr$eta_names[j], eta_at[, j], ev)
  } else for (e in pr$eta_names) assign(e, eta_at, ev)
  for (k in pr$cn) assign(k, AA[, k], ev)
  for (k in pr$dn) assign(k, cell[[k]], ev)
  # A study declares a covariate one of two ways: as a DISTRIBUTION to
  # marginalise over, or as the VALUE it is CONDITIONED at. Only the first is an
  # integral, so only the first appears in the design -- but a conditioned
  # covariate in the SAME assignment still has to be in scope, or the probe
  # cannot evaluate it and the whole study is refused. That is a mixed study:
  # marginalising over weight while sitting in a reported age stratum.
  for (k in names(pr$cov_fixed))
    if (!(k %in% pr$cn) && !(k %in% pr$dn)) assign(k, pr$cov_fixed[[k]], ev)
  out <- vector("list", length(pr$hit)); j <- 0L
  for (ii in seq_along(pr$lst)) {
    e <- pr$lst[[ii]]
    if (!isTRUE(pr$is_asgn[ii])) next
    # An assignment that will not evaluate in R -- cp <- linCmt(), an ODE line,
    # anything reaching the solver -- is SKIPPED rather than fatal. It simply
    # does not get defined, and if a covariate-reading assignment needed it,
    # THAT one fails and is caught. Bailing on the first unevaluable line
    # refused every model with a linCmt(), which is most of them.
    v <- tryCatch(eval(e[[3L]], ev), error = function(e) NULL)
    if (is.null(v)) {
      if (ii %in% pr$hit) return(NULL)
      next
    }
    assign(as.character(e[[2L]]), v, ev)
    if (ii %in% pr$hit) {
      j <- j + 1L
      if (length(v) != nrw || !all(is.finite(v))) return(NULL)
      out[[j]] <- as.numeric(v)
    }
  }
  out <- Filter(Negate(is.null), out)
  if (length(out) != length(pr$hit)) return(NULL)
  do.call(cbind, out)
}

# The direction each covariate-reading assignment depends on the latent normal
# through: affine where that holds, single index otherwise.
#
# `routes` replays a decision already made instead of re-deciding it. Two
# reasons, and the second is the important one:
#
#   - COST. This runs on every objective call now, and the affinity test is an
#     lm.fit plus a residual norm per column. Replaying the chosen route is one
#     lm.fit and no test.
#   - CONTINUITY. The test is a threshold. A column sitting near it could be
#     read as affine on one call and as a single index on the next, and the two
#     do not agree to machine precision -- so the objective would step, for no
#     reason the optimizer can see. Which route a column takes is a property of
#     the MODEL, so it is settled once, at admission.
.admCovLoadings <- function(P, Z, pc, routes = NULL) {
  if (is.null(P)) return(NULL)
  np <- nrow(Z)
  W1 <- rep(1 / np, np)
  B  <- matrix(0, pc, ncol(P))
  rt <- vector("list", ncol(P))
  for (k in seq_len(ncol(P))) {
    pk <- P[, k]
    if (stats::sd(pk) <= 0) { rt[[k]] <- "const"; next }
    if (!is.null(routes)) {
      r <- routes[[k]]
      if (identical(r, "const")) next
      b <- if (identical(r, "index")) .admIndexDir(pk, Z) else {
        y  <- if (identical(r, "affine_log")) log(pk) else pk
        cf <- tryCatch(stats::lm.fit(cbind(1, Z), y)$coefficients,
                       error = function(e) NULL)
        if (is.null(cf) || !all(is.finite(cf))) NULL else cf[-1L]
      }
      if (is.null(b)) return(NULL)
      B[, k] <- b
      next
    }
    b <- NULL
    # AFFINE first, on the LOG and then the raw scale. Exact where it holds,
    # and it holds for the multiplicative and allometric forms that dominate.
    for (yi in seq_len(2L)) {
      y <- if (yi == 1L) { if (all(pk > 0)) log(pk) else NULL } else pk
      if (is.null(y)) next
      if (.admShiftAffineResid(matrix(y, ncol = 1L), W1, Z) <
          .ADM_SHIFT_GAUSS_TOL) {
        cf <- tryCatch(stats::lm.fit(cbind(1, Z), y)$coefficients,
                       error = function(e) NULL)
        if (!is.null(cf) && all(is.finite(cf))) {
          b <- cf[-1L]
          rt[[k]] <- if (yi == 1L) "affine_log" else "affine_raw"
          break
        }
      }
    }
    # SINGLE INDEX otherwise. Affine is far stronger than the construction
    # needs: the design places Gauss-Hermite nodes in z, so it is enough that
    # the parameter be SOME function of one linear combination u. Then u is
    # normal, a GH rule integrates the composition exactly to degree 2n-1, and
    # the preimage is unchanged. Affine is the special case of an identity link
    # -- and requiring it refuses, for one, an Emax link on a product of
    # lognormal covariates, whose INDEX is affine and whose LINK is not.
    if (is.null(b)) { b <- .admIndexDir(pk, Z); rt[[k]] <- "index" }
    if (is.null(b)) return(NULL)
    B[, k] <- b
  }
  attr(B, "routes") <- rt
  B
}

# Re-aim the collapsed design at the CURRENT structural thetas.
#
# The rotation depends on them. A covariate coefficient is an ESTIMATED
# parameter, and moving it turns the direction the covariates reach the model
# through -- so a design built once at admission pins the covariates to the line
# the STARTING values implied, and the variation orthogonal to that line is
# missed entirely. Measured against a 15^3 product-grid reference: exact to
# 1e-5 for every non-covariate parameter, and 53 to 163 -2LL units out for a 0.1
# move in ONE coefficient. Every unit and moment test passed throughout, because
# they all evaluate at the initial point, where the cached design is correct by
# construction. It took a real fit to see it.
#
# .adghGrid already recomputes the shift path's Delta from pars$struct on every
# objective call for exactly this reason; this is the collapse's version of it.
# What stays fixed at admission is everything STRUCTURAL -- the rank, the node
# counts, the certificate -- none of which a coefficient's VALUE can change.
# The rank is set by how many assignments read covariates, not by how strongly.
#
# Costs no solves: a 128-point probe, an SVD of a pc x m matrix, and the design
# build. Against an rxSolve at ~11 ms this does not register.
.admCovRefresh <- function(co, st) {
  if (is.null(co) || is.null(co$pr) || is.null(st)) return(co)
  P <- .admCovProbeAt(co$pr, st, 0, co$cell_list[[1L]], co$Ap)
  B <- .admCovLoadings(P, co$Zp, co$pc, co$pr$routes)
  if (is.null(B)) return(co)
  sv <- tryCatch(svd(B), error = function(e) NULL)
  if (is.null(sv) || length(sv$d) < co$r) return(co)
  U  <- sv$u[, seq_len(co$r), drop = FALSE]
  Sr <- t(U) %*% co$Rc %*% U
  Lr <- tryCatch(chol(Sr), error = function(e) NULL)
  if (is.null(Lr)) return(co)
  gl <- lapply(co$nv, function(m) .adghNodes1(m))
  Xg <- as.matrix(expand.grid(lapply(gl, function(g) g$x)))
  Wg <- as.numeric(apply(expand.grid(lapply(gl, function(g) g$w)), 1L, prod))
  dimnames(Xg) <- NULL
  Zc <- Xg %*% Lr %*% t(U)
  Xc <- vapply(seq_len(co$pc), function(k)
    .admCovQuantile(co$cd[[co$cn[k]]],
                    pmin(pmax(stats::pnorm(Zc[, k]), .Machine$double.eps),
                         1 - .Machine$double.eps)), numeric(nrow(Zc)))
  if (!is.matrix(Xc)) Xc <- matrix(Xc, nrow(Zc), co$pc)
  colnames(Xc) <- co$cn
  # A refresh that cannot evaluate leaves the previous design in place. It can
  # only happen where the parameter assignments themselves fail -- a log of a
  # negative, an overflow -- and the solve rejects that region anyway.
  if (!all(is.finite(Xc))) return(co)
  Wc <- Wg / sum(Wg)
  nq <- nrow(Xc); nl <- max(nrow(co$cells), 1L)
  ix <- rep(seq_len(nq), times = nl)
  Xf <- Xc[ix, , drop = FALSE]
  Wf <- Wc[ix] * rep(co$pcell, each = nq)
  if (length(co$dn))
    Xf <- cbind(Xf, co$cells[rep(seq_len(nl), each = nq), , drop = FALSE])
  Xf <- Xf[, co$nms, drop = FALSE]
  co$X <- Xf; co$W <- Wf / sum(Wf); co$z <- Zc[ix, , drop = FALSE]
  co$U <- U; co$Lr <- Lr
  co
}

# How many nodes ONE COLLAPSED DIRECTION deserves.
#
# It is not cov_nodes. A collapsed direction carries the COMBINED spread of the
# pc covariate axes it replaced, so it is wider than any one of them and needs
# proportionally more resolution -- the same reasoning .adghGrid records for the
# shift path's n_u, which is min(101, 4 * nn0) rather than cov_nodes for exactly
# this reason ("fixing n_u at cov_nodes left the shift path ~10x LESS accurate
# than the grid it replaces").
#
# Measured, three lognormal covariates collapsing to one direction, against a
# 21^3 product-grid reference, at cov_nodes = 7:
#
#            design pts   at start   b +0.3    b +0.6    b +1.0
#   grid 7^3        343    3.1e-08   4.5e-05   1.2e-03   8.1e-02
#   collapse 7        7    1.6e-06   2.3e-02   3.6e-01   1.8e+01
#   collapse 21      21    2.4e-10   9.4e-10   5.0e-06   2.6e-02
#
# At the cap it is BETTER than the grid at every point, on 16x fewer design
# points. At cov_nodes it is worse everywhere except the starting values -- and
# the starting values are where every moment test evaluates, which is why this
# survived until a real fit walked away from them.
#
# pc/r is how many axes each direction absorbs on average, so r == pc recovers
# cov_nodes exactly and no collapse claims more than it merged.
.admCovDirNodes <- function(n_nodes, pc, r)
  min(101L, as.integer(ceiling(as.numeric(n_nodes) * pc / max(r, 1L))))

.admCovCollapse <- function(ui, pinfo, cov_dist, n_nodes, n_probe = 128L,
                            max_rows = 20000L, n_ver = 8192L,
                            cov_fixed = NULL) {
  cd  <- .admCovDistCanon(cov_dist)
  nms <- .admCovSpecNames(cd)
  p   <- length(nms)
  if (p < 2L) return(NULL)
  dsc <- vapply(nms, function(n) !is.null(cd[[n]][["values"]]), logical(1))
  cn  <- nms[!dsc]                     # CONTINUOUS: what collapses
  dn  <- nms[dsc]                      # DISCRETE: enumerated, as strata
  pc  <- length(cn)
  # one continuous covariate is already a one-dimensional integral
  if (pc < 2L) return(NULL)
  R <- cd[["latentR"]]
  # An opaque user `joint` publishes no latent structure to project along. The
  # `cor` sampler admixr2 builds itself DOES -- it is a Gaussian copula and
  # records latentR -- so a CORRELATED covariate set is workable, not refused.
  if (is.function(cd[["joint"]]) && is.null(R)) return(NULL)
  # latentR is indexed POSITIONALLY -- it carries no dimnames, and indexing it
  # by covariate name fails outright rather than silently.
  ic <- match(cn, nms); id <- match(dn, nms)
  # A discrete covariate DEPENDENT on a continuous one makes each level a
  # TRUNCATION of the latent normal rather than a point, so the continuous
  # conditional law differs cell to cell and one shared design would be the
  # wrong design in every cell. The same refusal .admCovTaylorDesign carries.
  if (length(dn) && !is.null(R) &&
      any(abs(R[id, ic, drop = FALSE]) > 0)) return(NULL)
  Rc <- if (is.null(R)) diag(1, pc) else R[ic, ic, drop = FALSE]
  Lc <- tryCatch(chol(Rc), error = function(e) NULL)
  if (is.null(Lc)) return(NULL)
  # discrete cells: exact product enumeration with their own probabilities
  if (length(dn)) {
    lev <- lapply(dn, function(n) as.numeric(cd[[n]][["values"]]))
    prb <- lapply(dn, function(n) { sp <- cd[[n]]
      pr <- sp[["probs"]] %||% rep(1 / length(sp[["values"]]),
                                   length(sp[["values"]]))
      pr / sum(pr) })
    cells <- as.matrix(expand.grid(lev, KEEP.OUT.ATTRS = FALSE))
    colnames(cells) <- dn
    pcell <- apply(as.matrix(expand.grid(prb, KEEP.OUT.ATTRS = FALSE)), 1L, prod)
  } else {
    cells <- matrix(numeric(0), 1L, 0L); pcell <- 1
  }
  cell_list <- lapply(seq_len(nrow(cells)), function(i)
    if (!length(dn)) list() else as.list(stats::setNames(cells[i, ], dn)))
  lst <- tryCatch(ui$lstExpr, error = function(e) NULL)
  if (is.null(lst)) return(NULL)
  # DIRECT readers only. An assignment reading an INTERMEDIATE that reads a
  # covariate (cl <- exp(tcl + eta) * wtf) depends on the covariates only
  # through that intermediate, so it adds no direction the span does not
  # already carry -- and including it would double-count the same loading.
  is_asgn <- vapply(lst, function(e) is.call(e) && length(e) == 3L &&
                      (identical(e[[1L]], quote(`<-`)) ||
                       identical(e[[1L]], quote(`=`))), logical(1))
  hit <- which(is_asgn & vapply(lst, function(e)
    is.call(e) && length(e) == 3L &&
      length(intersect(all.vars(e[[3L]]), cn)) > 0L, logical(1)))
  if (!length(hit)) return(NULL)

  # deterministic latent probe, and the covariate values it maps to
  Z <- tryCatch(suppressWarnings(
         stats::qnorm(randtoolbox::sobol(n_probe, dim = pc, seed = 7L))),
       error = function(e) NULL)
  if (is.null(Z) || !is.matrix(Z) || !all(is.finite(Z))) return(NULL)
  Z <- Z %*% Lc
  A <- vapply(seq_len(pc), function(k)
    .admCovQuantile(cd[[cn[k]]], stats::pnorm(Z[, k])), numeric(n_probe))
  colnames(A) <- cn
  if (!all(is.finite(A))) return(NULL)
  # A SEPARATE, much larger probe for the VERIFICATION reference. The loadings
  # come off the small one -- an average derivative needs no precision -- but
  # the reference the design is judged against must be more accurate than the
  # design, and at 128 points its own moments are only good to ~1e-2, which is
  # looser than the tolerance. Costs no solves: these are R evaluations of the
  # parameter assignment.
  Zv <- tryCatch(suppressWarnings(
          stats::qnorm(randtoolbox::sobol(n_ver, dim = pc, seed = 11L))),
        error = function(e) NULL)
  if (is.null(Zv) || !is.matrix(Zv) || !all(is.finite(Zv))) return(NULL)
  Zv <- Zv %*% Lc
  Av <- vapply(seq_len(pc), function(k)
    .admCovQuantile(cd[[cn[k]]], stats::pnorm(Zv[, k])), numeric(n_ver))
  colnames(Av) <- cn
  if (!all(is.finite(Av))) return(NULL)

  # Evaluate the assignments IN ORDER, so an intermediate is defined before the
  # assignment that reads it, and collect the direct readers' values.
  st <- .admShiftStruct(pinfo)
  pr <- list(lst = lst, is_asgn = is_asgn, hit = hit, cn = cn, dn = dn,
             eta_names = pinfo$eta_col_names, cov_fixed = cov_fixed)
  # evaluate the covariate-reading assignments at an ARBITRARY covariate matrix
  # -- the probe uses A, the verification uses the design points Xc
  probe_gen <- function(eta_at, cell, AA, st_use = st) {
    .admCovProbeAt(pr, st_use, eta_at, cell, AA)
  }
  probe    <- function(eta_at, cell) probe_gen(eta_at, cell, A)
  probe_at <- function(AA, cell)     probe_gen(0, cell, AA)
  loadings <- function(P) .admCovLoadings(P, Z, pc)
  B <- loadings(probe(0, cell_list[[1L]]))
  if (is.null(B)) return(NULL)
  pr$routes <- attr(B, "routes")   # settled here, replayed on every refresh
  # THE LOADING MUST NOT DEPEND ON THE RANDOM EFFECT. A covariate-by-eta
  # interaction (cl <- exp(tcl + b * WT * eta.cl)) has a direction that moves
  # with eta, and the probe at eta = 0 would report b = 0 -- a collapse onto
  # the wrong subspace, silently. Re-probe away from zero and require the same
  # loadings.
  # ... or with the STRATUM: a covariate-by-SEX interaction has a direction that
  # differs cell to cell, and one shared design would be wrong in all but one.
  chk <- list()
  if (length(pinfo$eta_col_names)) chk <- c(chk, list(list(0.5, cell_list[[1L]])))
  if (length(cell_list) > 1L)
    chk <- c(chk, lapply(cell_list[-1L], function(cl) list(0, cl)))
  for (cc in chk) {
    B2 <- loadings(probe(cc[[1L]], cc[[2L]]))
    if (is.null(B2) || !isTRUE(all.equal(B, B2, tolerance = 1e-6)))
      return(NULL)
  }

  sv <- tryCatch(svd(B), error = function(e) NULL)
  if (is.null(sv) || !length(sv$d)) return(NULL)
  r <- sum(sv$d > max(sv$d) * 1e-8)
  # r == pc is refused: no rank reduction to make, and with the node search gone
  # there is nothing else on this path to gain. The rotation alone buys nothing
  # there -- with B diagonal, U is a permutation and redistributes nothing.
  if (!is.finite(r) || r < 1L || r >= pc) return(NULL)   # no reduction to make
  U  <- sv$u[, seq_len(r), drop = FALSE]                # pc x r, orthonormal
  nn <- as.integer(n_nodes)
  if (nn^r * max(nrow(cells), 1L) > max_rows) return(NULL)

  # w = t(U) z ~ N(0, t(U) Rc U). INDEPENDENT margins give the identity and the
  # ordinary r-dimensional grid; a CORRELATED set factors it instead, which
  # costs one Cholesky and not a single extra point.
  Sr <- t(U) %*% Rc %*% U
  Lr <- tryCatch(chol(Sr), error = function(e) NULL)
  if (is.null(Lr)) return(NULL)

  # The design for one PER-DIRECTION node count.
  build <- function(nv) {
    gl <- lapply(nv, function(m) .adghNodes1(m))
    Xg <- as.matrix(expand.grid(lapply(gl, function(g) g$x)))
    Wg <- as.numeric(apply(expand.grid(lapply(gl, function(g) g$w)), 1L, prod))
    dimnames(Xg) <- NULL
    Zg <- Xg %*% Lr %*% t(U)                            # preimage z = U w
    # CLAMP before a margin quantile function sees a probability. pnorm()
    # returns exactly 1 from |z| >= 8.3, and the ROTATION reaches further than
    # the one-dimensional node range does -- an r-direction corner sits at
    # sqrt(r) times it, so cov_nodes = 15 already saturates at r = 2 and an
    # unbounded quantile comes back infinite. .admCovNodesFor carries the same
    # guard for the same reason. Here the finite check below caught it, so it
    # was a silent loss of the collapse at raised cov_nodes rather than a wrong
    # answer -- but losing a feature silently is still the wrong outcome.
    Xg2 <- vapply(seq_len(pc), function(k)
      .admCovQuantile(cd[[cn[k]]],
                      pmin(pmax(stats::pnorm(Zg[, k]), .Machine$double.eps),
                           1 - .Machine$double.eps)), numeric(nrow(Zg)))
    if (!is.matrix(Xg2)) Xg2 <- matrix(Xg2, nrow(Zg), pc)
    colnames(Xg2) <- cn
    if (!all(is.finite(Xg2))) return(NULL)
    list(X = Xg2, W = Wg / sum(Wg), z = Zg)
  }

  # Nodes per direction: the CAP, uniform. A search that reduced each direction
  # to where its own moments stopped moving was tried and reverted -- it is a
  # measurement made at the ADMISSION thetas, and a covariate coefficient is
  # estimated, so a direction that looks converged at the starting values is not
  # converged where the optimizer goes. Measured, it shaved one node off a
  # strongly loaded direction for 14% fewer rows and 5-10x the error once b
  # moved. The saving that survives is the rank reduction, which is structural.
  nv <- rep(.admCovDirNodes(nn, pc, r), r)
  if (prod(nv) * max(nrow(cells), 1L) > max_rows) return(NULL)
  dd <- build(nv)
  if (is.null(dd)) return(NULL)
  Xc <- dd$X; Wc <- dd$W; Zc <- dd$z

  # cross the collapsed continuous design with the EXACT discrete enumeration
  nq <- nrow(Xc); nl <- max(nrow(cells), 1L)
  ix <- rep(seq_len(nq), times = nl)
  Xf <- Xc[ix, , drop = FALSE]
  Wf <- Wc[ix] * rep(pcell, each = nq)
  if (length(dn))
    Xf <- cbind(Xf, cells[rep(seq_len(nl), each = nq), , drop = FALSE])
  Xf <- Xf[, nms, drop = FALSE]
  Wf <- Wf / sum(Wf)

  # VERIFY THE DESIGN, NOT THE CERTIFICATE. What the collapse needs is that the
  # reduced design reproduce the LAW of every covariate-reading assignment --
  # that is the property, and every affinity or single-index test is only a
  # proxy for it. Proxies were tried and are not sharp enough: a within-bin
  # spread reports 0.18 for an exact identity link, and a spline residual
  # separates a genuine index from a sum of separate nonlinearities by only a
  # factor of five.
  #
  # So evaluate the assignments at the DESIGN points and compare their weighted
  # moments against a large probe, which is the truth here. Costs no solves --
  # these are R evaluations of the parameter assignment. This subsumes the
  # affine test rather than replacing it: an affine case passes trivially.
  ver <- function(cell, wcell) {
    Pd <- probe_at(Xc, cell)
    if (is.null(Pd)) return(FALSE)
    Pr <- probe_at(Av, cell)
    if (is.null(Pr)) return(FALSE)
    for (k in seq_len(ncol(Pr))) {
      tgt <- Pr[, k]; got <- Pd[, k]
      sc  <- max(stats::sd(tgt), abs(mean(tgt)), .Machine$double.xmin)
      # first two moments, and the reciprocal where it is defined: a PK
      # parameter enters the prediction as both v and 1/v
      mm <- list(function(x) x, function(x) x^2)
      if (all(tgt > 0) && all(got > 0)) mm <- c(mm, list(function(x) 1 / x))
      for (f in mm) {
        a1 <- mean(f(tgt)); a2 <- sum(wcell * f(got))
        if (!is.finite(a1) || !is.finite(a2)) return(FALSE)
        # 5e-3 sits an order of magnitude above the reference own error
        # (~3e-4 at 8192 points) and two orders below the discrepancy a
        # genuine failure produces, which is of order 1.
        if (abs(a2 - a1) / max(abs(a1), sc) > 5e-3) return(FALSE)
      }
    }
    TRUE
  }
  for (i in seq_along(cell_list))
    if (!ver(cell_list[[i]], Wc)) return(NULL)

  # U and Lr are published so a SAMPLER can use the same subspace: admc draws
  # sobol(n, dim = n_eta + p) and QMC error grows with dimension, so drawing w
  # in r dimensions and mapping z = U Lr w gives the identical law of the
  # parameter from a lower-dimensional sequence. Measured on a rank-1 collapse
  # of three covariates, Frobenius error against a large reference: the
  # covariance improves ~2x at every sample size, and the covariance is what
  # log|V| + tr(V^-1 V_obs) leans on.
  list(X = Xf, W = Wf, z = Zc[ix, , drop = FALSE],
       collapsed = TRUE, r = r, p = p, pc = pc, m = ncol(B), n_cell = nl,
       nv = nv,
       U = U, Lr = Lr, cn = cn, dn = dn, cd = cd, nms = nms,
       cells = cells, pcell = pcell,
       # everything .admCovRefresh() needs to redo the rotation at the CURRENT
       # structural thetas. The probe ingredients, not a closure: a closure
       # captures its whole defining environment and has to survive being
       # stored on the study and shipped to a daemon.
       pr = pr, st0 = st, Zp = Z, Ap = A, Rc = Rc, cell_list = cell_list)
}

# =============================================================================
# JOINT COLLAPSE -- etas and covariates in ONE latent space
# =============================================================================
#
# The design crosses two things that are both latent normal directions: the eta
# grid (n_nodes^n_eta) and the covariate design. They are crossed as if
# independent, and frequently they are not -- an eta and a covariate index reach
# the model through the SAME sum:
#
#   cl <- exp(tcl + eta.cl) * (W1/70)^b1 * (W2/70)^b2 * (W3/70)^b3
#
# In the joint latent space that is ONE direction, not four. And the bound is
# structural: with xi = (eta_std, z) standard normal in both blocks,
#
#   log theta = a + B' xi        B's eta rows scale with L = chol(Omega)
#
# so rank(B) <= the number of PARAMETERS the latents reach, however many etas
# and covariates there are. A two-parameter model never needs more than a
# two-dimensional design.
#
# Measured on 2 etas + 3 covariates on the eta'd parameter, against a 15^2 x
# 15^3 product reference and against the shipping design, at six parameter
# points: rank 2 of 5, 324 design points against 525, and a covariance error of
# 1e-12 to 3e-09 against the shipping design's 6e-07 to 1.7e-05.
#
# THE ROTATION MOVES WITH THE FIT, for two reasons rather than the covariate
# collapse's one: the covariate coefficients are estimated, AND Omega enters B
# through L. So it is re-aimed on every objective call, and nothing about the
# direction is cached. What is fixed at admission is structural only -- the
# rank, the node count, the certificate.
.admJointCollapse <- function(ui, pinfo, cov_dist, n_nodes, s, out_var,
                              n_probe = 512L, max_rows = 20000L,
                              n_ver = 8192L, cov_fixed = NULL) {
  ne <- pinfo$n_eta
  if (!is.finite(ne) || ne < 1L) return(NULL)
  cd  <- .admCovDistCanon(cov_dist)
  nms <- .admCovSpecNames(cd)
  dsc <- vapply(nms, function(n) !is.null(cd[[n]][["values"]]), logical(1))
  cn  <- nms[!dsc]                    # CONTINUOUS: what rotates
  dn  <- nms[dsc]                     # DISCRETE: enumerated, as strata
  pc  <- length(cn)
  # A discrete covariate is a stratum, not a direction -- it has no latent
  # normal to rotate into. It is crossed with the continuous design exactly as
  # .admCovCollapse crosses it, at its declared levels and probabilities, so
  # the two constructions stack instead of one disqualifying the other. Sex,
  # genotype and formulation are about as common as covariates get, and they
  # used to turn the whole joint path off.
  if (pc < 1L) return(NULL)
  R <- cd[["latentR"]]
  if (is.function(cd[["joint"]]) && is.null(R)) return(NULL)
  ic <- match(cn, nms); id <- match(dn, nms)
  # A discrete covariate DEPENDENT on a continuous one makes each level a
  # TRUNCATION of the latent normal rather than a point, so the continuous
  # conditional law differs cell to cell and one shared rotation would be the
  # wrong rotation in every cell but one. The same refusal .admCovCollapse and
  # .admCovTaylorDesign carry.
  if (length(dn) && !is.null(R) &&
      any(abs(R[id, ic, drop = FALSE]) > 0)) return(NULL)
  Rc <- if (is.null(R)) diag(1, pc) else R[ic, ic, drop = FALSE]
  Lc <- tryCatch(chol(Rc), error = function(e) NULL)
  if (is.null(Lc)) return(NULL)
  # exact product enumeration of the discrete levels, with their probabilities
  if (length(dn)) {
    lev <- lapply(dn, function(n) as.numeric(cd[[n]][["values"]]))
    prb <- lapply(dn, function(n) { sp <- cd[[n]]
      pv <- sp[["probs"]] %||% rep(1 / length(sp[["values"]]),
                                   length(sp[["values"]]))
      pv / sum(pv) })
    cells <- as.matrix(expand.grid(lev, KEEP.OUT.ATTRS = FALSE))
    colnames(cells) <- dn
    pcell <- apply(as.matrix(expand.grid(prb, KEEP.OUT.ATTRS = FALSE)), 1L, prod)
  } else {
    cells <- matrix(numeric(0), 1L, 0L); pcell <- 1
  }
  cell_list <- lapply(seq_len(nrow(cells)), function(i)
    if (!length(dn)) list() else as.list(stats::setNames(cells[i, ], dn)))
  lst <- tryCatch(ui$lstExpr, error = function(e) NULL)
  if (is.null(lst)) return(NULL)
  is_asgn <- vapply(lst, function(e) is.call(e) && length(e) == 3L &&
                      (identical(e[[1L]], quote(`<-`)) ||
                       identical(e[[1L]], quote(`=`))), logical(1))
  # DIRECT readers of a covariate OR an eta -- the joint space spans both, so a
  # parameter carrying only an eta is as much a direction as one carrying only
  # covariates.
  lat <- c(cn, dn, pinfo$eta_col_names)
  hit <- which(is_asgn & vapply(lst, function(e)
    is.call(e) && length(e) == 3L &&
      length(intersect(all.vars(e[[3L]]), lat)) > 0L, logical(1)))
  if (!length(hit)) return(NULL)
  pr <- list(lst = lst, is_asgn = is_asgn, hit = hit, cn = cn,
             dn = dn, eta_names = pinfo$eta_col_names,
             cov_fixed = cov_fixed)
  nl <- ne + pc
  mkXi <- function(n, seed) {
    Z <- tryCatch(suppressWarnings(
           stats::qnorm(randtoolbox::sobol(n, dim = nl, seed = seed))),
         error = function(e) NULL)
    if (is.null(Z) || !is.matrix(Z) || !all(is.finite(Z))) return(NULL)
    Z
  }
  Xi <- mkXi(n_probe, 13L); Xv <- mkXi(n_ver, 17L)
  if (is.null(Xi) || is.null(Xv)) return(NULL)
  # r, m and routes are settled by .admJointAdmit(), but they are declared HERE,
  # holding NULL, and that is load-bearing rather than tidiness.
  #
  # `$` PARTIAL-MATCHES on lists. While `m` was absent, jc$m resolved to
  # jc$max_rows -- 20000 -- so the row cap compared m^r against itself and
  # rejected every design, silently and at every parameter point. Declaring the
  # field means `$` always finds an exact match and can never fall through to a
  # prefix. `list(m = NULL)` does create the element; it is not dropped.
  #
  # The [[ ]] reads downstream stay as a second line of defence, and
  # test-covariate.R runs these paths under warnPartialMatchDollar so a field
  # added later cannot quietly reintroduce it.
  list(pr = pr, cn = cn, cd = cd, nms = nms, Rc = Rc, Lc = Lc, ne = ne, pc = pc,
       nl = nl, Xi = Xi, Xv = Xv, out_var = out_var,
       n_nodes = as.integer(n_nodes), max_rows = max_rows, joint = TRUE,
       dn = dn, nms = nms, cells = cells, pcell = pcell,
       cell_list = cell_list, n_cell = max(nrow(cells), 1L),
       r = NULL, m = NULL, routes = NULL)
}

# The covariate values a latent block maps to, with the same clamp the collapse
# and the sampler use: pnorm() saturates at exactly 1 from |z| >= 8.3, and a
# rotation reaches further than the one-dimensional node range.
.admJointCov <- function(jc, Zc) {
  X <- vapply(seq_len(jc$pc), function(k)
    .admCovQuantile(jc$cd[[jc$cn[k]]],
                    pmin(pmax(stats::pnorm(Zc[, k]), .Machine$double.eps),
                         1 - .Machine$double.eps)), numeric(nrow(Zc)))
  if (!is.matrix(X)) X <- matrix(X, nrow(Zc), jc$pc)
  colnames(X) <- jc$cn
  X
}

# The joint loading matrix at a given (struct, L). Etas enter scaled by L, so
# this moves with Omega as well as with the structural thetas.
.admJointB <- function(jc, st, L, Xi, routes = NULL, cell = NULL) {
  Zc <- Xi[, jc$ne + seq_len(jc$pc), drop = FALSE] %*% jc$Lc
  Et <- Xi[, seq_len(jc$ne), drop = FALSE] %*% t(L)
  AA <- .admJointCov(jc, Zc)
  P  <- .admCovProbeAt(jc$pr, st, Et, cell %||% jc$cell_list[[1L]], AA)
  if (is.null(P)) return(NULL)
  .admCovLoadings(P, Xi, jc$nl, routes)
}

# Re-aim the joint design at the CURRENT parameters, and build it.
.admJointDesign <- function(jc, st, L) {
  B <- .admJointB(jc, st, L, jc$Xi, jc[["routes"]])
  if (is.null(B)) return(NULL)
  sv <- tryCatch(svd(B), error = function(e) NULL)
  if (is.null(sv) || !length(sv$d) || max(sv$d) <= 0) return(NULL)
  # [[ ]] not $: `$` PARTIAL-MATCHES on lists, so jc$m silently resolved to
  # jc$max_rows (20000) and the row cap then rejected every design. Both of
  # these are deliberately absent until admission fixes them, which is exactly
  # the case partial matching turns into a wrong answer instead of a NULL.
  r <- jc[["r"]] %||% sum(sv$d > max(sv$d) * 1e-8)
  if (!is.finite(r) || r < 1L || r > jc$nl) return(NULL)
  U <- sv$u[, seq_len(r), drop = FALSE]
  # the cap lesson from .admCovDirNodes, over the joint space: a direction
  # absorbs (n_eta + pc)/r axes, so it needs that much more resolution than one
  m <- jc[["m"]] %||% .admCovDirNodes(jc$n_nodes, jc$nl, r)
  nl_c <- jc$n_cell %||% 1L
  if (m^r * nl_c > jc$max_rows) return(NULL)
  g  <- .adghNodeGrid(m, r)
  Xz <- g$X %*% t(U)                                   # preimage xi = U w
  Xe <- Xz[, seq_len(jc$ne), drop = FALSE]
  eta <- Xe %*% t(L)
  colnames(eta) <- jc$pr$eta_names
  Zc <- Xz[, jc$ne + seq_len(jc$pc), drop = FALSE] %*% jc$Lc
  X  <- .admJointCov(jc, Zc)
  if (!all(is.finite(X)) || !all(is.finite(eta))) return(NULL)
  Wg <- g$W / sum(g$W)
  # CROSS the rotated continuous design with the exact discrete enumeration.
  # Same stride .admCovCollapse uses -- the continuous block cycles fastest --
  # so eta, X and cov_rows stay aligned row for row with the weights.
  if (length(jc$dn)) {
    nq <- nrow(X); ix <- rep(seq_len(nq), times = nl_c)
    ic <- rep(seq_len(nl_c), each = nq)
    X   <- cbind(X[ix, , drop = FALSE],
                 jc$cells[ic, , drop = FALSE])[, jc$nms, drop = FALSE]
    eta <- eta[ix, , drop = FALSE]; colnames(eta) <- jc$pr$eta_names
    Xe  <- Xe[ix, , drop = FALSE]
    Wg  <- Wg[ix] * jc$pcell[ic]
    Wg  <- Wg / sum(Wg)
  }
  # Xe is handed back as the node matrix the omega chain rule differentiates.
  # eta = Xe L' has exactly the form the ordinary grid has (eta = X L'), so
  # d(eta[q,])/d(L_ij) = Xe[q,j] * e_i and .adghGrad needs no new branch. What
  # it does NOT carry is the rotation's own dependence on Omega -- U moves with
  # L too. That term is the quadrature re-choosing itself: the integral is the
  # same for any U spanning B's column space, so it vanishes to the accuracy the
  # design is verified to. FD-checked rather than argued.
  list(eta = eta, X = Xe, cov_rows = X, W = Wg, r = r, m = m,
       U = U, d = sv$d, routes = attr(B, "routes"))
}

# Settle everything STRUCTURAL about the joint design, once, and verify it.
#
# Rank, node count and the per-column loading route are fixed here and replayed
# on every refresh. None of them may be re-derived per call: rank and node count
# would change the NUMBER of design points mid-fit and step the objective, and a
# route is a threshold decision that a borderline column could flip. The
# DIRECTION is the only thing that moves, and it has to.
#
# Verification is the same instrument the covariate collapse uses: score the
# reader assignments at the design points against a much larger probe, on the
# first two moments and the reciprocal. Costs no solves.
.admJointAdmit <- function(jc, st, L, tol = 5e-3) {
  if (is.null(jc)) return(NULL)
  jd <- .admJointDesign(jc, st, L)
  if (is.null(jd)) return(NULL)
  jc$r <- jd$r; jc$m <- jd$m; jc$routes <- jd$routes
  cl_list <- jc$cell_list %||% list(list())
  # THE ROTATION MUST NOT DIFFER BETWEEN STRATA. A covariate-by-stratum
  # interaction -- (WT/70)^(b + c*SEX) -- has a direction that changes cell to
  # cell, so a single shared design would be the right design in one cell and
  # the wrong one in all the others. Re-probe in each and require the same
  # loadings, as .admCovCollapse does.
  if (length(cl_list) > 1L) {
    B0 <- .admJointB(jc, st, L, jc$Xi, jc$routes, cl_list[[1L]])
    if (is.null(B0)) return(NULL)
    for (cc in cl_list[-1L]) {
      Bk <- .admJointB(jc, st, L, jc$Xi, jc$routes, cc)
      if (is.null(Bk) || !isTRUE(all.equal(B0, Bk, tolerance = 1e-6,
                                           check.attributes = FALSE)))
        return(NULL)
    }
  }
  # the truth to score against: a large probe in the SAME latent space
  Zv <- jc$Xv[, jc$ne + seq_len(jc$pc), drop = FALSE] %*% jc$Lc
  Ev <- jc$Xv[, seq_len(jc$ne), drop = FALSE] %*% t(L)
  Av <- .admJointCov(jc, Zv)
  nq <- nrow(jd$eta) %/% length(cl_list)
  for (ci in seq_along(cl_list)) {
    cell <- cl_list[[ci]]
    Pv <- .admCovProbeAt(jc$pr, st, Ev, cell, Av)
    if (is.null(Pv)) return(NULL)
    # this cell's slice of the design, and its weights renormalised within it
    ix <- (ci - 1L) * nq + seq_len(nq)
    Wc <- jd$W[ix]; sw <- sum(Wc)
    if (!is.finite(sw) || sw <= 0) return(NULL)
    Wc <- Wc / sw
    Pd <- .admCovProbeAt(jc$pr, st, jd$eta[ix, , drop = FALSE], cell,
                         jd$cov_rows[ix, , drop = FALSE])
    if (is.null(Pd) || ncol(Pd) != ncol(Pv)) return(NULL)
    for (k in seq_len(ncol(Pv))) {
      tgt <- Pv[, k]; got <- Pd[, k]
      sc  <- max(stats::sd(tgt), abs(mean(tgt)), .Machine$double.xmin)
      mm  <- list(function(x) x, function(x) x^2)
      if (all(tgt > 0) && all(got > 0)) mm <- c(mm, list(function(x) 1 / x))
      for (f in mm) {
        a1 <- mean(f(tgt)); a2 <- sum(Wc * f(got))
        if (!is.finite(a1) || !is.finite(a2)) return(NULL)
        if (abs(a2 - a1) / max(abs(a1), sc) > tol) return(NULL)
      }
    }
  }
  jc
}
