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
# that is where the observed V carries it. Three paths compute those moments --
# collapse, u-quantile and the general per-row path -- chosen per study by
# .admCheckCovariates(); see each section below.
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
    # tryCatch into an Inf objective everywhere. It also made .admCovDelta's
    # probe -- which passes cov_s = NULL with a probe matrix -- unable to supply
    # a second covariate, which is what killed the u-quantile path for any model
    # with more than one covariate.
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

# =============================================================================
# Covariate collapse (exact marginal moments, no nodes)
# =============================================================================

# Map covariate -> (coefficient theta, eta) from rxode2's OWN mu-reference
# metadata, so nothing has to be declared by hand. muRefCovariateDataFrame gives
# (theta, covariate, covariateParameter); muRefDataFrame gives (theta, eta).
# Joining on `theta` says which eta's mu-referenced argument each covariate
# enters, and with which coefficient.
.admCovMap <- function(ui) {
  cd <- tryCatch(ui$muRefCovariateDataFrame, error = function(e) NULL)
  md <- tryCatch(ui$muRefDataFrame,          error = function(e) NULL)
  if (is.null(cd) || is.null(md)) return(NULL)
  if (!NROW(cd) || !NROW(md)) return(NULL)
  m <- merge(cd, md, by = "theta")
  if (!NROW(m)) return(NULL)
  data.frame(covariate = as.character(m$covariate),
             coef      = as.character(m$covariateParameter),
             eta       = as.character(m$eta),
             stringsAsFactors = FALSE)
}

# Sigma_a for the named covariates, in the order given.
.admCovDistSigma <- function(cov_dist, nms) {
  d   <- length(nms)
  # [["sd"]] not $sd: `$` PARTIAL-MATCHES, so a lognormal spec's `sdlog` was
  # returned as if it were a natural-scale `sd` -- one routing change away from
  # folding a log-scale spread into Omega.
  sds <- vapply(nms, function(n) {
    v <- cov_dist[[n]][["sd"]]
    if (is.null(v)) NA_real_ else as.numeric(v)
  }, numeric(1))
  if (anyNA(sds))
    stop(".admCovDistSigma: covariate(s) ",
         paste(sQuote(nms[is.na(sds)]), collapse = ", "),
         " have no natural-scale `sd`; the collapse is defined for normal ",
         "covariates only.", call. = FALSE)
  if (!is.null(cov_dist$Sigma)) return(cov_dist$Sigma[nms, nms, drop = FALSE])
  rho <- cov_dist$rho %||% 0
  S <- diag(sds^2, d)
  if (d >= 2L && abs(rho) > 1e-12)
    for (i in seq_len(d - 1L)) for (j in seq(i + 1L, d))
      S[i, j] <- S[j, i] <- rho * sds[i] * sds[j]
  S
}

# Effective Cholesky for a study whose subjects SPAN a covariate distribution.
#
# With mu-referencing the covariate and the random effect enter the SAME
# argument:  param <- g(theta + theta_cov*a + eta).  So for a ~ N(mu_a, Sigma_a),
#
#     theta + theta_cov*a + eta  ==  theta + theta_cov*mu_a + s,
#     s ~ N(0, Omega + J Sigma_a J'),      J[eta_k, k] = theta_cov_k
#
# EXACTLY. Solving at the covariate MEAN and inflating Omega therefore reproduces
# the marginal moments with one solve set -- no quadrature nodes, no importance
# weights, and the same cost as a fit with no covariate at all.
#
# Exact only for a LINEAR covariate effect on the mu-referenced scale and a
# normal covariate. Returns NULL whenever that cannot be established, so every
# caller falls back to the plain Cholesky rather than silently approximating.
.admCovInflateL <- function(pars, pinfo, s) {
  cd <- s$cov_dist
  if (is.null(cd) || pinfo$n_eta == 0L || is.null(pars$L)) return(NULL)
  cmap <- pinfo$cov_map
  if (is.null(cmap) || !NROW(cmap)) return(NULL)
  eta_nms <- pinfo$eta_col_names
  has_sd  <- vapply(names(cd), function(n)
    is.list(cd[[n]]) && !is.null(cd[[n]]$sd), logical(1))
  keep <- cmap$covariate %in% names(cd)[has_sd] & cmap$eta %in% eta_nms
  use  <- cmap[keep, , drop = FALSE]
  if (!NROW(use)) return(NULL)
  J <- matrix(0, pinfo$n_eta, NROW(use))
  for (k in seq_len(NROW(use))) {
    # pars$struct holds ESTIMATED thetas only. `[[` with an unmatched name on a
    # LIST errors rather than returning NULL, so the guard below was dead and a
    # fix()ed covariate coefficient killed the fit at the first objective
    # evaluation (.admNLL calls .admStudyL with no tryCatch). Look the value up
    # by name, falling back to the model's fixed value.
    cf <- if (use$coef[k] %in% names(pars$struct))
      pars$struct[[use$coef[k]]] else pinfo$fixed_theta[[use$coef[k]]]
    if (is.null(cf) || !is.finite(cf)) return(NULL)
    J[match(use$eta[k], eta_nms), k] <- cf
  }
  Sig <- .admCovDistSigma(cd, use$covariate)
  Om  <- tcrossprod(pars$L) + J %*% Sig %*% t(J)
  tryCatch(t(chol(Om)), error = function(e) NULL)
}

# The Cholesky a given study should be simulated with: inflated when the study
# declares a covariate DISTRIBUTION, the plain one otherwise.
#
# Returns NULL -- never a silent fall back to pars$L -- when a study DOES declare
# cov_dist and the inflation cannot be built. Falling back there would solve at
# the covariate mean and understate the variance, which is a plausible-looking
# fit of the wrong model. Callers treat NULL as "cannot evaluate" (Inf / NA).
# .admCheckCovariates() rejects the reachable causes up front, so in practice
# this only fires on a pathological Omega.
.admStudyL <- function(pars, pinfo, s) {
  if (!identical(s$.adm_cov_path, "collapse")) return(pars$L)
  .admCovInflateL(pars, pinfo, s)
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
# study unchanged on the collapse path (where the covariate is held at its mean
# and its variance lives in Omega) and when no distribution is declared.
.admStudyCovRows <- function(s, pinfo, n_row) {
  if (!identical(s$.adm_cov_path, "rows")) return(s)
  s$cov_rows <- .admCovRowsFor(s$cov_dist, n_row, pinfo$n_eta)
  s
}

# Up-front validation for every study that declares a covariate distribution.
#
# Scope, deliberately narrow and enforced rather than assumed: the collapse is
# EXACT only for a covariate that enters as a bare `theta * COV` product inside a
# mu-referenced expression, appears nowhere else in the model, and is normally
# distributed. rxode2 records the first of those in muRefCovariateDataFrame --
# and records it ONLY for that bare product. Measured on 6 covariate models:
#
#   cl <- exp(tcl + tcov*WT + eta.cl)              -> 1 row   (collapsible)
#   cl <- exp(tcl + tcov*log(WT/70) + eta.cl)      -> 0 rows
#   cl <- exp(tcl + eta.cl)*(WT/70)^tcov           -> 0 rows
#   cl <- exp(tcl + eta.cl)*exp(tcov*WT)           -> 0 rows
#   cl <- exp(tcl + eta.cl)*(1 + emax*WT/(ec50+WT))-> 0 rows
#   cl <- exp(tcl + tcov*SEX + eta.cl)             -> 1 row   (but SEX is not normal)
#
# The allometric form is mathematically collapsible (it is linear in log(WT)) and
# is still refused here, because rxode2's metadata cannot distinguish it from the
# Emax form. Refusing is the safe direction: the alternative is a fit that runs,
# converges, and reports an omega that has quietly eaten the covariate spread.
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

.admCheckCovariates <- function(.ui, pinfo, studies, grad, est = NULL) {
  has <- vapply(studies, function(s) !is.null(s$cov_dist), logical(1))
  if (!any(has)) return(studies)
  bad <- function(...) stop("admixr2: ", ..., call. = FALSE)

  # Expand the friendly grammar ONCE, before anything routes on it: `cor`
  # becomes a `joint` sampler, and `joint` is what adgh refuses.
  for (nm in names(studies)[has])
    studies[[nm]]$cov_dist <- .admCovDistCanon(studies[[nm]]$cov_dist)

  cmap <- pinfo$cov_map
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
    ok_collapse <- pinfo$n_eta > 0L
    ok_uq       <- pinfo$n_eta > 0L
    uq          <- list()

    # `rho` and `Sigma` are metadata siblings of the covariate specs, not
    # covariates. Looping over them made a correlated spec fail with
    # "declares cov_dist for 'rho', which the model never reads" -- a message
    # pointing at the wrong thing, and the reason .admCovDistSigma()'s
    # correlated branch was unreachable.
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

      once <- .admNameOccurrence(.ui, cv)[[cv]] == 1L
      pe   <- if (once) .admCovParamEta(.ui, cv, pinfo$eta_col_names) else NULL

      # COLLAPSE: exact closed form. Needs a normal covariate entering as a
      # single bare `theta * COV` term.
      # The collapse is exact only when the model is solved AT the covariate
      # mean: Omega + J Sigma_a J' reproduces the marginal moments about mu_a and
      # nothing else. `cov` is filled from the distribution above when absent,
      # but a user-supplied value passes straight through -- and a reference
      # weight supplied alongside a distribution then biases the whole study
      # silently (measured: 13-19% on the mean, 10-14% on the covariance for a
      # value 0.2 away from the mean). The "rows" path ignores `cov` entirely, so
      # the two paths would also disagree on identical input.
      .mu_a <- .admCovMeanOf(sp)
      .at_mean <- !is.null(.mu_a) && is.finite(.mu_a) &&
                  isTRUE(all.equal(as.numeric(studies[[nm]][["cov"]][[cv]]),
                                   as.numeric(.mu_a), tolerance = 1e-8))
      if (!normal || is.null(cmap) || !cv %in% cmap$covariate || !once ||
          !.at_mean)
        ok_collapse <- FALSE
      # u-QUANTILE: the covariate's WHOLE effect has to fit in one eta column, as
      # u = Delta(a) + eta. That substitution is only valid under conditions the
      # path cannot check for itself, so it is granted only when EVERY one of
      # them is positively established, and refused by default:
      #
      #  * exactly once in the model, sharing a parameter assignment with an eta,
      #    with a reference value to measure Delta against (the original test);
      #  * a SMOOTH, UNIMODAL covariate distribution. .admCovUQuantile solves the
      #    mixture CDF by clamped Newton with no bracketing: on a discrete or
      #    multi-modal Delta the iterate oscillates between the clamps, exits at
      #    maxit WITHOUT warning, and the garbage quantiles go straight into the
      #    eta column. Measured on a two-point covariate: 13-20% on the mean and
      #    3-5x on the covariance, and SMALLER omega makes it worse.
      #  * a Gauss-Hermite-representable spec. A user `quantile` is integrated by
      #    32 equal-weight midpoints, which understates the spread of Delta
      #    (measured: sd 0.44117 against a truth of 0.45000) and biases the
      #    coefficient upward by ~7%.
      #  * ONE covariate per eta. .admCovUQEta writes u into eta_mat[, j] per
      #    covariate, so two covariates resolving to the same j silently destroy
      #    each other -- the result is bit-for-bit the second covariate alone.
      #
      # Delta is also measured as a LOG shift (.admCovDelta), so the covariate
      # must act multiplicatively on a log-mu-referenced parameter. An additive
      # effect, an identity mu-reference or an expit link all break it (16-21% on
      # the covariance). That one is checked by .admCovDeltaValid() below.
      #
      # "rows" costs more and assumes nothing; it was exact in every
      # configuration tested. When in doubt this must fall back to it.
      # ... and the honest answer is that admixr2 cannot establish them. Each of
      # the four conditions above is a property of the MODEL TEXT, and the
      # syntactic checks that would decide them are exactly as fragile as the
      # assumption they are meant to protect: a check on `exp(` and the covariate
      # name accepts `cl <- exp(tcl + eta.cl) + tcov*WT`, where the covariate is
      # outside the exp and the shift does not hold. So the path is NOT ROUTED TO.
      #
      # It is an optimisation, not a capability: "rows" computes the same moments
      # without any of these preconditions, and was exact in every configuration
      # an independent audit could construct (residual 2e-5..1.6e-3, i.e. Monte
      # Carlo noise). uq is also only ever reachable with grad = "none", which no
      # estimator defaults to. So what it bought was speed in one narrow case, at
      # the cost of four measured silent-wrong-answer paths:
      #   - a discrete or multi-modal covariate breaks .admCovUQuantile's
      #     unbracketed Newton solve: 13-20% on the mean, 3-5x on the covariance,
      #     WORSE for smaller omega, with no warning;
      #   - a non-log link or an additive effect breaks .admCovDelta's log-shift:
      #     16-21% on the covariance;
      #   - two covariates on one eta silently overwrite each other, leaving the
      #     result bit-for-bit identical to the second covariate alone;
      #   - a user `quantile` is integrated by 32 equal-weight midpoints, which
      #     understates the spread and biases the coefficient ~7% high.
      #
      # .admCovUQuantile/.admCovDelta/.admCovUQEta and their tests are kept: the
      # quantile solve is correct and well-tested where Delta is smooth and
      # unimodal, and the u-distribution framing is the basis for handling
      # dependent covariates. Reinstating the route needs a POSITIVE check of the
      # shift property, and the way to do that is numerically -- verify
      # f(a, eta) == f(a_ref, eta + Delta(a)) at a couple of (a, eta) points --
      # not by reading the model text.
      ok_uq <- FALSE
    }

    # A GRADIENT is only carried through the general path. On "rows" the
    # covariate is data -- a per-row params column -- so the existing sensitivity
    # directions differentiate exactly the function the NLL evaluates, with no
    # new chain rule. "collapse" moves Omega itself (Omega + J Sigma_a J', with J
    # carrying the covariate coefficient) and "uq" replaces an eta column with
    # u = F_u^-1(Phi(z)); neither derivative exists yet.
    #
    # So the gradient mode is part of what makes a path VALID, not a reason to
    # refuse the fit. Erroring here would make a covariate model fail out of the
    # box, since every estimator defaults to a gradient.
    if (!identical(grad, "none")) { ok_collapse <- FALSE; ok_uq <- FALSE }

    # A JOINT sampler describes a DEPENDENT covariate distribution, which the
    # COLLAPSE still cannot represent: it assembles Sigma_a from per-covariate
    # variances plus an optional scalar rho, so it is refused and the general
    # path takes over.
    #
    # adgh no longer refuses one. Its grid is a tensor product over the
    # SAMPLER'S UNIFORMS, not over the covariates, and a copula's uniforms are
    # independent by construction -- that is what a copula is -- so the product
    # rule is exact there whatever the dependence. See .admCovGrid(). The
    # earlier refusal read the product grid as "assumes independence", which is
    # true of a grid over covariate margins and false of a grid over u.
    if (is.function(cd[["joint"]])) ok_collapse <- FALSE


    # Most efficient VALID path wins. "rows" assumes nothing at all: every
    # simulated subject carries its own covariate value, so rxode2 evaluates the
    # whole model -- covariate on several parameters, on a parameter with no eta,
    # or interacting with one. It is the only path with no structural
    # precondition, hence the fallback.
    studies[[nm]]$.adm_cov_path <- if (ok_collapse) "collapse"
                                   else if (ok_uq)  "uq" else "rows"
    studies[[nm]]$.adm_cov_collapse <- ok_collapse
    studies[[nm]]$.adm_cov_uq <- if (ok_uq) uq else NULL

    # cov_integration = "shift": if the covariates act only as a shift of a
    # mu-referenced argument, the whole covariate dimension leaves the solve.
    # ADMITTED BY NUMERICAL VERIFICATION, not by reading the model text -- the
    # earlier `uq` route inferred this property from syntax and had four
    # measured silent-wrong-answer modes as a result.
    if (identical(studies[[nm]]$.adm_cov_path, "rows") &&
        identical(pinfo$cov_integration %||% "quadrature", "shift")) {
      .cn  <- .admCovSpecNames(cd)
      .sp  <- .admShiftSpec(.ui, .cn, pinfo$eta_col_names)
      .ar  <- if (is.null(.sp)) NULL else .admShiftRef(cd, .cn)
      .why <- NULL
      if (is.null(.sp)) .why <- paste0(
        "the covariates do not all enter one model assignment together with ",
        "exactly one random effect, so no single shifted column can carry them")
      else if (is.null(.ar)) .why <- "a covariate has no finite mean to shift against"
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
      if (!is.null(.why))
        bad("study '", nm, "': cov_integration = \"shift\" is not applicable -- ",
            .why, ". Use cov_integration = \"quadrature\".")
      # The shift replaces ONE eta column. If that eta is correlated with
      # another, its marginal is no longer sqrt(Omega_jj) and the column cannot
      # be substituted independently, so refuse rather than quietly use the
      # diagonal.
      .j <- match(.sp$eta, pinfo$eta_col_names)
      .off <- which(!pinfo$chol_diag &
                    (pinfo$chol_i %in% .j | pinfo$chol_j %in% .j))
      if (length(.off))
        bad("study '", nm, "': cov_integration = \"shift\" is not applicable -- ",
            "the random effect(s) ", paste(sQuote(.sp$eta), collapse = ", "),
            " that would carry the covariate are CORRELATED (an off-diagonal ",
            "Omega element is estimated), ",
            "element is estimated), so its marginal is not sqrt(Omega_jj) and ",
            "the column cannot be substituted on its own. Use ",
            "cov_integration = \"quadrature\".")
      .g <- .admCovGrid(cd, pinfo$cov_nodes %||% 7L)
      studies[[nm]]$.adm_cov_shift <- list(
        spec = .sp, aref = .ar, cov_names = .cn,
        eta_idx = match(.sp$eta, pinfo$eta_col_names),
        m = length(.sp$eta), X = .g$X, W = .g$W)
      studies[[nm]]$.adm_cov_path <- "shift"
    }

    # cov_integration = "taylor" replaces the product grid on the "rows" path
    # (and only there -- the collapse is exact and cheaper still). Its refusals
    # are properties of `cov_dist` alone, so build the design once HERE, where
    # the message can name the study, rather than letting the first objective
    # evaluation error out of the middle of a fit.
    if (identical(studies[[nm]]$.adm_cov_path, "rows") &&
        identical(pinfo$cov_integration %||% "quadrature", "taylor")) {
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
# The collapse above is exact but narrow: it needs the covariate to enter as a
# bare `theta * COV` term so its variance can be folded into Omega. Everything
# else -- (WT/70)^theta, exp(theta*WT), allometric-in-log, Emax, if/else,
# categorical -- goes through here instead.
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
  # `cor` is the user-facing one; `rho`/`Sigma` are the collapse's older
  # normal-scale form. Until this was unified, a spec written with `rho` built
  # its Gaussian copula for the collapse and a DIAGONAL grid for everything
  # else -- the correlation silently present in one path and absent in the
  # other. `rho`/`Sigma` are left in place afterwards, because
  # .admCovDistSigma() still reads them directly.
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
    # scalar rho applies to EVERY pair -- .admCovDistSigma's own convention
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
  # `rho`/`Sigma` are metadata for the collapse; `joint` is the sampler itself.
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
# it IS the omega'^2 = omega^2 + beta^2 Var(a) inflation the collapse computes
# in closed form, appearing here as a rank-p addition to V. Dropping it leaves
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
       cont = cn, disc = dn, n_cell = L, p_cell = pcell, n_cpt = n_cpt,
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
# coordinates, which are independent by construction whatever the dependence.
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
    gq  <- .adghNodes1(n_nodes)
    ixg <- as.matrix(expand.grid(rep(list(seq_len(n_nodes)), length(iS)),
                                 KEEP.OUT.ATTRS = FALSE))
    K   <- nrow(ixg)
    wk  <- apply(matrix(gq$w[ixg], K, length(iS)), 1L, prod); wk <- wk / sum(wk)
    zS  <- matrix(gq$x[ixg], K, length(iS))
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
      return(list(X = X, W = rep(1 / nrow(X), nrow(X))))
    }
    g  <- .adghNodes1(n_nodes)
    ix <- as.matrix(expand.grid(rep(list(seq_len(n_nodes)), d),
                                KEEP.OUT.ATTRS = FALSE))
    z  <- matrix(g$x[ix], nrow(ix), d)
    W  <- apply(matrix(g$w[ix], nrow(ix), d), 1L, prod)
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
    return(list(X = X, W = W / sum(W)))
  }
  one <- lapply(nms, function(n) .admCovNodesFor(cov_dist[[n]], n_nodes))
  X   <- as.matrix(expand.grid(lapply(one, `[[`, "x"), KEEP.OUT.ATTRS = FALSE))
  W   <- Reduce(`*`, lapply(seq_along(one), function(k)
           rep(rep(one[[k]]$w, each = prod(vapply(one[seq_len(k - 1L)],
                 function(o) length(o$x), integer(1)))),
               length.out = nrow(X))))
  colnames(X) <- nms
  list(X = X, W = W / sum(W))
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
.admShiftNodes <- function(Delta, W, om, n_u, tol = 1e-12, maxit = 30L) {
  g  <- .adghNodes1(n_u)
  tg <- stats::pnorm(g$x)
  mD <- sum(W * Delta); vD <- max(sum(W * Delta^2) - mD^2, 0)
  u  <- mD + sqrt(vD + om^2) * g$x
  for (it in seq_len(maxit)) {
    Z  <- outer(u, Delta, "-") / om
    Fu <- as.numeric(stats::pnorm(Z) %*% W)
    fu <- as.numeric(stats::dnorm(Z) %*% W) / om
    st <- (Fu - tg) / pmax(fu, 1e-300)
    # safeguarded: an unbounded Newton step on a mixture CDF can leap past the
    # support of a neighbouring component and stall in a flat tail
    u  <- u - sign(st) * pmin(abs(st), 2 * om)
    if (max(abs(st)) < tol) break
  }
  list(u = u, w = g$w)
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
.admShiftDu <- function(spec, struct, X, aref, Delta, W, om, u, h = 1e-6) {
  # Scalar shift only. With m > 1 the later coordinates' nodes also move through
  # the POSTERIOR weights of the Rosenblatt recursion, which is a second chain
  # this does not carry -- .adghGrad finite-differences that case instead.
  Delta <- as.numeric(as.matrix(Delta)[, 1L])
  Z   <- outer(u, Delta, "-") / om
  Ph  <- stats::dnorm(Z)
  den <- pmax(as.numeric(Ph %*% W), 1e-300)
  du_dom <- as.numeric((Ph * Z) %*% W) / den
  nm <- names(struct)
  dth <- vapply(nm, function(k) {
    s1 <- struct; s1[[k]] <- s1[[k]] + h
    s2 <- struct; s2[[k]] <- s2[[k]] - h
    d1 <- .admShiftDelta(spec, s1, X, aref)
    d2 <- .admShiftDelta(spec, s2, X, aref)
    if (is.null(d1) || is.null(d2)) return(rep(NA_real_, length(u)))
    dD <- (as.matrix(d1)[, 1L] - as.matrix(d2)[, 1L]) / (2 * h)
    as.numeric((Ph %*% (W * dD)) / den)
  }, numeric(length(u)))
  if (!is.matrix(dth)) dth <- matrix(dth, length(u), length(nm))
  dimnames(dth) <- list(NULL, nm)
  list(du_dtheta = dth, du_domega = du_dom)
}

# Quadrature over a VECTOR shift, u = Delta(a) + eta in R^m.
#
# u's law is the mixture sum_j W_j N(Delta_j, Omega_S). The affected etas are
# required to be mutually uncorrelated (.admCheckCovariates refuses otherwise),
# so Omega_S is diagonal and the components factor GIVEN the mixture index --
# but not marginally, because the components share that index. The Rosenblatt
# construction handles exactly that:
#
#   u_1        ~ sum_j W_j N(Delta_j1, s_1^2)                      invert its CDF
#   u_2 | u_1  ~ sum_j W_j^(u_1) N(Delta_j2, s_2^2),
#                W_j^(v) proportional to W_j phi((v - Delta_j1)/s_1)
#
# i.e. each conditional is again a ONE-dimensional mixture, so the same Newton
# inversion serves every level. Cost is n_u^m nodes and n_u^(m-1) inversions,
# against n_cov^p * n_node^m for the product grid -- still independent of the
# number of covariates, which is the whole point.
.admShiftNodesMulti <- function(D, W, om, n_u) {
  D <- as.matrix(D)
  m <- ncol(D)
  if (m == 1L) {                       # identical to the scalar path
    un <- .admShiftNodes(D[, 1L], W, om[1L], n_u)
    return(list(u = matrix(un$u, ncol = 1L), w = un$w))
  }
  rec <- function(k, Wc) {
    un <- .admShiftNodes(D[, k], Wc, om[k], n_u)
    if (k == m)
      return(list(u = matrix(un$u, ncol = 1L), w = un$w))
    out_u <- list(); out_w <- numeric(0)
    for (i in seq_along(un$u)) {
      # reweight the components by how well each explains the value just drawn
      pw <- Wc * stats::dnorm((un$u[i] - D[, k]) / om[k])
      sp <- sum(pw)
      pw <- if (sp > 0) pw / sp else Wc
      sub <- rec(k + 1L, pw)
      out_u[[i]] <- cbind(un$u[i], sub$u)
      out_w <- c(out_w, un$w[i] * sub$w)
    }
    list(u = do.call(rbind, out_u), w = out_w)
  }
  r <- rec(1L, W)
  list(u = r$u, w = r$w / sum(r$w))
}

# Does the shift identity actually hold for THIS model? Compiled model, a few
# rows, once per fit. Returns the worst relative discrepancy.
.admShiftVerify <- function(spec, ui, rxMod, pinfo, s, out_var, cov_names,
                            aref, cores = 1L) {
  if (is.null(out_var)) out_var <- tryCatch(.admOutputVar(ui), error = function(e) NULL)
  if (is.null(rxMod))   rxMod   <- tryCatch(.admLoadModel(ui), error = function(e) NULL)
  if (is.null(rxMod) || is.null(out_var)) return(NA_real_)
  cd <- s[["cov_dist"]]
  hi <- vapply(cov_names, function(n) {
    m <- .admCovMeanOf(cd[[n]]); sd <- .admCovSdOf(cd[[n]])
    if (is.null(m) || !is.finite(m)) return(NA_real_)
    if (is.null(sd) || !is.finite(sd)) m else m + sd
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
# u-quantile marginalisation -- ANY functional form
# =============================================================================
#
# The model sees the covariate and the random effect only through their sum in
# the mu-referenced argument, u = Delta(a) + eta. Delta(a) is arithmetic, so u's
# distribution can be pinned down exactly BEFORE any solve, and the whole solve
# budget spent representing it:
#
#     F_u(t) = E_a[ Phi( (t - Delta(a)) / sd_eta ) ]
#
# Nothing here inspects the functional form. Delta is measured from the model
# itself (see .admCovDelta), so power, exponential, allometric, Emax, if/else
# and categorical covariate effects are all handled identically.

# Which parameter assignment does the covariate enter, and which eta shares it?
# The shift structure the whole method rests on is exactly "covariate and eta in
# the same mu-referenced argument", so this also decides whether it applies.
.admCovParamEta <- function(ui, cov, eta_names) {
  lst <- tryCatch(ui$lstExpr, error = function(e) NULL)
  if (is.null(lst)) return(NULL)
  for (e in lst) {
    if (!is.call(e) || length(e) < 3L) next
    if (!as.character(e[[1L]])[1L] %in% c("<-", "=")) next
    v <- all.vars(e[[3L]])
    if (!cov %in% v) next
    et <- intersect(eta_names, v)
    if (length(et) != 1L) next
    return(list(param = as.character(e[[2L]])[1L], eta = et))
  }
  NULL
}

# Delta(a) MEASURED from the model, not parsed out of it.
#
# One rxSolve at eta = 0 over a covariate grid; the mu-referenced parameter comes
# back as an lhs column, and Delta(a) = log(param(a) / param(a_ref)). Exact to
# machine precision (measured 8e-17 on an allometric model) and indifferent to
# how the covariate effect is written.
#
# Costs one extra rxSolve call (~11 ms) per objective evaluation. Delta depends
# on the covariate coefficients, which move, so it cannot be cached across
# iterations -- but 11 ms against a solve budget in seconds is noise.
# ngrid = 32: Gauss-Hermite over a smooth Delta converges long before this.
# Measured against a 200-node reference, F_u agreed to 3.8e-15 at 16 nodes and
# 2.2e-15 at 24, over covariate spreads from sdlog 0.15 to 0.45. The node count
# costs twice -- once in the probe's rxSolve rows, once in every Newton
# iteration's n_sim x n_node matrices -- so 64 was paying for nothing.
.admCovDelta <- function(rxMod, pars, pinfo, s, cov, param_nm, cores,
                         ngrid = 32L) {
  spec  <- s$cov_dist[[cov]]
  a_ref <- s[["cov"]][[cov]]
  if (is.null(a_ref)) return(NULL)
  an <- .admCovANodes(spec, ngrid)
  a  <- an$x
  aa <- c(a, a_ref)
  pm <- .admMakeParamsList(length(aa), pinfo, 1L)[[1L]]
  for (nm in pinfo$struct_names) pm[, nm] <- pars$struct[[nm]]
  pm <- .admCovCols(pm, rxMod$params, NULL,
                    matrix(aa, ncol = 1L, dimnames = list(NULL, cov)))
  out <- tryCatch(rxode2::rxSolve(rxMod, params = as.data.frame(pm),
                                  events = s$ev_full, cores = cores,
                                  nDisplayProgress = .Machine$integer.max),
                  error = function(e) NULL)
  if (is.null(out) || is.null(out[[param_nm]])) return(NULL)
  v <- out[[param_nm]][out[["time"]] == s$times[[1L]]]
  if (length(v) != length(aa) || any(!is.finite(v)) || any(v <= 0)) return(NULL)
  list(a = a, delta = log(v[seq_along(a)]) - log(v[length(aa)]), w = an$w)
}

# The covariate value the model is SOLVED at when the study does not name one.
#
# Every path needs it: the collapse solves at the covariate mean so the model
# itself produces the mean shift, and u-quantile measures Delta relative to it.
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
  if (is.function(spec$quantile))
    return(stats::sd(spec$quantile((seq_len(1024L) - 0.5) / 1024L)))
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
  has <- vapply(studies, function(s) !is.null(s[["cov_dist"]]), logical(1))
  if (!any(has)) return(invisible(NULL))
  cds <- lapply(studies[has], `[[`, "cov_dist")
  for (cv in unique(unlist(lapply(cds, names)))) {
    # only covariates that share an argument with an eta sit on the ridge
    if (is.null(.admCovParamEta(.ui, cv, pinfo$eta_col_names))) next
    sp <- Filter(Negate(is.null), lapply(cds, `[[`, cv))
    if (length(sp) == 0L) next
    mu <- vapply(sp, function(x) .admCovMeanOf(x) %||% NA_real_, numeric(1))
    sd <- vapply(sp, function(x) .admCovSdOf(x)  %||% NA_real_, numeric(1))
    varies <- (length(unique(signif(mu, 10))) > 1L) ||
              (length(unique(signif(sd, 10))) > 1L)
    if (!varies)
      warning("admixr2: the coefficient on covariate '", cv, "' is not ",
              "identifiable from these data. It enters the same argument as a ",
              "random effect, and every study declaring it has the SAME ",
              "covariate distribution, so the likelihood is exactly flat along ",
              "a trade-off between that coefficient, the corresponding fixed ",
              "effect and omega. Identification needs studies whose covariate ",
              "MEANS or SPREADS differ.", call. = FALSE)
  }
  invisible(NULL)
}

# Nodes + weights for the a-integral inside F_u.
#
# GAUSS-HERMITE, not equal-weight quantiles. Here we are INTEGRATING over the
# covariate, and GH carries the tails through its weights; N equal-weight
# quantiles truncate the covariate at the (0.5/N) and (1-0.5/N) points and drop
# the tail mass entirely. That understates the spread of Delta, and the fit then
# inflates the covariate coefficient to make up the missing variance -- measured
# as tcov 0.80 against a truth of 0.75 on an allometric model, while omega stayed
# correct. The mirror-image mistake (using GH to REPRESENT u's distribution) is
# equally wrong; see .admCovUQuantile, which uses the closed-form CDF.
.admCovANodes <- function(spec, ngrid) {
  if (!is.null(spec$values)) {
    pr <- spec$probs %||% rep(1, length(spec$values))
    return(list(x = as.numeric(spec$values), w = pr / sum(pr)))
  }
  if (is.function(spec$quantile)) {           # arbitrary user CDF: no GH rule
    p <- (seq_len(ngrid) - 0.5) / ngrid
    return(list(x = spec$quantile(p), w = rep(1 / ngrid, ngrid)))
  }
  g <- .adghNodes1(ngrid)                     # standard-normal nodes/weights
  if (!is.null(spec$meanlog))
    list(x = exp(spec$meanlog + spec$sdlog * g$x), w = g$w)
  else
    list(x = spec$mu + spec$sd * g$x, w = g$w)
}

# Inverse CDF of u = Delta(a) + eta, evaluated at the uniforms `p`.
# Phi carries the tails analytically, which a grid of summed draws cannot.
# Inverse CDF of u = Delta(a) + eta at the uniforms `p`, by NEWTON on the exact
# mixture CDF rather than interpolation on a grid.
#
# F_u(t) = sum_j w_j Phi((t - Delta_j)/s) is smooth and strictly increasing, and
# its density is available in closed form, so Newton converges in a handful of
# steps from a moment-matched normal start.
#
# Why not the grid: interpolating F_u on a grid makes u PIECEWISE LINEAR in the
# parameters, because Delta and s move the grid's knots. The objective then has
# small kinks, which (a) is the exact mismatch that stops an analytic gradient
# being consistent with the objective it differentiates, and (b) defeats a
# trust-region derivative-free search, which builds quadratic models. Newton
# gives the true implicit function, so d(u)/d(theta) = -(dF_u/d(theta)) / f_u(u)
# is the derivative of what is actually computed.
#
# Bracketed and clamped: f_u underflows in the far tails, so an undamped Newton
# step there can throw the iterate out of the support entirely.
.admCovUQuantile <- function(dl, sd_eta, p, tol = 1e-11, maxit = 50L) {
  d <- dl$delta; w <- dl$w
  lo <- min(d) - 12 * sd_eta
  hi <- max(d) + 12 * sd_eta
  # moment-matched normal start: exact when Delta is degenerate, close otherwise
  m  <- sum(w * d)
  sdu <- sqrt(sum(w * (d - m)^2) + sd_eta^2)
  u  <- pmin(pmax(m + sdu * stats::qnorm(p), lo), hi)
  for (it in seq_len(maxit)) {
    z  <- outer(u, d, "-") / sd_eta
    r  <- as.numeric(stats::pnorm(z) %*% w) - p
    if (max(abs(r)) < tol) break
    fu <- as.numeric(stats::dnorm(z) %*% w) / sd_eta
    u  <- pmin(pmax(u - r / pmax(fu, 1e-300), lo), hi)
  }
  u
}
# Replace the covariate-affected eta column with draws of u = Delta(a) + eta.
#
# The substitution is by INVERSE TRANSFORM of the very dimension that would have
# produced that eta: u = F_u^-1(Phi(z_j)). So the low-discrepancy quality of the
# existing draws carries over, u gets exactly the right marginal, and its
# independence from the other eta columns is preserved by construction -- which
# a sorted set of quantiles pasted in would have destroyed.
#
# Requires Omega DIAGONAL in the affected row (or a single eta): with a
# correlated Omega, column j is a mixture of several z columns and replacing it
# would break the joint distribution. Enforced by .admCheckCovariates().
.admCovUQEta <- function(eta_mat, z, pars, pinfo, s, rxMod, cores) {
  info <- s$.adm_cov_uq
  if (is.null(info)) return(eta_mat)
  for (k in seq_along(info)) {
    it <- info[[k]]
    dl <- .admCovDelta(rxMod, pars, pinfo, s, it$cov, it$param, cores)
    if (is.null(dl)) return(NULL)                 # never silently fall back
    j  <- match(it$eta, pinfo$eta_col_names)
    sd_eta <- sqrt(tcrossprod(pars$L)[j, j])
    u <- .admCovUQuantile(dl, sd_eta, stats::pnorm(z[, j]))
    if (any(!is.finite(u))) return(NULL)
    eta_mat[, j] <- u
  }
  eta_mat
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
  # at the default 11 nodes -- so any covariate with a CV above about 0.19 gets
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
