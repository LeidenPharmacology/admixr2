# =============================================================================
# Covariate marginalisation for aggregate data modelling (admc estimator)
# =============================================================================
#
# Aggregate data are often stratified by a covariate (e.g. body weight): each
# study/arm reports a mean vector + covariance matrix for a sub-population with
# a known covariate value. To recover the population marginal model we treat
# each stratum as a study carrying its covariate value (`cov`) and a quadrature
# `weight`, then combine the per-study -2LL contributions:
#
#   * "gl"/"gh"  -- weighted sum   sum_k w_k * NLL_k
#   * "taylor"   -- 2nd-order (Laplace) correction of the central node NLL
#
# `admBuildQuadrature()` produces the nodes + weights for a covariate
# distribution; `admBuildCovStudies()` pairs those nodes with per-node
# aggregate data to build a `studies` list ready for `admControl(studies=)`.
# Pass the quadrature object via `admControl(quadrature=)` so the estimator
# knows how to combine the per-study contributions.
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
    if (length(nms) == 0L) return(mat)
    if (nrow(cov_rows) != nrow(mat))
      stop(".admCovCols: cov_rows has ", nrow(cov_rows), " rows but the params ",
           "frame has ", nrow(mat), ". Recycling here would hand subjects the ",
           "wrong covariate values silently.", call. = FALSE)
    add <- as.matrix(cov_rows[, nms, drop = FALSE])
    return(if (is.data.frame(mat)) cbind(mat, as.data.frame(add, check.names = FALSE))
           else                    cbind(mat, add))
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

# Covariance matrix of the covariate distribution from a `covariate_dist` spec.
# Entries with $mu and $sd are covariates; "rho" and "Sigma" are metadata.
# A single scalar rho is applied to all pairs; an explicit `Sigma` overrides.
.admGetSigmaCov <- function(covariate_dist) {
  is_cov  <- vapply(covariate_dist,
                    function(x) is.list(x) && !is.null(x[["mu"]]), logical(1))
  cov_ent <- covariate_dist[is_cov]
  d       <- length(cov_ent)
  sds     <- vapply(cov_ent, `[[`, numeric(1), "sd")
  if (!is.null(covariate_dist$Sigma)) return(covariate_dist$Sigma)
  rho <- covariate_dist$rho %||% 0
  S   <- diag(sds^2, d)
  if (d >= 2L && abs(rho) > 1e-12)
    for (i in seq_len(d - 1L)) for (j in seq(i + 1L, d))
      S[i, j] <- S[j, i] <- rho * sds[i] * sds[j]
  S
}

#' Build quadrature nodes and weights for covariate marginalisation
#'
#' Computes the covariate node values and integration weights used to
#' marginalise an aggregate-data model over a (multivariate) normal covariate
#' distribution. Pass the result to [admBuildCovStudies()] to build a study
#' list, and to `admControl(quadrature = )` so the `admc` estimator combines
#' the per-node -2LL contributions correctly.
#'
#' @param covariate_dist Named list describing the covariate distribution. Each
#'   covariate is a list with `mu` and `sd`, e.g.
#'   `list(wt = list(mu = 70, sd = 10))`. Optional `rho` (scalar correlation
#'   applied to all pairs) or `Sigma` (explicit covariance matrix) add
#'   correlation for the Taylor method.
#' @param method `"gl"` (Gauss-Legendre over a truncated normal), `"gh"`
#'   (Gauss-Hermite, exact normal) or `"taylor"` (Laplace, 3 or 5 nodes
#'   per covariate depending on `order`, plus cross nodes when correlated).
#' @param n_nodes Number of nodes per covariate for `"gl"`/`"gh"`.
#' @param truncation_sd Half-width of the Gauss-Legendre integration interval,
#'   in covariate SDs.
#' @param h Finite-difference step (in covariate units) for the `"taylor"`
#'   Hessian approximation. Scalar or one value per covariate dimension. Sets
#'   the spacing between evaluation nodes; the 2nd-order curvature is estimated
#'   as `(NLL₊ − 2·NLL₀ + NLL₋) / h²` and the 4th-order curvature as
#'   `(NLL₊₂ − 4NLL₊ + 6NLL₀ − 4NLL₋ + NLL₋₂) / h⁴`. A step of 1–5
#'   covariate units is typically adequate. Ignored for `"gl"` and `"gh"`.
#' @param order Approximation order for `"taylor"`: `2` (default, 3 nodes per
#'   covariate) or `4` (5 nodes per covariate — adds the
#'   \eqn{\sigma^4\,\mathrm{NLL}^{(4)}/8} correction term). Ignored for
#'   `"gl"` and `"gh"`.
#'
#' @return A list with the quadrature specification: `method`, `cov_names`,
#'   `d` (number of covariates), `mu`, `sd`, `Sigma_cov`, `wt_nodes`
#'   (numeric vector for 1 covariate, `n_node x d` matrix otherwise),
#'   `weights` (`NULL` for Taylor), `h`, `order`, `node_signs` (Taylor only)
#'   and `is_correlated`.
#'
#' @examples
#' q_gl <- admBuildQuadrature(list(wt = list(mu = 70, sd = 10)), "gl", n_nodes = 9L)
#' length(q_gl$wt_nodes)
#'
#' @export
admBuildQuadrature <- function(covariate_dist, method = c("gl", "gh", "taylor"),
                               n_nodes = 9L, truncation_sd = 3.5, h = 2.0,
                               order = 2L) {
  method <- match.arg(method)
  order  <- as.integer(order)
  if (method == "taylor" && !order %in% c(2L, 4L))
    stop("admBuildQuadrature: order must be 2 or 4 for method = 'taylor'.",
         call. = FALSE)
  if (method %in% c("gl", "gh") && !requireNamespace("statmod", quietly = TRUE))
    stop("admBuildQuadrature(method='", method,
         "') requires the 'statmod' package.", call. = FALSE)

  is_cov  <- vapply(covariate_dist,
                    function(x) is.list(x) && !is.null(x[["mu"]]), logical(1))
  cov_ent <- covariate_dist[is_cov]
  d       <- length(cov_ent)
  if (d == 0L) stop("admBuildQuadrature: no covariates with $mu found.", call. = FALSE)
  nms     <- names(cov_ent)
  mu      <- vapply(cov_ent, `[[`, numeric(1), "mu")
  sd_vec  <- vapply(cov_ent, `[[`, numeric(1), "sd")
  Sigma_cov <- .admGetSigmaCov(covariate_dist)

  if (d == 1L) {
    mu1 <- mu[[1]]; sd1 <- sd_vec[[1]]
    h1  <- h[[1L]]
    if (method == "gl") {
      gl      <- statmod::gauss.quad(n_nodes, kind = "legendre")
      wt_lo   <- mu1 - truncation_sd * sd1
      wt_hi   <- mu1 + truncation_sd * sd1
      scale   <- (wt_hi - wt_lo) / 2
      nodes   <- wt_lo + scale * (gl$nodes + 1)
      w_raw   <- scale * gl$weights * dnorm(nodes, mu1, sd1)
      weights <- w_raw / sum(w_raw)
    } else if (method == "gh") {
      gh      <- statmod::gauss.quad(n_nodes, kind = "hermite")
      nodes   <- mu1 + sqrt(2) * sd1 * gh$nodes
      weights <- gh$weights / sqrt(pi)
    } else {
      if (order == 4L) {
        nodes      <- c(mu1 - 2*h1, mu1 - h1, mu1, mu1 + h1, mu1 + 2*h1)
        node_signs <- matrix(c(-2L, -1L, 0L, 1L, 2L), 5L, 1L)
      } else {
        nodes      <- c(mu1 - h1, mu1, mu1 + h1)
        node_signs <- matrix(c(-1L, 0L, 1L), 3L, 1L)
      }
      weights <- NULL
    }
    return(list(method        = method,
                cov_name      = nms,            # backward-compat alias
                cov_names     = nms,
                d             = 1L,
                mu            = mu1, sd = sd1,
                Sigma_cov     = Sigma_cov,
                wt_nodes      = nodes,          # numeric vector for 1-d
                weights       = weights,
                h             = h1,
                order         = if (method == "taylor") order else 2L,
                node_signs    = if (method == "taylor") node_signs else NULL,
                is_correlated = FALSE))
  }

  # ---- Multi-covariate -------------------------------------------------------
  if (method %in% c("gl", "gh")) {
    one_d  <- lapply(seq_len(d), function(j)
      admBuildQuadrature(cov_ent[j], method, n_nodes, truncation_sd, h))
    args_n <- setNames(lapply(one_d, `[[`, "wt_nodes"), nms)
    args_w <- setNames(lapply(one_d, `[[`, "weights"),  nms)
    grid_n <- do.call(expand.grid, args_n)
    grid_w <- do.call(expand.grid, args_w)
    return(list(method        = method,
                cov_names     = nms, d = d,
                mu = mu, sd = sd_vec, Sigma_cov = Sigma_cov,
                wt_nodes      = as.matrix(grid_n),               # n_total x d
                weights       = apply(as.matrix(grid_w), 1L, prod),
                h             = h,
                order         = 2L,
                node_signs    = NULL, is_correlated = FALSE))
  }

  # Taylor (multi-d): 1 central node + axis nodes (±1 per dim, always) +
  # ±2 axis nodes per dim when order = 4 + cross nodes (±1,±1) for correlated
  # pairs. Mixed 4th-order cross terms omitted (would require (±1,±1,±2) etc.).
  h_vec <- if (length(h) == 1L) rep(h, d) else h
  is_correlated <- any(abs(Sigma_cov[lower.tri(Sigma_cov)]) > 1e-12)

  signs <- matrix(0L, 1L, d)
  for (j in seq_len(d)) {
    ep1 <- integer(d); ep1[j] <- +1L
    em1 <- integer(d); em1[j] <- -1L
    signs <- rbind(signs, ep1, em1)
    if (order == 4L) {
      ep2 <- integer(d); ep2[j] <- +2L
      em2 <- integer(d); em2[j] <- -2L
      signs <- rbind(signs, ep2, em2)
    }
  }
  if (is_correlated) {
    for (i in seq_len(d - 1L)) for (j in seq(i + 1L, d)) {
      if (abs(Sigma_cov[i, j]) < 1e-12) next
      for (si in c(+1L, -1L)) for (sj in c(+1L, -1L)) {
        s <- integer(d); s[i] <- si; s[j] <- sj
        signs <- rbind(signs, s)
      }
    }
  }
  rownames(signs) <- NULL
  wt_mat <- sweep(sweep(signs, 2L, h_vec, `*`), 2L, mu, `+`)
  colnames(wt_mat) <- nms

  list(method        = "taylor",
       cov_names     = nms, d = d,
       mu = mu, sd = sd_vec, Sigma_cov = Sigma_cov,
       wt_nodes      = wt_mat,           # n_node x d
       weights       = NULL,
       h             = h_vec,
       order         = order,
       node_signs    = signs,
       is_correlated = is_correlated)
}

#' Build a covariate-stratified study list from quadrature nodes
#'
#' Pairs each quadrature node (covariate value + weight) from
#' [admBuildQuadrature()] with the corresponding per-node aggregate data to
#' produce a `studies` list for `admControl(studies = )`.
#'
#' @param agg_list List of per-node aggregate data, one element per quadrature
#'   node (matching `quad$wt_nodes` order). Each element is a list with `E`
#'   (mean vector) and `V` (covariance matrix or variance vector).
#' @param quad Quadrature object from [admBuildQuadrature()].
#' @param ev `rxode2::et()` dosing event table (shared across nodes).
#' @param times Numeric vector of observation times (shared across nodes).
#' @param n Sample size used for -2LL scaling (per node).
#' @param prefix Study-name prefix; nodes are named `<prefix>01`, `<prefix>02`, ...
#'
#' @return A named list of normalised study specifications, each carrying `cov`
#'   (covariate values at that node) and `weight` (quadrature weight).
#'
#' @export
admBuildCovStudies <- function(agg_list, quad, ev, times, n, prefix = "study") {
  multi_d <- is.matrix(quad$wt_nodes)
  nms     <- if (multi_d) quad$cov_names else (quad$cov_names %||% quad$cov_name)

  n_nodes <- if (multi_d) nrow(quad$wt_nodes) else length(quad$wt_nodes)

  if (length(agg_list) != n_nodes)
    stop(sprintf(
      "admBuildCovStudies: length(agg_list) (%d) must match the number of quadrature nodes (%d).",
      length(agg_list), n_nodes), call. = FALSE)

  # The COMBINATION coefficient, not the raw quadrature weight. For gl/gh they
  # are the same thing; for taylor the quadrature has no weights at all and the
  # coefficients are the Laplace stencil (1 - sigma^2/h^2 in the centre, so
  # negative). Taking `quad$weights[k] %||% 1` here gave every Taylor node a
  # weight of 1, i.e. summed the stencil's evaluation points as if they were a
  # quadrature -- a different objective, silently.
  coefs <- .admCovNodeCoefs(quad)

  out <- lapply(seq_along(agg_list), function(k) {
    a <- agg_list[[k]]
    node_vals <- if (multi_d)
      setNames(quad$wt_nodes[k, , drop = TRUE], nms)
    else
      setNames(quad$wt_nodes[[k]], nms)
    nm <- sprintf("%s%02d", prefix, k)
    s  <- .admNormaliseStudy(
      list(E = a$E, V = a$V, n = n, times = times, ev = ev,
           cov = node_vals, weight = coefs[[k]]),
      nm)
    s
  })
  setNames(out, vapply(seq_along(agg_list),
                       function(k) sprintf("%s%02d", prefix, k), character(1)))
}

# -- Combine per-study contributions -------------------------------------------

# Weighted sum (gl/gh or no quadrature) or Laplace approximation (Taylor) of
# the per-study -2LL values.
# Taylor order=2: NLL_0 + 0.5 * tr(Sigma_cov * H_wt)
#   H_wt[j,j]  <- (NLL_+1 - 2*NLL_0 + NLL_-1) / h^2   (2nd FD)
#   H_wt[i,j]  <- cross-node FD when covariates are correlated
# Taylor order=4: adds sigma_j^4/8 * (NLL_+2 - 4*NLL_+1 + 6*NLL_0 - 4*NLL_-1 + NLL_-2) / h^4
#   (diagonal 4th-derivative terms; mixed cross-4th terms omitted)
.adm_combine_nll <- function(nll_vec, studies, quad) {
  if (!is.null(quad) && identical(quad$method, "taylor")) {
    ord   <- quad$order %||% 2L
    ci    <- which(rowSums(abs(quad$node_signs)) == 0L)
    nll_0 <- nll_vec[[ci]]
    d     <- quad$d
    h_vec <- rep_len(quad$h, d)

    # 2nd-order diagonal: sigma_j^2 * NLL''_j
    hess2 <- 0
    for (j in seq_len(d)) {
      other <- if (d > 1L) rowSums(abs(quad$node_signs[, -j, drop = FALSE])) == 0L
               else        rep(TRUE, nrow(quad$node_signs))
      pi <- which(quad$node_signs[, j] == +1L & other)
      mi <- which(quad$node_signs[, j] == -1L & other)
      hess2 <- hess2 + quad$Sigma_cov[j, j] *
                 (nll_vec[[pi]] - 2 * nll_0 + nll_vec[[mi]]) / h_vec[j]^2
    }
    # 2nd-order cross-terms (correlated covariates)
    if (isTRUE(quad$is_correlated) && d > 1L) {
      for (i in seq_len(d - 1L)) for (j in seq(i + 1L, d)) {
        if (abs(quad$Sigma_cov[i, j]) < 1e-12) next
        pp <- which(quad$node_signs[, i] == +1L & quad$node_signs[, j] == +1L)
        pm <- which(quad$node_signs[, i] == +1L & quad$node_signs[, j] == -1L)
        mp <- which(quad$node_signs[, i] == -1L & quad$node_signs[, j] == +1L)
        mm <- which(quad$node_signs[, i] == -1L & quad$node_signs[, j] == -1L)
        hess2 <- hess2 + quad$Sigma_cov[i, j] *
                   (nll_vec[[pp]] - nll_vec[[pm]] - nll_vec[[mp]] + nll_vec[[mm]]) /
                   (h_vec[i] * h_vec[j])
      }
    }
    result <- nll_0 + 0.5 * hess2

    # 4th-order diagonal: sigma_j^4/8 * NLL''''_j (diagonal terms only)
    if (ord >= 4L) {
      hess4 <- 0
      for (j in seq_len(d)) {
        other <- if (d > 1L) rowSums(abs(quad$node_signs[, -j, drop = FALSE])) == 0L
                 else        rep(TRUE, nrow(quad$node_signs))
        p2i <- which(quad$node_signs[, j] == +2L & other)
        m2i <- which(quad$node_signs[, j] == -2L & other)
        p1i <- which(quad$node_signs[, j] == +1L & other)
        m1i <- which(quad$node_signs[, j] == -1L & other)
        hess4 <- hess4 + quad$Sigma_cov[j, j]^2 / 8 *
                   (nll_vec[[p2i]] - 4*nll_vec[[p1i]] + 6*nll_0 -
                    4*nll_vec[[m1i]] + nll_vec[[m2i]]) / h_vec[j]^4
      }
      result <- result + hess4
    }
    return(result)
  }
  w <- vapply(studies, function(s) s$weight %||% 1, numeric(1))
  sum(w * nll_vec)
}

# Combine per-study gradient vectors consistently with .adm_combine_nll.
.adm_combine_grad <- function(grad_list, nll_vec, studies, quad) {
  if (!is.null(quad) && identical(quad$method, "taylor")) {
    ord   <- quad$order %||% 2L
    ci    <- which(rowSums(abs(quad$node_signs)) == 0L)
    g_0   <- grad_list[[ci]]
    d     <- quad$d
    h_vec <- rep_len(quad$h, d)

    # 2nd-order diagonal
    hess2_g <- numeric(length(g_0))
    for (j in seq_len(d)) {
      other <- if (d > 1L) rowSums(abs(quad$node_signs[, -j, drop = FALSE])) == 0L
               else        rep(TRUE, nrow(quad$node_signs))
      pi <- which(quad$node_signs[, j] == +1L & other)
      mi <- which(quad$node_signs[, j] == -1L & other)
      hess2_g <- hess2_g + quad$Sigma_cov[j, j] *
                   (grad_list[[pi]] - 2 * g_0 + grad_list[[mi]]) / h_vec[j]^2
    }
    # 2nd-order cross-terms (correlated covariates)
    if (isTRUE(quad$is_correlated) && d > 1L) {
      for (i in seq_len(d - 1L)) for (j in seq(i + 1L, d)) {
        if (abs(quad$Sigma_cov[i, j]) < 1e-12) next
        pp <- which(quad$node_signs[, i] == +1L & quad$node_signs[, j] == +1L)
        pm <- which(quad$node_signs[, i] == +1L & quad$node_signs[, j] == -1L)
        mp <- which(quad$node_signs[, i] == -1L & quad$node_signs[, j] == +1L)
        mm <- which(quad$node_signs[, i] == -1L & quad$node_signs[, j] == -1L)
        hess2_g <- hess2_g + quad$Sigma_cov[i, j] *
                     (grad_list[[pp]] - grad_list[[pm]] - grad_list[[mp]] + grad_list[[mm]]) /
                     (h_vec[i] * h_vec[j])
      }
    }
    result_g <- g_0 + 0.5 * hess2_g

    # 4th-order diagonal
    if (ord >= 4L) {
      hess4_g <- numeric(length(g_0))
      for (j in seq_len(d)) {
        other <- if (d > 1L) rowSums(abs(quad$node_signs[, -j, drop = FALSE])) == 0L
                 else        rep(TRUE, nrow(quad$node_signs))
        p2i <- which(quad$node_signs[, j] == +2L & other)
        m2i <- which(quad$node_signs[, j] == -2L & other)
        p1i <- which(quad$node_signs[, j] == +1L & other)
        m1i <- which(quad$node_signs[, j] == -1L & other)
        hess4_g <- hess4_g + quad$Sigma_cov[j, j]^2 / 8 *
                     (grad_list[[p2i]] - 4*grad_list[[p1i]] + 6*g_0 -
                      4*grad_list[[m1i]] + grad_list[[m2i]]) / h_vec[j]^4
      }
      result_g <- result_g + hess4_g
    }
    return(result_g)
  }
  w <- vapply(studies, function(s) s$weight %||% 1, numeric(1))
  Reduce(`+`, Map(function(g, wi) wi * g, grad_list, w))
}

# =============================================================================
# Node batching
# =============================================================================
#
# Covariate-marginalisation studies differ only in their covariate value, so
# every node can ride in ONE rxSolve rather than one per node -- an rxSolve call
# costs ~11 ms before it integrates anything, while a marginal row costs
# ~0.015 ms, so stacking nodes as ROWS is close to free.
#
# NOT currently called by any estimator: admc's NLL/gradient were rewritten
# around the residual-row arrays and the 8-argument C++ kernels, and these
# helpers predate that. They are kept, unit-tested and corrected, as the basis
# for rewiring covariate marginalisation into the current admc.

# Covariate column names to add to a stacked params frame: only names the node
# studies actually DECLARE in `cov` and the model actually reads. Deliberately
# not setdiff(mod$params, colnames(mat)) -- see .admCovCols() for the two bugs
# that pattern caused.
.admCovNodeNames <- function(mat, mod_params, studies) {
  nms <- unique(unlist(lapply(studies, function(s) names(s[["cov"]]))))
  setdiff(intersect(nms, mod_params), colnames(mat))
}

# One node's covariate values in `nms` order; 0 for a name it does not declare.
.admCovNodeVals <- function(cov_s, nms)
  vapply(nms, function(nm)
    if (!is.null(cov_s) && nm %in% names(cov_s)) as.numeric(cov_s[[nm]]) else 0,
    numeric(1))

# Studies are "node-batchable" when they share identical events and observation
# times, differing only in covariate values (and possibly their eta/z draws).
# The covariate-marginalisation studies built by datagen(covariate=)
# satisfy this condition, so all nodes can be stacked into a single rxSolve call.
# Each node may carry independent eta draws (GL) or shared CRN draws (Taylor);
# both are handled by .adm_sim_nodes_batched / .adm_grad_presim via eta_mat_list.
.adm_nodes_batchable <- function(studies) {
  n <- length(studies)
  if (n < 2L) return(FALSE)
  t1  <- studies[[1L]]$times
  ev1 <- tryCatch(as.data.frame(studies[[1L]]$ev_full), error = function(e) NULL)
  if (is.null(ev1)) return(FALSE)
  for (i in 2:n) {
    if (!isTRUE(all.equal(studies[[i]]$times, t1))) return(FALSE)
    evi <- tryCatch(as.data.frame(studies[[i]]$ev_full), error = function(e) NULL)
    if (is.null(evi) || !isTRUE(all.equal(evi, ev1))) return(FALSE)
  }
  TRUE
}

# Stack node-batchable studies into one rxSolve (per-node eta draws, per-node
# covariate). `eta_mat_list` is a list of n_sim x n_eta matrices, one per node
# (identical matrices for CRN/Taylor, independent for GL). Returns a list of
# n_sim x n_t prediction matrices, one per study, or NULL on failure.
# `solve_fn` is injectable for testing.
.adm_sim_nodes_batched <- function(rxMod, struct, sigma_names, eta_mat_list, studies,
                                   output_var, cores, solve_fn = rxode2::rxSolve) {
  n_nodes <- length(studies)
  n_sim   <- nrow(eta_mat_list[[1L]])
  n_t     <- length(studies[[1L]]$times)
  eta_col <- colnames(eta_mat_list[[1L]])
  base_cols <- c(names(struct), eta_col, sigma_names, paste0("rxerr.", output_var))
  M <- matrix(0, n_nodes * n_sim, length(base_cols), dimnames = list(NULL, base_cols))
  M[, paste0("rxerr.", output_var)] <- 1L
  for (nm in names(struct)) M[, nm] <- struct[nm]
  if (length(eta_col) > 0L)
    for (i in seq_len(n_nodes))
      M[(i - 1L) * n_sim + seq_len(n_sim), eta_col] <- eta_mat_list[[i]]

  extra <- .admCovNodeNames(M, rxMod$params, studies)
  if (length(extra) > 0L) {
    cm <- matrix(0, n_nodes * n_sim, length(extra), dimnames = list(NULL, extra))
    for (i in seq_len(n_nodes))
      cm[(i - 1L) * n_sim + seq_len(n_sim), ] <-
        matrix(rep(.admCovNodeVals(studies[[i]][["cov"]], extra), each = n_sim), n_sim, length(extra))
    M <- cbind(M, cm)
  }

  out <- tryCatch(
    solve_fn(rxMod, params = as.data.frame(M), events = studies[[1L]]$ev_full,
             cores = cores, nDisplayProgress = .Machine$integer.max),
    error = function(e) NULL)
  if (is.null(out)) return(NULL)
  keep <- out[["time"]] %in% studies[[1L]]$times
  vals <- out[[output_var]][keep]
  if (is.null(vals)) vals <- out[["ipredSim"]][keep]
  lapply(seq_len(n_nodes), function(i) {
    idx <- (i - 1L) * n_sim * n_t + seq_len(n_sim * n_t)
    matrix(vals[idx], nrow = n_sim, ncol = n_t, byrow = TRUE)
  })
}

# Batched forward simulation for the gradient across node-batchable studies:
# one rxSolve for the base predictions + eta sensitivities (sensitivity model or
# FD) and one/two more for the unpaired-theta perturbations -- all nodes share
# struct/eta and differ only by covariate. Returns the per-node inputs the
# per-node gradient assembly consumes (cp_mat, dpred_list, batched_hi,
# batched_lo), or NULL on any failure so the caller falls back to per-node.
# Mirrors the stacking/slicing of .admGradBatch (configs -> nodes). `solve_fn`
# is injectable for testing.
.adm_grad_presim <- function(studies, eta_mat_list, pars, pinfo, rxMod, sensModel,
                             output_var, cores, h, use_central, n_eta, n_sim,
                             eta_col_names, unpaired_k, solve_fn = rxode2::rxSolve) {
  n_nodes <- length(studies)
  n_t     <- length(studies[[1L]]$times)
  n_unp   <- length(unpaired_k)
  ev_full <- studies[[1L]]$ev_full
  obs_t   <- studies[[1L]]$times

  cp_list    <- vector("list", n_nodes)
  dpred_list <- vector("list", n_nodes)

  # --- main solve: base predictions + d(pred)/d(eta) ------------------------
  use_sens <- !is.null(sensModel) && n_eta > 0L
  if (use_sens) {
    rmap      <- sensModel$rename_map
    all_src   <- c(pinfo$struct_names, pinfo$sigma_names, eta_col_names)
    inner_nms <- rmap[all_src]; inner_nms <- inner_nms[!is.na(inner_nms)]
    inner_df  <- as.data.frame(matrix(0, n_nodes * n_sim, length(inner_nms),
                                      dimnames = list(NULL, unname(inner_nms))),
                               check.names = FALSE)
    for (nm in pinfo$struct_names) {
      mp <- rmap[nm]; if (!is.na(mp)) inner_df[[mp]][] <- pars$struct[nm]
    }
    for (j in seq_along(eta_col_names)) {
      mp <- rmap[eta_col_names[j]]
      if (!is.na(mp)) {
        vals <- numeric(n_nodes * n_sim)
        for (i in seq_len(n_nodes))
          vals[(i - 1L) * n_sim + seq_len(n_sim)] <- eta_mat_list[[i]][, j]
        inner_df[[mp]] <- vals
      }
    }
    extra_s <- .admCovNodeNames(inner_df, sensModel$mod$params, studies)
    if (length(extra_s) > 0L) {
      cm <- matrix(0, n_nodes * n_sim, length(extra_s), dimnames = list(NULL, extra_s))
      for (i in seq_len(n_nodes))
        cm[(i - 1L) * n_sim + seq_len(n_sim), ] <-
          matrix(rep(.admCovNodeVals(studies[[i]][["cov"]], extra_s), each = n_sim), n_sim, length(extra_s))
      inner_df <- cbind(inner_df, as.data.frame(cm, check.names = FALSE))
    }
    out <- tryCatch(suppressWarnings(
      solve_fn(sensModel$mod, params = inner_df, events = ev_full, cores = cores,
               nDisplayProgress = .Machine$integer.max)), error = function(e) NULL)
    if (is.null(out) || !all(sensModel$sens_cols %in% names(out))) {
      use_sens <- FALSE
    } else {
      keep      <- out[["time"]] %in% obs_t
      vals_pred <- out[["rx_pred_"]][keep]
      vals_sens <- lapply(sensModel$sens_cols, function(col) out[[col]][keep])
      for (i in seq_len(n_nodes)) {
        idx <- (i - 1L) * n_sim * n_t + seq_len(n_sim * n_t)
        cp_list[[i]]    <- matrix(vals_pred[idx], n_sim, n_t, byrow = TRUE)
        dpred_list[[i]] <- lapply(vals_sens, function(vs)
          matrix(vs[idx], n_sim, n_t, byrow = TRUE))
      }
    }
  }
  if (!use_sens) {
    # FD: per node a block of [base, eta perturbations] (central or forward).
    n_blk <- if (n_eta > 0L) (if (use_central) 1L + 2L * n_eta else 1L + n_eta) else 1L
    base_cols <- c(pinfo$struct_names, eta_col_names, pinfo$sigma_names, paste0("rxerr.", output_var))
    M <- matrix(0, n_nodes * n_blk * n_sim, length(base_cols),
                dimnames = list(NULL, base_cols))
    M[, paste0("rxerr.", output_var)] <- 1L
    for (nm in names(pars$struct)) M[, nm] <- pars$struct[nm]
    if (n_eta > 0L) for (i in seq_len(n_nodes)) {
      cfg    <- (i - 1L) * n_blk * n_sim
      em_i   <- eta_mat_list[[i]]
      M[cfg + seq_len(n_sim), eta_col_names] <- em_i
      if (use_central) {
        for (j in seq_len(n_eta)) {
          rh <- cfg + n_sim * (2L*j - 1L) + seq_len(n_sim)
          rl <- cfg + n_sim * (2L*j)      + seq_len(n_sim)
          eh <- em_i; eh[, j] <- eh[, j] + h
          el <- em_i; el[, j] <- el[, j] - h
          M[rh, eta_col_names] <- eh; M[rl, eta_col_names] <- el
        }
      } else {
        for (j in seq_len(n_eta)) {
          rh <- cfg + n_sim * j + seq_len(n_sim)
          eh <- em_i; eh[, j] <- eh[, j] + h
          M[rh, eta_col_names] <- eh
        }
      }
    }
    extra_m <- .admCovNodeNames(M, rxMod$params, studies)
    if (length(extra_m) > 0L) {
      cm <- matrix(0, nrow(M), length(extra_m), dimnames = list(NULL, extra_m))
      for (i in seq_len(n_nodes)) {
        rows <- (i - 1L) * n_blk * n_sim + seq_len(n_blk * n_sim)
        cm[rows, ] <- matrix(rep(.admCovNodeVals(studies[[i]][["cov"]], extra_m), each = n_blk * n_sim),
                             n_blk * n_sim, length(extra_m))
      }
      M <- cbind(M, cm)
    }
    out <- tryCatch(solve_fn(rxMod, params = as.data.frame(M), events = ev_full,
                    cores = cores, nDisplayProgress = .Machine$integer.max),
                    error = function(e) NULL)
    if (is.null(out)) return(NULL)
    keep <- out[["time"]] %in% obs_t
    vals <- out[[output_var]][keep]; if (is.null(vals)) vals <- out[["ipredSim"]][keep]
    for (i in seq_len(n_nodes)) {
      cfg <- (i - 1L) * n_blk * n_sim * n_t
      cp_list[[i]] <- matrix(vals[cfg + seq_len(n_sim * n_t)], n_sim, n_t, byrow = TRUE)
      dpred_list[[i]] <- if (n_eta > 0L) {
        if (use_central) lapply(seq_len(n_eta), function(j) {
          oh <- cfg + n_sim * n_t * (2L*j - 1L); ol <- cfg + n_sim * n_t * (2L*j)
          (matrix(vals[oh + seq_len(n_sim*n_t)], n_sim, n_t, byrow = TRUE) -
           matrix(vals[ol + seq_len(n_sim*n_t)], n_sim, n_t, byrow = TRUE)) / (2 * h)
        }) else lapply(seq_len(n_eta), function(j) {
          oh <- cfg + n_sim * n_t * j
          (matrix(vals[oh + seq_len(n_sim*n_t)], n_sim, n_t, byrow = TRUE) -
           cp_list[[i]]) / h
        })
      } else list()
    }
  }
  if (any(vapply(cp_list, function(m) is.null(m) || anyNA(m), logical(1)))) return(NULL)

  # --- unpaired-theta perturbation solve(s): all node x unpaired stacked ----
  hi_list <- vector("list", n_nodes)
  lo_list <- if (use_central) vector("list", n_nodes) else NULL
  for (i in seq_len(n_nodes)) {
    hi_list[[i]] <- vector("list", n_unp)
    if (use_central) lo_list[[i]] <- vector("list", n_unp)
  }
  if (n_unp > 0L) {
    base_cols <- c(pinfo$struct_names, eta_col_names, pinfo$sigma_names, paste0("rxerr.", output_var))
    n_cu      <- n_nodes * n_unp
    .solve_unp <- function(sign_h) {
      M <- matrix(0, n_cu * n_sim, length(base_cols), dimnames = list(NULL, base_cols))
      M[, paste0("rxerr.", output_var)] <- 1L
      for (nm in names(pars$struct)) M[, nm] <- pars$struct[nm]
      cu <- 0L
      for (i in seq_len(n_nodes)) for (bi in seq_len(n_unp)) {
        cu <- cu + 1L; rows <- (cu - 1L) * n_sim + seq_len(n_sim)
        if (n_eta > 0L) M[rows, eta_col_names] <- eta_mat_list[[i]]
        nm_u <- pinfo$struct_names[unpaired_k[bi]]
        M[rows, nm_u] <- pars$struct[nm_u] + sign_h * h
      }
      extra <- .admCovNodeNames(M, rxMod$params, studies)
      if (length(extra) > 0L) {
        cm <- matrix(0, nrow(M), length(extra), dimnames = list(NULL, extra)); cu <- 0L
        for (i in seq_len(n_nodes)) for (bi in seq_len(n_unp)) {
          cu <- cu + 1L; rows <- (cu - 1L) * n_sim + seq_len(n_sim)
          cm[rows, ] <- matrix(rep(.admCovNodeVals(studies[[i]][["cov"]], extra), each = n_sim),
                               n_sim, length(extra))
        }
        M <- cbind(M, cm)
      }
      out <- tryCatch(solve_fn(rxMod, params = as.data.frame(M), events = ev_full,
                      cores = cores, nDisplayProgress = .Machine$integer.max),
                      error = function(e) NULL)
      if (is.null(out)) return(NULL)
      keep <- out[["time"]] %in% obs_t
      v <- out[[output_var]][keep]; if (is.null(v)) v <- out[["ipredSim"]][keep]
      v
    }
    vhi <- .solve_unp(+1); if (is.null(vhi)) return(NULL)
    vlo <- if (use_central) { x <- .solve_unp(-1); if (is.null(x)) return(NULL); x } else NULL
    cu <- 0L
    for (i in seq_len(n_nodes)) for (bi in seq_len(n_unp)) {
      cu <- cu + 1L; idx <- (cu - 1L) * n_sim * n_t + seq_len(n_sim * n_t)
      hi_list[[i]][[bi]] <- matrix(vhi[idx], n_sim, n_t, byrow = TRUE)
      if (use_central) lo_list[[i]][[bi]] <- matrix(vlo[idx], n_sim, n_t, byrow = TRUE)
    }
  }

  list(cp = cp_list, dpred = dpred_list, hi = hi_list, lo = lo_list)
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
  sds <- vapply(nms, function(n) as.numeric(cov_dist[[n]]$sd), numeric(1))
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
    cf <- pars$struct[[use$coef[k]]]
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
  # A NODE study carries no `cov_dist` -- its covariate is a fixed value -- so
  # the cov_dist test below does not see it. What makes it a node is its
  # combination coefficient, and only the estimators that call
  # .admCovApplyNodeWeights() apply one. Here the coefficient would simply be
  # ignored and the nodes summed with weight 1: for gl/gh that inflates the
  # objective by ~1/w, for taylor it sums three stencil points as though they
  # were a quadrature. Both are finite, plausible and wrong.
  wt <- vapply(studies, function(s) {
    x <- s[["weight"]] %||% s[[".adm_node_c"]]
    if (is.null(x)) 1 else as.numeric(x)[[1L]]
  }, numeric(1))
  if (any(wt != 1))
    stop("admixr2: `", est, "` does not support node-quadrature covariate ",
         "studies (gl/gh/taylor). Stud",
         if (sum(wt != 1) > 1L) "ies " else "y ",
         paste(sQuote(names(studies)[wt != 1]), collapse = ", "),
         " carr", if (sum(wt != 1) > 1L) "y " else "ies ",
         "a combination coefficient that only `admc` and `adgh` apply; ",
         "summing the nodes here would give a different objective, not a ",
         "coarser one.", call. = FALSE)
  has <- vapply(studies, function(s) !is.null(s[["cov_dist"]]), logical(1))
  if (!any(has)) return(invisible(NULL))
  stop("admixr2: `", est, "` does not support covariate marginalisation. ",
       "Stud", if (sum(has) > 1L) "ies " else "y ",
       paste(sQuote(names(studies)[has]), collapse = ", "),
       " declare(s) `cov_dist`. Use `admc`, which marginalises over the ",
       "covariate distribution; running here would silently solve at the ",
       "covariate mean and inflate omega instead.", call. = FALSE)
}

.admCheckCovariates <- function(.ui, pinfo, studies, grad, est = NULL) {
  has <- vapply(studies, function(s) !is.null(s$cov_dist), logical(1))
  if (!any(has)) return(studies)
  bad <- function(...) stop("admixr2: ", ..., call. = FALSE)

  cmap <- pinfo$cov_map
  covs <- tryCatch(.ui$allCovs, error = function(e) character(0))
  for (nm in names(studies)[has]) {
    cd <- studies[[nm]]$cov_dist
    ok_collapse <- pinfo$n_eta > 0L
    ok_uq       <- pinfo$n_eta > 0L
    uq          <- list()

    for (cv in names(cd)) {
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
      if (!normal || is.null(cmap) || !cv %in% cmap$covariate || !once)
        ok_collapse <- FALSE
      # u-QUANTILE: any functional form, but the covariate's WHOLE effect has to
      # fit in one eta column -- so it must appear exactly once, in a parameter
      # assignment that also carries an eta, with a reference value to measure
      # Delta against.
      if (!once || is.null(pe) || is.null(studies[[nm]][["cov"]][[cv]])) ok_uq <- FALSE
      else uq[[length(uq) + 1L]] <- list(cov = cv, param = pe$param, eta = pe$eta,
                                         idx = match(pe$eta, pinfo$eta_col_names))
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


    # Most efficient VALID path wins. "rows" assumes nothing at all: every
    # simulated subject carries its own covariate value, so rxode2 evaluates the
    # whole model -- covariate on several parameters, on a parameter with no eta,
    # or interacting with one. It is the only path with no structural
    # precondition, hence the fallback.
    studies[[nm]]$.adm_cov_path <- if (ok_collapse) "collapse"
                                   else if (ok_uq)  "uq" else "rows"
    studies[[nm]]$.adm_cov_collapse <- ok_collapse
    studies[[nm]]$.adm_cov_uq <- if (ok_uq) uq else NULL
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
.admCovRowsFor <- function(cov_dist, n, n_eta) {
  nms <- names(cov_dist)
  d   <- length(nms)
  u   <- randtoolbox::sobol(n, dim = n_eta + d)
  if (!is.matrix(u)) u <- matrix(u, nrow = n)
  u   <- u[, n_eta + seq_len(d), drop = FALSE]
  # sobol emits an exact 0 in its first row; qnorm(0) is -Inf.
  u   <- pmin(pmax(u, .Machine$double.eps), 1 - .Machine$double.eps)
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
  if (!is.null(spec$meanlog)) list(x = exp(spec$meanlog + spec$sdlog * g$x), w = g$w)
  else                        list(x = spec$mu + spec$sd * g$x,              w = g$w)
}

# Product grid over several covariates: every combination, weights multiplied.
.admCovGrid <- function(cov_dist, n_nodes) {
  nms <- names(cov_dist)
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


# =============================================================================
# Node objectives: gl / gh / taylor, on the data a study actually reports
# =============================================================================
#
# These are the published alternatives to marginalisation. Each scores ONE
# study's aggregate (E, V) at a set of fixed covariate values and combines the
# per-node -2LL values linearly:
#
#     -2LL_study = sum_k c_k * NLL( E_obs, V_obs ; prediction at covariate a_k )
#
# All three are that sum; only the coefficients c_k differ. Taylor's stencil,
#
#     NLL(mu) + 0.5*sigma^2 * [NLL(mu+h) - 2 NLL(mu) + NLL(mu-h)] / h^2 ,
#
# is a linear combination too, so it folds into the same shape instead of
# needing a combiner of its own -- and the gradient of a linear combination is
# that combination of the gradients, which is why no gradient path has to know
# any of this exists.
#
# IMPLEMENTATION. c_k is folded into the node study's `n`. Every C++ kernel
# takes `n` as a double and uses it as a linear multiplier on
# (log|V| + tr(V^-1 V_obs) + r'V^-1 r), so n_k = c_k * n IS c_k * NLL_k --
# exactly, for the NLL, for the analytic gradient and for both batch paths,
# with no accumulation site aware that nodes exist. Taylor's central
# coefficient is 1 - sigma^2/h^2 and so is normally NEGATIVE: that `n` is a
# combination coefficient, not a subject count, which is why .admCovExpandNodes
# writes it after .admNormaliseStudy has done its validation.
#
# CONTRAST with admBuildCovStudies(), which pairs each node with its OWN
# aggregate data. That construction recovers truth exactly -- see
# validation/covariate-matched-conditional.R -- but it needs one (E, V) per
# node, 17 datasets for GH-17, which no published study reports. Scoring every
# node against the single (E, V) a study does report is what these methods
# amount to in practice, and is what this builds.

.ADM_COV_NODE_METHODS <- c("gl", "gh", "taylor")

# Linear combination coefficients, one per node, such that the study objective
# is sum_k c_k * NLL_k. For gl/gh these ARE the quadrature weights; for taylor
# they are the FD stencil, assembled to match .adm_combine_nll() term for term.
.admCovNodeCoefs <- function(quad) {
  if (!identical(quad$method, "taylor")) {
    w <- quad$weights
    if (is.null(w))
      stop(".admCovNodeCoefs: method '", quad$method, "' produced no weights.",
           call. = FALSE)
    return(as.numeric(w))
  }
  ns <- quad$node_signs
  K  <- nrow(ns); d <- quad$d; h <- rep_len(quad$h, d)
  cf <- numeric(K)
  ci <- which(rowSums(abs(ns)) == 0L)
  cf[ci] <- 1
  axis_only <- function(j)
    if (d > 1L) rowSums(abs(ns[, -j, drop = FALSE])) == 0L else rep(TRUE, K)

  # 2nd-order diagonal: 0.5 * sigma_jj * (f_+1 - 2 f_0 + f_-1) / h_j^2
  for (j in seq_len(d)) {
    o  <- axis_only(j)
    a  <- 0.5 * quad$Sigma_cov[j, j] / h[j]^2
    p1 <- which(ns[, j] == +1L & o); m1 <- which(ns[, j] == -1L & o)
    cf[p1] <- cf[p1] + a; cf[m1] <- cf[m1] + a; cf[ci] <- cf[ci] - 2 * a
  }
  # 2nd-order cross terms (correlated covariates)
  if (isTRUE(quad$is_correlated) && d > 1L) {
    for (i in seq_len(d - 1L)) for (j in seq(i + 1L, d)) {
      if (abs(quad$Sigma_cov[i, j]) < 1e-12) next
      a  <- 0.5 * quad$Sigma_cov[i, j] / (h[i] * h[j])
      pp <- which(ns[, i] == +1L & ns[, j] == +1L)
      pm <- which(ns[, i] == +1L & ns[, j] == -1L)
      mp <- which(ns[, i] == -1L & ns[, j] == +1L)
      mm <- which(ns[, i] == -1L & ns[, j] == -1L)
      cf[pp] <- cf[pp] + a; cf[mm] <- cf[mm] + a
      cf[pm] <- cf[pm] - a; cf[mp] <- cf[mp] - a
    }
  }
  # 4th-order diagonal
  if ((quad$order %||% 2L) >= 4L) {
    for (j in seq_len(d)) {
      o  <- axis_only(j)
      a  <- quad$Sigma_cov[j, j]^2 / 8 / h[j]^4
      p2 <- which(ns[, j] == +2L & o); m2 <- which(ns[, j] == -2L & o)
      p1 <- which(ns[, j] == +1L & o); m1 <- which(ns[, j] == -1L & o)
      cf[p2] <- cf[p2] + a;     cf[m2] <- cf[m2] + a
      cf[p1] <- cf[p1] - 4 * a; cf[m1] <- cf[m1] - 4 * a
      cf[ci] <- cf[ci] + 6 * a
    }
  }
  cf
}

# (mu, sd, link) per covariate for the node builders.
#
# gl/gh/taylor are all defined for a NORMAL covariate. A lognormal spec is
# handled on the log scale and the nodes exponentiated -- which is exact
# quadrature for the lognormal, and the only stable choice: moment-matching a
# normal to WT ~ lognormal(log 20.2, 0.45) puts the 3.5-SD Gauss-Legendre
# lower limit at a NEGATIVE weight, and WT^clwt is not defined there.
# Everything else (discrete `values`, a bare `quantile` function) has no mean
# and spread these methods can expand around, and is refused rather than
# silently moment-matched.
.admCovNodeSpec <- function(cd, method) {
  cd  <- cd[!names(cd) %in% c("rho", "Sigma")]
  out <- lapply(names(cd), function(nm) {
    sp <- cd[[nm]]
    if (!is.null(sp$mu) && !is.null(sp$sd))
      return(list(mu = sp$mu, sd = sp$sd, link = "identity"))
    if (!is.null(sp$meanlog) && !is.null(sp$sdlog))
      return(list(mu = sp$meanlog, sd = sp$sdlog, link = "exp"))
    stop("cov_method = '", method, "' needs a covariate with a mean and a ",
         "spread to expand around; covariate ", sQuote(nm), " is specified as ",
         paste(sQuote(names(sp)), collapse = "/"), ".\n",
         "  gl, gh and taylor are all defined for a normal covariate (a ",
         "lognormal is handled on the log scale).\n",
         "  Use cov_method = 'marginal', which integrates over the ",
         "distribution as given.", call. = FALSE)
  })
  stats::setNames(out, names(cd))
}

# Node values and combination coefficients for one study's covariate spec.
# Returns `values` (n_node x d, on the covariate's own scale) and `coefs`.
.admCovNodes <- function(cd, method, ctl = list()) {
  sp <- .admCovNodeSpec(cd, method)
  d  <- length(sp)
  if (d == 0L)
    stop("cov_method = '", method, "' but the study declares no covariate.",
         call. = FALSE)
  links <- vapply(sp, `[[`, character(1), "link")
  if (d > 1L && any(links == "exp"))
    stop("cov_method = '", method, "' with a lognormal covariate is supported ",
         "for a single covariate only (the nodes are built on the log scale, ",
         "and mixing scales in one product grid is not something these methods ",
         "define).", call. = FALSE)
  # admBuildQuadrature() works in (mu, sd); hand it the linked scale.
  lin <- stats::setNames(lapply(sp, function(x) list(mu = x$mu, sd = x$sd)),
                         names(sp))
  if (!is.null(cd$rho))   lin$rho   <- cd$rho
  if (!is.null(cd$Sigma)) lin$Sigma <- cd$Sigma
  quad <- admBuildQuadrature(
    lin,
    method        = method,
    n_nodes       = ctl$n_nodes       %||% 9L,
    truncation_sd = ctl$truncation_sd %||% 3.5,
    h             = ctl$h             %||% (vapply(sp, `[[`, numeric(1), "sd") / 2),
    order         = ctl$order         %||% 2L)
  vals <- quad$wt_nodes
  if (!is.matrix(vals))
    vals <- matrix(vals, ncol = 1L, dimnames = list(NULL, names(sp)))
  colnames(vals) <- names(sp)
  for (j in seq_len(d)) if (links[[j]] == "exp") vals[, j] <- exp(vals[, j])
  list(values = vals, coefs = .admCovNodeCoefs(quad), quad = quad)
}

# Expand every node-method study into its nodes. Marginal studies pass through.
#
# A node study is an ORDINARY fixed-covariate study: cov_dist is dropped, so
# nothing downstream marginalises and `.adm_cov_path` is never set. The
# combination coefficient rides in `n` (see the section header).
.admCovExpandNodes <- function(studies) {
  meth <- vapply(studies, function(s) s[["cov_method"]] %||% "marginal",
                 character(1))
  if (all(meth == "marginal")) return(studies)
  bad <- setdiff(unique(meth), c("marginal", .ADM_COV_NODE_METHODS))
  if (length(bad))
    stop("Unknown cov_method ", paste(sQuote(bad), collapse = ", "), ". ",
         "Use one of: ",
         paste(sQuote(c("marginal", .ADM_COV_NODE_METHODS)), collapse = ", "),
         ".", call. = FALSE)
  snames <- names(studies) %||% paste0("study", seq_along(studies))
  out <- list(); onames <- character(0)
  for (i in seq_along(studies)) {
    s <- studies[[i]]; m <- meth[[i]]
    if (identical(m, "marginal")) {
      out <- c(out, list(s)); onames <- c(onames, snames[[i]]); next
    }
    if (isTRUE(s$is_joint))
      stop("cov_method = '", m, "' is not supported for a joint (same-subject) ",
           "unit.", call. = FALSE)
    cd <- s[["cov_dist"]]
    if (is.null(cd))
      stop("Study ", sQuote(snames[[i]]), " sets cov_method = '", m,
           "' but declares no `cov_dist` to integrate over.", call. = FALSE)
    nd <- .admCovNodes(cd, m, s[["cov_control"]] %||% list())
    for (k in seq_len(nrow(nd$values))) {
      sk <- s
      sk[["cov"]] <- as.list(stats::setNames(nd$values[k, ],
                                             colnames(nd$values)))
      sk[["cov_dist"]]    <- NULL
      sk[["cov_rows"]]    <- NULL
      sk[["cov_method"]]  <- NULL
      sk[["cov_control"]] <- NULL
      # The combination coefficient, folded into the likelihood's multiplier.
      sk$n            <- nd$coefs[[k]] * s$n
      sk$.adm_node_of <- snames[[i]]
      sk$.adm_node_k  <- k
      sk$.adm_node_c  <- nd$coefs[[k]]
      sk$.adm_node_n  <- s$n
      out    <- c(out, list(sk))
      onames <- c(onames, sprintf("%s.node%02d", snames[[i]], k))
    }
  }
  stats::setNames(out, onames)
}

# Fold each node study's combination coefficient into its `n`.
#
# This is the route the development workflow used, and the one datagen(covariate=)
# generates: each node carries its OWN aggregate (E, V), simulated at that node's
# covariate value, and the study objective is sum_k c_k * NLL_k. A publication
# reporting summaries BY COVARIATE STRATUM gives exactly this shape.
#
# Same mechanism as .admCovExpandNodes(): `n` is a linear multiplier on the -2LL
# in every kernel, so n_k = c_k * n IS c_k * NLL_k -- for the NLL, the analytic
# gradient and both batch paths alike.
#
# `weight` is set by admBuildCovStudies() from .admCovNodeCoefs(), so Taylor's
# negative central coefficient arrives here already correct. A study with no
# weight is an ordinary study and passes through untouched.
.admCovApplyNodeWeights <- function(studies) {
  w <- vapply(studies, function(s) {
    x <- s[["weight"]]
    if (is.null(x)) 1 else as.numeric(x)[[1L]]
  }, numeric(1))
  if (all(w == 1)) return(studies)
  if (anyNA(w))
    stop("A study carries a non-numeric `weight`; it must be the node's ",
         "combination coefficient.", call. = FALSE)
  for (i in seq_along(studies)) {
    if (w[[i]] == 1) next
    studies[[i]]$.adm_node_n <- studies[[i]]$n
    studies[[i]]$.adm_node_c <- w[[i]]
    studies[[i]]$n           <- w[[i]] * studies[[i]]$n
    studies[[i]]$weight      <- NULL
  }
  studies
}
