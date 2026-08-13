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
.admCovCols <- function(mat, mod_params, cov_s) {
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

  out <- lapply(seq_along(agg_list), function(k) {
    a <- agg_list[[k]]
    node_vals <- if (multi_d)
      setNames(quad$wt_nodes[k, , drop = TRUE], nms)
    else
      setNames(quad$wt_nodes[[k]], nms)
    nm <- sprintf("%s%02d", prefix, k)
    s  <- .admNormaliseStudy(
      list(E = a$E, V = a$V, n = n, times = times, ev = ev,
           cov = node_vals, weight = (quad$weights[k]) %||% 1),
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
  nms <- unique(unlist(lapply(studies, function(s) names(s$cov))))
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
        matrix(rep(.admCovNodeVals(studies[[i]]$cov, extra), each = n_sim), n_sim, length(extra))
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
          matrix(rep(.admCovNodeVals(studies[[i]]$cov, extra_s), each = n_sim), n_sim, length(extra_s))
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
        cm[rows, ] <- matrix(rep(.admCovNodeVals(studies[[i]]$cov, extra_m), each = n_blk * n_sim),
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
          cm[rows, ] <- matrix(rep(.admCovNodeVals(studies[[i]]$cov, extra), each = n_sim),
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
