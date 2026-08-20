# -- adgh: aggregate Gauss-Hermite quadrature estimator -------------------------
# Computes population moments E[f] and Cov[f] for eta ~ N(0, Omega) by
# deterministic Gauss-Hermite quadrature over the random-effects distribution,
# then plugs them into the same aggregate MVN -2LL as adfo/admc.
#
# Structurally this is admc with a small fixed deterministic node grid in place
# of n_sim random draws, and quadrature weights in place of uniform ones.
# The objective is noise-free -> clean gradient/Hessian, fast reproducible opt.
#
# The measure is prior N(0, Omega) (prior-predictive population moments, not
# data-conditional posterior), so plain (non-adaptive) GH is exactly right.

# -- Node grid -----------------------------------------------------------------

# Probabilists' GH nodes/weights for E_{N(0,1)}[g] = sum_i w_i g(x_i).
# Golub-Welsch via symmetric tridiagonal eigendecomposition. No external deps.
# sum(w) = 1, sum(w * x^2) = 1.
# Memoised: the nodes depend on nothing but `m`, and .admTBSMoments/.admTBSMomentsD
# ask for the 81-node set on EVERY residual evaluation of a transformed endpoint --
# an eigen() of an 81x81 matrix each time, ~0.9 ms, inside the objective's inner
# loop. The cache is keyed by m, holds a handful of tiny numeric vectors, and lives
# in the package namespace (created at the top of R/zzz.R).
.adghNodes1 <- function(m) {
  if (m < 1L) stop("n_nodes must be >= 1")
  if (m == 1L) return(list(x = 0, w = 1))
  # The cache env is a package-level binding, and this function runs inside mirai
  # restart workers, where assignInNamespace() cannot ADD a binding to the locked
  # installed namespace. Degrade to recomputing rather than erroring if it is absent.
  .env <- tryCatch(get(".adm_node_env", envir = asNamespace("admixr2")),
                   error = function(e) NULL)
  if (is.null(.env)) {
    i <- seq_len(m - 1L)
    J <- matrix(0, m, m)
    J[cbind(i, i + 1L)] <- sqrt(i); J[cbind(i + 1L, i)] <- sqrt(i)
    e <- eigen(J, symmetric = TRUE)
    return(list(x = e$values, w = (e$vectors[1L, ])^2))
  }
  .key <- paste0("gh_", m)
  .hit <- tryCatch(get(.key, envir = .env, inherits = FALSE), error = function(e) NULL)
  if (!is.null(.hit)) return(.hit)
  i <- seq_len(m - 1L)
  J <- matrix(0, m, m)
  J[cbind(i, i + 1L)] <- sqrt(i)
  J[cbind(i + 1L, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE)
  out <- list(x = e$values, w = (e$vectors[1L, ])^2)
  assign(.key, out, envir = .env)
  out
}

# Tensor-product GH grid for n_eta dimensions.
# Returns X (n_node x n_eta standard-normal nodes) and W (length n_node weights).
.adghNodeGrid <- function(m, n_eta) {
  if (n_eta == 0L) return(list(X = matrix(0, 1L, 0L), W = 1))
  g <- .adghNodes1(m)
  X <- as.matrix(expand.grid(rep(list(g$x), n_eta)))
  W <- as.numeric(apply(expand.grid(rep(list(g$w), n_eta)), 1L, prod))
  dimnames(X) <- NULL
  list(X = X, W = W)
}

# -- Moments -------------------------------------------------------------------

# Population moments (E, V) for one study via GH quadrature.
# One batched .admSimulate over the node grid; weighted ML mean/cov; residual
# error added to the diagonal exactly as adfo/admc.
# Quadrature grid (eta nodes + weights) for the current Omega. Shared by the
# single and batched moment paths.
.adghGrid <- function(pars, pinfo, grid, s = NULL) {
  if (pinfo$n_eta > 0L) {
    eta <- grid$X %*% t(pars$L)
    colnames(eta) <- pinfo$eta_col_names
    g <- list(eta = eta, W = grid$W, X = grid$X, cov_rows = NULL)
  } else {
    g <- list(eta = matrix(0, 1L, 0L), W = 1, X = grid$X, cov_rows = NULL)
  }
  # General path, adgh's analogue of admc's per-row covariate draws: a PRODUCT
  # GRID over the covariate quadrature and the eta grid. Deterministic, so adgh
  # stays noise-free, and it is still ONE rxSolve -- n_cov x n_node rows rather
  # than n_node. The eta block cycles fastest, so the weights are
  # rep(W_eta, times = n_cov) * rep(W_cov, each = n_eta), which is what
  # as.numeric(outer(W_eta, W_cov)) produces column-major.
  # SHIFT: the covariate never reaches the solver. The affected eta column is
  # replaced by quantiles of u = Delta(a) + eta and the covariates are held at
  # their reference, so the solve costs n_u * (nodes for the OTHER etas) rows --
  # CONSTANT in the number of covariates, against n_node^n_eta * n_cov^p.
  sh <- if (!is.null(s)) s[[".adm_cov_shift"]] else NULL
  if (!is.null(sh) && identical(s$.adm_cov_path, "shift") && pinfo$n_eta > 0L) {
    j  <- sh$eta_idx
    D  <- .admShiftDelta(sh$spec, .admShiftStruct(pinfo, pars$struct),
                         sh$X, sh$aref)
    if (!is.null(D)) {
      D <- as.matrix(D)
      om <- sqrt(pmax(diag(as.matrix(pars$omega))[j], .Machine$double.eps))
      # HOW MANY NODES u DESERVES. u = Delta(a) + eta carries variance
      # Var(Delta) + omega^2, which is WIDER than the omega^2 that the eta
      # column carries in the product grid -- the covariate's spread has been
      # folded into this one dimension. The integrand is explored over a
      # correspondingly wider range, so resolving it as well as n_nodes resolves
      # eta needs n_nodes scaled by the ratio of standard deviations. Fixing
      # n_u at cov_nodes instead left the shift path ~10x LESS accurate than the
      # grid it replaces, which defeats the point of it.
      #
      # Cost is linear in n_u and CONSTANT in the number of covariates, so this
      # is cheap; the cap keeps a pathological covariate spread from blowing the
      # node count up without bound.
      nn0 <- as.integer(round(nrow(grid$X)^(1 / max(pinfo$n_eta, 1L))))
      # n_u MUST NOT depend on the current parameters. It used to scale with
      # sqrt((Var(Delta) + omega^2)/omega^2), so the node count changed as the
      # optimizer moved omega and the objective stepped discontinuously across
      # each switch -- 0.078 -2LL units, enough to send an FD Hessian entry from
      # -3169 to -256961. Fixed at admission instead (sh$n_u), from cov_dist,
      # which is data.
      n_u <- sh$n_u %||% min(101L, 4L * nn0)
      # A correlated Omega takes the absorption instead: the covariate becomes
      # Omega + P and the ORDINARY eta grid carries it, which is the only route
      # that keeps the off-diagonals the column substitution would drop.
      if (isTRUE(sh$absorb)) {
        ab <- .admShiftAbsorb(D, sh$W, sh$z, pars$omega, j, n_u, nn0)
        if (!is.null(ab)) {
          colnames(ab$eta) <- pinfo$eta_col_names
          cr <- matrix(rep(unlist(sh$aref[sh$cov_names]), each = nrow(ab$eta)),
                       nrow(ab$eta), length(sh$cov_names),
                       dimnames = list(NULL, sh$cov_names))
          # X is the STANDARD normal node matrix. Under the absorption
          # eta = mu + X chol(Omega + P)', so d(eta)/d(L_ab) is a Cholesky
          # differential rather than a single column -- carried in `shift`,
          # which .adghGrad applies through .admAbsorbBase.
          dv <- .admShiftAbsorbDeriv(sh$spec, .admShiftStruct(pinfo, pars$struct),
                                     sh$X, sh$aref, ab, pinfo$n_eta)
          return(list(eta = ab$eta, W = ab$W, X = ab$X, cov_rows = cr,
                      shift = if (is.null(dv)) list(degraded = TRUE) else
                        list(absorb = TRUE, Lt = ab$Lt, dmu = dv$dmu,
                             dP = dv$dP)))
        }
      }
      # A correlated Omega that did NOT absorb: condition instead of dropping.
      # The column substitution below rebuilds the unaffected etas from the
      # DIAGONAL, so it cannot carry an off-diagonal; conditioning can. See
      # .admCondShiftParts() for the construction and why it costs nothing --
      # w's law is free of eta_O, so the mixture inversion still runs once, and
      # the node count is unchanged.
      if (isTRUE(sh$cond)) {
        cp <- .admCondShiftParts(pars$omega, j)
        if (!is.null(cp)) {
          m_s <- ncol(D); Oc <- cp$O
          # rotate: w | node ~ N(Ls^-1 Delta, I), a UNIT-covariance mixture,
          # which is what the existing inversion already handles. om is 1 in
          # every direction after the rotation, so every direction below has
          # dom = 0 and only dD varies.
          Dw  <- t(solve(cp$Ls, t(D)))
          .st0 <- .admShiftStruct(pinfo, pars$struct)
          .dD0 <- .admShiftDDelta(sh$spec, .st0, sh$X, sh$aref)
          # struct thetas move Delta, hence Dw, with Ls held; omega parameters
          # move Ls, K and Lo, hence Dw through -Ls^-1 dLs Dw as well as eta
          # directly. Both are node-only: no solve, no extra rxSolve row.
          .om_d <- lapply(seq_along(pinfo$omega_par), function(rr) {
            a <- pinfo$chol_i[rr]; b <- pinfo$chol_j[rr]
            E <- matrix(0, pinfo$n_eta, pinfo$n_eta); E[a, b] <- 1
            .admCondShiftDeriv(cp, pars$omega, j,
                               E %*% t(pars$L) + pars$L %*% t(E), Dw)
          })
          .z0 <- matrix(0, nrow(Dw), m_s)
          dirs <- c(
            lapply(names(.st0), function(k)
              if (is.null(.dD0[[k]])) NULL else
                list(dD = t(solve(cp$Ls, t(as.matrix(.dD0[[k]])))),
                     dom = numeric(m_s))),
            lapply(.om_d, function(d) list(dD = d$dDw, dom = numeric(m_s))))
          if (any(vapply(dirs, is.null, logical(1)))) dirs <- NULL
          unc <- .admShiftNodesStrat(Dw, sh$W, rep(1, m_s), n_u, sh$strata,
                                     dirs)
          if (!is.null(unc)) {
            g1c  <- .adghNodes1(nn0)
            lstc <- c(list(seq_len(nrow(unc$u))),
                      lapply(Oc, function(k) seq_along(g1c$x)))
            ixc <- as.matrix(expand.grid(lstc, KEEP.OUT.ATTRS = FALSE))
            # eta_O = X_O Lo', on the ORDINARY grid, so Cov(eta_O) = Omega_OO
            # exactly -- including between two etas the covariate never touches.
            XO <- matrix(0, nrow(ixc), length(Oc))
            Wc <- unc$w[ixc[, 1L]]
            for (kk in seq_along(Oc)) {
              XO[, kk] <- g1c$x[ixc[, kk + 1L]]
              Wc <- Wc * g1c$w[ixc[, kk + 1L]]
            }
            wN   <- unc$u[ixc[, 1L], , drop = FALSE]
            eO   <- if (length(Oc)) XO %*% t(cp$Lo) else NULL
            etac <- matrix(0, nrow(ixc), pinfo$n_eta)
            if (length(Oc)) etac[, Oc] <- eO
            etac[, j] <- wN %*% t(cp$Ls) +
              (if (length(Oc)) eO %*% t(cp$K) else 0)
            colnames(etac) <- pinfo$eta_col_names
            crc <- matrix(rep(unlist(sh$aref[sh$cov_names]), each = nrow(etac)),
                          nrow(etac), length(sh$cov_names),
                          dimnames = list(NULL, sh$cov_names))
            shc <- list(degraded = TRUE)
            if (!is.null(dirs) && !is.null(unc$du)) {
              duN <- unc$du[ixc[, 1L], , , drop = FALSE]
              n_th <- length(.st0)
              # d(eta) per direction, assembled once. Struct directions move
              # only w; omega directions move Lo, K and Ls as well, so every
              # column responds and the contraction is taken in full.
              mkE <- function(d, om_i) {
                dE <- matrix(0, nrow(etac), pinfo$n_eta)
                dw <- duN[, , d, drop = FALSE]; dim(dw) <- dim(duN)[1:2]
                dS <- dw %*% t(cp$Ls)
                if (!is.null(om_i)) {
                  o <- .om_d[[om_i]]
                  dS <- dS + wN %*% t(o$dLs)
                  if (length(Oc)) {
                    dEO <- XO %*% t(o$dLo)
                    dE[, Oc] <- dEO
                    dS <- dS + dEO %*% t(cp$K) + eO %*% t(o$dK)
                  }
                }
                dE[, j] <- dS
                dE
              }
              shc <- list(cond = TRUE, th_names = names(.st0),
                          dEta_th = lapply(seq_len(n_th), function(k)
                            mkE(k, NULL)),
                          dEta_om = lapply(seq_along(pinfo$omega_par),
                            function(rr) mkE(n_th + rr, rr)))
            }
            # X is read ONLY by the omega chain, which assumes eta = X L'. That
            # does not hold here, so it is zero and dEta_om carries the whole
            # omega path instead.
            return(list(eta = etac, W = Wc / sum(Wc),
                        X = matrix(0, nrow(etac), pinfo$n_eta),
                        cov_rows = crc, shift = shc))
          }
        }
      }
      # A VECTOR shift that did not absorb takes the Rosenblatt recursion, and
      # carries its derivatives with it: every u_k moves both because Delta does
      # and because the posterior weights conditioning level k do. Directions
      # are the structural thetas (through d(Delta)/d(theta), two vectorised
      # evaluations each and no solve) and the shifted etas' own scales.
      .mdu <- NULL
      if (ncol(D) > 1L) {
        .stn <- pinfo$struct_names
        .st0 <- .admShiftStruct(pinfo, pars$struct)
        .dD0 <- .admShiftDDelta(sh$spec, .st0, sh$X, sh$aref)
        .dirs <- c(
          lapply(names(.st0), function(k)
            if (is.null(.dD0[[k]])) NULL else
              list(dD = .dD0[[k]], dom = numeric(ncol(D)))),
          lapply(seq_len(ncol(D)), function(a) {
            e <- numeric(ncol(D)); e[a] <- 1
            list(dD = matrix(0, nrow(D), ncol(D)), dom = e)
          }))
        if (any(vapply(.dirs, is.null, logical(1)))) .dirs <- NULL
        un <- .admShiftNodesStrat(D, sh$W, om, n_u, sh$strata, .dirs)
        if (!is.null(un) && !is.null(.dirs))
          .mdu <- list(n_th = length(.st0), th_names = names(.st0))
      } else if (!is.null(sh$strata)) {
        # A DISCRETE covariate: condition on its exactly-enumerated levels, so
        # each cell is the mild sub-mixture the quadrature resolves well. The
        # derivatives come from the SAME construction -- .admShiftDu answers for
        # a single mixture and would disagree with a stratified node set, which
        # is the objective-and-gradient split this file exists to avoid.
        .st1 <- .admShiftStruct(pinfo, pars$struct)
        .dD1 <- .admShiftDDelta(sh$spec, .st1, sh$X, sh$aref)
        .d1  <- c(lapply(names(.st1), function(k)
                    if (is.null(.dD1[[k]])) NULL else
                      list(dD = as.matrix(.dD1[[k]]), dom = 0)),
                  list(list(dD = matrix(0, nrow(D), 1L), dom = 1)))
        if (any(vapply(.d1, is.null, logical(1)))) .d1 <- NULL
        un <- .admShiftNodesStrat(D, sh$W, om, n_u, sh$strata, .d1)
        if (!is.null(un) && !is.null(.d1) && !is.null(un$du))
          .mdu <- list(n_th = length(.st1), th_names = names(.st1))
      } else {
        un0 <- .admShiftNodes(D[, 1L], sh$W, om[1L], n_u, z = sh$z)
        un  <- if (is.null(un0)) NULL else
          list(u = matrix(un0$u, ncol = 1L), w = un0$w)
      }
      # n_nodes per eta, recovered from the grid: nrow = n_nodes^n_eta. round(),
      # not a bare fractional power -- 343^(1/3) is 6.999999999999999.
      g1 <- .adghNodes1(nn0)
      other <- setdiff(seq_len(pinfo$n_eta), j)
      # the shifted columns move together (one index over the u node SET), the
      # remaining etas keep their own product grid
      lst <- c(list(seq_len(nrow(un$u))),
               lapply(other, function(k) seq_along(g1$x)))
      ix <- as.matrix(expand.grid(lst, KEEP.OUT.ATTRS = FALSE))
      eta <- matrix(0, nrow(ix), pinfo$n_eta)
      eta[, j] <- un$u[ix[, 1L], , drop = FALSE]
      W <- un$w[ix[, 1L]]
      for (kk in seq_along(other)) {
        k <- other[kk]
        eta[, k] <- sqrt(pars$omega[k, k]) * g1$x[ix[, kk + 1L]]
        W <- W * g1$w[ix[, kk + 1L]] }
      colnames(eta) <- pinfo$eta_col_names
      cr <- matrix(rep(unlist(sh$aref[sh$cov_names]), each = nrow(eta)),
                   nrow(eta), length(sh$cov_names),
                   dimnames = list(NULL, sh$cov_names))
      # THE OMEGA CHAIN NEEDS NO SPECIAL CASE. It forms d(f)/d(L_ab) as
      # Jl[[a]] * X[, b], and for the affected column eta_j = u with
      # du/dp = (du/domega)(domega/dp) = (du/domega) * L_jj/2 -- exactly the
      # shape the existing loop applies. So putting du/domega into X[, j] makes
      # that loop correct as written; the OTHER columns keep their standard
      # normal nodes, which is what they are.
      Xz <- matrix(0, nrow(eta), pinfo$n_eta)
      for (kk in seq_along(other))
        Xz[, other[kk]] <- g1$x[ix[, kk + 1L]]
      shinfo <- NULL
      # A STRATIFIED node set joins the vector shift in needing its derivatives
      # from the same construction: .admShiftDu below answers for a SINGLE
      # mixture, so letting the m == 1 case fall through to it would pair
      # stratified nodes with unstratified derivatives.
      if ((ncol(D) > 1L || !is.null(sh$strata)) &&
          (is.null(.mdu) || is.null(un$du))) {
        # The node derivatives could not be built. Say so, rather than return a
        # grid with no `shift`: the gradient would then simply omit this study's
        # shift chain -- finite, plausible and a direction the objective does
        # not follow, which is the failure this file keeps meeting.
        shinfo <- list(degraded = TRUE)
      } else if (!is.null(.mdu) && !is.null(un$du)) {
        # Every shifted coordinate responds to every direction, so the omega
        # chain cannot be folded into an X column the way the scalar case can:
        # d(om_1) moves u_2 through the posterior weights. .adghGrad forms the
        # full sum instead.
        shinfo <- list(multi = TRUE, eta_idx = j, th_names = .mdu$th_names,
                       du = un$du[ix[, 1L], , , drop = FALSE],
                       n_th = .mdu$n_th)
      } else if (length(j) == 1L) {
        du <- .admShiftDu(sh$spec, .admShiftStruct(pinfo, pars$struct), sh$X,
                          sh$aref, D, sh$W, om, un$u[, 1L], z = sh$z)
        Xz[, j] <- du$du_domega[ix[, 1L]]
        shinfo <- list(eta_idx = j,
                       du_dtheta = du$du_dtheta[ix[, 1L], , drop = FALSE])
      }
      return(list(eta = eta, W = W / sum(W), X = Xz, cov_rows = cr,
                  shift = shinfo))
    }
  }
  # JOINT COLLAPSE: one design over the etas AND the covariates together, where
  # they reach the model through the same directions. It replaces the eta grid
  # as well as the covariate design, so it returns before either is built.
  #
  # X is the node matrix the omega chain rule differentiates. eta = X L' holds
  # here exactly as it does for the ordinary grid -- the joint preimage's eta
  # block IS that matrix -- so .adghGrad needs no branch of its own. What it
  # does not carry is the rotation's own dependence on Omega; that term is the
  # quadrature re-choosing itself within the same column space, and vanishes to
  # the accuracy the design is verified to.
  .jc <- if (!is.null(s)) s[[".adm_cov_joint"]] else NULL
  if (!is.null(.jc) && !identical(s$.adm_cov_path, "shift")) {
    jd <- .admJointDesign(.jc, .admShiftStruct(pinfo, pars$struct), pars$L)
    if (!is.null(jd))
      return(list(eta = jd$eta, W = jd$W, X = jd$X, cov_rows = jd$cov_rows))
  }
  if (!is.null(s) && !identical(s$.adm_cov_path, "shift") &&
      !is.null(s[["cov_dist"]])) {
    nq <- max(nrow(g$eta), 1L)
    # cov_integration = "taylor": 1 + 2p design points in place of the product
    # grid's n_nodes^p, with SIGNED combination weights. Every expansion below
    # uses the same stride as the quadrature branch, so everything downstream --
    # the omega chain's X, the per-row covariate columns, the params frame --
    # is laid out identically and only the weights and the centring differ.
    .taylor <- identical(pinfo$cov_integration %||% "quadrature", "taylor")
    if (.taylor && isTRUE(s$is_joint))
      stop("admixr2: covariate marginalisation is not supported for a JOINT ",
           "(same-subject, multi-output) unit. The shared-eta joint solve has ",
           "no per-row covariate path, so this would silently solve at the ",
           "covariate mean.", call. = FALSE)
    # The design is cached on the study by .admCheckCovariates -- it depends
    # only on `cov_dist` and `cov_taylor_h`, both data, and rebuilding it here
    # costs a Sobol pass through the joint sampler on EVERY objective call.
    # The %||% keeps a hand-built study (Tier-1 mocks, direct .adghMoments
    # calls) working, at the old cost.
    cg <- if (.taylor) s[[".adm_cov_taylor"]] %||%
                       .admCovTaylorDesign(s[["cov_dist"]],
                                           pinfo$cov_taylor_h %||% 1)
          # The COLLAPSED design when the covariates reach the model through a
          # single scalar: the same integral in the dimension it actually has,
          # so this is not an approximation the grid would beat. Cached at
          # admission (.admCovCollapse costs a probe, no solves); the %||% keeps
          # a hand-built study working at the old cost.
          # RE-AIMED at the current thetas, not read from admission. The
          # rotation depends on the covariate coefficients, which are estimated,
          # so a design cached at the starting values integrates over the wrong
          # line in latent space as soon as the optimizer moves them -- measured
          # at 53 to 163 -2LL units for a 0.1 move in one coefficient. This is
          # the same thing the shift branch above does with .admShiftDelta.
          else         .admCovRefresh(s[[".adm_cov_collapse"]],
                                      .admShiftStruct(pinfo, pars$struct)) %||%
                       .admCovGrid(s[["cov_dist"]], pinfo$cov_nodes %||% 7L)
    nc <- nrow(cg$X)
    g$eta      <- g$eta[rep(seq_len(nq), times = nc), , drop = FALSE]
    colnames(g$eta) <- pinfo$eta_col_names
    # The node matrix has to be expanded with the SAME stride: the omega chain
    # differentiates eta = X %*% t(L) row by row, so a gradient using the
    # unexpanded X against an expanded eta would be silently misaligned.
    g$X        <- g$X[rep(seq_len(nq), times = nc), , drop = FALSE]
    g$cov_rows <- cg$X[rep(seq_len(nc), each = nq), , drop = FALSE]
    if (.taylor) {
      g$taylor <- .admCovTaylorRows(cg, g$W, nq)
      g$W      <- as.numeric(outer(g$W, cg$c))
    } else {
      g$W      <- as.numeric(outer(g$W, cg$W))
    }
  }
  g
}

# Structural moments (before any residual) from an already-solved node matrix.
#
# The quadrature path is one weighted mean and one weighted crossproduct about
# it. The Taylor path is THE SAME TWO EXPRESSIONS -- with signed weights
# c_k * w_q, and each design point's rows centred at their OWN conditional mean
# rather than at the pooled one -- plus the rank-p term:
#
#   crossprod(W, cp)                 = sum_k c_k E_k               = E_marg
#   crossprod(cpc, W * cpc)          = sum_k c_k Vc_k              (block-centred)
#   crossprod(dE, var * dE)          = sum_j v_j g'_j g'_j'        = Cov_a(g(a))
#
# `dE` is returned because the gradient needs it: the rank-p term is quadratic
# in the conditional means, so it is the one part of V whose derivative is not
# already carried by the weighted-crossproduct contraction.
.adghStructMoments <- function(cp, W, tay = NULL) {
  mu <- as.numeric(crossprod(W, cp))
  if (is.null(tay)) {
    cpc <- sweep(cp, 2L, mu)
    return(list(mu = mu, cpc = cpc, V = crossprod(cpc, W * cpc), dE = NULL))
  }
  cpc <- cp
  for (idx in tay$rows) {
    ek <- as.numeric(crossprod(tay$w, cp[idx, , drop = FALSE]))
    cpc[idx, ] <- sweep(cp[idx, , drop = FALSE], 2L, ek)
  }
  dE <- crossprod(tay$Dw, cp)                    # d x n_t, g'_j(mu_a)
  list(mu = mu, cpc = cpc,
       V = crossprod(cpc, W * cpc) + crossprod(dE, tay$var * dE),
       dE = dE)
}

# Attach the grid's per-row covariate values to a study, so .admSimulate writes
# them into the params frame. Returns the study untouched when there are none.
.adghStudyCov <- function(study, g) {
  if (!is.null(g$cov_rows)) study$cov_rows <- g$cov_rows
  study
}

# Weighted moments + residual error from an already-solved quadrature matrix.
# Split out of .adghMoments so the solve and the assembly can be driven
# independently: the assembly depends on sigma, but the SOLVE does not (sigma is
# zeroed into it and re-added analytically here), so a set of configurations
# that share a solve can each be assembled cheaply.
.adghMomentsFromCp <- function(cp, W, pars, pinfo, out_var, times = NULL,
                               tay = NULL) {
  sm  <- .adghStructMoments(cp, W, tay)
  mu  <- sm$mu
  V   <- sm$V

  # Restrict residual error to this output's sigma(s) (no-op single-output).
  arr <- .admUnitResidRows(pinfo, out_var, pars$sigma_var, length(mu),
                           phi = attr(cp, "phi"))   # beta precision (SOLVED)
  m <- .admResidMoments(mu, diag(V), arr, V, times)
  list(E = m$mu, V = m$V)
}

.adghMoments <- function(pars, pinfo, study, rxMod, out_var, grid, cores) {
  g  <- .adghGrid(pars, pinfo, grid, study)
  study <- .adghStudyCov(study, g)
  pm <- .admMakeParamsList(nrow(g$eta), pinfo, 1L)[[1L]]
  cp <- .admSimulate(rxMod, pars$struct, pinfo$sigma_names, g$eta, study,
                     out_var, pm, cores, pinfo$nDisplayProgress,
                             pinfo$sigdig)
  .adghMomentsFromCp(cp, g$W, pars, pinfo, out_var, study$times, g$taylor)
}

# Moments for a SET of structural-theta configurations in ONE rxSolve.
# The node grid and Omega are identical across configurations -- only the
# structural thetas move -- so the n_cfg quadrature solves stack into one call
# of n_cfg * n_node subjects instead of n_cfg calls of n_node.
.adghMomentsBatch <- function(struct_mat, pars, pinfo, study, rxMod, out_var, grid, cores) {
  g     <- .adghGrid(pars, pinfo, grid, study)
  study <- .adghStudyCov(study, g)
  Q     <- nrow(g$eta)
  n_cfg <- nrow(struct_mat)

  sm_big  <- struct_mat[rep(seq_len(n_cfg), each = Q), , drop = FALSE]
  # The frame stacks n_cfg blocks of Q rows, so the grid's covariate rows have to
  # be tiled to match -- every configuration must see the SAME covariate nodes,
  # or the struct-theta differences stop comparing like with like.
  if (!is.null(g$cov_rows))
    study$cov_rows <- g$cov_rows[rep(seq_len(Q), times = n_cfg), , drop = FALSE]
  # ONE node grid for every configuration. It used to be rebuilt per
  # configuration whenever the study declared a covariate distribution, because
  # the retired "collapse" path made Omega itself a function of the covariate
  # COEFFICIENT -- a structural theta, so it moved between configurations. With
  # collapse gone the grid depends only on Omega, which does not move here.
  #
  # The SHIFT path does still move its grid with the structural thetas (Delta,
  # and hence the u nodes AND their weights, are functions of the covariate
  # coefficient), and a per-configuration weight vector is something this
  # function cannot express -- it assembles every block against the single
  # `g$W`. .adghGrad therefore never routes a shift study here; see the guard on
  # its FD block.
  eta_big <- g$eta[rep(seq_len(Q), times = n_cfg), , drop = FALSE]
  colnames(eta_big) <- colnames(g$eta)

  pm     <- .admMakeParamsList(n_cfg * Q, pinfo, 1L)[[1L]]
  cp_all <- .admSimulateRows(rxMod, sm_big, pinfo$sigma_names, eta_big, study,
                             out_var, pm, cores, pinfo$nDisplayProgress,
                             pinfo$sigdig)

  # beta: phi = b1 + b2 comes back as one row per SOLVED row, and subsetting a
  # matrix drops attributes -- so each configuration's block has to carry its own
  # forward. phi is eta-independent (checked at parse) but not theta-independent,
  # and each block holds a different theta vector, so it is the block's own first
  # row that is representative, not the whole matrix's.
  .phi_all <- attr(cp_all, "phi")
  lapply(seq_len(n_cfg), function(k) {
    .cp <- cp_all[(k - 1L) * Q + seq_len(Q), , drop = FALSE]
    if (!is.null(.phi_all)) attr(.cp, "phi") <- .phi_all[(k - 1L) * Q + 1L, ]
    .adghMomentsFromCp(.cp, g$W, pars, pinfo, out_var, study$times, g$taylor)
  })
}

# GH-quadrature joint moments for a same-subject unit: one shared-eta node grid
# gives every output; stacked weighted mean/cov (full, cross blocks included) +
# per-output residual. Returns list(E = mu_sigma, V).
.adghMomentsJoint <- function(pars, pinfo, unit, rxMod, grid, cores) {
  n_eta <- pinfo$n_eta
  if (n_eta > 0L) {
    if (!is.null(unit[["cov_dist"]]))
      stop("admixr2: covariate marginalisation is not supported for a JOINT ",
           "(same-subject, multi-output) unit. The shared-eta joint solve has ",
           "no per-row covariate path, so this would silently solve at the ",
           "covariate mean.", call. = FALSE)
    eta <- grid$X %*% t(pars$L)
    colnames(eta) <- pinfo$eta_col_names; W <- grid$W
  } else { eta <- matrix(0, 1L, 0L); W <- 1 }
  pm <- .admMakeParamsList(nrow(eta), pinfo, 1L)[[1L]]
  cp <- .admSimulateJoint(rxMod, pars$struct, pinfo$sigma_names, eta, unit, pm, cores,
                          pinfo$nDisplayProgress, pinfo$sigdig)
  mu  <- as.numeric(crossprod(W, cp))
  cpc <- sweep(cp, 2L, mu)
  V   <- crossprod(cpc, W * cpc)
  jr  <- .admJointResidual(mu, V, unit, pinfo, pars$sigma_var)
  list(E = jr$mu, V = jr$V)
}

# -- NLL -----------------------------------------------------------------------

#' @noRd
.adghNLL <- function(p, pinfo, studies, rxMod, out_var, grid, cores) {
  pars <- tryCatch(.admUnpack(p, pinfo), error = function(e) NULL)
  if (is.null(pars)) return(Inf)
  # Reject non-finite parameters before the solve; see .admParsFinite().
  if (!.admParsFinite(pars, pinfo)) return(Inf)
  total <- 0
  for (s in studies) {
    if (isTRUE(s$is_joint)) {
      m   <- .adghMomentsJoint(pars, pinfo, s, rxMod, grid, cores)
      nll <- nll_cov_cpp(s$E, s$V, m$E, m$V, s$n)
    } else {
      m <- .adghMoments(pars, pinfo, s, rxMod, s$output %||% out_var, grid, cores)
      nll <- if (identical(s$method, "var"))
        nll_var_cpp(s$E, s$v_diag, m$E, diag(m$V), s$n)
      else
        nll_cov_cpp(s$E, s$V, m$E, m$V, s$n)
    }
    if (!is.finite(nll)) return(Inf)
    total <- total + nll
  }
  total
}

# -- Analytic gradient ---------------------------------------------------------

# Analytic gradient of the GH NLL w.r.t. optimizer vector p.
# One batched sensitivity solve per study over the node grid; closed-form
# contractions for struct thetas, omega Cholesky, sigma (add/prop/lnorm).
# Unpaired struct thetas: central FD of .adghNLL (like admc).
#
# Var-method studies use a diagonal derivative path; cov-method uses full B.
# For lnorm sigma: Jl scaled by exp(sv/2) for mean path; analytical sigma grad.
.adghGradNLL <- function(p, pinfo, studies, sensModel, rxMod, out_var, grid, cores,
                       grad_h = 1e-4) {
  pars  <- .admUnpack(p, pinfo)
  # Non-finite parameters never reach rxSolve -- see the note in .adghNLL. The
  # gradient unpacks `p` itself, so the NLL's guard does not cover this entry.
  if (!.admParsFinite(pars, pinfo))
    return(list(grad = stats::setNames(rep(NA_real_, length(p)), names(p)),
                nll = Inf))
  L     <- pars$L
  n_eta <- pinfo$n_eta
  n_s   <- length(pinfo$struct_names)
  n_e   <- length(pinfo$sigma_names)
  X     <- grid$X
  W     <- grid$W
  grad  <- numeric(length(p)); names(grad) <- names(p)

  # (#5) The moments this gradient is built from are exactly the moments the NLL
  # needs, so returning the NLL alongside costs nothing and lets the driver skip
  # a whole second solve per iterate (see .adghFusedFns). Returned as a list
  # rather than an attribute on the gradient: an attribute travels through
  # unname() and into expect_equal(), where it broke a gradient-vs-FD comparison
  # that had nothing to do with the NLL. Callers wanting just the gradient use
  # the .adghGrad wrapper below. nll = NULL on the FD-fallback returns -- those
  # never form these moments.
  nll_total <- 0

  # Which struct thetas are unpaired (no mu-referencing eta)?
  # struct_has_eta is struct-indexed (length n_s); struct_eta_idx is eta-indexed
  # (length n_eta) so is.na() on it never flags unpaired struct thetas.
  unpaired_k <- if (!is.null(pinfo$struct_has_eta))
    which(!pinfo$struct_has_eta) else integer(0)

  # An unpaired theta shifts the quadrature moments exactly like an eta does, so
  # given d(pred)/d(theta) from the augmented sens model it goes through the SAME
  # contrib() + sigma-V-coupling the paired thetas use -- no FD, no step size.
  # Accumulated separately: if ANY study fails to return theta columns, the whole
  # theta gradient falls back to the FD block below (mixing the two across
  # studies would double-count the studies already accumulated here).
  theta_sens_ok <- length(unpaired_k) > 0L
  g_theta       <- numeric(length(p))

  for (s in studies) {
    # Per study, through the SAME helper the objective uses. Deriving the grid
    # here from pars$L instead is what made adgh's analytical gradient blind to
    # the covariate product grid -- it differentiated a different function than
    # .adghNLL evaluated. X, W and eta must all come from one place.
    # A VECTOR shift (a covariate on more than one mu-referenced parameter) moves
    # the later coordinates' nodes through the Rosenblatt posterior weights as
    # well, a second chain .admShiftDu does not carry. Finite-difference the
    # objective instead of silently using the scalar chain -- still far cheaper
    # than the product grid, because the objective is what got cheap.
    .gS <- .adghGrid(pars, pinfo, grid, s)
    # A shift whose node derivatives could not be built degrades the WHOLE
    # gradient to finite differences, the same way a failed sensitivity solve
    # does. Continuing would drop this study's shift chain silently.
    if (isTRUE(.gS$shift$degraded))
      return(list(grad = .adghFDGrad(p, pinfo, studies, rxMod, out_var, grid,
                                     cores, grad_h), nll = NULL))
    X   <- .gS$X
    W   <- .gS$W
    ty  <- .gS$taylor          # NULL unless cov_integration = "taylor"
    s   <- .adghStudyCov(s, .gS)
    eta <- .gS$eta
    colnames(eta) <- pinfo$eta_col_names

    # --- Joint (same-subject) analytical quadrature gradient -----------------
    # Stacked weighted moments over all outputs (shared-eta node grid); paired
    # struct + omega + sigma analytical on the joint covariance, per output rows.
    # Unpaired struct thetas are ALSO analytical here (via js$dtheta_list, same path
    # as the paired ones); they fall back to the FD block only when the augmented
    # sens columns are unavailable (theta_sens_ok = FALSE).
    if (isTRUE(s$is_joint)) {
      js <- .admSimulateJointSens(sensModel, pars$struct, pinfo$sigma_names, eta, s, cores,
                                  pinfo$nDisplayProgress, pars$sigma_var, pinfo$sigdig)
      # A failed sens solve used to `next`, which SILENTLY DROPPED this study's
      # entire contribution -- the optimizer then walked a gradient that was
      # missing whole studies, with no error and no warning. Degrade the whole
      # gradient to finite differences instead (what admc/adfo already do).
      if (is.null(js))
        return(list(grad = .adghFDGrad(p, pinfo, studies, rxMod, out_var, grid,
                                       cores, grad_h), nll = NULL))
      f  <- js$cp_mat; Jl <- js$dpred_list
      mu  <- as.numeric(crossprod(W, f))
      cpc <- sweep(f, 2L, mu)
      # per-row residual: mean scaling (lnorm), then the residual-adjusted moments
      arr    <- .admResidRows(pinfo, .admRowOutput(s, s$n_total), pars$sigma_var, s$n_total)
      V_str  <- crossprod(cpc, W * cpc)
      var_f  <- diag(V_str)                 # Var_eta(f), pre-residual
      jr <- .admJointResidual(mu, V_str, s, pinfo, pars$sigma_var)
      mu_sigma <- jr$mu; V <- jr$V
      nll_total <- nll_total + nll_cov_cpp(s$E, s$V, mu_sigma, V, s$n)
      r  <- as.numeric(s$E) - mu_sigma
      G  <- tryCatch(chol2inv(chol(V)),
                     error = function(e) tryCatch(solve(V), error = function(e2) NULL))
      # A singular predicted V (tiny omega, near-duplicate observation times) used
      # to `next` -- returning a gradient that silently OMITTED this study. It is
      # finite and looks valid, so nloptr steps along a direction that is not a
      # descent direction for the true objective. Degrade to FD, as below.
      if (is.null(G))
        return(list(grad = .adghFDGrad(p, pinfo, studies, rxMod, out_var, grid,
                                       cores, grad_h), nll = NULL))
      B    <- s$n * (G - G %*% (s$V + tcrossprod(r)) %*% G)
      dNLL_dmu_sig <- as.numeric(-2 * s$n * (G %*% r))
      # Bt is chained to the STRUCTURAL covariance (ms_i*ms_j off-diagonal,
      # dv_dv0 on it), so contrib_j takes the RAW Jacobian -- see the
      # single-output branch below for why pre-scaling the caller's gmat is wrong.
      Bdiag <- diag(B)
      # ONE moment tail for this unit: .admResidDeriv, the V_pred -> V_struct
      # chain, the TBS mean-from-covariance diagonal fold, the sigma contraction.
      ch     <- .admResidChain(mu, var_f, arr, pinfo, dNLL_dmu_sig, Bdiag, B,
                               V_str, .admRowTimes(s, length(mu)))
      dres   <- ch$deriv
      ls_vec <- ch$dmu_df              # d(mu_pred)/df -- see the single-output branch
      Bsj    <- ch$dV
      Bt <- cpc %*% Bsj
      contrib_j <- function(graw) {              # graw = d(f)/dpsi rows, RAW
        dmu <- as.numeric(crossprod(W, graw))
        sum(dNLL_dmu_sig * dmu * ls_vec) + 2 * sum(W * rowSums(graw * Bt))
      }
      # V-path of the mean: the residual variance itself depends on mu, so a
      # parameter that moves mu also moves diag(V). dv_df = d(var)/d(mu).
      # ms = m'(f) (TBS) also reaches V's off-diagonal -- see the single-output branch.
      ms_off_j <- numeric(length(mu))
      # na.rm: dms_df is NaN when f crosses the transform bound (a line-search trial
      # point); any(NaN != 0) is NA and `if (NA)` throws instead of the optimizer
      # backing off the already-Inf objective. Bit-identical when dms_df is finite.
      if (!is.null(dres$dms_df) && any(dres$dms_df != 0, na.rm = TRUE))
        ms_off_j <- 2 * dres$dms_df * .admMsOffDiag(B, V_str, dres$ms)
      sig_V_extra <- function(dmu_raw)            # dmu_raw = d(mu)/dpsi (pre-lnorm)
        sum((Bdiag * dres$dv_df + ms_off_j) * dmu_raw)
      # paired struct thetas
      for (k in seq_len(n_s)) {
        ei <- which(pinfo$struct_eta_idx == k); if (length(ei) == 0L) next; ei <- ei[[1L]]
        dmu_raw <- as.numeric(crossprod(W, Jl[[ei]]))
        grad[k] <- grad[k] + contrib_j(Jl[[ei]]) + sig_V_extra(dmu_raw)
      }
      # unpaired struct thetas (augmented sens model): same path as the paired ones
      if (length(unpaired_k) > 0L) {
        if (is.null(js$dtheta_list)) {
          theta_sens_ok <- FALSE
        } else for (k in unpaired_k) {
          Dt      <- js$dtheta_list[[pinfo$struct_names[k]]]
          dmu_raw <- as.numeric(crossprod(W, Dt))
          g_theta[k] <- g_theta[k] + contrib_j(Dt) + sig_V_extra(dmu_raw)
        }
      }
      # omega Cholesky
      if (n_eta > 0L) for (rr in seq_along(pinfo$omega_par)) {
        i <- pinfo$chol_i[rr]; j <- pinfo$chol_j[rr]
        base    <- Jl[[i]] * X[, j]
        dmu_raw <- as.numeric(crossprod(W, base))
        dL  <- contrib_j(base) + sig_V_extra(dmu_raw)
        pos <- n_s + n_e + rr
        grad[pos] <- grad[pos] + if (pinfo$chol_diag[rr]) dL * L[i, i] / 2 else dL
      }
      # sigma (each row's own endpoint; other endpoints' derivatives are zero)
      grad[n_s + seq_len(n_e)] <- grad[n_s + seq_len(n_e)] +
        ch$sigma_grad()
      next
    }

    ov <- s$output %||% out_var

    res <- .admSimulateSens(sensModel, pars$struct, pinfo$sigma_names, eta, s, cores,
                            pinfo$nDisplayProgress, pars$sigma_var, pinfo$sigdig)
    # .admSimulateSens returns NULL when the solve fails. `next` skipped the study
    # -- i.e. returned a gradient that silently omitted it. Degrade the whole
    # gradient to finite differences instead (what admc/adfo already do).
    if (is.null(res))
      return(list(grad = .adghFDGrad(p, pinfo, studies, rxMod, out_var, grid,
                                     cores, grad_h), nll = NULL))
    f   <- res$cp_mat     # Q x n_t
    Jl  <- res$dpred_list # list n_eta of Q x n_t

    sm  <- .adghStructMoments(f, W, ty)
    mu  <- sm$mu
    cpc <- sm$cpc
    V   <- sm$V
    cov_f <- V                            # STRUCTURAL Cov_eta(f), before any residual

    # Residual error (and its lnorm scaling of the mean) -- this output only
    arr   <- .admResidRows(pinfo, ov, pars$sigma_var, length(mu))
    var_f <- diag(V)                      # Var_eta(f), pre-residual
    pm <- .admResidMoments(mu, var_f, arr, cov_f, s$times)
    V  <- pm$V
    mu_sigma <- pm$mu

    # The sens Jacobians give d(f)/d(psi). `contrib()` below takes them RAW and
    # applies the two chains itself:
    #   mean : d(mu_sigma)/dpsi = ms * d(E[f])/dpsi
    #   cov  : d(NLL)/d(V_struct) = d(NLL)/d(V_pred) o .admResidVChain()
    # Previously the caller pre-multiplied the Jacobian by ms and fed that to BOTH
    # terms, which gave the covariance path a single factor of ms where it needs
    # one per index (ms_i*ms_j) and no dv_dv0 on the diagonal at all.
    dres        <- .admResidDeriv(mu, var_f, arr, pinfo)
    # d(mu_pred)/d(f), which is ap$ms for every form EXCEPT TBS -- there the mean
    # carries a curvature term of its own. ap$ms stays the COVARIANCE scale.
    lnorm_scale <- dres$dmu_df

    r <- as.numeric(s$E) - mu_sigma

    is_var <- identical(s$method, "var")

    nll_total <- nll_total + if (is_var)
      nll_var_cpp(s$E, s$v_diag, mu_sigma, diag(V), s$n)
    else
      nll_cov_cpp(s$E, s$V, mu_sigma, V, s$n)

    if (is_var) {
      # ------ Var method: diagonal derivative path ----------------------------
      V_diag        <- diag(V)
      dNLL_dmu_sig  <- as.numeric(-2 * s$n * r / V_diag)  # d(NLL)/d(mu_sigma)
      dNLL_dV_diag  <- s$n * (1/V_diag - (s$v_diag + r^2) / V_diag^2)
      # + the mean's dependence on Var_eta(f) (TBS only; see .admResidVChain).
      # The moment tail, once: the V_pred -> V_struct chain plus the TBS
      # mean-from-covariance diagonal. NULL dNLL_dV/cov_f: diagonal path.
      ch            <- .admResidChain(mu, var_f, arr, pinfo, dNLL_dmu_sig,
                                      dNLL_dV_diag, NULL, NULL, s$times,
                                      deriv = dres)
      dNLL_dV_dg_s  <- ch$dV_diag                        # -> d(NLL)/d(var_f)

      contrib <- function(graw) {
        # graw: Q x n_t, RAW derivative of the structural f w.r.t. psi
        dmu     <- as.numeric(crossprod(W, graw))
        dV_diag <- 2 * colSums(W * cpc * graw)
        # Taylor: V also carries sum_j v_j g'_j g'_j', which is QUADRATIC in the
        # conditional means and so contributes a term the weighted crossproduct
        # above does not. d(g'_j)/dpsi is the same central difference applied to
        # the sensitivity columns.
        if (!is.null(ty))
          dV_diag <- dV_diag +
            2 * colSums(ty$var * sm$dE * crossprod(ty$Dw, graw))
        sum(dNLL_dmu_sig * dmu * lnorm_scale) + sum(dNLL_dV_dg_s * dV_diag)
      }

    } else {
      # ------ Cov method: full-matrix derivative path -------------------------
      G    <- tryCatch(chol2inv(chol(V)),
                       error = function(e) tryCatch(solve(V), error = function(e2) NULL))
      # Singular predicted V -- see the joint branch above. `next` silently dropped
      # the study from the gradient; degrade the whole gradient to FD instead.
      if (is.null(G))
        return(list(grad = .adghFDGrad(p, pinfo, studies, rxMod, out_var, grid,
                                       cores, grad_h), nll = NULL))
      Vhat      <- s$V + tcrossprod(r)
      B         <- s$n * (G - G %*% Vhat %*% G)
      dNLL_dmu_sig <- as.numeric(-2 * s$n * (G %*% r))  # d(NLL)/d(mu_sigma)
      Bdiag     <- diag(B)
      ch        <- .admResidChain(mu, var_f, arr, pinfo, dNLL_dmu_sig, Bdiag, B,
                                  cov_f, s$times, deriv = dres)
      Bs        <- ch$dV                  # mean-from-covariance fold included
      Bt        <- cpc %*% Bs             # Q x n_t; chained to V_struct

      # Bs is symmetric (B and the vchain both are), so the rank-p term's
      # contraction sum_st Bs_st d(v_j g'_j g'_j')_st collapses to
      # 2 * v_j * g'_j' Bs d(g'_j)/dpsi. Precompute the left half once.
      BsdE <- if (is.null(ty)) NULL else sm$dE %*% Bs

      contrib <- function(graw) {
        # graw: RAW derivative of the structural f w.r.t. psi
        dmu      <- as.numeric(crossprod(W, graw))
        term_mu  <- sum(dNLL_dmu_sig * dmu * lnorm_scale)
        term_cov <- 2 * sum(W * rowSums(graw * Bt))
        if (!is.null(ty))
          term_cov <- term_cov +
            2 * sum(ty$var * rowSums(BsdE * crossprod(ty$Dw, graw)))
        term_mu + term_cov
      }
    }

    # V-path of the mean: the residual variance depends on mu, so a parameter that
    # moves mu also moves diag(V) by dv_df = d(var)/d(mu).
    Bvec <- if (!is_var) Bdiag else dNLL_dV_diag  # length n_t

    # For a TBS endpoint the mean scale ms = m'(f) itself depends on f, so moving
    # the structural mean also moves the OFF-diagonal of V_pred (= ms_i ms_j cov_ij).
    # Row k gains 2*m''(f_k)*(A ms)_k with A = dNLL_dV o Cov_eta(f), zero diagonal --
    # the same contraction .admResidMuCoupling() applies for admc. Identically zero
    # unless ms varies with f, so every other error model is untouched.
    ms_off <- numeric(length(mu))
    if (!is_var && !is.null(dres$dms_df) && any(dres$dms_df != 0, na.rm = TRUE))
      ms_off <- 2 * dres$dms_df * .admMsOffDiag(B, cov_f, dres$ms)

    .sigma_V_extra <- function(dmu_raw) sum((Bvec * dres$dv_df + ms_off) * dmu_raw)

    # Struct thetas paired with an eta: reuse the eta's sensitivity column
    # (d(pred)/d(theta) == d(pred)/d(eta) for a mu-referenced theta).
    # struct_eta_idx is eta-indexed (value = struct paired with each eta), so the
    # eta for struct k is which(struct_eta_idx == k).
    for (k in seq_len(n_s)) {
      if (!is.null(pinfo$struct_has_eta) && !pinfo$struct_has_eta[k]) next  # unpaired
      ei <- which(pinfo$struct_eta_idx == k)[1L]  # struct k -> its eta dim
      if (is.na(ei)) next  # nocov -- defensive; ei always found when struct_has_eta[k]
      dmu_raw <- as.numeric(crossprod(W, Jl[[ei]]))  # d(mu_t)/d(psi) before lnorm scaling
      grad[k] <- grad[k] + contrib(Jl[[ei]]) + .sigma_V_extra(dmu_raw)
    }

    # Unpaired struct thetas: their own sensitivity column from the augmented
    # sens model, through the identical formula. Missing -> FD block below.
    if (length(unpaired_k) > 0L) {
      if (is.null(res$dtheta_list)) {
        theta_sens_ok <- FALSE
      } else for (k in unpaired_k) {
        Dt      <- res$dtheta_list[[pinfo$struct_names[k]]]
        dmu_raw <- as.numeric(crossprod(W, Dt))
        g_theta[k] <- g_theta[k] + contrib(Dt) + .sigma_V_extra(dmu_raw)
      }
    }

    # SHIFT path: every structural theta also moves the u nodes, because u's law
    # is sum_j W_j N(Delta_j, omega^2) and Delta depends on the thetas. The
    # chain factor is du/dtheta = E[dDelta/dtheta | u] (see .admShiftDu), and
    # the derivative it multiplies is the affected eta's own sensitivity column,
    # so the contribution is that column scaled ROW-WISE.
    #
    # A theta with no covariate in its Delta -- including every mu-referenced
    # TYPICAL VALUE, which cancels out of a difference of the same expression at
    # two covariate values -- gets a zero column here and is unaffected.
    .sh <- .gS$shift
    if (isTRUE(.sh$cond)) {
      # Conditioned shift: eta is not X L', so every eta column responds to
      # every direction and d(eta) is carried whole (see .admShiftCondBase).
      for (k in seq_len(n_s)) {
        kk <- match(pinfo$struct_names[k], .sh$th_names)
        if (is.na(kk)) next
        base <- .admShiftCondBase(Jl, .sh$dEta_th[[kk]])
        if (is.null(base)) next
        dmu_raw <- as.numeric(crossprod(W, base))
        grad[k] <- grad[k] + contrib(base) + .sigma_V_extra(dmu_raw)
      }
    } else if (isTRUE(.sh$multi)) {
      # Every shifted coordinate moves with every theta, through Delta and
      # through the posterior weights that condition the later levels.
      for (k in seq_len(n_s)) {
        nmk <- pinfo$struct_names[k]
        kk  <- match(nmk, .sh$th_names)
        if (is.na(kk)) next
        base <- .admShiftBase(Jl, .sh$eta_idx, .sh$du[, , kk, drop = FALSE])
        if (is.null(base)) next
        dmu_raw <- as.numeric(crossprod(W, base))
        grad[k] <- grad[k] + contrib(base) + .sigma_V_extra(dmu_raw)
      }
    } else if (isTRUE(.sh$absorb)) {
      # theta moves eta through BOTH mu and Omega + B B'; every eta dimension
      # responds, so the contribution is the full .admAbsorbBase sum.
      for (k in seq_len(n_s)) {
        nmk <- pinfo$struct_names[k]
        if (!nmk %in% colnames(.sh$dmu)) next
        dLt  <- .admCholDiff(.sh$Lt, .sh$dP[[nmk]])
        base <- .admAbsorbBase(Jl, X, dLt, .sh$dmu[, nmk])
        if (is.null(base)) next
        dmu_raw <- as.numeric(crossprod(W, base))
        grad[k] <- grad[k] + contrib(base) + .sigma_V_extra(dmu_raw)
      }
    } else if (!is.null(.sh)) {
      Jsh <- Jl[[.sh$eta_idx]]
      for (k in seq_len(n_s)) {
        nmk <- pinfo$struct_names[k]
        if (!nmk %in% colnames(.sh$du_dtheta)) next
        dk <- .sh$du_dtheta[, nmk]
        if (all(dk == 0)) next
        base    <- Jsh * dk
        dmu_raw <- as.numeric(crossprod(W, base))
        grad[k] <- grad[k] + contrib(base) + .sigma_V_extra(dmu_raw)
      }
    }

    # Omega Cholesky L: d(eta[q,])/d(L_ij) = x[q,j] * e_i (unit vector eta dim i)
    # So d(f[q,])/d(L_ij) = Jl[[i]][q,] * X[q,j]
    # Chain: L_ii stored as log(Omega_ii) -> d(L_ii)/dp = L_ii/2.
    if (n_eta > 0L) for (rr in seq_along(pinfo$omega_par)) {
      i <- pinfo$chol_i[rr]; j <- pinfo$chol_j[rr]
      base <- if (isTRUE(.sh$cond)) {
        # X is zero under conditioning, so Jl[[i]] * X[, j] would contribute
        # nothing at all: the whole omega path is in dEta_om.
        .admShiftCondBase(Jl, .sh$dEta_om[[rr]])
      } else if (isTRUE(.sh$multi)) {
        # A shifted eta's own scale moves EVERY shifted coordinate, so it cannot
        # be folded into X[, j]; an unaffected eta keeps the standard column.
        aa <- match(i, .sh$eta_idx)
        if (i == j && !is.na(aa))
          .admShiftBase(Jl, .sh$eta_idx,
                        .sh$du[, , .sh$n_th + aa, drop = FALSE])
        else Jl[[i]] * X[, j]
      } else if (isTRUE(.sh$absorb)) {
        # Omega enters through chol(Omega + P). dOmega/d(L_ij) = E_ij L' + L
        # E_ij', and P does not depend on the omega parameters at all.
        E <- matrix(0, n_eta, n_eta); E[i, j] <- 1
        .admAbsorbBase(Jl, X, .admCholDiff(.sh$Lt, E %*% t(L) + L %*% t(E)),
                       numeric(n_eta))
      } else Jl[[i]] * X[, j]
      if (is.null(base)) next
      dmu_raw <- as.numeric(crossprod(W, base))
      dL      <- contrib(base) + .sigma_V_extra(dmu_raw)
      pos <- n_s + n_e + rr
      grad[pos] <- grad[pos] + if (pinfo$chol_diag[rr]) dL * L[i, i] / 2 else dL
    }

    # Sigma. Only this output's residual parameters have a nonzero derivative.
    grad[n_s + seq_len(n_e)] <- grad[n_s + seq_len(n_e)] +
      ch$sigma_grad()
  }

  # Unpaired struct thetas: the sens path above already has them exactly.
  if (length(unpaired_k) > 0L && theta_sens_ok) {
    grad[unpaired_k] <- grad[unpaired_k] + g_theta[unpaired_k]
    return(list(grad = grad, nll = if (is.finite(nll_total)) nll_total else Inf))
  }

  # Otherwise CENTRAL FD of .adghNLL.
  #
  # Every perturbed configuration differs from the others only in its structural
  # thetas -- same node grid, same Omega, same sigma -- so they all share one
  # solve per study (.adghMomentsBatch).
  #
  # CENTRAL, not forward, and the STEP is the reason: `grad_h` arrives here as
  # Shi21's measured per-parameter step (.admShi21GradH, probed over exactly this
  # parameter set), and Shi21 minimises the error of a CENTRAL difference,
  # h* = (3 eps_f/|f'''|)^(1/3). The forward optimum is a square root and far
  # coarser, so a forward difference at the central step lands where the eps_f/h
  # noise term dominates. The baseline configuration that used to ride along as
  # configuration 1 is gone with it -- a central difference never evaluates the
  # centre -- so this is 2*n_u configurations against the old n_u + 1, still ONE
  # rxSolve per study, which is the cost that matters (~11 ms per CALL).
  #
  # Joint units keep the per-configuration path (their solve is per output block)
  # and do pay: 2*n_u .adghNLL calls against the old n_u + 1.
  if (length(unpaired_k) > 0L) {
    n_u <- length(unpaired_k)
    hs  <- pmax(abs(p[unpaired_k]), 0.1) * .admGH(grad_h, unpaired_k)
    # 1..n_u are p + h_i; n_u+1..2*n_u are p - h_i, SAME order, so the two halves
    # pair off by index.
    p_pert <- c(
      lapply(seq_len(n_u), function(i) {
        pp <- p; pp[unpaired_k[i]] <- p[unpaired_k[i]] + hs[i]; pp
      }),
      lapply(seq_len(n_u), function(i) {
        pp <- p; pp[unpaired_k[i]] <- p[unpaired_k[i]] - hs[i]; pp
      }))
    n_cfg <- length(p_pert)

    # A SHIFT study joins the joint units on the per-configuration route. Its
    # node grid is a function of the structural thetas (Delta moves the u nodes
    # and their weights), and .adghMomentsBatch assembles every configuration
    # against ONE weight vector, so batching it would score each perturbed
    # configuration on the unperturbed grid -- a finite, plausible, wrong
    # gradient. The shift solve is small enough that 2*n_u of them is cheap.
    if (any(vapply(studies, function(u) isTRUE(u$is_joint) ||
                     identical(u$.adm_cov_path, "shift"), logical(1)))) {
      for (i in seq_len(n_u))
        grad[unpaired_k[i]] <-
          (.adghNLL(p_pert[[i]], pinfo, studies, rxMod, out_var, grid, cores) -
             .adghNLL(p_pert[[n_u + i]], pinfo, studies, rxMod, out_var, grid, cores)) /
          (2 * hs[i])
    } else {
      struct_mat <- do.call(rbind,
        lapply(p_pert, function(pp) .admUnpack(pp, pinfo)$struct))
      colnames(struct_mat) <- names(pars$struct)

      nll_cfg <- numeric(n_cfg)
      for (s in studies) {
        ovs <- s$output %||% out_var
        ms  <- .adghMomentsBatch(struct_mat, pars, pinfo, s, rxMod, ovs, grid, cores)
        for (cfg in seq_len(n_cfg)) {
          if (!is.finite(nll_cfg[cfg])) next
          m     <- ms[[cfg]]
          nll_c <- if (identical(s$method, "var"))
            nll_var_cpp(s$E, s$v_diag, m$E, diag(m$V), s$n)
          else
            nll_cov_cpp(s$E, s$V, m$E, m$V, s$n)
          nll_cfg[cfg] <- if (is.finite(nll_c)) nll_cfg[cfg] + nll_c else Inf
        }
      }
      gk <- (nll_cfg[seq_len(n_u)] - nll_cfg[n_u + seq_len(n_u)]) / (2 * hs)
      # Inf - Inf is NaN; see the same guard in .adfoGrad's Pass 2.
      gk[is.nan(gk)] <- Inf
      grad[unpaired_k] <- gk
    }
  }

  list(grad = grad, nll = if (is.finite(nll_total)) nll_total else Inf)
}

# The gradient alone, with exactly the contract it has always had: a plain named
# numeric and nothing else. Callers that compare or subset it (.adghCalcCov, the
# FD cross-checks in the integration tests) must not have to know that the NLL
# rides along -- an earlier revision returned it as an attribute on this vector
# and it leaked into an expect_equal() of gradient values (#113). The NLL is
# computed either way; forming it from moments that already exist is free.
.adghGrad <- function(p, pinfo, studies, sensModel, rxMod, out_var, grid, cores,
                      grad_h = 1e-4)
  .adghGradNLL(p, pinfo, studies, sensModel, rxMod, out_var, grid, cores,
               grad_h)$grad

# (#5) Pair the objective and the gradient onto ONE solve.
#
# nloptr asks for them as two separate calls, but LBFGS always asks at the same
# p (measured: every gradient call shares its p with an NLL call), and
# .adghGradNLL already forms the moments the NLL needs -- so .adghNLL's solve was
# pure duplicate work. Memoising on p collapses the pair to one solve: ~2x on a
# gradient-mode fit (8.13s -> 3.77s, 58 -> 23 rxSolve on a 3-cmt/5-eta/40-time
# fit; estimates unchanged to 6 dp).
#
# Consequence worth knowing: the objective now comes from the SENSITIVITY solve
# rather than the plain one. Both integrate the same base ODEs, but the
# augmented system makes rxode2's adaptive stepper land ~1e-6 apart (4.6e-11
# relative on the NLL) -- so this is NOT bit-identical to the pre-fusion
# objective, though both sit at the solver's own rtol. It also makes f and
# grad-f self-consistent (one trajectory), where before they came from two.
#
# Only for the analytical-sens path; grad = "fd"/"none" keep the old route.
.adghFusedFns <- function(pinfo, studies, sensModel, rxMod, out_var, grid, cores,
                          grad_h) {
  memo <- new.env(parent = emptyenv())
  memo$key <- NULL; memo$nll <- NULL; memo$grad <- NULL
  ev <- function(p) {
    if (!is.null(memo$key) && identical(memo$key, p)) return(invisible(NULL))
    gn <- .adghGradNLL(p, pinfo, studies, sensModel, rxMod, out_var, grid, cores,
                       grad_h)
    n <- gn$nll
    # nll = NULL => .adghGradNLL degraded to .adghFDGrad, which never formed the
    # moments. Fall back to the ordinary NLL solve rather than invent a value.
    if (is.null(n)) n <- .adghNLL(p, pinfo, studies, rxMod, out_var, grid, cores)
    memo$key <- p; memo$grad <- gn$grad; memo$nll <- n
    invisible(NULL)
  }
  list(nll_fn  = function(p) { ev(p); memo$nll },
       grad_fn = function(p) { ev(p); memo$grad })
}

# -- FD gradient ---------------------------------------------------------------

.adghFDGrad <- function(p, pinfo, studies, rxMod, out_var, grid, cores,
                          grad_h = 1e-4) {
  g <- numeric(length(p)); names(g) <- names(p)
  for (k in seq_along(p)) {
    hk <- pmax(abs(p[k]), 0.1) * .admGH(grad_h, k)
    pp <- p; pp[k] <- p[k] + hk
    pm <- p; pm[k] <- p[k] - hk
    g[k] <- (.adghNLL(pp, pinfo, studies, rxMod, out_var, grid, cores) -
             .adghNLL(pm, pinfo, studies, rxMod, out_var, grid, cores)) / (2 * hk)
  }
  g
}

# -- Covariance ----------------------------------------------------------------

# Post-fit covariance via numerical Hessian over struct + sigma + omega
# (falls back to the struct+sigma sub-block if the full Hessian is not PD).
# Noise-free GH surface -> use tighter eps^(1/4) default step vs admc's eps^(1/5).
# use_grad=TRUE: central FD of gradient (2*np grad evals).
# use_grad=FALSE: full NLL-FD quadratic form (1+2*np+4*n_off NLL evals).
.adghCalcCov <- function(p_hat, pinfo, studies, sensModel, rxMod, out_var,
                           grid, cores,
                           use_grad = TRUE, grad_h = 1e-3,
                           cov_h_outer = .Machine$double.eps^(1/4),
                           sandwich = FALSE
                           ) {
  n_s     <- length(pinfo$struct_names)
  n_e     <- length(pinfo$sigma_names)
  n_o     <- length(pinfo$omega_par)
  # The Hessian now spans struct + sigma + OMEGA. Excluding omega does not just
  # forgo omega's own SEs -- it makes the STRUCTURAL SEs too small, because a theta
  # that carries an eta is correlated with that eta's variance and profiling it out
  # is not the same as fixing it. Measured against the empirical sampling SD over
  # 40 simulated datasets: SE(tcl) rose 8.8% on prop and 8.6% on lnorm when omega
  # was included (it was that much too small), while a purely additive model was
  # unaffected (+0.01%) -- exactly the models where the residual and the IIV
  # compete to explain the same spread. An omega SE from the full Hessian was
  # accurate to about +-20% of the empirical SD.
  #
  # The omega Cholesky is more weakly identified than struct/sigma, so the full
  # Hessian can be non-PD where the struct+sigma block is fine. That is handled by
  # falling back to the struct+sigma sub-block rather than returning nothing.
  n_sub   <- n_s + n_e
  cov_idx <- seq_len(n_sub + n_o)
  np_cov  <- length(cov_idx)
  nms_cov <- names(p_hat)[cov_idx]

  nll_fn  <- function(p)
    suppressMessages(.adghNLL(p, pinfo, studies, rxMod, out_var, grid, cores))
  grad_fn <- function(p)
    suppressMessages(.adghGrad(p, pinfo, studies, sensModel, rxMod, out_var,
                                grid, cores, grad_h = grad_h))

  nll0 <- nll_fn(p_hat)
  if (!is.finite(nll0)) {
    warning("adghCalcCov: NLL not finite at p_hat -- covariance not computed",
            call. = FALSE)
    return(NULL)
  }

  H <- matrix(0, np_cov, np_cov, dimnames = list(nms_cov, nms_cov))

  # Covariate studies: differentiate the NLL, not the gradient.
  #
  # The original reason given here was that .adghGradNLL builds its quadrature
  # from pars$L rather than through .adghGrid(), so it could not carry the
  # covariate grid. That is NOT true and has not been for some time -- it calls
  # .adghGrid() per study and then .adghStudyCov(), so it is covariate-aware.
  # Measured on a 1-cmt lognormal-covariate model at its optimum, the two forms
  # agree to 2.7e-05 on the covariance with identical standard errors, and the
  # gradient form is ~1.5x faster.
  #
  # The guard is kept anyway, deliberately. It buys one Hessian per fit, that
  # measurement covers a single model, and the failure it would expose --
  # standard errors computed from a different objective than the estimates --
  # is silent and severe. Removing it wants a broader comparison than one
  # model, not a rewritten comment.
  if (isTRUE(use_grad) &&
      any(vapply(studies, function(s) !is.null(s[["cov_dist"]]), logical(1))))
    use_grad <- FALSE

  if (use_grad) {
    # CENTRAL difference of the gradient -- see .adfoCalcCov() for the reasoning
    # and the cost (2*np_cov gradient evaluations against the old np_cov+1).
    h_c <- pmax(abs(p_hat[cov_idx]), 0.1) * cov_h_outer
    for (jj in seq_len(np_cov)) {
      ph      <- p_hat; ph[cov_idx[jj]] <- ph[cov_idx[jj]] + h_c[jj]
      pm      <- p_hat; pm[cov_idx[jj]] <- pm[cov_idx[jj]] - h_c[jj]
      gp      <- grad_fn(ph)[cov_idx]
      gm      <- grad_fn(pm)[cov_idx]
      H[, jj] <- if (anyNA(gp) || anyNA(gm)) 0 else (gp - gm) / (2 * h_c[jj])
    }
    H <- (H + t(H)) / 2
  } else {
    # Step selection. `pmax(abs(p), 0.1) * cov_h_outer` is a guess about how much
    # noise the objective carries, applied identically to every parameter -- and
    # it is the guess behind the "Hessian not positive definite ... try
    # increasing cov_h_outer" warning below. Gill83 measures instead: it probes
    # THIS objective and returns the step where condition error and truncation
    # error balance, per parameter. Exact fit here, since the function it probes
    # is the one being differenced.
    h_fd <- .admHessSteps(nll_fn, p_hat, cov_idx, cov_h_outer,
                            .var.name = "adghCalcCov")
    for (k in seq_len(np_cov)) {
      ki <- cov_idx[k]; hk <- h_fd[k]
      p_p <- p_hat; p_p[ki] <- p_p[ki] + hk
      p_m <- p_hat; p_m[ki] <- p_m[ki] - hk
      H[k, k] <- (nll_fn(p_p) - 2 * nll0 + nll_fn(p_m)) / hk^2
    }
    for (i in seq_len(np_cov - 1L)) {
      for (j in seq(i + 1L, np_cov)) {
        ii <- cov_idx[i]; ji <- cov_idx[j]
        hi <- h_fd[i];  hj <- h_fd[j]
        p_pp <- p_hat; p_pp[ii] <- p_pp[ii] + hi; p_pp[ji] <- p_pp[ji] + hj
        p_pm <- p_hat; p_pm[ii] <- p_pm[ii] + hi; p_pm[ji] <- p_pm[ji] - hj
        p_mp <- p_hat; p_mp[ii] <- p_mp[ii] - hi; p_mp[ji] <- p_mp[ji] + hj
        p_mm <- p_hat; p_mm[ii] <- p_mm[ii] - hi; p_mm[ji] <- p_mm[ji] - hj
        H[i, j] <- H[j, i] <-
          (nll_fn(p_pp) - nll_fn(p_pm) - nll_fn(p_mp) + nll_fn(p_mm)) / (4 * hi * hj)
      }
    }
  }

  if (!all(is.finite(H))) {
    warning("adghCalcCov: Hessian has non-finite entries -- covariance not computed",
            call. = FALSE)
    return(NULL)
  }

  eig_dec <- tryCatch(eigen(H, symmetric = TRUE), error = function(e) NULL)
  H_eigs  <- if (!is.null(eig_dec)) eig_dec$values else rep(NA_real_, np_cov)

  # If the weakly-identified omega Cholesky makes the full Hessian indefinite, drop
  # back to the struct+sigma sub-block (which is what this returned before omega was
  # included) rather than reporting nothing -- same threshold and retained rows as
  # adfo/admc, via the shared helper. .invert() then does adgh's own inversion.
  .red <- .admReduceNpdOmega(H, H_eigs, eig_dec, nms_cov, n_o, n_sub)
  if (.red$reduced)
    warning("adghCalcCov: the full Hessian including omega was not positive ",
            "definite or was numerically singular; reporting structural and ",
            "sigma standard errors only.", call. = FALSE)
  H <- .red$H; nms_cov <- .red$nms_cov; eig_dec <- .red$eig_dec; H_eigs <- .red$H_eigs

  .invert <- function(M) {
    e <- tryCatch(eigen(M, symmetric = TRUE, only.values = TRUE),
                  error = function(e) NULL)
    if (is.null(e) || min(e$values) < 0) return(NULL)
    tryCatch(chol2inv(chol(M)),
             error = function(e) tryCatch(solve(M), error = function(e2) NULL))
  }
  Hinv <- .invert(H)
  if (is.null(Hinv)) {
    warning(sprintf(
      "adghCalcCov: Hessian not positive definite or not invertible (min eigenvalue %.3e). Covariance not computed. Try increasing cov_h_outer (currently %.3e), e.g. cov_h_outer = %.3e.",
      if (length(H_eigs)) min(H_eigs) else NA_real_, cov_h_outer, cov_h_outer * 4),
      call. = FALSE)
    return(NULL)
  }

  cov_full <- (2 * Hinv + t(2 * Hinv)) / 2
  # covMethod = "r,s": replace the 2H^-1 filling with H^-1 J H^-1, built on the
  # SAME H so the two cannot disagree about the half they share. Under correct
  # specification J = 2H and this returns what "r" would have. Anything the
  # sandwich cannot supply -- a residual outside the conditionally-normal family,
  # a joint unit, a singular ingredient -- degrades to "r" and says so, rather
  # than reporting a number of unknown provenance.
  sw_used <- FALSE; sw_HJ <- NULL
  if (isTRUE(sandwich)) {
    sw <- tryCatch(
      .admSandwichCov(p_hat, pinfo, studies, rxMod, out_var, grid, cores,
                      H = H, keep = match(nms_cov, names(p_hat)),
                      sensModel = sensModel),
      error = function(e) NULL)
    ok <- !is.null(sw) && all(is.finite(sw$cov)) && all(diag(sw$cov) > 0)
    if (ok) {
      cov_full <- (sw$cov + t(sw$cov)) / 2
      sw_used  <- TRUE
      # H and J travel with the covariance because two more things are built
      # from exactly this pair: the TIC penalty tr(H^-1 J), and the eigenvalue
      # weights admCompare() rescales dOFV by. Recomputing them later would mean
      # re-solving, and would let them drift from the SE actually reported.
      sw_HJ    <- list(H = sw$H, J = sw$J, par_names = nms_cov,
                       Om = sw$Om, study_names = names(studies))
    } else {
      warning("adghCalcCov: the sandwich correction could not be computed; ",
              "reporting the covMethod = \"r\" covariance instead.", call. = FALSE)
    }
  }
  dimnames(cov_full) <- list(nms_cov, nms_cov)
  # Rotate onto the reported scale (residual delta factors + omega Jacobian). One
  # shared implementation for all three estimators -- see .admScaleReportedCov().
  out <- .admScaleReportedCov(cov_full, p_hat, pinfo, n_s, n_e, n_o, n_sub)
  attr(out, "sandwich") <- sw_used
  attr(out, "sandwich_HJ") <- sw_HJ
  out
}

# -- Restart worker ------------------------------------------------------------

# Self-contained GH optimization run (one restart); serializable to a worker.
# Signature mirrors .adfoRestartWorker: same base_args from .admRunRestarts().
# n_sim, sampling accepted for interface compatibility but not used.
.adghRestartWorker <- function(restart_id, p_init, ui_lstExpr, pinfo,
                                ov_lower, ov_upper, scale_c = NULL, studies, n_sim,
                                seed, n_nodes, algorithm, ftol_rel, maxeval,
                                use_grad, grad_h, grad_bounds,
                                output_var = "cp",
                                sampling = "sobol",
                                use_pure_fd = FALSE,
                                print_progress = TRUE, print = 10L,
                                cores = NULL, no_lock = FALSE,
                                sens_cache_file = NULL, sens_cols = NULL,
                                sens_rename = NULL,
                                rxMod_direct = NULL, sensModel_direct = NULL) {
  library(admixr2)
  tryCatch(.admPatchDevNamespace(), error = function(e) NULL)

  m <- .admWorkerLoadModels(ui_lstExpr, rxMod_direct, cores,
                            sens_cache_file, sens_cols, sens_rename, sensModel_direct,
                            pinfo)

  grid <- .adghNodeGrid(n_nodes, pinfo$n_eta)
  set.seed(seed + restart_id)

  # (#5) Same objective/gradient fusion the single-fit driver uses, so a restart
  # is not twice the solves of the fit it restarts.
  .fz <- if (!use_pure_fd && !is.null(m$sensModel))
    .adghFusedFns(pinfo, studies, m$sensModel, m$rxMod, output_var, grid,
                  m$cores_w, grad_h) else NULL

  nll_fn <- if (!is.null(.fz)) .fz$nll_fn else
    function(p) .adghNLL(p, pinfo, studies, m$rxMod, output_var, grid, m$cores_w)

  grad_fn <- if (use_pure_fd) {
    function(p) .adghFDGrad(p, pinfo, studies, m$rxMod, output_var, grid, m$cores_w,
                            grad_h)
  } else if (!is.null(.fz)) {
    .fz$grad_fn
  } else {
    function(p) .adghGrad(p, pinfo, studies, m$sensModel, m$rxMod, output_var,
                          grid, m$cores_w, grad_h)
  }

  # adgh loads its own model in-process and does not lock (single-nloptr path).
  .res <- .admScaledOptimize(restart_id, p_init, ov_lower, ov_upper, scale_c,
                     use_grad, grad_bounds, algorithm, ftol_rel, maxeval,
                     nll_fn, grad_fn, pinfo, print_progress, print,
                     lock_rxMod = NULL)
  # Carried back so .admRunRestarts() can report a worker that silently
  # dropped to a finite-difference gradient -- a daemon's own warning is
  # swallowed by mirai. NOT a new worker ARGUMENT: the signatures must stay
  # stable for a daemon resolving them from the installed namespace.
  .res$sens_fallback <- m$sens_fallback
  .res
}

# -- Control object ------------------------------------------------------------

#' Control settings for the Gauss-Hermite (GH) quadrature estimator
#'
#' Creates a control object for `nlmixr2(est = "adgh")`. The GH estimator
#' integrates model predictions against the random-effects prior
#' \eqn{\eta \sim N(0, \Omega)} using a deterministic tensor-product
#' Gauss-Hermite quadrature grid. It is unbiased at any IIV magnitude (unlike
#' FO), noise-free (unlike MC), and much faster than MC for models with up to
#' ~4 etas.
#'
#' @param studies Named list of study specifications (same format as
#'   [admControl()]: `E`, `V`, `n`, `times`, `ev`, optional `method`; or an
#'   `observations` list for multi-compartment fits -- see [admControl()]).
#' @param cov_nodes Gauss-Hermite nodes per covariate used to integrate the
#'   COVARIATE distribution when a study declares `cov_dist` (default 7). This is
#'   a separate dial from `n_nodes`, which refines the random-effect dimensions
#'   only: raising `n_nodes` alone leaves the covariate integration exactly where
#'   it was. Measured on a two-compartment model with an allometric weight effect
#'   and a lognormal weight distribution, 7 nodes place the marginal moments
#'   within 2e-06 (mean) and 2e-05 (covariance) of an exact reference, and the
#'   remaining error is the ODE solver's rather than the quadrature's. A wider or
#'   more skewed covariate distribution, or a more strongly non-linear covariate
#'   effect, warrants more. Measured against an exact reference on a
#'   two-compartment model with an allometric weight effect and a lognormal
#'   weight distribution: 3 nodes give 7.3e-04 / 8.2e-03 (mean / covariance),
#'   5 give 2.8e-05 / 3.7e-04, 7 give 2.2e-06 / 2.4e-05, and 9 onwards sit at
#'   ~1.2e-06 / ~1.0e-06, which is the ODE solver's accuracy rather than the
#'   quadrature's. The default is set past that knee, and raising it further
#'   buys nothing: against a per-subject reference the accuracy is identical at
#'   5, 9 and 15 nodes. Ignored when `cov_integration = "taylor"`.
#' @param cov_integration How a study's covariate distribution is integrated.
#'   `"quadrature"` (default) evaluates the model on a product Gauss-Hermite
#'   grid of `cov_nodes` points per covariate and forms the marginal moments
#'   from the whole grid; it is the accurate route and the one every existing
#'   fit uses. `"taylor"` instead expands the marginal MOMENTS to second order
#'   about the covariate mean,
#'   \eqn{E \approx g(m) + \frac{1}{2}\sum_j v_j g''_j(m)} and
#'   \eqn{V \approx V_c(m) + \frac{1}{2}\sum_j v_j V_{c,j}''(m) + \sum_j v_j
#'   g'_j(m) g'_j(m)^T} for \eqn{g(a) = E_\eta[f(a,\eta)]}, and evaluates the
#'   likelihood once at those approximate moments. That costs `1 + 2p` covariate
#'   points for `p` covariates instead of `cov_nodes^p`, so it is a speed lever
#'   for models with several covariates or a wide grid.
#'
#'   It is an APPROXIMATION, and how good depends on how far the covariate
#'   pushes the model relative to the random effects. With
#'   `ratio = theta_cov * sd_a / omega`, the measured relative error of the
#'   moments against the exact marginal is 1e-07 / 1e-05 (mean / covariance) at
#'   `ratio = 0.1`, 5e-05 / 5e-03 at 0.5, 7e-04 / 5e-02 at 1, and 3e-02 / 4e-01
#'   at 3. Check an important fit against `"quadrature"`.
#'
#'   DEPENDENT covariates (`cor`, `rho`, `Sigma`, or a `joint` sampler) are
#'   supported by both routes. The quadrature grid is a product rule over the
#'   sampler's UNIFORMS rather than over the covariate margins, and a copula
#'   maps independent uniforms to dependent values, so the product rule stays
#'   exact whatever the dependence. The expansion differences along the
#'   EIGENVECTORS of the covariate covariance, which keeps it at `1 + 2p` points
#'   at any correlation --- a coordinate-basis version would need the mixed
#'   partials explicitly, at `1 + 2p + p(p-1)` points for the same accuracy.
#'   (This is the third-degree spherical-radial cubature rule, i.e. the
#'   unscented transform's sigma points.)
#'
#'   Refused rather than approximated: a discrete `values` covariate (the
#'   differencing step would land between the levels), and (near-)perfectly
#'   collinear covariates (the step along the null direction collapses, so the
#'   second difference is cancellation rather than a derivative).
#'
#'   `"shift"` removes the covariate from the solve entirely. When the
#'   covariates act on the model only through a mu-referenced argument they are
#'   a pure shift of that argument's random effect,
#'   \eqn{f(a,\eta) = f(a_{ref}, \eta + \Delta(a))}, so the two-dimensional
#'   \eqn{(a, \eta)} integral collapses onto one over \eqn{u = \Delta(a) +
#'   \eta}. The solve cost is then CONSTANT in the number of covariates, against
#'   `cov_nodes^p` for the grid, and it holds for any covariate distribution ---
#'   discrete, skewed, dependent alike. The precondition is checked NUMERICALLY
#'   against the compiled model, never read off the model text (measured
#'   separation on the probe set: valid forms 1e-12--1e-14, invalid forms
#'   3e-02--7e-01), and a model that fails it is an ERROR naming the reason.
#'
#'   `"auto"` tries `"shift"` and falls back to `"quadrature"` with a message
#'   when the identity does not hold, so it is the fast route wherever the fast
#'   route is valid and the accurate route everywhere else. It is not the
#'   default only because switching it on moves an existing covariate fit's
#'   numbers at the ~1e-5 level (the tolerance at which shift and grid agree);
#'   for a new fit it is the recommended setting.
#' @param cov_taylor_h Radius of the design points for
#'   `cov_integration = "taylor"`, as a multiple of the moment-matched radius
#'   \eqn{\sqrt{3\lambda_k}} along each expansion direction (default 1). The
#'   \eqn{1+2p} design is a cubature rule rather than a finite-difference
#'   stencil, so the radius is chosen to integrate moments exactly, not to be
#'   small: at the default the design coincides with 3-point Gauss--Hermite in
#'   each direction and is exact through degree 5, where a smaller radius
#'   matches only through degree 3. For independent covariates the directions
#'   are the covariates themselves; when they are dependent they are the
#'   eigenvectors of the covariate covariance. Values below 1 pull the points
#'   back toward the covariate mean -- the scaled unscented transform -- which
#'   costs the fourth-moment match but keeps the design inside the range a
#'   strongly skewed margin actually spans. Ignored when
#'   `cov_integration = "quadrature"`.
#' @param resid_nodes Gauss-Hermite nodes used to integrate the RESIDUAL for a
#'   transform-both-sides endpoint (`boxCox`, `yeoJohnson`, `logitNorm`,
#'   `probitNorm`), where `y = g(h(f) + sigma*eps)` has no closed-form mean and
#'   variance. Ignored by every other error model, which has closed forms. Default
#'   81. Measured worst-case relative error against an independent quadrature, over
#'   all four transforms and residual SD of 0.5, 1, 2 and 3: n = 15 gives 5.7e-2,
#'   31 gives 4.5e-3, 81 gives 5.0e-5. The error is dominated by large residual SD;
#'   at SD <= 1, n = 31 already gives 1e-7 or better.
#'
#'   This is an ACCURACY dial, not a speed one. The quadrature is linear in
#'   `resid_nodes` in isolation (~50 us at 15, 300 us at 81 for an 8-row study) but
#'   negligible beside the ODE solve: a full NLL evaluation measured 0.750 s per 60
#'   evaluations at BOTH 31 and 81 nodes. Raise it if you have a saturating endpoint
#'   with a large residual SD; there is little to gain by lowering it.
#' @param n_nodes Number of quadrature nodes per eta dimension (default 5).
#'   Total nodes = `n_nodes^n_eta`. `n_nodes = 5` achieves near-exact covariance
#'   moments for IIV SD up to ~0.5; `n_nodes = 7` extends coverage to SD ~0.7.
#'   For models with >= 5 etas the node count grows steeply; consider reducing
#'   `n_nodes` or using a different estimator.
#' @param grad Gradient mode. `"analytical"` (default) uses closed-form
#'   contractions through the sensitivity equations -- cheapest and exact.
#'   `"fd"` uses central finite differences (forward differencing was removed
#'   in 0.4.1; see `adfoControl()`).
#'   `"none"` uses derivative-free BOBYQA.
#' @param algorithm nloptr algorithm, or `NULL` (default) to pick the default
#'   that matches `grad`: `"NLOPT_LD_LBFGS"` with a gradient, `"NLOPT_LN_BOBYQA"`
#'   when `grad = "none"`. Any algorithm reported by
#'   [nloptr::nloptr.print.options()] is accepted. An explicit algorithm is
#'   reconciled with `grad`: when `grad = "none"` a gradient-based algorithm
#'   (`NLOPT_LD_*` / `NLOPT_GD_*`) falls back to `"NLOPT_LN_BOBYQA"`; when a
#'   gradient is requested a derivative-free algorithm (`NLOPT_LN_*` /
#'   `NLOPT_GN_*`) turns the gradient off. Both emit a message.
#' @param maxeval Maximum function evaluations (default 500).
#' @param ftol_rel Relative tolerance (default `sqrt(.Machine$double.eps)`).
#' @param print Print-frequency for live progress (0 = silent).
#' @param seed Random seed (used for restarts).
#' @param cores OpenMP threads for `rxSolve()`. Defaults to
#'   `rxode2::rxCores()`. When `workers > 1` it is a *total* budget, split
#'   across the workers.
#' @param nDisplayProgress Passed to `rxSolve()`: show the solver's text
#'   progress bar only once a single solve exceeds this many subjects. The
#'   default (`.Machine$integer.max`) keeps it off for clean script/vignette
#'   output; lower it (e.g. `1000L`) to see progress during long fits.
#' @param grad_h Finite-difference step for unpaired struct theta gradient and
#'   FD Jacobian fallback.
#' @param grad_bounds Box-constraint half-width when using gradients: the fit is
#'   confined to `p0 +/- grad_bounds` on the optimizer scale, which for a
#'   log-scale parameter is a factor of `exp(grad_bounds)` (~148 at the default
#'   5). This bound is admixr2's, not the model's -- an unbounded parameter has
#'   no other -- and nloptr reports normal convergence at a box corner, so a
#'   warning is emitted if an estimate finishes on it.
#' @param cov_h Inner FD step for the gradient-based Hessian (only used when
#'   `covMethod = "r"` and `grad != "none"`).
#' @param cov_h_outer Outer step scale for numerical Hessian. Default
#'   `eps^(1/4)` (tighter than admc's `eps^(1/5)` because the GH surface is
#'   noise-free).
#' @param covMethod `"r"` computes covariance via a numerical Hessian over the
#'   structural, residual-error and omega parameters; `"none"` skips it. Omega is
#'   included because excluding it also biases the STRUCTURAL standard errors
#'   downward -- a theta carrying an eta is correlated with that eta's variance.
#'   If the weakly-identified omega Cholesky makes the Hessian non-positive
#'   definite, the structural + residual sub-block is reported with a warning.
#'
#'   `"r,s"` adds a sandwich correction, `H^-1 J H^-1`, on the same Hessian.
#'   The aggregate objective scores the reported mean and covariance as though
#'   the subjects behind them were multivariate normal; they are not, because the
#'   model is nonlinear in the random effects, so the sampling law of `(E, V)` is
#'   not the one the objective assumes. `"r,s"` scores that law from the model
#'   instead. Point estimates are untouched -- only the reported uncertainty
#'   changes -- and under correct specification it reduces to `"r"` exactly.
#'   Available for the conditionally-normal residual family (`add`, `prop`,
#'   `pow`, `combined1`, `combined2`, `lnorm`) and for non-joint studies;
#'   anything else falls back to `"r"` with a warning.
#'
#'   All three blocks are reported on the scale the ESTIMATES are printed on, as
#'   `nlmixr2est` does: structural thetas on the log/optimizer scale, residual
#'   error as an SD, and omega as the variance/covariance entries (named
#'   `om.<eta>` and `cov.<eta_i>.<eta_j>`). The omega block is rotated by the
#'   full Jacobian of Omega with respect to the log-Cholesky, which is not
#'   diagonal once omega is correlated.
#' @param n_restarts Number of optimizer restarts (1 = no multi-start).
#' @param restart_sd SD of random perturbations of initial struct thetas at
#'   each restart.
#' @param workers Number of parallel workers (mirai daemons) for multi-restart
#'   (default 1 = sequential). Requires the `mirai` package.
#' @param rxControl `rxode2::rxControl()` object. Created automatically when `NULL`.
#' @param sigdig Significant digits asked of the ODE solver, or `NULL` (the
#'   default) to leave rxode2's own solver tolerances alone. When set, it is
#'   passed to `rxode2::rxSolve()`'s own `sigdig` argument for every solve the
#'   estimator issues -- rxode2 owns the mapping to `atol`/`rtol` and has changed
#'   it between releases, which is why the digits, not the tolerances, are what
#'   travels -- and to `nlmixr2est::foceiControl()` for the post-fit tables.
#'
#'   It is a speed lever, and an opt-in one because it is not free. The
#'   estimators finite-difference the solve with steps of the same order:
#'   `grad_h` (1e-4), `cov_h` (1e-3) and `cov_h_outer` (~2.5e-3), while
#'   `sigdig = 4` maps to a relative tolerance of ~1e-4 on current rxode2.
#'   Differencing a solution whose own noise is 1e-4 with a 1e-4 step returns
#'   noise, and it surfaces as a moved objective and an indefinite covariance
#'   Hessian (every `SE` reported `NA`) rather than as an error. Most worthwhile
#'   where the gradient is fully analytic and nothing differences the solve --
#'   `adfoControl(grad = "analytical")` measured ~4.8x faster at `sigdig = 4`
#'   with standard errors unchanged to 4 significant figures. Elsewhere, compare
#'   the objective and the standard errors against `NULL` before relying on it.
#'   Table formatting is unaffected either way: `sigdigTable` defaults to 4
#'   regardless.
#' @param calcTables,compress,ci,sigdigTable,optExpression,sumProd,literalFix
#'   Passed to `nlmixr2est::foceiControl()` for the table/output machinery.
#' @param addProp How combined additive+proportional error is parameterised in
#'   the nlmixr2 output tables: `"combined2"` (default) or `"combined1"`.
#' @param returnAdmr If `TRUE`, return a plain list instead of the full
#'   nlmixr2 fit object.
#' @param ... Unused arguments (trigger an error).
#'
#' @return An `adghControl` object (a named list).
#'
#' @seealso [admControl()], [adfoControl()], [adirmcControl()]
#'
#' @examples
#' ctl <- adghControl()
#' ctl$n_nodes
#' ctl$grad
#'
#' # More nodes for large IIV, analytical gradient
#' ctl2 <- adghControl(n_nodes = 7L, grad = "analytical", maxeval = 300L)
#'
#' \donttest{
#' library(rxode2)
#' library(nlmixr2)
#'
#' data("examplomycin")
#' obs    <- examplomycin[examplomycin$EVID == 0, ]
#' obs    <- obs[order(obs$ID, obs$TIME), ]
#' times  <- sort(unique(obs$TIME))
#' ids    <- unique(obs$ID)
#' dv_mat <- do.call(rbind, lapply(ids, function(i) {
#'   sub <- obs[obs$ID == i, ]; sub$DV[order(sub$TIME)]
#' }))
#' E <- colMeans(dv_mat)
#' V <- cov.wt(dv_mat, method = "ML")$cov
#'
#' pk_model <- function() {
#'   ini({
#'     tcl <- log(5); tv <- log(30)
#'     prop.sd <- c(0, 0.2)
#'     eta.cl ~ 0.09; eta.v ~ 0.04
#'   })
#'   model({
#'     cl <- exp(tcl + eta.cl)
#'     v  <- exp(tv  + eta.v)
#'     d/dt(central) <- -(cl/v) * central
#'     cp <- central / v
#'     cp ~ prop(prop.sd)
#'   })
#' }
#'
#' fit <- nlmixr2(
#'   pk_model, admData(), est = "adgh",
#'   control = adghControl(
#'     studies = list(study1 = list(E = E, V = V, n = length(ids),
#'                                  times = times, ev = et(amt = 100)))
#'   )
#' )
#' }
#'
#' @export
adghControl <- function(
    studies     = list(),
    n_nodes     = 5L,
    grad        = c("analytical", "fd", "none"),
    algorithm   = NULL,
    maxeval     = 500L,
    ftol_rel    = .Machine$double.eps^(1/2),
    print       = 10L,
    seed        = 12345L,
    cores       = rxode2::rxCores(),
    nDisplayProgress = .Machine$integer.max,
    grad_h      = 1e-4,
    grad_bounds = 5,
    cov_h       = 1e-3,
    cov_h_outer = .Machine$double.eps^(1/4),
    covMethod   = c("r", "r,s", "none"),
    n_restarts  = 1L,
    restart_sd  = 0.5,
    workers     = 1L,
    rxControl     = NULL,
    calcTables    = FALSE,
    compress      = TRUE,
    ci            = 0.95,
    sigdig        = NULL,
    sigdigTable   = NULL,
    addProp       = c("combined2", "combined1"),
    optExpression = TRUE,
    sumProd       = FALSE,
    literalFix    = TRUE,
    returnAdmr    = FALSE,
    # LAST on purpose: inserting an argument mid-signature silently rebinds every
    # positional call -- adghControl(studies, 7L) used to set n_nodes = 7.
    resid_nodes   = 81L,
    # LAST on purpose: a new argument inserted mid-signature silently rebinds
    # every positional call. See the resid_nodes note in CLAUDE.md.
    cov_nodes     = 7L,
    # LAST on purpose, as above. These two are the covariate-integration pair:
    # cov_integration selects the method, cov_taylor_h the differencing step it
    # uses. Appended together so the tail of this signature reads as one addition.
    cov_integration = c("quadrature", "auto", "taylor", "shift"),
    cov_taylor_h    = 1,
    ...) {

  .xtra <- list(...)
  if (length(.xtra) > 0L)
    stop("adghControl: unused argument(s): ",
         paste(paste0("'", names(.xtra), "'"), collapse = ", "), call. = FALSE)

  addProp   <- match.arg(addProp)
  grad      <- match.arg(grad)
  covMethod <- match.arg(covMethod)
  cov_integration <- match.arg(cov_integration)

  checkmate::assertList(studies)
  checkmate::assertIntegerish(n_nodes,     lower = 1L, len = 1)
  # A residual quadrature needs a real grid. .adghNodes1() refuses m < 1, but it
  # accepts 1..4 happily and returns a rule that integrates nothing usefully --
  # the measured error at 5 nodes is already 3.3e-1. Refuse here, where the
  # message can name the argument, rather than silently scoring a wrong NLL.
  checkmate::assertIntegerish(resid_nodes, lower = 5L, len = 1)
  checkmate::assertIntegerish(cov_nodes, lower = 1L, len = 1)
  # cov_taylor_h is a MULTIPLIER on the moment-matched radius sqrt(3*lambda),
  # not a raw step: 1 puts the design points exactly where 3-point
  # Gauss-Hermite does, which is what makes the rule exact through degree 5.
  # Below 1 pulls them back toward the covariate mean (the scaled unscented
  # transform), trading the fourth-moment match for design points that stay
  # nearer the covariate range a skewed margin actually spans. Only positivity
  # is enforced -- h = 0 is a division by zero in the second difference.
  checkmate::assertNumeric(cov_taylor_h, lower = .Machine$double.eps, len = 1,
                           finite = TRUE)
  # NOT assertString(algorithm) here: NULL is now the default and means "pick the
  # one that matches grad". .admResolveAlgorithm() asserts the string and checks
  # it against the installed nloptr, which is more than this line ever did.
  checkmate::assertIntegerish(maxeval,     lower = 1L, len = 1)
  checkmate::assertNumeric(ftol_rel,       lower = 0,  len = 1)
  checkmate::assertIntegerish(print,       lower = 0L, len = 1)
  checkmate::assertIntegerish(seed,                    len = 1)
  checkmate::assertIntegerish(cores,       lower = 1L, len = 1)
  checkmate::assertIntegerish(nDisplayProgress, lower = 1L, len = 1,
                              .var.name = "nDisplayProgress")
  checkmate::assertNumeric(grad_h,         lower = 0,  len = 1)
  checkmate::assertNumeric(grad_bounds,    lower = 0,  len = 1)
  checkmate::assertNumeric(cov_h,          lower = 0,  len = 1)
  checkmate::assertNumeric(cov_h_outer,    lower = 0,  len = 1)
  checkmate::assertIntegerish(n_restarts,  lower = 1L, len = 1)
  checkmate::assertNumeric(restart_sd,     lower = 0,  len = 1)
  checkmate::assertIntegerish(workers,     lower = 1L, len = 1)
  checkmate::assertNumeric(ci, lower = 0, upper = 1,   len = 1)
  if (!is.null(sigdig))
  checkmate::assertIntegerish(sigdig,      lower = 1L, len = 1)
  checkmate::assertLogical(returnAdmr,                 len = 1)

  # .admResolveAlgorithm(), like the other three controls -- adgh had its own
  # two-line rule instead, and it was one-directional and unvalidated:
  #   * adghControl(algorithm = "NOT_AN_ALGO") was ACCEPTED and carried all the
  #     way to nloptr, where it surfaces as a cryptic error mid-fit;
  #   * adghControl(grad = "analytical", algorithm = "NLOPT_LN_NELDERMEAD") kept
  #     BOTH -- a derivative-free algorithm with the gradient still switched on,
  #     so every iteration paid for a gradient nloptr discards. The other three
  #     turn the gradient off here and say so.
  # The four behaviours test-adgh-nodes.R pins are unchanged: NULL + grad "fd" /
  # "fd" / "analytical" -> LBFGS, grad "none" -> BOBYQA, and an explicit
  # gradient-based algorithm is kept as given.
  .alg <- .admResolveAlgorithm(algorithm, grad,
                               .var.name = "adghControl: algorithm")
  algorithm <- .alg$algorithm
  grad      <- .alg$grad

  # sigdig = NULL (the DEFAULT) means "leave rxode2's own solver defaults alone".
  # It is the one setting whose meaning does not move under an rxode2 upgrade,
  # and it is the default because a looser solve is not free: this release is
  # what first routed sigdig into the estimators' own rxSolve calls, and every
  # finite-difference step that consumes those solves (grad_h 1e-4, cov_h 1e-3,
  # cov_h_outer ~2.5e-3) is the same order as the tolerance sigdig = 4 asks for.
  # rxode2 5.1.5 maps sigdig = 4 to rtol = 1e-4 (5.1.4 mapped it to 5e-7 -- 200x
  # tighter for the same request), so differencing with a 1e-4 step differences
  # noise: a moved objective and an indefinite Hessian, not an error. Shipping it
  # on by default would have changed the numerics of every existing script
  # silently, for a knob that looked like table formatting before this release.
  #
  # NULL is also the only way back: the sigdig -> tolerance map is
  # one-dimensional while rxode2's defaults are not (atol 1e-8 vs rtol 1e-6), so
  # no sigdig value reproduces them. The tables still need a number, so they fall
  # back to 4 -- i.e. sigdigTable is unchanged whichever way sigdig is set.
  if (is.null(rxControl))   rxControl   <- if (is.null(sigdig))
    rxode2::rxControl() else rxode2::rxControl(sigdig = sigdig)
  if (is.null(sigdigTable)) sigdigTable <- if (is.null(sigdig)) 4L else
    max(round(sigdig), 3L)

  .ret <- list(
    studies       = studies,
    resid_nodes   = as.integer(resid_nodes),
    cov_nodes     = as.integer(cov_nodes),
    cov_integration = cov_integration,
    cov_taylor_h    = cov_taylor_h,
    n_nodes       = as.integer(n_nodes),
    n_sim         = 1L,       # interface compat with .admRunRestarts()
    sampling      = "sobol",  # idem
    grad          = grad,
    algorithm     = algorithm,
    maxeval       = as.integer(maxeval),
    ftol_rel      = ftol_rel,
    print         = as.integer(print),
    seed          = as.integer(seed),
    cores         = as.integer(cores),
    nDisplayProgress = as.integer(nDisplayProgress),
    grad_h        = grad_h,
    grad_bounds   = grad_bounds,
    cov_h         = cov_h,
    cov_h_outer   = cov_h_outer,
    covMethod     = covMethod,
    n_restarts    = as.integer(n_restarts),
    restart_sd    = restart_sd,
    workers       = as.integer(workers),
    rxControl     = rxControl,
    calcTables    = calcTables,
    compress      = compress,
    ci            = ci,
    sigdig        = sigdig,
    sigdigTable   = as.integer(sigdigTable),
    addProp       = addProp,
    optExpression = optExpression,
    sumProd       = sumProd,
    literalFix    = literalFix,
    returnAdmr    = returnAdmr
  )
  class(.ret) <- "adghControl"
  .ret
}

# -- nlmixr2 S3 hooks ----------------------------------------------------------

#' @noRd
getValidNlmixrCtl.adgh <- function(control) {
  if (inherits(control, "adghControl")) return(control)
  .ctl <- control[[1]]
  if (inherits(.ctl, "adghControl")) return(.ctl)
  if (is.list(.ctl) && "studies" %in% names(.ctl))
    return(do.call(adghControl, .ctl[intersect(names(.ctl), names(formals(adghControl)))]))
  if (is.list(control) && length(names(control)) > 0L)
    return(do.call(adghControl, control[intersect(names(control), names(formals(adghControl)))]))
  adghControl()
}

#' @noRd
nmObjHandleControlObject.adghControl <- function(control, env) {
  assign("adghControl", control, envir = env)
}

#' @noRd
nmObjGetControl.adgh <- function(x, ...) {
  .env <- x[[1]]
  for (.nm in c("adghControl", "control")) {
    if (exists(.nm, .env)) {
      .ctl <- get(.nm, .env)
      if (inherits(.ctl, "adghControl")) return(.ctl)
    }
  }
  stop("cannot find adgh control object", call. = FALSE)
}

# -- Main estimation entry point -----------------------------------------------

#' Fit an aggregate data model via Gauss-Hermite quadrature
#'
#' Called automatically by `nlmixr2(model, admData(), est = "adgh",
#' control = adghControl(...))`. Not typically called directly.
#'
#' @param env nlmixr2 environment containing `ui` and `control`.
#' @param ... Unused.
#'
#' @return An `admFit` nlmixr2 fit object.
#'
#' @method nlmixr2Est adgh
#' @importFrom nlmixr2est nlmixr2Est
#' @export
nlmixr2Est.adgh <- function(env, ...) {
  .ui  <- env$ui
  .ctl <- env$control

  if (!inherits(.ctl, "adghControl")) .ctl <- getValidNlmixrCtl.adgh(.ctl)
  if (!inherits(.ctl, "adghControl"))
    stop("Could not recover adghControl", call. = FALSE)
  assign("control", .ctl, envir = .ui)

  .ds     <- .admDriverStudies(.ui, .ctl, "adgh")
  studies <- .ds$studies
  pinfo   <- .ds$pinfo
  output_var <- .admOutputVar(.ui)
  n_nodes    <- .ctl$n_nodes

  .u         <- .admDriverUnits(studies, .ui, output_var)
  studies    <- .u$studies
  multi_out  <- .u$multi_out
  any_joint  <- .u$any_joint


  # RETURNS the studies, annotated with which covariate path each takes.
  # Discarding the value silently disables covariate handling entirely.
  studies <- .admCheckCovariates(.ui, pinfo, studies)
  .admCheckAR(pinfo, studies)
  .admCheckOrdinal(pinfo, studies)
  .admCheckMixedEndpoints(.ui)

  # A beta endpoint's prediction is derived from TWO solved columns; the pair
  # travels on each study so the solve paths can combine them (see .admSimulate).
  .bpair <- .admBetaPair(.ui)
  if (!is.null(.bpair)) {
    studies <- lapply(studies, function(u) { u$out_pair <- .bpair; u })
    # ... and it is fitted DERIVATIVE-FREE, for the reason spelled out in
    # nlmixr2Est.admc(): beta's conditional variance is mu(1-mu)/(1+phi) with the
    # precision phi SOLVED from the structural model, so a theta reaches the
    # objective through phi as well as through mu, and every gradient path here
    # chains through mu alone. BOBYQA differences the objective itself.
    if (.ctl$grad != "none") {
      # Name the ALGORITHM change too. This is the one place an algorithm is
      # chosen outside .admResolveAlgorithm(), and it overrides whatever the user
      # asked for -- reporting only the grad change left an explicit
      # algorithm = "NLOPT_LD_SLSQP" silently replaced by BOBYQA.
      message("adghControl: a beta() endpoint is fitted derivative-free ",
              "(grad = \"none\", algorithm = \"", .admDefaultAlgorithm("none"),
              "\"): its precision is solved from the structural model, and the ",
              "gradient paths carry only d(prediction)/d(theta).")
      .ctl$grad      <- "none"
      .ctl$algorithm <- .admDefaultAlgorithm("none")
    }
  }

  want_grad    <- .ctl$grad != "none"
  want_sens    <- .ctl$grad == "analytical"
  use_pure_fd  <- .ctl$grad == "fd"
  # Joint (same-subject) fits keep the analytical quadrature gradient: .adghGrad's
  # joint branch computes the stacked-MVN gradient from shared-eta per-output
  # sensitivities (grad = "analytical"). grad = "fd" uses the FD gradient.

  if (pinfo$n_eta > 0L) {
    n_total <- n_nodes^pinfo$n_eta
    if (n_total > 5000L)
      message(sprintf(
        "adgh: n_nodes=%d x n_eta=%d = %d nodes. Consider reducing n_nodes or using est='admc'.",
        n_nodes, pinfo$n_eta, n_total))
  }


  # ORDERING INVARIANT: .admLoadSensModel() before .admLoadModel().
  sensModel <- if (want_sens) {
    sm <- tryCatch(.admLoadSensModel(.ui), error = function(e) NULL)
    if (is.null(sm)) {
      warning("adghControl(grad='analytical'): sensitivity model unavailable -- falling back to FD")
      want_sens   <- FALSE
      use_pure_fd <- TRUE
    }
    sm
  } else NULL

  # Unpaired (non-mu-referenced) struct thetas: the sens model carries an explicit
  # THETA_j_ direction for each (.admBuildThetaSens), so their sensitivities come from the same
  # solve as the etas'. Without those columns they fall back to FD of .adghNLL.
  .unpaired <- if (!is.null(pinfo$struct_has_eta))
    names(pinfo$struct_has_eta)[!pinfo$struct_has_eta] else character(0)
  .theta_sens <- want_sens && !is.null(sensModel) &&
    !is.null(sensModel$theta_sens_cols) &&
    all(.unpaired %in% names(sensModel$theta_sens_cols))
  if (length(.unpaired)) {
    message(sprintf("adgh: struct theta(s) without mu-referencing: %s. %s",
                    paste(.unpaired, collapse = ", "),
                    if (.theta_sens) "Sens model carries their sensitivities (no FD)."
                    else "FD for these parameters."))
  }

  rxMod <- .admLoadModel(.ui)
  rxode2::rxLock(rxMod)

  # Free the models this fit registered with rxode2's own idiom (the same
  # gc(); rxUnloadAll() nlmixr2est runs), so many fits in a session stay bounded.
  on.exit({ rxode2::rxUnlock(rxMod); rxode2::rxSolveFree(); gc(FALSE); rxode2::rxUnloadAll() },
          add = TRUE)

  # Node grid: fixed in standard-normal space; L applied per-eval in .adghMoments.
  grid  <- .adghNodeGrid(n_nodes, pinfo$n_eta)

  ov    <- .admBuildOptVec(pinfo)
  cores <- .ctl$cores
  .iter <- 0L

  .nll_trace <- numeric(0)
  .par_trace <- NULL
  .best_nll  <- Inf

  # Measure the gradient's FD steps ONCE, here, and let every later
  # difference reuse them (the mechanism FOCEI's numericGrad uses at nF == 1).
  # Only the parameters actually finite-differenced: all of them under
  # grad = "fd", otherwise the unpaired struct thetas -- but ONLY when the
  # sens model carries no THETA_j_ column for them (`.theta_sens`).
  #
  # An unpaired theta is NOT automatically an FD theta. .admBuildThetaSens emits
  # a direction per unpaired theta, so .adghGrad reads their gradient off the
  # same solve as the etas' and returns BEFORE its FD block; that block is the
  # fallback for when those columns could not be built. Probing them regardless
  # would spend ten NLL evaluations each choosing a step nothing uses.
  # `.ctl$grad_h` stays a scalar unless this fires.
  .fd_idx <- if (!want_grad) integer(0)
    else if (use_pure_fd) seq_along(ov$p0)
    else if (length(.unpaired) && !.theta_sens)
      which(pinfo$struct_names %in% .unpaired)
    else integer(0)
  grad_h_v <- if (length(.fd_idx))
    .admShi21GradH(function(p) .adghNLL(p, pinfo, studies, rxMod, output_var,
                                       grid, cores),
                  ov$p0, .fd_idx, .ctl$grad_h, scaled = TRUE,
                  .var.name = "adgh gradient")
  else .ctl$grad_h

  # (#5) One solve serves both the objective and the gradient -- see .adghFusedFns.
  .fz <- if (want_sens && !is.null(sensModel))
    .adghFusedFns(pinfo, studies, sensModel, rxMod, output_var, grid, cores,
                  grad_h_v) else NULL

  eval_f <- function(p) {
    .iter <<- .iter + 1L
    val <- if (!is.null(.fz)) .fz$nll_fn(p)
           else .adghNLL(p, pinfo, studies, rxMod, output_var, grid, cores)
    if (is.finite(val) && val < .best_nll) {
      .best_nll  <<- val
      .nll_trace <<- c(.nll_trace, val)
      .par_trace <<- rbind(.par_trace, p)
    }
    if (.ctl$print > 0L && .iter %% .ctl$print == 0L) {
      row <- .admProgressRow(sprintf("%04d", .iter), val, p, pinfo)
      if (!is.null(row)) message(row)
    }
    val
  }

  eval_grad_f <- if (!want_grad) {
    NULL
  } else if (use_pure_fd) {
    function(p) .adghFDGrad(p, pinfo, studies, rxMod, output_var, grid, cores,
                              grad_h_v)
  } else if (!is.null(.fz)) {
    .fz$grad_fn
  } else {
    function(p) .adghGrad(p, pinfo, studies, sensModel, rxMod, output_var,
                           grid, cores, grad_h_v)
  }

  grad_label <- if (!want_grad) "none"
                else if (!is.null(sensModel)) "Analytical"
                else "CFD"
  n_total_nodes <- if (pinfo$n_eta > 0L) n_nodes^pinfo$n_eta else 1L
  message("=== admixr2: Aggregate Data Modeling (GH) ===")
  message(sprintf("  Obs units: %d | Params: %d | Nodes: %d^%d=%d | Cores: %d | Grad: %s | Restarts: %d",
                  length(studies), length(ov$p0),
                  n_nodes, pinfo$n_eta, n_total_nodes,
                  cores, grad_label, .ctl$n_restarts))
  t0 <- proc.time()

  lb <- if (want_grad) pmax(ov$lower, ov$p0 - .ctl$grad_bounds) else ov$lower
  ub <- if (want_grad) pmin(ov$upper, ov$p0 + .ctl$grad_bounds) else ov$upper

  sc           <- ov$scale_c
  p0_sc        <- ov$p0 / sc
  lb_sc        <- lb  / sc
  ub_sc        <- ub  / sc
  eval_f_sc    <- function(p_s) eval_f(p_s * sc)
  eval_grad_sc <- if (!is.null(eval_grad_f)) {
    function(p_s) eval_grad_f(p_s * sc) * sc
  } else NULL

  if (.ctl$n_restarts == 1L) {
    message(.admProgressHeader(pinfo))
    opt_raw <- nlmixr2est::nlmixrWithTiming("adgh", {
      nloptr::nloptr(x0 = p0_sc, eval_f = eval_f_sc,
                     eval_grad_f = eval_grad_sc,
                     lb = lb_sc, ub = ub_sc,
                     opts = list(algorithm = .ctl$algorithm,
                                 ftol_rel  = .ctl$ftol_rel,
                                 maxeval   = .ctl$maxeval))
    })
    opt <- list(objective  = opt_raw$objective,
                solution   = opt_raw$solution * sc,
                message    = opt_raw$message,
                all_traces = list(list(restart_id = 1L,
                                       nll_trace  = .nll_trace,
                                       par_trace  = .par_trace)))
    if (.ctl$print > 0L) {
      row <- .admProgressRow(sprintf("%04d \u2713", .iter), opt$objective, opt$solution, pinfo)
      if (!is.null(row)) message(paste0(row, "\n",
        .admProgressTimingRow((proc.time() - t0)["elapsed"], pinfo)))
    }
  } else {
    .admSetupDaemons(.ctl, .ctl$n_restarts)
    on.exit(.admStopDaemons(), add = TRUE)
    opt <- .admRunRestarts(
      worker_fn  = .adghRestartWorker,
      p0         = ov$p0, ov = ov, pinfo = pinfo,
      .ctl       = .ctl, ui = .ui, studies = studies,
      extra_args = list(
        n_nodes          = n_nodes,
        algorithm        = .ctl$algorithm,
        ftol_rel         = .ctl$ftol_rel,
        maxeval          = .ctl$maxeval,
        use_grad         = want_grad,
        use_pure_fd      = use_pure_fd,
        # The measured steps, not the constant -- restarts must difference the
        # same way the sequential path does.
        grad_h           = grad_h_v,
        grad_bounds      = .ctl$grad_bounds,
        output_var       = output_var,
        print_progress   = TRUE,
        print            = .ctl$print,
        cores            = .ctl$cores,
        rxMod_direct     = rxMod,
        sensModel_direct = sensModel
      )
    )
    .admStopDaemons()
    .iter <- opt$n_iter
  }

  t_opt  <- (proc.time() - t0)["elapsed"]
  # A gradient fit is confined to <box centre> +/- grad_bounds; say so if it
  # stopped there rather than at an interior optimum. The centre is the winning
  # RESTART's own init where there was one -- see .admScaledOptimize()'s box_centre.
  if (want_grad)
    .admWarnOnBounds(opt$solution, opt$box_centre %||% ov$p0, ov,
                     .ctl$grad_bounds, pinfo)
  final  <- .admUnpack(opt$solution, pinfo)
  fullTheta <- .admFullTheta(final, pinfo)
  p_hat  <- setNames(opt$solution, names(ov$p0))

  t0_cov <- proc.time()
  .want_cov <- .ctl$covMethod %in% c("r", "r,s")
  .cov <- if (.want_cov) {
    # struct + sigma + OMEGA: the Hessian spans all three, so the evaluation
    # count must too.
    np_cov    <- length(pinfo$struct_names) + length(pinfo$sigma_names) +
                 length(pinfo$omega_par)
    use_grad_cov <- want_grad && !is.null(sensModel)
    n_evals   <- if (use_grad_cov) np_cov + 1L
                 else { n_off <- np_cov * (np_cov - 1L) / 2L; 1L + 2L * np_cov + 4L * n_off }
    evals_lbl <- if (use_grad_cov) "gradient evaluations" else "NLL evaluations"
    hess_lbl  <- if (!use_grad_cov) "" else if (!is.null(sensModel)) ", Analytical-Hessian" else ", FD-Hessian"
    sw_lbl    <- if (.ctl$covMethod == "r,s") ", sandwich" else ""
    message(sprintf("  Computing covariance (R method%s%s, %d %s)",
                    hess_lbl, sw_lbl, n_evals, evals_lbl))
    tryCatch(
      .adghCalcCov(p_hat, pinfo, studies, sensModel, rxMod, output_var, grid, cores,
                   use_grad    = use_grad_cov,
                   grad_h      = .ctl$cov_h,
                   cov_h_outer = .ctl$cov_h_outer,
                   sandwich    = .ctl$covMethod == "r,s"),
      error = function(e) { warning("adghCalcCov failed: ", conditionMessage(e)); NULL })
  } else NULL
  # A NULL covariance used to be completely silent: no warning reached the user,
  # `warnings()` was empty, covMethod came back "" and every SE was NA with no
  # indication why. Say so once, from the driver, where it cannot be swallowed.
  if (.want_cov && is.null(.cov))
    warning("covariance could not be computed (the Hessian was singular or ",
            "non-finite); standard errors are unavailable for this fit.",
            call. = FALSE)
  # The attribute records what the covariance IS, not what was asked for: a
  # requested sandwich that degraded must not be reported as one.
  .cov_lbl  <- if (isTRUE(attr(.cov, "sandwich"))) "r,s" else "r"
  .sw_HJ    <- attr(.cov, "sandwich_HJ")
  # iniDf order first (nlmixr2est maps SEs positionally), then snapshot the names
  # BEFORE nlmixr2est sees it -- .admCovThetaOrder()/.admRestoreCovNames().
  .cov      <- .admCovThetaOrder(.cov, .ui)
  .cov_nms  <- .admCovNames(.cov)
  t_cov     <- (proc.time() - t0_cov)["elapsed"]
  t_elapsed <- t_opt + t_cov

  if (.ctl$returnAdmr) {
    return(list(objective = opt$objective, fullTheta = fullTheta,
                struct = final$struct, sigma_var = final$sigma_var,
                omega = final$omega, L = final$L, nloptr = opt,
                cov = .cov))
  }

  .ret            <- new.env(parent = emptyenv())
  .ret$table      <- env$table
  .ret$ui         <- .ui
  .ret$fullTheta  <- fullTheta
  .ret$objective  <- opt$objective
  .ret$est        <- "adgh"
  .ret$ofvType    <- "adgh"
  .ret$adjObf     <- FALSE
  .ret$covMethod  <- if (!is.null(.cov)) .cov_lbl else ""
  .ret$cov        <- .cov
  .ret$message    <- opt$message
  .ret$extra      <- ""
  .ret$origData   <- studies

  .ret$admExtra <- list(sandwich = .sw_HJ,
                        struct         = final$struct,
                        sigma_var      = final$sigma_var,
                        sigma_is_prop  = pinfo$sigma_is_prop,
                        sigma_is_lnorm = pinfo$sigma_is_lnorm,
                        # the TBS residual quadrature the FIT used -- see adfo.R
                        resid_nodes    = pinfo$resid_nodes,
                        # ... and the solver tolerance the FIT used, for the
                        # same reason: plot.admFit() re-solves the model to build
                        # the predicted mean and covariance panels, and a fit run
                        # at a looser sigdig diagnosed against rxode2's stock
                        # tolerances shows standardised-residual structure the
                        # objective was never minimised on.
                        sigdig         = pinfo$sigdig,
                        omega          = final$omega,
                        L              = final$L,
                        eta_col_names  = pinfo$eta_col_names,
                        par_names      = names(ov$p0),
                        npar           = length(ov$p0),
                        nloptr         = opt,
                        nll_trace      = .nll_trace,
                        par_trace      = .par_trace,
                        all_traces     = opt$all_traces,
                        n_iter         = .iter,
                        time           = t_elapsed,
                        t_opt          = t_opt,
                        t_cov          = t_cov,
                        studies        = studies,
                        n_nodes        = n_nodes,
                        n_sim          = 5000L,
                        sampling       = "sobol",
                        n_gh           = n_total_nodes)

  .admFinaliseFit(.ret, .ui, .ctl, est = "adgh", objective = opt$objective,
                  ov = ov, studies = studies, cov = .cov,
                  cov_nms = .cov_nms, multi_out = multi_out,
                  extra_field = "admExtra",
                  handle_ctl = nmObjHandleControlObject.adghControl,
                  t_opt = t_opt, t_cov = t_cov, t_elapsed = t_elapsed,
                  pinfo = pinfo)
}
