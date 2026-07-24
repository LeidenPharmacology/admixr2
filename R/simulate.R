# Single rxSolve pass for one study given pre-computed eta_mat (n_sim x n_eta).
# Returns n_sim x n_times matrix of predicted concentrations.
# params_mat is a named numeric matrix (from .admMakeParamsList); converted to
# data.frame only at the rxSolve call to avoid repeated list COW copies.
.admSimulate <- function(rxMod, struct_theta, sigma_names, eta_mat, study,
                         output_var, params_mat, cores,
                         ndp = .Machine$integer.max) {
  eta_cols <- colnames(eta_mat)
  for (nm in names(struct_theta)) params_mat[, nm] <- struct_theta[nm]
  if (length(eta_cols) > 0L)      params_mat[, eta_cols] <- eta_mat
  for (nm in sigma_names)         params_mat[, nm] <- 0
  # Only the parameters we vary are supplied; rxSolve fills the rest (rxerr.*,
  # CMT, hard-coded model constants) from the model's own defaults.
  out  <- rxode2::rxSolve(rxMod, params = as.data.frame(params_mat),
                          events = study$ev_full, cores = cores,
                          nDisplayProgress = ndp)
  keep <- out[["time"]] %in% study$times
  # A beta endpoint is defined by TWO solved columns, so its prediction is
  # derived: mu = b1/(b1+b2). study$out_pair carries their names; every other
  # endpoint reads a single column exactly as before. Inlined rather than
  # factored out -- see the dev-mode daemon note at the top of this file.
  .phi <- NULL
  vals <- if (!is.null(study$out_pair)) {
    .b1 <- out[[study$out_pair[[1L]]]]; .b2 <- out[[study$out_pair[[2L]]]]
    .phi <- .b1 + .b2                      # precision; needed for the variance
    # Guard the denominator exactly as the sibling solve paths do (.admSimulateRows,
    # .admSimulateJoint, and the admc inlined copies): if phi = b1 + b2 hits 0 at a
    # draw, an unguarded b1/phi is 0/0 = NaN and poisons the whole moment/objective.
    .b1 / { .d <- .phi; .d[.d == 0] <- .Machine$double.eps; .d }
  } else {
    .v <- out[[output_var]]                # linCmt yields "ipredSim", not rx_pred_
    if (is.null(.v)) out[["ipredSim"]] else .v
  }
  m <- matrix(vals[keep],
              nrow = nrow(eta_mat), ncol = length(study$times), byrow = TRUE)
  # beta's Var(y|eta) = mu(1-mu)/(1+phi) needs phi, which is SOLVED rather than a
  # residual parameter -- so it rides back as an attribute. A matrix with an extra
  # attribute is still a matrix, so every existing consumer is unaffected; the one
  # caller that needs it reads it immediately after this returns. phi must be
  # eta-independent for the aggregate variance to factor; .admBetaPhiConst()
  # VERIFIES that across the draws and returns the representative row, so the
  # assumption and its use stay together (it used to be asserted by a comment only).
  if (!is.null(.phi))
    attr(m, "phi") <- .admBetaPhiConst(
      matrix(.phi[keep], nrow = nrow(eta_mat),
             ncol = length(study$times), byrow = TRUE))
  m
}

# Row-varying variants of .admSimulate / .admSimulateSens.
#
# .admSimulate broadcasts ONE structural-theta vector across every row. These
# take a struct_mat (n_row x n_struct, natural scale, colnames = theta names)
# so each row can carry its own thetas. That is what lets the FO/GH estimators
# put a whole set of finite-difference directions into a SINGLE rxSolve: an
# rxSolve call costs ~11 ms before it does any work, and FO's solves are one
# subject each, so the call overhead -- not the integration -- was the cost.
#
# eta_mat is n_row x n_eta and lines up row-for-row with struct_mat.
.admSimulateRows <- function(rxMod, struct_mat, sigma_names, eta_mat, study,
                             output_var, params_mat, cores,
                             ndp = .Machine$integer.max) {
  eta_cols <- colnames(eta_mat)
  for (nm in colnames(struct_mat)) params_mat[, nm] <- struct_mat[, nm]
  if (length(eta_cols) > 0L)       params_mat[, eta_cols] <- eta_mat
  for (nm in sigma_names)          params_mat[, nm] <- 0
  out  <- rxode2::rxSolve(rxMod, params = as.data.frame(params_mat),
                          events = study$ev_full, cores = cores,
                          nDisplayProgress = ndp)
  keep <- out[["time"]] %in% study$times
  # beta: derived prediction mu = b1/(b1+b2) -- see .admSimulate
  .phi <- NULL
  vals <- if (!is.null(study$out_pair)) {
    .b1 <- out[[study$out_pair[[1L]]]]; .b2 <- out[[study$out_pair[[2L]]]]
    .phi <- .b1 + .b2
    .b1 / { .d <- .phi; .d[.d == 0] <- .Machine$double.eps; .d }
  } else {
    .v <- out[[output_var]]
    if (is.null(.v)) out[["ipredSim"]] else .v
  }
  m <- matrix(vals[keep], nrow = nrow(struct_mat), ncol = length(study$times),
              byrow = TRUE)
  # phi rides back the same way it does from .admSimulate: beta's
  # Var(y | eta) = mu(1 - mu)/(1 + phi) needs it, and it is SOLVED rather than
  # fitted, so a caller that only has the prediction matrix cannot recover it.
  # Without this the row-varying paths left arr$phi at NA and every entry of the
  # predicted covariance came back NA. phi is eta-independent (verified by
  # .admBetaPhiConst() on the eta-draw paths) but NOT theta-independent, so each
  # ROW keeps its own -- these paths exist precisely to put different thetas in
  # different rows, so the rows here are CONFIGURATIONS, not draws, and must not
  # be collapsed or passed through the eta-independence check.
  if (!is.null(.phi))
    attr(m, "phi") <- matrix(.phi[keep], nrow = nrow(struct_mat),
                             ncol = length(study$times), byrow = TRUE)
  m
}

# NOTE ON THE FIXED-THETA FILL BELOW (repeated inline in three solve paths rather
# than factored into a helper -- deliberately):
#
# A FIXED theta is not an estimated parameter, so nothing in pinfo (and nothing in
# the solve paths) writes its THETA[k] column -- but the sens model still HAS that
# slot and rxSolve REQUIRES every model parameter. Left unset, the sens solve
# errors and returns NULL: admc/adfo then silently drop to a finite-difference
# gradient, and .adghGrad silently skipped the study altogether.
#
# These three loops run INSIDE mirai daemons. utils::assignInNamespace() can
# REPLACE a binding in the locked installed namespace but cannot ADD one, so a new
# helper called from here would be missing in a dev-mode daemon (devtools::load_all
# with workers > 1 and no prior devtools::install) -- "could not find function" --
# and abort the whole fit. Same hazard CLAUDE.md records for .admRestartWorker's
# signature. Three copies of two lines is the cheaper price.

.admSimulateSensRows <- function(sensModel, struct_mat, sigma_names, eta_mat, study,
                                 cores, ndp = .Machine$integer.max,
                                 sigma_var = NULL) {
  rmap      <- sensModel$rename_map
  n_row     <- nrow(struct_mat)
  theta_nms <- colnames(struct_mat)
  eta_cols  <- colnames(eta_mat)

  all_src   <- c(theta_nms, sigma_names, eta_cols)
  inner_nms <- rmap[all_src]
  inner_nms <- inner_nms[!is.na(inner_nms)]

  inner_df  <- as.data.frame(matrix(0, nrow = n_row, ncol = length(inner_nms),
                                    dimnames = list(NULL, unname(inner_nms))),
                              check.names = FALSE)
  for (nm in theta_nms) {
    mapped <- rmap[nm]; if (!is.na(mapped)) inner_df[[mapped]][] <- struct_mat[, nm]
  }
  for (j in seq_along(eta_cols)) {
    mapped <- rmap[eta_cols[j]]; if (!is.na(mapped)) inner_df[[mapped]][] <- eta_mat[, j]
  }
  # fixed thetas: constants the estimated-parameter loops above never write (see
  # the note at the top of this file -- inlined on purpose, no helper)
  for (nm in names(sensModel$fixed_theta))
    inner_df[[nm]] <- rep(unname(sensModel$fixed_theta[[nm]]), nrow(inner_df))

  # An ESTIMATED boxCox/yeoJohnson lambda is a SIGMA name, and the zero-fill above
  # therefore handed the solve lambda = 0 -- so rx_pred_ came back as a plain log
  # transform while the back-transform below inverted with the model's STARTING
  # lambda. Two different transforms, hence a sens gradient that was ~60x wrong for
  # boxCox and NaN for yeoJohnson. Write the current lambda into the solve and
  # invert with that same number, so the two agree by construction.
  # Inlined in each solve path on purpose -- see the dev-mode daemon note above.
  .tb <- sensModel$pred_tbs
  .lam <- if (is.null(.tb)) NA_real_ else .tb$lam
  if (!is.null(.tb) && !is.na(.tb$lam_name %||% NA_character_) &&
      !is.null(sigma_var) && .tb$lam_name %in% names(sigma_var)) {
    .lam   <- unname(sigma_var[[.tb$lam_name]])
    .mapped <- rmap[.tb$lam_name]
    if (!is.na(.mapped) && .mapped %in% names(inner_df)) inner_df[[.mapped]][] <- .lam
  }

  # do.call + sensModel$solve_args: a DDE model's sensitivity solve is forced onto
  # pure dop853 (see .admLoadSensModel). solve_args is NULL for every ordinary
  # model, and c(list(...), NULL) is the original list, so nothing else changes.
  # Spliced inline rather than through a helper -- see the note at the top of this
  # file about dev-mode daemons and new functions.
  out <- tryCatch(
    suppressWarnings(
      do.call(rxode2::rxSolve,
              c(list(sensModel$mod, params = inner_df,
                     events = study$ev_full, cores = cores,
                     nDisplayProgress = ndp), sensModel$solve_args))),
    error = function(e) NULL)
  if (is.null(out)) return(NULL)
  if (!all(sensModel$sens_cols %in% names(out))) return(NULL)

  keep  <- out[["time"]] %in% study$times
  n_t   <- length(study$times)
  n_eta <- ncol(eta_mat)

  cp_mat     <- matrix(out[["rx_pred_"]][keep], nrow = n_row, ncol = n_t, byrow = TRUE)
  dpred_list <- lapply(seq_len(n_eta), function(j)
    matrix(out[[sensModel$sens_cols[j]]][keep], nrow = n_row, ncol = n_t, byrow = TRUE))
  dtheta_list <- .admThetaSens(sensModel, out, keep, n_row, n_t)

  # lnorm endpoint: rx_pred_ is log(f). Back-transform to the natural scale the
  # NLL works on, chaining every sensitivity by d(exp(g))/dp = exp(g)*dg/dp.
  # Inlined rather than factored out -- see the dev-mode daemon note at the top.
  # Transformed endpoint: rx_pred_ is on the MODELLING scale. Back-transform and
  # chain every sensitivity by g'(z) -- d(g(z))/dp = g'(z) * dz/dp. Covers lnorm
  # (yj = 0, lambda = 0) and logit/probit/boxCox/yeoJohnson alike.
  if (!is.null(.tb)) {
    .gp        <- .admTBSid(cp_mat, .lam, .tb$yj, .tb$lo, .tb$hi)
    cp_mat     <- .admTBSi(cp_mat, .lam, .tb$yj, .tb$lo, .tb$hi)
    dpred_list <- lapply(dpred_list, function(D) D * .gp)
    if (!is.null(dtheta_list)) dtheta_list <- lapply(dtheta_list, function(D) D * .gp)
  }

  list(cp_mat = cp_mat, dpred_list = dpred_list, dtheta_list = dtheta_list)
}

# Joint (same-subject) simulation: one rxSolve with SHARED eta produces every
# observed output; each block is extracted by output name at its own times and
# stacked column-wise into an n_sim x n_total matrix (columns in block/row
# order). Used for joint units where the compartments share random effects.
.admSimulateJoint <- function(rxMod, struct_theta, sigma_names, eta_mat, unit,
                              params_mat, cores, ndp = .Machine$integer.max) {
  eta_cols <- colnames(eta_mat)
  for (nm in names(struct_theta)) params_mat[, nm] <- struct_theta[nm]
  if (length(eta_cols) > 0L)      params_mat[, eta_cols] <- eta_mat
  for (nm in sigma_names)         params_mat[, nm] <- 0
  out  <- rxode2::rxSolve(rxMod, params = as.data.frame(params_mat),
                          events = unit$ev_full, cores = cores,
                          nDisplayProgress = ndp)
  n_sim <- nrow(eta_mat)
  cp    <- matrix(0, nrow = n_sim, ncol = unit$n_total)
  time  <- out[["time"]]
  for (blk in unit$blocks) {
    vals <- if (!is.null(blk$out_pair)) {
      .b1 <- out[[blk$out_pair[[1L]]]]; .b2 <- out[[blk$out_pair[[2L]]]]
      .b1 / { .d <- .b1 + .b2; .d[.d == 0] <- .Machine$double.eps; .d }
    } else out[[blk$output]]
    if (is.null(vals)) vals <- out[["ipredSim"]]
    keep <- time %in% blk$times
    cp[, blk$rows] <- matrix(vals[keep], nrow = n_sim,
                             ncol = length(blk$times), byrow = TRUE)
  }
  cp
}

# Joint sensitivity simulation: for each observed output, one sens solve with the
# SHARED eta draws (obs tagged with that output's cmt) gives its prediction and
# d(pred)/d(eta_j); stacked column-wise into n_sim x n_total matrices. Returns
# list(cp_mat, dpred_list) or NULL if any block's sens solve fails (caller then
# falls back to FD). Enables the analytical gradient of a joint (same-subject)
# unit's stacked MVN.
.admSimulateJointSens <- function(sensModel, struct, sigma_names, eta_mat, unit,
                                  cores, ndp = .Machine$integer.max,
                                  sigma_var = NULL) {
  n_sim <- nrow(eta_mat); n_eta <- ncol(eta_mat)
  th_nms <- names(sensModel$theta_sens_cols)
  cp_mat     <- matrix(0, n_sim, unit$n_total)
  dpred_list <- lapply(seq_len(n_eta), function(j) matrix(0, n_sim, unit$n_total))
  dtheta_list <- if (length(th_nms))
    stats::setNames(lapply(th_nms, function(nm) matrix(0, n_sim, unit$n_total)), th_nms)
  else NULL
  for (blk in unit$blocks) {
    bs  <- list(ev_full = blk$ev_full, times = blk$times)
    res <- .admSimulateSens(sensModel, struct, sigma_names, eta_mat, bs, cores, ndp,
                            sigma_var)
    if (is.null(res)) return(NULL)
    cp_mat[, blk$rows] <- res$cp_mat
    for (j in seq_len(n_eta)) dpred_list[[j]][, blk$rows] <- res$dpred_list[[j]]
    # any block missing its theta columns disables the theta path for the whole
    # joint unit (a half-filled stacked derivative would be silently wrong)
    if (!is.null(dtheta_list)) {
      if (is.null(res$dtheta_list)) dtheta_list <- NULL
      else for (nm in th_nms) dtheta_list[[nm]][, blk$rows] <- res$dtheta_list[[nm]]
    }
  }
  list(cp_mat = cp_mat, dpred_list = dpred_list, dtheta_list = dtheta_list)
}

# Extract d(pred)/d(theta) for the unpaired thetas from a sens solve.
# NULL when the sens model has no theta directions (the nlmixr2est-inner
# fallback) or the solve did not return them -- the caller then falls back to
# finite differences.
.admThetaSens <- function(sensModel, out, keep, n_row, n_t) {
  tsc <- sensModel$theta_sens_cols
  if (is.null(tsc) || length(tsc) == 0L) return(NULL)
  if (!all(tsc %in% names(out))) return(NULL)
  stats::setNames(
    lapply(tsc, function(col)
      matrix(out[[col]][keep], nrow = n_row, ncol = n_t, byrow = TRUE)),
    names(tsc))
}

# Single pass on the sensitivity model returning predictions, d(pred)/d(eta_j)
# and (augmented model only) d(pred)/d(theta_k) for the unpaired struct thetas.
# Returns list(cp_mat, dpred_list, dtheta_list) or NULL on failure (caller falls
# back to FD). dtheta_list is NULL when the model carries no theta directions.
.admSimulateSens <- function(sensModel, struct_theta, sigma_names,
                             eta_mat, study, cores,
                             ndp = .Machine$integer.max, sigma_var = NULL) {
  eta_cols  <- colnames(eta_mat)
  rmap      <- sensModel$rename_map
  n_sim     <- nrow(eta_mat)
  theta_nms <- names(struct_theta)

  all_src   <- c(theta_nms, sigma_names, eta_cols)
  inner_nms <- rmap[all_src]
  inner_nms <- inner_nms[!is.na(inner_nms)]

  # check.names=FALSE: preserve THETA[1]/ETA[1] bracket notation so column
  # assignments below find existing columns rather than creating duplicates.
  inner_df  <- as.data.frame(matrix(0, nrow = n_sim, ncol = length(inner_nms),
                                    dimnames = list(NULL, unname(inner_nms))),
                              check.names = FALSE)
  for (nm in theta_nms) {
    mapped <- rmap[nm]; if (!is.na(mapped)) inner_df[[mapped]][] <- struct_theta[nm]
  }
  for (j in seq_along(eta_cols)) {
    mapped <- rmap[eta_cols[j]]; if (!is.na(mapped)) inner_df[[mapped]][] <- eta_mat[, j]
  }
  # fixed thetas: constants the estimated-parameter loops above never write (see
  # the note at the top of this file -- inlined on purpose, no helper)
  for (nm in names(sensModel$fixed_theta))
    inner_df[[nm]] <- rep(unname(sensModel$fixed_theta[[nm]]), nrow(inner_df))

  # An ESTIMATED boxCox/yeoJohnson lambda is a SIGMA name, and the zero-fill above
  # therefore handed the solve lambda = 0 -- so rx_pred_ came back as a plain log
  # transform while the back-transform below inverted with the model's STARTING
  # lambda. Two different transforms, hence a sens gradient that was ~60x wrong for
  # boxCox and NaN for yeoJohnson. Write the current lambda into the solve and
  # invert with that same number, so the two agree by construction.
  # Inlined in each solve path on purpose -- see the dev-mode daemon note above.
  .tb <- sensModel$pred_tbs
  .lam <- if (is.null(.tb)) NA_real_ else .tb$lam
  if (!is.null(.tb) && !is.na(.tb$lam_name %||% NA_character_) &&
      !is.null(sigma_var) && .tb$lam_name %in% names(sigma_var)) {
    .lam   <- unname(sigma_var[[.tb$lam_name]])
    .mapped <- rmap[.tb$lam_name]
    if (!is.na(.mapped) && .mapped %in% names(inner_df)) inner_df[[.mapped]][] <- .lam
  }

  # do.call + sensModel$solve_args: DDE sensitivity solves are forced onto pure
  # dop853 (see .admLoadSensModel); NULL, hence a no-op, for every other model.
  out <- tryCatch(
    suppressWarnings(
      do.call(rxode2::rxSolve,
              c(list(sensModel$mod, params = inner_df,
                     events = study$ev_full, cores = cores,
                     nDisplayProgress = ndp), sensModel$solve_args))),
    error = function(e) NULL)
  if (is.null(out)) return(NULL)

  out_cols <- names(out)
  if (!all(sensModel$sens_cols %in% out_cols)) return(NULL)

  keep  <- out[["time"]] %in% study$times
  n_t   <- length(study$times)
  n_eta <- ncol(eta_mat)

  cp_mat     <- matrix(out[["rx_pred_"]][keep], nrow = n_sim, ncol = n_t, byrow = TRUE)
  dpred_list <- lapply(seq_len(n_eta), function(j)
    matrix(out[[sensModel$sens_cols[j]]][keep], nrow = n_sim, ncol = n_t, byrow = TRUE))
  dtheta_list <- .admThetaSens(sensModel, out, keep, n_sim, n_t)

  # lnorm endpoint: rx_pred_ is log(f) -- back-transform with the chain rule so
  # the gradient differentiates the same quantity the NLL scores.
  # Transformed endpoint: rx_pred_ is on the MODELLING scale. Back-transform and
  # chain every sensitivity by g'(z) -- d(g(z))/dp = g'(z) * dz/dp. Covers lnorm
  # (yj = 0, lambda = 0) and logit/probit/boxCox/yeoJohnson alike.
  if (!is.null(.tb)) {
    .gp        <- .admTBSid(cp_mat, .lam, .tb$yj, .tb$lo, .tb$hi)
    cp_mat     <- .admTBSi(cp_mat, .lam, .tb$yj, .tb$lo, .tb$hi)
    dpred_list <- lapply(dpred_list, function(D) D * .gp)
    if (!is.null(dtheta_list)) dtheta_list <- lapply(dtheta_list, function(D) D * .gp)
  }

  list(cp_mat = cp_mat, dpred_list = dpred_list, dtheta_list = dtheta_list)
}
