# `sigdig` passes pinfo$sigdig to rxSolve's `sigdig` (NULL = rxode2 defaults).
# Last formal with NULL default preserves backward compatibility for callers.
# Single rxSolve pass for one study given pre-computed eta_mat (n_sim x n_eta).
# Returns n_sim x n_times matrix of predicted concentrations.
.admSimulate <- function(rxMod, struct_theta, sigma_names, eta_mat, study,
                         output_var, params_mat, cores,
                         ndp = .Machine$integer.max, sigdig = NULL) {
  eta_cols <- colnames(eta_mat)
  for (nm in names(struct_theta)) params_mat[, nm] <- struct_theta[nm]
  if (length(eta_cols) > 0L)      params_mat[, eta_cols] <- eta_mat
  for (nm in sigma_names)         params_mat[, nm] <- 0
  # Supply varied parameters; rxSolve fills model defaults. Add study covariates.
  params_mat <- .admCovCols(params_mat, rxMod$params, study[["cov"]], study[["cov_rows"]])
  out  <- rxode2::rxSolve(rxMod, params = as.data.frame(params_mat),
                          events = study$ev_full, cores = cores,
                          nDisplayProgress = ndp,
                          sigdig = sigdig)
  keep <- out[["time"]] %in% study$times
  # Beta endpoints derive prediction as mu = b1 / (b1 + b2).
  .phi <- NULL
  vals <- if (!is.null(study$out_pair)) {
    .b1 <- out[[study$out_pair[[1L]]]]; .b2 <- out[[study$out_pair[[2L]]]]
    .phi <- .b1 + .b2                      # precision; needed for the variance
    # Guard against division by zero when b1 + b2 == 0.
    .b1 / { .d <- .phi; .d[.d == 0] <- .Machine$double.eps; .d }
  } else {
    .v <- out[[output_var]]                # linCmt yields "ipredSim", not rx_pred_
    if (is.null(.v)) out[["ipredSim"]] else .v
  }
  m <- matrix(vals[keep],
              nrow = nrow(eta_mat), ncol = length(study$times), byrow = TRUE)
  # Attach solved beta precision phi = b1 + b2 as attribute (verified eta-independent).
  if (!is.null(.phi))
    attr(m, "phi") <- .admBetaPhiConst(
      matrix(.phi[keep], nrow = nrow(eta_mat),
             ncol = length(study$times), byrow = TRUE))
  m
}

# Row-varying variants taking struct_mat (n_row x n_struct) and eta_mat (n_row x n_eta)
# to evaluate multiple parameter configurations in a single rxSolve call.
.admSimulateRows <- function(rxMod, struct_mat, sigma_names, eta_mat, study,
                             output_var, params_mat, cores,
                             ndp = .Machine$integer.max, sigdig = NULL) {
  eta_cols <- colnames(eta_mat)
  for (nm in colnames(struct_mat)) params_mat[, nm] <- struct_mat[, nm]
  if (length(eta_cols) > 0L)       params_mat[, eta_cols] <- eta_mat
  for (nm in sigma_names)          params_mat[, nm] <- 0
  # Populate study covariates into parameter matrix.
  params_mat <- .admCovCols(params_mat, rxMod$params, study[["cov"]],
                            study[["cov_rows"]])
  out  <- rxode2::rxSolve(rxMod, params = as.data.frame(params_mat),
                          events = study$ev_full, cores = cores,
                          nDisplayProgress = ndp,
                          sigdig = sigdig)
  keep <- out[["time"]] %in% study$times
  # beta: derived prediction mu = b1/(b1+b2)
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
  # Attach solved beta precision matrix (one row per configuration).
  if (!is.null(.phi))
    attr(m, "phi") <- matrix(.phi[keep], nrow = nrow(struct_mat),
                             ncol = length(study$times), byrow = TRUE)
  m
}

# Fixed-theta values are filled inline across daemon-executed solve loops
# to avoid adding unexported helper functions in worker namespaces.

.admSimulateSensRows <- function(sensModel, struct_mat, sigma_names, eta_mat, study,
                                 cores, ndp = .Machine$integer.max,
                                 sigma_var = NULL, sigdig = NULL) {
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
  # fixed thetas: constants the estimated-parameter loops above never write
  for (nm in names(sensModel$fixed_theta))
    inner_df[[nm]] <- rep(unname(sensModel$fixed_theta[[nm]]), nrow(inner_df))

  # Supply current lambda estimate for TBS sensitivity solves.
  .tb <- sensModel$pred_tbs
  .lam <- if (is.null(.tb)) NA_real_ else .tb$lam
  if (!is.null(.tb) && !is.na(.tb$lam_name %||% NA_character_) &&
      !is.null(sigma_var) && .tb$lam_name %in% names(sigma_var)) {
    .lam   <- unname(sigma_var[[.tb$lam_name]])
    .mapped <- rmap[.tb$lam_name]
    if (!is.na(.mapped) && .mapped %in% names(inner_df)) inner_df[[.mapped]][] <- .lam
  }

  # Populate model covariates for sensitivity solve.
  inner_df <- .admCovCols(inner_df, sensModel$mod$params, study[["cov"]],
                          study[["cov_rows"]])
  # Forward solve_args (e.g. forced dop853 for DDE sensitivity models).
  out <- tryCatch(
    suppressWarnings(
      do.call(rxode2::rxSolve,
              c(list(sensModel$mod, params = inner_df,
                     events = study$ev_full, cores = cores,
                     nDisplayProgress = ndp,
                     sigdig = sigdig), sensModel$solve_args))),
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

  # Second-order cross block d2(pred)/(d eta_i d dir) for order-2 sens models;
  # dropped for transformed endpoints which use finite differences.
  d2_list <- NULL
  if (!is.null(sensModel$d2_cols) && all(sensModel$d2_cols %in% names(out))) {
    d2_list <- lapply(seq_len(ncol(sensModel$d2_cols)), function(b)
      lapply(seq_len(nrow(sensModel$d2_cols)), function(i)
        matrix(out[[sensModel$d2_cols[i, b]]][keep], nrow = n_row, ncol = n_t,
               byrow = TRUE)))
    names(d2_list) <- colnames(sensModel$d2_cols)
  }

  # Back-transform TBS/lnorm predictions and chain sensitivities by g'(z).
  if (!is.null(.tb)) {
    .gp        <- .admTBSid(cp_mat, .lam, .tb$yj, .tb$lo, .tb$hi)
    cp_mat     <- .admTBSi(cp_mat, .lam, .tb$yj, .tb$lo, .tb$hi)
    dpred_list <- lapply(dpred_list, function(D) D * .gp)
    if (!is.null(dtheta_list)) dtheta_list <- lapply(dtheta_list, function(D) D * .gp)
    d2_list    <- NULL                      # see the note above
  }

  list(cp_mat = cp_mat, dpred_list = dpred_list, dtheta_list = dtheta_list,
       d2_list = d2_list)
}

# Joint (same-subject) simulation: one rxSolve with shared eta produces every
# observed output; columns stacked in block/row order.
.admSimulateJoint <- function(rxMod, struct_theta, sigma_names, eta_mat, unit,
                              params_mat, cores, ndp = .Machine$integer.max,
                              sigdig = NULL) {
  eta_cols <- colnames(eta_mat)
  for (nm in names(struct_theta)) params_mat[, nm] <- struct_theta[nm]
  if (length(eta_cols) > 0L)      params_mat[, eta_cols] <- eta_mat
  for (nm in sigma_names)         params_mat[, nm] <- 0
  # Forward covariates to joint solve.
  params_mat <- .admCovCols(params_mat, rxMod$params, unit[["cov"]],
                            unit[["cov_rows"]])
  out  <- rxode2::rxSolve(rxMod, params = as.data.frame(params_mat),
                          events = unit$ev_full, cores = cores,
                          nDisplayProgress = ndp,
                          sigdig = sigdig)
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
                                  sigma_var = NULL, sigdig = NULL) {
  n_sim <- nrow(eta_mat); n_eta <- ncol(eta_mat)
  th_nms <- names(sensModel$theta_sens_cols)
  cp_mat     <- matrix(0, n_sim, unit$n_total)
  dpred_list <- lapply(seq_len(n_eta), function(j) matrix(0, n_sim, unit$n_total))
  dtheta_list <- if (length(th_nms))
    stats::setNames(lapply(th_nms, function(nm) matrix(0, n_sim, unit$n_total)), th_nms)
  else NULL
  for (blk in unit$blocks) {
    # Forward covariates to per-block sensitivity solve.
    bs  <- list(ev_full = blk$ev_full, times = blk$times,
                cov = unit[["cov"]], cov_rows = unit[["cov_rows"]])
    res <- .admSimulateSens(sensModel, struct, sigma_names, eta_mat, bs, cores, ndp,
                            sigma_var, sigdig)
    if (is.null(res)) return(NULL)
    cp_mat[, blk$rows] <- res$cp_mat
    for (j in seq_len(n_eta)) dpred_list[[j]][, blk$rows] <- res$dpred_list[[j]]
    # Disable theta path for whole unit if any block lacks theta sensitivities.
    if (!is.null(dtheta_list)) {
      if (is.null(res$dtheta_list)) dtheta_list <- NULL
      else for (nm in th_nms) dtheta_list[[nm]][, blk$rows] <- res$dtheta_list[[nm]]
    }
  }
  list(cp_mat = cp_mat, dpred_list = dpred_list, dtheta_list = dtheta_list)
}

# Extract d(pred)/d(theta) for unpaired thetas from a sens solve (NULL on failure/absence).
.admThetaSens <- function(sensModel, out, keep, n_row, n_t) {
  tsc <- sensModel$theta_sens_cols
  if (is.null(tsc) || length(tsc) == 0L) return(NULL)
  if (!all(tsc %in% names(out))) return(NULL)
  stats::setNames(
    lapply(tsc, function(col)
      matrix(out[[col]][keep], nrow = n_row, ncol = n_t, byrow = TRUE)),
    names(tsc))
}

# Single pass on sensitivity model returning predictions and sensitivities.
.admSimulateSens <- function(sensModel, struct_theta, sigma_names,
                             eta_mat, study, cores,
                             ndp = .Machine$integer.max, sigma_var = NULL,
                             sigdig = NULL) {
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
  # fixed thetas: constants the estimated-parameter loops above never write
  for (nm in names(sensModel$fixed_theta))
    inner_df[[nm]] <- rep(unname(sensModel$fixed_theta[[nm]]), nrow(inner_df))

  # Supply current lambda estimate for TBS sensitivity solves.
  .tb <- sensModel$pred_tbs
  .lam <- if (is.null(.tb)) NA_real_ else .tb$lam
  if (!is.null(.tb) && !is.na(.tb$lam_name %||% NA_character_) &&
      !is.null(sigma_var) && .tb$lam_name %in% names(sigma_var)) {
    .lam   <- unname(sigma_var[[.tb$lam_name]])
    .mapped <- rmap[.tb$lam_name]
    if (!is.na(.mapped) && .mapped %in% names(inner_df)) inner_df[[.mapped]][] <- .lam
  }

  # Model covariates for this study (only names in study$cov -- see .admCovCols).
  inner_df <- .admCovCols(inner_df, sensModel$mod$params, study[["cov"]],
                          study[["cov_rows"]])

  # Forward solve_args (e.g. forced dop853 for DDE sensitivity models).
  out <- tryCatch(
    suppressWarnings(
      do.call(rxode2::rxSolve,
              c(list(sensModel$mod, params = inner_df,
                     events = study$ev_full, cores = cores,
                     nDisplayProgress = ndp,
                     sigdig = sigdig), sensModel$solve_args))),
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

  # Back-transform TBS/lnorm predictions and chain sensitivities by g'(z).
  if (!is.null(.tb)) {
    .gp        <- .admTBSid(cp_mat, .lam, .tb$yj, .tb$lo, .tb$hi)
    cp_mat     <- .admTBSi(cp_mat, .lam, .tb$yj, .tb$lo, .tb$hi)
    dpred_list <- lapply(dpred_list, function(D) D * .gp)
    if (!is.null(dtheta_list)) dtheta_list <- lapply(dtheta_list, function(D) D * .gp)
  }

  list(cp_mat = cp_mat, dpred_list = dpred_list, dtheta_list = dtheta_list)
}
