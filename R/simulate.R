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
  # linCmt simulationModel outputs "ipredSim" rather than "rx_pred_"
  vals <- out[[output_var]]
  if (is.null(vals)) vals <- out[["ipredSim"]]
  matrix(vals[keep],
         nrow = nrow(eta_mat), ncol = length(study$times), byrow = TRUE)
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
  vals <- out[[output_var]]
  if (is.null(vals)) vals <- out[["ipredSim"]]
  matrix(vals[keep], nrow = nrow(struct_mat), ncol = length(study$times),
         byrow = TRUE)
}

.admSimulateSensRows <- function(sensModel, struct_mat, sigma_names, eta_mat, study,
                                 cores, ndp = .Machine$integer.max) {
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

  out <- tryCatch(
    suppressWarnings(
      rxode2::rxSolve(sensModel$mod, params = inner_df,
                      events = study$ev_full, cores = cores,
                      nDisplayProgress = ndp)),
    error = function(e) NULL)
  if (is.null(out)) return(NULL)
  if (!all(sensModel$sens_cols %in% names(out))) return(NULL)

  keep  <- out[["time"]] %in% study$times
  n_t   <- length(study$times)
  n_eta <- ncol(eta_mat)

  cp_mat     <- matrix(out[["rx_pred_"]][keep], nrow = n_row, ncol = n_t, byrow = TRUE)
  dpred_list <- lapply(seq_len(n_eta), function(j)
    matrix(out[[sensModel$sens_cols[j]]][keep], nrow = n_row, ncol = n_t, byrow = TRUE))

  list(cp_mat = cp_mat, dpred_list = dpred_list)
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
    vals <- out[[blk$output]]
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
                                  cores, ndp = .Machine$integer.max) {
  n_sim <- nrow(eta_mat); n_eta <- ncol(eta_mat)
  cp_mat     <- matrix(0, n_sim, unit$n_total)
  dpred_list <- lapply(seq_len(n_eta), function(j) matrix(0, n_sim, unit$n_total))
  for (blk in unit$blocks) {
    bs  <- list(ev_full = blk$ev_full, times = blk$times)
    res <- .admSimulateSens(sensModel, struct, sigma_names, eta_mat, bs, cores, ndp)
    if (is.null(res)) return(NULL)
    cp_mat[, blk$rows] <- res$cp_mat
    for (j in seq_len(n_eta)) dpred_list[[j]][, blk$rows] <- res$dpred_list[[j]]
  }
  list(cp_mat = cp_mat, dpred_list = dpred_list)
}

# Single pass on the sensitivity model returning predictions + d(pred)/d(eta_j).
# Returns list(cp_mat, dpred_list) or NULL on failure (caller falls back to FD).
.admSimulateSens <- function(sensModel, struct_theta, sigma_names,
                             eta_mat, study, cores,
                             ndp = .Machine$integer.max) {
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

  out <- tryCatch(
    suppressWarnings(
      rxode2::rxSolve(sensModel$mod, params = inner_df,
                      events = study$ev_full, cores = cores,
                      nDisplayProgress = ndp)),
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

  list(cp_mat = cp_mat, dpred_list = dpred_list)
}
