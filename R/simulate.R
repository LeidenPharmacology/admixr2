# Single rxSolve pass for one study given pre-computed eta_mat (n_sim x n_eta).
# Returns n_sim x n_times matrix of predicted concentrations.
# params_mat is a named numeric matrix (from .admMakeParamsList); converted to
# data.frame only at the rxSolve call to avoid repeated list COW copies.
.admSimulate <- function(rxMod, struct_theta, sigma_names, eta_mat, study,
                         output_var, params_mat, cores) {
  eta_cols <- colnames(eta_mat)
  for (nm in names(struct_theta)) params_mat[, nm] <- struct_theta[nm]
  if (length(eta_cols) > 0L)      params_mat[, eta_cols] <- eta_mat
  for (nm in sigma_names)         params_mat[, nm] <- 0
  # linCmt simulationModel requires rxerr.rxLinCmt; add any rxMod params not
  # already in params_mat as zeros so rxSolve receives a complete params frame.
  extra <- setdiff(rxMod$params, colnames(params_mat))
  if (length(extra) > 0L)
    params_mat <- cbind(params_mat,
                        matrix(0, nrow(params_mat), length(extra),
                               dimnames = list(NULL, extra)))
  out  <- rxode2::rxSolve(rxMod, params = as.data.frame(params_mat),
                          events = study$ev_full, cores = cores,
                          nDisplayProgress = .Machine$integer.max)
  keep <- out[["time"]] %in% study$times
  # linCmt simulationModel outputs "ipredSim" rather than "rx_pred_"
  vals <- out[[output_var]]
  if (is.null(vals)) vals <- out[["ipredSim"]]
  matrix(vals[keep],
         nrow = nrow(eta_mat), ncol = length(study$times), byrow = TRUE)
}

# Single pass on the sensitivity model returning predictions + d(pred)/d(eta_j).
# Returns list(cp_mat, dpred_list) or NULL on failure (caller falls back to FD).
.admSimulateSens <- function(sensModel, struct_theta, sigma_names,
                             eta_mat, study, cores) {
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
                      nDisplayProgress = .Machine$integer.max)),
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
