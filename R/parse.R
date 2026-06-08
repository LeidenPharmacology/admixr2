# Back-transform one struct-theta from optimizer scale to natural scale.
.admBackTransform <- function(p, tr) {
  if (is.null(tr)) return(exp(p))
  switch(tr$curEval,
    exp = , log = exp(p),
    expit = , logit = rxode2::expit(p, tr$low, tr$hi),
    probitInv = , probit = tr$low + (tr$hi - tr$low) * pnorm(p),
    p)
}

# log(back(p)) -- used in IRMC gradient chain rule for paired struct thetas.
.admLogBackTransform <- function(p, tr) {
  if (is.null(tr)) return(p)
  switch(tr$curEval,
    exp = , log = p,
    expit = , logit = log(rxode2::expit(p, tr$low, tr$hi)),
    probitInv = , probit = log(tr$low + (tr$hi - tr$low) * pnorm(p)),
    log(p))
}

# Extract parameter structure from ui$iniDf.
# Builds: struct/sigma/omega rows, Cholesky index vectors, per-theta transform metadata.
.admParseIniDf <- function(iniDf, ui = NULL) {
  theta_rows  <- iniDf[is.na(iniDf$neta1), , drop = FALSE]
  is_err      <- !is.na(theta_rows$err)
  struct_rows <- theta_rows[!is_err & !theta_rows$fix, , drop = FALSE]
  sigma_rows  <- theta_rows[ is_err & !theta_rows$fix, , drop = FALSE]

  eta_rows  <- iniDf[!is.na(iniDf$neta1) & !iniDf$fix, , drop = FALSE]
  diag_rows <- eta_rows[eta_rows$neta1 == eta_rows$neta2, , drop = FALSE]
  eta_names <- diag_rows$name
  n_eta     <- length(eta_names)

  omega_init <- matrix(0, n_eta, n_eta, dimnames = list(eta_names, eta_names))
  omega_par <- numeric(0); omega_par_names <- character(0)
  chol_i <- integer(0); chol_j <- integer(0); chol_diag <- logical(0)

  if (n_eta > 0) {
    for (r in seq_len(nrow(eta_rows))) {
      i <- eta_rows$neta1[r]; j <- eta_rows$neta2[r]
      omega_init[i, j] <- eta_rows$est[r]
      omega_init[j, i] <- eta_rows$est[r]
    }
    L_init <- t(chol(omega_init))
    for (r in seq_len(nrow(eta_rows))) {
      i <- eta_rows$neta1[r]; j <- eta_rows$neta2[r]
      if (i == j) {
        # Store log(Omega_ii) = 2*log(L_ii), NOT log(L_ii).
        # With log(L_ii), a unit optimizer step changes Omega_ii by 2*Omega_ii (chain rule x2),
        # making IS weights 2x more sensitive per LBFGS step -> fast IS degeneracy in IRMC.
        # log(Omega_ii) gives unit step -> Omega_ii change of Omega_ii, matching adm reference behavior.
        omega_par       <- c(omega_par, 2 * log(L_init[i, i]))
        omega_par_names <- c(omega_par_names, paste0("logchol_", eta_names[i]))
        chol_i <- c(chol_i, i); chol_j <- c(chol_j, j); chol_diag <- c(chol_diag, TRUE)
      } else if (i > j) {
        omega_par       <- c(omega_par, L_init[i, j])
        omega_par_names <- c(omega_par_names,
                             paste0("chol_", eta_names[i], "_", eta_names[j]))
        chol_i <- c(chol_i, i); chol_j <- c(chol_j, j); chol_diag <- c(chol_diag, FALSE)
      }
    }
  }
  names(omega_par) <- omega_par_names

  # has_kappa: TRUE when some struct thetas are not mu-referenced to an eta.
  has_kappa <- if (!is.null(ui) && !is.null(ui$muRefDataFrame) &&
                   "theta" %in% names(ui$muRefDataFrame)) {
    any(!(struct_rows$name %in% ui$muRefDataFrame$theta))
  } else {
    n_eta < nrow(struct_rows)
  }

  # Per-struct-theta transform metadata from ui$muRefCurEval.
  # Stores only curEval/low/hi -- no closures.  Use .admBackTransform() /
  # .admLogBackTransform() to evaluate transforms at runtime.
  .ce <- if (!is.null(ui)) tryCatch(ui$muRefCurEval, error = function(e) NULL) else NULL

  struct_transforms <- setNames(lapply(struct_rows$name, function(nm) {
    .w <- if (!is.null(.ce)) which(.ce$parameter == nm) else integer(0)
    if (length(.w) != 1L)
      return(list(curEval = "exp", low = NA_real_, hi = NA_real_))
    list(curEval = .ce$curEval[.w],
         low     = if ("low" %in% names(.ce)) .ce$low[.w] else NA_real_,
         hi      = if ("hi"  %in% names(.ce)) .ce$hi[.w]  else NA_real_)
  }), struct_rows$name)

  .err_vals  <- sigma_rows$err
  .prop_approx <- unique(.err_vals[!is.na(.err_vals) & .err_vals %in% c("propT", "propF")])
  if (length(.prop_approx) > 0L)
    warning("Residual error type(s) ", paste(.prop_approx, collapse = ", "),
            " modelled as proportional (prop). Transform-aware scaling ignored.",
            call. = FALSE)
  .add_approx <- unique(.err_vals[!is.na(.err_vals) & .err_vals %in% c("norm", "dnorm")])
  if (length(.add_approx) > 0L)
    warning("Residual error type(s) ", paste(.add_approx, collapse = ", "),
            " modelled as additive (add). Likelihood-path distinction ignored.",
            call. = FALSE)
  .lnorm_approx <- unique(.err_vals[!is.na(.err_vals) & .err_vals %in% c("dlnorm", "logn", "dlogn")])
  if (length(.lnorm_approx) > 0L)
    warning("Residual error type(s) ", paste(.lnorm_approx, collapse = ", "),
            " modelled as lognormal (lnorm). Likelihood-path distinction ignored.",
            call. = FALSE)
  .supported <- c("add", "norm", "dnorm", "prop", "propT", "propF",
                  "lnorm", "dlnorm", "logn", "dlogn")
  .unsupported <- unique(.err_vals[!is.na(.err_vals) & !.err_vals %in% .supported])
  if (length(.unsupported) > 0L)
    warning("Unsupported residual error type(s) detected: ",
            paste(.unsupported, collapse = ", "),
            ". Treated as additive. Supported: add/norm, prop, lnorm.",
            call. = FALSE)
  sigma_names <- sigma_rows$name
  sigma_is_prop  <- setNames(.err_vals %in% c("prop", "propT", "propF"), sigma_names)
  sigma_is_lnorm <- setNames(.err_vals %in% c("lnorm", "dlnorm", "logn", "dlogn"), sigma_names)

  list(struct_names    = struct_rows$name,
       struct_init     = setNames(struct_rows$est,   struct_rows$name),
       struct_lower    = setNames(struct_rows$lower, struct_rows$name),
       struct_upper    = setNames(struct_rows$upper, struct_rows$name),
       sigma_names     = sigma_names,
       sigma_init      = setNames(2 * log(sigma_rows$est), sigma_rows$name),
       sigma_lower     = setNames(sigma_rows$lower, sigma_rows$name),
       sigma_upper     = setNames(sigma_rows$upper, sigma_rows$name),
       sigma_is_prop   = sigma_is_prop,
       sigma_is_lnorm  = sigma_is_lnorm,
       eta_names       = eta_names, n_eta = n_eta,
       eta_col_names   = paste0("eta.", gsub("^eta\\.", "", eta_names)),
       omega_par       = omega_par,
       omega_par_names = omega_par_names,
       chol_i          = chol_i, chol_j = chol_j, chol_diag = chol_diag,
       iniDf              = iniDf,
       eta_rows_df        = eta_rows,
       has_kappa          = has_kappa,
       struct_has_eta     = setNames(
         struct_rows$name %in% if (!is.null(ui) && !is.null(ui$muRefDataFrame) &&
                                   "theta" %in% names(ui$muRefDataFrame))
           ui$muRefDataFrame$theta else eta_names,
         struct_rows$name),
       struct_eta_idx     = {
         mrd <- if (!is.null(ui) && !is.null(ui$muRefDataFrame) &&
                    "theta" %in% names(ui$muRefDataFrame))
           ui$muRefDataFrame else NULL
         if (!is.null(mrd)) {
           eta_col <- if ("eta" %in% names(mrd)) "eta" else names(mrd)[2]
           vapply(eta_names, function(en) {
             theta_nm <- mrd$theta[mrd[[eta_col]] == en]
             if (length(theta_nm) == 0L) {
               en_bare  <- gsub("^eta\\.", "", en)
               theta_nm <- mrd$theta[gsub("^eta\\.", "", mrd[[eta_col]]) == en_bare]
             }
             if (length(theta_nm) == 0L) return(NA_integer_)
             idx <- match(theta_nm[1L], struct_rows$name)
             if (is.na(idx)) NA_integer_ else idx
           }, integer(1))
         } else {
           seq_len(n_eta)
         }
       },
       struct_transforms  = struct_transforms)
}

# Convert pinfo to initial optimizer vector with bounds.
.admBuildOptVec <- function(pinfo) {
  p  <- c(pinfo$struct_init, pinfo$sigma_init, pinfo$omega_par)
  nm <- c(pinfo$struct_names, pinfo$sigma_names, pinfo$omega_par_names)
  names(p) <- nm
  sig_lb <- ifelse(is.finite(pinfo$sigma_lower) & pinfo$sigma_lower > 0,
                   2 * log(pinfo$sigma_lower), -Inf)
  sig_ub <- ifelse(is.finite(pinfo$sigma_upper),
                   2 * log(pinfo$sigma_upper),  Inf)
  lb <- c(pinfo$struct_lower, sig_lb, rep(-Inf, length(pinfo$omega_par)))
  ub <- c(pinfo$struct_upper, sig_ub, rep( Inf, length(pinfo$omega_par)))
  list(p0 = p, lower = lb, upper = ub, names = nm,
       scale_c = .admComputeScaleC(pinfo))
}

# Per-parameter optimizer pre-conditioning scales.
# Optimizer sees p_scaled = p_real / scale_c; gradients rescaled by scale_c.
# - struct thetas (exp): 1 (log-scale already normalized).
# - struct thetas (expit/probitInv): derivative-based magnitude at init point.
# - sigma: 1 (log(sigma^2) encoding self-normalizing).
# - omega diagonal: 1 (log(Omega_ii) encoding self-normalizing).
# - omega off-diagonal: pmax(|L_ij_init|, 0.1) (raw L values need magnitude scaling).
.admComputeScaleC <- function(pinfo) {
  struct_sc <- vapply(pinfo$struct_names, function(nm) {
    tr <- pinfo$struct_transforms[[nm]]
    if (is.null(tr)) return(1.0)
    p  <- pinfo$struct_init[[nm]]
    switch(tr$curEval,
      exp = ,
      log = 1.0,
      expit = ,
      logit = {
        a <- tr$low; b <- tr$hi
        pmax(exp(p) * (1 + exp(-p))^2 * (a + (b - a) / (1 + exp(-p))) / (b - a), 0.01)
      },
      probitInv = ,
      probit = {
        a <- tr$low; b <- tr$hi
        pmax(sqrt(2) * exp(0.5 * p^2) * sqrt(pi) *
               (a + 0.5 * (b - a) * (1 + rxode2::erf(p / sqrt(2)))) / (b - a), 0.01)
      },
      1.0)
  }, double(1))

  sigma_sc <- rep(1.0, length(pinfo$sigma_names))

  n_o <- length(pinfo$omega_par)
  omega_sc <- rep(1.0, n_o)
  if (n_o > 0L && any(!pinfo$chol_diag)) {
    off <- !pinfo$chol_diag
    omega_sc[off] <- pmax(abs(pinfo$omega_par[off]), 0.1)
  }

  setNames(c(struct_sc, sigma_sc, omega_sc),
           c(pinfo$struct_names, pinfo$sigma_names, pinfo$omega_par_names))
}

# Unpack optimizer vector p into named parameter lists.
.admUnpack <- function(p, pinfo) {
  n_s   <- length(pinfo$struct_names)
  n_e   <- length(pinfo$sigma_names)
  n_o   <- length(pinfo$omega_par)
  n_eta <- pinfo$n_eta

  struct    <- setNames(p[seq_len(n_s)], pinfo$struct_names)
  sigma_var <- setNames(exp(p[n_s + seq_len(n_e)]), pinfo$sigma_names)

  if (n_eta > 0) {
    om_p <- p[n_s + n_e + seq_len(n_o)]
    L    <- matrix(0, n_eta, n_eta)
    d <- pinfo$chol_diag; nd <- !d
    # Diagonal: p = log(Omega_ii), so L_ii = sqrt(Omega_ii) = exp(p/2).
    L[cbind(pinfo$chol_i[d],  pinfo$chol_j[d])]  <- exp(om_p[d] / 2)
    # Off-diagonal: p = L_ij directly (no transform).
    L[cbind(pinfo$chol_i[nd], pinfo$chol_j[nd])] <- om_p[nd]
    omega <- L %*% t(L)
  } else {
    L <- NULL; omega <- matrix(0, 0, 0)
  }

  list(struct = struct, sigma_var = sigma_var, omega = omega, L = L)
}
