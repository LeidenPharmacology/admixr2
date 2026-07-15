# Emit a warning at most once per R session, keyed by `key`.
# Used to avoid repeated identical warnings when the same model is parsed
# once per study in multi-study fits.
.adm_warn_once <- function(key, msg) {
  if (is.null(.adm_warn_env[[key]])) {
    .adm_warn_env[[key]] <- TRUE
    warning(msg, call. = FALSE)
  }
}

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

# How many times each name appears in the model expressions.
#
# No model text (a hand-built mock ui in the Tier-1 tests) -> 1L, i.e. the
# shared-eta guard below is simply not applied and muRefDataFrame is trusted as
# before. A real rxUi always carries lstExpr, so this only affects mocks. (This
# is the opposite of nlmixr2est's conservative 2L default: there an unknown case
# costs one extra direction, here it would silently change struct_has_eta for
# every mock-based unit test.)
.admNameOccurrence <- function(ui, nms) {
  if (length(nms) == 0L) return(setNames(integer(0), character(0)))
  lst <- tryCatch(ui$lstExpr, error = function(e) NULL)
  if (is.null(lst)) return(setNames(rep(1L, length(nms)), nms))
  # unique = FALSE: count repeated uses within one line too
  syms <- unlist(lapply(lst, function(e) all.vars(e, unique = FALSE)))
  setNames(vapply(nms, function(nm) sum(syms == nm), integer(1)), nms)
}

# The theta <-> eta mu-reference map, with SHARED etas removed.
#
# A mu-referenced theta may reuse its eta's sensitivity column, because theta and
# eta enter the parameter identically: d(pred)/d(theta) == d(pred)/d(eta). That
# identity FAILS if the eta also appears in another parameter (`eta.cl` in both
# `cl` and `v`): d(pred)/d(eta.cl) then collects a path through `v` that
# d(pred)/d(tcl) does not have. Such a theta must not reuse the eta column -- it
# is treated as unpaired and gets its own sensitivity direction (a dummy eta in
# the augmented sens model), which is always exact.
#
# In practice rxode2 already declines to mu-reference a shared eta (it drops the
# rows from muRefDataFrame), so this is belt-and-braces -- it mirrors the
# equivalent guard in nlmixr2est (.foceiEtaOccurrence > 1 in foceiCovAnalytic.R).
# A false positive costs one extra sensitivity direction and stays exact.
# NULL means "no mu-reference information at all" (no ui, or no muRefDataFrame) --
# the callers then keep their historical fallbacks. A ZERO-ROW frame is different:
# it means the information exists and says nothing is paired (a non-mu-referenced
# model, or every pair dropped by the shared-eta guard). Conflating the two makes
# struct_eta_idx fall back to identity pairing (eta j <-> struct j), which would
# add the eta-path gradient on top of the theta column -- double counting.
.admMuRefPairs <- function(ui) {
  if (is.null(ui)) return(NULL)
  mrd <- tryCatch(ui$muRefDataFrame, error = function(e) NULL)
  if (is.null(mrd) || !is.data.frame(mrd) || !("theta" %in% names(mrd))) return(NULL)
  if (nrow(mrd) == 0L) return(mrd)
  eta_col <- if ("eta" %in% names(mrd)) "eta" else names(mrd)[2]
  occ  <- .admNameOccurrence(ui, unique(as.character(mrd[[eta_col]])))
  keep <- occ[as.character(mrd[[eta_col]])] <= 1L
  mrd[keep, , drop = FALSE]
}

# Extract parameter structure from ui$iniDf.
# Builds: struct/sigma/omega rows, Cholesky index vectors, per-theta transform metadata.
.admParseIniDf <- function(iniDf, ui = NULL) {
  theta_rows  <- iniDf[is.na(iniDf$neta1), , drop = FALSE]
  is_err      <- !is.na(theta_rows$err)
  struct_rows <- theta_rows[!is_err & !theta_rows$fix, , drop = FALSE]
  sigma_rows  <- theta_rows[ is_err & !theta_rows$fix, , drop = FALSE]

  # A FIXED omega is not supported. It is dropped from eta_rows here, so n_eta and
  # omega_init no longer cover every eta index -- omega_init[neta1, neta1] then
  # runs off the end of the matrix and the parse dies with a bare "subscript out of
  # bounds". Worse, if it did not, the eta's variance would be silently EXCLUDED
  # from the model rather than held at its fixed value, so the fit would quietly
  # ignore that source of between-subject variability. Fail with something the user
  # can act on instead. (Proper support means: keep the eta, hold its variance at
  # the fixed value, and exclude only its Cholesky entries from the optimizer.)
  .eta_all <- iniDf[!is.na(iniDf$neta1), , drop = FALSE]
  if (any(.eta_all$fix)) {
    .fx <- unique(.eta_all$name[.eta_all$fix])
    stop("admixr2 does not support FIXED omega entries (", paste(.fx, collapse = ", "),
         "). Remove fix() from the random effect, or drop the eta from the model ",
         "and fold its parameter into a structural theta.", call. = FALSE)
  }

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

  .err_vals   <- sigma_rows$err
  sigma_names <- sigma_rows$name

  # `pow(b, c)` emits TWO iniDf rows: the coefficient (err "pow") and the
  # EXPONENT (err "pow2"). The exponent is not a variance -- it must not be
  # squared, floored at zero, or reported as an SD -- so it carries its own
  # optimizer role and an identity transform. See .admSigmaRole()/.admSigmaNat().
  sigma_role <- setNames(
    ifelse(.err_vals %in% .ADM_ERR_POW_EXP, "pow_exp", "var"), sigma_names)

  sigma_is_prop  <- setNames(.err_vals %in% .ADM_ERR_PROP, sigma_names)
  sigma_is_lnorm <- setNames(.err_vals %in% .ADM_ERR_LNORM, sigma_names)

  # Refuse anything we cannot represent. This runs BEFORE .admBuildResidSpecs(),
  # so it -- not the predDf gates below it -- is what a user with a boxCox or
  # Student-t endpoint actually sees; it therefore has to give the same quality of
  # explanation. .ADM_ERR_WHY supplies the per-type reason.
  .unsupported <- unique(.err_vals[!is.na(.err_vals) & !.err_vals %in% .ADM_ERR_KNOWN])
  if (length(.unsupported) > 0L) {
    .bad <- sort(.unsupported)[1L]
    .ep  <- if ("condition" %in% names(sigma_rows)) {
      .w <- which(.err_vals == .bad)[1L]
      as.character(sigma_rows$condition[.w])
    } else NA_character_
    .admStopErrModel(
      .ep,
      paste0(.bad, "()",
             if (length(.unsupported) > 1L)
               paste0(" (also: ", paste(sort(.unsupported)[-1L], collapse = ", "), ")")
             else ""),
      .ADM_ERR_WHY[[.bad]] %||%
        paste0("'", .bad, "' has no aggregate mean/variance admixr2 can score"))
  }

  # Output variable each residual-error parameter belongs to. In a real nlmixr2
  # iniDf the error rows carry a `condition` column naming the endpoint (e.g.
  # "cp"); with several observed outputs this maps each sigma to its output so a
  # given output's likelihood uses only its own residual error. NA when the
  # column is absent (Tier-1 mock iniDf, or single-output models) -- callers then
  # treat every sigma as belonging to the single output. See .admSigmaSel().
  sigma_output <- if ("condition" %in% names(sigma_rows))
    as.character(sigma_rows$condition) else rep(NA_character_, length(sigma_names))

  # Residual parameters enter the optimizer on their role's scale: variances as
  # log(sigma^2), a pow() exponent as itself.
  .is_var    <- sigma_role == "var"
  sigma_init <- setNames(ifelse(.is_var, 2 * log(sigma_rows$est), sigma_rows$est),
                         sigma_names)

  list(struct_names    = struct_rows$name,
       struct_init     = setNames(struct_rows$est,   struct_rows$name),
       struct_lower    = setNames(struct_rows$lower, struct_rows$name),
       struct_upper    = setNames(struct_rows$upper, struct_rows$name),
       sigma_names     = sigma_names,
       sigma_init      = sigma_init,
       sigma_lower     = setNames(sigma_rows$lower, sigma_rows$name),
       sigma_upper     = setNames(sigma_rows$upper, sigma_rows$name),
       sigma_role      = sigma_role,
       resid           = .admBuildResidSpecs(ui, sigma_rows, sigma_names),
       sigma_is_prop   = sigma_is_prop,
       sigma_is_lnorm  = sigma_is_lnorm,
       sigma_output    = sigma_output,
       # rxSolve progress-bar threshold; the estimators override this from the
       # control, but default it here so pinfo built directly (e.g. in tests that
       # call the NLL/gradient internals) always carries a valid integer.
       nDisplayProgress = .Machine$integer.max,
       eta_names       = eta_names, n_eta = n_eta,
       eta_col_names   = paste0("eta.", gsub("^eta\\.", "", eta_names)),
       omega_par       = omega_par,
       omega_par_names = omega_par_names,
       chol_i          = chol_i, chol_j = chol_j, chol_diag = chol_diag,
       iniDf              = iniDf,
       eta_rows_df        = eta_rows,
       has_kappa          = has_kappa,
       struct_has_eta     = setNames(
         struct_rows$name %in% {
           .mrd <- .admMuRefPairs(ui)
           if (!is.null(.mrd)) .mrd$theta else eta_names
         },
         struct_rows$name),
       struct_eta_idx     = {
         mrd <- .admMuRefPairs(ui)
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
  # Variance params are bounded on the log-variance scale; a pow() exponent is
  # unconstrained and its bounds pass through untransformed.
  .is_var <- .admSigmaRole(pinfo) == "var"
  sig_lb <- ifelse(!.is_var, pinfo$sigma_lower,
                   ifelse(is.finite(pinfo$sigma_lower) & pinfo$sigma_lower > 0,
                          2 * log(pinfo$sigma_lower), -Inf))
  sig_ub <- ifelse(!.is_var, pinfo$sigma_upper,
                   ifelse(is.finite(pinfo$sigma_upper),
                          2 * log(pinfo$sigma_upper),  Inf))
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
  # Natural scale: variances back-transform as exp(p); a pow() exponent is
  # already on its natural scale.
  sigma_var <- .admSigmaNat(p[n_s + seq_len(n_e)], pinfo)

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
