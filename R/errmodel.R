# Residual error models -------------------------------------------------------
#
# Every estimator needs the same three things from the residual error model:
# the mean it induces, the variance it adds to diag(V), and the derivative of
# both with respect to the residual parameters. Those used to be spelled out as
# an `if (prop) ... else if (lnorm) ... else ...` chain repeated at ~8 sites in
# R and 3 more in C++, which is why the supported set never grew past
# add/prop/lnorm. This file is the single place that knows the error models;
# everything else consumes the row arrays it builds.
#
# The general form, for a prediction f and an endpoint's residual parameters:
#
#   combined2   var = a^2 + b^2 * f^(2c)          (variance-additive)
#   combined1   var = (a + b * f^c)^2             (SD-additive)
#   lnorm       mu  = f * exp(s/2)
#               var = mu^2 * (exp(s) - 1)         (moment-matched lognormal)
#
# with a the additive SD, b the proportional/power SD and c the power exponent
# (c = 1 recovers prop, b = 0 recovers add). combined2 with c = 1 is exactly the
# old independent per-sigma addition, so add/prop/lnorm fits are unchanged.
#
# Encoding of `form`: 0 = combined2, 1 = combined1, 2 = lnorm. The C++ kernels
# take the same four row arrays (form, a2, b2, cc) and so are error-model
# agnostic -- they evaluate the residual against whatever mu they compute.

.ADM_RESID_COMBINED2 <- 0L
.ADM_RESID_COMBINED1 <- 1L
.ADM_RESID_LNORM     <- 2L

# iniDf$err vocabularies (rxode2 spells the same model several ways).
.ADM_ERR_ADD     <- c("add", "norm", "dnorm")
.ADM_ERR_PROP    <- c("prop", "propT", "propF")
.ADM_ERR_POW     <- c("pow", "powT", "powF")
.ADM_ERR_POW_EXP <- c("pow2", "powT2", "powF2")   # the EXPONENT row of pow(b, c)
.ADM_ERR_LNORM   <- c("lnorm", "dlnorm", "logn", "dlogn")
.ADM_ERR_KNOWN   <- c(.ADM_ERR_ADD, .ADM_ERR_PROP, .ADM_ERR_POW,
                      .ADM_ERR_POW_EXP, .ADM_ERR_LNORM)

# Build one residual spec per endpoint from ui$predDf.
#
# predDf is authoritative and iniDf$err alone is not: `errType` tells us which
# terms are present, `addProp` picks combined1 vs combined2, `transform` and
# `distribution` say whether the endpoint is representable at all, and `errTypeF`
# says which scale the proportional/power term multiplies.
#
# Returns NULL when there is no predDf (Tier-1 mock iniDf / no ui); callers then
# fall back to the legacy sigma_is_prop/sigma_is_lnorm flags.
.admBuildResidSpecs <- function(ui, sigma_rows, sigma_names) {
  predDf <- if (!is.null(ui)) tryCatch(as.data.frame(ui$predDf), error = function(e) NULL) else NULL
  if (is.null(predDf) || nrow(predDf) == 0L) return(NULL)
  if (length(sigma_names) == 0L) return(NULL)

  err  <- sigma_rows$err
  cond <- if ("condition" %in% names(sigma_rows))
    as.character(sigma_rows$condition) else rep(NA_character_, length(sigma_names))

  # nlmixr2's addProp default resolves to combined2.
  add_prop_default <- getOption("rxode2.addProp", "combined2")

  specs <- list(); keys <- character(0)
  for (i in seq_len(nrow(predDf))) {
    ep   <- as.character(predDf$cond[i])
    tr   <- as.character(predDf$transform[i] %||% "untransformed")
    etf  <- as.character(predDf$errTypeF[i]  %||% "untransformed")
    dist <- as.character(predDf$distribution[i] %||% "norm")

    if (!dist %in% c("norm", "dnorm", "t", "cauchy") || dist %in% c("t", "cauchy"))
      stop("Endpoint '", ep, "': distribution '", dist,
           "' is not supported. admixr2 models aggregate data as multivariate ",
           "normal; only normally-distributed residuals can be represented.",
           call. = FALSE)
    if (!tr %in% c("untransformed", "lnorm"))
      stop("Endpoint '", ep, "': residual transform '", tr, "' is not supported. ",
           "Supported transforms: untransformed, lnorm.", call. = FALSE)
    if (identical(etf, "f"))
      stop("Endpoint '", ep, "': propF()/powF() scale the residual by a ",
           "user-supplied model variable, which admixr2 cannot resolve from the ",
           "aggregate mean. Use prop()/pow() instead.", call. = FALSE)
    if (isTRUE(predDf$variance[i]))
      stop("Endpoint '", ep, "': variance-parameterised residual error ",
           "(variance = TRUE) is not supported; parameterise in SD units.",
           call. = FALSE)

    # Rows of iniDf belonging to this endpoint.
    idx    <- which(!is.na(cond) & cond == ep)
    if (length(idx) == 0L) idx <- which(is.na(cond))   # single-endpoint / no condition
    e      <- err[idx]
    k_add  <- idx[e %in% c(.ADM_ERR_ADD, .ADM_ERR_LNORM)]
    k_prop <- idx[e %in% c(.ADM_ERR_PROP, .ADM_ERR_POW)]
    k_pow  <- idx[e %in% .ADM_ERR_POW_EXP]

    # lnorm: the "add" parameter is the SD on the log scale; it cannot be
    # combined with a proportional or power term.
    if (identical(tr, "lnorm")) {
      if (length(k_prop) > 0L)
        stop("Endpoint '", ep, "': lnorm() combined with a proportional/power ",
             "term is not supported.", call. = FALSE)
      form <- .ADM_RESID_LNORM
    } else {
      ap <- as.character(predDf$addProp[i] %||% "default")
      if (identical(ap, "default")) ap <- add_prop_default
      form <- if (identical(ap, "combined1")) .ADM_RESID_COMBINED1 else .ADM_RESID_COMBINED2
    }

    key <- .admOutputColName(ep)
    specs[[key]] <- list(
      output = ep, form = form,
      k_add  = if (length(k_add)  > 0L) k_add[1L]  else NA_integer_,
      k_prop = if (length(k_prop) > 0L) k_prop[1L] else NA_integer_,
      k_pow  = if (length(k_pow)  > 0L) k_pow[1L]  else NA_integer_)
    keys <- c(keys, key)
  }
  if (length(specs) == 0L) return(NULL)
  specs
}

# Optimizer-scale role of each residual parameter.
#   "var"     log-variance: natural value is a VARIANCE, exp(p).  (a^2, b^2)
#   "pow_exp" a power exponent, not a variance: natural value is p itself.
# The distinction matters because a pow() exponent must not be squared, bounded
# below at zero, or reported as an SD. Absent (legacy/hand-built pinfo) => all "var".
.admSigmaRole <- function(pinfo) {
  r <- pinfo$sigma_role
  if (is.null(r)) rep("var", length(pinfo$sigma_names)) else r
}

# Optimizer vector -> natural-scale residual parameters.
.admSigmaNat <- function(p_sigma, pinfo) {
  role <- .admSigmaRole(pinfo)
  out  <- p_sigma
  iv   <- role == "var"
  out[iv] <- exp(p_sigma[iv])   # log-variance -> variance
  out[!iv] <- p_sigma[!iv]      # exponent: identity
  setNames(out, pinfo$sigma_names)
}

# Per-endpoint residual specs.
#
# `pinfo$resid` is built by .admParseIniDf from ui$predDf. When it is absent --
# a hand-built pinfo (unit tests) or a Tier-1 mock iniDf with no ui -- fall back
# to the legacy sigma_is_prop / sigma_is_lnorm flags, which describe exactly the
# combined2/c=1 subset. This keeps every existing caller working unchanged.
.admResidSpecs <- function(pinfo) {
  if (!is.null(pinfo$resid)) return(pinfo$resid)

  n <- length(pinfo$sigma_names)
  if (n == 0L) return(list())
  is_prop  <- pinfo$sigma_is_prop
  is_lnorm <- pinfo$sigma_is_lnorm
  so       <- pinfo$sigma_output
  if (is.null(so)) so <- rep(NA_character_, n)

  # One spec per distinct output; a legacy sigma set is at most {add, prop} plus
  # possibly a standalone lnorm, all with exponent 1.
  keys <- unique(so)
  specs <- lapply(keys, function(key) {
    sel <- if (is.na(key)) rep(TRUE, n) else !is.na(so) & so == key
    idx <- which(sel)
    k_lnorm <- idx[vapply(idx, function(k) isTRUE(is_lnorm[[k]]), logical(1))]
    k_prop  <- idx[vapply(idx, function(k) isTRUE(is_prop[[k]]),  logical(1))]
    k_add   <- setdiff(idx, c(k_lnorm, k_prop))
    if (length(k_lnorm) > 0L)
      return(list(output = key, form = .ADM_RESID_LNORM,
                  k_add = k_lnorm[1L], k_prop = NA_integer_, k_pow = NA_integer_))
    list(output = key, form = .ADM_RESID_COMBINED2,
         k_add  = if (length(k_add)  > 0L) k_add[1L]  else NA_integer_,
         k_prop = if (length(k_prop) > 0L) k_prop[1L] else NA_integer_,
         k_pow  = NA_integer_)
  })
  names(specs) <- as.character(keys)
  specs
}

# Look up the spec governing `output`. NA/NULL output (single-output model, or a
# pinfo with no sigma->output map) selects the sole spec.
.admResidSpecFor <- function(pinfo, output) {
  specs <- .admResidSpecs(pinfo)
  if (length(specs) == 0L) return(NULL)
  if (length(specs) == 1L) return(specs[[1L]])
  if (is.null(output) || is.na(output)) return(specs[[1L]])
  for (sp in specs) {
    if (is.na(sp$output)) return(sp)
    if (.admOutputColName(sp$output) == output) return(sp)
  }
  NULL
}

# Row arrays for a stacked mean vector.
#
# `row_output` gives each row's output column name (length n_t); pass a single
# value, or NULL, for a single-output vector. Returns the four parallel arrays
# the R and C++ residual evaluators both consume, plus `k_*` row->sigma maps used
# by the gradient.
.admResidRows <- function(pinfo, row_output, sigma_nat, n_t) {
  form <- integer(n_t)
  a2   <- numeric(n_t)   # additive VARIANCE  (a^2), or the lnorm log-variance
  b2   <- numeric(n_t)   # prop/pow VARIANCE  (b^2)
  cc   <- rep(1.0, n_t)  # power exponent
  k_add <- k_prop <- k_pow <- rep(NA_integer_, n_t)

  if (length(pinfo$sigma_names) == 0L)
    return(list(form = form, a2 = a2, b2 = b2, cc = cc,
                k_add = k_add, k_prop = k_prop, k_pow = k_pow))

  if (is.null(row_output)) row_output <- rep(NA_character_, n_t)
  if (length(row_output) == 1L) row_output <- rep(row_output, n_t)

  for (ov in unique(row_output)) {
    rows <- which(row_output == ov | (is.na(row_output) & is.na(ov)))
    sp   <- .admResidSpecFor(pinfo, ov)
    if (is.null(sp)) next
    form[rows] <- sp$form
    if (!is.na(sp$k_add))  { a2[rows] <- sigma_nat[[sp$k_add]];  k_add[rows]  <- sp$k_add }
    if (!is.na(sp$k_prop)) { b2[rows] <- sigma_nat[[sp$k_prop]]; k_prop[rows] <- sp$k_prop }
    if (!is.na(sp$k_pow))  { cc[rows] <- sigma_nat[[sp$k_pow]];  k_pow[rows]  <- sp$k_pow }
  }
  list(form = form, a2 = a2, b2 = b2, cc = cc,
       k_add = k_add, k_prop = k_prop, k_pow = k_pow)
}

# f^c, with an exact fast path at c == 1 so that prop/add models reproduce the
# old `mu^2` arithmetic bit-for-bit (and never hit pow()'s NaN for negative f).
.admPowF <- function(f, cc) {
  if (all(cc == 1)) return(f)
  ifelse(cc == 1, f, f^cc)
}

# Apply the residual to a structural mean/variance.
# Returns the residual-adjusted mean (lnorm shifts it) and diagonal variance.
.admResidApply <- function(mu_struct, dv, arr) {
  mu <- mu_struct
  ln <- arr$form == .ADM_RESID_LNORM
  nz <- !ln

  if (any(nz)) {
    fe <- .admPowF(mu_struct[nz], arr$cc[nz])
    c1 <- arr$form[nz] == .ADM_RESID_COMBINED1
    # combined2: accumulate add then prop, in that order, so add/prop models keep
    # the exact floating-point sequence the old per-sigma loop produced.
    add <- arr$a2[nz]; prp <- arr$b2[nz] * fe^2
    v   <- add + prp
    if (any(c1)) {
      s <- sqrt(arr$a2[nz][c1]) + sqrt(arr$b2[nz][c1]) * fe[c1]
      v[c1] <- s^2
    }
    dv[nz] <- dv[nz] + v
  }
  if (any(ln)) {
    sv <- arr$a2[ln]
    mu[ln] <- mu_struct[ln] * exp(sv / 2)
    dv[ln] <- dv[ln] + mu[ln]^2 * (exp(sv) - 1)
  }
  list(mu = mu, dv = dv)
}

# Derivatives of the residual w.r.t. the residual parameters (optimizer scale)
# and w.r.t. the structural prediction f.
#
# Returns:
#   dmu   n_t x n_sigma   d(mu)/d(p_k)         (nonzero only for lnorm)
#   dvar  n_t x n_sigma   d(var)/d(p_k)
#   dv_df n_t             d(var)/d(f)          (the V-path of the struct-theta chain)
#
# For a "var"-role parameter the optimizer holds p = log(sv), so d(sv)/dp = sv;
# for a "pow_exp" parameter p is the exponent itself and d(c)/dp = 1.
.admResidDeriv <- function(mu_struct, arr, pinfo) {
  n_t   <- length(mu_struct)
  n_sig <- length(pinfo$sigma_names)
  dmu   <- matrix(0, n_t, n_sig)
  dvar  <- matrix(0, n_t, n_sig)
  dv_df <- numeric(n_t)
  if (n_sig == 0L) return(list(dmu = dmu, dvar = dvar, dv_df = dv_df))

  for (t in seq_len(n_t)) {
    f  <- mu_struct[t]
    ka <- arr$k_add[t]; kb <- arr$k_prop[t]; kc <- arr$k_pow[t]

    if (arr$form[t] == .ADM_RESID_LNORM) {
      sv    <- arr$a2[t]
      mu_s  <- f * exp(sv / 2)
      # d(mu)/dp = mu * sv/2;  d(var)/dp = sv * mu^2 * (2*exp(sv) - 1)
      if (!is.na(ka)) {
        dmu[t, ka]  <- mu_s * sv / 2
        dvar[t, ka] <- sv * mu_s^2 * (2 * exp(sv) - 1)
      }
      dv_df[t] <- 2 * mu_s * exp(sv / 2) * (exp(sv) - 1)
      next
    }

    cval <- arr$cc[t]
    fe   <- if (cval == 1) f else f^cval
    sa   <- arr$a2[t]; sb <- arr$b2[t]

    if (arr$form[t] == .ADM_RESID_COMBINED1) {
      # var = (a + b*f^c)^2, a = sqrt(sa), b = sqrt(sb); da/dp_a = a/2, db/dp_b = b/2
      a <- sqrt(sa); b <- sqrt(sb)
      s <- a + b * fe
      if (!is.na(ka)) dvar[t, ka] <- s * a
      if (!is.na(kb)) dvar[t, kb] <- s * b * fe
      if (!is.na(kc) && f > 0) dvar[t, kc] <- 2 * s * b * fe * log(f)
      dv_df[t] <- if (cval == 1) 2 * s * b else 2 * s * b * cval * f^(cval - 1)
    } else {
      # var = sa + sb * f^(2c)
      fe2 <- fe^2
      if (!is.na(ka)) dvar[t, ka] <- sa
      if (!is.na(kb)) dvar[t, kb] <- sb * fe2
      if (!is.na(kc) && f > 0) dvar[t, kc] <- sb * fe2 * 2 * log(f)
      dv_df[t] <- if (cval == 1) 2 * sb * f else 2 * sb * cval * f^(2 * cval - 1)
    }
  }
  list(dmu = dmu, dvar = dvar, dv_df = dv_df)
}

# d(mu)/d(f): the multiplicative mean scaling the residual applies. exp(s/2) on
# lnorm rows, 1 everywhere else. (adgh's `ls_vec`.)
.admResidMuScale <- function(arr) {
  ms <- rep(1, length(arr$form))
  ln <- arr$form == .ADM_RESID_LNORM
  if (any(ln)) ms[ln] <- exp(arr$a2[ln] / 2)
  ms
}

# Residual-parameter gradient, shared by every estimator.
#
#   d(-2LL)/d(p_k) = sum_t dNLL_dvar[t] * d(var[t])/d(p_k)
#                  + sum_t dNLL_dmu[t]  * d(mu[t]) /d(p_k)
#
# The mu-path is nonzero only for lnorm (the only form that shifts the mean).
# `mu_struct` is the STRUCTURAL mean (before any lnorm scaling).
.admSigmaGrad <- function(mu_struct, arr, pinfo, dNLL_dvar, dNLL_dmu) {
  d <- .admResidDeriv(mu_struct, arr, pinfo)
  drop(crossprod(d$dvar, dNLL_dvar) + crossprod(d$dmu, dNLL_dmu))
}

# How the residual couples a change in the STRUCTURAL mean into the objective,
# over and above the plain dNLL_dmu path the eta/omega kernels already carry:
#
#   * lnorm rescales the mean, so d(mu)/d(mu_struct) = exp(s/2), not 1;
#   * the residual VARIANCE depends on mu (prop/pow/combined), so moving mu also
#     moves diag(V) by dv_df.
#
# Returns a length-n_t vector to be added to dNLL_dmu before the eta/omega chain
# rule. Zero on purely additive rows. (admc's `sigma_mu_scale`.)
.admResidMuCoupling <- function(mu_struct, arr, pinfo, dNLL_dvar, dNLL_dmu) {
  d <- .admResidDeriv(mu_struct, arr, pinfo)
  dNLL_dmu * (.admResidMuScale(arr) - 1) + dNLL_dvar * d$dv_df
}
