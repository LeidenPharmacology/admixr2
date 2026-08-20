#' @importFrom stats cov dnorm pnorm qnorm rnorm runif setNames
#' @importFrom utils assignInNamespace head
#' @importFrom Rcpp sourceCpp
#' @useDynLib admixr2, .registration = TRUE
NULL

# Suppress R CMD check notes for ggplot2 NSE column names used in plot.admFit.
utils::globalVariables(c(
  "time", "pred_lo", "pred_hi", "pred_mean",
  "obs_lo", "obs_hi", "obs_mean",
  "t_col", "t_row", "value",
  "nll", "restart", "iter",
  "z", "z_label", "z_vjust",
  "resid", "lo", "hi"
))

`%||%` <- function(x, y) if (is.null(x)) y else x


# Assemble fullTheta vector in iniDf row order (thetas + sigma SDs + omega entries).
.admFullTheta <- function(pars, pinfo) {
  .ini <- pinfo$iniDf
  setNames(vapply(seq_len(nrow(.ini)), function(i) {
    nm <- .ini$name[i]
    if (.ini$fix[i])                      return(.ini$est[i])
    if (nm %in% names(pars$struct))       return(unname(pars$struct[nm]))
    # Residual params report on their iniDf scale: variances as an SD; a pow()
    # exponent and a t() degrees-of-freedom as themselves (.admSigmaNat has
    # already mapped log(nu - 2) back to nu).
    if (nm %in% names(pars$sigma_var)) {
      .k <- match(nm, pinfo$sigma_names)
      .v <- unname(pars$sigma_var[nm])
      return(if (.admSigmaRole(pinfo)[.k] == "var") sqrt(.v) else .v)
    }
    if (!is.na(.ini$neta1[i]))
      return(pars$omega[.ini$neta1[i], .ini$neta2[i]])
    .ini$est[i]
  }, double(1)), .ini$name)
}

# Compute objective function statistics for a completed fit.
#
# nobs = sum(n_subjects * n_times) -- observation RECORDS, which is nlmixr2's
# convention (nlmixr2est R/ofv.R builds BIC as objf + log(nobs) * df, and rxode2
# fills nobs by counting observation records). Keep it: users compare BIC across
# papers and against other software, and silently redefining it breaks that.
#
# An earlier revision changed this to sum(n) on the argument that the aggregate
# likelihood has SUBJECTS as its unit of independence. That argument is sound as
# far as it goes, but the fix was wrong: the field is split (nlme and SPSS use
# records; Monolix, saemix and SAS NLMIXED use subjects; NONMEM reports no BIC at
# all), so neither extreme is "the" convention, and unilaterally moving admixr2
# to one of them just makes its BIC incomparable in a different direction.
#
# The principled form is neither -- see .admBICh() below.
#
# extra_nobs: observations contributing to `objective` that are NOT in `studies`
# -- individual-level records in a joint fit (#60). The objective already counts
# them, so leaving them out of nobs made BIC penalise with the wrong sample size
# and stamped logLik with an nobs that did not match the value it carried.
.admCalcObjStats <- function(objective, npar, studies, extra_nobs = 0L,
                             pinfo = NULL) {
  nobs <- sum(vapply(studies, function(s)
    as.integer(s$n) * (s$n_total %||% length(s$times)), integer(1))) +
    as.integer(extra_nobs)
  ll   <- -objective / 2
  attr(ll, "df")   <- npar
  attr(ll, "nobs") <- nobs
  class(ll) <- "logLik"
  # AIC/BIC from stats' own methods on the logLik built two lines up, rather than
  # rewriting `objective + 2*npar` / `objective + log(nobs)*npar` here: `ll`
  # already carries the df and nobs those formulas need, so stats::AIC()/BIC()
  # have every input and the numbers admixr2 PRINTS are then the same ones a user
  # calling AIC(fit)/BIC(fit) gets, by construction rather than by agreement.
  # Bit-identical to the hand-written forms -- `-2 * (-objective/2)` is exact in
  # binary floating point -- checked with identical() over objective in
  # {0, 1e-16, 50, 1234.5678901234, 1e12, -3.25} x npar {1..11} x nobs {1..120000}.
  objDf <- data.frame(
    OBJF             = objective,
    AIC              = stats::AIC(ll),
    BIC              = stats::BIC(ll),
    "Log-likelihood" = as.numeric(ll),
    check.names      = FALSE
  )
  bh <- .admBICh(objective, studies, pinfo, extra_nobs)
  if (!is.null(bh)) objDf$BIC_h <- bh
  list(ll = ll, nobs = nobs, npar = npar, objDf = objDf, BIC_h = bh)
}

# Delattre, Lavielle & Poursat (2014) EJS 8:456-475, eq. (2.6): the BIC whose
# sample size is derived rather than chosen.
#
#   BIC_h = -2 log p(y|theta) + dim(theta_R) log N + dim(theta_F) log ntot
#
# theta_R -- parameters of the RANDOM-EFFECT distribution -> log N (subjects)
# theta_F -- everything else                               -> log ntot (records)
#
# The two familiar criteria are its endpoints: all-random gives BIC_N, all-fixed
# gives BIC_ntot. Their reasoning maps onto ADM directly -- an aggregate study is
# n_s independent random m-vectors compressed into sufficient statistics, and the
# objective is their exact density -- which is why the effective sample size for
# anything describing the random-effect distribution is the subject count.
#
# Returns NULL rather than a guess when the split cannot be determined.
# pinfo$struct_has_eta is NULL-vs-zero-row sensitive: NULL means "no mu-ref
# information at all", where a zero-row frame means "the information exists and
# says nothing is paired". Conflating them sends the split to an identity
# fallback -- see the note in CLAUDE.md.
.admBICh <- function(objective, studies, pinfo, extra_nobs = 0L) {
  if (is.null(pinfo) || is.null(pinfo$struct_has_eta)) return(NULL)
  n_subj <- sum(vapply(studies, function(s) as.integer(s$n), integer(1)))
  n_tot  <- sum(vapply(studies, function(s)
    as.integer(s$n) * (s$n_total %||% length(s$times)), integer(1))) +
    as.integer(extra_nobs)
  if (!is.finite(n_subj) || n_subj < 1L || !is.finite(n_tot) || n_tot < 1L)
    return(NULL)
  has_eta <- as.logical(pinfo$struct_has_eta)
  # theta_R: every omega entry, plus the structural thetas that carry an eta.
  dim_R <- length(pinfo$omega_par) + sum(has_eta, na.rm = TRUE)
  # theta_F: the sigmas, plus structural thetas with no eta (covariate
  # coefficients on a parameter without one land here, which is correct --
  # they are estimated at the observation rate, not the subject rate).
  dim_F <- length(pinfo$sigma_names) + sum(!has_eta, na.rm = TRUE)
  objective + dim_R * log(n_subj) + dim_F * log(n_tot)
}

# Bridge admControl/adirmcControl fields into foceiControl for nlmixr2 table machinery.
#
# `skip_cov` says which thetas nlmixr2est should NOT expect a standard error for.
# It fills its SE column by walking the thetas in iniDf order and taking the next
# entry of sqrt(diag(cov)) for each one it is not skipping, so this vector and the
# row order of the covariance are two halves of the same contract -- see
# .admCovSkip()/.admCovThetaOrder(). Left to nlmixr2est's own default when NULL.
.admToFoceiControl <- function(ctl, skip_cov = NULL) {
  .args <- list(
    rxControl          = ctl$rxControl,
    maxOuterIterations = 0L,
    maxInnerIterations = 0L,
    covMethod          = 0L,
    sumProd            = ctl$sumProd,
    optExpression      = ctl$optExpression,
    literalFix         = ctl$literalFix,
    scaleTo            = 0,
    calcTables         = ctl$calcTables,
    addProp            = ctl$addProp,
    interaction        = 0L,
    compress           = ctl$compress,
    ci                 = ctl$ci,
    sigdigTable        = ctl$sigdigTable)
  if (!is.null(skip_cov)) .args$skipCov <- skip_cov
  do.call(nlmixr2est::foceiControl, .args)
}

# LHS: one sample per stratum per dimension, independently permuted.
.lhsSample <- function(n, d) {
  m <- matrix(0, nrow = n, ncol = d)
  for (j in seq_len(d))
    m[, j] <- (sample.int(n) - 1L + runif(n)) / n
  m
}

# Generate z matrices -- one per study (seed must be set by caller).
# sampling: "sobol" (quasi-random), "lhs" (Latin hypercube), "rnorm" (iid normal).
.admMakeZ <- function(n_sim, pinfo, n_studies, sampling = "sobol") {
  replicate(n_studies, {
    if (pinfo$n_eta == 0L)
      return(matrix(0, nrow = n_sim, ncol = 1L))
    switch(sampling,
      sobol  = {
        z <- qnorm(randtoolbox::sobol(n = n_sim, dim = pinfo$n_eta))
        if (!is.matrix(z)) z <- matrix(z, ncol = 1L)
        z
      },
      halton = {
        z <- qnorm(randtoolbox::halton(n = n_sim, dim = pinfo$n_eta))
        if (!is.matrix(z)) z <- matrix(z, ncol = 1L)
        z
      },
      torus  = {
        z <- qnorm(randtoolbox::torus(n = n_sim, dim = pinfo$n_eta))
        if (!is.matrix(z)) z <- matrix(z, ncol = 1L)
        z
      },
      lhs    = qnorm(.lhsSample(n_sim, pinfo$n_eta)),
      rnorm  = matrix(rnorm(n_sim * pinfo$n_eta), nrow = n_sim),
      stop("admMakeZ: unknown sampling method '", sampling, "'", call. = FALSE)
    )
  }, simplify = FALSE)
}

# Pre-allocate params matrix list -- one per study.
# Matrix avoids data.frame list COW overhead: first col-write copies once,
# subsequent writes modify in-place. as.data.frame() wraps at rxSolve call site.
#
# Includes one residual-error placeholder column per observed output
# (rxerr.<output>). rxSolve requires every endpoint's rxerr parameter to be
# present in the params frame (its value is immaterial -- admixr2 adds residual
# error analytically and reads the structural prediction). There is exactly one
# rxerr per prediction line regardless of the error model (add/prop/lnorm or a
# combined add+prop), so the endpoint set is `unique(sigma_output)`; the mapping
# holds for named endpoints and linCmt (`rxerr.rxLinCmt`). All other model
# parameters (CMT, hard-coded constants) are left for rxSolve to fill from the
# model's own defaults -- do not add them here.
.admMakeParamsList <- function(n_sim, pinfo, n_studies = 1L) {
  so    <- unique(pinfo$sigma_output[!is.na(pinfo$sigma_output)])
  rxerr <- if (length(so)) paste0("rxerr.", so) else "rxerr.cp"
  col_names <- c(pinfo$struct_names, pinfo$eta_col_names,
                 pinfo$sigma_names, rxerr)
  replicate(n_studies, {
    m <- matrix(0, nrow = n_sim, ncol = length(col_names),
                dimnames = list(NULL, col_names))
    m[, rxerr] <- 1L
    m
  }, simplify = FALSE)
}

# -- FOCEI-style aligned progress output ---------------------------------------

# Column names: -2LL, back-transformed struct thetas, sigma SDs, omega diagonal variances.
.admProgressNames <- function(pinfo) {
  omega_diag_nms <- if (pinfo$n_eta > 0L)
    pinfo$eta_names[pinfo$chol_i[pinfo$chol_diag]]
  else character(0)
  c("-2LL", pinfo$struct_names, pinfo$sigma_names, omega_diag_nms)
}

# Bordered header block (separator + header row + separator). iter_w sets the
# label column width; data columns are max(name_width, 10). bottom=FALSE omits
# the closing separator so a phase divider can follow immediately.
.admProgressHeader <- function(pinfo, iter_w = 8L, bottom = TRUE) {
  nms <- .admProgressNames(pinfo)
  cws <- pmax(nchar(nms), 8L)
  sep <- paste0("+", strrep("-", iter_w + 2L), "+",
                paste(vapply(cws, function(w) strrep("-", w + 2L), character(1)),
                      collapse = "+"), "+")
  hdr <- paste0("| ", formatC("", width = iter_w), " | ",
                paste(mapply(formatC, nms, width = cws), collapse = " | "), " |")
  if (bottom) paste0(sep, "\n", hdr, "\n", sep) else paste0(sep, "\n", hdr)
}

# Full-width section divider matching .admProgressHeader outer width; pads label right with dashes.
.admProgressDivider <- function(label, pinfo, iter_w = 8L) {
  nms     <- .admProgressNames(pinfo)
  cws     <- pmax(nchar(nms), 8L)
  inner_w <- iter_w + 2L + sum(cws) + 3L * length(cws)
  n_dash  <- max(inner_w - 2L - nchar(label), 0L)
  paste0("+--", label, strrep("-", n_dash), "+")
}

.admProgressPhase <- function(phase_idx, phase_name, ph_step, pinfo, iter_w = 8L)
  .admProgressDivider(sprintf(" Phase %d: %s (+/-%.2f) ", phase_idx, phase_name, ph_step), pinfo, iter_w)

.admProgressRestart <- function(r, n_r, pinfo, iter_w = 8L)
  .admProgressDivider(sprintf(" Restart %d / %d ", r, n_r), pinfo, iter_w)

# One bordered data row aligned to the same column widths as .admProgressHeader.
.admProgressRow <- function(iter_label, nll, p, pinfo, iter_w = 8L) {
  pars <- tryCatch(.admUnpack(p, pinfo), error = function(e) NULL)
  if (is.null(pars)) return(NULL)
  nms <- .admProgressNames(pinfo)
  cws <- pmax(nchar(nms), 8L)
  struct_vals <- vapply(pinfo$struct_names, function(nm)
    .admBackTransform(pars$struct[[nm]], pinfo$struct_transforms[[nm]]), double(1))
  sigma_vals  <- sqrt(pars$sigma_var)
  omega_vals  <- if (pinfo$n_eta > 0L)
    diag(pars$omega)[pinfo$chol_i[pinfo$chol_diag]]
  else numeric(0)
  nll_str  <- local({
    s <- formatC(nll, format = "f", digits = 2)
    if (nchar(s) <= cws[1L]) formatC(nll, format = "f", digits = 2, width = cws[1L])
    else                     formatC(nll, format = "e", digits = 1, width = cws[1L])
  })
  par_strs <- mapply(function(x, w) formatC(x, format = "g", digits = 4, width = w),
                     c(struct_vals, sigma_vals, omega_vals), cws[-1L])
  val_strs <- c(nll_str, par_strs)
  paste0("| ", formatC(iter_label, width = iter_w, flag = "-"), " | ",
         paste(val_strs, collapse = " | "), " |")
}

# Timing row: label column shows elapsed time, all data columns blank.
.admProgressTimingRow <- function(sec, pinfo, iter_w = 8L) {
  nms    <- .admProgressNames(pinfo)
  cws    <- pmax(nchar(nms), 8L)
  blanks <- vapply(cws, function(w) formatC("", width = w), character(1))
  time_label <- if (sec >= 60) sprintf("%.1f min", sec / 60) else sprintf("%.1f sec", sec)
  paste0("| ", formatC(time_label, width = iter_w, flag = "-"), " | ",
         paste(blanks, collapse = " | "), " |")
}
