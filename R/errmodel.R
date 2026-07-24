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
.ADM_RESID_TBS       <- 3L   # logitNorm/probitNorm/boxCox/yeoJohnson (quadrature)
.ADM_RESID_POIS      <- 4L   # y ~ pois(f)        E = f,   Var = f
.ADM_RESID_BINOM     <- 5L   # y ~ binom(N, p)    E = N p, Var = N p (1-p)
.ADM_RESID_NBINOM    <- 6L   # y ~ nbinomMu(k, m) E = m,   Var = m + m^2/k
.ADM_RESID_BETA      <- 7L   # y ~ beta(b1, b2)   E = mu,  Var = mu(1-mu)/(1+phi)
.ADM_RESID_ORDINAL   <- 8L   # y ~ c(p1,..,pK-1)  multinomial category indicators

# Ordinal endpoints -----------------------------------------------------------
#
# `y ~ c(p1, p2)` is a MULTINOMIAL over K = length(args) + 1 categories, with the
# listed probabilities MARGINAL (not cumulative, not sequential) and the last
# category taking the remainder. Verified against rxord() directly:
#   pa=0.6 pb=0.3 -> 0.6014 0.2989 0.0998
#   pa=0.2 pb=0.5 -> 0.1992 0.5017 0.2991
#
# The aggregate observation is the vector of category INDICATORS, so a study
# contributes the observed proportion in each category. Per subject:
#
#   E[1_k]            = p_k
#   Var(1_k)          = p_k (1 - p_k)
#   Cov(1_j, 1_k)     = -p_j p_k        (j != k, SAME time point)
#
# i.e. diag(p) - p p' within a time, and zero residual covariance across times
# (given eta). That is exactly the shape .admResidApply already supports: a
# diagonal plus an off-diagonal `rmat`, the machinery added for ar(). The K
# probabilities stack like a joint same-subject unit, and each p_k is an ordinary
# derived model expression, so they get analytic sensitivities like any other.
.admOrdinalSpec <- function(ui, var) {
  a <- .admDistArgs(ui, var)
  if (is.null(a) || length(a) == 0L) return(NULL)
  if (!all(vapply(a, is.name, logical(1)))) return(NULL)
  vapply(a, function(x) paste(deparse(x), collapse = ""), character(1))
}

# Beta endpoints --------------------------------------------------------------
#
# `y ~ beta(b1, b2)` is the only endpoint whose distribution is defined by TWO
# model quantities rather than one prediction. In the mean/precision
# parameterisation mu = b1/(b1+b2), phi = b1+b2:
#
#   E[y|eta]   = mu
#   Var(y|eta) = mu (1 - mu) / (1 + phi)
#
# which is exactly binom's shape with N replaced by 1/(1+phi). So once the solve
# can hand back mu, everything downstream is the machinery already in place --
# E_eta[Var] = (mu_bar - mu_bar^2 - var_f)/(1 + phi) via the usual E[f^2] identity.
#
# The one restriction: PHI MUST NOT DEPEND ON ETA. With an eta-dependent
# precision, E_eta[mu(1-mu)/(1+phi)] no longer factors and would need the joint
# distribution of (mu, phi), which admixr2's (mu, var_f) summary cannot carry.
# That matches how beta regression is written in practice (a varying mean and a
# scalar precision); .admBetaPhiConst() checks it rather than assuming it.

# Verify that the solved beta precision phi = b1 + b2 really is eta-independent,
# and return the representative row.
#
# `phi_mat` is (n_draw x n_time) for ONE parameter configuration, i.e. its ROWS are
# random-effect draws. If phi does not depend on eta every row is identical and the
# first is representative -- which is what the solve paths assume when they collapse
# the matrix with [1L, ]. This function performs the check AND the collapse together,
# so the assumption cannot drift away from its use (for a long time the check was
# only ever promised in a comment; an eta-dependent precision silently populated the
# variance from whichever draw happened to land first, biasing the objective, the
# estimates and every SE with no error and no warning).
#
# Rows that legitimately differ (e.g. one row per STRUCTURAL configuration, as in
# .admSimulateRows) are NOT eta draws and must not be passed here.
.admBetaPhiConst <- function(phi_mat, tol = 1e-6) {
  if (is.null(phi_mat)) return(NULL)
  if (!is.matrix(phi_mat)) return(phi_mat)
  ref <- phi_mat[1L, , drop = TRUE]
  if (nrow(phi_mat) > 1L) {
    .sc  <- pmax(abs(ref), 1e-8)
    .dev <- abs(phi_mat - rep(ref, each = nrow(phi_mat))) /
      rep(.sc, each = nrow(phi_mat))
    .rel <- suppressWarnings(max(.dev[is.finite(.dev)], -Inf))
    if (is.finite(.rel) && .rel > tol)
      stop("admixr2: the beta precision phi = b1 + b2 varies across random-effect ",
           "draws (max relative variation ", format(.rel, digits = 3), ").\n",
           "  admixr2 matches E_eta[Var(y|eta)] = (mu - mu^2 - var_f)/(1 + phi), ",
           "which factors ONLY when phi is independent of the random effects; with ",
           "an eta-dependent precision it would need the joint distribution of ",
           "(mu, phi), which the (mu, var_f) summary cannot carry.\n",
           "  Write the endpoint with a scalar precision, e.g.\n",
           "    phi <- exp(tphi); b1 <- mu * phi; b2 <- (1 - mu) * phi\n",
           "  (a varying MEAN is fine -- only phi must not carry an eta).",
           call. = FALSE)
  }
  ref
}

.admBetaSpec <- function(ui, var) {
  a <- .admDistArgs(ui, var)
  if (is.null(a) || length(a) < 2L) return(NULL)
  nm <- vapply(a[1:2], function(x) paste(deparse(x), collapse = ""), character(1))
  if (!all(vapply(a[1:2], is.name, logical(1)))) return(NULL)
  list(b1 = nm[[1L]], b2 = nm[[2L]])
}

# Count endpoints -------------------------------------------------------------
#
# nlmixr2 gives a count endpoint a DIFFERENT shape from a residual-error one: the
# mean is the distribution's ARGUMENT, not the endpoint, so predDf$var names the
# DV (`y`) while the quantity admixr2 must solve for is `cp` in `y ~ pois(cp)`.
# There are no iniDf error rows at all -- which is exactly why these models used
# to sail through every gate and fit with zero residual variance.
#
# Their conditional moments are closed form, so once the right column is solved
# they drop straight into the law of total variance already in place:
#
#   pois(f)        E[y|eta] = f          Var(y|eta) = f
#   binom(N, p)    E[y|eta] = N p        Var(y|eta) = N p (1 - p)
#   nbinomMu(k, m) E[y|eta] = m          Var(y|eta) = m + m^2 / k
#
# and E_eta[.] uses the same E[f^2] = mu^2 + var_f identity as every other form.
# N / k must be constants (literal or fix()ed): an ESTIMATED size enters the
# objective only through the variance, never through the ODE solve, so its
# gradient would have to travel a path the struct-theta machinery does not have.

# Arguments of the distribution call on an endpoint's model line, e.g.
# `y ~ binom(20, cp)` -> list(20, cp). NULL when the line cannot be found.
.admDistArgs <- function(ui, var) {
  lst <- tryCatch(ui$lstExpr, error = function(e) NULL)
  if (is.null(lst)) return(NULL)
  for (e in lst) {
    if (is.call(e) && length(e) == 3L && identical(as.character(e[[1L]]), "~") &&
        identical(paste(deparse(e[[2L]]), collapse = ""), var)) {
      rhs <- e[[3L]]
      if (is.call(rhs)) return(as.list(rhs)[-1L])
    }
  }
  NULL
}

# Resolve a count endpoint's MEAN argument to a model variable name, and any
# size/N argument to a constant. Returns NULL when the shape is not usable.
.admCountSpec <- function(ui, var, dist) {
  a <- .admDistArgs(ui, var)
  if (is.null(a)) return(NULL)
  .const <- function(x) {
    # A literal is stored in the AST as an actual R number: `binom(20, p)` gives the
    # double 20, `binom(20L, p)` the integer 20L. deparse(20L) is "20L", and
    # as.numeric("20L") is NA -- so an integer-suffixed size was misread as
    # non-constant and refused with advice to fix() a parameter that does not exist.
    # Read a bare numeric literal (double OR integer) directly, before deparsing.
    if (is.numeric(x) && length(x) == 1L) return(as.numeric(x))
    v <- suppressWarnings(as.numeric(paste(deparse(x), collapse = "")))
    if (!is.na(v)) return(v)
    nm  <- paste(deparse(x), collapse = "")
    idf <- tryCatch(ui$iniDf, error = function(e) NULL)
    if (!is.null(idf) && nm %in% idf$name) {
      r <- idf[idf$name == nm, , drop = FALSE][1L, ]
      if (isTRUE(r$fix)) return(as.numeric(r$est))
    }
    # A constant written in the MODEL BLOCK -- `nt <- 20; y ~ binom(nt, p)` -- is
    # every bit as constant as the literal, and is how a number of trials is
    # usually written. Refusing it sent the user to fix() a parameter that does not
    # exist. Only a bare numeric assignment counts: anything computed could depend
    # on a theta, which is the case that genuinely has no gradient path.
    lst <- tryCatch(ui$lstExpr, error = function(e) NULL)
    for (e in lst %||% list()) {
      if (is.call(e) && length(e) == 3L &&
          as.character(e[[1L]]) %in% c("<-", "=") &&
          identical(paste(deparse(e[[2L]]), collapse = ""), nm)) {
        lit <- suppressWarnings(as.numeric(paste(deparse(e[[3L]]), collapse = "")))
        if (!is.na(lit)) return(lit)
      }
    }
    NA_real_
  }
  if (dist %in% c("pois", "dpois")) {
    if (length(a) < 1L) return(NULL)
    return(list(mean_var = paste(deparse(a[[1L]]), collapse = ""), size = NA_real_))
  }
  if (length(a) < 2L) return(NULL)
  list(size = .const(a[[1L]]),
       mean_var = paste(deparse(a[[2L]]), collapse = ""))
}

# iniDf$err vocabularies (rxode2 spells the same model several ways).
.ADM_ERR_ADD     <- c("add", "norm", "dnorm")
.ADM_ERR_PROP    <- c("prop", "propT", "propF")
.ADM_ERR_POW     <- c("pow", "powT", "powF")
.ADM_ERR_POW_EXP <- c("pow2", "powT2", "powF2")   # the EXPONENT row of pow(b, c)
.ADM_ERR_LNORM   <- c("lnorm", "dlnorm", "logn", "dlogn")
.ADM_ERR_T       <- c("t", "dt")                  # the DEGREES-OF-FREEDOM row of t(nu)
.ADM_ERR_AR      <- "ar"                          # the CORRELATION row of ar(rho)
.ADM_ERR_TBS_SD  <- c("logitNorm", "probitNorm")  # the SD row of a bounded transform
.ADM_ERR_TBS_LAM <- c("boxCox", "tbs", "yeoJohnson", "tbsYj")   # the LAMBDA row
.ADM_ERR_NB_SIZE <- c("nbinom", "dnbinom", "nbinomMu", "dnbinomMu")  # the SIZE row
# Discrete endpoints put their OWN name in iniDf$err when a distribution argument
# is estimated (an estimated binomial N emits err = "binom"). Those are supported
# forms, so they must not trip the generic gate in .admParseIniDf() -- which fires
# BEFORE .admBuildResidSpecs() and printed a .ADM_ERR_WHY reason that is no longer
# true ("binomial data have no separate residual parameter... admixr2 would fit
# with no residual variance"). Listing them here lets the accurate, case-specific
# refusal in .admBuildResidSpecs() (a non-constant size/N has no gradient path) be
# the message the user actually sees.
.ADM_ERR_COUNT   <- c("pois", "dpois", "binom", "dbinom", "nbinomMu", "dnbinomMu",
                      "beta", "dbeta", "ordinal", "dordinal")
.ADM_ERR_KNOWN   <- c(.ADM_ERR_ADD, .ADM_ERR_PROP, .ADM_ERR_POW,
                      .ADM_ERR_POW_EXP, .ADM_ERR_LNORM, .ADM_ERR_T, .ADM_ERR_AR,
                      .ADM_ERR_TBS_SD, .ADM_ERR_TBS_LAM, .ADM_ERR_NB_SIZE,
                      .ADM_ERR_COUNT)

# AR(1) residual correlation matrix over the observation times. rxode2 implements
# ar(rho) as a stationary CONTINUOUS-time AR(1), so the correlation between two
# observations depends on their time GAP, not their index.
.admARCor <- function(rho, times) {
  if (is.na(rho) || rho <= 0) return(NULL)
  d <- abs(outer(times, times, "-"))
  r <- rho^d
  diag(r) <- 1
  r
}

# Student-t variance multiplier: Var(scale * T_nu) = scale^2 * nu/(nu-2).
.admTMult <- function(nu) nu / (nu - 2)

# Which ordinal rows are observations of the SAME time point?
#
# A joint ordinal unit stacks one block per category, so the categories of one
# time point are rows in different blocks and the row-time vector is what
# identifies them. Three places need that grouping and must agree exactly:
# .admResidApply (which emits the -p_j*p_k cross term), .admResidVChain (which
# zeroes d(V_pred)/d(V_struct) for those same entries, because the structural
# covariance cancels out of them) and .admResidMuCoupling (the mu-coupling term).
#
# ONE definition, because two of them disagreeing is worse than both being wrong:
# the objective then carries a cross term the gradient does not know about, and
# the optimizer descends a direction the function does not follow. That is exactly
# what happened when the tolerance grouping below was first added at one site only.
#
# Grouped by TOLERANCE. The row times come from the per-category blocks, i.e. from
# independent user inputs: `seq(0.1, 0.7, by = 0.2)` and `c(0.1, 0.3, 0.5, 0.7)`
# are the same grid to a reader and differ in the last bit to match(), which put
# the two categories in different groups and silently dropped the cross term.
# NA row times (a hand-built block with no `times`) must NOT collapse into one
# group -- match(NA, ...) matches, which would invent cross terms.
#
# Package-level, not inlined: this is the "one place that knows" the rule. A
# dev-mode mirai daemon cannot see a NEW binding in the installed namespace, so
# run devtools::install() before testing parallel restarts (the documented rule;
# it applies to any new function).
.admOrdTimeGroup <- function(times, or_) {
  out <- rep(NA_integer_, length(times))
  keep <- which(or_ & !is.na(times))
  if (!length(keep)) return(out)
  k <- 0L; ref <- NA_real_
  for (i in keep[order(times[keep])]) {
    if (is.na(ref) || abs(times[i] - ref) > 1e-8 * max(1, abs(times[i]), abs(ref))) {
      k <- k + 1L; ref <- times[i]      # compare against the group's FIRST member,
    }                                   # so groups cannot chain
    out[i] <- k
  }
  out
}

# rxode2's _eps (`#define _eps sqrt(DBL_EPSILON)`, rxode2.h). Every clamp in the
# transform code below uses THIS, so admixr2 and the solve agree on where a
# transform runs out of support. NOT the same constant as safeLog/safeZero/safePow,
# which use DBL_EPSILON itself -- do not conflate them.
.ADM_EPS <- sqrt(.Machine$double.eps)

# Ceiling on a second-order delta-expansion term before it is treated as diverged.
# The expansion E[f^k] ~ mu^k + k(k-1)/2 * mu^(k-2) * var is only meaningful while
# the correction is small against the leading term; for pow() with c < 1 near a zero
# prediction it is not, and an uncapped term produced a negative variance.
.ADM_MOM_CAP <- 1 / sqrt(.Machine$double.eps)

# rxode2's transform selector codes (`rx_yj_` in the emitted simulation model).
.ADM_TBS_YJ <- c(boxCox = 0L, tbs = 0L, yeoJohnson = 1L, tbsYj = 1L,
                 untransformed = 2L, logit = 4L, probit = 6L)

# Forward / inverse transform: rxode2's OWN kernel, not a re-implementation.
#
# `rxode2::.rxTransform(x, lambda, low, high, transform, inverse)` is exported and
# is what every one of rxode2's own boxCox()/yeoJohnson()/logit()/probit() (and
# their inverses) calls; it bottoms out in `.Call(_rxode2_powerD, ...)`, i.e. the
# very C routine the SOLVE transforms with. Its `transform` codes are the ones
# admixr2 already uses (0 boxCox, 1 yeoJohnson, 2 untransformed, 4 logit,
# 6 probit), because admixr2 took them from there.
#
# These used to be a careful line-by-line port of `_powerD`/`_powerDi` -- roughly
# ninety lines of branch order, clamps and short-circuits, each one a documented
# gotcha (lambda == 1 short-circuits BEFORE the clamp; log1p not log(x+1); an
# out-of-bounds logit argument is NaN, not +-Inf). It agreed with the kernel
# exactly (0 mismatches over every code x lambda x bounds combination) but it
# could only ever agree by re-deriving, which is what test-transform-vs-rxode2.R
# was written to police. Calling the kernel makes the agreement structural.
#
# The two things the kernel does NOT do for us, and the only reason these are
# wrappers rather than direct calls:
#
#   1. It DROPS dim(). The solve paths hand these an n_sim x n_t MATRIX, and a
#      dropped dim turned cp_mat into a flat vector so every downstream
#      colMeans()/sweep() produced NA.
#   2. Its scalar-argument fast path asserts `any.missing = FALSE` on `lambda`,
#      `low` and `high` (not on `x` -- non-finite x is handled and returns NA).
#
# Not routed through the public boxCox()/yeoJohnson()/logit()/probit() wrappers:
# those add their own `checkmate::assertNumeric(x, lower = 0, any.missing = FALSE)`
# and so ERROR on the NA/Inf that the quadrature's +-12 SD tail nodes legitimately
# produce, where the kernel returns NA.
.admTBSxf <- function(x, lam, yj, lo, hi, inverse) {
  d   <- dim(x)
  out <- rxode2::.rxTransform(x, lam, lo, hi, as.integer(yj), inverse)
  dim(out) <- d
  out
}

.admTBS  <- function(y, lam, yj, lo, hi) .admTBSxf(y, lam, yj, lo, hi, FALSE)

# Box-Cox inverse with rxode2's out-of-support rule. (lam*z + 1)^(1/lam) is only
# defined while lam*z + 1 > 0 -- the transform has bounded support -- and an
# 81-node quadrature grid reaches +-12 SD, so tail nodes fall outside it and
# return NaN, which propagates into the whole sigma gradient. From _powerDi()
# (rxode2.h, case 0):
#
#     x0 = x*lambda + 1.0;
#     if (x0 <= _eps) return _eps;                 // the RESULT is _eps
#     ret = pow(x0, 1.0/lambda);
#     if (ISNA(ret)) return _eps;
#
# i.e. clamp the RESULT to _eps, not the power's base to some floor. (An earlier
# version here clamped the base to 1e-12 and then raised it to 1/lambda, giving
# 1e-24 at lambda = 0.5 and 1e-6 at lambda = 2: two different "zeros", neither
# one rxode2's.) `pw` is 1/lambda for the value and 1/lambda - 1 for its
# derivative; rxode2 guards both the same way.
.admTBSp <- function(base, pw) {
  out <- base^pw
  # rxode2 tests `ISNA(ret)`, NOT !is.finite: an overflow to +Inf is RETURNED.
  # Clamping it turned the largest representable magnitude into the smallest.
  out[is.na(out) | base <= .ADM_EPS] <- .ADM_EPS
  out
}

# Inverse transform: modelling scale -> natural scale. See .admTBSxf() above.
.admTBSi <- function(z, lam, yj, lo, hi) .admTBSxf(z, lam, yj, lo, hi, TRUE)

# Derivative of the INVERSE transform, g'(z). Everything else follows from it:
# h'(f) = 1 / g'(h(f)) by the inverse-function theorem, so no separate forward
# derivative is needed and the two can never drift apart.
.admTBSid <- function(z, lam, yj, lo, hi) {
  if (!all(is.finite(z))) {                       # _powerDD's R_finite guard
    out <- rep(NA_real_, length(z)); dim(out) <- dim(z)
    ok  <- is.finite(z)
    if (any(ok)) out[ok] <- .admTBSid(z[ok], lam, yj, lo, hi)
    return(out)
  }
  # lambda == 1 is the identity transform, so g'(z) == 1 EXACTLY. Routing it
  # through .admTBSp computed base^0 = 1 and then CLAMPED it to _eps for z <= -1 --
  # a derivative wrong by a factor of 6.7e7.
  if (yj == 2L || ((yj == 0L || yj == 1L) && lam == 1)) {
    out <- rep(1, length(z)); dim(out) <- dim(z); return(out)
  }
  if (yj == 0L) return(if (lam == 0) exp(z) else .admTBSp(lam * z + 1, 1 / lam - 1))
  if (yj == 4L) { e <- 1 / (1 + exp(-z)); return((hi - lo) * e * (1 - e)) }
  if (yj == 6L) return((hi - lo) * stats::dnorm(z))
  if (yj == 1L) {
    out <- numeric(length(z)); p <- z >= 0
    out[p]  <- if (lam == 0) exp(z[p]) else (lam * z[p] + 1)^(1 / lam - 1)
    # NOTE: rxode2's _powerDD case 1 returns -1/(1-x) here, which is a SIGN ERROR:
    # h(x) = -log1p(-x) gives h'(x) = +1/(1-x), and a finite difference of rxode2's
    # OWN yeoJohnson() confirms the positive sign (+0.1667 at x = -5 vs _powerDD's
    # -0.1667). admixr2 keeps the correct sign deliberately -- matching rxode2
    # everywhere else does not extend to reproducing a defect.
    out[!p] <- if (lam == 2) exp(-z[!p])
               else (1 - (2 - lam) * z[!p])^(1 / (2 - lam) - 1)
    dim(out) <- dim(z)                  # see .admTBS(): a dropped dim reads as NA
    return(out)
  }
  stop("unsupported transform code ", yj)                                # nocov
}

# Conditional moments of y = g(h(f) + sd*eps), eps ~ N(0,1), by Gauss-Hermite
# quadrature (the same standard-normal nodes adgh uses for the eta integral).
# Returns m = E[y|eta] and v = Var(y|eta), vectorised over f (and over sd, which
# is f-dependent whenever the endpoint carries a prop()/pow() term).
#
# Node count, measured against stats::integrate() at rel.tol 1e-13 (an independent
# rule, not a bigger version of this one). Worst case over boxCox/yeoJohnson/
# logit/probit x residual sd in {0.5, 1, 2, 3}:
#
#   n     5      15       31       61       81      121
#   err  3.3e-1  5.7e-2   4.5e-3   6.5e-5   5.0e-5  2.8e-5
#
# Cost is linear in n in isolation -- ~50 us (n=15), 150 (31), 300 (81), 500 (121)
# per call for an 8-row study -- but negligible beside the ODE solve: a full NLL
# evaluation measured 0.750 s per 60 evaluations at BOTH 31 and 81 nodes. So
# `resid_nodes` is an ACCURACY dial, not a speed one.
#
# The worst case is dominated ENTIRELY by sd = 3, and there by boxCox/yeoJohnson,
# whose inverse hits its bounded support and is clamped (see .admTBSp) -- a kink
# Gauss-Hermite converges on slowly, which is why the error plateaus around 3e-5
# rather than continuing down. It is not quadrature error at that point. At sd <= 1,
# a realistic residual on a transformed scale, n = 31 already gives 1e-7 or better
# and logit/probit reach 1e-13.
#
# (An earlier note here claimed "81 keeps even an extreme sd = 3 below 1e-6". That
# was measured on logit only -- 3.6e-08 -- and does not hold for boxCox, which is
# 5.0e-05. Corrected rather than dropped, because the number was load-bearing for
# the choice of default.)
#
# 81 stays the DEFAULT because it is safe across the whole grid above; it is now
# configurable per fit via the `resid_nodes` control argument, since a user with a
# small residual SD can halve the cost and one with a saturating endpoint and a
# large SD may want more.
.ADM_TBS_NODES <- 81L

# Moments PLUS their exact partials w.r.t. f and the residual sd, by
# differentiating the quadrature term by term rather than differencing its
# result. Nesting finite differences here (an outer FD over an inner FD-derived
# m'/m'') amplified the inner step's error by ~1e6 and made the boxCox sigma
# gradient wrong by 300%, so the chain is analytic wherever a closed form exists:
#
#   z_q      = h(f) + sd*x_q
#   dm/df    = h'(f) * sum_q w_q g'(z_q)          h'(f) = 1 / g'(h(f))
#   dm/dsd   =          sum_q w_q g'(z_q) x_q
#   dv/df    = h'(f) * sum_q w_q 2 g g' - 2 m dm/df
#   dv/dsd   =          sum_q w_q 2 g g' x_q - 2 m dm/dsd
.admTBSMomentsD <- function(f, sd, lam, yj, lo, hi, nodes = .ADM_TBS_NODES) {
  gq <- .adghNodes1(nodes)
  hz <- .admTBS(f, lam, yj, lo, hi)
  hp <- 1 / .admTBSid(hz, lam, yj, lo, hi)          # h'(f)
  n  <- length(f)
  # The transform is evaluated for the WHOLE node grid in one call rather than
  # once per node. `.admTBSi` bottoms out in rxode2's C kernel, whose R-level
  # preamble (argument checks + a non-finite scan) costs far more than the kernel
  # itself on a length-n vector: measured 81 short calls = 2.48 ms against one
  # stacked call = 0.05 ms for n = 8 at 81 nodes, i.e. the loop, not the maths,
  # was the cost. z_all[i, q] = h(f_i) + sd_i * x_q (hz recycles down columns).
  #
  # The per-node ACCUMULATION below is deliberately left as a loop over columns:
  # summing in a different order would move the moments by an ulp and the NLL with
  # them, and this stays bit-identical to the previous implementation.
  z_all  <- hz + outer(sd, gq$x)
  yq_all <- .admTBSi(z_all, lam, yj, lo, hi)
  gp_all <- .admTBSid(z_all, lam, yj, lo, hi)
  # The +-12 SD tail nodes can overflow the inverse transform to +-Inf/NaN --
  # yeoJohnson at large lambda over a low prediction is the case that bites, where
  # m2 - m*m becomes Inf - Inf = NaN. Such a node's GH weight is ~1e-30, so dropping
  # it from the weighted sum leaves the moments unchanged to machine precision; NOT
  # dropping it turns the whole moment (and thus the NLL AND the gradient) NaN. The
  # optimizer tolerates a NaN objective -- it just rejects the step -- but nloptr
  # errors on a NaN gradient ("missing value where TRUE/FALSE needed"), so the fit
  # crashed instead of the line search backing off. Zeroing the offending node makes
  # both consistent and finite. For a well-behaved endpoint no node is non-finite,
  # so `any(.bad)` is FALSE and this is a bit-identical no-op.
  .bad <- !is.finite(yq_all) | !is.finite(gp_all)
  if (any(.bad)) { yq_all[.bad] <- 0; gp_all[.bad] <- 0 }
  m <- m2 <- dmf <- dms <- dvf0 <- dvs0 <- numeric(n)
  for (q in seq_along(gq$x)) {
    yq <- yq_all[, q]
    gp <- gp_all[, q]
    w  <- gq$w[q]
    m    <- m    + w * yq
    m2   <- m2   + w * yq * yq
    dmf  <- dmf  + w * gp
    dms  <- dms  + w * gp * gq$x[q]
    dvf0 <- dvf0 + w * 2 * yq * gp
    dvs0 <- dvs0 + w * 2 * yq * gp * gq$x[q]
  }
  dm_df <- hp * dmf
  dm_ds <- dms
  list(m = m, v = pmax(m2 - m * m, 0),
       dm_df = dm_df, dm_ds = dm_ds,
       dv_df = hp * dvf0 - 2 * m * dm_df,
       dv_ds = dvs0      - 2 * m * dm_ds)
}

# Why a residual error model is REFUSED.
#
# A list(), not c(): `[[` on a named CHARACTER VECTOR errors for a missing name,
# so the `%||%` fallback at every call site only behaves if a miss yields NULL.
#
# Only genuinely-refused types belong here. This table used to carry entries for
# pois, binom, nbinomMu, beta, ordinal, logitNorm, probitNorm, boxCox, tbs,
# yeoJohnson, tbsYj and ar -- every one of which is now implemented, so those
# reasons were stale documentation being printed to users as if it were current.
# An entry here is a promise that admixr2 will not fit the model; keep it that way.
.ADM_ERR_WHY <- list(
  cauchy     = paste0(
    "the Cauchy is a STABLE distribution: the mean of n draws has the SAME ",
    "distribution as one draw,\nso a study mean carries no more information than ",
    "a single observation and has no finite\nvariance. No aggregate moment exists ",
    "to match -- this is not a limitation of admixr2's\napproximation"),
  dcauchy    = "see cauchy()",
  propF      = paste0(
    "propF()/powF() scale the residual by a user-supplied model variable, which is ",
    "an\nindividual-level quantity admixr2 cannot recover from the aggregate mean"),
  powF       = "see propF()",
  # Distributions with no natural-scale mean/variance admixr2 can match.
  chisq      = "a chi-squared endpoint has no residual-error parameterisation admixr2 can score",
  dexp       = "an exponential endpoint has no residual-error parameterisation admixr2 can score",
  f          = "an F-distributed endpoint has no residual-error parameterisation admixr2 can score",
  geom       = "a geometric endpoint has no residual-error parameterisation admixr2 can score",
  unif       = "a uniform endpoint has no residual-error parameterisation admixr2 can score",
  weibull    = "a Weibull endpoint has no residual-error parameterisation admixr2 can score",
  dgamma     = "a gamma endpoint has no residual-error parameterisation admixr2 can score"
)

# The supported-model reference block, shown on every refusal.
.admSupportedErrText <- function() paste(
  "Supported residual error models (f = the model prediction):",
  "  add(a)              var = a^2",
  "  prop(b)             var = (b*f)^2",
  "  pow(b, c)           var = (b*f^c)^2",
  "  lnorm(a)            lognormal, moment-matched",
  "  add(a) + prop(b)    var = a^2 + (b*f)^2   [combined2, the default]",
  "                      var = (a + b*f)^2     [combined1, via combined1()]",
  "  add(a) + pow(b, c)  either combined form",
  "  ... + t(nu)         any of the above with Student-t residuals (nu > 2):",
  "                      the scale family, var = <above> * nu/(nu-2)",
  sep = "\n")

# One consistent refusal message.
#   what   what the model asked for, e.g. "logitNorm()"
#   why    the reason, as a clause ("residuals are normal on the LOGIT scale, ...")
#   fix    optional concrete suggestion
.admStopErrModel <- function(endpoint, what, why, fix = NULL) {
  .para <- function(txt, prefix = "  ")
    paste(strwrap(txt, width = 76, prefix = prefix), collapse = "\n")

  ep <- if (is.null(endpoint) || is.na(endpoint)) ""
        else paste0(" for endpoint '", endpoint, "'")

  stop(
    "Unsupported residual error model", ep, ": ", what, ".\n\n",
    .para(paste0("Why: ", why, ".")), "\n\n",
    .para(paste0("admixr2 fits AGGREGATE data -- each study contributes a mean and a ",
                 "covariance, scored as a multivariate normal -- so the residual model ",
                 "must reduce to a mean and a variance on the natural scale.")), "\n\n",
    if (!is.null(fix)) paste0(.para(paste0("Fix: ", fix)), "\n\n") else "",
    .admSupportedErrText(), "\n\n",
    .para(paste0("Note: earlier versions of admixr2 accepted this model with a warning ",
                 "and then fitted it as ADDITIVE error. Any results carried over from ",
                 "that are not the model you specified.")),
    call. = FALSE)
}

# Validate every endpoint's DISTRIBUTION, independently of whether it has any
# residual-error parameters.
#
# This gate has to be separate from .admBuildResidSpecs(), which returns early
# when `sigma_names` is empty -- and a count/categorical endpoint has NO iniDf
# error rows at all (`y ~ pois(cp)` puts the rate in the distribution argument;
# `cp ~ c(...)` emits only fixed probability rows). So those models sailed
# through every check and produced a pinfo with sigma_names = character(0),
# i.e. a converged fit with ZERO residual variance and no warning. That is the
# same silent-wrong-model failure as the historical pow() bug, so it is refused
# here on the authoritative field (predDf$distribution) before anything else.
#
# `t` is supported (a scale family, see .ADM_ERR_T); everything else non-normal
# is not, and .ADM_ERR_WHY supplies the per-distribution reason.
.admCheckEndpointDist <- function(ui) {
  predDf <- if (!is.null(ui)) tryCatch(as.data.frame(ui$predDf), error = function(e) NULL) else NULL
  if (is.null(predDf) || nrow(predDf) == 0L) return(invisible(NULL))
  if (!"distribution" %in% names(predDf)) return(invisible(NULL))
  for (i in seq_len(nrow(predDf))) {
    dist <- as.character(predDf$distribution[i] %||% "norm")
    if (dist %in% c("norm", "dnorm") || dist %in% .ADM_ERR_T) next
    if (dist %in% c("pois", "dpois", "binom", "dbinom",
                    "nbinomMu", "dnbinomMu", "beta", "dbeta",
                    "ordinal", "dordinal")) next   # closed-form moments
    ep <- as.character(predDf$cond[i] %||% predDf$var[i] %||% NA_character_)
    .admStopErrModel(
      ep, paste0("a '", dist, "' endpoint"),
      .ADM_ERR_WHY[[dist]] %||%
        paste0("'", dist, "' is not a distribution admixr2 can reduce to an aggregate ",
               "mean and variance"))
  }
  invisible(NULL)
}

# Build one residual spec per endpoint from ui$predDf.
#
# predDf is authoritative and iniDf$err alone is not: `errType` tells us which
# terms are present, `addProp` picks combined1 vs combined2, `transform` and
# `distribution` say whether the endpoint is representable at all, and `errTypeF`
# says which scale the proportional/power term multiplies.
#
# Returns NULL when there is no predDf (Tier-1 mock iniDf / no ui); callers then
# fall back to the legacy sigma_is_prop/sigma_is_lnorm flags.
# Constructor for a residual-error SPEC: the complete field set and its defaults
# in ONE place. `...` overrides any default.
#
# The four branches that build specs (ordinal, beta, count, and the
# residual-parameter endpoint) each used to spell out the whole ~18-field literal
# by hand, and three of the four already disagreed about which fields they set --
# only ordinal set add_fixed/prop_fixed/pow_fixed, only beta set out_pair, only
# ordinal set ord_p -- which is precisely why .admResidRows() has to read several
# of them through `%||% NA_real_`. Adding a field meant a four-site edit, and
# missing one let the endpoint fall back to the legacy spec (form 0, i.e. NO
# residual variance at all) -- the recurring bug this file's own comment above
# describes, of which beta was the fifth instance.
.admNewSpec <- function(output, form, ...) {
  .d <- list(
    output     = output,      form       = form,
    k_add      = NA_integer_, k_prop     = NA_integer_, k_pow     = NA_integer_,
    k_tdf      = NA_integer_, tdf_fixed  = NA_real_,
    k_ar       = NA_integer_, ar_fixed   = NA_real_,
    k_lam      = NA_integer_, lam_fixed  = NA_real_,
    yj         = 2L,          tr_lo      = 0,          tr_hi     = 1,
    csize      = NA_real_,    k_size     = NA_integer_,
    add_fixed  = NA_real_,    prop_fixed = NA_real_,   pow_fixed = NA_real_,
    # TBS only: form .ADM_RESID_TBS loses the combined1/combined2 distinction and
    # the errTypeF, so they ride alongside. FALSE matches what the other branches
    # got implicitly, since consumers read them through isTRUE().
    tbs_c1     = FALSE,       tbs_ftr    = FALSE,
    out_pair   = NULL,        ord_p      = NULL,       dv_name   = NA_character_)
  .o <- list(...)
  if (length(.o)) .d[names(.o)] <- .o
  .d
}

.admBuildResidSpecs <- function(ui, sigma_rows, sigma_names) {
  predDf <- if (!is.null(ui)) tryCatch(as.data.frame(ui$predDf), error = function(e) NULL) else NULL
  if (is.null(predDf) || nrow(predDf) == 0L) return(NULL)
  # NOT `if (length(sigma_names) == 0L) return(NULL)`: pois/binom have no residual
  # parameter at all, so bailing here left them with no spec -- which is exactly
  # how they ended up being fitted with zero residual variance.
  # Any endpoint whose moments come from the model rather than from residual
  # parameters. Missing one here is the recurring bug in this file: the endpoint
  # silently falls back to the legacy spec (form 0, no residual at all). beta was
  # the fifth instance.
  .no_sigma_ok <- "distribution" %in% names(predDf) &&
    any(as.character(predDf$distribution) %in%
          c("pois", "dpois", "binom", "dbinom", "nbinomMu", "dnbinomMu",
            "beta", "dbeta", "ordinal", "dordinal"))
  # ... and so is an endpoint whose residual parameters are ALL fix()ed: it has no
  # optimizer slots, but it certainly has a residual variance. Bailing here is what
  # made `cp ~ add(a)` with `a <- fix(0.7)` fit with none.
  .all_fixed_resid <- {
    .fi <- tryCatch(as.data.frame(ui$iniDf), error = function(e) NULL)
    !is.null(.fi) && all(c("err", "fix") %in% names(.fi)) &&
      any(!is.na(.fi$err) & .fi$fix & .fi$err %in% .ADM_ERR_KNOWN)
  }
  if (length(sigma_names) == 0L && !.no_sigma_ok && !.all_fixed_resid) return(NULL)

  err  <- sigma_rows$err
  cond <- if ("condition" %in% names(sigma_rows))
    as.character(sigma_rows$condition) else rep(NA_character_, length(sigma_names))

  # The FULL iniDf, FIXED error rows included. sigma_rows has had them removed
  # (.admParseIniDf keeps only estimated parameters), but a fixed t() df still has
  # to scale the variance -- see the tdf_fixed branch below.
  full_ini <- tryCatch(as.data.frame(ui$iniDf), error = function(e) NULL)
  if (is.null(full_ini) || !all(c("err", "fix", "est") %in% names(full_ini)))
    full_ini <- data.frame(err = character(0), fix = logical(0),
                           est = numeric(0), condition = character(0),
                           stringsAsFactors = FALSE)
  if (!"condition" %in% names(full_ini))
    full_ini$condition <- rep(NA_character_, nrow(full_ini))
  full_ini$condition <- as.character(full_ini$condition)

  # nlmixr2's addProp default resolves to combined2.
  add_prop_default <- getOption("rxode2.addProp", "combined2")

  specs <- list(); keys <- character(0)
  for (i in seq_len(nrow(predDf))) {
    ep   <- as.character(predDf$cond[i])
    tr   <- as.character(predDf$transform[i] %||% "untransformed")
    etf  <- as.character(predDf$errTypeF[i]  %||% "untransformed")
    dist <- as.character(predDf$distribution[i] %||% "norm")

    # "t" is supported as a scale family (see the .ADM_ERR_T block above); every
    # other non-normal distribution is refused.
    is_t <- dist %in% .ADM_ERR_T
    # Count endpoints: closed-form moments, but the mean is the distribution's
    # argument (see the .ADM_RESID_POIS block).
    cnt_form <- switch(dist, pois = , dpois = .ADM_RESID_POIS,
                             binom = , dbinom = .ADM_RESID_BINOM,
                             nbinomMu = , dnbinomMu = .ADM_RESID_NBINOM, NA_integer_)
    if (dist %in% c("ordinal", "dordinal")) {
      op <- .admOrdinalSpec(ui, as.character(predDf$var[i]))
      if (is.null(op))
        .admStopErrModel(
          ep, "an 'ordinal' endpoint whose probabilities are not plain model variables",
          "admixr2 reads the K-1 category probabilities from the c(...) arguments, so each must be a model variable it can solve",
          fix = "Assign them first, e.g. `p1 <- ...; p2 <- ...; y ~ c(p1, p2)`.")
      # An ordinal observation is the VECTOR of K-1 category indicators at each
      # time, so a study contributes one stacked mean per (time, category). Each
      # probability p_k is an ordinary model variable, so the study is expressed as
      # a JOINT same-subject unit with one observation block per category -- the
      # machinery multi-compartment fits already use. The spec is registered under
      # EVERY p_k, not just the first: .admResidSpecFor() looks a row's spec up by
      # its output name, so registering only op[[1]] left every other category with
      # no spec at all (form 0, zero residual variance).
      sp_ord <- .admNewSpec(
        output = op[[1L]], form = .ADM_RESID_ORDINAL,
        ord_p = op, dv_name = as.character(predDf$var[i]))
      for (.pk in op) {
        .key <- .admOutputColName(.pk)
        sp_ord$output <- .pk
        specs[[.key]] <- sp_ord
        keys <- c(keys, .key)
      }
      next
    }
    if (dist %in% c("beta", "dbeta")) {
      bs <- .admBetaSpec(ui, as.character(predDf$var[i]))
      if (is.null(bs))
        .admStopErrModel(
          ep, "a 'beta' endpoint whose shapes are not plain model variables",
          "admixr2 derives the mean as b1/(b1+b2) from the two SOLVED shape columns, so both arguments must be model variables",
          fix = "Assign the shapes first, e.g. `b1 <- mu*phi; b2 <- (1-mu)*phi; y ~ beta(b1, b2)`.")
      key <- .admOutputColName(bs$b1)
      specs[[key]] <- .admNewSpec(
        output = bs$b1, form = .ADM_RESID_BETA,
        out_pair = c(bs$b1, bs$b2), dv_name = as.character(predDf$var[i]))
      keys <- c(keys, key)
      next
    }
    if (!is.na(cnt_form)) {
      cs <- .admCountSpec(ui, as.character(predDf$var[i]), dist)
      if (is.null(cs) || !nzchar(cs$mean_var %||% ""))
        .admStopErrModel(
          ep, paste0("a '", dist, "' endpoint whose mean cannot be resolved"),
          "admixr2 reads the mean from the distribution's argument, and this model's argument is not a plain model variable",
          fix = "Write the mean as a model variable, e.g. `cp <- ...; y ~ pois(cp)`.")
      # idx/e are computed further down for the residual-error path; the count
      # branch returns before that, so derive them here.
      .ix <- which(!is.na(cond) & cond == ep)
      if (length(.ix) == 0L) .ix <- which(is.na(cond))
      k_size <- .ix[err[.ix] %in% .ADM_ERR_NB_SIZE]
      # nbinomMu's size IS an iniDf error row (err = "nbinomMu"), so unlike
      # binom's N it has a real optimizer slot and can be estimated: it enters
      # the variance as m^2/k, which is a gradient path admixr2 does have.
      if (cnt_form == .ADM_RESID_NBINOM && length(k_size) > 0L) cs$size <- NA_real_
      if (cnt_form == .ADM_RESID_BINOM && !is.finite(cs$size))
        .admStopErrModel(
          ep, paste0("a '", dist, "' endpoint with a non-constant size/N"),
          paste0("the size (or number of trials) enters the objective ONLY through the ",
                 "variance, never through
the ODE solve, so an ESTIMATED one has no ",
                 "gradient path in admixr2"),
          fix = "Give it as a literal or fix() it, e.g. `y ~ binom(20, p)` or `sz <- fix(3)`.")
      key <- .admOutputColName(cs$mean_var)
      specs[[key]] <- .admNewSpec(
        output = cs$mean_var, form = cnt_form,
        csize = cs$size, dv_name = as.character(predDf$var[i]),
        k_size = if (length(k_size) > 0L) k_size[1L] else NA_integer_)
      keys <- c(keys, key)
      next
    }
    if (!is_t && !dist %in% c("norm", "dnorm"))
      .admStopErrModel(ep, paste0("a '", dist, "' error distribution"),
                       .ADM_ERR_WHY[[dist]] %||%
                         paste0("'", dist, "' is not a normal distribution"))
    # Transform-both-sides endpoints (logit/probit/boxCox/yeoJohnson) are handled
    # by quadrature -- see .ADM_RESID_TBS. `lnorm` keeps its exact closed form.
    yj_code <- if (tr %in% names(.ADM_TBS_YJ)) unname(.ADM_TBS_YJ[[tr]]) else NA_integer_
    is_tbs  <- !is.na(yj_code) && !tr %in% c("untransformed", "lnorm")
    if (!is_tbs && !tr %in% c("untransformed", "lnorm"))
      .admStopErrModel(ep, paste0("the '", tr, "' transform"),
                       .ADM_ERR_WHY[[tr]] %||%
                         paste0("residuals under '", tr,
                                "' are normal on the transformed scale, not the natural one"))
    if (identical(etf, "f"))
      .admStopErrModel(
        ep, "propF() / powF()",
        paste0("propF()/powF() scale the residual by a user-supplied model variable,\n",
               "which is an individual-level quantity admixr2 cannot recover from the ",
               "aggregate mean"),
        fix = "Use prop() or pow(), which scale by the prediction itself.")
    # dv(): the residual scales by the OBSERVED DV rather than the prediction. That
    # is an individual-level quantity -- an aggregate study contributes a mean and a
    # covariance, and the individual observations that dv() multiplies are exactly
    # what admixr2 never sees. Note rxode2's SIMULATION ignores dv() (it emits the
    # same rx_r_ as without it), so this was silently fitting the prediction-scaled
    # model instead: a converged fit of a model the user did not write. Same
    # reasoning as propF()/powF() above.
    if (isTRUE(predDf$dv[i]))
      .admStopErrModel(
        ep, "dv()",
        paste0("dv() scales the residual by the OBSERVED DV, an individual-level ",
               "quantity that an
aggregate mean and covariance cannot recover"),
        fix = "Use prop()/pow(), which scale by the prediction.")
    if (isTRUE(predDf$variance[i]))
      .admStopErrModel(
        ep, "variance-parameterised residual error",
        "admixr2 parameterises residual error in SD units, as nlmixr2 does by default",
        fix = "Drop `variance = TRUE` and give the error parameter(s) as standard deviations.")

    # Rows of iniDf belonging to this endpoint.
    idx    <- which(!is.na(cond) & cond == ep)
    if (length(idx) == 0L) idx <- which(is.na(cond))   # single-endpoint / no condition
    e      <- err[idx]
    k_add  <- idx[e %in% c(.ADM_ERR_ADD, .ADM_ERR_LNORM, .ADM_ERR_TBS_SD)]
    k_lam  <- idx[e %in% .ADM_ERR_TBS_LAM]
    k_prop <- idx[e %in% c(.ADM_ERR_PROP, .ADM_ERR_POW)]
    k_pow  <- idx[e %in% .ADM_ERR_POW_EXP]
    k_tdf  <- idx[e %in% .ADM_ERR_T]
    k_ar   <- idx[e %in% .ADM_ERR_AR]

    # ar(): a FIXED rho never reaches sigma_rows, same trap as a fixed nu.
    ar_fixed <- NA_real_
    if (length(k_ar) == 0L) {
      .fa <- full_ini[!is.na(full_ini$err) & full_ini$err %in% .ADM_ERR_AR &
                        full_ini$fix &
                        (is.na(full_ini$condition) | full_ini$condition == ep), , drop = FALSE]
      if (nrow(.fa) > 0L) ar_fixed <- as.numeric(.fa$est[1L])
      if (is.finite(ar_fixed) && (ar_fixed <= 0 || ar_fixed >= 1))
        .admStopErrModel(
          ep, sprintf("ar() fixed at rho = %g", ar_fixed),
          "an AR(1) residual correlation must lie strictly inside (0, 1)",
          fix = "Fix rho in (0, 1), or drop ar() for independent residuals.")
    }
    # ar(): admixr2 implements the STATIONARY continuous-time AR(1), which is what
    # rxode2's ESTIMATION lines mean. rxode2 has two separate ar() emitters and they
    # do not agree with each other:
    #
    #   .rxArSimLines   (simulation)  rx.arRes <- phi*lag0(rx.arRes,1) +
    #                                             sqrt(rx_r_*(1-phi^2))*innovation
    #   .rxArEstNormLines (estimation) rx_pred_ ~ rx_pred_ + phi*rx_arEp
    #                                  rx_r_    ~ variance * (1 - phi^2)
    #
    # The estimation form is the prediction-error decomposition, conditional on the
    # PREVIOUS OBSERVED residual; its implied MARGINAL variance is exactly `variance`,
    # i.e. the stationary process -- so the covariance admixr2 scores,
    # sqrt(v_i v_j) rho^|t_i - t_j| with diag = v, is precisely nlmixr2's own
    # estimation semantics. (admixr2 cannot use the conditional form: it needs each
    # subject's previous residual, and aggregate data has no individual observations.)
    #
    # rxode2's SIMULATION, by contrast, is not stationary when a dose record precedes
    # the first observation. Measured with plain rxode2 (no admixr2), omega ~ 0,
    # a = 0.5, rho = 0.6, N = 1.5e5, MC se 0.0032 -- V/a^2 at t = 0.5, 1, 2, 4, 8:
    #
    #   observations only, no dose record : 0.998 1.004 1.001 0.999 1.004   correct
    #   dose at t = 0                     : 1.975 1.591 1.214 1.027 1.004
    #   ZERO-amount dose (amt = 0)        : 1.975 1.591 1.214 1.027 1.004
    #   dose + an observation AT t = 0    : 1.001 0.998 1.002 0.995 0.996   correct
    #   plain add(a) with a dose (control): 0.998 1.008 0.993 1.000 1.005   correct
    #
    # A dose of amt = 0 doubles the residual variance, and plain add() is unaffected,
    # so this is record-driven and specific to ar(). It is an upstream defect, not a
    # modelling choice: nlmixr2's OWN focei cannot recover rho from rxode2's OWN
    # simulation either -- fitting individual-level simulated data (N = 400, 8
    # observations) returns rho = 0.4617 against a truth of 0.60, improving to 0.5398
    # when an observation is placed at the dose time. admixr2 is not involved in that
    # test at all.
    #
    # Consequence to be aware of: simulating an ar() fit through rxode2 will not
    # reproduce the covariance admixr2 fitted until that is fixed (49% of max|V| on a
    # typical design), while every other error model here round-trips to within
    # Monte-Carlo noise. Validated against a genuinely stationary AR(1) truth,
    # admixr2 recovers rho to +1.0%. Mimicking the simulator would put admixr2 at
    # odds with nlmixr2's estimator and would break when rxode2 is fixed.

    if ((length(k_ar) > 0L || is.finite(ar_fixed)) && identical(tr, "lnorm"))
      .admStopErrModel(
        ep, "lnorm() combined with ar()",
        paste0("ar() correlates residuals on the NATURAL scale while lnorm() is a ",
               "multiplicative
log-scale residual; composing them has no single ",
               "stationary covariance admixr2 can score"),
        fix = "Use add(a) + ar(rho), or lnorm(a) alone.")

    # A FIXED nu never reaches sigma_rows (.admParseIniDf drops fixed error rows),
    # so it has no optimizer slot to index -- but it still has to scale the
    # variance. Without this the multiplier would be silently DROPPED and the
    # endpoint fitted as a plain normal: a converged fit of a model the user did
    # not write, exactly the failure mode pow() had. And since the identifiability
    # warning below tells users to fix() nu, this is the RECOMMENDED path, not an
    # edge case. Read the value straight off the full iniDf instead.
    tdf_fixed <- NA_real_
    if (is_t && length(k_tdf) == 0L) {
      .fx <- full_ini[!is.na(full_ini$err) & full_ini$err %in% .ADM_ERR_T &
                        full_ini$fix &
                        (is.na(full_ini$condition) | full_ini$condition == ep), , drop = FALSE]
      if (nrow(.fx) > 0L) tdf_fixed <- as.numeric(.fx$est[1L])
      if (is.finite(tdf_fixed) && tdf_fixed <= 2)
        .admStopErrModel(
          ep, sprintf("t() fixed at nu = %g", tdf_fixed),
          "a Student-t has no finite variance at nu <= 2, so there is no aggregate variance to match",
          fix = "Fix nu above 2, or drop t() and use a normal residual.")
    }

    # A t() endpoint must actually have a scale to multiply: `cp ~ t(nu)` alone
    # leaves errType "none" and no add/prop row, so there is no residual magnitude
    # to estimate and the multiplier would scale nothing. Count FIXED scale rows
    # too -- add(a) with a fixed `a` is still a scale.
    .has_scale <- length(k_add) > 0L || length(k_prop) > 0L ||
      any(!is.na(full_ini$err) &
            full_ini$err %in% c(.ADM_ERR_ADD, .ADM_ERR_PROP, .ADM_ERR_POW) &
            (is.na(full_ini$condition) | full_ini$condition == ep))
    if (is_t && !.has_scale)
      .admStopErrModel(
        ep, "t(nu) with no scale parameter",
        paste0("nlmixr2 writes Student-t residuals as a SCALE FAMILY -- residual = ",
               "scale * T_nu -- so t(nu)\nneeds an add()/prop()/pow() term to supply ",
               "that scale. On its own it fixes the scale at 1,\nwhich admixr2 cannot ",
               "tell apart from a model with no residual error at all"),
        fix = "Write cp ~ add(a) + t(nu) (or prop(b) + t(nu)).")

    # nu is NOT IDENTIFIABLE from aggregate data, and this is structural, not a
    # small-sample problem: nu enters the aggregate moments ONLY through the
    # multiplier m = nu/(nu-2), and the scale enters only as a^2*m (and b^2*m).
    # So the data see one number per coefficient and the pair (a, nu) traces a
    # flat ridge -- an optimizer will return whatever the starting value drifts
    # to. Measured: truth (a=0.8, nu=5) came back as (0.915, 8.27) with the
    # PRODUCT a^2*nu/(nu-2) accurate to 3.6%. Anyone reading the reported nu as an
    # estimate of tail weight would be reading a starting-value artefact.
    #
    # Deliberately emitted AFTER the structural refusals above: a model with no
    # scale at all should get the scale error, not advice about being aliased with
    # a scale that does not exist.
    if (is_t && length(k_tdf) > 0L)
      warning(
        "Student-t degrees of freedom (", sigma_names[k_tdf[1L]], ") cannot be ",
        "estimated from aggregate data.\n",
        "  admixr2 matches the residual MEAN and VARIANCE, and nu enters only ",
        "through the variance\n  multiplier nu/(nu-2) -- so nu and the scale (",
        paste(sigma_names[c(k_add, k_prop)], collapse = ", "),
        ") are aliased: only a^2*nu/(nu-2)\n  is identified. The reported nu will ",
        "reflect its starting value, not the data.\n",
        "  Fix nu instead, e.g. ", sigma_names[k_tdf[1L]], " <- fix(5).",
        call. = FALSE)

    if (is_t && identical(tr, "lnorm"))
      .admStopErrModel(
        ep, "lnorm() combined with t()",
        paste0("lnorm() moment-matches a LOGNORMAL residual; t() is a scale family on ",
               "the NATURAL\nscale. Composing them describes no single distribution with ",
               "a defined aggregate mean\nand variance"),
        fix = "Use lnorm(a) alone, or add(a) + t(nu) on the natural scale.")

    # Transform-both-sides: resolve lambda (estimated row, fixed row, or the
    # transform's fixed default) and the (lo, hi) bounds. A FIXED lambda never
    # reaches sigma_rows, same trap as fixed nu / fixed rho.
    lam_fixed <- NA_real_
    if (is_tbs && length(k_lam) == 0L) {
      .fl <- full_ini[!is.na(full_ini$err) & full_ini$err %in% .ADM_ERR_TBS_LAM &
                        full_ini$fix &
                        (is.na(full_ini$condition) | full_ini$condition == ep), , drop = FALSE]
      lam_fixed <- if (nrow(.fl) > 0L) as.numeric(.fl$est[1L]) else 1.0
    }
    # A FIXED add/prop/pow parameter never reaches sigma_rows either -- exactly the
    # trap tdf_fixed / ar_fixed / lam_fixed above exist to close, but the SCALE
    # parameters had no equivalent, so they were silently dropped from the residual:
    #   add(a),  a <- fix(0.7)          -> fitted with NO residual variance at all
    #   add(a) + prop(b), b <- fix(0.2) -> the proportional term vanished
    #   pow(b, c), c <- fix(0.75)       -> reverted to prop(), i.e. c = 1
    # The last is precisely the historical pow() failure this file documents:
    # a converged fit, plausible numbers, a model the user did not write. And
    # fix()ing a residual parameter is routine -- it is what this file's own t()
    # advice tells users to do for nu.
    # est is the SD/coefficient, so the stored VARIANCE is est^2; a pow exponent
    # is stored as itself.
    .fixval <- function(codes, sq) {
      .f <- full_ini[!is.na(full_ini$err) & full_ini$err %in% codes & full_ini$fix &
                       (is.na(full_ini$condition) | full_ini$condition == ep), , drop = FALSE]
      if (nrow(.f) == 0L) return(NA_real_)
      v <- as.numeric(.f$est[1L])
      if (sq) v * v else v
    }
    add_fixed  <- if (length(k_add)  == 0L)
      .fixval(c(.ADM_ERR_ADD, .ADM_ERR_LNORM, .ADM_ERR_TBS_SD), TRUE)  else NA_real_
    prop_fixed <- if (length(k_prop) == 0L)
      .fixval(c(.ADM_ERR_PROP, .ADM_ERR_POW), TRUE)                    else NA_real_
    pow_fixed  <- if (length(k_pow)  == 0L)
      .fixval(.ADM_ERR_POW_EXP, FALSE)                                 else NA_real_

    # ar() is only stationary when the residual variance is CONSTANT in time.
    # rxode2 emits
    #   rx.arRes = phi*lag0(rx.arRes,1) + sqrt(rx_r_*(1 - phi^2))*rxerr,
    # so the marginal variance obeys V_i = phi_i^2 V_{i-1} + r_i (1 - phi_i^2).
    # With r constant that telescopes to V_i = r and admixr2's
    # "diagonal unchanged, off-diagonal = sqrt(ev_i ev_j) rho^|dt|" is exact --
    # which is why add(a) + ar(rho) matches the simulator to within MC noise.
    # With a prediction-dependent r (prop/pow/combined) it does NOT: the process
    # is non-stationary, V_i depends on the whole history, and it is additionally
    # seeded from the first record -- usually the DOSE row, where f = 0 so r = 0.
    # Measured against rxode2's own simulation, admixr2's variance came out
    # 2.4-12x too high and its correlations 3-11x too high, silently.
    #
    # Refused rather than approximated: reproducing it would mean encoding
    # rxode2's recursion AND its dose-row initialisation into an aggregate
    # likelihood, which is an implementation detail rather than a stated model.
    if ((length(k_ar) > 0L || is.finite(ar_fixed)) &&
        (length(k_prop) > 0L || is.finite(prop_fixed) ||
         length(k_pow)  > 0L || is.finite(pow_fixed)))
      .admStopErrModel(
        ep, "ar() combined with a proportional or power residual",
        paste0("rxode2 scales the AR(1) innovation by sqrt(rx_r_*(1 - phi^2)), which ",
               "leaves the marginal\nvariance equal to rx_r_ only when rx_r_ is ",
               "CONSTANT. With prop()/pow() the residual variance\nvaries with the ",
               "prediction, so the process is non-stationary and its covariance ",
               "depends on\nthe whole observation history -- there is no ",
               "sqrt(v_i v_j)*rho^|dt| form for admixr2 to score"),
        fix = "Use add(a) + ar(rho), or drop ar() and keep prop()/pow().")

    tr_lo <- suppressWarnings(as.numeric(predDf$trLow[i] %||% NA))
    tr_hi <- suppressWarnings(as.numeric(predDf$trHi[i]  %||% NA))
    if (is_tbs) {
      if (length(k_add) == 0L && !is.finite(add_fixed))
        .admStopErrModel(
          ep, paste0(tr, " with no residual SD"),
          paste0("a transform-both-sides residual is normal on the ", tr, " scale with some ",
                 "SD;
without that parameter there is no residual to integrate"),
          fix = "Give the endpoint an add() term (or use logitNorm(sd, lo, hi)).")
      if (yj_code %in% c(4L, 6L) && (!is.finite(tr_lo) || !is.finite(tr_hi) || tr_hi <= tr_lo))
        .admStopErrModel(
          ep, paste0(tr, " with unusable bounds"),
          "a logit/probit endpoint needs finite bounds lo < hi to transform through",
          fix = "Write cp ~ logitNorm(sd, lo, hi) with numeric lo < hi.")
      if (!is.na(k_ar[1L]) || is.finite(ar_fixed))
        .admStopErrModel(
          ep, paste0(tr, " combined with ar()"),
          "ar() correlates residuals on the NATURAL scale; a transformed residual is normal on the TRANSFORMED scale, so the two have no common stationary covariance",
          fix = "Use add(a) + ar(rho) on the natural scale.")
    }

    # lnorm: the "add" parameter is the SD on the log scale; it cannot be
    # combined with a proportional or power term.
    if (identical(tr, "lnorm")) {
      if (length(k_prop) > 0L)
        .admStopErrModel(
          ep, "lnorm() combined with a proportional or power term",
          paste0("lnorm()'s parameter is the SD on the LOG scale, which already makes the\n",
                 "residual proportional to the prediction; adding prop()/pow() on top has no\n",
                 "single well-defined aggregate variance"),
          fix = "Use lnorm(a) alone, or add(a) + prop(b) on the natural scale.")
      form <- .ADM_RESID_LNORM
    } else if (is_tbs) {
      form <- .ADM_RESID_TBS
    } else {
      ap <- as.character(predDf$addProp[i] %||% "default")
      if (identical(ap, "default")) ap <- add_prop_default
      form <- if (identical(ap, "combined1")) .ADM_RESID_COMBINED1 else .ADM_RESID_COMBINED2
    }

    key <- .admOutputColName(ep)
    specs[[key]] <- .admNewSpec(
      output = ep, form = form,
      k_add  = if (length(k_add)  > 0L) k_add[1L]  else NA_integer_,
      k_prop = if (length(k_prop) > 0L) k_prop[1L] else NA_integer_,
      k_pow  = if (length(k_pow)  > 0L) k_pow[1L]  else NA_integer_,
      k_tdf  = if (length(k_tdf)  > 0L) k_tdf[1L]  else NA_integer_,
      tdf_fixed = tdf_fixed,
      k_ar   = if (length(k_ar)   > 0L) k_ar[1L]   else NA_integer_,
      ar_fixed = ar_fixed,
      k_lam  = if (length(k_lam)  > 0L) k_lam[1L]  else NA_integer_,
      lam_fixed = lam_fixed, yj = yj_code, tr_lo = tr_lo, tr_hi = tr_hi,
      add_fixed = add_fixed, prop_fixed = prop_fixed, pow_fixed = pow_fixed,
      # A TBS endpoint's residual is normal on the TRANSFORMED scale with variance
      # rx_r_ -- which rxode2 builds from the SAME add/prop/pow/combined machinery
      # as any other endpoint (`rx_r_ ~ (a)^2 + (rx_pred_f_)^2*(b)^2`). form is
      # .ADM_RESID_TBS, which loses the combined1/combined2 distinction and the
      # errTypeF, so carry both: without them a prop() term on a transformed
      # endpoint contributed NOTHING to the objective and got an exactly-zero
      # gradient. errTypeF "transformed" (propT) scales by rx_pred_ = h(f) rather
      # than rx_pred_f_ = f.
      tbs_c1  = identical(if (identical(as.character(predDf$addProp[i] %||% "default"),
                                        "default")) add_prop_default
                          else as.character(predDf$addProp[i]), "combined1"),
      tbs_ftr = identical(etf, "transformed"))
    keys <- c(keys, key)
  }
  if (length(specs) == 0L) return(NULL)
  specs
}

# Optimizer-scale role of each residual parameter.
#   "var"     log-variance: natural value is a VARIANCE, exp(p).  (a^2, b^2)
#   "pow_exp" a power exponent, not a variance: natural value is p itself.
#   "t_df"    Student-t degrees of freedom: natural value is 2 + exp(p).
#   "tbs_lam" a boxCox/yeoJohnson lambda: identity, unconstrained -- it may be
#             zero (which IS the log transform) or negative.
#   "nb_size" a negative-binomial size: natural value is exp(p) (strictly > 0).
#   "ar_cor"  ar() correlation: natural value is expit(p), so it stays in (0,1)
#             -- rxode2 itself validates rho in [0,1), and rho >= 1 is a
#             non-stationary process with no stationary covariance to match.
# The distinction matters because a pow() exponent must not be squared, bounded
# below at zero, or reported as an SD. Absent (legacy/hand-built pinfo) => all "var".
#
# Why t_df is estimated as log(nu - 2) rather than unconstrained: the variance
# multiplier nu/(nu-2) has a pole at nu = 2 and is NEGATIVE below it, so an
# optimizer step to nu <= 2 would hand the MVN kernel a negative variance. The
# shifted log keeps nu > 2 for every real p, so the constraint cannot be violated
# by a line search instead of merely being checked at the start.
.admSigmaRole <- function(pinfo) {
  r <- pinfo$sigma_role
  if (is.null(r)) rep("var", length(pinfo$sigma_names)) else r
}

# Optimizer vector -> natural-scale residual parameters.
.admSigmaNat <- function(p_sigma, pinfo) {
  role <- .admSigmaRole(pinfo)
  out  <- p_sigma
  iv   <- role == "var"
  it   <- role == "t_df"
  ia <- role == "ar_cor"
  ib <- role == "nb_size"
  out[iv] <- exp(p_sigma[iv])                    # log-variance -> variance
  out[it] <- 2 + exp(p_sigma[it])                # log(nu - 2)  -> nu (> 2)
  out[ia] <- 1 / (1 + exp(-p_sigma[ia]))         # logit        -> rho in (0,1)
  out[ib] <- exp(p_sigma[ib])                    # log(size)    -> size (> 0)
  ie <- !iv & !it & !ia & !ib
  out[ie] <- p_sigma[ie]               # exponent: identity
  setNames(out, pinfo$sigma_names)
}

# The value .admFullTheta() REPORTS for ONE residual parameter, as a function of
# its optimizer value.
#
# plot.R traces parameters one at a time and so needs the map per name rather
# than for the whole sigma vector. It is built out of .admSigmaNat() instead of
# re-deriving the roles, because a trace plotted on a different scale from the
# one print(fit) reports is a silent disagreement: an ar() correlation of 0.6
# (optimizer value 0.405) came out as 1.22 under the generic sigma rule -- past
# the top of its own support -- and a Box-Cox lambda of 0.5 as 1.28.
.admSigmaReportFn <- function(pinfo, nm) {
  .k <- match(nm, pinfo$sigma_names)
  if (is.na(.k)) return(function(v) exp(v / 2))
  # unname: pinfo$sigma_role is a NAMED character vector, and a named "var" is not
  # identical() to "var" -- which silently sent every residual SD back as a
  # VARIANCE (0.2 plotted as 0.04).
  .p1 <- list(sigma_names = nm, sigma_role = unname(.admSigmaRole(pinfo)[[.k]]))
  .sd <- identical(.p1$sigma_role, "var")
  function(v) vapply(v, function(x) {
    nat <- unname(.admSigmaNat(x, .p1))
    if (.sd) sqrt(nat) else nat            # a "var" role is reported as an SD
  }, double(1))
}

# Delta-method factor d(REPORTED value)/d(optimizer p), one per residual parameter.
#
# The post-fit Hessian is taken w.r.t. the optimizer's parameterisation, but
# `fit$cov` sits beside `Estimate` in nlmixr2est's parFixed table, so the two must
# be on the same scale or `Estimate +- 1.96*SE` is wrong by this factor. It lives
# here, next to .admSigmaNat() and .admSigmaRole(), because it is the derivative of
# exactly that map -- the three CalcCov functions used to each carry their own
# copy of the switch(), which is three places to forget when a role is added.
#
# `.admFullTheta()` reports a "var" role as an SD, so the factor is d(sd)/dp with
# p = log(sd^2): sd = exp(p/2), d(sd)/dp = sd/2. The identity roles (pow_exp,
# tbs_lam) get 1 because they are reported exactly as held.
.admSigmaReportJac <- function(p_sigma, pinfo) {
  if (length(p_sigma) == 0L) return(numeric(0))
  role <- .admSigmaRole(pinfo)
  nat  <- .admSigmaNat(p_sigma, pinfo)
  vapply(seq_along(p_sigma), function(k)
    switch(role[k],
           var     = sqrt(nat[[k]]) / 2,             # reported sd   = exp(p/2)
           t_df    = nat[[k]] - 2,                   # reported nu   = 2 + exp(p)
           ar_cor  = nat[[k]] * (1 - nat[[k]]),      # reported rho  = expit(p)
           nb_size = nat[[k]],                       # reported size = exp(p)
           1), double(1))                            # pow_exp / tbs_lam: identity
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
# Residual row array for a unit, with the SOLVED beta precision folded in.
#
# .admResidRows() knows the error model from pinfo but NOT phi = b1 + b2, which for
# a beta endpoint is derived from the SOLVED shapes and rides back as an attribute
# on the simulate matrix -- it cannot come from pinfo. Threading it was a two-line
# `.ph <- attr(cp, "phi"); if (!is.null(.ph)) arr$phi <- .ph` copied at every
# estimator moment path, datagen and plot; forgetting it left arr$phi NA and every
# predicted-covariance entry NaN, silently (a bug the call-site comments record more
# than once). Folding the assignment into the build makes "residual rows for a unit,
# with its solved phi" one call. The phi VALUE stays the caller's -- it legitimately
# varies (attr(cp_mat), a per-config .phi_all row, or a diagnostic argument) -- and
# is passed as `phi`; NULL leaves the NA field untouched (every non-beta endpoint).
.admUnitResidRows <- function(pinfo, row_output, sigma_nat, n_t, phi = NULL) {
  arr <- .admResidRows(pinfo, row_output, sigma_nat, n_t)
  if (!is.null(phi)) arr$phi <- phi
  arr
}

.admResidRows <- function(pinfo, row_output, sigma_nat, n_t) {
  form <- integer(n_t)
  a2   <- numeric(n_t)   # additive VARIANCE  (a^2), or the lnorm log-variance
  b2   <- numeric(n_t)   # prop/pow VARIANCE  (b^2)
  cc   <- rep(1.0, n_t)  # power exponent
  vmul <- rep(1.0, n_t)  # Student-t variance multiplier nu/(nu-2); 1 when normal
  rho  <- rep(NA_real_, n_t)   # ar() residual correlation; NA when independent
  lam  <- rep(1.0, n_t); yj <- rep(2L, n_t)       # TBS transform + parameter
  csz  <- rep(NA_real_, n_t)                     # binom N / nbinom size
  # Residual-quadrature node count, carried on `arr` because .admResidApply()
  # takes `arr` but not `pinfo`. Set from the control's `resid_nodes` by the
  # drivers; the default keeps every existing fit bit-identical.
  nodes   <- pinfo$resid_nodes %||% .ADM_TBS_NODES
  tbs_c1  <- rep(FALSE, n_t)   # TBS endpoint: combined1 rather than combined2
  tbs_ftr <- rep(FALSE, n_t)   # TBS endpoint: propT/powT (scale by the TRANSFORMED pred)
  phi  <- rep(NA_real_, n_t)                     # beta precision (SOLVED, set by caller)
  ord_grp <- rep(NA_integer_, n_t)               # ordinal: which time group a row is in
  tlo  <- rep(0.0, n_t); thi <- rep(1.0, n_t)
  k_add <- k_prop <- k_pow <- k_tdf <- k_ar <- k_lam <- k_size <-
    rep(NA_integer_, n_t)

  # Again: a count endpoint has no sigma rows but DOES have a spec, so only bail
  # when there is genuinely nothing to apply.
  if (length(pinfo$sigma_names) == 0L && length(.admResidSpecs(pinfo)) == 0L)
    return(list(form = form, a2 = a2, b2 = b2, cc = cc, vmul = vmul, rho = rho,
                lam = lam, yj = yj, tlo = tlo, thi = thi, csz = csz, phi = phi,
                ord_grp = ord_grp, tbs_c1 = tbs_c1, tbs_ftr = tbs_ftr,
                nodes = nodes,
                k_add = k_add, k_prop = k_prop, k_pow = k_pow, k_tdf = k_tdf,
                k_ar = k_ar, k_lam = k_lam, k_size = k_size))

  if (is.null(row_output)) row_output <- rep(NA_character_, n_t)
  if (length(row_output) == 1L) row_output <- rep(row_output, n_t)

  for (ov in unique(row_output)) {
    rows <- which(row_output == ov | (is.na(row_output) & is.na(ov)))
    sp   <- .admResidSpecFor(pinfo, ov)
    if (is.null(sp)) next
    form[rows] <- sp$form
    # Estimated parameters index sigma_nat; fix()ed ones have NO optimizer slot and
    # come from the spec constants instead (add_fixed/prop_fixed/pow_fixed), exactly
    # as a fixed nu/rho/lambda does below. Only the estimated ones get a k_* gradient
    # slot -- a fixed parameter must contribute to the variance but not to the gradient.
    if (!is.na(sp$k_add))  { a2[rows] <- sigma_nat[[sp$k_add]];  k_add[rows]  <- sp$k_add }
    else if (is.finite(sp$add_fixed  %||% NA_real_)) a2[rows] <- sp$add_fixed
    if (!is.na(sp$k_prop)) { b2[rows] <- sigma_nat[[sp$k_prop]]; k_prop[rows] <- sp$k_prop }
    else if (is.finite(sp$prop_fixed %||% NA_real_)) b2[rows] <- sp$prop_fixed
    if (!is.na(sp$k_pow))  { cc[rows] <- sigma_nat[[sp$k_pow]];  k_pow[rows]  <- sp$k_pow }
    else if (is.finite(sp$pow_fixed  %||% NA_real_)) cc[rows] <- sp$pow_fixed
    # Student-t: fold nu/(nu-2) into the variance coefficients. Exact for both
    # combined forms (see the .ADM_ERR_T block), so no downstream consumer -- R or
    # C++ -- has to know that this endpoint is not normal. nu comes either from an
    # estimated slot (k_tdf) or, when the user fix()ed it, from the spec constant.
    .nu <- if (!is.null(sp$k_tdf) && !is.na(sp$k_tdf)) sigma_nat[[sp$k_tdf]]
           else if (!is.null(sp$tdf_fixed) && is.finite(sp$tdf_fixed)) sp$tdf_fixed
           else NA_real_
    if (!is.na(.nu)) {
      m <- .admTMult(.nu)
      vmul[rows] <- m
      a2[rows]   <- a2[rows] * m
      b2[rows]   <- b2[rows] * m
      # only an ESTIMATED nu gets a gradient slot; a fixed one has none
      if (!is.null(sp$k_tdf) && !is.na(sp$k_tdf)) k_tdf[rows] <- sp$k_tdf
    }
    if (!is.null(sp$csize)) csz[rows] <- sp$csize
    if (!is.null(sp$k_size) && !is.na(sp$k_size)) {
      csz[rows]    <- sigma_nat[[sp$k_size]]
      k_size[rows] <- sp$k_size
    }
    # TBS: transform code, bounds and lambda (estimated slot or fixed constant)
    if (identical(sp$form, .ADM_RESID_TBS)) {
      yj[rows]  <- sp$yj
      tlo[rows] <- if (is.finite(sp$tr_lo)) sp$tr_lo else 0
      thi[rows] <- if (is.finite(sp$tr_hi)) sp$tr_hi else 1
      lam[rows] <- if (!is.null(sp$k_lam) && !is.na(sp$k_lam)) sigma_nat[[sp$k_lam]]
                   else if (!is.null(sp$lam_fixed) && is.finite(sp$lam_fixed)) sp$lam_fixed
                   else 1
      if (!is.null(sp$k_lam) && !is.na(sp$k_lam)) k_lam[rows] <- sp$k_lam
      tbs_c1[rows]  <- isTRUE(sp$tbs_c1)
      tbs_ftr[rows] <- isTRUE(sp$tbs_ftr)
    }
    .r <- if (!is.null(sp$k_ar) && !is.na(sp$k_ar)) sigma_nat[[sp$k_ar]]
          else if (!is.null(sp$ar_fixed) && is.finite(sp$ar_fixed)) sp$ar_fixed
          else NA_real_
    if (!is.na(.r)) {
      rho[rows] <- .r
      if (!is.null(sp$k_ar) && !is.na(sp$k_ar)) k_ar[rows] <- sp$k_ar
    }
  }
  # phi and ord_grp MUST be here too. They were present only in the early-return
  # branch above (and duplicated there), so on the normal path arr$ord_grp was NULL
  # and .admResidApply's `if (any(or_) && !is.null(arr$ord_grp))` never fired --
  # an ordinal endpoint silently lost its -p_j p_k cross-category covariance, and
  # beta's phi had nowhere to be patched in.
  list(form = form, a2 = a2, b2 = b2, cc = cc, vmul = vmul, rho = rho,
       lam = lam, yj = yj, tlo = tlo, thi = thi, csz = csz, phi = phi,
       ord_grp = ord_grp, tbs_c1 = tbs_c1, tbs_ftr = tbs_ftr,
       nodes = nodes,
       k_add = k_add, k_prop = k_prop, k_pow = k_pow, k_tdf = k_tdf, k_ar = k_ar,
       k_lam = k_lam, k_size = k_size)
}

# E[f^k] over the random-effects distribution, from the first two moments of f.
#
# Exact for k = 1 and k = 2 (the k(k-1)/2 factor vanishes at k = 1, and at k = 2
# gives E[f^2] = mu^2 + var exactly); a second-order delta expansion otherwise,
# which is only reached by pow() with c != 1 and c != 0.5. A caller holding the
# actual sample/node ensemble should pass exact moments via `mom` instead.
.admMomF <- function(mu, var_f, k) {
  k <- as.numeric(k)
  out <- mu^k
  nz  <- k != 1
  if (any(nz)) {
    kk <- k[nz]
    mm <- mu[nz]
    # mu^(k-2) at or NEAR mu == 0. TWO different situations, needing opposite answers:
    #
    #  k == 2 exactly (add/prop/combined1/combined2): the coefficient k(k-1)/2 * ...
    #    is not zero but mu^0 = 1, so nothing diverges -- E[f^2] = mu^2 + var, exact.
    #  k <  2 (pow with c < 1): mu^(k-2) is a genuine POLE. The second-order delta
    #    expansion simply does not exist at mu = 0.
    #
    # rxode2's safePow (substitute DBL_EPSILON for a zero base) is the right rule
    # for a SOLVE evaluating x^y, but the WRONG one inside a Taylor expansion:
    # eps^(k-2) is astronomically large and gets multiplied by a non-zero
    # coefficient. Cap the correction against the leading term instead.
    kk2  <- kk - 2
    lead <- out[nz]
    corr <- (kk * (kk - 1) / 2) * mm^kk2 * var_f[nz]
    cap  <- kk2 < 0 & is.finite(lead)
    if (any(cap))
      corr[cap] <- sign(corr[cap]) * pmin(abs(corr[cap]), abs(lead[cap]))
    corr[!is.finite(corr)] <- 0
    out[nz] <- lead + corr
  }
  out
}

# E[f^k] together with its partials w.r.t. mu, var_f and k. Same expansion as
# .admMomF (so the gradient differentiates exactly what the NLL evaluates, term
# for term). `dk` is only consumed by the pow() exponent, whose optimizer scale is
# the exponent itself; note k = 2c there, so the caller multiplies dk by 2.
.admMomFd <- function(mu, v0, k) {
  g   <- k * (k - 1) / 2
  mk  <- mu^k
  # rxode2's safeLog (_safe_log_, rxode2_model_shared.h): log(x) for x > 0, else
  # log(DBL_EPSILON). Was ifelse(mu > 0, log(mu), 0) -- a third convention.
  lg  <- ifelse(mu > 0, log(mu), log(.Machine$double.eps))
  # 0^(negative exponent). A structural prediction of exactly 0 is not exotic -- it
  # is what a dose-at-t=0 model returns for an observation at t = 0, which real
  # aggregate datasets have constantly (Theophylline does). Every synthetic
  # gradient test here starts at t > 0, which is why this reached the real-data fit
  # before it was caught: one NaN row poisoned dv_df and, through it, every
  # structural and omega gradient.
  #
  # The rule is rxode2's, NOT one invented here. rxode2 solves `x^y` through
  # Rx_pow_() with safePow = TRUE (its rxSolve default), which substitutes
  # x = DBL_EPSILON when x == 0 and y < 0 (rxode2_model_shared.h; documented in
  # rxode2 NEWS under #775 alongside safeZero/safeLog). Mirroring it means the
  # moment expansion degrades exactly the way the ODE solve of the SAME model
  # already does, instead of two different conventions inside one fit.
  #
  # For the exponents admixr2 actually uses this is not even an approximation:
  # c = 1 (add/prop/combined1/combined2) gives k = 2, where the coefficient
  # g*(k-2) is EXACTLY zero, so the substituted term contributes a clean 0 and
  # E[f^2] = mu^2 + var stays exact. Only pow() with c < 1 -- where the
  # second-order expansion genuinely diverges at f = 0 -- sees the eps value.
  #
  # THE CAP IS .admMomF's, NOT A SECOND ONE. A negative exponent is a genuine pole
  # of the expansion near mu = 0, and .admMomF (and adm_mom_f in src/nll.cpp)
  # handle it by capping the correction against the LEADING term. This function
  # used to zero mu^(k-2) once it passed .ADM_MOM_CAP instead -- a different rule,
  # so for pow(b, c) with c < 1 near a zero prediction the objective used a
  # correction of size mu^(2c) while the gradient used the uncapped expansion, and
  # the optimizer was handed a direction that does not descend the function it is
  # minimising. The two must agree term for term, so the derivatives below are the
  # derivatives of the CAPPED expression, piecewise.
  #
  # Where the cap binds, corr = sign(corr)*|lead| and lead = mu^k > 0, so
  # m = lead*(1 + s) with s = sign(corr): the correction no longer depends on
  # var_f at all (dv0 = 0) and both remaining partials just scale the leading term.
  lead <- mk
  pk2  <- mu^(k - 2)
  corr <- g * pk2 * v0
  bind <- (k - 2) < 0 & is.finite(lead) & is.finite(corr) & abs(corr) > abs(lead)
  bind[is.na(bind)] <- FALSE
  s    <- sign(corr)
  # non-finite corrections contribute nothing, exactly as in .admMomF
  nf   <- !is.finite(corr)
  corr[nf] <- 0
  pk1  <- mu^(k - 1)
  pk3  <- mu^(k - 3)
  pk1[!is.finite(pk1)] <- 0
  pk3[!is.finite(pk3)] <- 0
  pk2z <- pk2; pk2z[!is.finite(pk2z)] <- 0

  m   <- ifelse(bind, lead * (1 + s), lead + corr)
  dmu <- ifelse(bind, (1 + s) * k * pk1,
                k * pk1 + g * (k - 2) * pk3 * v0)
  dv0 <- ifelse(bind, 0, g * pk2z)
  dk  <- ifelse(bind, (1 + s) * lg * mk,
                lg * mk + v0 * pk2z * ((2 * k - 1) / 2 + g * lg))
  list(m = m, dmu = dmu, dv0 = dv0, dk = dk)
}

# Compose a full structural covariance with an .admResidApply() result: the
# lnorm/TBS off-diagonal mean-scale (ms_i ms_j), the composed diagonal (ap$dv) and
# any ar() correlation matrix (ap$rmat). This three-line tail was hand-copied at
# ~11 sites -- every estimator's moment/objective path, plot.R, datagen.R and
# .admJointResidual -- so an added off-diagonal residual channel meant editing all
# of them, and missing one silently dropped that endpoint's off-diagonal predicted
# covariance on that path (the exact hazard CLAUDE.md flags).
#
# The na.rm guard is load-bearing, not cosmetic: .admTBSi() returns NaN outside a
# transform's support, which a line search inside grad_bounds can reach, and a bare
# `any(ap$ms != 1)` is then NA -- so `if (NA)` ABORTS the whole fit instead of the
# optimizer rejecting the point. With na.rm the multiply is skipped; a NaN in ap$dv
# still yields a non-finite objective, so the point is rejected either way. For a
# constant scale (ms == 1: add/prop/count) tcrossprod(ms) is all ones and the
# multiply is a no-op, so add() models stay bit-identical.
.admApplyResidTail <- function(V, ap) {
  if (any(ap$ms != 1, na.rm = TRUE)) V <- V * tcrossprod(ap$ms)
  diag(V) <- ap$dv
  if (!is.null(ap$rmat)) V <- V + ap$rmat
  V
}

# Apply the residual to a structural mean/variance -- the LAW OF TOTAL VARIANCE.
#
#   mu_pred = ms * E[f]
#   V_pred  = ms^2 * Cov_eta(f)  +  diag( E_eta[ Var(y | eta) ] )
#
# `dv` comes IN as diag(Cov_eta(f)) -- the STRUCTURAL variance, which is what all
# five estimator call sites already pass -- and goes OUT as the full V_pred
# diagonal. `ms` is returned so the caller can scale the OFF-diagonals too:
# lnorm's conditional mean is f*exp(s/2), so the whole covariance is scaled, not
# just its diagonal. A caller holding a full matrix passes the result to
# .admApplyResidTail(V_struct, ap) (which does the ms-scale, diag and ar() rmat).
#
# Why this is not the old `dv + v(mu)`: the residual variance of a prop/pow/lnorm
# model depends on f, so E_eta[Var(y|eta)] != Var(y | eta = mean). Evaluating at
# the population mean is the NONMEM "no eta-eps interaction" convention; it
# understates V_pred by b^2*Var_eta(f) for prop, and for lnorm it additionally
# drops an exp(s) factor from every off-diagonal. Both are systematic and both
# were measured against individual-level simulation (see the oracle tests).
#
# The moments are ALWAYS taken from (mu, var_f) via .admMomF -- never from the
# caller's sample ensemble, even where one exists. That is deliberate: the
# gradient chains analytically through (mu, var_f), so using sample moments in the
# NLL and the closed form in the gradient would make the two disagree and the
# optimizer chase a discontinuity. One formula, used by every estimator and by
# both C++ kernels, keeps NLL and gradient consistent by construction.
# `times` is needed only by ar(): the residual correlation depends on the OBSERVATION
# TIMES (rho^|t_i - t_j|), which no other form cares about. Callers holding a full
# covariance pass them and add the returned `rmat` to the off-diagonals; callers on
# the diagonal-only (`method = "var"`) path pass NULL, and an ar() model is refused
# there because a diagonal V carries no information about rho at all.
.admResidApply <- function(mu_struct, dv, arr, times = NULL, cov_f = NULL) {
  mu <- mu_struct
  ms <- .admResidMuScale(arr)
  ln <- arr$form == .ADM_RESID_LNORM
  nz <- arr$form == .ADM_RESID_COMBINED1 | arr$form == .ADM_RESID_COMBINED2
  vf <- dv                      # Var_eta(f): the structural diagonal
  ev <- numeric(length(mu_struct))

  if (any(nz)) {
    cc  <- arr$cc[nz]
    m2c <- .admMomF(mu_struct[nz], vf[nz], 2 * cc)
    a2  <- arr$a2[nz]; b2 <- arr$b2[nz]
    c1  <- arr$form[nz] == .ADM_RESID_COMBINED1
    # combined2: add then prop, in that order, so a purely ADDITIVE model (b2 = 0,
    # ms = 1) keeps the exact floating-point result it always had.
    v <- a2 + b2 * m2c
    if (any(c1)) {
      mc <- .admMomF(mu_struct[nz][c1], vf[nz][c1], cc[c1])
      # E[(a + b f^c)^2] = a^2 + 2ab E[f^c] + b^2 E[f^2c]
      v[c1] <- a2[c1] + 2 * sqrt(a2[c1] * b2[c1]) * mc + b2[c1] * m2c[c1]
    }
    ev[nz] <- v
  }
  # Transform-both-sides. E[y|eta] = m(f) is NONLINEAR in f, so unlike every other
  # form the mean scale is not a constant: ms becomes m'(f) (used to scale the
  # covariance) and mu_pred picks up a curvature term. Both reduce exactly to the
  # linear case when m(f) = ms*f. m/v come from Gauss-Hermite quadrature over the
  # residual; m'/m'' by central difference on that same quadrature, so the NLL and
  # its derivatives differentiate one consistent function.
  # Count endpoints. E[y|eta] is LINEAR in the mean argument for all three, so ms
  # is a constant (1, or N for binom) and only ev differs. E_eta[.] uses the same
  # E[f^2] = mu^2 + var_f identity as every other form.
  cp_ <- arr$form == .ADM_RESID_POIS
  if (any(cp_)) ev[cp_] <- mu_struct[cp_]                       # E[f]
  cb_ <- arr$form == .ADM_RESID_BINOM
  if (any(cb_)) {
    N <- arr$csz[cb_]; p <- mu_struct[cb_]
    ms[cb_] <- N                                                # E[y|eta] = N p
    mu[cb_] <- N * p
    ev[cb_] <- N * (p - (p * p + vf[cb_]))                      # N(E[p] - E[p^2])
  }
  cn_ <- arr$form == .ADM_RESID_NBINOM
  if (any(cn_)) {
    k <- arr$csz[cn_]; m <- mu_struct[cn_]
    ev[cn_] <- m + (m * m + vf[cn_]) / k                        # m + E[m^2]/k
  }
  # Ordinal: Var(1_k) = p_k(1-p_k) on the diagonal; the -p_j p_k cross terms are
  # off-diagonal and are emitted with rmat below (same channel as ar()).
  or_ <- arr$form == .ADM_RESID_ORDINAL
  if (any(or_)) {
    pk <- mu_struct[or_]
    ev[or_] <- pk - (pk * pk + vf[or_])          # E[p] - E[p^2]
  }
  bt_ <- arr$form == .ADM_RESID_BETA
  if (any(bt_)) {
    # E[y|eta] = mu (the derived output IS mu), Var = mu(1-mu)/(1+phi):
    # binom's shape with N -> 1/(1+phi).
    ph <- arr$phi
    if (is.null(ph)) ph <- rep(NA_real_, length(mu_struct))
    if (length(ph) == 1L) ph <- rep(ph, length(mu_struct))
    m_ <- mu_struct[bt_]
    ev[bt_] <- (m_ - (m_ * m_ + vf[bt_])) / (1 + ph[bt_])
  }
  tb <- arr$form == .ADM_RESID_TBS
  if (any(tb)) {
    for (i in which(tb)) {
      # ONE definition of this assembly, shared with .admResidDeriv's gradient
      # path -- see .admTBSRow(). m'/v' are ANALYTIC (.admTBSMomentsD); only the
      # second-order curvature differences the analytic first derivative.
      .r <- .admTBSRow(mu_struct[i], vf[i], arr$a2[i], arr$b2[i], arr$cc[i],
                       arr$lam[i], arr$yj[i], arr$tlo[i], arr$thi[i],
                       !is.null(arr$tbs_ftr) && isTRUE(arr$tbs_ftr[i]),
                       !is.null(arr$tbs_c1)  && isTRUE(arr$tbs_c1[i]),
                       arr$nodes %||% .ADM_TBS_NODES)
      ms[i] <- .r$ms                                          # m'(f), exact
      mu[i] <- .r$mu                                          # E_eta[m(f)]
      ev[i] <- .r$ev                                          # E_eta[Var(y|eta)]
    }
  }
  if (any(ln)) {
    sv  <- arr$a2[ln]
    mu[ln] <- mu_struct[ln] * exp(sv / 2)
    m2  <- vf[ln] + mu_struct[ln]^2
    # E[Var(y|eta)] = E[f^2] * exp(s) * (exp(s) - 1); the Var(E[y|eta]) part is
    # carried by the ms^2 scaling below, so do NOT also write an exp(2s)*vf term
    # -- that double-counts var_f (it cost an hour once; the two equivalent forms
    # are exp(s)*vf + E[f^2]*exp(s)(exp(s)-1) and exp(2s)*vf + mu^2*exp(s)(exp(s)-1)).
    ev[ln] <- m2 * exp(sv) * (exp(sv) - 1)
  }
  dvo <- ms^2 * vf + ev

  # ar(): pure off-diagonal. rxode2 scales the AR innovation so the MARGINAL
  # variance is unchanged, so the diagonal above is already right and the only
  # new term is the correlation between distinct observation times. The residual
  # variance that gets correlated is E[Var(y|eta)] (`ev`), NOT the total dvo --
  # the eta-driven part of the covariance is already carried by ms^2 (x) Cov_eta(f).
  rmat <- NULL
  # Ordinal cross-category covariance: -p_j p_k for categories at the SAME time.
  # arr$ord_grp labels each row's time group (NA for non-ordinal rows), so rows in
  # different groups stay uncorrelated.
  if (any(or_) && !is.null(times) && length(times) == length(mu_struct)) {
    # Two categories observed at the SAME time. By the law of total covariance
    #   Cov(1_j, 1_k) = E[Cov(1_j,1_k | eta)] + Cov_eta(p_j, p_k)
    #                 = E[-p_j p_k]           + Cov_eta(p_j, p_k)
    #                 = -(E p_j E p_k + Cov_eta) + Cov_eta
    #                 = -E[p_j] E[p_k],
    # i.e. the STRUCTURAL covariance cancels exactly. So this entry must REPLACE
    # V_struct, not add to it -- and since every caller does `V <- V + rmat`, the
    # emitted value carries the cancellation itself (-mu_j mu_k - cov_f). Getting
    # this wrong leaves the off-diagonal too large by exactly Cov_eta(p_j, p_k).
    # Rows at the same TIME are the same category group; ordinal rows for one time
    # are stacked across the joint unit's per-category blocks, so the time vector
    # identifies the group without any extra bookkeeping.
    g   <- .admOrdTimeGroup(times, or_)
    rm0 <- matrix(0, length(mu_struct), length(mu_struct))
    same <- outer(g, g, function(a, b) !is.na(a) & !is.na(b) & a == b)
    diag(same) <- FALSE
    cross <- -outer(mu, mu)
    if (!is.null(cov_f)) cross <- cross - cov_f
    rm0[same] <- cross[same]
    rmat <- rm0
  }
  .rho <- if (is.null(arr$rho)) NA_real_ else arr$rho[1L]
  if (!is.na(.rho) && .rho > 0 && !is.null(times) && length(times) == length(mu_struct)) {
    R <- .admARCor(.rho, times)
    if (!is.null(R)) {
      sd_e <- sqrt(pmax(ev, 0))
      rmat <- outer(sd_e, sd_e) * R
      diag(rmat) <- 0                 # diagonal already in dvo
    }
  }
  list(mu = mu, dv = dvo, ms = ms, ev = ev, rmat = rmat)
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
# for a "pow_exp" parameter p is the exponent itself and d(c)/dp = 1; for a
# "t_df" parameter p = log(nu - 2), so nu = 2 + exp(p).
#
# Student-t needs NO change to the a/b/c derivatives below. .admResidRows() has
# already folded the multiplier m into a2/b2, and every formula here is written in
# terms of those (scaled) coefficients, so each one differentiates the scaled
# variance correctly. Worked through for combined2: var = A + B f^2c with
# A = m*exp(p_a), so d(var)/d(p_a) = m*exp(p_a) = A = arr$a2 -- which is exactly
# what the unscaled code already returns. Same for combined1, where
# a = sqrt(A) and d(a)/d(p_a) = a/2 regardless of m.
#
# The ONLY new term is d(var)/d(p_nu). With m = nu/(nu-2) and nu = 2 + exp(p):
#   m = 1 + 2*exp(-p)  =>  dm/dp = -2*exp(-p) = -(m - 1)
# and since var = m * V0 (V0 the unscaled variance = var/m),
#   d(var)/d(p_nu) = V0 * dm/dp = -var * (m - 1) / m.
#   dv_dv0 n_t            d(V_pred_diag)/d(Var_eta(f))
#
# dv_dv0 is NEW and load-bearing: V_pred now depends on the STRUCTURAL variance
# through more than the identity (ms^2 on every row, plus dE[v]/dvar_f wherever
# the residual is f-dependent). It is the diagonal factor the estimators must
# apply to dNLL/dV_pred before handing it to the eta/omega kernels, which
# differentiate the structural covariance. Off-diagonals take ms_i*ms_j instead.
# `.admResidVChain()` packages both.
.admResidDeriv <- function(mu_struct, var_f, arr, pinfo) {
  n_t   <- length(mu_struct)
  n_sig <- length(pinfo$sigma_names)
  dmu   <- matrix(0, n_t, n_sig)
  dvar  <- matrix(0, n_t, n_sig)
  dv_df <- numeric(n_t)
  dv_dv0 <- rep(1.0, n_t)
  ev_resid <- numeric(n_t)      # E_eta[Var(y|eta)] per row; ar() correlates THIS
  # The mean scale reaches the OFF-diagonal of V_pred (V_pred_ij = ms_i ms_j cov_ij),
  # so anything ms depends on has a gradient path the row-indexed dvar/dv_df cannot
  # carry. dms is d(ms)/d(sigma_k) and dms_df is d(ms)/df = m''(f); .admSigmaGrad()
  # and .admResidMuCoupling() chain them over the off-diagonals. Both are zero
  # wherever ms is a constant (everything except lnorm and TBS), so additive,
  # proportional and count models are bit-identical.
  dms    <- matrix(0, n_t, n_sig)
  dms_df <- numeric(n_t)
  # d(mu_pred)/d(Var_eta(f)). Nonzero ONLY for TBS, where E_eta[m(f)] carries the
  # curvature term 0.5*m''(f)*var_f -- i.e. the PREDICTED MEAN depends on omega.
  # No kernel routes a mean-from-covariance path, so estimators fold this into the
  # diagonal of d(NLL)/d(V_struct) (see the `dmu_dv0` attribute of .admResidVChain).
  dmu_dv0 <- numeric(n_t)
  # d(mu_pred)/df beyond the plain mean scale. For TBS mu = m(f) + 0.5*m''(f)*var_f,
  # so the mean moves with f by m'(f) + 0.5*m'''(f)*var_f -- and that is NOT the same
  # number as `ms`, which scales the COVARIANCE and must stay m'(f) exactly. Keeping
  # one vector for both is what left the TBS struct-theta gradient ~1e-3 out.
  dmu_extra <- numeric(n_t)
  # The mean scale as .admResidApply() actually uses it. .admResidMuScale() only
  # knows the CONSTANT scales (1, or exp(s/2) for lnorm); for a TBS endpoint the
  # scale is m'(f), which depends on f and so cannot come from arr alone. Letting
  # the chain and the mu-coupling keep using .admResidMuScale() meant they applied
  # ms = 1 while the NLL applied m'(f) -- a 115% gradient error, constant across
  # every FD step size (which is what proved it a formula bug, not noise).
  ms_out <- .admResidMuScale(arr)
  if (is.null(var_f)) var_f <- numeric(n_t)
  # Third instance of the same trap (see .admBuildResidSpecs / .admResidRows):
  # pois/binom have NO residual parameters, but they still have f- and var_f-paths
  # (Var = f, and N p(1-p) respectively). Returning the defaults here left dv_df at
  # 0 and dv_dv0 at 1, so their structural and omega gradients were simply absent.
  # The SAME trap reopens for a continuous residual whose only parameter is fix()ed:
  # `.all_fixed_resid` keeps the spec alive with n_sig == 0, yet a fixed prop/pow/
  # lnorm/TBS residual is still prediction-dependent (Var(y|eta) moves with f). The
  # early return must fire ONLY when there is genuinely no f- or var_f-path -- i.e.
  # a purely additive residual (combined form with b^2 == 0, no count/lnorm/TBS),
  # where dv_df = 0 and dv_dv0 = 1 are the correct answers. `arr$b2 > 0` covers a
  # fixed prop/pow/combined coefficient; lnorm/TBS are f-dependent with any params.
  .f_dep <- any(arr$form %in% c(.ADM_RESID_POIS, .ADM_RESID_BINOM,
                                .ADM_RESID_NBINOM, .ADM_RESID_BETA,
                                .ADM_RESID_ORDINAL, .ADM_RESID_LNORM,
                                .ADM_RESID_TBS)) ||
            (!is.null(arr$b2) && any(is.finite(arr$b2) & arr$b2 > 0))
  if (n_sig == 0L && !.f_dep)
    return(list(dmu = dmu, dvar = dvar, dv_df = dv_df, dv_dv0 = dv_dv0,
                ev_resid = ev_resid, ms = ms_out, dms = dms, dms_df = dms_df,
                dmu_dv0 = dmu_dv0, dmu_df = ms_out))

  # A hand-built arr (unit tests predating t() support) carries neither field.
  vmul  <- if (is.null(arr$vmul))  rep(1.0, n_t)        else arr$vmul
  k_tdf <- if (is.null(arr$k_tdf)) rep(NA_integer_, n_t) else arr$k_tdf

  for (t in seq_len(n_t)) {
    f  <- mu_struct[t]
    v0 <- var_f[t]
    ka <- arr$k_add[t]; kb <- arr$k_prop[t]; kc <- arr$k_pow[t]

    if (arr$form[t] == .ADM_RESID_POIS) {
      # V = var_f + mu  ->  d/df = 1, d/d(var_f) = 1. No residual parameter at all.
      dv_df[t] <- 1; dv_dv0[t] <- 1; ev_resid[t] <- f
      next
    }
    if (arr$form[t] == .ADM_RESID_BINOM) {
      # mu_pred = N p, V = N^2 var_f + N(p - p^2 - var_f)
      N <- arr$csz[t]
      ms_out[t] <- N                       # E[y|eta] = N p
      dv_df[t]  <- N * (1 - 2 * f)
      dv_dv0[t] <- N * N - N
      ev_resid[t] <- N * (f - (f * f + v0))
      next
    }
    if (arr$form[t] == .ADM_RESID_NBINOM) {
      # V = var_f + m + (m^2 + var_f)/k
      k <- arr$csz[t]
      dv_df[t]  <- 1 + 2 * f / k
      dv_dv0[t] <- 1 + 1 / k
      ev_resid[t] <- f + (f * f + v0) / k
      # d(V)/d(k) = -(m^2 + var_f)/k^2; optimizer holds p = log(k), dk/dp = k.
      ks <- if (is.null(arr$k_size)) NA_integer_ else arr$k_size[t]
      if (!is.na(ks)) dvar[t, ks] <- -(f * f + v0) / k
      next
    }
    if (arr$form[t] == .ADM_RESID_ORDINAL) {
      # V_diag = var_f + p - (p^2 + var_f) = p - p^2  (the var_f cancels)
      dv_df[t]  <- 1 - 2 * f
      dv_dv0[t] <- 0
      ev_resid[t] <- f - (f * f + v0)
      next
    }
    if (arr$form[t] == .ADM_RESID_BETA) {
      ph <- if (is.null(arr$phi)) NA_real_ else arr$phi[t]
      dv_df[t]  <- (1 - 2 * f) / (1 + ph)
      dv_dv0[t] <- 1 - 1 / (1 + ph)
      ev_resid[t] <- (f - (f * f + v0)) / (1 + ph)
      next
    }
    if (arr$form[t] == .ADM_RESID_TBS) {
      # f-derivatives stay ANALYTIC (the quadrature's own partials, plus the sd
      # path now that sd depends on f); the four parameter directions (a, b, c,
      # lambda) are central differences of the SAME assembly .admResidApply()
      # builds, so analytic and objective agree by construction. Differencing a
      # smooth closed form: lambda measured 7.7e-11 against a reference FD.
      lam0 <- arr$lam[t]; yjc <- arr$yj[t]; lo <- arr$tlo[t]; hi <- arr$thi[t]
      f    <- mu_struct[t]
      hstep <- max(abs(f), 1) * 1e-4
      .ftr <- !is.null(arr$tbs_ftr) && isTRUE(arr$tbs_ftr[t])
      .c1i <- !is.null(arr$tbs_c1)  && isTRUE(arr$tbs_c1[t])

      # Everything .admResidApply assembles for THIS row, plus the pieces the
      # chain needs -- the SAME .admTBSRow() the objective calls, so analytic and
      # objective agree by construction rather than by two copies staying in step.
      # One row only, so a parameter direction costs O(1).
      .asm1 <- function(a2_, b2_, cc_, lam_)
        .admTBSRow(f, v0, a2_, b2_, cc_, lam_, yjc, lo, hi, .ftr, .c1i,
                   arr$nodes %||% .ADM_TBS_NODES)

      a20 <- max(arr$a2[t], 0); b20 <- max(arr$b2[t], 0); cc0 <- arr$cc[t]
      b0  <- .asm1(a20, b20, cc0, lam0)
      ms_out[t]   <- b0$ms
      dms_df[t]   <- b0$d2m                      # d(ms)/df = m''(f)
      dmu_dv0[t]  <- 0.5 * b0$d2m                # mu carries 0.5*m''*var_f
      dv_dv0[t]   <- b0$ms^2 + 0.5 * b0$d2v      # exact d(...)/d(var_f)
      # ev_raw, NOT ev: this path carries the UNCORRECTED q$v[2] (the objective's
      # `ev` adds the 0.5*v''(f)*v0 curvature term). ar() correlates rows by
      # sqrt(ev_i ev_j), and that is the quantity it has always used here.
      ev_resid[t] <- b0$ev_raw
      # One order beyond m''/v'': second differences of the analytic TOTAL first
      # derivatives, at the eps^(1/4) optimum for a second difference.
      m3 <- (b0$dm_t[3L] - 2 * b0$dm_t[2L] + b0$dm_t[1L]) / hstep^2
      v3 <- (b0$dv_t[3L] - 2 * b0$dv_t[2L] + b0$dv_t[1L]) / hstep^2
      dmu_extra[t] <- 0.5 * m3 * v0
      dv_df[t]     <- 2 * b0$ms * b0$d2m * v0 + b0$dv_t[2L] + 0.5 * v3 * v0

      # Parameter directions. The optimizer holds log(a^2) and log(b^2), so
      # d(a2)/dp = a2 and d(b2)/dp = b2; a pow exponent and lambda are identity.
      .cdiff <- function(hi_, lo_, hh, chain, k) {
        dmu[t, k]  <<- (hi_$mu - lo_$mu) / (2 * hh) * chain
        dvar[t, k] <<- (hi_$dv - lo_$dv) / (2 * hh) * chain
        dms[t, k]  <<- (hi_$ms - lo_$ms) / (2 * hh) * chain
      }
      if (!is.na(ka) && a20 > 0) {
        hh <- a20 * 1e-5
        .cdiff(.asm1(a20 + hh, b20, cc0, lam0), .asm1(a20 - hh, b20, cc0, lam0),
               hh, a20, ka)
      }
      kbp <- arr$k_prop[t]
      if (!is.na(kbp) && b20 > 0) {
        hh <- b20 * 1e-5
        .cdiff(.asm1(a20, b20 + hh, cc0, lam0), .asm1(a20, b20 - hh, cc0, lam0),
               hh, b20, kbp)
      }
      kcp <- arr$k_pow[t]
      if (!is.na(kcp)) {
        hh <- max(abs(cc0), 1) * 1e-5
        .cdiff(.asm1(a20, b20, cc0 + hh, lam0), .asm1(a20, b20, cc0 - hh, lam0),
               hh, 1, kcp)
      }
      klm <- if (is.null(arr$k_lam)) NA_integer_ else arr$k_lam[t]
      if (!is.na(klm)) {
        hh <- max(abs(lam0), 1) * 1e-5
        .cdiff(.asm1(a20, b20, cc0, lam0 + hh), .asm1(a20, b20, cc0, lam0 - hh),
               hh, 1, klm)
      }
      # Student-t degrees of freedom. nu reaches the objective ONLY through the
      # multiplier m = nu/(nu-2) that .admResidRows() already folded into a2/b2
      # (it does so for EVERY form, TBS included), so differentiate along that
      # scaling rather than re-deriving the quadrature. With p = log(nu-2),
      # dm/dp = -2/(nu-2) = -(m-1); for the closed forms below vt is proportional
      # to m and this same chain reduces to their -vt*(m-1)/m, so the two branches
      # agree by construction.
      #
      # Without this the TBS branch fell through `next` before the closed-form nu
      # block at the end of the loop, so `cp ~ boxCox(lam) + add(a) + t(nu)` had an
      # identically ZERO nu gradient while the NLL genuinely moved with nu -- the
      # optimizer left nu at its start value and drove every other parameter along
      # a direction the objective does not follow. Nothing refuses TBS + t() (only
      # lnorm + t() is refused), so the combination is reachable.
      ktd <- k_tdf[t]
      if (!is.na(ktd)) {
        m0 <- vmul[t]
        if (is.finite(m0) && m0 > 0) {
          hh <- m0 * 1e-5
          .cdiff(.asm1(a20 * (m0 + hh) / m0, b20 * (m0 + hh) / m0, cc0, lam0),
                 .asm1(a20 * (m0 - hh) / m0, b20 * (m0 - hh) / m0, cc0, lam0),
                 hh, -(m0 - 1), ktd)
        }
      }
      next
    }
    if (arr$form[t] == .ADM_RESID_LNORM) {
      sv   <- arr$a2[t]
      es   <- exp(sv)
      mu_s <- f * exp(sv / 2)
      # V = exp(s)*v0 + (v0 + f^2)*exp(s)*(exp(s)-1)
      if (!is.na(ka)) {
        dmu[t, ka]  <- mu_s * sv / 2
        # d(V)/ds = exp(s)*v0 + (v0+f^2)*exp(s)*(2exp(s)-1); optimizer holds log(s)
        dvar[t, ka] <- sv * (es * v0 + (v0 + f * f) * es * (2 * es - 1))
        # ms = exp(s/2) also multiplies every OFF-diagonal of V_pred; d(ms)/ds =
        # ms/2, times sv for the log chain. Omitting this made the lnorm sigma
        # gradient 3.5% wrong under method = "cov" (exactly right under "var",
        # where there is no off-diagonal -- which is why it hid for so long).
        dms[t, ka]  <- exp(sv / 2) * sv / 2
      }
      dv_df[t]  <- 2 * f * es * (es - 1)
      dv_dv0[t] <- es * es                      # exp(2s)
      ev_resid[t] <- (v0 + f * f) * es * (es - 1)
      next
    }

    cval <- arr$cc[t]
    sa   <- arr$a2[t]; sb <- arr$b2[t]
    M2   <- .admMomFd(f, v0, 2 * cval)          # E[f^2c]

    if (arr$form[t] == .ADM_RESID_COMBINED1) {
      # V = v0 + a^2 + 2ab*E[f^c] + b^2*E[f^2c]
      M1 <- .admMomFd(f, v0, cval)              # E[f^c]
      ab <- sqrt(sa * sb)
      if (!is.na(ka)) dvar[t, ka] <- sa + ab * M1$m
      if (!is.na(kb)) dvar[t, kb] <- ab * M1$m + sb * M2$m
      if (!is.na(kc)) dvar[t, kc] <- 2 * ab * M1$dk + sb * M2$dk * 2
      dv_df[t]  <- 2 * ab * M1$dmu + sb * M2$dmu
      dv_dv0[t] <- 1 + 2 * ab * M1$dv0 + sb * M2$dv0
      vt <- sa + 2 * ab * M1$m + sb * M2$m
      ev_resid[t] <- vt
    } else {
      # V = v0 + a^2 + b^2*E[f^2c]
      if (!is.na(ka)) dvar[t, ka] <- sa
      if (!is.na(kb)) dvar[t, kb] <- sb * M2$m
      if (!is.na(kc)) dvar[t, kc] <- sb * M2$dk * 2      # k = 2c => dk/dc = 2
      dv_df[t]  <- sb * M2$dmu
      dv_dv0[t] <- 1 + sb * M2$dv0
      vt <- sa + sb * M2$m
      ev_resid[t] <- vt
    }

    # Student-t degrees of freedom. The multiplier scales only the CONDITIONAL
    # variance E[Var(y|eta)] (= vt), never the ms^2*v0 term -- and lnorm+t is
    # refused, so ms == 1 whenever this fires.
    kt <- k_tdf[t]
    if (!is.na(kt)) {
      m <- vmul[t]
      dvar[t, kt] <- -vt * (m - 1) / m
    }
  }
  list(dmu = dmu, dvar = dvar, dv_df = dv_df, dv_dv0 = dv_dv0,
       ev_resid = ev_resid, ms = ms_out, dms = dms, dms_df = dms_df,
       dmu_dv0 = dmu_dv0, dmu_df = ms_out + dmu_extra)
}

# Does this observation unit observe any of `outs`? A joint (same-subject) unit
# stacks several outputs and names them on its blocks; an ordinary unit has one.
# Used by the ar()/ordinal guards, which are about a particular ENDPOINT and must
# not judge a unit that observes a different one.
.admUnitTouches <- function(u, outs) {
  if (!length(outs)) return(FALSE)
  got <- if (isTRUE(u$is_joint) && length(u$blocks %||% list()))
    vapply(u$blocks, function(b) as.character(b$output %||% NA_character_), character(1))
  else as.character(u$output %||% NA_character_)
  any(got %in% outs)
}

# ar() needs a FULL observed covariance. A `method = "var"` study contributes only
# diag(V), which contains no information about a residual correlation at all -- the
# fit would run, report whatever rho started at, and mean nothing. Refuse instead,
# naming the study so the user knows which one to fix.
# Ordinal endpoints: a diagonal V carries no cross-category information, and the
# same-time -p_j p_k structure is the whole point of the model. Refuse `var`
# studies for the same reason ar() does. Also refuse an ordinal endpoint that is
# NOT expressed as a joint unit -- one block per category probability -- because a
# non-joint study would score a single category and drop the rest.
.admCheckOrdinal <- function(pinfo, studies) {
  sp <- .admResidSpecs(pinfo)
  if (length(sp) == 0L ||
      !any(vapply(sp, function(x) identical(x$form, .ADM_RESID_ORDINAL), logical(1))))
    return(invisible(NULL))
  .ord <- vapply(sp, function(x) identical(x$form, .ADM_RESID_ORDINAL), logical(1))
  n_cat <- length(sp[[which(.ord)[1L]]]$ord_p)
  # ONLY the units that actually observe an ordinal category. Deciding from a
  # model-level scan and then rejecting every unit made a PK + ordinal model
  # (`cp ~ add(a); y ~ c(p1, p2)`) unfittable: the ordinary `cp` study is neither
  # joint nor supplies n_cat blocks, so it tripped a check about an endpoint it
  # has nothing to do with.
  studies <- Filter(function(u) .admUnitTouches(u, names(sp)[.ord]), studies)
  if (length(studies) == 0L) return(invisible(NULL))
  bad_v <- names(studies)[vapply(studies, function(u) identical(u$method, "var"), logical(1))]
  if (length(bad_v) > 0L)
    stop("An ordinal endpoint needs a full observed covariance.
",
         "  Stud", if (length(bad_v) > 1L) "ies" else "y", " ",
         paste(sQuote(bad_v), collapse = ", "), " supplied only variances ",
         "(method = \"var\"),
",
         "  but the -p_j*p_k covariance between categories observed at the SAME time
",
         "  is exactly what a diagonal V throws away. Supply the full V.", call. = FALSE)
  bad_j <- names(studies)[vapply(studies, function(u)
    !isTRUE(u$is_joint) || length(u$blocks %||% list()) != n_cat, logical(1))]
  if (length(bad_j) > 0L)
    stop("An ordinal endpoint must be given one observation block per category.
",
         "  Stud", if (length(bad_j) > 1L) "ies" else "y", " ",
         paste(sQuote(bad_j), collapse = ", "), " do not supply ", n_cat,
         " blocks.
",
         "  Each observation of an ordinal endpoint is the VECTOR of the ", n_cat,
         " category
  probabilities, so write e.g.
",
         "    studies = list(s = list(observations = list(
",
         "      c1 = list(output = \"p1\", times = tt, E = E1, V = V1, n = n),
",
         "      c2 = list(output = \"p2\", times = tt, E = E2, V = V2, n = n)),
",
         "      V = Vjoint, n = n, ev = ev))", call. = FALSE)
  invisible(NULL)
}

.admCheckAR <- function(pinfo, studies) {
  sp <- .admResidSpecs(pinfo)
  .isar <- vapply(sp, function(x)
    (!is.null(x$k_ar) && !is.na(x$k_ar)) ||
      (!is.null(x$ar_fixed) && is.finite(x$ar_fixed)), logical(1))
  if (length(sp) == 0L || !any(.isar)) return(invisible(NULL))

  # ONLY the units that actually observe an ar() endpoint -- see .admCheckOrdinal.
  # `cp ~ add(a) + ar(rho); ct ~ add(a2)` used to refuse a `ct` study whose V
  # happened to be diagonal (which .admNormaliseStudy auto-detects as
  # method = "var"), for an autocorrelation ct does not have.
  studies <- Filter(function(u) .admUnitTouches(u, names(sp)[.isar]), studies)
  if (length(studies) == 0L) return(invisible(NULL))

  # ar() inside a JOINT (same-subject, multi-output) unit is not representable.
  # The unit's rows stack several outputs, so its row times REPEAT, and
  # .admARCor's rho^|t_i - t_j| yields rho^0 = 1 between two DIFFERENT outputs
  # observed at the same time -- a perfect correlation nobody asked for. On top of
  # that .admResidApply reads arr$rho[1L], so one output's rho would be applied to
  # every row even when the other endpoint has no ar() term at all. ar() describes
  # autocorrelation WITHIN one output's series; across outputs it has no meaning,
  # so refuse rather than invent one.
  jnt <- names(studies)[vapply(studies, function(u) isTRUE(u$is_joint), logical(1))]
  if (length(jnt) > 0L)
    stop("ar() cannot be combined with a joint (same-subject, multi-output) study.\n",
         "  Stud", if (length(jnt) > 1L) "ies" else "y", " ",
         paste(sQuote(jnt), collapse = ", "), " stack several outputs, so the row\n",
         "  times repeat and rho^|dt| would correlate DIFFERENT outputs observed at\n",
         "  the same time. ar() describes autocorrelation within one output's series.\n",
         "  Fit the outputs separately, or drop ar().", call. = FALSE)

  bad <- names(studies)[vapply(studies, function(u)
    identical(u$method, "var"), logical(1))]
  if (length(bad) == 0L) return(invisible(NULL))
  stop("ar() residual correlation needs a full observed covariance.
",
       "  Stud", if (length(bad) > 1L) "ies" else "y", " ",
       paste(sQuote(bad), collapse = ", "), " supplied only variances (method = \"var\"),
",
       "  and a diagonal V carries no information about a residual correlation --
",
       "  rho would be reported at whatever value it started from.
",
       "  Supply the full covariance matrix V for th", 
       if (length(bad) > 1L) "ose studies" else "at study",
       ", or drop ar().",
       call. = FALSE)
}

# Can the fused C++ kernels score this residual, or must the moments be assembled
# in R first?
#
# adm_apply_residual() in src/nll.cpp implements forms 0/1/2 only (combined2,
# combined1, lnorm) and has no channel for an off-diagonal contribution. Every
# form added since -- TBS(3), pois(4), binom(5), nbinom(6), beta(7), ordinal(8) --
# would fall into its `else` branch and be scored as COMBINED2, and ar()'s
# correlation would be dropped entirely. That is silent: the NLL simply computes a
# different model than the gradient, which is how it was found (a 1.5e8 gradient
# mismatch for ar in the FD audit, with the sens and FD gradients agreeing with
# each other and both disagreeing with the NLL).
#
# Rather than duplicate quadrature and multinomial covariance in C++, such a study
# assembles its moments through .admResidApply() in R and then calls the plain
# nll_cov_cpp/nll_var_cpp kernels -- which take an already-built (mu, V) and are
# therefore form-agnostic. The fused path stays for the forms it genuinely covers.
.admResidCppOK <- function(arr) {
  all(arr$form <= .ADM_RESID_LNORM) &&
    (is.null(arr$rho) || all(is.na(arr$rho)))
}

# Chain factors converting d(NLL)/d(V_pred) into d(NLL)/d(V_struct), which is what
# the eta/omega kernels consume. V_pred = ms^2 (x) V_struct + diag(E[v]), so
#
#   off-diagonal (i != j):  ms_i * ms_j
#   diagonal     (i == i):  dv_dv0_i   (already contains ms_i^2)
#
# Returns the full n_t x n_t multiplier matrix; elementwise-multiply dNLL_dV by it.
.admResidVChain <- function(mu_struct, var_f, arr, pinfo, times = NULL,
                            deriv = NULL) {
  # `deriv` (LAST, optional): a precomputed .admResidDeriv(mu_struct, var_f, arr,
  # pinfo) the caller already holds. The three consumers of .admResidDeriv --
  # .admResidVChain, .admSigmaGrad, .admResidMuCoupling -- are each called once per
  # study/unit on the SAME (mu_struct, var_f, arr, pinfo), so computing it once in
  # the estimator and threading it here avoids two of the three quadratures a TBS/
  # logitNorm/probitNorm endpoint (resid_nodes, default 81) does per gradient eval.
  # NULL -> compute it, so every legacy/hand-built caller is unchanged.
  d  <- deriv %||% .admResidDeriv(mu_struct, var_f, arr, pinfo)
  M  <- tcrossprod(d$ms)          # NOT .admResidMuScale(): see .admResidDeriv$ms
  diag(M) <- d$dv_dv0
  # Ordinal: a same-time cross-category entry of V_pred is -E[p_j]E[p_k], which does
  # NOT depend on the structural covariance at all (it cancelled -- see
  # .admResidApply). d(V_pred)/d(V_struct) is therefore 0 there, not 1.
  or_ <- arr$form == .ADM_RESID_ORDINAL
  if (any(or_) && !is.null(times) && length(times) == length(mu_struct)) {
    g <- .admOrdTimeGroup(times, or_)
    same <- outer(g, g, function(a, b) !is.na(a) & !is.na(b) & a == b)
    diag(same) <- FALSE
    M[same] <- 0
  }
  # Carried alongside, not folded in, because it is an ADDITIVE term on a
  # DIFFERENT derivative: for TBS the predicted MEAN depends on Var_eta(f)
  # (mu = m(f) + 0.5*m''(f)*var_f), so omega reaches the objective through the mean
  # as well as through the covariance. No kernel has a mean-from-covariance input,
  # but adding dNLL_dmu * dmu_dv0 to the DIAGONAL of d(NLL)/d(V_struct) routes it
  # exactly: that diagonal is what every kernel multiplies by d(var_f)/d(param).
  # An attribute rather than a new return shape so existing callers are untouched.
  attr(M, "dmu_dv0") <- d$dmu_dv0
  M
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
# `dNLL_dV`/`times` are needed only by ar(), whose parameter enters the objective
# ONLY through off-diagonal entries -- so, unlike every other residual parameter,
# its gradient cannot be read off the diagonal `dNLL_dvar`. Omitting them simply
# leaves rho out of the analytic gradient (a var-branch study, where rho is
# unidentifiable anyway and is refused up front).
# ONE row of the TBS moment assembly -- the single definition shared by the
# OBJECTIVE (.admResidApply) and the GRADIENT (.admResidDeriv's `.asm1`).
#
# Both need the same thing: the residual SD on the TRANSFORMED scale and its d/df
# (rxode2's rx_r_), the quadrature moments over a 3-point f-stencil, and the
# second-order curvature terms differenced from the ANALYTIC first derivatives
# (one FD level, not two -- nesting FD here cost 300% accuracy). They were written
# out twice; if the two copies drift by a term the optimizer descends a direction
# the NLL does not follow, which is the worst failure class in this file.
#
# Returns BOTH conventions for E[Var(y|eta)] because the callers genuinely differ:
#   ev_raw = q$v[2]                      -- what the gradient path carries as
#                                           ev_resid (the ar() correlation uses it)
#   ev     = q$v[2] + 0.5*v''(f)*v0      -- the objective's curvature-corrected
#                                           E_eta[Var(y|eta)]
# `dv` is the full predicted variance for the row, ms^2*v0 + ev.
.admTBSRow <- function(f, v0, a2, b2, cc, lam, yj, lo, hi, ftr, c1, nodes) {
  hstep <- max(abs(f), 1) * 1e-4
  fv    <- c(f - hstep, f, f + hstep)
  a2 <- max(a2, 0); b2 <- max(b2, 0)
  # Residual SD on the TRANSFORMED scale, exactly as rxode2 builds rx_r_:
  #   combined2 : var = a^2 + (x^c * b)^2      combined1 : sd = a + x^c * b
  # where x is rx_pred_f_ = f, or rx_pred_ = h(f) when errTypeF is "transformed".
  if (ftr) {
    xb  <- .admTBS(fv, lam, yj, lo, hi)
    xbd <- 1 / .admTBSid(xb, lam, yj, lo, hi)          # dh/df
  } else { xb <- fv; xbd <- rep(1, length(fv)) }
  if (cc == 1) { pw <- xb; pwd <- rep(1, length(xb)) }
  else         { pw <- .admTBSp(xb, cc); pwd <- cc * .admTBSp(xb, cc - 1) }
  if (b2 == 0) {
    sdv <- rep(sqrt(a2), length(fv)); dsd <- numeric(length(fv))
  } else if (c1) {
    sdv <- sqrt(a2) + sqrt(b2) * pw
    dsd <- sqrt(b2) * pwd * xbd
  } else {
    .vv <- a2 + b2 * pw * pw
    sdv <- sqrt(.vv)
    dsd <- (b2 * pw * pwd * xbd) / pmax(sdv, .Machine$double.xmin)
  }
  q <- .admTBSMomentsD(fv, sdv, lam, yj, lo, hi, nodes)
  # sd depends on f, so the TOTAL derivative picks up the sd path (both partials
  # analytic): dm/df = @m/@f + @m/@sd * dsd/df
  dm_t <- q$dm_df + q$dm_ds * dsd
  dv_t <- q$dv_df + q$dv_ds * dsd
  d2m  <- (dm_t[3L] - dm_t[1L]) / (2 * hstep)
  d2v  <- (dv_t[3L] - dv_t[1L]) / (2 * hstep)
  .ms  <- dm_t[2L]
  .ev  <- q$v[2L] + 0.5 * d2v * v0
  list(mu = q$m[2L] + 0.5 * d2m * v0, dv = .ms^2 * v0 + .ev,
       ms = .ms, ev = .ev, ev_raw = q$v[2L], d2m = d2m, d2v = d2v,
       dm_t = dm_t, dv_t = dv_t, hstep = hstep)
}

# The OFF-diagonal mean-scale factor, shared by every consumer of dms/dms_df.
#
# V_pred = ms (x) ms o Cov_eta(f) + diag(E[Var(y|eta)]), so any parameter moving the
# mean scale `ms` reaches every off-diagonal entry:
#   d/dp sum_{i!=j} dNLL_dV_ij ms_i ms_j cov_ij = 2 * (dms/dp)' (A ms),
#   A = dNLL_dV o cov_f  with a ZERO diagonal (the diagonal belongs to dvar/dv_df).
# This returns the common `A %*% ms` factor; each caller applies its own
# contraction -- elementwise by dms_df for a per-ROW term, crossprod by the dms
# matrix for a per-SIGMA one.
#
# This is TBS/lnorm gradient maths, which errmodel.R owns, but it had been copied
# inline four times: .admSigmaGrad(), .admResidMuCoupling(), and TWICE inside
# .adghGrad(). Adding a form whose ms varies and missing the adgh copies would have
# given adgh a gradient that disagrees with its own objective -- the exact class of
# bug this file's notes blame for the recurring silent-wrong-model failures.
# Returns NULL when the term does not apply, so callers can skip it cheaply.
.admMsOffDiag <- function(dNLL_dV, cov_f, ms) {
  if (is.null(dNLL_dV) || is.null(cov_f) || is.null(ms)) return(NULL)
  A <- dNLL_dV * cov_f
  diag(A) <- 0
  drop(A %*% ms)
}

.admSigmaGrad <- function(mu_struct, arr, pinfo, dNLL_dvar, dNLL_dmu, var_f = NULL,
                          dNLL_dV = NULL, times = NULL, cov_f = NULL,
                          deriv = NULL) {
  d <- deriv %||% .admResidDeriv(mu_struct, var_f, arr, pinfo)
  g <- drop(crossprod(d$dvar, dNLL_dvar) + crossprod(d$dmu, dNLL_dmu))

  # The OFF-diagonal ms path. V_pred = ms(x)ms o Cov_eta(f) + diag(E[v]), and for
  # lnorm/TBS the scale ms depends on a residual parameter, so sigma reaches every
  # off-diagonal entry as well as the diagonal that dvar already covers:
  #   d/dp sum_{i!=j} dNLL_dV_ij ms_i ms_j cov_ij = 2 * dms' (A ms),  A = dNLL_dV o cov_f
  # with A's diagonal zeroed (the diagonal is dvar's). Needs the STRUCTURAL
  # covariance, which the caller has before it applies the ms scaling; without
  # cov_f the term is skipped, which is exact for every form with a constant ms.
  if (!is.null(cov_f) && !is.null(dNLL_dV) && !is.null(d$dms) && any(d$dms != 0, na.rm = TRUE))
    g <- g + 2 * drop(crossprod(d$dms, .admMsOffDiag(dNLL_dV, cov_f, d$ms)))

  # ar(): rmat[i,j] = sqrt(ev_i ev_j) * rho^d_ij  (i != j), rho = expit(p), so
  #   d(rmat)/dp = sqrt(ev_i ev_j) * d_ij * rho^(d_ij - 1) * rho(1 - rho)
  # and the chain is over the OFF-diagonal entries of dNLL_dV.
  kar <- if (is.null(arr$k_ar)) NA_integer_ else arr$k_ar[1L]
  rho <- if (is.null(arr$rho))  NA_real_    else arr$rho[1L]
  if (!is.na(kar) && !is.na(rho) && rho > 0 &&
      !is.null(dNLL_dV) && !is.null(times) && length(times) == length(mu_struct)) {
    ev   <- d$ev_resid
    sd_e <- sqrt(pmax(ev, 0))
    dmat <- abs(outer(times, times, "-"))
    drho <- outer(sd_e, sd_e) * dmat * rho^(dmat - 1) * (rho * (1 - rho))
    diag(drho) <- 0
    g[kar] <- g[kar] + sum(dNLL_dV * drho)

    # Every OTHER residual parameter also reaches the off-diagonal, through the
    # sqrt(ev_i*ev_j) scale of rmat -- missing this leaves the add/prop gradient
    # wrong by ~25% whenever ar() is active (measured). With
    #   rmat_ij = sqrt(ev_i ev_j) rho^d_ij
    #   d/dp    = rho^d_ij * (dev_i * ev_j + ev_i * dev_j) / (2 sqrt(ev_i ev_j))
    # d(ev)/dp is d$dvar: ar() is refused with lnorm, so ms == 1 and the ms^2*v0
    # part of dv carries no sigma dependence, leaving dvar == d(ev)/dp exactly.
    rpow <- rho^dmat; diag(rpow) <- 0
    # rxode2's safeZero (_div0_): a zero denominator becomes DBL_EPSILON. Was Inf
    # here, which silently ZEROED the term instead of letting it blow up -- a
    # different answer, and not one rxode2 would give.
    den  <- outer(sd_e, sd_e); den[den == 0] <- .Machine$double.eps
    for (k in seq_len(ncol(d$dvar))) {
      if (k == kar) next
      dk <- d$dvar[, k]
      if (all(dk == 0)) next
      g[k] <- g[k] + sum(dNLL_dV * (rpow * (outer(dk, ev) + outer(ev, dk)) / (2 * den)))
    }
  }
  g
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
.admResidMuCoupling <- function(mu_struct, arr, pinfo, dNLL_dvar, dNLL_dmu,
                                var_f = NULL, dNLL_dV = NULL, cov_f = NULL,
                                times = NULL, deriv = NULL) {
  d <- deriv %||% .admResidDeriv(mu_struct, var_f, arr, pinfo)
  # d$dmu_df, NOT d$ms: they differ for TBS (see .admResidDeriv$dmu_extra). Callers
  # already carry the plain dNLL_dmu, so only the excess over 1 belongs here.
  out <- dNLL_dmu * (d$dmu_df - 1) + dNLL_dvar * d$dv_df
  # For a TBS endpoint ms = m'(f) depends on f, so moving the structural mean also
  # moves the OFF-diagonal of V_pred (V_pred_ij = ms_i ms_j cov_ij). Row k picks up
  #   2 * m''(f_k) * (A ms)_k,   A = dNLL_dV o cov_f with a zero diagonal,
  # which is exactly the same contraction .admSigmaGrad() applies for the sigma
  # path. Zero unless ms varies with f (i.e. TBS only), so nothing else moves.
  if (!is.null(cov_f) && !is.null(dNLL_dV) && !is.null(d$dms_df) && any(d$dms_df != 0, na.rm = TRUE))
    out <- out + 2 * d$dms_df * .admMsOffDiag(dNLL_dV, cov_f, d$ms)
  # ORDINAL: V_pred[i,j] = -E[p_i]E[p_j] for i != j at the same time (the structural
  # covariance cancels -- see .admResidApply), so mu_i reaches those entries directly
  # with d(V_pred[i,j])/d(mu_i) = -mu_j. Nothing else carries that: dv_df is the
  # diagonal and dms_df is zero for ordinal. Omitting it left this term with the
  # WRONG SIGN and ~16x relative error against a numeric reference.
  #
  # Currently unreachable -- .admLoadSensModel() returns NULL for ordinal, which puts
  # every estimator on a finite-difference gradient -- but the omission was silent
  # and depended on that invariant holding, so carry the term rather than the risk.
  or_ <- arr$form == .ADM_RESID_ORDINAL
  if (any(or_) && !is.null(dNLL_dV) && !is.null(times) &&
      length(times) == length(mu_struct)) {
    g <- .admOrdTimeGroup(times, or_)
    same <- outer(g, g, function(a, b) !is.na(a) & !is.na(b) & a == b)
    diag(same) <- FALSE
    out <- out + 2 * drop((dNLL_dV * same) %*% (-mu_struct))
  }
  out
}
