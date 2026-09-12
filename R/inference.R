# =============================================================================
# anova() on nested admFits
# =============================================================================
#
# The ORDINARY likelihood-ratio test: the objective difference against a
# chi-squared reference with as many degrees of freedom as the larger model has
# extra parameters. No reweighting, no sandwich, and no dependence on which
# `covMethod` the fits used -- anova(full, reduced) works on what they already
# carry.
#
# WHAT THIS DELIBERATELY DOES NOT DO. Under misspecification the plain LRT is
# not exactly chi-squared: H and J disagree and the exact reference becomes a
# weighted sum of chi-squares. Correcting for that needs H and J from both fits
# and a series evaluation of the weighted reference, and it constrains which
# covMethod a fit must have used. That is the wrong default for the ordinary
# question "is this extra parameter worth it".

# The optimizer-scale parameter names a fit was built on.
.admFitParNames <- function(fit) {
  e <- tryCatch(fit$env, error = function(e) NULL)
  if (is.null(e)) return(NULL)
  (e$admExtra %||% e$adirmcExtra)$par_names
}

# One nested comparison: dOFV, its degrees of freedom, and the p-value.
.admLRT <- function(full, reduced) {
  nm_f <- .admFitParNames(full)
  nm_r <- .admFitParNames(reduced)
  if (is.null(nm_f) || is.null(nm_r))
    stop("anova(): could not read the fits' parameter names.", call. = FALSE)
  # NESTING IS A PRECONDITION, not something to approximate around. A
  # non-nested pair is a different problem (Vuong) and must not come back with
  # a p-value.
  if (!all(nm_r %in% nm_f))
    stop("anova(): the models are not nested -- the smaller one has parameters ",
         "the larger does not (", paste(setdiff(nm_r, nm_f), collapse = ", "),
         "). A non-nested comparison needs a Vuong test, not this.",
         call. = FALSE)
  gamma <- setdiff(nm_f, nm_r)
  if (!length(gamma))
    stop("anova(): the two fits have the same parameters; nothing to test.",
         call. = FALSE)
  # THE SAME ARGUMENT, FOR THE ESTIMATOR AND FOR THE GRID. A -2LL is only a
  # likelihood ratio against another -2LL computed the same way, and the
  # estimators do not compute it the same way: adfo scores FO-LINEARISED
  # moments, adgh a quadrature, admc a Monte Carlo sample. anova(adfo_fit,
  # adgh_fit) differenced two numbers on different scales and returned a
  # perfectly finite p. Two adgh fits on different node counts have the same
  # problem: the objective moves with the grid, which is why .adghGrid refuses
  # to change the point count mid-fit in the first place.
  .m1 <- tryCatch(full$env$method, error = function(e) NULL)
  .m2 <- tryCatch(reduced$env$method, error = function(e) NULL)
  if (!is.null(.m1) && !is.null(.m2) && !identical(.m1, .m2))
    stop("anova(): these fits were made by different estimators (", .m1,
         " and ", .m2, "). Each scores its own approximation to the same ",
         "likelihood -- FO-linearised, quadrature or Monte Carlo -- so the ",
         "difference of their objectives is not a likelihood ratio. Refit ",
         "both with the same `est` before comparing.", call. = FALSE)
  .n1 <- tryCatch(full$env$nNodes, error = function(e) NULL)
  .n2 <- tryCatch(reduced$env$nNodes, error = function(e) NULL)
  if (!is.null(.n1) && !is.null(.n2) && !identical(.n1, .n2))
    stop("anova(): these fits used different node counts (", .n1, " and ",
         .n2, "). The objective moves with the grid, so the difference is not ",
         "a likelihood ratio -- refit both with the same `n_nodes`.",
         call. = FALSE)
  # SAME ARGUMENT, FOR admc/adirmc: their objective is a Monte Carlo average
  # over `n_sim` draws, so it moves with `n_sim` exactly as adgh's moves with
  # `n_nodes`. NULL for adfo/adgh, so this is a no-op there.
  .s1 <- tryCatch(full$env$nSim, error = function(e) NULL)
  .s2 <- tryCatch(reduced$env$nSim, error = function(e) NULL)
  if (!is.null(.s1) && !is.null(.s2) && !identical(.s1, .s2))
    stop("anova(): these fits used different `n_sim` (", .s1, " and ",
         .s2, "). The objective is a Monte Carlo average, so the difference ",
         "is not a likelihood ratio -- refit both with the same `n_sim`.",
         call. = FALSE)
  o_f <- as.numeric(full$objective)
  o_r <- as.numeric(reduced$objective)
  if (!is.finite(o_f) || !is.finite(o_r))
    stop("anova(): one of these fits has no finite objective.", call. = FALSE)
  d  <- o_r - o_f
  df <- length(gamma)
  # A NEGATIVE dOFV IS REPORTED, NOT CLAMPED. The larger model cannot fit worse
  # at its own optimum, so a negative difference says one of the two did not
  # converge -- which the user needs to see, not have rounded to zero.
  list(dOFV = d, df = df, gamma = gamma,
       p = if (d > 0) stats::pchisq(d, df, lower.tail = FALSE) else NA_real_)
}

#' Compare nested admixr2 fits by a likelihood-ratio test
#'
#' The ordinary LRT: the objective difference against a chi-squared reference
#' with `Df` equal to the number of parameters the larger model adds.
#'
#' Both fits must come from the same estimator and, for the quadrature
#' estimators, the same node count. Each scores its own approximation to the
#' likelihood, so objectives from different ones are not comparable and the
#' comparison is refused rather than reported.
#'
#' Testing a variance AT ZERO puts the null on the boundary of the parameter
#' space, where the exact reference is a chi-bar-squared mixture rather than a
#' chi-squared. The p-value reported there is CONSERVATIVE -- too large -- so a
#' significant result stays significant, but treat a borderline one with care.
#'
#' @param object An `admFit`.
#' @param ... Further `admFit`s to compare it with.
#' @return A data frame of class `anova.admFit`, smallest model first.
#' @export
anova.admFit <- function(object, ...) {
  fits <- c(list(object), list(...))
  ok   <- vapply(fits, function(z) inherits(z, "admFit"), logical(1))
  if (!all(ok))
    stop("anova(): every argument must be an admFit; got ",
         paste(vapply(fits[!ok], function(z) class(z)[1L], character(1)),
               collapse = ", "), ".", call. = FALSE)
  if (length(fits) < 2L)
    stop("anova() on a single admFit has nothing to compare it with. Pass the ",
         "nested pair, e.g. anova(full, reduced).", call. = FALSE)
  npar <- vapply(fits, function(z) length(.admFitParNames(z)), integer(1))
  # Smallest first, the way anova.lme and anova.lm order their tables, so each
  # test row reads "the model above, plus these parameters".
  ord  <- order(npar)
  fits <- fits[ord]; npar <- npar[ord]
  nms  <- make.unique(vapply(fits, function(z) z$env$method %||% "adm",
                             character(1)), sep = "_")
  n   <- length(fits)
  num <- function(nm) vapply(fits, function(z)
    as.numeric(z$env[[nm]] %||% NA_real_), numeric(1))
  out <- data.frame(
    Npar = npar,
    OBJF = vapply(fits, function(z) as.numeric(z$objective), numeric(1)),
    AIC  = num("AIC"),
    BIC  = num("BIC"),
    Test = c(NA_character_, sprintf("%d vs %d", seq_len(n - 1L), 2:n)),
    dOFV = NA_real_,
    Df   = NA_integer_,
    p    = NA_real_,
    row.names        = nms,
    stringsAsFactors = FALSE,
    check.names      = FALSE)
  for (i in seq_len(n - 1L)) {
    r <- .admLRT(fits[[i + 1L]], fits[[i]])
    out$dOFV[i + 1L] <- r$dOFV
    out$Df[i + 1L]   <- r$df
    out$p[i + 1L]    <- r$p
  }
  class(out) <- c("anova.admFit", "data.frame")
  attr(out, "heading") <- "Likelihood-ratio test"
  out
}

#' @export
print.anova.admFit <- function(x, ...) {
  # NOT stats::print.anova. It routes through printCoefmat, which calls
  # data.matrix() on the frame -- and data.matrix() turns a CHARACTER column
  # into its factor CODES, so `Test` would print as 1, 2, ... instead of
  # "1 vs 2". Printing the frame directly keeps it as written.
  h <- attr(x, "heading")
  if (!is.null(h)) cat(paste(h, collapse = "\n"), "\n", sep = "")
  y <- x
  attr(y, "heading") <- NULL
  class(y) <- "data.frame"
  print(format(y, digits = max(3L, getOption("digits") - 3L)), ...)
  invisible(x)
}
