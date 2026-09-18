# =============================================================================
# anova() on nested admFits
# =============================================================================
# Ordinary LRT vs a chi-squared reference (df = added parameters). Deliberately
# does not correct for H/J disagreement under misspecification (needs H and J
# from both fits, see R/adfweight.R) -- wrong default for "is this parameter
# worth it".

# The optimizer-scale parameter names a fit was built on.
.admFitParNames <- function(fit) {
  e <- tryCatch(fit$env, error = function(e) NULL)
  if (is.null(e)) return(NULL)
  (e$admExtra %||% e$adirmcExtra)$par_names
}

# One nested comparison: dOFV, its degrees of freedom, and the p-value.
.admFitHasModelSource <- function(fit) {
  e <- tryCatch(fit$env, error = function(e) NULL)
  !is.null(e) && isTRUE((e$admExtra %||% e$adirmcExtra)$has_model_source)
}

.admLRT <- function(full, reduced) {
  if (.admFitHasModelSource(full) || .admFitHasModelSource(reduced))
    stop("anova(): unavailable for fits containing a published model source because its parameter sampling law is unknown.", call. = FALSE)
  nm_f <- .admFitParNames(full)
  nm_r <- .admFitParNames(reduced)
  if (is.null(nm_f) || is.null(nm_r))
    stop("anova(): could not read the fits' parameter names.", call. = FALSE)
  # Nesting check: smaller model parameters must be a subset of the larger model
  if (!all(nm_r %in% nm_f))
    stop("anova(): the models are not nested -- the smaller one has parameters ",
         "the larger does not (", paste(setdiff(nm_r, nm_f), collapse = ", "),
         "). A non-nested comparison needs a Vuong test, not this.",
         call. = FALSE)
  gamma <- setdiff(nm_f, nm_r)
  if (!length(gamma))
    stop("anova(): the two fits have the same parameters; nothing to test.",
         call. = FALSE)
  # Verify identical estimation method and quadrature grid
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
  # Verify identical Monte Carlo sample size (admc/adirmc)
  .s1 <- tryCatch(full$env$nSim, error = function(e) NULL)
  .s2 <- tryCatch(reduced$env$nSim, error = function(e) NULL)
  if (!is.null(.s1) && !is.null(.s2) && !identical(.s1, .s2))
    stop("anova(): these fits used different `n_sim` (", .s1, " and ",
         .s2, "). The objective is a Monte Carlo average, so the difference ",
         "is not a likelihood ratio -- refit both with the same `n_sim`.",
         call. = FALSE)
  # STRATUM RESOLUTION, per covariate, compared only where the two fits OVERLAP.
  # A covariate only one of them reads cannot put them on different scales: if a
  # model's prediction does not move across a source's nodes, the mixture those
  # nodes collapse to is a sufficient statistic for it, so its objective is the
  # same at either resolution (measured: 4.7e-05 on a dOFV of 409). Comparing
  # one number per fit refused exactly the nested pair a covariate test is made
  # of, because the null model had dropped the term and with it the nodes.
  .j1 <- tryCatch(full$env$strataNodes, error = function(e) NULL)
  .j2 <- tryCatch(reduced$env$strataNodes, error = function(e) NULL)
  # AN UNSET STAMP IS NOT A DISAGREEMENT: a model reading no covariate records
  # nothing, which is the null of a single-covariate test. BOTH sides must be
  # named to compare per covariate -- with `||`, one named side entered that
  # branch, where an unnamed stamp intersects to nothing and passes, so a
  # pre-rename `5` was accepted against `c(WT = 9)`. as.character() throughout,
  # since the stamp carries "3/9" now and a saved integer must still match.
  .nm1 <- names(.j1); .nm2 <- names(.j2)
  .named <- !is.null(.nm1) && !is.null(.nm2)
  .sh <- if (.named) intersect(.nm1, .nm2) else NULL
  .bad <- if (!length(.j1) || !length(.j2)) FALSE
          else if (.named) length(.sh) > 0L &&
            !identical(as.character(.j1[.sh]), as.character(.j2[.sh]))
          else !identical(as.character(unname(.j1)),
                          as.character(unname(.j2)))
  if (isTRUE(.bad)) {
    .w <- if (!.named) "" else {
      .d <- .sh[as.character(.j1[.sh]) != as.character(.j2[.sh])]
      paste0(" on ", paste(sQuote(.d), collapse = ", "))
    }
    # ", " between covariates: "/" separates node counts within one stamp.
    stop("anova(): these fits were built at different stratum resolutions",
         .w, " (", paste(.j1, collapse = ", "), " and ",
         paste(.j2, collapse = ", "), "). The objective is J-dependent, so ",
         "their difference is not a likelihood ratio -- refit both at the ",
         "same resolution.", call. = FALSE)
  }
  o_f <- as.numeric(full$objective)
  o_r <- as.numeric(reduced$objective)
  if (!is.finite(o_f) || !is.finite(o_r))
    stop("anova(): one of these fits has no finite objective.", call. = FALSE)
  d  <- o_r - o_f
  df <- length(gamma)
  # Report negative dOFV directly without clamping (indicates convergence failure)
  list(dOFV = d, df = df, gamma = gamma,
       p = if (d > 0) stats::pchisq(d, df, lower.tail = FALSE) else NA_real_)
}

#' Compare nested admixr2 fits by a likelihood-ratio test
#'
#' The ordinary LRT: the objective difference against a chi-squared reference
#' with `Df` equal to the number of parameters the larger model adds.
#'
#' Both fits must come from the same estimator and, for the quadrature
#' estimators, the same node count (`n_nodes`) or, for the Monte Carlo ones
#' (`admc`, `adirmc`), the same sample size (`n_sim`). Each scores its own
#' approximation to the likelihood, so objectives from different ones are not
#' comparable and the comparison is refused rather than reported.
#'
#' The **stratum** resolution is checked per covariate, and only where the two
#' fits overlap. A covariate a model does not read cannot put the two on
#' different scales: if its prediction does not move across a source's nodes,
#' the mixture those nodes collapse to is a sufficient statistic for it, so its
#' objective is the same at either resolution. That is what makes the nested
#' pair of a covariate test comparable --- the null model drops the term, so its
#' sources are not cut along it. Two fits that both read a covariate and cut it
#' differently, including one cutting it and the other integrating over it
#' whole, are refused.
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
  # Order models smallest first (by parameter count)
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
  # Print data.frame directly to prevent printCoefmat from converting character 'Test' to factors
  h <- attr(x, "heading")
  if (!is.null(h)) cat(paste(h, collapse = "\n"), "\n", sep = "")
  y <- x
  attr(y, "heading") <- NULL
  class(y) <- "data.frame"
  print(format(y, digits = max(3L, getOption("digits") - 3L)), ...)
  invisible(x)
}
