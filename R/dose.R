#' A dosing table by covariate band
#'
#' Transcribes a paper's dosing table, e.g. "CrCl below 30 mL/min: 100 mg;
#' otherwise 200 mg", for the `dose` argument of [admStudy()]:
#' `dose = doseBands(CRCL, breaks = 30, amt = c(100, 200))`.
#'
#' `breaks` are the cut points between bands and `amt` the amount for each band,
#' lowest first, so there is always one more amount than cut points and the
#' bands cover the whole range with no gaps. A value equal to a cut point falls
#' in the band above it. Anything else about how each dose is given ---
#' `ii`, `addl`, `dur`, `rate`, `cmt`, ... --- is written as in [rxode2::et()]
#' and applies to every band.
#'
#' The study is expanded into one study per band: the covariate restricted to
#' the band (and to the study's enrolled `range`, if given), `n` split by the
#' band's share of the population, and that band's dose, so the analysis model
#' needs no dose rule. Covariates correlated with the banded one are conditioned
#' on the band.
#'
#' @param covariate The covariate the bands are cut on, bare or quoted
#'   (`CRCL` or `"CRCL"`); a variable holding the name also works. It must be declared, as a continuous covariate, in
#'   the study's `population`.
#' @param breaks Increasing cut points between bands.
#' @param amt Amount for each band, lowest band first; one more than `breaks`.
#' @param ... How each dose is given, as [rxode2::et()] arguments: `time`,
#'   `cmt`, `ii`, `addl`, `dur`, `rate`, `ss`, `until`. Without them each band
#'   gets one bolus at time 0.
#' @return An `admDoseBands` object, to pass as `dose` to [admStudy()].
#' @seealso [admStudy()]; a function of covariates as `dose` for dosing that is
#'   proportional to (or a formula of) a covariate, e.g. mg/kg.
#' @examples
#' doseBands(CRCL, breaks = 30, amt = c(100, 200))
#' # a 1-h infusion every 12 h, seven doses, in every band
#' doseBands(CRCL, breaks = 30, amt = c(500, 1000), ii = 12, addl = 6, dur = 1)
#' @export
doseBands <- function(covariate, breaks, amt, ...) {
  cv <- substitute(covariate)
  if (is.name(cv)) {
    # a bare name is the covariate (doseBands(CRCL, ...)), unless it is a
    # variable holding a covariate name (cv <- "CRCL"; doseBands(cv, ...))
    v <- tryCatch(covariate, error = function(e) NULL)
    cv <- if (is.character(v) && length(v) == 1L) v else as.character(cv)
  } else cv <- covariate
  bad <- function(...) stop("admixr2: doseBands(): ", ..., call. = FALSE)
  if (!is.character(cv) || length(cv) != 1L || is.na(cv) || !nzchar(cv))
    bad("`covariate` must name one covariate, e.g. doseBands(CRCL, ...).")
  if (!is.numeric(breaks) || !length(breaks) || any(!is.finite(breaks)) ||
      is.unsorted(breaks, strictly = TRUE))
    bad("`breaks` must be finite, strictly increasing cut points, e.g. ",
        "breaks = 30 or breaks = c(30, 60).")
  if (!is.numeric(amt) || length(amt) != length(breaks) + 1L)
    bad("`amt` needs one amount per band: ", length(breaks), " cut point",
        if (length(breaks) > 1L) "s make " else " makes ",
        length(breaks) + 1L, " bands, but ", length(amt),
        if (length(amt) == 1L) " amount was" else " amounts were", " given.")
  if (any(!is.finite(amt)) || any(amt < 0))
    bad("each `amt` must be a finite, non-negative amount.")
  timing <- list(...)
  if (length(timing)) {
    tn <- names(timing)
    if (is.null(tn) || any(!nzchar(tn)))
      bad("the dosing details in `...` must be named, e.g. ii = 12, addl = 6.")
    # rxode2::et() IGNORES arguments it does not know, so a typo (`adl = 6`)
    # would silently turn a regimen into a single dose.
    unk <- setdiff(tn, .ADM_ET_TIMING)
    if (length(unk))
      bad(paste(sQuote(unk), collapse = ", "), " is not an rxode2::et() ",
          "dosing argument this applies to every band; use ",
          paste(.ADM_ET_TIMING, collapse = ", "), ".")
    tryCatch(do.call(rxode2::et, c(list(amt = 1), timing)), error = function(e)
      bad("these dosing details do not build an event table: ",
          conditionMessage(e)))
  }
  structure(list(covariate = cv, breaks = as.numeric(breaks),
                 amt = as.numeric(amt), timing = timing),
            class = "admDoseBands")
}

#' @export
format.admDoseBands <- function(x, ...) {
  b <- x$breaks; a <- x$amt; k <- length(a); cv <- x$covariate
  lab <- character(k)
  lab[1L] <- sprintf("%s < %g: %g", cv, b[1L], a[1L])
  if (k > 2L) for (j in 2:(k - 1L))
    lab[j] <- sprintf("%g to < %g: %g", b[j - 1L], b[j], a[j])
  lab[k] <- sprintf("%s >= %g: %g", cv, b[k - 1L], a[k])
  out <- paste(lab, collapse = "; ")
  if (length(x$timing))
    out <- paste0(out, "  (", paste(sprintf("%s = %s", names(x$timing),
      vapply(x$timing, function(v) paste(format(v), collapse = ","),
             character(1))), collapse = ", "), ")")
  out
}

#' @export
print.admDoseBands <- function(x, ...) {
  cat("dose by ", x$covariate, " band: ", format(x), "\n", sep = "")
  invisible(x)
}

# The rxode2::et() arguments that say WHEN and HOW a dose is given, not how
# much: what doseBands() passes on to every band.
.ADM_ET_TIMING <- c("time", "cmt", "ii", "addl", "dur", "rate", "ss", "until")

# What kind of dose rule `dose` is: a fixed amount, a full event table, a
# banded table, or a function of covariates. NA for anything else.
.admDoseKind <- function(d) {
  if (inherits(d, "admDoseBands")) return("bands")
  if (is.function(d)) return("fun")
  if (inherits(d, "rxEt")) return("ev")
  if (is.numeric(d)) return("amt")
  NA_character_
}

# The covariates a dose rule reads: the banded one, or a function's arguments.
.admDoseCovs <- function(d) {
  if (inherits(d, "admDoseBands")) return(d$covariate)
  if (is.function(d)) return(names(formals(d)))
  character(0)
}

# The bands as a lo/hi/amt table (open at both ends).
.admDoseBandTable <- function(d)
  data.frame(lo = c(-Inf, d$breaks), hi = c(d$breaks, Inf), amt = d$amt)

# One band's event table: its amount under the table's shared dosing details.
.admDoseBandEv <- function(d, k)
  do.call(rxode2::et, c(list(amt = d$amt[k]), d$timing))

# A dose function's event table at one node, checked. The function returns an
# amount (one bolus at time 0) or an rxode2::et() table, exactly as an nlmixr2
# user writes a subject's dosing; every amount must be finite and
# non-negative. The error names the node, since that is what the user can act on.
.admDoseEvAt <- function(f, cov) {
  cv  <- names(formals(f))
  arg <- as.list(cov)[cv]
  at  <- paste(sprintf("%s = %s", cv, vapply(arg, function(v)
    format(signif(as.numeric(v), 4)), character(1))), collapse = ", ")
  r <- tryCatch(do.call(f, arg), error = function(e)
    stop("the dose function failed at ", at, ": ", conditionMessage(e),
         call. = FALSE))
  ev <- if (inherits(r, "rxEt")) r
        else if (is.numeric(r) && length(r) == 1L) NULL
        else stop("the dose function must return one amount or an ",
                  "rxode2::et() event table; at ", at, " it returned ",
                  class(r)[1L], if (is.numeric(r)) paste0(" of length ", length(r)),
                  ".", call. = FALSE)
  a <- if (is.null(ev)) r else {
    d <- as.data.frame(ev)
    # only dosing records carry an amount: observations (0), other-type (2)
    # and reset (3) rows legitimately have amt = NA
    d$amt[!is.na(d$evid) & !(d$evid %in% c(0L, 2L, 3L))]
  }
  if (!length(a) || any(!is.finite(a)) || any(a < 0))
    stop("the dose function gives ", paste(format(a), collapse = ", "),
         " at ", at, "; every amount must be finite and non-negative. A ",
         "covariate whose margin reaches zero or below (a NORMAL weight) can ",
         "be bounded with `range`, or declared with dist = \"lnorm\".",
         call. = FALSE)
  ev %||% rxode2::et(amt = as.numeric(r))
}

# One line describing a dose rule, for print.admStudy().
.admDoseLabel <- function(d) {
  switch(.admDoseKind(d),
    amt   = paste0("dose ", format(d)),
    ev    = "event table",
    bands = paste0("dose by ", d$covariate, " band"),
    fun   = paste0("dose = f(", paste(.admDoseCovs(d), collapse = ", "), ")"))
}
