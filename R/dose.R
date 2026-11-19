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
#' `ii`, `addl`, `dur`, `rate`, `cmt`, ... --- is written as in [rxode2::et()].
#' Each is one value, shared by every band, or one value per band, so the
#' regimen itself can depend on the covariate (`ii = c(24, 12)`).
#'
#' Several doses in a band (an oral dose plus an infusion, or an IV phase
#' followed by an oral phase): give `amt` as a list and the per-dose details as
#' vectors inside lists (a length-1 list is shared by all bands), e.g.
#' `amt = list(c(500, 100), c(1000, 200)), cmt = list(c("depot", "centr")),
#' dur = list(c(0, 1))`.
#'
#' **A band's dose may itself depend on a covariate**, e.g. "CrCl below 30:
#' 5 mg/kg, otherwise 10 mg/kg":
#' * `per = WT` makes every amount a dose per unit of `WT`:
#'   `doseBands(CRCL, breaks = 30, amt = c(5, 10), per = WT)`;
#' * for any other formula, give a band a function of the covariates it reads,
#'   returning that band's amount(s): `amt = list(function(WT) 5 * WT,
#'   function(WT) pmin(10 * WT, 1000))`. Bands may mix numbers and functions.
#'
#' The band edges stay exact either way; only the covariate the dose is
#' proportional to is cut into `strata_nodes` nodes inside each band. Write a
#' step in the banded covariate this way rather than as `ifelse()` in a dose
#' function, which nodes can only approximate.
#'
#' The study is expanded into one study per band (and per node, where a band's
#' dose reads a covariate): the covariate restricted to the band (and to the
#' study's enrolled `range`, if given), `n` split by the band's share of the
#' population, and that band's dose, so the analysis model needs no dose rule.
#' Covariates correlated with the banded one are conditioned on the band.
#'
#' @param covariate The covariate the bands are cut on, bare or quoted
#'   (`CRCL` or `"CRCL"`); a variable holding the name also works. It must be declared, as a continuous covariate, in
#'   the study's `population`.
#' @param breaks Increasing cut points between bands.
#' @param amt Amount for each band, lowest band first; one more than `breaks`.
#'   A band that gets several doses at once (an oral dose plus an infusion) takes
#'   a list, one element per band, each holding that band's amounts. A list
#'   element may instead be a function of covariates returning the band's
#'   amount(s), e.g. `function(WT) 5 * WT`.
#' @param ... How each dose is given, as [rxode2::et()] arguments: `time`,
#'   `cmt`, `ii`, `addl`, `dur`, `rate`, `ss`, `until`; length 1 (all bands) or
#'   one per band. Without them each band gets one bolus at time 0.
#' @param per Optional covariate every amount is given per unit of, bare or
#'   quoted: `per = WT` reads `amt = c(5, 10)` as 5 and 10 mg/kg. It must be
#'   declared in the study's `population`. For a cap or any other formula use
#'   functions in `amt` instead.
#' @return An `admDoseBands` object, to pass as `dose` to [admStudy()].
#' @seealso [admStudy()]; a function of covariates as `dose` for dosing that is
#'   a smooth formula of covariates with no step, e.g. mg/kg for everyone.
#' @examples
#' doseBands(CRCL, breaks = 30, amt = c(100, 200))
#' # a 1-h infusion every 12 h, seven doses, in every band
#' doseBands(CRCL, breaks = 30, amt = c(500, 1000), ii = 12, addl = 6, dur = 1)
#' # the interval depends on the band too: daily below 30, twice daily above
#' doseBands(CRCL, breaks = 30, amt = c(500, 1000), ii = c(24, 12), addl = c(2, 6))
#' # every regimen detail depends on the band: dose, interval, duration, route
#' doseBands(CRCL, breaks = c(25, 40, 60), amt = c(450, 450, 450, 900),
#'           ii = c(48, 24, 12, 12), addl = 16, dur = c(0, 0, 1, 1),
#'           cmt = c("depot", "depot", "centr", "centr"))
#' # an oral dose AND a 1-h infusion at the same time, in every band
#' doseBands(CRCL, breaks = 30, amt = list(c(500, 100), c(1000, 200)),
#'           cmt = list(c("depot", "centr")), dur = list(c(0, 1)))
#' # 5 days of IV (daily 1-h infusions) followed by 16 days oral, per band
#' doseBands(CRCL, breaks = 30, amt = list(c(500, 250), c(1000, 500)),
#'           time = list(c(0, 120)), cmt = list(c("centr", "depot")),
#'           dur = list(c(1, 0)), ii = 24, addl = list(c(4, 15)))
#' # mg/kg, with the mg/kg depending on the band
#' doseBands(CRCL, breaks = 30, amt = c(5, 10), per = WT)
#' # any formula per band: 5 mg/kg below 30, 10 mg/kg capped at 1 g above
#' doseBands(CRCL, breaks = 30,
#'           amt = list(function(WT) 5 * WT, function(WT) pmin(10 * WT, 1000)))
#' # a formula with no step needs no bands: give `dose` a function, and any
#' # et() argument may read the covariates
#' f <- function(WT, CRCL) {
#'   rxode2::et(amt = 5 * WT, cmt = "depot") |>
#'     rxode2::et(amt = 2 * WT, cmt = "centr", dur = 1, ii = 12 + CRCL / 10, addl = 3)
#' }
#' # admStudy(model, est, n, times, population, dose = f)
#' @export
doseBands <- function(covariate, breaks, amt, ..., per = NULL) {
  # a bare name is the covariate (doseBands(CRCL, ...)), unless it is a
  # variable holding a covariate name (cv <- "CRCL"; doseBands(cv, ...))
  cv  <- .admBareName(substitute(covariate), parent.frame())
  per <- .admBareName(substitute(per), parent.frame())
  bad <- function(...) stop("admixr2: doseBands(): ", ..., call. = FALSE)
  if (!is.character(cv) || length(cv) != 1L || is.na(cv) || !nzchar(cv))
    bad("`covariate` must name one covariate, e.g. doseBands(CRCL, ...).")
  if (!is.null(per) &&
      (!is.character(per) || length(per) != 1L || is.na(per) || !nzchar(per)))
    bad("`per` must name one covariate, e.g. per = WT for a dose in mg/kg.")
  if (!is.numeric(breaks) || !length(breaks) || any(!is.finite(breaks)) ||
      is.unsorted(breaks, strictly = TRUE))
    bad("`breaks` must be finite, strictly increasing cut points, e.g. ",
        "breaks = 30 or breaks = c(30, 60).")
  if (length(amt) != length(breaks) + 1L)
    bad("`amt` needs one amount per band: ", length(breaks), " cut point",
        if (length(breaks) > 1L) "s make " else " makes ",
        length(breaks) + 1L, " bands, but ", length(amt),
        if (length(amt) == 1L) " amount was" else " amounts were", " given.")
  # a list holds several simultaneous doses per band, e.g. oral plus infusion,
  # or a band's function of covariates
  isfn <- if (is.list(amt)) vapply(amt, is.function, logical(1))
          else rep(FALSE, length(amt))
  av <- if (is.list(amt)) unlist(amt[!isfn]) else amt
  if (length(av) && (!is.numeric(av) || any(!is.finite(av)) || any(av < 0)))
    bad("each `amt` must be a finite, non-negative amount",
        if (is.list(amt)) ", or a function of covariates returning one", ".")
  if (any(isfn) && !is.null(per))
    bad("`per` multiplies numeric amounts; with a function in `amt`, write ",
        "the covariate into the function instead, e.g. function(", per, ") 5 * ",
        per, ".")
  for (j in which(isfn)) {
    fa <- names(formals(amt[[j]]))
    if (!length(fa) || "..." %in% fa)
      bad("the function for band ", j, " must take the covariates it reads as ",
          "named arguments, e.g. function(WT) 5 * WT",
          if (!length(fa)) "; a dose that reads no covariate is just a number.",
          ".")
  }
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
    k <- length(amt)
    badlen <- tn[!vapply(timing, function(v) length(v) %in% c(1L, k), logical(1))]
    if (length(badlen))
      bad(paste(sQuote(badlen), collapse = ", "), " must have length 1 (the ",
          "same in every band) or ", k, " (one per band).")
    # a band's function is evaluated at its nodes later; a stand-in amount
    # checks the timing builds an event table now
    probe <- if (is.list(amt)) lapply(amt, function(a) if (is.function(a)) 1 else a)
             else amt
    for (j in seq_len(k))
      tryCatch(.admDoseBandEv(list(amt = probe, timing = timing), j),
               error = function(e)
        bad("these dosing details do not build an event table: ",
            conditionMessage(e)))
  }
  structure(list(covariate = cv, breaks = as.numeric(breaks),
                 amt = if (is.list(amt))
                         lapply(amt, function(a) if (is.function(a)) a else as.numeric(a))
                       else as.numeric(amt),
                 timing = timing, per = per),
            class = "admDoseBands")
}

# A covariate named in an argument, bare (`CRCL`) or quoted (`"CRCL"`), or a
# variable holding the name. NULL stays NULL.
.admBareName <- function(sx, env) {
  if (is.null(sx)) return(NULL)
  if (is.name(sx)) {
    v <- tryCatch(eval(sx, env), error = function(e) NULL)
    return(if (is.character(v) && length(v) == 1L) v else as.character(sx))
  }
  eval(sx, env)
}

#' @export
format.admDoseBands <- function(x, ...) {
  k <- length(x$amt)
  a <- vapply(x$amt, function(v) {
    if (is.function(v))
      return(paste0("f(", paste(names(formals(v)), collapse = ", "), ")"))
    s <- vapply(v, format, "")
    if (!is.null(x$per)) s <- paste0(s, " x ", x$per)
    paste(s, collapse = "+")
  }, "")
  lab <- paste0(vapply(seq_len(k), function(j) .admDoseBandWhere(x, j), ""),
                ": ", a)
  out <- paste(lab, collapse = "; ")
  if (length(x$timing))
    out <- paste0(out, "  (", paste(sprintf("%s = %s", names(x$timing),
      vapply(x$timing, function(v) paste(format(v, trim = TRUE), collapse = ","),
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

# The covariates a dose rule reads: the banded one plus whatever a band's dose
# is proportional to (`per`) or a function of; or a function's arguments.
.admDoseCovs <- function(d) {
  if (inherits(d, "admDoseBands"))
    return(unique(c(d$covariate, unlist(lapply(seq_along(d$amt),
                                               .admDoseBandAmtCovs, d = d)))))
  if (is.function(d)) return(names(formals(d)))
  character(0)
}

# The bands as a lo/hi/amt table (open at both ends).
.admDoseBandTable <- function(d)
  data.frame(lo = c(-Inf, d$breaks), hi = c(d$breaks, Inf))

# Band j as a reader would say it: "CRCL < 30", "30 to < 60", "CRCL >= 60".
.admDoseBandWhere <- function(d, j) {
  b <- d$breaks; k <- length(b) + 1L; cv <- d$covariate
  if (j == 1L) sprintf("%s < %g", cv, b[1L])
  else if (j == k) sprintf("%s >= %g", cv, b[k - 1L])
  else sprintf("%g to < %g", b[j - 1L], b[j])
}

# The covariates band j's AMOUNT reads (not the banded one, which only picks the
# band): its function's arguments, or `per`. Empty for a plain amount.
.admDoseBandAmtCovs <- function(d, j) {
  a <- d$amt[[j]]
  if (is.function(a)) names(formals(a))
  else if (!is.null(d$per)) d$per
  else character(0)
}

# Band j's amount(s) at one node `cov` (a named list of covariate values), checked.
.admDoseBandAmt <- function(d, j, cov = NULL) {
  a  <- d$amt[[j]]
  cv <- .admDoseBandAmtCovs(d, j)
  if (!length(cv)) return(a)
  if (is.null(cov) || !all(cv %in% names(cov)))
    stop("admixr2: internal: band ", j, "'s dose reads ",
         paste(cv, collapse = ", "), ", but no value was given.", call. = FALSE)
  at <- paste(sprintf("%s = %s", cv, vapply(cov[cv], function(v)
    format(signif(as.numeric(v), 4)), character(1))), collapse = ", ")
  where <- paste0("in band ", sQuote(.admDoseBandWhere(d, j)), " at ", at)
  r <- if (is.function(a))
         tryCatch(do.call(a, as.list(cov)[cv]), error = function(e)
           stop("admixr2: doseBands(): the dose function failed ", where, ": ",
                conditionMessage(e), call. = FALSE))
       else a * as.numeric(cov[[d$per]])
  if (inherits(r, "rxEt"))
    stop("admixr2: doseBands(): a band's function returns the band's AMOUNT(S), ",
         "not an event table; give the timing (ii, addl, dur, ...) to doseBands() ",
         "itself, e.g. doseBands(", d$covariate, ", ..., ii = 12, addl = 6).",
         call. = FALSE)
  if (!is.numeric(r) || !length(r) || any(!is.finite(r)) || any(r < 0))
    stop("admixr2: doseBands(): the dose ", where, " is ",
         if (is.numeric(r)) paste(format(r), collapse = ", ") else class(r)[1L],
         "; every amount must be a finite, non-negative number. A covariate ",
         "whose margin reaches zero or below (a NORMAL weight) can be bounded ",
         "with `range`, or declared with dist = \"lnorm\".", call. = FALSE)
  as.numeric(r)
}

# Band j's value of one dosing argument: length-1 is shared, length-k (or a
# list of k) is per band. A band's value may itself hold several doses.
.admDoseAt <- function(v, j, k) {
  if (is.list(v)) v[[if (length(v) == 1L) 1L else j]] else if (length(v) == k && k > 1L) v[[j]] else v
}

# One band's event table, at node `cov` when its amount reads a covariate.
# When a band's arguments hold several values they are several simultaneous
# doses (recycled to the longest), e.g. amt = list(c(500, 100), 1000),
# cmt = c("depot", "centr") for an oral dose plus an infusion.
.admDoseBandEv <- function(d, j, cov = NULL) {
  k <- length(d$amt)
  a <- c(list(amt = .admDoseBandAmt(d, j, cov)),
         lapply(d$timing, .admDoseAt, j = j, k = k))
  m <- max(lengths(a))
  rows <- lapply(seq_len(m), function(i)
    lapply(a, function(v) v[[(i - 1L) %% length(v) + 1L]]))
  Reduce(function(ev, r) do.call(rxode2::et, c(list(ev), r)),
         rows[-1L], do.call(rxode2::et, rows[[1L]]))
}

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

# The amount a dose record set gives in total, for the smoothness check.
.admDoseTotal <- function(ev) {
  d <- as.data.frame(ev)
  sum(d$amt[!is.na(d$evid) & !(d$evid %in% c(0L, 2L, 3L))])
}

# Band j's span [lo, hi] inside the enrolled range `rng` on the banded
# covariate; NULL when the band lies outside it.
.admDoseBandSpan <- function(d, j, rng = NULL) {
  bt  <- .admDoseBandTable(d)
  rng <- sort(as.numeric(rng %||% c(-Inf, Inf)))
  lo  <- max(bt$lo[j], rng[1L]); hi <- min(bt$hi[j], rng[2L])
  if (lo >= hi) NULL else c(lo, hi)
}

# The covariates band j's nodes are cut over, given the canonical population
# `cd`: what its amount reads, plus the banded covariate itself when the amount
# reads it or it is latently correlated with what the amount reads -- the nodes
# then condition the amount's covariate on the band. NULL for a plain amount.
.admDoseBandNodeSet <- function(cd, d, j) {
  acv <- .admDoseBandAmtCovs(d, j)
  if (!length(acv)) return(NULL)
  bcv <- d$covariate
  nms <- .admCovSpecNames(cd)
  R   <- cd[["latentR"]]
  oth <- match(setdiff(acv, bcv), nms)
  corr <- !is.null(R) && length(oth) &&
    any(abs(R[match(bcv, nms), oth]) > 1e-12)
  unique(c(if (bcv %in% acv || corr) bcv, acv))
}

# Nodes over the dosing covariates S, inside `cov_range`, as .admCovStrata()
# returns them -- but CONDITIONED SEQUENTIALLY when a range-truncated covariate
# is latently correlated with another one in S.
#
# .admCovStrata() lays a joint grid over S rotated by chol(R_SS), and a
# truncated margin enters it through its OWN uniform, so the correlation acts on
# the band's score, not on the full-margin latent the copula couples through:
# with CRCL < 30 and rho = 0.6, the WT nodes came out at E[WT] = 75 in both
# bands against 63.9 / 76.6 from the full law. The marginalised block maps
# through `.band` and is exact, so cut the truncated covariate FIRST (the rest
# conditioned exactly on each node), then cut the rest from each node's
# conditional law (the incoming latentMu/latentScale path). A product of exact
# conditional rules is a rule for the joint law. Otherwise the joint grid is
# kept, so every design that was right stays bit-identical.
.admDoseStrata <- function(cd, S, n_nodes, cov_range = NULL) {
  cd  <- .admCovDistCanon(cd)
  nms <- .admCovSpecNames(cd)
  R   <- cd[["latentR"]]
  cont <- S[vapply(S, function(v) is.null(cd[[v]][["values"]]), logical(1))]
  rngd <- intersect(names(cov_range), cont)
  coupled <- function(v) {
    o <- match(setdiff(cont, v), nms)
    !is.null(R) && length(o) && any(abs(R[match(v, nms), o]) > 1e-12)
  }
  first <- Filter(coupled, rngd)
  if (length(S) < 2L || !length(first))
    return(.admCovStrata(cd, S, n_nodes, warn_range = FALSE,
                         cov_range = cov_range))
  first <- first[1L]
  rest  <- setdiff(S, first)
  st1 <- .admCovStrata(cd, first, n_nodes, warn_range = FALSE,
                       cov_range = cov_range)
  # the remaining margins come back already truncated, so the ranges are spent
  out <- list()
  for (nd in st1) {
    st2 <- .admDoseStrata(nd$cov_dist, rest, n_nodes, NULL)
    for (nd2 in st2)
      out[[length(out) + 1L]] <- list(cov = c(nd$cov, nd2$cov)[S],
                                      cov_dist = nd2$cov_dist,
                                      weight = nd$weight * nd2$weight,
                                      n_pool_cell = NA_integer_)
  }
  out
}

# One spec per node of `stl` (from .admCovStrata over covariates S), built from
# the group spec `sp`: the node's share of `sp$n`, its covariates pinned, the
# rest kept CONDITIONAL on it. S's ranges were spent cutting the nodes and are
# dropped -- re-applied to the POINT spec datagen() dies in approx() -- and S
# is taken out of `stratify`, or datagen() would cut a pinned covariate again.
# Re-applying the OTHER ranges is a no-op (F(a) = 0, F(b) = 1 already).
.admDoseNodeSpecs <- function(sp, stl, S) {
  lapply(stl, function(nd) {
    spk <- sp
    spk$n <- sp$n * nd$weight
    spk[["cov"]] <- utils::modifyList(as.list(sp[["cov"]] %||% list()), nd$cov)
    spk$cov_dist <- nd$cov_dist
    for (.c in S) spk$cov_range[[.c]] <- NULL
    if (!length(spk$cov_range)) spk$cov_range <- NULL
    if (!identical(spk$stratify, FALSE)) {
      spk$stratify <- setdiff(spk$stratify, S)
      if (!length(spk$stratify)) spk$stratify <- FALSE
    }
    spk
  })
}

# WHERE A DOSE RULE IS NOT SMOOTH. Quadrature nodes cannot place a jump (a step,
# "CrCl < 30: half dose") or a sharp bend (a cap, pmin(15 * WT, 2000)), so dose
# groups cut from nodes only approximate either -- and the error does not go
# away smoothly as nodes are added (measured for a CrCl step: 0.125, 0.124,
# 0.139 at 9, 15, 25 nodes against an exact share of 0.129).
#
# Probed along each continuous covariate the rule reads, on a fine grid over
# its central 99% (inside `cov_range`), with the others at their quartiles
# (discrete ones at every level). Checking at the NODES instead only worked for
# one covariate: with two, the nodes do not lie on a line, so a step in CRCL
# under a dose of WT and CRCL passed without a word.
#
# `f` is the amount rule, called with the covariates in `covs` as named
# arguments; `total(cov)` is the total amount at one named list of values.
# Returns NULL when smooth, else list(covariate, kind = "step"/"kink", at).
#
# Most rules are vectorised (2 * WT, ifelse(), pmin()), so every probe point
# goes through ONE call to `f` first. Building an et() per point instead cost
# 3-6 s at construction; a rule that does not vectorise (it returns an et(),
# or several doses) falls back to point-by-point on a coarser probe.
.admDoseRough <- function(f, total, cd, covs, cov_range = NULL) {
  cd <- .admCovDistCanon(cd)
  if (length(cov_range))
    cd <- .admCovApplyRange(cd, cov_range[intersect(names(cov_range),
                                                    .admCovSpecNames(cd))])
  isdisc <- function(cv) !is.null(cd[[cv]][["values"]])
  levels_of <- function(cv, mid_only = FALSE)
    if (isdisc(cv)) as.numeric(cd[[cv]][["values"]])
    else .admCovQuantile(cd[[cv]], if (mid_only) 0.5 else c(0.25, 0.5, 0.75))
  # every probe line: one continuous covariate on a grid, the others fixed
  build <- function(np, mid_only) {
    lines <- list()
    for (cv in covs[!vapply(covs, isdisc, logical(1))]) {
      lim <- .admCovQuantile(cd[[cv]], c(0.005, 0.995))
      if (!all(is.finite(lim)) || lim[1L] >= lim[2L]) next
      x <- seq(lim[1L], lim[2L], length.out = np)
      oth <- setdiff(covs, cv)
      grid <- if (length(oth))
        expand.grid(lapply(stats::setNames(oth, oth), levels_of,
                           mid_only = mid_only), KEEP.OUT.ATTRS = FALSE)
      if (!is.null(grid) && nrow(grid) > 9L)
        grid <- expand.grid(lapply(stats::setNames(oth, oth), levels_of,
                                   mid_only = TRUE), KEEP.OUT.ATTRS = FALSE)
      for (r in seq_len(if (is.null(grid)) 1L else nrow(grid))) {
        pts <- stats::setNames(list(x), cv)
        for (o in oth) pts[[o]] <- rep(grid[[o]][r], np)
        lines[[length(lines) + 1L]] <- list(cv = cv, x = x, pts = pts[covs])
      }
    }
    lines
  }
  # a_of(i) gives the amounts along line i
  check <- function(lines, a_of) {
    for (i in seq_along(lines)) {
      ln  <- lines[[i]]
      hit <- .admDoseRoughAlong(ln$x, a_of(i))
      if (is.null(hit)) next
      # a step's grid interval holds the cut point; bisect it, so the advice
      # reads breaks = 30 rather than the grid's 29.6
      if (hit$kind == "step") {
        at_x <- function(xv) {
          p <- lapply(ln$pts, `[`, 1L); p[[ln$cv]] <- xv
          tryCatch(total(p), error = function(e) NA_real_)
        }
        lo <- ln$x[hit$j]; hi <- ln$x[hit$j + 1L]
        alo <- at_x(lo); ahi <- at_x(hi)
        if (is.finite(alo) && is.finite(ahi)) {
          for (it in seq_len(40L)) {
            mid <- (lo + hi) / 2; am <- at_x(mid)
            if (!is.finite(am)) break
            if (abs(am - alo) <= abs(am - ahi)) lo <- mid else hi <- mid
          }
          hit$at <- (lo + hi) / 2
        }
      }
      hit$j <- NULL
      return(c(list(covariate = ln$cv), hit))
    }
    NULL
  }
  lines <- build(101L, FALSE)
  if (!length(lines)) return(NULL)
  allp <- lapply(stats::setNames(covs, covs), function(v)
    unlist(lapply(lines, function(ln) ln$pts[[v]]), use.names = FALSE))
  va <- if (is.function(f))
    tryCatch(suppressWarnings(do.call(f, allp)), error = function(e) NULL)
  if (is.numeric(va) && !inherits(va, "rxEt") &&
      length(va) == length(allp[[1L]])) {
    off <- cumsum(c(0L, vapply(lines, function(ln) length(ln$x), 1L)))
    return(check(lines, function(i) va[(off[i] + 1L):off[i + 1L]]))
  }
  lines <- build(41L, TRUE)
  check(lines, function(i) {
    ln <- lines[[i]]
    vapply(seq_along(ln$x), function(k)
      tryCatch(total(lapply(ln$pts, `[`, k)), error = function(e) NA_real_),
      numeric(1))
  })
}

# One probe line: a JUMP is one grid step that moves the dose far more than its
# neighbours do; a KINK is a slope that changes by a quarter of its largest
# value between two neighbouring steps. A smooth rule (mg/kg, allometric) drifts
# slowly on a 101-point grid and trips neither.
.admDoseRoughAlong <- function(x, a) {
  if (length(a) < 4L || any(!is.finite(a))) return(NULL)
  da <- diff(a); span <- diff(range(a))
  if (span <= 0) return(NULL)
  ad <- abs(da)
  nb <- pmax(c(0, ad[-length(ad)]), c(ad[-1L], 0))
  j  <- which(ad > 10 * nb & ad > 1e-3 * span)
  if (length(j)) return(list(kind = "step", at = mean(x[j[1L] + 0:1]), j = j[1L]))
  s  <- da / diff(x); ds <- abs(diff(s)); m <- max(abs(s))
  if (m > 0 && max(ds) > 0.25 * m)
    return(list(kind = "kink", at = x[which.max(ds) + 1L]))
  NULL
}

# The warning for a rough dose rule, saying what was found, where, and what to
# write instead. `what` names the rule ("the dose function", "the dose in band
# 'CRCL >= 30'"); `nodes` is the strata_nodes in force.
.admDoseRoughWarn <- function(hit, what, nodes, band_of = NULL) {
  at <- format(signif(hit$at, 3))
  msg <- if (hit$kind == "step")
    paste0("admixr2: ", what, " is not smooth: it jumps where ", hit$covariate,
           " is about ", at, " (a step?). Doses are evaluated at quadrature ",
           "nodes, which cannot place a jump, so the share of patients on ",
           "each side is only approximate, and adding nodes does not fix it. ",
           if (is.null(band_of))
             paste0("Write the step as a dosing table so the split is exact: ",
                    "dose = doseBands(", hit$covariate, ", breaks = ", at,
                    ", amt = ...), with `per = <covariate>` for a mg/kg dose or ",
                    "a function per band for any other formula.")
           else paste0("Add a cut point there: doseBands(", band_of,
                       ", breaks = c(..., ", at, ", ...), ...), giving each side ",
                       "its own amount, or band on ", hit$covariate, " instead."))
  else
    paste0("admixr2: ", what, " is not smooth: it bends sharply where ",
           hit$covariate, " is about ", at, " (a cap?). Doses are evaluated at ",
           nodes, " quadrature nodes, which only approximate a bend; raise ",
           "`strata_nodes` to resolve it better.")
  warning(msg, call. = FALSE)
}

# One line describing a dose rule, for print.admStudy().
.admDoseLabel <- function(d) {
  switch(.admDoseKind(d),
    amt   = paste0("dose ", format(d)),
    ev    = "event table",
    bands = paste0("dose by ", d$covariate, " band"),
    fun   = paste0("dose = f(", paste(.admDoseCovs(d), collapse = ", "), ")"))
}
