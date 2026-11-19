#' A dosing table by covariate band
#'
#' Transcribes a paper's dosing table, e.g. "CrCl below 30 mL/min: 100 mg;
#' otherwise 200 mg", for the `dose` argument of [admStudy()]:
#' `dose = doseBands(CRCL, breaks = 30, amt = c(100, 200))`.
#'
#' `breaks` are the cut points between bands and `amt` the amount for each band,
#' lowest first, so there is always one more amount than cut points and the
#' bands cover the whole range with no gaps. A value equal to a cut point falls
#' in the band above it.
#'
#' **Everything else is [rxode2::et()].** Each band's doses are built by
#' `rxode2::et()` itself, and every `et()` argument may be given in `...`
#' (`ii`, `addl`, `until`, `dur`, `rate`, `cmt`, `ss`, `evid`, `time`, units,
#' ...): one value shared by every band, or one value per band, so the regimen
#' itself can depend on the covariate (`ii = c(24, 12)`). The accepted names are
#' read from rxode2 at run time, so an argument rxode2 adds is accepted with no
#' change here; a misspelt one is refused, since `et()` would silently ignore
#' it. A band may also be a whole `et()` table, `amt = list(rxode2::et(...),
#' rxode2::et(...))`, for anything one call cannot express.
#'
#' Without `cmt` a dose goes to the model's default (first) compartment --- the
#' depot in an oral model --- so give an IV dose its compartment
#' (`cmt = "centr"`); [admStudy()] warns about an infusion that does not.
#'
#' **Phases** (5 days IV, then 16 days oral) are one `doseBands()` per phase,
#' combined exactly as rxode2 combines event tables: `c(iv, oral)` (or
#' `seq(iv, oral)`) runs one after the other per band, through
#' [rxode2::etSeq()], so the oral phase starts where the IV phase ends;
#' `rbind(a, b)` merges them by time, through [rxode2::etRbind()]. The phases
#' may use different cut points (each band takes the doses of the phase band it
#' falls in) but must band on the same covariate. An `et()` table among them is
#' a phase every band gets, e.g. a loading dose. Put a `doseBands()` first:
#' with an `et()` table first, rxode2 handles the call and drops the bands
#' (see [combine-doseBands]).
#'
#' **A band's dose may itself depend on a covariate**, e.g. "CrCl below 30:
#' 5 mg/kg, otherwise 10 mg/kg":
#' * `per = WT` makes every amount a dose per unit of `WT`:
#'   `doseBands(CRCL, breaks = 30, amt = c(5, 10), per = WT)`;
#' * for any other formula, give a band a function of the covariates it reads,
#'   returning that band's amount(s) or a whole `et()` table:
#'   `amt = list(function(WT) 5 * WT, function(WT) pmin(10 * WT, 1000))`.
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
#'   (`CRCL` or `"CRCL"`); a variable holding the name also works. It must be
#'   declared, as a continuous covariate, in the study's `population`.
#' @param breaks Increasing cut points between bands.
#' @param amt Amount for each band, lowest band first; one more than `breaks`.
#'   As a list, one element per band, each element is that band's amounts
#'   (several doses at once, e.g. an oral dose plus an infusion), a function of
#'   covariates returning them or an `et()` table, or an `rxode2::et()` table.
#' @param ... [rxode2::et()] arguments for every band's doses; length 1 (all
#'   bands) or one per band. A list gives one element per band, each holding
#'   that band's per-dose values. Not used with a band given as an `et()`
#'   table, which carries its own.
#' @param per Optional covariate every amount is given per unit of, bare or
#'   quoted: `per = WT` reads `amt = c(5, 10)` as 5 and 10 mg/kg. It must be
#'   declared in the study's `population`. For a cap or any other formula use
#'   functions in `amt` instead.
#' @return An `admDoseBands` object, to pass as `dose` to [admStudy()]. It
#'   prints as the dosing table: each band's event table as `rxode2::et()`
#'   built it.
#' @seealso [admStudy()]; a function of covariates as `dose` for dosing that is
#'   a smooth formula of covariates with no step, e.g. mg/kg for everyone.
#' @examples
#' doseBands(CRCL, breaks = 30, amt = c(100, 200))
#' # a 1-h infusion every 12 h, seven doses, in every band
#' doseBands(CRCL, breaks = 30, amt = c(500, 1000), ii = 12, addl = 6, dur = 1)
#' # the interval depends on the band too: daily below 30, twice daily above
#' doseBands(CRCL, breaks = 30, amt = c(500, 1000), ii = c(24, 12), addl = c(2, 6))
#' # a paper's table: mg/kg as a 1-h IV infusion, interval by CrCl, for 7 days
#' doseBands(CRCL, breaks = c(25, 50, 70), amt = c(1.25, 2.5, 2.5, 5),
#'           per = WT, cmt = "centr", dur = 1, ii = c(24, 24, 12, 12),
#'           addl = c(6, 6, 13, 13))
#' # 5 days IV, then 16 days oral: one doseBands() per phase, in sequence
#' # (name the IV compartment: without `cmt` a dose goes to the model's FIRST
#' # compartment, which in an oral model is the depot)
#' iv   <- doseBands(CRCL, breaks = 30, amt = c(250, 500), cmt = "centr",
#'                   dur = 1, ii = 24, addl = 4)
#' oral <- doseBands(CRCL, breaks = 30, amt = c(250, 500), cmt = "depot",
#'                   ii = 24, addl = 15)
#' c(iv, oral)
#' # an oral dose AND a 1-h infusion at the same time: merge by time
#' rbind(doseBands(CRCL, breaks = 30, amt = c(500, 1000), cmt = "depot"),
#'       doseBands(CRCL, breaks = 30, amt = c(100, 200), cmt = "centr", dur = 1))
#' # ... or in one call, two doses per band
#' doseBands(CRCL, breaks = 30, amt = list(c(500, 100), c(1000, 200)),
#'           cmt = list(c("depot", "centr")), dur = list(c(0, 1)))
#' # any et() argument passes through, e.g. a reset-and-dose record
#' doseBands(CRCL, breaks = 30, amt = c(100, 200), evid = 4)
#' # a band given as its own et() table
#' doseBands(CRCL, breaks = 30, amt = list(
#'   rxode2::et(amt = 100, ii = 24, addl = 6),
#'   rxode2::et(amt = 400, dur = 1) |> rxode2::et(amt = 200, time = 12, ii = 12, addl = 12)))
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
#' @importFrom rxode2 et
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
  if (inherits(amt, "rxEt") || is.function(amt) || length(amt) != length(breaks) + 1L)
    bad("`amt` needs one amount per band: ", length(breaks), " cut point",
        if (length(breaks) > 1L) "s make " else " makes ",
        length(breaks) + 1L, " bands, but ",
        if (inherits(amt, "rxEt") || is.function(amt)) "one rule was given; put one per band in a list."
        else paste0(length(amt), if (length(amt) == 1L) " amount was" else " amounts were",
                    " given."))
  # each band: amount(s), a function of covariates, or a whole et() table
  kind <- if (is.list(amt)) vapply(amt, function(a)
    if (is.function(a)) "fun" else if (inherits(a, "rxEt")) "et" else "num", "")
    else rep("num", length(amt))
  av <- if (is.list(amt)) unlist(amt[kind == "num"]) else amt
  if (length(av) && (!is.numeric(av) || any(!is.finite(av)) || any(av < 0)))
    bad("each `amt` must be a finite, non-negative amount",
        if (is.list(amt)) ", a function of covariates, or an rxode2::et() table", ".")
  if (any(kind != "num") && !is.null(per))
    bad("`per` multiplies numeric amounts; with a function or an et() table in ",
        "`amt`, write the covariate into it instead, e.g. function(", per,
        ") 5 * ", per, ".")
  for (j in which(kind == "fun")) {
    fa <- names(formals(amt[[j]]))
    if (!length(fa) || "..." %in% fa)
      bad("the function for band ", j, " must take the covariates it reads as ",
          "named arguments, e.g. function(WT) 5 * WT",
          if (!length(fa)) "; a dose that reads no covariate is just a number.",
          ".")
  }
  timing <- list(...)
  tn <- names(timing)
  if (length(timing) && (is.null(tn) || any(!nzchar(tn))))
    bad("the dosing details in `...` must be named, e.g. ii = 12, addl = 6.")
  # rxode2::et() IGNORES arguments it does not know, so a typo (`adl = 6`)
  # would silently turn a regimen into a single dose. The names are rxode2's own.
  ok  <- .admEtArgs()
  unk <- if (!is.null(ok)) setdiff(tn, ok)
  if (length(unk))
    bad(paste(sQuote(unk), collapse = ", "), " is not an rxode2::et() argument; ",
        "et() takes ", paste(ok, collapse = ", "), ".")
  if (length(timing) && any(kind == "et"))
    bad("band ", which(kind == "et")[1L], " is a whole et() table, which carries ",
        "its own timing; give the details for the other bands inside et() too, ",
        "or use no et() table.")
  k <- length(amt)
  badlen <- tn[!vapply(timing, function(v) length(v) %in% c(1L, k), logical(1))]
  if (length(badlen))
    bad(paste(sQuote(badlen), collapse = ", "), " must have length 1 (the ",
        "same in every band) or ", k, " (one per band).")
  d <- .admDoseNew(covariate = cv, breaks = as.numeric(breaks),
                   amt = if (is.list(amt)) lapply(amt, function(a)
                           if (is.numeric(a)) as.numeric(a) else a)
                         else as.numeric(amt),
                   timing = timing, per = per)
  # build every band once now, through et() itself, with stand-in amounts where
  # they need covariates: a detail et() rejects fails here, not at fit time
  for (j in seq_len(k))
    tryCatch(.admDoseBandEv(d, j, labels = .admDoseLabels()),
             error = function(e)
      bad("band ", sQuote(.admDoseBandWhere(d, j)), " does not build an ",
          "event table: ", sub("^admixr2: (doseBands\\(\\): )?", "",
                               conditionMessage(e))))
  d
}

# The rxode2::et() argument names, read from rxode2's own method at run time,
# minus the event table itself (`x`), `...`, `amt` (given as doseBands()'s own)
# and `envir`. NULL when rxode2 cannot be asked, and then nothing is refused.
.admEtArgs <- function() {
  m <- tryCatch(utils::getS3method("et", "default"), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  setdiff(names(formals(m)), c("x", "...", "amt", "envir"))
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

#' Combine dosing phases, as rxode2 combines event tables
#'
#' `c(iv, oral)` and `seq(iv, oral)` give each band the first phase and then
#' the next, through [rxode2::etSeq()], so a phase starts where the previous one
#' ends --- exactly what `c()` does for two `rxode2::et()` tables. `rbind(a, b)`
#' merges them by time, through [rxode2::etRbind()], e.g. an oral dose and an
#' infusion given together. The phases may use different cut points (every
#' band of the result gets the doses of the phase band it falls in) but must
#' band on the same covariate. An `rxode2::et()` table among them is a phase
#' every band gets, e.g. a loading dose; a number in `c()`/`seq()` is what it is
#' to `etSeq()`, a wait in hours between the phases (`c(iv, 12, oral)`).
#'
#' **Put a `doseBands()` FIRST.** R sends `c()`, `seq()` and `rbind()` to the
#' method of their first argument, so `c(rxode2::et(...), bands)` or
#' `rbind(rxode2::et(...), bands)` is handled by rxode2, which keeps only the
#' event table and silently drops the bands --- the study would then be fitted
#' with the loading dose alone, and admixr2 never sees the bands to object.
#' Write `rbind(bands, rxode2::et(...))`.
#'
#' @param ... [doseBands()] objects, [rxode2::et()] tables, or (in `c()`/
#'   `seq()`) waits in hours; the FIRST must be a `doseBands()`.
#' @param deparse.level Not used; part of the `rbind()` generic.
#' @return An `admDoseBands` object, to pass as `dose` to [admStudy()].
#' @examples
#' iv   <- doseBands(CRCL, breaks = c(25, 50, 70), amt = c(1.25, 2.5, 2.5, 5),
#'                   per = WT, cmt = "centr", dur = 1, ii = c(24, 24, 12, 12),
#'                   addl = c(4, 4, 9, 9))
#' oral <- doseBands(CRCL, breaks = c(25, 50, 70), amt = c(5, 5, 10, 10),
#'                   per = WT, cmt = "depot", ii = 24, addl = 15)
#' c(iv, oral)
#' # a 1000 mg loading dose for everyone, then the banded maintenance dose from
#' # 12 h: merged by time (c() would start it one etSeq() gap after the load)
#' rbind(doseBands(CRCL, breaks = 30, amt = c(250, 500), dur = 1, time = 12,
#'                 ii = 12, addl = 12),
#'       rxode2::et(amt = 1000, dur = 1))
#' @name combine-doseBands
#' @aliases c.admDoseBands seq.admDoseBands rbind.admDoseBands
NULL

#' @rdname combine-doseBands
#' @export
c.admDoseBands <- function(...) .admDoseCombine(list(...), "seq")

#' @rdname combine-doseBands
#' @export
seq.admDoseBands <- function(...) .admDoseCombine(list(...), "seq")

#' @rdname combine-doseBands
#' @export
rbind.admDoseBands <- function(..., deparse.level = 1)
  .admDoseCombine(list(...), "rbind")

# `how` is "seq" (rxode2::etSeq) or "rbind" (rxode2::etRbind), applied per
# band. A combination is a tree: its `parts` are doseBands() objects, other
# combinations, fixed et() phases, or -- in a sequence -- a NUMBER, which is
# what it is to rxode2::etSeq(): a wait between the phases, not a dose.
.admDoseCombine <- function(xs, how) {
  parts <- lapply(xs, function(x) {
    if (inherits(x, "admDoseBands")) return(x)
    if (inherits(x, "rxEt")) return(.admDoseFixedPhase(x))
    if (is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 0) {
      if (how == "seq") return(as.numeric(x))
      stop("admixr2: rbind() merges event tables by time, and a number is not ",
           "one (rxode2::etRbind() would drop it). A dose for everyone is ",
           "rxode2::et(amt = ", x, "); a wait between phases goes in c() or ",
           "seq().", call. = FALSE)
    }
    stop("admixr2: only doseBands(), an rxode2::et() table",
         if (how == "seq") " or a wait (a number of hours)",
         " can be combined with doseBands(); got ", class(x)[1L], ".",
         call. = FALSE)
  })
  cvs <- unique(unlist(lapply(parts, .admDoseCovOf)))
  if (length(cvs) != 1L)
    stop("admixr2: dosing phases combined together must band on the SAME ",
         "covariate; these band on ", paste(sQuote(cvs), collapse = " and "),
         ". For a dose that depends on both, give one phase per ", cvs[1L],
         " band with a function of ", cvs[2L], " as its amount.", call. = FALSE)
  d <- .admDoseNew(covariate = cvs,
                   breaks = sort(unique(unlist(lapply(parts, .admDoseBreaksOf)))),
                   parts = parts, how = how)
  # rxode2::etSeq() spaces phases with no `ii` 24 h apart and warns about it on
  # EVERY build -- once per node and probe point, hundreds per study. Say it
  # once here; the builds below run with that warning muffled.
  if (how == "seq") {
    # build each band's parts, then call etSeq() here directly: the builds
    # themselves muffle its warning, so this is the one place it can be heard.
    # EVERY band, since whether a phase sets `ii` can differ band to band.
    gap <- FALSE
    for (j in seq_len(.admDoseNBands(d))) {
      evs <- tryCatch(lapply(parts, .admDoseNodeEv, lo = .admDoseLo(d, j),
                             cov = NULL,
                             where = function() .admDoseBandWhere(d, j),
                             labels = .admDoseLabels()),
                      error = function(e) NULL)
      if (length(evs) > 1L)
        withCallingHandlers(
          tryCatch(do.call(rxode2::etSeq, evs), error = function(e) NULL),
          warning = function(w) if (grepl(.ADM_ETSEQ_NOII, conditionMessage(w))) {
            gap <<- TRUE; invokeRestart("muffleWarning")
          })
      if (gap) break
    }
    if (gap)
      warning("admixr2: in at least one band none of these phases sets a dosing ",
              "interval (`ii`), so rxode2::etSeq() starts each phase 24 h after ",
              "the previous one's last dose. Give `ii`, or a wait in hours ",
              "between the phases, e.g. c(a, 12, b).", call. = FALSE)
  }
  d
}
.ADM_ETSEQ_NOII <- "No inter-dose interval found in event tables"

# The banded covariate and the cut points of one part of a combination; a
# fixed phase or a wait has none.
.admDoseCovOf <- function(p) if (inherits(p, "admDoseBands")) p[["covariate"]]
.admDoseBreaksOf <- function(p)
  if (inherits(p, "admDoseBands")) p[["breaks"]] else numeric(0)

# Every admDoseBands object, whichever kind -- a doseBands() table (`amt`,
# `timing`, `per`), a combination (`parts`, `how`) or a fixed phase (`fixed`,
# `ev`) -- carries every field, NULL where it does not apply, and is read with
# [[ ]]: `$` partial-matches, so a field absent on one kind would silently
# resolve to a longer one (see CLAUDE.md, "`$` PARTIAL-MATCHES ON LISTS").
.admDoseNew <- function(covariate = NULL, breaks = numeric(0), amt = NULL,
                        timing = list(), per = NULL, parts = NULL, how = NULL,
                        fixed = NULL, ev = NULL)
  structure(list(covariate = covariate, breaks = breaks, amt = amt,
                 timing = timing, per = per, parts = parts, how = how,
                 fixed = fixed, ev = ev),
            class = "admDoseBands")

# A phase every band gets: an event table, combined with bands, and checked
# like any band's (doses only, one subject, a non-negative amount).
.admDoseFixedPhase <- function(x)
  .admDoseNew(fixed = TRUE,
              ev = .admDoseEvChecked(x, function() NULL))

#' @export
format.admDoseBands <- function(x, ...) {
  if (!is.null(x[["parts"]]))
    return(paste(vapply(x[["parts"]], function(p)
                   if (is.numeric(p)) paste0("wait ", format(p), " h") else format(p), ""),
                 collapse = if (identical(x[["how"]], "rbind")) "  and  " else "  then  "))
  if (isTRUE(x[["fixed"]]))
    return(paste0("everyone: ", paste(format(as.numeric(.admEvDoses(x[["ev"]])$amt)),
                                      collapse = "+")))
  a <- vapply(x[["amt"]], function(v) {
    if (is.function(v))
      return(paste0("f(", paste(names(formals(v)), collapse = ", "), ")"))
    if (inherits(v, "rxEt")) return("et()")
    s <- paste(vapply(v, format, ""), collapse = "+")
    if (!is.null(x[["per"]]))
      s <- paste0(if (length(v) > 1L) paste0("(", s, ")") else s, " x ", x[["per"]])
    s
  }, "")
  lab <- paste0(vapply(seq_along(x[["amt"]]), function(j) .admDoseBandWhere(x, j), ""),
                ": ", a)
  out <- paste(lab, collapse = "; ")
  if (length(x[["timing"]]))
    out <- paste0(out, "  (", paste(sprintf("%s = %s", names(x[["timing"]]),
      vapply(x[["timing"]], function(v)
        if (is.list(v)) paste(vapply(v, function(e)
          paste(format(e, trim = TRUE), collapse = ","), ""), collapse = " | ")
        else paste(format(v, trim = TRUE), collapse = ","),
             character(1))), collapse = ", "), ")")
  out
}

#' @export
print.admDoseBands <- function(x, ...) {
  if (isTRUE(x[["fixed"]])) {
    cat("Dose for everyone:", format(x), "\n")
    return(invisible(x))
  }
  cat("Dosing by ", x[["covariate"]], " band:\n", sep = "")
  tab <- .admDoseTable(x)
  print(tab, row.names = FALSE, right = FALSE)
  if (any(tab == "?"))
    cat("? : set by that band's function, which returns its own et(); print()",
        "of the study shows it at the population's medians\n")
  invisible(x)
}

# THE TABLE IS WHAT rxode2::et() BUILT, not a re-creation of it: each band's
# event table is built and printed through as.data.frame(), columns and order
# as rxode2 gives them, so a new et() argument shows up with no change here.
# Amounts that need covariates (`per`, a function) are built as stand-in
# numbers -- `.ADM_DOSE_TAG` + an index into `labels` -- and shown as their
# rule ("5 x WT", "f(WT)"). An unset `cmt` is the model's default dosing
# compartment, named when `default_cmt` is known, and an "IV" dose without one
# goes to an oral depot, so it is shown rather than left blank.
.admDoseTable <- function(x, default_cmt = NULL, cov = NULL) {
  dflt <- if (is.null(default_cmt)) "(default)" else paste0(default_cmt, " (default)")
  built_at <- FALSE
  bands <- lapply(seq_len(.admDoseNBands(x)), function(j) {
    lab <- .admDoseLabels()
    d   <- as.data.frame(.admDoseBandEv(x, j, cov = cov, labels = lab))
    built_at <<- built_at || lab$at
    d$id <- NULL
    if (!is.null(d$evid) && all(d$evid %in% 1L)) d$evid <- NULL
    a  <- suppressWarnings(as.numeric(d[["amt"]]))
    tg <- !is.na(a) & a > .ADM_DOSE_TAG
    # element by element: a column formatted whole takes the stand-ins' width
    # and scientific notation (1000 printed as 1e+03)
    out <- lapply(d, function(v) {
      s <- vapply(seq_along(v), function(r) format(v[r], trim = TRUE, digits = 4), "")
      s[is.na(v) | s == "NA"] <- "-"
      s
    })
    # a stand-in shows its rule, keeping any amount unit rxode2 printed
    unit <- sub("^[^[]*", "", out$amt[tg])
    out$amt[tg] <- paste0(lab$lab[a[tg] - .ADM_DOSE_TAG],
                          ifelse(nzchar(unit), paste0(" ", unit), ""))
    if (is.null(out$cmt)) {
      # rxode2 leaves cmt out when no dose sets it; show it next to time
      out$cmt <- rep(dflt, nrow(d))
      at <- match("time", names(out))
      if (!is.na(at))
        out <- out[append(setdiff(names(out), "cmt"), "cmt", after = at)]
    }
    out$cmt[out$cmt %in% c("(default)", "-")] <- dflt
    # a function returning its own et(): nothing but its rule is known here
    wh <- tg
    wh[tg] <- lab$whole[a[tg] - .ADM_DOSE_TAG]
    for (cn in setdiff(names(out), "amt")) out[[cn]][wh] <- "?"
    band <- rep("", nrow(d))
    band[1L] <- sub(paste0("^", x[["covariate"]], " "), "", .admDoseBandWhere(x, j))
    as.data.frame(c(list(band = band), out), stringsAsFactors = FALSE,
                  check.names = FALSE)
  })
  cols <- unique(unlist(lapply(bands, names)))
  df <- do.call(rbind, lapply(bands, function(b) {
    for (cn in setdiff(cols, names(b))) b[[cn]] <- "-"
    b[cols]
  }))
  names(df)[1L] <- x[["covariate"]]
  # TRUE when a band was a function returning et(), built at `cov`
  attr(df, "built_at") <- built_at
  df
}

# Stand-in amounts for the table: a value above .ADM_DOSE_TAG indexes `lab`;
# `whole` marks a function returning its own et(), whose timing is unknown
# until it is called; `at` records that one was built at given covariates.
.ADM_DOSE_TAG <- 1e15
.admDoseLabels <- function() {
  e <- new.env(parent = emptyenv())
  e$lab <- character(0); e$whole <- logical(0); e$at <- FALSE
  e
}

# Whether a band's function returns a whole et() table: called once with every
# argument at 1, which only reads the class. NA when it cannot be called there.
.admDoseFunIsEt <- function(f)
  tryCatch(inherits(suppressWarnings(do.call(f, stats::setNames(
    rep(list(1), length(formals(f))), names(formals(f))))), "rxEt"),
    error = function(e) NA)

# A model's default dosing compartment: rxode2 doses cmt = 1, the first state
# (`depot` for an ODE or oral linCmt() model, `central` for an IV linCmt()).
.admDefaultCmt <- function(ui)
  tryCatch({ s <- ui$state; if (length(s)) as.character(s[[1L]]) },
           error = function(e) NULL)

# THE DOSING RECORDS of an event table: rxode2's evid 1 (dose) and 4 (reset,
# then dose), and their classic 101/104 spellings. evid 2/3 are other/reset
# records, 5/6 replace/multiply a compartment's amount and 7 is a phantom dose
# -- none of them puts an amount in. The ONE place admixr2 decides this; it
# used to be decided in five, which disagreed about 5/6/7.
.admEvDoses <- function(ev) {
  d <- as.data.frame(ev)
  if (is.null(d$evid)) return(d[0L, , drop = FALSE])
  d[!is.na(d$evid) & d$evid %in% c(1, 4, 101, 104), , drop = FALSE]
}

# The total amount a dose record set gives, addl doses included (expanded by
# rxode2::etExpand()), for the step/cap check.
.admDoseTotal <- function(ev)
  sum(as.numeric(.admEvDoses(rxode2::etExpand(ev))$amt))

# Warn when an infusion goes to a default compartment that is the depot -- the
# name rxode2 itself gives an absorption compartment (linCmt() with ka) and the
# conventional one in an ODE model: an "IV" phase written without
# cmt = "centr" doses it, and every number downstream still looks plausible.
# Each band's doses are built by et() at the population's median covariates.
.admDoseCheckCmt <- function(dose, ui, nm, population = NULL) {
  dc <- .admDefaultCmt(ui)
  if (is.null(dc) || !identical(tolower(dc), "depot")) return(invisible())
  med <- .admDoseMedians(dose, population)
  # each band on its own: one band that cannot be built at the medians must not
  # switch the check off for all the others
  quiet <- function(expr) tryCatch(expr, error = function(e) NULL)
  evs <- switch(.admDoseKind(dose),
    ev    = list(dose),
    bands = lapply(seq_len(.admDoseNBands(dose)), function(j)
              quiet(.admDoseBandEv(dose, j, med))),
    fun   = list(quiet(.admDoseEvAt(dose, med))),
    NULL)
  inf <- any(vapply(Filter(Negate(is.null), evs), function(ev) {
    d <- .admEvDoses(ev)
    if (!nrow(d)) return(FALSE)
    num <- function(cn) if (is.null(d[[cn]])) rep(0, nrow(d))
                        else { v <- suppressWarnings(as.numeric(d[[cn]])); v[is.na(v)] <- 0; v }
    cm  <- if (is.null(d[["cmt"]])) rep("(default)", nrow(d)) else as.character(d[["cmt"]])
    # unset, or the default by NUMBER (cmt = 1): the depot without saying so.
    # A dose that NAMES the depot is a choice (zero-order absorption), not a
    # slip, and so is rate/dur = -1/-2, which hands the rate or duration to
    # the model's own rate(depot)/dur(depot): only a POSITIVE one is an infusion.
    dflt <- is.na(cm) | cm %in% c("(default)", "1")
    any((num("dur") > 0 | num("rate") > 0) & dflt)
  }, logical(1)))
  if (inf) {
    oth <- setdiff(ui$state, dc)
    warning("admixr2: study ", sQuote(nm), ": a dose given as an infusion ",
            "(dur or rate) does not name its compartment, so it goes to the ",
            "model's default dosing compartment, ", sQuote(dc), ", which ",
            "absorbs orally. ",
            if (length(oth))
              paste0("An IV dose belongs in the central compartment: give cmt = \"",
                     oth[1L], "\". This model's compartments: ",
                     paste(ui$state, collapse = ", "), ".")
            else paste0("This model has no other compartment, so it has no ",
                        "route for an IV dose."), call. = FALSE)
  }
  invisible()
}

# The population's median of every covariate a dose rule reads, as a named
# list: where a rule is shown or checked at one representative patient.
.admDoseMedians <- function(dose, population) {
  covs <- .admDoseCovs(dose)
  if (!length(covs) || is.null(population)) return(NULL)
  cd <- .admCovDistCanon(population)
  covs <- intersect(covs, .admCovSpecNames(cd))
  stats::setNames(lapply(covs, function(v) .admCovQuantile(cd[[v]], 0.5)), covs)
}

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
    return(unique(c(d[["covariate"]], unlist(lapply(seq_len(.admDoseNBands(d)),
                                               .admDoseBandAmtCovs, d = d)))))
  if (is.function(d)) return(names(formals(d)))
  character(0)
}

# The bands, on the (merged) cut points, as a lo/hi table, open at both ends.
.admDoseNBands <- function(d) length(d[["breaks"]]) + 1L
.admDoseBandTable <- function(d)
  data.frame(lo = c(-Inf, d[["breaks"]]), hi = c(d[["breaks"]], Inf))

# Band j as a reader would say it: "CRCL < 30", "30 to < 60", "CRCL >= 60".
.admDoseBandWhere <- function(d, j) {
  b <- d[["breaks"]]; k <- length(b) + 1L; cv <- d[["covariate"]]
  if (j == 1L) sprintf("%s < %g", cv, b[1L])
  else if (j == k) sprintf("%s >= %g", cv, b[k - 1L])
  else sprintf("%g to < %g", b[j - 1L], b[j])
}

# A band of a combination is named by its lower edge `lo`: each part takes the
# band of its OWN cut points that `lo` falls in (a value on a cut point is in
# the band above). For a plain doseBands() that is band j itself.
.admDoseLo <- function(d, j) c(-Inf, d[["breaks"]])[j]
.admDoseLeafIdx <- function(p, lo) findInterval(lo, p[["breaks"]]) + 1L

# The covariates band j's AMOUNT reads (not the banded one, which only picks the
# band): a function's arguments or `per`, over every part. Empty otherwise.
.admDoseBandAmtCovs <- function(d, j) .admDoseNodeAmtCovs(d, .admDoseLo(d, j))
.admDoseNodeAmtCovs <- function(node, lo) {
  if (is.numeric(node)) return(character(0))       # a wait in a sequence
  if (!is.null(node[["parts"]]))
    return(unique(unlist(lapply(node[["parts"]], .admDoseNodeAmtCovs, lo = lo))))
  if (isTRUE(node[["fixed"]])) return(character(0))
  a <- node[["amt"]][[.admDoseLeafIdx(node, lo)]]
  if (is.function(a)) names(formals(a))
  else if (inherits(a, "rxEt")) character(0)
  else if (!is.null(node[["per"]])) node[["per"]]
  else character(0)
}

# Band j's event table, at node `cov` (a named list of covariate values) when
# its amount reads a covariate. A combination is built part by part and joined
# with rxode2::etSeq() or rxode2::etRbind(). With `labels`, amounts that need
# covariates are stand-ins for the table (see .admDoseTable()).
.admDoseBandEv <- function(d, j, cov = NULL, labels = NULL)
  .admDoseNodeEv(d, .admDoseLo(d, j), cov, function() .admDoseBandWhere(d, j),
                 labels)
.admDoseNodeEv <- function(node, lo, cov, where, labels) {
  # a number is a wait, handed to etSeq() as it is
  if (is.numeric(node)) return(node)
  if (!is.null(node[["parts"]])) {
    evs <- lapply(node[["parts"]], .admDoseNodeEv, lo = lo, cov = cov,
                  where = where, labels = labels)
    if (length(evs) == 1L) return(evs[[1L]])
    if (identical(node[["how"]], "rbind"))
      return(do.call(rxode2::etRbind, evs))
    # the no-`ii` 24-h spacing is announced once, by .admDoseCombine()
    return(withCallingHandlers(do.call(rxode2::etSeq, evs),
      warning = function(w) if (grepl(.ADM_ETSEQ_NOII, conditionMessage(w)))
        invokeRestart("muffleWarning")))
  }
  if (isTRUE(node[["fixed"]])) return(node[["ev"]])
  .admDoseLeafEv(node, .admDoseLeafIdx(node, lo), cov, where, labels)
}

# One doseBands() band's event table. Its amounts -- numbers, `per`, or a
# function's result -- go to rxode2::et() with the band's `...` values; one
# et() call per dose, chained, since et() does not vectorise across its
# arguments. A band given as an et() table, or a function returning one, is
# used as it is.
.admDoseLeafEv <- function(p, i, cov, where, labels = NULL) {
  k <- length(p[["amt"]])
  a <- p[["amt"]][[i]]
  if (inherits(a, "rxEt")) return(.admDoseEvChecked(a, where))
  tim <- lapply(p[["timing"]], .admDoseAt, j = i, k = k)
  if (!is.null(labels)) {
    whole <- FALSE
    if (is.function(a)) {
      f  <- paste0("f(", paste(names(formals(a)), collapse = ", "), ")")
      # a function that returns a WHOLE et() table carries its own timing, which
      # a stand-in cannot show: build it for real where the covariate values
      # are known (print() of a study passes the medians), and otherwise say
      # what it is rather than print one bolus at time 0 that it is not
      if (!length(tim) && isTRUE(.admDoseFunIsEt(a))) {
        if (!is.null(cov) && all(names(formals(a)) %in% names(cov))) {
          labels$at <- TRUE
          return(.admDoseEvChecked(.admDoseLeafAmt(p, i, cov, where), where))
        }
        lab <- paste(f, "-> et()"); whole <- TRUE
      } else {
        m   <- max(c(1L, lengths(tim)))
        lab <- if (m > 1L) paste0(f, " [", seq_len(m), "]") else f
      }
    } else {
      s   <- vapply(a, format, "")
      lab <- if (!is.null(p[["per"]])) paste(s, "x", p[["per"]]) else s
    }
    amt <- .ADM_DOSE_TAG + length(labels$lab) + seq_along(lab)
    labels$lab   <- c(labels$lab, lab)
    labels$whole <- c(labels$whole, rep(whole, length(lab)))
    return(.admDoseEvChecked(.admDoseChain(c(list(amt = amt), tim)), where))
  }
  amt <- .admDoseLeafAmt(p, i, cov, where)
  if (inherits(amt, "rxEt")) {
    if (length(tim))
      stop("admixr2: doseBands(): the function for band ", sQuote(where()),
           " returns a whole et() table, which carries its own timing; give ",
           "the details inside that et() rather than to doseBands().",
           call. = FALSE)
    return(.admDoseEvChecked(amt, where))
  }
  .admDoseEvChecked(.admDoseChain(c(list(amt = amt), tim)), where)
}

# One et() call per dose, chained: values that hold several entries are several
# doses (recycled to the longest), e.g. amt = c(500, 100), cmt = c("depot",
# "centr") for an oral dose plus an infusion.
.admDoseChain <- function(a) {
  m <- max(lengths(a))
  rows <- lapply(seq_len(m), function(r)
    lapply(a, function(v) v[[(r - 1L) %% length(v) + 1L]]))
  Reduce(function(ev, r) do.call(rxode2::et, c(list(ev), r)),
         rows[-1L], do.call(rxode2::et, rows[[1L]]))
}

# Band i's amount(s) at node `cov`, checked: numbers, `per` times the
# covariate, or a function's result (amounts, or an et() table returned whole).
.admDoseLeafAmt <- function(p, i, cov, where) {
  a  <- p[["amt"]][[i]]
  cv <- if (is.function(a)) names(formals(a)) else if (!is.null(p[["per"]])) p[["per"]]
  if (!length(cv)) return(a)
  if (is.null(cov) || !all(cv %in% names(cov)))
    stop("admixr2: internal: band ", i, "'s dose reads ",
         paste(cv, collapse = ", "), ", but no value was given.", call. = FALSE)
  at <- paste(sprintf("%s = %s", cv, vapply(cov[cv], function(v)
    format(signif(as.numeric(v), 4)), character(1))), collapse = ", ")
  w <- paste0("in band ", sQuote(where()), " at ", at)
  r <- if (is.function(a))
         tryCatch(do.call(a, as.list(cov)[cv]), error = function(e)
           stop("admixr2: doseBands(): the dose function failed ", w, ": ",
                conditionMessage(e), call. = FALSE))
       else a * as.numeric(cov[[p[["per"]]]])
  if (inherits(r, "rxEt")) return(r)
  if (!is.numeric(r) || !length(r) || any(!is.finite(r)) || any(r < 0))
    stop("admixr2: doseBands(): the dose ", w, " is ",
         if (is.numeric(r)) paste(format(r), collapse = ", ") else class(r)[1L],
         "; every amount must be a finite, non-negative number. A covariate ",
         "whose margin reaches zero or below (a NORMAL weight) can be bounded ",
         "with `range`, or declared with dist = \"lnorm\".", call. = FALSE)
  as.numeric(r)
}

# What rxode2 built, checked for what a DOSE RULE cannot be: observation
# records (the times come from admStudy(times =); et()'s sampling arguments
# would add them, or leave no dose at all), several subjects (`id`), no dose,
# or a negative amount. Checked on the result, not by argument name, so it
# holds for any et() argument.
.admDoseEvChecked <- function(ev, where) {
  d <- as.data.frame(ev)
  # `where()` names the band; NULL for an et() table added for every band
  w <- function() { b <- where()
    if (is.null(b)) "added for every band" else paste0("in band ", sQuote(b)) }
  if (!is.null(d$id) && length(unique(d$id)) > 1L)
    stop("admixr2: doseBands(): the doses ", w(), " are for several subjects ",
         "(`id`); an aggregate study's dose is one regimen.", call. = FALSE)
  if (!is.null(d$evid) && any(d$evid %in% 0L))
    stop("admixr2: doseBands(): the event table ", w(), " has observation ",
         "records; a dose rule holds doses only, and the observation times ",
         "come from admStudy(times =).", call. = FALSE)
  dz <- .admEvDoses(ev)
  if (!nrow(dz))
    stop("admixr2: doseBands(): the event table ", w(), " holds no dose.",
         call. = FALSE)
  a <- as.numeric(dz$amt)
  if (any(!is.finite(a)) || any(a < 0))
    stop("admixr2: doseBands(): the dose ", w(), " is ",
         paste(format(a), collapse = ", "), "; every amount must be finite and ",
         "non-negative.", call. = FALSE)
  ev
}

# Band j's value of one dosing argument: length-1 is shared, length-k (or a
# list of k) is per band. A band's value may itself hold several doses.
.admDoseAt <- function(v, j, k) {
  if (is.list(v)) v[[if (length(v) == 1L) 1L else j]] else if (length(v) == k && k > 1L) v[[j]] else v
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
  a <- if (is.null(ev)) r else as.numeric(.admEvDoses(ev)$amt)
  if (!length(a) || any(!is.finite(a)) || any(a < 0))
    stop("the dose function gives ", paste(format(a), collapse = ", "),
         " at ", at, "; every amount must be finite and non-negative. A ",
         "covariate whose margin reaches zero or below (a NORMAL weight) can ",
         "be bounded with `range`, or declared with dist = \"lnorm\".",
         call. = FALSE)
  ev %||% rxode2::et(amt = as.numeric(r))
}

# Band j's span [lo, hi] inside the enrolled range `rng` on the banded
# covariate; NULL when the band lies outside it.
.admDoseBandSpan <- function(d, j, rng = NULL) {
  bt  <- .admDoseBandTable(d)
  rng <- sort(as.numeric(rng %||% c(-Inf, Inf)))
  lo  <- max(bt$lo[j], rng[1L]); hi <- min(bt$hi[j], rng[2L])
  if (lo >= hi) NULL else c(lo, hi)
}

# Each band's SHARE of the enrolled patients when an enrolled range on another
# covariate is latently correlated with the banded one. The band's mass is then
# P(banded covariate in band, other in its range), not P(banded in band): with
# WT enrolled at 40-70 kg and rho = 0.6, CrCl < 30 holds 27.5% of the enrolled
# patients against 12.9% of the population. On the latent scale that is
#   int_band phi(z) [Phi((b2 - rho z) / s) - Phi((b1 - rho z) / s)] dz,
# s = sqrt(1 - rho^2) -- one-dimensional quadrature. NULL when no such range
# exists, and the caller keeps the banded covariate's own shares, bit for bit.
.admDoseBandShares <- function(cd, bcv, spans, cov_range) {
  nms <- .admCovSpecNames(cd); R <- cd[["latentR"]]
  oth <- setdiff(intersect(names(cov_range), nms), bcv)
  oth <- oth[vapply(oth, function(v) is.null(cd[[v]][["values"]]), logical(1))]
  if (is.null(R) || !length(oth)) return(NULL)
  ia <- match(bcv, nms)
  cp <- oth[abs(R[ia, match(oth, nms)]) > 1e-12]
  if (!length(cp)) return(NULL)
  Fa <- .admCovCdfOf(cd[[bcv]])
  if (length(cp) > 1L) {
    # several correlated ranges: P(band) x P(all of them | band), the second
    # from the sequential conditioner's UNNORMALISED weights (the band covariate
    # first), whose total is exactly that conditional probability
    mass <- vapply(spans, function(sp) {
      cr  <- c(stats::setNames(list(sp), bcv), cov_range[cp])
      res <- .admCovSeqCore(cd, c(bcv, cp), .ADM_STRATA_NODES, 32768L, cr)
      diff(Fa(sp)) * sum(vapply(res$nodes, function(z) z$weight, numeric(1)))
    }, numeric(1))
    if (!(sum(mass) > 0))
      stop("admixr2: the bands on ", sQuote(bcv), " and the enrolled ranges of ",
           paste(sQuote(cp), collapse = " and "), " leave no patients under the ",
           "declared correlation.", call. = FALSE)
    return(mass / sum(mass))
  }
  rho <- R[ia, match(cp, nms)]; s <- sqrt(1 - rho^2)
  zb  <- stats::qnorm(.admCovCdfOf(cd[[cp]])(sort(as.numeric(cov_range[[cp]]))))
  mass <- vapply(spans, function(sp) {
    za <- stats::qnorm(Fa(sp))
    stats::integrate(function(z) stats::dnorm(z) *
                       (stats::pnorm((zb[2L] - rho * z) / s) -
                        stats::pnorm((zb[1L] - rho * z) / s)),
                     za[1L], za[2L], rel.tol = 1e-10)$value
  }, numeric(1))
  if (!(sum(mass) > 0))
    stop("admixr2: the bands on ", sQuote(bcv), " and the enrolled range of ",
         sQuote(cp), " leave no patients under the declared correlation.",
         call. = FALSE)
  mass / sum(mass)
}

# The covariates band j's nodes are cut over, given the canonical population
# `cd`: what its amount reads, plus the banded covariate itself when the amount
# reads it or it is latently correlated with what the amount reads -- the nodes
# then condition the amount's covariate on the band. NULL for a plain amount.
.admDoseBandNodeSet <- function(cd, d, j) {
  acv <- .admDoseBandAmtCovs(d, j)
  if (!length(acv)) return(NULL)
  bcv <- d[["covariate"]]
  nms <- .admCovSpecNames(cd)
  R   <- cd[["latentR"]]
  oth <- match(setdiff(acv, bcv), nms)
  corr <- !is.null(R) && length(oth) &&
    any(abs(R[match(bcv, nms), oth]) > 1e-12)
  unique(c(if (bcv %in% acv || corr) bcv, acv))
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
# neighbours do; a KINK is a slope that CHANGES AND STAYS CHANGED -- the mean
# slope over the `w` steps after a point differs from the mean over the `w`
# before it by a quarter of the largest such mean. A smooth rule (mg/kg,
# allometric) drifts slowly and trips neither. Comparing NEIGHBOURING steps
# instead read rounding as a cap: round(5 * WT) alternates 3- and 4-mg steps on
# the grid, a slope change of a quarter at every other step, while its window
# means stay at 5 mg/kg.
.admDoseRoughAlong <- function(x, a, w = 5L) {
  if (length(a) < 4L || any(!is.finite(a))) return(NULL)
  da <- diff(a); span <- diff(range(a))
  if (span <= 0) return(NULL)
  ad <- abs(da)
  nb <- pmax(c(0, ad[-length(ad)]), c(ad[-1L], 0))
  j  <- which(ad > 10 * nb & ad > 1e-3 * span)
  if (length(j)) return(list(kind = "step", at = mean(x[j[1L] + 0:1]), j = j[1L]))
  s <- da / diff(x); n <- length(s)
  if (n < 2L * w) return(NULL)
  i  <- w:(n - w)
  cs <- c(0, cumsum(s))
  sl <- (cs[i + 1L] - cs[i - w + 1L]) / w           # mean slope over s[(i-w+1):i]
  sr <- (cs[i + w + 1L] - cs[i + 1L]) / w           # mean slope over s[(i+1):(i+w)]
  ds <- abs(sr - sl); m <- max(abs(c(sl, sr)))
  if (!(m > 0)) return(NULL)
  # a kink is a change CONCENTRATED at one point: big against the slope, and
  # big against the change two windows either side. A saturating curve,
  # 100 * x / (x + 10), changes its slope a lot near its foot but gradually,
  # so it fails the second test; a cap, pmin(15 * x, 900), passes both.
  L  <- length(ds); off <- 2L * w
  nb <- vapply(seq_len(L), function(k) {
    v <- c(if (k - off >= 1L) ds[k - off], if (k + off <= L) ds[k + off])
    if (length(v)) max(v) else NA_real_
  }, numeric(1))
  ok <- ds > 0.25 * m & !is.na(nb) & ds > 3 * nb
  if (any(ok)) {
    k <- which(ok)[which.max(ds[ok])]
    return(list(kind = "kink", at = x[i[k] + 1L]))
  }
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
           else if (identical(hit$covariate, band_of))
             paste0("Add a cut point there: doseBands(", band_of,
                    ", breaks = c(..., ", at, ", ...), ...), giving each side ",
                    "its own amount.")
           # a step in ANOTHER covariate cannot become a cut on this one: the
           # bands split exactly on one covariate only
           else paste0("doseBands() splits exactly on ", band_of, " only, so a step ",
                       "in ", hit$covariate, " inside a band stays approximate; if ",
                       "the rule's main step is in ", hit$covariate, ", band on ",
                       hit$covariate, " instead and put the ", band_of,
                       " dependence into each band's function."))
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
    bands = paste0("dose by ", d[["covariate"]], " band"),
    fun   = paste0("dose = f(", paste(.admDoseCovs(d), collapse = ", "), ")"))
}
