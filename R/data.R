#' Dummy data frame for nlmixr2 dispatch
#'
#' Returns a minimal NONMEM-style data frame that satisfies nlmixr2's data
#' argument requirement. The single observation row carries a non-`NA`
#' placeholder `DV` (`1`); the dose row keeps `DV = NA`. The placeholder is
#' purely for dispatch and output construction and never enters the reported
#' objective -- each estimator overwrites `fit$env$objective` (and
#' OBJF/logLik/AIC/BIC) with its own aggregate -2LL. A non-`NA` observation is
#' required because nlmixr2's post-fit output construction
#' (`nlmixr2CreateOutputFromUi`) solves the model over this frame, and rxode2's
#' event-table translation rejects a dataset with no non-`NA` observation rows
#' ("no rows in event table or input data"), the same reason the multi-endpoint
#' frame below uses a placeholder `DV`.
#'
#' @param outputs Optional character vector of observed output (endpoint) names
#'   for a multi-compartment model with several prediction lines (e.g.
#'   `c("cp", "cCSF")`). One observation row is emitted per endpoint, keyed by
#'   name in the `DVID` column (nlmixr2's endpoint identifier; `CMT` is `NA` on
#'   those rows), so nlmixr2's data translation recognises every endpoint. These
#'   rows carry a non-`NA` placeholder `DV` (`1`) because nlmixr2's
#'   multi-endpoint translator rejects an all-`NA`-DV dataset; the placeholder is
#'   purely for dispatch and never enters the reported objective (each estimator
#'   overwrites it with its own aggregate -2LL). When `NULL` (default) the
#'   single-endpoint dummy frame is returned unchanged.
#'
#' @param occasions Optional named number of occasion levels for an inter-occasion
#'   variability term (`iov.cl ~ 0.1 | occ`), e.g. `c(occ = 3)`. Adds that column
#'   with the levels `1:3`: nlmixr2est reads the levels of an occasion variable
#'   from the data it is given, and an aggregate fit has no per-record column.
#'   Which times belong to which occasion is given on each study, with
#'   `admStudy(occasions = )`. Needs at least two levels. Supported by `adgh` and
#'   `admc`.
#'
#' @return A data frame with columns `ID`, `TIME`, `DV`, `AMT`, `EVID`, `CMT`
#'   (single-endpoint), plus a `DVID` endpoint column when `outputs` is given.
#'
#' @examples
#' admData()
#' admData(c("cp", "cCSF"))
#' admData(occasions = c(occ = 2))
#'
#' @export
admData <- function(outputs = NULL, occasions = NULL) {
  out <- .admData(outputs)
  # nlmixr2est reads the occasion LEVELS of an inter-occasion term
  # (`iov.cl ~ 0.1 | occ`) off this column; the aggregate study supplies the
  # times, so only the levels matter here.
  if (!is.null(occasions)) {
    if (!is.numeric(occasions) || is.null(names(occasions)) ||
        any(!nzchar(names(occasions))))
      stop("`occasions` must be a named numeric, e.g. occasions = c(occ = 3).",
           call. = FALSE)
    for (v in names(occasions)) {
      k <- occasions[[v]]
      if (length(k) != 1L || !is.finite(k) || k < 2 || k != round(k))
        stop("`occasions[[\"", v, "\"]]` must be a whole number of occasions >= 2.",
             call. = FALSE)
      out <- .admWithLevels(out, v, seq_len(k))
    }
  }
  out
}

# `var` as a column holding every one of `levels`: the frame is padded with copies
# of its last row when it has fewer rows than levels, so no level is lost.
.admWithLevels <- function(d, var, levels) {
  while (nrow(d) < length(levels)) {
    r <- d[nrow(d), , drop = FALSE]
    if ("TIME" %in% names(r)) r$TIME <- r$TIME + 1
    d <- rbind(d, r)
  }
  d[[var]] <- rep_len(levels, nrow(d))
  d
}

.admData <- function(outputs = NULL) {
  if (is.null(outputs))
    return(data.frame(ID   = c(1L, 1L),
                      TIME = c(0, 1),
                      DV   = c(NA_real_, 1),
                      AMT  = c(100, 0),
                      EVID = c(101L, 0L),
                      CMT  = c(1L, 2L)))

  outputs <- as.character(outputs)
  n_o     <- length(outputs)
  # DVID keys each endpoint (cf. nlmixr2data::warfarin); dose row targets CMT = 1.
  # Real dosing/observation times come from each study's `ev`.
  data.frame(ID   = 1L,
             TIME = c(0, seq_len(n_o)),
             DV   = c(NA_real_, rep(1, n_o)),
             AMT  = c(100, rep(0, n_o)),
             EVID = c(1L, rep(0L, n_o)),
             CMT  = c(1L, rep(NA_integer_, n_o)),
             DVID = c(NA_character_, outputs),
             stringsAsFactors = FALSE)
}
