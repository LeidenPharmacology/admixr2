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
#' @return A data frame with columns `ID`, `TIME`, `DV`, `AMT`, `EVID`, `CMT`
#'   (single-endpoint), plus a `DVID` endpoint column when `outputs` is given.
#'
#' @examples
#' admData()
#' admData(c("cp", "cCSF"))
#'
#' @export
admData <- function(outputs = NULL) {
  if (is.null(outputs))
    return(data.frame(ID   = c(1L, 1L),
                      TIME = c(0, 1),
                      DV   = c(NA_real_, 1),
                      AMT  = c(100, 0),
                      EVID = c(101L, 0L),
                      CMT  = c(1L, 2L)))

  outputs <- as.character(outputs)
  n_o     <- length(outputs)
  # Multi-endpoint dummy: one observation per endpoint, keyed by name in a `DVID`
  # column (nlmixr2's endpoint identifier, cf. nlmixr2data::warfarin). The dose
  # row targets a real dosing compartment (`CMT = 1`) rather than an endpoint.
  # As in the single-endpoint frame the observation DV values are a non-NA
  # placeholder (1): nlmixr2's multi-endpoint data translation rejects an
  # all-NA-DV dataset. This is purely for dispatch -- each estimator overwrites
  # `fit$env$objective` (and OBJF/logLik/AIC/BIC) with its own aggregate -2LL,
  # so the placeholder never enters the reported objective; actual dosing and
  # observation times come from each study's `ev`.
  data.frame(ID   = 1L,
             TIME = c(0, seq_len(n_o)),
             DV   = c(NA_real_, rep(1, n_o)),
             AMT  = c(100, rep(0, n_o)),
             EVID = c(1L, rep(0L, n_o)),
             CMT  = c(1L, rep(NA_integer_, n_o)),
             DVID = c(NA_character_, outputs),
             stringsAsFactors = FALSE)
}
