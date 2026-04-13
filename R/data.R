#' Dummy data frame for nlmixr2 dispatch
#'
#' Returns a minimal NONMEM-style data frame that satisfies nlmixr2's data
#' argument requirement. All `DV` values are `NA` so nlmixr2 adds zero
#' log(2pi) constants to OBJF, keeping `fit$objective == our -2LL` exactly.
#'
#' @return A data frame with columns `ID`, `TIME`, `DV`, `AMT`, `EVID`, `CMT`.
#'
#' @examples
#' admData()
#'
#' @export
admData <- function() {
  data.frame(ID   = c(1L, 1L),
             TIME = c(0, 1),
             DV   = c(NA_real_, NA_real_),
             AMT  = c(100, 0),
             EVID = c(101L, 0L),
             CMT  = c(1L, 2L))
}
