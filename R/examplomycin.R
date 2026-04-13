#' Examplomycin dataset
#'
#' A simulated pharmacokinetic dataset for the fictional drug examplomycin,
#' intended as a worked example for aggregate data modelling with `admixr2`.
#' The dataset contains 500 subjects, each with 9 observation time points,
#' generated from a two-compartment model with first-order absorption.
#'
#' **True population parameters:**
#'
#' | Parameter | Value |
#' |-----------|-------|
#' | CL (L/hr) | 5 |
#' | V1 (L) | 10 |
#' | V2 (L) | 30 |
#' | Q (L/hr) | 10 |
#' | ka (1/hr) | 1 |
#' | IIV (all, SD on log scale) | 0.3 |
#' | Proportional error (SD) | 0.2 |
#'
#' Single oral dose of 100 mg; sampling at 0.1, 0.25, 0.5, 1, 2, 3, 5, 8,
#' and 12 hours post-dose.
#'
#' @format A data frame with 5000 rows and 6 columns:
#'   - `ID`: Subject identifier (integer, 1--500).
#'   - `TIME`: Time after dose (hours).
#'   - `DV`: Observed plasma concentration (mg/L).
#'   - `AMT`: Dose amount (mg); 100 for dosing records, 0 otherwise.
#'   - `EVID`: Event type (101 = dose, 0 = observation).
#'   - `CMT`: Compartment (1 = depot, 2 = central).
#'
#' @examples
#' data("examplomycin")
#' head(examplomycin)
#'
#' # Compute aggregate statistics
#' obs <- examplomycin[examplomycin$EVID == 0, ]
#' obs <- obs[order(obs$ID, obs$TIME), ]
#' times <- sort(unique(obs$TIME))
#' E <- sapply(times, function(t) mean(obs$DV[obs$TIME == t]))
#' round(E, 3)
#'
#' @source Generated from a two-compartment PK model using `rxode2::rxSolve()`.
#'   See `vignette("admixr2")` for a full modelling example.
"examplomycin"
