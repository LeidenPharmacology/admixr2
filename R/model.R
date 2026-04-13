# Load (or compile + cache) the rxode2 simulation model.
# Compiled DLL is cached to disk via qs2, keyed by model digest.
.admLoadModel <- function(ui) {
  .model_key <- digest::digest(ui$lstExpr)
  .cacheFile <- file.path(
    rxode2::rxTempDir(),
    paste0("adm-sim-", .model_key, ".qs2")
  )
  if (file.exists(.cacheFile)) {
    mod <- tryCatch(qs2::qs_read(.cacheFile), error = function(e) NULL)
    load_ok <- !is.null(mod) &&
      tryCatch({ rxode2::rxLoad(mod); TRUE }, error = function(e) FALSE)
    if (load_ok) {
      .old_wd <- tryCatch(getwd(), error = function(e) NULL)
      on.exit(if (!is.null(.old_wd)) setwd(.old_wd), add = TRUE)
      setwd(rxode2::rxTempDir())
      suppressMessages(tryCatch(rxode2::rxode2(ui), error = function(e) NULL))
      return(mod)
    }
    tryCatch(file.remove(.cacheFile), error = function(e) NULL)
  }
  # rxode2 compilation calls setwd() internally -- save/restore to avoid
  # "cannot change working directory" error on first compile (Windows).
  .old_wd <- tryCatch(getwd(), error = function(e) NULL)
  on.exit(if (!is.null(.old_wd)) setwd(.old_wd), add = TRUE)
  setwd(rxode2::rxTempDir())
  mod <- rxode2::rxode2(ui)$simulationModel
  tryCatch(suppressWarnings(qs2::qs_save(mod, .cacheFile)), error = function(e) NULL)
  rxode2::rxLoad(mod)
  mod
}

# Load the rxode2 sensitivity model (ui$foceiModel$inner) if available.
#
# Returns list(type="ode", mod, sens_cols, rename_map, is_lincmt) or NULL.
# Works for both ODE and linCmt models; ui$foceiModel$inner is non-NULL for
# both after compilation.
#
# Pinning: ui$foceiModel creates companion objects ($outer, $predOnly,
# $predNoLhs) with live C++ DLL pointers. rxUi is a locked environment so
# assign(..., envir = ui) fails silently. Instead we pin the full foceiModel
# result in .adm_pin_env (package-level, always writable), keyed by model
# digest. This keeps companions alive for the session and prevents Windows GC
# finalizer heap corruption (STATUS_HEAP_CORRUPTION / -1073740940).
.admLoadSensModel <- function(ui) {
  .model_key <- digest::digest(ui$lstExpr)
  .sens_key  <- paste0("sens_",  .model_key)
  .pin_key   <- paste0("focei_", .model_key)

  # In-memory cache: avoids disk read and rxLoad on repeat calls within a session.
  .cached <- tryCatch(
    get(.sens_key, envir = .adm_pin_env, inherits = FALSE),
    error = function(e) NULL
  )
  if (!is.null(.cached)) return(.cached)

  .focei_model <- tryCatch(ui$foceiModel, error = function(e) NULL)
  # Pin the full foceiModel to keep companion objects ($outer, $predOnly,
  # $predNoLhs) alive. See pinning note in function header.
  tryCatch(assign(.pin_key, .focei_model, envir = .adm_pin_env), error = function(e) NULL)
  inner <- .focei_model$inner
  if (is.null(inner)) return(NULL)

  lhs <- tryCatch(inner$lhs, error = function(e) NULL)
  if (is.null(lhs)) return(NULL)

  ini_df <- tryCatch(ui$iniDf, error = function(e) NULL)
  if (is.null(ini_df)) return(NULL)
  eta_rows    <- ini_df[!is.na(ini_df$neta1) & ini_df$neta1 == ini_df$neta2 & !ini_df$fix, ]
  struct_rows <- ini_df[is.na(ini_df$neta1) & !ini_df$fix, ]
  n_eta       <- nrow(eta_rows)

  rename_map <- c(
    setNames(paste0("THETA[", seq_len(nrow(struct_rows)), "]"), struct_rows$name),
    setNames(paste0("ETA[",   seq_len(n_eta),             "]"),
             paste0("eta.", gsub("^eta\\.", "", eta_rows$name)))
  )

  sens_cols <- lhs[grepl("sens_rx_pred.*ETA|sens.*pred.*BY.*ETA", lhs, ignore.case = TRUE)]
  if (length(sens_cols) == 0L) return(NULL)

  eta_idx <- suppressWarnings(as.integer(regmatches(sens_cols, regexpr("[0-9]+", sens_cols))))
  if (any(is.na(eta_idx))) return(NULL)
  sens_cols <- sens_cols[order(eta_idx)]
  if (length(sens_cols) != n_eta) return(NULL)

  .cacheFile <- file.path(
    rxode2::rxTempDir(),
    paste0("adm-sens-", digest::digest(inner), ".qs2")
  )

  .old_wd <- tryCatch(getwd(), error = function(e) NULL)
  on.exit(if (!is.null(.old_wd)) setwd(.old_wd), add = TRUE)
  setwd(rxode2::rxTempDir())

  if (file.exists(.cacheFile)) {
    mod <- tryCatch({ m <- qs2::qs_read(.cacheFile); rxode2::rxLoad(m); m },
                    error = function(e) NULL)
    if (!is.null(mod)) {
      result <- list(type = "ode", mod = mod, sens_cols = sens_cols,
                     rename_map = rename_map, is_lincmt = FALSE)
      tryCatch(assign(.sens_key, result, envir = .adm_pin_env), error = function(e) NULL)
      return(result)
    }
  }

  # inner is already a compiled "rxode2" object -- load its DLL directly.
  mod <- tryCatch({ rxode2::rxLoad(inner); inner }, error = function(e) NULL)
  # Fallback: re-compile if load fails (e.g., stale DLL path after clean session).
  if (is.null(mod))
    mod <- tryCatch({ m <- rxode2::rxode2(inner); rxode2::rxLoad(m); m },
                    error = function(e) NULL)
  if (is.null(mod)) return(NULL)

  mvars     <- tryCatch(rxode2::rxModelVars(mod), error = function(e) NULL)
  is_lincmt <- if (!is.null(mvars))
    any(grepl("linCmtB", mvars$model, fixed = TRUE)) else FALSE

  result <- list(type = "ode", mod = mod, sens_cols = sens_cols,
                 rename_map = rename_map, is_lincmt = is_lincmt)
  tryCatch(qs2::qs_save(result, .cacheFile), error = function(e) NULL)
  tryCatch(assign(.sens_key, result, envir = .adm_pin_env), error = function(e) NULL)
  result
}
