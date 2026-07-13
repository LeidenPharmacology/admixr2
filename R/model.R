# Load (or compile + cache) the rxode2 simulation model.
# Compiled DLL is cached to disk via qs2, keyed by model digest.
.admLoadModel <- function(ui) {
  # Accessing $simulationModel (below) caches the compiled model in
  # ui$meta$.simModelBase as a side effect -- a live, self-referential rxode2
  # object that breaks nlmixr2's ui-cloning during fit assembly. Drop it (and any
  # sibling artifacts) on every exit so the ui stays in the canonical state
  # nlmixr2 expects; see .admDropSimModelMeta() for the full rationale.
  on.exit(.admDropSimModelMeta(ui), add = TRUE)
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

# Remove transient rxode2 model objects that $simulationModel / $foceiModel leave
# behind in ui$meta.
#
# nlmixr2's output machinery (nlmixr2CreateOutputFromUi -> ... -> nmObjGet.*)
# deep-clones the ui with nlmixr2est's internal .cloneEnv(), which recurses into
# every environment-valued member and has no cycle detection. rxode2's compiled
# model objects hold a back-reference to the global .rxModels registry
# (registry -> model -> .rx -> .rxModels -> registry ...), so cloning one loops
# forever -- surfacing as "evaluation nested too deeply: infinite recursion"
# (interactive) or "node stack overflow" (batch). A normal nlmixr2 fit never
# hits this because its estimators do not populate ui$meta with these objects;
# admixr2 does, because it simulates via $simulationModel. Keeping our ui clean
# is the in-framework fix: no wrapping of nlmixr2's code, we just do not feed it
# a ui it was never designed to clone. Safe because admixr2 simulates via its
# own cached model (the return value of .admLoadModel), and rxode2 regenerates
# these lazily if any downstream method needs them.
.admDropSimModelMeta <- function(ui) {
  .meta <- ui$meta
  if (!is.environment(.meta)) return(invisible())
  for (.nm in ls(.meta, all.names = TRUE)) {
    .v <- get(.nm, envir = .meta, inherits = FALSE)
    if (is.environment(.v) && inherits(.v, "rxode2"))
      rm(list = .nm, envir = .meta)
  }
  invisible()
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

  # The cache holds the whole result list, so the compiled model to rxLoad() is
  # m$mod -- rxLoad(m) errors ("need an rxode2-type object"), which silently sent
  # every cache hit down the recompile path.
  if (file.exists(.cacheFile)) {
    result <- tryCatch({ m <- qs2::qs_read(.cacheFile); rxode2::rxLoad(m$mod); m },
                       error = function(e) NULL)
    if (!is.null(result)) {
      # Caches written by older versions predate these fields.
      result$cache_file <- .cacheFile
      if (is.null(result$sens_cols))  result$sens_cols  <- sens_cols
      if (is.null(result$rename_map)) result$rename_map <- rename_map
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

  # cache_file travels with the result: parallel workers reload the sens model
  # from exactly the file that was written here. Re-deriving the key from a later
  # ui$foceiModel$inner does not work -- that access yields a different object
  # than the one digested above, so the lookup always missed.
  result <- list(type = "ode", mod = mod, sens_cols = sens_cols,
                 rename_map = rename_map, is_lincmt = is_lincmt,
                 cache_file = .cacheFile)
  tryCatch(qs2::qs_save(result, .cacheFile), error = function(e) NULL)
  tryCatch(assign(.sens_key, result, envir = .adm_pin_env), error = function(e) NULL)
  result
}
