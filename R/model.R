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
# digest, to prevent Windows GC finalizer heap corruption
# (STATUS_HEAP_CORRUPTION / -1073740940). This pin is Windows-only
# (.admPinCompanions()): off Windows the finalizers are safe and pinning every
# model's companions just accumulates native memory. Both the companion pin and
# the sens-result cache go through .adm_pin_env's LRU (.adm_pin_limit()), so a
# long multi-model session does not grow without bound. See R/zzz.R.
.admLoadSensModel <- function(ui) {
  .model_key <- digest::digest(ui$lstExpr)
  .sens_key  <- paste0("sens_",  .model_key)
  .pin_key   <- paste0("focei_", .model_key)

  # In-memory cache (LRU): avoids disk read and rxLoad on repeat calls within a
  # session. .admPinGet touches the entry so it survives eviction while in use.
  .cached <- .admPinGet(.sens_key)
  if (!is.null(.cached)) return(.cached)

  .focei_model <- tryCatch(ui$foceiModel, error = function(e) NULL)
  # Pin the full foceiModel to keep companion objects ($outer, $predOnly,
  # $predNoLhs) off the GC. See pinning note in function header. Windows-only:
  # off Windows the finalizers are safe and pinning only accumulates native
  # memory (.admPinCompanions()), so we let .focei_model fall out of scope here
  # and rxode2 reclaim its companion DLLs.
  if (.admPinCompanions())
    tryCatch(.admPinSet(.pin_key, .focei_model), error = function(e) NULL)
  inner <- .focei_model$inner
  if (is.null(inner)) return(NULL)

  lhs <- tryCatch(inner$lhs, error = function(e) NULL)
  if (is.null(lhs)) return(NULL)

  ini_df <- tryCatch(ui$iniDf, error = function(e) NULL)
  if (is.null(ini_df)) return(NULL)
  eta_rows <- ini_df[!is.na(ini_df$neta1) & ini_df$neta1 == ini_df$neta2 & !ini_df$fix, ]
  th_rows  <- ini_df[!is.na(ini_df$ntheta), ]
  n_eta    <- nrow(eta_rows)

  # Index by ntheta, NOT by position among the non-fixed thetas. The sens model's
  # THETA[k] is numbered by ntheta and INCLUDES fixed thetas, so dropping a fixed
  # theta from the map shifted every later theta into the wrong slot.
  rename_map <- c(
    setNames(paste0("THETA[", th_rows$ntheta, "]"), th_rows$name),
    setNames(paste0("ETA[",   seq_len(n_eta), "]"),
             paste0("eta.", gsub("^eta\\.", "", eta_rows$name)))
  )

  # A FIXED theta is not an estimated parameter, so it never reaches the solve
  # paths (pinfo carries only the estimated ones) -- but the sens model still has
  # a THETA[k] slot for it and rxSolve REQUIRES every parameter. Left unset the
  # sens solve errors and returns NULL, which silently drops admc/adfo to a
  # finite-difference gradient and, worse, made .adghGrad skip the study entirely.
  # Carry the fixed values so the solve paths can fill those columns.
  fix_rows <- th_rows[th_rows$fix, , drop = FALSE]
  fixed_theta <- if (nrow(fix_rows) > 0L)
    setNames(as.numeric(fix_rows$est), paste0("THETA[", fix_rows$ntheta, "]"))
  else numeric(0)

  sens_cols <- lhs[grepl("sens_rx_pred.*ETA|sens.*pred.*BY.*ETA", lhs, ignore.case = TRUE)]
  if (length(sens_cols) == 0L) return(NULL)

  eta_idx <- suppressWarnings(as.integer(regmatches(sens_cols, regexpr("[0-9]+", sens_cols))))
  if (any(is.na(eta_idx))) return(NULL)
  sens_cols <- sens_cols[order(eta_idx)]
  if (length(sens_cols) != n_eta) return(NULL)

  # Cache key: the MODEL (ui$lstExpr) plus a schema tag for the fields we store
  # alongside the compiled model.
  #
  # NOT digest(inner): ui$foceiModel$inner returns a DIFFERENT object on its first
  # access than on later ones (the same foceiModel instability documented for the
  # pin), so digesting it gives an unstable key -- the path written on the first
  # (compile) load does not match the path recomputed on a later load, so the
  # disk-cache branch below never fires in-session and the sens model is silently
  # RECOMPILED on every reload. ui$lstExpr is the model source: stable, and already
  # the key used by the in-memory pin (.sens_key) and by .admLoadModel's own cache.
  #
  # The schema tag makes a cache written BEFORE the fixed-theta fix a miss: without
  # it a parallel worker -- which reads this file directly and cannot re-derive --
  # would keep the old position-indexed rename_map and a NULL fixed_theta, so the
  # fit would silently disagree between workers = 1 and workers > 1.
  .cacheFile <- file.path(
    rxode2::rxTempDir(),
    paste0("adm-sens-",
           digest::digest(list(ui$lstExpr, "ntheta-map+fixed-theta")), ".qs2")
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
      # Overwrite the derived fields from the parent's fresh derivation rather than
      # trusting the file. The cache key carries a schema tag, so a hit is always
      # current-schema -- but a parallel WORKER reads this same file and cannot
      # re-derive, so what the parent writes here is what the worker gets. Trusting
      # a stale position-indexed rename_map would put a theta's value in the wrong
      # THETA[k] slot and silently diverge the parallel fit from the sequential one.
      result$cache_file  <- .cacheFile
      result$rename_map  <- rename_map
      result$fixed_theta <- fixed_theta
      result$sens_cols   <- sens_cols
      tryCatch(.admPinSet(.sens_key, result), error = function(e) NULL)
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
                 rename_map = rename_map, fixed_theta = fixed_theta,
                 is_lincmt = is_lincmt, cache_file = .cacheFile)
  tryCatch(qs2::qs_save(result, .cacheFile), error = function(e) NULL)
  tryCatch(.admPinSet(.sens_key, result), error = function(e) NULL)
  result
}
