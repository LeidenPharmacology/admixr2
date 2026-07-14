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

# Structural thetas with no usable mu-referenced eta ("unpaired"): the ones whose
# gradient cannot come from an eta sensitivity column and would otherwise be
# finite-differenced. These are the thetas that get a dummy eta in the augmented
# sens model.
#
# Uses .admMuRefPairs() -- the SAME map pinfo$struct_has_eta is built from -- so
# the set of thetas the estimators route through the theta columns and the set the
# sens model actually builds columns for cannot drift apart. That includes the
# shared-eta guard: a theta whose eta appears in another parameter is unpaired.
.admUnpairedThetas <- function(ui) {
  ini <- tryCatch(ui$iniDf, error = function(e) NULL)
  if (is.null(ini)) return(character(0))
  struct <- ini[is.na(ini$neta1) & is.na(ini$err) & !ini$fix, , drop = FALSE]
  mrd    <- .admMuRefPairs(ui)
  paired <- if (!is.null(mrd)) as.character(mrd$theta) else character(0)
  setdiff(struct$name, paired)
}

# Substitute symbols throughout an expression tree: map is a named list of
# replacement expressions keyed by symbol name.
.admSubstSym <- function(e, map) {
  if (is.symbol(e)) {
    r <- map[[as.character(e)]]
    return(if (is.null(r)) e else r)
  }
  if (is.call(e)) {
    # index 1 is the function itself -- never substituted (a theta is never a call head)
    for (i in seq_along(e)[-1L]) e[[i]] <- .admSubstSym(e[[i]], map)
  }
  e
}

# Build a SENSITIVITY-ONLY copy of the model in which every unpaired structural
# theta gains a dummy eta: `tka` is rewritten to `(tka + eta.admSens.tka)`
# wherever it appears, and `eta.admSens.tka ~ <var>` is appended to ini().
#
# Why: rxode2's focei inner model only emits d(pred)/d(eta_j). With the dummy eta
# pinned at 0 the prediction is unchanged (verified to ~1e-10), while the chain
# rule makes d(pred)/d(eta.admSens.tka) == d(pred)/d(tka) EXACTLY -- so an extra
# sensitivity direction replaces the finite-difference solve the estimators used
# to run for these thetas. Substituting the SYMBOL (not the parameter line) means
# it does not matter how the parameter is written: `exp(tka)`, `exp(tka)*exp(eta)`
# and a theta used in several equations all give the correct total derivative.
#
# This works for linCmt() too (the derivative goes through rxode2's linCmtB chain
# rule), where nlmixr2est's own augmented outer model cannot go -- that one needs
# SECOND-order state sensitivities, which linCmtB does not provide.
#
# The dummy omega variance is never used: this ui is only ever a source of
# sensitivity columns, is never fitted, and is always solved at eta_dummy = 0.
# It must be non-zero, though, or nlmixr2's zero-omega pre-processing drops the eta.
#
# Returns the augmented rxUi, or NULL if anything about the rewrite fails (the
# caller then falls back to the plain sens model + FD).
.admBuildSensUi <- function(ui, unpaired, dummy_var = 0.1) {
  if (length(unpaired) == 0L) return(NULL)
  fn <- tryCatch(as.function(ui), error = function(e) NULL)
  if (is.null(fn)) return(NULL)
  b  <- body(fn)
  idx_ini <- idx_mod <- NA_integer_
  for (i in seq_along(b)) {
    if (!is.call(b[[i]])) next
    head_i <- as.character(b[[i]][[1L]])
    if (identical(head_i, "ini"))   idx_ini <- i
    if (identical(head_i, "model")) idx_mod <- i
  }
  if (is.na(idx_ini) || is.na(idx_mod)) return(NULL)

  eta_nms <- paste0("eta.admSens.", unpaired)
  map <- stats::setNames(
    lapply(seq_along(unpaired), function(i)
      bquote((.(as.symbol(unpaired[i])) + .(as.symbol(eta_nms[i]))))),
    unpaired)

  mblk <- b[[idx_mod]][[2L]]
  for (i in seq_along(mblk)[-1L]) mblk[[i]] <- .admSubstSym(mblk[[i]], map)
  b[[idx_mod]][[2L]] <- mblk

  iblk <- b[[idx_ini]][[2L]]
  for (nm in eta_nms) iblk[[length(iblk) + 1L]] <- bquote(.(as.symbol(nm)) ~ .(dummy_var))
  b[[idx_ini]][[2L]] <- iblk

  body(fn) <- b
  tryCatch(rxode2::rxode2(fn), error = function(e) NULL)
}

# Load the rxode2 sensitivity model (ui$foceiModel$inner) if available.
#
# Returns list(type="ode", mod, sens_cols, theta_sens_cols, dummy_eta_inner,
# rename_map, is_lincmt, cache_file) or NULL. Works for both ODE and linCmt
# models; ui$foceiModel$inner is non-NULL for both after compilation.
#
# The model compiled here is the DUMMY-ETA augmented one (.admBuildSensUi) when
# the model has unpaired structural thetas, so the solve also returns
# d(pred)/d(theta) for those (theta_sens_cols); `dummy_eta_inner` names the extra
# ETA[k] columns the solve paths must pin at 0. The first n_eta sensitivity
# columns are the real etas either way, so every existing consumer of `sens_cols`
# is unaffected. Any failure in the augmented build falls back to the plain
# inner model with theta_sens_cols = NULL, which routes the estimators back to
# their finite-difference path for those thetas.
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

  # In-memory cache: avoids disk read and rxLoad on repeat calls within a session.
  .cached <- tryCatch(
    get(.sens_key, envir = .adm_pin_env, inherits = FALSE),
    error = function(e) NULL
  )
  if (!is.null(.cached)) return(.cached)

  ini_df <- tryCatch(ui$iniDf, error = function(e) NULL)
  if (is.null(ini_df)) return(NULL)
  n_eta_real <- nrow(ini_df[!is.na(ini_df$neta1) & ini_df$neta1 == ini_df$neta2 &
                              !ini_df$fix, , drop = FALSE])

  # Augmented (theta-sensitivity) model first; plain model as the fallback.
  unpaired <- if (n_eta_real > 0L) .admUnpairedThetas(ui) else character(0)
  sens_ui  <- if (length(unpaired) > 0L)
    .admBuildSensUi(ui, unpaired) else NULL

  result <- NULL
  if (!is.null(sens_ui))
    result <- tryCatch(.admSensFromUi(sens_ui, ui, unpaired),
                       error = function(e) NULL)
  if (is.null(result))
    result <- tryCatch(.admSensFromUi(ui, ui, character(0)),
                       error = function(e) NULL)
  if (is.null(result)) return(NULL)

  tryCatch(assign(.sens_key, result, envir = .adm_pin_env), error = function(e) NULL)
  result
}

# Compile/load the sens model from `sens_ui` (augmented or plain) and describe it
# relative to the ORIGINAL `ui` (whose parameter names the estimators speak).
# `unpaired` is empty for the plain model. Returns NULL on any failure.
.admSensFromUi <- function(sens_ui, ui, unpaired) {
  .pin_key <- paste0("focei_", digest::digest(sens_ui$lstExpr))

  .focei_model <- tryCatch(sens_ui$foceiModel, error = function(e) NULL)
  # Pin the full foceiModel to keep companion objects ($outer, $predOnly,
  # $predNoLhs) alive. See pinning note in .admLoadSensModel's header.
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
  # The augmented model appends one dummy eta per unpaired theta, in `unpaired`
  # order, AFTER the real etas -- so its sens columns split n_eta | n_unpaired.
  n_unp <- length(unpaired)
  if (length(sens_cols) != n_eta + n_unp) return(NULL)

  theta_sens_cols <- NULL
  dummy_eta_inner <- character(0)
  if (n_unp > 0L) {
    theta_sens_cols <- setNames(sens_cols[n_eta + seq_len(n_unp)], unpaired)
    dummy_eta_inner <- paste0("ETA[", n_eta + seq_len(n_unp), "]")
    sens_cols       <- sens_cols[seq_len(n_eta)]
  }

  # Cache key: the inner model AND everything that changes the model we compile
  # from it. "jump" = the eventSens flag below; without it in the key, a cache
  # written before that fix would be reloaded and silently serve a model whose
  # dose-parameter sensitivities are zero.
  .cacheFile <- file.path(
    rxode2::rxTempDir(),
    paste0("adm-sens-", digest::digest(list(inner, "jump")), ".qs2")
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
      # Caches written by older versions predate these fields. The theta-sens
      # fields are re-derived rather than trusted: the cache key is the inner
      # model, so a cache written before this feature has the right model but no
      # theta columns.
      result$cache_file      <- .cacheFile
      result$theta_sens_cols <- theta_sens_cols
      result$dummy_eta_inner <- dummy_eta_inner
      if (is.null(result$sens_cols))  result$sens_cols  <- sens_cols
      if (is.null(result$rename_map)) result$rename_map <- rename_map
      return(result)
    }
  }

  # Recompile the inner model with eventSens = "jump".
  #
  # WITHOUT this, a parameter that enters a DOSING MODIFIER -- bioavailability
  # f(), lag()/alag(), rate(), dur() -- has a sensitivity of exactly ZERO in the
  # solve, silently. nlmixr2est's inner model does not carry event/dose-parameter
  # sensitivities (FOCEI computes those separately, by finite differences, via its
  # `predNoLhs` "events FD model"), but admixr2 reads the inner model's columns
  # directly. So `grad = "sens"` produced a zero gradient for every eta -- and now
  # every theta -- driving a modelled F/lag/rate/dur. Verified: analytic column 0
  # vs a true derivative of ~0.4-0.7.
  #
  # eventSens = "jump" attaches rxode2's analytic forward variational jumps at
  # dose times for the rx__sens_* compartments; this is the same flag nlmixr2est
  # passes when it builds its augmented outer model, for exactly this reason.
  # Predictions are bit-identical (verified 0.0e+00) and ordinary (non-dose)
  # sensitivities are unchanged, so this is a strict fix. Falls back to the inner
  # model as-is if the recompile fails.
  .normMod <- tryCatch(rxode2::rxModelVars(inner)$model[["normModel"]],
                       error = function(e) NULL)
  mod <- if (!is.null(.normMod))
    tryCatch({ m <- rxode2::rxode2(.normMod, eventSens = "jump"); rxode2::rxLoad(m); m },
             error = function(e) NULL)
  else NULL
  # inner is already a compiled "rxode2" object -- load its DLL directly.
  if (is.null(mod))
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
                 theta_sens_cols = theta_sens_cols,
                 dummy_eta_inner = dummy_eta_inner,
                 rename_map = rename_map, is_lincmt = is_lincmt,
                 cache_file = .cacheFile)
  # suppressWarnings: the sens model is now compiled inside admixr2 (for
  # eventSens="jump"), so its environment chain references the package namespace
  # and serialising it warns "'package:admixr2' may not be available when
  # loading". Harmless -- the worker reloads the DLL via rxLoad(), not from the
  # serialised env. .admLoadModel() suppresses the same warning on its own cache.
  tryCatch(suppressWarnings(qs2::qs_save(result, .cacheFile)), error = function(e) NULL)
  result
}
