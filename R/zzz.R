# Session-scoped cache for once-per-session warnings.
# Keys are error-type strings; presence of a key means the warning was already emitted.
.adm_warn_env <- new.env(parent = emptyenv())

# Gauss-Hermite node cache (see .adghNodes1).
.adm_node_env <- new.env(parent = emptyenv())

# Session-scoped registry of compiled rxode2 models, keyed by model digest.
# Holds compiled models permanently reachable to prevent GC finalizers from
# unloading shared libraries while native code is executing. Bounded and
# cleared wholesale at .ADM_MODEL_CACHE_MAX entries (matching nlmixr2est).
.ADM_MODEL_CACHE_MAX <- 64L
.adm_model_env <- new.env(parent = emptyenv())

# Session cache for compiled sensitivity models (separate env so wipes do not cross).
.adm_sens_env <- new.env(parent = emptyenv())

# Bound-and-wipe cache assignment in one place.
.admCacheAssign <- function(key, value, envir) {
  if (length(ls(envir, all.names = TRUE)) >= .ADM_MODEL_CACHE_MAX)
    rm(list = ls(envir, all.names = TRUE), envir = envir)
  assign(key, value, envir = envir)
  value
}

# Write a compiled-model cache entry: warn once on failure and carry on with
# the model in memory (non-fatal optimisation).
.admCacheWrite <- function(object, file, what) {
  # Write to unique tempfile in same directory, then rename atomically to prevent
  # concurrent readers (e.g. parallel workers) seeing 0-byte or truncated files.
  # tempfile() uses an internal counter + pid without perturbing R's RNG stream.
  .tmp <- tryCatch(
    tempfile(pattern = paste0(basename(file), ".tmp"), tmpdir = dirname(file)),
    error = function(e) NULL)
  .fail <- function(e) {
    .adm_warn_once(paste0("cache_write:", file), sprintf(
      paste0("admixr2: could not write the %s cache to '%s' (%s). ",
             "The fit continues, but the model will be recompiled ",
             "for every fit in this session, and parallel restarts ",
             "(workers > 1) cannot read it."),
      what, file, conditionMessage(e)))
    FALSE
  }
  ok <- if (is.null(.tmp)) .fail(simpleError("no writable temporary name")) else
    tryCatch({ suppressWarnings(saveRDS(object, .tmp)); TRUE }, error = .fail)
  if (ok) {
    # file.rename() returns FALSE/warns on Windows when target is open for reading.
    # Leaving existing content-addressed entry in place is safe.
    ok <- isTRUE(tryCatch(suppressWarnings(file.rename(.tmp, file)),
                          error = function(e) FALSE))
    if (!ok)
      .adm_warn_once(paste0("cache_write:", file), sprintf(
        paste0("admixr2: could not publish the %s cache entry '%s' ",
               "(another process is most likely reading it). The fit continues ",
               "using the model already in hand."), what, file))
  }
  # Clean up temp file
  if (!is.null(.tmp) && file.exists(.tmp))
    tryCatch(suppressWarnings(file.remove(.tmp)), error = function(e) NULL)
  invisible(ok)
}

# Update the version stamp in rxTempDir() on package load.
# Does not delete files on mismatch (avoids wiping models mid-run or across workers).
.admResetCacheIfNeeded <- function() {
  .wd <- rxode2::rxTempDir()
  if (.wd != "") {
    .verFile <- file.path(.wd, "admixr2.version")
    .ver <- as.character(utils::packageVersion("admixr2"))
    if (file.exists(.verFile)) {
      # Validate stamp format; any read failure or non-scalar counts as mismatch
      .stamp <- tryCatch(readLines(.verFile, n = 1L, warn = FALSE),
                         error = function(e) character(0))
      if (length(.stamp) != 1L || !identical(.stamp, .ver)) {
        # Deliberately avoid rxClean() which wipes shared rxTempDir() under live models
        writeLines(.ver, .verFile)
      }
    } else {
      writeLines(.ver, .verFile)
    }
  }
}

.onLoad <- function(libname, pkgname) {
  tryCatch(.admResetCacheIfNeeded(), error = function(e) NULL)
  tryCatch(.register_adm(),  error = function(e)
    warning("admixr2: admc registration failed (", conditionMessage(e), ")", call. = FALSE))
  tryCatch(.register_adirmc(), error = function(e)
    warning("admixr2: adirmc registration failed (", conditionMessage(e), ")", call. = FALSE))
  tryCatch(.register_adfo(), error = function(e)
    warning("admixr2: adfo registration failed (", conditionMessage(e), ")", call. = FALSE))
  tryCatch(.register_adgh(), error = function(e)
    warning("admixr2: adgh registration failed (", conditionMessage(e), ")", call. = FALSE))
  # Register knit_print methods into knitr's namespace (knitr is in Suggests).
  # If knitr loads after admixr2 the setHook fires and registers then.
  tryCatch(.register_knit_print(), error = function(e) NULL)
  setHook(packageEvent("knitr", "onLoad"),
          function(...) tryCatch(.register_knit_print(), error = function(e) NULL))
}

.register_adm <- function() {
  ns <- asNamespace("nlmixr2est")
  registerS3method("nlmixr2Est",              "admc",        nlmixr2Est.admc,                    envir = ns)
  registerS3method("getValidNlmixrCtl",       "admc",        getValidNlmixrCtl.admc,             envir = ns)
  registerS3method("nmObjGetControl",         "admc",        nmObjGetControl.admc,               envir = ns)
  registerS3method("nmObjHandleControlObject","admControl",  nmObjHandleControlObject.admControl, envir = ns)
}

.register_knit_print <- function() {
  if (!isNamespaceLoaded("knitr")) return(invisible(NULL))
  ns <- asNamespace("knitr")
  registerS3method("knit_print", "admFit", .admKnitPrint, envir = ns)
}

.register_adirmc <- function() {
  ns <- asNamespace("nlmixr2est")
  registerS3method("nlmixr2Est",              "adirmc",        nlmixr2Est.adirmc,                        envir = ns)
  registerS3method("getValidNlmixrCtl",       "adirmc",        getValidNlmixrCtl.adirmc,                 envir = ns)
  registerS3method("nmObjGetControl",         "adirmc",        nmObjGetControl.adirmc,                   envir = ns)
  registerS3method("nmObjHandleControlObject","adirmcControl", nmObjHandleControlObject.adirmcControl,   envir = ns)
}

.register_adfo <- function() {
  ns <- asNamespace("nlmixr2est")
  registerS3method("nlmixr2Est",              "adfo",        nlmixr2Est.adfo,                      envir = ns)
  registerS3method("getValidNlmixrCtl",       "adfo",        getValidNlmixrCtl.adfo,               envir = ns)
  registerS3method("nmObjGetControl",         "adfo",        nmObjGetControl.adfo,                 envir = ns)
  registerS3method("nmObjHandleControlObject","adfoControl", nmObjHandleControlObject.adfoControl, envir = ns)
}

.register_adgh <- function() {
  ns <- asNamespace("nlmixr2est")
  registerS3method("nlmixr2Est",              "adgh",        nlmixr2Est.adgh,                      envir = ns)
  registerS3method("getValidNlmixrCtl",       "adgh",        getValidNlmixrCtl.adgh,               envir = ns)
  registerS3method("nmObjGetControl",         "adgh",        nmObjGetControl.adgh,                 envir = ns)
  registerS3method("nmObjHandleControlObject","adghControl", nmObjHandleControlObject.adghControl, envir = ns)
}
