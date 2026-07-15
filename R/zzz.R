# Session-scoped environment for pinning foceiModel companion objects and
# caching sens model results across calls within a session.
#
# rxode2 companion objects ($outer, $predOnly, $predNoLhs) hold live C++ DLL
# pointers. If they become GC-eligible while another allocation is in progress,
# their finalizers unload DLLs mid-allocation -> STATUS_HEAP_CORRUPTION on
# Windows. We cannot pin them to `ui` because rxode2 locks the rxUi environment
# (assign() fails silently there). A package-level env is always writable.
#
# Keys: paste0("focei_", digest(ui$lstExpr))  -> foceiModel list (companion pin)
#       paste0("sens_",  digest(ui$lstExpr))  -> sens model result (in-memory cache)
.adm_pin_env <- new.env(parent = emptyenv())

# Whether to pin the foceiModel companion objects ($outer/$predOnly/$predNoLhs).
# They are NEVER read back -- the pin exists only to keep them off the GC so their
# rxode2 finalizers cannot unload a DLL mid-allocation (STATUS_HEAP_CORRUPTION,
# a Windows-only failure). On other platforms the finalizers are safe, so pinning
# them is pure overhead: it holds every distinct model's companion DLLs resident
# for the whole session AND blocks rxode2's own DLL unloading (.rxShouldUnload
# sees a live reference), so native memory climbs without bound across many models
# -- the ubuntu-devel test-suite hang. Skip it off Windows and let rxode2 reclaim.
.admPinCompanions <- function() identical(.Platform$OS.type, "windows")

# Upper bound on entries kept in .adm_pin_env, LRU-evicted beyond it -- the same
# shape as rxode2's rxSolveCacheLimit / rxUnloadAll orphan ceiling, which is why
# rxode2 and nlmixr2est stay flat while fitting many models. Two entries per model
# on Windows (focei_ + sens_), one elsewhere. Generous default so a normal
# multi-model session never evicts; caps pathological long-running sessions.
# Override with options(admixr2.pin_limit = N).
.adm_pin_limit <- function() {
  .v <- suppressWarnings(as.integer(getOption("admixr2.pin_limit", 16L)))
  if (length(.v) != 1L || is.na(.v) || .v < 1L) 16L else .v
}

# LRU bookkeeping lives in .adm_pin_env$.order (most-recent first). It is not a
# model key, so it is never itself a Get/Set target or an eviction candidate.
.admPinTouch <- function(key) {
  .ord <- .adm_pin_env$.order
  .adm_pin_env$.order <- c(key, .ord[.ord != key])
}

.admPinGet <- function(key) {
  if (!exists(key, envir = .adm_pin_env, inherits = FALSE)) return(NULL)
  .admPinTouch(key)
  get(key, envir = .adm_pin_env, inherits = FALSE)
}

.admPinSet <- function(key, value) {
  assign(key, value, envir = .adm_pin_env)
  .admPinTouch(key)
  .ord <- .adm_pin_env$.order
  .lim <- .adm_pin_limit()
  if (length(.ord) > .lim) {
    .drop <- .ord[-seq_len(.lim)]
    if (length(.drop)) {
      rm(list = .drop, envir = .adm_pin_env)
      # Drop the evicted models' native DLLs at this quiescent point (no solve in
      # flight) so any rxode2 finalizer runs here, controlled, rather than firing
      # mid-allocation later -- the Windows heap-corruption trigger the pin guards.
      gc(FALSE)
    }
    .adm_pin_env$.order <- .ord[seq_len(.lim)]
  }
  invisible(value)
}

# Session-scoped cache for once-per-session warnings.
# Keys are error-type strings; presence of a key means the warning was already emitted.
.adm_warn_env <- new.env(parent = emptyenv())

#' Clear the admixr2 model cache
#'
#' Removes all cached simulation and sensitivity models from both the
#' session-level in-memory cache and the qs2 disk files written to
#' `rxode2::rxTempDir()`. Call this in long-running sessions to free memory
#' and disk space after fitting many distinct models.
#'
#' @return Invisibly returns the number of in-memory objects removed.
#' @export
admClearCache <- function() {
  nms <- ls(envir = .adm_pin_env, all.names = TRUE)
  rm(list = nms, envir = .adm_pin_env)   # also drops the ".order" LRU bookkeeping
  qs2_files <- list.files(rxode2::rxTempDir(),
                          pattern = "^adm-.*\\.qs2$", full.names = TRUE)
  unlink(qs2_files)
  invisible(length(setdiff(nms, ".order")))   # count only real cache entries
}

.onLoad <- function(libname, pkgname) {
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
