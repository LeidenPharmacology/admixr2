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
  rm(list = nms, envir = .adm_pin_env)
  qs2_files <- list.files(rxode2::rxTempDir(),
                          pattern = "^adm-.*\\.qs2$", full.names = TRUE)
  unlink(qs2_files)
  invisible(length(nms))
}

.onLoad <- function(libname, pkgname) {
  tryCatch(.register_adm(),  error = function(e)
    warning("admixr2: admc registration failed (", conditionMessage(e), ")", call. = FALSE))
  tryCatch(.register_adirmc(), error = function(e)
    warning("admixr2: adirmc registration failed (", conditionMessage(e), ")", call. = FALSE))
  tryCatch(.register_adfo(), error = function(e)
    warning("admixr2: adfo registration failed (", conditionMessage(e), ")", call. = FALSE))
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
