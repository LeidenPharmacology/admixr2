# Session-scoped cache for once-per-session warnings.
# Keys are error-type strings; presence of a key means the warning was already emitted.
.adm_warn_env <- new.env(parent = emptyenv())

# Gauss-Hermite node cache (see .adghNodes1). The nodes depend on nothing but the
# node count, and a transformed endpoint asks for the 81-node set on every residual
# evaluation -- an 81x81 eigen() inside the objective's inner loop otherwise.
.adm_node_env <- new.env(parent = emptyenv())

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
