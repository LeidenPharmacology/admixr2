# Session-scoped cache for once-per-session warnings.
# Keys are error-type strings; presence of a key means the warning was already emitted.
.adm_warn_env <- new.env(parent = emptyenv())

# Gauss-Hermite node cache (see .adghNodes1). The nodes depend on nothing but the
# node count, and a transformed endpoint asks for the 81-node set on every residual
# evaluation -- an 81x81 eigen() inside the objective's inner loop otherwise.
.adm_node_env <- new.env(parent = emptyenv())

# Session-scoped registry of compiled rxode2 models, keyed by model digest.
#
# Two reasons, and the second is the load-bearing one:
#
#  1. .admLoadModel() used to deserialise and re-load the model from disk on
#     EVERY call, even on a cache hit -- readRDS of a compiled model plus a
#     dyn.load, repeated for every fit and every test.
#  2. Each of those calls minted a NEW rxode2 object wrapping the same DLL, and
#     every one of them is a finalizer waiting to run. rxode2's finalizers unload
#     the model's shared library, so a gc() at an unlucky moment can unload a
#     library while native code is still using it. Holding the objects here makes
#     them permanently reachable, so they never become collectable at all.
#
# Bounded and cleared WHOLESALE at .ADM_MODEL_CACHE_MAX entries, which is the
# mechanism nlmixr2est uses for its own compiled-model session cache
# (.foceiAnalyticAugCache in foceiCovAnalytic.R: get0()/assign() on an
# emptyenv-parented env, a composite key covering everything that changes the
# emitted model, and `if (length(ls(...)) >= 64L) rm(list = ls(...), ...)` to
# bound retained compiled models). Copied rather than invented so the two
# packages behave the same way under the same pressure.
#
# NOTE the trade-off that bound carries here: a wipe makes the evicted models
# collectable again, so a session fitting more than this many DISTINCT models
# re-enters the finalizer window this cache exists to avoid. 64 is upstream's
# number and a suite that registered 51 stayed under it.
.ADM_MODEL_CACHE_MAX <- 64L
.adm_model_env <- new.env(parent = emptyenv())

# Session cache for compiled SENSITIVITY models, same mechanism. Separate env
# per purpose, as upstream keeps .foceiAnalyticAugCache and .foceiAnalyticEtCache
# apart, so one being wiped does not take the other with it.
.adm_sens_env <- new.env(parent = emptyenv())

# Bound-and-wipe, upstream's shape, in one place so the two caches cannot drift.
.admCacheAssign <- function(key, value, envir) {
  if (length(ls(envir, all.names = TRUE)) >= .ADM_MODEL_CACHE_MAX)
    rm(list = ls(envir, all.names = TRUE), envir = envir)
  assign(key, value, envir = envir)
  value
}

# Write a compiled-model cache entry: neither swallowed nor fatal.
#
# SWALLOWED (`tryCatch(saveRDS(...), error = function(e) NULL)`, which is what
# 0.4.0 did) is wrong because the parallel restart workers FIND these models by
# reading exactly these files. A silent failure here surfaced much later, and
# somewhere else, as every restart failing to read the cache.
#
# FATAL (a bare saveRDS, matching upstream's bare `qs2::qs_save`) is wrong too,
# and is what this replaces. The cache is an OPTIMISATION -- by the time it is
# written the model is already compiled and loaded -- so on an unwritable or full
# rxTempDir() (a locked-down HPC home, a cache directory owned by another user) a
# fit that used to work, slowly, died with `cannot open the connection` from
# inside model loading. The sens model was worse: its callers wrap the whole build
# in tryCatch(error = function(e) NULL), so a WRITE failure discarded a
# successfully compiled model and silently dropped adfo from the order-2 analytic
# structural gradient to forward FD.
#
# So: warn once per file, and carry on with the model already in hand. The next
# fit recompiles, which is exactly the pre-cache behaviour.
#
# suppressWarnings: these models are compiled inside admixr2, so their environment
# chain references the package namespace and serialising warns "'package:admixr2'
# may not be available when loading". Harmless -- a worker reloads the DLL via
# rxLoad(), not from the serialised environment.
.admCacheWrite <- function(object, file, what) {
  # Sampled before the write: afterwards a failed-at-open call and a
  # failed-halfway call look identical, and only the second one made the mess.
  .pre_existing <- file.exists(file)
  ok <- tryCatch({ suppressWarnings(saveRDS(object, file)); TRUE },
                 error = function(e) {
                   .adm_warn_once(paste0("cache_write:", file), sprintf(
                     paste0("admixr2: could not write the %s cache to '%s' (%s). ",
                            "The fit continues, but the model will be recompiled ",
                            "for every fit in this session, and parallel restarts ",
                            "(workers > 1) cannot read it."),
                     what, file, conditionMessage(e)))
                   FALSE
                 })
  # A half-written file is worse than none: a later session would read it back as
  # a corrupt entry rather than a miss.
  #
  # But only clean up a file THIS call created. `ok` is FALSE for any saveRDS
  # failure, including one that fails at open time and so leaves a complete,
  # valid pre-existing entry untouched -- and the caches are shared, keyed only
  # by model content, across sessions and processes. Session A writes
  # adm-sim-<hash>.rds and starts daemons about to read it; session B fits the
  # same model, its write fails (disk full, or the file momentarily open on
  # Windows), and an unconditional remove deletes A's entry. A's workers take the
  # deliberately fallback-free `.cacheFile <- pinfo$sim_cache_file` branch, find
  # nothing and stop() -- surfacing as "parallel restart 1 failed" on a fit that
  # was running correctly.
  #
  # existed-before is sampled BEFORE the write, because after a partial write the
  # two cases are indistinguishable from the file alone.
  # suppressWarnings, not just tryCatch(error=): file.remove() signals a WARNING
  # on failure, not an error, so an error handler alone lets "cannot remove file
  # ..., reason 'Permission denied'" leak out right behind the warning above --
  # two warnings for one event, the second of them noise.
  if (!ok && !.pre_existing && file.exists(file))
    tryCatch(suppressWarnings(file.remove(file)), error = function(e) NULL)
  invisible(ok)
}

# Drop the disk cache when the package version changes.
#
# nlmixr2est::.resetCacheIfNeeded(), with admixr2's version in place of its
# build-time md5 constant. Same mechanism, same place in the load sequence:
# stamp our identity into rxTempDir() and act when the stored stamp differs.
#
# IT DELETES NOTHING, and the comment here used to claim otherwise ("admixr2 had
# no equivalent, so superseded compiled models accumulated there indefinitely").
# They still do. Upstream can honestly say there is nothing to sweep -- its
# focei-*.rds are keyed on model identity alone -- whereas .admPkgKey() puts a
# SOURCE DIGEST in the sensitivity key, so every edit to an emitter orphans a
# fresh adm-sens-*.rds.
#
# A targeted sweep of our own files on a stamp mismatch was considered and
# REJECTED. The accumulation it would recover is bounded for a user: the digest
# is constant within an installed version, so it is one set of files per version,
# and the stamp changes exactly when those become dead. It is unbounded only
# while DEVELOPING the emitters, where rxode2::rxClean() is already at hand. The
# risk is not symmetric: a mirai daemon runs library(admixr2), so this hook fires
# in the WORKER, and under the documented dev-mode version skew (daemon on a
# different installed build) a mismatching daemon would delete the .rds the
# parent wrote seconds earlier and its siblings are about to read. That trades
# unbounded disk for a fit that fails -- a bad trade for a hygiene fix. Sweep
# manually with rxode2::rxClean() instead.
.admResetCacheIfNeeded <- function() {
  .wd <- rxode2::rxTempDir()
  if (.wd != "") {
    .verFile <- file.path(.wd, "admixr2.version")
    .ver <- as.character(utils::packageVersion("admixr2"))
    if (file.exists(.verFile)) {
      ## n = 1L / warn = FALSE, and a length check, because `readLines(f) != .ver`
      ## is not a scalar condition: an EMPTY stamp gives character(0) -> logical(0)
      ## -> "argument is of length zero", and a multi-line one gives a vector,
      ## which R >= 4.2 also errors on. .onLoad() wraps this in tryCatch(), so
      ## either would be SWALLOWED -- the stamp would then never be refreshed and
      ## the whole check would silently do nothing on every subsequent load. A
      ## truncated or hand-edited stamp is exactly the state where it most needs
      ## to work, so anything unexpected counts as a mismatch and is overwritten.
      .stamp <- tryCatch(readLines(.verFile, n = 1L, warn = FALSE),
                         error = function(e) character(0))
      if (length(.stamp) != 1L || !identical(.stamp, .ver)) {
        ## Deliberately NOT rxClean() here -- upstream removed exactly this call
        ## (nlmixr2est fef5be69). It wipes the whole SHARED rxTempDir(), and
        ## deleting a generated model's compiled artifact out from under a live
        ## model object makes rxode2's deferred-compile thunk rebuild it,
        ## emitting different code than the build it replaces (nlmixr2/rxode2#1171).
        ## For admixr2 it is worse still: a mirai daemon runs library(admixr2),
        ## so this hook fires in the WORKER and can wipe the adm-sim-*.rds the
        ## parent wrote seconds earlier and the sibling daemons are about to read.
        ##
        ## Rewriting the stamp is also what stops the comparison failing forever:
        ## the mismatch branch never refreshed it, so every load for the rest of
        ## the installation's life would clean again, not just the first after an
        ## upgrade.
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
