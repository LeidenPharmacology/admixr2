# Tier 1 unit tests for .admDropSimModelMeta() (the #81 recursion fix).
# No rxode2 required: a ui is mocked as a list with a plain `meta` environment.

test_that(".admDropSimModelMeta drops rxode2 envs from ui$meta, keeps everything else", {
  meta <- new.env(parent = emptyenv())

  rx_model <- new.env(parent = emptyenv())
  class(rx_model) <- c("rxode2tos", "rxode2")   # the cyclic artifact we must remove
  plain_env <- new.env(parent = emptyenv())      # unrelated env: must be kept
  class(plain_env) <- "someOtherClass"

  assign(".simModelBase", rx_model,  envir = meta)
  assign(".keepEnv",      plain_env, envir = meta)
  assign(".keepScalar",   42L,       envir = meta)

  admixr2:::.admDropSimModelMeta(list(meta = meta))

  expect_false(exists(".simModelBase", envir = meta, inherits = FALSE))
  expect_true(exists(".keepEnv",    envir = meta, inherits = FALSE))
  expect_true(exists(".keepScalar", envir = meta, inherits = FALSE))
})

test_that(".admDropSimModelMeta is an invisible no-op when ui$meta is not an environment", {
  expect_invisible(admixr2:::.admDropSimModelMeta(list(meta = NULL)))
  expect_invisible(admixr2:::.admDropSimModelMeta(list()))
})

# ---- cache keys (Tier 1: a ui is mocked as a list holding an iniDf) ----------
#
# The sensitivity cache carries a fix()ed theta's VALUE to the solve as data
# (sensModel$fixed_theta), and a parallel worker reads that file without a `ui`
# to re-derive from. So two models differing only in a fixed value must not
# share a key -- otherwise every restart solves at the other fit's fixed value,
# silently, and across sessions because the cache directory is persistent.

.mock_ini_ui <- function(est, names = c("tka", "tcl", "add.err"),
                         fix = c(TRUE, FALSE, FALSE),
                         err = c(NA, NA, "add")) {
  list(iniDf = data.frame(name = names, est = est, fix = fix, err = err,
                          stringsAsFactors = FALSE))
}

test_that(".admIniKey separates models that differ only in a FIXED value", {
  a <- .mock_ini_ui(c(0.5, log(5), 0.3))
  b <- .mock_ini_ui(c(0.9, log(5), 0.3))
  expect_false(identical(admixr2:::.admIniKey(a), admixr2:::.admIniKey(b)))
})

test_that(".admIniKey ignores a changed STARTING value", {
  # A starting value is optimizer state, not model text: invalidating on it
  # would force a recompile for every tweak and buy nothing.
  a <- .mock_ini_ui(c(0.5, log(5), 0.3))
  b <- .mock_ini_ui(c(0.5, log(7), 0.3))
  expect_identical(admixr2:::.admIniKey(a), admixr2:::.admIniKey(b))
})

test_that(".admIniKey separates a reordered ini, a changed fix flag, and a changed err", {
  base <- .mock_ini_ui(c(0.5, log(5), 0.3))
  expect_false(identical(
    admixr2:::.admIniKey(base),
    admixr2:::.admIniKey(.mock_ini_ui(c(log(5), 0.5, 0.3),
                                      names = c("tcl", "tka", "add.err")))))
  expect_false(identical(
    admixr2:::.admIniKey(base),
    admixr2:::.admIniKey(.mock_ini_ui(c(0.5, log(5), 0.3),
                                      fix = c(TRUE, TRUE, FALSE)))))
  expect_false(identical(
    admixr2:::.admIniKey(base),
    admixr2:::.admIniKey(.mock_ini_ui(c(0.5, log(5), 0.3),
                                      err = c(NA, NA, "prop")))))
})

test_that(".admIniKey survives a ui with no iniDf", {
  # A mock ui (Tier-1 tests, a hand-built pinfo) has no iniDf. It must return a
  # single stable string rather than erroring or leaking a stray token.
  k <- admixr2:::.admIniKey(list())
  expect_type(k, "character")
  expect_length(k, 1L)
  expect_false(grepl("NULL", k, fixed = TRUE))
  expect_identical(k, admixr2:::.admIniKey(list()))
})

test_that(".admPkgKey carries the package version AND the emitter's source", {
  # Keying the version is what replaced a hand-maintained schema-tag string,
  # which had to be edited whenever the emitted model changed -- and was not.
  # The version alone only moves at RELEASE, so every commit in between shares
  # one key against a cache that persists across sessions: editing what the
  # order-2 build emits and re-running would hit the previously compiled model
  # and contract its columns against the new direction map. Digesting the
  # emitter's own body closes that with nothing to remember.
  k <- admixr2:::.admPkgKey()
  expect_true(nzchar(k))
  expect_identical(k, admixr2:::.admPkgKey())        # stable within a session
  parts <- strsplit(k, "/", fixed = TRUE)[[1L]]
  expect_length(parts, 2L)
  expect_identical(parts[[1L]], as.character(utils::packageVersion("admixr2")))
  # ... and the second component really is a digest of the emitter, not a
  # constant: recomputing it from a DIFFERENT body must give something else.
  .other <- digest::digest(list(deparse(body(admixr2:::.admIniKey)),
                                deparse(body(admixr2:::.admLoadSensModel))))
  expect_false(identical(parts[[2L]], .other))
})

# ---- simulation-model cache key ---------------------------------------------
#
# The same fix()ed-value collision, on the OTHER cache. A fixed theta is not an
# estimated parameter, so .admMakeParamsList() builds no column for it and the
# value the solve uses is the one rxode2 baked into $simulationModel. Keyed on
# the model({}) block alone, two models differing only in `tka <- fix(0.5)` vs
# `fix(0.9)` therefore shared one compiled model and the second silently solved
# at the first's fixed value.

.mock_sim_ui <- function(lst, est, fix = c(TRUE, FALSE, FALSE)) {
  c(.mock_ini_ui(est, fix = fix), list(lstExpr = lst))
}

test_that(".admModelCacheFile separates models that differ only in a FIXED value", {
  skip_if_not_installed("rxode2")
  lst <- list(quote(cl <- exp(tcl)), quote(cp ~ add(add.err)))
  a <- admixr2:::.admModelCacheFile(.mock_sim_ui(lst, c(0.5, log(5), 0.3)))
  b <- admixr2:::.admModelCacheFile(.mock_sim_ui(lst, c(0.9, log(5), 0.3)))
  expect_false(identical(a, b))
  # Both still land in the rxode2 temp dir with the expected prefix, so the
  # worker's fallback formula and this one name files of the same shape.
  expect_true(all(grepl("^adm-sim-.*[.]rds$", basename(c(a, b)))))
})

test_that(".admModelCacheFile is stable and ignores a changed STARTING value", {
  skip_if_not_installed("rxode2")
  lst <- list(quote(cl <- exp(tcl)), quote(cp ~ add(add.err)))
  a <- admixr2:::.admModelCacheFile(.mock_sim_ui(lst, c(0.5, log(5), 0.3)))
  expect_identical(a, admixr2:::.admModelCacheFile(.mock_sim_ui(lst, c(0.5, log(5), 0.3))))
  # A starting value must NOT force a recompile -- same rule as .admIniKey.
  expect_identical(a, admixr2:::.admModelCacheFile(.mock_sim_ui(lst, c(0.5, log(7), 0.3))))
  # ... but a different model({}) block must.
  expect_false(identical(a, admixr2:::.admModelCacheFile(
    .mock_sim_ui(list(quote(cl <- exp(tcl) * 2), quote(cp ~ add(add.err))),
                 c(0.5, log(5), 0.3)))))
})

# --- .admSameDir(): the one guard whose failure mode is endless recompiles -----
#
# .admRxLoadAll() rejects a cached model whose DLL sits in another session's
# build directory. If .admSameDir() ever answered FALSE for two spellings of the
# SAME directory, every load would be rejected and every fit would recompile
# forever -- slow, not wrong, and so easy to miss. Windows is where that would
# happen: rxode2 hands back 8.3 short components ("RT6EC4~1") while tempdir()
# reports the long form, and the two are the same directory.

test_that(".admSameDir equates spellings of one directory and separates two", {
  d <- normalizePath(tempdir(), winslash = "/", mustWork = FALSE)

  expect_true(admixr2:::.admSameDir(d, d))
  expect_true(admixr2:::.admSameDir(paste0(d, "/"), d))   # trailing separator

  # The backslash spelling and the 8.3 short form are WINDOWS-ONLY spellings of a
  # path. On Linux and macOS a backslash is an ordinary character in a file name,
  # so "\tmp\RtmpX" names a different (non-existent) thing and .admSameDir is
  # RIGHT to say FALSE -- asserting otherwise fails there and only there, which
  # is how this first reached CI green on Windows and red on the other two.
  if (.Platform$OS.type == "windows") {
    expect_true(admixr2:::.admSameDir(gsub("/", "\\\\", d), d))
    short <- tryCatch(utils::shortPathName(d), error = function(e) d)
    expect_true(admixr2:::.admSameDir(short, d))
  }

  # ... and it must still say FALSE for a genuinely different directory, or the
  # guard stops guarding.
  expect_false(admixr2:::.admSameDir(file.path(d, "someOtherSession"), d))
  expect_false(admixr2:::.admSameDir("/no/such/path/at/all", d))
})

test_that(".admSameDir on a vanished directory is FALSE, not an error", {
  # The case the guard exists for: a cached DLL under a killed session's tempdir.
  # normalizePath(mustWork = FALSE) must not throw, and must not accidentally
  # compare equal to ours.
  gone <- file.path(tempdir(), "adm-no-such-dir-12345")
  expect_false(admixr2:::.admSameDir(gone, tempdir()))
})

# --- .admPkgKey(): the emitter set must stay in step with the code ------------

test_that(".admPkgKey is a version and a real digest, stable within a session", {
  k <- admixr2:::.admPkgKey()
  expect_type(k, "character")
  expect_length(k, 1L)
  parts <- strsplit(k, "/", fixed = TRUE)[[1L]]
  expect_length(parts, 2L)
  expect_identical(parts[[1L]], as.character(utils::packageVersion("admixr2")))
  expect_true(nzchar(parts[[2L]]))
  # "NA" is the sentinel the digest degrades to when nothing could be resolved.
  # Seeing it here would mean cache invalidation is silently switched off.
  expect_false(identical(parts[[2L]], "NA"))
  expect_identical(k, admixr2:::.admPkgKey())
})

test_that("editing ANY emitter changes .admPkgKey", {
  # The property, asserted directly rather than by inspecting a list of names --
  # the point of the key is that an edit to the code that decides what gets
  # cached invalidates the cache, with nothing to remember. Digesting only the
  # two entry points left .admSensFromInner/.admLinCmtToOde/.admRxode2/
  # .admModName/.admJumpCovers free to change what is emitted without changing
  # the key, which is a cache HIT on a model built by superseded code: a finite,
  # plausible, silently wrong gradient under a normal-looking objective.
  ns <- asNamespace("admixr2")
  base <- admixr2:::.admPkgKey()
  emitters <- c(".admBuildThetaSens", ".admLoadSensModel", ".admSensNameMaps",
                ".admSensFromInner", ".admLinCmtToOde", ".admRxode2",
                ".admModName", ".admJumpCovers", ".admToRx", ".admEndpointVar",
                ".admCountSpec", ".admBetaSpec", ".admBetaPair", ".admDistArgs",
                ".admIsLinCmtMod")
  for (nm in emitters) {
    orig <- get(nm, envir = ns)
    patched <- orig
    # prepend a no-op: different source text, identical behaviour
    body(patched) <- as.call(c(as.name("{"), quote(.adm_probe <- 1L),
                               as.list(body(orig))[-1L]))
    utils::assignInNamespace(nm, patched, ns = "admixr2")
    moved <- !identical(admixr2:::.admPkgKey(), base)
    utils::assignInNamespace(nm, orig, ns = "admixr2")
    expect_true(moved, info = nm)
  }
  # ... and everything is restored, so the key is what it was
  expect_identical(admixr2:::.admPkgKey(), base)
})

test_that("no function an emitter CALLS escapes the key unaccounted for", {
  # The test above only checks the names it is handed, which is precisely how the
  # count/beta helpers were missed: they were added to the emitter PATH long after
  # the list was written, and an edit to .admEndpointVar was demonstrated serving
  # a cached model whose rx_pred_ was 3x the prediction the objective scored --
  # same session, no error.
  #
  # So walk the call graph instead. Every admixr2 function reachable from a listed
  # emitter must be either listed itself, or on the allowlist below with a stated
  # reason. A new helper in the emitter path fails this until someone decides
  # which it is, which is the decision that was skipped.
  ns <- asNamespace("admixr2")
  emitters <- c(".admBuildThetaSens", ".admLoadSensModel", ".admSensNameMaps",
                ".admSensFromInner", ".admLinCmtToOde", ".admRxode2",
                ".admModName", ".admJumpCovers", ".admToRx", ".admEndpointVar",
                ".admCountSpec", ".admBetaSpec", ".admBetaPair", ".admDistArgs",
                ".admIsLinCmtMod")

  # Reachable but deliberately NOT digested:
  allow <- c(
    # covered BY VALUE -- their result is itself a component of the cache key, so
    # changing the body moves the key through the value
    ".admIniKey", ".admUnpairedThetas", ".admMuRefPairs", ".admNameOccurrence",
    # the key itself; it cannot digest its own body
    ".admPkgKey",
    # cache PLUMBING -- decides where a payload is stored or whether a stored one
    # is usable, never what the payload contains
    ".admCacheAssign", ".admCacheWrite", ".admRxLoadAll", ".admModDir",
    ".admModelCacheFile", ".admNormPath", ".admSameDir", ".admUnderTemp",
    ".admDropSimModelMeta", ".adm_warn_once")

  called <- function(nm) {
    o <- tryCatch(get(nm, envir = ns), error = function(e) NULL)
    if (!is.function(o)) return(character(0))
    b <- paste(deparse(body(o)), collapse = " ")
    u <- unique(unlist(regmatches(b, gregexpr("\\.adm[A-Za-z0-9_]*", b))))
    u[vapply(u, function(x)
      exists(x, envir = ns, inherits = FALSE) && is.function(get(x, envir = ns)),
      logical(1))]
  }
  reach <- unique(unlist(lapply(emitters, called)))
  escaped <- setdiff(reach, c(emitters, allow))
  expect_identical(escaped, character(0),
                   info = paste("reachable from an emitter but neither digested",
                                "nor allowlisted:", paste(escaped, collapse = ", ")))
})

test_that(".admResetCacheIfNeeded survives a malformed version stamp", {
  # `readLines(f) != .ver` is not a scalar condition: an EMPTY stamp yields
  # logical(0) ("argument is of length zero") and a multi-line one a vector,
  # which R >= 4.2 also errors on. .onLoad() wraps this in tryCatch(), so either
  # would be swallowed and the stamp never refreshed -- the check would then
  # silently never run again. Anything unexpected must count as a mismatch.
  wd <- rxode2::rxTempDir()
  skip_if(!nzchar(wd) || !dir.exists(wd), "no rxode2 temp dir")
  f <- file.path(wd, "admixr2.version")
  keep <- if (file.exists(f)) readLines(f, warn = FALSE) else NULL
  on.exit({
    if (is.null(keep)) unlink(f) else writeLines(keep, f)
  }, add = TRUE)

  for (content in list(character(0), c("a", "b"), "")) {
    writeLines(content, f)
    expect_silent(admixr2:::.admResetCacheIfNeeded())
    # ... and the stamp is repaired to this version, so the next load is a match
    expect_identical(readLines(f, n = 1L, warn = FALSE),
                     as.character(utils::packageVersion("admixr2")))
  }
})

# --- artifact naming guards (Tier 1: pure functions, no rxode2 needed) --------

test_that(".admModName refuses a name it cannot make unique", {
  # Returning NULL is not a failure: .admRxode2() then falls back to rxode2's own
  # anonymous naming, which is still built in OUR directory. A malformed name
  # would be worse than no name, since two models could share it.
  expect_null(admixr2:::.admModName("d/dt(a) = -a", role = NULL))
  expect_null(admixr2:::.admModName("d/dt(a) = -a", role = ""))
})

test_that(".admModName folds the role and eventSens in beside the md5", {
  a <- admixr2:::.admModName("d/dt(a) = -a", "admSens", eventSens = "jump")
  b <- admixr2:::.admModName("d/dt(a) = -a", "admSens", eventSens = "fd")
  c3 <- admixr2:::.admModName("d/dt(a) = -a", "admSensInner", eventSens = "jump")
  d <- admixr2:::.admModName("d/dt(a) = -2*a", "admSens", eventSens = "jump")
  # eventSens is the whole point (nlmixr2/rxode2#1171): two builds of one text
  # that differ only there must not share an artifact.
  expect_false(identical(a, b))
  expect_false(identical(a, c3))   # ... nor two roles
  expect_false(identical(a, d))    # ... nor two model texts
  expect_true(all(nzchar(c(a, b, c3, d))))
})

test_that(".admRxLoadAll is TRUE for a payload holding no compiled model", {
  # The cached object is a list; anything in it that is not an rxode2 model is
  # not this function's business, and a non-list payload is vacuously loadable.
  expect_true(admixr2:::.admRxLoadAll(42))
  expect_true(admixr2:::.admRxLoadAll("not a model"))
  expect_true(admixr2:::.admRxLoadAll(list(a = 1, b = "x")))
})
