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

# --- .admPkgKey(): the emitter list must stay in step with the code ------------

test_that(".admPkgKey digests every emitter it names, and they all exist", {
  # A name that cannot be resolved degrades to the name itself rather than
  # collapsing the whole digest -- but it must not come to that, so pin the list.
  for (nm in admixr2:::.ADM_SENS_EMITTERS)
    expect_true(is.function(get(nm, envir = asNamespace("admixr2"))),
                info = nm)

  k <- admixr2:::.admPkgKey()
  expect_type(k, "character")
  expect_length(k, 1L)
  # version/source-digest, both non-empty
  parts <- strsplit(k, "/", fixed = TRUE)[[1L]]
  expect_length(parts, 2L)
  expect_true(all(nzchar(parts)))
  # the digest half must NOT be the "everything failed" sentinel
  expect_false(identical(parts[2L], "NA"))
  expect_identical(k, admixr2:::.admPkgKey())   # stable within a session
})

test_that(".admPkgKey covers the emitters that shape the cached payload", {
  # Regression on the gap this list closed: digesting only the two entry points
  # left .admSensFromInner/.admLinCmtToOde/.admRxode2/.admModName/.admJumpCovers
  # able to change what is cached without changing the key.
  expect_true(all(c(".admBuildThetaSens", ".admLoadSensModel", ".admSensFromInner",
                    ".admLinCmtToOde", ".admRxode2", ".admModName",
                    ".admJumpCovers") %in% admixr2:::.ADM_SENS_EMITTERS))
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
