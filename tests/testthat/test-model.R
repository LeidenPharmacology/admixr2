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
