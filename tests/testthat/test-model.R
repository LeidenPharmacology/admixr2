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

test_that(".admPkgKey tracks the installed package version", {
  # Keying the version is what replaced a hand-maintained schema-tag string,
  # which had to be edited whenever the emitted model changed -- and was not.
  expect_true(nzchar(admixr2:::.admPkgKey()))
})
