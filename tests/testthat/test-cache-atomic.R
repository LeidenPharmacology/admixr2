# .admCacheWrite() publishes a cache entry atomically.
#
# These entries are content-addressed files in a SHARED, persistent rxTempDir(),
# so their readers are other processes: another session fitting the same model
# derives the same path, and a parallel fit's own daemons read the entry the
# parent just wrote. An in-place saveRDS() publishes the path at ZERO BYTES from
# the moment it opens the connection, and a reader in that window gets a
# truncated payload -- which is how the intermittent
#   "a parallel worker could not read the compiled-model cache"
# failure was produced. Truncation at any fraction breaks readRDS(), so there is
# no "mostly written is good enough" here.

.ca_tmpfiles <- function(dir, base)
  list.files(dir, pattern = paste0("^", base, "\\.tmp"), all.files = TRUE)

test_that("a successful cache write publishes the object and leaves no temp file", {
  d <- file.path(tempdir(), "ca-ok"); dir.create(d, showWarnings = FALSE)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  f <- file.path(d, "adm-sim-ok.rds")

  expect_true(admixr2:::.admCacheWrite(list(v = "NEW", pad = seq_len(1000L)),
                                       f, "simulation model"))
  expect_identical(readRDS(f)$v, "NEW")
  # the sibling temp is ours alone and must never survive the call
  expect_identical(.ca_tmpfiles(d, "adm-sim-ok\\.rds"), character(0))
})

test_that("a cache write that cannot happen warns, reports FALSE, and creates nothing", {
  # A directory that does not exist: the temp write fails at open, exactly as an
  # unwritable or full rxTempDir() would. The fit must continue -- the model is
  # already compiled and in hand by the time this is called -- so this is a
  # warning and a FALSE, never an error.
  missing_dir <- file.path(tempdir(), "ca-absent-dir")
  unlink(missing_dir, recursive = TRUE)
  f <- file.path(missing_dir, "adm-sim-nowhere.rds")

  expect_warning(
    expect_false(admixr2:::.admCacheWrite(list(v = 1), f, "simulation model")),
    "could not write the simulation model cache")
  expect_false(file.exists(f))
})

test_that("the target is never opened for writing, so a reader's entry survives", {
  # THE regression test. Windows refuses to rename over a file that is open for
  # reading ("Access is denied"), which is exactly a concurrent reader -- so this
  # exercises the failure branch and asserts the property that matters: the
  # previous COMPLETE entry is still there, byte for byte.
  #
  # With the in-place saveRDS() this replaced, the same call overwrote the open
  # file's contents (measured) -- and a reader mid-readRDS() saw the prefix of a
  # different object.
  #
  # POSIX allows the rename and publishes the new entry atomically instead, so
  # the branch is unreachable there; the invariant it guards ("a reader sees one
  # complete entry or the other, never a prefix") is checked for both below.
  skip_if(.Platform$OS.type != "windows", "rename-over-open is Windows-specific")

  d <- file.path(tempdir(), "ca-open"); dir.create(d, showWarnings = FALSE)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  f <- file.path(d, "adm-sim-open.rds")
  saveRDS(list(v = "OLD", pad = seq_len(1000L)), f)
  before <- readBin(f, "raw", file.size(f))

  h <- file(f, "rb")
  res <- withCallingHandlers(
    admixr2:::.admCacheWrite(list(v = "NEW"), f, "simulation model"),
    warning = function(w) invokeRestart("muffleWarning"))
  close(h)

  expect_false(res)
  expect_identical(readBin(f, "raw", file.size(f)), before)
  expect_identical(readRDS(f)$v, "OLD")
  expect_identical(.ca_tmpfiles(d, "adm-sim-open\\.rds"), character(0))
})

test_that("after a cache write the entry is always a complete, readable object", {
  # The cross-platform statement of the same invariant: whatever the outcome, the
  # published path holds a whole object. Windows keeps OLD (the rename is
  # refused); POSIX publishes NEW. Never a prefix of either.
  d <- file.path(tempdir(), "ca-complete"); dir.create(d, showWarnings = FALSE)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  f <- file.path(d, "adm-sim-complete.rds")
  saveRDS(list(v = "OLD"), f)

  h <- file(f, "rb")
  suppressWarnings(admixr2:::.admCacheWrite(list(v = "NEW"), f, "simulation model"))
  close(h)

  expect_true(readRDS(f)$v %in% c("OLD", "NEW"))
})

test_that("writing a cache entry does not disturb the RNG stream", {
  # The temp name must not come from sample()/runif(): every estimator draws its
  # Sobol/rnorm z matrices from this stream, so consuming from it here would move
  # the objective of any fit that happens to compile a model.
  d <- file.path(tempdir(), "ca-rng"); dir.create(d, showWarnings = FALSE)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  set.seed(42); expected <- runif(3)
  set.seed(42)
  admixr2:::.admCacheWrite(list(v = 1), file.path(d, "adm-sim-rng.rds"), "simulation model")
  expect_identical(runif(3), expected)
})
