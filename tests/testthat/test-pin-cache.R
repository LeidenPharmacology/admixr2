# Tier 1: LRU model-pin cache in .adm_pin_env (R/zzz.R). Pure R, no rxode2.
# Guards the ubuntu-devel native-memory fix: the pin must round-trip, evict the
# least-recently-used entries past the limit, and keep a touched entry alive.

test_that(".admPinGet returns NULL for a missing key and round-trips a set", {
  on.exit(admixr2::admClearCache(), add = TRUE)
  admixr2::admClearCache()

  expect_null(admixr2:::.admPinGet("sens_absent"))

  admixr2:::.admPinSet("sens_a", list(x = 1L))
  expect_equal(admixr2:::.admPinGet("sens_a"), list(x = 1L))
})

test_that(".adm_pin_env LRU-evicts the oldest entries past the limit", {
  old <- options(admixr2.pin_limit = 3L)
  on.exit({ options(old); admixr2::admClearCache() }, add = TRUE)
  admixr2::admClearCache()

  expect_equal(admixr2:::.adm_pin_limit(), 3L)

  for (k in paste0("m", 1:5)) admixr2:::.admPinSet(k, k)   # insert m1..m5

  # Only the 3 most-recently-set survive; m1/m2 are evicted.
  expect_null(admixr2:::.admPinGet("m1"))
  expect_null(admixr2:::.admPinGet("m2"))
  expect_equal(admixr2:::.admPinGet("m3"), "m3")
  expect_equal(admixr2:::.admPinGet("m4"), "m4")
  expect_equal(admixr2:::.admPinGet("m5"), "m5")
})

test_that("a Get touch protects an entry from eviction (true LRU, not FIFO)", {
  old <- options(admixr2.pin_limit = 2L)
  on.exit({ options(old); admixr2::admClearCache() }, add = TRUE)
  admixr2::admClearCache()

  admixr2:::.admPinSet("a", 1L)
  admixr2:::.admPinSet("b", 2L)
  admixr2:::.admPinGet("a")            # touch a -> a is now most-recent
  admixr2:::.admPinSet("c", 3L)        # inserting c must evict b, not a

  expect_equal(admixr2:::.admPinGet("a"), 1L)
  expect_null(admixr2:::.admPinGet("b"))
  expect_equal(admixr2:::.admPinGet("c"), 3L)
})

test_that(".adm_pin_limit falls back to 16 on a bad option", {
  old <- options(admixr2.pin_limit = "nonsense")
  on.exit(options(old), add = TRUE)
  expect_equal(admixr2:::.adm_pin_limit(), 16L)

  options(admixr2.pin_limit = 0L)
  expect_equal(admixr2:::.adm_pin_limit(), 16L)
})

test_that(".admPinCompanions tracks the platform (Windows-only pin)", {
  expect_equal(admixr2:::.admPinCompanions(),
               identical(.Platform$OS.type, "windows"))
})

test_that("admClearCache resets the LRU order and counts only real entries", {
  admixr2::admClearCache()
  admixr2:::.admPinSet("sens_x", 1L)
  admixr2:::.admPinSet("sens_y", 2L)

  # Two real entries removed; the ".order" bookkeeping key is not counted.
  expect_equal(admixr2::admClearCache(), 2L)
  expect_null(admixr2:::.admPinGet("sens_x"))
})
