# Tier 1: per-fit teardown (.admFitTeardown, R/utils.R). Pure R -- no rxode2 fit.
# The heavy end-to-end effect (bounded RSS across the integration suite) is
# exercised by the Tier-2 integration tests running many fits in one process;
# here we check the unit behaviour.

# Helper: drop any test keys we injected into the rxode2 registry / tracker.
.tw_cleanup <- function(reg, keys) {
  suppressWarnings(rm(list = intersect(keys, ls(reg, all.names = TRUE)), envir = reg))
  suppressWarnings(rm(list = intersect(keys, ls(admixr2:::.adm_model_keys, all.names = TRUE)),
                      envir = admixr2:::.adm_model_keys))
}

test_that(".admClearPins clears the in-memory model pins", {
  on.exit(admixr2::admClearCache(), add = TRUE)
  assign("sens_dummy",  list(x = 1L), envir = admixr2:::.adm_pin_env)
  assign("focei_dummy", list(y = 2L), envir = admixr2:::.adm_pin_env)
  expect_gte(length(ls(admixr2:::.adm_pin_env, all.names = TRUE)), 2L)

  admixr2:::.admClearPins()

  expect_equal(length(ls(admixr2:::.adm_pin_env, all.names = TRUE)), 0L)
})

test_that(".admFitTeardown removes this fit's registry additions, keeps user models", {
  reg <- admixr2:::.admRxRegistry()
  skip_if(is.null(reg), "rxode2 model registry unreachable")
  keys <- c("adm_user_model", "adm_fit_model")
  on.exit(.tw_cleanup(reg, keys), add = TRUE)
  options(admixr2.fit_teardown = TRUE)

  assign("adm_user_model", 0L, envir = reg)        # a pre-existing (user) model
  before <- admixr2:::.admRegistrySnapshot()
  assign("adm_fit_model", 0L, envir = reg)         # "registered by this fit"

  admixr2:::.admFitTeardown(before)

  expect_true(exists("adm_user_model", envir = reg, inherits = FALSE))   # preserved
  expect_false(exists("adm_fit_model", envir = reg, inherits = FALSE))   # reclaimed
})

test_that(".admFitTeardown reclaims TRACKED models even without a fit snapshot", {
  reg <- admixr2:::.admRxRegistry()
  skip_if(is.null(reg), "rxode2 model registry unreachable")
  keys <- c("adm_user_model2", "adm_tracked_model")
  on.exit(.tw_cleanup(reg, keys), add = TRUE)
  options(admixr2.fit_teardown = TRUE)

  assign("adm_user_model2", 0L, envir = reg)        # user model (never tracked)
  before <- admixr2:::.admRegistrySnapshot()
  assign("adm_tracked_model", 0L, envir = reg)      # loaded "by admixr2" outside a fit
  admixr2:::.admTrackRegistry(before)               # record it as admixr2's

  admixr2:::.admFitTeardown(NULL)                    # no fit delta -> clears tracked only

  expect_true(exists("adm_user_model2",  envir = reg, inherits = FALSE))   # preserved
  expect_false(exists("adm_tracked_model", envir = reg, inherits = FALSE)) # reclaimed
})

test_that(".admFitTeardown honours the disable toggle", {
  on.exit({ options(admixr2.fit_teardown = TRUE); admixr2::admClearCache() }, add = TRUE)
  assign("sens_keep", list(x = 1L), envir = admixr2:::.adm_pin_env)
  options(admixr2.fit_teardown = FALSE)

  admixr2:::.admFitTeardown(character(0))

  expect_true(exists("sens_keep", envir = admixr2:::.adm_pin_env, inherits = FALSE))
})

test_that(".admFitTeardown clears pins off Windows, keeps them on Windows", {
  on.exit(admixr2::admClearCache(), add = TRUE)
  options(admixr2.fit_teardown = TRUE)
  assign("sens_dummy", list(x = 1L), envir = admixr2:::.adm_pin_env)

  admixr2:::.admFitTeardown(admixr2:::.admRegistrySnapshot())

  kept <- exists("sens_dummy", envir = admixr2:::.adm_pin_env, inherits = FALSE)
  if (identical(.Platform$OS.type, "windows"))
    expect_true(kept)    # pins anchor GC finalizers on Windows -> not cleared here
  else
    expect_false(kept)   # off Windows the pins are dropped for R-heap hygiene
})
