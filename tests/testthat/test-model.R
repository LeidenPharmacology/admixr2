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
