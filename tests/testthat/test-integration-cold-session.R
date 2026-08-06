# Regression test for the #81 first-fit recursion.
#
# Accessing ui$simulationModel during model compilation leaves a self-referential
# rxode2 object in ui$meta$.simModelBase. With covMethod = "r", nlmixr2's
# fit-assembly deep-clone of the ui (.cloneEnv, no cycle detection) recursed
# forever -- "node stack overflow" / "evaluation nested too deeply" -- aborting
# the fit. It only triggers when THIS fit compiles the model (a cache miss, so
# .simModelBase is (re)created); a fit that hits the compiled-model cache never
# touches $simulationModel and is safe.
#
# The regular pipeline (test-integration-pipeline.R) builds several
# covMethod = "none" fits first, which cache the model, so its later
# covMethod = "r" fit hits the cache and would NOT reproduce the crash. Dropping
# THIS model's simulation-cache entry forces the recompile that repopulates
# .simModelBase, so the test reproduces the exact trigger independent of test
# execution order.
#
# It used to do that with rxode2::rxClean(), which does not scope to this test:
# it empties the whole shared rxTempDir(), so every model any other test compiled
# is gone and every later file pays to rebuild it. Harmless when files run one
# after another, fatal when they do not -- under a parallel testthat run this
# deleted the cache entry a mirai daemon in test-integration-parallel.R had
# already been promised, and that fit died with "a parallel worker could not read
# the compiled-model cache". Deleting one content-addressed entry, for the model
# this test actually fits, reproduces the same cache miss and touches nobody
# else's: test-integration-parallel.R fits one_cmt_kappa_fn, a different digest.

test_that("covMethod='r' after a cache clear does not recurse in fit assembly (#81)", {
  skip_on_cran()
  skip_if_not_installed("rxode2")
  skip_if_not_installed("nlmixr2est")
  nlmixr2 <- nlmixr2est::nlmixr2

  times  <- c(0.5, 1, 2, 4)
  E_true <- .one_cmt_mean(5, 20, 100, times)
  study1 <- list(E = E_true, V = diag((0.3 * E_true)^2), n = 200L,
                 times = times, ev = rxode2::et(amt = 100))

  # Force .admLoadModel() to recompile -> repopulates ui$meta$.simModelBase.
  # Scoped to one_cmt_fn's own entry; see the header for why not rxClean().
  .ui <- suppressMessages(tryCatch(rxode2::rxode2(one_cmt_fn), error = function(e) NULL))
  skip_if(is.null(.ui), "model would not parse")
  .cache <- admixr2:::.admModelCacheFile(.ui)
  unlink(.cache)
  # The whole test rests on this fit being the one that COMPILES. A cache hit
  # exercises a different path and would pass while proving nothing, silently --
  # which is precisely how the trigger was lost before. Assert the miss.
  expect_false(file.exists(.cache))

  fit <- suppressMessages(nlmixr2(
    one_cmt_fn, admData(), est = "admc",
    control = admControl(studies = list(s1 = study1), n_sim = 300L,
                         maxeval = 12L, seed = 1L, grad = "sens",
                         covMethod = "r", cov_n_sim = 2000L)))

  expect_s3_class(fit, "admFit")
  expect_true(is.finite(fit$objective))
  expect_identical(fit$env$covMethod, "r")

  # The transient artifact must not be left on the returned fit's ui either.
  .meta <- fit$env$ui$meta
  if (is.environment(.meta)) {
    .rx_left <- vapply(ls(.meta, all.names = TRUE), function(nm) {
      v <- tryCatch(get(nm, envir = .meta, inherits = FALSE), error = function(e) NULL)
      is.environment(v) && inherits(v, "rxode2")
    }, logical(1))
    expect_false(any(.rx_left))
  }
})
