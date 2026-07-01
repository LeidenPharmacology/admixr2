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
# covMethod = "r" fit hits the cache and would NOT reproduce the crash.
# admClearCache() forces the recompile that repopulates .simModelBase, so this
# test reproduces the exact trigger independent of test execution order.

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
  admClearCache()

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
