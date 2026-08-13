# A daemon's own startup can delete the model cache it is about to read.
#
# `library(admixr2)` in a daemon loads nlmixr2est, and the INSTALLED 6.2.0's
# .resetCacheIfNeeded() calls rxode2::rxClean() whenever the nlmixr2est.md5 stamp
# in rxTempDir() does not match its build -- which wipes the WHOLE shared cache
# directory, including the adm-sim-*.rds the parent wrote seconds earlier. The
# mismatch branch never rewrites the stamp, so it is not self-limiting: it
# re-fires on every load, in every daemon. Having two nlmixr2est builds in play
# (a pkgload::load_all() source tree next to the installed package) is enough,
# and that is an ordinary state during upstream development.
#
# .admWarmDaemons() is the answer: get every daemon's startup -- and therefore
# every daemon-side rxClean() -- out of the way FIRST, then re-derive anything
# that went missing before a single restart is dispatched.
skip_on_cran()
skip_if_not_installed("rxode2")
skip_if_not_installed("mirai")

test_that(".admWarmDaemons rebuilds a model cache deleted after the parent wrote it", {
  env <- .int_grad_setup()
  ui  <- env$ui
  pinfo <- env$pinfo
  pinfo$sim_cache_file <- admixr2:::.admModelCacheFile(ui)

  admixr2:::.admLoadSensModel(ui)          # INVARIANT: sens before sim
  admixr2:::.admLoadModel(ui)
  skip_if(!file.exists(pinfo$sim_cache_file), "simulation cache not written")

  # Stand in for the daemon-side rxClean(): the entry the workers were told to
  # read is simply gone by the time they look.
  unlink(pinfo$sim_cache_file)
  expect_false(file.exists(pinfo$sim_cache_file))

  suppressMessages(admixr2:::.admWarmDaemons(ui, pinfo))

  # Rebuilt, and a real compiled model rather than any old file -- the worker
  # asserts exactly this shape before using it.
  expect_true(file.exists(pinfo$sim_cache_file))
  expect_s3_class(readRDS(pinfo$sim_cache_file), "rxode2")
})

test_that(".admWarmDaemons leaves an intact cache alone", {
  # The common case must stay free: no recompile, no rewrite, when nothing was
  # cleared. Pinned by mtime rather than by content, since a rebuild would
  # republish an equivalent payload and be invisible otherwise.
  env <- .int_grad_setup()
  ui  <- env$ui
  pinfo <- env$pinfo
  pinfo$sim_cache_file <- admixr2:::.admModelCacheFile(ui)
  admixr2:::.admLoadSensModel(ui)
  admixr2:::.admLoadModel(ui)
  skip_if(!file.exists(pinfo$sim_cache_file), "simulation cache not written")

  before <- file.mtime(pinfo$sim_cache_file)
  Sys.sleep(1.1)                            # coarse mtime resolution on some filesystems
  suppressMessages(admixr2:::.admWarmDaemons(ui, pinfo))
  expect_identical(file.mtime(pinfo$sim_cache_file), before)
})
