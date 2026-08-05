# .admWorkerLoadModels(): the loader a mirai daemon uses.
#
# In a normal fit this branch is unreachable: the parent hands the worker
# `rxMod_direct`/`sensModel_direct` in-process, so the cache-reading path only
# ever executes inside a daemon -- a separate process that resolves admixr2 from
# its OWN .libPaths(), which means neither the test suite nor covr observes a
# single line of it. It is directly callable, though, so it is tested here in
# process with rxMod_direct = NULL, which is the state a daemon is really in.
skip_on_cran()
skip_if_not_installed("rxode2")

.wl_model <- function() {
  ini({ tcl <- log(4); tv <- log(70); prop.err <- 0.1
        eta.cl ~ 0.09; eta.v ~ 0.04 })
  model({ cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
          d/dt(central) = -(cl / v) * central
          cp <- central / v
          cp ~ prop(prop.err) })
}

.wl_setup <- function() {
  ui <- rxode2::rxode2(.wl_model)
  sm <- admixr2:::.admLoadSensModel(ui)          # must run BEFORE .admLoadModel
  rx <- admixr2:::.admLoadModel(ui)
  pinfo <- admixr2:::.admParseIniDf(ui$iniDf, ui)
  pinfo$sim_cache_file <- admixr2:::.admModelCacheFile(ui)
  list(ui = ui, sm = sm, rx = rx, pinfo = pinfo)
}

test_that("the worker loads both models from the cache the parent wrote", {
  s <- .wl_setup()
  skip_if(is.null(s$sm), "sensitivity model unavailable")
  skip_if(!file.exists(s$pinfo$sim_cache_file), "simulation cache not written")

  m <- admixr2:::.admWorkerLoadModels(
    ui_lstExpr = s$ui$lstExpr, rxMod_direct = NULL, cores = 1L,
    sens_cache_file = s$sm$cache_file, sens_cols = s$sm$sens_cols,
    sens_rename = s$sm$rename_map, sensModel_direct = NULL, pinfo = s$pinfo)

  expect_type(m, "list")
  expect_s3_class(m$rxMod, "rxode2")
  expect_true(is.numeric(m$cores_w) && m$cores_w >= 1L)
  # the sens model came back, and carries the PARENT's derived fields rather
  # than whatever the cached file happened to hold
  expect_false(is.null(m$sensModel))
  expect_identical(m$sensModel$sens_cols, s$sm$sens_cols)
})

test_that("the worker refuses a simulation cache it cannot read", {
  s <- .wl_setup()
  bad <- s$pinfo
  bad$sim_cache_file <- file.path(tempdir(), "adm-sim-does-not-exist.rds")
  expect_error(
    admixr2:::.admWorkerLoadModels(
      ui_lstExpr = s$ui$lstExpr, rxMod_direct = NULL, cores = 1L,
      sens_cache_file = NULL, sens_cols = NULL, sens_rename = NULL,
      sensModel_direct = NULL, pinfo = bad),
    "could not read the compiled-model cache", fixed = TRUE)
})

test_that("the worker refuses a NULL cache path rather than guessing one", {
  # It used to fall back to the pre-.admIniKey formula, which cannot separate
  # fix(0.5) from fix(0.9) -- so a false hit solved at another model's fixed
  # value. Finding *a* file is the bug; the legible error is the fix. The message
  # must survive a NULL path, since `file.exists(NULL)` is logical(0) and would
  # otherwise replace it with "argument is of length zero".
  s <- .wl_setup()
  msg <- tryCatch({
    admixr2:::.admWorkerLoadModels(
      ui_lstExpr = s$ui$lstExpr, rxMod_direct = NULL, cores = 1L,
      sens_cache_file = NULL, sens_cols = NULL, sens_rename = NULL,
      sensModel_direct = NULL, pinfo = list())
    ""
  }, error = function(e) conditionMessage(e))
  expect_match(msg, "could not read the compiled-model cache")
  expect_false(grepl("length zero", msg, fixed = TRUE))
})

test_that("a worker with no sens model at all is quiet", {
  # sens_cache_file = NULL means the parent has no sensitivity model either, so
  # there is no divergence to report. Warning here would fire on every restart of
  # every gradient-free fit.
  s <- .wl_setup()
  skip_if(!file.exists(s$pinfo$sim_cache_file), "simulation cache not written")
  expect_no_warning(
    m <- admixr2:::.admWorkerLoadModels(
      ui_lstExpr = s$ui$lstExpr, rxMod_direct = NULL, cores = 1L,
      sens_cache_file = NULL, sens_cols = NULL, sens_rename = NULL,
      sensModel_direct = NULL, pinfo = s$pinfo))
  expect_null(m$sensModel)
})

test_that("a worker degrades to FD, and REPORTS it, when the sens cache is gone", {
  # Not an error: a worker without a sensitivity model still fits, by finite
  # differences. But it must SAY so -- it is then computing a different gradient
  # from the parent, which is invisible in the objective.
  #
  # Reported on the RETURN VALUE, not by warning(). This branch is reachable only
  # inside a mirai daemon (the sequential path passes sensModel_direct and
  # short-circuits), and mirai does not relay a daemon's conditions to the parent
  # -- so the warning this used to assert was raised exactly where nothing could
  # hear it. .admRunRestarts() collects the field and warns in the parent.
  s <- .wl_setup()
  skip_if(!file.exists(s$pinfo$sim_cache_file), "simulation cache not written")
  m <- admixr2:::.admWorkerLoadModels(
    ui_lstExpr = s$ui$lstExpr, rxMod_direct = NULL, cores = 1L,
    sens_cache_file = file.path(tempdir(), "adm-sens-gone.rds"),
    sens_cols = NULL, sens_rename = NULL, sensModel_direct = NULL,
    pinfo = s$pinfo)
  expect_null(m$sensModel)
  expect_s3_class(m$rxMod, "rxode2")
  expect_true(is.character(m$sens_fallback) && nzchar(m$sens_fallback))
  expect_match(m$sens_fallback, "sensitivity model cache")
})

test_that("rxMod_direct short-circuits the cache read entirely", {
  # The sequential path: the parent already holds the models, so no file is
  # touched. Passing a deliberately bogus cache path proves it is not read.
  s <- .wl_setup()
  bad <- s$pinfo
  bad$sim_cache_file <- file.path(tempdir(), "adm-sim-never-read.rds")
  m <- admixr2:::.admWorkerLoadModels(
    ui_lstExpr = s$ui$lstExpr, rxMod_direct = s$rx, cores = 1L,
    sens_cache_file = NULL, sens_cols = NULL, sens_rename = NULL,
    sensModel_direct = s$sm, pinfo = bad)
  expect_identical(m$rxMod, s$rx)
  expect_false(is.null(m$sensModel))
})

test_that("the worker refuses a cached model whose shared library is gone", {
  # rxLoad() does not error on a vanished DLL -- it silently re-runs the deferred
  # compile -- so the worker checks file.exists(rxDll()) itself and hard-fails.
  # Without that it would walk on with a live-looking model over an unloaded
  # library, and the first rxSolve() dereferences a dead pointer: an opaque
  # rxode2 error deep in the restart, or STATUS_HEAP_CORRUPTION on Windows.
  s <- .wl_setup()
  skip_if(!file.exists(s$pinfo$sim_cache_file), "simulation cache not written")

  fake <- file.path(tempdir(), "adm-sim-dll-gone.rds")
  mod  <- readRDS(s$pinfo$sim_cache_file)
  dll  <- tryCatch(rxode2::rxDll(if (inherits(mod, "rxode2")) mod else mod[[1L]]),
                   error = function(e) NA_character_)
  skip_if(is.na(dll), "could not resolve the cached model's DLL")
  saveRDS(mod, fake)
  on.exit(unlink(fake), add = TRUE)

  # Point the cached object at a DLL that does not exist, by copying it to a name
  # nothing built. The object still deserialises; only the file is missing.
  gone <- s$pinfo
  gone$sim_cache_file <- fake
  moved <- paste0(dll, ".moved")
  file.rename(dll, moved)
  on.exit(file.rename(moved, dll), add = TRUE)

  expect_error(
    admixr2:::.admWorkerLoadModels(
      ui_lstExpr = s$ui$lstExpr, rxMod_direct = NULL, cores = 1L,
      sens_cache_file = NULL, sens_cols = NULL, sens_rename = NULL,
      sensModel_direct = NULL, pinfo = gone),
    "could not (read|load)")
})

test_that("the worker REPORTS and degrades when the sens cache holds an unusable model", {
  # A file that exists but whose $mod is not an rxode2 object: the load must fail
  # INSIDE the tryCatch and become an FD gradient plus a recorded reason, not an
  # error that kills the restart.
  #
  # Recorded on the return value rather than warned -- see the cache-is-gone test
  # above for why a daemon's warning cannot reach the parent.
  s <- .wl_setup()
  skip_if(!file.exists(s$pinfo$sim_cache_file), "simulation cache not written")
  bad <- file.path(tempdir(), "adm-sens-not-a-model.rds")
  saveRDS(list(mod = 42, sens_cols = "rx_f1_ETA_1_"), bad)
  on.exit(unlink(bad), add = TRUE)

  m <- admixr2:::.admWorkerLoadModels(
    ui_lstExpr = s$ui$lstExpr, rxMod_direct = NULL, cores = 1L,
    sens_cache_file = bad, sens_cols = NULL, sens_rename = NULL,
    sensModel_direct = NULL, pinfo = s$pinfo)
  expect_null(m$sensModel)
  expect_s3_class(m$rxMod, "rxode2")
  expect_true(is.character(m$sens_fallback) && nzchar(m$sens_fallback))
  expect_match(m$sens_fallback, "could not load the sensitivity model")
})

test_that("the worker REFUSES a cache payload holding no compiled model", {
  # This test used to assert the opposite -- that shape is "not this function's
  # job" -- and that was wrong. .load_all() is upstream's iterate-the-container
  # step and is vacuously TRUE for anything not rxode2-classed, so a readable
  # .rds holding something else was not merely accepted, it was handed on AS
  # `rxMod`: the first rxSolve() then dereferences a non-model, which is an
  # opaque error deep in the restart or, on Windows, the STATUS_HEAP_CORRUPTION
  # crash the surrounding comment warns about.
  #
  # The parent asserts the same shape (.admLoadModel) and recovers by deleting
  # and rebuilding. A worker cannot rebuild -- it has no `ui`, which is the whole
  # reason the cache exists -- so its only correct move is the legible stop().
  s <- .wl_setup()
  odd <- file.path(tempdir(), "adm-sim-not-a-list.rds")
  saveRDS(42, odd)
  on.exit(unlink(odd), add = TRUE)
  p <- s$pinfo
  p$sim_cache_file <- odd
  expect_error(
    admixr2:::.admWorkerLoadModels(
      ui_lstExpr = s$ui$lstExpr, rxMod_direct = NULL, cores = 1L,
      sens_cache_file = NULL, sens_cols = NULL, sens_rename = NULL,
      sensModel_direct = NULL, pinfo = p),
    "could not read the compiled-model cache")
})
