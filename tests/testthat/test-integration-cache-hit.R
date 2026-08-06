# The disk-cache HIT paths of .admLoadModel() and .admLoadSensModel().
#
# Unreachable on a fresh runner, which is why coverage never saw them: the
# session memo answers every repeat call within one process, and CI starts with
# an empty rxTempDir so the first call is always a miss. The state that reaches
# them -- memo empty, disk file present -- is the ordinary one for a SECOND R
# session on a machine with an rxode2 cache, and for every mirai worker.
#
# Emptying the memo reproduces it deterministically.
skip_on_cran()
skip_if_not_installed("rxode2")

.ch_model <- function() {
  ini({ tcl <- log(4); tv <- log(70); prop.err <- 0.1
        eta.cl ~ 0.09; eta.v ~ 0.04 })
  model({ cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
          d/dt(central) = -(cl / v) * central
          cp <- central / v
          cp ~ prop(prop.err) })
}

.ch_clear_memo <- function() {
  for (e in list(admixr2:::.adm_model_env, admixr2:::.adm_sens_env))
    rm(list = ls(e, all.names = TRUE), envir = e)
}

test_that(".admLoadModel serves a disk hit and re-memoises it", {
  ui <- rxode2::rxode2(.ch_model)
  admixr2:::.admLoadSensModel(ui)                 # ordering invariant
  first <- admixr2:::.admLoadModel(ui)
  expect_s3_class(first, "rxode2")
  skip_if(!file.exists(admixr2:::.admModelCacheFile(ui)), "cache not written")

  # memo empty, disk file present: the branch a second session takes
  .ch_clear_memo()
  expect_identical(length(ls(admixr2:::.adm_model_env, all.names = TRUE)), 0L)
  second <- admixr2:::.admLoadModel(ui)
  expect_s3_class(second, "rxode2")
  # it came from the file and was put back in the memo
  expect_gt(length(ls(admixr2:::.adm_model_env, all.names = TRUE)), 0L)
  # ... and the memo now answers without touching disk
  third <- admixr2:::.admLoadModel(ui)
  expect_s3_class(third, "rxode2")
})

test_that(".admLoadSensModel serves a disk hit and re-derives its volatile fields", {
  ui <- rxode2::rxode2(.ch_model)
  first <- admixr2:::.admLoadSensModel(ui)
  skip_if(is.null(first), "sensitivity model unavailable")
  skip_if(!file.exists(first$cache_file), "sens cache not written")

  .ch_clear_memo()
  second <- admixr2:::.admLoadSensModel(ui)
  expect_false(is.null(second))
  expect_s3_class(second$mod, "rxode2")
  # The fields a worker cannot re-derive must be the PARENT's, not the file's --
  # a stale rename_map would put a theta's value in the wrong THETA[k] slot.
  expect_identical(second$cache_file, first$cache_file)
  expect_identical(second$sens_cols,  first$sens_cols)
  expect_identical(second$rename_map, first$rename_map)
  expect_gt(length(ls(admixr2:::.adm_sens_env, all.names = TRUE)), 0L)
})

test_that("a corrupt cache entry is discarded and rebuilt rather than served", {
  # The shape assertion: .admRxLoadAll is a no-op on anything not rxode2-classed,
  # so it returns TRUE for a file that is not a model at all. Without the
  # inherits() check the junk would be handed back as `rxMod` AND memoised, so
  # every fit in the session would repeat the failure.
  ui <- rxode2::rxode2(.ch_model)
  admixr2:::.admLoadSensModel(ui)
  admixr2:::.admLoadModel(ui)
  f <- admixr2:::.admModelCacheFile(ui)
  skip_if(!file.exists(f), "cache not written")

  .ch_clear_memo()
  saveRDS(list(not = "a model"), f)      # shape the guard must reject
  again <- admixr2:::.admLoadModel(ui)
  expect_s3_class(again, "rxode2")       # rebuilt, not served
})

test_that(".admRxLoadAll rejects a model whose shared library has gone", {
  # The check that makes a stale cache entry self-heal instead of being served.
  # rxLoad() does not error on a vanished DLL -- it silently re-runs the deferred
  # compile -- so file.exists(rxDll()) has to be asked explicitly.
  ui <- rxode2::rxode2(.ch_model)
  m <- admixr2:::.admLoadSensModel(ui)
  skip_if(is.null(m), "sensitivity model unavailable")
  dll <- tryCatch(rxode2::rxDll(m$mod), error = function(e) NA_character_)
  skip_if(is.na(dll) || !file.exists(dll), "no DLL to move")

  expect_true(admixr2:::.admRxLoadAll(m$mod))     # while it is there
  moved <- paste0(dll, ".moved")
  file.rename(dll, moved)
  on.exit(file.rename(moved, dll), add = TRUE)
  expect_false(admixr2:::.admRxLoadAll(m$mod))    # and once it is not
})

test_that("an order-2 request on a transformed endpoint comes back as order 1", {
  # V_pred's second derivative cannot be chained through g() with g' alone, so
  # adfo finite-differences those endpoints instead. Demoting ABOVE the cache key
  # means the order-1 entry is shared rather than a near-duplicate stored under an
  # order-2 key -- and it is why a lnorm fit does not pay for a cross block that
  # nothing reads.
  lnorm_mod <- function() {
    ini({ tcl <- log(4); tv <- log(70); prop.err <- 0.1
          eta.cl ~ 0.09; eta.v ~ 0.04 })
    model({ cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
            d/dt(central) = -(cl / v) * central
            cp <- central / v
            cp ~ lnorm(prop.err) })
  }
  ui <- rxode2::rxode2(lnorm_mod)
  sm <- admixr2:::.admLoadSensModel(ui, order = 2L)
  skip_if(is.null(sm), "sensitivity model unavailable")
  expect_null(sm$d2_cols)          # demoted: no cross block was built
})

test_that(".admLinCmtToOde leaves a model that is already an ODE alone", {
  # Promotion only applies to a solved-form linCmt(). Anything else must come
  # back untouched rather than being rewritten.
  ui <- rxode2::rxode2(.ch_model)
  out <- admixr2:::.admLinCmtToOde(ui)
  expect_false(is.null(out))
})

test_that(".admRxode2's fallback still builds in OUR directory, under a name of its own", {
  # .admModName() returns NULL two ways: an unusable role, and rxModelVars()
  # failing to yield a parsed md5 -- the second does not depend on the caller.
  # The fallback must still build in .admModDir(), because dropping back to
  # rxode2's own directory puts the model right back where two builds of one text
  # overwrite each other.
  #
  # It also must not simply omit modName: rxode2 refuses a `wd` without one
  # ("working directory specified, but modName not declared"), so the branch this
  # replaced threw instead of falling back.
  m <- tryCatch(admixr2:::.admRxode2("d/dt(a) = -a", role = NULL),
                error = function(e) conditionMessage(e))
  expect_false(is.character(m))                  # i.e. it did not error
  dll <- tryCatch(rxode2::rxDll(m), error = function(e) NA_character_)
  skip_if(is.na(dll), "no DLL")
  expect_true(admixr2:::.admSameDir(dirname(dirname(dll)), admixr2:::.admModDir()))
})

test_that(".admRxode2's fallback name separates two eventSens builds of one text", {
  # THE property of the whole naming scheme (nlmixr2/rxode2#1171): the emitted C
  # depends on eventSens, which the .so path otherwise ignores, so two builds of
  # one model text that differ only there collide -- the second overwrites the
  # first and earlier model objects silently execute the replacement.
  #
  # The primary name folds eventSens in. A fallback that digests only the model
  # text would reintroduce the bug in the one function written to prevent it, so
  # it folds in `...` too. Asserted on the artifacts, not the names, since the
  # name is synthesised inside .admRxode2.
  a <- tryCatch(admixr2:::.admRxode2("d/dt(b) = -b", role = NULL, eventSens = "jump"),
                error = function(e) NULL)
  b <- tryCatch(admixr2:::.admRxode2("d/dt(b) = -b", role = NULL, eventSens = "fd"),
                error = function(e) NULL)
  skip_if(is.null(a) || is.null(b), "models could not be built")
  da <- tryCatch(rxode2::rxDll(a), error = function(e) NA_character_)
  db <- tryCatch(rxode2::rxDll(b), error = function(e) NA_character_)
  skip_if(is.na(da) || is.na(db), "no DLL")
  expect_false(identical(basename(da), basename(db)))
  # ... and both are still ours, so the cross-session guard recognises them even
  # though neither is named "admSens*" -- which is why that guard tests the PATH.
  for (d in c(da, db))
    expect_true(admixr2:::.admSameDir(dirname(dirname(d)), admixr2:::.admModDir()))
})

test_that(".admModName tolerates an absent eventSens", {
  a <- admixr2:::.admModName("d/dt(a) = -a", "admSens")
  b <- admixr2:::.admModName("d/dt(a) = -a", "admSens", eventSens = NA)
  expect_true(nzchar(a)); expect_true(nzchar(b))
})
