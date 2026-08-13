# The cross-session cache guard.
#
# `.rds` cache entries live in rxTempDir(), which can persist; the compiled
# artifacts they name are built in .admModDir(), which is under tempdir() and is
# session-local. A cache entry written by another session therefore names an
# artifact this session must not use -- and R removes its temp directory only on
# a CLEAN exit, so a killed session leaves one behind whose DLL passes
# file.exists() indefinitely.
#
# Serving such an entry is silent: rxLoad() does not error, and the model solves
# with its prediction frozen at the t = 0 value and NA structural gradients.
#
# This cannot be produced by ordinary single-session test code, so the foreign
# session is simulated by building into a different directory and then restoring
# .admModDir() -- which leaves exactly the artifact-exists-but-is-not-ours state.
skip_on_cran()
skip_if_not_installed("rxode2")

.xs_model <- function() {
  ini({ tcl <- log(4); tv <- log(70); prop.err <- 0.1
        eta.cl ~ 0.09; eta.v ~ 0.04 })
  model({ cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
          d/dt(central) = -(cl / v) * central
          cp <- central / v
          cp ~ prop(prop.err) })
}

test_that(".admRxLoadAll rejects an artifact built by another session", {
  ui <- rxode2::rxode2(.xs_model)

  # Build one model as if by a different session: same *Sens directory NAME (so
  # the discriminator has something to match on) under a different SESSION.
  #
  # OUTSIDE tempdir(), not merely outside .admModDir(). A real foreign session's
  # build directory is a sibling of ours under the system temp root -- it cannot
  # be nested inside our own tempdir, which is created for us and removed on our
  # exit. Putting the fixture at <tempdir>/xsOtherSession/admixr2Sens instead made
  # it foreign only to a discriminator that compares against .admModDir(), and
  # that comparison is precisely what could never accept nlmixr2est's own
  # <tempdir>/nlmixr2estSens model either. The guard now asks the question that
  # actually identifies a session -- see .admUnderTemp() -- so the fixture has to
  # ask it honestly.
  foreign_root <- file.path(dirname(tempdir()),
                            paste0("admxsOther", Sys.getpid()))
  foreign_dir  <- file.path(foreign_root, "admixr2Sens")
  dir.create(foreign_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(foreign_root, recursive = TRUE, force = TRUE), add = TRUE)
  skip_if(!dir.exists(foreign_dir), "could not create a foreign build directory")
  orig <- admixr2:::.admModDir
  utils::assignInNamespace(".admModDir", function() foreign_dir, ns = "admixr2")
  foreign <- tryCatch(admixr2:::.admBuildThetaSens(ui, character(0)),
                      error = function(e) NULL)
  utils::assignInNamespace(".admModDir", orig, ns = "admixr2")
  skip_if(is.null(foreign), "sensitivity model could not be built")

  dll <- rxode2::rxDll(foreign$mod)
  # The artifact must genuinely EXIST -- otherwise file.exists() alone would
  # reject it and the directory half of the guard would go untested.
  expect_true(file.exists(dll))
  expect_false(admixr2:::.admSameDir(dirname(dirname(dll)),
                                     admixr2:::.admModDir()))
  expect_false(admixr2:::.admUnderTemp(dll))

  # The guard must refuse it.
  expect_false(admixr2:::.admRxLoadAll(foreign$mod))

  # ... and must still accept a model this session built.
  ours <- tryCatch(admixr2:::.admBuildThetaSens(ui, character(0)),
                   error = function(e) NULL)
  skip_if(is.null(ours), "sensitivity model could not be rebuilt")
  expect_true(admixr2:::.admRxLoadAll(ours$mod))
})

test_that("a session-local build directory that is not ours is still accepted", {
  # The other half of the invariant, and the one the previous discriminator got
  # wrong: "session-local *Sens directory" is not the same set as ".admModDir()".
  # nlmixr2est 7.x builds its own inner model in <tempdir>/nlmixr2estSens, which
  # matches the name pattern and can NEVER equal <tempdir>/admixr2Sens -- so a
  # cached .admSensFromInner() result was rejected on every call, in the very
  # session that built it, and recompiled (~3 s) for every fit.
  #
  # Asserted on the predicate rather than by standing up a real nlmixr2est model:
  # this must hold on 6.x (where the inner model is elsewhere and the pattern does
  # not fire) as well as on 7.x.
  nlm <- file.path(tempdir(), "nlmixr2estSens", "mod.d", "mod_x64.dll")
  expect_true(admixr2:::.admUnderTemp(nlm))
  expect_false(admixr2:::.admSameDir(dirname(dirname(nlm)),
                                     admixr2:::.admModDir()))
  # ... and our own directory is of course still ours.
  expect_true(admixr2:::.admUnderTemp(admixr2:::.admModDir()))
})

test_that("the discriminator survives Windows 8.3 short paths", {
  # THE regression. rxDll() returns a path whose DIRECTORY components are in 8.3
  # short form ("...\\Temp\\RT4F27~1\\ADMIXR~1\\admSens_...dll"), so matching the
  # raw string against "admixr2Sens" never fires and the guard silently does
  # nothing at all. Caught only by a genuine two-process test; pinned here.
  d <- admixr2:::.admModDir()
  short <- tryCatch(utils::shortPathName(d), error = function(e) d)
  skip_if(identical(short, d), "no 8.3 short name on this platform")

  expect_false(grepl("admixr2Sens", short, fixed = TRUE))   # the trap
  expect_true(grepl("admixr2Sens",
                    normalizePath(short, winslash = "/", mustWork = FALSE),
                    fixed = TRUE))                           # the fix
  # and .admSameDir must equate the two spellings
  expect_true(admixr2:::.admSameDir(short, d))
})
