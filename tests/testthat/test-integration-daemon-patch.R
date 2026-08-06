# Dev-mode patching of a daemon's namespace.
#
# Under devtools::load_all() the parent serialises its dev functions and the
# daemon patches them into the admixr2 it loaded from its OWN .libPaths() -- the
# INSTALLED one, whose namespace is locked. utils::assignInNamespace() can REPLACE
# a binding there but cannot ADD one, and .admDaemonRestart() swallows the error,
# so a function the branch INTRODUCES is simply absent in the worker while
# everything looks healthy.
#
# The failure mode is not always a clean error: .ADM_SENS_EMITTERS was a new
# top-level constant, and its absence sent .admPkgKey() down its tryCatch so every
# daemon derived a DIFFERENT sens cache key from its parent -- a silent divergence
# between the sequential and parallel results.
#
# This has to run in a real daemon. In-process the namespace is the dev one and is
# unlocked, so assignInNamespace() adds the binding happily and the test would
# pass whether or not the mechanism works.
skip_on_cran()
skip_if_not_installed("mirai")
skip_if_not_installed("rxode2")

test_that("a daemon resolves a function the installed package does not have", {
  # Two names that cannot exist in any released admixr2: one standing in for a
  # restart worker, one for a helper it calls. The helper is the point -- patching
  # only the worker was never the problem.
  fn_list <- list(
    .admDaemonPatchProbeWorker = function(restart_id, p_init, cores)
      .admDaemonPatchProbeHelper(p_init) + cores,
    .admDaemonPatchProbeHelper = function(x) x * 2)

  expect_false(exists(".admDaemonPatchProbeHelper",
                      envir = asNamespace("admixr2"), inherits = FALSE))

  .compute <- "admixr2DaemonPatchProbe"
  mirai::daemons(1L, .compute = .compute)
  on.exit(try(mirai::daemons(0L, .compute = .compute), silent = TRUE), add = TRUE)

  res <- mirai::mirai_map(
    1L, admixr2:::.admDaemonRestart,
    .args = list(worker_fn_name    = ".admDaemonPatchProbeWorker",
                 fn_list           = fn_list,
                 inits             = list(20),
                 all_args          = list(),
                 cores_vec         = 2L,
                 effective_workers = 1L),
    .compute = .compute)[]

  out <- res[[1L]]
  # A miraiError here IS the regression: "could not find function" for either name.
  expect_false(inherits(out, "miraiError") || inherits(out, "errorValue"),
               info = paste(as.character(out), collapse = " "))
  expect_identical(out, 42)          # 20 * 2 + 2
})

test_that("the dev payload does not ship the parent's model caches", {
  # The name filter that collects dev functions (^\\.(adm|adfo|adirmc|adgh|...))
  # also catches the memo ENVIRONMENTS, and an environment serialises by VALUE --
  # so every dev-mode dispatch was copying the parent's compiled rxode2 models to
  # every daemon. Heavy, and pointless: a deserialised model carries a dead pointer
  # and .admRxLoadAll() rejects it on the DLL check, so the worker rebuilds anyway.
  ns <- asNamespace("admixr2")
  nms <- ls(ns, all.names = TRUE)
  nms <- nms[grepl("^\\.(adm|adfo|adirmc|adgh|softmax|logdmvnorm)", nms)]
  envs <- nms[vapply(nms, function(x) is.environment(get(x, envir = ns)), logical(1))]
  # The filter has something to catch -- otherwise this test proves nothing.
  expect_true(length(envs) > 0L)
  expect_true(".adm_model_env" %in% envs)

  payload <- setNames(lapply(nms, get, envir = ns), nms)
  payload <- payload[!vapply(payload, is.environment, logical(1))]
  expect_false(any(envs %in% names(payload)))
  # ... and the functions are all still there.
  expect_true(".admGrad" %in% names(payload))
})
