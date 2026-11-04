# Run the Tier 2 integration tests under coverage.
#
# covr::package_coverage() executes the test suite in a *separate* R session
# (via callr) whose environment is built from rcmd_safe_env() and therefore
# does NOT inherit the NOT_CRAN variable set at the CI job level. The
# integration tests guard on skip_on_cran(), so under covr they all skip and
# their code (the estimator pipelines, gradients, .adirmcPhaseLoop, the
# plot.admFit panels, ...) is reported as uncovered. That is why coverage sat
# at ~77% even after NOT_CRAN was added to test-coverage.yaml in #75 — the env
# var was set on the covr driver process, not the session that runs the tests.
#
# covr DOES export R_COVR=true into that session, so detect it here and turn
# NOT_CRAN on from inside the process where the tests actually run. This only
# affects coverage runs:
#   * local devtools::test() already sets NOT_CRAN -> integration tests run
#   * R CMD check runs them via rcmdcheck's default NOT_CRAN=true
#   * plain testthat::test_dir() (R_COVR unset) is unaffected -> Tier 2 skips
if (identical(Sys.getenv("R_COVR"), "true") &&
    !identical(Sys.getenv("NOT_CRAN"), "true")) {
  Sys.setenv(NOT_CRAN = "true")
}

# Stop admixr2's own worker processes before the session ends, under coverage.
#
# covr instruments the INSTALLED copy and merges every trace any process
# loading it wrote. admixr2's restart workers are mirai daemons -- separate R
# processes that each dump a trace as they exit -- so one still writing while
# covr walks the directory leaves a half-written file, and the merge dies with
#   Error in readRDS(f) : error reading from connection
# which is how this job has been failing. The workflow already opted out of
# parallel testthat for the same reason; these are the OTHER processes, which
# that setting does not reach. Shut them down deterministically and let the
# writes land.
if (identical(Sys.getenv("R_COVR"), "true")) {
  withr::defer({
    try(suppressMessages(admixr2::admStopWorkers()), silent = TRUE)
    Sys.sleep(3)
  }, testthat::teardown_env())
}
