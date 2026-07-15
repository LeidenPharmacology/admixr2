# Temporary CI diagnostic, split so we can SEE where things go:
#   stage=heavy    -> run every integration file EXCEPT parallel, per-file RSS+peak
#                     (fast; answers "does the per-fit teardown fix cap memory?")
#   stage=parallel -> run ONLY test-integration-parallel.R with a per-test_that
#                     reporter logging RSS + wall-time per block (answers "does the
#                     parallel/mirai test hang, and is it memory or daemons?")
# ADM_TEARDOWN_OFF=1 disables the fix for a baseline.
stage <- commandArgs(trailingOnly = TRUE)[1]
Sys.setenv(NOT_CRAN = "true")
suppressMessages({library(testthat); library(nlmixr2); library(admixr2); library(rxode2); library(R6)})
options(admixr2.teardown_unload = !identical(Sys.getenv("ADM_TEARDOWN_OFF"), "1"))
rss <- function() { v <- tryCatch(as.numeric(sub("[^0-9]*([0-9]+).*","\\1",
  grep("VmRSS", readLines("/proc/self/status"), value=TRUE))), error=function(e) NA); round(v/1024) }
st  <- function(...) { cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep=""); flush(stdout()) }
st("teardown_unload = ", getOption("admixr2.teardown_unload"), "  stage=", stage)
setwd("tests/testthat")

if (identical(stage, "heavy")) {
  files <- setdiff(sort(list.files(".", pattern = "^test-integration.*\\.R$")),
                   "test-integration-parallel.R")
  peak <- rss(); totfail <- 0
  for (f in files) {
    r <- tryCatch(as.data.frame(testthat::test_file(f, reporter = "silent")), error = function(e) e)
    fail <- if (is.data.frame(r)) sum(r$failed) + sum(r$error) else NA
    if (!is.na(fail)) totfail <- totfail + fail
    peak <- max(peak, rss())
    st(sprintf("ran %-42s RSS=%6d peak=%6d fails=%s", f, rss(), peak, fail))
  }
  st("HEAVY COMPLETE  peak RSS=", peak, " MiB  total fails=", totfail)

} else if (identical(stage, "parallel")) {
  rep <- R6::R6Class("R", inherit = testthat::Reporter, public = list(
    t0 = NULL,
    start_test = function(context, test) { self$t0 <- Sys.time() },
    end_test = function(context, test)
      st(sprintf("  PAR block done in %4.1fs  RSS=%6d   %s",
         as.numeric(difftime(Sys.time(), self$t0, units="secs")), rss(), substr(test,1,50)))))$new()
  st("start RSS=", rss(), " -- running test-integration-parallel.R")
  r <- tryCatch(as.data.frame(testthat::test_file("test-integration-parallel.R", reporter = rep)),
                error = function(e) { st("ERR: ", conditionMessage(e)); NULL })
  if (!is.null(r)) st("PARALLEL COMPLETE fails=", sum(r$failed)+sum(r$error), " RSS=", rss())
}
st("DONE stage=", stage)
