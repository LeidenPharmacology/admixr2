# Probe A -- the MECHANISM, isolated from admixr2 entirely.
#
# Compiling and loading N distinct rxode2 models is claimed to grow native
# memory that gc() cannot reach. This asks three things at once:
#
#   1. does it, and at what MB per model?
#   2. is that growth ANONYMOUS (malloc/arena) or FILE-BACKED (dlopen'd .so)?
#      -- the split that picks the root cause; see meminfo.R
#   3. does R-devel grow it faster than R-release, same machine, same rxode2?
#
# Run twice per job, with and without MALLOC_ARENA_MAX=2. That pair is the test
# of the glibc-arena hypothesis: if the memory is freed-but-retained, capping
# arenas at 2 must flatten the curve; if it is genuinely still allocated, or is
# file-backed, capping arenas changes nothing.
#
# One flushed row per model, so an OOM kill still leaves the curve up to it.
# Run from the repository root.

source(".github/mem-probe/meminfo.R")

args <- commandArgs(trailingOnly = TRUE)
tag  <- if (length(args)) args[[1L]] else "default"
N    <- as.integer(Sys.getenv("PROBE_N", "40"))

.memenv()
suppressMessages(library(rxode2))
.memhdr("PROBE_A", c("tag", "step", "dll_kb"))
.memrow("PROBE_A", c(tag, "baseline", "0"), .meminfo())

ev <- rxode2::eventTable()
ev$add.dosing(dose = 100)
ev$add.sampling(seq(0, 24, by = 4))

dll_tot <- 0
for (i in seq_len(N)) {
  # A distinct constant per model -> a distinct DLL, i.e. no cache hits.
  txt <- sprintf(paste0("ka = 1.0;\ncl = %0.6f;\nv = 20;\n",
                        "d/dt(depot) = -ka*depot;\n",
                        "d/dt(central) = ka*depot - (cl/v)*central;\n",
                        "cp = central/v;\n"), 1 + i / 1000)
  m <- rxode2::rxode2(txt)
  invisible(rxode2::rxSolve(m, ev))
  # Size of the .so rxode2 just built and dlopen'd. If the growth is
  # file-backed, this is the per-model quantum it should track.
  kb <- tryCatch(file.size(rxode2::rxDll(m)) / 1024, error = function(e) NA_real_)
  if (!is.na(kb)) dll_tot <- dll_tot + kb
  rm(m)
  .memrow("PROBE_A", c(tag, sprintf("model_%03d", i),
                       if (is.na(kb)) "NA" else sprintf("%.0f", kb)), .meminfo())
}
cat(sprintf("# total model .so built: %.1f MB\n", dll_tot / 1024))
flush(stdout())

# Does anything come back? Three escalating attempts, each measured. n_so is the
# one to read: if it does not fall, the mappings were never unmapped.
gc(full = TRUE)
.memrow("PROBE_A", c(tag, "after_gc_full", "0"), .meminfo())
tryCatch({ rxode2::rxUnloadAll(); gc(full = TRUE) }, error = function(e) NULL)
.memrow("PROBE_A", c(tag, "after_rxUnloadAll", "0"), .meminfo())
tryCatch({ rxode2::rxClean(); gc(full = TRUE) }, error = function(e) NULL)
.memrow("PROBE_A", c(tag, "after_rxClean", "0"), .meminfo())
