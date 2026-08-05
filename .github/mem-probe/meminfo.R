# Shared instrument for both probes.
#
# The point of this file: "native memory grew" is a symptom, not a cause. The
# split that discriminates between the candidate causes is ANONYMOUS versus
# FILE-BACKED resident memory, and /proc/self/smaps_rollup gives it directly:
#
#   Rss - Anonymous  = file-backed resident  -> dlopen'd .so text/data, i.e.
#                      rxode2 compiling a model per fit and never dlclose()ing
#   Anonymous        -> malloc'd and never freed, or freed but retained by
#                      glibc in per-thread arenas (the MALLOC_ARENA_MAX theory)
#
# If devel and release differ, ONE of those two columns is where it shows, and
# that single fact picks the root cause. /proc/self/maps additionally counts how
# many distinct shared objects are mapped, which tests the dlopen story head on:
# a leak of mappings must show a rising .so count.

.meminfo <- function() {
  .kv <- function(lines, key) {
    v <- grep(paste0("^", key, ":"), lines, value = TRUE)[1L]
    if (is.na(v)) return(NA_real_)
    as.numeric(sub("[^0-9]*([0-9]+).*", "\\1", v)) / 1024      # kB -> MB
  }

  st <- readLines("/proc/self/status", warn = FALSE)
  out <- c(rss = .kv(st, "VmRSS"), hwm = .kv(st, "VmHWM"),
           vsz = .kv(st, "VmSize"), thr = NA_real_)
  tv <- grep("^Threads:", st, value = TRUE)[1L]
  if (!is.na(tv)) out[["thr"]] <- as.numeric(sub("[^0-9]*([0-9]+).*", "\\1", tv))

  ro <- tryCatch(readLines("/proc/self/smaps_rollup", warn = FALSE),
                 error = function(e) character(0))
  anon <- if (length(ro)) .kv(ro, "Anonymous") else NA_real_
  pdirty <- if (length(ro)) .kv(ro, "Private_Dirty") else NA_real_
  sclean <- if (length(ro)) .kv(ro, "Shared_Clean") else NA_real_

  mp <- tryCatch(readLines("/proc/self/maps", warn = FALSE),
                 error = function(e) character(0))
  paths <- sub("^\\S+\\s+\\S+\\s+\\S+\\s+\\S+\\s+\\S+\\s*", "", mp)
  so <- unique(paths[grepl("\\.so($|\\.)", paths)])

  gcv <- gc(verbose = FALSE)
  c(out,
    heap = sum(gcv[, 2L]),
    anon = anon,
    filebacked = if (is.na(anon)) NA_real_ else out[["rss"]] - anon,
    pdirty = pdirty, sclean = sclean,
    n_maps = length(mp), n_so = length(so))
}

.MEM_COLS <- c("rss", "hwm", "vsz", "thr", "heap", "anon", "filebacked",
               "pdirty", "sclean", "n_maps", "n_so")

.memhdr <- function(prefix, extra = character(0)) {
  cat(paste0(prefix, "\t", paste(c(extra, .MEM_COLS), collapse = "\t"), "\n"))
  flush(stdout())
}

.memrow <- function(prefix, extra, m) {
  vals <- vapply(.MEM_COLS, function(k) {
    v <- unname(m[[k]])
    if (is.na(v)) "NA" else if (k %in% c("n_maps", "n_so", "thr"))
      sprintf("%d", as.integer(v)) else sprintf("%.1f", v)
  }, character(1))
  cat(paste0(prefix, "\t", paste(c(extra, vals), collapse = "\t"), "\n"))
  flush(stdout())
}

.memenv <- function() {
  cat(sprintf("# R              : %s\n", R.version.string))
  cat(sprintf("# platform       : %s\n", R.version$platform))
  cat(sprintf("# MALLOC_ARENA_MAX: %s\n", Sys.getenv("MALLOC_ARENA_MAX", "<unset>")))
  cat(sprintf("# glibc          : %s\n",
              tryCatch(system("ldd --version 2>/dev/null | head -1", intern = TRUE)[1L],
                       error = function(e) "?")))
  for (p in c("rxode2", "nlmixr2est", "RcppEigen", "RcppParallel")) {
    v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA)
    if (!is.na(v)) {
      # Source-built .so keep their symbol/debug sections; P3M binaries are
      # stripped. A size difference here is itself a candidate explanation.
      so <- file.path(find.package(p), "libs", paste0(p, ".so"))
      sz <- if (file.exists(so)) sprintf("%.1f MB", file.size(so) / 1024^2) else "-"
      cat(sprintf("# %-14s : %-10s .so=%s\n", p, v, sz))
    }
  }
  flush(stdout())
}
