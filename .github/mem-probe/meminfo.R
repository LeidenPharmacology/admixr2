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

  # SYSTEM-WIDE, because the parent process is not the whole story: the suite
  # spawns mirai daemons (test-integration-parallel), and a daemon's RSS is
  # invisible to /proc/self. What kills a runner is total machine memory, so
  # measure that too -- rss_all is every R process on the box, and sys_used is
  # MemTotal - MemAvailable, i.e. what the OOM killer actually sees.
  mi <- tryCatch(readLines("/proc/meminfo", warn = FALSE),
                 error = function(e) character(0))
  memtot <- if (length(mi)) .kv(mi, "MemTotal") else NA_real_
  memav  <- if (length(mi)) .kv(mi, "MemAvailable") else NA_real_

  rss_all <- NA_real_; n_r <- NA_real_
  pids <- suppressWarnings(as.integer(list.files("/proc")))
  pids <- pids[!is.na(pids)]
  if (length(pids)) {
    tot <- 0; k <- 0L
    for (p in pids) {
      st2 <- tryCatch(readLines(file.path("/proc", p, "status"), warn = FALSE),
                      error = function(e) character(0))
      if (!length(st2)) next
      nm <- sub("^Name:\\s*", "", grep("^Name:", st2, value = TRUE)[1L])
      if (is.na(nm) || !grepl("^(R|exec/R|Rscript)$", nm)) next
      v <- .kv(st2, "VmRSS")
      if (!is.na(v)) { tot <- tot + v; k <- k + 1L }
    }
    rss_all <- tot; n_r <- k
  }

  # Ncells vs Vcells splits the heap by KIND, which object.size() cannot do:
  # Ncells are cons cells -- language objects, closures, environments, srcrefs;
  # Vcells are vector payloads -- numeric data. Which one grows says whether
  # the retained material is CODE/structure or DATA.
  gcv <- gc(verbose = FALSE)
  c(out,
    heap = sum(gcv[, 2L]),
    ncells = gcv[1L, 2L], vcells = gcv[2L, 2L],
    anon = anon,
    filebacked = if (is.na(anon)) NA_real_ else out[["rss"]] - anon,
    pdirty = pdirty, sclean = sclean,
    n_maps = length(mp), n_so = length(so),
    rss_all = rss_all, n_rproc = n_r,
    sys_used = if (is.na(memtot) || is.na(memav)) NA_real_ else memtot - memav,
    sys_tot = memtot)
}

.MEM_COLS <- c("rss", "hwm", "vsz", "thr", "heap", "ncells", "vcells",
               "anon", "filebacked", "pdirty", "sclean", "n_maps", "n_so",
               "rss_all", "n_rproc", "sys_used", "sys_tot")

.memhdr <- function(prefix, extra = character(0)) {
  cat(paste0(prefix, "\t", paste(c(extra, .MEM_COLS), collapse = "\t"), "\n"))
  flush(stdout())
}

.memrow <- function(prefix, extra, m) {
  vals <- vapply(.MEM_COLS, function(k) {
    v <- if (k %in% names(m)) unname(m[[k]]) else NA_real_
    if (is.na(v)) "NA" else if (k %in% c("n_maps", "n_so", "thr", "n_rproc"))
      sprintf("%d", as.integer(v)) else sprintf("%.1f", v)
  }, character(1))
  cat(paste0(prefix, "\t", paste(c(extra, vals), collapse = "\t"), "\n"))
  flush(stdout())
}

.memenv <- function() {
  cat(sprintf("# R              : %s\n", R.version.string))
  cat(sprintf("# platform       : %s\n", R.version$platform))
  cat(sprintf("# MALLOC_ARENA_MAX: %s\n", Sys.getenv("MALLOC_ARENA_MAX", "<unset>")))
  # srcref state -- the leading hypothesis for the devel/release heap gap.
  # R_KEEP_PKG_SOURCE governs whether an INSTALLED package keeps a srcref per
  # function; keep.source governs code parsed at RUNTIME. R-devel installs its
  # whole stack from source on the runner, so the first one actually bites there
  # while release gets pre-built binaries that carry none.
  cat(sprintf("# R_KEEP_PKG_SOURCE: %s | keep.source: %s | keep.source.pkgs: %s\n",
              Sys.getenv("R_KEEP_PKG_SOURCE", "<unset>"),
              isTRUE(getOption("keep.source")),
              isTRUE(getOption("keep.source.pkgs"))))
  # Direct evidence: how many functions in the dependency carry a srcref.
  for (p in c("rxode2", "nlmixr2est", "admixr2")) {
    ns <- tryCatch(asNamespace(p), error = function(e) NULL)
    if (is.null(ns)) next
    nms <- tryCatch(ls(ns, all.names = TRUE), error = function(e) character(0))
    n <- 0L; k <- 0L
    for (x in nms) {
      o <- tryCatch(get(x, envir = ns, inherits = FALSE), error = function(e) NULL)
      if (!is.function(o)) next
      k <- k + 1L
      if (!is.null(attr(o, "srcref"))) n <- n + 1L
    }
    cat(sprintf("#   %-11s functions=%4d with srcref=%4d\n", p, k, n))
  }
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
