# Inter-occasion variability (IOV) for adgh and admc.
#
# nlmixr2est owns the hard part. With `attr(nlmixr2Est.<est>, "iov") <- TRUE` it
# rewrites `iov.cl ~ 0.1 | occ` BEFORE the estimator runs -- one estimated
# magnitude theta (`iov.cl`) plus one fixed unit-variance eta per occasion level,
# `rx.iov.cl <- sqrt(iov.cl^2) * (rx.iov.cl.1 * (occ == 1) + ...)` -- and restores
# the user's model on the fit afterwards. The occasion levels are read from the
# `occ` column of the data handed to nlmixr2(), hence admData(occasions =).
#
# What is left for admixr2: an aggregate study has no per-record `occ` column, so
# `occ` is derived from time and from the study's own occasion start times, which
# arrive as the covariates `iov_t2`, `iov_t3`, ... (see admStudy(occasions =)).

# The occasion variable and its levels, read off nlmixr2est's rewritten model
# (`<sym> == <level>` comparisons on the right of a `rx.*` line). NULL when the
# model has no occasion expansion.
.admIovOcc <- function(ui) {
  vars <- list()
  walk <- function(e) {
    if (is.call(e)) {
      if (identical(e[[1L]], as.name("==")) && length(e) == 3L &&
          is.name(e[[2L]]) && is.numeric(e[[3L]])) {
        v <- as.character(e[[2L]])
        vars[[v]] <<- c(vars[[v]], e[[3L]])
      }
      for (a in as.list(e)[-1L]) walk(a)
    }
  }
  for (ex in ui$lstExpr) {
    lhs <- tryCatch(as.character(ex[[2L]]), error = function(e) "")
    if (length(lhs) == 1L && startsWith(lhs, "rx.")) walk(ex[[3L]])
  }
  vars <- lapply(vars[names(vars) %in% ui$allCovs], function(l) sort(unique(l)))
  if (!length(vars)) return(NULL)
  if (length(vars) > 1L)
    stop("admixr2: only one occasion variable is supported (the model has ",
         paste(sQuote(names(vars)), collapse = ", "), ").", call. = FALSE)
  list(var = names(vars), levels = vars[[1L]])
}

# The model with `occ` DEFINED from time instead of read from a data column:
#   occ <- L1*(time < iov_t2) + L2*(time >= iov_t2 & time < iov_t3) + ...
.admIovOccUi <- function(ui, occ) {
  L <- occ$levels; K <- length(L)
  if (K < 2L)
    stop("admixr2: the occasion variable ", sQuote(occ$var), " has one level in ",
         "the data given to nlmixr2(); IOV needs at least two. Give the levels ",
         "with admData(occasions = c(", occ$var, " = K)).", call. = FALSE)
  nm <- paste0("iov_t", seq_len(K))
  part <- vapply(seq_len(K), function(k) {
    lo <- if (k > 1L) paste0("time >= ", nm[k])
    hi <- if (k < K) paste0("time < ", nm[k + 1L])
    paste0(format(L[k], digits = 15), " * (", paste(c(lo, hi), collapse = " & "), ")")
  }, "")
  line <- str2lang(paste0(occ$var, " <- ", paste(part, collapse = " + ")))
  # piping returns a compressed rxUi; the drivers assign into the model's environment
  rxode2::rxUiDecompress(suppressWarnings(eval(bquote(rxode2::model(.(ui), .(line), append = 0)))))
}

# Dose times of an et() event table (addl/ii expanded); NULL for anything else.
.admEtDoseTimes <- function(ev) {
  if (!inherits(ev, "rxEt")) return(NULL)
  d <- as.data.frame(ev)
  if (!"evid" %in% names(d)) return(NULL)
  d <- d[d$evid %in% c(1, 4, 101, 104), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  ii <- if ("ii" %in% names(d)) d$ii else 0
  ad <- if ("addl" %in% names(d)) d$addl else 0
  ii[is.na(ii)] <- 0; ad[is.na(ad)] <- 0
  sort(unique(unlist(lapply(seq_len(nrow(d)), function(r) d$time[r] + ii[r] * seq(0, ad[r])))))
}

# Where each model occasion level starts for ONE study. Level k is active on
# [t_k, t_(k+1)); `labels` is the occasion of each observation time, as an `occ`
# data column would give it. A level starts at the dose that opens it, which the
# labels identify only when exactly ONE dose falls between the previous occasion's
# last observation and this one's first (no dose: the first observation). Several
# candidate doses are ambiguous and refused: say where the occasion starts with
# `starts` (one time per observed occasion after the first). A level the study
# never observes gets an empty interval.
.admOccStarts <- function(times, labels, levels, dose_times = NULL, starts = NULL) {
  K <- length(levels); never <- 1e10
  t <- rep(never, K)
  if (is.null(labels)) { t[-1L] <- never; return(t) }
  o <- order(times); times <- times[o]; labels <- labels[o]
  if (is.unsorted(match(labels, levels), na.rm = FALSE) || anyNA(match(labels, levels)))
    stop("admixr2: `occasions` must name the model's occasion levels (",
         paste(levels, collapse = ", "), ") and not go back to an earlier occasion.",
         call. = FALSE)
  present <- which(levels %in% labels)
  first   <- present[1L]
  if (!is.null(starts) && length(starts) != length(present) - 1L)
    stop("admixr2: `occasion_starts` needs one time per observed occasion after the first (",
         length(present) - 1L, ").", call. = FALSE)
  for (k in rev(seq_len(K))) {
    nxt <- if (k < K) t[k + 1L] else never
    if (!k %in% present) { t[k] <- if (k < first) -never else nxt; next }
    ok <- labels == levels[k]
    s  <- min(times[ok])
    if (k > first) {
      prev_end <- max(times[match(labels, levels) < k])
      j <- match(k, present[-1L])
      if (!is.null(starts)) {
        if (!(starts[j] > prev_end && starts[j] <= s))
          stop("admixr2: `occasion_starts` for occasion ", levels[k], " (", starts[j],
               ") must fall after the last time of the earlier occasion (", prev_end,
               ") and not after its own first time (", s, ").", call. = FALSE)
        s <- starts[j]
      } else {
        d <- dose_times[dose_times <= s & dose_times > prev_end]
        if (length(d) > 1L)
          stop("admixr2: occasion ", levels[k], " could have opened at any of the doses at ",
               paste(d, collapse = ", "), " (between times ", prev_end, " and ", s,
               "). Say which with `occasion_starts = ", format(d[length(d)]), "`.",
               call. = FALSE)
        if (length(d)) s <- d
      }
    } else s <- -never
    t[k] <- s
  }
  t
}

# Which studies can tell between-occasion variance from between-subject variance: those
# that observe two or more occasions AND report their covariance (a model source
# generates it; a digitised study needs a full matrix V, not just SDs). With none, the
# two variances share one flat ridge -- only their SUM is identified.
.admIovInformative <- function(studies) {
  vapply(studies, function(s) {
    lab <- s[["occasions"]]
    if (length(unique(lab)) < 2L) return(FALSE)
    has_model <- !is.null(s[["ui"]]) || !is.null(s[["model"]])
    has_model || is.matrix(s[["V"]])
  }, logical(1))
}

# Give every study the covariates the rewritten model reads: iov_t2 ... iov_tK.
.admIovPad <- function(studies, occ) {
  if (!any(.admIovInformative(studies)))
    warning("admixr2: no study observes two or more occasions with a covariance (a ",
            "published model, or a full matrix `V`), so the between-occasion variance ",
            "cannot be separated from the between-subject variance: only their sum is ",
            "identified, and the split stays wherever the optimizer starts. Give ",
            "`occasions` to a study that reports its covariance.", call. = FALSE)
  K <- length(occ$levels)
  nm <- paste0("iov_t", seq_len(K))[-1L]
  one <- function(s) {
    raw <- !inherits(s, "admStudy")
    lab <- s[["occasions"]]
    dt  <- .admEtDoseTimes(if (raw) s[["ev"]] else s[["dose"]])
    t   <- .admOccStarts(s[["times"]], lab, occ$levels, dt, s[["occasion_starts"]])[-1L]
    if (!raw) s[[".iov_levels"]] <- occ$levels
    add <- as.list(stats::setNames(t, nm))
    if (raw) s[["cov"]] <- utils::modifyList(add, as.list(s[["cov"]]))
    else     s$at <- utils::modifyList(add, as.list(s$at))
    s
  }
  studies[] <- lapply(studies, one)
  studies
}

# A published model WRITTEN with an IOV term (`iov.cl ~ 0.05 | occ`) as a source: the
# variance is a published number, so there is nothing for nlmixr2est's estimation
# rewrite to do. Each occasion `level` gets its own eta carrying that variance, and
# the term is replaced in the model by the one active at `time` (the study's own
# occasion start times arrive as the covariates iov_t2, iov_t3, ...). NULL when the
# model has no IOV term.
.admIovSourceUi <- function(ui, levels) {
  d <- ui$iniDf
  cond <- as.character(d$condition)
  r <- d[!is.na(d$neta1) & d$neta1 == d$neta2 & !is.na(cond) & nzchar(cond) & cond != "id", ,
         drop = FALSE]
  if (!nrow(r)) return(NULL)
  if (length(unique(cond[match(r$name, d$name)])) > 1L)
    stop("admixr2: only one occasion variable is supported.", call. = FALSE)
  K <- length(levels); nm <- paste0("iov_s", seq_len(K))
  swap <- function(e, v, with) {
    if (is.name(e) && identical(as.character(e), v)) return(with)
    if (is.call(e)) for (i in seq_along(e)[-1L]) e[[i]] <- swap(e[[i]], v, with)
    e
  }
  etas <- list()
  u <- ui
  for (v in r$name) {
    en <- paste0(v, ".occ", seq_len(K))
    part <- vapply(seq_len(K), function(k) {
      lo <- if (k > 1L) paste0("time >= ", nm[k])
      hi <- if (k < K) paste0("time < ", nm[k + 1L])
      cnd <- paste(c(lo, hi), collapse = " & ")
      if (K == 1L) en[k] else paste0(en[k], " * (", cnd, ")")
    }, "")
    with <- str2lang(paste0("(", paste(part, collapse = " + "), ")"))
    for (ex in u$lstExpr) {
      if (!any(all.vars(ex) == v)) next
      line <- swap(ex, v, with)
      u <- rxode2::rxUiDecompress(suppressWarnings(suppressMessages(
        eval(bquote(rxode2::model(.(u), .(line)))))))
    }
    etas[[v]] <- list(names = en, est = r$est[r$name == v])
  }
  for (v in names(etas)) for (en in etas[[v]]$names)
    u <- rxode2::rxUiDecompress(suppressWarnings(suppressMessages(
      eval(bquote(rxode2::ini(.(u), .(as.name(en)) ~ .(etas[[v]]$est)))))))
  u
}

# TRUE when the model declares an inter-occasion term (a diagonal eta on a level other than `id`).
.admHasIov <- function(ui) {
  d <- ui$iniDf; cond <- as.character(d$condition)
  any(!is.na(d$neta1) & d$neta1 == d$neta2 & !is.na(cond) & nzchar(cond) & cond != "id")
}

# `iov_t2..` for ONE concrete study (its final `ev`), overriding any earlier value: the
# dose function / dose bands of an admStudy() only become event tables when expanded.
.admIovCov <- function(cov, times, occasions, starts, levels, ev, prefix = "iov_t") {
  t <- .admOccStarts(times, occasions, levels, .admEtDoseTimes(ev), starts)[-1L]
  utils::modifyList(as.list(cov %||% list()),
                    as.list(stats::setNames(t, paste0(prefix, seq_along(levels))[-1L])))
}
