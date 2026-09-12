# Writing down a published study
# Reserved data columns cannot be covariates; every other column is one.
.ADM_DATA_COLS <- c("ID", "TIME", "DV", "AMT", "EVID", "CMT", "DVID", "MDV",
                    "RATE", "DUR", "SS", "II", "ADDL", "CENS", "LIMIT")

.admPopFromData <- function(data, dist, given) {
  if (!is.data.frame(data)) data <- as.data.frame(data)
  given_nms <- names(given)
  nms <- setdiff(names(data), given_nms) # anything stated in `...` wins outright
  bad <- function(nm, ...) stop("admixr2: covariate '", nm, "' ", ...,
                                call. = FALSE)
  binary <- function(p) list(values = c(0, 1), probs = c(1 - p, p))

  # Overridden margins still contribute cohort correlations.
  cont <- Filter(function(nm) {
    v <- data[[nm]]
    is.numeric(v) && all(is.finite(v)) &&
      !all(unique(v) %in% c(0, 1)) &&
      stats::sd(v) >= .Machine$double.eps^0.5 * max(1, abs(mean(v))) &&
      (!identical(dist, "lnorm") || all(v > 0)) &&
      is.null(tryCatch(.admPopSpec(given[[nm]], nm, dist)[["values"]],
                       error = function(e) TRUE))
  }, intersect(given_nms, names(data)))
  specs <- list()
  for (nm in nms) {
    v <- data[[nm]]
    if (is.factor(v) || is.character(v)) {
      lv <- sort(unique(as.character(v)))
      if (length(lv) != 2L)
        bad(nm, "has ", length(lv), " levels. A data frame supplies BINARY ",
            "covariates; give a multi-level one yourself as `", nm,
            " = list(values = ..., probs = ...)`.")
      # The second sorted level maps to 1.
      specs[[nm]] <- binary(mean(as.character(v) == lv[2L]))
      next
    }
    v <- as.numeric(v)
    if (!all(is.finite(v)))
      bad(nm, "has missing or non-finite values. A margin has to describe ",
          "every subject the study reports on -- drop the column, or state ",
          "it yourself.")
    if (all(unique(v) %in% c(0, 1))) { specs[[nm]] <- binary(mean(v)); next }
    # Constant columns are pinned values, not distributions.
    if (stats::sd(v) < .Machine$double.eps^0.5 * max(1, abs(mean(v))))
      bad(nm, "is CONSTANT at ", signif(v[1L], 6), " in this data frame, so ",
          "it has no distribution to describe. A covariate every subject ",
          "shares is a value the study was reported AT, not a population it ",
          "was reported OVER -- pin it with `at = list(", nm, " = ",
          signif(v[1L], 6), ")` and leave it out of `population`.")
    if (identical(dist, "lnorm") && any(v <= 0))
      bad(nm, "has values at or below zero, which a lognormal margin cannot ",
          "describe. Use `dist = \"normal\"`, or state that margin yourself.")
    specs[[nm]] <- c(mean = mean(v), sd = stats::sd(v))
    cont <- c(cont, nm)
  }

  # Copula correlations use each resolved margin's latent scale.
  resolved <- c(given, specs)
  resolved <- stats::setNames(
    lapply(names(resolved), function(nm) .admPopSpec(resolved[[nm]], nm, dist)),
    names(resolved))
  tf <- function(nm) {
    v <- data[[nm]]
    if (!is.null(resolved[[nm]][["meanlog"]])) log(v) else v
  }
  rho <- numeric(0)
  if (length(cont) > 1L)
    for (p in asplit(utils::combn(cont, 2L), 2L)) {
      r <- stats::cor(tf(p[1L]), tf(p[2L]))
      if (is.finite(r) && abs(r) > 1e-8) rho[paste(p, collapse = ".")] <- r
    }

  # Report unsupported discrete-continuous dependence above a practical floor
  # and three sampling SEs.
  disc <- names(resolved)[vapply(resolved, function(s) !is.null(s[["values"]]),
                                 logical(1))]
  if (length(disc) && length(cont)) {
    lim  <- max(0.1, 3 / sqrt(nrow(data)))
    hit  <- function(d, c) {
      v <- data[[d]]
      if (!is.numeric(v)) v <- as.numeric(factor(v))
      r <- suppressWarnings(stats::cor(v, data[[c]]))
      is.finite(r) && abs(r) > lim
    }
    drop <- outer(disc, cont, Vectorize(hit))
    if (any(drop))
      message("admixr2: `data` shows ",
              paste(paste0(disc[row(drop)[drop]], "/", cont[col(drop)[drop]]),
                    collapse = ", "),
              " correlated, and that correlation is being DROPPED -- the ",
              "levels are enumerated exactly and taken as independent of the ",
              "rest. Band on the discrete one (`stratify`) if the association ",
              "carries information you need.")
  }
  list(specs = specs, cor = rho)
}

#' Describe the population a study enrolled
#'
#' Written the way a baseline demographics table reads. Each covariate takes
#' whichever summary the paper printed --- `mean`/`sd`, `median`/`iqr`,
#' a `cv` as a percent, or a proportion for a binary one --- and correlations
#' are given for the PAIRS that were reported, everything else being
#' independent.
#'
#' @param ... Named covariates. A continuous one takes a named vector, e.g.
#'   `WT = c(mean = 75, sd = 16)` or `CRCL = c(median = 92, iqr = c(62, 118))`.
#'   A binary one takes a single named proportion, e.g. `SEX = c(male = 0.55)`,
#'   which becomes levels `0`/`1` with that probability on `1`.
#' @param cor Correlations between covariate PAIRS, named `A.B`, e.g.
#'   `cor = c(WT.CRCL = 0.45)`. Pairs not named are independent, so a partial
#'   table needs no identity padding. A full matrix is accepted too.
#' @param dist `"lnorm"` (default) or `"normal"`, for the continuous margins.
#'   Lognormal is the usual choice for a positive covariate --- a normal margin
#'   wide enough to matter puts mass at or below zero, which is `NaN` inside any
#'   power or log term.
#'
#' @param data A data frame of individual covariates to derive the table
#'   FROM, instead of typing it out --- a digitised baseline listing, or the
#'   cohort itself in a simulation study. Each numeric column becomes a margin
#'   (a 0/1 column becomes a proportion, everything else a continuous margin
#'   matching that column's mean and SD), and every continuous PAIR gets its
#'   correlation --- taken on the LATENT scale, so on the logs for a lognormal
#'   margin, which is the conversion easiest to get wrong by hand. Anything
#'   named in `...` or `cor` overrides what the data would have given, so a
#'   column you would rather state yourself simply gets stated.
#' @return A covariate specification, as [covDist()] returns.
#' @seealso [admStudy()], which takes one; [covDraw()] to inspect it.
#' @export
admPopulation <- function(..., cor = NULL, dist = c("lnorm", "normal"),
                          data = NULL) {
  dist <- match.arg(dist)
  a <- list(...)
  if (!is.null(data)) {
    d <- .admPopFromData(data, dist, a)
    a <- c(a, d$specs)
    # Explicit correlations override cohort values pairwise.
    if (is.null(cor)) { if (length(d$cor)) cor <- d$cor }
    else if (!is.matrix(cor)) {
      # Canonicalise unordered pairs without splitting dotted covariate names.
      pair_key <- function(x) {
        if (length(a) < 2L) return(x)
        pairs <- utils::combn(names(a), 2L, simplify = FALSE)
        hit <- vapply(pairs, function(p)
          identical(x, paste(p, collapse = ".")) ||
          identical(x, paste(rev(p), collapse = ".")), logical(1))
        if (sum(hit) == 1L) paste(sort(pairs[[which(hit)]]), collapse = "\r") else x
      }
      user_keys <- vapply(names(cor), pair_key, character(1))
      data_keys <- vapply(names(d$cor), pair_key, character(1))
      cor <- c(cor, d$cor[!data_keys %in% user_keys])
    }
  }
  if (!length(a) || is.null(names(a)) || any(!nzchar(names(a))))
    stop("admixr2: `admPopulation()` needs NAMED covariates, e.g. ",
         "admPopulation(WT = c(mean = 75, sd = 16)), or a `data` frame to ",
         "derive them from.", call. = FALSE)
  nms <- names(a)
  specs <- stats::setNames(
    lapply(seq_along(a), function(i) .admPopSpec(a[[i]], nms[i], dist)), nms)
  lv <- vapply(specs, function(s) !is.null(s[["values"]]), logical(1))
  specs <- lapply(specs, function(s) { s$.level <- NULL; s })
  # Unspecified pairs remain independent.
  R <- NULL
  if (!is.null(cor)) {
    if (is.matrix(cor)) {
      R <- cor
      if (is.null(rownames(R))) dimnames(R) <- list(nms, nms)
    } else {
      R <- diag(length(nms)); dimnames(R) <- list(nms, nms)
      cn <- names(cor)
      if (is.null(cn) || any(!nzchar(cn)))
        stop("admixr2: `cor` must name the PAIR it applies to, e.g. ",
             "cor = c(WT.CRCL = 0.45).", call. = FALSE)
      for (k in seq_along(cor)) {
        pr <- strsplit(cn[k], ".", fixed = TRUE)[[1L]]
        # Match declared names because covariate names may contain dots.
        if (length(pr) != 2L || !all(pr %in% nms)) {
          hit <- nms[vapply(nms, function(x) startsWith(cn[k], paste0(x, ".")),
                            logical(1))]
          pr <- if (length(hit) == 1L)
            c(hit, substring(cn[k], nchar(hit) + 2L)) else pr
        }
        if (length(pr) != 2L || !all(pr %in% nms))
          stop("admixr2: `cor` entry '", cn[k], "' does not name two declared ",
               "covariates. Declared: ", paste(sQuote(nms), collapse = ", "),
               "; write the pair as 'A.B'.", call. = FALSE)
        if (any(lv[pr]))
          stop("admixr2: `cor` entry '", cn[k], "' correlates ",
               paste(sQuote(pr[lv[pr]]), collapse = " and "),
               ", which is a DISCRETE covariate. A level would then be a ",
               "truncation of the latent normal rather than a point, which ",
               "admixr2 refuses -- see covStrata(). Declare it independent.",
               call. = FALSE)
        R[pr[1L], pr[2L]] <- R[pr[2L], pr[1L]] <- cor[[k]]
      }
    }
  }
  do.call(covDist, c(specs, list(cor = R, dist = dist)))
}

#' Write down a published study
#'
#' One study, transcribed from one paper. It holds what the paper reported and
#' nothing else, is checked as you build it, and prints so the transcription can
#' be read back against the source before anything is fitted.
#'
#' A study contributes in one of two currencies, and this takes either:
#'
#' * **a published model** --- `model`, with the paper's parameter table as
#'   `est`;
#' * **digitised aggregate data** --- `E` with `sd` (or `sem`, or `V`).
#'
#' A study generated from a model is not a sample: its mean and covariance are
#' exact functions of that model's parameters, and `n` sets its RELATIVE WEIGHT
#' against the other studies rather than its precision. No standard error is
#' reported for a fit that includes one --- see [admStudies()].
#'
#' Nothing is solved here. The study is generated when it reaches the fit, so
#' building one is cheap and a mistake surfaces on `print()` rather than after a
#' long run.
#'
#' @param model The published model, as an nlmixr2-style function (or a parsed
#'   `rxUi`). Omit for a digitised study.
#' @param est Named vector of the paper's parameter estimates, on the scale the
#'   model's `ini()` is written on --- so `tcl = log(5.2)` for a clearance
#'   reported as 5.2. Omitted parameters keep the model's own `ini()` value.
#' @param v_denom Which denominator the supplied spread uses: `"unbiased"`
#'   (`n - 1`) or `"ml"` (`n`). Usually leave it unset --- the currency you
#'   wrote the study in already says which it is, and `admStudy()` resolves it
#'   and shows the answer when the study is printed.
#'
#'   `sd`/`sem` default to `"unbiased"`, because a **published** spread is the
#'   `n - 1` one; `V` defaults to `"ml"`, because handing over a covariance
#'   matrix is the deliberate act of someone who computed it, and the
#'   likelihood wants the ML denominator. Set it explicitly when your source
#'   breaks that pattern --- a `sd` you computed yourself with the ML
#'   denominator, say. A published `model` is always `"ml"`: its moments are
#'   generated at that denominator, so `"unbiased"` is refused rather than
#'   applied to a `V` that is already right.
#' @param E,V,sd,sem Digitised aggregate data: the reported mean profile `E`,
#'   with its spread as `sd` (per timepoint), `sem` (converted using `n`), or a
#'   full covariance `V`.
#' @param n Number of subjects the study reports on. Taken from `population`
#'   when that is the cohort's data frame.
#' @param times Observation times.
#' @param dose Dose amount, as shorthand for a single-dose `ev`.
#' @param ev A dosing event table from [rxode2::et()], for anything `dose`
#'   cannot express.
#' @param population The population the study enrolled --- see
#'   [admPopulation()]. Covariates the model reads are integrated over it.
#'   A **data frame** of enrolled subjects is accepted directly and read into
#'   a table, using the columns this model actually reads and taking `n` from
#'   its row count.
#' @param at Named list pinning a covariate at a value, for a study reported in
#'   one subgroup, e.g. `at = list(SEX = 1)`. A pinned covariate must be omitted
#'   from `population`.
#' @param by Covariate the paper reports results SEPARATELY by, e.g.
#'   `by = "SEX"`. Expands a published model into one study per level; digitised
#'   subgroup profiles must be supplied as separate studies.
#' @param stratify Band the source into strata, so a covariate it fitted
#'   contributes a CONTRAST rather than one pooled number. `TRUE` bands over
#'   every covariate the source's own model ESTIMATED a coefficient for, and
#'   leaves the rest marginalised --- derived from the model, so nothing has to
#'   be restated. A covariate the model merely READS is not banded: weight at a
#'   fixed allometric exponent carries no fitted effect to recover, and banding
#'   on it would buy strata and no evidence. Name a covariate explicitly to
#'   override that judgement. Available for published model sources only.
#'   Measured over 720 replicates: coverage 0.933
#'   with one source banded and 0.925 with all of them, against 0.817 with
#'   none --- one banded source is as good as three. See [covStrata()].
#' @param strata_nodes,range Resolution and the enrolled range for `stratify`.
#' @param label Optional display name; otherwise taken from the argument name in
#'   [admStudies()].
#'
#' @return An `admStudy` object.
#' @seealso [admStudies()] to collect several, [admPopulation()] for the
#'   baseline table.
#' @export
admStudy <- function(model = NULL, est = NULL,
                     E = NULL, V = NULL, sd = NULL, sem = NULL,
                     n = NULL, times = NULL, dose = NULL, ev = NULL,
                     population = NULL, at = NULL, by = NULL,
                     stratify = NULL, strata_nodes = NULL, range = NULL,
                     label = NULL, v_denom = NULL) {
  # Record the supplied spread before deriving V or replacing sd from sem.
  .from_spread <- is.null(V) && (!is.null(sd) || !is.null(sem))
  # Use the model expression only in construction-time messages.
  nm  <- label %||% tryCatch({
    .e <- substitute(model)
    if (is.name(.e)) as.character(.e) else "study"
  }, error = function(e) "study")
  bad <- function(...) stop("admixr2: study '", nm, "': ", ..., call. = FALSE)
  has_model <- !is.null(model)
  has_data  <- !is.null(E)
  if (has_model && has_data)
    bad("has BOTH a `model` and digitised `E`. A study contributes in one ",
        "currency: the model it published, or the aggregate data it printed.")
  if (!has_model && !has_data)
    bad("needs either a `model` (with `est`) or digitised `E` (with ",
        "`sd`, `sem` or `V`).")
  if (!has_model && (!is.null(by) ||
      (!is.null(stratify) && !identical(stratify, FALSE)) ||
      !is.null(strata_nodes) || !is.null(range)))
    bad("digitised data cannot be expanded with `by` or `stratify`: one ",
        "reported mean/spread profile contains no separate subgroup profiles. ",
        "Create one admStudy(..., at = ...) per reported subgroup instead.")
  # A cohort data frame supplies its own sample size.
  if (is.data.frame(population) && is.null(n)) n <- nrow(population)
  if (is.null(n) || !is.finite(n) || n <= 0)
    bad("needs a positive `n` -- the number of subjects it reports on.")
  if (is.null(times) || !length(times)) bad("needs `times`.")
  if (is.null(ev) && is.null(dose))
    bad("needs a `dose`, or an `ev` event table for anything a single dose ",
        "cannot express.")
  if (!is.null(ev) && !is.null(dose))
    bad("has both `dose` and `ev`; `dose` is only shorthand for one.")

  ui <- NULL
  if (has_model) {
    ui <- tryCatch(suppressMessages(rxode2::rxode2(model)),
                   error = function(e)
                     bad("`model` could not be parsed: ", conditionMessage(e)))
    known <- ui$iniDf$name
    chk <- function(x, what) {
      if (is.null(x)) return(invisible())
      if (is.null(names(x)) || any(!nzchar(names(x))))
        bad("`", what, "` must be a NAMED vector, e.g. ", what,
            " = c(tcl = ...).")
      unk <- setdiff(names(x), known)
      if (length(unk))
        bad("`", what, "` names ", paste(sQuote(unk), collapse = ", "),
            ", which the model's `ini()` does not declare. Declared: ",
            paste(sQuote(known), collapse = ", "), ".")
    }
    chk(est, "est")
    # Put published estimates directly into the source model.
    if (length(est)) {
      d <- ui$iniDf
      d$est[match(names(est), d$name)] <- unname(est)
      ui$iniDf <- d
    }
  } else {
    if (!is.null(sem)) {
      if (!is.null(sd)) bad("has both `sd` and `sem`; give one.")
      # Convert SEM covariance back to per-subject covariance.
      sd <- as.numeric(sem) * sqrt(n)
    }
    if (is.null(V) && is.null(sd))
      bad("digitised data needs a spread: `sd` per timepoint, `sem`, or a ",
          "full `V`.")
    if (is.null(V)) V <- as.numeric(sd)^2
    if (length(E) != length(times))
      bad("`E` has ", length(E), " values but `times` has ", length(times), ".")
  }
  # Published sd/sem uses n-1; supplied V and model moments use ML.
  if (!is.null(v_denom)) {
    if (!is.character(v_denom) || length(v_denom) != 1L ||
        !v_denom %in% c("ml", "unbiased"))
      bad("`v_denom` must be \"ml\" or \"unbiased\".")
    if (has_model && identical(v_denom, "unbiased"))
      bad("`v_denom = \"unbiased\"` does not apply to a published `model`: ",
          "its moments are generated at the ML denominator, so there is ",
          "nothing to convert.")
  } else {
    v_denom <- if (.from_spread) "unbiased" else "ml"
  }
  # Keep every non-reserved cohort column; cross-study covariates may identify
  # effects absent from the source model.
  if (is.data.frame(population)) {
    drop <- names(population)[toupper(names(population)) %in% .ADM_DATA_COLS]
    keep <- setdiff(names(population), drop)
    if (!length(keep))
      bad("`population` is a data frame with no covariate columns -- it has ",
          "only ", paste(sQuote(names(population)), collapse = ", "),
          ", which are data columns rather than covariates.")
    if (length(drop))
      message("admixr2: study ", sQuote(nm), ": reading ",
              paste(sQuote(keep), collapse = ", "), " from `population`; ",
              paste(sQuote(drop), collapse = ", "),
              if (length(drop) == 1L) " is a data column, not a covariate."
              else " are data columns, not covariates.")
    population <- admPopulation(data = population[, keep, drop = FALSE])
  }
  overlap <- intersect(names(at), .admCovSpecNames(population))
  if (length(overlap))
    bad("`at` pins ", paste(sQuote(overlap), collapse = ", "),
        ", but `population` also gives ",
        if (length(overlap) == 1L) "it a distribution" else "them distributions",
        ". Remove the pinned covariate from `population`; conditioning a ",
        "dependent population requires its conditional distribution.")
  structure(list(
    ui = ui, model = model,
    E = E, V = V, n = as.numeric(n), times = as.numeric(times),
    ev = ev, dose = dose, population = population, at = at, by = by,
    stratify = stratify, strata_nodes = strata_nodes, range = range,
    label = label, v_denom = v_denom), class = "admStudy")
}

#' @export
print.admStudy <- function(x, ...) {
  nm <- x$label %||% "study"
  kind <- if (is.null(x$ui)) "digitised DATA" else "published MODEL"
  cat("admixr2 study '", nm, "'  -- ", kind, "\n", sep = "")
  dz <- if (!is.null(x$dose)) paste0("dose ", x$dose) else "ev supplied"
  cat(sprintf("  design    n = %s, %s, %d times (%g - %g)\n",
              format(x$n), dz, length(x$times), min(x$times), max(x$times)))
  # Show resolved defaults in the collection pre-flight, not per study.
  if (is.null(x$ui))
    cat(sprintf("  V scale   %s denominator%s\n", x$v_denom %||% "ml",
                if (identical(x$v_denom, "unbiased"))
                  "  (published spread; converted to ML for the fit)" else ""))
  if (!is.null(x$ui)) {
    ini <- x$ui$iniDf
    cvs <- tryCatch(x$ui$allCovs, error = function(e) character(0))
    cat(sprintf("  model     %d estimated parameter%s%s\n",
                sum(!ini$fix), if (sum(!ini$fix) == 1L) "" else "s",
                if (length(cvs)) paste0("; reads ", paste(cvs, collapse = ", "))
                else ""))
    cat("  reported  a published model; no standard error is available for a ",
        "fit that includes one\n", sep = "")
  } else {
    cat(sprintf("  reported  mean profile, %s\n",
                if (is.matrix(x$V)) "full covariance" else "per-time spread"))
  }
  if (!is.null(x$population)) {
    pn <- .admCovSpecNames(x$population)
    cat("  population", paste(pn, collapse = ", "), "\n")
    if (!is.null(x$ui)) {
      cvs <- tryCatch(x$ui$allCovs, error = function(e) character(0))
      marg <- setdiff(pn, cvs)
      if (length(marg))
        cat("            ", paste(marg, collapse = ", "),
            " not in this model -> marginalised (no contrast from this source)\n",
            sep = "")
    }
  }
  if (!is.null(x$at)) cat("  pinned at ",
    paste(sprintf("%s = %s", names(x$at), unlist(x$at)), collapse = ", "), "\n")
  if (!is.null(x$by))       cat("  reported by", x$by, "-> one study per level\n")
  if (!is.null(x$stratify)) cat("  banded on ", x$stratify, "\n")
  invisible(x)
}

#' Collect the studies a meta-analysis draws on
#'
#' Takes [admStudy()] objects and names them. Names come from the argument names
#' where given, so `admStudies(smith2019, jones2021)` labels them by the objects
#' they were built as.
#'
#' Nothing is generated here either. The studies are materialised once, inside
#' the fit, so the whole specification stays cheap to build and inspect.
#'
#' @param ... [admStudy()] objects, optionally named.
#' @return An `admStudies` object, to pass as `studies` to any of
#'   [admControl()], [adghControl()], [adfoControl()] or [adirmcControl()].
#' @seealso [admStudy()], [admPopulation()].
#' @export
admStudies <- function(...) {
  a <- list(...)
  if (length(a) == 1L && inherits(a[[1L]], "admStudies")) return(a[[1L]])
  if (length(a) == 1L && is.list(a[[1L]]) && !inherits(a[[1L]], "admStudy"))
    a <- a[[1L]]
  ok <- vapply(a, inherits, logical(1), "admStudy")
  if (!all(ok))
    stop("admixr2: `admStudies()` takes admStudy() objects; element",
         if (sum(!ok) > 1L) "s " else " ",
         paste(which(!ok), collapse = ", "), " ",
         if (sum(!ok) > 1L) "are" else "is", " not one.", call. = FALSE)
  nms <- names(a)
  # Derive labels only when `...` directly contains the studies.
  auto <- if (length(a) == length(unwrapped <- as.list(substitute(list(...)))[-1L]))
    vapply(unwrapped, function(e) if (is.name(e)) as.character(e) else "",
           character(1))
  else rep("", length(a))
  if (is.null(nms)) nms <- rep("", length(a))
  nms[is.na(nms)] <- ""
  for (i in seq_along(a)) {
    if (!nzchar(nms[i])) nms[i] <- a[[i]]$label %||% auto[i]
    if (is.na(nms[i]) || !nzchar(nms[i])) nms[i] <- paste0("study", i)
    a[[i]]$label <- nms[i]
  }
  if (anyDuplicated(nms))
    stop("admixr2: study names must be unique; repeated: ",
         paste(sQuote(unique(nms[duplicated(nms)])), collapse = ", "), ".",
         call. = FALSE)
  structure(stats::setNames(a, nms), class = "admStudies")
}

#' @export
print.admStudies <- function(x, ...) {
  n_src <- sum(vapply(x, function(s) !is.null(s$ui), logical(1)))
  plural <- function(n, one, many) paste0(n, if (n == 1L) one else many)
  cat("admixr2 studies: ", length(x), "  (",
      plural(n_src, " published model, ", " published models, "),
      length(x) - n_src, " digitised)
", sep = "")
  for (nm in names(x)) {
    s <- x[[nm]]
    cat(sprintf("  %-14s %s  n = %-6s %d times%s
", nm,
                if (is.null(s$ui)) "data " else "model", format(s$n),
                length(s$times), ""))
  }

  # Identifiability is a property of the complete study set.
  role <- function(s, cv) {
    band <- .admStudyBandNames(s)
    if (cv %in% c(names(s[["at"]]), as.character(s[["by"]]))) "conditioned"
    else if (cv %in% band)                                    "banded"
    else if (cv %in% .admCovSpecNames(s[["population"]]))     "marginal"
    else                                                      "-"
  }
  allcov <- unique(unlist(lapply(x, function(s)
    .admCovSpecNames(s[["population"]]))))
  if (length(allcov)) {
    cat("
covariate      ",
        paste(sprintf("%-10s", substr(names(x), 1L, 10L)), collapse = ""),
        "
", sep = "")
    flat <- character(0); fit_by <- character(0)
    for (cv in allcov) {
      r <- vapply(x, role, character(1), cv = cv)
      if (all(r %in% c("marginal", "-"))) {
        flat <- c(flat, cv)
        if (any(vapply(x, function(s)
              length(.admCovCoefThetas(s$ui, cv, s[["population"]])) > 0L,
              logical(1))))
          fit_by <- c(fit_by, cv)
      }
      cat(sprintf("  %-13s", cv), paste(sprintf("%-10s", r), collapse = ""),
          "
", sep = "")
    }
    if (length(flat))
      cat("
NOTE: ", paste(sQuote(flat), collapse = ", "),
          if (length(flat) == 1L)
            " is marginal in every source, so its effect is identified ONLY"
          else
            " are marginal in every source, so their effects are identified ONLY",
          " by the contrast BETWEEN sources -- no source",
          " carries a within-source contrast to separate it from a random",
          " effect on the same parameter. Legitimate when the sources really",
          " do differ; a flat ridge when they do not, and a flat ridge",
          " converges to a confident wrong number rather than failing.",
          # Suggest only bands .admExpandStrata() can construct.
          if (length(fit_by))
            paste0(" Band the source that fitted it (`stratify = \"",
                   fit_by[1L], "\"`).")
          else
            paste0(" No source fitted ",
                   if (length(flat) == 1L) "it" else "any of them",
                   ", so banding is not available and this is the",
                   " between-source contrast and nothing more."),
          "
", sep = "")
  }
  src <- names(x)[vapply(x, function(s) !is.null(s$ui), logical(1))]
  if (length(src))
    cat("
NOTE: ", paste(sQuote(src), collapse = ", "),
        if (length(src) == 1L) " contributes as a published MODEL, weighted"
        else " contribute as published MODELS, weighted",
        " as if `n` patients had been sampled -- `n` sets RELATIVE WEIGHT",
        " against the other studies, not precision. No standard error is",
        " reported for a fit that includes one.
",
        sep = "")
  cat("
print() a single study to check its transcription.
")
  invisible(x)
}

# Materialise lazy specs once at the shared driver entry point.
.admMaterialise <- function(studies) {
  if (inherits(studies, "admStudies")) studies <- unclass(studies)
  if (!is.list(studies)) return(studies)
  spec <- vapply(studies, inherits, logical(1), "admStudy")
  if (!any(spec)) return(studies)
  out <- list()
  add <- function(nm, value) {
    if (nm %in% names(out))
      stop("admixr2: materialising the studies produced duplicate name ",
           sQuote(nm), ". Rename the study whose `by` or `stratify` expansion ",
           "collides with it.", call. = FALSE)
    out[[nm]] <<- value
  }
  for (nm in names(studies)) {
    s <- studies[[nm]]
    if (!inherits(s, "admStudy")) { add(nm, s); next }
    ev <- s$ev %||% rxode2::et(amt = s$dose)
    if (is.null(s$ui)) {
      # Digitised studies need no generation.
      g <- list(E = as.numeric(s$E), V = s$V, n = s$n, times = s$times,
                ev = ev, v_denom = s$v_denom)
      if (!is.null(s[["population"]])) g[["cov_dist"]] <- s[["population"]]
      if (!is.null(s[["at"]]))         g[["cov"]]      <- s[["at"]]
      add(nm, g)
      next
    }
    # Use exact indexing: `$cov` would partially match `cov_dist`.
    sp <- list(times = s$times, ev = ev, n = s$n)
    if (!is.null(s[["population"]])) sp[["cov_dist"]] <- s[["population"]]
    if (!is.null(s[["at"]]))         sp[["cov"]]      <- s[["at"]]
    if (!is.null(s$stratify)) {
      sp$stratify <- s$stratify
      if (!is.null(s$strata_nodes)) sp$strata_nodes <- s$strata_nodes
      # Resolve `stratify = TRUE` before assigning an unnamed range.
      if (!is.null(s$range)) {
        sp$cov_range <- if (is.list(s$range) && !is.null(names(s$range)))
          s$range
        else {
          .bn <- .admStudyBandNames(s)
          if (length(.bn) != 1L)
            stop("admixr2: study '", nm, "': `range` does not say which ",
                 "covariate it is the enrolled range of, and this source ",
                 "bands ", if (!length(.bn)) "none" else length(.bn),
                 ". Give a named list, e.g. range = list(WT = c(52, 118)).",
                 call. = FALSE)
          stats::setNames(list(s$range), .bn)
        }
      }
    }
    # `by` creates one ordinary pinned study per reported level.
    if (!is.null(s$by)) {
      lv <- s$population[[s$by]][["values"]]
      if (is.null(lv))
        stop("admixr2: study '", nm, "': `by = \"", s$by, "\"` needs that ",
             "covariate declared with levels in `population`, e.g. ",
             s$by, " = c(male = 0.55).", call. = FALSE)
      pr <- s$population[[s$by]][["probs"]]
      pr <- if (is.null(pr)) rep(1 / length(lv), length(lv)) else pr / sum(pr)
      for (k in seq_along(lv)) {
        spk <- sp
        spk$n   <- s$n * pr[k]
        spk[["cov"]] <- c(sp[["cov"]], stats::setNames(list(lv[k]), s$by))
        spk$cov_dist <- .admCovDropMargin(s$population, s$by)
        g <- datagen(stats::setNames(list(spk), paste0(nm, "_", s$by, lv[k])),
                     model = s$ui, control = datagenControl(method = "gh"))
        # `by` plus `stratify` may yield several bands per level.
        for (kk in names(g)) add(kk, g[[kk]])
      }
      next
    }
    g <- datagen(stats::setNames(list(sp), nm), model = s$ui,
                 control = datagenControl(method = "gh"))
    for (k in names(g)) add(k, g[[k]])
  }
  out
}

# Covariates `stratify` can band. `TRUE` means effects estimated by the source,
# not covariates merely read by it.
.admStudyBandNames <- function(s) {
  st <- s[["stratify"]]
  if (is.null(st) || identical(st, FALSE)) return(character(0))
  if (!isTRUE(st)) return(as.character(st))
  Filter(function(cv)
           length(.admCovCoefThetas(s$ui, cv, s[["population"]])) > 0L,
         intersect(.admCovSpecNames(s[["population"]]),
                   tryCatch(s$ui$allCovs, error = function(e) character(0))))
}
