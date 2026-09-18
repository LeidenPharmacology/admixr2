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
  missing_msg <- function(nm)
    bad(nm, "has missing or non-finite values. A margin has to describe ",
        "every subject the study reports on -- drop the column, or state ",
        "it yourself.")
  for (nm in nms) {
    v <- data[[nm]]
    # Reject NAs before sort() or factor operations.
    if (anyNA(v)) missing_msg(nm)
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
    if (!all(is.finite(v))) missing_msg(nm)
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

  # Report unsupported discrete-continuous correlation in cohort data.
  disc <- intersect(
    names(resolved)[vapply(resolved, function(s) !is.null(s[["values"]]),
                           logical(1))],
    names(data))
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
              "rest. Condition on the discrete one (`stratify`) if the association ",
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
  pair_key <- function(x) {
    if (length(a) < 2L) return(x)
    pairs <- utils::combn(names(a), 2L, simplify = FALSE)
    hit <- vapply(pairs, function(p)
      identical(x, paste(p, collapse = ".")) ||
      identical(x, paste(rev(p), collapse = ".")), logical(1))
    if (sum(hit) == 1L) paste(sort(pairs[[which(hit)]]), collapse = "\r") else x
  }
  if (!is.null(data)) {
    d <- .admPopFromData(data, dist, a)
    a <- c(a, d$specs)
    # Explicit correlations override cohort values pairwise.
    if (is.null(cor)) { if (length(d$cor)) cor <- d$cor }
    else if (!is.matrix(cor)) {
      # Canonicalise unordered pairs without splitting dotted covariate names.
      user_keys <- vapply(names(cor), pair_key, character(1))
      data_keys <- vapply(names(d$cor), pair_key, character(1))
      cor <- c(cor, d$cor[!data_keys %in% user_keys])
    } else if (length(d$cor)) {
      # Full matrix overrides all cohort correlations.
      message("admixr2: `cor` is a full matrix, so the correlations the ",
              "cohort itself shows (",
              paste(sprintf("%s = %.2f", names(d$cor), d$cor),
                    collapse = ", "),
              ") are REPLACED rather than merged. Give `cor` as a named ",
              "vector of the pairs you mean to override if that is not what ",
              "you want.")
    }
  }
  if (!length(a) || is.null(names(a)) || any(!nzchar(names(a))))
    stop("admixr2: `admPopulation()` needs NAMED covariates, e.g. ",
         "admPopulation(WT = c(mean = 75, sd = 16)), or a `data` frame to ",
         "derive them from.", call. = FALSE)
  if (anyDuplicated(names(a)))
    stop("admixr2: `admPopulation()` covariate names must be unique; repeated: ",
         paste(sQuote(unique(names(a)[duplicated(names(a))])), collapse = ", "),
         ".", call. = FALSE)
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
      # One set of names, however they were spelled.
      nmR <- rownames(R) %||% colnames(R)
      if (!is.null(nmR)) dimnames(R) <- list(nmR, nmR)
      if (is.null(nmR)) {
        # Require dimnames when matrix cor is used alongside data.
        if (!is.null(data))
          stop("admixr2: a matrix `cor` given alongside `data` must have ",
               "dimnames -- the covariates are ordered as declared and then ",
               "as derived from the data frame (",
               paste(sQuote(nms), collapse = ", "),
               "), which is not an order to rely on positionally.",
               call. = FALSE)
        dimnames(R) <- list(nms, nms)
      } else {
        if (!all(nms %in% nmR))
          stop("admixr2: matrix `cor` is missing ",
               paste(sQuote(setdiff(nms, nmR)), collapse = ", "),
               ". Name every declared covariate. Declared: ",
               paste(sQuote(nms), collapse = ", "), ".", call. = FALSE)
        R <- R[nms, nms, drop = FALSE]
      }
      # Refuse correlated discrete margins in matrix cor.
      off <- nms[lv][vapply(nms[lv], function(x)
        any(abs(R[x, setdiff(nms, x)]) > 0), logical(1))]
      if (length(off))
        stop("admixr2: matrix `cor` correlates ",
             paste(sQuote(off), collapse = " and "),
             ", which is a DISCRETE covariate. A level would then be a ",
             "truncation of the latent normal rather than a point, which ",
             "admixr2 refuses -- see covStrata(). Declare it independent.",
             call. = FALSE)
    } else {
      R <- diag(length(nms)); dimnames(R) <- list(nms, nms)
      cn <- names(cor)
      if (is.null(cn) || any(!nzchar(cn)))
        stop("admixr2: `cor` must name the PAIR it applies to, e.g. ",
             "cor = c(WT.CRCL = 0.45).", call. = FALSE)
      keys <- vapply(cn, pair_key, character(1))
      if (anyDuplicated(keys))
        stop("admixr2: `cor` names the same covariate pair more than once: ",
             paste(sQuote(cn[duplicated(keys)]), collapse = ", "), ".",
             call. = FALSE)
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
#'   a table, taking `n` from its row count. EVERY column that is not a
#'   reserved data column (`ID`, `TIME`, `DV`, `AMT`, ...) becomes a covariate,
#'   including ones this model does not read --- a covariate another source
#'   fitted is still evidence here, and each one adds a quadrature axis, so
#'   drop the columns you do not want integrated over.
#' @param at Named list pinning a covariate at a value, for a study reported in
#'   one subgroup, e.g. `at = list(SEX = 1)`. A pinned covariate must be omitted
#'   from `population`.
#' @param by Covariate the paper reports results SEPARATELY by, e.g.
#'   `by = "SEX"`. Expands a published model into one study per level; digitised
#'   subgroup profiles must be supplied as separate studies.
#' @param strata_nodes Quadrature resolution for a conditional CONTINUOUS
#'   covariate, i.e. how many nodes it is cut into. Default 9. Raise it for
#'   precision, at a multiplicative cost: a source conditional on two continuous
#'   covariates becomes `strata_nodes^2` studies, and each is solved. Discrete
#'   covariates are unaffected --- their levels are exact and cost one stratum
#'   each. Note that the objective VALUE depends on this, so two fits are only
#'   comparable at the same resolution; admixr2 records it per study and
#'   [anova()] checks it.
#' @param range Optional named list giving the covariate span the source
#'   ENROLLED, e.g. `range = list(WT = c(52, 118))`. The declared distribution
#'   is truncated to it, whether that covariate ends up conditional (strata are
#'   then cut from the truncated margin) or marginal (the truncated margin is
#'   integrated over). Needed when that distribution is
#'   wider than the enrolment --- a `mean +/- SD` transcribed from a baseline
#'   table has tails where nobody was --- and *not* wanted when `population` is
#'   the patients themselves, since the margins are then fitted to exactly who
#'   was enrolled and truncating would score the paper against a sub-population
#'   of its own. An unnamed value is allowed only when the source itself is
#'   conditional on exactly one covariate --- which covariate an unnamed range
#'   belongs to is a question about the source, so the answer does not change
#'   with the model being fitted to it.
#' @param label Optional display name; otherwise taken from the argument name in
#'   [admStudies()].
#'
#' @section Conditional and marginal covariates are derived, not declared:
#'
#' Whether a covariate is **conditional** or **marginal** for a source is a
#' property of that source's own model, so admixr2 works it out and there is
#' nothing to set.
#'
#' A covariate is **conditional** when the source's model ESTIMATED its
#' coefficient. That paper reports a contrast along it, so the source is cut
#' into nodes and the contrast is carried into the fit. A covariate the model
#' merely READS is not conditional --- weight at a fixed allometric exponent carries
#' no fitted effect to recover, and conditioning on it would buy nodes and no
#' evidence.
#'
#' A covariate is **marginal** when the model did not estimate it. There is no
#' contrast in that paper to condition on, and splitting it on a covariate its model
#' never saw would manufacture evidence, so admixr2 integrates over the
#' population instead. That covers both the covariate the model never mentions
#' and the one it reads at an asserted coefficient; [print()] names them
#' separately, because only the second is easy to mistake for conditional.
#'
#' Which covariates a source is conditional on comes from that source. Which of
#' them admixr2 has to cut into nodes is narrowed to the ones the **analysis**
#' model reads, and that narrowing is exact: if the model's prediction does not
#' move across a source's nodes, the mixture those nodes collapse to is a
#' sufficient statistic for it, so the objective is unchanged to eight decimal
#' places while the study count falls. Measured on a source conditional on two
#' covariates with an analysis model reading one, collapsing the unread one
#' moved the objective by 0.00004 and collapsing the read one by 73.6.
#'
#' This matters and is why it is not left to the caller: a covariate every
#' source marginalises is not identified against a random effect on the same
#' parameter. Measured over 720 replicates, coverage was 0.933 with one source
#' conditional and 0.925 with all of them conditional, against **0.817 with
#' none**.
#'
#' A source given as digitised `E`/`V` has no model to derive from and is
#' marginal over whatever it declares. Use one `admStudy(..., at = ...)` per
#' reported subgroup to enter a paper that published by subgroup.
#'
#' `stratify` was removed; [covStrata()] still takes it for working with a
#' covariate distribution directly. `strata_nodes` and `range` remain: the
#' first is a precision setting and the second is a fact about the source,
#' neither of which is a statement about which covariates are conditional.
#'
#' @return An `admStudy` object.
#' @seealso [admStudies()] to collect several, [admPopulation()] for the
#'   baseline table.
#' @export
admStudy <- function(model = NULL, est = NULL,
                     E = NULL, V = NULL, sd = NULL, sem = NULL,
                     n = NULL, times = NULL, dose = NULL, ev = NULL,
                     population = NULL, at = NULL, by = NULL,
                     strata_nodes = NULL, range = NULL,
                     label = NULL, v_denom = NULL) {
  # WHETHER A COVARIATE IS CONDITIONAL OR MARGINAL IS NOT A USER CHOICE. It is a
  # property of the source's own model -- conditional when that model ESTIMATED
  # the covariate's coefficient, marginal when it never mentions the covariate
  # or reads it at an asserted one -- so admixr2 derives it and there is nothing
  # to declare.
  # `strata_nodes` and `range` STAY, because neither is that choice: one is a
  # precision setting and the other is a fact about the source. It is the span the source
  # ENROLLED, and it is needed exactly when the declared distribution is wider
  # than the enrolment -- a `mean +/- SD` transcribed from a baseline table has
  # tails where nobody was, and strata cut from it credit the paper there.
  # Measured on a lognormal declared at meanlog log(76), sdlog 0.30: cut from
  # the full distribution the nodes span 20-294 kg and the mixture reproduces
  # E[(WT/70)^0.75] = 1.09089, the declared population's value; cut from the
  # 10th-90th they span 52-112 and reproduce 1.07533, the truncated
  # population's. Both are exact quadratures -- of DIFFERENT populations, and
  # only one of them is the one the paper reported moments for.
  #
  # `stratify` IS NOT A FIELD HERE. It stays an internal field of the plain
  # study specs .admMaterialise() builds, where covStrata(), the internal
  # callers and a deliberately uncut reference fit set it; an `admStudy`
  # object carries the model the derivation reads instead. It used to be
  # carried as a hardcoded NULL, which made every reader of it on this object
  # -- print()'s "conditional on nothing" line among them -- unreachable code that
  # looked live.
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
  if (has_model && any(!vapply(list(E, V, sd, sem), is.null, logical(1))))
    bad("has BOTH a `model` and digitised summary fields (`E`, `V`, `sd` or ",
        "`sem`). A study contributes in one ",
        "currency: the model it published, or the aggregate data it printed.")
  if (!has_model && !has_data)
    bad("needs either a `model` (with `est`) or digitised `E` (with ",
        "`sd`, `sem` or `V`).")
  if (!has_model && !is.null(est))
    bad("digitised data cannot have `est`; estimates belong to a published ",
        "`model` source.")
  if (!has_model && !is.null(by))
    bad("digitised data cannot be expanded with `by`: one reported mean/spread ",
        "profile contains no separate subgroup profiles. Create one ",
        "admStudy(..., at = ...) per reported subgroup instead.")
  if (!has_model && (!is.null(range) || !is.null(strata_nodes)))
    bad("digitised data cannot take ",
        paste(sQuote(c("range", "strata_nodes")[
          c(!is.null(range), !is.null(strata_nodes))]), collapse = " or "),
        ": both bear only on the strata a published model is cut into, and ",
        "there is no model here to condition.")
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
  if (length(at)) {
    if (is.null(names(at)) || any(!nzchar(names(at))) || anyDuplicated(names(at)))
      bad("`at` must have unique, non-empty covariate names, e.g. ",
          "at = list(SEX = 1).")
    if (!all(vapply(as.list(at), function(x)
      length(x) == 1L && (is.numeric(x) || is.logical(x)) && is.finite(x),
      logical(1))))
      bad("every `at` value must be one finite number.")
  }

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
    if (anyDuplicated(names(est)))
      bad("`est` parameter names must be unique; repeated: ",
          paste(sQuote(unique(names(est)[duplicated(names(est))])),
                collapse = ", "), ".")
    if (length(est) && (!is.numeric(est) || any(!is.finite(est))))
      bad("`est` must contain finite numeric parameter estimates.")
    # Put published estimates directly into the source model.
    if (length(est)) {
      d <- ui$iniDf
      d$est[match(names(est), d$name)] <- unname(est)
      ui$iniDf <- d
    }
  } else {
    spread <- if (is.null(V)) sem %||% sd else NULL
    if (!is.null(spread) &&
        (!is.numeric(spread) || any(!is.finite(spread)) || any(spread < 0)))
      bad("`sd` and `sem` must contain finite, non-negative numbers.")
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
  # Extract non-reserved columns as candidate covariates.
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
  # Canonicalise population so plain lists match admPopulation() structures.
  if (!is.null(population))
    population <- tryCatch(.admCovDistCanon(population),
                           error = function(e) bad("`population` is not a ",
                                                   "valid covariate ",
                                                   "specification: ",
                                                   conditionMessage(e)))
  # THE SECOND DOOR. covDist() refuses `joint`, but `population` also takes a
  # plain list, and the canon lets one through untouched. `jointOwn` is what
  # separates the sampler admixr2 BUILT from `cor` -- which is the supported
  # path and arrives here on every correlated population -- from one the caller
  # wrote.
  if (is.function(population[["joint"]]) &&
      !isTRUE(population[["jointOwn"]]))
    bad("`population` carries its own `joint` sampler, which is not accepted ",
        "yet. Declare the margins and give `cor` instead; an arbitrary ",
        "sampler -- a vine copula among them -- is planned, not released.")
  if (!is.null(by)) {
    if (!is.character(by) || length(by) != 1L || is.na(by) || !nzchar(by))
      bad("`by` must be one non-empty covariate name.")
    declared <- .admCovSpecNames(population)
    if (!by %in% declared)
      bad("`by` names ", sQuote(by), ", which `population` does not declare. ",
          "Declared: ", paste(sQuote(declared), collapse = ", "), ".")
    if (is.null(population[[by]][["values"]]))
      bad("`by` must name a discrete population margin with declared levels.")
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
    strata_nodes = strata_nodes, range = range,
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
  # What the fit will condition on, resolved the same way .admMaterialise() resolves it.
  .bn <- .admStudyBandNames(x)
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
      # READ AT AN ASSERTED COEFFICIENT is the other way to be marginal, and
      # the one a reader would otherwise have to work out: the model does
      # mention the covariate, so the line above says nothing about it.
      .fx <- setdiff(intersect(pn, cvs), .bn)
      if (length(.fx))
        cat("            ", paste(.fx, collapse = ", "),
            " read at an asserted coefficient -> marginalised (no fitted ",
            "effect to recover)\n", sep = "")
    }
  }
  if (!is.null(x$at)) cat("  pinned at ",
    paste(sprintf("%s = %s", names(x$at), unlist(x$at)), collapse = ", "), "\n")
  if (!is.null(x$by))       cat("  reported by", x$by, "-> one study per level\n")
  if (length(.bn)) cat("  conditional on ", paste(.bn, collapse = ", "), "\n")
  else if (!is.null(x$ui) && !is.null(x$population))
    cat("  conditional on  nothing -- this model estimated no coefficient for a ",
        "covariate this source describes\n", sep = "")
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
    cond <- .admStudyBandNames(s)
    if (cv %in% c(names(s[["at"]]), as.character(s[["by"]]))) "conditioned"
    else if (cv %in% cond)                                    "conditional"
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
          # Conditioning is derived, so there is nothing for the caller to set --
          # this reports what the derivation found rather than suggesting a
          # flag. A source that FITTED the covariate is conditional on it already.
          if (length(fit_by))
            paste0(" ", sQuote(fit_by[1L]),
                   " did fit it, so that source is conditional on it and the",
                   " contrast is carried; this warns because the others are",
                   " not.")
          else
            paste0(" No source fitted ",
                   if (length(flat) == 1L) "it" else "any of them",
                   ", so conditioning is not available and this is the",
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
## WHICH COVARIATES A SOURCE IS CONDITIONAL ON is a property of that source: its
## own model estimated their coefficients, so its paper reports a contrast along
## them. WHICH OF THEM ADMIXR2 HAS TO CUT INTO NODES is a different question,
## and the answer is: only the ones the ANALYSIS model can be moved by.
##
## `analysis_covs` is that narrowing, and it is exact rather than a saving with
## a cost. Measured at identical parameters, on a source conditional on CRCL and
## WT with an analysis model reading WT only:
##
##   both conditional        50 studies   OFV -3041.72602426
##   only WT conditional     10 studies   OFV -3041.72606427   diff -0.00004
##   only CRCL conditional   10 studies   OFV -3115.32302414   diff -73.59700
##
## Collapsing the direction the analysis CANNOT see costs four decimal places;
## collapsing one it can see costs 73.6. The invariance is exact rather than
## approximate -- against a fully marginal reference it is -0.158 at 3 nodes,
## +0.00003 at 5 and -0.00000001 at 9 -- because if a model's prediction does
## not move across a source's nodes, the mixture collapse those nodes reduce to
## is a sufficient statistic for it. That is the same law .admMixMoments()
## applies, so this is the exactly removable part of the work and nothing else.
.admMaterialise <- function(studies, analysis_covs = NULL) {
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
    # RESOLVED HERE, and nowhere else: what reaches datagen() is the set of
    # covariates this SOURCE estimated a coefficient for, narrowed to the ones
    # the ANALYSIS model reads -- see the note above for why the narrowing is
    # exact.
    # An explicit `stratify` on the SPEC still wins. It is not an `admStudy()`
    # argument and not a field of the object -- see the note in admStudy() --
    # but it stays an internal field that covStrata(), the internal callers and
    # a deliberately UNBANDED reference fit can set, and `FALSE` has to reach
    # the spec or the opt-out silently becomes its opposite.
    # TWO SETS, and the difference between them is the whole point of the
    # narrowing. `.src` is what this SOURCE estimated a coefficient for -- a
    # property of the paper, true whatever is fitted to it. `.bn` is that set
    # narrowed to what the ANALYSIS model reads, which is what gets cut into
    # nodes. Questions about the source are answered from `.src`; only the
    # cutting uses `.bn`.
    .src <- if (!is.null(s[["stratify"]])) s[["stratify"]]
            else .admStudyBandNames(s)
    .bn <- if (identical(.src, FALSE) || is.null(analysis_covs)) .src
           else intersect(.src, analysis_covs)
    # THE ANSWER REACHES THE SPEC EITHER WAY, `FALSE` included. An empty
    # conditional set used to record nothing, and .admExpandStrata() then
    # derived one of its own -- so a source that resolved to nothing came back
    # cut into nodes
    # anyway, and `stratify = FALSE` silently became its opposite.
    sp$stratify <- if (identical(.bn, FALSE) || !length(.bn)) FALSE else .bn
    # OUTSIDE the branch, so the spec carries what the caller asked for however
    # the derivation came out. Copied only when something was conditional, a study
    # that resolved to nothing took the default 9 while its neighbour took the
    # caller's setting -- and `.adm_strata_nodes` varying between studies is
    # exactly what anova() refuses to compare across.
    if (!is.null(s$strata_nodes)) sp$strata_nodes <- s$strata_nodes
    # The enrolled span truncates the declared distribution WHATEVER becomes of
    # the covariate afterwards: strata cut from a conditional one are cut from
    # the truncated margin, and a marginal one is integrated over the truncated
    # margin. `range` says which patients the source had, and that is true
    # before admixr2 decides how to use them -- so this is deliberately outside
    # the conditioning branch above. An unnamed value still needs exactly one
    # covariate to be the range OF; otherwise there is nothing to attach it to.
    #
    # FROM `.src`, NOT `.bn`: which covariate an unnamed range belongs to is a
    # question about the SOURCE, and it is asked before the analysis model has
    # any say. Read from the narrowed set, one `studies` object gave a source
    # conditional on a single covariate a clean fit under a model that reads it
    # and a hard error under the null that drops it -- which is precisely the
    # nested pair this PR is built around -- and, for a source conditional on
    # two, silently attached a weight range to whichever one happened to
    # survive. .admCovSourceRange() already asks the un-narrowed question, so
    # the fit and the plot disagreed about the same range.
    if (!is.null(s$range))
      sp$cov_range <- if (is.list(s$range) && !is.null(names(s$range)))
        s$range
      else {
        .one <- if (identical(.src, FALSE)) character(0) else .src
        if (length(.one) != 1L)
          stop("admixr2: study '", nm, "': `range` does not say which ",
               "covariate it is the enrolled range of, and this source conditions on ",
               length(.one),
               if (length(.one)) paste0(" (", paste(sQuote(.one),
                                                    collapse = ", "), ")"),
               ". Give a named list, e.g. range = list(WT = c(52, 118)).",
               call. = FALSE)
        stats::setNames(list(s$range), .one)
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
        # `by` plus `stratify` may yield several nodes per level.
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

# Covariates this source's model ESTIMATED a coefficient for, and whose
# distribution it also declares -- the ones conditioning can extract a contrast from.
#
# Read by .admMaterialise() at fit time and by print.admStudy() beforehand, so
# what a study prints is what the fit does. Without the second reader a study
# printed every covariate as marginal and then came back conditional.
.admStudyBandNames <- function(s) {
  st <- s[["stratify"]]
  if (identical(st, FALSE)) return(character(0))
  if (!is.null(st) && !isTRUE(st)) return(as.character(st))
  if (is.null(s$ui)) return(character(0))
  # `by` IS RESOLVED HERE, not at the call site, so print() cannot claim a
  # covariate the fit will not cut on: `by` pins its covariate one level per
  # study, each level's `population` stops declaring it, and there is nothing
  # left to cut. Resolved in .admMaterialise() alone, a source with
  # `by = "SEX"` whose model estimates a sex effect printed "conditional on
  # SEX" and then came back from the fit with SEX pinned instead.
  # ESTIMATED, not merely READ, which is the rule the man page states and the
  # one `stratify = TRUE` has always enforced. A covariate the model reads at an
  # ASSERTED coefficient -- weight at a fixed allometric exponent -- carries no
  # fitted effect to recover: conditioning on it buys nodes and no evidence, and
  # credits the source with information it never earned. One the population
  # declares and the model never mentions is marginal; one the model estimates
  # but the population never described has no distribution to condition on.
  cvs <- intersect(.admCovSpecNames(s[["population"]]),
                   tryCatch(s$ui$allCovs, error = function(e) character(0)))
  if (!length(cvs)) return(character(0))
  cvs <- setdiff(cvs, s[["by"]] %||% character(0))
  if (!length(cvs)) return(character(0))
  cvs[vapply(cvs, function(cv) length(tryCatch(
    .admCovCoefThetas(s$ui, cv, s[["population"]]),
    error = function(e) character(0))) > 0L, logical(1))]
}
