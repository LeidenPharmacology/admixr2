# =============================================================================
# Writing down a published study
# =============================================================================
#
# What a user actually has is a PAPER: a parameter table with %RSE, a baseline
# demographics table, and a design. Not an nlmixr2 fit object. Everything here
# is shaped so that transcribing those three tables is the whole job, and so
# that the transcription can be CHECKED before anything is solved.
#
# The alternative -- raw nested lists with magic field names, and a model whose
# ini() reads global variables that a helper reassigns with `<<-` -- is what the
# covariates vignette used to do. It has no autocomplete, no validation until a
# fit fails, and it makes the reader hold five moving parts at once.

# Column names rxode2/nlmixr2 reserve for the data itself, so they can never
# be covariates. Everything else in a supplied data frame is one.
.ADM_DATA_COLS <- c("ID", "TIME", "DV", "AMT", "EVID", "CMT", "DVID", "MDV",
                    "RATE", "DUR", "SS", "II", "ADDL", "CENS", "LIMIT")

.admPopFromData <- function(data, dist, given) {
  if (!is.data.frame(data)) data <- as.data.frame(data)
  nms <- setdiff(names(data), given)     # anything stated in `...` wins outright
  if (!length(nms))
    stop("admixr2: `data` adds no covariate that is not already given by name.",
         call. = FALSE)
  bad <- function(nm, ...) stop("admixr2: covariate '", nm, "' ", ...,
                                call. = FALSE)
  binary <- function(p) list(values = c(0, 1), probs = c(1 - p, p))

  specs <- list(); cont <- character(0)
  for (nm in nms) {
    v <- data[[nm]]
    if (is.factor(v) || is.character(v)) {
      lv <- sort(unique(as.character(v)))
      if (length(lv) != 2L)
        bad(nm, "has ", length(lv), " levels. A data frame supplies BINARY ",
            "covariates; give a multi-level one yourself as `", nm,
            " = list(values = ..., probs = ...)`.")
      # probability goes on the SECOND level in sort order, matching the
      # `SEX = c(male = 0.55)` convention: levels 0/1, the name on 1
      specs[[nm]] <- binary(mean(as.character(v) == lv[2L]))
      next
    }
    v <- as.numeric(v)
    if (!all(is.finite(v)))
      bad(nm, "has missing or non-finite values. A margin has to describe ",
          "every subject the study reports on -- drop the column, or state ",
          "it yourself.")
    if (all(unique(v) %in% c(0, 1))) { specs[[nm]] <- binary(mean(v)); next }
    # A CONSTANT COLUMN IS NOT A DISTRIBUTION, it is a pinned value. Left to
    # run, mean/sd gives sd = 0 -> sdlog = 0, which no margin accepts, and the
    # refusal surfaces much later inside datagen() as "not a supported
    # distribution" -- far from the column that caused it.
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

  # THE CORRELATION IS ON THE LATENT SCALE. For a lognormal margin the copula
  # runs over log(x), so the number wanted is cor(log(WT), log(CRCL)) and NOT
  # cor(WT, CRCL) -- close enough to look right, wrong enough to matter, and
  # the reason this helper earns its place.
  tf  <- if (identical(dist, "lnorm")) log else identity
  rho <- numeric(0)
  if (length(cont) > 1L)
    for (p in asplit(utils::combn(cont, 2L), 2L)) {
      r <- stats::cor(tf(data[[p[1L]]]), tf(data[[p[2L]]]))
      if (is.finite(r) && abs(r) > 1e-8) rho[paste(p, collapse = ".")] <- r
    }

  # A discrete margin correlated with a continuous one makes each level a
  # truncation of the latent normal rather than a point, which admPopulation()
  # refuses -- so those pairs are dropped, and said out loud where the data
  # actually show one.
  #
  # THE BAR IS SAMPLING NOISE, not a fixed number: under independence a sample
  # correlation has SE ~ 1/sqrt(n), so a flat 0.1 fires on a perfectly
  # independent cohort of 60 more often than not. Three SE, because this
  # message proposes changing the design; 0.1 stays as the floor so a huge
  # cohort does not report an association too small to matter.
  disc <- setdiff(names(specs), cont)
  if (length(disc) && length(cont)) {
    lim  <- max(0.1, 3 / sqrt(nrow(data)))
    hit  <- function(d, c) {
      r <- suppressWarnings(stats::cor(as.numeric(data[[d]]), data[[c]]))
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
  # A BASELINE TABLE READ OFF A COHORT, rather than transcribed. The three
  # things this absorbs are the three that were written out by hand at every
  # source: the lognormal mean/sd, the sex proportion, and -- the one with a
  # wrong answer rather than a tedious one -- the correlation, which is on the
  # LATENT scale and so runs over log(WT) and log(CRCL), not WT and CRCL.
  if (!is.null(data)) {
    d <- .admPopFromData(data, dist, names(a))
    a <- c(a, d$specs)
    # `cor` the user gave WINS, pair by pair, so a stated correlation is never
    # silently replaced by the sample one.
    # A single continuous margin leaves NO pair, and a zero-length `cor` is
    # not an empty correlation table -- it has no names, so the pair check
    # below rejects it. Absent means NULL.
    if (is.null(cor)) { if (length(d$cor)) cor <- d$cor }
    else if (!is.matrix(cor))
      cor <- c(cor, d$cor[setdiff(names(d$cor), names(cor))])
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
  # PAIRWISE correlations, so a partial table needs no identity padding. The
  # old idiom -- rbind/cbind an identity row on for every independent covariate
  # -- is easy to get wrong and impossible to read.
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
        # a covariate name may itself contain a dot, so match against the
        # declared names rather than assuming the split is two-way
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
#'   one subgroup, e.g. `at = list(SEX = 1)`.
#' @param by Covariate the paper reports results SEPARATELY by, e.g.
#'   `by = "SEX"`. Expands into one study per level.
#' @param stratify Band the source into strata, so a covariate it fitted
#'   contributes a CONTRAST rather than one pooled number. `TRUE` bands over
#'   every covariate the source's own model ESTIMATED a coefficient for, and
#'   leaves the rest marginalised --- derived from the model, so nothing has to
#'   be restated. A covariate the model merely READS is not banded: weight at a
#'   fixed allometric exponent carries no fitted effect to recover, and banding
#'   on it would buy strata and no evidence. Name a covariate explicitly to
#'   override that judgement. Measured over 720 replicates: coverage 0.933
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
  # Recorded BEFORE `sd` is overwritten from `sem` and before `V` is derived
  # from either, because the default below turns on which currency the study
  # was actually written in -- and a study given both `V` and `sd` is using the
  # `V`.
  .from_spread <- is.null(V) && (!is.null(sd) || !is.null(sem))
  # A study has no name until admStudies() gives it one, so an unlabelled
  # study used to report itself as `study 'study'`. The symbol the model was
  # passed as is the one thing on hand that the user will recognise -- used for
  # MESSAGES ONLY, never stored, so it cannot leak into the study names
  # admStudies() derives from its own argument expressions.
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
  # A COHORT KNOWS ITS OWN SIZE. When the population is handed over as the
  # data frame of enrolled subjects, `n` is its row count and asking for it
  # again is asking the user to restate what they just supplied. print() shows
  # the number that was taken, so a subset handed over by mistake is visible
  # before anything is fitted.
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
    # The paper's numbers go INTO the model, so nothing downstream has to carry
    # them separately -- and no global variable is involved. This is what the
    # `CLp <<- p$CL` idiom was working around.
    if (length(est)) {
      d <- ui$iniDf
      d$est[match(names(est), d$name)] <- unname(est)
      ui$iniDf <- d
    }
  } else {
    if (!is.null(sem)) {
      if (!is.null(sd)) bad("has both `sd` and `sem`; give one.")
      # V is a PER-SUBJECT covariance, so a standard error of the mean has to be
      # scaled back up by n. Reading a SEM as an SD understates the spread by
      # sqrt(n), which is the commonest error in digitising a published figure.
      sd <- as.numeric(sem) * sqrt(n)
    }
    if (is.null(V) && is.null(sd))
      bad("digitised data needs a spread: `sd` per timepoint, `sem`, or a ",
          "full `V`.")
    if (is.null(V)) V <- as.numeric(sd)^2
    if (length(E) != length(times))
      bad("`E` has ", length(E), " values but `times` has ", length(times), ".")
  }
  # WHICH DENOMINATOR THE SUPPLIED V USES -- and the currency already says.
  #
  # The likelihood scores V as the sample covariance of n iid subjects, which
  # is the ML (divide by n) one. A PUBLISHED spread is the unbiased (n - 1)
  # one: that is what every paper and every stats package reports. So a study
  # written as `sd =` or `sem =` is in the n - 1 currency by construction, and
  # defaulting it to "ml" silently scores a V that is n/(n-1) too large --
  # 1.7% at n = 60, more as n falls, and it inflates the variance parameters
  # rather than the estimates, which is the half of the fit nobody eyeballs.
  #
  # `V =` is left at "ml": handing over a covariance MATRIX is the deliberate
  # act of someone who computed it, and the documentation tells them to use the
  # ML denominator. Guessing "unbiased" there would corrupt a correct V.
  #
  # A MODEL SOURCE HAS NO CHOICE. Its moments come from datagen(), which is ML
  # by construction, so "unbiased" is not a convention to declare but a claim
  # that would shrink a V that is already right -- refused rather than honoured.
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
  # A DATA FRAME IS A POPULATION, so accept one directly.
  #
  # EVERY column is a covariate except the ones that structurally cannot be.
  # Restricting to what the SOURCE model reads was tried and is WRONG: a
  # covariate the source never fitted is exactly the case this package exists
  # for -- three trials declaring CRCL that no published model conditions on is
  # what identifies a renal effect at all -- so dropping it would delete the
  # evidence silently and leave a fit that converges. The filter is therefore a
  # denylist of names rxode2/nlmixr2 reserve for the data itself.
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
  # PRINTED, not messaged. A resolved default the user never sees is the thing
  # this package keeps getting bitten by -- but one message per study is noise
  # in a ten-study meta-analysis, so it goes where they are already looking.
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
  # the argument EXPRESSION is the natural label -- admStudies(smith2019, ...)
  # should not need the name typed twice.
  #
  # ONLY when `...` IS the studies. The single-list branch above replaced `a`
  # with the list's contents, so the expressions in `...` describe the list,
  # not its elements -- `auto` then has length 1 against k studies and every
  # element after the first indexes out of bounds to NA_character_.
  #
  # And NA is not caught by the `!nzchar()` fallbacks, because
  # nzchar(NA_character_) is TRUE. So the name stayed NA, and
  # .admMaterialise() dropped that study silently -- studies[[NA_character_]]
  # is NULL, not an error, so a two-study fit quietly became a one-study fit.
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

  # THE PRE-FLIGHT. Whether a covariate is IDENTIFIED is a property of the
  # study SET, not of any one study, so this is the only place it can be said
  # -- and it has to be said BEFORE the fit, because the failure does not look
  # like one afterwards. A covariate marginalised identically everywhere enters
  # only through the mixture it induces, which is what a random effect on the
  # same parameter does; the profile is flat (0.019 units across its whole
  # range, measured) and a deterministic optimiser stops in the same place
  # every run, reading as converged at a wrong number.
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
          # only suggest banding where some source could actually be banded --
          # .admExpandStrata() refuses a covariate whose coefficient every
          # source asserted, so naming one here would send the reader into an
          # error rather than a fix
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

# Turn specs into generated studies, once, inside the fit.
#
# LAZY ON PURPOSE. Building a study solves nothing, so a transcription error
# surfaces on print() rather than after a long generate, and a specification can
# be assembled, inspected and edited freely. Everything a spec needs to generate
# itself is already on it, so this needs no arguments beyond the list.
#
# Called from .admDriverUnits(), which is the one place all four estimators
# pass through -- a step wired into three drivers of four is a silent
# divergence, not an error.
.admMaterialise <- function(studies) {
  if (inherits(studies, "admStudies")) studies <- unclass(studies)
  if (!is.list(studies)) return(studies)
  spec <- vapply(studies, inherits, logical(1), "admStudy")
  if (!any(spec)) return(studies)
  out <- list()
  for (nm in names(studies)) {
    s <- studies[[nm]]
    if (!inherits(s, "admStudy")) { out[[nm]] <- s; next }
    ev <- s$ev %||% rxode2::et(amt = s$dose)
    if (is.null(s$ui)) {
      # digitised: it IS the data, nothing to generate
      g <- list(E = as.numeric(s$E), V = s$V, n = s$n, times = s$times, ev = ev)
      if (!is.null(s[["population"]])) g[["cov_dist"]] <- s[["population"]]
      if (!is.null(s[["at"]]))         g[["cov"]]      <- s[["at"]]
      out[[nm]] <- g
      next
    }
    # `[[ ]]` THROUGHOUT on the datagen spec. `sp` carries `cov_dist`, so with
    # no `cov` element `sp$cov` PARTIAL-MATCHES to it and hands back the whole
    # population spec -- which then got concatenated with the pinned level
    # below, producing a `cov` with two SEX entries and a fit that returned no
    # standard errors at all. Third instance of this trap in one day; the house
    # rule exists for a reason.
    sp <- list(times = s$times, ev = ev, n = s$n)
    if (!is.null(s[["population"]])) sp[["cov_dist"]] <- s[["population"]]
    if (!is.null(s[["at"]]))         sp[["cov"]]      <- s[["at"]]
    if (!is.null(s$stratify)) {
      sp$stratify <- s$stratify
      if (!is.null(s$strata_nodes)) sp$strata_nodes <- s$strata_nodes
      # `range` names the covariate it belongs to. With `stratify = TRUE` the
      # banded set is derived, not typed, so setNames(list(range), TRUE) built
      # a cov_range for a covariate called "TRUE" -- unusable with the very
      # spelling the documentation recommends. Resolve the band the same way
      # the pre-flight does, and ask for a named list when it is ambiguous.
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
    # `by`: the paper reported results separately per level, which is one
    # ordinary study per level with that covariate PINNED -- not a stratified
    # one, because the source really did report each subgroup.
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
        # every element: `by` combined with `stratify` expands each level into
        # J bands, and taking only the first silently dropped J-1 of them
        for (kk in names(g)) out[[kk]] <- g[[kk]]
      }
      next
    }
    g <- datagen(stats::setNames(list(sp), nm), model = s$ui,
                 control = datagenControl(method = "gh"))
    for (k in names(g)) out[[k]] <- g[[k]]
  }
  out
}

# =============================================================================
# The standard error that would otherwise be a silent wrong answer
# =============================================================================

# Does any study contribute as a published MODEL rather than as a sample?
#
# Three shapes reach here: an admStudy() spec, a raw datagen spec written by
# hand, and an already-generated study. `[[ ]]` throughout -- `$cov` partial-
# matches `cov_dist` on a datagen spec, which is the trap .admMaterialise()
# records.
# The covariates `stratify` actually bands, for one study.
#
# `TRUE` bands what the source ESTIMATED, not what it merely reads: weight at a
# fixed allometric exponent carries no fitted effect to recover. Shared by the
# pre-flight print and the expansion below, so the table cannot promise a
# banding the fit will not perform.
.admStudyBandNames <- function(s) {
  st <- s[["stratify"]]
  if (is.null(st) || identical(st, FALSE)) return(character(0))
  if (!isTRUE(st)) return(as.character(st))
  Filter(function(cv)
           length(.admCovCoefThetas(s$ui, cv, s[["population"]])) > 0L,
         intersect(.admCovSpecNames(s[["population"]]),
                   tryCatch(s$ui$allCovs, error = function(e) character(0))))
}
