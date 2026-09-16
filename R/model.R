# Cache-key material shared by the simulation and sensitivity models.
#
# Names, fix flags, error types, and the fixed thetas' VALUES: a fix()ed
# theta travels to the solve as DATA, not as an optimizer parameter, so
# `tka <- fix(0.5)` and `tka <- fix(0.9)` must key differently -- otherwise
# a parallel worker solves every restart at the other fit's fixed value
# (rxTempDir() persists across sessions).
#
# Only FIXED rows: folding in `est` would invalidate the cache on every
# starting-value change, which is optimizer state, not model text.
.admIniKey <- function(ui) {
  tryCatch({
    .i <- ui$iniDf
    .fx <- !is.na(.i$fix) & .i$fix
    paste(paste(.i$name, collapse = "|"),
          paste(as.integer(.i$fix), collapse = "|"),
          paste(.i$err, collapse = "|"),
          # format(NULL) is the literal string "NULL"; keep the empty case empty
          # so a ui without an iniDf gets a clean key rather than a stray token.
          if (any(.fx)) paste(format(.i$est[.fx], digits = 17), collapse = "|") else "",
          sep = "//")
  }, error = function(e) "")
}

# Identity of the code that EMITS a cached model, replacing a hand-maintained
# schema tag that was forgotten once and served a stale entry.
#
# The version ALONE is not enough: `Version:` moves only at release, so a
# whole dev cycle shares one key, and editing .g2, the direction set, or
# the f2 naming would then hit a stale compiled model and produce a
# finite, plausible, silently wrong gradient. So the key also digests the
# BODIES of the functions that decide what gets emitted.
.admPkgKey <- function() {
  .ver <- tryCatch(as.character(utils::packageVersion("admixr2")),
                   error = function(e) "dev")
  # deparse(), not the closure: the package namespace digests differently
  # between load_all() and an installed build. NAME LIST, not just the two
  # entry points, since a helper like .admJumpCovers() also shapes the
  # payload; digested BY VALUE, per-name tryCatch so one unresolvable name
  # doesn't disable invalidation entirely. LOCAL, not package-level:
  # .admDaemonRestart() patches a dev body into the stale installed
  # namespace via assignInNamespace(), which can REPLACE a binding but not
  # ADD one. RULE: belongs here if changing its BODY changes the cached
  # payload.
  .emitters <- c(
    ".admBuildThetaSens",   # emits the direction set, the chains and the f2 block
    ".admLoadSensModel",    # assembles the cached list and its fallbacks
    ".admSensNameMaps",     # THETA[k]/ETA[i] maps that fill the emitted columns
    ".admSensFromInner",    # builds the whole type = "inner" payload
    ".admLinCmtToOde",      # emits the promoted ODE an order-2 linCmt differentiates
    ".admRxode2",           # artifact name + wd + eventSens handed to the compiler
    ".admModName",          # ... and the name itself
    ".admJumpCovers",       # decides whether a model is cached at all
    ".admToRx",             # rxFromSE wrapper -- emits EVERY line of the model text
    ".admEndpointVar",      # which variable an endpoint's sensitivities follow
    ".admCountSpec",        # ... and the two it delegates to for a discrete or
    ".admBetaSpec",         #     bounded endpoint, which shape the same pred_expr
    ".admBetaPair",         # the beta shape pair the emitted pred_expr is built from
    ".admDistArgs",         #     ... and the argument parser those two locate it with
    ".admIsLinCmtMod")      # fills the cached is_lincmt flag the consumers branch on
  .src <- tryCatch(
    digest::digest(c(
      lapply(.emitters,
             function(.n) tryCatch(deparse(body(get(.n))), error = function(e) .n)),
      list(.admDoseModRe, .ADM_TBS_YJ))),
    error = function(e) "NA")
  paste(.ver, .src, sep = "/")
}


# Compile a generated model under a stable, role-tagged artifact name.
#
# Ported from nlmixr2est's
# .nlmixr2estRxode2()/.nlmixr2estModName()/.nlmixr2estModDir()
# (main@fef5be69, post-6.2.0), with admixr2's own roles.
#
# rxode2 names a model's .c/.so from the parsed text ALONE for an ANONYMOUS
# model -- it can't see inputs like the event-sensitivity code injected
# afterward. Two builds of one text differing only there land on a single
# .so; the later build wins, and since entry points resolve BY NAME an
# object bound to the earlier build silently executes the replacement
# (nlmixr2/rxode2#1171). admixr2 hits this hardest: .admSensFromInner()
# recompiles nlmixr2est's OWN inner model text with a different eventSens,
# same parsed md5 -- and rxTempDir() persists, so the collision survives restarts.
.admRxode2 <- function(model, role, ...) {
  .nm <- .admModName(model, role, ...)
  .wd <- .admModDir()
  # The fallback still builds in OUR directory and SYNTHESISES a name rather
  # than omitting one -- rxode2's own naming+directory would overwrite as
  # before, and `rxode2(model, wd=)` errors without a modName. The
  # synthesised name folds in `...` (carries eventSens, the #1171 mechanism)
  # and is `admMod_*`, not `admSens*` -- see .admRxLoadAll()'s PATH check.
  if (is.null(.nm))
    .nm <- paste0("admMod_", digest::digest(list(model, role, list(...))))
  rxode2::rxode2(model, modName = .nm, wd = .wd, ...)
}

# Stable artifact name: role + everything that changes the emitted code.
# `eventSens` is folded in since a "jump" and an "fd" build of one model text
# would otherwise share an artifact. The parsed md5 stays in the name so two
# genuinely different models still never collide.
.admModName <- function(model, role, ...) {
  if (is.null(role) || !nzchar(role)) return(NULL)
  .md5 <- tryCatch(rxode2::rxModelVars(model)$md5[["parsed_md5"]], error = function(e) NULL)
  if (is.null(.md5) || !nzchar(.md5)) return(NULL)
  .dots <- list(...)
  # eventSens can arrive as an un-evaluated match.arg default, i.e. c("jump","fd"):
  # take the first element, as match.arg would, so the name stays length 1. A
  # zero-length value must fall back to "" rather than produce character(0), which
  # would sail past an is.null() check and reach rxode2 as an empty modName.
  .es <- .dots$eventSens
  .es <- if (is.null(.es) || length(.es) == 0L) "" else gsub("\\W", "", as.character(.es)[1L])
  if (is.na(.es)) .es <- ""
  # "_" separates the parts: without it role "rxA" + es "bc" and role "rxAb" + es
  # "c" would produce one name. The md5 is last and fixed-width.
  .nm <- paste0(role, "_", .es, "_", .md5)
  if (length(.nm) != 1L || is.na(.nm) || !nzchar(.nm)) return(NULL)
  .nm
}

# Where generated models are built.
#
# NOT getwd(): scatters <modName>.d directories through the user's working
# directory. NOT rxTempDir(): upstream measured a generated model's .so
# replaced mid-run by a build emitting different event-sensitivity code, and
# for admixr2 it's also a PERSISTENT user cache. R's session temp directory
# is neither -- cost is that artifacts aren't shared across sessions, worth
# it while nlmixr2/rxode2#1171 is open.
#
# The .rds caches still live in rxTempDir(), so a cross-session hit may
# reference a DLL this session no longer has; rxLoad() doesn't reliably
# error on that, it quietly binds to whatever shares the entry-point name
# and solves to garbage. .admRxLoadAll() guards this explicitly. admixr2
# crosses artifact and cache lifetimes for the SENSITIVITY model only
# (session-local artifact, persistent adm-sens-*.rds cache); the
# simulation model stays on upstream's disk pairing -- chosen over a
# session-only cache (breaks parallel grad = "sens") for the cold-start win.
#
# Disappears entirely once nlmixr2/rxode2#1171 is fixed.
.admModDir <- function() {
  .d <- file.path(tempdir(), "admixr2Sens")
  if (!dir.exists(.d)) {
    dir.create(.d, recursive = TRUE, showWarnings = FALSE)
    # dir.create() is quiet about failure here (another process may have won the
    # race, which is fine), so confirm rather than hand back a path that does not
    # exist -- rxode2 would then fail deep in the compile with a confusing error.
    if (!dir.exists(.d)) return(tempdir())
  }
  .d
}

# Re-load every compiled model inside a cached object.
#
# nlmixr2est's own load step, verbatim (rxUiGet.foceiModel). A deserialised
# rxode2 object carries a dead pointer until rxLoad() re-attaches its shared
# library, and the CONTAINER is what's cached -- so every rxode2-classed
# element must be re-loaded, not just the first one a caller reads.
#
# TRUE if everything loaded; FALSE treats it as a stale cache entry.
# Canonical spelling of a path, for comparing two of them.
#
# normalizePath() resolves Windows 8.3 short components against the long
# form tempdir() reports, and warns (suppressed) for a nonexistent path.
# Case-folded on Windows ONLY, since only that filesystem is
# case-insensitive -- folding elsewhere would false-match two real directories.
.admNormPath <- function(p) {
  .p <- tryCatch(normalizePath(p, winslash = "/", mustWork = FALSE),
                 error = function(e) p, warning = function(w) p)
  if (.Platform$OS.type == "windows") tolower(.p) else .p
}

.admSameDir <- function(a, b) identical(.admNormPath(a), .admNormPath(b))

# Does `path` lie inside THIS session's temporary directory?
#
# More than one build directory counts as "this session" (admixr2's
# .admModDir() and nlmixr2est 7.x's <tempdir>/nlmixr2estSens), so this
# prefix-tests against the tempdir rather than checking .admModDir() alone
# -- which once rejected the second, recompiling the model every fit.
#
# BOTH SPELLINGS OF tempdir() are needed: on macOS tempdir() is a symlink
# that only resolves for a path that EXISTS, so a not-yet-built artifact
# path under it fails a prefix test against the resolved root.
.admUnderTemp <- function(path) {
  # Collapse repeated separators before comparing. R's tempdir() on macOS is
  # commonly ".../T//RtmpXXXX", and an unresolved path keeps that doubled slash
  # while a resolved one loses it. A LEADING "//" is preserved, being a UNC root.
  .sq <- function(x) gsub("(?<=.)/{2,}", "/", x, perl = TRUE)
  .p     <- .sq(.admNormPath(path))
  # The raw spelling still has to match .admNormPath()'s conventions -- forward
  # slashes, and case-folded on Windows -- or it is not a candidate at all, just
  # a string that can never match.
  .raw   <- gsub("\\\\", "/", tempdir())
  if (.Platform$OS.type == "windows") .raw <- tolower(.raw)
  .cands <- unique(.sq(c(.admNormPath(tempdir()), .raw)))
  any(vapply(.cands, function(.t) {
    .t <- sub("/+$", "", .t)
    identical(.p, .t) || startsWith(.p, paste0(.t, "/"))
  }, logical(1)))
}

.admRxLoadAll <- function(x) {
  .one <- function(e) {
    if (!inherits(e, "rxode2")) return(TRUE)
    # rxLoad() does NOT error on a cached model whose shared object has
    # gone -- it silently recompiles (~2.9 s), and N daemons reading one
    # .rds all recompile concurrently into the SAME path. So the DLL is
    # checked EXPLICITLY first.
    .dll <- tryCatch(rxode2::rxDll(e), error = function(err) NA_character_)
    if (is.na(.dll) || !nzchar(.dll) || !file.exists(.dll)) return(FALSE)
    # file.exists() alone isn't enough: R only removes its temp directory on
    # a CLEAN exit, so a killed session's build directory lingers and a
    # concurrent LIVE session's is equally readable -- so isolation must be
    # enforced, not assumed. Discriminated on the PATH (a basename test
    # misses .admRxode2()'s anonymous fallback and nlmixr2est's inner
    # model), tested against the SESSION TEMPDIR via .admUnderTemp(), after
    # NORMALISING: rxDll() hands back Windows 8.3 short form, so a literal
    # grepl("admixr2Sens", .) would silently never match.
    .dllN <- tryCatch(normalizePath(.dll, winslash = "/", mustWork = FALSE),
                      error = function(err) .dll)
    if (grepl("(admixr2Sens|nlmixr2estSens)", .dllN) &&
        !.admUnderTemp(.dll)) return(FALSE)
    tryCatch({ rxode2::rxLoad(e); TRUE }, error = function(err) FALSE)
  }
  if (inherits(x, "rxode2")) return(.one(x))
  if (!is.list(x)) return(TRUE)
  all(vapply(x, .one, logical(1)))
}

# Disk-cache path for the compiled simulation model.
#
# The key covers the model({}) block AND .admIniKey(ui) -- names, fix()
# flags, error types, and the VALUES of the fixed thetas. The lstExpr digest
# alone isn't enough: `tka <- fix(0.5)` and `tka <- fix(0.9)` digest
# identically, and since a fixed theta never reaches the optimizer, the
# second fit would silently solve at the first's baked-in value. Same
# collision .admIniKey() already closes for the SENSITIVITY cache.
#
# Split out as its own function: a parallel worker has no `ui` to recompute
# it, so the parent stores this path on pinfo for .admWorkerLoadModels().
.admModelCacheFile <- function(ui) {
  # The rxode2 VERSION is part of the key, as it already is for the sens cache.
  # rxTempDir() is a persistent user cache and nothing in admixr2 or rxode2 sweeps
  # it. rxode2 itself treats its own binary as part of a model's identity, so after
  # an rxode2-only upgrade a fresh compile takes a NEW artifact name while this
  # entry still points at the old one. In practice nlmixr2est's
  # .resetCacheIfNeeded() calls rxClean() on ITS version change, which is why this
  # has not bitten -- but that is someone else's hook doing our invalidation.
  .rx_ver <- tryCatch(as.character(utils::packageVersion("rxode2")),
                      error = function(e) "NA")
  file.path(rxode2::rxTempDir(),
            paste0("adm-sim-",
                   digest::digest(list(ui$lstExpr, .admIniKey(ui), .rx_ver)), ".rds"))
}

# Load (or compile + cache) the rxode2 simulation model.
# Compiled DLL is cached to disk with saveRDS(), keyed by model digest -- the
# same shape nlmixr2est uses for its own compiled models (rxTempDir(), a
# digest-named .rds, saveRDS/readRDS); see .admRxLoadAll for the load step.
.admLoadModel <- function(ui) {
  # Accessing $simulationModel (below) caches the compiled model in
  # ui$meta$.simModelBase as a side effect -- a live, self-referential rxode2
  # object that breaks nlmixr2's ui-cloning during fit assembly. Drop it (and any
  # sibling artifacts) on every exit so the ui stays in the canonical state
  # nlmixr2 expects; see .admDropSimModelMeta() for the full rationale.
  on.exit(.admDropSimModelMeta(ui), add = TRUE)
  .cacheFile <- .admModelCacheFile(ui)
  .model_key <- sub("\\.rds$", "", basename(.cacheFile))
  # In-session registry (.adm_model_env): a disk cache HIT still costs a
  # readRDS + dyn.load and a finalizer-bearing object. Gated on the disk
  # cache STILL EXISTING, so a cleared cache still forces a recompile
  # instead of being served from the memo.
  .memo <- get0(.model_key, envir = .adm_model_env, inherits = FALSE)
  if (!is.null(.memo) && file.exists(.cacheFile) && .admRxLoadAll(.memo))
    return(.memo)
  if (file.exists(.cacheFile)) {
    mod <- tryCatch(readRDS(.cacheFile), error = function(e) NULL)
    # inherits() FIRST, then load. .admRxLoadAll mirrors nlmixr2est's load step
    # exactly, and that step is a no-op on anything not rxode2-classed -- so on its
    # own it reports TRUE for a file whose content is not a model at all. Asserting
    # the payload's SHAPE here means the file is deleted and recompiled rather than
    # handed back and memoised, in which case every fit in the session repeats it.
    load_ok <- inherits(mod, "rxode2") && .admRxLoadAll(mod)
    if (load_ok) {
      return(.admCacheAssign(.model_key, mod, .adm_model_env))
    }
    # suppressWarnings: file.remove() signals a WARNING on failure, not an error, so
    # tryCatch(error=) alone lets it through -- a stale entry whose DLL is still
    # loaded cannot be unlinked on Windows. Recompiling is the correct recovery
    # either way; the warning is noise.
    tryCatch(suppressWarnings(file.remove(.cacheFile)), error = function(e) NULL)
  }
  # rxode2 compilation calls setwd() internally -- save/restore to avoid
  # "cannot change working directory" error on first compile (Windows).
  #
  # STAYS rxTempDir(): setwd()-ing here instead breaks multi-endpoint
  # models, since `rxode2(ui)$simulationModel` compiles companion models
  # that resolve against the working directory -- moving it made every
  # multi-output fit return an all-NA structural gradient. The build-directory
  # change applies only where admixr2 passes `wd =` explicitly (.admRxode2()).
  .old_wd <- tryCatch(getwd(), error = function(e) NULL)
  on.exit(if (!is.null(.old_wd)) setwd(.old_wd), add = TRUE)
  setwd(rxode2::rxTempDir())
  mod <- rxode2::rxode2(ui)$simulationModel
  .admCacheWrite(mod, .cacheFile, "simulation model")
  rxode2::rxLoad(mod)
  .admCacheAssign(.model_key, mod, .adm_model_env)
}

# Remove transient rxode2 model objects that $simulationModel /
# $foceiModel leave behind in ui$meta.
#
# nlmixr2's deep-clone (.cloneEnv()) recurses into every environment-valued
# member with no cycle detection, and rxode2's compiled model objects hold a
# back-reference to the global .rxModels registry, so cloning one loops
# forever ("node stack overflow"). admixr2 populates ui$meta with these
# (simulating via $simulationModel); safe to drop since rxode2 regenerates
# them lazily if needed.
.admDropSimModelMeta <- function(ui) {
  .meta <- ui$meta
  if (!is.environment(.meta)) return(invisible())
  for (.nm in ls(.meta, all.names = TRUE)) {
    .v <- get(.nm, envir = .meta, inherits = FALSE)
    if (is.environment(.v) && inherits(.v, "rxode2"))
      rm(list = .nm, envir = .meta)
  }
  invisible()
}

# Structural thetas with no usable mu-referenced eta ("unpaired") -- these
# get their OWN sensitivity direction (THETA_j_) in the sens model rather
# than reusing an eta's gradient column.
#
# Uses .admMuRefPairs() -- the SAME map pinfo$struct_has_eta is built from --
# so the estimators' theta-column routing and the sens model's built columns
# cannot drift apart. Includes the shared-eta guard: a theta whose eta
# appears in another parameter is unpaired.
.admUnpairedThetas <- function(ui) {
  ini <- tryCatch(ui$iniDf, error = function(e) NULL)
  if (is.null(ini)) return(character(0))
  struct <- ini[is.na(ini$neta1) & is.na(ini$err) & !ini$fix, , drop = FALSE]
  mrd    <- .admMuRefPairs(ui)
  paired <- if (!is.null(mrd)) as.character(mrd$theta) else character(0)
  setdiff(struct$name, paired)
}

# rxode2::rxFromSE() substitutes its argument, so it MUST be called through a
# wrapper -- calling it directly on an inline expression emits the literal call
# text instead of the model code. (nlmixr2est's aug builder has the same wrapper.)
.admToRx <- function(l) rxode2::rxFromSE(l)

# Promote a linCmt() model to the equivalent explicit ODE system, for
# SECOND-order sensitivities: linCmt() has no second derivative
# (nlmixr2est refuses linCmt outright), but ordinary state sensitivities do.
#
# `rxode2::linToOde()` names the prediction `rxLinCmt`, an `rx`-prefixed lhs
# that is RESERVED and dies with "syntax errors" in ui$loadPruneSens --
# renaming it is the whole fix. Reproduces the analytic prediction to
# 1.8e-08 relative on a 1-cmt oral model.
#
# NULL on failure -- caller falls back to first-order. Only called for
# order-2; admc/adgh stay on the SOLVED form.
.admLinCmtToOde <- function(ui) {
  tryCatch({
    .u   <- rxode2::linToOde(ui)
    .txt <- paste(deparse(.u$fun), collapse = "\n")
    if (!grepl("rxLinCmt", .txt, fixed = TRUE)) return(.u)
    # Pick a name that is neither reserved (no `rx` prefix) nor already in the
    # model. A collision would silently merge two different quantities.
    .nm <- "admLinCmtOut"
    while (grepl(.nm, .txt, fixed = TRUE)) .nm <- paste0(.nm, "X")
    .txt <- gsub("rxLinCmt", .nm, .txt, fixed = TRUE)
    suppressMessages(rxode2::rxode2(eval(parse(text = .txt))))
  }, error = function(e) NULL)
}

# The dosing-modifier variables rxode2 emits into the pruned sens env
# as rx_<mod>_<state>_ (f/lag/rate/dur; lag() stored as alag()). ONE
# regex used everywhere a dose modifier is found or extracted, so a
# rxode2 naming change is a single edit. Group 1 = modifier, 2 = state.
.admDoseModRe <- "^rx_(f|lag|alag|rate|dur)_(.+)_$"

# Can this rxode2 build differentiate every dosing modifier one of OUR
# directions actually feeds?
#
# A direction entering f()/lag()/rate()/dur() only has a sensitivity if
# rxode2 attaches analytic variational jumps at dose times (eventSens =
# "jump"); otherwise the column is silently ZERO. FEATURE-DETECT via
# eventSensInfo$derivs rather than version-compare.
#
# Tested per DIRECTION, not "the model has a lag()": nlmixr2est's inner
# model has no theta directions, so `alag(depot) = exp(tlag)`'s derivs$lag
# is legitimately empty.
#
# FALSE -> caller falls back to a finite-difference gradient rather than a
# silently zero component.
.admJumpCovers <- function(mod, s, dirs) {
  vars <- grep(.admDoseModRe, ls(envir = s, all.names = TRUE), value = TRUE)
  if (length(vars) == 0L || length(dirs) == 0L) return(TRUE)

  need <- character(0)
  for (v in vars) {
    ex <- tryCatch(get(v, envir = s), error = function(e) NULL)
    if (is.null(ex)) next
    depends <- any(vapply(dirs, function(p)
      !identical(tryCatch(.admToRx(symengine::D(ex, symengine::S(p))),
                          error = function(e) "0"), "0"),
      logical(1)))
    if (!depends) next
    key <- sub(.admDoseModRe, "\\1", v)
    need <- c(need, if (identical(key, "alag")) "lag" else key)
  }
  need <- unique(need)
  if (length(need) == 0L) return(TRUE)

  info <- tryCatch(mod$eventSensInfo, error = function(e) NULL)
  if (is.null(info) || !identical(info$mode, "jump")) return(FALSE)
  d <- info$derivs
  if (!is.list(d)) return(FALSE)
  all(vapply(need, function(m) is.data.frame(d[[m]]) && nrow(d[[m]]) > 0L, logical(1)))
}

# Build the sensitivity model over an explicit DIRECTION SET:
#
#   dirs = ETA_1_ .. ETA_n_        (one per random effect)
#        + THETA_j_                (one per UNPAIRED structural theta)
#
# A mu-referenced theta reuses its eta's column (d(pred)/d(theta) ==
# d(pred)/d(eta)) rather than getting its own direction. Sigmas get none.
# Same scheme as nlmixr2est's fast-focei (.foceiAnalyticDirections).
#
# `order = 2L` additionally emits the CROSS second-order block
#   rx_f2_<eta_i>_<dir_b> = d2(pred)/(d eta_i d dir_b)
# which only adfo needs (its objective needs dJ/d(theta) = d2f/(d eta d
# theta)); admc/adgh stay at order 1. Deliberately ASYMMETRIC -- eta
# directions x ALL directions, not nlmixr2est's full symmetric triangle --
# since adfo differentiates only ONCE per theta and skips the chains that
# dominate FOCEI's order-2 compile time.
#
# An ODE model gets rxode2::.rxSens() variational (state-sensitivity)
# compartments; a linCmt model has no states to augment, so D(pred, dir)
# resolves directly through linCmtB. Compiled with eventSens = "jump" so a
# dosing modifier (f/lag/rate/dur) gets analytic variational jumps at dose
# times, else its sensitivity is silently ZERO.
#
# Returns list(mod, dirs, sens_cols, theta_sens_cols) or NULL (caller falls
# back to nlmixr2est's inner model + FD). `pred_expr` defaults to
# `rx_pred_`, but a likelihood-form endpoint has the LOG-LIKELIHOOD there
# instead of the MEAN admixr2 needs, so such an endpoint passes its own
# derived expression.
.admBuildThetaSens <- function(ui, unpaired, pred_expr = NULL, order = 1L) {
  order <- as.integer(order)
  s <- tryCatch(ui$loadPruneSens, error = function(e) NULL)
  if (is.null(s)) return(NULL)
  st <- tryCatch(rxode2::rxStateOde(s), error = function(e) NULL)
  if (is.null(st)) return(NULL)

  ini      <- tryCatch(ui$iniDf, error = function(e) NULL)
  if (is.null(ini)) return(NULL)
  eta_rows <- ini[!is.na(ini$neta1) & ini$neta1 == ini$neta2 & !ini$fix, , drop = FALSE]
  eta_rows <- eta_rows[order(eta_rows$neta1), , drop = FALSE]
  th_rows  <- ini[!is.na(ini$ntheta), , drop = FALSE]

  # linCmt() carries no SECOND derivative, so an order-2 request promotes
  # the model to its ODE form first (NULL on failure retries at order 1).
  #
  # DETECTED WITH rxode2::testRxLinCmt(), not predDf$linCmt: on rxode2 5.1.4
  # that column is FALSE for a genuine `cp <- linCmt()` model, so the gate
  # never fired and EVERY order-2 request on a linCmt model silently dropped
  # to the FD pass it was meant to replace. A PROMOTED solved-form linCmt is
  # caught instead by the linCmtB text check below.
  if (order >= 2L) {
    .lin <- tryCatch(isTRUE(rxode2::testRxLinCmt(ui)), error = function(e) NULL)
    if (is.null(.lin))   # older rxode2 without the predicate
      .lin <- tryCatch(isTRUE(any(as.logical(ui$predDf$linCmt), na.rm = TRUE)),
                       error = function(e) FALSE)
    if (isTRUE(.lin)) {
      ui <- .admLinCmtToOde(ui)
      if (is.null(ui)) return(NULL)
      s  <- tryCatch(ui$loadPruneSens, error = function(e) NULL)
      if (is.null(s)) return(NULL)
      st <- tryCatch(rxode2::rxStateOde(s), error = function(e) NULL)
      if (is.null(st) || length(st) == 0L) return(NULL)
      ini <- tryCatch(ui$iniDf, error = function(e) NULL)
      if (is.null(ini)) return(NULL)
      # RE-DERIVE eta_rows/th_rows from the PROMOTED model's iniDf. These drive the
      # direction set, so taking them from the pre-promotion ui means that if
      # linToOde() yields ANY iniDf difference -- a renumbered ntheta, an added row,
      # a different eta ordering -- the emitted rx_f1_THETA_k_ differentiates a
      # different parameter than the caller thinks, and adfo descends a structural
      # gradient computed for the wrong theta with no error. linToOde() does
      # preserve the iniDf on the models measured here, so this was latent rather
      # than firing; deriving from the model actually being differentiated makes
      # that an outcome rather than an assumption.
      eta_rows <- ini[!is.na(ini$neta1) & ini$neta1 == ini$neta2 & !ini$fix, , drop = FALSE]
      eta_rows <- eta_rows[order(eta_rows$neta1), , drop = FALSE]
      th_rows  <- ini[!is.na(ini$ntheta), , drop = FALSE]
    }
  }

  eta_dirs <- paste0("ETA_", seq_len(nrow(eta_rows)), "_")
  # NB: paste0("THETA_", integer(0), "_") is "THETA__", not character(0) -- R
  # recycles the zero-length argument to "". Guard, or a model with no unpaired
  # theta gets a phantom direction.
  theta_dirs <- character(0)
  if (length(unpaired) > 0L) {
    theta_idx <- th_rows$ntheta[match(unpaired, th_rows$name)]
    if (anyNA(theta_idx)) return(NULL)
    theta_dirs <- paste0("THETA_", theta_idx, "_")
  }
  dirs <- c(eta_dirs, theta_dirs)
  if (length(dirs) == 0L) return(NULL)

  # matExp() / indLin(): rxStateOde() can return the states REVERSED (an indLin
  # state parses as compartment 1), so emitting the ODEs in that order would put
  # the dose in the wrong compartment. nlmixr2est fixes this with an internal
  # reorder (.rxMatExpStateOrder); rather than reimplement it, bail out and let
  # the caller fall back to nlmixr2est's inner model + FD -- correct, just slower.
  .mv <- tryCatch(rxode2::rxModelVars(s), error = function(e) NULL)
  if (!is.null(.mv) && is.list(.mv$indLin) && length(.mv$indLin) == 4L) return(NULL)

  res <- tryCatch({
    sens_lines <- character(0)
    sens2_lines <- character(0)
    if (length(st) > 0L) {
      rxode2::.rxJacobian(s, c(st, dirs))
      sens_lines <- rxode2::.rxSens(s, dirs)
      if (length(sens_lines) == 0L) return(NULL)
      # Cross block only: rows = eta directions, columns = every direction.
      # .rxSens assigns rx__sens_*_BY_*_BY_*__ into `s` as a side effect, so
      # this must run BEFORE the chains are built (.g2 below resolves them).
      #
      # Requested one eta ROW at a time, against only the directions at or
      # after it: the eta x eta half is SYMMETRIC but rxExpandSens2_ doesn't
      # know it, and the full rectangle would emit both d2/(d eta_1 d eta_2)
      # and its mirror as separate compartments -- 18 of 93 extra states on
      # a 3-state/4-eta/2-theta model for no information.
      if (order >= 2L) {
        for (.i in seq_along(eta_dirs)) {
          .cols <- c(eta_dirs[.i:length(eta_dirs)], theta_dirs)
          .l <- rxode2::.rxSens(s, eta_dirs[.i], .cols)
          if (length(.l) == 0L) return(NULL)
          sens2_lines <- c(sens2_lines, .l)
        }
      }
    }
    pred <- if (!is.null(pred_expr)) pred_expr else get("rx_pred_", envir = s)
    .Dn  <- function(e, v) symengine::D(e, symengine::S(v))
    # Variadic, matching the compartment naming rxExpandSens2_ emits:
    # one direction  -> rx__sens_<state>_BY_<p>__
    # two directions -> rx__sens_<state>_BY_<p>_BY_<q>__
    .sn1 <- function(j, ...)
      symengine::S(paste0("rx__sens_", j, "_BY_",
                          paste(c(...), collapse = "_BY_"), "__"))
    # linCmt: st is empty, so the state sum drops out and D(pred, dir) alone
    # resolves through the linCmtB derivative rules.
    .g1 <- function(ex, p) {
      e <- .Dn(ex, p)
      for (j in st) e <- e + .Dn(ex, j) * .sn1(j, p)
      e
    }
    # Second-order chain. Differentiating .g1(ex, q) w.r.t. p picks up three terms:
    # the direct partial, the first-order state paths of the already-chained
    # expression, and the second-order state sensitivities themselves. Same
    # construction as nlmixr2est's .g2. Getting any one term wrong yields a
    # plausible-but-wrong Jacobian derivative, so it is FD-checked in
    # test-integration-sens2.R rather than trusted.
    .g2 <- function(ex, p, q) {
      gq <- .g1(ex, q)
      e  <- .Dn(gq, p)
      for (k in st) e <- e + .Dn(gq, k) * .sn1(k, p)
      for (j in st) e <- e + .Dn(ex, j) * .sn1(j, p, q)
      e
    }

    base_ode <- if (length(st))
      vapply(st, function(x)
        paste0("d/dt(", x, ")=", .admToRx(get(paste0("rx__d_dt_", x, "__"), envir = s))),
        character(1)) else character(0)

    # dosing modifiers (bioavailability, lag, rate, duration) live in the pruned
    # env as rx_<mod>_<state>_ and are NOT part of rx__d_dt_*; rxode2 stores lag()
    # as alag().
    dos_vars <- grep(.admDoseModRe, ls(envir = s, all.names = TRUE), value = TRUE)
    dose <- vapply(dos_vars, function(v) {
      m   <- regmatches(v, regexec(.admDoseModRe, v))[[1L]]
      fun <- if (identical(m[2L], "lag")) "alag" else m[2L]
      paste0(fun, "(", m[3L], ")=", .admToRx(get(v, envir = s)))
    }, character(1))

    # state ICs + their direction derivatives. The IC is evaluated at t = 0, before
    # integration, so its direction derivative is a direct partial (no state chain).
    # Skip any compartment whose IC .rxSens already emitted.
    ic_done <- trimws(sub("\\(0\\)=.*$", "",
                          grep("\\(0\\)=", unlist(strsplit(c(sens_lines, sens2_lines), "\n")),
                               value = TRUE)))
    ic <- character(0)
    for (x in st) {
      x0 <- tryCatch(get(paste0("rx_", x, "_ini_0__"), envir = s), error = function(e) NULL)
      if (is.null(x0)) next
      if (!(x %in% ic_done)) ic <- c(ic, paste0(x, "(0)=", .admToRx(x0)))
      for (p in dirs) {
        cmt <- paste0("rx__sens_", x, "_BY_", p, "__")
        d   <- .admToRx(.Dn(x0, p))
        if (!identical(d, "0") && !(cmt %in% ic_done))
          ic <- c(ic, paste0(cmt, "(0)=", d))
      }
      # Second-order ICs: an IC depending on two directions leaves the cross
      # sensitivity compartment starting at 0 unless d2(x0)/(dp dq) is emitted.
      # Same argument as the first-order block above (the IC is evaluated before
      # integration, so this is a plain double partial, no state chain).
      # Iterated over the SAME canonical pairs the compartments were emitted for
      # -- an IC naming a mirrored pair would declare a compartment that no
      # longer exists.
      if (order >= 2L) for (.i in seq_along(eta_dirs)) {
        p <- eta_dirs[.i]
        for (q in c(eta_dirs[.i:length(eta_dirs)], theta_dirs)) {
          cmt <- paste0("rx__sens_", x, "_BY_", p, "_BY_", q, "__")
          d   <- .admToRx(.Dn(.Dn(x0, p), q))
          if (!identical(d, "0") && !(cmt %in% ic_done))
            ic <- c(ic, paste0(cmt, "(0)=", d))
        }
      }
    }

    # DDE pre-history. A non-constant delay() needs `past(state, tau) <- expr`
    # lines plus the per-sensitivity-compartment histories that .rxSens()
    # accumulates as a side effect (rxode2's .rxDelaySensAugment). NULL for an
    # ordinary model, and for a CONSTANT delay -- but omitting them when they do
    # exist would silently give a wrong sensitivity, so emit them where
    # nlmixr2est's own augmented builder does: after the ODEs/ICs, before the
    # prediction.
    past_lines <- tryCatch(s$..pastLines, error = function(e) NULL)
    if (is.null(past_lines)) past_lines <- character(0)

    f1 <- vapply(dirs, function(p) paste0("rx_f1_", p, "=", .admToRx(.g1(pred, p))),
                 character(1))

    # Cross second-order block: rows = eta directions, columns = all directions.
    # d2_cols[i, b] names the solve column holding d2(pred)/(d eta_i d dir_b).
    #
    # Only the CANONICAL pairs are emitted -- eta row i against directions at or
    # after it -- and the matrix mirrors the eta x eta half onto them, so
    # d2_cols["ETA_2_", "ETA_1_"] and d2_cols["ETA_1_", "ETA_2_"] are the SAME column
    # name. Exact, not an approximation (mixed partials commute). The consumer reads
    # d2_cols purely as a name lookup, so a repeated name needs no handling there.
    f2 <- character(0); d2_cols <- NULL
    if (order >= 2L) {
      n_eta_d <- length(eta_dirs)
      P2 <- do.call(rbind, lapply(seq_len(n_eta_d), function(i)
        cbind(i = i, j = c(i:n_eta_d, if (length(theta_dirs))
                             n_eta_d + seq_along(theta_dirs) else integer(0)))))
      nm2 <- paste0("rx_f2_", eta_dirs[P2[, "i"]], "_", dirs[P2[, "j"]])
      f2  <- vapply(seq_len(nrow(P2)), function(r)
                    paste0(nm2[r], "=",
                           .admToRx(.g2(pred, eta_dirs[P2[r, "i"]], dirs[P2[r, "j"]]))),
                    character(1))
      key <- stats::setNames(nm2, paste(P2[, "i"], P2[, "j"], sep = "|"))
      d2_cols <- matrix(NA_character_, nrow = n_eta_d, ncol = length(dirs),
                        dimnames = list(eta_dirs, dirs))
      for (i in seq_len(n_eta_d)) for (b in seq_along(dirs)) {
        cn <- if (b <= n_eta_d) c(min(i, b), max(i, b)) else c(i, b)
        d2_cols[i, b] <- key[[paste(cn[1L], cn[2L], sep = "|")]]
      }
      if (anyNA(d2_cols)) return(NULL)
    }

    # Endpoint routing for a MULTI-ENDPOINT model: rx_pred_ is
    # CMT-conditional, so the pseudo-compartments need `cmt(cp); cmt(ct);
    # dvid(3,4);` lines as nlmixr2est's inner model does (`..stateInfo` isn't
    # populated on ui$loadPruneSens, so built from ui$predDf instead).
    # Single-endpoint models get no routing lines.
    #
    # BASE states are also declared up front to pin them to compartments
    # 1..n_base -- declared implicitly instead, the rx__sens_* compartments
    # interleave and dvid() resolves to the wrong ones.
    outs <- tryCatch(as.character(ui$predDf$var), error = function(e) character(0))
    multi <- length(outs) > 1L
    # The base count must include linCmt's implicit `central` compartment:
    # rxStateOde() lists only d/dt states (empty for a pure linCmt model) but
    # rxState() reports the linCmt compartment too. Using rxStateOde() numbered a
    # multi-endpoint linCmt model's endpoints one too low, mis-routing the
    # CMT-conditional columns. For a pure-ODE model the two agree.
    st_all <- tryCatch(rxode2::rxState(s), error = function(e) st)
    n_base <- length(st_all)
    # The one case we still cannot number reliably: a model mixing linCmt with
    # EXPLICIT ODE states. The ordering of the implicit linCmt central versus the
    # declared ODE-state cmt() lines is not reproducible from here, so bail to the
    # inner model + FD, as matExp()/indLin() do above.
    if (multi && length(st) > 0L && n_base > length(st)) return(NULL)
    head_lines <- if (multi && length(st)) paste0("cmt(", st, ")") else character(0)
    tail_lines <- if (multi)
      c(paste0("cmt(", outs, ")"),
        paste0("dvid(", paste(n_base + seq_along(outs), collapse = ","), ")"))
    else character(0)

    txt <- paste(c(head_lines, base_ode, dose, sens_lines, sens2_lines, ic, past_lines,
                   paste0("rx_pred_=", .admToRx(pred)), f1, f2, tail_lines), collapse = "\n")
    # Second-order backstop for a PROMOTED solved-form linCmt: rxStateOde() is
    # non-empty and ui$predDf$linCmt is cleared, so neither gate above fires, yet
    # rx_pred_ still resolves through linCmtB -- whose second derivative rxFromSE
    # cannot emit. Refuse before compiling rather than risk a silently wrong f2
    # column; order 1 is unaffected.
    if (order >= 2L && any(grepl("linCmtB", c(f2, sens2_lines), fixed = TRUE)))
      return(NULL)

    txt <- tryCatch(rxode2::rxOptExpr(txt, "admixr2 sensitivity model"),
                    error = function(e) txt)
    # Role-tagged and built outside rxTempDir() -- see .admRxode2(). This model's
    # text is admixr2's own emission, but it shares a parsed md5 with any other
    # build of the same text, and it is compiled with eventSens = "jump".
    mod <- .admRxode2(txt, "admSens", eventSens = "jump")
    rxode2::rxLoad(mod)
    # This rxode2 cannot differentiate a dosing modifier one of our directions
    # feeds -> that column would be silently zero. Refuse the sens model entirely;
    # the caller falls back to a finite-difference gradient.
    if (!.admJumpCovers(mod, s, dirs)) return(NULL)
    list(mod = mod, dirs = dirs, order = order,
         # The iniDf these directions were NUMBERED from -- the promoted one for an
         # order-2 linCmt build. The caller derives its rename_map from this rather
         # than from its own `ui`, so the map that FILLS THETA[k] and the derivative
         # that differentiates it cannot come from different frames.
         ini_used = ini,
         sens_cols = paste0("rx_f1_", eta_dirs),
         theta_sens_cols = if (length(unpaired))
           stats::setNames(paste0("rx_f1_", theta_dirs), unpaired) else NULL,
         # NULL at order 1; the consumer (adfo) treats that as "no dJ available"
         # and keeps its finite-difference pass.
         d2_cols = d2_cols,
         # The theta each direction belongs to, so a consumer can go from a struct
         # theta to its column of d2_cols without re-deriving the pairing: a
         # mu-referenced theta uses its ETA direction, an unpaired one its own.
         theta_dirs = if (length(unpaired))
           stats::setNames(theta_dirs, unpaired) else NULL,
         eta_dirs = eta_dirs)
  }, error = function(e) NULL)
  res
}


# Load (or compile + cache) the sensitivity model.
#
# Returns list(type, mod, sens_cols, theta_sens_cols, rename_map, is_lincmt,
# cache_file) or NULL.
#
#   sens_cols       -- one column per eta, in eta order:      d(pred)/d(eta_i)
#   theta_sens_cols -- named by theta, for the UNPAIRED ones: d(pred)/d(theta_k)
#                      (NULL when the model has none, or when the emitted model
#                       couldn't be built and we fell back to nlmixr2est's inner
#                       model -- the estimators then finite-difference those thetas)
#
# Parameter names the estimators speak -> the model's THETA[j] / ETA[i], for
# one iniDf. Returns list(rename_map, fixed_theta, n_eta), or NULL when the
# frame has no estimated eta.
#
# ONE function because there are TWO frames it can legitimately be asked
# about, and they must not diverge: an order-2 request on a linCmt() model
# promotes the model to explicit ODE form, and .admBuildThetaSens numbers
# its emitted directions from the PROMOTED iniDf, while the rename_map
# that FILLS those THETA[k] columns at solve time used to build from the
# original -- any iniDf difference across the promotion silently mis-slots
# each theta's value against the derivative that differentiates it.
.admSensNameMaps <- function(ini_df) {
  eta_rows <- ini_df[!is.na(ini_df$neta1) & ini_df$neta1 == ini_df$neta2 &
                       !ini_df$fix, , drop = FALSE]
  # Order by neta1 so rename_map's ETA[i] labels line up with
  # .admBuildThetaSens's ETA_i_ directions (numbered after order(neta1));
  # otherwise, for an iniDf whose eta rows are out of neta1 order, sens_cols[i]
  # would report d(pred)/d(eta) for a different eta than rename_map fills ETA[i].
  eta_rows <- eta_rows[order(eta_rows$neta1), , drop = FALSE]
  n_eta    <- nrow(eta_rows)
  if (n_eta == 0L) return(NULL)

  # Indexed by ntheta / neta1, NOT by position among the non-fixed thetas: the
  # sens model's THETA[k] is numbered by ntheta and INCLUDES fixed thetas, so a
  # position-indexed map would put every theta after a fixed one in the wrong slot.
  th_rows    <- ini_df[!is.na(ini_df$ntheta), , drop = FALSE]
  rename_map <- c(
    stats::setNames(paste0("THETA[", th_rows$ntheta, "]"), th_rows$name),
    stats::setNames(paste0("ETA[", seq_len(n_eta), "]"),
                    paste0("eta.", gsub("^eta\\.", "", eta_rows$name))))

  # A FIXED theta is not an estimated parameter, so it never reaches the solve
  # paths -- but the EMITTED sens model still has a THETA[k] slot for it and
  # rxSolve REQUIRES every parameter. Left unset the sens solve errors and
  # returns NULL, silently dropping admc/adfo to a finite-difference gradient.
  # Carry the fixed values so the solve paths can fill those columns
  # (.admFillFixedTheta in simulate.R).
  fix_rows <- th_rows[th_rows$fix, , drop = FALSE]
  fixed_theta <- if (nrow(fix_rows) > 0L)
    stats::setNames(as.numeric(fix_rows$est), paste0("THETA[", fix_rows$ntheta, "]"))
  else numeric(0)

  list(rename_map = rename_map, fixed_theta = fixed_theta, n_eta = n_eta)
}

# Preferred model: admixr2's own direction-set model (.admBuildThetaSens),
# a direction per eta plus one per unpaired theta, eventSens = "jump".
#
# Fallback: nlmixr2est's `ui$foceiModel$inner` (eta columns only),
# recompiled with eventSens = "jump" -- WITHOUT that flag a dosing modifier
# (f/lag/rate/dur) has a sensitivity of exactly ZERO, silently.
#
# `order`: 1L (default) emits d(pred)/d(dir) only -- what admc/adgh need. 2L
# adds the cross block d2(pred)/(d eta d dir) that adfo needs for
# dJ/dtheta; a failed order-2 build falls back to order 1.
# Would .admLoadSensModel() return NULL for this model BY DESIGN?
#
# Distinguishes correct, permanent refusals (below) from genuine failures
# (a compile error, or an unwritable rxTempDir() that discards a model
# that actually compiled) -- the latter would otherwise silently degrade a
# user to a coarser struct-theta gradient. The drivers ask this first:
# by-design NULL stays a message, an unexplained NULL becomes a warning.
#
# A SEPARATE predicate rather than a reason code from .admLoadSensModel():
# that function's body is digested into the sens cache key (.admPkgKey),
# so adding plumbing there would invalidate every cached model for nothing.
.admSensNullByDesign <- function(ui, pinfo = NULL) {
  # 1. No random effects: there is nothing to take a sensitivity with respect to.
  .n_eta <- tryCatch(pinfo$n_eta, error = function(e) NULL)
  if (is.null(.n_eta))
    .n_eta <- tryCatch(sum(!is.na(ui$iniDf$neta1)), error = function(e) NA_integer_)
  if (isTRUE(.n_eta == 0L)) return(TRUE)
  # 2. Ordinal: rx_pred_ is the ordinal log-likelihood, not a category
  #    probability, so its columns differentiate a different function.
  .d <- tryCatch(as.character(ui$predDf$distribution), error = function(e) character(0))
  if (length(.d) > 0L && any(.d %in% c("ordinal", "dordinal"))) return(TRUE)
  # 3. Mixed transformed/untransformed endpoints, or several endpoints that are
  #    transformed DIFFERENTLY: rx_pred_ carries different scales per row with no
  #    per-row map to undo it.
  .tr <- tryCatch(as.character(ui$predDf$transform), error = function(e) character(0))
  .ln <- .tr %in% c("lnorm", "logit", "probit", "boxCox", "tbs",
                    "yeoJohnson", "tbsYj")
  if (length(.ln) > 0L && any(.ln) && !all(.ln)) return(TRUE)
  if (length(.ln) > 0L && all(.ln) && length(.tr) > 1L) {
    .bnd <- tryCatch(
      paste(suppressWarnings(as.numeric(ui$predDf$trLow)),
            suppressWarnings(as.numeric(ui$predDf$trHi))),
      error = function(e) rep("", length(.tr)))
    if (length(unique(.tr)) > 1L || length(unique(.bnd)) > 1L) return(TRUE)
  }
  FALSE
}

.admLoadSensModel <- function(ui, order = 1L) {
  order <- as.integer(order)
  ini_df <- tryCatch(ui$iniDf, error = function(e) NULL)
  if (is.null(ini_df)) return(NULL)
  # ORDINAL endpoints get no sensitivity model. rx_pred_ for `y ~ c(p1, p2)` is the
  # ordinal LOG-LIKELIHOOD, not any one category probability, so its columns
  # differentiate a different function than admixr2 scores. Returning NULL here is
  # the single lever that routes every estimator onto the finite-difference path,
  # rather than gating grad in four drivers separately.
  .d <- tryCatch(as.character(ui$predDf$distribution), error = function(e) character(0))
  if (length(.d) > 0L && any(.d %in% c("ordinal", "dordinal"))) return(NULL)
  # A TRANSFORMED endpoint cannot use the second-order block, so do not build it.
  # .admSimulateSensRows() extracts d2_list and then drops it for a transformed
  # endpoint, deliberately: chaining a second derivative through g() needs
  # g''(z) z_p z_q + g'(z) z_pq. But nothing STOPPED the order-2 build, so
  # `cp ~ lnorm(sd)` compiled and then integrated the cross compartments on every
  # solve only to throw them away -- ~2.5x the integrated system for no gain, on the
  # default path for every lnorm/TBS adfo fit.
  #
  # Demoted HERE, above the cache key, so the order-1 model is shared with the
  # order-1 key rather than duplicated under an "order2" one.
  if (order >= 2L) {
    .tr0 <- tryCatch(as.character(ui$predDf$transform), error = function(e) character(0))
    if (any(.tr0 %in% c("lnorm", "logit", "probit", "boxCox", "tbs",
                        "yeoJohnson", "tbsYj"))) order <- 1L
  }
  .maps <- .admSensNameMaps(ini_df)
  if (is.null(.maps)) return(NULL)
  n_eta       <- .maps$n_eta
  rename_map  <- .maps$rename_map
  fixed_theta <- .maps$fixed_theta

  unpaired <- .admUnpairedThetas(ui)

  # Cache key: the MODEL (ui$lstExpr), the DIRECTION SET (unpaired -- so a model
  # cached before a theta gained its own direction is a miss), .admPkgKey(), and the
  # rxode2 VERSION. NOT digest(inner): ui$foceiModel$inner returns a DIFFERENT object
  # on its first access than on later ones, so digesting it gives an unstable key.
  #
  # The rxode2 version keys the transition where a dosing modifier's jump derivative
  # becomes available (lag()/rate()/dur() gain jumps in 5.1.3): a model fitted on the
  # older rxode2 caches the FD-fallback sens model, and without the version in the
  # key that stale fallback could be served after the upgrade.
  .rx_ver <- tryCatch(as.character(utils::packageVersion("rxode2")),
                      error = function(e) "NA")
  # The key MUST include the iniDf parameter ORDER, not just the model({}) block.
  # `rename_map` numbers THETA[i] by iniDf row order, and `theta_sens_cols` is served
  # straight from the cache. Two models with an identical model({}) block and a
  # reordered ini({}) therefore collided: the second was handed the first's column
  # map and read the wrong sensitivity column -- measured end-to-end as an objective
  # stuck at the starting value with every SE NA and no warning, surviving restarts
  # because rxTempDir() persists. Folding the order into the digest fixes every
  # cache-served field at once.
  .ini_key <- .admIniKey(ui)
  .cacheFile <- file.path(
    rxode2::rxTempDir(),
    paste0("adm-sens-",
           digest::digest(list(ui$lstExpr, unpaired, .ini_key,
                               .admPkgKey(), .rx_ver,
                               paste0("order", order))),
    # NOTE the ORDER-1 FALLBACK is cached under the ORDER-2 key. That is what we want
    # at runtime (an order-2 build that cannot succeed must not be retried on every
    # gradient call), but it means a change to what the order-2 build EMITS would be
    # invisible to an existing entry -- hence .admPkgKey(), which digests the
    # emitter's own source as well as the package version.
           ".rds"))

  # STAYS rxTempDir() -- see the same note in .admLoadModel(). The build
  # directory is applied through .admRxode2()'s `wd =` argument, on the models
  # admixr2 emits itself; setwd()ing the whole load path elsewhere breaks
  # multi-endpoint models.
  .old_wd <- tryCatch(getwd(), error = function(e) NULL)
  on.exit(if (!is.null(.old_wd)) setwd(.old_wd), add = TRUE)
  setwd(rxode2::rxTempDir())

  # pred_tbs is derived BEFORE the cache read, because it must also be applied on a
  # cache HIT. The key digests ui$lstExpr -- the model({}) block only -- but lambda's
  # starting value and its fix() status live in ini({}), so `lam <- fix(0.5)` and
  # `lam <- 0.5` COLLIDE. pred_tbs is what tells .admSimulateSens which lambda to
  # write into the solve and which to invert with, so serving a stale one produced
  # gradients wrong by 1e2-1e4x while the NLL stayed bit-identical. Pure metadata off
  # `ui`, so re-deriving costs nothing.
  .tr <- tryCatch(as.character(ui$predDf$transform), error = function(e) character(0))
  .ln <- .tr %in% c("lnorm", "logit", "probit", "boxCox", "tbs",
                    "yeoJohnson", "tbsYj")
  if (length(.ln) > 0L && any(.ln) && !all(.ln)) {
    # Mixed transformed/untransformed endpoints: rx_pred_ then carries DIFFERENT
    # scales in different rows and the solve paths have no per-row map to undo it.
    # Refuse the sens model so the estimators finite-difference instead -- correct,
    # just slower, and the alternative is a silently wrong gradient.
    return(NULL)
  }
  # ... and equally: transformed endpoints that are not transformed the SAME WAY.
  # pred_tbs below is ONE spec, derived from predDf row 1, and .admSimulateSens()
  # inverts the whole stacked rx_pred_ with it -- so `cp ~ lnorm(a); ct ~ boxCox(b,
  # lam)` applied exp() to ct's Box-Cox rows, and two logitNorm endpoints with
  # different bounds applied endpoint 1's to endpoint 2's rows. The residual path is
  # already per-endpoint, so the gradient described a different function than the NLL
  # scored and the second endpoint converged to the wrong estimate silently.
  #
  # Refuse rather than build a per-row spec: the solve paths would each need a row
  # map, and finite differences are correct today.
  if (length(.ln) > 0L && all(.ln) && length(.tr) > 1L) {
    .bnd <- tryCatch(
      paste(suppressWarnings(as.numeric(ui$predDf$trLow)),
            suppressWarnings(as.numeric(ui$predDf$trHi))),
      error = function(e) rep("", length(.tr)))
    .n_lam <- tryCatch(nrow(ui$iniDf[!is.na(ui$iniDf$err) &
                                       ui$iniDf$err %in% .ADM_ERR_TBS_LAM, ,
                                     drop = FALSE]), error = function(e) 0L)
    if (length(unique(.tr)) > 1L || length(unique(.bnd)) > 1L || .n_lam > 1L)
      return(NULL)
  }

  .pred_tbs <- NULL
  if (length(.ln) > 0L && all(.ln)) {
    .t1 <- .tr[[1L]]
    .yj <- if (identical(.t1, "lnorm")) 0L else unname(.ADM_TBS_YJ[[.t1]])
    .lm <- 0
    .lnm <- NA_character_
    if (.yj %in% c(0L, 1L) && !identical(.t1, "lnorm")) {
      .lr <- ui$iniDf[!is.na(ui$iniDf$err) &
                        ui$iniDf$err %in% .ADM_ERR_TBS_LAM, , drop = FALSE]
      .lm <- if (nrow(.lr) > 0L) as.numeric(.lr$est[1L]) else 1
      # An ESTIMATED lambda moves; this `lam` is only its starting value. The solve
      # paths must use the CURRENT one -- both to fill lambda's parameter column (it
      # is a sigma name, and .admSimulateSens zero-fills those) and to invert with
      # the matching lambda. `lam_name` is how they look it up in pars$sigma_var; NA
      # when lambda is fixed, where the frozen value is already correct.
      if (nrow(.lr) > 0L && !isTRUE(.lr$fix[1L])) .lnm <- as.character(.lr$name[1L])
    }
    .pred_tbs <- list(
      lam = .lm, yj = .yj, lam_name = .lnm,
      lo = suppressWarnings(as.numeric(ui$predDf$trLow[1L]  %||% 0)),
      hi = suppressWarnings(as.numeric(ui$predDf$trHi[1L]   %||% 1)))
    if (!is.finite(.pred_tbs$lo)) .pred_tbs$lo <- 0
    if (!is.finite(.pred_tbs$hi)) .pred_tbs$hi <- 1
  }

  # Session cache, same reasons as .admLoadModel; gated on the disk file so
  # a cache clear still forces a rebuild.
  #
  # The four fields below are re-derived on a hit, NOT trusted from the
  # cached object: the key covers ini NAMES/fix flags but not VALUES, so
  # `lam <- 0.5` vs `lam <- 0.7` would collide and serve a stale pred_tbs
  # (gradient 1e2-1e4x wrong, objective unchanged). Name maps are re-derived
  # from the frame the CACHED build numbered directions from, not this
  # call's `ini_df` -- they differ only across an order-2 linCmt promotion.
  # Falls back to the caller's maps when the entry predates ini_used.
  .fromCached <- function(m) {
    if (is.null(m$ini_used) || identical(m$ini_used, ini_df))
      return(list(rename_map = rename_map, fixed_theta = fixed_theta))
    .mm <- .admSensNameMaps(m$ini_used)
    if (is.null(.mm)) return(list(rename_map = rename_map, fixed_theta = fixed_theta))
    list(rename_map = .mm$rename_map, fixed_theta = .mm$fixed_theta)
  }

  .smemo <- get0(basename(.cacheFile), envir = .adm_sens_env, inherits = FALSE)
  if (!is.null(.smemo) && file.exists(.cacheFile) && .admRxLoadAll(.smemo)) {
    .mp <- .fromCached(.smemo)
    .smemo$cache_file  <- .cacheFile
    .smemo$rename_map  <- .mp$rename_map
    .smemo$fixed_theta <- .mp$fixed_theta
    .smemo$pred_tbs    <- .pred_tbs
    return(.smemo)
  }

  if (file.exists(.cacheFile)) {
    result <- tryCatch({
      m <- readRDS(.cacheFile)
      # Assert the payload's SHAPE before trusting the load. .admRxLoadAll is
      # upstream's load step and is a no-op on anything not rxode2-classed, so it
      # returns TRUE for a file that is not a sens result at all; see the same
      # guard in .admLoadModel. Every branch below stores the model in $mod.
      if (!inherits(m$mod, "rxode2") || !.admRxLoadAll(m)) NULL else m
    }, error = function(e) NULL)
    if (!is.null(result)) {
      # Overwrite the worker-inherited fields from the parent's fresh derivation
      # rather than trusting the file. A parallel WORKER reads this same file and
      # cannot re-derive, so what the parent writes here is what the worker gets.
      # (sens_cols / dirs are NOT re-derived: they are keyed by `unpaired` in the
      # cache path, so a hit has the same direction set by construction.)
      .mp <- .fromCached(result)           # ... from the frame THIS entry was built on
      result$cache_file  <- .cacheFile
      result$rename_map  <- .mp$rename_map
      result$fixed_theta <- .mp$fixed_theta
      result$pred_tbs    <- .pred_tbs      # see the derivation above -- key collision
      return(.admCacheAssign(basename(.cacheFile), result, .adm_sens_env))
    }
  }

  # A beta endpoint's rx_pred_ is llikBeta(DV, b1, b2) -- the LOG-LIKELIHOOD, which
  # is what FOCEI maximises but NOT what admixr2 moment-matches. Emit sensitivities
  # of the derived mean mu = b1/(b1+b2) instead. The state-sensitivity chain is
  # shared across the system, so this costs nothing extra (see .admBuildThetaSens).
  .dist <- tryCatch(as.character(ui$predDf$distribution), error = function(e) character(0))
  .pred_expr <- NULL

  # A COUNT endpoint has the same shape as beta: `y ~ pois(cp)` puts the
  # LOG-LIKELIHOOD in rx_pred_ and differentiates that rather than the mean, and both
  # need DV, which an aggregate fit does not have -- so the solve returned NULL. admc
  # coped (it falls back to FD) but .adghGrad returned all-NA, killing adgh at
  # iteration 0, and .admGradBatch returned all-NA, silently giving admc a ZERO
  # Hessian and no standard errors. Emit sensitivities of the count MEAN instead.
  if (any(.dist %in% c("pois", "dpois", "binom", "dbinom", "nbinomMu", "dnbinomMu"))) {
    .mv <- tryCatch(.admEndpointVar(ui, which(.dist %in% c("pois", "dpois", "binom",
                                                           "dbinom", "nbinomMu",
                                                           "dnbinomMu"))[1L]),
                    error = function(e) NULL)
    .se <- tryCatch(ui$loadPruneSens, error = function(e) NULL)
    if (is.null(.mv) || is.null(.se)) return(NULL)
    .pred_expr <- tryCatch(.se[[.mv]], error = function(e) NULL)
    if (is.null(.pred_expr)) return(NULL)
  }

  if (any(.dist %in% c("beta", "dbeta"))) {
    .bp <- tryCatch(.admBetaPair(ui), error = function(e) NULL)
    .se <- tryCatch(ui$loadPruneSens, error = function(e) NULL)
    if (is.null(.bp) || is.null(.se)) return(NULL)
    .pred_expr <- tryCatch({
      .e1 <- .se[[.bp[[1L]]]]; .e2 <- .se[[.bp[[2L]]]]
      if (is.null(.e1) || is.null(.e2)) NULL else .e1 / (.e1 + .e2)
    }, error = function(e) NULL)
    if (is.null(.pred_expr)) return(NULL)
  }

  built <- .admBuildThetaSens(ui, unpaired, .pred_expr, order = order)
  # An order-2 build can fail where order 1 succeeds -- linCmt has no second
  # linCmtB derivative, and a model can simply be too large. Retry at order 1
  # rather than dropping the caller to finite differences for BOTH orders: adfo
  # then gets exactly the model it had before, and its FD struct-theta pass.
  if (is.null(built) && order >= 2L)
    built <- .admBuildThetaSens(ui, unpaired, .pred_expr, order = 1L)
  if (!is.null(built)) {
    # Adopt the frame the directions were actually numbered from. Identical to
    # ini_df except across an order-2 linCmt promotion -- see .admSensNameMaps().
    if (!is.null(built$ini_used) && !identical(built$ini_used, ini_df)) {
      .m2 <- .admSensNameMaps(built$ini_used)
      if (is.null(.m2)) return(NULL)
      rename_map <- .m2$rename_map; fixed_theta <- .m2$fixed_theta
      n_eta      <- .m2$n_eta
    }
    result <- list(type = "dirs", mod = built$mod,
                   sens_cols = built$sens_cols,
                   theta_sens_cols = built$theta_sens_cols,
                   d2_cols = built$d2_cols,
                   theta_dirs = built$theta_dirs,
                   eta_dirs = built$eta_dirs,
                   order = built$order,
                   dirs = built$dirs,
                   rename_map = rename_map,
                   fixed_theta = fixed_theta,
                   # The frame the directions were NUMBERED from -- the promoted one
                   # across an order-2 linCmt build. Stored so the cache-hit paths
                   # can re-derive the maps from it too: they re-derive correctly but
                   # had only `ini_df`, the UNPROMOTED frame, so a warm hit rebuilt
                   # the map from a different iniDf than the cold build used.
                   ini_used = built$ini_used,
                   is_lincmt = .admIsLinCmtMod(built$mod),
                   cache_file = .cacheFile)
  } else if (!is.null(.pred_expr)) {
    # The count/beta branches above exist BECAUSE nlmixr2est's inner model puts the
    # log-likelihood in rx_pred_ and differentiates that, not the mean. Falling back
    # to it here would hand the estimators exactly the object those branches were
    # written to avoid. NULL routes every estimator onto finite differences instead,
    # the same lever the ordinal guard at the top of this function pulls.
    return(NULL)
  } else {
    result <- .admSensFromInner(ui, rename_map, fixed_theta, n_eta, .cacheFile)
    if (is.null(result)) return(NULL)
  }

  # rx_pred_ is on the MODELLING scale (e.g. log for lnorm), but .admSimulate
  # always reads the natural-scale column -- so .admGrad and .admNLL would
  # disagree on scale for any transformed endpoint (logit/probit/boxCox/
  # yeoJohnson, not just lnorm) unless back-transformed via the chain rule.
  #
  # Back-transform spec: (lambda, yj, lo, hi) for .admTBSi()/.admTBSid().
  # NULL for an untransformed endpoint, which leaves every solve path
  # byte-identical.
  result$pred_tbs <- .pred_tbs

  # DDE: force pure dop853 for the SENSITIVITY solve. A delayed model's
  # variational compartments trip rxode2's hasDelay AutoSwitch into ros4,
  # whose dense delay history drifts once delay() reads recorded (not
  # pre-)history; dop853's dense output reproduces the base solve exactly.
  # Mirrors nlmixr2est's ed03b8dfc fix for the same failure.
  #
  # Stored on the result and folded into the cache key (.admPkgKey) since a
  # parallel worker reads the file directly. NULL for an ordinary model.
  result$solve_args <- if (isTRUE(tryCatch(
        rxode2::rxModelVars(result$mod)$flags[["hasDelay"]] == 1L,
        error = function(e) FALSE)))
    list(method = "dop853", stiff2 = 0L, dense = TRUE) else NULL

  .admCacheWrite(result, .cacheFile, "sensitivity model")
  .admCacheAssign(basename(.cacheFile), result, .adm_sens_env)
}

.admIsLinCmtMod <- function(mod) {
  mv <- tryCatch(rxode2::rxModelVars(mod), error = function(e) NULL)
  if (is.null(mv)) FALSE else any(grepl("linCmtB", mv$model, fixed = TRUE))
}

# Fallback sens model: nlmixr2est's `ui$foceiModel$inner`. Eta columns only --
# its sensitivity block is keyed on etas, so there are no theta columns and the
# estimators finite-difference the unpaired thetas, as they always did.
# Recompiled with eventSens = "jump" (see .admLoadSensModel's header).
.admSensFromInner <- function(ui, rename_map, fixed_theta, n_eta, cacheFile) {
  # .admLoadSensModel already pinned $foceiModel (the Windows finalizer guard).
  .focei_model <- tryCatch(ui$foceiModel, error = function(e) NULL)
  inner <- .focei_model$inner
  if (is.null(inner)) return(NULL)

  lhs <- tryCatch(inner$lhs, error = function(e) NULL)
  if (is.null(lhs)) return(NULL)
  sens_cols <- lhs[grepl("sens_rx_pred.*ETA|sens.*pred.*BY.*ETA", lhs, ignore.case = TRUE)]
  if (length(sens_cols) == 0L) return(NULL)
  eta_idx <- suppressWarnings(as.integer(regmatches(sens_cols, regexpr("[0-9]+", sens_cols))))
  if (anyNA(eta_idx)) return(NULL)
  sens_cols <- sens_cols[order(eta_idx)]
  if (length(sens_cols) != n_eta) return(NULL)

  .normMod <- tryCatch(rxode2::rxModelVars(inner)$model[["normModel"]],
                       error = function(e) NULL)
  # THE collision this package was most exposed to: `.normMod` is nlmixr2est's
  # OWN inner model text (rxModelVars(inner)$model[["normModel"]]), and admixr2
  # rebuilds it with eventSens = "jump" where nlmixr2est built it with a
  # different one. Anonymous, that is the same parsed md5 in the same directory
  # emitting different C -- the later build wins for BOTH packages. Role-tagged
  # and built in admixr2's own directory it cannot collide with either
  # nlmixr2est's build or an anonymous one. See .admRxode2().
  mod <- if (!is.null(.normMod))
    tryCatch({ m <- .admRxode2(.normMod, "admSensInner", eventSens = "jump")
               rxode2::rxLoad(m); m },
             error = function(e) NULL)
  else NULL
  if (is.null(mod))
    mod <- tryCatch({ rxode2::rxLoad(inner); inner }, error = function(e) NULL)
  if (is.null(mod))
    mod <- tryCatch({ m <- .admRxode2(inner, "admSensInnerRaw"); rxode2::rxLoad(m); m },
                    error = function(e) NULL)
  if (is.null(mod)) return(NULL)

  # Same guard as the emitter, over the inner model's directions (etas only): if
  # this rxode2 cannot differentiate a dosing modifier an ETA feeds, that eta's
  # column is identically zero. Refuse the sens model so the estimators use FD.
  .s <- tryCatch(ui$loadPruneSens, error = function(e) NULL)
  if (!is.null(.s) &&
      !.admJumpCovers(mod, .s, paste0("ETA_", seq_len(n_eta), "_"))) return(NULL)

  # theta_sens_cols = NULL: the inner model has no theta directions, so the
  # estimators finite-difference the unpaired thetas. fixed_theta still travels so
  # the solve paths fill a fixed theta's THETA[k] (the inner model needs it too).
  list(type = "inner", mod = mod, sens_cols = sens_cols,
       theta_sens_cols = NULL, rename_map = rename_map, fixed_theta = fixed_theta,
       is_lincmt = .admIsLinCmtMod(mod), cache_file = cacheFile)
}
