# Load (or compile + cache) the rxode2 simulation model.
# Compiled DLL is cached to disk via qs2, keyed by model digest.
.admLoadModel <- function(ui) {
  # Accessing $simulationModel (below) caches the compiled model in
  # ui$meta$.simModelBase as a side effect -- a live, self-referential rxode2
  # object that breaks nlmixr2's ui-cloning during fit assembly. Drop it (and any
  # sibling artifacts) on every exit so the ui stays in the canonical state
  # nlmixr2 expects; see .admDropSimModelMeta() for the full rationale.
  on.exit(.admDropSimModelMeta(ui), add = TRUE)
  # Record the rxode2 model(s) this load registers so .admFitTeardown can reclaim
  # them later, even when called outside an estimator fit (test setup, datagen).
  .before_reg <- .admRegistrySnapshot()
  on.exit(.admTrackRegistry(.before_reg), add = TRUE)
  .model_key <- digest::digest(ui$lstExpr)
  .cacheFile <- file.path(
    rxode2::rxTempDir(),
    paste0("adm-sim-", .model_key, ".qs2")
  )
  if (file.exists(.cacheFile)) {
    mod <- tryCatch(qs2::qs_read(.cacheFile), error = function(e) NULL)
    load_ok <- !is.null(mod) &&
      tryCatch({ rxode2::rxLoad(mod); TRUE }, error = function(e) FALSE)
    if (load_ok) {
      return(mod)
    }
    tryCatch(file.remove(.cacheFile), error = function(e) NULL)
  }
  # rxode2 compilation calls setwd() internally -- save/restore to avoid
  # "cannot change working directory" error on first compile (Windows).
  .old_wd <- tryCatch(getwd(), error = function(e) NULL)
  on.exit(if (!is.null(.old_wd)) setwd(.old_wd), add = TRUE)
  setwd(rxode2::rxTempDir())
  mod <- rxode2::rxode2(ui)$simulationModel
  tryCatch(suppressWarnings(qs2::qs_save(mod, .cacheFile)), error = function(e) NULL)
  rxode2::rxLoad(mod)
  mod
}

# Remove transient rxode2 model objects that $simulationModel / $foceiModel leave
# behind in ui$meta.
#
# nlmixr2's output machinery (nlmixr2CreateOutputFromUi -> ... -> nmObjGet.*)
# deep-clones the ui with nlmixr2est's internal .cloneEnv(), which recurses into
# every environment-valued member and has no cycle detection. rxode2's compiled
# model objects hold a back-reference to the global .rxModels registry
# (registry -> model -> .rx -> .rxModels -> registry ...), so cloning one loops
# forever -- surfacing as "evaluation nested too deeply: infinite recursion"
# (interactive) or "node stack overflow" (batch). A normal nlmixr2 fit never
# hits this because its estimators do not populate ui$meta with these objects;
# admixr2 does, because it simulates via $simulationModel. Keeping our ui clean
# is the in-framework fix: no wrapping of nlmixr2's code, we just do not feed it
# a ui it was never designed to clone. Safe because admixr2 simulates via its
# own cached model (the return value of .admLoadModel), and rxode2 regenerates
# these lazily if any downstream method needs them.
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

# Structural thetas with no usable mu-referenced eta ("unpaired"): the ones whose
# gradient cannot come from an eta sensitivity column and would otherwise be
# finite-differenced. These are the thetas that get their OWN sensitivity
# direction (THETA_j_) in the sens model.
#
# Uses .admMuRefPairs() -- the SAME map pinfo$struct_has_eta is built from -- so
# the set of thetas the estimators route through the theta columns and the set the
# sens model actually builds columns for cannot drift apart. That includes the
# shared-eta guard: a theta whose eta appears in another parameter is unpaired.
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

# Which dosing modifiers does the model actually use? (f/lag/rate/dur live in the
# pruned env as rx_<mod>_<state>_; rxode2 stores lag() as alag().)
.admDoseMods <- function(s) {
  v <- grep("^rx_(f|lag|alag|rate|dur)_.+_$", ls(envir = s, all.names = TRUE), value = TRUE)
  m <- sub("^rx_(f|lag|alag|rate|dur)_.+_$", "\\1", v)
  unique(ifelse(m == "alag", "lag", m))
}

# Can this rxode2 build differentiate every dosing modifier that one of OUR
# directions actually feeds?
#
# A direction entering f()/lag()/rate()/dur() only has a sensitivity if rxode2
# attaches its analytic variational jumps at dose times (eventSens = "jump");
# otherwise its column is silently ZERO. Support is version-dependent -- rxode2
# 5.1.2 has no lag() jumps, 5.1.3 does -- so FEATURE-DETECT rather than
# version-compare: the compiled model carries eventSensInfo$derivs, one table per
# modifier, and an unsupported one has zero rows.
#
# The test must be per DIRECTION, not merely "the model has a lag()": nlmixr2est's
# inner model has no theta directions, so for `alag(depot) = exp(tlag)` its
# derivs$lag is legitimately empty -- nothing depends on lag there, so nothing can
# be wrong. Only a modifier some direction differentiates to non-zero needs cover.
#
# FALSE -> the caller returns NULL for the whole sens model and the estimators
# fall back to a finite-difference gradient (correct, if slower). Far better than
# the alternative: an identically-zero gradient component, silently.
.admJumpCovers <- function(mod, s, dirs) {
  vars <- grep("^rx_(f|lag|alag|rate|dur)_.+_$", ls(envir = s, all.names = TRUE), value = TRUE)
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
    key <- sub("^rx_(f|lag|alag|rate|dur)_.+_$", "\\1", v)
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
# A mu-referenced theta needs no direction of its own -- d(pred)/d(theta) ==
# d(pred)/d(eta) -- so it reuses its eta's column for free. Only a theta with no
# usable eta (eta-less, non-mu-referenced, or one whose eta is shared across
# parameters) gets its own direction. Sigmas get none (they never enter the
# prediction). This is the same direction/linking scheme nlmixr2est's fast-focei
# uses (.foceiAnalyticDirections), FIRST-ORDER only: admixr2's MC/quadrature
# moments need d(pred)/d(dir) and nothing higher, so the O(ndir^2) second-order
# tier that FOCEI's Laplace term requires is skipped entirely.
#
# Two branches:
#   * ODE    -- rxode2::.rxSens() augments the system with the variational
#               (state-sensitivity) compartments for each direction; the emitted
#               prediction chain is
#                 rx_f1_<dir> = d(pred)/d(dir) + sum_states d(pred)/d(state)
#                                                 * d(state)/d(dir)
#   * linCmt -- there are no states to augment (.rxSens errors), so the state sum
#               drops out and D(pred, dir) alone is emitted: symengine resolves it
#               through rxode2's linCmtB derivative rules (.rxD$linCmtB), which
#               give d(linCmt)/d(micro parameter) in closed form. This is why the
#               direction set works for linCmt at first order even though
#               nlmixr2est's (second-order) augmented outer model cannot build it.
#
# Compiled with eventSens = "jump" so a parameter entering a dosing modifier
# (f/lag/rate/dur) gets rxode2's analytic variational jumps at dose times --
# without it such a sensitivity is silently ZERO. State initial conditions and
# their direction derivatives are emitted too (a parameter-dependent IC otherwise
# starts every sensitivity compartment at 0).
#
# Returns list(mod, dirs, sens_cols, theta_sens_cols) or NULL on any failure (the
# caller then falls back to nlmixr2est's inner model + FD for the thetas).
.admBuildThetaSens <- function(ui, unpaired) {
  s <- tryCatch(ui$loadPruneSens, error = function(e) NULL)
  if (is.null(s)) return(NULL)
  st <- tryCatch(rxode2::rxStateOde(s), error = function(e) NULL)
  if (is.null(st)) return(NULL)

  ini      <- tryCatch(ui$iniDf, error = function(e) NULL)
  if (is.null(ini)) return(NULL)
  eta_rows <- ini[!is.na(ini$neta1) & ini$neta1 == ini$neta2 & !ini$fix, , drop = FALSE]
  eta_rows <- eta_rows[order(eta_rows$neta1), , drop = FALSE]
  th_rows  <- ini[!is.na(ini$ntheta), , drop = FALSE]

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
    if (length(st) > 0L) {
      rxode2::.rxJacobian(s, c(st, dirs))
      sens_lines <- rxode2::.rxSens(s, dirs)
      if (length(sens_lines) == 0L) return(NULL)
    }
    pred <- get("rx_pred_", envir = s)
    .Dn  <- function(e, v) symengine::D(e, symengine::S(v))
    .sn1 <- function(j, p) symengine::S(paste0("rx__sens_", j, "_BY_", p, "__"))
    # linCmt: st is empty, so the state sum drops out and D(pred, dir) alone
    # resolves through the linCmtB derivative rules.
    .g1 <- function(ex, p) {
      e <- .Dn(ex, p)
      for (j in st) e <- e + .Dn(ex, j) * .sn1(j, p)
      e
    }

    base_ode <- if (length(st))
      vapply(st, function(x)
        paste0("d/dt(", x, ")=", .admToRx(get(paste0("rx__d_dt_", x, "__"), envir = s))),
        character(1)) else character(0)

    # dosing modifiers (bioavailability, lag, rate, duration) live in the pruned
    # env as rx_<mod>_<state>_ and are NOT part of rx__d_dt_*; rxode2 stores lag()
    # as alag().
    dos_vars <- grep("^rx_(f|lag|alag|rate|dur)_.+_$", ls(envir = s, all.names = TRUE),
                     value = TRUE)
    dose <- vapply(dos_vars, function(v) {
      m   <- regmatches(v, regexec("^rx_(f|lag|alag|rate|dur)_(.+)_$", v))[[1L]]
      fun <- if (identical(m[2L], "lag")) "alag" else m[2L]
      paste0(fun, "(", m[3L], ")=", .admToRx(get(v, envir = s)))
    }, character(1))

    # state ICs + their direction derivatives. The IC is evaluated at t = 0, before
    # integration, so its direction derivative is a direct partial (no state chain).
    # Skip any compartment whose IC .rxSens already emitted.
    ic_done <- trimws(sub("\\(0\\)=.*$", "",
                          grep("\\(0\\)=", unlist(strsplit(sens_lines, "\n")), value = TRUE)))
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

    # Endpoint routing for a MULTI-ENDPOINT model. Its rx_pred_ is a CMT-conditional
    # expression (`CMT==3 ? ... : ...`), so the solve needs the endpoint
    # pseudo-compartments declared and mapped -- nlmixr2est's inner model ends with
    #   cmt(cp); cmt(ct); dvid(3,4);
    # admixr2 tags each unit's observations with its output's cmt
    # (.admBuildEvFull(tag_cmt = TRUE)) precisely so the solve can disambiguate.
    #
    # The dvid indices are (number of BASE states) + endpoint position: the
    # rx__sens_* variational compartments are not counted. `..stateInfo` (which
    # nlmixr2est uses) is not populated on ui$loadPruneSens, so build the lines from
    # ui$predDf instead. Single-endpoint models need no routing -- every observation
    # is that endpoint -- and get no lines, which is what they already did.
    # The BASE states must also be declared up front (`cmt(central); cmt(periph);`),
    # as nlmixr2est's inner model does: that pins them to compartments 1..n_base so
    # the endpoint numbering is n_base + i. Declared implicitly (by d/dt alone) the
    # rx__sens_* compartments get interleaved and dvid() resolves to the wrong ones.
    outs <- tryCatch(as.character(ui$predDf$var), error = function(e) character(0))
    multi <- length(outs) > 1L
    # Endpoint pseudo-compartments are numbered AFTER the base compartments, so the
    # dvid() indices are (number of base compartments) + endpoint position. The base
    # count must include linCmt's implicit `central` compartment: rxStateOde() lists
    # only d/dt states (empty for a pure linCmt model) but rxState() reports the
    # linCmt compartment too. Using rxStateOde() here numbered a multi-endpoint
    # linCmt model's endpoints one too low (dvid(1,2) instead of nlmixr2est's
    # dvid(2,3)), mis-routing the CMT-conditional rx_pred_/rx_f1_ columns; rxState()
    # matches nlmixr2est's inner model exactly. For a pure-ODE model rxState() ==
    # rxStateOde(), so ODE numbering is unchanged.
    st_all <- tryCatch(rxode2::rxState(s), error = function(e) st)
    n_base <- length(st_all)
    # The one case we still cannot number reliably: a model mixing linCmt with
    # EXPLICIT ODE states (n_base > number of d/dt states). The ordering of the
    # implicit linCmt central versus the declared ODE-state cmt() lines is not
    # reproducible from here, so bail to the inner model + FD (correct, just slower),
    # as matExp()/indLin() do above.
    if (multi && length(st) > 0L && n_base > length(st)) return(NULL)
    head_lines <- if (multi && length(st)) paste0("cmt(", st, ")") else character(0)
    tail_lines <- if (multi)
      c(paste0("cmt(", outs, ")"),
        paste0("dvid(", paste(n_base + seq_along(outs), collapse = ","), ")"))
    else character(0)

    txt <- paste(c(head_lines, base_ode, dose, sens_lines, ic, past_lines,
                   paste0("rx_pred_=", .admToRx(pred)), f1, tail_lines), collapse = "\n")
    txt <- tryCatch(rxode2::rxOptExpr(txt, "admixr2 sensitivity model"),
                    error = function(e) txt)
    mod <- rxode2::rxode2(txt, eventSens = "jump")
    rxode2::rxLoad(mod)
    # This rxode2 cannot differentiate a dosing modifier one of our directions
    # feeds -> that column would be silently zero. Refuse the sens model entirely;
    # the caller falls back to a finite-difference gradient.
    if (!.admJumpCovers(mod, s, dirs)) return(NULL)
    list(mod = mod, dirs = dirs,
         sens_cols = paste0("rx_f1_", eta_dirs),
         theta_sens_cols = if (length(unpaired))
           stats::setNames(paste0("rx_f1_", theta_dirs), unpaired) else NULL)
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
#                       could not be built and we fell back to nlmixr2est's inner
#                       model -- the estimators then finite-difference those thetas)
#
# Preferred model: admixr2's own direction-set model (.admBuildThetaSens), which
# carries a direction per eta plus one per unpaired theta, and is compiled with
# eventSens = "jump".
#
# Fallback: nlmixr2est's `ui$foceiModel$inner`, which only ever emits eta columns
# (its sensitivity block is keyed on etas), recompiled with eventSens = "jump" --
# WITHOUT that flag a parameter entering a dosing modifier (f/lag/rate/dur) has a
# sensitivity of exactly ZERO, silently, because FOCEI computes event/dose
# sensitivities separately (its `predNoLhs` FD model) and admixr2 reads the inner
# model's columns directly.
#
# Pinning: ui$foceiModel creates companion objects ($outer, $predOnly,
# $predNoLhs) with live C++ DLL pointers. rxUi is a locked environment so
# assign(..., envir = ui) fails silently. Instead we pin the full foceiModel
# result in .adm_pin_env (package-level, always writable), keyed by model
# digest. This keeps companions alive for the session and prevents Windows GC
# finalizer heap corruption (STATUS_HEAP_CORRUPTION / -1073740940).
.admLoadSensModel <- function(ui) {
  .model_key <- digest::digest(ui$lstExpr)
  .sens_key  <- paste0("sens_", .model_key)

  # Record the rxode2 models this load registers (the foceiModel companions read
  # via ui$foceiModel, plus the inner sens model) so .admFitTeardown can reclaim
  # them -- including when loaded outside an estimator fit (test setup, datagen).
  .before_reg <- .admRegistrySnapshot()
  on.exit(.admTrackRegistry(.before_reg), add = TRUE)

  # In-memory cache: avoids disk read and rxLoad on repeat calls within a session.
  .cached <- tryCatch(get(.sens_key, envir = .adm_pin_env, inherits = FALSE),
                      error = function(e) NULL)
  if (!is.null(.cached)) return(.cached)

  ini_df <- tryCatch(ui$iniDf, error = function(e) NULL)
  if (is.null(ini_df)) return(NULL)
  eta_rows <- ini_df[!is.na(ini_df$neta1) & ini_df$neta1 == ini_df$neta2 &
                       !ini_df$fix, , drop = FALSE]
  n_eta    <- nrow(eta_rows)
  if (n_eta == 0L) return(NULL)

  unpaired <- .admUnpairedThetas(ui)

  # Parameter names the estimators speak -> the model's THETA[j] / ETA[i]. Indexed
  # by ntheta / neta1, NOT by position among the non-fixed thetas: the sens model's
  # THETA[k] is numbered by ntheta and INCLUDES fixed thetas, so a position-indexed
  # map would put every theta after a fixed one in the wrong slot.
  th_rows    <- ini_df[!is.na(ini_df$ntheta), , drop = FALSE]
  rename_map <- c(
    stats::setNames(paste0("THETA[", th_rows$ntheta, "]"), th_rows$name),
    stats::setNames(paste0("ETA[", seq_len(n_eta), "]"),
                    paste0("eta.", gsub("^eta\\.", "", eta_rows$name))))

  # A FIXED theta is not an estimated parameter, so it never reaches the solve
  # paths (pinfo carries only the estimated ones) -- but the EMITTED sens model
  # still has a THETA[k] slot for it (the model text references every theta) and
  # rxSolve REQUIRES every parameter. Left unset the sens solve errors and returns
  # NULL, which silently drops admc/adfo to a finite-difference gradient and, worse,
  # made .adghGrad skip the study entirely. Carry the fixed values so the solve
  # paths can fill those columns (.admFillFixedTheta in simulate.R).
  fix_rows <- th_rows[th_rows$fix, , drop = FALSE]
  fixed_theta <- if (nrow(fix_rows) > 0L)
    stats::setNames(as.numeric(fix_rows$est), paste0("THETA[", fix_rows$ntheta, "]"))
  else numeric(0)

  # Cache key: the MODEL (ui$lstExpr), the DIRECTION SET (unpaired -- so a model
  # cached before a theta gained its own direction is a miss), a schema tag, and
  # the rxode2 VERSION.
  # NOT digest(inner): ui$foceiModel$inner returns a DIFFERENT object on its first
  # access than on later ones, so digesting it gives an unstable key. The schema
  # tag ("+fixed-theta") makes a cache written before the fixed-theta fix a miss:
  # a parallel worker reads this file directly and cannot re-derive, so it would
  # otherwise inherit a stale rename_map / NULL fixed_theta and silently diverge
  # from the sequential fit.
  # The rxode2 version keys the transition where a dosing-modifier's jump
  # derivative becomes available (e.g. lag()/rate()/dur() gain jumps in 5.1.3):
  # a model fitted on the older rxode2 caches the FD-fallback sens model
  # (theta_sens_cols = NULL), and without the version in the key that stale
  # fallback could be served after the upgrade instead of rebuilding the now-full
  # jump model. rxTempDir() is session-scoped and compiled models are
  # version-stamped, so this is belt-and-braces for a pinned persistent tempdir --
  # but it makes the 5.1.2 -> 5.1.3 handoff automatic regardless.
  .rx_ver <- tryCatch(as.character(utils::packageVersion("rxode2")),
                      error = function(e) "NA")
  .cacheFile <- file.path(
    rxode2::rxTempDir(),
    paste0("adm-sens-",
           digest::digest(list(ui$lstExpr, unpaired, "dirs-jump+fixed-theta", .rx_ver)),
           ".qs2"))

  .old_wd <- tryCatch(getwd(), error = function(e) NULL)
  on.exit(if (!is.null(.old_wd)) setwd(.old_wd), add = TRUE)
  setwd(rxode2::rxTempDir())

  if (file.exists(.cacheFile)) {
    result <- tryCatch({ m <- qs2::qs_read(.cacheFile); rxode2::rxLoad(m$mod); m },
                       error = function(e) NULL)
    if (!is.null(result)) {
      # Overwrite the worker-inherited fields from the parent's fresh derivation
      # rather than trusting the file. A parallel WORKER reads this same file and
      # cannot re-derive, so what the parent writes here is what the worker gets;
      # a stale position-indexed rename_map or a NULL fixed_theta would silently
      # diverge the parallel fit from the sequential one. (sens_cols / dirs are NOT
      # re-derived: they are keyed by `unpaired` in the cache path, so a hit is
      # guaranteed to have the same direction set.)
      result$cache_file  <- .cacheFile
      result$rename_map  <- rename_map
      result$fixed_theta <- fixed_theta
      tryCatch(assign(.sens_key, result, envir = .adm_pin_env), error = function(e) NULL)
      return(result)
    }
  }

  # Pin ui$foceiModel even when the emitted model does not need it. The estimator
  # drivers access $foceiModel anyway (nlmixr2CreateOutputFromUi / double-compile
  # prevention), which creates companion objects ($outer, $predOnly, $predNoLhs)
  # holding live C++ DLL pointers. Unpinned they are immediately GC-eligible and
  # their finalizers unload DLLs mid-allocation -> Windows heap corruption
  # (STATUS_HEAP_CORRUPTION / -1073740940). rxUi is a locked env, so they are
  # pinned in .adm_pin_env (package-level, always writable) instead.
  tryCatch(assign(paste0("focei_", .model_key),
                  suppressMessages(ui$foceiModel), envir = .adm_pin_env),
           error = function(e) NULL)

  built <- .admBuildThetaSens(ui, unpaired)
  if (!is.null(built)) {
    result <- list(type = "dirs", mod = built$mod,
                   sens_cols = built$sens_cols,
                   theta_sens_cols = built$theta_sens_cols,
                   dirs = built$dirs,
                   rename_map = rename_map,
                   fixed_theta = fixed_theta,
                   is_lincmt = .admIsLinCmtMod(built$mod),
                   cache_file = .cacheFile)
  } else {
    result <- .admSensFromInner(ui, rename_map, fixed_theta, n_eta, .cacheFile)
    if (is.null(result)) return(NULL)
  }

  # suppressWarnings: the sens model is compiled inside admixr2, so its environment
  # chain references the package namespace and serialising it warns "'package:
  # admixr2' may not be available when loading". Harmless -- a worker reloads the
  # DLL via rxLoad(), not from the serialised env (.admLoadModel does the same).
  tryCatch(suppressWarnings(qs2::qs_save(result, .cacheFile)), error = function(e) NULL)
  tryCatch(assign(.sens_key, result, envir = .adm_pin_env), error = function(e) NULL)
  result
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
  mod <- if (!is.null(.normMod))
    tryCatch({ m <- rxode2::rxode2(.normMod, eventSens = "jump"); rxode2::rxLoad(m); m },
             error = function(e) NULL)
  else NULL
  if (is.null(mod))
    mod <- tryCatch({ rxode2::rxLoad(inner); inner }, error = function(e) NULL)
  if (is.null(mod))
    mod <- tryCatch({ m <- rxode2::rxode2(inner); rxode2::rxLoad(m); m },
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
