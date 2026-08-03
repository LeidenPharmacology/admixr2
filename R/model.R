# Load (or compile + cache) the rxode2 simulation model.
# Compiled DLL is cached to disk with saveRDS(), keyed by model digest.
.admLoadModel <- function(ui) {
  # Accessing $simulationModel (below) caches the compiled model in
  # ui$meta$.simModelBase as a side effect -- a live, self-referential rxode2
  # object that breaks nlmixr2's ui-cloning during fit assembly. Drop it (and any
  # sibling artifacts) on every exit so the ui stays in the canonical state
  # nlmixr2 expects; see .admDropSimModelMeta() for the full rationale.
  on.exit(.admDropSimModelMeta(ui), add = TRUE)
  .model_key <- digest::digest(ui$lstExpr)
  .cacheFile <- file.path(
    rxode2::rxTempDir(),
    paste0("adm-sim-", .model_key, ".rds")
  )
  # In-session registry (see .adm_model_env): a disk cache HIT still costs a
  # readRDS plus a dyn.load, and hands back a fresh object whose finalizer will
  # eventually unload the shared library this one is using.
  #
  # Gated on the disk cache STILL EXISTING, so the registry is a fast path for a
  # valid cache rather than a replacement for one. Without that check, clearing
  # the cache (rxode2::rxClean(), or the cold-session regression test that uses
  # it to force a recompile) would no longer force anything: the memo would be
  # served and the recompile path silently stop being exercised. The extra
  # file.exists() is nothing against the readRDS it avoids.
  .memo <- get0(.model_key, envir = .adm_model_env, inherits = FALSE)
  if (!is.null(.memo) && file.exists(.cacheFile) &&
      tryCatch({ rxode2::rxLoad(.memo); TRUE }, error = function(e) FALSE))
    return(.memo)
  if (file.exists(.cacheFile)) {
    mod <- tryCatch(readRDS(.cacheFile), error = function(e) NULL)
    load_ok <- !is.null(mod) &&
      tryCatch({ rxode2::rxLoad(mod); TRUE }, error = function(e) FALSE)
    if (load_ok) {
      return(.admCacheAssign(.model_key, mod, .adm_model_env))
    }
    tryCatch(file.remove(.cacheFile), error = function(e) NULL)
  }
  # rxode2 compilation calls setwd() internally -- save/restore to avoid
  # "cannot change working directory" error on first compile (Windows).
  .old_wd <- tryCatch(getwd(), error = function(e) NULL)
  on.exit(if (!is.null(.old_wd)) setwd(.old_wd), add = TRUE)
  setwd(rxode2::rxTempDir())
  mod <- rxode2::rxode2(ui)$simulationModel
  tryCatch(suppressWarnings(saveRDS(mod, .cacheFile)), error = function(e) NULL)
  rxode2::rxLoad(mod)
  .admCacheAssign(.model_key, mod, .adm_model_env)
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

# Promote a linCmt() model to the equivalent explicit ODE system, so it can carry
# SECOND-order sensitivities.
#
# Why this exists: linCmt()'s first derivative resolves through rxode2's linCmtB
# rules, but there is no second one -- rxFromSE cannot emit the nested linCmtB
# derivative, so d2f/(d eta d theta) is unavailable and adfo would be stuck with
# its finite-difference struct-theta pass. nlmixr2est hits the identical wall and
# refuses linCmt outright for its analytic gradient and covariance. Promotion
# sidesteps it: the ODE form has ordinary state sensitivities, so .rxSens()
# expands it to any order.
#
# `rxode2::linToOde()` (exported) does the translation, but its output cannot be
# used as-is: it names the prediction `rxLinCmt`, and an `rx`-prefixed lhs is
# RESERVED -- ui$loadPruneSens then dies with "syntax errors" on the model rxode2
# itself just generated, at every order. Renaming that one variable is the whole
# fix. Measured on a 1-cmt oral model: the promoted solve reproduces the analytic
# linCmt prediction to 1.8e-08 relative, and its second-order block agrees with a
# central difference of its own first-order columns to 1.8e-06 (the FD floor).
#
# NULL on any failure -- the caller then falls back to the first-order model, i.e.
# exactly the behaviour before this existed.
#
# Only ever called for an order-2 request. admc/adgh stay on the SOLVED form,
# which is faster and needs nothing higher than first order.
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

# The dosing-modifier variables rxode2 emits into the pruned sens env as
# rx_<mod>_<state>_ (f/lag/rate/dur; lag() is stored as alag()). ONE regex, used
# everywhere a dose modifier is found or its name extracted, so a change to
# rxode2's naming is a single edit rather than four. Group 1 = modifier, 2 = state.
.admDoseModRe <- "^rx_(f|lag|alag|rate|dur)_(.+)_$"

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
# A mu-referenced theta needs no direction of its own -- d(pred)/d(theta) ==
# d(pred)/d(eta) -- so it reuses its eta's column for free. Only a theta with no
# usable eta (eta-less, non-mu-referenced, or one whose eta is shared across
# parameters) gets its own direction. Sigmas get none (they never enter the
# prediction). This is the same direction/linking scheme nlmixr2est's fast-focei
# uses (.foceiAnalyticDirections).
#
# `order = 2L` additionally emits the CROSS second-order block
#   rx_f2_<eta_i>_<dir_b> = d2(pred)/(d eta_i d dir_b)
# which is what adfo -- and only adfo -- needs. Its objective is built on
# V_pred = J Omega J' + resid(mu, diag(J Omega J')) with J = df/d(eta)|_0, so a
# structural theta's gradient needs dJ/d(theta) = d2f/(d eta d theta): a second
# derivative, which is why adfo alone still finite-differenced its struct thetas
# while admc/adgh read theirs off the first-order columns.
#
# The block is deliberately ASYMMETRIC -- eta directions x ALL directions -- and
# that is the whole economy of it. nlmixr2est's builder expands the full
# symmetric triangle over every direction (it needs the Laplace Hessian);
# adfo needs no theta x theta pair at all, because the objective is only ever
# differentiated ONCE with respect to a theta. rxode2's rxExpandSens2_() accepts
# two different direction sets, so the cross block can be requested directly:
# for 2 states / 2 etas / 1 unpaired theta that is 12 second-order compartments
# instead of 18. admixr2 also needs no rx_rvar*/rx_rsig* variance chains at all
# (errmodel.R derives the residual analytically from (mu, var_f)), which is the
# larger part of what FOCEI's order-2 build spends its compile time on.
#
# admc/adgh stay at order 1: their moments need d(pred)/d(dir) and nothing
# higher, so they must not pay for the second-order compartments.
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
# `pred_expr`: the model expression whose sensitivities to emit, as a symengine
# object. Defaults to `rx_pred_`, which is the prediction for every ordinary
# endpoint -- but NOT for a likelihood-form endpoint, where rxode2 puts the
# LOG-LIKELIHOOD there (`llikBeta(DV, b1, b2)`). FOCEI wants exactly that; admixr2
# moment-matches and needs the MEAN, so such an endpoint passes its own derived
# expression (beta: b1/(b1+b2)).
#
# This works because .rxSens() builds the variational compartments for the WHOLE
# ODE SYSTEM -- d(state)/d(dir), once -- not for a particular target. .g1() then
# applies the chain rule to any expression on top of them. A derived prediction is
# therefore no harder than rx_pred_; only the direct partial differs.
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

  # linCmt() carries no SECOND derivative (see .admLinCmtToOde), so an order-2
  # request promotes the model to its ODE form and builds from that instead. The
  # rest of this function then runs on an ordinary ODE model and needs no special
  # case. Promotion failing is not fatal -- NULL here means the caller retries at
  # order 1, where linCmt works exactly as it always has.
  #
  # ui$predDf only (NOT ui$predDfFocei): reaching for the focei-side frame pulls on
  # nlmixr2est's model machinery, which is the documented Windows GC/finalizer
  # hazard. A PROMOTED solved-form linCmt (predDf flag cleared, real ODE states)
  # is caught instead by the linCmtB text check on the emitted model, below.
  if (order >= 2L) {
    .lin <- tryCatch(as.logical(ui$predDf$linCmt), error = function(e) NULL)
    if (isTRUE(any(.lin, na.rm = TRUE))) {
      ui <- .admLinCmtToOde(ui)
      if (is.null(ui)) return(NULL)
      s  <- tryCatch(ui$loadPruneSens, error = function(e) NULL)
      if (is.null(s)) return(NULL)
      st <- tryCatch(rxode2::rxStateOde(s), error = function(e) NULL)
      if (is.null(st) || length(st) == 0L) return(NULL)
      ini <- tryCatch(ui$iniDf, error = function(e) NULL)
      if (is.null(ini)) return(NULL)
      # eta_rows/th_rows are re-derived below from `ini`; the promoted model has
      # the SAME ini({}) block, so every direction index and the caller's
      # rename_map (built from the original ui) stay valid.
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
      # The cross block only: vars = eta directions, vars2 = every direction.
      # rxExpandSens2_ takes the two sets independently, so no theta x theta
      # compartment is generated. .rxSens also assigns the rx__sens_*_BY_*_BY_*__
      # symbols into `s` as a side effect, which is what makes them resolvable in
      # .g2 below -- so this must run BEFORE the chains are built.
      if (order >= 2L) {
        sens2_lines <- rxode2::.rxSens(s, eta_dirs, dirs)
        if (length(sens2_lines) == 0L) return(NULL)
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
    # Second-order chain. Differentiating .g1(ex, q) w.r.t. p picks up three
    # terms: the direct partial, the first-order state paths of the ALREADY
    # chained expression, and the second-order state sensitivities themselves.
    # Same construction as nlmixr2est's .g2 -- it is simply the product rule
    # applied to the first-order chain, and getting any one term wrong yields a
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
      if (order >= 2L) for (p in eta_dirs) for (q in dirs) {
        cmt <- paste0("rx__sens_", x, "_BY_", p, "_BY_", q, "__")
        d   <- .admToRx(.Dn(.Dn(x0, p), q))
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

    # Cross second-order block: rows = eta directions, columns = all directions.
    # expand.grid varies `i` fastest, so filling the name matrix column-major puts
    # eta i in row i and direction b in column b -- i.e. d2_cols[i, b] names the
    # column holding d2(pred)/(d eta_i d dir_b), which is column i of dJ/d(dir_b).
    f2 <- character(0); d2_cols <- NULL
    if (order >= 2L) {
      P2  <- expand.grid(i = eta_dirs, j = dirs, stringsAsFactors = FALSE)
      nm2 <- paste0("rx_f2_", P2$i, "_", P2$j)
      f2  <- vapply(seq_len(nrow(P2)), function(r)
                    paste0(nm2[r], "=", .admToRx(.g2(pred, P2$i[r], P2$j[r]))),
                    character(1))
      d2_cols <- matrix(nm2, nrow = length(eta_dirs), ncol = length(dirs),
                        dimnames = list(eta_dirs, dirs))
    }

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

    txt <- paste(c(head_lines, base_ode, dose, sens_lines, sens2_lines, ic, past_lines,
                   paste0("rx_pred_=", .admToRx(pred)), f1, f2, tail_lines), collapse = "\n")
    # Second-order backstop for a PROMOTED solved-form linCmt: rxStateOde() is
    # non-empty and ui$predDf$linCmt is cleared, so neither gate above fires, yet
    # rx_pred_ still resolves through linCmtB -- whose second derivative rxFromSE
    # cannot emit (nlmixr2est documents the same failure). Refuse before compiling
    # rather than risk a silently wrong f2 column; order 1 is unaffected and keeps
    # linCmt working exactly as it does today.
    if (order >= 2L && any(grepl("linCmtB", c(f2, sens2_lines), fixed = TRUE)))
      return(NULL)

    txt <- tryCatch(rxode2::rxOptExpr(txt, "admixr2 sensitivity model"),
                    error = function(e) txt)
    mod <- rxode2::rxode2(txt, eventSens = "jump")
    rxode2::rxLoad(mod)
    # This rxode2 cannot differentiate a dosing modifier one of our directions
    # feeds -> that column would be silently zero. Refuse the sens model entirely;
    # the caller falls back to a finite-difference gradient.
    if (!.admJumpCovers(mod, s, dirs)) return(NULL)
    list(mod = mod, dirs = dirs, order = order,
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
# `order`: 1L (default) emits d(pred)/d(dir) only -- what admc/adgh need. 2L adds
# the cross second-order block d2(pred)/(d eta d dir) that adfo needs for dJ/dtheta;
# it costs extra state compartments, so only adfo asks for it, and a failed order-2
# build falls back to order 1 (adfo then keeps its finite-difference pass).
.admLoadSensModel <- function(ui, order = 1L) {
  order <- as.integer(order)
  ini_df <- tryCatch(ui$iniDf, error = function(e) NULL)
  if (is.null(ini_df)) return(NULL)
  # ORDINAL endpoints get no sensitivity model. rx_pred_ for `y ~ c(p1, p2)` is the
  # ordinal LOG-LIKELIHOOD, not any one category probability, so its sensitivity
  # columns differentiate a different function than admixr2 scores -- the same
  # class of mismatch that made lnorm's gradient ~200x wrong. Returning NULL here
  # is the single lever that routes every estimator onto the finite-difference
  # path (audited at ~5e-06), rather than gating grad in four drivers separately.
  .d <- tryCatch(as.character(ui$predDf$distribution), error = function(e) character(0))
  if (length(.d) > 0L && any(.d %in% c("ordinal", "dordinal"))) return(NULL)
  eta_rows <- ini_df[!is.na(ini_df$neta1) & ini_df$neta1 == ini_df$neta2 &
                       !ini_df$fix, , drop = FALSE]
  # Order by neta1 so rename_map's ETA[i] labels below line up with
  # .admBuildThetaSens's ETA_i_ directions (which it numbers after order(neta1));
  # otherwise, for an iniDf whose eta rows are out of neta1 order, sens_cols[i]
  # would report d(pred)/d(eta) for a different eta than rename_map fills ETA[i].
  eta_rows <- eta_rows[order(eta_rows$neta1), , drop = FALSE]
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
  # The key MUST include the iniDf parameter ORDER, not just the model({}) block.
  # `rename_map` numbers THETA[i] by iniDf row order, and `theta_sens_cols` -- which
  # names the emitted rx_f1_THETA_j_ columns -- is served straight from the cache.
  # Two models with an identical model({}) block and a reordered ini({}) therefore
  # collided: the second was handed the first's column map and read the wrong
  # sensitivity column. Measured end-to-end (adgh, same model, ini order swapped):
  # objective 1081.08 with a clean cache vs 2355.77 when served from the other
  # model's entry -- i.e. stuck at the starting value, every SE NA, and no warning
  # of any kind. rxTempDir() persists ACROSS SESSIONS, so this survived restarts.
  #
  # This is the same failure class the pred_tbs block below documents and fixes by
  # re-deriving on a cache hit; the reasoning had simply not been carried across to
  # the field that names the columns. Folding the order into the digest fixes every
  # cache-served field at once rather than one at a time.
  .ini_key <- tryCatch({
    .i <- ui$iniDf
    paste(paste(.i$name, collapse = "|"),
          paste(as.integer(.i$fix), collapse = "|"),
          paste(.i$err, collapse = "|"), sep = "//")
  }, error = function(e) "")
  .cacheFile <- file.path(
    rxode2::rxTempDir(),
    paste0("adm-sens-",
           digest::digest(list(ui$lstExpr, unpaired, .ini_key,
                               "dirs-jump+fixed-theta+dde+predtbs+derivpred+tbslam+countpred+inikey+order2lin", .rx_ver,
                               paste0("order", order))),
    # NOTE the ORDER-1 FALLBACK is cached under the ORDER-2 key. That is what we
    # want at runtime (an order-2 build that cannot succeed must not be retried on
    # every gradient call), but it means a change to what the order-2 build EMITS
    # is invisible until the schema tag above is bumped -- a stale "linCmt cannot
    # do order 2" entry outlived exactly that change once already.
           ".rds"))

  .old_wd <- tryCatch(getwd(), error = function(e) NULL)
  on.exit(if (!is.null(.old_wd)) setwd(.old_wd), add = TRUE)
  setwd(rxode2::rxTempDir())

  # pred_tbs is derived BEFORE the cache read, because it must also be applied on a
  # cache HIT. The cache key digests ui$lstExpr -- the model({}) block only -- but
  # lambda's starting value and its fix() status live in ini({}), so
  # `lam <- fix(0.5)` and `lam <- 0.5` COLLIDE on one key. pred_tbs is what tells
  # .admSimulateSens which lambda to write into the solve and which to invert with,
  # so serving a stale one produced gradients wrong by 1e2-1e4x (one component with
  # the wrong sign) while the NLL stayed bit-identical -- nothing warned, and the
  # optimizer simply stalled near its starting values. Pure metadata off `ui`, so
  # re-deriving costs nothing; same reason rename_map/fixed_theta are re-derived.
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
  # inverts the whole stacked rx_pred_ with it. So `cp ~ lnorm(a); ct ~ boxCox(b,
  # lam)` -- which passes the mixed-vs-untransformed guard above, since both are
  # transformed -- applied exp() to ct's Box-Cox rows; two logitNorm endpoints with
  # different (trLow, trHi) applied endpoint 1's bounds to endpoint 2's rows; and
  # two boxCox endpoints with separate lambdas used endpoint 1's lambda for both.
  # The residual path is already per-endpoint (errmodel.R reads predDf$trLow[i] and
  # carries per-row lam/yj), so under the default grad = "sens" the gradient
  # described a different function than the NLL scored and the second endpoint
  # converged to the wrong estimate with no error and no warning.
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
      # paths must use the CURRENT one -- both to fill lambda's parameter column
      # (it is a sigma name, and .admSimulateSens zero-fills those, so rx_pred_ was
      # built with lambda = 0, i.e. a plain log transform) and to invert with the
      # matching lambda. `lam_name` is how they look it up in pars$sigma_var; NA
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

  # Session cache, ahead of the disk cache and on the same key (.cacheFile already
  # digests the model, the direction set, the ini order, the schema tag, the
  # rxode2 version and the ORDER -- upstream's "composite key covering everything
  # that changes the emitted model"). Same three reasons as .admLoadModel: skip a
  # readRDS + dyn.load, stop minting a finalizer-bearing object per call, and keep
  # the compiled model reachable. Gated on the disk file so a cache clear still
  # forces a rebuild.
  #
  # The four fields below are re-derived on a hit exactly as the disk path does --
  # NOT trusted from the cached object. `.cacheFile`'s key includes the ini NAMES
  # and fix flags but not their VALUES, so `lam <- 0.5` and `lam <- 0.7` collide;
  # serving a stale pred_tbs left the objective bit-identical and the gradient
  # 1e2-1e4x wrong. Sharing that hazard with the disk path is the point.
  .smemo <- get0(basename(.cacheFile), envir = .adm_sens_env, inherits = FALSE)
  if (!is.null(.smemo) && file.exists(.cacheFile) &&
      tryCatch({ rxode2::rxLoad(.smemo$mod); TRUE }, error = function(e) FALSE)) {
    .smemo$cache_file  <- .cacheFile
    .smemo$rename_map  <- rename_map
    .smemo$fixed_theta <- fixed_theta
    .smemo$pred_tbs    <- .pred_tbs
    return(.smemo)
  }

  if (file.exists(.cacheFile)) {
    result <- tryCatch({ m <- readRDS(.cacheFile); rxode2::rxLoad(m$mod); m },
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

  # A COUNT endpoint has the same shape as beta: `y ~ pois(cp)` emits
  #   rx_pred_     = llikPois(DV, cp)
  #   rx_f1_ETA_1_ = ... llikPoisDlambda(DV, cp) * d(central)/d(eta)
  # i.e. rx_pred_ is the LOG-LIKELIHOOD and the sensitivity columns differentiate
  # it, not the mean -- and both need DV, which an aggregate fit does not have, so
  # the solve returned NULL. admc coped (it falls back to FD) but .adghGrad returned
  # all-NA, which killed adgh at iteration 0 with the default grad = "analytical",
  # and .admGradBatch returned all-NA, which silently gave admc a ZERO Hessian and
  # therefore no standard errors at all. Emit sensitivities of the count MEAN -- the
  # distribution's argument, which is an ordinary model variable -- exactly as the
  # beta branch below emits them for b1/(b1+b2).
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
                   is_lincmt = .admIsLinCmtMod(built$mod),
                   cache_file = .cacheFile)
  } else if (!is.null(.pred_expr)) {
    # The count/beta branches above exist BECAUSE nlmixr2est's inner model puts the
    # log-likelihood in rx_pred_ (llikPois(DV, cp)) and differentiates that, not the
    # mean -- and needs a DV an aggregate fit does not have. Falling back to it here
    # would hand the estimators exactly the object those branches were written to
    # avoid: .adghGrad returns all-NA and .admGradBatch a zero Hessian. NULL routes
    # every estimator onto finite differences instead, the same lever the ordinal
    # guard at the top of this function pulls.
    return(NULL)
  } else {
    result <- .admSensFromInner(ui, rename_map, fixed_theta, n_eta, .cacheFile)
    if (is.null(result)) return(NULL)
  }

  # The sensitivity model's rx_pred_ is on the endpoint's MODELLING scale, which
  # for an lnorm endpoint is the LOG scale:
  #
  #   cp ~ add(a)   ->  rx_pred_ = <cp>              (natural)
  #   cp ~ lnorm(a) ->  rx_pred_ = log(<cp>)         (log!)
  #
  # .admSimulate (the NLL path) always reads the natural-scale output column, so
  # an lnorm model had .admGrad differentiating log(f) while .admNLL scored f --
  # a gradient of a different function entirely. It went unnoticed because lnorm
  # appears in no gradient test. The solve paths back-transform with the chain
  # rule (d(exp(g))/dp = exp(g)*dg/dp) when this flag is set.
  # EVERY transformed endpoint puts rx_pred_ on the MODELLING scale, not just
  # lnorm: logit/probit/boxCox/yeoJohnson all emit rx_pred_ = rxTBS(f, ...). The
  # original fix handled only the log case, so the other four had a sens model
  # that predicted the transformed scale -- which is why their analytic gradient
  # was unusable (the FD audit showed them as NA). Generalise to the full inverse
  # transform; lnorm is exactly yj = 0 with lambda = 0 (Box-Cox's log branch), so
  # it stays bit-identical to what pred_log did. `.tr`/`.ln`/`.pred_tbs` are all
  # derived ABOVE the cache read -- see there for why.
  # Back-transform spec: (lambda, yj, lo, hi) for .admTBSi()/.admTBSid(). NULL for
  # an untransformed endpoint, which leaves every solve path byte-identical.
  result$pred_tbs <- .pred_tbs

  # DDE: force pure dop853 for the SENSITIVITY solve.
  #
  # A delay() model's sensitivity system is the base ODEs plus one variational
  # compartment per state per direction, all of them delayed -- stiff enough to trip
  # rxode2's hasDelay AutoSwitch composite (dop853+ros4) into its ros4 leg, whose
  # dense delay-history is inaccurate for this system. The symptom is not an error:
  # the augmented prediction agrees with the fit's own solve for the first
  # observations and then drifts, once delay() starts reading the RECORDED (solved)
  # history rather than the pre-history. dop853's 8th-order dense output reproduces
  # the base solve exactly, so forcing it also keeps the gradient's predictions
  # consistent with .admSimulate's (which is NOT augmented, does not trip, and is
  # deliberately left alone).
  #
  # This mirrors nlmixr2est's ed03b8dfc, which found and fixed the same failure in
  # its own augmented-sensitivity solve. Stored on the result -- and folded into the
  # cache schema tag above -- because a parallel worker reads the cache file directly
  # and cannot re-derive it. NULL for an ordinary model, which leaves every existing
  # solve call byte-for-byte as it was.
  result$solve_args <- if (isTRUE(tryCatch(
        rxode2::rxModelVars(result$mod)$flags[["hasDelay"]] == 1L,
        error = function(e) FALSE)))
    list(method = "dop853", stiff2 = 0L, dense = TRUE) else NULL

  # suppressWarnings: the sens model is compiled inside admixr2, so its environment
  # chain references the package namespace and serialising it warns "'package:
  # admixr2' may not be available when loading". Harmless -- a worker reloads the
  # DLL via rxLoad(), not from the serialised env (.admLoadModel does the same).
  tryCatch(suppressWarnings(saveRDS(result, .cacheFile)), error = function(e) NULL)
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
