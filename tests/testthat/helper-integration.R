# Integration test helpers.
# Sourced automatically before any test file.
# All rxSolve calls happen inside memoised setup functions; individual tests
# operate only on pre-computed scalars/matrices — no rxSolve during test
# execution.

# ---- Model: 1-cmt with 2 etas -----------------------------------------------
# Two etas (cl, v) ensure sobol(dim=2) returns a matrix (not a vector).
# Analytic: C(t) = (D/V)*exp(-CL/V*t). True: CL=5, V=20, D=100.

one_cmt_fn <- function() {
  ini({
    tcl     <- log(5)  ; label("Log CL")
    tv      <- log(20) ; label("Log V")
    add.err <- 0.1     ; label("Additive SD")
    eta.cl  ~ 0.09
    eta.v   ~ 0.04
  })
  model({
    cl <- exp(tcl + eta.cl)
    v  <- exp(tv  + eta.v)
    d/dt(central) <- -(cl / v) * central
    cp <- central / v
    cp ~ add(add.err)
  })
}

# ---- Analytic mean for a bolus dose -----------------------------------------
.one_cmt_mean <- function(cl, v, dose, times)
  (dose / v) * exp(-(cl / v) * times)

# ---- Memoised setup caches ---------------------------------------------------
# All caches live here (helper env) so they persist across test files.
# If defined inside a test file, the binding is cleaned up when that file's
# eval-env is GC'd — which unloads the rxMod DLL and crashes later rxSolve calls.
.int_cache        <- NULL
.int_grad_cache   <- NULL
.int_grad_lin_kappa_cache <- NULL
.int_irmc_exact_kappa_cache <- NULL
.int_lincmt_cache <- NULL
.int_irmc_cache   <- NULL
.int_plot_cache   <- NULL
.int_cov_cache    <- NULL
.int_adfo_cache        <- NULL
.int_adfo_kappa_cache  <- NULL
.int_adfo_cov_cache    <- NULL

# ---- linCmt model (1-cmt, 2 etas) -------------------------------------------

one_cmt_lincmt_fn <- function() {
  ini({
    tcl     <- log(5)  ; label("Log CL")
    tv      <- log(20) ; label("Log V")
    add.err <- 0.1     ; label("Additive SD")
    eta.cl  ~ 0.09
    eta.v   ~ 0.04
  })
  model({
    cl <- exp(tcl + eta.cl)
    v  <- exp(tv  + eta.v)
    linCmt() ~ add(add.err)
  })
}

# ---- Model: 1-cmt with 2 etas + 1 unpaired struct theta ----------------------

one_cmt_kappa_fn <- function() {
  ini({
    tcl     <- log(5)  ; label("Log CL")
    tv      <- log(20) ; label("Log V")
    # Intentionally unpaired from any eta so IRMC takes the kappa path.
    tsc     <- log(1)  ; label("Log scale")
    add.err <- 0.1     ; label("Additive SD")
    eta.cl  ~ 0.09
    eta.v   ~ 0.04
  })
  model({
    cl <- exp(tcl + eta.cl)
    v  <- exp(tv  + eta.v)
    # Neutral scale keeps the same analytic mean while exercising has_kappa=TRUE.
    sc <- exp(tsc)
    d/dt(central) <- -(cl / v) * central
    cp <- sc * central / v
    cp ~ add(add.err)
  })
}

# ---- Grad setup: passes ui so struct_has_eta = TRUE for struct thetas --------
# Separate from .int_setup() which calls .admParseIniDf(iniDf) without ui,
# routing struct thetas through unpaired FD (2x gradient error).

.int_grad_setup <- function(n_sim = 500L, seed = 42L) {
  if (!is.null(.int_grad_cache) &&
      .int_grad_cache$n_sim == n_sim && .int_grad_cache$seed == seed)
    return(.int_grad_cache)

  skip_if_not_installed("rxode2")

  ui <- suppressMessages(tryCatch(
    rxode2::rxode2(one_cmt_fn),
    error = function(e) NULL
  ))
  if (is.null(ui)) skip("rxode2 model parse failed")

  pinfo      <- admixr2:::.admParseIniDf(ui$iniDf, ui)
  output_var <- "cp"

  rxMod <- tryCatch(admixr2:::.admLoadModel(ui), error = function(e) NULL)
  if (is.null(rxMod)) skip("Model compilation failed")

  times  <- c(0.5, 1, 2, 4)
  E_true <- .one_cmt_mean(5, 20, 100, times)
  V_true <- diag((0.3 * E_true)^2)

  study <- list(E = E_true, V = V_true, n = 200L, times = times,
                ev = rxode2::et(amt = 100))
  study <- admixr2:::.admNormaliseStudy(study, "s")
  study$ev_full <- study$ev |> rxode2::et(study$times)
  studies <- list(s = study)

  z_list      <- admixr2:::.admMakeZ(n_sim, pinfo, 1L, "sobol")
  params_list <- admixr2:::.admMakeParamsList(n_sim, pinfo, 1L)
  vec         <- admixr2:::.admBuildOptVec(pinfo)

  p0   <- vec$p0
  h_fd <- 1e-3

  # Sens model (ODE always has one); re-rxLoad rxMod after to restore DLL state.
  sensModel <- tryCatch(admixr2:::.admLoadSensModel(ui), error = function(e) NULL)
  rxode2::rxLoad(rxMod)

  # ---- MC gradient -----------------------------------------------------------
  g_ana_p0 <- admixr2:::.admGrad(p0, pinfo, studies, z_list, rxMod, output_var,
                                 params_list, cores = 1L, h = h_fd,
                                 sensModel = NULL, use_central = TRUE)

  # ---- IRMC inner gradient ---------------------------------------------------
  pars0 <- admixr2:::.admUnpack(p0, pinfo)
  irmc_proposals <- lapply(seq_along(studies), function(si)
    admixr2:::.adirmcProposal(
      rxMod, pars0$struct, pinfo$sigma_names,
      pinfo$sigma_is_prop, pinfo$sigma_is_lnorm,
      pars0$omega, omega_expansion = 2,
      studies[[si]], z_list[[si]], output_var,
      params_list[[si]], cores = 1L,
      pinfo$eta_col_names,
      has_kappa         = pinfo$has_kappa,
      struct_transforms = pinfo$struct_transforms,
      struct_eta_idx    = pinfo$struct_eta_idx,
      use_grad          = TRUE
    )
  )
  proposals_ok <- !any(vapply(irmc_proposals, is.null, logical(1)))

  g_irmc_ana <- NULL
  g_irmc_fd  <- NULL
  if (proposals_ok) {
    g_irmc_ana <- admixr2:::.adirmcInnerGrad(p0, pinfo, studies, irmc_proposals)
    names(g_irmc_ana) <- names(p0)

    h_irmc <- 1e-4
    g_irmc_fd <- vapply(seq_along(p0), function(k) {
      ph <- p0; ph[k] <- ph[k] + h_irmc
      pl <- p0; pl[k] <- pl[k] - h_irmc
      nh <- admixr2:::.adirmcNLL(ph, pinfo, studies, irmc_proposals)
      nl <- admixr2:::.adirmcNLL(pl, pinfo, studies, irmc_proposals)
      (nh - nl) / (2 * h_irmc)
    }, double(1))
    names(g_irmc_fd) <- names(p0)
  }

  .int_grad_cache <<- list(
    ui = ui, pinfo = pinfo, rxMod = rxMod, sensModel = sensModel,
    studies = studies, z_list = z_list,
    params_list = params_list, vec = vec,
    output_var = output_var, times = times,
    E_true = E_true, n_sim = n_sim, seed = seed,
    g_ana_p0 = g_ana_p0,
    irmc_proposals = irmc_proposals, proposals_ok = proposals_ok,
    g_irmc_ana = g_irmc_ana, g_irmc_fd = g_irmc_fd
  )
  .int_grad_cache
}

# ---- IRMC linearized-kappa grad setup ----------------------------------------

.int_irmc_kappa_setup_impl <- function(kappa_method, n_sim, seed) {
  skip_if_not_installed("rxode2")

  ui <- suppressMessages(tryCatch(
    rxode2::rxode2(one_cmt_kappa_fn),
    error = function(e) NULL
  ))
  if (is.null(ui)) skip("rxode2 model parse failed")

  pinfo      <- admixr2:::.admParseIniDf(ui$iniDf, ui)
  output_var <- "cp"

  rxMod <- tryCatch(admixr2:::.admLoadModel(ui), error = function(e) NULL)
  if (is.null(rxMod)) skip("Model compilation failed")

  times  <- c(0.5, 1, 2, 4)
  E_true <- .one_cmt_mean(5, 20, 100, times)
  V_true <- diag((0.3 * E_true)^2)

  study <- list(E = E_true, V = V_true, n = 200L, times = times,
                ev = rxode2::et(amt = 100))
  study <- admixr2:::.admNormaliseStudy(study, "s")
  study$ev_full <- study$ev |> rxode2::et(study$times)
  studies <- list(s = study)

  z_list      <- admixr2:::.admMakeZ(n_sim, pinfo, 1L, "sobol")
  params_list <- admixr2:::.admMakeParamsList(n_sim, pinfo, 1L)
  vec         <- admixr2:::.admBuildOptVec(pinfo)

  p0    <- vec$p0
  pars0 <- admixr2:::.admUnpack(p0, pinfo)
  irmc_proposals <- lapply(seq_along(studies), function(si)
    admixr2:::.adirmcProposal(
      rxMod, pars0$struct, pinfo$sigma_names,
      pinfo$sigma_is_prop, pinfo$sigma_is_lnorm,
      pars0$omega, omega_expansion = 2,
      studies[[si]], z_list[[si]], output_var,
      params_list[[si]], cores = 1L,
      pinfo$eta_col_names,
      has_kappa         = pinfo$has_kappa,
      kappa_method      = kappa_method,
      struct_transforms = pinfo$struct_transforms,
      struct_eta_idx    = pinfo$struct_eta_idx,
      use_grad          = TRUE
    )
  )
  proposals_ok <- !any(vapply(irmc_proposals, is.null, logical(1)))

  g_irmc_ana <- NULL
  g_irmc_fd  <- NULL
  if (proposals_ok) {
    g_irmc_ana <- admixr2:::.adirmcInnerGrad(p0, pinfo, studies, irmc_proposals)
    names(g_irmc_ana) <- names(p0)

    h_irmc <- 1e-4
    g_irmc_fd <- vapply(seq_along(p0), function(k) {
      ph <- p0; ph[k] <- ph[k] + h_irmc
      pl <- p0; pl[k] <- pl[k] - h_irmc
      nh <- admixr2:::.adirmcNLL(ph, pinfo, studies, irmc_proposals)
      nl <- admixr2:::.adirmcNLL(pl, pinfo, studies, irmc_proposals)
      (nh - nl) / (2 * h_irmc)
    }, double(1))
    names(g_irmc_fd) <- names(p0)
  }

  list(
    ui = ui, pinfo = pinfo, rxMod = rxMod,
    studies = studies, z_list = z_list,
    params_list = params_list, vec = vec,
    output_var = output_var, times = times,
    E_true = E_true, n_sim = n_sim, seed = seed,
    irmc_proposals = irmc_proposals, proposals_ok = proposals_ok,
    g_irmc_ana = g_irmc_ana, g_irmc_fd = g_irmc_fd
  )
}

.int_grad_lin_kappa_setup <- function(n_sim = 500L, seed = 42L) {
  if (!is.null(.int_grad_lin_kappa_cache) &&
      .int_grad_lin_kappa_cache$n_sim == n_sim &&
      .int_grad_lin_kappa_cache$seed == seed)
    return(.int_grad_lin_kappa_cache)
  .int_grad_lin_kappa_cache <<- .int_irmc_kappa_setup_impl("linearized", n_sim, seed)
  .int_grad_lin_kappa_cache
}

# ---- linCmt setup ------------------------------------------------------------

.int_lincmt_setup <- function(n_sim = 500L, seed = 42L) {
  if (!is.null(.int_lincmt_cache) &&
      .int_lincmt_cache$n_sim == n_sim && .int_lincmt_cache$seed == seed)
    return(.int_lincmt_cache)

  skip_if_not_installed("rxode2")

  ui_lc <- suppressMessages(tryCatch(
    rxode2::rxode2(one_cmt_lincmt_fn),
    error = function(e) NULL
  ))
  if (is.null(ui_lc)) skip("rxode2 linCmt model parse failed")

  pinfo      <- admixr2:::.admParseIniDf(ui_lc$iniDf, ui_lc)
  output_var <- "rx_pred_"

  sensModel <- tryCatch(admixr2:::.admLoadSensModel(ui_lc), error = function(e) NULL)
  rxMod     <- tryCatch(admixr2:::.admLoadModel(ui_lc),     error = function(e) NULL)
  if (is.null(rxMod)) skip("linCmt model compilation failed")

  times  <- c(0.5, 1, 2, 4)
  E_true <- .one_cmt_mean(5, 20, 100, times)
  V_true <- diag((0.3 * E_true)^2)

  study <- list(E = E_true, V = V_true, n = 200L, times = times,
                ev = rxode2::et(amt = 100))
  study <- admixr2:::.admNormaliseStudy(study, "s")
  study$ev_full <- study$ev |> rxode2::et(study$times)
  studies <- list(s = study)

  z_list      <- admixr2:::.admMakeZ(n_sim, pinfo, 1L, "sobol")
  params_list <- admixr2:::.admMakeParamsList(n_sim, pinfo, 1L)
  vec         <- admixr2:::.admBuildOptVec(pinfo)
  p0          <- vec$p0

  nll_p0 <- admixr2:::.admNLL(p0, pinfo, studies, z_list, rxMod, output_var, params_list, 1L)

  .int_lincmt_cache <<- list(
    ui = ui_lc, pinfo = pinfo, rxMod = rxMod, sensModel = sensModel,
    studies = studies, z_list = z_list,
    params_list = params_list, vec = vec,
    output_var = output_var, times = times,
    E_true = E_true, n_sim = n_sim, seed = seed,
    nll_p0 = nll_p0
  )
  .int_lincmt_cache
}

# ---- NLL setup ---------------------------------------------------------------

.int_setup <- function(n_sim = 300L, seed = 7L) {
  if (!is.null(.int_cache) &&
      .int_cache$n_sim == n_sim && .int_cache$seed == seed)
    return(.int_cache)

  skip_if_not_installed("rxode2")

  ui <- suppressMessages(tryCatch(
    rxode2::rxode2(one_cmt_fn),
    error = function(e) NULL
  ))
  if (is.null(ui)) skip("rxode2 model parse failed")

  iniDf <- ui$iniDf
  if (!all(c("neta1", "err", "fix") %in% names(iniDf)))
    skip("iniDf format unexpected -- rxode2 API may have changed")

  pinfo      <- admixr2:::.admParseIniDf(iniDf)
  output_var <- "cp"

  rxMod <- tryCatch(admixr2:::.admLoadModel(ui), error = function(e) NULL)
  if (is.null(rxMod)) skip("Model compilation failed")

  times  <- c(0.5, 1, 2, 4)
  E_true <- .one_cmt_mean(5, 20, 100, times)
  V_true <- diag((0.3 * E_true)^2)

  s <- list(E = E_true, V = V_true, n = 200L, times = times,
            ev = rxode2::et(amt = 100))
  s <- admixr2:::.admNormaliseStudy(s, "s")
  s$ev_full <- s$ev |> rxode2::et(s$times)
  studies <- list(s = s)

  z_list      <- admixr2:::.admMakeZ(n_sim, pinfo, 1L, "sobol")
  params_list <- admixr2:::.admMakeParamsList(n_sim, pinfo, 1L)
  vec         <- admixr2:::.admBuildOptVec(pinfo)

  p0      <- vec$p0
  p_bad   <- p0; p_bad["tcl"] <- p_bad["tcl"] + 1.0
  p_nonpd <- p0; p_nonpd[pinfo$omega_par_names[pinfo$chol_diag][1]] <- -1e10

  nll <- function(p)
    admixr2:::.admNLL(p, pinfo, studies, z_list, rxMod, output_var, params_list, 1L)

  nll_warnings <- character(0)
  nll_p0 <- withCallingHandlers(nll(p0), warning = function(w) {
    nll_warnings <<- c(nll_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  nll_p_bad   <- nll(p_bad)
  nll_p_nonpd <- nll(p_nonpd)

  .int_cache <<- list(
    ui = ui, pinfo = pinfo, rxMod = rxMod,
    studies = studies,
    z_list = z_list, params_list = params_list,
    vec = vec, output_var = output_var,
    times = times, E_true = E_true, V_true = V_true,
    n_sim = n_sim, seed = seed,
    p_bad = p_bad, p_nonpd = p_nonpd,
    nll_p0 = nll_p0, nll_p_bad = nll_p_bad, nll_p_nonpd = nll_p_nonpd,
    nll_warnings = nll_warnings
  )
  .int_cache
}

# ---- Plot integration setup --------------------------------------------------
# Builds a real admFit-like env from .int_grad_setup() components:
# real rxMod + real iniDf + real parameters. Avoids calling nlmixr2().
# The "mean" panel will run real rxSolve; "nll"/"par" traces are representative.

.int_plot_setup <- function() {
  if (!is.null(.int_plot_cache)) return(.int_plot_cache)

  env <- .int_grad_setup()   # handles skip_if_not_installed("rxode2") internally
  if (is.null(env$rxMod)) skip("rxMod unavailable from grad setup")

  tryCatch(rxode2::rxLoad(env$rxMod), error = function(e) NULL)

  pars      <- admixr2:::.admUnpack(env$vec$p0, env$pinfo)
  sigma_var <- setNames(pars$sigma_var, env$pinfo$sigma_names)

  n_trace <- 5L
  all_traces <- list(list(
    restart_id = 1L,
    nll_trace  = seq(300, 250, length.out = n_trace),
    par_trace  = matrix(rep(env$vec$p0, n_trace), nrow = n_trace, byrow = TRUE)
  ))

  fit_env <- new.env(parent = emptyenv())
  fit_env$adirmcExtra <- list(
    studies        = env$studies,
    n_sim          = 50L,
    omega          = pars$omega,
    L              = pars$L,
    sigma_var      = sigma_var,
    sigma_is_prop  = env$pinfo$sigma_is_prop,
    sigma_is_lnorm = env$pinfo$sigma_is_lnorm,
    eta_col_names  = env$pinfo$eta_col_names,
    struct         = pars$struct,
    sampling       = "sobol",
    all_traces     = all_traces,
    par_names      = names(env$vec$p0)
  )
  fit_env$ui <- list(
    simulationModel = env$rxMod,
    iniDf           = env$ui$iniDf
  )

  fit <- structure(list(env = fit_env), class = c("admFit", "list"))
  .int_plot_cache <<- list(fit = fit)
  .int_plot_cache
}

# ---- Covariance setup --------------------------------------------------------
# Uses 5000 Sobol points pre-built here (not inside .admCalcCov()) so the draw
# happens at a fixed position in the Sobol sequence. Passes z_cov directly with
# cov_n_sim = NULL to skip internal regeneration which would draw at an
# unpredictable sequence position mid-computation.
# Cached once per session — avoids 5x repeated expensive computation in tests.

.int_cov_setup <- function() {
  if (!is.null(.int_cov_cache)) return(.int_cov_cache)

  env <- .int_grad_setup()   # handles skip_if_not_installed("rxode2") internally

  # At the true params (add.err_sd = 0.1), sigma contributes only 0.01 variance
  # while IIV dominates (~1.7). The Hessian w.r.t. sigma is near-zero → non-PD
  # regardless of step size or n_sim. Shift sigma to sigma_sd = 1 (log_sigma_var = 0)
  # so it contributes ~1.0 variance (comparable to IIV). The result is no longer
  # at the ML optimum but the Hessian is well-conditioned for structural tests.
  pinfo   <- env$pinfo
  n_s     <- length(pinfo$struct_names)
  n_e     <- length(pinfo$sigma_names)
  p_cov   <- env$vec$p0
  p_cov[n_s + seq_len(n_e)] <- 0  # 2*log(sigma_sd) = 0 → sigma_sd = 1

  # Build a 5000-point Sobol z_list here, outside .admCalcCov(), so the draw
  # happens at a fixed position in the Sobol sequence (always right after
  # .int_grad_setup()'s 500-point draw). Passing z_cov directly with
  # cov_n_sim = NULL avoids the internal regeneration inside .admCalcCov()
  # which would draw at an unpredictable sequence position.
  z_cov      <- admixr2:::.admMakeZ(5000L, pinfo, 1L, "sobol")
  params_cov <- admixr2:::.admMakeParamsList(5000L, pinfo, 1L)

  result_nll <- suppressWarnings(admixr2:::.admCalcCov(
    p_cov, pinfo, env$studies, z_cov,
    env$rxMod, env$output_var, params_cov, 1L,
    cov_n_sim = NULL, use_grad = FALSE
  ))

  result_grad <- suppressWarnings(admixr2:::.admCalcCov(
    p_cov, pinfo, env$studies, z_cov,
    env$rxMod, env$output_var, params_cov, 1L,
    cov_n_sim = NULL, use_grad = TRUE, sensModel = NULL
  ))

  .int_cov_cache <<- list(
    env         = env,
    p_cov       = p_cov,
    n_struct    = n_s,
    struct_names = pinfo$struct_names,
    result_nll  = result_nll,
    result_grad = result_grad,
    z_cov       = z_cov,
    params_cov  = params_cov
  )
  .int_cov_cache
}

# ---- FO setup ----------------------------------------------------------------
# Reuses sens+rxMod from .int_grad_setup() (ordering invariant satisfied).
# Kept in helper so it is visible to test-integration-cov.R as well.

.int_adfo_setup <- function() {
  if (!is.null(.int_adfo_cache)) return(.int_adfo_cache)

  skip_on_cran()
  skip_if_not_installed("rxode2")

  env <- .int_grad_setup()
  if (is.null(env$rxMod)) skip("rxMod unavailable from grad setup")

  pinfo      <- env$pinfo
  studies    <- env$studies
  output_var <- env$output_var
  p0         <- env$vec$p0

  params_list <- admixr2:::.admMakeParamsList(1L, pinfo, length(studies))

  nll_p0 <- admixr2:::.adfoNLL(p0, pinfo, studies, env$sensModel, env$rxMod,
                                output_var, params_list, cores = 1L)

  p_bad     <- p0; p_bad["tcl"] <- p_bad["tcl"] + 0.5
  nll_p_bad <- admixr2:::.adfoNLL(p_bad, pinfo, studies, env$sensModel, env$rxMod,
                                   output_var, params_list, cores = 1L)

  h_fd  <- 1e-4
  g_ana <- admixr2:::.adfoGrad(p0, pinfo, studies, env$sensModel, env$rxMod,
                                output_var, params_list, cores = 1L, grad_h = h_fd)

  g_fd <- vapply(seq_along(p0), function(k) {
    ph <- p0; ph[k] <- ph[k] + h_fd
    pl <- p0; pl[k] <- pl[k] - h_fd
    nh <- admixr2:::.adfoNLL(ph, pinfo, studies, env$sensModel, env$rxMod,
                              output_var, params_list, 1L)
    nl <- admixr2:::.adfoNLL(pl, pinfo, studies, env$sensModel, env$rxMod,
                              output_var, params_list, 1L)
    (nh - nl) / (2 * h_fd)
  }, double(1))
  names(g_fd) <- names(p0)

  .int_adfo_cache <<- list(
    pinfo = pinfo, studies = studies,
    rxMod = env$rxMod, sensModel = env$sensModel,
    output_var = output_var, params_list = params_list,
    p0 = p0, p_bad = p_bad,
    nll_p0 = nll_p0, nll_p_bad = nll_p_bad,
    g_ana = g_ana, g_fd = g_fd, h_fd = h_fd
  )
  .int_adfo_cache
}

# ---- FO kappa setup ----------------------------------------------------------
# Exercises the unpaired-struct-theta gradient path in .adfoGrad().
# tsc is intentionally not mu-referenced by any eta.

.int_adfo_kappa_setup <- function() {
  if (!is.null(.int_adfo_kappa_cache)) return(.int_adfo_kappa_cache)

  skip_on_cran()
  skip_if_not_installed("rxode2")

  ui <- suppressMessages(tryCatch(
    rxode2::rxode2(one_cmt_kappa_fn),
    error = function(e) NULL
  ))
  if (is.null(ui)) skip("rxode2 model parse failed")

  pinfo      <- admixr2:::.admParseIniDf(ui$iniDf, ui)
  output_var <- "cp"

  sensModel <- tryCatch(admixr2:::.admLoadSensModel(ui), error = function(e) NULL)
  rxMod     <- tryCatch(admixr2:::.admLoadModel(ui),     error = function(e) NULL)
  if (is.null(rxMod))     skip("Model compilation failed")
  if (is.null(sensModel)) skip("Sensitivity model unavailable")

  times  <- c(0.5, 1, 2, 4)
  E_true <- .one_cmt_mean(5, 20, 100, times)
  V_true <- diag((0.3 * E_true)^2)

  study <- list(E = E_true, V = V_true, n = 200L, times = times,
                ev = rxode2::et(amt = 100))
  study <- admixr2:::.admNormaliseStudy(study, "s")
  study$ev_full <- study$ev |> rxode2::et(study$times)
  studies <- list(s = study)

  params_list <- admixr2:::.admMakeParamsList(1L, pinfo, length(studies))
  p0          <- admixr2:::.admBuildOptVec(pinfo)$p0
  h_fd        <- 1e-4

  g_ana <- admixr2:::.adfoGrad(p0, pinfo, studies, sensModel, rxMod,
                               output_var, params_list, cores = 1L, grad_h = h_fd)

  # Central FD reference using the same scaled step as .adfoGrad() so that
  # the perturbation matches for each parameter (pmax(|p[k]|, 0.1) * h_fd).
  g_fd <- vapply(seq_along(p0), function(k) {
    hk <- pmax(abs(p0[k]), 0.1) * h_fd
    ph <- p0; ph[k] <- ph[k] + hk
    pl <- p0; pl[k] <- pl[k] - hk
    nh <- admixr2:::.adfoNLL(ph, pinfo, studies, sensModel, rxMod,
                              output_var, params_list, 1L)
    nl <- admixr2:::.adfoNLL(pl, pinfo, studies, sensModel, rxMod,
                              output_var, params_list, 1L)
    (nh - nl) / (2 * hk)
  }, double(1))
  names(g_fd) <- names(p0)

  .int_adfo_kappa_cache <<- list(
    pinfo = pinfo, studies = studies,
    rxMod = rxMod, sensModel = sensModel,
    output_var = output_var, params_list = params_list,
    p0 = p0, g_ana = g_ana, g_fd = g_fd, h_fd = h_fd
  )
  .int_adfo_kappa_cache
}

# ---- FO covariance setup -----------------------------------------------------
# Mirrors .int_cov_setup() but calls .adfoCalcCov() instead of .admCalcCov().
# Shifts sigma to sd=1 (log_sigma_var=0) so the Hessian is well-conditioned.

.int_adfo_cov_setup <- function() {
  if (!is.null(.int_adfo_cov_cache)) return(.int_adfo_cov_cache)

  env <- .int_adfo_setup()   # handles skip guards internally

  pinfo <- env$pinfo
  n_s   <- length(pinfo$struct_names)
  n_e   <- length(pinfo$sigma_names)
  p_cov <- env$p0
  p_cov[n_s + seq_len(n_e)] <- 0  # 2*log(sigma_sd) = 0 → sigma_sd = 1

  result_nll <- suppressWarnings(admixr2:::.adfoCalcCov(
    p_cov, pinfo, env$studies, env$sensModel, env$rxMod, env$output_var,
    env$params_list, 1L,
    use_grad = FALSE
  ))

  result_grad <- suppressWarnings(admixr2:::.adfoCalcCov(
    p_cov, pinfo, env$studies, env$sensModel, env$rxMod, env$output_var,
    env$params_list, 1L,
    use_grad = TRUE
  ))

  .int_adfo_cov_cache <<- list(
    env          = env,
    p_cov        = p_cov,
    n_struct     = n_s,
    struct_names = pinfo$struct_names,
    result_nll   = result_nll,
    result_grad  = result_grad
  )
  .int_adfo_cov_cache
}

# ---- IRMC exact-kappa gradient setup -----------------------------------------
# Same as .int_grad_lin_kappa_setup() but kappa_method = "exact".
# Uses one_cmt_kappa_fn so has_kappa = TRUE (tsc is unpaired).

.int_irmc_exact_kappa_setup <- function(n_sim = 500L, seed = 42L) {
  if (!is.null(.int_irmc_exact_kappa_cache) &&
      .int_irmc_exact_kappa_cache$n_sim == n_sim &&
      .int_irmc_exact_kappa_cache$seed == seed)
    return(.int_irmc_exact_kappa_cache)
  .int_irmc_exact_kappa_cache <<- .int_irmc_kappa_setup_impl("exact", n_sim, seed)
  .int_irmc_exact_kappa_cache
}
