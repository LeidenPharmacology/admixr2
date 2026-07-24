# Tier-2 integration: multi-compartment fitting (multiple observed outputs).
# Fits a 2-compartment model that observes BOTH plasma (cp) and tissue (ct) --
# each observed compartment is an independent likelihood block with its own ev,
# times, E and V. Verifies adfo/admc/adgh recover near-truth estimates, that
# adirmc reports the not-yet-supported guard, and that plot() produces one panel
# set per observed compartment.

skip_on_cran()
skip_if_not_installed("rxode2")

# ---- Two-output model + deterministic aggregate-data generator ---------------

.mo_model <- function() {
  ini({
    tcl <- log(1.0); tv1 <- log(10); tq <- log(2.0); tv2 <- log(20)
    prop.cp <- 0.10
    add.ct  <- 0.05
    eta.cl ~ 0.09
    eta.v1 ~ 0.04
  })
  model({
    cl <- exp(tcl + eta.cl); v1 <- exp(tv1 + eta.v1)
    q  <- exp(tq);           v2 <- exp(tv2)
    d/dt(central) <- -(cl/v1)*central - (q/v1)*central + (q/v2)*periph
    d/dt(periph)  <-  (q/v1)*central - (q/v2)*periph
    cp <- central / v1
    ct <- periph  / v2
    cp ~ prop(prop.cp)
    ct ~ add(add.ct)
  })
}

# Build studies + fits once; keep references alive for the whole file so rxode2
# model DLLs are not GC-unloaded mid-run (Windows finalizer safety).
.mo_cache <- NULL
.mo_setup <- function() {
  if (!is.null(.mo_cache)) return(.mo_cache)

  set.seed(7)
  Omega <- diag(c(0.09, 0.04))
  simstudy <- function(outcol, times, n, resid_sd) {
    sim <- rxode2::rxode2({
      cl <- 1.0 * exp(eta.cl); v1 <- 10 * exp(eta.v1); q <- 2.0; v2 <- 20
      d/dt(central) <- -(cl/v1)*central - (q/v1)*central + (q/v2)*periph
      d/dt(periph)  <-  (q/v1)*central - (q/v2)*periph
      cp <- central / v1; ct <- periph / v2
    })
    ev  <- rxode2::et(amt = 100, cmt = "central") |> rxode2::et(times)
    df  <- as.data.frame(MASS::mvrnorm(n, c(0, 0), Omega))
    names(df) <- c("eta.cl", "eta.v1")
    out <- rxode2::rxSolve(sim, params = df, events = ev, returnType = "data.frame")
    M   <- matrix(out[[outcol]][out$time %in% times], nrow = n, byrow = TRUE)
    M   <- M + matrix(rnorm(length(M), 0, resid_sd), nrow = n)
    list(E = colMeans(M), V = cov.wt(M, method = "ML")$cov, n = n, times = times)
  }
  ev <- rxode2::et(amt = 100, cmt = "central")
  pl <- simstudy("cp", c(0.5, 1, 2, 4, 8), 60L, 0.1)
  ct <- simstudy("ct", c(1, 4, 8),          20L, 0.05)
  study <- list(observations = list(
    plasma = list(output = "cp", ev = ev, times = pl$times, E = pl$E, V = pl$V, n = pl$n),
    tissue = list(output = "ct", ev = ev, times = ct$times, E = ct$E, V = ct$V, n = ct$n)))

  dat <- admData(c("cp", "ct"))
  # Exercise the analytical / sensitivity gradient paths (not FD fallback).
  fo <- suppressMessages(nlmixr2est::nlmixr2(.mo_model, dat, est = "adfo",
    control = adfoControl(studies = list(rat = study), grad = "analytical",
                          maxeval = 200L, covMethod = "none")))
  mc <- suppressMessages(nlmixr2est::nlmixr2(.mo_model, dat, est = "admc",
    control = admControl(studies = list(rat = study), n_sim = 400L, maxeval = 80L,
                         grad = "sens", covMethod = "none", seed = 1L)))
  gh <- suppressMessages(nlmixr2est::nlmixr2(.mo_model, dat, est = "adgh",
    control = adghControl(studies = list(rat = study), n_nodes = 5L, maxeval = 200L,
                          grad = "analytical", covMethod = "none")))

  .mo_cache <<- list(study = study, fo = fo, mc = mc, gh = gh)
  .mo_cache
}

truth <- c(tcl = 0, tv1 = log(10), tq = log(2), tv2 = log(20))

test_that("adfo fits a multi-compartment (2 observed outputs) model near truth", {
  s <- .mo_setup()
  expect_true(is.finite(s$fo$objective))
  expect_s3_class(s$fo, "admFit")
  est <- unlist(s$fo$env$admExtra$struct)
  expect_equal(est[["tcl"]], truth[["tcl"]], tolerance = 0.5)
  expect_equal(est[["tv1"]], truth[["tv1"]], tolerance = 0.3)
})

test_that("admc fits a multi-compartment model near truth (sens gradient path)", {
  s <- .mo_setup()
  expect_true(is.finite(s$mc$objective))
  est <- unlist(s$mc$env$admExtra$struct)
  expect_equal(est[["tv1"]], truth[["tv1"]], tolerance = 0.3)
  expect_equal(est[["tv2"]], truth[["tv2"]], tolerance = 0.4)
})

test_that("multi-output analytical/sens gradients match finite differences", {
  s <- .mo_setup()
  ui    <- rxode2::rxode2(.mo_model)
  pinfo <- admixr2:::.admParseIniDf(ui$iniDf, ui)
  ov    <- admixr2:::.admOutputVar(ui)
  st    <- admixr2:::.admNormaliseStudy(s$study, "rat", ov)
  units <- admixr2:::.admBuildEvFull(admixr2:::.admFlattenStudies(list(rat = st)),
                                     tag_cmt = TRUE)
  rxMod <- admixr2:::.admLoadModel(ui)
  sm    <- admixr2:::.admLoadSensModel(ui)
  p     <- admixr2:::.admBuildOptVec(pinfo)$p0 +
           c(0.1, -0.1, 0.05, -0.05, 0.2, 0.15, 0.1, -0.1)

  # adfo analytical vs central FD (deterministic)
  pl1  <- admixr2:::.admMakeParamsList(1L, pinfo, length(units))
  g_an <- admixr2:::.adfoGrad(p, pinfo, units, sm, rxMod, ov, pl1, 1L, 1e-4)
  h    <- 3e-5
  g_fd <- vapply(seq_along(p), function(k) {
    pp <- p; pp[k] <- p[k] + h; pm <- p; pm[k] <- p[k] - h
    (admixr2:::.adfoNLL(pp, pinfo, units, sm, rxMod, ov, pl1, 1L) -
     admixr2:::.adfoNLL(pm, pinfo, units, sm, rxMod, ov, pl1, 1L)) / (2 * h)
  }, double(1))
  expect_equal(unname(g_an), unname(g_fd), tolerance = 0.02)
})

test_that("adgh fits a multi-compartment model near truth", {
  s <- .mo_setup()
  expect_true(is.finite(s$gh$objective))
  est <- unlist(s$gh$env$admExtra$struct)
  expect_equal(est[["tv1"]], truth[["tv1"]], tolerance = 0.3)
})

test_that("adirmc reports the multi-output guard", {
  s <- .mo_setup()
  expect_error(
    suppressMessages(nlmixr2est::nlmixr2(.mo_model, admData(c("cp", "ct")), est = "adirmc",
      control = adirmcControl(studies = list(rat = s$study)))),
    regexp = "does not yet support multiple observed outputs")
})

test_that("plot() yields one mean/cov panel set per observed compartment", {
  s <- .mo_setup()
  skip_if_not_installed("ggplot2")
  ps <- plot(s$fo, which = c("mean", "cov"))
  expect_true(all(c("mean_rat.plasma", "mean_rat.tissue",
                    "cov_rat.plasma",  "cov_rat.tissue") %in% names(ps)))
})

test_that("plot() renders all four panel types for a multi-output fit", {
  s <- .mo_setup()
  skip_if_not_installed("ggplot2")
  ps <- plot(s$mc, which = c("mean", "cov", "nll", "par"))
  expect_true(all(c("mean_rat.plasma", "mean_rat.tissue",
                    "cov_rat.plasma",  "cov_rat.tissue",
                    "nll_trace", "par_trace") %in% names(ps)))
  # every panel is a ggplot object
  expect_true(all(vapply(ps, function(p) inherits(p, "ggplot"), logical(1))))
})

test_that("datagen produces per-compartment aggregate data that round-trips", {
  skip_if_not_installed("MASS")
  ev <- rxode2::et(amt = 100, cmt = "central")
  dstudies <- list(rat = list(n = 50L, observations = list(
    plasma = list(output = "cp", times = c(0.5, 1, 2, 4, 8), ev = ev),
    tissue = list(output = "ct", times = c(1, 4, 8),          ev = ev))))
  d <- datagen(dstudies, model = .mo_model,
               control = datagenControl(method = "fo"))
  obs <- d$rat$observations
  expect_named(obs, c("plasma", "tissue"))
  expect_equal(obs$plasma$output, "cp")
  expect_equal(obs$tissue$output, "ct")
  expect_length(obs$plasma$E, 5L)
  expect_equal(dim(obs$tissue$V), c(3L, 3L))

  # Round-trip: fit the generated data and recover truth (FO data is exact).
  study_fit <- list(observations = list(
    plasma = list(output = "cp", ev = ev, times = obs$plasma$times,
                  E = obs$plasma$E, V = obs$plasma$V, n = 50L),
    tissue = list(output = "ct", ev = ev, times = obs$tissue$times,
                  E = obs$tissue$E, V = obs$tissue$V, n = 50L)))
  f <- suppressMessages(nlmixr2est::nlmixr2(.mo_model, admData(c("cp", "ct")),
    est = "adfo", control = adfoControl(studies = list(rat = study_fit),
                                        grad = "analytical", maxeval = 200L,
                                        covMethod = "none")))
  est <- unlist(f$env$admExtra$struct)
  expect_equal(est[["tcl"]], truth[["tcl"]], tolerance = 0.05)
  expect_equal(est[["tv1"]], truth[["tv1"]], tolerance = 0.05)
  expect_equal(est[["tv2"]], truth[["tv2"]], tolerance = 0.05)
})

test_that("objective equals the standalone aggregate -2LL at the estimate", {
  s <- .mo_setup()
  ui    <- rxode2::rxode2(.mo_model)
  pinfo <- admixr2:::.admParseIniDf(ui$iniDf, ui)
  ov    <- admixr2:::.admOutputVar(ui)
  st    <- admixr2:::.admNormaliseStudy(s$study, "rat", ov)
  units <- admixr2:::.admBuildEvFull(admixr2:::.admFlattenStudies(list(rat = st)),
                                     tag_cmt = TRUE)
  sm    <- admixr2:::.admLoadSensModel(ui)   # match the fit's analytical Jacobian
  rxMod <- admixr2:::.admLoadModel(ui)       # ordering invariant: sens before model
  pl    <- admixr2:::.admMakeParamsList(1L, pinfo, length(units))
  p_hat <- s$fo$env$admExtra$nloptr$solution
  nll   <- admixr2:::.adfoNLL(p_hat, pinfo, units, sm, rxMod, ov, pl, 1L)
  expect_equal(nll, s$fo$objective, tolerance = 1e-4)
})

# ---- Joint (same-subject) fits: cross-covariance across compartments ----------

.mo_joint_cache <- NULL
.mo_joint_setup <- function() {
  if (!is.null(.mo_joint_cache)) return(.mo_joint_cache)
  skip_if_not_installed("MASS")
  set.seed(11)
  Omega <- diag(c(0.09, 0.04))
  n <- 80L; tp <- c(0.5, 1, 2, 4, 8); tc <- c(1, 4, 8)
  sim <- rxode2::rxode2({
    cl <- 1.0*exp(eta.cl); v1 <- 10*exp(eta.v1); q <- 2.0; v2 <- 20
    d/dt(central) <- -(cl/v1)*central - (q/v1)*central + (q/v2)*periph
    d/dt(periph)  <-  (q/v1)*central - (q/v2)*periph
    cp <- central/v1; ct <- periph/v2 })
  ut <- sort(unique(c(tp, tc)))
  ev <- rxode2::et(amt = 100, cmt = "central") |> rxode2::et(ut)
  df <- as.data.frame(MASS::mvrnorm(n, c(0, 0), Omega)); names(df) <- c("eta.cl", "eta.v1")
  out <- rxode2::rxSolve(sim, params = df, events = ev, returnType = "data.frame")
  CP <- matrix(out$cp[out$time %in% ut], nrow = n, byrow = TRUE); colnames(CP) <- ut
  CT <- matrix(out$ct[out$time %in% ut], nrow = n, byrow = TRUE); colnames(CT) <- ut
  M  <- cbind(CP[, as.character(tp)] + matrix(rnorm(n*length(tp), 0, 0.1), n),
              CT[, as.character(tc)] + matrix(rnorm(n*length(tc), 0, 0.05), n))
  Ej <- colMeans(M); Vj <- cov.wt(M, method = "ML")$cov
  ip <- seq_along(tp); ic <- length(tp) + seq_along(tc)
  evd <- rxode2::et(amt = 100, cmt = "central")
  study_cross <- list(n = n, ev = evd, observations = list(
    plasma = list(output = "cp", times = tp, E = Ej[ip], V = Vj[ip, ip]),
    brain  = list(output = "ct", times = tc, E = Ej[ic], V = Vj[ic, ic])),
    cross = list("plasma:brain" = Vj[ip, ic]))
  study_fullV <- list(n = n, ev = evd, observations = list(
    plasma = list(output = "cp", times = tp, E = Ej[ip]),
    brain  = list(output = "ct", times = tc, E = Ej[ic])),
    V = Vj)
  fmc <- suppressMessages(nlmixr2est::nlmixr2(.mo_model, admData(c("cp","ct")), est = "admc",
    control = admControl(studies = list(rat = study_cross), n_sim = 400L, maxeval = 80L,
                         grad = "sens", covMethod = "none", seed = 1L)))
  .mo_joint_cache <<- list(cross = study_cross, fullV = study_fullV,
                           Vj = Vj, ip = ip, ic = ic, fmc = fmc)
  .mo_joint_cache
}

test_that("joint (same-subject) admc fit uses cross-covariance and recovers truth", {
  s <- .mo_joint_setup()
  expect_true(is.finite(s$fmc$objective))
  # cross-covariance is non-trivial (else 'joint' would be pointless here)
  expect_gt(max(abs(s$Vj[s$ip, s$ic])), 0.05)
  est <- unlist(s$fmc$env$admExtra$struct)
  expect_equal(est[["tv1"]], truth[["tv1"]], tolerance = 0.3)
  expect_equal(est[["tv2"]], truth[["tv2"]], tolerance = 0.4)
})

test_that("cross-pairs and full-V joint specs give the same fit", {
  s <- .mo_joint_setup()
  f_full <- suppressMessages(nlmixr2est::nlmixr2(.mo_model, admData(c("cp","ct")), est = "adfo",
    control = adfoControl(studies = list(rat = s$fullV), grad = "analytical",
                          maxeval = 150L, covMethod = "none")))
  f_cross <- suppressMessages(nlmixr2est::nlmixr2(.mo_model, admData(c("cp","ct")), est = "adfo",
    control = adfoControl(studies = list(rat = s$cross), grad = "analytical",
                          maxeval = 150L, covMethod = "none")))
  expect_equal(f_full$objective, f_cross$objective, tolerance = 1e-6)
})

test_that("adfo/adgh fit the joint model and joint covMethod='r' works", {
  s <- .mo_joint_setup()
  f_fo <- suppressMessages(nlmixr2est::nlmixr2(.mo_model, admData(c("cp","ct")), est = "adfo",
    control = adfoControl(studies = list(rat = s$cross), grad = "analytical",
                          maxeval = 150L, covMethod = "r")))
  expect_true(is.finite(f_fo$objective))
  expect_false(is.null(f_fo$env$cov))                 # NLL-FD Hessian on a joint fit
  f_gh <- suppressMessages(nlmixr2est::nlmixr2(.mo_model, admData(c("cp","ct")), est = "adgh",
    control = adghControl(studies = list(rat = s$cross), n_nodes = 5L,
                          grad = "analytical", maxeval = 150L, covMethod = "none")))
  expect_true(is.finite(f_gh$objective))
})

test_that("adirmc rejects a joint (same-subject) study", {
  s <- .mo_joint_setup()
  expect_error(
    suppressMessages(nlmixr2est::nlmixr2(.mo_model, admData(c("cp","ct")), est = "adirmc",
      control = adirmcControl(studies = list(rat = s$cross)))),
    regexp = "does not yet support multiple observed outputs")
})

test_that("the multi-output sens model carries theta columns for tq/tv2", {
  # .mo_model has two thetas with no eta (tq, tv2), so its sens model is augmented
  # with a dummy eta for each and the gradients below take the THETA-COLUMN path,
  # not the FD one. Asserted explicitly: without this, a regression that silently
  # dropped the theta columns would still pass the gradient tests (they would just
  # quietly fall back to FD, which is also correct -- only slower and coarser).
  ui <- rxode2::rxode2(.mo_model)
  expect_equal(admixr2:::.admUnpairedThetas(ui), c("tq", "tv2"))
  sm <- admixr2:::.admLoadSensModel(ui)
  expect_false(is.null(sm$theta_sens_cols))
  expect_equal(names(sm$theta_sens_cols), c("tq", "tv2"))
  # tq/tv2 are THETA[3]/THETA[4]; the two mu-referenced thetas reuse their etas
  expect_equal(sm$dirs, c("ETA_1_", "ETA_2_", "THETA_3_", "THETA_4_"))
})

test_that("joint + multi-output sens solves return d(pred)/d(theta) for every block", {
  # The joint path stacks each output's rows; a block that failed to return its
  # theta columns must disable the theta path for the WHOLE unit (a half-filled
  # stacked derivative would be silently wrong), so check the stacked shape.
  s     <- .mo_joint_setup()
  ui    <- rxode2::rxode2(.mo_model)
  pinfo <- admixr2:::.admParseIniDf(ui$iniDf, ui)
  ov    <- admixr2:::.admOutputVar(ui)
  st    <- admixr2:::.admNormaliseStudy(s$cross, "rat", ov)
  units <- admixr2:::.admBuildEvFull(admixr2:::.admFlattenStudies(list(rat = st)),
                                     tag_cmt = TRUE)
  sm    <- admixr2:::.admLoadSensModel(ui)
  admixr2:::.admLoadModel(ui)
  unit  <- units[[1]]
  expect_true(isTRUE(unit$is_joint))

  pars    <- admixr2:::.admUnpack(admixr2:::.admBuildOptVec(pinfo)$p0, pinfo)
  n_sim   <- 16L
  z       <- admixr2:::.admMakeZ(n_sim, pinfo, 1L, "sobol")[[1]]
  eta_mat <- z %*% t(pars$L)
  colnames(eta_mat) <- pinfo$eta_col_names

  js <- admixr2:::.admSimulateJointSens(sm, pars$struct, pinfo$sigma_names,
                                        eta_mat, unit, 1L)
  expect_false(is.null(js$dtheta_list))
  expect_equal(names(js$dtheta_list), c("tq", "tv2"))
  for (nm in c("tq", "tv2")) {
    # one column per stacked row across BOTH outputs, all filled (no zero block)
    expect_equal(dim(js$dtheta_list[[nm]]), c(n_sim, unit$n_total))
    expect_true(all(is.finite(js$dtheta_list[[nm]])))
    for (blk in unit$blocks)
      expect_false(all(js$dtheta_list[[nm]][, blk$rows] == 0),
                   info = paste(nm, blk$output))
  }
})

test_that("analytical joint (stacked-MVN) gradient matches finite differences", {
  s     <- .mo_joint_setup()
  ui    <- rxode2::rxode2(.mo_model)
  pinfo <- admixr2:::.admParseIniDf(ui$iniDf, ui)
  ov    <- admixr2:::.admOutputVar(ui)
  st    <- admixr2:::.admNormaliseStudy(s$cross, "rat", ov)
  units <- admixr2:::.admBuildEvFull(admixr2:::.admFlattenStudies(list(rat = st)),
                                     tag_cmt = TRUE)
  sm    <- admixr2:::.admLoadSensModel(ui)
  rxMod <- admixr2:::.admLoadModel(ui)
  set.seed(42); z_list <- admixr2:::.admMakeZ(20000L, pinfo, 1L, "sobol")
  pl    <- admixr2:::.admMakeParamsList(20000L, pinfo, 1L)
  p     <- admixr2:::.admBuildOptVec(pinfo)$p0 +
           c(0.1, -0.1, 0.05, -0.05, 0.2, 0.15, 0.1, -0.1)
  g_an <- admixr2:::.admGrad(p, pinfo, units, z_list, rxMod, ov, pl, 1L, 1e-4, sm)
  h    <- 1e-4
  g_fd <- vapply(seq_along(p), function(k) {
    pp <- p; pp[k] <- p[k] + h; pm <- p; pm[k] <- p[k] - h
    (admixr2:::.admNLL(pp, pinfo, units, z_list, rxMod, ov, pl, 1L) -
     admixr2:::.admNLL(pm, pinfo, units, z_list, rxMod, ov, pl, 1L)) / (2 * h)
  }, double(1))
  expect_equal(unname(g_an), unname(g_fd), tolerance = 0.02)
})

test_that("analytical joint gradients (adfo FO + adgh quadrature) match FD", {
  s     <- .mo_joint_setup()
  ui    <- rxode2::rxode2(.mo_model)
  pinfo <- admixr2:::.admParseIniDf(ui$iniDf, ui)
  ov    <- admixr2:::.admOutputVar(ui)
  st    <- admixr2:::.admNormaliseStudy(s$cross, "rat", ov)
  units <- admixr2:::.admBuildEvFull(admixr2:::.admFlattenStudies(list(rat = st)),
                                     tag_cmt = TRUE)
  sm    <- admixr2:::.admLoadSensModel(ui)
  rxMod <- admixr2:::.admLoadModel(ui)
  p     <- admixr2:::.admBuildOptVec(pinfo)$p0 +
           c(0.1, -0.1, 0.05, -0.05, 0.2, 0.15, 0.1, -0.1)
  h     <- 3e-5

  # adfo joint analytical vs central FD of the joint .adfoNLL
  pl1  <- admixr2:::.admMakeParamsList(1L, pinfo, 1L)
  g_fo <- admixr2:::.adfoGrad(p, pinfo, units, sm, rxMod, ov, pl1, 1L, 1e-4)
  fd_fo <- vapply(seq_along(p), function(k) {
    pp <- p; pp[k] <- p[k] + h; pm <- p; pm[k] <- p[k] - h
    (admixr2:::.adfoNLL(pp, pinfo, units, sm, rxMod, ov, pl1, 1L) -
     admixr2:::.adfoNLL(pm, pinfo, units, sm, rxMod, ov, pl1, 1L)) / (2 * h)
  }, double(1))
  expect_equal(unname(g_fo), unname(fd_fo), tolerance = 0.02)

  # adgh joint analytical vs central FD of the joint .adghNLL
  grid <- admixr2:::.adghNodeGrid(7L, pinfo$n_eta)
  g_gh <- admixr2:::.adghGrad(p, pinfo, units, sm, rxMod, ov, grid, 1L, 1e-4)
  fd_gh <- vapply(seq_along(p), function(k) {
    pp <- p; pp[k] <- p[k] + h; pm <- p; pm[k] <- p[k] - h
    (admixr2:::.adghNLL(pp, pinfo, units, rxMod, ov, grid, 1L) -
     admixr2:::.adghNLL(pm, pinfo, units, rxMod, ov, grid, 1L)) / (2 * h)
  }, double(1))
  expect_equal(unname(g_gh), unname(fd_gh), tolerance = 0.02)
})

test_that(".adghGrad degrades to FD on a JOINT unit (singular V, and failed sens solve)", {
  # The joint branch of .adghGrad gained the same two guards as the non-joint one:
  # a failed joint sens solve (`is.null(js)`) and a singular predicted V
  # (`is.null(G)`) each used to `next`, silently dropping the unit from the
  # gradient. Both must now degrade the whole gradient to .adghFDGrad.
  s     <- .mo_joint_setup()
  ui    <- rxode2::rxode2(.mo_model)
  pinfo <- admixr2:::.admParseIniDf(ui$iniDf, ui)
  ov    <- admixr2:::.admOutputVar(ui)
  st    <- admixr2:::.admNormaliseStudy(s$cross, "rat", ov)
  units <- admixr2:::.admBuildEvFull(admixr2:::.admFlattenStudies(list(rat = st)),
                                     tag_cmt = TRUE)
  sm    <- admixr2:::.admLoadSensModel(ui)
  rxMod <- admixr2:::.admLoadModel(ui)
  grid  <- admixr2:::.adghNodeGrid(5L, pinfo$n_eta)

  # (a) singular predicted V: omega and sigma driven to exactly 0
  p_sing <- admixr2:::.admBuildOptVec(pinfo)$p0
  p_sing[grep("^logchol", names(p_sing))] <- -1500
  p_sing[pinfo$sigma_names]               <- -1500
  g_sing <- admixr2:::.adghGrad(p_sing, pinfo, units, sm, rxMod, ov, grid, 1L, 1e-3)
  g_fd_s <- admixr2:::.adghFDGrad(p_sing, pinfo, units, rxMod, ov, grid, 1L, 1e-3)
  expect_equal(unname(g_sing), unname(g_fd_s))

  # (b) failed joint sens solve: a sens model whose columns do not exist forces
  # .admSimulateJointSens to return NULL
  p0     <- admixr2:::.admBuildOptVec(pinfo)$p0
  broken <- sm; broken$sens_cols <- paste0("no_such_", seq_along(broken$sens_cols))
  g_brk  <- admixr2:::.adghGrad(p0, pinfo, units, broken, rxMod, ov, grid, 1L, 1e-3)
  g_fd_b <- admixr2:::.adghFDGrad(p0, pinfo, units, rxMod, ov, grid, 1L, 1e-3)
  expect_true(all(is.finite(g_brk)))
  expect_equal(unname(g_brk), unname(g_fd_b), tolerance = 1e-8)
})

test_that("a joint fit gets per-endpoint aggregate diagnostics", {
  # .admAggData() reconstructs the residual itself rather than reusing what the
  # estimator computed, and it did so with ONE output name for the whole unit --
  # so a joint unit got the first endpoint's residual spec applied to every
  # stacked row, and its single-output solve returned plasma's trajectory for the
  # brain rows. In practice it then died on the dimnames (length(times) is 5 while
  # the stacked V is 8x8), and since .admAttachAggData() guards with tryCatch the
  # slot was simply never set: no diagnostics at all, silently.
  #
  # This model is the case that matters: cp ~ prop() and ct ~ add() are different
  # residual FORMS, so a shared spec cannot be right for both.
  s   <- .mo_joint_setup()
  fit <- s$fmc
  ad  <- fit$env$aggData
  expect_false(is.null(ad))
  e <- ad[[1L]]
  n_tot <- length(s$ip) + length(s$ic)
  expect_length(e$pred$E, n_tot)
  expect_identical(dim(e$pred$V), c(n_tot, n_tot))

  # labelled by the endpoint each row belongs to -- the times repeat across
  # blocks, so bare times would be duplicate labels
  expect_true(all(grepl("^(cp|ct)@", names(e$pred$E))))
  expect_identical(sum(grepl("^cp@", names(e$pred$E))), length(s$ip))
  expect_identical(sum(grepl("^ct@", names(e$pred$E))), length(s$ic))

  # and BOTH blocks are predicted, not just the first: a shared residual spec
  # (or plasma's trajectory in the brain rows) shows up here
  expect_true(all(is.finite(e$pred$V)))
  expect_true(all(diag(e$pred$V) > 0))
  r <- diag(e$pred$V) / diag(e$obs$V)
  expect_true(all(r > 0.4 & r < 2.5), info = paste(round(r, 3), collapse = " "))
  expect_equal(unname(e$pred$E), unname(e$obs$E), tolerance = 0.35)
})
