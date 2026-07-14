# Tier 2: models with a FIXED structural theta.
#
# A fixed theta is not an estimated parameter, so pinfo does not carry it and
# nothing wrote its THETA[k] column -- but the sensitivity model still HAS that
# slot and rxSolve requires every parameter. The sens solve therefore errored and
# returned NULL, which:
#   * dropped admc/adfo to a finite-difference gradient, silently; and
#   * made .adghGrad `next` -- returning a gradient that silently OMITTED the
#     study altogether.
# The map was also indexed by position among the non-fixed thetas, so every theta
# after a fixed one pointed at the wrong THETA[k].

skip_on_cran()
skip_if_not_installed("rxode2")

# tv is FIXED and sits BETWEEN two estimated thetas, so position (2) != ntheta (3)
# for tka -- the case the old position-indexed map got wrong.
fixed_theta_fn <- function() {
  ini({
    tcl     <- log(5)
    tv      <- fix(log(20))
    tka     <- log(1)
    add.err <- 0.1
    eta.cl  ~ 0.09
    eta.ka  ~ 0.04
  })
  model({
    cl <- exp(tcl + eta.cl)
    v  <- exp(tv)
    ka <- exp(tka + eta.ka)
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - (cl / v) * central
    cp <- central / v
    cp ~ add(add.err)
  })
}

.ft_setup <- function(n_sim = 64L) {
  ui    <- suppressMessages(rxode2::rxode2(fixed_theta_fn))
  pinfo <- admixr2:::.admParseIniDf(ui$iniDf, ui)
  sens  <- suppressMessages(admixr2:::.admLoadSensModel(ui))
  if (is.null(sens)) skip("sens model unavailable")
  rxMod <- admixr2:::.admLoadModel(ui)
  rxode2::rxLoad(rxMod)

  times  <- c(0.5, 1, 2, 4, 8)
  E_true <- c(1.8, 2.7, 3.1, 2.3, 0.9)
  study  <- admixr2:::.admNormaliseStudy(
    list(E = E_true, V = diag((0.2 * E_true)^2), n = 100L, times = times,
         ev = rxode2::et(amt = 100)), "s")
  study$ev_full <- study$ev |> rxode2::et(study$times)

  pars    <- admixr2:::.admUnpack(admixr2:::.admBuildOptVec(pinfo)$p0, pinfo)
  z       <- admixr2:::.admMakeZ(n_sim, pinfo, 1L, "sobol")[[1]]
  eta_mat <- z %*% t(pars$L)
  colnames(eta_mat) <- pinfo$eta_col_names

  list(ui = ui, pinfo = pinfo, sens = sens, rxMod = rxMod, study = study,
       studies = list(s = study), pars = pars, eta_mat = eta_mat, n_sim = n_sim,
       plist = admixr2:::.admMakeParamsList(n_sim, pinfo, 1L),
       p0 = admixr2:::.admBuildOptVec(pinfo)$p0)
}

test_that("the sens model carries the fixed theta's value and maps by ntheta", {
  e <- .ft_setup()
  # tv is THETA[2] and must NOT shift tka out of THETA[3]
  expect_equal(unname(e$sens$rename_map[["tcl"]]), "THETA[1]")
  expect_equal(unname(e$sens$rename_map[["tv"]]),  "THETA[2]")
  expect_equal(unname(e$sens$rename_map[["tka"]]), "THETA[3]")

  # the fixed value travels with the sens model so the solve paths can fill it
  expect_equal(e$sens$fixed_theta, c(`THETA[2]` = log(20)))

  # ... and pinfo (correctly) does not carry it
  expect_false("tv" %in% e$pinfo$struct_names)
})

test_that("the sens solve works with a fixed theta (was NULL)", {
  e <- .ft_setup()
  cp_sim <- admixr2:::.admSimulate(e$rxMod, e$pars$struct, e$pinfo$sigma_names,
                                   e$eta_mat, e$study, "cp", e$plist[[1]], 1L)
  res <- admixr2:::.admSimulateSens(e$sens, e$pars$struct, e$pinfo$sigma_names,
                                    e$eta_mat, e$study, 1L)
  expect_false(is.null(res))
  # tolerance 1e-4, not exact: rxode2's adaptive solver takes different steps for
  # a system carrying sensitivity compartments (pre-existing, documented).
  expect_equal(res$cp_mat, cp_sim, tolerance = 1e-4)
  expect_length(res$dpred_list, e$pinfo$n_eta)
})

test_that("admc gradient is analytical (not FD) and matches FD for a fixed-theta model", {
  e <- .ft_setup(n_sim = 500L)
  z_list <- admixr2:::.admMakeZ(500L, e$pinfo, 1L, "sobol")
  plist  <- admixr2:::.admMakeParamsList(500L, e$pinfo, 1L)
  h <- 1e-3
  g_an <- admixr2:::.admGrad(e$p0, e$pinfo, e$studies, z_list, e$rxMod, "cp",
                             plist, 1L, h, e$sens, use_central = TRUE)
  g_fd <- vapply(seq_along(e$p0), function(k) {
    ph <- e$p0; ph[k] <- ph[k] + h
    pl <- e$p0; pl[k] <- pl[k] - h
    (admixr2:::.admNLL(ph, e$pinfo, e$studies, z_list, e$rxMod, "cp", plist, 1L) -
     admixr2:::.admNLL(pl, e$pinfo, e$studies, z_list, e$rxMod, "cp", plist, 1L)) / (2 * h)
  }, numeric(1))
  ok <- abs(g_fd) > 1e-6
  expect_true(all(abs(g_an[ok] / g_fd[ok] - 1) < 0.05))
})

test_that(".adghGrad matches FD for a fixed-theta model (was silently dropping the study)", {
  e    <- .ft_setup()
  grid <- admixr2:::.adghNodeGrid(5L, e$pinfo$n_eta)
  g_an <- admixr2:::.adghGrad(e$p0, e$pinfo, e$studies, e$sens, e$rxMod, "cp",
                              grid, 1L, 1e-3)
  h <- 1e-4
  g_fd <- vapply(seq_along(e$p0), function(k) {
    ph <- e$p0; ph[k] <- ph[k] + h
    pl <- e$p0; pl[k] <- pl[k] - h
    (admixr2:::.adghNLL(ph, e$pinfo, e$studies, e$rxMod, "cp", grid, 1L) -
     admixr2:::.adghNLL(pl, e$pinfo, e$studies, e$rxMod, "cp", grid, 1L)) / (2 * h)
  }, numeric(1))
  ok <- abs(g_fd) > 1e-6
  expect_true(all(abs(g_an[ok] / g_fd[ok] - 1) < 0.05),
    info = sprintf("max ratio deviation %.4f",
                   max(abs(g_an[ok] / g_fd[ok] - 1))))
  expect_false(any(g_an == 0))     # a dropped study would leave zeros behind
})

test_that(".adghGrad degrades to FD (not a dropped study) when the sens solve fails", {
  # Force every sens solve to fail by handing .adghGrad a sens model whose columns
  # do not exist. The gradient must still be the (FD) gradient of the objective --
  # NOT a gradient with the study silently missing.
  e    <- .ft_setup()
  grid <- admixr2:::.adghNodeGrid(5L, e$pinfo$n_eta)
  broken <- e$sens
  broken$sens_cols <- paste0("no_such_column_", seq_along(broken$sens_cols))

  g_broken <- admixr2:::.adghGrad(e$p0, e$pinfo, e$studies, broken, e$rxMod, "cp",
                                  grid, 1L, 1e-3)
  g_fd <- admixr2:::.adghFDGrad(e$p0, e$pinfo, e$studies, e$rxMod, "cp", grid, 1L, 1e-3)

  expect_true(all(is.finite(g_broken)))
  expect_false(all(g_broken == 0))
  expect_equal(unname(g_broken), unname(g_fd), tolerance = 1e-8)
})
