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

test_that("the sens model re-derives its fields on a disk-cache hit (worker path)", {
  # The disk-cache branch (.admLoadSensModel, file.exists) runs on the SECOND load:
  # the first load writes the qs2 file, the reload reads it and OVERWRITES
  # rename_map / fixed_theta from the parent's fresh derivation rather than trusting
  # the file (sens_cols is keyed by `unpaired` in the cache path, so a hit is
  # guaranteed the same direction set and it is trusted). These are the fields a
  # parallel WORKER inherits -- it reads the same file and cannot re-derive -- so a
  # stale position-indexed map here would silently diverge the parallel fit from the
  # sequential one.
  ui    <- suppressMessages(rxode2::rxode2(fixed_theta_fn))
  first <- suppressMessages(admixr2:::.admLoadSensModel(ui))
  if (is.null(first)) skip("sens model unavailable")
  expect_true(file.exists(first$cache_file))

  # seed the file with a STALE (position-indexed) map + missing fixed value, then
  # force a disk read: the reload must return the corrected fields, not the stale
  # ones, proving the overwrite (not a fill-if-missing) actually runs.
  stale <- first
  stale$rename_map[["tka"]] <- "THETA[2]"    # the pre-fix bug
  stale$fixed_theta <- numeric(0)
  qs2::qs_save(stale, first$cache_file)

  back <- suppressMessages(admixr2:::.admLoadSensModel(ui))
  expect_equal(unname(back$rename_map[["tka"]]), "THETA[3]")     # corrected, not "THETA[2]"
  expect_equal(back$fixed_theta, c(`THETA[2]` = log(20)))        # re-derived, not empty
  # sens_cols / the direction set are keyed by `unpaired` in the cache path, so a
  # hit is guaranteed the same direction set -- trusted from the cache (a worker
  # also gets sens_cols passed to it directly), so it round-trips unchanged.
  expect_equal(back$sens_cols, first$sens_cols)

  # the disk-hit branch corrects the RETURNED result but does not rewrite the file,
  # so the stale seed we wrote is still on disk -- wipe it so it cannot pollute a
  # later test that reads this cache path.
  suppressWarnings(file.remove(first$cache_file))
})

test_that("adfo Jacobian (.admSimulateSensRows) fills the fixed theta", {
  # .adfoGetMuJ -> .adfoGetMuJBatch -> .admSimulateSensRows: the batched sens solve
  # the FO estimator uses. It has its own fixed-theta fill; without it the FO
  # linearisation of a fixed-theta model would solve a mis-parameterised model.
  e   <- .ft_setup()
  ov  <- admixr2:::.admOutputVar(e$ui)
  pm  <- admixr2:::.admMakeParamsList(1L, e$pinfo, 1L)[[1L]]
  mj  <- admixr2:::.adfoGetMuJ(e$pars, e$pinfo, e$study, e$sens, e$rxMod, ov, pm, 1L)
  expect_false(is.null(mj))
  expect_length(mj$mu, length(e$study$times))
  expect_true(all(is.finite(mj$mu)))
  expect_equal(dim(mj$J), c(length(e$study$times), e$pinfo$n_eta))
  # the FO mean must match the simulation model at eta = 0
  eta0   <- matrix(0, 1L, e$pinfo$n_eta, dimnames = list(NULL, e$pinfo$eta_col_names))
  cp_sim <- admixr2:::.admSimulate(e$rxMod, e$pars$struct, e$pinfo$sigma_names, eta0,
                                   e$study, ov, pm, 1L)
  expect_equal(mj$mu, as.numeric(cp_sim), tolerance = 1e-4)
})

test_that(".admGradBatch fills the fixed theta (grad-FD Hessian path)", {
  # .admGradBatch is the batched gradient behind .admCalcCov. Its sens block has
  # its own inlined fixed-theta fill; exercise it directly on a fixed-theta model.
  e      <- .ft_setup(n_sim = 200L)
  z_list <- admixr2:::.admMakeZ(200L, e$pinfo, 1L, "sobol")
  plist  <- admixr2:::.admMakeParamsList(200L, e$pinfo, 1L)
  p1 <- e$p0
  p2 <- e$p0; p2[["tcl"]] <- p2[["tcl"]] + 1e-3
  G  <- admixr2:::.admGradBatch(list(p1, p2), e$pinfo, e$studies, z_list, e$rxMod,
                                "cp", plist, 1L, 1e-3, e$sens)
  expect_equal(dim(G), c(2L, length(e$p0)))
  expect_true(all(is.finite(G)))
  expect_false(any(apply(G, 1, function(row) all(row == 0))))   # not a dropped/empty row
})

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

test_that("the sens cache key carries a schema tag (pre-fix caches miss)", {
  # A cache written by an older admixr2 carries the position-indexed rename_map and
  # no fixed_theta. The key now includes a schema tag, so such a file is not a hit;
  # a worker cannot then inherit a stale map. (What happens ON a hit -- the fields
  # are re-derived rather than trusted -- is the "disk-cache hit" test above; note
  # the disk file itself may legitimately hold stale data, which is exactly why the
  # reload overwrites rather than trusts it, so this test does not read its
  # contents.)
  e <- .ft_setup()

  no_tag <- file.path(rxode2::rxTempDir(),
                      paste0("adm-sens-", digest::digest(e$ui$lstExpr), ".qs2"))
  expect_false(identical(normalizePath(e$sens$cache_file, mustWork = FALSE),
                         normalizePath(no_tag, mustWork = FALSE)))
  expect_match(basename(e$sens$cache_file), "^adm-sens-.*\\.qs2$")
})

test_that("a FIXED omega is rejected with an actionable error (not a bare subscript error)", {
  # Dropping a fixed omega from eta_rows leaves omega_init too small, so the parse
  # died with "subscript out of bounds"; and had it not, the eta's variance would
  # have been silently excluded from the model instead of held fixed.
  fixed_omega_fn <- function() {
    ini({
      tcl <- log(5); tv <- log(20); add.err <- 0.1
      eta.cl ~ fix(0.09)
      eta.v  ~ 0.04
    })
    model({
      cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      d/dt(central) <- -(cl / v) * central
      cp <- central / v
      cp ~ add(add.err)
    })
  }
  ui <- suppressMessages(rxode2::rxode2(fixed_omega_fn))
  expect_error(admixr2:::.admParseIniDf(ui$iniDf, ui),
               regexp = "does not support FIXED omega")
})

test_that(".adghGrad degrades to FD when the predicted V is singular (was a dropped study)", {
  # A singular predicted V makes `G` NULL. That used to `next`, so the study was
  # dropped from the gradient -- and with a single study the caller got a FINITE,
  # ALL-ZERO gradient, which nloptr reads as a stationary point. It now degrades to
  # FD, which reports NaN here (the objective itself is not finite at this point) --
  # a failed evaluation the optimizer can see, rather than a fake optimum.
  #
  # Forced by driving omega and sigma to exactly 0 (exp(-1500)), so every quadrature
  # node predicts the same value at every time and V collapses to the zero matrix.
  e    <- .ft_setup()
  grid <- admixr2:::.adghNodeGrid(5L, e$pinfo$n_eta)

  # `G` only exists on the COV branch -- a diagonal V auto-detects as method="var",
  # which inverts the diagonal directly and never builds G. Use a full V.
  times  <- c(0.5, 1, 2, 4, 8)
  E_true <- c(1.8, 2.7, 3.1, 2.3, 0.9)
  A      <- diag((0.2 * E_true)^2)
  A[1, 2] <- A[2, 1] <- 0.01                      # off-diagonal -> method = "cov"
  study  <- admixr2:::.admNormaliseStudy(
    list(E = E_true, V = A, n = 100L, times = times, ev = rxode2::et(amt = 100)), "s")
  study$ev_full <- study$ev |> rxode2::et(study$times)
  studies <- list(s = study)
  expect_equal(study$method, "cov")

  p_sing <- e$p0
  p_sing[grep("^logchol", names(p_sing))] <- -1500
  p_sing[e$pinfo$sigma_names]             <- -1500
  expect_true(all(diag(admixr2:::.admUnpack(p_sing, e$pinfo)$omega) == 0))  # V -> 0

  g <- admixr2:::.adghGrad(p_sing, e$pinfo, studies, e$sens, e$rxMod, "cp",
                           grid, 1L, 1e-3)

  # the old (dropped-study) behaviour: a FINITE gradient of exactly zeros, which
  # nloptr would read as a stationary point
  expect_false(isTRUE(all(g == 0)))

  # what it must be instead: the FD gradient of the WHOLE objective
  g_fd <- admixr2:::.adghFDGrad(p_sing, e$pinfo, studies, e$rxMod, "cp",
                                grid, 1L, 1e-3)
  expect_equal(unname(g), unname(g_fd))
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
