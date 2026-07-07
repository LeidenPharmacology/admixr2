# Tier 2 integration: adgh gradient for an eta-less (unpaired) structural theta.
#
# Regression for the bug where .adghGrad detected "unpaired" struct thetas from
# the eta-indexed struct_eta_idx (which is empty whenever every eta is paired),
# so a structural parameter WITHOUT a random effect received a hard zero gradient
# and the optimizer never moved it off its initial value. The fix derives the
# unpaired set from the struct-indexed struct_has_eta.
#
# Here ka has no eta.ka, so d(NLL)/d(tka) must be finite, non-zero at an
# off-optimum point, and equal to an independent finite difference of the adgh
# NLL. Pre-fix it was exactly 0.

.adgh_unpaired_setup <- function() {
  skip_on_cran()
  skip_if_not_installed("rxode2")

  # 1-cmt oral; tka is unpaired (no eta.ka), tcl/tv are mu-referenced.
  m <- function() {
    ini({
      tka <- log(4)   # deliberately off the data-generating ka = 1.5
      tcl <- log(5)
      tv  <- log(20)
      add.sd <- c(0, 0.3)
      eta.cl ~ 0.09
      eta.v  ~ 0.04
    })
    model({
      ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      d/dt(depot)  <- -ka * depot
      d/dt(center) <-  ka * depot - (cl / v) * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }
  ui <- suppressMessages(tryCatch(rxode2::rxode2(m), error = function(e) NULL))
  if (is.null(ui)) skip("rxode2 model parse failed")

  pinfo <- admixr2:::.admParseIniDf(ui$iniDf, ui)
  # sanity: parser must see tka as unpaired and tcl/tv as paired
  if (is.null(pinfo$struct_has_eta) ||
      !identical(unname(pinfo$struct_has_eta), c(FALSE, TRUE, TRUE)))
    skip("unexpected struct_has_eta layout")
  out_var <- "cp"

  # data-generating population mean (analytic 1-cmt oral) at ka = 1.5
  oral_mean <- function(ka, cl, v, D, t) {
    ke <- cl / v; (D / v) * (ka / (ka - ke)) * (exp(-ke * t) - exp(-ka * t)) }
  times  <- c(0.5, 1, 2, 4, 8)
  E_true <- oral_mean(1.5, 5, 20, 100, times)
  V_true <- diag((0.25 * E_true)^2)

  study <- list(E = E_true, V = V_true, n = 200L, times = times,
                ev = rxode2::et(amt = 100))
  study <- admixr2:::.admNormaliseStudy(study, "s")
  study$ev_full <- study$ev |> rxode2::et(study$times)
  studies <- list(s = study)

  # ORDERING INVARIANT: sens model before sim model
  sensModel <- tryCatch(admixr2:::.admLoadSensModel(ui), error = function(e) NULL)
  if (is.null(sensModel)) skip("sensitivity model unavailable")
  rxMod <- tryCatch(admixr2:::.admLoadModel(ui), error = function(e) NULL)
  if (is.null(rxMod)) skip("model compilation failed")

  grid <- admixr2:::.adghNodeGrid(5L, pinfo$n_eta)
  p0   <- admixr2:::.admBuildOptVec(pinfo)$p0

  g <- admixr2:::.adghGrad(p0, pinfo, studies, sensModel, rxMod, out_var, grid, cores = 1L)

  # independent central FD of the adgh NLL w.r.t. tka
  h  <- max(abs(p0[["tka"]]), 0.1) * 1e-4
  pp <- p0; pp[["tka"]] <- pp[["tka"]] + h
  pm <- p0; pm[["tka"]] <- pm[["tka"]] - h
  fd_tka <- (admixr2:::.adghNLL(pp, pinfo, studies, rxMod, out_var, grid, 1L) -
             admixr2:::.adghNLL(pm, pinfo, studies, rxMod, out_var, grid, 1L)) / (2 * h)

  list(grad = g, fd_tka = fd_tka)
}

test_that("adgh gives an eta-less structural parameter a finite, non-zero gradient", {
  s <- .adgh_unpaired_setup()
  g <- s$grad
  expect_true(all(is.finite(g)))
  # the regression: tka (no eta.ka) must NOT be frozen at zero
  expect_gt(abs(unname(g[["tka"]])), 1e-2)
  expect_gt(abs(unname(g[["tcl"]])), 1e-2)   # paired thetas obviously non-zero
  expect_gt(abs(unname(g[["tv"]])),  1e-2)
})

test_that("adgh eta-less gradient matches a finite difference of the NLL", {
  s <- .adgh_unpaired_setup()
  expect_equal(unname(s$grad[["tka"]]), s$fd_tka, tolerance = 0.1)
})
