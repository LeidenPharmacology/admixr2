# .admGradBatch drives the gradient-FD Hessian (covMethod = "r" with a gradient).
# It must agree with .admGrad, which is the reference gradient. On multi-output
# fits it did not: it read the model-level `output_var` instead of each unit's
# own output, looped every sigma rather than the unit's own, and accumulated the
# sigma gradient at a RUNNING slot (n_s + 1, n_s + 2, ...) instead of the global
# sigma index -- so a second endpoint's sigma gradient landed in the first
# endpoint's slot and its own stayed zero. Observed before the fix on this model:
# the tq gradient came out as -1.08e6 against a true +523.7, and add.csf as 0.

test_that(".admGradBatch matches .admGrad on a multi-output model", {
  skip_on_cran()
  skip_if_not_installed("rxode2")
  skip_if_not_installed("nlmixr2")

  two_out_fn <- function() {
    ini({
      tcl     <- log(5)
      tv      <- log(20)
      tq      <- log(2)
      add.err <- 0.1
      add.csf <- 0.05
      eta.cl  ~ 0.09
      eta.v   ~ 0.04
    })
    model({
      cl <- exp(tcl + eta.cl)
      v  <- exp(tv  + eta.v)
      q  <- exp(tq)
      d/dt(central) <- -(cl / v) * central - q * central + q * csfc
      d/dt(csfc)    <-  q * central - q * csfc
      cp   <- central / v
      cCSF <- csfc / v
      cp   ~ add(add.err)
      cCSF ~ add(add.csf)
    })
  }

  ui        <- rxode2::rxode2(two_out_fn)
  sensModel <- admixr2:::.admLoadSensModel(ui)
  rxMod     <- admixr2:::.admLoadModel(ui)
  pinfo     <- admixr2:::.admParseIniDf(ui$iniDf, ui)
  ov        <- admixr2:::.admBuildOptVec(pinfo)

  tp <- c(0.5, 1, 2, 4)
  tc <- c(1, 4, 8)
  Ep <- (100 / 20) * exp(-(5 / 20) * tp)
  Ec <- 0.3 * (100 / 20) * exp(-(5 / 20) * tc)

  raw <- list(observations = list(
    plasma = list(output = "cp",   ev = rxode2::et(amt = 100), times = tp,
                  E = Ep, V = diag((0.3 * Ep)^2), n = 40L),
    csf    = list(output = "cCSF", ev = rxode2::et(amt = 100), times = tc,
                  E = Ec, V = diag((0.3 * Ec)^2), n = 40L)))

  st      <- list(m1 = admixr2:::.admNormaliseStudy(raw, "m1", "cp"))
  studies <- admixr2:::.admBuildEvFull(admixr2:::.admFlattenStudies(st), tag_cmt = TRUE)

  n_sim <- 400L
  z     <- admixr2:::.admMakeZ(n_sim, pinfo, length(studies), "sobol")
  pl    <- admixr2:::.admMakeParamsList(n_sim, pinfo, length(studies))

  g_ref <- admixr2:::.admGrad(ov$p0, pinfo, studies, z, rxMod, "cp", pl, 1L,
                              1e-4, sensModel)
  g_bat <- as.numeric(
    admixr2:::.admGradBatch(list(ov$p0), pinfo, studies, z, rxMod, "cp", pl, 1L,
                            1e-4, sensModel)[1L, ])

  expect_equal(g_bat, unname(g_ref), tolerance = 1e-10)

  # every sigma must receive a non-zero gradient: the second endpoint's used to be
  # identically zero because its contribution was written into the first's slot
  n_s   <- length(pinfo$struct_names)
  n_sig <- length(pinfo$sigma_names)
  expect_true(all(abs(g_bat[n_s + seq_len(n_sig)]) > 0))
})
