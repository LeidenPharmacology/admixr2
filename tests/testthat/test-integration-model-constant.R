# Regression: hard-coded numeric constants in a model's `model({})` block (common
# in PBPK/CNS models -- a fixed brain volume, flow, plasma fraction, ...) must keep
# their value. admixr2 used to hand-fill every model parameter it didn't set with
# 0, zeroing such a constant (e.g. `qout / vb` -> divide-by-zero -> NA objective).
# It now supplies only the parameters it varies and lets rxSolve fill the rest
# from the model's own defaults.

skip_on_cran()
skip_if_not_installed("rxode2")

.const_model <- function() {
  ini({
    tcl <- log(1); tv1 <- log(10); tqin <- log(3); tqout <- log(6)
    prop.cp <- 0.05
    eta.cl ~ 0.09
    eta.v1 ~ 0.04
  })
  model({
    cl <- exp(tcl + eta.cl); v1 <- exp(tv1 + eta.v1)
    qin <- exp(tqin); qout <- exp(tqout)
    vb <- 5                                   # hard-coded constant (was zeroed)
    d/dt(central) <- -(cl/v1)*central - (qin/v1)*central + (qout/vb)*brain
    d/dt(brain)   <-  (qin/v1)*central - (qout/vb)*brain
    cp <- central / v1
    cp ~ prop(prop.cp)
  })
}

test_that("a hard-coded model constant keeps its value (finite objective)", {
  study <- list(n = 40L, ev = rxode2::et(amt = 100, cmt = "central"),
                times = c(0.5, 1, 2, 4, 8),
                E = c(8.6, 7.5, 6.1, 4.9, 3.7),
                V = c(1.6, 1.3, 0.9, 0.7, 0.6)^2)
  f <- suppressMessages(nlmixr2est::nlmixr2(.const_model, admData(), est = "admc",
    control = admControl(studies = list(rat = study), n_sim = 200L, grad = "sens",
                         seed = 1L, covMethod = "none", maxeval = 20L)))
  expect_true(is.finite(f$objective))
})

test_that("the model constant is used, not zeroed (Kp = Qin/Qout recovered)", {
  # cp/cb ratio identifies Qin/Qout; with vb zeroed the solve was NA. Fit both
  # outputs and check the brain:plasma partition is recovered near truth (0.5).
  cns <- function() {
    ini({ tcl <- log(1); tv1 <- log(10); tqin <- log(3); tqout <- log(6)
          prop.cp <- 0.05; add.cb <- 0.02; eta.cl ~ 0.09; eta.v1 ~ 0.04 })
    model({ cl <- exp(tcl+eta.cl); v1 <- exp(tv1+eta.v1)
      qin <- exp(tqin); qout <- exp(tqout); vb <- 5
      d/dt(central) <- -(cl/v1)*central - (qin/v1)*central + (qout/vb)*brain
      d/dt(brain)   <-  (qin/v1)*central - (qout/vb)*brain
      cp <- central/v1; cb <- brain/vb; cp ~ prop(prop.cp); cb ~ add(add.cb) }) }
  ev <- rxode2::et(amt = 100, cmt = "central")
  study <- list(n = 60L, ev = ev, observations = list(
    plasma = list(output = "cp", times = c(0.5,1,2,4,8,12),
                  E = c(8.712,7.89,6.858,5.752,4.201,3.053),
                  V = c(1.46,1.126,0.857,0.68,0.754,0.745)^2),
    brain  = list(output = "cb", times = c(1,2,4,8,12),
                  E = c(3.016,3.42,3.069,2.238,1.644),
                  V = c(0.461,0.451,0.361,0.34,0.351)^2)))
  f <- suppressMessages(nlmixr2est::nlmixr2(cns, admData(c("cp","cb")), est = "admc",
    control = admControl(studies = list(rat = study), n_sim = 300L, grad = "sens",
                         seed = 1L, covMethod = "none", maxeval = 120L)))
  est <- f$env$admExtra$struct
  expect_equal(unname(exp(est[["tqin"]] - est[["tqout"]])), 0.5, tolerance = 0.1)
})
