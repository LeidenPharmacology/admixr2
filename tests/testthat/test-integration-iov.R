# IOV end to end: simulate aggregate data from a model with a TRUE per-occasion
# random effect, then recover its variance with `iov.cl ~ v | occ` (nlmixr2est's
# own rewrite) under adgh and admc.

skip_on_cran()
skip_if_not_installed("rxode2")

.iov_truth <- function() {
  ini({ tcl <- log(2.7); tv <- log(30); eta.cl ~ 0.09; eta.k1 ~ 0.05; eta.k2 ~ 0.05
        add.sd <- 0.3 })
  model({ iov <- eta.k1 * (time < 96) + eta.k2 * (time >= 96)
    cl <- exp(tcl + eta.cl + iov); v <- exp(tv)
    d/dt(centr) = -cl / v * centr; cp = centr / v; cp ~ add(add.sd) })
}
.iov_fit_model <- function() {
  ini({ tcl <- log(2.5); tv <- log(28); eta.cl ~ 0.06; iov.cl ~ 0.1 | occ
        add.sd <- 0.25 })
  model({ cl <- exp(tcl + eta.cl + iov.cl); v <- exp(tv)
    d/dt(centr) = -cl / v * centr; cp = centr / v; cp ~ add(add.sd) })
}
.iov_studies <- function(N) suppressWarnings(suppressMessages({
  tt <- c(2, 6, 12, 26, 50, 98, 102, 110, 122, 146)
  ev <- rxode2::et(amt = 100, ii = 24, addl = 6)
  sim <- admixr2::datagen(admixr2::admStudies(s = admixr2::admStudy(
    model = .iov_truth, est = c(tcl = log(2.7), tv = log(30)), n = N, times = tt,
    dose = ev, simulated = TRUE)))[[1]]
  admixr2::admStudies(s = admixr2::admStudy(
    E = sim$E, V = sim$V, n = N, times = tt, dose = ev,
    occasions = c(rep(1, 5), rep(2, 5))))
}))

test_that("adgh recovers a per-occasion variance and keeps it apart from BSV", {
  st <- .iov_studies(1000)
  f <- suppressMessages(suppressWarnings(nlmixr2est::nlmixr2(
    .iov_fit_model, admixr2::admData(occasions = c(occ = 2)), est = "adgh",
    control = admixr2::adghControl(studies = st, cores = 2L, covMethod = "none"))))
  expect_equal(f$omega$occ[["iov.cl", "iov.cl"]], 0.05, tolerance = 0.05)
  expect_equal(f$omega$id[["eta.cl", "eta.cl"]], 0.09, tolerance = 0.05)
  expect_equal(unname(f$fixef["tcl"]), log(2.7), tolerance = 0.01)
  expect_true(all(c("ID", "occ", "iov.cl") %in% names(f$iov$occ)))
})

test_that("admc recovers a per-occasion variance", {
  st <- .iov_studies(1000)
  f <- suppressMessages(suppressWarnings(nlmixr2est::nlmixr2(
    .iov_fit_model, admixr2::admData(occasions = c(occ = 2)), est = "admc",
    control = admixr2::admControl(studies = st, cores = 2L, covMethod = "none",
                                   n_sim = 3000L))))
  expect_equal(f$omega$occ[["iov.cl", "iov.cl"]], 0.05, tolerance = 0.25)
  expect_equal(f$omega$id[["eta.cl", "eta.cl"]], 0.09, tolerance = 0.25)
})

test_that("adfo refuses an IOV model instead of fitting it as a plain eta", {
  st <- .iov_studies(200)
  expect_error(suppressMessages(suppressWarnings(nlmixr2est::nlmixr2(
    .iov_fit_model, admixr2::admData(occasions = c(occ = 2)), est = "adfo",
    control = admixr2::adfoControl(studies = st, cores = 2L)))),
    "adgh and admc only")
})
