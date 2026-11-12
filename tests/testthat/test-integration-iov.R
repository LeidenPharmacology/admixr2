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
    occasions = c(rep(1, 5), rep(2, 5)), occasion_starts = 96))
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
  # post-fit aggregate predictions survive the occasion rewrite
  expect_gt(nrow(admixr2::admMoments(f)), 0L)
  expect_gt(length(suppressMessages(suppressWarnings(plot(f, which = "mean")))), 0L)
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

test_that("a published IOV model as a source generates what its explicit form does", {
  pub <- function() {
    ini({ tcl <- log(2.7); tv <- log(30); eta.cl ~ 0.09; iov.cl ~ 0.05 | occ; add.sd <- 0.3 })
    model({ cl <- exp(tcl + eta.cl + iov.cl); v <- exp(tv)
      d/dt(centr) = -cl / v * centr; cp = centr / v; cp ~ add(add.sd) })
  }
  tt <- c(2, 6, 12, 26, 50, 98, 102, 110, 122, 146)
  ev <- rxode2::et(amt = 100, ii = 24, addl = 6)
  est <- c(tcl = log(2.7), tv = log(30))
  a <- suppressWarnings(suppressMessages(admixr2::datagen(admixr2::admStudies(
    s = admixr2::admStudy(model = .iov_truth, est = est, n = 2000, times = tt,
                          dose = ev, simulated = TRUE)))))[[1]]
  sp <- suppressWarnings(suppressMessages(admixr2::admStudy(
    model = pub, est = est, n = 2000, times = tt, dose = ev, simulated = TRUE,
    occasions = c(rep(1, 5), rep(2, 5)), occasion_starts = 96)))
  expect_equal(sp$at$iov_s2, 96)                    # the dose that opens occasion 2
  b <- suppressWarnings(suppressMessages(admixr2::datagen(admixr2::admStudies(s = sp))))[[1]]
  # the labels survive datagen(), so fitting what it returns places the occasions the same way
  expect_equal(b$occasions, c(rep(1, 5), rep(2, 5)))
  expect_equal(b$cov$iov_t2, 96)
  rp <- admixr2:::.admIovPad(list(s = b), list(var = "occ", levels = c(1, 2)))
  expect_equal(rp$s$cov$iov_t2, 96)
  expect_equal(b$E, a$E, tolerance = 1e-6)
  expect_equal(b$V, a$V, tolerance = 1e-6)
  # an IOV model without occasions, and occasions without an IOV model, are refused
  expect_error(suppressWarnings(suppressMessages(admixr2::admStudy(
    model = pub, est = est, n = 10, times = tt, dose = ev))), "needs `occasions`")
})

test_that("a model without an IOV term can still be given occasions (it asserts none)", {
  plain <- function() {
    ini({ tcl <- log(2.7); tv <- log(30); eta.cl ~ 0.14; add.sd <- 0.3 })
    model({ cl <- exp(tcl + eta.cl); v <- exp(tv)
      d/dt(centr) = -cl / v * centr; cp = centr / v; cp ~ add(add.sd) })
  }
  tt <- c(2, 6, 12, 26, 50, 98, 102, 110, 122, 146)
  ev <- rxode2::et(amt = 100, ii = 24, addl = 6)
  est <- c(tcl = log(2.7), tv = log(30))
  lab <- c(rep(1, 5), rep(2, 5))
  mk <- function(...) suppressWarnings(suppressMessages(admixr2::admStudy(
    model = plain, est = est, n = 500, times = tt, dose = ev, simulated = TRUE, ...)))
  with <- mk(occasions = lab, occasion_starts = 96); without <- mk()
  expect_equal(with$occasions, lab)
  g <- function(s) suppressWarnings(suppressMessages(
    admixr2::datagen(admixr2::admStudies(s = s))))[[1]]
  expect_equal(g(with)$V, g(without)$V)             # generated as written
  # in the fitted model the labels split its times into occasions
  padded <- admixr2:::.admIovPad(list(s = with), list(var = "occ", levels = c(1, 2)))
  expect_equal(padded$s$at$iov_t2, 96)
})

test_that(".admIovCov places occasion starts from the concrete event table", {
  ev <- rxode2::et(amt = 100, ii = 24, addl = 6)
  tt <- c(2, 26, 98, 102)
  cv <- admixr2:::.admIovCov(list(WT = 70), tt, c(1, 1, 2, 2), 96, c(1, 2), ev)
  expect_equal(cv$iov_t2, 96)                      # the dose, not the first sample at 98
  expect_equal(cv$WT, 70)
})

test_that("materialising does not leak one study's occasions onto another", {
  ev <- rxode2::et(amt = 100, ii = 24, addl = 6)
  a <- admixr2::admStudy(E = c(1, 2, 3, 4), sd = rep(1, 4), n = 10, times = c(2, 26, 98, 102),
                         dose = ev, occasions = c(1, 1, 2, 2), occasion_starts = 96)
  raw <- list(E = c(1, 2), V = c(1, 1), n = 10, times = c(2, 26), ev = ev)
  for (ord in list(list("a", "b"), list("b", "a"))) {
    st <- list(a = a, b = raw)[unlist(ord)]
    m <- admixr2:::.admMaterialise(st)
    expect_null(m$b$occasions)
    expect_null(m$b$occasion_starts)
    expect_equal(m$a$occasions, c(1, 1, 2, 2))
  }
})

test_that("a published IOV source keeps its own boundaries when the analysis has more levels", {
  pub <- function() {
    ini({ tcl <- log(2.7); tv <- log(30); eta.cl ~ 0.09; iov.cl ~ 0.05 | occ; add.sd <- 0.3 })
    model({ cl <- exp(tcl + eta.cl + iov.cl); v <- exp(tv)
      d/dt(centr) = -cl / v * centr; cp = centr / v; cp ~ add(add.sd) })
  }
  mk <- function() suppressWarnings(suppressMessages(admixr2::admStudy(
    model = pub, est = c(tcl = log(2.7), tv = log(30)), n = 500,
    times = c(102, 110, 202, 206), occasions = c(2, 2, 3, 3), occasion_starts = 192,
    dose = rxode2::et(amt = 100, ii = 24, addl = 9), simulated = TRUE)))
  g <- function(s) suppressWarnings(suppressMessages(
    admixr2::datagen(admixr2::admStudies(s = s))))[[1]]
  alone <- g(mk())
  s3 <- mk(); s3[[".iov_levels"]] <- c(1, 2, 3)     # what an analysis model with 3 levels does
  glob <- g(s3)
  expect_equal(glob$E, alone$E, tolerance = 1e-6)
  expect_equal(glob$V, alone$V, tolerance = 1e-6)
  expect_equal(glob$cov$iov_t3, 192)                # global indexing for the analysis model
})
