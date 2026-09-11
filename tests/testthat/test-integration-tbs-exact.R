skip_if_not_installed("rxode2")
skip_on_cran()

# The exact node-wise composition for transform-both-sides endpoints.
#
# .admTBSNodeParts()/.admTBSAggregate() replace a second-order delta expansion
# about (mu_struct, var_f) that did NOT converge: measured against a 201-node
# evaluation of the defining integral, its relative error in V was a floor no
# node count removed -- 3.4e-03 (boxCox) to 3.1e-02 (probitNorm at omega 0.49),
# flat from 7 nodes to 25. Which also meant `n_nodes` bought a TBS fit nothing.
#
# The unit-level algebra is pinned in test-adfweight.R (moments vs a simulation
# of the conditional law, the weight vs simulated studies, S == V_pred, and the
# information equality). This file pins the paths that go THROUGH the package --
# every one of which was verified by hand while developing and is here so it
# stays verified.

.tbs_mk <- function(err, extra, struct = "cl <- exp(tcl + eta.cl); v <- exp(tv)")
  eval(parse(text = sprintf('function() {
    ini({ tcl <- log(5); tv <- log(20); %s })
    model({ %s; cp <- linCmt(); %s }) }', extra, struct, err)))

.tbs_times <- c(1, 2, 4, 8, 12)
.tbs_study <- function() {
  E0 <- 100 / 20 * exp(-0.25 * .tbs_times)
  list(s = list(E = E0, V = diag((0.3 * E0)^2), n = 200L, times = .tbs_times,
                ev = rxode2::et(amt = 100)))
}

test_that("only TBS takes the exact path; every other family is untouched", {
  # The guarantee that this change is inert outside TBS. .admTBSNodeParts()
  # returns NULL for any non-TBS form, so those units run byte-identical code --
  # by construction rather than by hope, which is what this checks.
  st <- .tbs_study()
  cases <- list(
    add    = c("cp ~ add(aa)",              "aa <- 0.5; eta.cl ~ 0.16"),
    prop   = c("cp ~ prop(bb)",             "bb <- 0.15; eta.cl ~ 0.16"),
    pow    = c("cp ~ pow(bb, pw)",          "bb <- 0.15; pw <- fix(0.8); eta.cl ~ 0.16"),
    lnorm  = c("cp ~ lnorm(aa)",            "aa <- 0.15; eta.cl ~ 0.16"),
    comb2  = c("cp ~ add(aa) + prop(bb)",   "aa <- 0.3; bb <- 0.1; eta.cl ~ 0.16"),
    boxCox = c("cp ~ add(aa) + boxCox(lam)","aa <- 0.3; lam <- fix(0.5); eta.cl ~ 0.16"))
  for (nm in names(cases)) {
    ui  <- suppressMessages(rxode2::rxode2(.tbs_mk(cases[[nm]][1], cases[[nm]][2])))
    ov  <- admixr2:::.admOutputVar(ui); rx <- admixr2:::.admLoadModel(ui)
    ctl <- adghControl(studies = st, n_nodes = 5L, print = 0L, covMethod = "none",
                       grad = "none")
    pin <- admixr2:::.admDriverPinfo(ui, ctl)
    u   <- admixr2:::.admDriverUnits(st, ui, ov)
    g   <- admixr2:::.adghNodeGrid(5L, pin$n_eta)
    pars <- admixr2:::.admUnpack(admixr2:::.admBuildOptVec(pin)$p0, pin)
    gg  <- admixr2:::.adghGrid(pars, pin, g)
    pm  <- admixr2:::.admMakeParamsList(nrow(gg$eta), pin, 1L)[[1L]]
    cp  <- admixr2:::.admSimulate(rx, pars$struct, pin$sigma_names, gg$eta,
                                  u$studies[[1L]], ov, pm, 1L, 1e6, pin$sigdig)
    arr <- admixr2:::.admUnitResidRows(pin, ov, pars$sigma_var,
                                       length(.tbs_times), phi = attr(cp, "phi"))
    np  <- admixr2:::.admTBSNodeParts(cp, arr)
    if (nm == "boxCox") expect_false(is.null(np), info = nm)
    else                expect_null(np, info = nm)
  }
})

test_that("n_nodes and resid_nodes both reach a TBS objective", {
  # Before this change neither did: the expansion depended on (mu_struct, var_f)
  # and on a fixed residual quadrature, so raising either returned the same
  # number. Both must now move the objective and settle.
  st <- .tbs_study()
  fn <- .tbs_mk("cp ~ add(aa) + boxCox(lam)", "aa <- 0.4; lam <- fix(0.5); eta.cl ~ 0.16")
  ui <- suppressMessages(rxode2::rxode2(fn))
  ov <- admixr2:::.admOutputVar(ui); rx <- admixr2:::.admLoadModel(ui)
  nll_at <- function(nq, rn) {
    ctl <- adghControl(studies = st, n_nodes = nq, print = 0L, covMethod = "none",
                       grad = "none", resid_nodes = rn)
    pin <- admixr2:::.admDriverPinfo(ui, ctl)
    u   <- admixr2:::.admDriverUnits(st, ui, ov)
    g   <- admixr2:::.adghNodeGrid(nq, pin$n_eta)
    admixr2:::.adghNLL(admixr2:::.admBuildOptVec(pin)$p0, pin, u$studies, rx, ov,
                       g, 1L)
  }
  n5 <- nll_at(5L, 81L); n9 <- nll_at(9L, 81L); n15 <- nll_at(15L, 81L)
  expect_false(isTRUE(all.equal(n5, n9)))          # n_nodes MOVES it ...
  expect_equal(n9, n15, tolerance = 1e-4)          # ... and it settles
  r21 <- nll_at(9L, 21L); r201 <- nll_at(9L, 201L)
  expect_false(isTRUE(all.equal(r21, n9)))         # resid_nodes reaches it too
  expect_equal(n9, r201, tolerance = 1e-4)
})

test_that("admc scores one TBS objective, batched or not", {
  # .admNLLBatch is the evaluator behind the batched gradient. If it composed
  # differently from .admNLL the gradient would difference a function the
  # objective never evaluates -- finite, plausible, and wrong.
  st <- .tbs_study()
  fn <- .tbs_mk("cp ~ add(aa) + boxCox(lam)", "aa <- 0.4; lam <- fix(0.5); eta.cl ~ 0.16")
  ui <- suppressMessages(rxode2::rxode2(fn))
  ov <- admixr2:::.admOutputVar(ui); rx <- admixr2:::.admLoadModel(ui)
  ctl <- admControl(studies = st, n_sim = 2000L, seed = 1L, print = 0L,
                    covMethod = "none", grad = "none")
  pin <- admixr2:::.admDriverPinfo(ui, ctl)
  u   <- admixr2:::.admDriverUnits(st, ui, ov)
  pl  <- admixr2:::.admMakeParamsList(2000L, pin, 1L)
  z   <- admixr2:::.admMakeZ(2000L, pin, length(u$studies), "sobol")
  p0  <- admixr2:::.admBuildOptVec(pin)$p0 + 0.03
  one <- admixr2:::.admNLL(p0, pin, u$studies, z, rx, ov, pl, 1L)
  bat <- admixr2:::.admNLLBatch(list(p0, p0 + 0.01), pin, u$studies, z, rx, ov, pl, 1L)
  expect_equal(one, bat[[1L]])
})

test_that("datagen composes a TBS study the way the estimators score it", {
  # datagen's `mc` branch used to apply the expansion analytically while `gh`
  # went through .adghMoments, so the two disagreed by the truncation. They must
  # now agree -- and this is also why datagen cannot serve as an INDEPENDENT
  # oracle for the composition: it shares the implementation.
  fn <- .tbs_mk("cp ~ add(aa) + boxCox(lam)", "aa <- 0.4; lam <- fix(0.5); eta.cl ~ 0.16")
  spec <- list(s1 = list(n = 200L, times = .tbs_times, ev = rxode2::et(amt = 100)))
  g_mc <- datagen(spec, fn, control = datagenControl(method = "mc", n_sim = 60000L,
                                                     seed = 3L))
  g_gh <- datagen(spec, fn, control = datagenControl(method = "gh", n_nodes = 15L))
  expect_equal(g_mc$s1$E, g_gh$s1$E, tolerance = 0.01)
  expect_equal(diag(g_mc$s1$V), diag(g_gh$s1$V), tolerance = 0.02)
})

test_that("t() with a transform is refused by the weight, not guessed at", {
  # t() folds nu/(nu-2) into the variance coefficients -- exact for the combined
  # forms, where only the residual's VARIANCE enters. Not under a transform:
  # composing integrates the inverse transform over the conditional law, and a t
  # error is not an inflated-sd normal one. The third and fourth moments the
  # weight would get are a normal's, so "r,s" must decline and report "r".
  skip_if_not_installed("nlmixr2est")
  st <- .tbs_study()
  fits <- function(fn) suppressMessages(suppressWarnings(nlmixr2est::nlmixr2(
    fn, admData(), est = "adgh",
    control = adghControl(studies = st, n_nodes = 5L, maxeval = 15L, seed = 1L,
                          grad = "analytical", covMethod = "r,s", print = 0L))))
  t_tbs <- fits(.tbs_mk("cp ~ add(aa) + boxCox(lam) + t(nu)",
                        "aa <- 0.4; lam <- fix(0.5); nu <- fix(6); eta.cl ~ 0.16"))
  plain <- fits(.tbs_mk("cp ~ add(aa) + boxCox(lam)",
                        "aa <- 0.4; lam <- fix(0.5); eta.cl ~ 0.16"))
  expect_identical(t_tbs$covMethod, "r")
  expect_identical(plain$covMethod, "r,s")
})

test_that("a TBS endpoint with no random effects still composes", {
  # The single-point ensemble: with no IIV the aggregation is one node, so V is
  # the conditional variance alone. Exercises .admTBSAggregate() at Q = 1.
  skip_if_not_installed("nlmixr2est")
  st <- .tbs_study()
  fn <- .tbs_mk("cp ~ add(aa) + boxCox(lam)", "aa <- 0.4; lam <- fix(0.5)",
                struct = "cl <- exp(tcl); v <- exp(tv)")
  f <- suppressMessages(suppressWarnings(nlmixr2est::nlmixr2(
    fn, admData(), est = "adgh",
    control = adghControl(studies = st, n_nodes = 5L, maxeval = 15L, seed = 1L,
                          grad = "analytical", covMethod = "r,s", print = 0L))))
  expect_identical(f$covMethod, "r,s")
  free <- !(rownames(f$parFixedDf) %in% f$ui$iniDf$name[f$ui$iniDf$fix %in% TRUE])
  expect_true(all(is.finite(f$parFixedDf$SE[free])))
})

test_that("plot() draws the predicted moments the fit was scored against", {
  # .add_sigma() takes the simulated matrix so a TBS diagnostic composes the same
  # way the objective does; without it the panel would show the expansion's
  # moments against a fit made on the exact ones.
  skip_if_not_installed("nlmixr2est")
  skip_if_not_installed("ggplot2")
  st <- .tbs_study()
  fn <- .tbs_mk("cp ~ add(aa) + boxCox(lam)", "aa <- 0.4; lam <- fix(0.5); eta.cl ~ 0.16")
  f <- suppressMessages(suppressWarnings(nlmixr2est::nlmixr2(
    fn, admData(), est = "adgh",
    control = adghControl(studies = st, n_nodes = 5L, maxeval = 15L, seed = 1L,
                          grad = "analytical", covMethod = "r", print = 0L))))
  expect_no_error(suppressWarnings(plot(f, which = "cov")))
  expect_no_error(suppressWarnings(plot(f, which = "mean")))
})
