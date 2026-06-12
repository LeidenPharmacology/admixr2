# Tier 1 unit tests for the nlmixr2 traceplot() bridge.
#
# admixr2 fits populate `fit$env$parHistData` so nlmixr2's traceplot() generic
# (traceplot.nlmixr2FitCore, which reads fit$parHistStacked) works natively.
# These tests exercise the two pure-R helpers behind that bridge without
# building a real rxode2 fit:
#   .admTraceDisplaySpec()  - shared display names + back-transforms
#   .admBuildParHistData()  - best-restart, natural-scale parameter history

# ---- .admTraceDisplaySpec ----------------------------------------------------

test_that("display spec: NULL pinfo or par_names returns NULL", {
  expect_null(.admTraceDisplaySpec(NULL, c("tcl")))
  expect_null(.admTraceDisplaySpec(list(), NULL))
})

test_that("display spec: omega diag labelled V(eta), off-diag labelled eta,eta", {
  pinfo     <- .admParseIniDf(make_inidf_2eta())
  par_names <- names(.admBuildOptVec(pinfo)$p0)
  spec      <- .admTraceDisplaySpec(pinfo, par_names, make_inidf_2eta())

  disp <- unlist(spec$disp_nms)
  expect_true("V(eta.cl)" %in% disp)
  expect_true("V(eta.v)"  %in% disp)
  expect_true("eta.v,eta.cl" %in% disp)   # off-diagonal Cholesky entry
  expect_true(all(c("tcl", "tv") %in% disp))
})

test_that("display spec: back-transforms map optimizer scale to natural scale", {
  pinfo     <- .admParseIniDf(make_inidf_1eta())
  par_names <- names(.admBuildOptVec(pinfo)$p0)
  spec      <- .admTraceDisplaySpec(pinfo, par_names, make_inidf_1eta())

  # struct theta: exp() back-transform
  tcl_fn <- spec$back_fns[["tcl"]]
  expect_equal(tcl_fn(log(5)), 5)

  # omega diagonal stored as log(Omega_ii) -> exp() recovers variance
  om_nm <- pinfo$omega_par_names[pinfo$chol_diag][1]
  expect_equal(spec$back_fns[[om_nm]](log(0.09)), 0.09)

  # sigma stored as log(sigma^2) -> exp(v/2) recovers SD
  sig_nm <- setdiff(par_names, c("tcl", pinfo$omega_par_names))[1]
  expect_equal(spec$back_fns[[sig_nm]](2 * log(0.2)), 0.2)
})

test_that("display spec: param_order follows iniDf row order", {
  pinfo     <- .admParseIniDf(make_inidf_2eta())
  par_names <- names(.admBuildOptVec(pinfo)$p0)
  spec      <- .admTraceDisplaySpec(pinfo, par_names, make_inidf_2eta())
  # struct thetas first, in iniDf order, then omega entries
  expect_equal(spec$param_order[1:3], c("tcl", "tv", "add.err"))
})

# ---- .admBuildParHistData ----------------------------------------------------

.mk_traces <- function(par_names, niter = 4L, seed = 1L) {
  set.seed(seed)
  np <- length(par_names)
  list(
    list(restart_id = 1L, nll_trace = c(100, 95, 93, 92),
         par_trace = matrix(rnorm(niter * np), niter, np)),
    list(restart_id = 2L, nll_trace = c(100, 90, 85, 80),   # best (lowest final)
         par_trace = matrix(rnorm(niter * np), niter, np))
  )
}

test_that("parHistData: NULL when no traces", {
  expect_null(.admBuildParHistData(NULL, c("tcl"), list(iniDf = make_inidf_1eta())))
  expect_null(.admBuildParHistData(list(), c("tcl"), list(iniDf = make_inidf_1eta())))
})

test_that("parHistData: NULL when all final NLLs are NA", {
  traces <- list(list(restart_id = 1L, nll_trace = numeric(0), par_trace = NULL))
  expect_null(.admBuildParHistData(traces, c("tcl"), list(iniDf = make_inidf_1eta())))
})

test_that("parHistData: has type/iter columns plus one column per parameter", {
  ini       <- make_inidf_2eta()
  pinfo     <- .admParseIniDf(ini)
  par_names <- names(.admBuildOptVec(pinfo)$p0)
  ph <- .admBuildParHistData(.mk_traces(par_names), par_names, list(iniDf = ini))

  expect_s3_class(ph, "data.frame")
  expect_true(all(c("type", "iter") %in% names(ph)))
  expect_equal(unique(ph$type), "Unscaled")
  expect_equal(ph$iter, 1:4)
  # one data column per parameter (display-named), plus type + iter
  expect_equal(ncol(ph), length(par_names) + 2L)
  expect_true("V(eta.cl)" %in% names(ph))
})

test_that("parHistData: selects the best restart (lowest final NLL)", {
  ini       <- make_inidf_1eta()
  pinfo     <- .admParseIniDf(ini)
  par_names <- names(.admBuildOptVec(pinfo)$p0)
  traces    <- .mk_traces(par_names)
  ph <- .admBuildParHistData(traces, par_names, list(iniDf = ini))

  # restart 2 is best; its tcl column back-transformed = exp(par_trace[,1])
  best_tcl_raw <- traces[[2]]$par_trace[, 1]
  expect_equal(ph[["tcl"]], exp(best_tcl_raw))
})

test_that("parHistData: columns follow iniDf facet order", {
  ini       <- make_inidf_2eta()
  pinfo     <- .admParseIniDf(ini)
  par_names <- names(.admBuildOptVec(pinfo)$p0)
  ph <- .admBuildParHistData(.mk_traces(par_names), par_names, list(iniDf = ini))
  data_cols <- setdiff(names(ph), c("type", "iter"))
  expect_equal(data_cols[1:3], c("tcl", "tv", "add.err"))
})

test_that("parHistData: NULL when par_trace col count disagrees with par_names", {
  ini       <- make_inidf_1eta()
  par_names <- names(.admBuildOptVec(.admParseIniDf(ini))$p0)
  bad <- list(list(restart_id = 1L, nll_trace = c(2, 1),
                   par_trace = matrix(0, 2L, length(par_names) + 1L)))
  expect_null(.admBuildParHistData(bad, par_names, list(iniDf = ini)))
})

test_that("parHistData feeds nlmixr2's parHistStacked contract (iter/par/val)", {
  # Mirror nlmixr2est:::.parHistCalc + nmObjGet.parHistStacked stacking so we
  # assert the shape traceplot.nlmixr2FitCore actually consumes. Runs even when
  # nlmixr2est is unavailable; the live integration check below complements it.
  ini       <- make_inidf_2eta()
  pinfo     <- .admParseIniDf(ini)
  par_names <- names(.admBuildOptVec(pinfo)$p0)
  ph <- .admBuildParHistData(.mk_traces(par_names), par_names, list(iniDf = ini))

  unscaled <- ph[ph$type == "Unscaled", names(ph) != "type"]
  stacked  <- data.frame(iter = unscaled$iter,
                         stack(unscaled[, names(unscaled) != "iter"]))
  names(stacked) <- sub("values", "val", sub("ind", "par", names(stacked)))

  expect_true(all(c("iter", "par", "val") %in% names(stacked)))
  expect_equal(nrow(stacked), 4L * length(par_names))
  expect_true(is.numeric(stacked$val))
})

test_that("parHistData is consumed by the real nlmixr2est parHistStacked getter", {
  # Live integration: feed our parHistData through nlmixr2est's actual
  # nmObjGet.parHistStacked (the path fit$parHistStacked -> traceplot uses)
  # rather than a re-implementation, so a contract drift in nlmixr2est is caught.
  skip_if_not_installed("nlmixr2est")
  getter <- tryCatch(getFromNamespace("nmObjGet.parHistStacked", "nlmixr2est"),
                     error = function(e) NULL)
  skip_if(is.null(getter), "nlmixr2est:::nmObjGet.parHistStacked unavailable")

  ini       <- make_inidf_2eta()
  pinfo     <- .admParseIniDf(ini)
  par_names <- names(.admBuildOptVec(pinfo)$p0)
  ph <- .admBuildParHistData(.mk_traces(par_names), par_names, list(iniDf = ini))

  # Minimal stand-in for a fit object: the getter reads x[[1]]$env$parHistData.
  e <- new.env(); e$parHistData <- ph
  stacked <- getter(list(list(env = e)))

  expect_s3_class(stacked, "data.frame")
  expect_true(all(c("iter", "par", "val") %in% names(stacked)))
  expect_equal(nrow(stacked), 4L * length(par_names))
  expect_true(is.numeric(stacked$val))
  # Display names survive the round-trip (e.g. omega diagonal V(eta.cl)).
  expect_true("V(eta.cl)" %in% as.character(stacked$par))
})
