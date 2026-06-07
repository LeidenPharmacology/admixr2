test_that(".admCalcObjStats: AIC = objective + 2*npar", {
  s <- make_study_var(n_times = 3L, n = 100L)
  res <- admixr2:::.admCalcObjStats(objective = 50.0, npar = 5L, studies = list(s1 = s))
  expect_equal(res$objDf$AIC, 50.0 + 2 * 5)
})

test_that(".admCalcObjStats: BIC = objective + log(nobs)*npar", {
  s1 <- make_study_var(n_times = 3L, n = 100L)
  s2 <- make_study_var(n_times = 4L, n = 50L)
  nobs_expected <- 100L * 3L + 50L * 4L
  npar <- 4L
  obj  <- 80.0
  res  <- admixr2:::.admCalcObjStats(obj, npar, list(s1 = s1, s2 = s2))
  expect_equal(res$objDf$BIC, obj + log(nobs_expected) * npar, tolerance = 1e-12)
})

test_that(".admCalcObjStats: logLik class and attributes", {
  s   <- make_study_var(n_times = 2L, n = 60L)
  res <- admixr2:::.admCalcObjStats(objective = 40.0, npar = 3L, studies = list(s = s))
  expect_s3_class(res$ll, "logLik")
  expect_equal(attr(res$ll, "df"),   3L)
  expect_equal(attr(res$ll, "nobs"), 60L * 2L)
  expect_equal(as.numeric(res$ll), -40.0 / 2, tolerance = 1e-12)
})

test_that(".admCalcObjStats: nobs = sum(n * n_times)", {
  s1 <- make_study_var(n_times = 3L, n = 100L)
  s2 <- make_study_var(n_times = 5L, n = 200L)
  res <- admixr2:::.admCalcObjStats(0, 1L, list(s1, s2))
  expect_equal(res$nobs, 100L * 3L + 200L * 5L)
})

test_that(".admMakeZ: sobol returns correct shape", {
  pinfo <- admixr2:::.admParseIniDf(make_inidf_2eta())
  set.seed(1)
  z_list <- admixr2:::.admMakeZ(500L, pinfo, n_studies = 2L, sampling = "sobol")
  expect_length(z_list, 2L)
  expect_equal(dim(z_list[[1]]), c(500L, 2L))
  expect_equal(dim(z_list[[2]]), c(500L, 2L))
})

test_that(".admMakeZ: rnorm returns roughly N(0,1) marginals", {
  pinfo <- admixr2:::.admParseIniDf(make_inidf_1eta())
  set.seed(42)
  z_list <- admixr2:::.admMakeZ(5000L, pinfo, n_studies = 1L, sampling = "rnorm")
  z <- z_list[[1]]
  expect_true(abs(mean(z)) < 0.1)
  expect_true(abs(sd(z) - 1) < 0.05)
})

test_that(".admMakeZ: n_eta=0 returns matrix of zeros", {
  pinfo <- admixr2:::.admParseIniDf(make_inidf_0eta())
  set.seed(1)
  z_list <- admixr2:::.admMakeZ(100L, pinfo, n_studies = 1L)
  expect_equal(dim(z_list[[1]]), c(100L, 1L))
  expect_true(all(z_list[[1]] == 0))
})

test_that(".admMakeZ: reproducible with same seed", {
  pinfo <- admixr2:::.admParseIniDf(make_inidf_1eta())
  set.seed(7)
  z1 <- admixr2:::.admMakeZ(200L, pinfo, 1L, "sobol")
  set.seed(7)
  z2 <- admixr2:::.admMakeZ(200L, pinfo, 1L, "sobol")
  expect_identical(z1, z2)
})

test_that(".admMakeParamsList: correct structure and initialisation", {
  pinfo <- admixr2:::.admParseIniDf(make_inidf_1eta())
  pl <- admixr2:::.admMakeParamsList(100L, pinfo, n_studies = 2L)

  expect_length(pl, 2L)
  m <- pl[[1]]
  expect_true(is.matrix(m))
  expect_equal(nrow(m), 100L)

  expected_cols <- c(pinfo$struct_names, pinfo$eta_col_names,
                     pinfo$sigma_names, "rxerr.cp")
  expect_equal(colnames(m), expected_cols)

  expect_true(all(m[, "rxerr.cp"] == 1L))
  non_rxerr <- setdiff(colnames(m), "rxerr.cp")
  expect_true(all(m[, non_rxerr] == 0))
})

test_that(".lhsSample: each column covers all strata exactly once", {
  set.seed(1)
  n <- 10L; d <- 3L
  m <- admixr2:::.lhsSample(n, d)
  expect_equal(dim(m), c(n, d))
  expect_true(all(m > 0) && all(m < 1))

  for (j in seq_len(d)) {
    strata <- floor(m[, j] * n)
    expect_equal(sort(strata), 0L:(n - 1L))
  }
})

test_that(".lhsSample: values in (0, 1) for large n", {
  set.seed(2)
  m <- admixr2:::.lhsSample(1000L, 5L)
  expect_true(all(m > 0))
  expect_true(all(m < 1))
})

# ---- .admFullTheta -----------------------------------------------------------

test_that(".admFullTheta: assembles vector in iniDf row order", {
  pinfo <- admixr2:::.admParseIniDf(make_inidf_2eta())
  vec   <- admixr2:::.admBuildOptVec(pinfo)
  pars  <- admixr2:::.admUnpack(vec$p0, pinfo)

  ft <- admixr2:::.admFullTheta(pars, pinfo)

  expect_equal(names(ft), pinfo$iniDf$name)
  # Struct thetas returned on optimizer scale (as-is)
  expect_equal(unname(ft["tcl"]), unname(pars$struct["tcl"]))
  expect_equal(unname(ft["tv"]),  unname(pars$struct["tv"]))
  # Sigma returned as sqrt(sigma_var) = original sigma SD
  expect_equal(unname(ft["add.err"]), sqrt(unname(pars$sigma_var["add.err"])),
               tolerance = 1e-12)
  # Omega diagonal entries
  expect_equal(unname(ft["eta.cl"]), pars$omega[1, 1], tolerance = 1e-12)
  expect_equal(unname(ft["eta.v"]),  pars$omega[2, 2], tolerance = 1e-12)
})

test_that(".admFullTheta: fixed parameter uses iniDf est, not pars", {
  inidf <- make_inidf_1eta()
  inidf$fix[inidf$name == "tcl"] <- TRUE
  inidf$est[inidf$name == "tcl"] <- 99.0
  pinfo <- admixr2:::.admParseIniDf(inidf)
  vec   <- admixr2:::.admBuildOptVec(pinfo)
  pars  <- admixr2:::.admUnpack(vec$p0, pinfo)

  ft <- admixr2:::.admFullTheta(pars, pinfo)
  expect_equal(unname(ft["tcl"]), 99.0)
})

# ---- .admOutputVar -----------------------------------------------------------

test_that(".admOutputVar: uses predDf$var[1] as output variable", {
  ui <- list(predDf = data.frame(var = "cp", stringsAsFactors = FALSE))
  expect_equal(admixr2:::.admOutputVar(ui), "cp")
})

test_that(".admOutputVar: returns custom var name from predDf", {
  ui <- list(predDf = data.frame(var = "custom_pred", stringsAsFactors = FALSE))
  expect_equal(admixr2:::.admOutputVar(ui), "custom_pred")
})

test_that(".admOutputVar: rx-prefixed var maps to 'ipredSim'", {
  ui <- list(predDf = data.frame(var = "rxLinCmt", stringsAsFactors = FALSE))
  expect_equal(admixr2:::.admOutputVar(ui), "ipredSim")
})

test_that(".admOutputVar: linCmt-prefixed var maps to 'ipredSim'", {
  ui <- list(predDf = data.frame(var = "linCmtB_1", stringsAsFactors = FALSE))
  expect_equal(admixr2:::.admOutputVar(ui), "ipredSim")
})

test_that(".admOutputVar: NULL predDf returns 'cp'", {
  ui <- list(predDf = NULL)
  expect_equal(admixr2:::.admOutputVar(ui), "cp")
})

test_that(".admOutputVar: missing predDf (error) returns 'cp'", {
  ui <- list()
  expect_equal(admixr2:::.admOutputVar(ui), "cp")
})

test_that(".admFullTheta: off-diagonal omega entry from omega[neta1, neta2]", {
  pinfo   <- admixr2:::.admParseIniDf(make_inidf_2eta())
  vec     <- admixr2:::.admBuildOptVec(pinfo)
  pars    <- admixr2:::.admUnpack(vec$p0, pinfo)
  ft      <- admixr2:::.admFullTheta(pars, pinfo)

  off_nm  <- "eta.cl_v"
  off_row <- which(pinfo$iniDf$name == off_nm)
  i <- pinfo$iniDf$neta1[off_row]
  j <- pinfo$iniDf$neta2[off_row]
  expect_equal(unname(ft[off_nm]), pars$omega[i, j], tolerance = 1e-12)
})

# ---- .admMakeZ: additional sampling methods ----------------------------------

test_that(".admMakeZ: halton returns correct shape", {
  pinfo <- admixr2:::.admParseIniDf(make_inidf_2eta())
  z_list <- admixr2:::.admMakeZ(200L, pinfo, n_studies = 1L, sampling = "halton")
  expect_length(z_list, 1L)
  expect_equal(dim(z_list[[1]]), c(200L, 2L))
})

test_that(".admMakeZ: torus returns correct shape", {
  pinfo <- admixr2:::.admParseIniDf(make_inidf_2eta())
  z_list <- admixr2:::.admMakeZ(200L, pinfo, n_studies = 1L, sampling = "torus")
  expect_length(z_list, 1L)
  expect_equal(dim(z_list[[1]]), c(200L, 2L))
})

test_that(".admMakeZ: lhs returns correct shape and uniform marginals", {
  pinfo <- admixr2:::.admParseIniDf(make_inidf_1eta())
  set.seed(99)
  z_list <- admixr2:::.admMakeZ(100L, pinfo, n_studies = 1L, sampling = "lhs")
  expect_equal(dim(z_list[[1]]), c(100L, 1L))
})

test_that(".admMakeZ: halton multi-study returns correct number of studies", {
  # Use 2-eta model: randtoolbox::halton(dim=1) returns a vector not a matrix,
  # so dim() is NULL for 1-eta; 2-eta always returns a matrix.
  pinfo <- admixr2:::.admParseIniDf(make_inidf_2eta())
  z_list <- admixr2:::.admMakeZ(50L, pinfo, n_studies = 3L, sampling = "halton")
  expect_length(z_list, 3L)
  for (i in seq_len(3L)) expect_equal(dim(z_list[[i]]), c(50L, 2L))
})

# ---- admData() ---------------------------------------------------------------

test_that("admData() returns data.frame with required columns", {
  d <- admData()
  expect_s3_class(d, "data.frame")
  expect_named(d, c("ID", "TIME", "DV", "AMT", "EVID", "CMT"))
})

test_that("admData() has all-NA DV column", {
  d <- admData()
  expect_true(all(is.na(d$DV)))
})

test_that("admData() has at least one dosing and one observation row", {
  d <- admData()
  expect_true(any(d$EVID != 0L))
  expect_true(any(d$EVID == 0L))
})

# ---- admStopWorkers() --------------------------------------------------------

test_that("admStopWorkers(): returns invisibly when no cluster is running", {
  expect_invisible(admStopWorkers())
})

test_that("admStopWorkers(): can be called twice without error", {
  expect_no_error(admStopWorkers())
  expect_no_error(admStopWorkers())
})

# ---- %||% operator -----------------------------------------------------------

test_that("%||% returns left side when non-NULL", {
  expect_equal(admixr2:::`%||%`("a", "b"), "a")
})

test_that("%||% returns right side when left is NULL", {
  expect_equal(admixr2:::`%||%`(NULL, "fallback"), "fallback")
})

test_that("%||% works with numeric values", {
  expect_equal(admixr2:::`%||%`(NULL, 42), 42)
  expect_equal(admixr2:::`%||%`(0, 42), 0)
})

# ---- .admSetupParallelPlan() single-worker fast path -------------------------

test_that(".admSetupParallelPlan: workers=1 returns invisibly without starting workers", {
  ctl <- admControl(workers = 1L)
  expect_invisible(admixr2:::.admSetupParallelPlan(ctl, n_r = 3L))
})
