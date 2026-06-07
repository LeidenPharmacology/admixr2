test_that(".admParseIniDf extracts correct structure for 1-eta model", {
  pinfo <- admixr2:::.admParseIniDf(make_inidf_1eta())

  expect_equal(pinfo$n_eta, 1L)
  expect_equal(pinfo$eta_names, "eta.cl")
  expect_equal(pinfo$struct_names, "tcl")
  expect_equal(pinfo$sigma_names, "prop.err")
  expect_equal(length(pinfo$omega_par), 1L)     # only one diagonal entry
  expect_equal(pinfo$chol_diag, TRUE)
  expect_equal(pinfo$chol_i, 1L)
  expect_equal(pinfo$chol_j, 1L)
})

test_that(".admParseIniDf extracts correct structure for 2-eta model", {
  pinfo <- admixr2:::.admParseIniDf(make_inidf_2eta())

  expect_equal(pinfo$n_eta, 2L)
  expect_equal(pinfo$struct_names, c("tcl", "tv"))
  expect_equal(pinfo$sigma_names, "add.err")
  # 2 diagonal + 1 off-diagonal
  expect_equal(length(pinfo$omega_par), 3L)
  expect_equal(sum(pinfo$chol_diag), 2L)
  expect_equal(sum(!pinfo$chol_diag), 1L)
  # off-diagonal entry has i > j
  off_idx <- which(!pinfo$chol_diag)
  expect_true(pinfo$chol_i[off_idx] > pinfo$chol_j[off_idx])
})

test_that(".admParseIniDf handles 0-eta model", {
  pinfo <- admixr2:::.admParseIniDf(make_inidf_0eta())

  expect_equal(pinfo$n_eta, 0L)
  expect_equal(length(pinfo$eta_names), 0L)
  expect_equal(length(pinfo$omega_par), 0L)
  expect_equal(pinfo$struct_names, "tcl")
})

test_that("Cholesky round-trip: admBuildOptVec -> admUnpack recovers omega_init", {
  inidf <- make_inidf_2eta()
  pinfo <- admixr2:::.admParseIniDf(inidf)
  vec   <- admixr2:::.admBuildOptVec(pinfo)
  pars  <- admixr2:::.admUnpack(vec$p0, pinfo)

  expect_equal(unname(pars$omega), unname(Omega_2x2), tolerance = 1e-10)
  # L is lower triangular
  expect_true(all(pars$L[upper.tri(pars$L)] == 0))
  # L L' = omega
  expect_equal(pars$L %*% t(pars$L), pars$omega, tolerance = 1e-12)
})

test_that("Diagonal omega parameterization: unit step multiplies Omega_ii by exp(1), not exp(2)", {
  # Regression for the log(Omega_ii) vs log(L_ii) bug documented in CLAUDE.md.
  # With log(Omega_ii) encoding: p -> p+1 should give Omega_ii * exp(1).
  # With the wrong log(L_ii) encoding it would give Omega_ii * exp(2).
  pinfo <- admixr2:::.admParseIniDf(make_inidf_1eta(omega = 0.09))
  vec   <- admixr2:::.admBuildOptVec(pinfo)

  om0 <- admixr2:::.admUnpack(vec$p0, pinfo)$omega[1, 1]

  p1 <- vec$p0
  # perturb only the diagonal omega parameter
  diag_nm <- pinfo$omega_par_names[pinfo$chol_diag][1]
  p1[diag_nm] <- p1[diag_nm] + 1

  om1 <- admixr2:::.admUnpack(p1, pinfo)$omega[1, 1]

  expect_equal(om1 / om0, exp(1), tolerance = 1e-10)   # NOT exp(2)
})

test_that("Diagonal omega initial encoding: omega_par = log(Omega_ii) = 2*log(L_ii)", {
  # Stored as 2*log(L_ii) = log(Omega_ii), so exp(par) = Omega_ii.
  omega_val <- 0.09
  pinfo <- admixr2:::.admParseIniDf(make_inidf_1eta(omega = omega_val))

  diag_par <- pinfo$omega_par[pinfo$chol_diag]
  expect_equal(unname(diag_par), log(omega_val), tolerance = 1e-12)  # = 2*log(sqrt(omega_val))
  # Back-transform: exp(par) = Omega_ii
  expect_equal(exp(unname(diag_par)), omega_val, tolerance = 1e-12)
})

test_that("Sigma encoding: sigma_init = 2*log(sigma_est); round-trip gives sigma_est^2", {
  sigma_est <- 0.2
  pinfo <- admixr2:::.admParseIniDf(make_inidf_1eta(sigma = sigma_est))
  vec   <- admixr2:::.admBuildOptVec(pinfo)
  pars  <- admixr2:::.admUnpack(vec$p0, pinfo)

  expect_equal(unname(pinfo$sigma_init), 2 * log(sigma_est), tolerance = 1e-12)
  expect_equal(unname(pars$sigma_var), sigma_est^2, tolerance = 1e-12)
  expect_equal(sqrt(unname(pars$sigma_var)), sigma_est, tolerance = 1e-12)
})

test_that("Sigma with finite lower bound: lb encoded as 2*log(sigma_lower)", {
  pinfo <- admixr2:::.admParseIniDf(make_inidf_1eta(sigma = 0.2))
  pinfo$sigma_lower <- c("prop.err" = 0.01)
  vec   <- admixr2:::.admBuildOptVec(pinfo)

  sig_lb_idx <- match(pinfo$sigma_names, vec$names)
  expect_equal(unname(vec$lower[sig_lb_idx]), 2 * log(0.01), tolerance = 1e-12)
})

test_that("Fixed parameters absent from opt vector", {
  inidf <- make_inidf_1eta()
  inidf$fix[inidf$name == "tcl"] <- TRUE
  pinfo <- admixr2:::.admParseIniDf(inidf)
  vec   <- admixr2:::.admBuildOptVec(pinfo)

  expect_false("tcl" %in% vec$names)
  expect_equal(pinfo$struct_names, character(0))
})

test_that("Omega parameters have unbounded box constraints", {
  pinfo <- admixr2:::.admParseIniDf(make_inidf_2eta())
  vec   <- admixr2:::.admBuildOptVec(pinfo)

  omega_idx <- match(pinfo$omega_par_names, vec$names)
  expect_true(all(vec$lower[omega_idx] == -Inf))
  expect_true(all(vec$upper[omega_idx] ==  Inf))
})

test_that(".admUnpack handles n_eta = 0 (no omega)", {
  pinfo <- admixr2:::.admParseIniDf(make_inidf_0eta())
  vec   <- admixr2:::.admBuildOptVec(pinfo)
  pars  <- admixr2:::.admUnpack(vec$p0, pinfo)

  expect_null(pars$L)
  expect_equal(dim(pars$omega), c(0L, 0L))
  expect_true("tcl" %in% names(pars$struct))
})

test_that("Struct theta round-trip on optimizer scale", {
  cl_init <- log(5)
  pinfo <- admixr2:::.admParseIniDf(make_inidf_1eta(cl = cl_init))
  vec   <- admixr2:::.admBuildOptVec(pinfo)
  pars  <- admixr2:::.admUnpack(vec$p0, pinfo)

  expect_equal(unname(pars$struct["tcl"]), cl_init, tolerance = 1e-12)
})

test_that("struct_transforms stores metadata only (curEval/low/hi, no closures)", {
  pinfo <- admixr2:::.admParseIniDf(make_inidf_1eta())
  tr <- pinfo$struct_transforms[["tcl"]]

  expect_equal(tr$curEval, "exp")
  expect_true(is.na(tr$low))
  expect_true(is.na(tr$hi))
  expect_null(tr$back_fn)
  expect_null(tr$log_back_fn)
})

test_that(".admBackTransform and .admLogBackTransform correct for exp/log", {
  pinfo <- admixr2:::.admParseIniDf(make_inidf_1eta())
  tr <- pinfo$struct_transforms[["tcl"]]

  expect_equal(admixr2:::.admBackTransform(0, tr), 1)
  expect_equal(admixr2:::.admLogBackTransform(0, tr), 0)
  expect_equal(admixr2:::.admBackTransform(log(4), tr), 4, tolerance = 1e-12)
})

test_that(".admBackTransform: NULL transform returns exp(p)", {
  expect_equal(admixr2:::.admBackTransform(log(7), NULL), 7, tolerance = 1e-12)
})

test_that(".admLogBackTransform: NULL transform returns p unchanged", {
  expect_equal(admixr2:::.admLogBackTransform(2.5, NULL), 2.5)
})

test_that(".admBackTransform: probitInv returns lo + (hi-lo)*pnorm(p)", {
  lo <- 0.1; hi <- 5.0; p <- 0.5
  tr <- list(curEval = "probitInv", low = lo, hi = hi)
  expect_equal(admixr2:::.admBackTransform(p, tr),
               lo + (hi - lo) * pnorm(p), tolerance = 1e-12)
})

test_that(".admLogBackTransform: probit returns log(lo + (hi-lo)*pnorm(p))", {
  lo <- 0.2; hi <- 4.0; p <- -0.3
  tr <- list(curEval = "probit", low = lo, hi = hi)
  expect_equal(admixr2:::.admLogBackTransform(p, tr),
               log(lo + (hi - lo) * pnorm(p)), tolerance = 1e-12)
})

test_that(".admBackTransform: expit returns lo + (hi-lo)*plogis(p)", {
  skip_if_not_installed("rxode2")
  lo <- 0.5; hi <- 5.0; p <- 1.0
  tr <- list(curEval = "expit", low = lo, hi = hi)
  expect_equal(admixr2:::.admBackTransform(p, tr),
               lo + (hi - lo) * plogis(p), tolerance = 1e-10)
})

test_that(".admLogBackTransform: logit returns log(lo + (hi-lo)*plogis(p))", {
  skip_if_not_installed("rxode2")
  lo <- 0.3; hi <- 3.0; p <- -0.5
  tr <- list(curEval = "logit", low = lo, hi = hi)
  expect_equal(admixr2:::.admLogBackTransform(p, tr),
               log(lo + (hi - lo) * plogis(p)), tolerance = 1e-10)
})

test_that(".admComputeScaleC: expit struct theta uses derivative-based scale (no rxode2)", {
  # Build pinfo with expit transform via mock ui$muRefCurEval.
  mock_ui <- list(muRefCurEval = data.frame(
    parameter = "tcl", curEval = "expit", low = 0.5, hi = 5.0,
    stringsAsFactors = FALSE
  ))
  pinfo <- admixr2:::.admParseIniDf(make_inidf_1eta(), mock_ui)
  sc    <- admixr2:::.admComputeScaleC(pinfo)

  p <- pinfo$struct_init[["tcl"]]
  a <- 0.5; b <- 5.0
  expected <- max(exp(p) * (1 + exp(-p))^2 * (a + (b - a) / (1 + exp(-p))) / (b - a), 0.01)
  expect_equal(unname(sc["tcl"]), expected, tolerance = 1e-12)
  expect_true(sc["tcl"] > 0.01)
})

test_that(".admComputeScaleC: probit struct theta uses derivative-based scale", {
  skip_if_not_installed("rxode2")
  mock_ui <- list(muRefCurEval = data.frame(
    parameter = "tcl", curEval = "probitInv", low = 0.1, hi = 5.0,
    stringsAsFactors = FALSE
  ))
  pinfo <- admixr2:::.admParseIniDf(make_inidf_1eta(), mock_ui)
  sc    <- admixr2:::.admComputeScaleC(pinfo)

  # erf(x) = 2*pnorm(x*sqrt(2)) - 1, so 1+erf(p/sqrt(2)) = 2*pnorm(p)
  p <- pinfo$struct_init[["tcl"]]
  a <- 0.1; b <- 5.0
  expected <- max(sqrt(2 * pi) * exp(0.5 * p^2) * (a + (b - a) * pnorm(p)) / (b - a), 0.01)
  expect_equal(unname(sc["tcl"]), expected, tolerance = 1e-10)
  expect_true(sc["tcl"] > 0.01)
})

test_that("struct_has_eta = FALSE for struct thetas when ui = NULL", {
  # Without ui, struct_has_eta = struct_name %in% eta_names.
  # "tcl" is not in eta_names ("eta.cl"), so struct_has_eta["tcl"] = FALSE.
  pinfo <- admixr2:::.admParseIniDf(make_inidf_1eta())
  expect_false(pinfo$struct_has_eta[["tcl"]])
})

test_that("struct_has_eta = TRUE when mock ui muRefDataFrame lists struct theta", {
  inidf   <- make_inidf_1eta()
  mock_ui <- list(muRefDataFrame = data.frame(theta = "tcl", eta = "eta.cl",
                                              stringsAsFactors = FALSE))
  pinfo   <- admixr2:::.admParseIniDf(inidf, mock_ui)
  expect_true(pinfo$struct_has_eta[["tcl"]])
})

test_that("has_kappa = FALSE when all struct thetas appear in ui$muRefDataFrame", {
  inidf   <- make_inidf_1eta()
  mock_ui <- list(muRefDataFrame = data.frame(theta = "tcl", eta = "eta.cl",
                                              stringsAsFactors = FALSE))
  pinfo   <- admixr2:::.admParseIniDf(inidf, mock_ui)
  expect_false(pinfo$has_kappa)
})

test_that("struct_eta_idx maps each eta to paired struct theta index with mock ui", {
  inidf   <- make_inidf_1eta()
  mock_ui <- list(muRefDataFrame = data.frame(theta = "tcl", eta = "eta.cl",
                                              stringsAsFactors = FALSE))
  pinfo   <- admixr2:::.admParseIniDf(inidf, mock_ui)
  # eta.cl pairs with tcl, which is struct index 1
  expect_equal(unname(pinfo$struct_eta_idx), 1L)
})

test_that("sigma_is_prop = TRUE for prop error, FALSE for add error", {
  # sigma_is_prop / sigma_is_lnorm are plain logical vectors (no names)
  pinfo_prop <- admixr2:::.admParseIniDf(make_inidf_1eta())   # err = "prop"
  expect_true(pinfo_prop$sigma_is_prop[[1]])
  expect_false(pinfo_prop$sigma_is_lnorm[[1]])

  pinfo_add <- admixr2:::.admParseIniDf(make_inidf_2eta())    # err = "add"
  expect_false(pinfo_add$sigma_is_prop[[1]])
  expect_false(pinfo_add$sigma_is_lnorm[[1]])
})

test_that("sigma_is_lnorm = TRUE for lnorm error type", {
  # "lnorm" is supported natively and does NOT trigger an approximation warning
  # (only "dlnorm"/"logn"/"dlogn" do)
  inidf <- make_inidf_1eta()
  inidf$err[inidf$name == "prop.err"] <- "lnorm"
  pinfo <- admixr2:::.admParseIniDf(inidf)
  expect_true(pinfo$sigma_is_lnorm[[1]])
  expect_false(pinfo$sigma_is_prop[[1]])
})

test_that(".admComputeScaleC: exp-transform struct and sigma give 1; off-diagonal L scaled by magnitude", {
  pinfo <- admixr2:::.admParseIniDf(make_inidf_2eta())
  sc    <- admixr2:::.admComputeScaleC(pinfo)

  # struct thetas (exp transform): scale = 1
  expect_equal(unname(sc["tcl"]), 1.0)
  expect_equal(unname(sc["tv"]),  1.0)

  # sigma (log-σ² encoding): scale = 1
  expect_equal(unname(sc["add.err"]), 1.0)

  # omega diagonal entries: log(Ω_ii) encoding, scale = 1
  expect_equal(unname(sc["logchol_eta.cl"]), 1.0)
  expect_equal(unname(sc["logchol_eta.v"]),  1.0)

  # off-diagonal: pmax(|L_21_init|, 0.1)
  off_nm <- "chol_eta.v_eta.cl"
  expect_equal(unname(sc[off_nm]), max(abs(pinfo$omega_par[off_nm]), 0.1))
})

test_that(".admBuildOptVec includes scale_c with correct length and names", {
  pinfo <- admixr2:::.admParseIniDf(make_inidf_2eta())
  vec   <- admixr2:::.admBuildOptVec(pinfo)

  expect_true("scale_c" %in% names(vec))
  expect_equal(length(vec$scale_c), length(vec$p0))
  expect_equal(names(vec$scale_c), vec$names)
  expect_true(all(vec$scale_c > 0))
})

test_that(".admComputeScaleC: 0-eta model returns only struct and sigma entries", {
  pinfo <- admixr2:::.admParseIniDf(make_inidf_0eta())
  sc    <- admixr2:::.admComputeScaleC(pinfo)

  expect_equal(length(sc), 2L)  # tcl + add.err
  expect_true(all(sc == 1.0))
})

# ---- .admBackTransform / .admLogBackTransform: identity/unknown ---------------

test_that(".admBackTransform: unknown curEval returns p unchanged", {
  tr <- list(curEval = "identity", low = NA_real_, hi = NA_real_)
  expect_equal(admixr2:::.admBackTransform(3.7, tr), 3.7)
})

test_that(".admLogBackTransform: unknown curEval returns log(p)", {
  tr <- list(curEval = "identity", low = NA_real_, hi = NA_real_)
  expect_equal(admixr2:::.admLogBackTransform(exp(2), tr), 2, tolerance = 1e-12)
})

test_that(".admLogBackTransform: log curEval returns p unchanged", {
  tr <- list(curEval = "log", low = NA_real_, hi = NA_real_)
  expect_equal(admixr2:::.admLogBackTransform(1.5, tr), 1.5)
})

# ---- .admParseIniDf: warning paths for residual error types ------------------

test_that(".admParseIniDf: propT/propF triggers proportional approximation warning", {
  inidf <- make_inidf_1eta()
  inidf$err[inidf$name == "prop.err"] <- "propT"
  expect_warning(
    admixr2:::.admParseIniDf(inidf),
    regexp = "modelled as proportional"
  )
})

test_that(".admParseIniDf: propF triggers proportional approximation warning", {
  inidf <- make_inidf_1eta()
  inidf$err[inidf$name == "prop.err"] <- "propF"
  expect_warning(
    admixr2:::.admParseIniDf(inidf),
    regexp = "modelled as proportional"
  )
})

test_that(".admParseIniDf: norm triggers additive approximation warning", {
  inidf <- make_inidf_1eta()
  inidf$err[inidf$name == "prop.err"] <- "norm"
  expect_warning(
    admixr2:::.admParseIniDf(inidf),
    regexp = "modelled as additive"
  )
})

test_that(".admParseIniDf: dnorm triggers additive approximation warning", {
  inidf <- make_inidf_1eta()
  inidf$err[inidf$name == "prop.err"] <- "dnorm"
  expect_warning(
    admixr2:::.admParseIniDf(inidf),
    regexp = "modelled as additive"
  )
})

test_that(".admParseIniDf: dlnorm triggers lognormal approximation warning", {
  inidf <- make_inidf_1eta()
  inidf$err[inidf$name == "prop.err"] <- "dlnorm"
  expect_warning(
    admixr2:::.admParseIniDf(inidf),
    regexp = "modelled as lognormal"
  )
})

test_that(".admParseIniDf: logn triggers lognormal approximation warning", {
  inidf <- make_inidf_1eta()
  inidf$err[inidf$name == "prop.err"] <- "logn"
  expect_warning(
    admixr2:::.admParseIniDf(inidf),
    regexp = "modelled as lognormal"
  )
})

test_that(".admParseIniDf: unsupported error type triggers warning", {
  inidf <- make_inidf_1eta()
  inidf$err[inidf$name == "prop.err"] <- "weird_err"
  expect_warning(
    admixr2:::.admParseIniDf(inidf),
    regexp = "Unsupported residual error"
  )
})

test_that(".admParseIniDf: propT sigma_is_prop = TRUE (approximated as prop)", {
  inidf <- make_inidf_1eta()
  inidf$err[inidf$name == "prop.err"] <- "propT"
  pinfo <- suppressWarnings(admixr2:::.admParseIniDf(inidf))
  expect_true(pinfo$sigma_is_prop[[1]])
})

test_that(".admParseIniDf: norm sigma_is_prop = FALSE, sigma_is_lnorm = FALSE", {
  inidf <- make_inidf_1eta()
  inidf$err[inidf$name == "prop.err"] <- "norm"
  pinfo <- suppressWarnings(admixr2:::.admParseIniDf(inidf))
  expect_false(pinfo$sigma_is_prop[[1]])
  expect_false(pinfo$sigma_is_lnorm[[1]])
})

test_that(".admParseIniDf: dlnorm sigma_is_lnorm = TRUE (approximated as lnorm)", {
  inidf <- make_inidf_1eta()
  inidf$err[inidf$name == "prop.err"] <- "dlnorm"
  pinfo <- suppressWarnings(admixr2:::.admParseIniDf(inidf))
  expect_true(pinfo$sigma_is_lnorm[[1]])
})
