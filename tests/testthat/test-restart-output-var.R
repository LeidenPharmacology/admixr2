test_that("restart workers expose output_var argument with cp default", {
  expect_true("output_var" %in% names(formals(admixr2:::.adfoRestartWorker)))
  expect_true("output_var" %in% names(formals(admixr2:::.admRestartWorker)))
  expect_true("output_var" %in% names(formals(admixr2:::.adirmcRestartWorker)))

  expect_identical(eval(formals(admixr2:::.adfoRestartWorker)$output_var), "cp")
  expect_identical(eval(formals(admixr2:::.admRestartWorker)$output_var), "cp")
  expect_identical(eval(formals(admixr2:::.adirmcRestartWorker)$output_var), "cp")
})

test_that("restart estimators pass detected output_var into restart workers", {
  adfo_txt <- paste(deparse(body(admixr2:::nlmixr2Est.adfo)), collapse = "\n")
  admc_txt <- paste(deparse(body(admixr2:::nlmixr2Est.admc)), collapse = "\n")
  irmc_txt <- paste(deparse(body(admixr2:::nlmixr2Est.adirmc)), collapse = "\n")

  pat <- "(?s)extra_args\\s*=\\s*list\\(.*?output_var\\s*=\\s*output_var"

  expect_true(grepl(pat, adfo_txt, perl = TRUE))
  expect_true(grepl(pat, admc_txt, perl = TRUE))
  expect_true(grepl(pat, irmc_txt, perl = TRUE))
})

test_that("restart workers run with non-cp output_var", {
  skip_on_cran()
  skip_if_not_installed("rxode2")
  skip_if_not_installed("nloptr")

  env <- .int_lincmt_setup(n_sim = 32L, seed = 123L)
  vec <- env$vec

  adfo_res <- admixr2:::.adfoRestartWorker(
    restart_id = 1L, p_init = vec$p0, ui_lstExpr = NULL, pinfo = env$pinfo,
    ov_lower = vec$lower, ov_upper = vec$upper, scale_c = vec$scale_c,
    studies = env$studies, n_sim = 1L, seed = 123L,
    algorithm = "NLOPT_LN_BOBYQA", ftol_rel = 1e-8, maxeval = 1L,
    use_grad = FALSE, grad_h = 1e-3, grad_bounds = 1,
    output_var = env$output_var, sampling = "sobol",
    use_central = FALSE, use_pure_fd = FALSE,
    print_progress = FALSE, print = 0L, cores = 1L, no_lock = TRUE,
    rxMod_direct = env$rxMod, sensModel_direct = env$sensModel
  )
  expect_true(is.finite(adfo_res$objective))

  admc_res <- admixr2:::.admRestartWorker(
    restart_id = 1L, p_init = vec$p0, ui_lstExpr = NULL, pinfo = env$pinfo,
    ov_lower = vec$lower, ov_upper = vec$upper, scale_c = vec$scale_c,
    studies = env$studies, n_sim = 32L, seed = 123L,
    algorithm = "NLOPT_LN_BOBYQA", ftol_rel = 1e-8, maxeval = 1L,
    use_grad = FALSE, grad_h = 1e-3, grad_bounds = 1,
    output_var = env$output_var, sampling = "sobol", use_central = FALSE,
    print_progress = FALSE, print = 0L, cores = 1L, no_lock = TRUE,
    rxMod_direct = env$rxMod, sensModel_direct = env$sensModel
  )
  expect_true(is.finite(admc_res$objective))

  irmc_res <- admixr2:::.adirmcRestartWorker(
    restart_id = 1L, p_init = vec$p0, ui_lstExpr = NULL, pinfo = env$pinfo,
    ov_lower = vec$lower, ov_upper = vec$upper,
    studies = env$studies, n_sim = 32L, seed = 123L,
    phases = 1L, outer_iter = 1L, maxeval = 1L, ftol_rel = 1e-8,
    algorithm = "NLOPT_LN_BOBYQA", omega_expansion = 2, convcrit = 1e-6,
    max_worse = 1L, grad_mode = "none",
    output_var = env$output_var, kappa_method = "exact", sampling = "sobol",
    print_progress = FALSE, print = 0L, cores = 1L, no_lock = TRUE,
    rxMod_direct = env$rxMod
  )
  expect_true(is.finite(irmc_res$objective))
})
