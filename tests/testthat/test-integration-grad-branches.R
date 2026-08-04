# Branches of the MC gradient that a default fit never enters.
#
# admc's default is grad = "sens" with a sensitivity model available, so the
# common-random-number finite-difference blocks -- the ones that perturb an eta
# or a struct theta inside the params frame -- only run as a FALLBACK, and their
# central-difference variants only when grad = "cfd". Neither is reached by any
# end-to-end fit in the suite, which is why the step re-keying done in this
# release (.admGH / .admGH0) landed there untested.
#
# These call the gradients directly with sensModel = NULL and use_central = TRUE
# to reach them, and assert what actually matters: a finite gradient, and the
# central difference agreeing with the forward one.
skip_on_cran()
skip_if_not_installed("rxode2")

test_that(".admGrad runs its CRN fallback under both forward and central FD", {
  env <- .int_grad_setup()
  p <- env$vec$p0
  h <- 1e-3

  fwd <- admixr2:::.admGrad(p, env$pinfo, env$studies, env$z_list, env$rxMod,
                            env$output_var, env$params_list, cores = 1L, h = h,
                            sensModel = NULL, use_central = FALSE)
  ctr <- admixr2:::.admGrad(p, env$pinfo, env$studies, env$z_list, env$rxMod,
                            env$output_var, env$params_list, cores = 1L, h = h,
                            sensModel = NULL, use_central = TRUE)
  expect_length(fwd, length(p))
  expect_true(all(is.finite(fwd)))
  expect_true(all(is.finite(ctr)))
  # Two finite-difference schemes for one derivative: they must agree to within
  # the step's own truncation error, not to machine precision.
  rel <- abs(ctr - fwd) / pmax(abs(ctr), 1e-6)
  expect_lt(max(rel), 0.25)
})

test_that(".admGradBatch runs its CRN fallback under both forward and central FD", {
  env <- .int_grad_setup()
  p <- env$vec$p0
  pl <- list(p, p + 1e-3)
  h <- 1e-3

  fwd <- admixr2:::.admGradBatch(pl, env$pinfo, env$studies, env$z_list, env$rxMod,
                                 env$output_var, env$params_list, cores = 1L, h = h,
                                 sensModel = NULL, use_central = FALSE)
  ctr <- admixr2:::.admGradBatch(pl, env$pinfo, env$studies, env$z_list, env$rxMod,
                                 env$output_var, env$params_list, cores = 1L, h = h,
                                 sensModel = NULL, use_central = TRUE)
  expect_equal(dim(fwd), c(length(pl), length(p)))
  expect_true(all(is.finite(fwd)))
  expect_true(all(is.finite(ctr)))
  # ... and the batch's first row must match the single-vector gradient it
  # duplicates, or the batching has drifted from what it batches.
  one <- admixr2:::.admGrad(p, env$pinfo, env$studies, env$z_list, env$rxMod,
                            env$output_var, env$params_list, cores = 1L, h = h,
                            sensModel = NULL, use_central = FALSE)
  expect_equal(as.numeric(fwd[1L, ]), unname(one), tolerance = 1e-8)
})

test_that(".admGradBatch rejects a step vector of the wrong length", {
  # Same guard as .admGrad's: anything but one step or one per parameter would
  # recycle into the n_sim-row eta perturbations without a warning.
  env <- .int_grad_setup()
  p <- env$vec$p0
  expect_error(
    admixr2:::.admGradBatch(list(p), env$pinfo, env$studies, env$z_list, env$rxMod,
                            env$output_var, env$params_list, cores = 1L,
                            h = rep(1e-3, length(p) - 1L)),
    "one step or one per parameter")
})

test_that(".admNLLGradFD supports central differences", {
  env <- .int_grad_setup()
  p <- env$vec$p0
  fwd <- admixr2:::.admNLLGradFD(p, env$pinfo, env$studies, env$z_list, env$rxMod,
                                 env$output_var, env$params_list, cores = 1L,
                                 h = 1e-3, use_central = FALSE)
  ctr <- admixr2:::.admNLLGradFD(p, env$pinfo, env$studies, env$z_list, env$rxMod,
                                 env$output_var, env$params_list, cores = 1L,
                                 h = 1e-3, use_central = TRUE)
  expect_true(all(is.finite(fwd)))
  expect_true(all(is.finite(ctr)))
  expect_equal(ctr, fwd, tolerance = 0.25)
})

test_that("a non-finite omega is refused before it reaches the solver", {
  # The covariance probe can perturb a parameter to exp(1e5/2) = Inf. The screen
  # tested only that the omega diagonal was positive, and Inf > 0 is TRUE, so Inf
  # went to rxSolve() which integrated garbage and emitted ~190k solver warnings
  # before the caller discarded the result. Both entry points now screen it: the
  # gradient unpacks the optimizer vector itself, so the objective's guard alone
  # did not cover it.
  env <- .int_grad_setup()
  p <- env$vec$p0
  k <- length(env$pinfo$struct_names) + length(env$pinfo$sigma_names) + 1L
  skip_if(k > length(p), "no omega parameter in this fixture")
  p[k] <- 1e5                       # exp(1e5/2) overflows to Inf on unpack

  expect_silent(
    nll <- admixr2:::.admNLL(p, env$pinfo, env$studies, env$z_list, env$rxMod,
                             env$output_var, env$params_list, cores = 1L))
  expect_identical(nll, Inf)
  expect_silent(
    g <- admixr2:::.admGrad(p, env$pinfo, env$studies, env$z_list, env$rxMod,
                            env$output_var, env$params_list, cores = 1L, h = 1e-3,
                            sensModel = env$sensModel))
  expect_true(all(is.na(g)))
})
