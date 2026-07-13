test_that("admControl: cores < workers emits informational message", {
  expect_message(
    admControl(studies = list(), cores = 1L, workers = 4L),
    regexp = "cores.*<.*workers",
    fixed  = FALSE
  )
})

test_that("parallel: remainder >= 0 when cores < effective_workers (no rep() crash)", {
  skip_on_cran()
  skip_if_not_installed("rxode2")
  skip_if_not_installed("nlmixr2")
  skip_if_not_installed("mirai")

  env <- .int_grad_setup()

  # cores = 1, workers = 2 used to crash with remainder = -1
  ctl <- admControl(
    studies    = env$studies,
    n_sim      = 100L,
    maxeval    = 3L,
    seed       = 1L,
    grad       = "sens",
    n_restarts = 2L,
    workers    = 2L,
    cores      = 1L
  )

  expect_no_error(suppressMessages(
    nlmixr2est::nlmixr2(one_cmt_fn, admData(), est = "admc", control = ctl)
  ))
})

# The daemon backend is the same code path on every platform (no fork/PSOCK
# split), so one comparison covers all of them.
test_that("parallel: fit NLL matches sequential (workers = 2, same seed)", {
  skip_on_cran()
  skip_if_not_installed("rxode2")
  skip_if_not_installed("nlmixr2")
  skip_if_not_installed("mirai")

  env <- .int_grad_setup()

  ctl <- admControl(
    studies    = env$studies,
    n_sim      = 200L,
    maxeval    = 5L,
    seed       = 1L,
    grad       = "sens",
    n_restarts = 2L
  )

  fit_seq <- suppressMessages(
    nlmixr2est::nlmixr2(one_cmt_fn, admData(), est = "admc",
                        control = modifyList(ctl, list(workers = 1L)))
  )
  fit_par <- suppressMessages(
    nlmixr2est::nlmixr2(one_cmt_fn, admData(), est = "admc",
                        control = modifyList(ctl, list(workers = 2L)))
  )

  expect_equal(fit_par$objective, fit_seq$objective, tolerance = 1e-2,
               label = "parallel NLL", expected.label = "sequential NLL")
})

test_that("parallel: daemon pool is shut down after the fit", {
  skip_on_cran()
  skip_if_not_installed("rxode2")
  skip_if_not_installed("nlmixr2")
  skip_if_not_installed("mirai")

  env <- .int_grad_setup()

  suppressMessages(
    nlmixr2est::nlmixr2(one_cmt_fn, admData(), est = "admc",
                        control = admControl(studies    = env$studies,
                                             n_sim      = 100L,
                                             maxeval    = 3L,
                                             seed       = 1L,
                                             grad       = "sens",
                                             n_restarts = 2L,
                                             workers    = 2L))
  )

  expect_equal(admixr2:::.adm_worker_env$n, 0L)
})
