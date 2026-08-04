# Second-order sensitivities: the cross block d2f/(d eta d dir) that lets adfo
# differentiate its structural thetas analytically instead of finite-differencing
# the whole NLL.
#
# Everything here is Tier 2: it compiles rxode2 models.

test_that("order 1 emits no second-order block (admc/adgh must not pay for it)", {
  env <- .int_sens2_setup()$ode
  skip_if(is.null(env) || is.null(env$sm1), "sens model unavailable")

  expect_null(env$sm1$d2_cols)
  expect_equal(env$sm1$order, 1L)
  # the first-order columns are unchanged by the order argument existing
  expect_equal(length(env$sm1$sens_cols), env$pinfo$n_eta)
  expect_true(all(c("tka") %in% names(env$sm1$theta_sens_cols)))
})

test_that("order 2 emits the eta x dir cross block, and only that", {
  env <- .int_sens2_setup()$ode
  skip_if(is.null(env) || is.null(env$sm2), "sens model unavailable")
  skip_if(is.null(env$sm2$d2_cols), "order-2 model unavailable on this rxode2")

  d2 <- env$sm2$d2_cols
  expect_equal(env$sm2$order, 2L)
  # rows = eta directions, columns = ALL directions (etas + unpaired thetas).
  # No theta x theta pair: adfo differentiates the objective once w.r.t. a theta,
  # so that block would be pure cost.
  expect_equal(nrow(d2), env$pinfo$n_eta)
  expect_equal(ncol(d2), length(env$sm2$dirs))
  expect_equal(rownames(d2), env$sm2$eta_dirs)
  expect_equal(colnames(d2), env$sm2$dirs)
  expect_true(all(grepl("^rx_f2_", d2)))
})

test_that("the second-order columns match a central difference of the first-order ones", {
  env <- .int_sens2_setup()$ode
  skip_if(is.null(env) || is.null(env$sm2$d2_cols), "order-2 model unavailable")
  # f1 comes from the variational compartments, f2 from a separate second-order
  # expansion -- so this compares two independent constructions. 1e-4 is loose
  # against the ~2e-6 observed; the reference is the noisy side.
  expect_lt(.sens2_cfd_check(env), 1e-4)
})

test_that("a solved-form linCmt model is DETECTED, whichever marker rxode2 uses", {
  # The gate that decides whether an order-2 request promotes the model. It read
  # `ui$predDf$linCmt`, which rxode2 5.1.4 leaves FALSE for a genuine
  # `cp <- linCmt()` (it marks `ui$.linCmtM` instead) -- so the gate never fired,
  # every order-2 linCmt build returned NULL, and adfo silently kept the
  # finite-difference struct-theta pass the order-2 block exists to replace.
  #
  # Asserted on the PREDICATE and on d2_cols, with NO skip: the test below skips
  # on `is.null(sm2$d2_cols)` to tolerate an old rxode2, and that guard is exactly
  # what hid this -- a feature that never ran looks identical to a feature the
  # platform cannot support. This one must fail rather than skip.
  env <- .int_sens2_setup()$lin
  skip_if(is.null(env) || is.null(env$ui), "linCmt fixture unavailable")

  expect_true(isTRUE(rxode2::testRxLinCmt(env$ui)))
  # ... and the promotion it gates actually happened.
  expect_false(is.null(env$sm2$d2_cols))
})

test_that("linCmt is promoted to ODE form for order 2 and stays solved-form at order 1", {
  env <- .int_sens2_setup()$lin
  skip_if(is.null(env) || is.null(env$sm1), "sens model unavailable")

  # order 1: the fast solved form, differentiated through linCmtB
  expect_null(env$sm1$d2_cols)
  expect_true(isTRUE(env$sm1$is_lincmt))

  skip_if(is.null(env$sm2$d2_cols), "order-2 promotion unavailable on this rxode2")
  # order 2: promoted, so no linCmtB left anywhere in the emitted model
  expect_false(isTRUE(env$sm2$is_lincmt))
  expect_equal(nrow(env$sm2$d2_cols), env$pinfo$n_eta)
  expect_lt(.sens2_cfd_check(env), 1e-4)
})

test_that("the promoted linCmt model reproduces the analytic linCmt prediction", {
  env <- .int_sens2_setup()$lin
  skip_if(is.null(env) || is.null(env$sm2$d2_cols), "order-2 promotion unavailable")

  tms <- .int_sens2_setup()$times
  ev  <- env$studies[[1L]]$ev_full
  r   <- env$sm2$rename_map
  p   <- stats::setNames(as.list(c(log(1.2), log(5), log(20), 0.3, 0, 0)),
                         c(r[["tka"]], r[["tcl"]], r[["tv"]], r[["add.err"]],
                           r[["eta.cl"]], r[["eta.v"]]))
  o2 <- rxode2::rxSolve(env$sm2$mod, params = as.data.frame(p, check.names = FALSE),
                        events = ev, nDisplayProgress = .Machine$integer.max)
  o2 <- o2[o2$time %in% tms, ]
  sl <- rxode2::rxSolve(env$rxMod,
                        data.frame(tka = log(1.2), tcl = log(5), tv = log(20),
                                   add.err = 0.3, eta.cl = 0, eta.v = 0),
                        ev, nDisplayProgress = .Machine$integer.max)
  vl <- sl$ipredSim[sl$time %in% tms]
  # An ODE solve reproducing a closed form: this is solver tolerance, not method
  # error. Observed 1.8e-08 relative.
  expect_lt(max(abs(vl - o2$rx_pred_)) / max(abs(vl)), 1e-5)
})

test_that("adfo's struct-theta gradient is analytic and beats the FD pass it replaces", {
  for (nm in c("ode", "lin")) {
    env <- .int_sens2_setup()[[nm]]
    skip_if(is.null(env) || is.null(env$sm2$d2_cols) || is.null(env$rxMod),
            "order-2 model unavailable")

    args <- list(env$pinfo, env$studies, NULL, env$rxMod, env$ov,
                 env$params_list, 1L)
    g_ana <- do.call(admixr2:::.adfoGrad,
                     c(list(env$p0), replace(args, 3L, list(env$sm2)), list(1e-4)))
    g_fd  <- do.call(admixr2:::.adfoGrad,
                     c(list(env$p0), replace(args, 3L, list(env$sm1)), list(1e-4)))
    h <- 1e-5
    g_ref <- vapply(seq_along(env$p0), function(k) {
      ph <- env$p0; ph[k] <- ph[k] + h
      pm <- env$p0; pm[k] <- pm[k] - h
      (do.call(admixr2:::.adfoNLL, c(list(ph), replace(args, 3L, list(env$sm2)))) -
       do.call(admixr2:::.adfoNLL, c(list(pm), replace(args, 3L, list(env$sm2))))) /
        (2 * h)
    }, double(1))

    ns   <- length(env$pinfo$struct_names)
    rel  <- function(g) max(abs(g[seq_len(ns)] - g_ref[seq_len(ns)]) /
                              pmax(1e-8, abs(g_ref[seq_len(ns)])))
    expect_true(all(is.finite(g_ana)), info = nm)
    # Analytic: observed 2e-07..2e-06. The FD pass it replaces: 8e-04..1e-02.
    expect_lt(rel(g_ana), 1e-4)
    expect_lt(rel(g_ana), rel(g_fd))
  }
})

test_that(".adfoGrad degrades to FD when the cached direction map does not cover a theta", {
  # .dir_of used `theta_dirs[[nm]]`, and `[[` with an unmatched name on an ATOMIC
  # vector throws "subscript out of bounds" rather than returning NULL -- so the
  # `%||% NA_character_` fallback beside it could never fire, and the anyNA()
  # check that turns use_d2 off was unreachable. .adfoGrad is not wrapped in a
  # tryCatch at that point, so the error propagated out of eval_grad_f and killed
  # the whole nloptr run with a bare "subscript out of bounds".
  #
  # Reachable whenever pinfo's unpaired set and the cached model's disagree --
  # they are computed from different objects at different times, and the sens
  # model can be served from a persistent disk cache keyed on the BUILD-time set.
  # Simulated here by emptying the map on a real order-2 model.
  for (nm in c("ode", "lin")) {
    env <- .int_sens2_setup()[[nm]]
    skip_if(is.null(env) || is.null(env$sm2$d2_cols) || is.null(env$rxMod),
            "order-2 model unavailable")

    args <- list(env$pinfo, env$studies, NULL, env$rxMod, env$ov,
                 env$params_list, 1L)
    for (broken in list(character(0), NULL, c(not.a.theta = "THETA_9_"))) {
      sm <- env$sm2
      sm$theta_dirs <- broken
      g <- do.call(admixr2:::.adfoGrad,
                   c(list(env$p0), replace(args, 3L, list(sm)), list(1e-4)))
      expect_true(all(is.finite(g)), info = paste(nm, length(broken)))
      expect_length(g, length(env$p0))
    }
  }
})
