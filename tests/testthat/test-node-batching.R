# Tier 1 -- node-batched forward/gradient simulation (no rxode2; mock solver).
# Confirms the stacked single-rxSolve path reproduces a per-node solve exactly.

.nb_setup <- function() {
  times   <- c(1, 2, 3); n_sim <- 4L
  set.seed(7)
  eta_mat <- matrix(rnorm(n_sim * 2), n_sim, 2,
                    dimnames = list(NULL, c("eta.cl", "eta.v")))
  struct  <- c(tcl = 1.0, tv = 2.0, tka = 0.5)
  rxMod   <- list(params = c("tcl","tv","tka","eta.cl","eta.v","wt","rxerr.cp"))
  ev_df   <- data.frame(time = times, amt = 0)
  studies <- lapply(c(60, 70, 80), function(w)
    list(times = times, cov = c(wt = w), ev_full = ev_df))
  z_list  <- replicate(3, matrix(0, n_sim, 2), simplify = FALSE)
  list(times = times, n_sim = n_sim, eta_mat = eta_mat, struct = struct,
       rxMod = rxMod, studies = studies, z_list = z_list,
       pinfo = list(struct_names = c("tcl","tv","tka"), sigma_names = character(0)))
}

# mock simulation: cp = tcl + 5*tka + 10*eta.cl + 3*eta.v + 0.01*wt + t (subject-major)
.nb_mock <- function(times) function(rxMod, params, events, cores, nDisplayProgress) {
  df <- params; N <- nrow(df); cp <- numeric(N * length(times)); k <- 0L
  for (r in seq_len(N)) for (t in times) { k <- k + 1L
    cp[k] <- df$tcl[r] + 5*df$tka[r] + 10*df$eta.cl[r] + 3*df$eta.v[r] + 0.01*df$wt[r] + t }
  data.frame(time = rep(times, N), cp = cp)
}
.nb_f <- function(s, ecl, ev, w, t) s["tcl"] + 5*s["tka"] + 10*ecl + 3*ev + 0.01*w + t
.nb_refcp <- function(d, w) outer(seq_len(d$n_sim), d$times,
  Vectorize(function(j, t) .nb_f(d$struct, d$eta_mat[j,1], d$eta_mat[j,2], w, t)))

test_that(".adm_nodes_batchable detects shared design", {
  d <- .nb_setup()
  expect_true(admixr2:::.adm_nodes_batchable(d$studies))
  # Different z draws no longer block batching (only ev/times matter)
  st <- d$studies; st[[2]]$times <- c(1, 2, 4)
  expect_false(admixr2:::.adm_nodes_batchable(st))                # different times
  expect_false(admixr2:::.adm_nodes_batchable(d$studies[1]))      # single study
})

test_that(".adm_sim_nodes_batched == per-node solve (NLL path)", {
  d <- .nb_setup(); mk <- .nb_mock(d$times)
  # Same eta_mat for all nodes (CRN / Taylor case)
  eta_list_crn <- replicate(3, d$eta_mat, simplify = FALSE)
  cps <- admixr2:::.adm_sim_nodes_batched(d$rxMod, d$struct, character(0), eta_list_crn,
                                          d$studies, "cp", 1L, solve_fn = mk)
  for (i in 1:3) expect_equal(unname(cps[[i]]), unname(.nb_refcp(d, c(60,70,80)[i])))
  expect_false(isTRUE(all.equal(cps[[1]], cps[[3]])))   # covariate actually injected

  # Per-node independent draws (GL case): each node uses its own eta_mat
  set.seed(42)
  eta_list_gl <- lapply(1:3, function(i)
    matrix(rnorm(d$n_sim * 2), d$n_sim, 2, dimnames = list(NULL, c("eta.cl","eta.v"))))
  cps_gl <- admixr2:::.adm_sim_nodes_batched(d$rxMod, d$struct, character(0), eta_list_gl,
                                              d$studies, "cp", 1L, solve_fn = mk)
  d2 <- d; d2$eta_mat <- eta_list_gl[[2]]
  expect_equal(unname(cps_gl[[2]]), unname(.nb_refcp(d2, 70)))   # node 2 uses its own eta
})

test_that(".adm_grad_presim FD path reproduces cp / dpred / unpaired hi", {
  d <- .nb_setup(); h <- 1e-3; mk <- .nb_mock(d$times)
  eta_list <- replicate(3, d$eta_mat, simplify = FALSE)
  ps <- admixr2:::.adm_grad_presim(
    d$studies, eta_list, list(struct = d$struct), d$pinfo, d$rxMod,
    sensModel = NULL, output_var = "cp", cores = 1L, h = h, use_central = FALSE,
    n_eta = 2L, n_sim = d$n_sim, eta_col_names = c("eta.cl","eta.v"),
    unpaired_k = 3L, solve_fn = mk)
  for (i in 1:3) {
    w <- c(60, 70, 80)[i]
    expect_equal(unname(ps$cp[[i]]), unname(.nb_refcp(d, w)))
    expect_equal(ps$dpred[[i]][[1]], matrix(10, d$n_sim, length(d$times)), tolerance = 1e-3)
    expect_equal(ps$dpred[[i]][[2]], matrix(3,  d$n_sim, length(d$times)), tolerance = 1e-3)
    expect_equal(unname(ps$hi[[i]][[1]]), unname(.nb_refcp(d, w) + 5*h), tolerance = 1e-9)
  }
  expect_null(ps$lo)   # forward FD

  # GL case: independent eta draws per node — each node uses its own eta_mat
  set.seed(42)
  eta_gl <- lapply(1:3, function(i)
    matrix(rnorm(d$n_sim * 2), d$n_sim, 2, dimnames = list(NULL, c("eta.cl","eta.v"))))
  ps2 <- admixr2:::.adm_grad_presim(
    d$studies, eta_gl, list(struct = d$struct), d$pinfo, d$rxMod,
    sensModel = NULL, output_var = "cp", cores = 1L, h = h, use_central = FALSE,
    n_eta = 2L, n_sim = d$n_sim, eta_col_names = c("eta.cl","eta.v"),
    unpaired_k = 3L, solve_fn = mk)
  d2 <- d; d2$eta_mat <- eta_gl[[2]]
  expect_equal(unname(ps2$cp[[2]]), unname(.nb_refcp(d2, 70)))
})

test_that(".adm_grad_presim central-FD path adds the lo perturbation", {
  d <- .nb_setup(); h <- 1e-3; mk <- .nb_mock(d$times)
  eta_list <- replicate(3, d$eta_mat, simplify = FALSE)
  ps <- admixr2:::.adm_grad_presim(
    d$studies, eta_list, list(struct = d$struct), d$pinfo, d$rxMod,
    sensModel = NULL, output_var = "cp", cores = 1L, h = h, use_central = TRUE,
    n_eta = 2L, n_sim = d$n_sim, eta_col_names = c("eta.cl","eta.v"),
    unpaired_k = 3L, solve_fn = mk)
  for (i in 1:3) {
    w <- c(60, 70, 80)[i]
    expect_equal(unname(ps$hi[[i]][[1]]), unname(.nb_refcp(d, w) + 5*h), tolerance = 1e-9)
    expect_equal(unname(ps$lo[[i]][[1]]), unname(.nb_refcp(d, w) - 5*h), tolerance = 1e-9)
  }
})
