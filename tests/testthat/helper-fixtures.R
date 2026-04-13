# Shared fixtures for admixr2 tests.
# Sourced automatically by testthat before any test file.

# ---- iniDf builders ----------------------------------------------------------

# Minimal iniDf row constructor.
.ini_row <- function(name, est, lower = -Inf, upper = Inf, fix = FALSE,
                     neta1 = NA_integer_, neta2 = NA_integer_, err = NA_character_) {
  data.frame(name  = name, est = est, lower = lower, upper = upper,
             fix   = fix,  neta1 = neta1, neta2 = neta2, err = err,
             stringsAsFactors = FALSE)
}

# 1 struct theta, 1 proportional sigma, 1 eta.
make_inidf_1eta <- function(cl = log(5), omega = 0.09, sigma = 0.2) {
  rbind(
    .ini_row("tcl",      cl,    -Inf, Inf),
    .ini_row("prop.err", sigma, 0,    Inf,  err = "prop"),
    .ini_row("eta.cl",   omega, NA,   NA,   neta1 = 1L, neta2 = 1L)
  )
}

# 1 struct theta, 1 additive sigma, no etas.
make_inidf_0eta <- function(cl = log(5), sigma = 0.1) {
  rbind(
    .ini_row("tcl",      cl,    -Inf, Inf),
    .ini_row("add.err",  sigma, 0,    Inf,  err = "add")
  )
}

# 2 struct thetas, 1 additive sigma, 2 etas with off-diagonal omega.
make_inidf_2eta <- function(cl = log(5), v = log(20),
                            omega11 = 0.09, omega22 = 0.04,
                            omega21 = 0.018,   # off-diagonal (covariance, not correlation)
                            sigma = 0.1) {
  rbind(
    .ini_row("tcl",      cl,    -Inf, Inf),
    .ini_row("tv",       v,     -Inf, Inf),
    .ini_row("add.err",  sigma, 0,    Inf,  err = "add"),
    .ini_row("eta.cl",   omega11, NA, NA,   neta1 = 1L, neta2 = 1L),
    .ini_row("eta.v",    omega22, NA, NA,   neta1 = 2L, neta2 = 2L),
    .ini_row("eta.cl_v", omega21, NA, NA,   neta1 = 2L, neta2 = 1L)
  )
}

# ---- Known Omega matrices ----------------------------------------------------

Omega_1x1 <- matrix(0.09)
Omega_2x2 <- matrix(c(0.09, 0.018, 0.018, 0.04), 2, 2,
                    dimnames = list(c("eta.cl", "eta.v"), c("eta.cl", "eta.v")))

# ---- Mini study specs (no rxode2 ev needed) ----------------------------------

make_study_var <- function(n_times = 3L, n = 100L) {
  v_diag <- rep(0.01, n_times)
  list(E = rep(1.0, n_times), V = diag(v_diag),
       v_diag = v_diag, n = n,
       times = seq_len(n_times), method = "var")
}

make_study_cov <- function(n_times = 3L, n = 100L) {
  V <- diag(0.01, n_times)
  V[upper.tri(V)] <- V[lower.tri(V)] <- 0.001
  list(E = rep(1.0, n_times), V = V, n = n,
       times = seq_len(n_times), method = "cov")
}

# ---- Synthetic cp_mat with known mean/cov ------------------------------------
# n_sim rows, n_times cols; rows are iid N(mu_true, sigma_true^2).
make_cp_mat <- function(n_sim = 200L, n_times = 3L, seed = 42L) {
  set.seed(seed)
  mu  <- c(1.0, 2.0, 1.5)
  sds <- c(0.1, 0.15, 0.12)
  m   <- matrix(0, n_sim, n_times)
  for (j in seq_len(n_times))
    m[, j] <- rnorm(n_sim, mean = mu[j], sd = sds[j])
  list(mat = m, mu_true = mu, sd_true = sds)
}
