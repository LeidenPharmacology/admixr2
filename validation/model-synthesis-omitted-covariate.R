## Combining models that DO NOT SHARE STRUCTURE, using the covariate JOINT.
##
## Two published models. NEITHER contains CRCL -- both report only a weight
## effect. Their populations differ in how strongly weight and CRCL are
## correlated (rho_A = 0.20, rho_B = 0.70), which is the sort of thing a baseline
## characteristics table gives you for free.
##
## Each source's weight coefficient is therefore a MARGINAL effect that has
## absorbed part of the CRCL effect through that correlation:
##
##   a_s = b_WT + b_CRCL * rho_s * sd_CRCL / sd_WT
##
## Two sources, two different rho, two equations -- so b_WT and b_CRCL SEPARATE,
## and the unified model can carry a covariate that no source model contains.
## The same omission also inflated each source's IIV by b_CRCL^2 sd^2 (1-rho^2),
## which the unified model should hand back.
##
## This works only if each source's blocks are projected over its own JOINT
## covariate distribution, so that within a weight stratum CRCL follows its
## CONDITIONAL distribution. The "naive" arm below draws CRCL from its
## unconditional margin -- destroying the dependence -- and is the control.
##
## Data are exact quadrature moments, no sampling noise, so a displaced argmin is
## a property of the construction.

TCL <- log(1); TV <- log(10); DOSE <- 100; TT <- c(0.5, 1, 2, 4, 8)
B_WT <- 0.75; B_CR <- 0.40; OM <- 0.30; ADD <- 0.20
S1 <- 0.30      # sd of log(WT/70)
S2 <- 0.35      # sd of log(CRCL/100)
SRC <- list(A = list(rho = 0.20, n = 120), B = list(rho = 0.70, n = 120))

gh <- function(k) {
  i <- seq_len(k - 1L)
  J <- matrix(0, k, k); J[cbind(i, i+1L)] <- sqrt(i); J[cbind(i+1L, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE); o <- order(e$values)
  list(x = e$values[o], w = e$vectors[1L, o]^2)
}
QE <- gh(15L)   # eta
QC <- gh(9L)    # CRCL given WT
conc <- function(lcl) DOSE / exp(TV) * exp(outer(-exp(lcl) / exp(TV), TT))

## what each source PUBLISHED: weight only, IIV inflated by the omission
src_par <- function(rho) list(
  a  = B_WT + B_CR * rho * S2 / S1,
  om = sqrt(OM^2 + B_CR^2 * S2^2 * (1 - rho^2)))

## weight nodes inside stratum k of K (equal-probability sub-slices)
wt_nodes <- function(k, K, m = 15L) {
  u <- ((k - 1L) + (seq_len(m) - 0.5) / m) / K
  list(x = S1 * stats::qnorm(u), w = rep(1/m, m))
}
moments <- function(lcl_of, QW) {          # lcl_of(x1) -> vector over eta nodes
  m1 <- 0; M2 <- 0
  for (i in seq_along(QW$x)) {
    Y <- conc(lcl_of(QW$x[i]))
    for (j in seq_along(QE$x)) {
      w <- QW$w[i] * QE$w[j]; y <- Y[j, ]
      m1 <- m1 + w * y; M2 <- M2 + w * outer(y, y)
    }
  }
  V <- M2 - outer(m1, m1); diag(V) <- diag(V) + ADD^2
  list(E = m1, V = V)
}
## SOURCE block: its own model (weight only), over stratum k of its population
src_block <- function(rho, k, K) {
  sp <- src_par(rho)
  moments(function(x1) TCL + sp$a * x1 + sp$om * QE$x, wt_nodes(k, K))
}
## UNIFIED prediction for that block. joint = TRUE integrates CRCL over its
## CONDITIONAL law given the stratum's weight; FALSE uses the unconditional one.
uni_block <- function(p, rho, k, K, joint = TRUE) {
  tcl <- p[1]; b1 <- p[2]; b2 <- p[3]; om <- exp(p[4])
  cm  <- if (joint) rho * S2 / S1 else 0
  cs  <- if (joint) S2 * sqrt(1 - rho^2) else S2
  QW  <- wt_nodes(k, K)
  m1 <- 0; M2 <- 0
  for (i in seq_along(QW$x)) for (c in seq_along(QC$x)) {
    x1 <- QW$x[i]; x2 <- cm * x1 + cs * QC$x[c]
    Y  <- conc(tcl + b1 * x1 + b2 * x2 + om * QE$x)
    for (j in seq_along(QE$x)) {
      w <- QW$w[i] * QC$w[c] * QE$w[j]; y <- Y[j, ]
      m1 <- m1 + w * y; M2 <- M2 + w * outer(y, y)
    }
  }
  V <- M2 - outer(m1, m1); diag(V) <- diag(V) + ADD^2
  list(E = m1, V = V)
}
nll2 <- function(obs, pred, n) {
  ch <- tryCatch(chol(pred$V), error = function(e) NULL)
  if (is.null(ch)) return(1e12)
  iv <- chol2inv(ch); r <- obs$E - pred$E
  n * (2*sum(log(diag(ch))) + sum(iv * obs$V) + as.numeric(t(r) %*% iv %*% r))
}
project <- function(K, joint = TRUE) {
  dat <- lapply(names(SRC), function(s)
    lapply(seq_len(K), function(k) src_block(SRC[[s]]$rho, k, K)))
  names(dat) <- names(SRC)
  obj <- function(p) {
    tot <- 0
    for (s in names(SRC)) for (k in seq_len(K))
      tot <- tot + nll2(dat[[s]][[k]],
                        uni_block(p, SRC[[s]]$rho, k, K, joint),
                        SRC[[s]]$n / K)
    if (!is.finite(tot)) 1e12 else tot
  }
  o <- optim(c(TCL, 0.5, 0.1, log(0.4)), obj, method = "Nelder-Mead",
             control = list(maxit = 4000, reltol = 1e-12))
  o <- optim(o$par, obj, method = "Nelder-Mead",
             control = list(maxit = 4000, reltol = 1e-12))
  c(b_WT = o$par[2], b_CRCL = o$par[3], om = exp(o$par[4]))
}

cat("Two published models, NEITHER containing CRCL.\n\n")
cat(sprintf("%-10s %10s %10s %10s\n", "source", "rho", "pub. b_WT", "pub. IIV sd"))
for (s in names(SRC)) {
  sp <- src_par(SRC[[s]]$rho)
  cat(sprintf("%-10s %10.2f %10.4f %10.4f\n", s, SRC[[s]]$rho, sp$a, sp$om))
}
cat(sprintf("\ntruth:      b_WT %.4f   b_CRCL %.4f   IIV sd %.4f\n\n",
            B_WT, B_CR, OM))
cat(sprintf("%-42s %9s %9s %9s\n", "projection", "b_WT", "b_CRCL", "IIV sd"))
for (K in c(2L, 4L, 8L)) {
  r <- project(K, joint = TRUE)
  cat(sprintf("%-42s %9.4f %9.4f %9.4f\n",
              sprintf("joint covariate law,  K = %d", K), r[1], r[2], r[3]))
}
r <- project(4L, joint = FALSE)
cat(sprintf("%-42s %9.4f %9.4f %9.4f   <- control\n",
            "CRCL drawn UNCONDITIONALLY, K = 4", r[1], r[2], r[3]))
