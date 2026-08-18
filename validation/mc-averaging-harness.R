## Harness for the model-averaging comparison.
##
## Reproduces adgh on the `rows` path EXACTLY: a product Gauss-Hermite grid over
## the covariate (.admCovGrid, cov_nodes) tensored with the eta grid, evaluated
## per row, then the law-of-total-variance composition. Node-for-node the same
## finite weighted sum admixr2 computes, so the two must agree to near machine
## precision -- that is the validation gate in mc-averaging-validate.R.
##
## Model: 1-cmt IV bolus.   x = log(WT/70) ~ N(mu_x, sd_x)
##   log CL = lcl + bCL*x + eta_CL      log V = lv + bV*x + eta_V
##   C(t) = D/V * exp(-CL/V * t)        y = f*(1+eps),  eps ~ N(0, sig^2)
PAR_NM <- c("lcl","lv","bCL","bV","oCL","oV","sig")

gh <- function(k) {                      # Golub-Welsch, probabilists' scaling
  i <- seq_len(k - 1L)
  J <- matrix(0, k, k); J[cbind(i, i+1L)] <- sqrt(i); J[cbind(i+1L, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE); o <- order(e$values)
  list(x = e$values[o], w = e$vectors[1L, o]^2)
}
GH_E <- gh(5L)      # eta nodes  (adghControl n_nodes)
GH_X <- gh(7L)      # covariate nodes (pinfo$cov_nodes default)

## Covariate nodes for stratum k of K, matching .admCovNodesFor on a `quantile`
## spec: GH nodes -> pnorm -> the stratum's own quantile function.
xnodes <- function(mu_x, sd_x, k = 1L, K = 1L)
  mu_x + sd_x * stats::qnorm(((k - 1L) + stats::pnorm(GH_X$x)) / K)

## The product grid: every covariate node x every eta node, weights multiplied.
grid_of <- function(mu_x, sd_x, k = 1L, K = 1L) {
  xs <- xnodes(mu_x, sd_x, k, K)
  g  <- expand.grid(ix = seq_along(xs), ie1 = seq_along(GH_E$x),
                    ie2 = seq_along(GH_E$x), KEEP.OUT.ATTRS = FALSE)
  w  <- GH_X$w[g$ix] * GH_E$w[g$ie1] * GH_E$w[g$ie2]
  list(x = xs[g$ix], z1 = GH_E$x[g$ie1], z2 = GH_E$x[g$ie2], w = w / sum(w))
}

## Aggregate moments for one block: E and V a study with this population and
## design would report under parameters `p`.
moments <- function(p, gr, times, dose) {
  cl <- exp(p[["lcl"]] + p[["bCL"]] * gr$x + p[["oCL"]] * gr$z1)
  v  <- exp(p[["lv"]]  + p[["bV"]]  * gr$x + p[["oV"]]  * gr$z2)
  f  <- outer(dose / v, times, function(a, b) a) * exp(-outer(cl / v, times))
  mu <- as.numeric(crossprod(gr$w, f))
  fc <- sweep(f, 2L, mu)
  V  <- crossprod(fc, fc * gr$w)                      # Cov_{x,eta}(f)
  ## law of total variance, proportional error: E[Var(y|.)] = sig^2 E[f^2]
  Ef2 <- as.numeric(crossprod(gr$w, f^2))
  diag(V) <- diag(V) + p[["sig"]]^2 * Ef2
  list(E = mu, V = V)
}

## admixr2's NLL: n * (log|V| + tr(V^-1 Vobs) + r' V^-1 r), summed over blocks.
nll_blocks <- function(p, blocks) {
  tot <- 0
  for (b in blocks) {
    pr <- moments(p, b$gr, b$times, b$dose)
    ch <- tryCatch(chol(pr$V), error = function(e) NULL)
    if (is.null(ch)) return(1e12)
    iv <- chol2inv(ch); r <- b$E - pr$E
    tot <- tot + b$n * (2*sum(log(diag(ch))) + sum(iv * b$V) +
                        as.numeric(t(r) %*% iv %*% r))
  }
  if (!is.finite(tot)) 1e12 else tot
}

## Blocks a source produces: K covariate strata, n/K each, at its own design.
make_blocks <- function(p, src, K) lapply(seq_len(K), function(k) {
  gr <- grid_of(src$mu_x, src$sd_x, k, K)
  m  <- moments(p, gr, src$times, src$dose)
  list(gr = gr, times = src$times, dose = src$dose, n = src$n / K,
       E = m$E, V = m$V)
})

## Fit a MODEL FORM (a subset of free parameters; the rest held at `fix`) to
## blocks. This is both the aggregate estimator and the inner problem of the
## induced map -- deliberately one function, so they cannot drift apart.
fit_form <- function(blocks, free, start, fix = NULL) {
  o <- stats::optim(start[free], function(v) {
    p <- start; p[free] <- v; if (length(fix)) p[names(fix)] <- fix
    p[["oCL"]] <- abs(p[["oCL"]]); p[["oV"]] <- abs(p[["oV"]])
    p[["sig"]] <- abs(p[["sig"]]); nll_blocks(p, blocks)
  }, method = "Nelder-Mead", control = list(maxit = 3000, reltol = 1e-11))
  o <- stats::optim(o$par, function(v) {
    p <- start; p[free] <- v; if (length(fix)) p[names(fix)] <- fix
    p[["oCL"]] <- abs(p[["oCL"]]); p[["oV"]] <- abs(p[["oV"]])
    p[["sig"]] <- abs(p[["sig"]]); nll_blocks(p, blocks)
  }, method = "Nelder-Mead", control = list(maxit = 3000, reltol = 1e-12))
  out <- start; out[free] <- o$par
  if (length(fix)) out[names(fix)] <- fix
  out[c("oCL","oV","sig")] <- abs(out[c("oCL","oV","sig")])
  out
}
