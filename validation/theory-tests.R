## THREE FALSIFIABLE PREDICTIONS OF THE THEORY, tested against exact quadrature.
## No sampling noise anywhere: data are exact moments at the true parameters, so a
## displaced argmin is a property of the OBJECTIVE, not a finite-sample effect.
##
## T1  DPI scaling.  If the marginal objective fails because marginalisation is
##     not a sufficient statistic, its bias must VANISH as Var(x) -> 0 and GROW
##     with Var(x). The conditional objective must be exact at every Var(x).
##     A bias that does not scale with covariate spread would refute the account.
##
## T2  Attenuation identity.  E_a[ r(a)' Vi r(a) ] = rbar' Vi rbar + tr(Vi Cov_a(mu)).
##     The second term carries no data -- it is a pure penalty on covariate-induced
##     mean spread, and is why the marginal/pooled route drives b -> 0.
##
## T3  Omitted-covariate binding function.  A covariate-FREE model fitted to data
##     from a covariate-CARRYING one must return  omega_s^2 = omega^2 + b^2 Var(x).
##     This is g_s in closed form, and it is the channel by which a publication
##     with no covariate term still carries covariate information.
say <- function(...) { cat(...,"\n"); utils::flush.console() }
gh <- function(k){ i <- seq_len(k-1L); J <- matrix(0,k,k)
  J[cbind(i,i+1L)] <- sqrt(i); J[cbind(i+1L,i)] <- sqrt(i)
  e <- eigen(J,symmetric=TRUE); o <- order(e$values); list(x=e$values[o],w=e$vectors[1L,o]^2) }
GE <- gh(15L); GX <- gh(15L)
TT <- c(0.5,1,2,4,8); DOSE <- 100; LCL <- log(1); LV <- log(10)
B <- 0.75; OM <- 0.30; SIG <- 0.10

conc <- function(cl) DOSE/exp(LV) * exp(-outer(cl/exp(LV), TT))
## moments at a FIXED covariate value (eta integrated only) -- the node "data"
mom_cond <- function(a, lcl, b, om, sig) {
  f <- conc(exp(lcl + b*a + om*GE$x))
  mu <- as.numeric(crossprod(GE$w, f)); fc <- sweep(f,2L,mu)
  V <- crossprod(fc, fc*GE$w)
  diag(V) <- diag(V) + sig^2*as.numeric(crossprod(GE$w, f^2)); list(E=mu, V=V)
}
## moments marginal over BOTH x and eta -- the pooled summary
mom_marg <- function(sx, lcl, b, om, sig) {
  m1 <- 0; M2 <- 0; Ef2 <- 0
  for (i in seq_along(GX$x)) for (j in seq_along(GE$x)) {
    w <- GX$w[i]*GE$w[j]
    y <- conc(exp(lcl + b*sx*GX$x[i] + om*GE$x[j]))[1,]
    m1 <- m1 + w*y; M2 <- M2 + w*outer(y,y); Ef2 <- Ef2 + w*y^2
  }
  V <- M2 - outer(m1,m1); diag(V) <- diag(V) + sig^2*Ef2; list(E=m1, V=V)
}
nll <- function(o, p, n) { ch <- tryCatch(chol(p$V), error=function(e) NULL)
  if (is.null(ch)) return(1e12); iv <- chol2inv(ch); r <- o$E-p$E
  n*(2*sum(log(diag(ch))) + sum(iv*o$V) + as.numeric(t(r)%*%iv%*%r)) }

fit_marg <- function(obs, sx, n=100) {                # marginal objective
  o <- optim(c(LCL,B,OM,SIG), function(v)
    nll(obs, mom_marg(sx, v[1], v[2], abs(v[3]), abs(v[4])), n),
    method="Nelder-Mead", control=list(maxit=3000, reltol=1e-12))
  c(lcl=o$par[1], b=o$par[2], om=abs(o$par[3]), sig=abs(o$par[4]))
}
fit_cond <- function(blocks, sx, n=100) {             # conditional (node) objective
  nd <- sx*GX$x; wt <- GX$w/sum(GX$w)
  o <- optim(c(LCL,B,OM,SIG), function(v) sum(vapply(seq_along(nd), function(k)
      nll(blocks[[k]], mom_cond(nd[k], v[1], v[2], abs(v[3]), abs(v[4])), wt[k]*n), 0)),
    method="Nelder-Mead", control=list(maxit=3000, reltol=1e-12))
  c(lcl=o$par[1], b=o$par[2], om=abs(o$par[3]), sig=abs(o$par[4]))
}
fit_nocov <- function(obs, sx, n=100) {               # covariate-FREE model, pooled data
  o <- optim(c(LCL,OM,SIG), function(v)
    nll(obs, mom_marg(0, v[1], 0, abs(v[2]), abs(v[3])), n),
    method="Nelder-Mead", control=list(maxit=3000, reltol=1e-12))
  c(lcl=o$par[1], om=abs(o$par[2]), sig=abs(o$par[3]))
}

## ---------------- T1 + T3 --------------------------------------------------
SX <- c(0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40)
say(sprintf("truth b = %.4f, omega = %.4f\n", B, OM))
say(sprintf("%6s | %9s %9s | %9s %9s | %9s %9s", "sd(x)",
            "b marginal", "rel err", "b nodes", "rel err", "om_free", "predicted"))
res <- NULL
for (sx in SX) {
  obs_m <- mom_marg(sx, LCL, B, OM, SIG)
  nd <- sx*GX$x
  obs_c <- lapply(nd, function(a) mom_cond(a, LCL, B, OM, SIG))
  fm <- fit_marg(obs_m, sx); fc <- fit_cond(obs_c, sx); fn <- fit_nocov(obs_m, sx)
  pred_om <- sqrt(OM^2 + B^2*sx^2)                    # T3 closed form
  say(sprintf("%6.2f | %9.4f %8.1f%% | %9.4f %8.1f%% | %9.4f %9.4f",
      sx, fm[["b"]], 100*(fm[["b"]]/B-1), fc[["b"]], 100*(fc[["b"]]/B-1),
      fn[["om"]], pred_om))
  res <- rbind(res, c(sx=sx, b_m=fm[["b"]], b_c=fc[["b"]], om_f=fn[["om"]], om_p=pred_om))
}
res <- as.data.frame(res)
say(sprintf("\nT1  marginal |rel err| vs sd(x): correlation = %.3f  (predicted: strongly +)",
            cor(res$sx, abs(res$b_m/B - 1))))
say(sprintf("T1  marginal err at smallest sd(x) = %.2f%%   at largest = %.2f%%",
            100*(res$b_m[1]/B-1), 100*(res$b_m[nrow(res)]/B-1)))
say(sprintf("T1  node-objective max |rel err| over all sd(x) = %.3f%%",
            100*max(abs(res$b_c/B - 1))))
say(sprintf("T3  max |om_free - sqrt(om^2+b^2 var(x))| = %.5f", max(abs(res$om_f-res$om_p))))

## ---------------- T2: the attenuation identity ------------------------------
sx <- 0.30; obs <- mom_marg(sx, LCL, B, OM, SIG)
nd <- sx*GX$x; wt <- GX$w/sum(GX$w)
MU <- vapply(nd, function(a) mom_cond(a, LCL, B, OM, SIG)$E, numeric(length(TT)))
Vi <- chol2inv(chol(obs$V)); rj <- obs$E - MU
pool  <- sum(wt * colSums(rj * (Vi %*% rj)))
mubar <- as.numeric(MU %*% wt); rb <- obs$E - mubar
Cova  <- (sweep(MU,1L,mubar) * rep(wt, each=length(TT))) %*% t(sweep(MU,1L,mubar))
say(sprintf("\nT2  E_a[r' Vi r]            = %12.8f", pool))
say(sprintf("T2  rbar' Vi rbar           = %12.8f   <- the only data-bearing term",
            as.numeric(t(rb)%*%Vi%*%rb)))
say(sprintf("T2  tr(Vi Cov_a(mu))        = %12.8f   <- pure covariate-spread penalty",
            sum(Vi*Cova)))
say(sprintf("T2  identity holds: %s",
    isTRUE(all.equal(pool, as.numeric(t(rb)%*%Vi%*%rb)+sum(Vi*Cova), tolerance=1e-9))))
saveRDS(res, "C:/Users/hidde/AppData/Local/Temp/claude/C--package-admixr2/3ff305c7-64a5-4cb0-ba61-1436e2e9b16e/scratchpad/theory-tests.rds")
