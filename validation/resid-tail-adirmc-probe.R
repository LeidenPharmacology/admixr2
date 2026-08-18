# adirmc half of the moment-tail bit-identity check ----------------------------
#
# .adirmcInnerNLL()/.adirmcInnerGrad() need a PROPOSAL object, which the general
# harness (validation/resid-moment-tail-bitidentity.R) does not build, so
# adirmc's chain conversion would otherwise be unmeasured. Same protocol:
#   Rscript validation/resid-tail-adirmc-probe.R before   (unconverted tree)
#   Rscript validation/resid-tail-adirmc-probe.R after    (converted tree)
#   Rscript validation/resid-tail-adirmc-probe.R compare

suppressMessages(devtools::load_all(".", quiet = TRUE))
mode <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(mode)) mode <- "before"
outf <- file.path("validation", sprintf("resid-tail-irmc-%s.rds", mode))

if (identical(mode, "compare")) {
  b <- readRDS(file.path("validation", "resid-tail-irmc-before.rds"))
  a <- readRDS(file.path("validation", "resid-tail-irmc-after.rds"))
  nm <- union(names(b), names(a)); tot <- 0L; eq <- 0L; mx <- 0
  for (k in nm) {
    x <- as.numeric(b[[k]]); y <- as.numeric(a[[k]])
    if (length(x) != length(y)) { cat(k, ": LENGTH MISMATCH\n"); next }
    bx <- vapply(x, function(v) paste(writeBin(v, raw()), collapse = ""), "")
    by <- vapply(y, function(v) paste(writeBin(v, raw()), collapse = ""), "")
    e <- bx == by
    tot <- tot + length(x); eq <- eq + sum(e)
    if (any(!e)) mx <- max(mx, max(abs(x[!e] - y[!e]) / pmax(abs(x[!e]), 1e-300)))
    cat(sprintf("%-28s n=%3d bit_eq=%3d\n", k, length(x), sum(e)))
  }
  cat(sprintf("\nadirmc: %d/%d bit-identical, max rel %.3e\n", eq, tot, mx))
  quit(save = "no")
}

.mk <- function(ini, err) suppressMessages(rxode2::rxode2(eval(parse(text = sprintf(
  'function() { ini({ tcl <- log(5); tv <- log(20); %s
                      eta.cl ~ 0.09; eta.v ~ 0.04 })
     model({ cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
             d/dt(central) <- -(cl/v)*central; cp <- central/v
             %s }) }', ini, err)))))

cases <- list(
  add   = c("a <- 0.5",            "cp ~ add(a)"),
  prop  = c("b <- 0.15",           "cp ~ prop(b)"),
  comb2 = c("a <- 0.5; b <- 0.15", "cp ~ add(a) + prop(b)"),
  comb1 = c("a <- 0.5; b <- 0.15", "cp ~ add(a) + prop(b) + combined1()"),
  lnorm = c("a <- 0.3",            "cp ~ lnorm(a)"))

tt <- c(0.5, 1, 2, 4); ev <- rxode2::et(amt = 100)
E0 <- 100 / 20 * exp(-(5 / 20) * tt)
Vcov <- outer(0.3 * E0, 0.3 * E0) * (0.4 ^ abs(outer(tt, tt, "-")))
Vvar <- diag((0.3 * E0) ^ 2)

res <- list()
for (nm in names(cases)) {
  cs <- cases[[nm]]
  ui <- .mk(cs[1], cs[2])
  pinfo <- suppressWarnings(.admParseIniDf(ui$iniDf, ui))
  rxMod <- .admLoadModel(ui)
  ovar  <- .admOutputVar(ui)
  for (meth in c("var", "cov")) {
    s <- .admNormaliseStudy(list(E = E0, V = if (meth == "cov") Vcov else Vvar,
                                 n = 200L, times = tt, ev = ev), "s")
    s$ev_full <- rxode2::et(s$ev, s$times)
    z  <- .admMakeZ(500L, pinfo, 1L, "sobol")
    pm <- .admMakeParamsList(500L, pinfo, 1L)
    p0 <- .admBuildOptVec(pinfo)$p0
    pars0 <- .admUnpack(p0, pinfo)
    prop <- .adirmcProposal(rxMod, pars0$struct, pinfo$sigma_names, pars0$omega,
                            omega_expansion = 2, s, z[[1L]], ovar, pm[[1L]],
                            cores = 1L, pinfo$eta_col_names,
                            has_kappa = pinfo$has_kappa,
                            struct_transforms = pinfo$struct_transforms,
                            struct_eta_idx = pinfo$struct_eta_idx,
                            use_grad = TRUE)
    if (is.null(prop)) { message("no proposal: ", nm, "/", meth); next }
    pp <- p0 * 1.03 + 0.01
    res[[sprintf("%s/%s/irmc.nll", nm, meth)]] <-
      .adirmcNLL(pp, pinfo, list(s), list(prop))
    res[[sprintf("%s/%s/irmc.grad", nm, meth)]] <-
      .adirmcInnerGrad(pp, pinfo, list(s), list(prop))
  }
  message("done: ", nm)
}
saveRDS(res, outf)
cat(sprintf("wrote %s (%d keys, %d doubles)\n", outf, length(res),
            sum(vapply(res, length, integer(1)))))
