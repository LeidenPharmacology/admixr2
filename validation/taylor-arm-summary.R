## Summarise the 2x2 study with the Taylor arm added.
##   AGG   stratified on ALL covariates (the mismatched-grid construction)
##   AGGm  per-covariate: stratify what the source fitted, marginalise the rest, by QUADRATURE
##   AGGt  same split, marginalised directions by the corrected 2nd-order TAYLOR expansion
##   AGGt0 the Taylor arm WITHOUT the rank-one s^2 dE dE^T term (the guard)
OUT <- "C:/Users/hidde/AppData/Local/Temp/claude/C--package-admixr2/3ff305c7-64a5-4cb0-ba61-1436e2e9b16e/scratchpad"
d <- readRDS(file.path(OUT, "aggregate-marginal-taylor.rds"))
R <- d$RES; TP <- d$TP; say <- function(...) cat(..., "\n")
cells <- vapply(R, `[[`, "", "cell"); CN <- c("baseline","struct","omit","both")
rmse <- function(a,b) sqrt(mean((a-b)^2)); mcse <- function(v) sd(v)/sqrt(length(v))
has <- function(nm) all(vapply(R, function(r) is.numeric(r[[nm]]) && length(r[[nm]])==length(TP), TRUE))
ARMS <- c("AGG","AGGm","AGGt")[c(TRUE,TRUE,has("AGGt"))]
if (has("AGGt0")) ARMS <- c(ARMS, "AGGt0")

say("=========== deviation from the INDIVIDUAL-DATA fit (primary) ===========")
say(sprintf("%-10s %5s %s", "cell","reps", paste(sprintf("%12s",ARMS), collapse="")))
for (cn in CN) { ix <- which(cells==cn); if(!length(ix)) next
  v <- sapply(ARMS, function(a) mean(vapply(ix, function(i) rmse(R[[i]][[a]],R[[i]]$GOLD),0)))
  say(sprintf("%-10s %5d %s", cn, length(ix), paste(sprintf("%12.4f",v), collapse=""))) }

say("\n=========== deviation from the simulation truth ===========")
say(sprintf("%-10s %5s %12s %s", "cell","reps","IPD (gold)", paste(sprintf("%12s",ARMS),collapse="")))
for (cn in CN) { ix <- which(cells==cn); if(!length(ix)) next
  g <- mean(vapply(ix, function(i) rmse(R[[i]]$GOLD,TP),0))
  v <- sapply(ARMS, function(a) mean(vapply(ix, function(i) rmse(R[[i]][[a]],TP),0)))
  say(sprintf("%-10s %5d %12.4f %s", cn, length(ix), g, paste(sprintf("%12.4f",v),collapse=""))) }

say(sprintf("\n=========== bCL2 (truth %.3f) ===========", TP[["bCL2"]]))
say(sprintf("%-10s %14s %s","cell","IPD (gold)", paste(sprintf("%14s",ARMS),collapse="")))
for (cn in CN) { ix <- which(cells==cn); if(!length(ix)) next
  f <- function(a){ x<-vapply(ix,function(i) R[[i]][[a]][["bCL2"]],0)
    sprintf("%6.3f%+7.3f", mean(x), mean(x)-TP[["bCL2"]]) }
  say(sprintf("%-10s %14s %s", cn, f("GOLD"), paste(sprintf("%14s",sapply(ARMS,f)),collapse=""))) }
say("(mean, and bias against truth)")

if (has("AGGt")) {
say("\n=========== PAIRED: Taylor vs quadrature on the marginalised directions ===========")
say(sprintf("%-10s %5s %14s %10s %10s %10s","cell","reps","mean d(bCL2)","se","AGGt wins","p"))
for (cn in CN) { ix <- which(cells==cn); if(!length(ix)) next
  dd <- vapply(ix, function(i) R[[i]]$AGGt[["bCL2"]]-R[[i]]$AGGm[["bCL2"]], 0)
  ae <- vapply(ix, function(i) abs(R[[i]]$AGGt[["bCL2"]]-TP[["bCL2"]]) -
                               abs(R[[i]]$AGGm[["bCL2"]]-TP[["bCL2"]]), 0)
  p <- if (sd(dd)>0) tryCatch(t.test(dd)$p.value, error=function(e) NA_real_) else NA_real_
  say(sprintf("%-10s %5d %14.5f %10.5f %10s %10.3g", cn, length(ix), mean(dd), mcse(dd),
              sprintf("%d/%d", sum(ae<0), length(ix)), p)) }
say("('AGGt wins' = Taylor closer to truth than quadrature on that replicate)")
say("\n=========== paired RMSE-to-gold difference (AGGt - AGGm) ===========")
for (cn in CN) { ix <- which(cells==cn); if(!length(ix)) next
  dd <- vapply(ix, function(i) rmse(R[[i]]$AGGt,R[[i]]$GOLD)-rmse(R[[i]]$AGGm,R[[i]]$GOLD),0)
  say(sprintf("  %-9s %+10.5f  (se %.5f)  AGGt better in %d/%d", cn, mean(dd), mcse(dd),
              sum(dd<0), length(ix))) } }

say("\n=========== moment accuracy vs a 21-node reference ===========")
say(sprintf("%-16s %12s %12s %12s %12s %10s","block","relE AGGm","relV AGGm","relE AGGt","relV AGGt","rank1/V"))
for (nm in names(d$MOMCHK)) { m <- d$MOMCHK[[nm]]
  say(sprintf("%-16s %12.2e %12.2e %12.2e %12.2e %10s", nm, m[1],m[2],m[3],m[4],
              if (is.na(m[5])) "-" else sprintf("%.3f", m[5]))) }

say("\n=========== cost per replicate ===========")
for (nm in c("tim","nev","rows")) {
  if (!all(vapply(R, function(r) is.numeric(r[[nm]]), TRUE))) next
  say(sprintf(" %s:", nm))
  for (cn in CN) { ix <- which(cells==cn); if(!length(ix)) next
    x <- R[[ix[1]]][[nm]]
    say(sprintf("   %-9s %s", cn, paste(sprintf("%s=%.4g", names(x), colMeans(
      do.call(rbind, lapply(ix, function(i) R[[i]][[nm]])))), collapse="  "))) } }
say(sprintf("\nhfrac = %.2f", d$HFRAC))
