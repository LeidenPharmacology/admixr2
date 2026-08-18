## Summarise the overnight 2x2 study, with the INDIVIDUAL-DATA fit as the reference.
##
## Why IPD is the gold standard and not the simulation truth. IPD-POOL fits the
## CORRECT structure -- 2-compartment, both covariates -- to the pooled individual
## records. It is therefore the only arm that can recover the data-generating model,
## and it is the ceiling for any summary-level method: no construction working from
## published summaries can do better than a fit that saw every observation.
##
## Reporting against truth alone is misleading, because in any one replicate the
## individual data do not point exactly at truth either. A summary method that lands
## nearer truth than IPD-POOL did has not beaten the data; it has been lucky. So the
## primary metric is deviation from IPD-POOL, with deviation from truth reported
## alongside for reference.
OUT <- "C:/Users/hidde/AppData/Local/Temp/claude/C--package-admixr2/3ff305c7-64a5-4cb0-ba61-1436e2e9b16e/scratchpad"
d <- readRDS(file.path(OUT,"overnight.rds")); RES <- d$RES; TP <- d$TP
UN <- names(TP); say <- function(...) cat(..., "\n")
cells <- vapply(RES, `[[`, "", "cell")
say(sprintf("%d replicates: %s\n", length(RES),
  paste(sprintf("%s %d", names(table(cells)), as.integer(table(cells))), collapse="  ")))
ARMS <- c(AGG="AGG", PLLs="PLLs", PLLr="PLLr")
rmse <- function(a,b) sqrt(mean((a-b)^2))
mcse <- function(v) stats::sd(v)/sqrt(length(v))

say("=========== deviation from the INDIVIDUAL-DATA fit (primary) ===========")
say(sprintf("%-10s %6s %18s %18s %18s", "cell","reps","AGG","PLL (computed W)","PLL (reported SE)"))
for (cn in unique(cells)) { ix <- which(cells==cn); if (!length(ix)) next
  v <- lapply(ARMS, function(a) vapply(ix, function(i) rmse(RES[[i]][[a]],RES[[i]]$GOLD), 0))
  say(sprintf("%-10s %6d %10.4f (%.4f) %10.4f (%.4f) %10.4f (%.4f)", cn, length(ix),
      mean(v$AGG),mcse(v$AGG), mean(v$PLLs),mcse(v$PLLs), mean(v$PLLr),mcse(v$PLLr))) }

say("\n=========== deviation from the simulation truth (reference) ===========")
say(sprintf("%-10s %6s %12s %12s %12s %12s", "cell","reps","IPD (gold)","AGG","PLLs","PLLr"))
for (cn in unique(cells)) { ix <- which(cells==cn); if (!length(ix)) next
  g <- vapply(ix, function(i) rmse(RES[[i]]$GOLD,TP), 0)
  v <- lapply(ARMS, function(a) vapply(ix, function(i) rmse(RES[[i]][[a]],TP), 0))
  say(sprintf("%-10s %6d %12.4f %12.4f %12.4f %12.4f", cn, length(ix),
      mean(g), mean(v$AGG), mean(v$PLLs), mean(v$PLLr))) }

say("\n=========== bCL2: the covariate ABSENT from source A in `omit`/`both` ===========")
say(sprintf("truth %.3f", TP[["bCL2"]]))
say(sprintf("%-10s %6s %14s %14s %14s %14s", "cell","reps","IPD (gold)","AGG","PLLs","PLLr"))
for (cn in unique(cells)) { ix <- which(cells==cn); if (!length(ix)) next
  f <- function(a) { v <- vapply(ix, function(i) RES[[i]][[a]][["bCL2"]], 0)
    sprintf("%6.3f%+7.3f", mean(v), mean(v)-TP[["bCL2"]]) }
  say(sprintf("%-10s %6d %14s %14s %14s %14s", cn, length(ix),
      f("GOLD"), f("AGG"), f("PLLs"), f("PLLr"))) }
say("\n(mean, and bias against truth)")

say("\n=========== does any summary method beat the individual data? ===========")
for (cn in unique(cells)) { ix <- which(cells==cn); if (!length(ix)) next
  g <- vapply(ix, function(i) rmse(RES[[i]]$GOLD,TP), 0)
  for (a in names(ARMS)) { v <- vapply(ix, function(i) rmse(RES[[i]][[a]],TP), 0)
    say(sprintf("  %-9s %-5s beats IPD in %2d/%2d replicates", cn, a, sum(v<g), length(ix))) } }
say("\n(a summary method beating IPD is luck, not skill -- IPD saw every observation)")
