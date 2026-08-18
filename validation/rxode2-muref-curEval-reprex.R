## ===========================================================================
## REPREX for an upstream rxode2 issue -- self-contained, rxode2 only.
##
## muRefCurEval loses the transform of a mu-referenced parameter when the
## mu-referenced expression uses a variable computed earlier in the model
## ("mu referencing 3.0" style, https://blog.nlmixr2.org/blog/2024-12-11-mu3/).
##
## The two models below are MATHEMATICALLY IDENTICAL -- bit-identical
## predictions -- and differ only in whether log(WT/70) is written inline or
## assigned to a variable first. rxode2 reports curEval = "exp" for the first
## and "" for the second, while muRefDataFrame correctly identifies tcl~eta.cl
## in both. So the parameter is known to be mu-referenced; only its transform
## is lost.
##
## Consequence for any consumer that back-transforms structural parameters via
## muRefCurEval: a log-scale estimate is reported untransformed (1.386 where
## the truth is 4).
##
## Paste into a GitHub issue as-is. No admixr2 dependency.
## ===========================================================================
suppressMessages(library(rxode2))
cat(sprintf("rxode2 %s\n\n", as.character(utils::packageVersion("rxode2"))))

inline <- function() {
  ini({tcl <- log(4); tv <- log(30); tcov <- 0.75; eta.cl ~ 0.09; a <- 0.1})
  model({cl <- exp(tcl + tcov * log(WT / 70) + eta.cl); v <- exp(tv)
         d/dt(centr) <- -cl / v * centr; cp <- centr / v; cp ~ add(a)})
}
precomputed <- function() {
  ini({tcl <- log(4); tv <- log(30); tcov <- 0.75; eta.cl ~ 0.09; a <- 0.1})
  model({wt70 <- log(WT / 70)
         cl <- exp(tcl + tcov * wt70 + eta.cl); v <- exp(tv)
         d/dt(centr) <- -cl / v * centr; cp <- centr / v; cp ~ add(a)})
}
ui1 <- suppressMessages(rxode2::rxode2(inline))
ui2 <- suppressMessages(rxode2::rxode2(precomputed))

## ---- 1. the two models really are the same ---------------------------------
ev   <- rxode2::et(rxode2::et(amt = 500), c(0.5, 1, 2, 4, 8, 12))
pars <- c(tcl = log(4), tv = log(30), tcov = 0.75, eta.cl = 0.2, WT = 92)
sim <- function(ui) {
  m <- rxode2::rxode2(ui$simulationModel)
  p <- pars[intersect(names(pars), m$params)]
  for (q in setdiff(m$params, names(p))) p[q] <- 0
  s <- rxode2::rxSolve(m, params = p, events = ev, returnType = "data.frame",
                       addDosing = FALSE, atol = 1e-12, rtol = 1e-12)
  s$cp[!is.na(s$cp)]
}
p1 <- sim(ui1); p2 <- sim(ui2)
cat(sprintf("predictions identical?              max rel diff %.3e\n",
            max(abs(p1 - p2) / abs(p1))))

## ---- 2. but the mu-reference metadata differs -------------------------------
ce <- function(ui, p) { d <- ui$muRefCurEval; d$curEval[d$parameter == p] }
cat(sprintf("muRefCurEval['tcl']    inline = %-6s  precomputed = %s\n",
            sQuote(ce(ui1, "tcl")), sQuote(ce(ui2, "tcl"))))
cat(sprintf("muRefCurEval['eta.cl'] inline = %-6s  precomputed = %s\n",
            sQuote(ce(ui1, "eta.cl")), sQuote(ce(ui2, "eta.cl"))))
cat(sprintf("muRefDataFrame         inline = %-6s  precomputed = %s\n",
            paste(ui1$muRefDataFrame$theta, ui1$muRefDataFrame$eta, sep = "~"),
            paste(ui2$muRefDataFrame$theta, ui2$muRefDataFrame$eta, sep = "~")))

## ---- 3. what that costs a consumer -----------------------------------------
bt <- function(cv, p) switch(cv, exp = , log = exp(p), p)
cat(sprintf("\nback-transform of tcl = log(4) = %.6f\n", log(4)))
cat(sprintf("  inline      curEval %-6s -> %.6f   (correct, exp(log 4) = 4)\n",
            sQuote(ce(ui1, "tcl")), bt(ce(ui1, "tcl"), log(4))))
cat(sprintf("  precomputed curEval %-6s -> %.6f   (reported on the log scale)\n",
            sQuote(ce(ui2, "tcl")), bt(ce(ui2, "tcl"), log(4))))
