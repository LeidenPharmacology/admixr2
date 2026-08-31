# NO SOURCE ESTIMATION AT ALL. Each "published model" is written down at the
# values its analyst would converge to in the limit, with the omega that
# exactly absorbs the CRCL variation their model omits. If bcrcl still comes
# back at ~0.55, the attenuation lives in the POOLED fit, not in noise
# propagating from the sources.
Sys.setenv(NOT_CRAN = "true"); options(warn = 1)
suppressMessages({library(nlmixr2); library(rxode2)
                  pkgload::load_all(".", quiet = TRUE)})
CL<-5; BCRCL<-0.6; BSEX<-0.18; OM<-0.05
SD_W <- 0.198; SD_C <- 0.217; RHO <- as.numeric(Sys.getenv("RHO", "0.45"))
v_C  <- SD_C^2
# what a source's omega must absorb: its model omits CRCL, so the CRCL-driven
# part of log-CL variance lands in the random effect
OM_SRC <- BCRCL^2 * v_C + OM
cat(sprintf("RHO = %.2f
", RHO)); cat(sprintf("omega a source must report: %.5f  (truth %.3f + CRCL %.5f)\n",
            OM_SRC, OM, BCRCL^2 * v_C))
Rc <- matrix(c(1, RHO, RHO, 1), 2, 2)
CRCL_MED <- c(38, 62, 95); WT_MED <- c(74, 76, 78)
OBS <- c(0.5, 1, 2, 4, 8, 12, 24)

cohort <- function(n, wt_med, crcl_med) {
  z <- matrix(rnorm(2*n), n, 2) %*% chol(Rc)
  data.frame(WT   = qlnorm(pnorm(z[,1]), log(wt_med),   SD_W),
             CRCL = qlnorm(pnorm(z[,2]), log(crcl_med), SD_C),
             SEX  = rbinom(n, 1L, .55))
}
# the published model: correct in every respect EXCEPT that it has no renal
# term, exactly as the three analysts' models did
src_model <- function(tcl_j) {
  f <- function() {
    ini({ tcl <- 1; tv <- log(50); bsex <- 0.18; add.err <- 0.10
          eta.cl ~ 0.067 })
    model({ cl <- exp(tcl + eta.cl)*(WT/70)^0.75*exp(bsex*SEX)
            v  <- exp(tv)*(WT/70); cp <- linCmt(); cp ~ add(add.err) })
  }
  u <- suppressMessages(rxode2::rxode2(f))
  d <- u$iniDf
  d$est[d$name == "tcl"]    <- tcl_j
  d$est[d$name == "eta.cl"] <- OM_SRC
  u$iniDf <- d
  u
}
pooled <- function() {
  ini({ tcl <- log(4.5); tv <- log(55); bcrcl <- 0.3; bsex <- 0.05
        add.err <- 0.1; eta.cl ~ 0.04 })
  model({ cl <- exp(tcl + eta.cl)*(WT/70)^0.75*(CRCL/90)^bcrcl*exp(bsex*SEX)
          v <- exp(tv)*(WT/70); cp <- linCmt(); cp ~ add(add.err) })
}
nm <- c("tcl","tv","bsex","add.err","eta.cl")
Csrc <- diag(c(.005,.005,.005,.0005,.001)^2); dimnames(Csrc) <- list(nm,nm)

res <- list()
for (.rep in 1:8) {
  set.seed(1300 + .rep)
  ok <- try({
    sts <- lapply(1:3, function(j) {
      coh <- cohort(150L, WT_MED[j], CRCL_MED[j])
      # the analyst's typical CL for THIS cohort, in the limit: the geometric
      # mean, which is exactly linear in log(median CRCL) with slope 0.6
      tcl_j <- log(CL) + BCRCL * log(CRCL_MED[j]/90)
      admStudy(model = src_model(tcl_j), cov = Csrc, population = coh,
               dose = 200, times = OBS, stratify = "SEX")
    })
    names(sts) <- paste0("src", 1:3)
    fit <- nlmixr2(pooled, admData(), est = "adgh",
                   control = adghControl(studies = do.call(admStudies, sts),
                                         print = 0L, cores = 2L))
    p <- fit$parFixedDf
    row <- data.frame(rep=.rep, bcrcl=p["bcrcl","Estimate"],
                      bcrcl_se=p["bcrcl","SE"], bsex=p["bsex","Estimate"],
                      om=as.numeric(fit$omega[1,1]), cov_method=fit$covMethod)
    res[[length(res)+1L]] <- row
    write.csv(do.call(rbind,res),
              Sys.getenv("CSV", "C:/Users/hidde/.claude/jobs/319faff2/tmp/exp_exact.csv"),
              row.names = FALSE)
    cat(sprintf("rep %d ok: bcrcl %.4f (se %.4f)  om %.4f  %s\n",
                .rep, row$bcrcl, row$bcrcl_se, row$om, row$cov_method))
  }, silent = TRUE)
  if (inherits(ok,"try-error")) cat("rep",.rep,"FAILED:",
      conditionMessage(attr(ok,"condition")), "\n")
}
cat("\n== done\n")
