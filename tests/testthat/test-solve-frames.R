# Every params frame handed to rxSolve must have had the study's covariates
# appended to it, IN THE SAME FUNCTION.
#
# This is a static check, and it exists because the failure it catches is
# invisible to every dynamic one. `.admCovCols()`/`.admCovColsTiled()` is spelled
# out at fifteen frame builders across simulate.R, admc.R, adirmc.R and plot.R,
# and `.admMakeParamsList()` does not take the study, so it cannot append them
# itself. Ten separate silent failures have come from a site that was missed:
#
#   simulate.R      the sens frame          -> whole fit aborted, untryCatch'd
#   simulate.R      the order-2 sens frame  -> struct-theta gradient silently FD
#   simulate.R      the joint frame         -> Inf objective, no diagnosis
#   adirmc.R        the three proposal solves
#   admc.R          .admNLLBatch            -> covMethod = "r" gave no covariance
#   admc.R          .admNLLBatch stride     -> died at the last Hessian step
#   admc.R          .admGradBatch (x2)      -> constant-ZERO Hessian row for the
#                                              covariate coefficient, valid TRUE
#
# Each was found by measurement, one at a time. The point of this test is that
# the ELEVENTH is found by the suite instead.
.frame_audit <- function() {
  files <- list.files(file.path("..", "..", "R"), pattern = "[.]R$",
                      full.names = TRUE)
  if (!length(files)) files <- list.files("../../R", pattern = "[.]R$",
                                          full.names = TRUE)
  out <- character(0)
  for (f in files) {
    src <- readLines(f, warn = FALSE)
    # top-level function definitions split the file into scopes
    starts <- grep("^[.A-Za-z][A-Za-z0-9._]* *<- *function", src)
    if (!length(starts)) next
    ends <- c(starts[-1L] - 1L, length(src))
    for (k in seq_along(starts)) {
      body <- src[starts[k]:ends[k]]
      fn   <- sub(" *<-.*$", "", src[starts[k]])
      # the variable each rxSolve is given as `params =`
      pl <- regmatches(body, regexpr(
        "params *= *(as[.]data[.]frame[(])?[A-Za-z.][A-Za-z0-9._]*", body))
      pl <- unique(sub("params *= *(as[.]data[.]frame[(])?", "", pl))
      pl <- setdiff(pl, c("params", "pm", "pl"))   # passed in, appended by caller
      if (!length(pl)) next
      # ... and the variables a covariate append writes to
      ap <- regmatches(body, regexpr(
        "[A-Za-z.][A-Za-z0-9._]* *<- *[.]admCovCols(Tiled)?[(]", body))
      ap <- unique(sub(" *<-.*$", "", ap))
      # A frame COPIED from an appended one carries the columns with it
      # (`pdf_lo <- pdf_hi` is the central-difference twin), so follow plain
      # variable-to-variable assignment until it stops adding anything.
      repeat {
        cp <- regmatches(body, regexpr(
          "^ *[A-Za-z.][A-Za-z0-9._]* *<- *[A-Za-z.][A-Za-z0-9._]* *$", body))
        add <- character(0)
        for (l in cp) {
          lhs <- trimws(sub(" *<-.*$", "", l))
          rhs <- trimws(sub("^.*<- *", "", l))
          if (rhs %in% ap && !lhs %in% ap) add <- c(add, lhs)
        }
        if (!length(add)) break
        ap <- c(ap, add)
      }
      miss <- setdiff(pl, ap)
      if (length(miss))
        out <- c(out, sprintf("%s in %s (params = %s)",
                              fn, basename(f), paste(miss, collapse = ", ")))
    }
  }
  out
}

test_that("every rxSolve params frame gets the study's covariate columns", {
  skip_on_cran()
  bad <- .frame_audit()
  # as.data.frame(mat) wrappers and frames appended by the caller are excluded
  # above; anything left is a builder that never appends.
  expect_identical(bad, character(0))
})
