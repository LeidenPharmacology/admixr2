# Shared knitr setup for every admixr2 vignette and pkgdown article.
#
# Sourced from each document's setup chunk as
#
#   source(Filter(file.exists, c("_common.R", "../_common.R", "vignettes/_common.R"))[1])
#
# which resolves whether the working directory is vignettes/ (R CMD build, and
# pkgdown for an ordinary vignette), vignettes/articles/ (pkgdown for a Blog
# article), or the package root.
#
# This block used to be copied verbatim into four vignettes.  The seven that did
# not have it rendered nlmixr2's cli-coloured fit print as raw escape codes on
# the pkgdown site -- visible on admixr2.html, datagen.html and
# multiple-studies.html, i.e. exactly the three that call print(fit).  Keeping it
# in one file is what stops a newly added article from missing it.

# nlmixr2's fit print method uses cli colour; when a build forces colour on (e.g.
# pkgdown/CI) the raw ANSI escapes leak into the HTML. Force plain-text output
# for this document only, saving the previous global state...
.adm_opts    <- options(cli.num_colors = 1, crayon.enabled = FALSE)
.adm_nocolor <- Sys.getenv("NO_COLOR", unset = NA)
Sys.setenv(NO_COLOR = "1")

# ...and restoring it once the document has finished knitting, so the setting
# does not leak into other articles when pkgdown builds them in one process.
knitr::knit_hooks$set(document = function(x) {
  options(.adm_opts)
  if (is.na(.adm_nocolor)) Sys.unsetenv("NO_COLOR") else Sys.setenv(NO_COLOR = .adm_nocolor)
  x
})

knitr::knit_hooks$set(output = local({
  .default <- knitr::knit_hooks$get("output")
  function(x, options) {
    # The estimators keep rxode2's progress bar off during fitting (the control's
    # nDisplayProgress argument), but a one-off bar still prints while each model
    # is *compiled*. Strip those animation frames from the rendered output.
    x <- gsub("[^\n]*\\[[=| ]+\\][^\n]*\n?", "", x, perl = TRUE)
    # Belt and braces: if anything still emits colour despite the options above,
    # drop the escape sequences rather than let them reach the HTML.  pkgdown
    # rewrites the ESC byte to U+2029, so a leak shows up as a stray "[1m" in the
    # rendered page rather than as anything recognisably an escape code.
    x <- gsub("\033\\[[0-9;]*m", "", x, perl = TRUE)
    .default(x, options)
  }
}))

# ---------------------------------------------------------------------------
# The examplomycin running example.
#
# Five vignettes (admixr2, multiple-studies, diagnostic-plots, advanced and
# estimator-comparison) open on the same thing: reshape the shipped individual
# records into a subjects x times matrix, take E and V off it, and fit the
# two-compartment model the data were simulated from.  Copied into each, the
# model had already drifted -- two copies carried label() calls and two did
# not -- so it lives here and each vignette shows only what it is about.

# Subjects x times matrix of observed DV, one row per subject.
admVignetteDvMatrix <- function() {
  obs   <- admixr2::examplomycin
  obs   <- obs[obs$EVID == 0, ]
  obs   <- obs[order(obs$ID, obs$TIME), ]
  times <- sort(unique(obs$TIME))
  ids   <- unique(obs$ID)
  dv    <- matrix(NA_real_, nrow = length(ids), ncol = length(times),
                  dimnames = list(NULL, times))
  for (i in seq_along(ids)) {
    sub     <- obs[obs$ID == ids[i], ]
    dv[i, ] <- sub$DV[order(sub$TIME)]
  }
  dv
}

# E, V and n for a set of rows of that matrix -- all 500 subjects by default.
# V uses the ML (denominator n) convention, which is what admixr2 assumes
# unless a study declares `v_denom = "unbiased"`.
admVignetteStats <- function(dv = admVignetteDvMatrix(), rows = seq_len(nrow(dv))) {
  dv <- dv[rows, , drop = FALSE]
  list(E = colMeans(dv),
       V = stats::cov.wt(dv, method = "ML")$cov,
       n = nrow(dv),
       times = as.numeric(colnames(dv)))
}

# The model examplomycin was simulated from: two compartments, first-order
# absorption, IIV on all five parameters, proportional residual error.
admVignetteModel <- function() {
  ini({
    tcl     <- log(5)  ; label("Log clearance (L/hr)")
    tv1     <- log(10) ; label("Log central volume (L)")
    tv2     <- log(30) ; label("Log peripheral volume (L)")
    tq      <- log(10) ; label("Log inter-compartmental CL (L/hr)")
    tka     <- log(1)  ; label("Log absorption rate constant (1/hr)")
    prop.sd <- c(0, 0.2); label("Proportional residual error SD")
    eta.cl ~ 0.09
    eta.v1 ~ 0.09
    eta.v2 ~ 0.09
    eta.q  ~ 0.09
    eta.ka ~ 0.09
  })
  model({
    cl <- exp(tcl + eta.cl)
    v1 <- exp(tv1 + eta.v1)
    v2 <- exp(tv2 + eta.v2)
    q  <- exp(tq  + eta.q)
    ka <- exp(tka + eta.ka)
    d/dt(depot)      <- -ka * depot
    d/dt(central)    <- ka * depot - (cl/v1 + q/v1) * central + (q/v2) * peripheral
    d/dt(peripheral) <- (q/v1) * central - (q/v2) * peripheral
    cp <- central / v1
    cp ~ prop(prop.sd)
  })
}
