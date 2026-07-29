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
