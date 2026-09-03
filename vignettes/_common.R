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
  # Progress bars are removed HERE, on the finished document, because
  # `collapse = TRUE` merges a chunk's source and output into one block and
  # never calls the `output` hook -- so a filter registered there is simply not
  # run, which is how nine compile bars reached covariates.html while this file
  # claimed to strip them. The per-hook filters below still catch the
  # uncollapsed case; this catches every case.
  .adm_scrub(x)
})

# The estimators keep rxode2's progress bar off during fitting (the control's
# nDisplayProgress argument), but a one-off bar still prints while each model is
# *compiled*. Strip those animation frames from the rendered output.
#
# Two things this has to get right, each of which let bars through before:
#
#  * Match the bar by its OPENING and its filler, never by how it ends. A
#    compile bar is flushed with no terminator at all, arriving as
#    "[====|====|====" with neither "]" nor "done".
#  * Register on message and warning as well as output. rxode2 writes the bar to
#    the message connection, so a hook on `output` alone never sees it -- which
#    is why nine bars survived into covariates.html while this file claimed to
#    remove them.
.adm_scrub <- function(x) {
  # Key on the FILLER RUN, not on any bracket. A bar redraws with \r, so what
  # survives into the document may have lost its "[" ("|====|====|====] 0:00:00")
  # or its tail ("[====|====|===="). Eight or more consecutive "=" / "|" occurs
  # in no real output and in no kable, whose separator rules use "-".
  x <- gsub("[^\n]*[=|]{8,}[^\n]*\n?", "", x, perl = TRUE)
  # Belt and braces: if anything still emits colour despite the options above,
  # drop the escape sequences rather than let them reach the HTML.  pkgdown
  # rewrites the ESC byte to U+2029, so a leak shows up as a stray "[1m" in the
  # rendered page rather than as anything recognisably an escape code.
  gsub("\033\\[[0-9;]*m", "", x, perl = TRUE)
}

for (.h in c("output", "message", "warning")) {
  local({
    .default <- knitr::knit_hooks$get(.h)
    knitr::knit_hooks$set(structure(
      list(function(x, options) {
        x <- .adm_scrub(x)
        if (!nzchar(gsub("[[:space:]]", "", x))) return("")
        .default(x, options)
      }), names = .h))
  })
}
