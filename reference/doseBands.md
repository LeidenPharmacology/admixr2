# A dosing table by covariate band

Transcribes a paper's dosing table, e.g. "CrCl below 30 mL/min: 100 mg;
otherwise 200 mg", for the `dose` argument of
[`admStudy()`](https://leidenpharmacology.github.io/admixr2/reference/admStudy.md):
`dose = doseBands(CRCL, breaks = 30, amt = c(100, 200))`.

## Usage

``` r
doseBands(covariate, breaks, amt, ...)
```

## Arguments

- covariate:

  The covariate the bands are cut on, bare or quoted (`CRCL` or
  `"CRCL"`); a variable holding the name also works. It must be
  declared, as a continuous covariate, in the study's `population`.

- breaks:

  Increasing cut points between bands.

- amt:

  Amount for each band, lowest band first; one more than `breaks`. A
  band that gets several doses at once (an oral dose plus an infusion)
  takes a list, one element per band, each holding that band's amounts.

- ...:

  How each dose is given, as
  [`rxode2::et()`](https://nlmixr2.github.io/rxode2/reference/et.html)
  arguments: `time`, `cmt`, `ii`, `addl`, `dur`, `rate`, `ss`, `until`;
  length 1 (all bands) or one per band. Without them each band gets one
  bolus at time 0.

## Value

An `admDoseBands` object, to pass as `dose` to
[`admStudy()`](https://leidenpharmacology.github.io/admixr2/reference/admStudy.md).

## Details

`breaks` are the cut points between bands and `amt` the amount for each
band, lowest first, so there is always one more amount than cut points
and the bands cover the whole range with no gaps. A value equal to a cut
point falls in the band above it. Anything else about how each dose is
given — `ii`, `addl`, `dur`, `rate`, `cmt`, ... — is written as in
[`rxode2::et()`](https://nlmixr2.github.io/rxode2/reference/et.html).
Each is one value, shared by every band, or one value per band, so the
regimen itself can depend on the covariate (`ii = c(24, 12)`). For
anything Several simultaneous doses in a band: give `amt` as a list and
the per-dose details as vectors inside lists (a length-1 list is shared
by all bands), e.g.
`amt = list(c(500, 100), c(1000, 200)), cmt = list(c("depot", "centr")), dur = list(c(0, 1))`.
For anything that is a formula of a covariate, pass a function as `dose`
that returns an
[`rxode2::et()`](https://nlmixr2.github.io/rxode2/reference/et.html)
table.

The study is expanded into one study per band: the covariate restricted
to the band (and to the study's enrolled `range`, if given), `n` split
by the band's share of the population, and that band's dose, so the
analysis model needs no dose rule. Covariates correlated with the banded
one are conditioned on the band.

## See also

[`admStudy()`](https://leidenpharmacology.github.io/admixr2/reference/admStudy.md);
a function of covariates as `dose` for dosing that is proportional to
(or a formula of) a covariate, e.g. mg/kg.

## Examples

``` r
doseBands(CRCL, breaks = 30, amt = c(100, 200))
#> dose by CRCL band: CRCL < 30: 100; CRCL >= 30: 200
# a 1-h infusion every 12 h, seven doses, in every band
doseBands(CRCL, breaks = 30, amt = c(500, 1000), ii = 12, addl = 6, dur = 1)
#> dose by CRCL band: CRCL < 30: 500; CRCL >= 30: 1000  (ii = 12, addl = 6, dur = 1)
# the interval depends on the band too: daily below 30, twice daily above
doseBands(CRCL, breaks = 30, amt = c(500, 1000), ii = c(24, 12), addl = c(2, 6))
#> dose by CRCL band: CRCL < 30: 500; CRCL >= 30: 1000  (ii = 24,12, addl = 2,6)
# every regimen detail depends on the band: dose, interval, duration, route
doseBands(CRCL, breaks = c(25, 40, 60), amt = c(450, 450, 450, 900),
          ii = c(48, 24, 12, 12), addl = 16, dur = c(0, 0, 1, 1),
          cmt = c("depot", "depot", "centr", "centr"))
#> dose by CRCL band: CRCL < 25: 450; 25 to < 40: 450; 40 to < 60: 450; CRCL >= 60: 900  (ii = 48,24,12,12, addl = 16, dur = 0,0,1,1, cmt = depot,depot,centr,centr)
# an oral dose AND a 1-h infusion at the same time, in every band
doseBands(CRCL, breaks = 30, amt = list(c(500, 100), c(1000, 200)),
          cmt = list(c("depot", "centr")), dur = list(c(0, 1)))
#> dose by CRCL band: CRCL < 30: 500+100; CRCL >= 30: 1000+200  (cmt = depot, centr, dur = 0, 1)
# a formula of the covariate instead: any et() argument may read it
f <- function(WT, CRCL) {
  rxode2::et(amt = 5 * WT, cmt = "depot") |>
    rxode2::et(amt = 2 * WT, cmt = "centr", dur = 1, ii = 12 + CRCL / 10, addl = 3)
}
# admStudy(model, est, n, times, population, dose = f)
```
