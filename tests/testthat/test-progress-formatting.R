# Helper: build a pinfo with known column names for deterministic width checks.
make_pinfo_simple <- function() admixr2:::.admParseIniDf(make_inidf_1eta())
make_pinfo_2eta   <- function() admixr2:::.admParseIniDf(make_inidf_2eta())

# ---- .admProgressHeader ------------------------------------------------------

test_that("Header has three lines: sep + header row + sep", {
  pinfo <- make_pinfo_simple()
  hdr   <- admixr2:::.admProgressHeader(pinfo)
  lines <- strsplit(hdr, "\n")[[1]]
  expect_equal(length(lines), 3L)
})

test_that("Header separator lines start and end with '+'", {
  pinfo <- make_pinfo_simple()
  lines <- strsplit(admixr2:::.admProgressHeader(pinfo), "\n")[[1]]
  expect_match(lines[1], "^\\+.*\\+$")
  expect_match(lines[3], "^\\+.*\\+$")
})

test_that("Header row starts and ends with '|'", {
  pinfo <- make_pinfo_simple()
  lines <- strsplit(admixr2:::.admProgressHeader(pinfo), "\n")[[1]]
  expect_match(lines[2], "^\\|.*\\|$")
})

test_that("Header contains '-2LL'", {
  pinfo <- make_pinfo_simple()
  expect_true(grepl("-2LL", admixr2:::.admProgressHeader(pinfo)))
})

test_that("Header contains all struct names", {
  pinfo <- make_pinfo_2eta()
  hdr   <- admixr2:::.admProgressHeader(pinfo)
  for (nm in pinfo$struct_names) expect_true(grepl(nm, hdr, fixed = TRUE))
})

test_that("Header contains sigma names", {
  pinfo <- make_pinfo_simple()
  hdr   <- admixr2:::.admProgressHeader(pinfo)
  for (nm in pinfo$sigma_names) expect_true(grepl(nm, hdr, fixed = TRUE))
})

test_that("Header separators have equal width", {
  pinfo <- make_pinfo_simple()
  lines <- strsplit(admixr2:::.admProgressHeader(pinfo), "\n")[[1]]
  expect_equal(nchar(lines[1]), nchar(lines[3]))
})

test_that("bottom=FALSE omits closing separator (2 lines instead of 3)", {
  pinfo <- make_pinfo_simple()
  hdr   <- admixr2:::.admProgressHeader(pinfo, bottom = FALSE)
  lines <- strsplit(hdr, "\n")[[1]]
  expect_equal(length(lines), 2L)
})

# ---- .admProgressRow ---------------------------------------------------------

test_that("Row starts and ends with '|'", {
  pinfo <- make_pinfo_simple()
  vec   <- admixr2:::.admBuildOptVec(pinfo)
  nll_val <- 123.45
  row   <- admixr2:::.admProgressRow("1", nll_val, vec$p0, pinfo)
  expect_match(row, "^\\|.*\\|$")
})

test_that("Row pipe count matches header pipe count", {
  pinfo <- make_pinfo_simple()
  vec   <- admixr2:::.admBuildOptVec(pinfo)
  hdr_lines <- strsplit(admixr2:::.admProgressHeader(pinfo), "\n")[[1]]
  row   <- admixr2:::.admProgressRow("1", 50.0, vec$p0, pinfo)

  n_pipes_hdr <- nchar(gsub("[^|]", "", hdr_lines[2]))
  n_pipes_row <- nchar(gsub("[^|]", "", row))
  expect_equal(n_pipes_row, n_pipes_hdr)
})

test_that("Row width matches header separator width", {
  pinfo <- make_pinfo_simple()
  vec   <- admixr2:::.admBuildOptVec(pinfo)
  hdr_lines <- strsplit(admixr2:::.admProgressHeader(pinfo), "\n")[[1]]
  row   <- admixr2:::.admProgressRow("eval 1", 50.0, vec$p0, pinfo)
  expect_equal(nchar(row), nchar(hdr_lines[1]))
})

test_that("Row contains NLL value formatted with 2 decimal places", {
  pinfo <- make_pinfo_simple()
  vec   <- admixr2:::.admBuildOptVec(pinfo)
  row   <- admixr2:::.admProgressRow("1", 99.12345, vec$p0, pinfo)
  expect_true(grepl("99.12", row, fixed = TRUE))
})

# ---- .admProgressDivider / .admProgressPhase / .admProgressRestart ----------

test_that("Divider starts and ends with '+'", {
  pinfo <- make_pinfo_simple()
  div   <- admixr2:::.admProgressDivider("label", pinfo)
  expect_match(div, "^\\+.*\\+$")
})

test_that("Divider width matches header separator width", {
  pinfo <- make_pinfo_simple()
  hdr_lines <- strsplit(admixr2:::.admProgressHeader(pinfo), "\n")[[1]]
  div   <- admixr2:::.admProgressDivider("test", pinfo)
  expect_equal(nchar(div), nchar(hdr_lines[1]))
})

test_that("Divider contains the label text", {
  pinfo <- make_pinfo_simple()
  div   <- admixr2:::.admProgressDivider("my label", pinfo)
  expect_true(grepl("my label", div, fixed = TRUE))
})

test_that(".admProgressRestart output contains restart indices", {
  pinfo <- make_pinfo_simple()
  div   <- admixr2:::.admProgressRestart(2L, 5L, pinfo)
  expect_true(grepl("2", div, fixed = TRUE))
  expect_true(grepl("5", div, fixed = TRUE))
})

test_that(".admProgressPhase output contains phase step", {
  pinfo <- make_pinfo_simple()
  div   <- admixr2:::.admProgressPhase(1L, "box", 0.5, pinfo)
  expect_true(grepl("0.50", div))
})

# ---- .admProgressTimingRow ---------------------------------------------------

test_that("Timing row starts and ends with '|'", {
  pinfo <- make_pinfo_simple()
  tr    <- admixr2:::.admProgressTimingRow(12.3, pinfo)
  expect_match(tr, "^\\|.*\\|$")
})

test_that("Timing row width matches header separator width", {
  pinfo <- make_pinfo_simple()
  hdr_lines <- strsplit(admixr2:::.admProgressHeader(pinfo), "\n")[[1]]
  tr <- admixr2:::.admProgressTimingRow(5.5, pinfo)
  expect_equal(nchar(tr), nchar(hdr_lines[1]))
})

test_that("Timing row contains the elapsed time", {
  pinfo <- make_pinfo_simple()
  tr <- admixr2:::.admProgressTimingRow(7.89, pinfo)
  expect_true(grepl("7.9", tr))
})

# ---- Consistency: 2-eta model ------------------------------------------------

test_that("All formatting functions produce consistent widths for 2-eta model", {
  pinfo <- make_pinfo_2eta()
  vec   <- admixr2:::.admBuildOptVec(pinfo)

  hdr_lines <- strsplit(admixr2:::.admProgressHeader(pinfo), "\n")[[1]]
  w <- nchar(hdr_lines[1])

  expect_equal(nchar(admixr2:::.admProgressDivider("x", pinfo)), w)
  expect_equal(nchar(admixr2:::.admProgressRow("1", 0.0, vec$p0, pinfo)), w)
  expect_equal(nchar(admixr2:::.admProgressTimingRow(1.0, pinfo)), w)
})

# ---- Large NLL switches to scientific notation --------------------------------

test_that("Row with large NLL uses scientific notation", {
  pinfo <- make_pinfo_simple()
  vec   <- admixr2:::.admBuildOptVec(pinfo)
  # Value too wide for 'f' format in 8-char column -> fall back to 'e'
  row   <- admixr2:::.admProgressRow("1", 1e15, vec$p0, pinfo)
  expect_true(grepl("e", row, ignore.case = TRUE))
})

# ---- Timing row: >= 60 seconds shows minutes ----------------------------------

test_that("Timing row >= 60s shows minutes", {
  pinfo <- make_pinfo_simple()
  tr    <- admixr2:::.admProgressTimingRow(90.0, pinfo)
  expect_true(grepl("min", tr, fixed = TRUE))
})

test_that("Timing row < 60s shows seconds", {
  pinfo <- make_pinfo_simple()
  tr    <- admixr2:::.admProgressTimingRow(45.0, pinfo)
  expect_true(grepl("sec", tr, fixed = TRUE))
})

# ---- 0-eta model (no omega diagonal) ------------------------------------------

test_that("Header and row work for 0-eta model", {
  pinfo <- admixr2:::.admParseIniDf(make_inidf_0eta())
  vec   <- admixr2:::.admBuildOptVec(pinfo)

  hdr_lines <- strsplit(admixr2:::.admProgressHeader(pinfo), "\n")[[1]]
  w   <- nchar(hdr_lines[1])
  row <- admixr2:::.admProgressRow("1", 50.0, vec$p0, pinfo)
  expect_equal(nchar(row), w)
})

# ---- .admProgressNames -------------------------------------------------------

test_that(".admProgressNames includes -2LL, struct, sigma, omega diagonal", {
  pinfo <- make_pinfo_2eta()
  nms   <- admixr2:::.admProgressNames(pinfo)
  expect_true("-2LL" %in% nms)
  for (nm in pinfo$struct_names) expect_true(nm %in% nms)
  for (nm in pinfo$sigma_names)  expect_true(nm %in% nms)
  diag_nms <- pinfo$eta_names[pinfo$chol_i[pinfo$chol_diag]]
  for (nm in diag_nms) expect_true(nm %in% nms)
})

test_that(".admProgressNames: 0-eta model has no omega entries", {
  pinfo <- admixr2:::.admParseIniDf(make_inidf_0eta())
  nms   <- admixr2:::.admProgressNames(pinfo)
  expect_equal(nms, c("-2LL", pinfo$struct_names, pinfo$sigma_names))
})
