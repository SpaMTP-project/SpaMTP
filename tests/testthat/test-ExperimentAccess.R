test_that("expression matrices are selected by exact character names", {
  object <- nativeFixture()
  SummarizedExperiment::assay(object, "1") <-
    SummarizedExperiment::assay(object, "counts") + 100
  expect_equal(.assayData(object, layer = "1"),
    SummarizedExperiment::assay(object, "1"))
  invalid <- list(1, NA_character_, character(), c("counts", "logcounts"),
    "", list("counts"), factor("counts"), matrix("counts"))
  for (name in invalid) {
    expect_error(.assayData(object, layer = name), "layer.*character name")
    expect_error(.setAssayData(object, .assayData(object), layer = name),
      "layer.*character name")
    expect_error(.experimentForAssay(object, name), "assay.*character name")
  }
  expect_error(normalizeSMData(object, slot = 1), "layer.*character name")
  expect_error(scaleSMData(object, slot = 1), "layer.*character name")
  expect_error(.assayData(object, layer = "count"), "not found")
  expect_identical(.experimentForAssay(object, NULL), object)
})

test_that("matrix replacement changes only the requested alternative experiment", {
  object <- nativeFixture()
  before <- object
  values <- log1p(.assayData(object, "transcriptome"))
  result <- .setAssayData(object, values, "transcriptome", "logcounts")
  expect_equal(SummarizedExperiment::assays(result),
    SummarizedExperiment::assays(before))
  expect_equal(.assayData(result, "transcriptome", "logcounts"), values)
  expect_equal(.assayData(result, "transcriptome", "counts"),
    .assayData(before, "transcriptome", "counts"))
  expect_identical(object, before)
  expect_true(methods::validObject(result))
})

test_that("merging uses common modality selectors without double-counting aliases", {
  object <- nativeFixture()
  SingleCellExperiment::mainExpName(object) <- "MSI"
  for (alias in c("main", "primary", "Spatial", "SPM", "counts", "MSI")) {
    merged <- createMergedModalityAssay(object, c(alias, "transcriptome"))
    expect_equal(nrow(SingleCellExperiment::altExp(merged, "merged")), 5L)
    expect_equal(SummarizedExperiment::assays(merged),
      SummarizedExperiment::assays(object))
  }
  for (alias in c("primary", "Spatial", "SPM", "counts", "MSI")) {
    expect_error(createMergedModalityAssay(object, c("main", alias)),
      "distinct modalities")
  }
  expect_error(createMergedModalityAssay(object, c("main", "missing")),
    "Unknown assay or modality")
  expect_error(createMergedModalityAssay(object, c(1, 2)), "distinct modalities")
  # An actual altExp name takes precedence over a legacy primary alias.
  SingleCellExperiment::altExp(object, "SPM") <-
    SingleCellExperiment::altExp(object, "transcriptome")
  merged <- createMergedModalityAssay(object, c("main", "SPM"))
  expect_equal(nrow(SingleCellExperiment::altExp(merged, "merged")), 5L)
})

test_that("native scaled values take precedence but legacy values remain readable", {
  object <- nativeFixture()
  old <- .assayData(object) + 10
  native <- .assayData(object) + 20
  SummarizedExperiment::assay(object, "scale.data") <- old
  SummarizedExperiment::assay(object, "scaled") <- native
  merged <- createMergedModalityAssay(object, c("main", "transcriptome"))
  values <- .assayData(merged, "merged", "scaled")
  expect_equal(unname(values[seq_len(nrow(object)), ]), unname(native))
  SummarizedExperiment::assay(object, "scaled") <- NULL
  merged <- createMergedModalityAssay(object, c("main", "transcriptome"))
  values <- .assayData(merged, "merged", "scaled")
  expect_equal(unname(values[seq_len(nrow(object)), ]), unname(old))
})

test_that("help pages distinguish expression names from S4 storage", {
  paths <- list.files(test_path("..", "..", "man"), pattern = "\\.Rd$",
    full.names = TRUE)
  docs <- if (length(paths)) {
    stats::setNames(lapply(paths, tools::parse_Rd), basename(paths))
  } else {
    tools::Rd_db("SpaMTP")
  }
  expect_true("experimentAccess.Rd" %in% names(docs))
  slotNames <- c("slot", "slots", "SM.slot", "ST.slot", "SM_slot", "ST_slot",
    "mz.slot", "pathway.slot")
  checked <- 0L
  for (name in names(docs)) {
    rd <- docs[[name]]
    text <- paste(unlist(rd), collapse = " ")
    expect_false(grepl("assay storage slot|assay slot|reduction slot|pixelData slot",
      text), info = name)
    sections <- Filter(function(x) identical(attr(x, "Rd_tag"), "\\arguments"), rd)
    for (section in sections) {
      items <- Filter(function(x) identical(attr(x, "Rd_tag"), "\\item"), section)
      for (item in items) {
        if (!trimws(paste(unlist(item[[1L]]), collapse = "")) %in% slotNames) next
        checked <- checked + 1L
        expect_match(paste(unlist(item[[2L]]), collapse = ""), "experimentAccess",
          info = name)
      }
    }
  }
  expect_gte(checked, 30L)
  expect_match(paste(unlist(docs[["plotPathways.Rd"]]), collapse = ""),
    'default = .*logcounts')
})
