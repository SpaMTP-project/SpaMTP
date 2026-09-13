rawCardinalFixture <- function() {
  Cardinal::MSImagingArrays(
    centroided = TRUE, continuous = FALSE,
    spectraData = list(
      intensity = list(c(1, 2, 3), c(4, 5, 6, 7)),
      mz = list(c(100, 101, 102), c(100, 100.25, 101, 104))),
    pixelData = Cardinal::PositionDataFrame(
      coord = data.frame(x = c(1, 2), y = c(1, 1)),
      run = factor(c("section", "section"))))
}

test_that("unaligned raw Cardinal spectra are binned before SPE conversion", {
  raw <- rawCardinalFixture()
  expect_s4_class(raw, "MSImagingArrays")
  expect_identical(asCardinal(raw), raw)
  expect_error(asSpatialExperiment(raw), "must be binned first")
  result <- binSpaMTP(raw, resolution = 0.25, units = "mz")
  expect_s4_class(result, "SpatialExperiment")
  expect_equal(ncol(result), 2L)
  expect_true(max(SummarizedExperiment::rowData(result)$mz) >= 104)
  expect_equal(as.numeric(Matrix::colSums(SummarizedExperiment::assay(result))),
                 c(6, 22))
  expect_error(binSpaMTP(raw, resolution = -1), "resolution")
})

test_that("imzML import stays raw unless binning is explicitly requested", {
  raw <- rawCardinalFixture()
  directory <- withr::local_tempdir()
  path <- file.path(directory, "tiny.imzML")
  Cardinal::writeImzML(raw, path, bundle = FALSE, verbose = FALSE)
  imported <- loadSM(path, verbose = FALSE)
  expect_s4_class(imported, "MSImagingArrays")
  expect_equal(lengths(Cardinal::mz(imported)), c(3L, 4L))
  binned <- loadSM(path, resolution = 0.25, units = "mz", verbose = FALSE)
  expect_s4_class(binned, "SpatialExperiment")
  expect_equal(as.numeric(Matrix::colSums(SummarizedExperiment::assay(binned))),
                 c(6, 22))
  expect_error(loadSM(path, returnType = "SpatialExperiment"), "resolution")
})

test_that("Cardinal normalization queues work that conversion actually executes", {
  raw <- rawCardinalFixture()
  normalized <- normalizeSMData(raw, verbose = FALSE)
  expect_s4_class(normalized, "MSImagingArrays")
  expected <- Cardinal::bin(Cardinal::process(
    Cardinal::normalize(raw, method = "tic", verbose = FALSE), verbose = FALSE),
    resolution = 0.25, units = "mz", mass.range = c(100, 104))
  result <- binSpaMTP(normalized, resolution = 0.25, units = "mz")
  expect_equal(as.matrix(SummarizedExperiment::assay(result)),
                 as.matrix(Cardinal::spectra(expected)), ignore_attr = TRUE)
  expect_error(normalizeSMData(raw, "RC"), "Cardinal input supports TIC")
  aligned <- asCardinal(nativeFixture())
  normalized <- normalizeSMData(aligned, verbose = FALSE)
  expected <- Cardinal::spectra(Cardinal::process(normalized, verbose = FALSE))
  expect_equal(as.matrix(SummarizedExperiment::assay(asSpatialExperiment(normalized))),
                 as.matrix(expected), ignore_attr = TRUE)
})
