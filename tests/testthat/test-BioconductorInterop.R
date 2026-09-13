makeCardinalFixture <- function(nFeatures = 8L, nPixels = 6L) {
  mz <- 100 + seq_len(nFeatures) / 100
  coordinates <- data.frame(
    x = seq_len(nPixels),
    y = rep(seq_len(2L), length.out = nPixels)
  )
  Cardinal::MSImagingExperiment(
    spectraData = matrix(
      seq_len(nFeatures * nPixels),
      nrow = nFeatures,
      ncol = nPixels
    ),
    featureData = Cardinal::MassDataFrame(
      mz = mz,
      annotation = paste0("feature", seq_len(nFeatures))
    ),
    pixelData = Cardinal::PositionDataFrame(
      coord = coordinates,
      run = factor(rep(c("run1", "run2"), length.out = nPixels)),
      tissue = rep(TRUE, nPixels)
    )
  )
}

test_that("Cardinal and SpatialExperiment use formal S4 coercions", {
  cardinal <- makeCardinalFixture()
  spatial <- asSpatialExperiment(cardinal)

  expect_s4_class(spatial, "SpatialExperiment")
  expect_equal(dim(spatial), c(8L, 6L))
  expect_equal(
    unname(SpatialExperiment::spatialCoords(spatial)[, "x"]),
    seq_len(6L)
  )
  expect_equal(
    SummarizedExperiment::rowData(spatial)$mz,
    Cardinal::mz(cardinal)
  )
  expect_true(methods::is(methods::as(cardinal, "SpatialExperiment"),
                          "SpatialExperiment"))

  roundTrip <- methods::as(spatial, "MSImagingExperiment")
  expect_s4_class(roundTrip, "MSImagingExperiment")
  expect_equal(Cardinal::mz(roundTrip), Cardinal::mz(cardinal))
  expect_equal(
    unname(as.matrix(Cardinal::spectra(roundTrip))),
    unname(as.matrix(Cardinal::spectra(cardinal)))
  )
})

test_that("matrix import returns SpatialExperiment by default", {
  path <- tempfile(fileext = ".csv")
  table <- data.frame(
    x = c(1, 2, 1),
    y = c(1, 1, 2),
    `100.1` = c(1, 2, 3),
    `101.2` = c(4, 5, 6),
    check.names = FALSE
  )
  utils::write.csv(table, path, row.names = FALSE)
  object <- readSMMatrix(path, verbose = FALSE)

  expect_s4_class(object, "SpatialExperiment")
  expect_equal(dim(object), c(2L, 3L))
  expect_equal(SummarizedExperiment::rowData(object)$mz, c(100.1, 101.2))
  expect_equal(
    unname(SpatialExperiment::spatialCoords(object)),
    unname(as.matrix(table[, c("x", "y")]))
  )
  expect_false(any(c("x", "y", "x_coord", "y_coord") %in%
    colnames(SummarizedExperiment::colData(object))))
  expect_identical(object$orig.ident, rep("SpaMTP", ncol(object)))
})

test_that("binning and normalization dispatch on SpatialExperiment", {
  spatial <- asSpatialExperiment(makeCardinalFixture())
  binned <- binSpaMTP(
    spatial,
    resolution = 0.02,
    units = "mz",
    method = "sum"
  )
  expect_s4_class(binned, "SpatialExperiment")
  expect_lt(nrow(binned), nrow(spatial))
  expect_equal(ncol(binned), ncol(spatial))

  normalized <- normalizeSMData(spatial, "LogNormalize", verbose = FALSE)
  expect_true("logcounts" %in% SummarizedExperiment::assayNames(normalized))
  expect_length(SingleCellExperiment::sizeFactors(normalized), ncol(spatial))
})

test_that("paired transcriptomes use altExp and Bioconductor reductions", {
  spatial <- normalizeSMData(
    asSpatialExperiment(makeCardinalFixture(12L, 8L)),
    "LogNormalize",
    verbose = FALSE
  )
  transcriptome <- matrix(
    seq_len(80L),
    nrow = 10L,
    dimnames = list(paste0("gene", seq_len(10L)), rev(colnames(spatial)))
  )
  spatial <- addTranscriptome(spatial, transcriptome)
  expect_identical(SingleCellExperiment::altExpNames(spatial), "transcriptome")
  expect_identical(
    colnames(SingleCellExperiment::altExp(spatial, "transcriptome")),
    colnames(spatial)
  )

  integrated <- suppressWarnings(multiOmicIntegration(
    spatial,
    reduction.list = list("spm.pca", "spt.pca"),
    dims.list = list(1:2, 1:2),
    verbose = FALSE
  ))
  expect_equal(
    dim(SingleCellExperiment::reducedDim(integrated, "integrated")),
    c(8L, 4L)
  )
})

test_that("images are stored through SpatialExperiment imgData", {
  spatial <- asSpatialExperiment(makeCardinalFixture())
  imagePath <- tempfile(fileext = ".png")
  png::writePNG(array(1, dim = c(2, 2, 3)), imagePath)
  spatial <- addSpatialImage(
    spatial,
    imageSource = imagePath,
    sampleId = unique(spatial$sample_id)[[1L]],
    imageId = "optical",
    load = FALSE
  )
  expect_equal(nrow(SpatialExperiment::imgData(spatial)), 1L)
  expect_identical(SpatialExperiment::imgData(spatial)$image_id, "optical")

  output <- tempfile()
  saveSpaMTPData(
    spatial,
    output,
    assay = "counts",
    slot = "counts",
    image = "optical",
    generate.h5 = FALSE,
    verbose = FALSE
  )
  expect_true(file.exists(file.path(output, "spatial", "tissue_lowres_image.png")))
  expect_true(file.exists(file.path(output, "spatial", "tissue_positions_list.csv")))
})

test_that("feature utilities and spatial plots use Bioconductor accessors", {
  spatial <- asSpatialExperiment(makeCardinalFixture())
  nearest <- findNearestMZ(spatial, 100.031)
  expect_identical(nearest, rownames(spatial)[[3L]])

  spatial <- binMetabolites(
    spatial,
    mzs = rownames(spatial)[1:2],
    assay = "counts",
    slot = "counts",
    bin_name = "combined"
  )
  expect_true("combined" %in% colnames(SummarizedExperiment::colData(spatial)))
  expect_s3_class(
    plotSpatialFeature(spatial, rownames(spatial)[[1L]]),
    "ggplot"
  )
  expect_s3_class(
    spatialMZPlot(spatial, mzs = 100.01, assay = "counts"),
    "ggplot"
  )
})

test_that("TMM and expression plotting reuse edgeR and scater", {
  spatial <- asSpatialExperiment(makeCardinalFixture(12L, 8L))
  spatial$sample <- rep(c("reference", "test"), each = 4L)
  normalized <- tmmNormalize(
    spatial,
    ident = "sample",
    refIdent = "reference",
    assay = "counts",
    slot = "counts"
  )
  expect_true("normcounts" %in% SummarizedExperiment::assayNames(normalized))
  expect_named(S4Vectors::metadata(normalized)$tmm$factors, c("reference", "test"))
  expect_s3_class(
    mzViolinPlot(
      normalized,
      group.by = "sample",
      mzs = rownames(normalized)[[1L]],
      assay = "counts",
      slot = "normcounts"
    ),
    "ggplot"
  )
})

test_that("merged modalities are represented as an altExp", {
  spatial <- normalizeSMData(
    asSpatialExperiment(makeCardinalFixture(12L, 8L)),
    "LogNormalize",
    verbose = FALSE
  )
  transcriptome <- matrix(
    seq_len(80L),
    nrow = 10L,
    dimnames = list(paste0("gene", seq_len(10L)), colnames(spatial))
  )
  spatial <- addTranscriptome(spatial, transcriptome)
  spatial <- createMergedModalityAssay(
    spatial,
    assays.to.merge = c("main", "transcriptome")
  )
  expect_true("merged" %in% SingleCellExperiment::altExpNames(spatial))
  expect_equal(
    nrow(SingleCellExperiment::altExp(spatial, "merged")),
    nrow(spatial) + nrow(transcriptome)
  )
})

test_that("pathway scores are represented as an altExp", {
  spatial <- asSpatialExperiment(makeCardinalFixture(6L, 6L))
  rownames(spatial) <- paste0("RAMP_C_", seq_len(nrow(spatial)))
  database <- list(
    analytehaspathway = data.frame(
      pathwayRampId = c("RAMP_P_1", "RAMP_P_1", "RAMP_P_2"),
      rampId = rownames(spatial)[1:3]
    ),
    pathway = data.frame(
      pathwayRampId = c("RAMP_P_1", "RAMP_P_2"),
      pathwayName = c("pathway one", "pathway two")
    )
  )
  spatial <- createPathwayObject(
    spatial,
    slot = "counts",
    database = database
  )
  expect_true("pathway" %in% SingleCellExperiment::altExpNames(spatial))
  pathway <- SingleCellExperiment::altExp(spatial, "pathway")
  expect_identical(SummarizedExperiment::assayNames(pathway), "pathwayScores")
  expect_equal(ncol(pathway), ncol(spatial))
})

test_that("precomputed alignment updates spatialCoords", {
  spatial <- asSpatialExperiment(makeCardinalFixture())
  original <- SpatialExperiment::spatialCoords(spatial)
  aligned <- data.frame(
    cell = colnames(spatial),
    x_aligned = original[, "x"] + 10,
    y_aligned = original[, "y"] - 5
  )
  output <- applySpatialAlignment(
    spatial,
    alignment = aligned,
    coordinate.columns = c("x_aligned", "y_aligned"),
    verbose = FALSE
  )
  expect_equal(
    unname(SpatialExperiment::spatialCoords(output)[, "x"]),
    unname(original[, "x"] + 10)
  )
  expect_true("spatial_alignment" %in% names(S4Vectors::metadata(output)))
})

test_that("dimensionality reduction uses reducedDims", {
  spatial <- asSpatialExperiment(makeCardinalFixture(8L, 8L))
  spatial <- runMetabolicPCA(
    spatial,
    npcs = 2L,
    assay = "counts",
    slot = "counts",
    reduction.name = "metabolicPCA",
    verbose = FALSE
  )
  expect_true("metabolicPCA" %in% SingleCellExperiment::reducedDimNames(spatial))

  spatial <- runSpatialGraphPCA(
    spatial,
    n_components = 2L,
    assay = "counts",
    slot = "counts",
    platform = "ST",
    n_neighbors = 2L,
    reduction_name = "spatialPCA",
    verbose = FALSE
  )
  expect_true("spatialPCA" %in% SingleCellExperiment::reducedDimNames(spatial))
  expect_true("SpatialKNN" %in% names(S4Vectors::metadata(spatial)$spatialGraphs))

  spatial <- getKmeanClusters(
    spatial,
    reduction = "spatialPCA",
    clusters = 2L,
    cluster.name = "region"
  )
  expect_true("region" %in% colnames(SummarizedExperiment::colData(spatial)))
})

test_that("Seurat interoperability is explicit and round-trips altExps", {
  skip_if_not_installed("SeuratObject")
  spatial <- asSpatialExperiment(makeCardinalFixture())
  transcriptome <- matrix(
    seq_len(30L),
    nrow = 5L,
    dimnames = list(paste0("gene", seq_len(5L)), colnames(spatial))
  )
  spatial <- addTranscriptome(spatial, transcriptome)
  seurat <- spatialExperimentToSeurat(spatial)
  restored <- seuratToSpatialExperiment(seurat)

  expect_s4_class(seurat, "Seurat")
  expect_s4_class(restored, "SpatialExperiment")
  expect_true("transcriptome" %in% SingleCellExperiment::altExpNames(restored))
  expected <- spatial[, colnames(restored), drop = FALSE]
  expect_equal(as.matrix(SummarizedExperiment::assay(restored)),
                 as.matrix(SummarizedExperiment::assay(expected)))
  expect_equal(unname(SpatialExperiment::spatialCoords(restored)),
                 unname(SpatialExperiment::spatialCoords(expected)))
  expect_equal(SummarizedExperiment::rowData(restored)$mz,
                 SummarizedExperiment::rowData(expected)$mz)

  # Seurat can keep its own cell order when subsetting. Compare by the
  # returned IDs, not by the requested index order.
  seurat <- seurat[, rev(colnames(seurat))]
  restored <- seuratToSpatialExperiment(seurat)
  expected <- spatial[, colnames(restored), drop = FALSE]
  expect_equal(as.matrix(SummarizedExperiment::assay(restored)),
                 as.matrix(SummarizedExperiment::assay(expected)))
  expect_equal(unname(SpatialExperiment::spatialCoords(restored)),
                 unname(SpatialExperiment::spatialCoords(expected)))
  expect_equal(as.matrix(SummarizedExperiment::assay(
    SingleCellExperiment::altExp(restored, "transcriptome"))),
    as.matrix(SummarizedExperiment::assay(
      SingleCellExperiment::altExp(expected, "transcriptome"))))
})
