test_that("Seurat namespace calls occur only in explicit adapters", {
  namespace <- asNamespace("SpaMTP")
  allowed <- c(".seuratLayer", ".seuratStoredMetadata", ".seuratCoordinates",
               "seuratToSingleCellExperiment", "spatialExperimentToSeurat")
  containsSeuratCall <- function(value) {
    tokens <- getParseData(parse(text = deparse(value), keep.source = TRUE))
    any(tokens$token == "SYMBOL_PACKAGE" &
          tokens$text %in% c("Seurat", "SeuratObject"))
  }
  names <- ls(namespace, all.names = TRUE)
  functions <- names[vapply(names, function(name) {
    value <- get(name, namespace)
    is.function(value) && containsSeuratCall(body(value))
  }, logical(1))]
  expect_setequal(functions, allowed)
  dependencies <- utils::packageDescription("SpaMTP")
  expect_false(grepl("\\bSeurat\\b", dependencies$Imports))
  expect_false(grepl("\\bSeurat\\b", dependencies$Suggests))
  expect_match(dependencies$Suggests, "SeuratObject")
})

test_that("core analysis rejects Seurat-shaped input at the boundary", {
  object <- structure(list(), class = "Seurat")
  calls <- list(
    function() normalizeSMData(object),
    function() binSpaMTP(object, 10),
    function() runMetabolicPCA(object),
    function() scaleSMData(object),
    function() runSpatialGraphPCA(object),
    function() multiOmicIntegration(object),
    function() subsetSPM(object),
    function() tmmNormalize(object, "region", "edge"),
    function() createPathwayObject(object),
    function() createPathwayAssay(object))
  for (call in calls) expect_error(call(), "convert|Convert")
  expect_error(readSMMatrix("missing.csv", returnType = "Seurat"), "arg")
  expect_error(loadSM("missing.imzML", returnType = "Seurat"), "arg")
  expect_error(loadMetaspace("not-requested", returnType = "Seurat"), "arg")
})

test_that("scaling preserves the input and acts on the selected modality", {
  object <- nativeFixture()
  original <- SummarizedExperiment::assays(object)
  result <- scaleSMData(object)
  expected <- t(base::scale(t(as.matrix(.assayData(object, layer = "logcounts")))))
  expected[!is.finite(expected)] <- 0
  expect_equal(unname(as.matrix(.assayData(result, layer = "scaled"))),
               unname(as.matrix(expected)), ignore_attr = TRUE)
  expect_equal(.assayData(result), .assayData(object))
  expect_equal(as.numeric(.assayData(result, layer = "scaled")[3, ]), rep(0, 6))
  result <- scaleSMData(object, assay = "transcriptome", slot = "counts")
  expect_equal(SummarizedExperiment::assays(result), original)
  expect_true("scaled" %in% SummarizedExperiment::assayNames(
    SingleCellExperiment::altExp(result, "transcriptome")))
  expect_error(scaleSMData(result, "transcriptome", "counts"), "not overwritten")
  expect_error(scaleSMData(object, center = NA), "TRUE or FALSE")
  expect_error(scaleSMData(object, slot = "absent"), "not found")
})

test_that("PCA and expression plots use the requested alternative experiment", {
  object <- nativeFixture()
  result <- suppressWarnings(runMetabolicPCA(object, assay = "transcriptome", npcs = 1))
  expected <- suppressWarnings(scater::runPCA(
    SingleCellExperiment::altExp(object, "transcriptome"),
    exprs_values = "counts", ncomponents = 1, name = "pca"))
  expect_equal(abs(SingleCellExperiment::reducedDim(result, "pca")),
               abs(SingleCellExperiment::reducedDim(expected, "pca")))
  expect_error(runMetabolicPCA(object, assay = "transcriptome", slot = "absent"),
               "not found")
  expect_error(runMetabolicPCA(object, npcs = 0), "npcs")
  result <- suppressWarnings(runMetabolicPCA(object, npcs = NULL,
                                               variance_explained_threshold = 0.8))
  expect_gte(sum(attr(SingleCellExperiment::reducedDim(result, "pca"), "percentVar")), 80)
  expect_s3_class(mzViolinPlot(object, assay = "transcriptome", mzs = "gene2",
                              group.by = "region"), "ggplot")
})

test_that("Seurat conversion preserves identities, metadata and named pixels", {
  skip_if_not_installed("SeuratObject", "5.0.0")
  source <- nativeFixture()
  SingleCellExperiment::colLabels(source) <- factor(source$region)
  S4Vectors::metadata(source)$db_3 <- data.frame(mz = 100, annotation = "test")
  seurat <- spatialExperimentToSeurat(source)
  sce <- seuratToSingleCellExperiment(seurat)
  expect_s4_class(sce, "SingleCellExperiment")
  expect_equal(as.character(SingleCellExperiment::colLabels(sce)),
               as.character(SingleCellExperiment::colLabels(source)))
  expect_equal(S4Vectors::metadata(sce)$db_3, S4Vectors::metadata(source)$db_3)
  xy <- SpatialExperiment::spatialCoords(source)[6:1, , drop = FALSE]
  rownames(xy) <- rev(colnames(source))
  result <- seuratToSpatialExperiment(seurat, coordinates = xy)
  expect_equal(unname(SpatialExperiment::spatialCoords(result)),
               unname(SpatialExperiment::spatialCoords(source)))
  expect_equal(.featureMetadata(result), .featureMetadata(source))
  expect_equal(as.matrix(.assayData(result)), as.matrix(.assayData(source)))
  expect_true("transcriptome" %in% SingleCellExperiment::altExpNames(result))
  expect_error(seuratToSpatialExperiment(seurat, coordinates = xy[-1, ]),
               "every selected pixel")
  expect_error(seuratToSpatialExperiment(seurat, image = "fov", coordinates = xy),
               "not both")
  xy[1, 1] <- Inf
  expect_error(seuratToSpatialExperiment(seurat, coordinates = xy), "finite")
})

test_that("conversion preserves non-syntactic metadata column names", {
  skip_if_not_installed("SeuratObject", "5.0.0")
  counts <- Matrix::Matrix(matrix(seq_len(12), 3,
    dimnames = list(paste0("mz-", 101:103), paste0("p", 1:4))), sparse = TRUE)
  cellName <- "2-Hydroxyglutric Aciduria (D And L Form)"
  featureName <- "mass error (ppm)"
  pixels <- data.frame(value = 1:4, row.names = colnames(counts))
  names(pixels) <- cellName
  features <- data.frame(value = c(0.1, 0.2, 0.3), row.names = rownames(counts))
  names(features) <- featureName
  seurat <- SeuratObject::CreateSeuratObject(counts, meta.data = pixels)
  seurat[["RNA"]] <- SeuratObject::AddMetaData(seurat[["RNA"]], features)
  sce <- seuratToSingleCellExperiment(seurat)
  expect_identical(SummarizedExperiment::colData(sce)[[cellName]], 1:4)
  expect_identical(SummarizedExperiment::rowData(sce)[[featureName]], c(0.1, 0.2, 0.3))
  coordinates <- cbind(x = 1:4, y = c(0, 1, 0, 1))
  rownames(coordinates) <- colnames(counts)
  spe <- seuratToSpatialExperiment(seurat, coordinates = coordinates)
  expect_identical(SummarizedExperiment::colData(spe)[[cellName]], 1:4)
  expect_identical(SummarizedExperiment::rowData(spe)[[featureName]], c(0.1, 0.2, 0.3))
  expect_identical(.cellMetadata(spe)[[cellName]], 1:4)
  expect_identical(.featureMetadata(spe)[[featureName]], c(0.1, 0.2, 0.3))
  expect_s3_class(plot3DFeature(spe, features = cellName, assays = "main"), "plotly")
  restored <- spatialExperimentToSeurat(spe)
  expect_identical(restored[[]][[cellName]], 1:4)
  expect_identical(restored[["Spatial"]][[]][[featureName]], c(0.1, 0.2, 0.3))
})

test_that("non-spatial and split-layer input cannot silently lose pixels", {
  skip_if_not_installed("SeuratObject", "5.0.0")
  counts <- Matrix::Matrix(matrix(seq_len(24), 4,
                   dimnames = list(paste0("g", 1:4), paste0("p", 1:6))), sparse = TRUE)
  seurat <- SeuratObject::CreateSeuratObject(counts)
  expect_s4_class(seuratToSingleCellExperiment(seurat), "SingleCellExperiment")
  expect_error(seuratToSpatialExperiment(seurat), "seuratToSingleCellExperiment")
  expect_error(seuratToSingleCellExperiment(seurat, layer = "count"), "exact")
  seurat[["RNA"]] <- SeuratObject::CreateAssay5Object(
    counts = list(a = counts[, 1:3], b = counts[, 4:6]))
  expect_error(seuratToSingleCellExperiment(seurat), "Join split layers")
  result <- seuratToSingleCellExperiment(seurat, layer = "counts.a")
  expect_equal(colnames(result), colnames(counts)[1:3])
  expect_equal(as.matrix(SummarizedExperiment::assay(result)), as.matrix(counts[, 1:3]))
})

test_that("the spatial accessor takes precedence over historical coordinate copies", {
  skip_if_not_installed("SeuratObject", "5.0.0")
  source <- nativeFixture()
  seurat <- spatialExperimentToSeurat(source)
  legacy <- data.frame(x_coord = 1000 + seq_len(ncol(source)),
    y_coord = 2000 + seq_len(ncol(source)), row.names = colnames(source))
  seurat <- SeuratObject::AddMetaData(seurat, legacy)
  before <- seurat
  result <- seuratToSpatialExperiment(seurat)
  expect_equal(unname(SpatialExperiment::spatialCoords(result)),
    unname(SpatialExperiment::spatialCoords(source)))
  expect_equal(result$x_coord, legacy$x_coord)
  expect_identical(seurat, before)
  provenance <- S4Vectors::metadata(result)$seurat_interoperability
  expect_identical(provenance$image, "fov")
  expect_identical(provenance$coordinate_source, "image")
  expect_identical(provenance$coordinate_columns, c("x", "y"))

  explicit <- legacy[, c("x_coord", "y_coord")]
  colnames(explicit) <- c("x", "y")
  result <- seuratToSpatialExperiment(seurat, coordinates = explicit)
  expect_equal(SpatialExperiment::spatialCoords(result), as.matrix(explicit))
  expect_identical(S4Vectors::metadata(result)$seurat_interoperability$coordinate_source,
    "coordinates")

  seurat[["fov"]] <- NULL
  result <- seuratToSpatialExperiment(seurat)
  expect_equal(SpatialExperiment::spatialCoords(result), as.matrix(explicit))
  expect_identical(S4Vectors::metadata(result)$seurat_interoperability$coordinate_source,
    "cell_metadata")
  expect_identical(S4Vectors::metadata(result)$seurat_interoperability$coordinate_columns,
    c("x_coord", "y_coord"))
})

test_that("metadata cannot silently disambiguate multiple spatial images", {
  skip_if_not_installed("SeuratObject", "5.0.0")
  seurat <- spatialExperimentToSeurat(nativeFixture())
  seurat$x <- seq_len(ncol(seurat))
  seurat$y <- seq_len(ncol(seurat))
  second <- seurat[["fov"]]
  SeuratObject::Key(second) <- "second_"
  seurat[["second"]] <- second
  expect_error(seuratToSpatialExperiment(seurat), "Multiple Seurat images")
  result <- seuratToSpatialExperiment(seurat, image = "second")
  expect_identical(S4Vectors::metadata(result)$seurat_interoperability$image, "second")
  result <- seuratToSpatialExperiment(seurat,
    coordinates = seurat[[]][, c("x", "y")])
  expect_equal(SpatialExperiment::spatialCoords(result),
    as.matrix(seurat[[]][, c("x", "y")]))
})

test_that("finite-radius centroids are read as one centre per pixel", {
  skip_if_not_installed("SeuratObject", "5.0.0")
  source <- nativeFixture()
  seurat <- spatialExperimentToSeurat(source)
  xy <- SpatialExperiment::spatialCoords(source)
  rownames(xy) <- colnames(source)
  centres <- SeuratObject::CreateCentroids(
    xy, nsides = 6L, radius = 0.4, theta = 0)
  seurat[["fov"]] <- SeuratObject::CreateFOV(
    coords = list(centroids = centres), assay = "Spatial")
  expanded <- SeuratObject::GetTissueCoordinates(seurat, image = "fov")
  expect_gt(nrow(expanded), ncol(source))
  result <- seuratToSpatialExperiment(seurat)
  expect_equal(unname(SpatialExperiment::spatialCoords(result)),
    unname(SpatialExperiment::spatialCoords(source)))
  expect_identical(colnames(result), colnames(source))
})

test_that("conversion reports unsupported metadata and derived alternative assays", {
  skip_if_not_installed("SeuratObject", "5.0.0")
  source <- nativeFixture()
  source <- createMergedModalityAssay(source, c("main", "transcriptome"))
  expect_warning(seurat <- spatialExperimentToSeurat(source), "Skipping.*merged")
  expect_false("merged" %in% SeuratObject::Assays(seurat))
  SeuratObject::Misc(seurat, slot = "nested") <- list(old = seurat)
  expect_warning(result <- seuratToSingleCellExperiment(seurat),
                  "Skipping Seurat-specific stored metadata: nested")
  expect_null(S4Vectors::metadata(result)$nested)
  expect_s4_class(result, "SingleCellExperiment")
})

test_that("native annotation and feature subsetting do not use Seurat assay indexing", {
  object <- nativeFixture()
  originalTranscriptome <- SingleCellExperiment::altExp(object, "transcriptome")
  selected <- subsetMZFeatures(object, c("constant", "a"), assay = "main")
  expect_identical(rownames(selected), c("constant", "a"))
  expect_equal(SingleCellExperiment::altExp(selected, "transcriptome"), originalTranscriptome)
  selected <- subsetMZFeatures(object, "gene2", assay = "transcriptome")
  expect_identical(rownames(selected), rownames(object))
  expect_identical(rownames(SingleCellExperiment::altExp(selected, "transcriptome")), "gene2")
  expect_error(subsetMZFeatures(object, "absent"), "selected experiment")
  custom <- data.frame(annotation = "target", mass = 100)
  selected <- addCustomMZAnnotations(object, custom, return.only.annotated = TRUE)
  expect_identical(rownames(selected), "a")
  expect_equal(SummarizedExperiment::rowData(selected)$all_IsomerNames, "target")
  expect_equal(SingleCellExperiment::altExp(selected, "transcriptome"), originalTranscriptome)
  target <- data.frame(formula = "C5H8O5", exactmass = 148.037173366,
                       id = "2HG", isomers_names = "D-2HG", max_exchangeable_protons = 2)
  index <- buildMZAnnotationIndex(target, adducts = "M+H")
  SummarizedExperiment::rowData(object)$mz <- c(index$expected_mz[1], 500, 600)
  result <- annotateSM(object, index = index, ppm_error = 5, verbose = FALSE)
  expect_identical(rownames(result), "a")
  expect_equal(SummarizedExperiment::rowData(result)$mz, index$expected_mz[1])
  expect_equal(SingleCellExperiment::altExp(result, "transcriptome"), originalTranscriptome)
  expect_true(is.list(S4Vectors::metadata(result)$mz_annotation))
  result <- annotateSM(object, index = index, ppm_error = 5,
                        return.only.annotated = FALSE, verbose = FALSE)
  expect_identical(rownames(result), rownames(object))
  expect_true(is.numeric(SummarizedExperiment::rowData(result)$mz))
  expect_equal(SummarizedExperiment::rowData(result)$all_IsomerNames[2:3],
               rep("No Annotation", 2))
})

test_that("curated FMP10 annotations retain panel fields and arbitrary feature IDs", {
  object <- nativeFixture()
  panel <- data.frame(mass = 100, annotation = "curated", Adduct = "FMP10",
    Formula = "C5H8O5", Isomers = "id1", Isomers_IDs = "hmdb:ID1",
    IsomerNames = "curated")
  database <- list(filtered_fmp10 = panel,
    chem_props = data.frame(ramp_id = "RAMP_C_1", chem_source_id = "hmdb:ID1"))
  result <- addFMP10Annotations(object, database = database,
    return.only.annotated = TRUE)
  expect_identical(rownames(result), "a")
  expect_identical(SummarizedExperiment::rowData(result)$all_Ramp_IDs, "RAMP_C_1")
  expect_equal(SummarizedExperiment::rowData(result)$mz, 100)
  expect_equal(S4Vectors::metadata(result)$db_3$observed_mz, 100)
})
