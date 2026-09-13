nativeRecipe <- function() {
  skip_if_not_installed("SpaMTPData", "0.99.2")
  skip_if_not_installed("SeuratObject", "5.0.0")
  path <- system.file("scripts", "native-resource-utils.R", package = "SpaMTPData")
  skip_if(!nzchar(path), "The native data preparation recipe is unavailable")
  environment <- new.env(parent = globalenv())
  sys.source(path, envir = environment)
  environment
}

nativeRecipeSpec <- function(alternatives = "transcriptome") {
  data.frame(resource = "fixture", primary_assay = "Spatial",
    alternative_assays = alternatives, image = "fov", modality = "metabolome",
    synthetic = TRUE, stringsAsFactors = FALSE)
}

test_that("one-time preparation keeps analysis data and omits legacy geometry", {
  recipe <- nativeRecipe()
  seurat <- spatialExperimentToSeurat(nativeFixture())
  data <- log1p(SeuratObject::LayerData(seurat, assay = "Spatial", layer = "counts"))
  SeuratObject::LayerData(seurat, assay = "Spatial", layer = "data") <- data
  SeuratObject::Misc(seurat, slot = "polygon_coordinates") <- list(p1 = matrix(1:8, 4))
  SeuratObject::Misc(seurat, slot = "db_3") <- data.frame(mz = 100, annotation = "kept")
  metadata <- data.frame(value = seq_len(ncol(seurat)), row.names = colnames(seurat))
  names(metadata) <- "2-Hydroxyglutarate intensity"
  seurat <- SeuratObject::AddMetaData(seurat, metadata)
  specification <- nativeRecipeSpec()
  object <- recipe$prepareNativeResource(seurat, specification)
  expect_s4_class(object, "SpatialExperiment")
  expect_true(recipe$validateNativeResource(object, seurat, specification))
  expect_identical(SummarizedExperiment::assay(object, "data"), data)
  expect_identical(SummarizedExperiment::colData(object)[[names(metadata)]],
    seq_len(ncol(seurat)))
  expect_null(S4Vectors::metadata(object)$polygon_coordinates)
  expect_identical(S4Vectors::metadata(object)$db_3$annotation, "kept")
  expect_identical(SingleCellExperiment::altExpNames(object), "transcriptome")
  expect_equal(ncol(SummarizedExperiment::colData(
    SingleCellExperiment::altExp(object, "transcriptome"))), 0L)
  expect_false("scaled" %in% SummarizedExperiment::assayNames(object))
  expect_true("polygon_coordinates" %in%
    S4Vectors::metadata(object)$native_conversion$omitted_metadata)
})

test_that("native preparation rejects incomplete layers and unpaired modalities", {
  recipe <- nativeRecipe()
  seurat <- spatialExperimentToSeurat(nativeFixture())
  counts <- SeuratObject::LayerData(seurat, assay = "Spatial", layer = "counts")
  SeuratObject::LayerData(seurat, assay = "Spatial", layer = "data") <-
    counts[1:2, , drop = FALSE]
  expect_error(recipe$prepareNativeResource(seurat, nativeRecipeSpec()),
    "complete counts matrix")
  seurat <- spatialExperimentToSeurat(nativeFixture())
  alternative <- SeuratObject::LayerData(seurat, assay = "transcriptome", layer = "counts")
  expect_warning(seurat[["transcriptome"]] <- SeuratObject::CreateAssay5Object(
    counts = alternative[, -1, drop = FALSE]), "Different cells")
  expect_error(recipe$prepareNativeResource(seurat, nativeRecipeSpec()),
    "not paired")
  expect_error(recipe$nativePortableMetadata(new.env()), "Non-portable")
})
