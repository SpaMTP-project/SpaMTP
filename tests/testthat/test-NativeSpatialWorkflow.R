test_that("metadata updates retain sample information and align by ID", {
  object <- nativeFixture()
  original <- SummarizedExperiment::colData(object)
  updates <- data.frame(score = 6:1, row.names = rev(colnames(object)))
  result <- .setCellMetadata(object, updates)
  expect_identical(result$sample_id, original$sample_id)
  expect_identical(result$region, original$region)
  expect_equal(result$score, 1:6)
  rownames(updates)[1] <- "unknown"
  expect_error(.setCellMetadata(object, updates), "row names must match")
  expect_identical(SummarizedExperiment::colData(object), original)
})

test_that("alternative experiments have independent feature and assay access", {
  object <- nativeFixture()
  expect_equal(as.numeric(.assayData(object, "transcriptome", "counts")["a", ]), 6:1)
  expect_equal(as.numeric(.assayData(object, "main", "counts")["a", ]), 1:6)
  expect_identical(rownames(.featureMetadata(object, "transcriptome")), c("a", "gene2"))
  updated <- .setAssayData(object, .assayData(object, "transcriptome") * 2,
                            "transcriptome")
  expect_equal(.assayData(updated, "main"), .assayData(object, "main"))
  expect_equal(.assayData(updated, "transcriptome"),
                 .assayData(object, "transcriptome") * 2)
  expect_error(.assayData(object, "missing"), "Unknown assay or modality")
})

test_that("correlation uses the requested modality and keeps undefined values", {
  object <- nativeFixture()
  result <- findCorrelatedFeatures(object, gene = "a", ST.assay = "transcriptome",
                                    nfeatures = NULL)
  expect_equal(result$correlation[result$features == "a" & result$modality == "gene"], 1)
  expect_equal(result$correlation[result$features == "a" & result$modality == "metabolite"], -1)
  expect_true(is.na(result$correlation[result$features == "constant"]))
  expect_true(all(is.na(result$mz[result$modality == "gene"])))
  groups <- findCorrelatedFeatures(object, ident = "region", nfeatures = NULL)
  expect_setequal(groups$ident, c("edge", "core"))
  expect_equal(groups$correlation[groups$features == "a" & groups$ident == "edge"],
                 stats::cor(1:6, as.numeric(object$region == "edge")))
  expect_error(findCorrelatedFeatures(object, mz = 100, gene = "a"), "exactly one")
})

test_that("Moran statistics match the equation and do not rank constant rows", {
  withr::local_seed(123)
  object <- nativeFixture()
  result <- findSpatiallyVariableMetabolites(object, max_spots = NULL, verbose = FALSE)
  values <- as.numeric(.assayData(object)["a", ])
  weights <- 1 / as.matrix(stats::dist(SpatialExperiment::spatialCoords(object)))^2
  diag(weights) <- 0
  weights <- weights / rowSums(weights)
  centered <- values - mean(values)
  expected <- length(values) / sum(weights) *
    as.numeric(crossprod(centered, weights %*% centered)) / sum(centered^2)
  metadata <- SummarizedExperiment::rowData(result)
  expect_equal(metadata$MoransI_observed[1], expected, tolerance = 1e-12)
  expect_equal(metadata$MoransI_p.value[1], ape::Moran.I(values, weights)$p.value)
  expect_true(is.na(metadata$MoransI_observed[3]))
  expect_false(metadata$moransi.spatially.variable[3])
  expect_false("constant" %in% getSpatiallyVariableMetabolites(result))
  expect_identical(metadata$annotation, SummarizedExperiment::rowData(object)$annotation)
  expect_identical(result$region, object$region)
  before <- .Random.seed
  sampled <- findSpatiallyVariableMetabolites(object, max_spots = 4, seed = 17,
                                              verbose = FALSE)
  expect_identical(.Random.seed, before)
  repeatResult <- findSpatiallyVariableMetabolites(object, max_spots = 4, seed = 17,
                                                   verbose = FALSE)
  expect_identical(S4Vectors::metadata(sampled)$moransi$pixels,
                     S4Vectors::metadata(repeatResult)$moransi$pixels)
  object <- nativeFixture(rep(c("s1", "s2"), each = 3))
  expect_error(findSpatiallyVariableMetabolites(object), "Choose sampleId")
})

test_that("pathway scores preserve metadata and plots work with SPE and altExp", {
  object <- nativeFixture()
  original <- SummarizedExperiment::colData(object)
  scored <- addGesecaScores(list(first = "a", second = "b"), object)
  expect_identical(scored$sample_id, object$sample_id)
  expect_identical(scored$region, object$region)
  expect_true(all(c("first", "second") %in% colnames(SummarizedExperiment::colData(scored))))
  expect_identical(SummarizedExperiment::colData(object), original)
  expect_error(addGesecaScores(list(absent = "unknown"), object), "no features")
  SingleCellExperiment::reducedDim(object, "PCA") <- cbind(1:6, c(1, 3, 2, 5, 4, 6))
  plot <- plotSinglePathway("a", object, assay = "transcriptome", slot = "counts")
  expect_s3_class(plot, "ggplot")
  expect_true(length(ggplot2::ggplot_build(plot)$data) > 0L)
  spatial <- plotSinglePathwaySpatially("a", object)
  expect_s3_class(spatial, "ggplot")
  expect_equal(nrow(ggplot2::ggplot_build(spatial)$data[[1]]), ncol(object))
})

test_that("mapping respects tissue identity and preserves target and source data", {
  source <- nativeFixture(rep(c("s1", "s2"), 3))[, c(1, 2), drop = FALSE]
  target <- nativeFixture(rep(c("s2", "s1", "s1"), 2))[, c(1, 2, 3), drop = FALSE]
  SpatialExperiment::spatialCoords(source) <- cbind(x = c(0, 0), y = c(0, 0))
  SpatialExperiment::spatialCoords(target) <- cbind(x = c(0, 0, 10), y = c(0, 0, 10))
  result <- mapSpatialOmics(source, target, ST.hires = TRUE, SM.pixel.width = 1,
                             map.data = TRUE)
  expect_equal(as.numeric(SummarizedExperiment::assay(result, "counts")["a", ]),
                 c(2, 1, 0))
  expect_identical(result$mapped, c(TRUE, TRUE, FALSE))
  expect_identical(result$region, target$region)
  expect_identical(colnames(result), colnames(target))
  expect_equal(SummarizedExperiment::assay(result, "logcounts")[, 1],
                 SummarizedExperiment::assay(source, "logcounts")[, 2], ignore_attr = TRUE)
  expect_equal(SummarizedExperiment::assay(
    SingleCellExperiment::altExp(result, "transcriptome"), "counts"),
    SummarizedExperiment::assay(target, "counts"))
  expect_identical(SummarizedExperiment::rowData(result)$annotation,
                     SummarizedExperiment::rowData(source)$annotation)
  expect_equal(ncol(mapSpatialOmics(source, target, ST.hires = TRUE,
    SM.pixel.width = 1, dropUnmapped = TRUE)), 2L)
})

test_that("area mapping averages intersecting MSI pixels and validates geometry", {
  source <- nativeFixture()[, 1:2, drop = FALSE]
  target <- nativeFixture()[, 1, drop = FALSE]
  SpatialExperiment::spatialCoords(source) <- cbind(x = c(-0.5, 0.5), y = 0)
  SpatialExperiment::spatialCoords(target) <- cbind(x = 0, y = 0)
  polygon <- sf::st_sf(cell = colnames(target), geometry = sf::st_sfc(
    sf::st_polygon(list(rbind(c(-1, -0.5), c(1, -0.5), c(1, 0.5),
                             c(-1, 0.5), c(-1, -0.5))))))
  result <- mapSpatialOmics(source, target, SM.pixel.width = 1,
                             ST.polygons = polygon, overlap.threshold = 0.49)
  expect_equal(as.numeric(SummarizedExperiment::assay(result)["a", ]), 1.5)
  expect_equal(result$n_source_pixels, 2L)
  expect_error(mapSpatialOmics(source, target, SM.pixel.width = 1), "ST.radius")
})

test_that("affine entry point transforms coordinates without changing assay data", {
  object <- nativeFixture()
  transform <- diag(3)
  transform[1:2, 3] <- c(10, -4)
  result <- alignSpatialOmics(object, object, transformation = transform)
  expected <- SpatialExperiment::spatialCoords(object)
  expected[, 1] <- expected[, 1] + 10
  expected[, 2] <- expected[, 2] - 4
  expect_equal(SpatialExperiment::spatialCoords(result), expected)
  expect_equal(SummarizedExperiment::assays(result), SummarizedExperiment::assays(object))
  expect_identical(result$sample_id, object$sample_id)
  expect_true("spatial_alignment" %in% names(S4Vectors::metadata(result)))
  expect_s3_class(checkAlignment(result, object), "ggplot")
})

test_that("image overlays honour the stored image scale factor", {
  object <- nativeFixture()
  path <- tempfile(fileext = ".png")
  on.exit(unlink(path))
  png::writePNG(array(1, c(4, 8, 3)), path)
  object <- addSMImage(path, object, scaleFactor = 2, interactive = FALSE)
  image <- .nativeImage(object, "optical", "section1")
  expect_equal(c(image$width, image$height), c(4, 2))
  plot <- plotSinglePathwaySpatially("a", object, images = "optical", crop = FALSE)
  expect_s3_class(plot, "ggplot")
  expect_true(length(ggplot2::ggplot_build(plot)$data) >= 2L)
  expect_equal(ncol(object), 6L)
})
