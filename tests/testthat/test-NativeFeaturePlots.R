test_that("m/z plots sum full windows and read the requested experiment", {
  object <- nativeFixture()
  SummarizedExperiment::rowData(object)$mz <- c(100, 100.1, 100.2)
  values <- .mzPlotValues(object, 100.1, "main", "counts", 0.1)
  expect_equal(as.numeric(values),
    as.numeric(Matrix::colSums(SummarizedExperiment::assay(object, "counts"))))
  alternative <- SingleCellExperiment::altExp(object, "transcriptome")
  SummarizedExperiment::rowData(alternative)$mz <- c(500, 600)
  SingleCellExperiment::altExp(object, "otherMSI") <- alternative
  values <- .mzPlotValues(object, 500, "otherMSI", "counts")
  expect_equal(as.numeric(values), as.numeric(SummarizedExperiment::assay(alternative)[1, ]))
  expect_error(imageMZPlot(object, 100, slot = "missing"), "was not found")
  expect_error(imageMZPlot(object, 100, blend = TRUE), "unsupported")
  expect_error(imageMZPlot(object, 100, plusminus = -1), "non-negative")
  expect_error(imageMZPlot(object, "unknown"), "Unknown")
})

test_that("annotation panels aggregate independently without double-counting peaks", {
  object <- nativeFixture()
  SummarizedExperiment::rowData(object)$mz <- c(100, 100.1, 100.2)
  SummarizedExperiment::rowData(object)$all_IsomerNames <-
    c("Glucose; Glutamine", "Glucose", "Glutamine")
  before <- object
  plots <- imageMZAnnotationPlot(object, c("Glucose", "Glutamine"), combine = FALSE)
  counts <- SummarizedExperiment::assay(object, "counts")
  expect_length(plots, 2)
  expect_equal(plots[[1]]$data$score, as.numeric(Matrix::colSums(counts[1:2, ])))
  expect_equal(plots[[2]]$data$score, as.numeric(Matrix::colSums(counts[c(1, 3), ])))
  values <- .mzPlotValues(object, NULL, "main", "counts", 0.1, "Glucose")
  expect_equal(as.numeric(values), as.numeric(Matrix::colSums(counts)))
  expect_error(imageMZAnnotationPlot(object, "Gluc"), "No matching")
  expect_s3_class(imageMZAnnotationPlot(object, "Gluc", plot.exact = FALSE), "ggplot")
  expect_s3_class(spatialMZAnnotationPlot(object, "Glucose"), "ggplot")
  expect_identical(object, before)
})

test_that("current annotation stores can label arbitrary native feature IDs", {
  object <- nativeFixture()
  SummarizedExperiment::rowData(object)$annotation <- NULL
  annotations <- data.frame(observed_mz = 100, Adduct = "M+H",
    Ramp_IDs = "RAMP_C_1", IsomerNames = "Glucose", Score = 1,
    MassScore = 1, ChemicalScore = 1, IsotopeScore = NA_real_,
    AdductNetworkScore = NA_real_)
  object <- .setStoredData(object, "mz_annotation",
    list(results = annotations, metadata = list(schema_version = 2L)))
  result <- searchAnnotations(object, "Glucose", search.exact = TRUE)
  expect_identical(result$mz_names, "a")
  expect_identical(rownames(result), "a")
  expect_equal(imageMZAnnotationPlot(object, "Glucose")$data$score, as.numeric(1:6))
})

test_that("sample, split and cutoff controls keep coordinates aligned", {
  object <- nativeFixture(rep(c("s1", "s2"), 3))[, c(6, 1, 4, 3, 2, 5)]
  plots <- imageMZPlot(object, c(100, 200), combine = FALSE, scale = "all")
  expect_length(plots, 4)
  for (plot in plots) {
    expect_length(unique(plot$data$sample_id), 1L)
    expect_equal(plot$data$x, SpatialExperiment::spatialCoords(object)[
      match(plot$data$cell, colnames(object)), "x"],
      ignore_attr = TRUE)
    expect_equal(plot$scales$get_scales("fill")$limits, c(1, 6))
  }
  selected <- imageMZPlot(object, 100, fov = "s1", cells = c("p1", "p3"))
  expect_setequal(selected$data$cell, c("p1", "p3"))
  expect_setequal(selected$data$score, c(1, 3))
  split <- imageMZPlot(nativeFixture(), 100, split.by = "region", combine = FALSE)
  expect_length(split, 2)
  clipped <- imageMZPlot(nativeFixture(), 100, min.cutoff = "q50", max.cutoff = 5)
  expect_equal(clipped$data$score, c(3.5, 3.5, 3.5, 4, 5, 5))
  expect_error(imageMZPlot(object, 100, min.cutoff = 8, max.cutoff = 1), "exceeds")
  expect_error(imageMZPlot(object, 100, min.cutoff = c(1, 2)), "one cutoff")
  expect_error(imageMZPlot(object, 100, min.cutoff = "q101"), "Cutoffs")
  expect_error(imageMZPlot(object, 100, max.cutoff = "typo"), "Cutoffs")
})

test_that("pathway-network payloads align native RNA and MSI with spatial pixels", {
  object <- nativeFixture()[, c(6, 2, 4, 1, 5, 3)]
  genes <- .pn_layer_data(object, "transcriptome", "counts")
  metabolites <- .pn_layer_data(object, "main", "counts")
  payload <- .pn_spatial_payload(object, "region", NULL,
    c("RAMP_G_1", "RAMP_C_1"),
    gene_de = data.frame(rampId = "RAMP_G_1", gene = "a", p_val_adj = 0.01),
    metabolite_de = data.frame(ramp_id = "RAMP_C_1", mz_name = "a", p_val_adj = 0.01),
    gene_matrix = genes, metabolite_matrix = metabolites)
  expect_identical(payload$clusters, as.character(object$region))
  expect_equal(payload$expression$rna$RAMP_G_1, as.numeric(genes["a", ]))
  expect_equal(payload$expression$mets$RAMP_C_1, as.numeric(metabolites["a", ]))
  expect_equal(payload$coordinates$x, c(1, 0.2, 0.6, 0, 0.8, 0.4))
  expect_error(.pn_layer_data(object, "transcriptome", "missing"), "was not found")
})

test_that("optical overlays use imgData scale and native interactive widgets", {
  object <- nativeFixture()
  path <- tempfile(fileext = ".png")
  on.exit(unlink(path))
  raster <- array(0, c(2, 3, 3))
  raster[, , 1] <- matrix(c(1, 0, 0, 0, 0, 0), 2)
  raster[2, 1, 2] <- 1
  raster[1, 3, 3] <- 1
  png::writePNG(raster, path)
  object <- addSpatialImage(object, path, imageId = "optical", scaleFactor = 2)
  image <- .nativeImage(object, "optical", unique(object$sample_id))
  faded <- .imageLayer(image, alpha = 0.4)$geom_params$raster
  expect_equal(dim(faded), dim(image$raster))
  expect_equal(grDevices::col2rgb(as.vector(faded)),
    grDevices::col2rgb(as.vector(image$raster)))
  expect_equal(as.numeric(grDevices::col2rgb(as.vector(faded), alpha = TRUE)[4, ]),
    rep(102, 6))
  plot <- spatialMZPlot(object, 100, images = "optical", crop = FALSE)
  expect_s3_class(plot, "ggplot")
  expect_equal(plot$layers[[1]]$geom_params$xmax, 1.5)
  expect_equal(plot$layers[[1]]$geom_params$ymax, 1)
  expect_equal(nrow(ggplot2::ggplot_build(plot)$data[[3]]), ncol(object))
  expect_error(spatialMZPlot(object, 100, images = "absent"), "not found")
  expect_s3_class(spatialMZPlot(object, 100, interactive = TRUE), "plotly")
  plots <- spatialMZPlot(object, c(100, 200), interactive = TRUE, combine = FALSE)
  expect_length(plots, 2)
  expect_true(all(vapply(plots, inherits, logical(1), "plotly")))
  widget <- plotly::plotly_build(plot3DFeature(object, "a", show.image = "optical"))
  traces <- widget$x$data
  expect_length(traces, 3)
  expect_equal(as.numeric(traces[[1]]$marker$color), as.numeric(6:1))
  expect_equal(as.numeric(traces[[2]]$marker$color), as.numeric(1:6))
  expect_equal(as.numeric(traces[[3]]$x), rep(c(0.25, 0.75, 1.25), 2))
  expect_equal(as.numeric(traces[[3]]$y), rep(c(0.25, 0.75), each = 3))
  expect_equal(grDevices::col2rgb(traces[[3]]$marker$color),
    grDevices::col2rgb(c("red", "black", "blue", "green", "black", "black")))
  multiple <- nativeFixture(rep(c("s1", "s2"), 3))
  expect_error(plot3DFeature(multiple, "a"), "Choose sampleId")
  expect_s3_class(plot3DFeature(multiple, "a", sampleId = "s2"), "plotly")
})

test_that("interactive mass windows use masses and the requested altExp", {
  object <- nativeFixture()
  alternative <- SingleCellExperiment::altExp(object, "transcriptome")
  SummarizedExperiment::rowData(alternative)$mz <- c(500, 501)
  SingleCellExperiment::altExp(object, "otherMSI") <- alternative
  app <- interactiveSpatialPlot(object, assay = "otherMSI")
  expect_s3_class(app, "shiny.appobj")
  suppressPackageStartupMessages(suppressWarnings(library(shiny)))
  shiny::testServer(app, {
    session$setInputs(curve_center = 500, curve_width = 0, bin_mode = "mz",
      spot_size = 1, x_min = 499, x_max = 502)
    expect_identical(selectedRows(), 1L)
    expect_equal(selectedValues(), as.numeric(6:1))
    session$setInputs(curve_width = 1)
    expect_identical(selectedRows(), 1:2)
    expect_equal(selectedValues(), as.numeric(6:1 + c(0, 1, 0, 1, 0, 1)))
    session$setInputs(curve_center = 700, curve_width = 0)
    expect_length(selectedRows(), 0)
  })
  expect_equal(.mzWindow(500, 20, "ppm"), c(499.99, 500.01))
  expect_error(.mzWindow(500, -1, "mz"), "non-negative")
})

test_that("density HTML exports native masses, empty peaks and safe annotations", {
  object <- nativeFixture()
  SummarizedExperiment::rowData(object)$annotation <- c("</script>", "b", "c")
  directory <- tempfile()
  on.exit(unlink(directory, recursive = TRUE))
  output <- densityMap(object, folder = directory)
  expect_true(file.exists(output))
  html <- paste(readLines(output), collapse = "\n")
  expect_true(grepl("\\u003c", html, fixed = TRUE))
  expect_match(html, '"mz":100')
  SummarizedExperiment::assay(object, "counts")[] <- 0
  expect_no_error(densityMap(object, folder = directory))
  multiple <- nativeFixture(rep(c("s1", "s2"), 3))
  expect_error(densityMap(multiple, folder = directory), "Choose sampleId")
})

test_that("mass spectra and pixel symbols do not infer masses or swap axes", {
  object <- nativeFixture()
  spectrum <- massIntensityPlot(object, group.by = "region", mz.labels = 100)
  expect_setequal(spectrum$data$mass, c(100, 200, 300))
  expected <- Matrix::rowMeans(SummarizedExperiment::assay(object, "counts")[, c(1, 3, 5)])
  expect_equal(spectrum$data$intensity[spectrum$data$variable == "edge"], as.numeric(expected))
  expect_no_error(ggplot2::ggplot_build(spectrum))
  plain <- massIntensityPlot(object)
  expect_equal(plain$data$intensity,
    as.numeric(Matrix::rowMeans(SummarizedExperiment::assay(object, "counts"))))
  expect_no_error(ggplot2::ggplot_build(massIntensityPlot(object, split.by = "region")))
  original <- imageMZPlot(object, 100)
  squared <- pixelPlot(original)
  expect_identical(squared$data, original$data)
  expect_equal(squared$layers[[1]]$aes_params$shape, 22)
  expect_equal(original$layers[[1]]$aes_params$shape, 21)
  expect_length(pixelPlot(list(original, original)), 2)
  binned <- binMetabolites(object, c("a", "a", "b"))
  expect_equal(binned$Binned_Metabolites,
    as.numeric(Matrix::colSums(SummarizedExperiment::assay(object, "counts")[1:2, ])))
  expect_error(binMetabolites(object, "a", bin_name = "sample_id"), "new colData")
})
