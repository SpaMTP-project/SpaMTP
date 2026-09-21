regionWorkflowFixture <- function() {
  withr::local_seed(83)
  E <- matrix(stats::rpois(7 * 16, 5), 7,
    dimnames = list(paste0("feature_", 1:7), paste0("observation_", 1:16)))
  x <- SpatialExperiment::SpatialExperiment(assays = list(counts = E),
    spatialCoords = cbind(x = rep(1:4, 4), y = rep(rep(1:2, each = 4), 2)),
    colData = S4Vectors::DataFrame(sample_id = rep(c("section_x", "section_y"), each = 8),
      condition = rep(c("control", "treated"), 8), donor = rep(paste0("donor", 1:8), each = 2),
      anatomy = rep(c("core", "rim"), each = 4, times = 2), row.names = colnames(E)))
  y <- E[1:5, ] + 2
  rownames(y) <- paste0("gene_", 1:5)
  SingleCellExperiment::altExp(x, "arbitrary_RNA_name") <-
    SingleCellExperiment::SingleCellExperiment(assays = list(counts = y))
  x
}

regionWorkflowSpecs <- function() list(main = list(type = "metabolomics", layer = "counts"),
  arbitrary_RNA_name = list(type = "transcriptomics", layer = "counts", species = "Mus musculus"))

test_that("regional reference excludes unlabelled observations and never invents P values", {
  E <- rbind(first = c(2, 4, 0, 0, 100, 999, 6, 0), second = c(-2, 0, 3, 1, 55, 88, 0, 1))
  r <- .wf_region_summary(E, c("A", "A", "B", "B", "", NA, "C", "C"), 2)
  expect_equal(r$unlabelled, 2)
  expect_equal(r$tables$A$effect, c(3 - 1.5, -1 - 1.25))
  expect_equal(r$tables$A$detected_reference, c(1/4, 3/4))
  expect_equal(r$tables$A$n_reference, c(4, 4))
  expect_equal(unname(r$means[, "A"]), c(3, -1))
  expect_false(any(c("P.Value", "FDR", "adj.P.Val") %in% names(r$tables$A)))
  empty <- .wf_region_summary(E, rep(NA_character_, ncol(E)), 2)
  expect_length(empty$tables, 0)
  expect_equal(dim(empty$means), c(2, 0))
  underpowered <- .wf_region_summary(E, c("tiny", rep("large", 7)), 2)
  expect_false(any(underpowered$groups$eligible))
})

test_that("preview sampling covers samples and regions without using only leading pixels", {
  md <- data.frame(sample_id = rep(c("s1", "s2"), each = 40), region = rep(rep(c("A", "B"), each = 20), 2))
  idx <- .wf_preview_indices(md, "region", 12)
  expect_length(idx, 12)
  expect_equal(as.integer(table(md$sample_id[idx], md$region[idx])), rep(3L, 4))
  expect_true(all(c(1, 20, 21, 40, 41, 60, 61, 80) %in% idx))
  expect_identical(.wf_preview_indices(md, "region", 100), seq_len(80))
  expect_length(.wf_preview_indices(md, "region", 3), 3)
})

test_that("native spatial results use recorded within-sample inputs and public functions", {
  x <- regionWorkflowFixture()
  r <- runSpaMTPWorkflow(x, regionWorkflowSpecs(), group = "condition", regions = "anatomy", npcs = 2,
    native = list(moran = TRUE, graph_pca = TRUE, max_points = 6, max_features = 4))
  expect_identical(r$region_analysis$field, "anatomy")
  expect_null(r$analysis$arbitrary_RNA_name$pathways)
  for (n in names(r$analysis)) for (s in c("section_x", "section_y")) {
    a <- r$native_analysis$modalities[[n]]$samples[[s]]
    cols <- match(a$pixels, colnames(r$object))
    expect_true(all(r$object$sample_id[cols] == s))
    expect_equal(nrow(a$moran), 4)
    expect_equal(a$moran$family_size, rep(4, 4))
    small <- .wf_modality_spatial(r$object, r$settings$assay_names[[n]], a$features, cols)
    expected <- findSpatiallyVariableMetabolites(small, slot = "workflow", max_spots = NULL, nfeatures = 4, verbose = FALSE)
    expect_equal(a$moran$I, SummarizedExperiment::rowData(expected)$MoransI_observed)
    expect_equal(a$moran$FDR, SummarizedExperiment::rowData(expected)$MoransI_p.adjust)
    expect_identical(rownames(a$graph_pca$scores), a$pixels)
    expect_equal(nrow(a$graph_pca$loadings), 4)
  }
  disabled <- analyzeSpaMTPRegions(r, native = list(moran = FALSE, graph_pca = FALSE))
  expect_identical(disabled$native_analysis$status, "not_requested")
  expect_identical(disabled$analysis, r$analysis)
  expect_error(analyzeSpaMTPRegions(r, native = list(TRUE)), "Unknown")
  expect_error(analyzeSpaMTPRegions(r, regions = "nonexistent"), "metadata field")
})

test_that("single-omic clustering produces browsable regions without claiming multi-omic integration", {
  x <- regionWorkflowFixture()
  E <- SummarizedExperiment::assay(x, "counts")
  E[1:3, 1:8] <- E[1:3, 1:8] + 50
  E[4:7, 9:16] <- E[4:7, 9:16] + 50
  SummarizedExperiment::assay(x, "counts") <- E
  r <- runSpaMTPWorkflow(x, regionWorkflowSpecs()[1], clusters = 2, npcs = 2,
    native = list(moran = FALSE, graph_pca = FALSE))
  expect_null(r$joint)
  expect_identical(r$clustering$input, "main / PCA")
  expect_equal(length(unique(r$object$workflow_cluster)), 2)
  expect_identical(r$region_analysis$field, "workflow_cluster")
  expect_length(r$region_analysis$modalities$main$tables, 2)
  expect_true(any(r$stages$stage == "clustering" & r$stages$status == "completed"))
  p <- .wf_explorer_payload(r, 10, 3, 4)
  expect_length(p$modalities, 1)
  expect_length(p$modalities[[1]]$modes, 2)
})

test_that("region extensions preserve replicate contrasts and are independent of study labels", {
  x <- regionWorkflowFixture()
  r <- runSpaMTPWorkflow(x, regionWorkflowSpecs(), group = "condition", replicate = "donor", regions = "anatomy",
    contrasts = list(treated_control = list(numerator = "treated", denominator = "control", paired = TRUE)),
    npcs = 2, native = list(moran = FALSE, graph_pca = FALSE))
  r$package_version <- "earlier_analysis_version"
  new <- analyzeSpaMTPRegions(r, regions = "anatomy", min_region_size = 2,
    native = list(moran = FALSE, graph_pca = FALSE))
  expect_identical(new$analysis, r$analysis)
  expect_identical(new$package_version, "earlier_analysis_version")
  expect_identical(new$region_analysis$package_version, as.character(utils::packageVersion("SpaMTP")))
  y <- x; y$anatomy <- ifelse(x$anatomy == "core", "甲 <script>", "unrelated tissue")
  renamed <- runSpaMTPWorkflow(y, regionWorkflowSpecs(), group = "condition", regions = "anatomy", npcs = 2,
    native = list(moran = FALSE, graph_pca = FALSE))
  for (n in names(r$analysis)) {
    actual <- renamed$region_analysis$modalities[[n]]$tables[["甲 <script>"]]
    expected <- r$region_analysis$modalities[[n]]$tables$core
    actual <- actual[match(expected$feature, actual$feature), ]
    fields <- setdiff(names(expected), "cluster")
    expect_equal(actual[, fields], expected[, fields], ignore_attr = TRUE)
  }
  p <- .wf_explorer_payload(r, points = 8, features = 3, rows = 4)
  expect_length(p$observations, 8)
  for (i in seq_along(p$modalities)) {
    m <- p$modalities[[i]]; n <- names(r$analysis)[i]
    E <- .assayData(r$object, r$settings$assay_names[[n]], "workflow")
    for (profile in m$profiles) expect_equal(unlist(profile$values),
      signif(as.numeric(E[profile$feature + 1L, match(unlist(p$observations), colnames(E))]), 6))
    test <- Filter(function(x) x$kind == "test", m$modes)[[1]]
    expected <- r$analysis[[n]]$comparisons$tests$treated_control$table
    ids <- match(unlist(m$features)[test$rows[, 1] + 1L], expected$feature)
    expect_equal(test$rows[, 7], expected$adj.P.Val[ids])
    expect_equal(test$rows[, 2], expected$logFC[ids])
    expect_equal(test$total, nrow(expected))
  }
  html <- .wf_explorer_html(renamed, points = 8, features = 3, rows = 4)
  expect_false(grepl("甲 <script>", html, fixed = TRUE))
  encoded <- sub('.*type="application/json">"([^\"]+)"</script>.*', '\\1', html)
  expect_identical(jsonlite::fromJSON(paste0('"', encoded, '"')), encoded)
  decoded <- jsonlite::fromJSON(rawToChar(jsonlite::base64_dec(encoded)))
  expect_true("甲 <script>" %in% decoded$regions)
  r$analysis$main$comparisons$tests$treated_control$table$adj.P.Val[1] <- 1.23456789e-12
  tiny <- .wf_explorer_html(r, points = 8, features = 3, rows = 4)
  encoded <- sub('.*type="application/json">"([^\"]+)"</script>.*', '\\1', tiny)
  wire <- jsonlite::fromJSON(rawToChar(jsonlite::base64_dec(encoded)), simplifyVector = FALSE)
  mode <- Filter(function(m) m$kind == "test", wire$modalities[[1]]$modes)[[1]]
  expect_equal(mode$rows[[1]][[7]], 1.23456789e-12, tolerance = 1e-20)
})
