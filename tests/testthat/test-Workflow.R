workflowFixture <- function() {
  withr::local_seed(11)
  ids <- paste0("obs", 1:18)
  E <- matrix(stats::rpois(8 * 18, 8), 8,
    dimnames = list(paste0("mass", 1:8), ids))
  x <- SpatialExperiment::SpatialExperiment(assays = list(counts = E),
    colData = S4Vectors::DataFrame(group = rep(c("A", "B"), 9),
      patient = rep(paste0("p", 1:9), each = 2), row.names = ids),
    spatialCoords = cbind(x = rep(1:6, each = 3), y = rep(1:3, 6)))
  G <- matrix(stats::rpois(6 * 18, 5), 6, dimnames = list(paste0("gene", 1:6), ids))
  SingleCellExperiment::altExp(x, "RNA") <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = G))
  P <- matrix(stats::rpois(3 * 18, 10), 3, dimnames = list(paste0("protein", 1:3), ids))
  SingleCellExperiment::altExp(x, "protein") <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = P))
  x
}

workflowSpecs <- function() list(main = list(type = "metabolomics", layer = "counts"),
  RNA = list(type = "transcriptomics", layer = "counts"),
  protein = list(type = "proteomics", layer = "counts"))

test_that("three-omic workflow preserves raw layers and integrates the selected normalization", {
  x <- workflowFixture(); withr::local_seed(27); previous <- .Random.seed
  r <- runSpaMTPWorkflow(x, workflowSpecs(), npcs = 2, clusters = 2, group = "group")
  expect_identical(.Random.seed, previous)
  expect_s3_class(r, "spamtp_workflow")
  expect_equal(SummarizedExperiment::assay(r$object, "counts"), SummarizedExperiment::assay(x, "counts"))
  E <- SummarizedExperiment::assay(x, "counts")
  expected <- log2(1 + sweep(E, 2, colSums(E), "/") * 10000)
  expect_equal(as.matrix(SummarizedExperiment::assay(r$object, "workflow")), expected)
  for (n in names(workflowSpecs())) {
    e <- r$analysis[[n]]$pca$scores
    expect_equal(unname(r$joint$embedding[, startsWith(colnames(r$joint$embedding), paste0(n, "_")), drop = FALSE]),
      unname(scale(e) / sqrt(ncol(e))), ignore_attr = TRUE)
  }
  expect_length(r$associations, 3)
  expect_length(r$analysis$RNA$comparisons$tests, 0)
  expect_equal(ncol(r$joint$display), 2)
  r2 <- runSpaMTPWorkflow(x, workflowSpecs(), npcs = 2, clusters = 2, group = "group")
  expect_equal(r$joint, r2$joint)
})

test_that("QC excludes joint observations but keeps measured constant features", {
  x <- workflowFixture()
  e <- SingleCellExperiment::altExp(x, "RNA")
  SummarizedExperiment::assay(e)[, 1] <- 0
  SummarizedExperiment::assay(e)[1, ] <- 0
  SingleCellExperiment::altExp(x, "RNA") <- e
  r <- runSpaMTPWorkflow(x, workflowSpecs(), npcs = 2)
  expect_equal(ncol(r$object), 17)
  expect_false(r$retained$included[1])
  expect_equal(nrow(SingleCellExperiment::altExp(r$object, "RNA")), 6)
  expect_false("gene1" %in% r$analysis$RNA$pca$features)
  expect_equal(ncol(x), 18)
})

test_that("independent modalities align to the chosen reference and audit unmatched locations", {
  x <- workflowFixture(); source <- x[, -18]
  transform <- diag(3); transform[1:2, 3] <- c(-10, 3)
  xy <- SpatialExperiment::spatialCoords(source)
  SpatialExperiment::spatialCoords(source) <- sweep(xy, 2, c(10, -3), "+")
  specs <- list(MS = list(type = "metabolomics", layer = "counts"),
    RNA = list(type = "transcriptomics", layer = "counts"),
    protein = list(type = "proteomics", layer = "counts"))
  r <- runSpaMTPWorkflow(list(MS = source, RNA = x, protein = x), specs,
    reference = "RNA", coordinate_units = "micrometres", npcs = 2,
    alignment = list(MS = list(method = "affine", alignment = transform), protein = list(method = "identity")),
    mapping = list(MS = list(method = "nearest", max_distance = 0.1),
      protein = list(method = "pixel", width = 0.5)))
  expect_equal(ncol(r$object), 17)
  expect_setequal(SingleCellExperiment::altExpNames(r$object), c("MS", "protein"))
  expect_equal(as.matrix(r$mapping$MS$weights), cbind(diag(17), 0), ignore_attr = TRUE)
  expect_false(tail(r$mapping$MS$observations$mapped, 1))
  expect_equal(as.matrix(SummarizedExperiment::assay(SingleCellExperiment::altExp(r$object, "MS"), "counts")),
    as.matrix(SummarizedExperiment::assay(source, "counts")), ignore_attr = TRUE)
})

test_that("nearest correspondence cannot cross samples or exceed the radius", {
  x <- workflowFixture(); y <- x
  y$sample_id <- rep("other", ncol(y))
  expect_error(runSpaMTPWorkflow(list(MS = x, RNA = y),
    list(MS = workflowSpecs()$main, RNA = workflowSpecs()$RNA),
    reference = "RNA", coordinate_units = "um", alignment = list(MS = list(method = "identity")),
    mapping = list(MS = list(method = "nearest", max_distance = 1))), "Fewer than three")
})

test_that("replicate contrasts match an independent paired limma design", {
  x <- workflowFixture()
  r <- runSpaMTPWorkflow(x, workflowSpecs(), npcs = 2, group = "group", replicate = "patient",
    contrasts = list(B_A = list(numerator = "B", denominator = "A", paired = TRUE)))
  t <- r$analysis$RNA$comparisons$tests$B_A
  expect_identical(t$status, "completed")
  E <- as.matrix(SummarizedExperiment::assay(SingleCellExperiment::altExp(r$object, "RNA"), "workflow"))
  design <- stats::model.matrix(~ factor(x$patient) + factor(x$group))
  f <- limma::eBayes(limma::lmFit(E, design), trend = TRUE)
  tab <- limma::topTable(f, coef = ncol(design), number = Inf, sort.by = "none")
  expect_equal(t$table$logFC, tab$logFC, tolerance = 1e-10)
  expect_equal(t$table$P.Value, tab$P.Value, tolerance = 1e-10)
  expect_error(runSpaMTPWorkflow(x, workflowSpecs(), group = "group", replicate = "patient",
    contrasts = list(B_A = list(numerator = "B", denominator = "A"))), "paired=TRUE")
  x$patient <- "one"
  skip <- runSpaMTPWorkflow(x, workflowSpecs(), group = "group", replicate = "patient",
    contrasts = list(B_A = list(numerator = "B", denominator = "A", paired = TRUE)))
  expect_identical(skip$analysis$RNA$comparisons$tests$B_A$status, "skipped")
})

test_that("reports embed figures, escape labels and preserve existing output", {
  x <- workflowFixture(); directory <- tempfile()
  on.exit(unlink(directory, recursive = TRUE))
  r <- runSpaMTPWorkflow(x, workflowSpecs()[1:2], output_dir = directory, npcs = 2,
    title = "<script>alert('bad')</script>")
  html <- paste(readLines(r$report), collapse = "\n")
  expect_match(html, "data:image/png;base64,")
  expect_false(grepl("<script>", html, fixed = TRUE))
  expect_match(html, "&lt;script&gt;", fixed = TRUE)
  saved <- readRDS(file.path(directory, "workflow.rds"))
  expect_identical(saved$report, r$report)
  manifest <- utils::read.csv(file.path(directory, "manifest.csv"))
  expect_equal(unname(tools::md5sum(file.path(directory, manifest$file))), manifest$md5)
  before <- tools::md5sum(r$report)
  expect_error(renderSpaMTPReport(r, directory), "new or empty")
  expect_equal(tools::md5sum(r$report), before)
})

test_that("bad specifications fail instead of selecting an unintended layer", {
  x <- workflowFixture(); s <- workflowSpecs()
  s$main$layer <- "absent"
  expect_error(runSpaMTPWorkflow(x, s), "absent")
  s <- workflowSpecs(); s$main$normalisation <- "none"
  expect_error(runSpaMTPWorkflow(x, s), "Unknown modality")
  expect_error(runSpaMTPWorkflow(x, list(main = list(type = "metabolomics"))), "type and layer")
  expect_error(runSpaMTPWorkflow(list(resource = "x", version = "latest"), workflowSpecs()), "Pin")
  expect_error(runSpaMTPWorkflow(x, workflowSpecs(), contrasts = list(test = list())), "group")
})

test_that("normalized non-spatial inputs are retained and pathway members share one index", {
  x <- workflowFixture()
  SpatialExperiment::spatialCoords(x) <- matrix(numeric(), ncol(x), 0)
  rna <- SingleCellExperiment::altExp(x, "RNA")
  rownames(rna) <- paste0("RAMP_G_", 1:6)
  SummarizedExperiment::assay(rna, "data") <- log2(1 + SummarizedExperiment::assay(rna)) - 2
  SingleCellExperiment::altExp(x, "RNA") <- rna
  db <- list(analytehaspathway = data.frame(rampId = paste0("RAMP_G_", c(1, 2, 3, 4, 5, 6)),
    pathwayRampId = rep(c("P1", "P2", "P3"), each = 2)),
    pathway = data.frame(pathwayRampId = c("P1", "P2", "P3"), pathwayName = c("one", "two", "three")))
  index <- buildPathwayIndex(db, gene_mapping = "ramp")
  specs <- workflowSpecs()[1:2]
  specs$RNA <- list(type = "transcriptomics", layer = "data", normalization = "none", species = "custom",
    pathways = list(index = index, min_size = 1, foreground = "RAMP_G_1"))
  r <- runSpaMTPWorkflow(x, specs, npcs = 2)
  expect_equal(as.matrix(SummarizedExperiment::assay(SingleCellExperiment::altExp(r$object, "RNA"), "workflow")),
    as.matrix(SummarizedExperiment::assay(rna, "data")))
  p <- r$analysis$RNA$pathways
  expect_equal(p$coverage$measured_size, rep(2L, 3))
  expect_equal(p$coverage$used_members, attr(p$enrichment, "pathway_coverage")$measured_members)
  expect_equal(ncol(SpatialExperiment::spatialCoords(r$object)), 0)
  expect_length(p$enrichment[[1]], 3)
  expect_identical(r$resources$RNA$pathway_index, index)
  specs$RNA$pathways$min_size <- 7
  skipped <- runSpaMTPWorkflow(x, specs, npcs = 2)
  expect_identical(skipped$analysis$RNA$pathways$status, "skipped")
  expect_equal(nrow(skipped$analysis$RNA$pathways$coverage), 3L)
})

test_that("resource acquisition is pinned and RDS inputs record content identity", {
  path <- tempfile(fileext = ".rds"); on.exit(unlink(path))
  saveRDS(workflowFixture(), path)
  r <- runSpaMTPWorkflow(path, workflowSpecs(), npcs = 2)
  expect_identical(r$origins$paired$md5, unname(tools::md5sum(path)))
  expect_identical(r$origins$paired$kind, "RDS")
})
