nativeWorkflowFixture <- function(n = 60L) {
  withr::local_seed(831)
  E <- matrix(stats::rnorm(12 * n), 12,
    dimnames = list(paste0("feature", seq_len(12)), paste0("obs", seq_len(n))))
  region <- rep(c("A", "B", "C"), length.out = n)
  E[1, region == "A"] <- E[1, region == "A"] + 2
  x <- SpatialExperiment::SpatialExperiment(assays = list(logcounts = E),
    colData = S4Vectors::DataFrame(region = region, row.names = colnames(E)),
    spatialCoords = cbind(x = seq_len(n) %% 10, y = seq_len(n) %/% 10))
  SingleCellExperiment::altExp(x, "RNA") <- SingleCellExperiment::SingleCellExperiment(
    assays = list(logcounts = E + matrix(stats::rnorm(length(E), sd = .5), nrow(E))))
  x
}

test_that("explicit native PCA feature selection cannot include high-variance outsiders", {
  x <- nativeWorkflowFixture()
  E <- SummarizedExperiment::assay(x); E[12, ] <- E[12, ] * 1000
  SummarizedExperiment::assay(x) <- E
  fit <- runMetabolicPCA(x, slot = "logcounts", npcs = 2,
    features = rownames(x)[1:4], scale = TRUE, verbose = FALSE)
  z <- SingleCellExperiment::reducedDim(fit, "pca")
  reference <- stats::prcomp(t(E[1:4, ]), scale. = TRUE, rank. = 2)$x[, 1:2]
  expect_equal(tcrossprod(z), tcrossprod(reference), tolerance = 1e-8, ignore_attr = TRUE)
  expect_setequal(rownames(attr(z, "rotation")), rownames(x)[1:4])
})

test_that("sparse Graph PCA agrees with the independent dense solution", {
  x <- scaleSMData(nativeWorkflowFixture(), slot = "logcounts")
  fit <- runSpatialGraphPCA(x, n_components = 3, lambda = .6, platform = "ST", n_neighbors = 4, fast = FALSE, verbose = FALSE)
  A <- as.matrix(S4Vectors::metadata(fit)$spatialGraphs$SpatialKNN)
  E <- t(as.matrix(SummarizedExperiment::assay(x, "scaled")))
  X <- solve(diag(nrow(A)) + .6 * (diag(rowSums(A)) - A), E)
  V <- eigen(crossprod(E, X), symmetric = TRUE)$vectors[, 1:3]
  expected <- X %*% V
  actual <- SingleCellExperiment::reducedDim(fit, "SpatialPCA")
  expect_equal(tcrossprod(actual), tcrossprod(expected), tolerance = 1e-8, ignore_attr = TRUE)
  expect_s4_class(S4Vectors::metadata(fit)$spatialGraphs$SpatialKNN, "sparseMatrix")
})

test_that("region markers agree with empirical pairwise concordance without invented P values", {
  x <- nativeWorkflowFixture()
  m <- findAllDEMs(x, ident = "region", assay = "main", slot = "logcounts", method = "markers", spatial_blocks = 4)
  a <- m$DEMs[m$DEMs$gene == "feature1" & m$DEMs$cluster == "A", ]
  E <- SummarizedExperiment::assay(x)
  empirical <- vapply(c("B", "C"), function(g) {
    difference <- outer(E[1, x$region == "A"], E[1, x$region == g], "-")
    mean(difference > 0) + .5 * mean(difference == 0)
  }, numeric(1))
  expect_equal(a$mean_auc, mean(empirical), tolerance = 1e-12)
  expect_equal(a$min_auc, min(empirical), tolerance = 1e-12)
  expect_false(any(c("P.Value", "FDR", "p_val_adj") %in% names(m$DEMs)))
  expect_true(a$stability_fits > 0)
  heat <- demsHeatmap(m, n = 2, only.pos = TRUE, logfc.threshold = 0, order.by = "mean_auc")
  expect_equal(heat$expression, m$expression[heat$selected_features, , drop = FALSE])
})

test_that("conditional native correlation removes a pure regional composition effect", {
  region <- rep(c("A", "B"), each = 8)
  shift <- rep(c(0, 100), each = 8)
  xvalue <- shift + rep(c(-1, -1, 1, 1), 4)
  yvalue <- shift + rep(c(-1, 1, -1, 1), 4)
  E <- rbind(protein = xvalue, constant_after_region = shift)
  colnames(E) <- paste0("p", 1:16)
  x <- SpatialExperiment::SpatialExperiment(assays = list(logcounts = E),
    colData = S4Vectors::DataFrame(region = region, row.names = colnames(E)),
    spatialCoords = cbind(x = 1:16, y = 0))
  SingleCellExperiment::altExp(x, "RNA") <- SingleCellExperiment::SingleCellExperiment(
    assays = list(logcounts = rbind(gene = stats::setNames(yvalue, colnames(E)))))
  r <- findCorrelatedFeatures(x, mz = "protein", ST.assay = "RNA", SM.slot = "logcounts", ST.slot = "logcounts",
    covariates = "region", nfeatures = NULL)
  gene <- r[r$modality == "gene", ]
  expect_gt(gene$correlation_raw, .99)
  expect_equal(gene$correlation, 0, tolerance = 1e-12)
  expect_true(is.na(r$correlation[r$features == "constant_after_region"]))
  expect_true(all(is.na(r$mz)))
  expect_s3_class(imageMZPlot(x, "protein", slot = "logcounts"), "ggplot")
})

test_that("the workflow feeds full-observation spatial representations into integration", {
  x <- nativeWorkflowFixture(610)
  r <- runSpaMTPWorkflow(x, list(main = list(type = "proteomics", layer = "logcounts"),
    RNA = list(type = "transcriptomics", layer = "logcounts")), npcs = 3, max_features = 10,
    clusters = 3, association_features = 0, regions = "region",
    structure = list(primary = "spatial", reference = "region", umap = FALSE, k_grid = c(2, 3), stability_blocks = 3))
  expect_equal(nrow(r$analysis$main$spatial_pca$scores), 610)
  expect_identical(r$analysis$main$pca$features, r$analysis$main$spatial_pca$features)
  expected <- do.call(cbind, lapply(r$analysis, function(a) {
    z <- a$spatial_pca$scores
    scale(z) / sqrt(ncol(z))
  }))
  expect_equal(r$joint$embedding, expected, tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(nrow(r$structure$metrics), 12)
  expect_true(all(r$structure$metrics$stability_omissions == 3))
})

test_that("native pathway displays use FDR and a dissimilarity distance", {
  db <- list(analytehaspathway = data.frame(pathwayRampId = c("p1", "p1", "p2", "p2", "p3"),
    rampId = c("RAMP_G_1", "RAMP_G_2", "RAMP_G_1", "RAMP_G_2", "RAMP_G_3")),
    pathway = data.frame(pathwayRampId = c("p1", "p2", "p3"), pathwayName = c("one", "two", "three")))
  t <- data.frame(pathwayRampId = c("p1", "p2", "p3"), pathwayName = c("one", "two", "three"),
    Cluster_id = "A", NES = c(1, -1, 2), size = c(2, 2, 1), pval = .001, padj = .8)
  p <- plotRegionalPathways(t, database = db)
  expect_equal(attr(p, "pathway_distance")["p1", "p2"], 0)
  expect_equal(attr(p, "pathway_distance")["p1", "p3"], 1)
  expect_true(all(p$data$.FDR == "Above cutoff / unavailable"))
  expect_s3_class(plotRegionalPathways(t[1, ], database = db), "ggplot")
  members <- data.frame(pathwayRampId = c("p1", "p2", "p3"),
    used_members = I(list(c("a", "b"), c("b", "a"), c("b", "c"))))
  attr(t, "pathway_coverage") <- list(A = members, B = members)
  groups <- .wf_pathway_equivalence(t)
  expect_equal(sort(groups$pathway_labels), c(1L, 2L))
  attr(t, "pathway_coverage")$B$used_members[[2]] <- c("b", "c")
  expect_equal(nrow(.wf_pathway_equivalence(t)), 3L)
})

test_that("complete regional rankings agree with direct competitive enrichment", {
  ids <- paste0("RAMP_G_", seq_len(40))
  sets <- list(P1 = ids[1:7], P2 = ids[12:22], P3 = ids[c(1, 13, 25, 37)])
  db <- list(source_df = data.frame(rampId = ids, commonName = ids),
    analytehaspathway = data.frame(rampId = unlist(sets, use.names = FALSE),
      pathwayRampId = rep(names(sets), lengths(sets))),
    pathway = data.frame(pathwayRampId = names(sets), pathwayName = names(sets)))
  index <- buildPathwayIndex(db, gene_mapping = "ramp", organism = "custom")
  withr::local_seed(323)
  ranks <- matrix(stats::rnorm(40), dimnames = list(ids, "A"))
  result <- withr::with_seed(12, findRegionalPathways(nativeWorkflowFixture(),
    ident = "region", ranks = list(genes = ranks), pathway_index = index,
    rank_scale = "none", min_path_size = 3, max_path_size = 15,
    nPermSimple = 1000, verbose = FALSE))
  direct <- withr::with_seed(12, as.data.frame(fgsea::fgseaMultilevel(
    pathways = sets, stats = sort(stats::setNames(ranks[, 1], ids), decreasing = TRUE),
    minSize = 3, maxSize = 15, nPermSimple = 1000, eps = 1e-10, nproc = 1)))
  direct <- direct[match(result$pathwayRampId, direct$pathway), ]
  expect_equal(result$NES, direct$NES, tolerance = 1e-12)
  expect_equal(result$pval, direct$pval, tolerance = 1e-12)
  expect_equal(result$padj, stats::p.adjust(direct$pval, "BH", n = length(sets)))
  expect_setequal(names(attr(result, "ranks")$A), ids)
  expect_true(all(result$tested_family == 3))
  expect_equal(attr(result, "pathway_coverage")$A$used_size, lengths(sets), ignore_attr = TRUE)
})

test_that("sparse compound assays retain exact weights and exclude ambiguous masses", {
  E <- Matrix::Matrix(matrix(c(0, 3, 8, 4, 2, 6), nrow = 3,
    dimnames = list(c("mz-10", "mz-20", "mz-30"), c("a", "b"))), sparse = TRUE)
  object <- SingleCellExperiment::SingleCellExperiment(assays = list(workflow = E),
    rowData = S4Vectors::DataFrame(mz = c(10, 20, 30)))
  annotations <- data.frame(observed_mz = c(10, 20, 30), Adduct = "M+H",
    Ramp_IDs = c("RAMP_C_1", "RAMP_C_1", "RAMP_C_2; RAMP_C_3"), Score = .9,
    MassScore = .99, ChemicalScore = 1, IsotopeScore = NA_real_, AdductNetworkScore = NA_real_)
  store <- list(results = annotations, metadata = list(engine = "indexed-chemical-v2", ramp_version = "fixture"))
  S4Vectors::metadata(object)$mz_annotation <- store
  db <- list(chem_props = data.frame(ramp_id = paste0("RAMP_C_", 1:3),
      common_name = letters[1:3], chem_source_id = paste0("hmdb:", 1:3)),
    source_df = data.frame(rampId = paste0("RAMP_C_", 1:3)))
  result <- createPathwayAssay(object, assay = "main", slot = "workflow", database = db,
    metabolite_ambiguity = "exclude", verbose = FALSE)
  mapped <- SingleCellExperiment::altExp(result, "pathway")
  expect_identical(SummarizedExperiment::assayNames(mapped), "workflow")
  expect_identical(rownames(mapped), "RAMP_C_1")
  expect_equal(as.numeric(SummarizedExperiment::assay(mapped)), as.numeric(Matrix::colMeans(E[1:2, ])))
  audit <- S4Vectors::metadata(mapped)$pathway_mapping
  expect_identical(audit$excluded_conflicts, "mz-30")
  expect_equal(as.numeric(audit$weights), c(.5, .5, 0))
  expect_equal(nrow(audit$inputs), 4)
  # A second MS modality must use its own stored candidates, even when the
  # feature names and masses happen to be identical to the main experiment.
  SingleCellExperiment::altExp(object, "MS2") <- SingleCellExperiment::SingleCellExperiment(
    assays = list(workflow = E), rowData = SummarizedExperiment::rowData(object))
  second <- store; second$results$Ramp_IDs <- "RAMP_C_3"
  S4Vectors::metadata(object)$modality_annotations <- list(main = list(mz_annotation = store),
    MS2 = list(mz_annotation = second))
  result <- createPathwayAssay(object, assay = "MS2", slot = "workflow", database = db,
    metabolite_ambiguity = "exclude", verbose = FALSE)
  expect_identical(rownames(SingleCellExperiment::altExp(result, "pathway")), "RAMP_C_3")
  expect_equal(as.numeric(.assayData(result, "pathway", "workflow")), as.numeric(Matrix::colMeans(E)))
})

test_that("the public network accepts audited compound identities and exports offline", {
  ids <- paste0("RAMP_C_", 1:3)
  db <- list(source_df = data.frame(rampId = ids, commonName = letters[1:3], sourceId = paste0("hmdb:", 1:3)),
    chem_props = data.frame(ramp_id = ids, common_name = letters[1:3], chem_source_id = paste0("hmdb:", 1:3)),
    pathway = data.frame(pathwayRampId = "P1", pathwayName = "Test pathway", sourceId = "WP1", type = "WikiPathways"),
    analytehaspathway = data.frame(rampId = ids, pathwayRampId = "P1"),
    ramp_db_metadata = list(ramp_version = "fixture"),
    ramp_wikipathway = list(`Test pathway` = list(id = "WP1", title = "Test pathway",
      mixedEdges = data.frame(src = ids[1:2], dest = ids[2:3], directed = 1L, reaction_type = 1L))))
  index <- buildPathwayIndex(db, gene_mapping = "ramp", organism = "custom")
  E <- matrix(1:18, 3, dimnames = list(ids, paste0("spot", 1:6)))
  object <- SpatialExperiment::SpatialExperiment(assays = list(workflow = E),
    colData = S4Vectors::DataFrame(region = rep(c("A", "B"), each = 3)),
    spatialCoords = cbind(x = 1:6, y = 1))
  S4Vectors::metadata(object)$pathway_mapping <- list(provenance = index$provenance,
    ambiguity = "exclude", original_features = list(RAMP_C_1 = "mz-1"))
  m <- findAllDEMs(object, ident = "region", method = "markers", assay = "main", slot = "workflow")
  regional <- cbind(db$pathway, Cluster_id = "A", NES = 1, padj = .5)
  regional$leadingEdge <- list(ids[1:2])
  attr(regional, "pathway_index") <- index$provenance
  folder <- withr::local_tempdir()
  path <- pathwayNetworkPlots(object, ident = "region", regpathway = regional,
    DE.list = list(metabolites = m$DEMs), SM_slot = "workflow", path = folder,
    analyte_types = "metabolites", database = db, pathway_index = index,
    organism = "custom", verbose = FALSE)
  html <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_false(grepl('<script src=', html, fixed = TRUE))
  expect_true(grepl('https://d3js.org v7.9.0', html, fixed = TRUE))
  payload <- sub('.*const payload = ([^\n]+);\n.*', '\\1', html)
  wire <- jsonlite::fromJSON(payload, simplifyVector = FALSE)
  payload <- jsonlite::fromJSON(payload)
  expect_false(payload$metadata$annotation_controls)
  expect_identical(payload$metadata$effect_label, "Regional effect (workflow units)")
  expect_equal(payload$metadata$pathway_coverage$measured_size, 3L)
  expect_equal(payload$metadata$pathway_coverage$used_size, 3L)
  expect_setequal(vapply(wire$networks[[1]][[1]]$nodes, function(n) n$id, character(1)), ids)
})
