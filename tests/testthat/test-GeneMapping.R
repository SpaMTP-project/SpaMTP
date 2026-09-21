gene_mapping_fixture <- function() {
  # Fictional genes and identifier relationships for deterministic offline tests.
  reference <- data.frame(
    hgnc_id = paste0("HGNC:", 1:4), symbol = c("GENEA", "GENEB", "GENEC", "GENED"),
    status = "Approved", entrez_id = as.character(101:104),
    ensembl_gene_id = paste0("ENSG000", 1:4), uniprot_ids = paste0("U", 1:4),
    prev_symbol = c("OLDA", "", "OLDC", ""),
    alias_symbol = c("SHARED|GENEB|ALIASA", "", "SHARED", "")
  )
  source <- data.frame(
    rampId = paste0("RAMP_G_", c(1, 2, 3, 4, 5, 5, 6, 7)),
    sourceId = c("entrez:101", "uniprot:U1", "gene_symbol:GENEB", "entrez:103",
                 "entrez:101", "entrez:103", "gene_symbol:OUTSIDE", "entrez:101"),
    commonName = c(NA, "OLDA", "GENEB", "GENEC", NA, NA, "OUTSIDE", "GENEC")
  )
  list(source_df = source, analyte = data.frame(rampId = unique(source$rampId)),
       analytehaspathway = data.frame(
         rampId = paste0("RAMP_G_", c(1, 2, 3, 2, 4, 5, 6)),
         pathwayRampId = c("P1", "P1", "P1", "P2", "P2", "conflict", "orphan")),
       pathway = data.frame(pathwayRampId = c("P1", "P2", "conflict", "orphan")),
       gene_reference = reference)
}

test_that("stable IDs recover missing symbols and merge records of one gene", {
  db <- gene_mapping_fixture()
  index <- buildGeneMappingIndex(db$source_df, db$gene_reference)
  out <- mapGeneIdentifiers(c("GENEA", "gene_symbol:olda", "101", "entrez:101",
                              "ENSG0001.12", "hgnc:1", "uniprot:U1", "RAMP_G_2"), index = index)
  expect_true(all(out$status == "mapped"))
  expect_true(all(out$gene_id == "HGNC:1"))
  expect_true(all(out$canonical_ramp_id == "RAMP_G_1"))
  expect_true(all(out$ramp_ids == "RAMP_G_1;RAMP_G_2;RAMP_G_7"))
  expect_equal(out$method[2], "previous_symbol")
  expect_equal(out$method[5], "ensembl_id")
  expect_equal(index$nodes$status[index$nodes$rampId == "RAMP_G_5"], "conflicting_record")
  expect_equal(index$nodes$gene_id[index$nodes$rampId == "RAMP_G_7"], "HGNC:1")
})

test_that("approved symbols win and ambiguous aliases never expand into genes", {
  db <- gene_mapping_fixture()
  index <- buildGeneMappingIndex(db$source_df, db$gene_reference)
  expect_equal(mapGeneIdentifiers("GENEB", index = index)$gene_id, "HGNC:2")
  expect_equal(mapGeneIdentifiers("ALIASA", index = index)$gene_id, "HGNC:1")
  expect_warning(out <- mapGeneIdentifiers(c("SHARED", "GENED", "MISSING", "RAMP_G_5"), index = index), "unresolved")
  expect_equal(out$status, c("ambiguous", "not_in_ramp", "unrecognized_identifier", "ambiguous"))
  expect_true(all(is.na(out$canonical_ramp_id)))
  expect_equal(out$candidate_gene_ids[c(1, 4)], c("HGNC:1;HGNC:3", "HGNC:1;HGNC:3"))
  expect_error(mapGeneIdentifiers("SHARED", index = index, ambiguous = "error"), "Ambiguous gene mapping")
  expect_equal(mapGeneIdentifiers("RAMP_G_6", index = index)$status, "mapped_ramp_only")
})

test_that("record conflicts cannot create transitive merges of distinct genes", {
  db <- gene_mapping_fixture()
  index <- buildGeneMappingIndex(db$source_df, db$gene_reference)
  out <- mapGeneIdentifiers(c("GENEA", "GENEC"), index = index)
  expect_equal(out$canonical_ramp_id, c("RAMP_G_1", "RAMP_G_4"))
  expect_false(any(grepl("RAMP_G_5", out$ramp_ids)))
  view <- .gene_pathway_view(db, db, gene_index = index)
  expect_false("RAMP_G_5" %in% view$resources$analytehaspathway$rampId)
  expect_equal(sum(view$resources$analytehaspathway$rampId == "RAMP_G_1" &
                     view$resources$analytehaspathway$pathwayRampId == "P1"), 1L)
})

test_that("a gene-specific stable ID disambiguates a shared protein accession", {
  db <- gene_mapping_fixture()
  db$gene_reference$uniprot_ids[3] <- "U1"
  db$source_df <- rbind(db$source_df,
                        data.frame(rampId = "RAMP_G_1", sourceId = "uniprot:U1", commonName = NA))
  index <- buildGeneMappingIndex(db$source_df, db$gene_reference)
  expect_equal(index$nodes$gene_id[index$nodes$rampId == "RAMP_G_1"], "HGNC:1")
  expect_equal(index$nodes$status[index$nodes$rampId == "RAMP_G_2"], "conflicting_record")
  expect_equal(index$nodes$status[index$nodes$rampId == "RAMP_G_5"], "conflicting_record")
  expect_warning(out <- mapGeneIdentifiers("uniprot:U1", index = index), "unresolved")
  expect_equal(out$status, "ambiguous")
})

test_that("differential features are never arbitrarily chosen after gene merging", {
  db <- gene_mapping_fixture()
  index <- buildGeneMappingIndex(db$source_df, db$gene_reference)
  de <- data.frame(gene = c("GENEA", "OLDA", "GENEB"), cluster = "N",
                    avg_log2FC = c(1, 1, -1), p_val_adj = 0.01)
  compact <- .gene_prepare_de(de, index)
  expect_equal(nrow(compact), 2L)
  expect_equal(compact$rampId, c("RAMP_G_1", "RAMP_G_3"))
  network <- .pn_prepare_gene_de(de, db$source_df, gene_index = index)
  expect_setequal(network$rampId, c("RAMP_G_1", "RAMP_G_2", "RAMP_G_7", "RAMP_G_3"))
  de$avg_log2FC[2] <- 2
  expect_error(.gene_prepare_de(de, index), "different differential statistics")
})

test_that("foreground and background collapse by gene before Fisher testing", {
  db <- gene_mapping_fixture()
  out <- fishersPathwayAnalysis(
    list(genes = c("GENEA", "OLDA", "RAMP_G_2")),
    universe = list(genes = c("101", "GENEB", "ENSG0003.5", "RAMP_G_7")),
    database = db, min_path_size = 1, pathway_all_info = TRUE, verbose = FALSE
  )
  expect_equal(out$pathwayRampId, c("P1", "P2"))
  expect_equal(out$foreground_analytes_number, c(1L, 1L))
  expect_equal(out$background_analytes_number, c(3L, 3L))
  expect_equal(out$total_in_pathways, c(2L, 2L))
  expect_equal(out$analytes_in_pathways, c(1L, 1L))
  expect_equal(out$p_val, c(2 / 3, 2 / 3))
  expect_true(all(out$gene_name_list == "GENEA"))
  audit <- attr(out, "enrichment")$gene_mapping
  expect_equal(audit$mode, "hgnc")
  expect_equal(audit$conflicting_records$rampId, "RAMP_G_5")
  expect_equal(length(unique(audit$groups$gene_id)), 3L)
  expect_equal(nrow(audit$inputs), 7L)
})

test_that("mapping status remains explicit for unresolved foreground and universe", {
  db <- gene_mapping_fixture()
  messages <- character()
  out <- withCallingHandlers(fishersPathwayAnalysis(
    list(genes = c("GENEA", "SHARED")), universe = list(genes = c("GENEA", "GENEB", "SHARED")),
    database = db, min_path_size = 1, verbose = FALSE
  ), warning = function(w) {
    messages <<- c(messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_length(messages, 2L)
  expect_true(all(grepl("unmapped", messages)))
  expect_equal(attr(out, "enrichment")$unmapped$Analyte$genes, "SHARED")
  expect_equal(attr(out, "enrichment")$gene_mapping$inputs$status, c("mapped", "ambiguous", "mapped"))
})

test_that("SpaMTPData paired experiments retain alignment and assays after mapping", {
  skip_if_not_installed("SpaMTPData")
  object <- SpaMTPData::spaMTPExampleData()
  paired <- SingleCellExperiment::altExp(object, "transcriptome")
  SummarizedExperiment::rowData(paired)$symbol <- c("GENEA", "OLDA", "GENEB", "SHARED", "GENED")
  SingleCellExperiment::altExp(object, "transcriptome") <- paired
  db <- gene_mapping_fixture()
  index <- buildGeneMappingIndex(db$source_df, db$gene_reference)
  expect_warning(mapped <- annotateGeneIdentifiers(object, assay = "transcriptome",
      id_column = "symbol", index = index), "2 gene input")
  target <- SingleCellExperiment::altExp(mapped, "transcriptome")
  expect_identical(SummarizedExperiment::assay(target), SummarizedExperiment::assay(paired))
  expect_identical(SummarizedExperiment::assay(mapped), SummarizedExperiment::assay(object))
  expect_identical(SpatialExperiment::spatialCoords(mapped), SpatialExperiment::spatialCoords(object))
  expect_identical(rownames(target), rownames(paired))
  expect_equal(SummarizedExperiment::rowData(target)$spamtp_gene_id[1:3], c("HGNC:1", "HGNC:1", "HGNC:2"))
  expect_equal(SummarizedExperiment::rowData(target)$spamtp_gene_status,
               c("mapped", "mapped", "mapped", "ambiguous", "not_in_ramp"))
  expect_equal(S4Vectors::metadata(target)$gene_mapping$provenance$organism, "Homo sapiens")
  expect_error(annotateGeneIdentifiers(object, assay = "transcriptome", id_column = "missing", index = index), "id_column")
  S4Vectors::metadata(object)$SpaMTPData <- list(organism = "Mus musculus")
  expect_error(annotateGeneIdentifiers(object, assay = "transcriptome", id_column = "symbol", index = index), "Experiment species")
})

test_that("mapping preserves row order, empty inputs and duplicate feature IDs", {
  db <- gene_mapping_fixture()
  index <- buildGeneMappingIndex(db$source_df, db$gene_reference)
  x <- c("GENEB", "GENEA", "GENEB")
  out <- mapGeneIdentifiers(x, index = index)
  expect_identical(out$input, x)
  expect_identical(out$gene_id, c("HGNC:2", "HGNC:1", "HGNC:2"))
  expect_equal(nrow(mapGeneIdentifiers(character(), index = index)), 0L)
  expect_error(mapGeneIdentifiers(NA_character_, index = index), "without missing")
  expect_error(mapGeneIdentifiers("", index = index), "blank")
  expect_error(buildGeneMappingIndex(db$source_df, db$gene_reference, organism = "Mus musculus"), "no orthology")
  expect_error(mapGeneIdentifiers("GENEA", index = index, organism = "Mus musculus"), "organism")
  duplicate <- rbind(db$gene_reference, db$gene_reference[1, ])
  expect_error(buildGeneMappingIndex(db$source_df, duplicate), "unique approved")
})

test_that("custom non-human and historical RaMP-only workflows stay explicit", {
  db <- gene_mapping_fixture()
  index <- buildGeneMappingIndex(db$source_df, db$gene_reference)
  expect_error(fishersPathwayAnalysis(list(genes = "GENEA"), database = db,
                                     gene_index = index, organism = "Mus musculus"), "organism")
  db$gene_reference <- NULL
  out <- fishersPathwayAnalysis(list(genes = "GENEB"), database = db,
                                min_path_size = 1, organism = "Mus musculus", verbose = FALSE)
  expect_equal(attr(out, "enrichment")$gene_mapping$mode, "ramp")
  bad <- index
  bad$nodes <- bad$nodes[-1, ]
  expect_error(fishersPathwayAnalysis(list(genes = "GENEA"), database = db,
                                     gene_index = bad), "does not match")
})

test_that("regional pathway ranks and memberships use the same gene identities", {
  skip_if_not_installed("SpaMTPData")
  object <- SpaMTPData::spaMTPExampleData()
  counts <- matrix(seq_len(36), 3, dimnames = list(c("GENEA", "GENEB", "GENEC"), colnames(object)))
  SingleCellExperiment::altExp(object, "transcriptome") <-
    SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
  db <- gene_mapping_fixture()
  db$chem_props <- data.frame()
  db$pathway$pathwayName <- c("first", "second", "conflict", "orphan")
  de <- data.frame(gene = rep(c("GENEA", "OLDA", "GENEB", "GENEC"), 2),
                    cluster = rep(c("edge", "core"), each = 4),
                    avg_log2FC = rep(c(2, 2, -1, 0), 2), p_val_adj = 0.01)
  seen <- list()
  testthat::local_mocked_bindings(fgsea = function(pathways, stats, ...) {
    seen[[length(seen) + 1L]] <<- list(pathways = pathways, stats = stats)
    data.table::data.table(pathway = character(), pval = numeric(), padj = numeric(),
                           NES = numeric(), leadingEdge = list())
  }, .package = "fgsea")
  result <- findRegionalPathways(object, ident = "region", DE.list = list(genes = de),
                                 analyte_types = "genes", database = db,
                                 min_path_size = 2, verbose = FALSE)
  expect_length(seen, 2L)
  expect_setequal(names(seen[[1]]$stats), c("RAMP_G_1", "RAMP_G_3", "RAMP_G_4"))
  expect_setequal(seen[[1]]$pathways$P1, c("RAMP_G_1", "RAMP_G_3"))
  expect_setequal(seen[[1]]$pathways$P2, c("RAMP_G_1", "RAMP_G_4"))
  expect_equal(nrow(result), 0L)
  expect_true("inputs" %in% names(attr(result, "gene_mapping")))
})

test_that("source fingerprints reject changed identities with unchanged RaMP IDs", {
  db <- gene_mapping_fixture()
  index <- buildGeneMappingIndex(db$source_df, db$gene_reference)
  for (column in c("sourceId", "commonName")) {
    changed <- db
    changed$source_df[[column]][changed$source_df$rampId == "RAMP_G_3"] <-
      if (column == "sourceId") "entrez:103" else "GENEC"
    expect_error(.gene_pathway_view(changed, changed, gene_index = index), "does not match")
    expect_error(fishersPathwayAnalysis(list(genes = "GENEB"), database = changed,
        gene_index = index, min_path_size = 1, verbose = FALSE), "does not match")
    expect_error(mapGeneIdentifiers("GENEB", database = changed, index = index), "does not match")
  }
  changed <- db
  changed$source_df <- changed$source_df[, c("rampId", "sourceId")]
  expect_error(.gene_pathway_view(changed, changed, gene_index = index), "does not match")
  old <- index
  old$provenance$schema_version <- 1L
  old$provenance$source_fingerprint <- NULL
  expect_error(.gene_pathway_view(db, db, gene_index = old), "Rebuild")
  expect_error(mapGeneIdentifiers("GENEA", index = old), "Rebuild")
})

test_that("source fingerprints ignore row order, duplicate evidence and unrelated columns", {
  db <- gene_mapping_fixture()
  index <- buildGeneMappingIndex(db$source_df, db$gene_reference)
  changed <- db
  changed$source_df <- rbind(changed$source_df[nrow(changed$source_df):1L, ], changed$source_df[1L, ])
  changed$source_df$note <- "irrelevant display metadata"
  changed$source_df$sourceId <- toupper(changed$source_df$sourceId)
  rownames(changed$source_df) <- NULL
  reordered <- buildGeneMappingIndex(changed$source_df, changed$gene_reference)
  expect_identical(index$provenance$source_fingerprint, reordered$provenance$source_fingerprint)
  expect_identical(index$nodes, reordered$nodes)
  expect_identical(.gene_pathway_view(changed, changed, gene_index = index)$index, index)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(index, path)
  expect_identical(.gene_pathway_view(changed, changed, gene_index = readRDS(path))$index, index)
})

test_that("shared crossreferences retain unresolved and conflicting competitors", {
  for (competitor in c("RAMP_G_6", "RAMP_G_5")) {
    db <- gene_mapping_fixture()
    db$source_df <- rbind(db$source_df, data.frame(
      rampId = c("RAMP_G_1", competitor), sourceId = "legacy:shared", commonName = NA))
    index <- buildGeneMappingIndex(db$source_df, db$gene_reference)
    expect_warning(out <- mapGeneIdentifiers("legacy:shared", index = index), "unresolved")
    expect_identical(out$status, "ambiguous")
    expect_true(is.na(out$gene_id))
    expect_true(is.na(out$canonical_ramp_id))
    expect_identical(out$ramp_ids, "")
    expect_identical(out$candidate_ramp_ids, paste(c("RAMP_G_1", competitor), collapse = ";"))
    expect_identical(out$candidate_gene_ids, if (competitor == "RAMP_G_6") "HGNC:1" else "HGNC:1;HGNC:3")
    expect_error(mapGeneIdentifiers("legacy:shared", index = index, ambiguous = "error"), "Ambiguous gene mapping")
    expect_warning(enriched <- fishersPathwayAnalysis(list(genes = "legacy:shared"),
        database = db, gene_index = index, min_path_size = 1, verbose = FALSE), "unmapped")
    expect_true(all(enriched$foreground_analytes_number == 0L))
    expect_identical(attr(enriched, "enrichment")$gene_mapping$inputs$status, "ambiguous")
  }
})

test_that("safe raw crossreferences retain their existing resolution", {
  db <- gene_mapping_fixture()
  db$source_df <- rbind(db$source_df, data.frame(
    rampId = c("RAMP_G_1", "RAMP_G_2", "RAMP_G_6"),
    sourceId = c("legacy:same_gene", "legacy:same_gene", "legacy:unique"), commonName = NA))
  index <- buildGeneMappingIndex(db$source_df, db$gene_reference)
  out <- mapGeneIdentifiers(c("legacy:same_gene", "legacy:unique"), index = index, ambiguous = "error")
  expect_identical(out$status, c("mapped", "mapped_ramp_only"))
  expect_identical(out$canonical_ramp_id, c("RAMP_G_1", "RAMP_G_6"))
})

test_that("conflict lookups and repeated queries preserve auditable statuses", {
  db <- gene_mapping_fixture()
  db$source_df <- db$source_df[db$source_df$rampId %in% c("RAMP_G_3", "RAMP_G_5", "RAMP_G_6"), ]
  index <- buildGeneMappingIndex(db$source_df, db$gene_reference)
  expect_setequal(index$conflict_gene_ids, c("HGNC:1", "HGNC:3"))
  genes <- c("GENEA", "GENED", "GENEB", "genea", "GENED", "RAMP_G_6", "outside")
  expect_warning(out <- mapGeneIdentifiers(genes, index = index), "4 gene input")
  expect_identical(out$input, genes)
  expect_identical(out$status, c("conflicting_ramp_records", "not_in_ramp", "mapped",
      "conflicting_ramp_records", "not_in_ramp", "mapped_ramp_only", "mapped_ramp_only"))
  expect_identical(rownames(out), as.character(seq_along(genes)))
  db$source_df <- db$source_df[db$source_df$rampId != "RAMP_G_5", ]
  expect_identical(buildGeneMappingIndex(db$source_df, db$gene_reference)$conflict_gene_ids, character())
})

test_that("detached SpaMTPData RNA retains species protection at the public boundary", {
  skip_if_not_installed("SpaMTPData", minimum_version = "0.99.5")
  path <- tempfile("gene-species-")
  dir.create(path)
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  object <- SpaMTPData::spaMTPExampleData()
  row <- SpaMTPData::spaMTPData("mouse_brain_dhb_striatum", metadata = TRUE)
  saveRDS(object, file.path(path, row$file_name))
  # Deliberate local fixture, not an official release payload.
  loaded <- SpaMTPData::spaMTPData("mouse_brain_dhb_striatum", local_dir = path,
      offline = TRUE, verify = FALSE)
  detached <- SingleCellExperiment::altExp(loaded, "transcriptome")
  db <- gene_mapping_fixture()
  index <- buildGeneMappingIndex(db$source_df, db$gene_reference)
  expect_error(annotateGeneIdentifiers(detached, index = index), "Experiment species.*Mus musculus")
})
