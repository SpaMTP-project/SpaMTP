fisher_fixture <- function() {
  ids <- paste0("RAMP_G_", seq_len(10))
  list(
    source_df = data.frame(rampId = ids, sourceId = paste0("gene_symbol:G", seq_len(10)),
                           commonName = paste0("G", seq_len(10))),
    analyte = data.frame(rampId = ids),
    analytehaspathway = data.frame(rampId = ids,
                                   pathwayRampId = rep(c("P1", "P2"), c(4, 6)),
                                   pathwaySource = "test"),
    pathway = data.frame(pathwayRampId = c("P1", "P2"), pathwayName = c("one", "two"),
                         sourceId = c("P1", "P2"), type = "test", pathwayCategory = "test")
  )
}

fisher_run <- function(db = fisher_fixture(), fg = c(1, 2, 5), ...) {
  fishersPathwayAnalysis(list(genes = paste0("gene_symbol:G", fg)),
                        database = db, min_path_size = 1, verbose = FALSE, ...)
}

test_that("Fisher regression tables use foreground minus overlap", {
  out <- fisher_run()
  expect_equal(out$pathway_id, c("P1", "P2"))
  expect_equal(out$p_val, c(1 / 3, 29 / 30), tolerance = 1e-12)
  expect_equal(out$fdr, p.adjust(c(1 / 3, 29 / 30), "BH"))
  expect_equal(out$analytes_in_pathways, c(2L, 1L))
  expect_equal(out$total_in_pathways, c(4L, 6L))
  expect_equal(out$foreground_analytes_number, c(3L, 3L))
  expect_equal(out$background_analytes_number, c(10L, 10L))
})

test_that("size filtering never changes the background or remaining raw p-values", {
  out <- fisher_run(max_path_size = 4)
  expect_equal(nrow(out), 1L)
  expect_equal(out$p_val, 1 / 3)
  expect_equal(out$background_analytes_number, 10L)
  expect_equal(attr(out, "enrichment")$n_tested, 1L)
})

test_that("measured universe controls all margins and retains non-pathway members", {
  db <- fisher_fixture()
  db$analytehaspathway <- db$analytehaspathway[1:8, ]
  out <- fisher_run(db, universe = list(genes = paste0("RAMP_G_", 1:9)))
  p1 <- out[out$pathway_id == "P1", ]
  expect_equal(p1$background_analytes_number, 9L)
  expect_equal(p1$p_val, fisher.test(matrix(c(2, 2, 1, 4), 2), alternative = "greater")$p.value)
  expect_equal(attr(out, "enrichment")$universe_source, "measured")
  expect_true("RAMP_G_9" %in% attr(out, "enrichment")$universe_ids)
  small <- fisher_run(db, universe = list(genes = paste0("RAMP_G_", 1:5)))
  expect_equal(small$total_in_pathways[match(c("P1", "P2"), small$pathway_id)], c(4L, 1L))
  expect_true(all(small$background_analytes_number == 5L))
  # A measured foreground ID without pathway annotation still contributes to K.
  with_unlinked <- fisher_run(db, fg = c(1, 2, 9),
                             universe = list(genes = paste0("RAMP_G_", 1:9)))
  expect_true(all(with_unlinked$foreground_analytes_number == 3L))
  expect_equal(with_unlinked$p_val[with_unlinked$pathway_id == "P1"], p1$p_val)
  expect_warning(default <- fisher_run(db, fg = c(1, 2, 9)), "without pathway membership")
  expect_equal(attr(default, "enrichment")$excluded_foreground_ids, "RAMP_G_9")
  expect_true(all(default$foreground_analytes_number == 2L))
})

test_that("all eligible zero-overlap pathways belong to the BH family", {
  db <- fisher_fixture()
  db$analytehaspathway <- data.frame(rampId = paste0("RAMP_G_", 1:8),
                                    pathwayRampId = rep(c("P1", "P2", "P3"), c(3, 2, 3)))
  # P3 deliberately lacks display metadata. It must still be tested.
  universe <- list(genes = paste0("RAMP_G_", 1:10))
  out <- fisher_run(db, fg = 1:3, universe = universe)
  expect_equal(nrow(out), 3L)
  expect_equal(out$p_val, c(1 / 120, 1, 1))
  expect_equal(out$fdr, c(3 / 120, 1, 1))
  expect_equal(out$analytes_in_pathways, c(3L, 0L, 0L))
  expect_equal(out$pathway_name[out$pathwayRampId == "P3"], "P3")
  expect_true(is.na(out$type[out$pathwayRampId == "P3"]))
  filtered <- fisher_run(db, fg = 1:3, universe = universe, pval_cutoff = 0.01)
  expect_equal(nrow(filtered), 0L)
  expect_equal(attr(filtered, "enrichment")$n_tested, 3L)
  selected <- fisher_run(db, fg = 1:3, universe = universe, pval_cutoff = 0.03)
  expect_equal(selected$fdr, 3 / 120)
  expect_equal(attr(selected, "enrichment")$n_tested, 3L)
  detailed <- fisher_run(db, fg = 1:3, universe = universe, pathway_all_info = TRUE)
  expect_equal(detailed[names(out)], out[names(out)])
  expect_equal(attr(detailed, "enrichment"), attr(out, "enrichment"))
  expect_equal(detailed$gene_id_list[2:3], c("", ""))
  expect_true(all(detailed$adduct_info == ""))
})

test_that("duplicate links, pathway metadata, aliases and repeated inputs do not inflate counts", {
  db <- fisher_fixture()
  alias <- db$source_df[1, ]
  alias$sourceId <- "entrez:1"
  db$source_df <- rbind(db$source_df, alias, db$source_df)
  duplicate <- db$analytehaspathway
  duplicate$pathwaySource <- "second provenance"
  db$analytehaspathway <- rbind(db$analytehaspathway, duplicate)
  db$pathway <- rbind(db$pathway, db$pathway)
  out <- fishersPathwayAnalysis(
    list(genes = c("GENE_SYMBOL:G1", "entrez:1", "g2", "ramp_g_5", "g2")),
    universe = list(genes = c(paste0("G", 1:10), "entrez:1")),
    database = db, min_path_size = 1, verbose = FALSE
  )
  expect_equal(out$p_val, c(1 / 3, 29 / 30))
  expect_equal(nrow(out), 2L)
  expect_equal(out$foreground_analytes_number, c(3L, 3L))
  expect_equal(out$background_analytes_number, c(10L, 10L))
})

test_that("empty foregrounds and empty eligible families have well-defined results", {
  for (alternative in c("greater", "less", "two.sided")) {
    out <- fishersPathwayAnalysis(list(genes = character()),
                                  database = fisher_fixture(), min_path_size = 1,
                                  verbose = FALSE, alternative = alternative)
    expect_equal(out$p_val, c(1, 1))
    expect_equal(out$fdr, c(1, 1))
    expect_equal(out$foreground_analytes_number, c(0L, 0L))
  }
  empty <- fisher_run(max_path_size = 2, pathway_all_info = TRUE)
  expect_equal(nrow(empty), 0L)
  expect_true(all(c("p_val", "fdr", "pathwayRampId", "gene_id_list") %in% names(empty)))
  expect_equal(attr(empty, "enrichment")$n_tested, 0L)
  expect_length(attr(empty, "enrichment")$universe_ids, 10L)
  full <- fisher_run(fg = 1:10)
  expect_equal(full$p_val, c(1, 1))
  db <- fisher_fixture()
  db$analytehaspathway <- db$analytehaspathway[1:4, ]
  unlinked <- fisher_run(db, fg = 9, universe = list(genes = c("G9", "G10")))
  expect_equal(nrow(unlinked), 0L)
  expect_true("pathwayRampId" %in% names(unlinked))
  expect_equal(attr(unlinked, "enrichment")$foreground_ids, "RAMP_G_9")
})

test_that("unmapped values warn and remain auditable", {
  expect_warning(out <- fisher_run(fg = c(1, 2, 5, 999)), "Analyte\\$genes: 1 unmapped")
  expect_equal(out$p_val, c(1 / 3, 29 / 30))
  expect_equal(attr(out, "enrichment")$unmapped$Analyte$genes, "gene_symbol:G999")
  expect_warning(out <- fisher_run(universe = list(genes = c(paste0("G", 1:10), "unknown"))),
                 "universe\\$genes: 1 unmapped")
  expect_equal(out$background_analytes_number, c(10L, 10L))
  expect_equal(attr(out, "enrichment")$unmapped$universe$genes, "unknown")
  expect_warning(empty_fg <- fisher_run(fg = 999), "unmapped")
  expect_equal(empty_fg$p_val, c(1, 1))
})

test_that("invalid universes and inputs fail before invalid tables can be constructed", {
  expect_error(fisher_run(universe = list(genes = c("G1", "G2"))), "outside universe: RAMP_G_5")
  expect_error(fisher_run(universe = list(genes = character())), "universe is empty")
  expect_error(fisher_run(universe = list(metabolites = "RAMP_C_1")), "same biological modalities")
  expect_error(fisher_run(max_path_size = 0), "size limits")
  expect_error(fisher_run(max_path_size = 1.5), "size limits")
  expect_error(fisher_run(pval_cutoff = NA_real_), "pval_cutoff")
  expect_error(fisher_run(pval_cutoff = 1.1), "pval_cutoff")
  expect_error(fisher_run(pathway_all_info = NA), "pathway_all_info")
  expect_error(fisher_run(alternative = "invalid"), "arg")
  for (invalid in list(list(gene = "G1"), list(genes = "G1", genes = "G2"),
                       list(genes = NA_character_), list(genes = " "),
                       list(mzs = "oops149.1"), list(mzs = -1), list(mzs = Inf))) {
    expect_error(fishersPathwayAnalysis(invalid, database = fisher_fixture()), "Analyte")
  }
})

test_that("mixed genes and compounds share one deduplicated measured universe", {
  db <- fisher_fixture()
  for (key in c("source_df", "analyte", "analytehaspathway")) {
    db[[key]]$rampId <- sub("RAMP_G_([6-9]|10)$", "RAMP_C_\\1", db[[key]]$rampId)
  }
  db$source_df$sourceId[6:10] <- paste0("hmdb:HMDB", 6:10)
  out <- fishersPathwayAnalysis(
    list(genes = c("G1", "G2"), metabolites = "hmdb:HMDB6"),
    universe = list(genes = paste0("G", 1:5), metabolites = paste0("hmdb:HMDB", 6:10)),
    database = db, min_path_size = 1, pathway_all_info = TRUE, verbose = FALSE
  )
  expect_equal(out$p_val, c(1 / 3, 29 / 30))
  expect_equal(out$gene_name_list[1], "G1;G2")
  expect_equal(out$metabolite_id_list[2], "hmdb:hmdb6")
  expect_true(all(out$background_analytes_number == 10L))
})

test_that("all alternatives agree with independently constructed membership tables", {
  set.seed(90513)
  db <- fisher_fixture()
  # Include all IDs through an anchor pathway, then arbitrary overlapping sets.
  sets <- c(list(anchor = paste0("RAMP_G_", 1:10)),
            setNames(replicate(30, sample(paste0("RAMP_G_", 1:10), sample(1:10, 1)),
                               simplify = FALSE), paste0("random", 1:30)))
  db$analytehaspathway <- data.frame(rampId = unlist(sets, use.names = FALSE),
                                    pathwayRampId = rep(names(sets), lengths(sets)))
  for (alternative in c("greater", "less", "two.sided")) {
    out <- fisher_run(db, alternative = alternative)
    expected <- vapply(sets, function(members) {
      ids <- paste0("RAMP_G_", 1:10)
      tab <- table(factor(ids %in% members, levels = c(TRUE, FALSE)),
                   factor(ids %in% paste0("RAMP_G_", c(1, 2, 5)), levels = c(TRUE, FALSE)))
      fisher.test(tab, alternative = alternative)$p.value
    }, numeric(1))
    expected <- expected[out$pathwayRampId]
    expect_equal(out$p_val, unname(expected), tolerance = 1e-12)
    expect_equal(out$fdr, unname(p.adjust(expected, "BH")), tolerance = 1e-12)
    if (alternative == "greater") {
      hyper <- with(out, phyper(analytes_in_pathways - 1, total_in_pathways,
                               background_analytes_number - total_in_pathways,
                               foreground_analytes_number, lower.tail = FALSE))
      expect_equal(out$p_val, hyper, tolerance = 1e-12)
    }
  }
})

test_that("m/z and metabolite universes use consistent current and legacy mappings", {
  db <- fisher_fixture()
  for (key in c("source_df", "analyte", "analytehaspathway")) {
    db[[key]]$rampId <- sub("RAMP_G_", "RAMP_C_", db[[key]]$rampId)
  }
  db$source_df$sourceId <- paste0("hmdb:HMDB", seq_len(10))
  chemical <- data.frame(
    ramp_id = paste0("RAMP_C_", 1:3), chem_source_id = paste0("hmdb:HMDB", 1:3),
    common_name = paste0("compound", 1:3), chem_data_source = "hmdb",
    mol_formula = c("C5H8O5", "C6H12O6", "C3H6O3"),
    monoisotop_mass = c(148.037173366, 180.063388104, 90.031694052)
  )
  index <- buildMZAnnotationIndex(chemical, adducts = "M+H", infer_structure = "never")
  mass <- chemical$monoisotop_mass + 1.007276466621
  out <- fishersPathwayAnalysis(
    list(mzs = c(mass[1], mass[1])), universe = list(mzs = mass),
    database = db, index = index, ppm_error = 0.001, min_path_size = 1,
    pathway_all_info = TRUE, verbose = FALSE
  )
  expect_equal(out$foreground_analytes_number, 1L)
  expect_equal(out$background_analytes_number, 3L)
  expect_equal(out$total_in_pathways, 3L)
  expect_match(out$adduct_info, "M\\+H")
  mixed <- fishersPathwayAnalysis(
    list(mzs = paste0("mz-", mass[1]), metabolites = "hmdb:HMDB1"),
    universe = list(metabolites = paste0("hmdb:HMDB", 1:10)),
    database = db, db = chemical, adducts = "M+H", ppm_error = 0.001,
    infer_structure = "never", min_path_size = 1, verbose = FALSE
  )
  expect_equal(mixed$p_val, c(0.4, 1))
  expect_true(all(mixed$foreground_analytes_number == 1L))
  legacy <- data.frame(formula = chemical$mol_formula,
                       exactmass = chemical$monoisotop_mass,
                       isomers = paste0("HMDB", 1:3), isomers_names = chemical$common_name)
  old <- fishersPathwayAnalysis(
    list(mzs = mass[1]), universe = list(metabolites = paste0("hmdb:HMDB", 1:10)),
    database = db, db = legacy, adducts = "M+H", ppm_error = 0.001,
    infer_structure = "never", min_path_size = 1, verbose = FALSE
  )
  expect_equal(old$p_val, mixed$p_val)
  expect_setequal(attr(old, "enrichment")$foreground_ids, "RAMP_C_1")
})
