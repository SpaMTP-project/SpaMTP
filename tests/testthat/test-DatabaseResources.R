test_that("custom database resources can be used without a Hub lookup", {
  example_database <- list(
    ramp_db_metadata = list(ramp_version = "example")
  )
  database <- loadSpaMTPDatabase(
    "ramp_db_metadata",
    database = example_database,
    refresh = TRUE
  )

  expect_named(database, "ramp_db_metadata")
  expect_type(database$ramp_db_metadata, "list")
  expect_identical(database$ramp_db_metadata$ramp_version, "example")
})

test_that("SpaMTPdb resources load from an offline staging directory", {
  skip_if_not_installed("SpaMTPdb")
  staging <- tempfile("spamtpdb-")
  dir.create(staging)
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)
  fixture <- list(ramp_version = "3.0.7")
  saveRDS(fixture, file.path(staging, "ramp_db_metadata.rds"))

  database <- loadSpaMTPDatabase(
    "ramp_db_metadata",
    source = "spamtpdb",
    local_dir = staging,
    offline = TRUE,
    refresh = TRUE,
    verify = FALSE
  )

  expect_identical(database$ramp_db_metadata$ramp_version, "3.0.7")
  expect_identical(
    attr(database$ramp_db_metadata, "spamtp_database")$source,
    "spamtpdb"
  )
})

test_that("custom resource bundles are validated and subset", {
  custom <- list(
    chem_props = data.frame(exactmass = 100),
    pathway = data.frame(pathwayRampId = "RAMP_P_1")
  )

  selected <- SpaMTP:::.spamtp_db_bundle("chem_props", database = custom)
  expect_named(selected, "chem_props")
  expect_equal(selected$chem_props$exactmass, 100)

  expect_error(
    SpaMTP:::.spamtp_db_bundle("source_df", database = custom),
    "missing resource"
  )
})

test_that("database registry reports canonical resource names", {
  skip_if_not_installed("SpaMTPdb")
  registry <- spaMTPDatabaseInfo()
  expect_s3_class(registry, "data.frame")
  expect_true("resource" %in% names(registry))
  expect_true(all(
    c("chem_props", "source_df", "analytehaspathway", "pathway") %in%
      registry$resource
  ))
})

test_that("cache keys follow configured directories and verification requirements", {
  paths <- c(tempfile("spamtpdb-a-"), tempfile("spamtpdb-b-"))
  for (path in paths) dir.create(path)
  on.exit(unlink(paths, recursive = TRUE), add = TRUE)
  for (i in seq_along(paths)) {
    saveRDS(list(label = paste0("synthetic-", i)),
      file.path(paths[i], "ramp_db_metadata.rds"))
  }
  withr::local_options(list(SpaMTPdb.resource_dir = paths[1]))
  first <- loadSpaMTPDatabase("ramp_db_metadata", offline = TRUE, verify = FALSE)
  options(SpaMTPdb.resource_dir = paths[2])
  second <- loadSpaMTPDatabase("ramp_db_metadata", offline = TRUE, verify = FALSE)
  expect_identical(first$ramp_db_metadata$label, "synthetic-1")
  expect_identical(second$ramp_db_metadata$label, "synthetic-2")
  expect_false(attr(second$ramp_db_metadata, "spamtp_database")$verify_local)
  expect_error(loadSpaMTPDatabase("ramp_db_metadata", offline = TRUE), "MD5")
  expect_error(loadSpaMTPDatabase("ramp_db_metadata", offline = TRUE, verify = NA),
    "TRUE or FALSE")
})

test_that("the native companion example works with the software annotation engine", {
  skip_if_not_installed("SpaMTPData", "0.99.1")
  object <- SpaMTPData::spaMTPExampleData()
  before <- intersect(c("Seurat", "SeuratObject"), loadedNamespaces())
  object <- normalizeSMData(object, "LogNormalize", verbose = FALSE)
  object <- runMetabolicPCA(object, npcs = 2, slot = "logcounts")
  expect_true("pca" %in% SingleCellExperiment::reducedDimNames(object))
  database <- data.frame(id = "demo-2HG", name = "2HG", formula = "C5H8O5",
    exactmass = 148.037173366)
  index <- buildMZAnnotationIndex(database, adducts = "M+H")
  result <- annotateSM(object, index = index, adducts = "M+H",
    return.only.annotated = FALSE, verbose = FALSE)
  expect_identical(dim(result), dim(object))
  expect_identical(SummarizedExperiment::rowData(result)$all_IsomerNames[1], "2HG")
  expect_identical(SingleCellExperiment::altExp(result, "transcriptome"),
    SingleCellExperiment::altExp(object, "transcriptome"))
  expect_identical(spaMTPDatabaseInfo(), SpaMTPdb::spaMTPdbResources())
  expect_identical(intersect(c("Seurat", "SeuratObject"), loadedNamespaces()), before)
})
