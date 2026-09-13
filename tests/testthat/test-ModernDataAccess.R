accessorViolations <- function(code) {
  tokens <- getParseData(parse(text = code, keep.source = TRUE))
  calls <- gsub("^`|`$", "", tokens$text[tokens$token == "SYMBOL_FUNCTION_CALL"])
  c(tokens$text[tokens$terminal & tokens$text == "@"],
    intersect(calls, c("slot", "slot<-", "GetAssayData", "SetAssayData",
                       "spatialData", "spatialDataNames")))
}

rVignetteChunks <- function(path) {
  lines <- readLines(path, warn = FALSE)
  chunks <- list()
  code <- character()
  inR <- FALSE
  for (line in lines) {
    if (!inR && grepl("^\\s*`{3,}\\s*\\{r([ ,}]|$)", line, perl = TRUE)) {
      inR <- TRUE
      code <- character()
    } else if (inR && grepl("^\\s*`{3,}\\s*$", line, perl = TRUE)) {
      chunks[[length(chunks) + 1L]] <- code
      inR <- FALSE
    } else if (inR) {
      code <- c(code, line)
    }
  }
  if (inR) stop("Unclosed R chunk in ", path)
  chunks
}

test_that("the accessor audit distinguishes code from comments and SMILES", {
  expect_length(accessorViolations(c(
    '# legacy object@meta.data',
    'smiles <- "O[C@H](CCC(O)=O)C(O)=O"',
    'SeuratObject::Misc(x, slot = "analysis")',
    'SummarizedExperiment::assay(x, "counts")')), 0L)
  expect_equal(accessorViolations('x@meta.data'), "@")
  expect_equal(accessorViolations('methods::slot(x, "assays")'), "slot")
  expect_equal(accessorViolations('methods::`slot<-`(x, "assays", value = y)'), "slot<-")
  expect_equal(accessorViolations('SeuratObject::GetAssayData(x)'), "GetAssayData")
  expect_equal(accessorViolations('SeuratObject::SetAssayData(x, new.data = y)'),
    "SetAssayData")
  expect_equal(accessorViolations('SpatialExperiment::spatialData(x)'), "spatialData")
})

test_that("package functions use accessors rather than internal or superseded APIs", {
  namespace <- asNamespace("SpaMTP")
  problems <- lapply(ls(namespace, all.names = TRUE), function(name) {
    value <- get(name, namespace)
    if (!is.function(value)) return(character())
    violations <- accessorViolations(deparse(body(value)))
    if (length(violations)) paste(name, violations, sep = ": ") else character()
  })
  expect_length(unlist(problems), 0L)
})

test_that("all shipped vignette R chunks use current data-access interfaces", {
  paths <- list.files(c(test_path("..", "..", "vignettes"),
    system.file("doc", package = "SpaMTP")), pattern = "\\.Rmd$", full.names = TRUE)
  skip_if(!length(paths), "Vignette sources are not available in this installation")
  for (path in paths) {
    chunks <- rVignetteChunks(path)
    expect_true(length(chunks) > 0L, info = basename(path))
    violations <- lapply(chunks, accessorViolations)
    expect_true(length(unlist(violations)) == 0L, info = basename(path))
  }
})

test_that("METASPACE import has one authoritative coordinate store", {
  fixture <- list(
    images = list("100.01" = matrix(1:6, nrow = 2),
                  "200.02" = matrix(7:12, nrow = 2)),
    annotations = data.frame(mz = c(100.01, 200.02), formula = c("A", "B")))
  local_mocked_bindings(getMetaspace = function(...) fixture, .package = "SpaMTP")
  object <- suppressMessages(loadMetaspace("demo", verbose = FALSE))
  expect_s4_class(object, "SpatialExperiment")
  expect_equal(unname(SpatialExperiment::spatialCoords(object)),
    unname(cbind(x = rep(1:3, each = 2), y = rep(1:2, 3))))
  expect_equal(unname(as.matrix(SummarizedExperiment::assay(object))),
    unname(rbind(1:6, 7:12)))
  expect_equal(SummarizedExperiment::rowData(object)$mz, c(100.01, 200.02))
  expect_false(any(c("x", "y", "x_coord", "y_coord") %in%
    colnames(SummarizedExperiment::colData(object))))
  expect_identical(object$metaspace_dataset, rep("demo", 6))
  shifted <- object
  SpatialExperiment::spatialCoords(shifted)[, "x"] <-
    SpatialExperiment::spatialCoords(shifted)[, "x"] + 10
  expect_equal(SummarizedExperiment::colData(shifted),
    SummarizedExperiment::colData(object))
})
