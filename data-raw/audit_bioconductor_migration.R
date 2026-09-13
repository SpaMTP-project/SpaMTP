# Run from the package root with Rscript data-raw/audit_bioconductor_migration.R.
# Static evidence only: namespace-call counts are not runtime dependency tests.

parseSources <- function(sources) {
  lapply(sources, function(source) {
    getParseData(parse(text = source, keep.source = TRUE))
  })
}

summarizeSources <- function(parsed, namespace) {
  tokens <- do.call(rbind, parsed)
  functions <- gsub("^`|`$", "", tokens$text[tokens$token == "SYMBOL_FUNCTION_CALL"])
  exports <- sub("^export\\((.*)\\)$", "\\1",
                 grep("^export\\(", namespace, value = TRUE))
  packages <- tokens$text[tokens$token == "SYMBOL_PACKAGE"]
  c(exports = length(exports),
    nonCamelExports = sum(!grepl("^[a-z][A-Za-z0-9]*$", exports)),
    directSlotOperators = sum(tokens$text == "@" & tokens$terminal),
    directSlotCalls = sum(functions %in% c("slot", "slot<-")),
    supersededDataCalls = sum(functions %in% c("GetAssayData", "SetAssayData",
                                               "spatialData", "spatialDataNames")),
    s4Generics = sum(functions == "setGeneric"),
    s4Methods = sum(functions == "setMethod"),
    registeredCoercions = sum(functions == "setAs"),
    qualifiedSeuratCalls = sum(packages %in% c("Seurat", "SeuratObject")))
}

files <- list.files("R", pattern = "\\.[Rr]$", full.names = TRUE)
current <- parseSources(lapply(files, readLines, warn = FALSE))
report <- list(workingTree = summarizeSources(current, readLines("NAMESPACE")))

if (nzchar(Sys.which("git"))) {
  tracked <- system2("git", c("ls-tree", "-r", "--name-only", "HEAD", "R"),
                     stdout = TRUE)
  tracked <- tracked[grepl("\\.[Rr]$", tracked)]
  if (length(tracked)) {
    previous <- parseSources(lapply(tracked, function(path) {
      system2("git", c("show", paste0("HEAD:", path)), stdout = TRUE)
    }))
    namespace <- system2("git", c("show", "HEAD:NAMESPACE"), stdout = TRUE)
    report <- c(list(HEAD = summarizeSources(previous, namespace)), report)
  }
}
print(do.call(cbind, report))

description <- read.dcf("DESCRIPTION")
dependencies <- function(field) {
  trimws(gsub("\\s*\\([^)]*\\)", "",
                strsplit(description[1, field], ",")[[1]]))
}
cat("\nDirect Imports:", length(dependencies("Imports")), "\n")
cat("Seurat packages in Imports:",
    any(c("Seurat", "SeuratObject") %in% dependencies("Imports")), "\n")
cat("SeuratObject is optional:", "SeuratObject" %in% dependencies("Suggests"), "\n")
cat("Seurat analysis package is declared:",
    "Seurat" %in% c(dependencies("Imports"), dependencies("Suggests")), "\n")
perFile <- vapply(current, function(tokens) {
  sum(tokens$token == "SYMBOL_PACKAGE" & tokens$text %in% c("Seurat", "SeuratObject"))
}, integer(1))
print(data.frame(file = files[perFile > 0L], calls = perFile[perFile > 0L]))
stopifnot(all(perFile[basename(files) != "ConvertingBetweenObjects.R"] == 0L))

# Parse actual R chunks (including eval=FALSE) rather than counting @ in text:
# SMILES chirality, roxygen tags and comments are not S4 slot operations.
vignettes <- list.files("vignettes", pattern = "\\.Rmd$", full.names = TRUE)
vignetteReport <- lapply(vignettes, function(path) {
  lines <- readLines(path, warn = FALSE)
  code <- character()
  chunks <- list()
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
  summary <- summarizeSources(parseSources(chunks), character())
  data.frame(file = path, chunks = length(chunks),
    directSlotOperators = summary[["directSlotOperators"]],
    directSlotCalls = summary[["directSlotCalls"]],
    supersededDataCalls = summary[["supersededDataCalls"]])
})
if (length(vignetteReport)) {
  vignetteReport <- do.call(rbind, vignetteReport)
  print(vignetteReport, row.names = FALSE)
  stopifnot(all(vignetteReport$directSlotOperators == 0L),
    all(vignetteReport$directSlotCalls == 0L),
    all(vignetteReport$supersededDataCalls == 0L))
}
stopifnot(all(report$workingTree[c("directSlotOperators", "directSlotCalls",
                                  "supersededDataCalls")] == 0L))
