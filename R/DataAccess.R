#' Access SpaMTP data through Bioconductor containers
#'
#' SpaMTP uses public S4 accessors, not internal slots. A SpatialExperiment is
#' a SingleCellExperiment and a SummarizedExperiment. The required container
#' class is stated on each function's help page; spatial operations need
#' coordinates, while paired modalities and reductions need SingleCellExperiment
#' infrastructure. Raw, unaligned MSI spectra remain in Cardinal until binned
#' or aligned to a common feature-by-pixel matrix.
#'
#' @section Experiments and expression assays:
#' In analysis functions, `assay` selects the primary experiment (`"main"`)
#' or a named `SingleCellExperiment::altExp()`. The compatibility arguments
#' `slot`, `slots`, `SM.slot`, `ST.slot`, `SM_slot`, `ST_slot`, `mz.slot` and
#' `pathway.slot` select expression assay names within those experiments.
#' They do not name S4 slots. Supply exact, non-empty character names, not
#' numeric matrix positions. For example, `assay = "transcriptome"` and
#' `slot = "logcounts"` read
#' `SummarizedExperiment::assay(SingleCellExperiment::altExp(x,
#' "transcriptome"), "logcounts")`.
#'
#' The expression assays within one experiment share feature and pixel axes;
#' `counts`, `normcounts`, `logcounts` and `scaled` represent different value
#' scales, not different modalities. Each alternative experiment has its own
#' feature rows and the same paired pixel columns. Inspect `assayNames()` and
#' `altExpNames()` before selecting data. No automatic alias maps `data` to
#' `logcounts` or `scale.data` to `scaled` during matrix access.
#'
#' For older analysis calls, `NULL`, `"primary"`, `"Spatial"`, `"SPM"`, the
#' primary experiment's `mainExpName()`, and its expression assay names can
#' also select the primary experiment. An exact alternative-experiment name
#' takes precedence over these aliases. Prefer `"main"` and explicit
#' `altExpNames()` in new code. An expression name supplied as `assay` does
#' not replace the separate `slot` selection. Constructors and explicit
#' converters have their own documented naming arguments.
#'
#' @section Metadata and spatial infrastructure:
#' Use `SummarizedExperiment::rowData()` for feature annotations and
#' `SummarizedExperiment::colData()` for pixel metadata. A primary experiment
#' and each `altExp()` have separate `rowData()` and `S4Vectors::metadata()`.
#' Parent pixel metadata are not automatically copied into an extracted
#' alternative experiment. Use `SingleCellExperiment::colLabels()` for groups,
#' `reducedDim()` for embeddings and `colPair()` for pixel-pair graphs.
#' SpatialExperiment coordinates and images are accessed with
#' `spatialCoords()`, `imgData()` and `imgRaster()`.
#'
#' Accessors have replacement forms. Standard column subsetting keeps paired
#' alternative experiments and spatial coordinates aligned. Row subsetting
#' changes only the selected experiment's feature space. SpaMTP does not
#' require users to define additional S4 slots or inspect internal storage.
#'
#' @section Optional Seurat conversion:
#' Core analysis does not use Seurat. SeuratObject is an optional dependency
#' for [seuratToSingleCellExperiment()], [seuratToSpatialExperiment()] and
#' [spatialExperimentToSeurat()]. Convert once at the boundary, then use the
#' native Bioconductor workflow. Seurat assay names become experiment names;
#' selected Seurat layers become expression assays as documented by the converter.
#'
#' @name experimentAccess
#' @seealso [asSpatialExperiment()], [addTranscriptome()], [normalizeSMData()],
#'   [scaleSMData()], [runMetabolicPCA()]
#' @examples
#' x <- SpatialExperiment::SpatialExperiment(
#'     assays = list(counts = rbind(a = 1:4, b = 4:1)),
#'     spatialCoords = cbind(x = 1:4, y = 0),
#'     colData = S4Vectors::DataFrame(row.names = paste0("p", 1:4)))
#' rna <- matrix(c(1, 3, 2, 4), nrow = 1,
#'     dimnames = list("gene1", colnames(x)))
#' x <- addTranscriptome(x, rna)
#' x <- normalizeSMData(x, "LogNormalize", assay = "transcriptome",
#'     slot = "counts", verbose = FALSE)
#' rna <- SingleCellExperiment::altExp(x, "transcriptome")
#' SummarizedExperiment::assayNames(rna)
#' SummarizedExperiment::assay(rna, "logcounts")
#' SummarizedExperiment::rowData(rna)$feature_label <- "example gene"
#' SingleCellExperiment::altExp(x, "transcriptome") <- rna
#' SummarizedExperiment::colData(x)$region <- c("A", "A", "B", "B")
#' selected <- x[, x$region == "A"]
#' stopifnot(identical(colnames(SingleCellExperiment::altExp(selected,
#'     "transcriptome")), colnames(selected)))
NULL
