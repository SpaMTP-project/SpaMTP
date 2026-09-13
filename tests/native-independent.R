## R CMD check runs this script in a fresh process, independently of the
## optional conversion tests. Neither Seurat namespace may be loaded.
library(SpaMTP)
stopifnot(!any(c("Seurat", "SeuratObject") %in% loadedNamespaces()))
values <- rbind(1:8, c(3, 1, 4, 8, 2, 5, 7, 6), 8:1, rep(c(1, 4), 4))
dimnames(values) <- list(paste0("mz-", 101:104), paste0("p", 1:8))
spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = values),
    rowData = S4Vectors::DataFrame(mz = 101:104),
    colData = S4Vectors::DataFrame(region = rep(c("a", "b"), 4)),
    spatialCoords = cbind(x = 1:8, y = rep(0:1, 4)))
spe <- normalizeSMData(spe, "LogNormalize", verbose = FALSE)
spe <- scaleSMData(spe)
spe <- suppressWarnings(runMetabolicPCA(spe, npcs = 2))
spe <- runSpatialGraphPCA(spe, n_components = 2, n_neighbors = 2,
                           platform = "ST", verbose = FALSE)
spe <- getKmeanClusters(spe, clusters = 2)
spe <- addTranscriptome(spe, values[1:3, ])
spe <- suppressWarnings(multiOmicIntegration(spe,
    dims.list = list(1:2, 1:2), verbose = FALSE))
subset <- subsetSPM(spe, subset = region == "a")
stopifnot(ncol(subset) == 4L,
          ncol(SingleCellExperiment::altExp(subset, "transcriptome")) == 4L,
          "scaled" %in% SummarizedExperiment::assayNames(spe),
          ncol(SingleCellExperiment::reducedDim(spe, "integrated")) == 4L,
          !any(c("Seurat", "SeuratObject") %in% loadedNamespaces()))
