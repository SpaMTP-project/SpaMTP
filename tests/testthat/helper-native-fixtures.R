nativeFixture <- function(sampleIds = rep("section1", 6)) {
  counts <- rbind(a = 1:6, b = c(2, 1, 5, 3, 6, 4), constant = rep(2, 6))
  colnames(counts) <- paste0("p", seq_len(ncol(counts)))
  object <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts, logcounts = log1p(counts)),
    rowData = S4Vectors::DataFrame(mz = c(100, 200, 300), annotation = letters[1:3]),
    colData = S4Vectors::DataFrame(region = rep(c("edge", "core"), 3)),
    spatialCoords = cbind(x = 1:6, y = c(0, 1, 0, 1, 0, 1)),
    sample_id = sampleIds)
  transcriptome <- rbind(a = 6:1, gene2 = c(0, 1, 0, 1, 0, 1))
  colnames(transcriptome) <- colnames(object)
  addTranscriptome(object, transcriptome)
}
