
#' Load imzML data into a Bioconductor-native spatial container
#'
#' Raw spectra remain in Cardinal's file-backed `MSImagingArrays`, without
#' asking the reader to infer a shared mass axis.
#' Once a resolution is supplied and Cardinal has
#' produced a common binned feature axis, the default return value is a
#' `SpatialExperiment`. Seurat conversion is explicit and optional.
#'
#' @param file Character string defining the directory path of the file and the file name.
#' @param mass.range Optional two-element mass range for binning. Requires
#'   resolution; raw spectral filtering should be performed in Cardinal.
#' @param resolution Numeric value defining the the accuracy to which the m/z values will be binned after reading. This value can be in either "ppm" or "mz" depending on the units type specified (default = 10).
#' @param units Character string defining the resolution value unit type, either c("ppm", "mz") (default = "ppm")
#' @param verbose Boolean indicating whether to show informative processing messages. If TRUE the message will be show, else the message will be suppressed (default = TRUE)
#' @param assay Label stored in mainExpName(); intensities are in the counts assay.
#' @param multi.run Retained for compatibility; Cardinal run IDs are preserved.
#' @param returnType Output container. `"auto"` returns Cardinal for raw data
#'   and SpatialExperiment after binning. SpatialExperiment output requires resolution.
#'
#' @return A Cardinal `MSImagingArrays`/`MSImagingExperiment`, or a `SpatialExperiment`.
#' @export
#'
#' @examples
#' utils::str(formals(loadSM))
#' # data <-loadSM(name = "run1", folder = "/Documents/SpaMTP_test_data/", mass.range = c(160,1500), resolution = 10, assay = "Spatial")
loadSM <- function(
    file,
    mass.range = NULL,
    resolution = NA,
    units = "ppm",
    verbose = TRUE,
    assay = "Spatial",
    multi.run = FALSE,
    returnType = c("auto", "Cardinal", "SpatialExperiment")
) {
  returnType <- match.arg(returnType)
  if (length(resolution) != 1L || (!is.na(resolution) &&
      (!is.numeric(resolution) || !is.finite(resolution) || resolution <= 0))) {
    stop("resolution must be NA (raw input) or one positive finite number.", call. = FALSE)
  }
  binned <- length(resolution) == 1L && !is.na(resolution)
  if (!binned && (!is.null(mass.range) || identical(returnType, "SpatialExperiment"))) {
    stop("Supply resolution for mass.range filtering or a rectangular output container.",
         call. = FALSE)
  }
  data <- Cardinal::readImzML(
    file = file,
    memory = FALSE,
    as = "MSImagingArrays",
    verbose = verbose
  )
  if (binned) {
    data <- .binCardinalData(
      data = data,
      mass.range = mass.range,
      resolution = resolution,
      units = units,
      method = "sum",
      verbose = verbose
    )
  }
  if (identical(returnType, "auto")) {
    returnType <- if (binned) "SpatialExperiment" else "Cardinal"
  }
  if (identical(returnType, "Cardinal")) {
    return(data)
  }
  data <- cardinalToSpatialExperiment(data)
  SingleCellExperiment::mainExpName(data) <- assay
  data
}



#' Read a binned spatial metabolomics matrix
#'
#' This function reads SM data stored in a table format, with one row per
#' pixel and two columns named 'x' and 'y' storing its coordinates.
#' All other columns should be the relative m/z values containing the intensity values of each pixel.
#' Coordinates are stored only in `spatialCoords()`, not duplicated in
#' pixel metadata where a later alignment could leave them out of date.
#'
#' @param mtx.file Character string defining the path of the spatial metabolomic image matrix .csv file.
#' @param assay Label stored in mainExpName(); intensities are in the counts assay.
#' @param verbose Boolean indicating whether to show the message. If TRUE, the message will be shown; else, it will be suppressed (default = TRUE).
#' @param feature.start.column Numeric value defining the start index containing the x, y, and m/z value columns within the table (default = 1).
#' @param mz.prefix Character string matching the prefix string in front of each m/z name (default = NULL).
#' @param project.name Character string defining the name of the sample to be assigned as `orig.ident` (default = "SpaMTP").
#' @param returnType Output container; only `"SpatialExperiment"` is supported.
#'
#' @details
#' **NOTE:** The input file must be in a format similar to the table below:
#'
#' ```r
#' A data.frame: 5 × 5
#'    x   y   mz1  mz2  mz3
#' 1  0   1    0    0   11
#' 2  0   2    0    0    0
#' 3  0   3    0    0    0
#' 4  0   4   20    0    0
#' 5  0   5    0    0    0
#' ```
#'
#' - The first two columns (`x`, `y`) contain the respective spatial coordinates.
#' - The subsequent columns contain the m/z values and their intensities for each spatial pixel.
#'
#' @return A `SpatialExperiment`.
#' @export
#'
#' @examples
#' utils::str(formals(readSMMatrix))
#' # msi_data <- readSMMatrix("~/Documents/msi_mtx.csv")
readSMMatrix <- function(
    mtx.file,
    assay = "Spatial",
    verbose = TRUE,
    feature.start.column = 1,
    mz.prefix = NULL,
    project.name = "SpaMTP",
    returnType = "SpatialExperiment"
) {
  returnType <- match.arg(returnType)

  verbose_message(message_text = "Reading mtx file.... ", verbose = verbose)

  data <- as.data.frame(data.table::fread(mtx.file))

  if (feature.start.column > 0){
    verbose_message(message_text = paste0("Spliting matrix data from column ", feature.start.column," onwards .... "), verbose = verbose)
    data <- data[,feature.start.column:dim(data)[2]]
  }

  if ("x" %in% colnames(data) && "y" %in% colnames(data)) {
    coords <- data[c("x", "y")]

    data <-  data[, !(names(data) %in% c("x", "y"))]

  } else {
    stop("X and Y Tissue coordinates not found. Expecting columns 'x' and 'y'")
  }

  barcodes <- paste0(coords$x, "_", coords$y)

  rownames(data) <- barcodes
  rownames(coords) <- barcodes

  rawNames <- colnames(data)
  if (!is.null(mz.prefix)) {
    rawNames <- sub(mz.prefix, "", rawNames)
  }
  rawNames <- sub("^mz-", "", rawNames)
  mz <- suppressWarnings(as.numeric(rawNames))
  if (any(!is.finite(mz))) {
    stop("All feature columns must resolve to numeric m/z values.", call. = FALSE)
  }
  featureNames <- make.unique(paste0("mz-", rawNames))

  if(is.null(project.name)){
    project.name <- "SpaMTP"
  }

  intensity <- Matrix::Matrix(t(as.matrix(data)), sparse = TRUE)
  rownames(intensity) <- featureNames
  colnames(intensity) <- barcodes
  featureData <- S4Vectors::DataFrame(
    mz = mz,
    raw_mz = mz,
    mz_names = featureNames,
    row.names = featureNames
  )
  pixelData <- S4Vectors::DataFrame(
    orig.ident = rep(project.name, length(barcodes)),
    row.names = barcodes
  )
  object <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = intensity),
    rowData = featureData,
    colData = pixelData,
    spatialCoords = as.matrix(coords[, c("x", "y"), drop = FALSE]),
    sample_id = project.name
  )
  S4Vectors::metadata(object)$spamtp_assay_name <- assay
  SingleCellExperiment::mainExpName(object) <- assay
  object
}
