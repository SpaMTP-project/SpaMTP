
#### SpaMTP Saving Data Objects ########################################################################################################################################################################################

#' Saves SpaMTP Object
#'
#' This function saves a Bioconductor experiment into a standard single-cell/spatial file format.
#' This includes a filtered_feature_bc_matrix folder containing files storing the features, barcode/pixels and intensity matrix.
#' Metadata and sapatial files (such as scale factors and hires/lowres images) are also stored.
#'
#' @param data A Spatial Metabolomic Bioconductor experiment being saved.
#' @param outdir Character string of the directory to save the mtx.mtx, barcode.tsv, features.tsv, barcode_metadata.csv and feature_metadata.csv in.
#' @param assay Character string defining the primary or alternative experiment that contains the m/z count data (default = "Spatial").
#' @param slot Character string defining the primary or alternative experiment slot that contains the m/z values directly (default = "counts").
#' @param image Image ID identifying exactly one imgData row; NULL omits image export.
#' @param annotations Boolean values defining if the Bioconductor experiment contains annotations to be saved (default = FALSE).
#' @param generate.h5 Boolean value indicating whether to generate a filtered_feature_bc_matrix.h5 file. Often used by data loading functions (e.g. scanpy.load_visium). If `FALSE`, only a filtered_feature_bc_matrix folder will be generated (default = TRUE).
#' @param verbose Boolean indicating whether to show informative processing messages. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#'
#' ### Details
#' * This can be used for saving data for transfer to python. Can be read in as Anndata using scanpy.read_10x_mtx().
#' * For saving in R saveRDS() is recommended.
#'
#' @return The output-directory path, invisibly.
#'
#' @export
#'
#' @examples
#' utils::str(formals(saveSpaMTPData))
saveSpaMTPData <- function(data, outdir, assay = "Spatial", slot = "counts", image = NULL, annotations = FALSE, generate.h5 = TRUE, verbose = TRUE){
  .requireExperiment(data)
  if (!is.null(image)) .requireExperiment(data, "SpatialExperiment")


  if (!dir.exists(outdir)) {
    verbose_message(message_text = paste0("Generating new directory to store output here: ", outdir), verbose = verbose)
    dir.create(outdir)
  } else {
    verbose_message(message_text = paste0("Directory already exists, storing output here: ", outdir), verbose = verbose)
  }

  verbose_message(message_text = paste0("Writing ", slot," slot to matrix.mtx, barcode.tsv, genes.tsv"), verbose = verbose)
  assayMatrix <- .assayData(data, assay, slot)
  if (!inherits(assayMatrix, "sparseMatrix")) {
    assayMatrix <- Matrix::Matrix(assayMatrix, sparse = TRUE)
  }
  DropletUtils::write10xCounts(assayMatrix, path = paste0(outdir,"/filtered_feature_bc_matrix/"), overwrite = TRUE)

  verbose_message(message_text = "Writing cell metadata to metadata.csv", verbose = verbose)
  data.table::fwrite(.cellMetadata(data), paste0(outdir,"/barcode_metadata.csv"))

  if(generate.h5){
    DropletUtils::write10xCounts(assayMatrix, path = paste0(outdir,"/filtered_feature_bc_matrix.h5"), type = "HDF5", overwrite = TRUE)
  }


  if (!is.null(image)){

    verbose_message(message_text ="Generating 'spatial' directory ... ", verbose = verbose)
    dir.create(paste0(outdir, "/spatial/"))

    imageData <- as.data.frame(SpatialExperiment::imgData(data))
    selected <- which(imageData$image_id == image)
    if (length(selected) != 1L) {
      stop("`image` must identify exactly one row of imgData(data).", call. = FALSE)
    }
    sampleId <- imageData$sample_id[[selected]]
    scaleFactor <- imageData$scaleFactor[[selected]]
    scaleFactors <- list(
      tissue_hires_scalef = scaleFactor,
      tissue_lowres_scalef = scaleFactor
    )
    sfJSON <- jsonlite::toJSON(
      rapply(
        scaleFactors,
        function(x) if (length(x) == 1L) jsonlite::unbox(x) else x,
        how = "replace"
      )
    )
    write(sfJSON, file = paste0(outdir, "/spatial/scalefactors_json.json"))

    imagePath <- paste0(outdir, "/spatial/tissue_lowres_image.png")
    source <- SpatialExperiment::imgSource(
      data,
      sample_id = sampleId,
      image_id = image
    )
    copied <- length(source) == 1L && !is.na(source) && file.exists(source) &&
      file.copy(source, imagePath, overwrite = TRUE)
    if (!isTRUE(copied)) {
      raster <- SpatialExperiment::imgRaster(
        data,
        sample_id = sampleId,
        image_id = image
      )
      magick::image_write(magick::image_read(raster), imagePath, format = "png")
    }

    keep <- as.character(data$sample_id) == sampleId
    coords <- as.data.frame(SpatialExperiment::spatialCoords(data)[keep, , drop = FALSE])
    coords$cell <- colnames(data)[keep]
    coords$in_tissue <- 1L
    coords$arrayrow <- match(coords$x, sort(unique(coords$x)))
    coords$arraycol <- match(coords$y, sort(unique(coords$y)))
    coords <- coords[c("cell", "in_tissue", "arrayrow", "arraycol", "x", "y")]
    colnames(coords) <- NULL
    data.table::fwrite(
      coords,
      paste0(outdir, "/spatial/tissue_positions_list.csv"),
      col.names = FALSE
    )

  }
  if (annotations){
    verbose_message(message_text = "Writing feature metadata annotations to feature_metadata.csv", verbose = verbose)
    data.table::fwrite(.featureMetadata(data, assay), paste0(outdir,"/feature_metadata.csv"))
  }

  invisible(normalizePath(outdir, mustWork = FALSE))

}

########################################################################################################################################################################################################################
