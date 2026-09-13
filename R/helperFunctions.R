#' Helper function for suppressing function progress messages
#'
#' @param message_text Character string containing the message being shown
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the messsage will be suppressed (default = TRUE).
#' @noRd
verbose_message <- function(message_text, verbose) {
  if (verbose) {
    message(message_text)
  }
}





#'@importFrom magrittr %>% %<>%
NULL

#' Subset features and pixels through Bioconductor containers
#'
#' Uses the container's standard bracket method so alternative experiments,
#' reduced dimensions and spatial edges remain synchronized with pixels.
#'
#' Function params/args:
#' @param object A SummarizedExperiment, including SingleCellExperiment or SpatialExperiment.
#' @param subset Logical expression indicating features/variables to keep
#' @param cells A vector of cells to keep; if \code{NULL}, defaults to all cells
#' @param idents A vector of identity classes to keep
#' @param features A vector of feature names or indices to keep
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = FALSE).
#' @param ... Reserved; additional arguments are rejected.
#'
#' @return A subset of the same Bioconductor container class.
#' @export
#'
#' @examples
#' utils::str(formals(subsetSPM))
#' # subsetSPM(spe, subset = region == "edge")
subsetSPM <- function(
    object = NULL,
    subset = NULL,
    cells = NULL,
    idents = NULL,
    features = NULL,
    verbose = FALSE,
    ...)
{
  .requireExperiment(object, "SummarizedExperiment")
  if (length(list(...))) stop("Unused arguments in ...", call. = FALSE)
  keep <- rep(TRUE, ncol(object))
  if (!missing(subset)) {
    selection <- eval(substitute(subset), .cellMetadata(object), parent.frame())
    if (!is.null(selection)) {
      if (!is.logical(selection) || length(selection) != ncol(object)) {
        stop("subset must evaluate to one logical value per pixel.", call. = FALSE)
      }
      keep <- keep & !is.na(selection) & selection
    }
  }
  if (!is.null(idents)) {
    labels <- if (methods::is(object, "SingleCellExperiment")) {
      SingleCellExperiment::colLabels(object)
    } else NULL
    if (is.null(labels)) {
      stop("Set colLabels(object), or use subset with a colData column.", call. = FALSE)
    }
    keep <- keep & labels %in% idents
  }
  indices <- if (is.null(cells)) seq_len(ncol(object)) else {
    if (is.character(cells)) match(cells, colnames(object)) else cells
  }
  if (anyNA(indices) || any(indices < 1 | indices > ncol(object)) ||
      any(indices != floor(indices))) {
    stop("cells must identify existing pixels.", call. = FALSE)
  }
  rows <- if (is.null(features)) seq_len(nrow(object)) else features
  return(object[rows, indices[keep[indices]], drop = FALSE])
}


#' Check Cardinal Version
#'
#' Checks if the Cardinal Package installed is >= version 3.6.0
#'
#' @return Boolean indicating if the version of Cardinal is >= 3.6.0
check_cardinal_version <- function(){
  if (numeric_version(utils::packageVersion("Cardinal")) >= numeric_version("3.6")) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


########################################################################
