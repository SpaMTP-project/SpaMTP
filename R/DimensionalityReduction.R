#' Compute PCA with scater on a selected experiment
#'
#' Uses scater::runPCA() and stores embeddings in reducedDim().
#' The user can provide a bin/resolution size to increase the bin size and reduce the dimensionality/noise of the SM dataset prior to calculating PCAs.
#'
#' @param SpaMTP Bioconductor experiment that contains spatial metabolic information.
#' @param npcs Positive number of PCs, capped by matrix dimensions. NULL uses the variance threshold.
#' @param variance_explained_threshold Cumulative variance fraction in (0, 1], used only when npcs is NULL.
#' @param assay Primary (`main`) or alternative experiment name.
#' @param slot Character string defining the assay slot containing the intensity values (default = "counts").
#' @param show_variance_plot Boolean indicating weather to display the variance plot output by this analysis (default = FALSE).
#' @param bin_resolution Numeric value defining the resolution to use for binning m/z peaks. If set to `NULL`, no binning will be performed (default = NULL).
#' @param resolution_units Character string specifying the units of the `bin_resolution`. Either 'ppm' or 'mz' can be provided. `bin_resolution` must be provided for this parameter to be implemented (default = "ppm").
#' @param bin_method Character string defining the method to use for binning respective m/z peaks that fall within a bin. Options for this parameter can be one of "sum", "mean", "max" or "min". `bin_resolution` must be provided for this parameter to be implemented (default = "sum").
#' @param reduction.name Character string indicating the name associated with the PCA results stored in the output Bioconductor experiment (default = "pca").
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#' @param features Optional exact feature identifiers. When supplied, all and
#'   only these rows enter PCA; the scater top-500 default is bypassed.
#' @param scale Scale features to unit variance before PCA.
#'
#'
#' @return The input with PCA embeddings in reducedDim().
#' @export
#'
#' @examples
#' counts <- matrix(seq_len(30), nrow = 5)
#' sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
#' sce <- runMetabolicPCA(sce, npcs = 2)
#' SingleCellExperiment::reducedDimNames(sce)
#' ## For running PCA on un-adjusted peak bin sizes
#' # spamtp_obj <- runMetabolicPCA(spamtp_obj, npcs = 50)
#'
#' ## For running PCA on increased peak bin sizes
#' # spamtp_obj <- runMetabolicPCA(spamtp_obj, npcs = 50, bin_resolution = 50)
runMetabolicPCA <- function(SpaMTP,
                            npcs = 30,
                            variance_explained_threshold = 0.9,
                            assay = "main",
                            slot = "counts",
                            show_variance_plot= FALSE,
                            bin_resolution = NULL,
                            resolution_units = "ppm",
                            bin_method = "sum",
                            reduction.name = "pca",
                            verbose = TRUE, features = NULL, scale = FALSE)
{
  .requireExperiment(SpaMTP, "SingleCellExperiment")
  analysisData <- .experimentForAssay(SpaMTP, assay)
  if (!is.null(bin_resolution)) {
    if (!identical(analysisData, SpaMTP)) {
      stop("Bin the primary MSI experiment, not an alternative modality.", call. = FALSE)
    }
    analysisData <- binSpaMTP(
      analysisData, resolution = bin_resolution, units = resolution_units,
      method = bin_method, assay = assay, slot = slot)
    slot <- "counts"
  }
  if (!methods::is(analysisData, "SingleCellExperiment")) {
    analysisData <- methods::as(analysisData, "SingleCellExperiment")
  }
  values <- .assayData(analysisData, layer = slot)
  if (!is.null(features)) {
    if (!is.character(features) || anyNA(features) || anyDuplicated(features) ||
        length(features) < 2L || !all(features %in% rownames(values)))
      stop("features must contain distinct, observed feature names.", call. = FALSE)
    values <- values[features, , drop = FALSE]
  }
  maximum <- min(dim(values)) - 1L
  if (maximum < 1L || any(!is.finite(values))) {
    stop("PCA requires finite values and at least two features and pixels.", call. = FALSE)
  }
  if (!is.null(npcs)) {
    .validateFeatureCount(npcs, "npcs")
  } else if (length(variance_explained_threshold) != 1L ||
             !is.finite(variance_explained_threshold) ||
             variance_explained_threshold <= 0 ||
             variance_explained_threshold > 1) {
    stop("variance_explained_threshold must be in (0, 1] when npcs is NULL.",
         call. = FALSE)
  }
  components <- if (is.null(npcs)) maximum else min(npcs, maximum)
  analysisData <- scater::runPCA(
    analysisData, exprs_values = slot,
    ncomponents = components, subset_row = features, scale = scale,
    BSPARAM = if (components >= min(dim(values)) / 2) BiocSingular::ExactParam() else BiocSingular::IrlbaParam(),
    name = reduction.name)
  embedding <- SingleCellExperiment::reducedDim(analysisData, reduction.name)
  percentVar <- attr(embedding, "percentVar")
  if (is.null(npcs)) {
    retained <- which(cumsum(percentVar) / 100 >= variance_explained_threshold)[1L]
    if (is.na(retained)) retained <- ncol(embedding)
    loadings <- attr(embedding, "rotation")
    embedding <- embedding[, seq_len(retained), drop = FALSE]
    attr(embedding, "percentVar") <- percentVar[seq_len(retained)]
    if (!is.null(loadings)) {
      attr(embedding, "rotation") <- loadings[, seq_len(retained), drop = FALSE]
    }
  }
  if (isTRUE(show_variance_plot)) {
    varianceData <- data.frame(
      component = seq_along(percentVar), cumulative = cumsum(percentVar) / 100)
    print(ggplot2::ggplot(varianceData,
                         ggplot2::aes(.data$component, .data$cumulative)) +
      ggplot2::geom_line() + ggplot2::geom_point() +
      ggplot2::labs(x = "Principal component", y = "Cumulative variance explained"))
  }
  SingleCellExperiment::reducedDim(SpaMTP, reduction.name) <- embedding
  S4Vectors::metadata(SpaMTP)$reduction_inputs[[reduction.name]] <- list(
    assay = assay, layer = slot, features = rownames(attr(embedding, "rotation")),
    scale = scale, method = "scater::runPCA")
  SpaMTP
}






#' Perform Dimensionality Reduction using Graph-Regularised PCA on Spatial Data
#'
#' Computes a graph-regularised PCA using spatial coordinates and scaled expression data. A k-nearest neighbour (k-NN) graph is computed using spatial locations and used to regularise the PCA decomposition via a graph Laplacian. The result is stored in reducedDim(), with spatial edges in colPairs().
#'
#' Note: This method has been adapted from the GraphPCA Python package
#' (\doi{10.1186/s13059-024-03429-x}).
#'
#' @param data A Bioconductor experiment containing spatial data (feature data and spatial coordinates).
#' @param n_components Integer specifying the number of principal components to compute (default = 50).
#' @param assay Primary (`main`) or alternative experiment name.
#' @param slot Character string defining the name of the slot to extract scaled data from (default = "scaled").
#' @param image Reserved; must be NULL. Subset pixels before running graph PCA.
#' @param platform Character string matching either `"Visium"` or `"ST"` to determine how the k-NN graph is constructed. If "Visium" k-nns will handle the hexagon spot arrangement, including setting `n_neighbors` = 6, else "ST" assignment will set `n_neighbors` = 4 unless a value is specifically provided (default = "Visium").
#' @param lambda Numeric value defining the regularisation parameter that controls the influence of the graph Laplacian (default = 0.5).
#' @param n_neighbors Integer value specifying the number of spatial neighbours to use. If `NULL`, will default of 6 for "Visium" data and 4 for "ST" platforms (default = NULL).
#' @param include_self Boolean logical value indicating whether to include self-connections in the graph (default = FALSE).
#' @param alg Character string specifying the algorithm to use for nearest neighbour search (passed to `FNN::get.knn()`) (default ="kd_tree").
#' @param fast Boolean logical value stating whether to use fast approximate eigendecomposition via `RSpectra::eigs_sym()`. For large datasets this is recommended (default =TRUE).
#' @param graph_name Character string used to store the computed spatial graph (default ="SpatialKNN").
#' @param reduction_name Character string used to store the dimensionality reduction (default ="SpatialPCA").
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#'
#' @return A Bioconductor experiment with a new graph and spatially-aware PCA reduction.
#' @export
#'
#' @rawNamespace import(Matrix, except = c(expand, head, pack, unpack))
#'
#' @examples
#' utils::str(formals(runSpatialGraphPCA))
#' # spamtp_obj <- runSpatialGraphPCA(spamtp_obj, platform = "Visium")
runSpatialGraphPCA <- function(data, n_components=50, assay = "main", slot = "scaled", image = NULL, platform="Visium", lambda=0.5, n_neighbors=NULL, include_self = FALSE, alg = "kd_tree", fast = TRUE, graph_name = "SpatialKNN", reduction_name = "SpatialPCA", verbose = TRUE){
  data <- .nativeSpatialObject(data)
  if (!is.null(image)) stop("Subset pixels first; image must be NULL.", call. = FALSE)


  if(!platform %in% c("Visium", "ST")){
    stop("Incorrect value for plantform! platform must be assigned either 'Visium' or 'ST'... If data is not Visium 'spot'-based then set platform to 'ST'")
  }

  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda < 0){
    stop("lambda must be one finite, non-negative number.", call. = FALSE)
  }
  .validateFeatureCount(n_components, "n_components")

  assayMatrix <- .assayData(data, assay, slot)
  if (min(dim(assayMatrix)) < 2L || any(!is.finite(assayMatrix))) {
    stop("Graph PCA needs finite values and at least two features and pixels.",
         call. = FALSE)
  }
  n_components <- min(n_components, nrow(assayMatrix), ncol(assayMatrix))
  Expr = t(assayMatrix)


  if(is.null(n_neighbors)){

    if (platform == "Visium"){
      n_neighbors <- 6
      verbose_message(message_text = "No nearest neighbour value provided. Based on 'Visium' platform type setting n-neighbours = 6 ..." , verbose = verbose)
    }else{
      verbose_message(message_text = "No nearest neighbour value provided. Based on 'ST' platform type setting n-neighbours = 4 ..." , verbose = verbose)
      n_neighbors <- 4
    }
  } else {
    n_neighbors <- as.integer(n_neighbors)
    verbose_message(message_text = paste0("Using a nearest neighbour value = ",n_neighbors,"! NOTE: the recomended values are -> 'Visium' data: n_neighbors = 6  / 'ST' data: n_neighbors = 4 ...." ), verbose = verbose)
  }

  location <- .nativeCoordinates(data)
  sampleIndices <- split(seq_len(nrow(location)), location$sample_id)
  graphs <- lapply(sampleIndices, function(index) {
    kNeighborsGraph(location[index, c("x", "y"), drop = FALSE],
      n_neighbors = n_neighbors, platform = platform,
      include_self = include_self, alg = alg)
  })
  graph <- Matrix::bdiag(graphs)
  originalOrder <- order(unlist(sampleIndices, use.names = FALSE))
  graph <- graph[originalOrder, originalOrder, drop = FALSE]
  graph <- 0.5 * (graph + t(graph))


  # Keep the Laplacian sparse: the previous conversion allocated an n-by-n
  # dense matrix and made full-tissue analysis unnecessarily expensive.
  graphL <- Matrix::Diagonal(x = Matrix::rowSums(graph)) - graph


  # Create identity matrix and add lambda * graphL
  n <- nrow(Expr)

  G <- Matrix::Diagonal(n) + (lambda * graphL)

  X <- Matrix::solve(G, as.matrix(Expr))


  rownames(graph) <- colnames(graph) <- colnames(assayMatrix)

  spatialGraphs <- S4Vectors::metadata(data)$spatialGraphs %||% list()
  spatialGraphs[[graph_name]] <- graph
  S4Vectors::metadata(data)$spatialGraphs <- spatialGraphs
  edges <- summary(methods::as(graph, "generalMatrix"))
  SingleCellExperiment::colPair(data, graph_name) <- S4Vectors::SelfHits(
    edges$i, edges$j, nnode = ncol(data), weight = edges$x)


  rm(G)
  rm(graph)
  rm(graphL)

  # 3. Eigendecomposition of C
  if (fast && n_components < ncol(Expr)){
    eig <- RSpectra::eigs_sym(crossprod(Expr, X), k = n_components)
  } else {
    eig <- eigen(t(Expr) %*% X, symmetric = TRUE)
  }

  W <- eig$vectors[, seq_len(n_components), drop = FALSE]
  Z <- X %*% W

  rownames(Z) <- colnames(assayMatrix)      # cells
  rownames(W) <- rownames(assayMatrix)      # genes (or features)
  colnames(Z) <- paste0("PC_", 1:ncol(Z))
  colnames(W) <- paste0("PC_", 1:ncol(W))

  # Create the PCA reduction
  SingleCellExperiment::reducedDim(data, reduction_name) <- as.matrix(Z)
  loadings <- S4Vectors::metadata(data)$reductionLoadings %||% list()
  loadings[[reduction_name]] <- W
  S4Vectors::metadata(data)$reductionLoadings <- loadings
  S4Vectors::metadata(data)$reduction_inputs[[reduction_name]] <- list(
    assay = assay, layer = slot, features = rownames(assayMatrix),
    observations = colnames(assayMatrix), lambda = lambda, neighbors = n_neighbors,
    graph = graph_name, sample_separated = TRUE, method = "runSpatialGraphPCA")


  return(data)

}


#' Construct a k-Nearest Neighbour Graph from Spatial Coordinates.
#'
#' Builds a sparse adjacency matrix representing the k-nearest neighbour relationships between spots or cells, using Euclidean distance between spatial coordinates.
#'
#' @param location A data frame or matrix with columns `x` and `y` representing spatial coordinates, and rownames match barcode/cell names.
#' @param n_neighbors Integer value defining the number of neighbours to consider.
#' @param platform Character string being either `"Visium"` or `"ST"`, specifying the platform of the dataset being analysed. Determines how distance thresholds and adjacency are calculated whereby "Visium" data is handled slightly differently to compensate for the hexagon shape of the spots.
#' @param include_self Boolean logical value stating whether to include self-connections (default = FALSE).
#' @param alg Character string specifying the algorithm to use for nearest neighbour search. Passed to `FNN::get.knn()` (default = "kd_tree").
#'
#' @return A sparse adjacency matrix (`dgCMatrix`) representing the k-NN graph.
#' @export
#'
#' @examples
#' utils::str(formals(kNeighborsGraph))
#' ### HELPER FUNCTION
kNeighborsGraph <- function(location, n_neighbors, platform, include_self = FALSE, alg = "kd_tree") {
  .validateFeatureCount(n_neighbors, "n_neighbors")
  if (any(!is.finite(as.matrix(location)))) {
    stop("Spatial coordinates must be finite.", call. = FALSE)
  }
  if (nrow(location) <= 1L) {
    return(Matrix::Diagonal(nrow(location), x = as.numeric(include_self)))
  }
  n_neighbors <- min(n_neighbors, nrow(location) - 1L)
  # Get the k-nearest neighbors using kd_tree (most similar to scikit-learn's default)
  knn_result <- FNN::get.knn(location, k = n_neighbors, algorithm = alg)

  if(platform == "Visium"){

    # Step 1: Flatten the distances and indices
    N <- nrow(knn_result$nn.index)

    # Reshape distances and column indices
    dists <- as.vector(t(knn_result$nn.dist))
    col_indices <- as.vector(t(knn_result$nn.index))

    # Create row indices (equivalent to np.repeat())
    row_indices <- rep(1:N, each = n_neighbors)

    # Apply neighbor correction if needed

    dist_cutoff <- median(dists) * 1.3  # Small amount of sway
    mask <- dists < dist_cutoff

    # Filter using the mask
    row_indices <- row_indices[mask]
    col_indices <- col_indices[mask]
    dists <- dists[mask]

    adjacency <- sparseMatrix(
      i = row_indices,
      j = col_indices,
      x = rep(1, length(row_indices)),
      dims = c(N, N),
      repr = "C"  # Ensures a "dgCMatrix" (like CSR in Python)
    )
  } else {

    # Create empty adjacency matrix
    n <- nrow(location)
    adjacency <- Matrix::sparseMatrix(i = rep(seq_len(n), each = n_neighbors),
      j = as.vector(t(knn_result$nn.index)), x = 1, dims = c(n, n))


  }

  # Include self-loops if requested (in this case, no)
  if (include_self) {
    diag(adjacency) <- 1
  }


  return(adjacency)
}



#' Perform K-means clustering on a specified reduction
#'
#' This function runs K-means clustering on a specified reduction in a Bioconductor experiment and adds the cluster assignments to the object metadata.
#'
#' @param data A Bioconductor experiment containing the results from `runSpatialGraphPCA()`.
#' @param reduction Character string stating the name of the reduction slot to use (default = "SpatialPCA").
#' @param cluster.name Character string of the name of the metadata column to store the cluster labels (default = "spatial_clusters").
#' @param clusters Integer defining the number of clusters to form (default = 8).
#' @param iter.max Integer defining the maximum number of iterations allowed (default = 10).
#' @param nstart Integer stating the number of random sets to choose (default = 1).
#' @param algorithm Character string defining the K-means algorithm to use. One of `"Hartigan-Wong"`, `"Lloyd"`, `"Forgy"`, or `"MacQueen"` (default = "Hartigan-Wong").
#' @param trace Logical boolean indicating whether to produce tracing information on the progress of the algorithm (default = FALSE).
#' @param seed Integer of the random seed to use for reproducibility (default = 888).
#'
#' @return A Bioconductor experiment with a new metadata column containing the K-means cluster assignments.
#' @export
#'
#' @examples
#' utils::str(formals(getKmeanClusters))
getKmeanClusters <- function(data, reduction = "SpatialPCA", cluster.name = "spatial_clusters", clusters = 8, iter.max = 10, nstart = 1, algorithm = c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"), trace = FALSE, seed = 888){
  .requireExperiment(data, "SingleCellExperiment")


  if (!reduction %in% SingleCellExperiment::reducedDimNames(data)) {
    stop("Reduction `", reduction, "` is not present in the object.", call. = FALSE)
  }
  embedding <- SingleCellExperiment::reducedDim(data, reduction)

  res <- withr::with_seed(
    seed,
    stats::kmeans(
      embedding,
      centers = clusters,
      iter.max = iter.max,
      nstart = nstart,
      algorithm = algorithm,
      trace = trace
    )
  )

  metadata <- .cellMetadata(data)
  metadata[[cluster.name]] <- res$cluster
  data <- .setCellMetadata(data, metadata)
  return(data)
}
