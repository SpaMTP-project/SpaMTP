
#### SpaMTP Differential Metabolite Analysis Functions ########################################################################################################################################################################################



#' Pools Bioconductor experiment into random pools for pseudo-bulking.
#'
#' Runs pooling of a SpaMTP dataset to generate pseudo-replicates for each unique identity provided.
#' This function is used by `findAllDEMs()`.
#' These random pools are technical partitions, not independent biological
#' replicates; resulting p-values must not be interpreted as biological
#' replication. For population inference aggregate independent samples instead.
#'
#' @param data.filt A Bioconductor experiment containing count values for pooling.
#' @param idents A character string defining the idents column to pool the data against.
#' @param n An integer defining the amount of pseudo-replicates to generate for each sample (default = 3).
#' @param assay Character string defining the assay where the mz count data and annotations are stored (default = "Spatial").
#' @param slot Character string defining the assay storage slot to pull the relative mz intensity values from (default = "counts").
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#' @param seed Numeric value used to set the seed for reproducible randomisation (default = 1234).
#'
#' @returns A SingleCellExperiment object which contains pooled (n)-pseudo-replicate counts data based on the Bioconductor experiment input
#' @export
#'
#' @examples
#' utils::str(formals(runPooling))
runPooling <- function(data.filt, idents, n, assay, slot, seed = 1234, verbose = TRUE) {
  .validateFeatureCount(n, "n")
  cell_metadata <- .cellMetadata(data.filt)
  if (length(idents) != 1L || !idents %in% colnames(cell_metadata) ||
      anyNA(cell_metadata[[idents]])) {
    stop("idents must name a colData column without missing values.", call. = FALSE)
  }
  samples <- unique(cell_metadata[[idents]])
  verbose_message(paste("Partitioning each group into", n, "technical pools."),
                  verbose = verbose)
  for(i in seq_along(samples)){
    wo<-which(cell_metadata[[idents]]== samples[i])
    if (length(wo) < n) {
      stop("Every group must have at least n pixels.", call. = FALSE)
    }
    pooled_ids <- withr::with_seed(
      seed + i,
      sample(rep(seq_len(n), length.out = length(wo)))
    )
    cell_metadata[wo,'orig.ident2']<-paste0("group", i, "_pool", pooled_ids)
  }
  expression <- .assayData(data.filt, assay, slot)
  if (!nrow(expression) || !ncol(expression) || any(!is.finite(expression)) ||
      any(expression < 0)) {
    stop("Pooling needs non-empty, finite, non-negative intensities.", call. = FALSE)
  }
  # edgeR materializes its input. Bound each call to about one million values
  # so sparse MSI input is not densified in its entirety.
  blockSize <- max(1L, floor(1e6 / ncol(expression)))
  blocks <- split(seq_len(nrow(expression)),
                  ceiling(seq_len(nrow(expression)) / blockSize))
  pooledCounts <- do.call(rbind, lapply(blocks, function(rows) {
    edgeR::sumTechReps(expression[rows, , drop = FALSE], ID = cell_metadata$orig.ident2)
  }))
  pools <- colnames(pooledCounts)
  pooledMetadata <- cell_metadata[match(pools, cell_metadata$orig.ident2), , drop = FALSE]
  for (column in colnames(cell_metadata)) {
    for (index in seq_along(pools)) {
      values <- cell_metadata[[column]][cell_metadata$orig.ident2 == pools[index]]
      if (length(unique(values)) != 1L) pooledMetadata[[column]][index] <- NA
    }
  }
  rownames(pooledMetadata) <- pools
  pooledMetadata$ncells <- as.integer(table(factor(cell_metadata$orig.ident2, levels = pools)))
  summed <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = pooledCounts),
    rowData = S4Vectors::DataFrame(.featureMetadata(data.filt, assay)),
    colData = S4Vectors::DataFrame(pooledMetadata))
  S4Vectors::metadata(summed)$spamtp_pooling <- list(
    group = idents, pools_per_group = n, seed = seed,
    replication = "technical partitions, not biological replicates")
  return(summed)
}




#' Runs EdgeR analysis for pooled data
#'
#' Worker function for calculating differentially abundant metabolites per pooling group.
#' This function is used by by `findAllDEMs()`.
#'
#' @param pooled_data A SingleCellExperiment object which contains the pooled pseudo-replicate data.
#' @param data A Bioconductor experiment containing the merged Xenium data being analysed (this is subset).
#' @param ident A character string defining the ident column to perform differential expression analysis against.
#' @param output_dir A character string defining the ident column to perform differential expression analysis against.
#' @param run_name A character string defining the title of this DE analysis (will be used when saving DEMs to .csv file).
#' @param n An integer that defines the number of pseudo-replicates per sample (default = 3).
#' @param logFC_threshold A numeric value indicating the logFC threshold to use for defining significant genes (default = 1.2).
#' @param annotation.column Character string defining the column where annotation information is stored in the assay metadata. This requires annotateSM() to be run where the default column to store annotations is "all_IsomerNames" (default = "None").
#' @param assay A character string defining the assay where the mz count data and annotations are stored (default = "Spatial").
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#' @param return.individual Boolean value defining whether to return a list of individual edgeR objects for each designated ident. If FALSE, one merged edgeR object will be returned (default = FALSE).
#'
#' @returns A modified edgeR object which contains the relative pseudo-bulking analysis outputs, including a DEMs data.frame with a list of differential expressed m/z metabolites
#' @export
#'
#'
#' @examples
#' utils::str(formals(runDE))
runDE <- function(pooled_data, data, ident, output_dir, run_name, n, logFC_threshold, annotation.column, assay, return.individual = FALSE, verbose = TRUE){

  cellMetadata <- .cellMetadata(data)
  verbose_message(message_text = paste("Running limma DE Analysis for ", run_name, " -> with samples [", paste(unique(unlist(cellMetadata[[ident]])), collapse = ", "), "]"), verbose = verbose)

  annotation_result <- list()

  for (condition in unique(cellMetadata[[ident]])) {

    # Create groups
    groups <- SingleCellExperiment::colData(pooled_data)[[ident]]
    groups <- ifelse(groups == condition, "Comp_A", "Comp_B")

    # Extract continuous expression data (e.g., intensity matrix)
    expression_data <- SingleCellExperiment::counts(pooled_data)  # Or the assay holding your continuous data

    y <- edgeR::DGEList(SingleCellExperiment::counts(pooled_data), samples=SingleCellExperiment::colData(pooled_data)$orig.ident2, group = groups)

    y$samples$condition <- groups
    y$samples$ident <- as.character(SummarizedExperiment::colData(pooled_data)[[ident]])


    # Optional: If your data is raw intensities, log-transform it here (add small offset if needed)
    expression_data <- log2(expression_data + 1)

    #keep <- rowMeans(expression_data) > some_threshold  # Define threshold based on your data
    #expression_data <- expression_data[keep, ]

    # Create design matrix
    design <- model.matrix(~groups)
    design[, 2] <- 1 - design[, 2]  # To match your original contrast logic

    # Fit linear model
    fit <- limma::lmFit(expression_data, design)

    # Empirical Bayes moderation
    fit <- limma::eBayes(fit, robust = TRUE)

    # Use treat for log fold change threshold testing if desired
    res <- limma::treat(fit, lfc = log2(logFC_threshold), robust = TRUE)

    decisions <- limma::decideTests(res)
    all_decisions <- stats::setNames(
      as.integer(decisions[, ncol(decisions)]), rownames(decisions))

    res_table <- limma::topTreat(res, coef = ncol(fit$design), n = nrow(expression_data))

    res_table$regulate <- dplyr::recode(
      as.character(all_decisions[rownames(res_table)]),
      "0" = "Normal",
      "1" = "Up",
      "-1" = "Down"
    )

    # Order by p-value or FDR
    de_group_limma <- res_table[order(res_table$adj.P.Val), ]

    # Rename adj.P.Val to FDR
    colnames(de_group_limma) <- ifelse(colnames(de_group_limma) == "adj.P.Val", "FDR", colnames(de_group_limma))

    # Add gene/metabolite names
    de_group_limma$gene <- rownames(de_group_limma)

    # Add annotations if requested
    if (!is.null(annotation.column)) {
      annotation.data <- .featureMetadata(data, assay)
      if (!(annotation.column %in% colnames(annotation.data))) {
        stop("The annotation column does not exist in the assay feature metadata.")
      } else {
        annotation.data_subset <- annotation.data[rownames(de_group_limma), ]
        de_group_limma$annotations <- annotation.data_subset[[annotation.column]]
      }
    }

    # Write CSV output if directory specified
    if (!is.null(output_dir)) {
      utils::write.csv(de_group_limma, file.path(output_dir, paste0(condition, "_", run_name, ".csv")))
    }

    # Store results
    y$DEMs <- de_group_limma
    annotation_result[[condition]] <- y

    verbose_message(message_text = paste("Analysis complete for condition:", condition), verbose = verbose)

  }



  if (return.individual){
    annotation_result <- lapply(names(annotation_result), function(x){
      annotation_result[[x]]$DEMs$cluster <- x
      annotation_result[[x]]
    })

    return(annotation_result)
  } else {

    edger <- edgeR::DGEList(
      counts = annotation_result[[1]]$counts,
      samples = annotation_result[[1]]$samples
    )
    edger$samples$group <- edger$samples$ident
    edger$samples$condition <- NULL

    dems <- lapply(names(annotation_result), function(x){
      annotation_result[[x]]$DEMs$cluster <- x
      rownames(annotation_result[[x]]$DEMs) <- NULL
      annotation_result[[x]]$DEMs
    })

    combined_dems <- do.call(rbind, dems)
    rownames(combined_dems) <- 1:length(combined_dems$cluster)

    edger$DEMs <- combined_dems
    return(edger)
  }

}


#' Finds differentially expressed m/z values/metabolites between all comparison groups.
#'
#' @param data A Bioconductor experiment containing mz values for differential expression analysis.
#' @param ident A character string defining the metadata column or groups to compare mz values between.
#' @param n An integer that defines the number of pseudo-replicates (pools) per sample (default = 3).
#' @param logFC_threshold A numeric value indicating the logFC threshold to use for defining significant genes (default = 1.2).
#' @param DE_output_dir A character string defining the directory path for all output files to be stored. This path must a new directory. Else, set to NULL as default.
#' @param run_name A character string defining the title of this DE analysis that will be used when saving DEMs to .csv file (default = 'findAllDEMs').
#' @param annotation.column Character string defining the column where annotation information is stored in the assay metadata. This requires annotateSM() to be run where the default column to store annotations is "all_IsomerNames" (default = "None").
#' @param assay A character string defining the assay where the mz count data and annotations are stored (default = "Spatial").
#' @param slot Character string defining the assay storage slot to pull the relative mz intensity values from. Note: EdgeR requires raw counts, all values must be positive (default = "counts").
#' @param return.individual Boolean value defining whether to return a list of individual edgeR objects for each designated ident. If FALSE, one merged edgeR object will be returned (default = FALSE).
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#' @param seed Numeric value used to set the seed for reproducible randomisation (default = 1234).
#'
#' @returns Returns an list() contains the EdgeR DE results. Pseudo-bulk counts are stored in $counts and DEMs are in $DEMs.
#' @export
#'
#' @examples
#' utils::str(formals(findAllDEMs))
findAllDEMs <- function(data, ident, n = 3, logFC_threshold = 1.2, DE_output_dir = NULL, run_name = "findAllDEMs", annotation.column = NULL, assay = "Spatial", slot = "counts", return.individual = FALSE, verbose = TRUE, seed = 1234){

  if (!(is.null(DE_output_dir))){
    if (dir.exists(DE_output_dir)){
      warning("Please supply a directory path that doesn't already exist")
      stop("dir.exists(DE_output_dir) = TRUE")
    } else{
      dir.create(DE_output_dir)
    }
  }

  if (!is.factor(.cellMetadata(data)[[ident]])){
    stop("ident provided is not a factor! Please convert the ident column using `factor()` ...")
  }

  #Step 1: Run Pooling to split each unique ident into 'n' number of pseudo-replicate pools
  pooled_data <- runPooling(data,ident, n = n, assay = assay, slot = slot, verbose = verbose, seed = seed)

  #Step 2: Run EdgeR to calculate differentially expressed m/z peaks
  DEM_results <- runDE(pooled_data, data, ident = ident, output_dir = DE_output_dir, run_name = run_name, n=n, logFC_threshold=logFC_threshold, annotation.column = annotation.column, assay = assay, verbose = verbose, return.individual = return.individual)

  # Returns an EDGEr object which contains the pseudo-bulk counts in $counts and DEMs in $DEMs
  return(DEM_results)

}




#' Heatmap of Differentially Expressed Metabolites
#'
#' Generates a heatmap of DEMs generated from edgeR analysis run using `findAllDEMs()`.
#' This function uses `pheatmap` to plot data.
#'
#' @param edgeR_output A list containing outputs from edgeR analysis (from findAllDEMs()). This includes pseudo-bulked counts and DEMs.
#' @param n A numeric integer that defines the number of UP and DOWN regulated peaks to plot (default = 25).
#' @param only.pos Boolean indicating if only positive markers should be returned (default = FALSE).
#' @param FDR.threshold Numeric value that defines the FDR threshold to use for defining most significant results (default = 0.05).
#' @param logfc.threshold Numeric value that defines the logFC threshold to use for filtering significant results (default = 0.5).
#' @param order.by Character string defining which parameter to order markers by, options are either 'FDR' or 'logFC' (default = "FDR").
#' @param scale A character string indicating if the values should be centered and scaled in either the row direction or the column direction, or none. Corresponding values are "row", "column" and "none"
#' @param color A vector of colors used in heatmap (default = grDevices::colorRampPalette(c("navy", "white", "red"))(50)).
#' @param cluster_cols Boolean value determining if columns should be clustered or hclust object (default = FALSE).
#' @param cluster_rows Boolean value determining if rows should be clustered or hclust object (default = TRUE).
#' @param fontsize_row A numeric value defining the fontsize of rownames (default = 15).
#' @param fontsize_col A numeric value defining the fontsize of colnames (default = 15).
#' @param cutree_cols A numeric value defining the number of clusters the columns are divided into, based on the hierarchical clustering(using cutree), if cols are not clustered, the argument is ignored (default = 9).
#' @param silent Boolean value indicating if the plot should not be draw (default = TRUE).
#' @param plot_annotations_column Character string indicating the column name that contains the metabolite annotations to plot. Annotations = TRUE must be used in findAllDEMs() for edgeR output to include annotations. If plot_annotations_column = NULL, m/z vaues will be plotted (default = NULL).
#' @param save_to_path Character string defining the full filepath and name of the plot to be saved as.
#' @param plot.save.width Integer value representing the width of the saved pdf plot (default = 20).
#' @param plot.save.height Integer value representing the height of the saved pdf plot (default = 20).
#' @param nlabels.to.show Numeric value defining the number of annotations to show per m/z (default = NULL).
#' @param annotation_colors List for specifying annotation_row and annotation_col track colors manually. Check pheatmap R-Package documentation for details. If set to 'NA', default coloring will be used (default = NA).
#'
#' @returns A heatmap plot of significantly differentially expressed metabolites defined in the edgeR ouput object.
#' @export
#'
#' @import dplyr
#'
#' @examples
#' utils::str(formals(demsHeatmap))
#'
#' # demsHeatmap(DEMs)
demsHeatmap <- function(edgeR_output,
                         n = 5,
                         only.pos = FALSE,
                         FDR.threshold = 0.05,
                         logfc.threshold = 0.5,
                         order.by = "FDR",
                         scale ="row",
                         color = grDevices::colorRampPalette(c("navy", "white", "red"))(50),
                         cluster_cols = FALSE,
                         cluster_rows = TRUE,
                         fontsize_row = 15,
                         fontsize_col = 15,
                         cutree_cols = 9,
                         silent = TRUE,
                         plot_annotations_column = NULL,
                         save_to_path = NULL,
                         plot.save.width = 20,
                         plot.save.height = 20,
                         nlabels.to.show = NULL,
                         annotation_colors = NULL){


  degs <- edgeR_output$DEMs
  degs <- subset(degs, FDR < FDR.threshold)

  if (order.by == "FDR"){

    grouped_pos<- degs %>%
      group_by(cluster) %>%
      filter( logFC > logfc.threshold) %>%
      arrange(desc(regulate)) %>%
      slice_head(n = n)


    if (only.pos) {
      grouped_neg <- NULL

    } else {
      grouped_neg <- degs %>%
        group_by(cluster) %>%
        filter(logFC < - logfc.threshold) %>%
        arrange(regulate) %>%
        slice_head(n = n)
    }
    df <- do.call(rbind, list(grouped_pos,grouped_neg))
    df <- df[order(df$cluster, dplyr::desc(df$regulate)), ]

  } else {
    if ( order.by != "logFC"){
      warning("order.by has invalid argument. Must be either 'FDR' or 'logFC'. Heatmap defaulting to order by logFC")
    }

    grouped_pos<- degs %>%
      group_by(cluster) %>%
      filter(logFC > logfc.threshold) %>%
      arrange(-logFC) %>%
      slice_head(n = n)


    if (only.pos) {
      grouped_neg <- NULL
    } else {
      grouped_neg <- degs %>%
        group_by(cluster) %>%
        filter(logFC < - logfc.threshold) %>%
        arrange(logFC) %>%
        slice_head(n = n)
    }
    df <- do.call(rbind, list(grouped_pos,grouped_neg))
    df <- df[order(df$cluster, -df$logFC), ]
  }



  col_annot <- data.frame(sample = edgeR_output$samples$ident)
  row.names(col_annot) <- colnames(as.data.frame(edgeR::cpm(edgeR_output,log=TRUE)))

  if (!is.null(annotation_colors)){
    annotation_colors <- list(sample = unlist(annotation_colors))
  } else {
    annotation_colors <- NA
  }

  mtx <- as.matrix(as.data.frame(edgeR::cpm(edgeR_output,log=TRUE))[unique(df$gene),])
  if (!(is.null(plot_annotations_column))){
    if (is.null(edgeR_output$DEMs[[plot_annotations_column]])){
      warning("There are no annotations present in the edgeR_output object. Run 'annotateSM()' prior to 'findAllDEMs' and set annotations = TRUE .....\n Heatmap will plot default m/z values ... ")
    } else{
      if (!is.null(nlabels.to.show)){
        df[[plot_annotations_column]] <- labels_to_show(df[[plot_annotations_column]], n = nlabels.to.show)
      }
      rownames(mtx) <- unique(df[[plot_annotations_column]])
    }
  }

  p <- pheatmap::pheatmap(mtx,scale=scale,color=color,cluster_cols = cluster_cols, annotation_col=col_annot, cluster_rows = cluster_rows,
                          fontsize_row = fontsize_row, fontsize_col = fontsize_col, cutree_cols = cutree_cols, silent = silent, annotation_colors = annotation_colors)

   if (!(is.null(save_to_path))){
     savePheatmapAsPDF(pheatmap = p, filename = save_to_path, width = plot.save.width, height = plot.save.height)
   }

  return(p)
}


#' Saves a demsHeatmap as a PDF
#'
#' @param pheatmap A pheatmap plot object that is being saved.
#' @param filename Character string defining the full filepath and name of the plot to be saved as.
#' @param width Integer value representing the width of the saved pdf plot (default = 20).
#' @param height Integer value representing the height of the saved pdf plot (default = 20).
#'
#' @return The generated PDF path, invisibly.
#'
#' @export
#'
#' @examples
#' utils::str(formals(savePheatmapAsPDF))
#' # savePheatmapAsPDF(pheatmap, filename = "/Documents/plots/pheatmap1")
savePheatmapAsPDF <- function(pheatmap, filename, width=20, height=20){

  output_file <- paste0(filename,".pdf")
  pdf(output_file, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(pheatmap$gtable)
  dev.off()
  invisible(output_file)
}





########################################################################################################################################################################################################################
