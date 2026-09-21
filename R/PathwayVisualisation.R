#' Visualise Significant Pathways
#'
#' Displays the pathway analysis results form running the 'fishersPathwayAnalysis()' function
#'
#' @param SpaMTP A SpatialExperiment used for fishersPathwayAnalysis.
#' @param pathway_df Dataframe containing the pathway enrichment results (output from SpaMTP::fishersPathwayAnalysis function).
#' @param assay Character string defining the SpaMTP assay that contains m/z values (default = "SPM").
#' @param slot Character string defining the assay slot contain the intensity values (default = "counts").
#' @param min_n Integer value specifying the minimum number of analytes required to be present in a pathway (default = 3).
#' @param p_val_threshold The p-val cutoff to keep the pathways generated from fisher exact test (default = "0.1").
#' @param method Character string defining the statistical method used to calculate hclust (default = "ward.D2").
#' @param verbose Boolean indicating whether to show informative messages. If FALSE these messages will be suppressed (default = TRUE).
#' @param database Optional named list of database resources, normally created
#'   by [loadSpaMTPDatabase()].
#' @param database_version SpaMTPdb/RaMP version used for pathway lookup.
#' @param database_source Database source; see [loadSpaMTPDatabase()].
#' @param database_local_dir Optional staged SpaMTPdb resource directory.
#' @param ... The arguments pass to stats::hclust
#'
#' @return A combined gg, ggplot object with pathway and dendrogram
#' @export
#'
#' @import grid
#'
#' @examples
#' utils::str(formals(visualisePathways))
#' #SpaMTP:::visualisePathways(SpaMTP =seurat,pathway_df = pathway_df,p_val_threshold = 0.1,assay = "Spatial",slot = "counts")
visualisePathways = function(SpaMTP,
                             pathway_df,
                             assay = "SPM",
                             slot = "counts",
                             min_n = 3,
                             p_val_threshold = 0.1,
                             method = "ward.D2",
                             verbose = TRUE,
                             database = NULL,
                             database_version = "latest",
                             database_source = c("auto", "spamtpdb"),
                             database_local_dir = NULL,
                             ...) {
  no_pathways_plot <- function(reason) {
    verbose_message(message_text = reason, verbose = verbose)
    ggplot2::ggplot() +
      ggplot2::annotate(
        "text",
        x = 0,
        y = 0,
        label = reason,
        size = 5
      ) +
      ggplot2::xlim(-1, 1) +
      ggplot2::ylim(-1, 1) +
      ggplot2::labs(title = "No pathways to visualise") +
      ggplot2::theme_void()
  }

  if (is.null(min_n)){
    stop("Incorrect minimum analyte number! `min_n` must be set to a value > 1. Please adjust this value accordingly ...")
  }
  pathway_df = pathway_df[which(pathway_df$analytes_in_pathways>=min_n),]
  if (nrow(pathway_df) == 0L) {
    return(no_pathways_plot(
      paste0("No pathways contain at least ", min_n, " detected analytes.")
    ))
  }

  pathway_df$duplicate_pathways <- rep(NA_character_, nrow(pathway_df))
  verbose_message(message_text = "Reducing synonymous pathways", verbose = verbose)
  index = seq_len(nrow(pathway_df))
  merged_pathways = data.frame()
  pb = txtProgressBar(
    min = 0,
    max = length(index),
    initial = 0,
    style = 3
  )
  while (length(index) != 0) {
    pattern = stringr::str_extract(pathway_df$pathway_name[index[1]], pattern = "[A-Z][A-Z]\\([a-z0-9]")
    if (length(pattern) != 0 & !is.na(pattern)) {
      pattern = sub("\\(", "\\\\(", pattern)
      name = strsplit(pathway_df$pathway_name[index[1]], pattern)[[1]][1]
      full_name = paste0(name, pattern)
      frst_ind = which(grepl(
        pattern = full_name,
        x = pathway_df$pathway_name
      ))
      all_pathways = pathway_df[frst_ind, ]
      second_ind = which(duplicated(all_pathways$p_val))
      index = index[-which(index %in% frst_ind)]
      if (length(second_ind)>0){
        name = stringr::str_trim(name, side = "right")

        unique_pathways <- all_pathways[-second_ind, ]
        for (i in seq_along(unique_pathways$p_val)) {
          # Get indices of duplicates for each unique p_val
          duplicates <- which(all_pathways$p_val == unique_pathways$p_val[i])
          duplicate_ids <- all_pathways$pathway_id[duplicates]
          unique_pathways$duplicate_pathways[i] <- paste(duplicate_ids, collapse = ", ")
        }

        unique_pathways$pathway_name <- name
        merged_pathways = rbind(merged_pathways, unique_pathways)
      } else{
        merged_pathways = rbind(merged_pathways, all_pathways)
      }

    } else{
      merged_pathways = rbind(merged_pathways, pathway_df[index[1], ])
      index = index[-1]
    }
    setTxtProgressBar(pb, nrow(pathway_df) - length(index))
  }
  close(pb)
  merged_pathways = merged_pathways %>% filter(p_val <= p_val_threshold) %>% mutate(signif_at_005level =   ifelse(p_val <= 0.05, "Significant", "Non-significant"))

  if (nrow(merged_pathways) == 0L) {
    return(no_pathways_plot(
      paste0("No pathways passed the p-value threshold (", p_val_threshold, ").")
    ))
  }

  retain_ind = seq_len(nrow(merged_pathways))

  gg_bar1 = with(
    merged_pathways,
    ggplot() +  geom_bar(
      aes(
        x = paste0(pathway_name, "(", pathway_id, ")"),
        y = total_in_pathways,
        color = "Total analytes in pathway"
      ),
      stat = "identity",
      fill = NA,
      position = "dodge"
    ) +
      geom_bar(
        aes(
          x = paste0(pathway_name, "(", pathway_id, ")"),
          y = analytes_in_pathways,
          colour = "Analytes detected in pathway"
        ),
        stat = "identity",
        fill = "blue",
        position = "dodge"
      ) +
      geom_point(
        aes(
          x = paste0(pathway_name, "(", pathway_id, ")"),
          y = 1.05 * max(total_in_pathways),
          size = p_val,
          fill = signif_at_005level
        ),
        color = "black",
        shape = 21,
        position = position_nudge(y = 0.5)
      ) +
      scale_size_area(max_size = 5, name = "p value") +
      scale_fill_manual(
        values = c(
          "Significant" = "green",
          "Non-significant" = "red"
        ),
        name = ""
      ) + scale_colour_manual(
        values = c(
          "Analytes detected in pathway" = "blue",
          "Total analytes in pathway" = "grey"
        ),
        name = ""
      ) +
      theme_minimal() + coord_flip() +
      theme(
        legend.position = "left",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )+labs(y= "Number of analytes in pathway", x = "Pathway names")
  )

  # Get the raster illustration for each pathway
  image_raster = list()
  verbose_message(message_text = "Generating images for visualising pathway enrichment across the sample ... ", verbose = verbose)
  pb = txtProgressBar(
    min = 0,
    max = nrow(merged_pathways),
    initial = 0,
    style = 3
  )
  mass_matrix = Matrix::t(.assayData(SpaMTP, assay, slot))

  for(z in seq_len(nrow(merged_pathways))){
    mzs = paste0("mz-",
                 stringr::str_extract_all(merged_pathways$adduct_info[z], "\\d+\\.\\d+")[[1]])
    if ("mz-" %in% mzs){
      image_raster[[z]] <- NULL
    } else{
      mat_ind = which(rownames(.assayData(SpaMTP, assay, slot)) %in% mzs)
      pca_result <- prcomp(mass_matrix[, mat_ind])
      nc = 3
      pca_re_df <- pca_result[["x"]][, 1:nc]
      pca_df_normalized <- as.data.frame(apply(
        pca_re_df,
        MARGIN = 2 ,
        FUN =  function(x)
          (x - min(x)) / (max(x) - min(x))
      ))
      coords <- as.data.frame(SpatialExperiment::spatialCoords(SpaMTP))[, c("x", "y")]
      # Convert the normalized UMAP result to an image matrix
      # 3 channels side by side

      rgb_m = array(dim = c(max(coords[, 1]), max(coords[, 2]), nc))
      for (j in 1:nc) {
        rgb_m[, , j] = matrix(pca_df_normalized[, j], nrow = max(coords[, 1]))
      }
      # Convert the image matrix to a raster object
      image_raster[[z]] <- as.raster(rgb_m)
      # # Plot the RGB image using ggplot2
      # ggplot() + annotation_custom(rasterGrob(image_raster, width = unit(1, "npc"), height = unit(1, "npc"))) +
      #   theme_void()
      setTxtProgressBar(pb, z)
    }
  }
  close(pb)
  if (length(image_raster)>0){
    for (k in 1:length(image_raster)) {
      if (!is.null(image_raster[[k]])){
        gg_bar1 =  gg_bar1 + ggplot2::annotation_custom(
          grid::rasterGrob(
            image_raster[[k]],
            width = unit(1, "npc"),
            height = unit(1, "npc")
          ),
          xmin = k - 0.5,
          xmax = k + 0.5,
          ymin = -3.5 ,
          ymax = -0.3
        )
      }
    }
  }
  gg_bar1 =  gg_bar1 + ylim(-2, max(merged_pathways$total_in_pathways) + 5)
  #gg_bar1

  if (nrow(merged_pathways) == 1L) {
    verbose_message(
      message_text = "\nOnly one pathway passed filtering; returning the pathway plot without a dendrogram.",
      verbose = verbose
    )
    return(gg_bar1)
  }

  database_resources <- .spamtp_db_bundle(
    c("analytehaspathway", "pathway"),
    database = database,
    version = database_version,
    source = match.arg(database_source),
    local_dir = database_local_dir
  )
  analytehaspathway <- database_resources$analytehaspathway
  pathway <- database_resources$pathway

  # Adding the dendrogram

  # data(pathway, package = "SpaMTP")
  # data(analytehaspathway, package = "SpaMTP")
  jaccard_matrix = matrix(nrow = nrow(merged_pathways),
                          ncol = nrow(merged_pathways))

  for (i in 1:(nrow(merged_pathways) - 1)) {
    pathway_id_i = merged_pathways$pathway_id[i]
    pathway_content_i = unique(analytehaspathway$rampId[which(analytehaspathway$pathwayRampId == pathway$pathwayRampId[which(pathway$sourceId == pathway_id_i)])])
    for (j in (i + 1):nrow(merged_pathways)) {
      pathway_id_j = merged_pathways$pathway_id[j]
      pathway_content_j = unique(analytehaspathway$rampId[which(analytehaspathway$pathwayRampId == pathway$pathwayRampId[which(pathway$sourceId == pathway_id_j)])])

      jc_simi = length(intersect(pathway_content_i, pathway_content_j)) / length(union(pathway_content_i, pathway_content_j))
      jaccard_matrix[i, j] = jaccard_matrix[j, i] = jc_simi
    }
  }
  diag(jaccard_matrix) = 1


  # Generate a dendrogram
  hc <- as.dendrogram(hclust(as.dist(jaccard_matrix), method = method, ...))
  segment_hc <- with(ggdendro::segment(ggdendro::dendro_data(hc)),
                     data.frame(
                       x = y,
                       y = x,
                       xend = yend,
                       yend = xend
                     ))
  pos_table <- with(ggdendro::dendro_data(hc)$labels,
                    data.frame(
                      y_center = x,
                      gene = as.character(label),
                      height = 1
                    ))
  axis_limits <- with(pos_table, c(min(y_center - 0.5 * height), max(y_center + 0.5 * height))) + 0.1 * c(-1, 1)

  plt_dendr <- ggplot(segment_hc) +
    geom_segment(aes(
      x = sqrt(x),
      y = y,
      xend = sqrt(xend),
      yend = yend
    )) +
    scale_x_continuous(expand = c(0, 0.5),
                       limits = c(0,max(segment_hc$xend))) +
    scale_y_continuous(
      expand = c(0, 0),
      breaks = pos_table$y_center,
      labels = pos_table$gene,
      limits = axis_limits,
      position = "right"
    ) +
    labs(
      x = "Jacard distance",
      y = "",
      colour = "",
      size = ""
    ) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())

  # # Combine the dendrogram and bar plot
  # combined_plot <- grid.arrange(gg_bar1,ggplotify::as.ggplot(dendro), widths = c(5, 1), nrow =1)

  combined_plot = cowplot::plot_grid(gg_bar1,
                                     plt_dendr,
                                     align = 'h',
                                     rel_widths = c(6, 1))
  verbose_message(message_text = "\nDone", verbose = verbose)
  # Show the plot
  return(combined_plot)
}



#' Plot significantly enriched pathways per region
#'
#' Visualisation of Set Enrichment Analysis Results from `SpaMTP::findRegionalPathways()`.
#'
#' @param regpathway A dataframe generated by the `SpaMTP::findRegionalPathways()` function, containing identified regional pathways.
#' @param ident.column A character string specifying the column name of the `regpathway` dataframe containing the idents or clusters to compare (default = "Cluster_id").
#' @param selected_pathways A character vector specifying the names or IDs of pathways to be included in the analysis (e.g., `c("Amino acid metabolism", "WP1902", "Aspartate and asparagine metabolism")`). This argument is not case-sensitive (default = NULL).
#' @param sig_cutoff Within-family BH-adjusted P-value cutoff for the point
#'   outline. NULL uses 0.05. Display selection does not recompute adjustment.
#' @param num_display An integer specifying the number of pathways to display in the plot. If set to null will plot the smaller of either 10 pathways or the number of unique pathways provided in 'regpathway' (default = NULL).
#' @param text_size A numeric value controlling the size of the text elements in the plot. If NULL the default text size is 12 (default = NULL).
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#' @param database Optional named list of database resources, normally created
#'   by [loadSpaMTPDatabase()].
#' @param database_version SpaMTPdb/RaMP version used for pathway lookup.
#' @param database_source Database source; see [loadSpaMTPDatabase()].
#' @param database_local_dir Optional staged SpaMTPdb resource directory.
#'
#' @return A `ggplot` object representing the set enrichment analysis results.
#' @export
#'
#' @examples
#' utils::str(formals(plotRegionalPathways))
#' # plotRegionalPathways(SpaMTP, ident = "clusters", regpathway = pathway_df)
plotRegionalPathways <- function(regpathway,
                                 ident.column = "Cluster_id",
                                 selected_pathways = NULL,
                                 sig_cutoff = NULL,
                                 num_display = NULL,
                                 text_size = NULL,
                                 verbose = TRUE,
                                 database = NULL,
                                 database_version = "latest",
                                 database_source = c("auto", "spamtpdb"),
                                 database_local_dir = NULL) {

  resources <- .spamtp_db_bundle(c("analytehaspathway", "pathway"),
    database = database, version = database_version, source = match.arg(database_source),
    local_dir = database_local_dir)
  required <- c("pathwayRampId", "pathwayName", ident.column, "NES", "size", "padj")
  if (!all(required %in% names(regpathway)))
    stop("Regional pathway plots require pathwayRampId, pathwayName, NES, size, padj and the region column.", call. = FALSE)
  if (!nrow(regpathway)) stop("No pathways to plot.", call. = FALSE)
  regpathway <- as.data.frame(regpathway)
  if (!is.null(selected_pathways)) regpathway <- regpathway[
    regpathway$pathwayRampId %in% selected_pathways | regpathway$pathwayName %in% selected_pathways |
      (!is.null(regpathway$sourceId) & regpathway$sourceId %in% selected_pathways), , drop = FALSE]
  importance <- tapply(abs(regpathway$NES), regpathway$pathwayRampId, function(z) max(z, na.rm = TRUE))
  ids <- head(names(sort(importance, decreasing = TRUE)), num_display %||% 10L)
  if (!length(ids)) stop("No selected pathways have estimable enrichment scores.", call. = FALSE)
  membership <- split(resources$analytehaspathway$rampId, resources$analytehaspathway$pathwayRampId)
  sets <- lapply(ids, function(id) unique(membership[[id]]))
  distance <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
  if (length(ids) > 1L) for (i in seq_len(length(ids) - 1L)) for (j in (i + 1L):length(ids)) {
    total <- length(union(sets[[i]], sets[[j]]))
    distance[i,j] <- distance[j,i] <- if (total) 1 - length(intersect(sets[[i]], sets[[j]])) / total else 1
  }
  ordering <- if (length(ids) > 1L) ids[stats::hclust(stats::as.dist(distance), method = "average")$order] else ids
  d <- regpathway[regpathway$pathwayRampId %in% ids, , drop = FALSE]
  labels <- stats::setNames(paste0(d$pathwayName, " [", d$pathwayRampId, "]"), d$pathwayRampId)
  d$.pathway <- factor(d$pathwayRampId, levels = rev(ordering))
  d$.region <- as.character(d[[ident.column]])
  d$.FDR <- factor(ifelse(!is.na(d$padj) & d$padj <= (sig_cutoff %||% .05), "Within-family BH FDR", "Above cutoff / unavailable"))
  limit <- max(abs(d$NES), na.rm = TRUE)
  p <- ggplot2::ggplot(d, ggplot2::aes(.data$.region, .data$.pathway)) +
    ggplot2::geom_point(ggplot2::aes(fill = .data$NES, size = .data$size, colour = .data$.FDR), shape = 21, stroke = 0.8) +
    ggplot2::scale_fill_gradient2(low = "#3b5b92", mid = "#f7f7f3", high = "#b73e45", limits = c(-limit, limit)) +
    ggplot2::scale_colour_manual(values = c("Above cutoff / unavailable" = "#d0d7dc", "Within-family BH FDR" = "#132f40")) +
    ggplot2::scale_size_area(max_size = 8) +
    ggplot2::scale_y_discrete(labels = labels[rev(ordering)]) +
    ggplot2::labs(x = ident.column, y = NULL, fill = "NES", size = "Tested members", colour = "Outline",
      subtitle = "Rows ordered by 1 - Jaccard member overlap; outlines use adjusted P values") +
    ggplot2::theme_minimal(base_size = text_size %||% 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "bottom")
  attr(p, "pathway_distance") <- distance
  attr(p, "selected_pathways") <- ids
  attr(p, "significance_column") <- "padj"
  p

}
