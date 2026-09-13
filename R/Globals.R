# Namespace declarations for base graphics/statistics helpers and non-standard
# evaluation columns used by dplyr, data.table, and ggplot2 expressions.

#' @importFrom data.table :=
#' @importFrom grDevices as.raster chull dev.off pdf
#' @importFrom graphics abline arrows par points text
#' @importFrom methods new
#' @importFrom purrr map_dfr
#' @importFrom stats as.dendrogram as.dist cor density dnorm hclust median model.matrix na.omit p.adjust prcomp quantile setNames
#' @rawNamespace import(ggplot2, except = last_plot)
#' @rawNamespace import(tidyr, except = c(expand, pack, unpack))
#' @noRd
NULL

utils::globalVariables(c(
  ".SD", "Adduct", "Error", "FDR", "Features", "Formula", "IsomerNames", "Isomers",
  "Isomers_IDs", "NES", "Ramp_IDs", "Score", "Significance", "abs_cor",
  "adduct_info", "all_IsomerNames", "analytes_in_pathways", "annotation",
  "bottom_cutoff", "cluster", "color", "commonName", "cor", "exactmass",
  "fdr", "formula", "gene", "gene_id_list", "gene_name_list", "grid_col",
  "grid_row", "group_importance", "isomers", "isomers_inchikey",
  "isomers_names", "logFC", "mass", "max_cor", "max_value",
  "metabolite_id_list", "metabolite_name_list",
  "moransi.spatially.variable.rank", "mz", "mz_name", "mz_names",
  "n_sig_path", "name", "observed_mz", "p_val", "p_val_adj", "pathnameid",
  "pathwayCategory", "pathwayName", "pathwayRampId", "pathwaySource.y",
  "pathway_id", "pathway_name", "pixel", "ppm_diff", "present", "pval",
  "pval_adj", "rampId", "rampId.y", "ramp_id", "ratio", "reaction_type",
  "regulate", "size", "sourceId", "top_cutoff", "total_in_pathways", "type",
  "intensity", "value", "var", "variable", "x", "xend", "y", "yend", "z_cor", "z_path",
  "z_score", "sample_id", "score", "dataset", "opacity"
))
