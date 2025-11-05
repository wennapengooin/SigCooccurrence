#' Plot a Signature Co-occurrence Heatmap
#'
#' This function visualizes the co-occurrence matrix from `calcCooccur()` as
#' a clustered heatmap using the `pheatmap` package.
#'
#' @details
#' The function uses `pheatmap::pheatmap()` to generate the plot.
#' Hierarchical clustering is applied by default to both rows and columns
#' to group similar signatures together, and dendrograms are displayed.
#'
#' A diverging color scale (blue-white-red) is used to represent the
#' correlation values, with -1 (blue) for negative correlations,
#' 1 (red) for positive correlations, and 0 (white) for no correlation.
#'
#' @param cooccur_matrix A square numeric matrix from `calcCooccur()`, with
#'   signatures in both rows and columns.
#' @param cluster_rows A logical value (default: `TRUE`) indicating whether to
#'   cluster the rows.
#' @param cluster_cols A logical value (default: `TRUE`) indicating whether to
#'   cluster the columns.
#' @param title An optional title for the plot (default: "Signature
#'   Co-occurrence Heatmap").
#' @param ... Additional arguments passed to `pheatmap::pheatmap()`.
#'
#' @return A `pheatmap` object representing the heatmap.
#'
#' @references
#' Kolde R (2025). pheatmap: Pretty Heatmaps. \emph{R package version 1.0.13}.
#' \href{https://github.com/raivokolde/pheatmap}{Link}
#'
#' @importFrom stats hclust dist as.dendrogram
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # --- Full Workflow -----------------
#'   muts_grl <- importMuts(vcf_files, genome = "hg19")
#'   snv_grl <- filterMuts(muts_grl, type = "SNV")
#'
#'   fit_res <- extractSigs(
#'     snv_grl,
#'     type = "SNV",
#'     genome = "hg19",
#'     mode = "fit_to_cosmic",
#'     cosmic_version = "v3.2"
#'   )
#'   exposures <- fit_res$contribution
#'   cooccur_matrix <- calcCooccur(exposures)
#'
#'   # Plot the heatmap
#'   p <- plotCooccurHeatmap(cooccur_matrix)
#'   # To display the plot (it prints automatically in interactive sessions):
#'   # print(p)
#' }
plotCooccurHeatmap <- function(cooccur_matrix,
                               cluster_rows = TRUE,
                               cluster_cols = TRUE,
                               title = "Signature Co-occurrence Heatmap",
                               ...) {

  # --- Input Validation -----------------
  if (!is.matrix(cooccur_matrix) ||
      !is.numeric(cooccur_matrix) ||
      nrow(cooccur_matrix) != ncol(cooccur_matrix)) {
    stop("`cooccur_matrix` must be a square numeric matrix.")
  }
  if (is.null(rownames(cooccur_matrix)) || is.null(colnames(cooccur_matrix))) {
    stop("`cooccur_matrix` must have both rownames and colnames.")
  }

  # --- Define Color Scale -----------------
  # create a diverging color palette: Blue -> White -> Red
  color_palette <- grDevices::colorRampPalette(c("blue", "white", "red"))(100)

  # define the breaks for the color scale, centered at 0
  max_val <- max(abs(cooccur_matrix), na.rm = TRUE)
  breaks <- seq(-max_val, max_val, length.out = 101)

  # handle the case where all correlations are 0
  if (max_val == 0) {
    breaks <- seq(-1, 1, length.out = 101)
  }

  # --- Create pheatmap -----------------
  message("Generating heatmap plot...")

  # replace NAs with 0 for plotting, otherwise clustering can fail
  plot_matrix <- cooccur_matrix
  plot_matrix[is.na(plot_matrix)] <- 0

  gg_heatmap <- pheatmap::pheatmap(
    plot_matrix,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    color = color_palette,
    breaks = breaks,
    main = title,
    fontsize = 8,
    ...
  )

  return(gg_heatmap)
}


# [END]
