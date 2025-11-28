#' Calculate Pairwise Signature Co-occurrence
#'
#' This function takes the signature exposure matrix (the `$contribution` matrix
#' from `extractSigs`) and calculates a pairwise correlation matrix between
#' all signatures.
#'
#' @details
#' The exposure matrix has signatures in rows and samples in columns. To calculate
#' the correlation between signatures *across samples*, this function
#' transposes the matrix and uses `stats::cor()`.
#'
#' The `method` argument allows you to choose the type of correlation:
#' \itemize{
#'   \item \strong{"spearman"} (default): Spearman's rank correlation. This is
#'     a non-parametric test that is robust to outliers and non-normal
#'     distributions, which is common in signature exposure data. It assesses
#'     monotonic relationships (i.e., as one signature's exposure increases,
#'     the other consistently increases or decreases, but not necessarily
#'     linearly).
#'   \item \strong{"pearson"}: Pearson's linear correlation. This assesses
#'     linear relationships.
#' }
#' The resulting square matrix is the core data for plotting heatmaps and
#' network graphs.
#'
#' @param exposure_matrix A numeric matrix of signature exposures
#'   (e.g., the `$contribution` output from `extractSigs`).
#'   Must have signatures in rows and samples in columns.
#' @param method A character string indicating which correlation coefficient
#'   is to be computed. One of "spearman" (default) or "pearson".
#' @param save_csv Logical (default: FALSE). If TRUE, the co-occurrence matrix
#'   will be saved as a CSV file.
#' @param csv_path Character string. The file path where the CSV should be
#'   saved (e.g., "results/cooccur.csv"). Required if `save_csv` is TRUE.
#'
#' @return A square numeric matrix where both rows and columns are
#'   signature names, and the values are the correlation coefficients
#'   (from -1 to 1).
#'
#' @importFrom stats cor var
#' @importFrom utils write.csv
#'
#' @export
#'
#' @examples
#' \dontrun{
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
#'
#'   # Get the exposure matrix
#'   exposures <- fit_res$contribution
#'
#'   # Calculate the co-occurrence matrix
#'   cooccur_matrix <- calcCooccur(exposures, method = "spearman")
#'
#'   # Calculate AND save to CSV
#'   cooccur_matrix <- calcCooccur(exposures,
#'                                 method = "spearman",
#'                                 save_csv = TRUE,
#'                                 csv_path = "cooccur_results.csv")
#'
#'   # View the top of the matrix
#'   print(round(cooccur_matrix[1:5, 1:5], 2))
#' }
calcCooccur <- function(exposure_matrix,
                        method = c("spearman", "pearson"),
                        save_csv = FALSE,
                        csv_path = NULL) {

  # --- Input Validation -----------------
  method <- match.arg(method)

  if (!is.matrix(exposure_matrix) || !is.numeric(exposure_matrix)) {
    stop("`exposure_matrix` must be a numeric matrix.")
  }

  if (is.null(rownames(exposure_matrix))) {
    stop("`exposure_matrix` must have rownames (signature names).")
  }

  if (nrow(exposure_matrix) < 2) {
    stop("`exposure_matrix` must have at least 2 signatures (rows) to correlate.")
  }

  if (ncol(exposure_matrix) < 2) {
    stop("`exposure_matrix` must have at least 2 samples (columns) to correlate.")
  }

  if (save_csv) {
    if (is.null(csv_path) || !is.character(csv_path)) {
      stop("`csv_path` must be provided as a string when `save_csv` is TRUE.")
    }
    # wheck if directory exists
    dir_name <- dirname(csv_path)
    if (dir_name != "." && !dir.exists(dir_name)) {
      stop("Directory '", dir_name, "' does not exist. Please create it first.")
    }
  }

  # warn if there are rows with no variance
  row_vars <- apply(exposure_matrix, 1, stats::var)
  if (any(row_vars == 0)) {
    zero_var_sigs <- rownames(exposure_matrix)[row_vars == 0]
    warning(
      "The following signatures have zero variance (e.g., all 0s) ",
      "across all samples and will result in 'NA' correlations:\n",
      paste(zero_var_sigs, collapse = ", ")
    )
  }

  # --- Calculate Correlation -----------------

  # transpose matrix to make signatures the columns.
  message(
    "Calculating pairwise '", method, "' correlation for ",
    nrow(exposure_matrix), " signatures across ",
    ncol(exposure_matrix), " samples..."
  )

  # suppress old warning since we have custom, more informative warning.
  cooccur_matrix <- suppressWarnings(
    stats::cor(
      t(exposure_matrix),
      method = method,
      use = "pairwise.complete.obs" # Handle potential NAs
    )
  )

  message("Co-occurrence matrix calculated.")

  # --- Save CSV -----------------
  if (save_csv) {
    message("Saving co-occurrence matrix to: ", csv_path)
    utils::write.csv(cooccur_matrix, file = csv_path, row.names = TRUE)
  }

  return(cooccur_matrix)
}


# [END]
