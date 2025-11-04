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
#'
#' @return A square numeric matrix where both rows and columns are
#'   signature names, and the values are the correlation coefficients
#'   (from -1 to 1).
#'
#' @importFrom stats cor var
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # --- Full Workflow ---
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
#'   # View the top of the matrix
#'   # print(round(cooccur_matrix[1:5, 1:5], 2))
#' }
calcCooccur <- function(exposure_matrix, method = c("spearman", "pearson")) {

  # --- 1. Input Validation ---
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

  # Warn if there are rows with no variance (e.g., all zeros)
  row_vars <- apply(exposure_matrix, 1, stats::var)
  if (any(row_vars == 0)) {
    zero_var_sigs <- rownames(exposure_matrix)[row_vars == 0]
    warning(
      "The following signatures have zero variance (e.g., all 0s) ",
      "across all samples and will result in 'NA' correlations:\n",
      paste(zero_var_sigs, collapse = ", ")
    )
  }

  # --- 2. Calculate Correlation ---

  # stats::cor() calculates correlations between *columns*.
  # Our matrix has signatures in *rows*.
  # We must transpose (t()) the matrix to make signatures the columns.
  message(
    "Calculating pairwise '", method, "' correlation for ",
    nrow(exposure_matrix), " signatures across ",
    ncol(exposure_matrix), " samples..."
  )

  # Suppress the "standard deviation is zero" warning from stats::cor,
  # because we already have our own custom, more informative warning.
  cooccur_matrix <- suppressWarnings(
    stats::cor(
      t(exposure_matrix),
      method = method,
      use = "pairwise.complete.obs" # Handle potential NAs
    )
  )

  message("Co-occurrence matrix calculated.")
  return(cooccur_matrix)
}

