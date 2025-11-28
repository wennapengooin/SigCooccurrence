#' Summarize Signature Co-occurrence Patterns
#'
#' Takes the co-occurrence matrix and produces a summary table of
#' the strongest relationships. It identifies the top co-occurring and mutually
#' exclusive signature pairs.
#'
#' @details
#' The function flattens the square correlation matrix into a long-format
#' table of unique pairs. It then sorts these pairs by the strength of their
#' correlation.
#'
#' @param cooccur_matrix A square numeric matrix from `calcCooccur()`.
#' @param top_n An integer (default: 10). The number of top positive and
#'   top negative pairs to return.
#'
#' @return A data.frame with the following columns:
#'   \itemize{
#'     \item \code{Signature_A}: The name of the first signature.
#'     \item \code{Signature_B}: The name of the second signature.
#'     \item \code{Correlation}: The correlation coefficient.
#'     \item \code{Type}: A categorical label ("Co-occurrence" or
#'           "Mutual Exclusivity") based on the sign of the correlation.
#'   }
#'
#' @importFrom reshape2 melt
#' @importFrom dplyr filter arrange desc mutate select
#' @importFrom magrittr `%>%`
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # Get the top 5 co-occurring and top 5 mutually exclusive pairs
#'   summary_table <- summarizeCooccur(cooccur_matrix, top_n = 5)
#'
#'   print(summary_table)
#' }
summarizeCooccur <- function(cooccur_matrix, top_n = 10) {

  if (!is.matrix(cooccur_matrix) || !is.numeric(cooccur_matrix)) {
    stop("`cooccur_matrix` must be a numeric matrix.")
  }
  if (is.null(rownames(cooccur_matrix)) || is.null(colnames(cooccur_matrix))) {
    stop("`cooccur_matrix` must have rownames and colnames.")
  }

  long_df <- reshape2::melt(cooccur_matrix, na.rm = TRUE)
  colnames(long_df) <- c("Signature_A", "Signature_B", "Correlation")

  # ensure columns are character
  long_df$Signature_A <- as.character(long_df$Signature_A)
  long_df$Signature_B <- as.character(long_df$Signature_B)

  unique_pairs <- long_df %>%
    dplyr::filter(Signature_A < Signature_B)

  # get top positive correlations
  top_positive <- unique_pairs %>%
    dplyr::arrange(dplyr::desc(Correlation)) %>%
    dplyr::filter(Correlation > 0) %>%
    head(top_n) %>%
    dplyr::mutate(Type = "Co-occurrence")

  # get top negative correlations
  top_negative <- unique_pairs %>%
    dplyr::arrange(Correlation) %>%
    dplyr::filter(Correlation < 0) %>%
    head(top_n) %>%
    dplyr::mutate(Type = "Mutual Exclusivity")

  summary_table <- rbind(top_positive, top_negative)

  summary_table <- summary_table %>%
    dplyr::arrange(dplyr::desc(abs(Correlation)))

  return(summary_table)
}


# [END]
