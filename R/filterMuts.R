#' Filter a GRangesList by Mutation Type
#'
#' Subsets a \code{GRangesList} object to include only mutations of a
#' specified type (SNV, DBS, MBS, or INDEL).
#'
#' @details
#' This function is a wrapper for \code{MutationalPatterns::get_mut_type}.
#' It takes the full \code{GRangesList} from \code{importMuts} and
#' returns a new \code{GRangesList} containing only the desired
#' mutation type.
#'
#' As per \code{MutationalPatterns}, DBS/MBSs are assumed to be
#' called as separate SNVs and will be merged into single variants.
#'
#' @param muts_grl A \code{GRangesList} object, typically from
#'   \code{importMuts}.
#' @param type A character string specifying the mutation type to
#'   keep. One of "SNV", "DBS", "MBS", or "INDEL". (Case-insensitive)
#'
#' @return A \code{GRangesList} containing only mutations of the
#'   specified type.
#'
#' @references
#' Manders, F., Brandsma, A.M., de Kanter, J. et al. (2022).
#' MutationalPatterns: the one stop shop for the analysis of mutational
#' processes. \emph{BMC Genomics}, 23, 134.
#' \href{https://doi.org/10.1186/s12864-022-08357-3}{Link}
#'
#' @importFrom MutationalPatterns get_mut_type
#' @importFrom GenomicRanges GRangesList
#' @importFrom S4Vectors elementNROWS
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # Import all mutations
#'   muts_grl <- importMuts(vcf_files, genome = "hg19")
#'
#'   # Filter for SNVs (for SBS analysis)
#'   snv_grl <- filterMuts(muts_grl, type = "SNV")
#'
#'   # Filter for DBSs (for DBS analysis)
#'   dbs_grl <- filterMuts(muts_grl, type = "DBS")
#'
#'   # Filter for INDELs (for ID analysis)
#'   indel_grl <- filterMuts(muts_grl, type = "INDEL")
#' }
filterMuts <- function(muts_grl, type = c("SNV", "DBS", "MBS", "INDEL")) {

  # --- 1. Input Validation ---
  if (!inherits(muts_grl, "GRangesList")) {
    stop("Input `muts_grl` must be a GRangesList object.")
  }

  # Match the user's input type
  type <- match.arg(toupper(type), c("SNV", "DBS", "MBS", "INDEL"))
  type_lower <- tolower(type)

  # --- 2. Filter using MutationalPatterns ---
  message(paste0("Filtering for '", type, "' mutations..."))

  # Suppress messages from the underlying function (e.g., about merging)
  filtered_grl <- suppressMessages({
    MutationalPatterns::get_mut_type(muts_grl, type = type_lower)
  })

  # --- 3. Check for empty samples and warn ---
  # Use S4Vectors::elementNROWS instead of BiocGenerics::lengths
  original_lengths <- S4Vectors::elementNROWS(muts_grl)
  filtered_lengths <- S4Vectors::elementNROWS(filtered_grl)

  empty_samples <- names(filtered_grl[filtered_lengths == 0 & original_lengths > 0])

  if (length(empty_samples) > 0) {
    warning(
      "The following samples had 0 mutations of type '", type, "':\n",
      paste(empty_samples, collapse = ", ")
    )
  }

  message(
    "Filtering complete. Output contains ",
    sum(filtered_lengths), " '", type, "' mutations."
  )

  return(filtered_grl)
}


# [END]
