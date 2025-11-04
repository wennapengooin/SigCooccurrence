#' Import Somatic Mutations from VCFs
#'
#' Reads somatic mutations from VCF files and prepares a
#' standardized GRangesList object for signature analysis.
#'
#' @details
#' This function serves as a standardized wrapper for reading somatic mutation
#' data from VCF files. It's designed to prepare the data for use in
#' subsequent filtering and analysis functions.
#'
#' It uses \code{\link[MutationalPatterns]{read_vcfs_as_granges}} and
#' accepts multiple VCF files, with each file treated as a separate sample.
#'
#' A reference genome package (e.g., "BSgenome.Hsapiens.UCSC.hg38")
#' must be installed for the specified \code{genome} string. This is
#' necessary to fetch the sequence context for mutations.
#'
#' @param vcf_files A character vector of file paths to VCF files
#'   (e.g., \code{c("sample1.vcf", "sample2.vcf")}).
#' @param genome A single string specifying the reference genome assembly
#'   (e.g., "hg38", "hg19").
#' @param sample_names (Optional) A character vector of sample names.
#'   If provided, it must be the same length as \code{vcf_files}.
#'   If \code{NULL} (default), sample names are derived from the VCF file names.
#'
#' @return A \code{GRangesList} object. Each element in the list is a
#'   \code{GRanges} object containing the mutations for a single sample.
#'
#' @importFrom MutationalPatterns read_vcfs_as_granges
#' @importFrom BSgenome getBSgenome
#' @importFrom tools file_ext file_path_sans_ext
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # --- Example with VCF files ---
#'
#'   # Get paths to example VCFs from the MutationalPatterns package
#'   vcf_files <- list.files(
#'     system.file("extdata", package = "MutationalPatterns"),
#'     pattern = "sample.vcf", full.names = TRUE
#'   )
#'
#'   # Define sample names
#'   sample_names <- c("sample_1", "sample_2", "sample_3")
#'
#'   # Load the mutations (assuming hg19 is installed)
#'   # BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
#'   muts_grl <- importMuts(vcf_files, genome = "hg19", sample_names = sample_names)
#'
#'   print(muts_grl)
#' }
importMuts <- function(vcf_files, genome, sample_names = NULL) {

  # --- 1. Input Validation ---
  if (!is.character(vcf_files) || length(vcf_files) == 0) {
    stop("`vcf_files` must be a character vector of one or more file paths.")
  }
  if (missing(genome) || !is.character(genome) || length(genome) != 1) {
    stop("`genome` must be a single string (e.g., 'hg38', 'hg19').")
  }

  # Check file extensions
  extensions <- unique(tools::file_ext(gsub("\\.gz$", "", vcf_files)))
  if (length(extensions) > 1 || !all(extensions == "vcf")) {
    stop("All files must be VCF files (with .vcf or .vcf.gz extension). ",
         "Found extensions: ", paste(extensions, collapse=", "))
  }

  # --- 2. Load Reference Genome ---
  # Dynamically build the BSgenome package name
  genome_pkg_name <- paste0("BSgenome.Hsapiens.UCSC.", genome)

  # Check if the genome package is installed
  if (!requireNamespace(genome_pkg_name, quietly = TRUE)) {
    stop(
      paste(
        "BSgenome package", genome_pkg_name, "is not installed.\n",
        "Please install it from Bioconductor using:\n",
        paste0("  BiocManager::install('", genome_pkg_name, "')")
      )
    )
  }

  # Load the genome object
  ref_genome <- BSgenome::getBSgenome(genome_pkg_name)

  # --- 3. Read VCF Data ---
  message("Detected VCF files. Reading mutations...")

  # If sample_names aren't provided, derive them from file names
  if (is.null(sample_names)) {
    message("`sample_names` not provided. Deriving from file names.")
    sample_names <- tools::file_path_sans_ext(basename(vcf_files))
    sample_names <- gsub("\\.vcf$", "", sample_names) # Handle .vcf.gz
  }

  # Final check on sample names and file paths
  if (length(sample_names) != length(vcf_files)) {
    stop("The number of `sample_names` must match the number of VCF files.")
  }

  # Use the MutationalPatterns function to read VCFs
  muts_grl <- MutationalPatterns::read_vcfs_as_granges(
    vcf_files = vcf_files,
    sample_names = sample_names,
    genome = ref_genome
  )

  message("Successfully imported mutations for ", length(muts_grl), " sample(s).")
  return(muts_grl)
}

