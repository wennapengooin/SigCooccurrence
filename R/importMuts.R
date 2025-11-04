#' Import Somatic Mutations from VCF Files
#'
#' Reads one or more VCF files into a \code{GRangesList} object, which is the
#' standard input for subsequent analysis functions.
#'
#' @details
#' This function serves as the primary data import step. It uses
#' \code{MutationalPatterns::read_vcfs_as_granges} to parse VCF files.
#'
#' It requires a reference genome (\code{BSgenome}) to fetch sequence context
#' and correctly classify mutations.
#'
#' The output \code{GRangesList} can be passed to \code{filterMuts} to
#' extract specific mutation types (SNVs, DBS, etc.).
#'
#' @param vcf_files A character vector of file paths to VCF files.
#' @param genome A string identifying the reference genome
#'   (e.g., "hg19", "hg38"). Must correspond to an installed
#'   \code{BSgenome} package (e.g., "BSgenome.Hsapiens.UCSC.hg19").
#' @param sample_names A character vector of sample names. If provided,
#'   its length must match the number of VCF files. If \code{NULL} (default),
#'   sample names are derived from the VCF file names.
#'
#' @return A \code{GRangesList} object where each element represents a sample.
#'
#' @importFrom MutationalPatterns read_vcfs_as_granges
#' @importFrom BSgenome getBSgenome
#' @importFrom tools file_path_sans_ext file_ext
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # --- Example with VCF files ---
#'
#'   # Get paths to example VCFs
#'   # (This requires the 'MutationalPatterns' package)
#'   vcf_files <- list.files(
#'     system.file("extdata", package = "MutationalPatterns"),
#'     pattern = "sample.vcf", full.names = TRUE
#'   )
#'
#'   # Define the reference genome
#'   # (Requires BiocManager::install("BSgenome.Hsapiens.UCSC.hg19"))
#'   genome <- "hg19"
#'
#'   # 1. Import mutations with derived sample names
#'   muts_grl <- importMuts(vcf_files, genome = genome)
#'
#'   # 2. Import mutations with custom sample names
#'   custom_names <- c("Sample_A", "Sample_B", "Sample_C")
#'   # (Assuming vcf_files has length 3)
#'   muts_grl_custom <- importMuts(vcf_files, genome, custom_names)
#' }
importMuts <- function(vcf_files, genome, sample_names = NULL) {

  # --- 1. Input Validation ---
  if (!is.character(vcf_files) || length(vcf_files) == 0) {
    stop("`vcf_files` must be a character vector of one or more file paths.")
  }
  if (any(!file.exists(vcf_files))) {
    stop("One or more VCF files do not exist.")
  }
  if (any(tools::file_ext(vcf_files) != "vcf")) {
    warning("Some input files do not have a .vcf extension.")
  }
  if (missing(genome) || !is.character(genome) || length(genome) != 1) {
    stop("`genome` must be a single string (e.g., 'hg19').")
  }
  if (!is.null(sample_names) && length(sample_names) != length(vcf_files)) {
    stop("The number of `sample_names` must match the number of VCF files.")
  }

  # --- 2. Resolve Sample Names ---
  if (is.null(sample_names)) {
    message("`sample_names` not provided. Deriving from file names.")
    sample_names <- tools::file_path_sans_ext(base::basename(vcf_files))
  }

  # --- 3. Load Genome and Read VCFs ---
  message("Detected VCF files. Reading mutations...")
  ref_genome_pkg <- paste0("BSgenome.Hsapiens.UCSC.", genome)

  if (!requireNamespace(ref_genome_pkg, quietly = TRUE)) {
    stop(
      "BSgenome package '", ref_genome_pkg, "' is not installed.\n",
      "Please install it using:\n",
      "  BiocManager::install('", ref_genome_pkg, "')"
    )
  }

  ref_genome <- BSgenome::getBSgenome(ref_genome_pkg)

  muts_grl <- MutationalPatterns::read_vcfs_as_granges(
    vcf_files = vcf_files,
    sample_names = sample_names,
    genome = ref_genome
  )

  # --- 4. Final Report ---
  message(
    "Successfully imported mutations for ",
    length(muts_grl), " sample(s)."
  )

  return(muts_grl)
}

