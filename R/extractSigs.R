#' Extract or Fit Mutational Signatures
#'
#' This function serves two modes:
#' 1.  **"de_novo"**: Extracts a specified number of *de novo* mutational
#'     signatures from a mutation count matrix using NMF.
#' 2.  **"fit_to_cosmic"**: Fits known COSMIC mutational signatures to a
#'     mutation count matrix to find the optimal exposure (contribution)
#'     of each signature in each sample.
#'
#' @details
#' This function is a wrapper for core functions from the `MutationalPatterns`
#' package. It first builds the appropriate mutation count matrix
#' (e.g., 96-channel for SNV) from the provided `muts_grl`. It then
#' performs one of the two modes.
#'
#' The output is a list containing the results. For signature co-occurrence
#' analysis, you will be most interested in the `contribution` matrix
#' (exposures) from the result.
#'
#' @param muts_grl A \code{GRangesList} object, **pre-filtered** by
#'   \code{filterMuts} for a single mutation type (e.g., "SNV").
#' @param type A character string specifying the mutation type. Must match
#'   the filtering used to create \code{muts_grl}. One of "SNV", "DBS",
#'   or "INDEL". (Case-insensitive)
#' @param genome A character string specifying the reference genome
#'   (e.g., "hg19" or "hg38"). Must match the one used in \code{importMuts}.
#' @param mode A character string specifying the analysis mode. One of:
#'   \itemize{
#'     \item \code{"de_novo"}: Extract \code{n_de_novo} new signatures.
#'     \item \code{"fit_to_cosmic"}: Fit known COSMIC signatures.
#'   }
#' @param n_de_novo An integer. Required only if \code{mode = "de_novo"}.
#'   Specifies the. number of signatures to extract (the NMF rank).
#' @param nmf_nrun An integer. Optional, for \code{mode = "de_novo"} only.
#'   Specifies the number of NMF runs. Defaults to 100 for stability.
#'   (Note: Test scripts may use `nmf_nrun = 1` for speed).
#' @param cosmic_version A character string specifying the COSMIC signature
#'   version to use. Required only if \code{mode = "fit_to_cosmic"}.
#'   \itemize{
#'     \item For "SNV": "v3.1", "v3.2" (default: "v3.2")
#'     \item For "DBS": "v3.2" (default: "v3.2")
#'     \item For "INDEL": "v3.1" (default: "v3.1")
#'   }
#'
#' @return A list containing the results of the NMF or fitting,
#'   as returned by \code{MutationalPatterns}. This list includes
#'   \code{signatures} (the signature profiles) and \code{contribution}
#'   (the exposure matrix, with signatures in rows and samples in columns).
#'
#' @importFrom MutationalPatterns mut_matrix count_dbs_contexts
#'   count_indel_contexts get_known_signatures extract_signatures
#'   fit_to_signatures
#' @importFrom BSgenome getBSgenome
#' @importFrom magrittr `%>%`
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # --- Setup ---
#'   muts_grl <- importMuts(vcf_files, genome = "hg19")
#'   snv_grl <- filterMuts(muts_grl, type = "SNV")
#'
#'   # --- Example 1: Fit to COSMIC v3.2 SBS signatures ---
#'   fit_res <- extractSigs(
#'     snv_grl,
#'     type = "SNV",
#'     genome = "hg19",
#'     mode = "fit_to_cosmic",
#'     cosmic_version = "v3.2"
#'   )
#'   # Get the exposure matrix for the next step
#'   exposures <- fit_res$contribution
#'
#'   # --- Example 2: Extract 5 de novo SBS signatures ---
#'   # This requires more samples and mutations to be stable
#'   de_novo_res <- extractSigs(
#'     snv_grl,
#'     type = "SNV",
#'     genome = "hg19",
#'     mode = "de_novo",
#'     n_de_novo = 5,
#'     nmf_nrun = 100 # Use more runs for a real analysis
#'   )
#'   # Get the de novo exposure matrix
#'   de_novo_exposures <- de_novo_res$contribution
#' }
extractSigs <- function(muts_grl,
                        type = c("SNV", "DBS", "INDEL"),
                        genome,
                        mode = c("de_novo", "fit_to_cosmic"),
                        n_de_novo = NULL,
                        nmf_nrun = 100,
                        cosmic_version = NULL) {

  # --- 1. Input Validation ---
  if (!inherits(muts_grl, "GRangesList")) {
    stop("Input `muts_grl` must be a GRangesList object.")
  }
  if (missing(genome)) {
    stop("Argument `genome` is required (e.g., 'hg19').")
  }

  type <- match.arg(toupper(type), c("SNV", "DBS", "INDEL"))
  mode <- match.arg(mode)

  # Validate mode-specific arguments
  if (mode == "de_novo" && is.null(n_de_novo)) {
    stop("`n_de_novo` (number of signatures) must be provided for 'de_novo' mode.")
  }
  if (mode == "fit_to_cosmic" && is.null(cosmic_version)) {
    # Set default versions if not provided
    cosmic_version <- switch(type,
                             "SNV" = "v3.2",
                             "DBS" = "v3.2",
                             "INDEL" = "v3.1"
    )
    message(
      "No `cosmic_version` provided. Using default: '", cosmic_version, "'."
    )
  }

  # Load the reference genome
  message("Loading reference genome '", genome, "'...")
  ref_genome <- BSgenome::getBSgenome(
    paste0("BSgenome.Hsapiens.UCSC.", genome)
  )

  # --- 2. Build Mutation Count Matrix ---
  message(paste0("Building ", type, " mutation count matrix..."))

  # Use the correct MutationalPatterns function based on type
  mut_mat <- switch(type,
                    "SNV" = MutationalPatterns::mut_matrix(
                      vcf_list = muts_grl,
                      ref_genome = ref_genome,
                      extension = 1 # For 96-channel context
                    ),
                    "DBS" = MutationalPatterns::count_dbs_contexts(
                      vcf_list = muts_grl
                    ),
                    "INDEL" = MutationalPatterns::count_indel_contexts(
                      vcf_list = muts_grl
                    )
  )

  # Add a small epsilon to avoid NMF errors with 0-count samples
  mut_mat <- mut_mat + 0.0001

  # --- 3. Execute Selected Mode ---

  if (mode == "de_novo") {

    # --- De Novo Extraction ---
    message(paste0("Extracting ", n_de_novo, " de novo signatures..."))

    # FIX: Explicitly load NMF library to fix `path.package` error in tests
    if (!requireNamespace("NMF", quietly = TRUE)) {
      stop("Package 'NMF' is required for 'de_novo' mode. Please install it.")
    }
    suppressPackageStartupMessages(library(NMF))

    # Run NMF
    signatures_res <- MutationalPatterns::extract_signatures(
      mut_mat,
      rank = n_de_novo,
      nrun = nmf_nrun
    )

    message("De novo extraction complete.")
    return(signatures_res)

  } else {

    # --- Fit to COSMIC ---

    # Map our 'type' to MutationalPatterns 'muttype'
    muttype_arg <- switch(type,
                          "SNV" = "snv",
                          "DBS" = "dbs",
                          "INDEL" = "indel"
    )

    # Map our 'genome' to MutationalPatterns 'genome'
    genome_arg <- switch(genome,
                         "hg19" = "GRCh37",
                         "hg38" = "GRCh38",
                         stop("Invalid `genome`. Must be 'hg19' or 'hg38'.")
    )

    # Map our 'cosmic_version' to MutationalPatterns 'source'
    source_arg <- paste0("COSMIC_", cosmic_version)

    message(paste0(
      "Loading COSMIC ",
      type,
      " signatures (",
      source_arg,
      ", ",
      genome_arg,
      ")..."
    ))

    # Get the signatures
    known_sigs <- tryCatch({
      MutationalPatterns::get_known_signatures(
        muttype = muttype_arg,
        source = source_arg,
        genome = genome_arg
      )
    }, error = function(e) {
      stop(
        "Could not load COSMIC signatures.\n",
        "Check `type` ('", type, "'), `cosmic_version` ('", cosmic_version,
        "'), and `genome` ('", genome, "').\n",
        "Original error: ", e$message
      )
    })

    message(paste0("Fitting to ", ncol(known_sigs), " known signatures..."))

    # Fit signatures
    fit_res <- MutationalPatterns::fit_to_signatures(
      mut_mat,
      signatures = known_sigs
    )

    message("Fitting complete.")
    return(fit_res)
  }
}


# [END]
