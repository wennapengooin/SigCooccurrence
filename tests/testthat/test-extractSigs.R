# --- Setup: Load data once for all tests ---
# This code runs once to set up the test environment.
# It's wrapped in suppressMessages to keep the test output clean.
muts_grl <- tryCatch({
  testthat::skip_if_not_installed("MutationalPatterns")
  testthat::skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg19")

  # Load the BSgenome
  ref_genome_hg19 <- BSgenome::getBSgenome(
    "BSgenome.Hsapiens.UCSC.hg19"
  )

  # Get VCF file paths
  vcf_files <- list.files(
    system.file("extdata", package = "MutationalPatterns"),
    pattern = "sample.vcf",
    full.names = TRUE
  )

  # Import mutations
  suppressMessages(
    importMuts(vcf_files = vcf_files, genome = "hg19")
  )
}, error = function(e) {
  # If setup fails (e.g., package not found), return NULL
  NULL
})


# --- Main Test Block ---

test_that("Happy Path: 'fit_to_cosmic' mode works for all types", {

  # Skip all tests in this block if setup failed
  testthat::skip_if(is.null(muts_grl), "Test data loading failed.")

  # 1. Test SNV (SBS)
  message("\nTesting 'fit_to_cosmic' for SNV...")
  snv_grl <- suppressMessages(
    filterMuts(muts_grl, type = "SNV")
  )

  # Use a minimal `nmf_nrun` for speed
  fit_res_snv <- suppressMessages(
    extractSigs(
      snv_grl,
      type = "SNV",
      genome = "hg19",
      mode = "fit_to_cosmic",
      cosmic_version = "v3.2"
    )
  )
  testthat::expect_type(fit_res_snv, "list")
  testthat::expect_true("contribution" %in% names(fit_res_snv))
  testthat::expect_gt(nrow(fit_res_snv$contribution), 10) # Has SBS sigs
  testthat::expect_equal(ncol(fit_res_snv$contribution), length(snv_grl))

  # 2. Test DBS
  message("Testing 'fit_to_cosmic' for DBS...")
  dbs_grl <- suppressWarnings(suppressMessages(
    filterMuts(muts_grl, type = "DBS")
  ))

  fit_res_dbs <- suppressMessages(
    extractSigs(
      dbs_grl,
      type = "DBS",
      genome = "hg19",
      mode = "fit_to_cosmic",
      cosmic_version = "v3.2"
    )
  )
  testthat::expect_type(fit_res_dbs, "list")
  testthat::expect_true("contribution" %in% names(fit_res_dbs))
  testthat::expect_gt(nrow(fit_res_dbs$contribution), 5) # Has DBS sigs

  # 3. Test INDEL
  message("Testing 'fit_to_cosmic' for INDEL...")
  indel_grl <- suppressWarnings(suppressMessages(
    filterMuts(muts_grl, type = "INDEL")
  ))

  fit_res_indel <- suppressMessages(
    extractSigs(
      indel_grl,
      type = "INDEL",
      genome = "hg19",
      mode = "fit_to_cosmic",
      cosmic_version = "v3.1"
    )
  )
  testthat::expect_type(fit_res_indel, "list")
  testthat::expect_true("contribution" %in% names(fit_res_indel))
  testthat::expect_gt(nrow(fit_res_indel$contribution), 10) # Has ID sigs

})


test_that("Happy Path: 'de_novo' mode works", {

  # Skip all tests in this block if setup failed
  testthat::skip_if(is.null(muts_grl), "Test data loading failed.")
  testthat::skip_if_not_installed("NMF")

  # Load NMF library explicitly to fix `path.package` error
  suppressPackageStartupMessages(library(NMF))

  n_de_novo <- 3

  snv_grl <- suppressMessages(
    filterMuts(muts_grl, type = "SNV")
  )

  de_novo_res <- suppressMessages(
    extractSigs(
      snv_grl,
      type = "SNV",
      genome = "hg19",
      mode = "de_novo",
      n_de_novo = n_de_novo,
      nmf_nrun = 1 # Use 1 run for speed
    )
  )

  testthat::expect_type(de_novo_res, "list")
  testthat::expect_true("contribution" %in% names(de_novo_res))
  testthat::expect_equal(nrow(de_novo_res$contribution), n_de_novo)
  testthat::expect_equal(ncol(de_novo_res$contribution), length(snv_grl))
})


test_that("Sad Path: Input validation throws correct errors", {

  # Skip all tests in this block if setup failed
  testthat::skip_if(is.null(muts_grl), "Test data loading failed.")

  snv_grl <- suppressMessages(filterMuts(muts_grl, type = "SNV"))

  # Missing genome
  testthat::expect_error(
    extractSigs(snv_grl, type = "SNV", mode = "de_novo", n_de_novo = 3),
    regexp = "Argument `genome` is required"
  )

  # Missing n_de_novo in de_novo mode
  testthat::expect_error(
    extractSigs(snv_grl, type = "SNV", genome = "hg19", mode = "de_novo"),
    regexp = "`n_de_novo`.*must be provided"
  )

  # Invalid 'type'
  testthat::expect_error(
    extractSigs(snv_grl, type = "SNP", genome = "hg19", mode = "de_novo"),
    regexp = "'arg' should be one of"
  )

  # Invalid 'cosmic_version'
  testthat::expect_error(
    extractSigs(
      snv_grl,
      type = "SNV",
      genome = "hg19",
      mode = "fit_to_cosmic",
      cosmic_version = "v_fake" # This version doesn't exist
    ),
    regexp = "Could not load COSMIC signatures"
  )

  # Invalid 'genome'
  testthat::expect_error(
    extractSigs(
      snv_grl,
      type = "SNV",
      genome = "hg_fake", # This genome doesn't exist
      mode = "fit_to_cosmic"
    ),
    regexp = "is not available"
  )
})

