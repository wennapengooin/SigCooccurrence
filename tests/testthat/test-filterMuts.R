# This file tests the "filterMuts" function

# --- Helper Data ---
# We need to run importMuts to get real data for filtering
# Skip all tests if the required packages are not installed.
testthat::skip_if_not_installed("MutationalPatterns")
testthat::skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg19")
testthat::skip_if_not_installed("GenomicRanges") # Added for mock object

# Get paths to the example VCF files
vcf_files <- list.files(
  system.file("extdata", package = "MutationalPatterns"),
  pattern = "sample.vcf",
  full.names = TRUE
)

# Load the genome
test_genome <- "hg19"

# Run the import function to get the test object
# We silence messages here as they are not relevant for the test
suppressMessages({
  muts_grl <- importMuts(vcf_files = vcf_files, genome = test_genome)
})

# --- Tests ---

test_that("Happy Path: filterMuts works for all types", {

  # Test SNV
  suppressMessages({
    snv_grl <- filterMuts(muts_grl, type = "SNV")
  })
  testthat::expect_s4_class(snv_grl, "GRangesList")
  # Check that the number of samples is the same
  testthat::expect_equal(length(snv_grl), length(muts_grl))

  # Test DBS
  # Suppress the expected warnings, since we only care that it runs.
  suppressMessages(suppressWarnings({
    dbs_grl <- filterMuts(muts_grl, type = "DBS")
  }))
  testthat::expect_s4_class(dbs_grl, "GRangesList")
  testthat::expect_equal(length(dbs_grl), length(muts_grl))

  # Test INDEL
  # FIX: Suppress the expected warnings, since we only care that it runs.
  suppressMessages(suppressWarnings({
    indel_grl <- filterMuts(muts_grl, type = "INDEL")
  }))
  testthat::expect_s4_class(indel_grl, "GRangesList")
  testthat::expect_equal(length(indel_grl), length(muts_grl))
})


test_that("Sad Path: filterMuts throws correct errors", {

  # Error 1: Wrong input class
  testthat::expect_error(
    filterMuts(list("a" = 1), type = "SNV"),
    regexp = "GRangesList object"
  )

  # Error 2: Invalid 'type' argument
  # The error message from match.arg is different.
  testthat::expect_error(
    filterMuts(muts_grl, type = "SNP"), # "SNP" is not a valid type
    regexp = "'arg' should be one of"
  )

  # Error 3: Case-insensitivity (should NOT throw an error)
  testthat::expect_no_error(
    suppressMessages(filterMuts(muts_grl, type = "snv"))
  )
})


test_that("Warning: filterMuts warns on 0-mutation samples", {

  # Create a mock GRangesList containing only sample 1,
  # which is known to have SNVs (original > 0) but no INDELs (filtered = 0).
  muts_grl_mock <- muts_grl[1]
  first_sample_name <- names(muts_grl_mock)[1]

  # Check that filtering for INDELs (which this sample doesn't have)
  # throws the expected warning.
  testthat::expect_warning(
    suppressMessages(filterMuts(muts_grl_mock, type = "INDEL")),
    regexp = first_sample_name
  )
})

