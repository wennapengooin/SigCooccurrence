# This file tests the "importMuts" function

# --- Helper variables ---

# CRITICAL: Skip this entire test file if MutationalPatterns (which provides
# the test data) is not installed.
testthat::skip_if_not_installed("MutationalPatterns")

# Get the paths to the example VCF files from MutationalPatterns
# We'll use these as our test data.
# The pattern "sample.vcf" finds the 9 "tissue-sample.vcf" files.
vcf_files <- list.files(
  system.file("extdata", package = "MutationalPatterns"),
  pattern = "sample.vcf", # Use the pattern from the function's example
  full.names = TRUE
)

# Define a set of correct sample names
sample_names_custom <- c("colon1", "colon2", "colon3",
                         "intestine1", "intestine2", "intestine3",
                         "liver1", "liver2", "liver3")

# Define the genome
test_genome <- "hg19"

# --- Tests ---

test_that("Happy Path: VCFs are read correctly with custom names", {

  # Skip this test if the hg19 BSgenome is not installed.
  testthat::skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg19")

  # Run the function
  muts_grl <- importMuts(
    vcf_files = vcf_files,
    genome = test_genome,
    sample_names = sample_names_custom
  )

  # Check 1: Is the output the correct class?
  # ***** THIS IS THE FIX *****
  # Use expect_s4_class for Bioconductor S4 objects
  testthat::expect_s4_class(muts_grl, "GRangesList")

  # Check 2: Is the length correct? (9 VCFs in = 9 samples out)
  testthat::expect_equal(length(muts_grl), 9)

  # Check 3: Are the names correct?
  testthat::expect_equal(names(muts_grl), sample_names_custom)
})


test_that("Happy Path: Sample names are derived correctly from filenames", {

  testthat::skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg19")

  # Run the function with sample_names = NULL
  muts_grl_derived <- importMuts(
    vcf_files = vcf_files,
    genome = test_genome,
    sample_names = NULL
  )

  # Check 1: Is the output the correct class? (Added for robustness)
  testthat::expect_s4_class(muts_grl_derived, "GRangesList")

  # The function's logic should derive these names from the filenames
  expected_derived_names <- c("colon1-sample", "colon2-sample", "colon3-sample",
                              "intestine1-sample", "intestine2-sample", "intestine3-sample",
                              "liver1-sample", "liver2-sample", "liver3-sample")

  # Check: Are the derived names correct?
  testthat::expect_equal(names(muts_grl_derived), expected_derived_names)
})


test_that("Sad Path: Function throws correct errors", {

  # Error 1: Missing 'genome' argument
  testthat::expect_error(
    importMuts(vcf_files = vcf_files),
    regexp = "genome" # Check that the error message mentions "genome"
  )

  # Error 2: Wrong file type (e.g., a text file)
  not_a_vcf <- tempfile(fileext = ".txt")
  writeLines("test", not_a_vcf)
  testthat::expect_error(
    importMuts(vcf_files = not_a_vcf, genome = test_genome),
    regexp = "All files must be VCF"
  )
  on.exit(unlink(not_a_vcf)) # Clean up the temp file

  # Error 3: Mismatch in length of sample_names and vcf_files
  testthat::expect_error(
    importMuts(
      vcf_files = vcf_files,
      genome = test_genome,
      sample_names = c("one", "two") # Only 2 names for 9 files
    ),
    regexp = "must match"
  )

  # Error 4: Missing BSgenome package
  testthat::expect_error(
    importMuts(vcf_files = vcf_files, genome = "fake_genome_123"),
    regexp = "is not installed"
  )
})

