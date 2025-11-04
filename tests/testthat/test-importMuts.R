# This file relies on data from MutationalPatterns
testthat::skip_if_not_installed("MutationalPatterns")
testthat::skip_if_not_installed("S4Vectors")
testthat::skip_if_not_installed("GenomicRanges")

# --- Setup: Load shared data ---

# Get VCF file paths from MutationalPatterns
vcf_files <- list.files(
  system.file("extdata", package = "MutationalPatterns"),
  pattern = "sample.vcf", full.names = TRUE
)

# Use only 3 for most tests to be faster
vcf_files_3 <- vcf_files[1:3]
test_genome <- "hg19"

# --- Tests ---

test_that("Happy Path: VCFs are read correctly with custom names", {

  # Skip if the required BSgenome is not installed
  testthat::skip_if_not_installed(
    paste0("BSgenome.Hsapiens.UCSC.", test_genome)
  )

  custom_names <- c("Sample_A", "Sample_B", "Sample_C")

  muts_grl <- suppressMessages(
    importMuts(
      vcf_files = vcf_files_3,
      genome = test_genome,
      sample_names = custom_names
    )
  )

  # 1. Check class
  testthat::expect_s4_class(muts_grl, "GRangesList")

  # 2. Check names
  testthat::expect_equal(names(muts_grl), custom_names)

  # 3. Check for content
  # Use S4Vectors::elementNROWS instead of BiocGenerics::lengths
  testthat::expect_gt(sum(S4Vectors::elementNROWS(muts_grl)), 0)
})


test_that("Happy Path: VCFs are read correctly with derived names", {

  # Skip if the required BSgenome is not installed
  testthat::skip_if_not_installed(
    paste0("BSgenome.Hsapiens.UCSC.", test_genome)
  )

  # Get expected names (e.g., "chr22_sample1", "chr22_sample2", ...)
  expected_derived_names <- tools::file_path_sans_ext(base::basename(vcf_files))

  muts_grl <- suppressMessages(
    importMuts(
      vcf_files = vcf_files,
      genome = test_genome
      # No sample_names provided
    )
  )

  # 1. Check class
  testthat::expect_s4_class(muts_grl, "GRangesList")

  # 2. Check names
  testthat::expect_equal(names(muts_grl), expected_derived_names)

  # 3. Check for content
  # Use S4Vectors::elementNROWS instead of BiocGenerics::lengths
  testthat::expect_gt(sum(S4Vectors::elementNROWS(muts_grl)), 0)
})


test_that("Sad Path: Function throws correct errors", {

  # Error on non-existent files
  testthat::expect_error(
    importMuts(vcf_files = "non_existent_file.vcf", genome = test_genome),
    regexp = "do not exist"
  )

  # Error on missing genome
  testthat::expect_error(
    importMuts(vcf_files = vcf_files_3), # Missing 'genome'
    regexp = "genome"
  )

  # Error on uninstalled genome
  testthat::expect_error(
    importMuts(vcf_files = vcf_files_3, genome = "hg_fake"),
    regexp = "is not installed"
  )

  # Error on sample_names mismatch
  testthat::expect_error(
    importMuts(
      vcf_files = vcf_files_3,
      genome = test_genome,
      sample_names = c("one", "two") # 2 names for 3 files
    ),
    regexp = "must match the number"
  )
})

