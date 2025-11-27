testthat::skip_if_not_installed("MutationalPatterns")
testthat::skip_if_not_installed("S4Vectors")
testthat::skip_if_not_installed("GenomicRanges")

# --- Setup -----------------

vcf_files <- list.files(
  system.file("extdata", package = "MutationalPatterns"),
  pattern = "sample.vcf", full.names = TRUE
)

vcf_files_3 <- vcf_files[1:3]
test_genome <- "hg19"

# --- Tests -----------------

test_that("Happy Path: VCFs are read correctly with custom names", {

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

  # check class
  testthat::expect_s4_class(muts_grl, "GRangesList")

  # check names
  testthat::expect_equal(names(muts_grl), custom_names)

  # check for content
  testthat::expect_gt(sum(S4Vectors::elementNROWS(muts_grl)), 0)
})


test_that("Happy Path: VCFs are read correctly with derived names", {

  testthat::skip_if_not_installed(
    paste0("BSgenome.Hsapiens.UCSC.", test_genome)
  )

  # get expected names
  expected_derived_names <- tools::file_path_sans_ext(base::basename(vcf_files))

  muts_grl <- suppressMessages(
    importMuts(
      vcf_files = vcf_files,
      genome = test_genome
    )
  )

  # check class
  testthat::expect_s4_class(muts_grl, "GRangesList")

  # check names
  testthat::expect_equal(names(muts_grl), expected_derived_names)

  # check for content
  testthat::expect_gt(sum(S4Vectors::elementNROWS(muts_grl)), 0)
})


test_that("Sad Path: Function throws correct errors", {

  # error on non-existent files
  testthat::expect_error(
    importMuts(vcf_files = "non_existent_file.vcf", genome = test_genome),
    regexp = "do not exist"
  )

  # error on missing genome
  testthat::expect_error(
    importMuts(vcf_files = vcf_files_3),
    regexp = "genome"
  )

  # error on uninstalled genome
  testthat::expect_error(
    importMuts(vcf_files = vcf_files_3, genome = "hg_fake"),
    regexp = "is not installed"
  )

  # error on sample_names mismatch
  testthat::expect_error(
    importMuts(
      vcf_files = vcf_files_3,
      genome = test_genome,
      sample_names = c("one", "two")
    ),
    regexp = "must match the number"
  )
})

