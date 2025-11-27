# --- Helper Data -----------------
testthat::skip_if_not_installed("MutationalPatterns")
testthat::skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg19")
testthat::skip_if_not_installed("GenomicRanges") # Added for mock object

# get paths to the example VCF files
vcf_files <- list.files(
  system.file("extdata", package = "MutationalPatterns"),
  pattern = "sample.vcf",
  full.names = TRUE
)

# load the genome
test_genome <- "hg19"

suppressMessages({
  muts_grl <- importMuts(vcf_files = vcf_files, genome = test_genome)
})

# --- Tests -----------------

test_that("Happy Path: filterMuts works for all types", {

  # test SNV
  suppressMessages({
    snv_grl <- filterMuts(muts_grl, type = "SNV")
  })
  testthat::expect_s4_class(snv_grl, "GRangesList")
  # check that the number of samples is the same
  testthat::expect_equal(length(snv_grl), length(muts_grl))

  # test DBS
  suppressMessages(suppressWarnings({
    dbs_grl <- filterMuts(muts_grl, type = "DBS")
  }))
  testthat::expect_s4_class(dbs_grl, "GRangesList")
  testthat::expect_equal(length(dbs_grl), length(muts_grl))

  # test INDEL
  suppressMessages(suppressWarnings({
    indel_grl <- filterMuts(muts_grl, type = "INDEL")
  }))
  testthat::expect_s4_class(indel_grl, "GRangesList")
  testthat::expect_equal(length(indel_grl), length(muts_grl))
})


test_that("Sad Path: filterMuts throws correct errors", {

  # wrong input class
  testthat::expect_error(
    filterMuts(list("a" = 1), type = "SNV"),
    regexp = "GRangesList object"
  )

  # invalid 'type' argument
  testthat::expect_error(
    filterMuts(muts_grl, type = "SNP"), # "SNP" is not a valid type
    regexp = "'arg' should be one of"
  )

  # case-insensitivity
  testthat::expect_no_error(
    suppressMessages(filterMuts(muts_grl, type = "snv"))
  )
})


test_that("Warning: filterMuts warns on 0-mutation samples", {

  # create a mock GRangesList containing only sample 1
  muts_grl_mock <- muts_grl[1]
  first_sample_name <- names(muts_grl_mock)[1]

  # check that filtering for INDELs throws the expected warning.
  testthat::expect_warning(
    suppressMessages(filterMuts(muts_grl_mock, type = "INDEL")),
    regexp = first_sample_name
  )
})

