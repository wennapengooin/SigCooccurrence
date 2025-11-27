# Test script for summarizeCooccur()

# --- Setup -----------------
# create mock correlation matrix
# SigA vs SigB:  0.9 (strong co-occurrence)
# SigA vs SigC: -0.8 (strong mutual exclusivity)
# SigA vs SigD:  0.1 (weak co-occurrence)
# SigB vs SigC: -0.5 (moderate mutual exclusivity)
# SigD vs SigE:  0.0 (no correlation)

mock_cor_matrix <- matrix(
  c( 1.0,  0.9, -0.8,  0.1,  0.0,
     0.9,  1.0, -0.5,  0.2,  0.0,
     -0.8, -0.5,  1.0, -0.1,  0.0,
     0.1,  0.2, -0.1,  1.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  1.0),
  nrow = 5,
  dimnames = list(
    c("SigA", "SigB", "SigC", "SigD", "SigE"),
    c("SigA", "SigB", "SigC", "SigD", "SigE")
  )
)

# --- Tests -----------------

test_that("Happy Path: Function returns correct structure", {

  summary_df <- summarizeCooccur(mock_cor_matrix)

  # check class
  testthat::expect_s3_class(summary_df, "data.frame")

  # Check columns
  expected_cols <- c("Signature_A", "Signature_B", "Correlation", "Type")
  testthat::expect_true(all(expected_cols %in% colnames(summary_df)))

  # Check that self-correlations (correlation = 1) are removed
  testthat::expect_equal(max(summary_df$Correlation), 0.9)
})

test_that("Happy Path: Identifies Co-occurrence and Mutual Exclusivity", {

  summary_df <- summarizeCooccur(mock_cor_matrix)

  # find the row for SigA vs SigB (0.9)
  row_pos <- summary_df[summary_df$Correlation == 0.9, ]
  testthat::expect_equal(row_pos$Type, "Co-occurrence")

  # find the row for SigA vs SigC (-0.8)
  row_neg <- summary_df[summary_df$Correlation == -0.8, ]
  testthat::expect_equal(row_neg$Type, "Mutual Exclusivity")
})

test_that("Happy Path: Sorting is correct (by absolute strength)", {

  summary_df <- summarizeCooccur(mock_cor_matrix)

  # the first row should be the strongest correlation
  testthat::expect_equal(summary_df$Correlation[1], 0.9)

  # the second strongest
  testthat::expect_equal(summary_df$Correlation[2], -0.8)
})

test_that("Happy Path: top_n parameter works", {

  # Request only top 1 positive and top 1 negative
  summary_small <- summarizeCooccur(mock_cor_matrix, top_n = 1)

  testthat::expect_equal(nrow(summary_small), 2)

  testthat::expect_true(0.9 %in% summary_small$Correlation)  # Strongest Pos
  testthat::expect_true(-0.8 %in% summary_small$Correlation) # Strongest Neg
})

test_that("Sad Path: Input validation throws correct errors", {

  # not a matrix
  testthat::expect_error(
    summarizeCooccur(as.data.frame(mock_cor_matrix)),
    regexp = "must be a numeric matrix"
  )

  # missing rownames/colnames
  mock_no_names <- mock_cor_matrix
  rownames(mock_no_names) <- NULL
  testthat::expect_error(
    summarizeCooccur(mock_no_names),
    regexp = "must have rownames and colnames"
  )
})
