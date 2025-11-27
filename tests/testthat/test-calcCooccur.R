# Test script for calcCooccur()

# --- Setup -----------------
mock_exposures <- matrix(
  c(
    1, 2, 3, 4, 5,
    2, 4, 6, 8, 10,
    5, 4, 3, 2, 1,
    1, 5, 2, 4, 3,
    10, 10, 10, 10, 10
  ),
  nrow = 5,
  byrow = TRUE,
  dimnames = list(
    c("SigA", "SigB", "SigC", "SigD", "SigE"),
    paste0("Sample", 1:5)
  )
)

# --- Tests -----------------

test_that("Happy Path: Spearman correlation works correctly", {

  # run function
  cooccur_mat <- suppressWarnings(
    calcCooccur(mock_exposures, method = "spearman")
  )

  # check object structure
  testthat::expect_true(is.matrix(cooccur_mat))
  testthat::expect_equal(nrow(cooccur_mat), 5)
  testthat::expect_equal(ncol(cooccur_mat), 5)
  testthat::expect_equal(rownames(cooccur_mat), rownames(mock_exposures))

  # check correlations
  testthat::expect_equal(cooccur_mat["SigA", "SigB"], 1)

  # should be -1
  testthat::expect_equal(cooccur_mat["SigA", "SigC"], -1)

  # should be NA
  testthat::expect_true(is.na(cooccur_mat["SigA", "SigE"]))
})

test_that("Happy Path: Pearson correlation works correctly", {

  cooccur_mat <- suppressWarnings(
    calcCooccur(mock_exposures, method = "pearson")
  )

  # should be 1
  testthat::expect_equal(cooccur_mat["SigA", "SigB"], 1)

  # should be -1
  testthat::expect_equal(cooccur_mat["SigA", "SigC"], -1)

  # should be 0.3
  expected_cor <- stats::cor(c(1,2,3,4,5), c(1,5,2,4,3))
  testthat::expect_equal(cooccur_mat["SigA", "SigD"], expected_cor)
})

test_that("Happy Path: Default method is Spearman", {

  # run without specifying method
  cooccur_default <- suppressWarnings(calcCooccur(mock_exposures))

  # run explicit Spearman
  cooccur_spearman <- suppressWarnings(
    calcCooccur(mock_exposures, method = "spearman")
  )

  # should be identical
  testthat::expect_equal(cooccur_default, cooccur_spearman)
})

test_that("Sad Path: Input validation throws correct errors", {

  # not a matrix
  testthat::expect_error(
    calcCooccur(as.data.frame(mock_exposures)),
    regexp = "`exposure_matrix` must be a numeric matrix"
  )

  # non-numeric matrix
  char_matrix <- matrix(as.character(mock_exposures), nrow=5)
  testthat::expect_error(
    calcCooccur(char_matrix),
    regexp = "`exposure_matrix` must be a numeric matrix"
  )

  # missing rownames
  unnamed_matrix <- mock_exposures
  rownames(unnamed_matrix) <- NULL
  testthat::expect_error(
    calcCooccur(unnamed_matrix),
    regexp = "must have rownames"
  )

  # not enough signatures (rows)
  testthat::expect_error(
    calcCooccur(mock_exposures[1, , drop=FALSE]),
    regexp = "must have at least 2 signatures"
  )

  # not enough samples (columns)
  testthat::expect_error(
    calcCooccur(mock_exposures[, 1, drop=FALSE]),
    regexp = "must have at least 2 samples"
  )

  # invalid method
  testthat::expect_error(
    calcCooccur(mock_exposures, method = "kendall"),
    regexp = "'arg' should be one of"
  )
})

test_that("Warning: Zero variance rows trigger a warning", {

  # expect a warning because zero variance
  testthat::expect_warning(
    calcCooccur(mock_exposures),
    regexp = "The following signatures have zero variance"
  )

  # warning should specifically mention "SigE"
  testthat::expect_warning(
    calcCooccur(mock_exposures),
    regexp = "SigE"
  )
})
