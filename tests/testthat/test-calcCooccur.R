# --- Setup -----------------
get_mock_exposures <- function() {
  matrix(
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
}

# --- Tests -----------------

test_that("Happy Path: Spearman correlation works correctly", {
  mock_exposures <- get_mock_exposures()

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
  testthat::expect_equal(cooccur_mat["SigA", "SigC"], -1)

  # actual value is 0.3
  testthat::expect_equal(round(cooccur_mat["SigA", "SigD"], 1), 0.3)

  testthat::expect_true(is.na(cooccur_mat["SigA", "SigE"]))
})

test_that("Happy Path: Pearson correlation works correctly", {
  mock_exposures <- get_mock_exposures()

  cooccur_mat <- suppressWarnings(
    calcCooccur(mock_exposures, method = "pearson")
  )

  testthat::expect_equal(cooccur_mat["SigA", "SigB"], 1)
  testthat::expect_equal(cooccur_mat["SigA", "SigC"], -1)

  expected_cor <- stats::cor(c(1,2,3,4,5), c(1,5,2,4,3))
  testthat::expect_equal(cooccur_mat["SigA", "SigD"], expected_cor)
})

test_that("Happy Path: Default method is Spearman", {
  mock_exposures <- get_mock_exposures()

  cooccur_default <- suppressWarnings(calcCooccur(mock_exposures))
  cooccur_spearman <- suppressWarnings(
    calcCooccur(mock_exposures, method = "spearman")
  )

  testthat::expect_equal(cooccur_default, cooccur_spearman)
})

test_that("Happy Path: CSV saving works", {
  mock_exposures <- get_mock_exposures()

  # create a temp file path
  tmp_csv <- tempfile(fileext = ".csv")

  # run with save_csv = TRUE
  suppressWarnings(
    suppressMessages(
      calcCooccur(mock_exposures, save_csv = TRUE, csv_path = tmp_csv)
    )
  )

  # check if file exists
  testthat::expect_true(file.exists(tmp_csv))

  # check if we can read it back
  saved_data <- read.csv(tmp_csv, row.names = 1)
  testthat::expect_equal(nrow(saved_data), 5)
  testthat::expect_equal(ncol(saved_data), 5)

  # cleanup
  unlink(tmp_csv)
})

test_that("Sad Path: Input validation throws correct errors", {
  mock_exposures <- get_mock_exposures()

  # not a matrix
  testthat::expect_error(
    calcCooccur(as.data.frame(mock_exposures)),
    regexp = "`exposure_matrix` must be a numeric matrix"
  )

  # missing rownames
  unnamed_matrix <- mock_exposures
  rownames(unnamed_matrix) <- NULL
  testthat::expect_error(
    calcCooccur(unnamed_matrix),
    regexp = "must have rownames"
  )

  # CSV Path missing
  testthat::expect_error(
    calcCooccur(mock_exposures, save_csv = TRUE),
    regexp = "`csv_path` must be provided"
  )

  # CSV Directory doesn't exist
  testthat::expect_error(
    calcCooccur(mock_exposures, save_csv = TRUE, csv_path = "fake_dir/test.csv"),
    regexp = "Directory 'fake_dir' does not exist"
  )
})

test_that("Warning: Zero variance rows trigger a warning", {
  mock_exposures <- get_mock_exposures()

  testthat::expect_warning(
    calcCooccur(mock_exposures),
    regexp = "The following signatures have zero variance"
  )
})
