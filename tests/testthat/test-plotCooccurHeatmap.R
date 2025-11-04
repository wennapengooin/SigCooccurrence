# Test script for plotCooccurHeatmap()

# --- 1. Setup ---
# Create a mock co-occurrence matrix
mock_cooccur <- matrix(
  c(
    1.0, 0.8, -0.2, 0.1,
    0.8, 1.0, -0.1, 0.0,
    -0.2, -0.1, 1.0, 0.9,
    0.1, 0.0, 0.9, 1.0
  ),
  nrow = 4,
  dimnames = list(c("SigA", "SigB", "SigC", "SigD"),
                  c("SigA", "SigB", "SigC", "SigD"))
)


# --- 2. Tests ---

test_that("Happy Path: Plot is generated correctly", {
  testthat::skip_if_not_installed("pheatmap")

  # Suppress the "Generating heatmap plot..." message
  # and capture the plot object
  p <- suppressMessages(
    plotCooccurHeatmap(mock_cooccur)
  )

  # Check 1: Is it a pheatmap object?
  testthat::expect_s3_class(p, "pheatmap")

  # Check 2: Does the tree_row (dendrogram) exist?
  testthat::expect_true(!is.null(p$tree_row))

  # Check 3: Does the tree_col (dendrogram) exist?
  testthat::expect_true(!is.null(p$tree_col))

  # Check 4: Is the title correct?
  testthat::expect_equal(p$gtable$grobs[[1]]$label, "Signature Co-occurrence Heatmap")
})


test_that("Happy Path: Clustering parameters work", {
  testthat::skip_if_not_installed("pheatmap")

  # Test cluster_rows = FALSE
  p_no_row <- suppressMessages(
    plotCooccurHeatmap(mock_cooccur, cluster_rows = FALSE, cluster_cols = TRUE)
  )
  # The row dendrogram should be NA (not NULL)
  testthat::expect_true(is.na(p_no_row$tree_row))

  # Test cluster_cols = FALSE
  p_no_col <- suppressMessages(
    plotCooccurHeatmap(mock_cooccur, cluster_rows = TRUE, cluster_cols = FALSE)
  )
  # The col dendrogram should be NA (not NULL)
  testthat::expect_true(is.na(p_no_col$tree_col))
})


test_that("Sad Path: Input validation throws correct errors", {
  testthat::skip_if_not_installed("pheatmap")

  # Test for non-matrix input
  testthat::expect_error(
    plotCooccurHeatmap(as.data.frame(mock_cooccur)),
    regexp = "`cooccur_matrix` must be a square numeric matrix."
  )

  # Test for non-square matrix
  testthat::expect_error(
    plotCooccurHeatmap(mock_cooccur[1:2, ]),
    regexp = "`cooccur_matrix` must be a square numeric matrix."
  )

  # Test for missing rownames
  mock_no_rownames <- mock_cooccur
  rownames(mock_no_rownames) <- NULL
  testthat::expect_error(
    plotCooccurHeatmap(mock_no_rownames),
    regexp = "`cooccur_matrix` must have both rownames and colnames."
  )

  # Test for missing colnames
  mock_no_colnames <- mock_cooccur
  colnames(mock_no_colnames) <- NULL
  testthat::expect_error(
    plotCooccurHeatmap(mock_no_colnames),
    regexp = "`cooccur_matrix` must have both rownames and colnames."
  )
})

