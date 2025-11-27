# --- Setup -----------------
# create a mock co-occurrence matrix
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


# --- Tests -----------------

test_that("Happy Path: Plot is generated correctly", {
  testthat::skip_if_not_installed("pheatmap")

  p <- suppressMessages(
    plotCooccurHeatmap(mock_cooccur)
  )

  testthat::expect_s3_class(p, "pheatmap")

  testthat::expect_true(!is.null(p$tree_row))

  testthat::expect_true(!is.null(p$tree_col))

  testthat::expect_equal(p$gtable$grobs[[1]]$label, "Signature Co-occurrence Heatmap")
})


test_that("Happy Path: Clustering parameters work", {
  testthat::skip_if_not_installed("pheatmap")

  # test cluster_rows = FALSE
  p_no_row <- suppressMessages(
    plotCooccurHeatmap(mock_cooccur, cluster_rows = FALSE, cluster_cols = TRUE)
  )
  # the row dendrogram should be NA
  testthat::expect_true(is.na(p_no_row$tree_row))

  # test cluster_cols = FALSE
  p_no_col <- suppressMessages(
    plotCooccurHeatmap(mock_cooccur, cluster_rows = TRUE, cluster_cols = FALSE)
  )
  # the col dendrogram should be NA
  testthat::expect_true(is.na(p_no_col$tree_col))
})


test_that("Sad Path: Input validation throws correct errors", {
  testthat::skip_if_not_installed("pheatmap")

  # test for non-matrix input
  testthat::expect_error(
    plotCooccurHeatmap(as.data.frame(mock_cooccur)),
    regexp = "`cooccur_matrix` must be a square numeric matrix."
  )

  # test for non-square matrix
  testthat::expect_error(
    plotCooccurHeatmap(mock_cooccur[1:2, ]),
    regexp = "`cooccur_matrix` must be a square numeric matrix."
  )

  # test for missing rownames
  mock_no_rownames <- mock_cooccur
  rownames(mock_no_rownames) <- NULL
  testthat::expect_error(
    plotCooccurHeatmap(mock_no_rownames),
    regexp = "`cooccur_matrix` must have both rownames and colnames."
  )

  # test for missing colnames
  mock_no_colnames <- mock_cooccur
  colnames(mock_no_colnames) <- NULL
  testthat::expect_error(
    plotCooccurHeatmap(mock_no_colnames),
    regexp = "`cooccur_matrix` must have both rownames and colnames."
  )
})

