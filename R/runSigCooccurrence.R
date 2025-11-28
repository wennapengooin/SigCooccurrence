#' Launch SigCooccurrence Shiny App
#'
#' Launches the interactive Shiny application for visualizing signature
#' co-occurrence.
#'
#' @details
#' The Shiny app allows users to upload a co-occurrence matrix (as a CSV file)
#' or load built-in demo data. Users can interactively cluster rows and
#' columns and download the resulting heatmap.
#'
#' @return None. The function opens a Shiny app window.
#'
#' @importFrom shiny runApp
#' @export
#'
#' @examples
#' \dontrun{
#'   runSigCooccurrence()
#' }
runSigCooccurrence <- function() {
  appDir <- system.file("shiny", package = "SigCooccurrence")

  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `SigCooccurrence`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}


# [END]
