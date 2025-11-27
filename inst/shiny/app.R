# Load required packages
library(shiny)
library(SigCooccurrence)
library(pheatmap)

ui <- fluidPage(

  titlePanel("SigCooccurrence Explorer"),

  sidebarLayout(

    sidebarPanel(

      helpText("Analyze and visualize mutational signature co-occurrence patterns."),

      # select file
      fileInput("file1", "Upload Co-occurrence Matrix (.csv)",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),

      actionButton("load_demo", "Load Demo Data"),

      tags$hr(),

      # clustering options
      checkboxInput("cluster_rows", "Cluster Rows", TRUE),
      checkboxInput("cluster_cols", "Cluster Columns", TRUE),

      tags$hr(),

      # download button
      downloadButton("downloadPlot", "Download Heatmap")
    ),

    # panel displaying outputs
    mainPanel(

      plotOutput("heatmapPlot", height = "600px"),

      tags$h4("Data Preview"),
      tableOutput("contents")
    )
  )
)

server <- function(input, output, session) {

  values <- reactiveValues(matrix_data = NULL)

  observeEvent(input$file1, {
    req(input$file1)
    tryCatch(
      {
        # read csv
        df <- read.csv(input$file1$datapath, row.names = 1, check.names = FALSE)
        values$matrix_data <- as.matrix(df)
      },
      error = function(e) {
        stop(safeError(e))
      }
    )
  })

  observeEvent(input$load_demo, {

    demo_path <- system.file("extdata", "demo_cooccur.csv", package = "SigCooccurrence")

    if (demo_path != "" && file.exists(demo_path)) {
      df <- read.csv(demo_path, row.names = 1, check.names = FALSE)
      values$matrix_data <- as.matrix(df)
      showNotification("Demo data loaded successfully!", type = "message")
    } else {
      local_path <- "../../../inst/extdata/demo_cooccur.csv"
      if(file.exists(local_path)) {
        df <- read.csv(local_path, row.names = 1, check.names = FALSE)
        values$matrix_data <- as.matrix(df)
        showNotification("Demo data loaded (local fallback).", type = "message")
      } else {
        showNotification("Demo data file not found.", type = "error")
      }
    }
  })

  # render heatmap
  output$heatmapPlot <- renderPlot({
    req(values$matrix_data)

    plotCooccurHeatmap(
      cooccur_matrix = values$matrix_data,
      cluster_rows = input$cluster_rows,
      cluster_cols = input$cluster_cols,
      title = "Signature Co-occurrence"
    )
  })

  output$contents <- renderTable({
    req(values$matrix_data)
    head(values$matrix_data[, 1:min(5, ncol(values$matrix_data))], 5)
  }, rownames = TRUE)

  output$downloadPlot <- downloadHandler(
    filename = function() { paste("cooccurrence_heatmap", ".png", sep = "") },
    content = function(file) {
      req(values$matrix_data)

      p <- plotCooccurHeatmap(
        cooccur_matrix = values$matrix_data,
        cluster_rows = input$cluster_rows,
        cluster_cols = input$cluster_cols,
        title = "Signature Co-occurrence"
      )

      # save as png
      png(file, width = 800, height = 700)
      print(p)
      dev.off()
    }
  )
}

shinyApp(ui = ui, server = server)
