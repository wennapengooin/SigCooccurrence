# Load required packages
library(shiny)
library(SigCooccurrence)
library(pheatmap)

# Define UI
ui <- fluidPage(

  # App title
  titlePanel("SigCooccurrence Explorer"),

  # Sidebar layout with input and output definitions
  sidebarLayout(

    # Sidebar panel for inputs
    sidebarPanel(

      # Introduction text
      helpText("Analyze and visualize mutational signature co-occurrence patterns."),

      # Input: Select a file
      fileInput("file1", "Upload Co-occurrence Matrix (.csv)",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),

      # Input: Demo data button
      actionButton("load_demo", "Load Demo Data"),

      # Horizontal line
      tags$hr(),

      # Input: Clustering options
      checkboxInput("cluster_rows", "Cluster Rows", TRUE),
      checkboxInput("cluster_cols", "Cluster Columns", TRUE),

      # Horizontal line
      tags$hr(),

      # Download button
      downloadButton("downloadPlot", "Download Heatmap")
    ),

    # Main panel for displaying outputs
    mainPanel(

      # Output: Plot
      plotOutput("heatmapPlot", height = "600px"),

      # Output: Data Table (preview)
      tags$h4("Data Preview"),
      tableOutput("contents")
    )
  )
)

# Define server logic
server <- function(input, output, session) {

  # Reactive value to store the data
  values <- reactiveValues(matrix_data = NULL)

  # Observer: Handle File Upload
  observeEvent(input$file1, {
    req(input$file1)
    tryCatch(
      {
        # Read the CSV file
        df <- read.csv(input$file1$datapath, row.names = 1, check.names = FALSE)
        values$matrix_data <- as.matrix(df)
      },
      error = function(e) {
        stop(safeError(e))
      }
    )
  })

  # Observer: Handle Demo Data Load
  observeEvent(input$load_demo, {
    # UPDATED: Locate demo data in inst/extdata
    # When installed, this maps to the package's 'extdata' folder
    demo_path <- system.file("extdata", "demo_cooccur.csv", package = "SigCooccurrence")

    if (demo_path != "" && file.exists(demo_path)) {
      df <- read.csv(demo_path, row.names = 1, check.names = FALSE)
      values$matrix_data <- as.matrix(df)
      showNotification("Demo data loaded successfully!", type = "message")
    } else {
      # Fallback for local testing if package isn't installed yet
      # Assuming app is run from inst/shiny/SigCooccurrenceApp/
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

  # Output: Render the Heatmap
  output$heatmapPlot <- renderPlot({
    req(values$matrix_data)

    # Use your package's function!
    plotCooccurHeatmap(
      cooccur_matrix = values$matrix_data,
      cluster_rows = input$cluster_rows,
      cluster_cols = input$cluster_cols,
      title = "Signature Co-occurrence"
    )
  })

  # Output: Data Preview
  output$contents <- renderTable({
    req(values$matrix_data)
    # Show first 5 rows/cols
    head(values$matrix_data[, 1:min(5, ncol(values$matrix_data))], 5)
  }, rownames = TRUE)

  # Handler: Download Plot
  output$downloadPlot <- downloadHandler(
    filename = function() { paste("cooccurrence_heatmap", ".png", sep = "") },
    content = function(file) {
      req(values$matrix_data)

      # Generate the pheatmap object
      p <- plotCooccurHeatmap(
        cooccur_matrix = values$matrix_data,
        cluster_rows = input$cluster_rows,
        cluster_cols = input$cluster_cols,
        title = "Signature Co-occurrence"
      )

      # Save using png (pheatmap prints directly)
      png(file, width = 800, height = 700)
      print(p)
      dev.off()
    }
  )
}

# Create Shiny app
shinyApp(ui = ui, server = server)
