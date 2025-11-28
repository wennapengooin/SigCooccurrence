# Load required packages
library(shiny)
library(SigCooccurrence)
library(pheatmap)
library(dplyr)

ui <- fluidPage(

  titlePanel("SigCooccurrence Explorer"),

  sidebarLayout(

    sidebarPanel(

      helpText("Analyze and visualize mutational signature co-occurrence patterns."),

      fileInput("file1", "Upload Co-occurrence Matrix (.csv)",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),

      actionButton("load_demo", "Load Demo Data"),
      tags$hr(),

      tags$h4("Heatmap Settings"),
      checkboxInput("cluster_rows", "Cluster Rows", TRUE),
      checkboxInput("cluster_cols", "Cluster Columns", TRUE),
      tags$hr(),

      tags$h4("Summary Settings"),
      numericInput("top_n", "Number of Top Pairs to Show:",
                   value = 5, min = 1, max = 50),
      tags$hr(),

      downloadButton("downloadPlot", "Download Heatmap")
    ),

    mainPanel(

      tabsetPanel(type = "tabs",
                  tabPanel("Heatmap",
                           plotOutput("heatmapPlot", height = "600px")
                  ),
                  # UPDATED: Three separate summary tabs
                  tabPanel("Top Co-occurring",
                           tags$h4("Top Co-occurring Signature Pairs"),
                           tableOutput("table_cooccur")
                  ),
                  tabPanel("Top Mutually Exclusive",
                           tags$h4("Top Mutually Exclusive Signature Pairs"),
                           tableOutput("table_mutex")
                  ),
                  tabPanel("No Correlation",
                           tags$h4("Pairs with Near-Zero Correlation (-0.1 to 0.1)"),
                           tableOutput("table_no_cor")
                  ),
                  tabPanel("Raw Data",
                           tags$h4("Input Matrix Preview"),
                           tableOutput("contents")
                  )
      )
    )
  )
)

server <- function(input, output, session) {

  values <- reactiveValues(matrix_data = NULL)

  observeEvent(input$file1, {
    req(input$file1)
    tryCatch(
      {
        df <- read.csv(input$file1$datapath, row.names = 1, check.names = FALSE)
        values$matrix_data <- as.matrix(df)
      },
      error = function(e) { stop(safeError(e)) }
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

  summary_data <- reactive({
    req(values$matrix_data)

    SigCooccurrence::summarizeCooccur(values$matrix_data, top_n = 10)
  })

  output$heatmapPlot <- renderPlot({
    req(values$matrix_data)
    plotCooccurHeatmap(
      cooccur_matrix = values$matrix_data,
      cluster_rows = input$cluster_rows,
      cluster_cols = input$cluster_cols,
      title = "Signature Co-occurrence"
    )
  })

  output$table_cooccur <- renderTable({
    req(summary_data())
    summary_data() %>%
      dplyr::filter(Type == "Co-occurrence") %>%
      head(input$top_n)
  })

  output$table_mutex <- renderTable({
    req(summary_data())
    summary_data() %>%
      dplyr::filter(Type == "Mutual Exclusivity") %>%
      head(input$top_n)
  })

  output$table_no_cor <- renderTable({
    req(summary_data())
    summary_data() %>%
      dplyr::filter(Correlation > -0.1 & Correlation < 0.1) %>%
      head(input$top_n) %>%
      dplyr::mutate(Type = "No Correlation")
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
      png(file, width = 800, height = 700)
      print(p)
      dev.off()
    }
  )
}

shinyApp(ui = ui, server = server)
