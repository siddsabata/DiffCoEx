# Jason Hyun (jasonhyu)
# Siddharth Sabata (ssabata)
# Darrick Lo (ddlo)
# Katie Wang (kcw2)

# Dec 1, 2024

# NOTE: Generative AI used to produce following code:

library(shiny)

ui <- fluidPage(
  titlePanel("Gene Expression Data Preprocessing and Clustering"),
  
  sidebarLayout(
    sidebarPanel(
      # Dataset selection
      radioButtons("dataset", "Choose Dataset:",
                   choices = list("Rat Data (GDS2901)" = "rat",
                                "Golub Data" = "golub"),
                   selected = "rat"),
      
      # File path inputs
      conditionalPanel(
        condition = "input.dataset == 'rat'",
        textInput("rat_path", "GDS2901.soft file path:", 
                  value = "data/GDS2901.soft")
      ),
      conditionalPanel(
        condition = "input.dataset == 'golub'",
        textInput("golub_path", "Golub.txt file path:", 
                  value = "data/golub.txt")
      ),
      
      # Process, Cluster, and Test buttons
      actionButton("process", "Process Data", 
                   class = "btn-primary"),
      br(), br(),
      actionButton("cluster", "Perform Clustering",
                   class = "btn-success"),
      br(), br(),
      actionButton("test", "Perform Significance Testing",
                   class = "btn-info"),
      
      # Add module selection and plot button
      conditionalPanel(
        condition = "input.tabset == 'Significance Testing'",  # Only show when in Significance Testing tab
        hr(),
        selectInput("module", "Select Module to Plot:", choices = NULL),
        actionButton("plot", "Generate Module Plots", class = "btn-info"),
        hr(),
        actionButton("heatmap", "Generate Correlation Heatmap", class = "btn-warning")
      )
    ),
    
    mainPanel(
      tabsetPanel(id = "tabset",  # Add id for conditional panel
        # Existing tabs
        tabPanel("Processing", 
          h4("Processing Status:"),
          verbatimTextOutput("status"),
          verbatimTextOutput("log")
        ),
        
        tabPanel("Clustering", 
          h4("Clustering Results:"),
          verbatimTextOutput("clustering_info"),
          h4("Clustering Visualizations:"),
          plotOutput("cluster_comparison_plot"),
          plotOutput("diffcoex_distribution_plot"),
          plotOutput("coxpress_distribution_plot"),
          plotOutput("cluster_overlap_plot")
        ),
        
        # New significance testing tab
        tabPanel("Significance Testing",
          h4("Significance Testing Results"),
          verbatimTextOutput("sig_testing_status"),
          hr(),
          h4("Module Correlation Results:"),
          tableOutput("module_correlation_table"),
          hr(),
          h4("Null Distribution Results:"),
          tableOutput("null_distribution_table"),
          hr(),
          h4("Module Distribution Plots:"),
          div(style = "text-align: center;",
              uiOutput("condition1_plot")),
          div(style = "text-align: center;",
              uiOutput("condition2_plot")),
          hr(),
          h4("Correlation Heatmap:"),
          div(style = "text-align: center;",
              uiOutput("heatmap_plot"))
        )
      )
    )
  )
)

server <- function(input, output, session) {
  # Create reactive values to store processing status
  values <- reactiveValues(
    status = "",
    log = "",
    preprocessing_complete = FALSE,
    clustering_complete = FALSE,
    testing_complete = FALSE,
    testing_status = "",
    plotting_status = "",
    clustering_plots = NULL,
    clustering_info = NULL,
    testing_results = NULL,
    module_correlation_data = NULL,
    null_distribution_data = NULL,
    available_modules = NULL,
    heatmap_status = ""  # Add new reactive value for heatmap status
  )
  
  # Update module choices after significance testing
  observe({
    req(values$module_correlation_data)
    modules <- unique(values$module_correlation_data$Module)
    updateSelectInput(session, "module",
                     choices = modules)
  })
  
  # Process button handler
  observeEvent(input$process, {
    values$status <- "Processing..."
    values$preprocessing_complete <- FALSE
    values$clustering_complete <- FALSE
    
    # Determine which file path to use
    file_path <- if(input$dataset == "rat") {
      input$rat_path
    } else {
      input$golub_path
    }
    
    # Run preprocessing
    tryCatch({
      result <- system2("./preprocess",
                        args = c(input$dataset, file_path),
                        stdout = TRUE,
                        stderr = TRUE,
                        wait = TRUE)
      
      values$status <- "Processing complete! Ready for clustering."
      values$log <- paste(result, collapse = "\n")
      values$preprocessing_complete <- TRUE
      
    }, error = function(e) {
      values$status <- "Error!"
      values$log <- paste("Error processing data:", e$message)
      values$preprocessing_complete <- FALSE
    })
  })
  
  # Cluster button handler
  observeEvent(input$cluster, {
    values$status <- "Clustering..."
    
    tryCatch({
      # Get the appropriate file paths based on dataset
      if (input$dataset == "rat") {
        data_paths <- list(
          diffcoex = list(
            condition1 = "output/diffcoex/rat_wild_types.csv",
            condition2 = "output/diffcoex/rat_eker_mutants.csv"
          ),
          coxpress = list(
            condition1 = "output/coxpress/rat_wild_types.csv",
            condition2 = "output/coxpress/rat_eker_mutants.csv"
          )
        )
      } else {
        data_paths <- list(
          diffcoex = list(
            condition1 = "output/diffcoex/golub_ALL_samples.csv",
            condition2 = "output/diffcoex/golub_AML_samples.csv"
          ),
          coxpress = list(
            condition1 = "output/coxpress/golub_ALL_samples.csv",
            condition2 = "output/coxpress/golub_AML_samples.csv"
          )
        )
      }
      
      # Source the clustering script and perform clustering
      source("clustering.R")
      results <- performClustering(data_paths)
      
      # Store results in reactive values
      values$clustering_complete <- TRUE
      values$clustering_plots <- results$plots
      values$clustering_info <- paste(
        "DiffCoEx Summary:",
        "Total Genes:", results$diffcoex_summary$Total_Genes,
        "\nTotal Clusters:", results$diffcoex_summary$Total_Clusters,
        "\nGenes Per Cluster Summary:",
        "\n  Min:", results$diffcoex_summary$Genes_Per_Cluster[1],
        "\n  1st Qu:", results$diffcoex_summary$Genes_Per_Cluster[2],
        "\n  Median:", results$diffcoex_summary$Genes_Per_Cluster[3],
        "\n  Mean:", round(results$diffcoex_summary$Genes_Per_Cluster[4], 2),
        "\n  3rd Qu:", results$diffcoex_summary$Genes_Per_Cluster[5],
        "\n  Max:", results$diffcoex_summary$Genes_Per_Cluster[6],
        "\n\nCoXpress Summary:",
        "\nTotal Genes:", results$coxpress_summary$Total_Genes,
        "\nTotal Clusters:", results$coxpress_summary$Total_Clusters,
        "\nGenes Per Cluster Summary:",
        "\n  Min:", results$coxpress_summary$Genes_Per_Cluster[1],
        "\n  1st Qu:", results$coxpress_summary$Genes_Per_Cluster[2],
        "\n  Median:", results$coxpress_summary$Genes_Per_Cluster[3],
        "\n  Mean:", round(results$coxpress_summary$Genes_Per_Cluster[4], 2),
        "\n  3rd Qu:", results$coxpress_summary$Genes_Per_Cluster[5],
        "\n  Max:", results$coxpress_summary$Genes_Per_Cluster[6],
        sep = " "
      )
      values$status <- "Clustering complete!"
      
    }, error = function(e) {
      values$status <- "Clustering Error!"
      values$log <- paste("Error during clustering:", e$message)
      values$clustering_complete <- FALSE
    })
  })
  
  # Significance testing handler
  observeEvent(input$test, {
    req(values$clustering_complete)
    
    # Immediately show that testing is starting
    values$testing_status <- "Running significance testing... (this may take a few minutes)"
    values$plotting_status <- ""  # Clear any previous plotting status
    
    # Force the UI to update immediately
    invalidateLater(0)
    
    tryCatch({
      if(input$dataset == "rat") {
        module_map <- file.path(getwd(), "data/rat/rat_diffcoex.csv")
        condition1_data <- file.path(getwd(), "output/diffcoex/rat_wild_types.csv")
        condition2_data <- file.path(getwd(), "output/diffcoex/rat_eker_mutants.csv")
      } else {
        module_map <- file.path(getwd(), "data/golub/golub_diffcoex.csv")
        condition1_data <- file.path(getwd(), "output/diffcoex/golub_ALL_samples.csv")
        condition2_data <- file.path(getwd(), "output/diffcoex/golub_AML_samples.csv")
      }
      
      executable_path <- "significanceTesting/significanceTesting"
      
      result <- system2(executable_path,
                       args = c(module_map, condition1_data, condition2_data),
                       stdout = TRUE,
                       stderr = TRUE,
                       wait = TRUE)
      
      # Read results after execution
      values$module_correlation_data <- read.csv("output/sigTesting/module_correlation_results.csv")
      values$null_distribution_data <- read.csv("output/sigTesting/null_distribution_results.csv")
      
      values$testing_complete <- TRUE
      values$testing_status <- "Significance testing complete!"
      
    }, error = function(e) {
      values$testing_status <- "Testing Error!"
      values$log <- paste("Error during significance testing:", e$message)
      values$testing_complete <- FALSE
    })
  })
  
  # Plot handler
  observeEvent(input$plot, {
    req(values$testing_complete, input$module)
    values$plotting_status <- paste("Generating plots for module", input$module, "...")
    
    tryCatch({
      if(input$dataset == "rat") {
        module_map <- file.path(getwd(), "data/rat/rat_diffcoex.csv")
        condition1_data <- file.path(getwd(), "output/diffcoex/rat_wild_types.csv")
        condition2_data <- file.path(getwd(), "output/diffcoex/rat_eker_mutants.csv")
      } else {
        module_map <- file.path(getwd(), "data/golub/golub_diffcoex.csv")
        condition1_data <- file.path(getwd(), "output/diffcoex/golub_ALL_samples.csv")
        condition2_data <- file.path(getwd(), "output/diffcoex/golub_AML_samples.csv")
      }
      
      # Run plotting executable
      executable_path <- "plotSignificanceTesting/plotSignificanceTesting"
      result <- system2(executable_path,
                       args = c(module_map, 
                              condition1_data, 
                              condition2_data,
                              input$module),
                       stdout = TRUE,
                       stderr = TRUE,
                       wait = TRUE)
      
      # Check if files were created
      condition1_file <- file.path("output", "plotting", 
                                 paste0(input$module, "_condition1_distribution.png"))
      condition2_file <- file.path("output", "plotting", 
                                 paste0(input$module, "_condition2_distribution.png"))
      
      if (file.exists(condition1_file) && file.exists(condition2_file)) {
        values$plotting_status <- paste("Plots generated for module", input$module)
      } else {
        values$plotting_status <- paste("Error: Plot files not found for module", input$module)
      }
      
      # Print debugging information
      print(paste("Result from executable:", paste(result, collapse = "\n")))
      print(paste("Condition 1 file exists:", file.exists(condition1_file)))
      print(paste("Condition 2 file exists:", file.exists(condition2_file)))
      
    }, error = function(e) {
      values$plotting_status <- paste("Error generating plots:", e$message)
    })
  })
  
  # Add heatmap handler
  observeEvent(input$heatmap, {
    req(values$testing_complete)
    values$heatmap_status <- "Generating correlation heatmap..."
    
    tryCatch({
      if(input$dataset == "rat") {
        condition1_data <- file.path(getwd(), "output/diffcoex/rat_wild_types.csv")
        condition2_data <- file.path(getwd(), "output/diffcoex/rat_eker_mutants.csv")
      } else {
        condition1_data <- file.path(getwd(), "output/diffcoex/golub_ALL_samples.csv")
        condition2_data <- file.path(getwd(), "output/diffcoex/golub_AML_samples.csv")
      }
      
      # Run heatmap executable
      executable_path <- "correlationHeatmap/correlationHeatmap"
      result <- system2(executable_path,
                       args = c(condition1_data, 
                              condition2_data),
                       stdout = TRUE,
                       stderr = TRUE,
                       wait = TRUE)
      
      values$heatmap_status <- "Correlation heatmap generated!"
      
    }, error = function(e) {
      values$heatmap_status <- paste("Error generating heatmap:", e$message)
    })
  })
  
  # Add heatmap output
  output$heatmap_plot <- renderUI({
    req(values$testing_complete, input$heatmap)
    
    # Get absolute path to the image
    filename <- normalizePath(file.path("output", "plotting", "heatmap_merged.png"),
                            mustWork = FALSE)
    
    # Check if file exists
    if (!file.exists(filename)) {
      return(p("Heatmap not available"))
    }
    
    # Add timestamp to force image refresh
    timestamp <- format(Sys.time(), "%Y%m%d%H%M%OS6")
    
    tags$img(
      src = paste0("plotting/heatmap_merged.png?t=", timestamp),
      style = "max-width: 100%; height: auto;",
      alt = "Correlation Heatmap"
    )
  })
  
  # Update combined status output to include heatmap status
  output$sig_testing_status <- renderText({
    status_text <- values$testing_status
    
    if (values$plotting_status != "") {
      if (status_text != "") {
        status_text <- paste(status_text, "\n\n", values$plotting_status)
      } else {
        status_text <- values$plotting_status
      }
    }
    
    if (values$heatmap_status != "") {
      if (status_text != "") {
        status_text <- paste(status_text, "\n\n", values$heatmap_status)
      } else {
        status_text <- values$heatmap_status
      }
    }
    
    status_text
  })
  
  # Display status and log
  output$status <- renderText({
    values$status
  })
  
  output$log <- renderText({
    values$log
  })
  
  # Display clustering results
  output$clustering_info <- renderText({
    req(values$clustering_complete)
    values$clustering_info
  })
  
  output$cluster_comparison_plot <- renderPlot({
    req(values$clustering_complete)
    values$clustering_plots$cluster_comparison
  })
  
  output$diffcoex_distribution_plot <- renderPlot({
    req(values$clustering_complete)
    ggplot(values$clustering_plots$diffcoex_distribution, 
           aes(x = Cluster, y = Gene_Count)) +
      geom_bar(stat = "identity", fill = "blue", alpha = 0.7) +
      labs(title = "Number of Genes Per Cluster - DiffCoEx", 
           x = "Cluster", y = "Number of Genes") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            text = element_text(size = 14))
  })
  
  output$coxpress_distribution_plot <- renderPlot({
    req(values$clustering_complete)
    ggplot(values$clustering_plots$coxpress_distribution, 
           aes(x = Cluster, y = Gene_Count)) +
      geom_bar(stat = "identity", fill = "orange", alpha = 0.7) +
      labs(title = "Number of Genes Per Cluster - CoXpress", 
           x = "Cluster", y = "Number of Genes") +
      theme_minimal() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            text = element_text(size = 14))
  })
  
  output$cluster_overlap_plot <- renderPlot({
    req(values$clustering_complete)
    values$clustering_plots$cluster_overlap
  })
  
  # Render the tables using native Shiny
  output$module_correlation_table <- renderTable({
    req(values$testing_complete)
    req(values$module_correlation_data)
    values$module_correlation_data
  }, digits = 4)  # Round numbers to 4 decimal places
  
  output$null_distribution_table <- renderTable({
    req(values$testing_complete)
    req(values$null_distribution_data)
    values$null_distribution_data
  }, digits = 4)  # Round numbers to 4 decimal places
  
  # Modified plot outputs
  output$condition1_plot <- renderUI({
    req(values$testing_complete, input$module, input$plot)
    
    # Get absolute path to the image
    filename <- normalizePath(file.path("output", "plotting", 
                            paste0(input$module, "_condition1_distribution.png")),
                            mustWork = FALSE)
    
    # Check if file exists
    if (!file.exists(filename)) {
      return(p("Plot not available"))
    }
    
    # Add timestamp to force image refresh
    timestamp <- format(Sys.time(), "%Y%m%d%H%M%OS6")
    
    tags$img(
      src = paste0("plotting/", basename(filename), "?t=", timestamp),
      style = "max-width: 100%; height: auto;",
      alt = paste("Distribution Plot for", input$module, "Condition 1")
    )
  })
  
  output$condition2_plot <- renderUI({
    req(values$testing_complete, input$module, input$plot)
    
    # Get absolute path to the image
    filename <- normalizePath(file.path("output", "plotting", 
                            paste0(input$module, "_condition2_distribution.png")),
                            mustWork = FALSE)
    
    # Check if file exists
    if (!file.exists(filename)) {
      return(p("Plot not available"))
    }
    
    # Add timestamp to force image refresh
    timestamp <- format(Sys.time(), "%Y%m%d%H%M%OS6")
    
    tags$img(
      src = paste0("plotting/", basename(filename), "?t=", timestamp),
      style = "max-width: 100%; height: auto;",
      alt = paste("Distribution Plot for", input$module, "Condition 2")
    )
  })
}

# Add this at the end of your app.R, before shinyApp()
addResourcePath("plotting", "output/plotting")

shinyApp(ui = ui, server = server)