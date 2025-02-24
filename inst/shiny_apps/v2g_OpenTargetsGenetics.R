library(shiny)
library(reactable)
library(dplyr)
library(purrr)
library(otargen)
library(crosstalk)
library(openxlsx)

# Define UI for application
ui <- fluidPage(
  titlePanel("Variant to Gene (V2G) Mapping"),
  sidebarLayout(
    sidebarPanel(
      textAreaInput(
        "rsids",
        "Enter variant rsids (one per line):",
        value = "",
        rows = 10,
        resize = "vertical",
        placeholder = "rs12913832"
      ),
      actionButton("retrieve", "Retrieve V2G"),
      downloadButton("downloadData", "Download as Excel"),
      uiOutput("filterUI"),
      radioButtons("outputType", "Select Output Type:",
                   choices = list("Table" = "table", "Variant per Row" = "variant_row"),
                   selected = "table")
    ),
    mainPanel(reactableOutput("v2gTable"), textOutput("errorMessage"))
  )
)

# Define server logic
server <- function(input, output) {
  possibly_genesForVariant <- purrr::possibly(otargen::genesForVariant, otherwise = list(tssd = list(gene.symbol = "Unknown")))

  v2g_data <- eventReactive(input$retrieve, {
    rsids <- strsplit(input$rsids, "\n")[[1]]
    rsids <- trimws(rsids)
    rsids <- rsids[rsids != ""]

    data <- rsids |>
      purrr::set_names() |>
      purrr::map(possibly_genesForVariant) |>
      purrr::map(\(x) x$v2g) |>
      dplyr::bind_rows(.id = "variant_input")

    if (nrow(data) == 0) {
      return(NULL)
    } else {
      return(data)
    }
  })

  shared_data <- reactive({
    data <- v2g_data()
    if (!is.null(data)) {
      SharedData$new(data)
    } else {
      NULL
    }
  })

  output$v2gTable <- renderReactable({
    data <- v2g_data()
    if (is.null(data)) {
      output$errorMessage <- renderText("No results found for the provided rsids.")
      return(NULL)
    } else {
      output$errorMessage <- renderText("")
      if (input$outputType == "table") {
        reactable(
          shared_data(),
          columns = list(
            gene.symbol = colDef(name = "Gene Symbol"),
            variant = colDef(name = "Variant"),
            overallScore = colDef(name = "Overall Score"),
            gene.id = colDef(name = "Gene ID")
          ),
          searchable = TRUE,
          resizable = TRUE,
          filterable = TRUE,
          showPageSizeOptions = TRUE,
          paginationType = "jump",
          defaultPageSize = 10,
          pageSizeOptions = c(10, 25, 50, 100, 200)
        )
      } else {
        variant_data <- data %>%
          group_by(variant_input) %>%
          summarize(genes = paste(gene.symbol[order(-overallScore)], collapse = ", "))
        reactable(
          variant_data,
          columns = list(
            variant_input = colDef(name = "Variant"),
            genes = colDef(name = "Genes (sorted by V2G score)")
          ),
          searchable = TRUE,
          resizable = TRUE,
          filterable = TRUE,
          showPageSizeOptions = TRUE,
          paginationType = "jump",
          defaultPageSize = 10,
          pageSizeOptions = c(10, 25, 50, 100, 200)
        )
      }
    }
  })

  output$filterUI <- renderUI({
    data <- shared_data()
    if (!is.null(data)) {
      tagList(
        filter_select(
          "filterRsidInput",
          "Filter by rsid (input)",
          data,
          ~ variant_input
        ),
        filter_select("filterRsid", "Filter by rsid", data, ~ variant),
        filter_select("filterGene", "Filter by gene", data, ~ gene.symbol)
      )
    }
  })

  output$downloadData <- downloadHandler(
    filename = function() {
      paste("v2g_data-", Sys.Date(), ".xlsx", sep = "")
    },
    content = function(file) {
      data <- v2g_data()
      if (!is.null(data)) {
        variant_data <- data %>%
          group_by(variant_input) %>%
          summarize(genes = paste(gene.symbol[order(-overallScore)], collapse = ", "))
        
        wb <- createWorkbook()
        addWorksheet(wb, "Full table")
        writeData(wb, "Full table", data)
        addWorksheet(wb, "Variant per row")
        writeData(wb, "Variant per row", variant_data)
        saveWorkbook(wb, file, overwrite = TRUE)
      }
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
