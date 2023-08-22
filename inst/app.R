library(shiny)
library(RGeneTools)
library(stringr)

ui <- fluidPage(textAreaInput("rsids", "RSIDs (one on each new line)"),
                tableOutput("nearest_genes"))

server <- function(input, output, session) {
  output$nearest_genes <-
    renderTable(OpenTargetFunc(rsids = str_split(input$rsids, "\\n")[[1]]))
}

shinyApp(ui, server)