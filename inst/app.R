library(shiny)
library(RGeneTools)
library(stringr)

ui <- fluidPage(textAreaInput("rsids", "RSIDs (one on each new line)", value = "rs12345\nrs1061170\nrs859705", height = '200px'),
                tableOutput("nearest_genes"))

server <- function(input, output, session) {
  output$nearest_genes <-
    renderTable(OpenTargetFunc(rsids = str_split(input$rsids, "\\n")[[1]]))
}

shinyApp(ui, server)