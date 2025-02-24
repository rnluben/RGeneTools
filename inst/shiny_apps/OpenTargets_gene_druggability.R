library(shiny)
library(RGeneTools)
library(dplyr)
library(purrr)
library(tibble)
library(readr)
library(flextable)
library(openxlsx)

default_druggable_categories <- c("approved_drug", "advanced_clinical", "phase_1_clinical", "high_quality_ligand", "high_quality_pocket", "druggable_family", "uni_prot_loc_high_conf", "go_cc_high_conf")

druggability_labels <- c(
  "Approved Drug" = "approved_drug",
  "Advanced Clinical" = "advanced_clinical",
  "Phase 1 Clinical" = "phase_1_clinical",
  "Structure with Ligand" = "structure_with_ligand",
  "High-Quality Ligand" = "high_quality_ligand",
  "High-Quality Pocket" = "high_quality_pocket",
  "Med-Quality Pocket" = "med_quality_pocket",
  "Druggable Family" = "druggable_family",
  "UniProt loc high conf" = "uni_prot_loc_high_conf",
  "GO CC high conf" = "go_cc_high_conf",
  "UniProt loc med conf" = "uni_prot_loc_med_conf",
  "UniProt SigP or TMHMM" = "uni_prot_sig_p_or_tmhmm",
  "GO CC med conf" = "go_cc_med_conf",
  "Human Protein Atlas loc" = "human_protein_atlas_loc",
  "Literature" = "literature",
  "UniProt Ubiquitination" = "uni_prot_ubiquitination",
  "Database Ubiquitination" = "database_ubiquitination",
  "Half-life Data" = "half_life_data",
  "Small Molecule Binder" = "small_molecule_binder"
)

# Define UI for application
ui <- fluidPage(
  titlePanel("Gene Druggability"),
  sidebarLayout(
    sidebarPanel(
      textAreaInput("genes_input", "Enter gene names (one per line):", "CFH\nABCA1", rows = 10, resize = "both", placeholder = "CFH\nABCA1"),
      selectizeInput("druggability_criteria", "Select Druggability Criteria:", 
                     choices = druggability_labels,
                     selected = default_druggable_categories,
                     multiple = TRUE),
      actionButton("analyze", "Analyze"),
      actionButton("reset", "Reset"),
      downloadButton("download_excel", "Download as Excel")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Druggable Genes", uiOutput("table1")),
        tabPanel("Druggable Genes Compounds", uiOutput("table2"))
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  observeEvent(input$reset, {
    updateSelectizeInput(session, "druggability_criteria", selected = default_druggable_categories)
  })
  
  analyze_genes <- eventReactive(input$analyze, {
    req(input$genes_input)
    
    showNotification("Retrieving Ensembl IDs...", type = "message", duration = NULL, id = "retrieve_ensemblids")
    genes <- strsplit(input$genes_input, "\n")[[1]] |> unique() |> RGeneTools::lookup_gene_symbols()
    
    if (nrow(genes) == 0) {
      removeNotification(id = "retrieve_ensemblids")
      return(NULL)
    }
    
    showNotification("Querying OpenTargets...", type = "message", duration = NULL, id = "query_opentargets")
    ot_data_all <- list(
      known_drugs = RGeneTools::ot_knowndrugs(ensgIds = na.omit(genes$id)),
      tractability = RGeneTools::ot_target_tractability(ensgIds = na.omit(genes$id))
    )
    removeNotification(id = "retrieve_ensemblids")
    removeNotification(id = "query_opentargets")
    
    gene_symbols_ids <- genes %>%
      dplyr::select(approvedSymbol = display_name, ensgId = id)
    
    selected_druggable_categories <- input$druggability_criteria
    
    tractability <- ot_data_all$tractability %>%
      dplyr::mutate(druggable_target_score = rowSums(dplyr::across(tidyselect::all_of(selected_druggable_categories)), na.rm = TRUE))
    
    tractable_targets <- tractability %>% 
      dplyr::filter(druggable_target_score > 0) %>% 
      dplyr::pull(approvedSymbol) %>% 
      unique()
    
    known_drugs <- ot_data_all$known_drugs %>%
      dplyr::left_join(gene_symbols_ids, by = "ensgId") %>%
      dplyr::relocate(approvedSymbol, .after = ensgId)
    
    if ("drug.name" %in% names(known_drugs)) {
      existing_drug <- known_drugs %>%
        dplyr::filter(!is.na(drug.name)) %>%
        dplyr::pull(approvedSymbol) %>%
        unique()
    } else {
      existing_drug <- character()
    }
    
    druggable_or_drugged <- list(
      existing_drug = sort(existing_drug),
      tractable_targets = sort(tractable_targets),
      union = sort(union(existing_drug, tractable_targets))
    )
    
    if ("drug.name" %in% names(known_drugs)) {
    tab_druggable_genes <- ot_data_all$known_drugs %>%
      dplyr::full_join(gene_symbols_ids, by = "ensgId") %>%
      dplyr::relocate(approvedSymbol, .after = ensgId) |>
      # dplyr::filter(!is.na(drug.name)) %>%
      dplyr::group_by(approvedSymbol, drug.name) %>%
      dplyr::mutate(phase = dplyr::case_when(is.na(phase) ~ -1L,
                                            TRUE ~ phase)) %>% 
      dplyr::mutate(max_phase = max(phase, na.rm = TRUE)) %>%
      dplyr::slice(which.max(phase)) %>%
      dplyr::ungroup() %>%
      dplyr::select(approvedSymbol, drug.name, max_phase, status, drugType, mechanismOfAction) %>%
      dplyr::arrange(approvedSymbol, drug.name, max_phase, status, drugType, mechanismOfAction) %>%
      dplyr::bind_rows(tibble::tibble(approvedSymbol = subset(druggable_or_drugged$tractable_targets, !druggable_or_drugged$tractable_targets %in% druggable_or_drugged$existing_drug))) |>
      dplyr::mutate(dplyr::across(dplyr::everything(), \(x) ifelse(is.na(x), "-", x))) |>
      dplyr::arrange(approvedSymbol) %>% 
      dplyr::mutate(max_phase = dplyr::case_when(max_phase == -1L ~ "-",
                                             TRUE ~ as.character(max_phase))) |>
      dplyr::distinct()
    
    tab_druggable_genes_compounds <- ot_data_all$known_drugs %>% 
      dplyr::left_join(gene_symbols_ids, by = "ensgId") |>
      dplyr::filter(approvedSymbol %in% gene_symbols_ids$approvedSymbol) %>% 
      dplyr::filter(!is.na(disease.name)) %>% 
      dplyr::group_by(drug.name) %>% 
      dplyr::summarise(target_disease = paste(sort(unique(disease.name)), sep = "", collapse = "; ")) %>% 
      dplyr::rename(`Compound name` = drug.name, `Target disease` = target_disease)
    
    result <- list(tab_druggable_genes = tab_druggable_genes,
                   tab_druggable_genes_compounds = tab_druggable_genes_compounds)
    } else {
      result <- NULL
      showNotification("No results found. Please check gene names and try again.", type = "error", duration = 5)
    }
    
    result
  })
  
  output$tables_ui <- renderUI({
    results <- analyze_genes()
    if (is.null(results)) {
      NULL
    }
    
    tagList(
      h4("Druggable Genes"),
      uiOutput("table1"),
      h4("Druggable Genes Compounds"),
      uiOutput("table2")
    )
  })
  
  output$table1 <- renderUI({
    results <- analyze_genes()
    req(results$tab_druggable_genes)
    results$tab_druggable_genes |> 
      flextable::as_grouped_data(groups = "approvedSymbol") |>
      flextable::flextable() |>
      bold(j = 1) |>
      flextable::autofit() |> 
      htmltools_value()
  })
  
  output$table2 <- renderUI({
    results <- analyze_genes()
    req(results$tab_druggable_genes_compounds)
    flextable::flextable(results$tab_druggable_genes_compounds) |> 
      flextable::theme_zebra() |>
      flextable::autofit() |> 
      htmltools_value()
  })
  
  output$download_excel <- downloadHandler(
    filename = function() {
      paste("gene_druggability_analysis", Sys.Date(), ".xlsx", sep = "")
    },
    content = function(file) {
      results <- analyze_genes()
      wb <- createWorkbook()
      addWorksheet(wb, "Druggable Genes")
      writeData(wb, "Druggable Genes", results$tab_druggable_genes)
      addWorksheet(wb, "Druggable Genes Compounds")
      writeData(wb, "Druggable Genes Compounds", results$tab_druggable_genes_compounds)
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)
