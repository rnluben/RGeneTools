#' Request info for a gene symbol
#'
#' Queries Ensembl REST API (build [38](https://rest.ensembl.org/) or
#' [37](https://grch37.rest.ensembl.org/)). Includes gene coordinates,
#' description and biotype (e.g. whether protein coding or not).
#'
#' @param gene_symbols A character vector of gene symbols e.g. "ESPN"
#' @param build Either 'GRCh38' (default) or 'GRCh37'.
#'
#' @return A tibble.
#' @export
#' @examples
#'  lookup_gene_symbols(c("ARMS2", "CFH"))
lookup_gene_symbols <- function(gene_symbols,
                                build = "GRCh38") {
  gene_symbols |>
    purrr::map(\(x) lookup_gene_symbol_single(gene_symbol = x,
                                              build = build)) |>
    dplyr::bind_rows()
}

#' Get VEP annotation for a variant
#'
#' Queries Ensembl REST API (build [38](https://rest.ensembl.org/) or
#' [37](https://grch37.rest.ensembl.org/)). Includes most severe consequence
#' across transcripts.
#'
#' @param variants Character vector of variants matching format specified by
#'   [Ensembl](https://rest.ensembl.org/documentation/info/vep_region_get), e.g.
#'   '9:22125503-22125502:1/C'.
#' @inheritParams lookup_gene_symbols
#'
#' @return A tibble
#' @export
#'
#' @examples
#' lookup_variant_vep(variants = "9:22125503-22125502:1/C")
lookup_variant_vep <- function(variants,
                               build = "GRCh38") {
  variants |>
    purrr::map(\(x) lookup_variant_vep_single(variant = x,
                                              build = build)) |>
    dplyr::bind_rows()
}

# PRIVATE -----------------------------------------------------------------

lookup_gene_symbol_single <- function(gene_symbol,
                                      build) {
  req <- req_build_ensembl(
    build = build,
    ext = paste0("lookup/symbol/homo_sapiens/", gene_symbol, "?")
  )
  
  tryCatch(
    expr = req |>
      httr2::req_perform() |>
      httr2::resp_body_json() |>
      tibble::as_tibble(),
    error = function(e)
      tibble::tibble(display_name = gene_symbol,
                     assembly_name = build)
  )
}

lookup_variant_vep_single <- function(build,
                                      variant) {
  req <- req_build_ensembl(build = build,
                           ext = paste0("/vep/human/region/", variant, "?"))
  
  
  tryCatch(
    expr = {
      result <- req |>
        httr2::req_perform() |>
        httr2::resp_body_json()
      
      c(purrr::discard(result[[1]], is.list),
        list(transcript_consequences = purrr::keep(result[[1]], is.list))) |>
        tibble::as_tibble_row()
    },
    error = function(e)
      tibble::tibble(input = variant,
                     assembly_name = build)
  )
}

req_build_ensembl <- function(build,
                              ext) {
  # validate args
  match.arg(build,
            c("GRCh38", "GRCh37"))
  
  server <- switch(build,
                   GRCh38 = "https://rest.ensembl.org/",
                   GRCh37 = "https://grch37.rest.ensembl.org/")
  
  # define REST query to get the gene ID from the gene name
  req <- httr2::request(paste0(server, ext))
  
  req |>
    httr2::req_headers("Accept" = "application/json")
}