#' Nearest gene for given variant
#'
#' @param rsids Character vector of rsids
#'
#' @return A data frame
#' @export
#'
#' @examples OpenTargetFunc("rs12345")
OpenTargetFunc <- function(rsids) {
  rsids |>
    purrr::map(OpenTargetFunc_single) |>
    dplyr::bind_rows()
}

OpenTargetFunc_single <- function(query_rsID) {
   # Build query string
   query_string = "
   query useSearchToConvertRSIDIntoIDFormat($query_rsID: String!) {
     search(queryString: $query_rsID) {
       variants {
         id
         rsId
         nearestGene {
           id
           start
           symbol
           tss
           description
           chromosome
           exons
         }
         nearestGeneDistance
       }
     }
   } 
   "
   
   # Set base URL of GraphQL API endpoint
   base_url <- "https://api.genetics.opentargets.org/graphql"
   
   # Set variables object of arguments to be passed to endpoint
   variables <- list("query_rsID" = query_rsID)
   
   # Construct POST request body object with query string and variables
   post_body <- list(query = query_string, variables = variables)
   
   # Perform POST request
   r <- httr::POST(url=base_url, body=post_body, encode='json')
   
   if(length(httr::content(r)$data$search$variants)==1) {
      NearestGene <- httr::content(r)$data$search$variants[[1]]$nearestGene$symbol
      Variant <- httr::content(r)$data$search$variants[[1]]$id
      rsid <- httr::content(r)$data$search$variants[[1]]$rsId
      GeneID <- httr::content(r)$data$search$variants[[1]]$nearestGene$id
      Position <- httr::content(r)$data$search$variants[[1]]$nearestGene$start
      Desc <- httr::content(r)$data$search$variants[[1]]$nearestGene$description
      Chr <- httr::content(r)$data$search$variants[[1]]$nearestGene$chromosome
      Result <- tibble::tibble(Chr=Chr,Position=Position,Variant=Variant,rsid=rsid,NearestGene=NearestGene,GeneID=GeneID,Desc=Desc)
   } else {
      Result <- tibble::tibble(Chr="",Position=NA_integer_,Variant="",rsid=query_rsID,NearestGene="",GeneID="",Desc="")
   }
   Result
}
