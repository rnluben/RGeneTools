#' Get drugs known to target a given EnsemblId product
#'
#' @param ensgIds Character vector of EnsemblIDs
#' @param size Maximum number of rows to return per EnsemblID queried
#'
#' @return A tibble
#' @export
#' @examples
#' 
#' ot_knowndrugs("ENSG00000157764")
#' 
ot_knowndrugs <- function(ensgIds,
                          size = 1000) {
  ensgIds |>
    purrr::set_names() |>
    purrr::map(\(x) ot_knowndrugs_single(ensgId = x, size = size)) |>
    dplyr::bind_rows(.id = "ensgId")
}

ot_knowndrugs_single <- function(ensgId,
                          size = 1000) {
  result <- submit_ot_query(
    variables = list("ensgId" = ensgId, "size" = size),
    query = "
      query KnownDrugsQuery(
        $ensgId: String!
        $cursor: String
        $freeTextQuery: String
        $size: Int = 10
      ) {
        target(ensemblId: $ensgId) {
          id
          knownDrugs(cursor: $cursor, freeTextQuery: $freeTextQuery, size: $size) {
            count
            cursor
            rows {
              phase
              status
              urls {
                name
                url
              }
              disease {
                id
                name
              }
              drug {
                id
                name
                mechanismsOfAction {
                  rows {
                    actionType
                    targets {
                      id
                    }
                  }
                }
              }
              drugType
              mechanismOfAction
            }
          }
        }
      }
    
      "
  )
  
  result$data$target$knownDrugs$rows |>
    purrr::map(\(x) {
      if (is.null(x$status)) {
        x$status <- "Unknown"
      }
      
      x
    }) |>
    purrr::map(as.data.frame) |>
    dplyr::bind_rows() |>
    tidyr::unite(col = "urls.name", 
                 tidyselect::starts_with("urls.name"), 
                 sep = "|", 
                 remove = TRUE, 
                 na.rm = TRUE) |>
    tidyr::unite(col = "urls.url", 
                 tidyselect::starts_with("urls.url"), 
                 sep = "|", 
                 remove = TRUE, 
                 na.rm = TRUE) |>
    tidyr::unite(col = "drug.mechanismsOfAction", 
                 tidyselect::starts_with("drug.mechanismsOfAction"), 
                 sep = "|", 
                 remove = TRUE, 
                 na.rm = TRUE) |>
    tibble::as_tibble()
}


ot_associatedDiseases <- function(ensemblId) {
  submit_ot_query(
    variables = list("ensemblId" = ensemblId),
    query = "
  query target($ensemblId: String!){
    target(ensemblId: $ensemblId){
      id
      approvedSymbol
      biotype
      associatedDiseases {
        count
        rows {
          disease {
            id
            name
          }
          datasourceScores {
            id
            score
          }
        }
      }
      geneticConstraint {
        constraintType
        exp
        obs
        score
        oe
        oeLower
        oeUpper
      }
      tractability {
        label
        modality
        value
      }
    }
  }
"
  )
}

#' Base function for submitting queries to OpenaTargets GraphQL APi
#'
#' See  [OpenTargets GraphQL API
#' documentation](https://platform-docs.opentargets.org/data-access/graphql-api),
#' the [API playground](https://platform.opentargets.org/api) and the online
#' browser, which includes option to download results or retrieve the API query
#' that will retrieve these
#' e.g.[ENSG00000169083](https://platform.opentargets.org/target/ENSG00000169083)
#'
#' @param query A query string
#' @param variables A named list of variables, corresponding to `query`
#'
#' @return A response in JSON format.
#' @noRd
submit_ot_query <- function(variables,
                            query) {
  # Set base URL of GraphQL API endpoint
  base_url <- "https://api.platform.opentargets.org/api/v4/graphql"
  
  # Construct POST request body object with query string and variables
  post_body <- list(query = query, variables = variables)
  
  # Perform POST request
  httr::POST(url = base_url,
             body = post_body,
             encode = 'json') |>
    httr::content()
}
