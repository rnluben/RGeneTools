#' Request features for a region
#'
#' Retrieves features (e.g. genes, transcripts, variants and more) that overlap
#' a given region.
#'
#' @param regions A character vector of genomic regions e.g.
#'   ""1:55516888-55516888"
#' @param build Either 'GRCh38' (default) or 'GRCh37'.
#'
#' @return A tibble.
#' @export
#' @examples
#'  ensembl_get_overlap_region(regions = "1:55516888-55516888", build = "GRCh37")
ensembl_get_overlap_region <- function(regions, build = "GRCh38") {
  regions |>
    purrr::map(\(x) ensembl_get_overlap_region_single(region = x, build = build)) |>
    dplyr::bind_rows()
}

# PRIVATE -----------------------------------------------------------------

ensembl_get_overlap_region_single <- function(region,
                                      build) {
  req <- req_build_ensembl(
    build = build,
    ext = paste0("overlap/region/homo_sapiens/", region, "?content-type=application/json;feature=variation")
  )
  
  tryCatch(
    expr = req |>
      httr2::req_perform() |>
      httr2::resp_body_json() |>
      purrr::pluck(1) |>
      purrr::map_if(.p = is.list, \(x) paste(x, sep = "", collapse = "_")) |>
      
      # alternatively...
      # purrr::map_if(.p = is.list, \(x) list(x)) |>
      tibble::as_tibble(),
    error = function(e)
      tibble::tibble(display_name = region,
                     assembly_name = build)
  )
}

