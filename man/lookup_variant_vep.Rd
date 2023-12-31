% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ensemble_lookup_gene_symbol.R
\name{lookup_variant_vep}
\alias{lookup_variant_vep}
\title{Get VEP annotation for a variant}
\usage{
lookup_variant_vep(variants, build = "GRCh38")
}
\arguments{
\item{variants}{Character vector of variants matching format specified by
\href{https://rest.ensembl.org/documentation/info/vep_region_get}{Ensembl}, e.g.
'9:22125503-22125502:1/C'.}

\item{build}{Either 'GRCh38' (default) or 'GRCh37'.}
}
\value{
A tibble
}
\description{
Queries Ensembl REST API (build \href{https://rest.ensembl.org/}{38} or
\href{https://grch37.rest.ensembl.org/}{37}). Includes most severe consequence
across transcripts.
}
\examples{
lookup_variant_vep(variants = "9:22125503-22125502:1/C")
}
