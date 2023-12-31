% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ensemble_lookup_gene_symbol.R
\name{lookup_gene_symbols}
\alias{lookup_gene_symbols}
\title{Request info for a gene symbol}
\usage{
lookup_gene_symbols(gene_symbols, build = "GRCh38")
}
\arguments{
\item{gene_symbols}{A character vector of gene symbols e.g. "ESPN"}

\item{build}{Either 'GRCh38' (default) or 'GRCh37'.}
}
\value{
A tibble.
}
\description{
Queries Ensembl REST API (build \href{https://rest.ensembl.org/}{38} or
\href{https://grch37.rest.ensembl.org/}{37}). Includes gene coordinates,
description and biotype (e.g. whether protein coding or not).
}
\examples{
 lookup_gene_symbols(c("ARMS2", "CFH"))
}
