% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_known_genes.R
\name{check_known_genes}
\alias{check_known_genes}
\title{check that all symbols are in the GENCODE data}
\usage{
check_known_genes(known_genes, known_genes_file)
}
\arguments{
\item{known_genes}{vector of known gene symbols}

\item{known_genes_file}{name of the known_genes file}
}
\value{
A file of variant-gene pair predictions, with associated scores, saved in the given output directory.
}
\description{
check that all symbols are in the GENCODE data
}
