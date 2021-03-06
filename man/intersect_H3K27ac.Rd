% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/intersect_H3K27ac.R
\name{intersect_H3K27ac}
\alias{intersect_H3K27ac}
\title{Intersect with the DHS annotations}
\usage{
intersect_H3K27ac(l, query, DHSs, H3K27ac, ...)
}
\arguments{
\item{query}{A query bed}

\item{H3K27ac}{A BED tibble (the H3K27ac-in-DHSs annotations) to be intersected by the query bed}

\item{...}{query columns to be retained in the output}
}
\value{
A tibble of intersected columns, with .query and .annotation suffixes
}
\description{
Intersect a query BED with the DHS annotations.
}
