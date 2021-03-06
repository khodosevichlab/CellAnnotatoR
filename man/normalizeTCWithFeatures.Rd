% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotation.R
\name{normalizeTCWithFeatures}
\alias{normalizeTCWithFeatures}
\title{Normalize Total Count with Features}
\usage{
normalizeTCWithFeatures(
  cm,
  max.quantile = 0.95,
  max.smooth = 1e-10,
  transpose = T
)
}
\arguments{
\item{cm}{gene count matrix with cells by columns and genes by rows. May be in raw, TC-normalized or tf-idf-normalized format (in case of tf-idf, `prenormalized` must be set to `T`)}

\item{max.quantile}{quantile to be used for max estimation during scaling}
}
\value{
Normalized matrix of the same shape as `cm`
}
\description{
Normalize `cm` matrix using total counts and then column-wise min-max scaling
}
