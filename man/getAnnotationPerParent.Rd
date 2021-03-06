% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{getAnnotationPerParent}
\alias{getAnnotationPerParent}
\title{Get Annotation Per Parent}
\usage{
getAnnotationPerParent(classification.tree, annotation)
}
\arguments{
\item{classification.tree}{cell type hierarchy represented by graph. Part of `clf_data` object from `getClassificationData`}

\item{annotation}{annotation for high-resolution. Cell type names must correspond to nodes in the `classification.tree`}
}
\value{
list of sub-annotations named by parent types
}
\description{
for each cell type get annotation of its subtypes on the next hierarchy level
}
\examples{
  ann_by_parent <- getAnnotationPerParent(clf_data$classification.tree, annotation)
}
