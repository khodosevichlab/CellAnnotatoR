getAnnotationConfidence <- function(annotation, scores) {
  mapply(function(i,j) scores[i, j], 1:nrow(scores), match(annotation, colnames(scores))) %>%
    setNames(rownames(scores)) %>% unlist()
}

#' Score Per Cell Uncertainty
#' @description Score uncertainty per each cell
#' @param annotation assigned annotation labels. Can be obtrained with `assignCellsByScores`.
#' @param scores.norm assignment scores. Can be obtrained with `assignCellsByScores`.
#' @param score.info list of marker scores per cell type. Can be obtained with `getMarkerScoreInfo`
#' @param cur.types subset of types used to measure uncertainty. Default: all values presented in annotation
#' @param coverage.max.quantile all coverage values above this quantille are winsorized
#' @return list of scores with elements:\itemize{
#'   \item{positive: uncertainty based on conflicting positive markers}
#'   \item{negative: uncertainty based on large expression of negative markers}
#'   \item{coverage: uncertainty based on low coverage}
#' }
scorePerCellUncertainty <- function(annotation, scores.norm, score.info, cur.types=unique(annotation), coverage.max.quantile=0.75) {
  positive <- getAnnotationConfidence(annotation, scores.norm[names(annotation),])

  cbs.per.type <- split(names(annotation), annotation)
  negative.mult <- cur.types %>%
    lapply(function(n) score.info[[n]]$mult[cbs.per.type[[n]]]) %>%
    Reduce(c, .)

  coverage <- cur.types %>%
    lapply(function(n) score.info[[n]]$scores[cbs.per.type[[n]]]) %>%
    Reduce(c, .) %>% pmin(quantile(., coverage.max.quantile))

  return(list(
    positive=(1 - positive),
    negative=(1 - negative.mult),
    coverage=(1 - coverage / max(coverage))
  ))
}

#' Score Cell Uncertainty Per Level
#' @description Estimate uncertainty scores per annotation level
#' @param ann.info.per.level annotation info per level. Result of `assignCellsByScores`.
#' @inheritParams scorePerCellUncertainty
#' @param verbose show progress bar for estimation
#' @inheritDotParams scorePerCellUncertainty cur.types coverage.max.quantile
#' @return Uncertainty per level for each cell
#' @examples
#'   score_info <- getMarkerScoreInfo(clf_data)
#'   unc_info <- scoreCellUncertaintyPerLevel(ann_by_level, score_info)
#'
#' @export
scoreCellUncertaintyPerLevel <- function(ann.info.per.level, score.info, verbose=F, ...) {
  names(ann.info.per.level[[1]]) %>% setNames(., .) %>% plapply(function(n)
    scorePerCellUncertainty(ann.info.per.level$annotation[[n]], ann.info.per.level$scores[[n]], score.info, ...),
    verbose=verbose)
}

scorePerClusterUncertainty <- function(unc.per.cell, clusters) {
  return(unc.per.cell %$% list(
    coverage = coverage %>% split(clusters[names(.)]) %>% sapply(function(x) mean(x > 0.99)),
    positive = positive %>% split(clusters[names(.)]) %>% sapply(mean),
    negative = negative %>% split(clusters[names(.)]) %>% sapply(mean)
  ))
}

#' Score Cluster Uncertainty Per Level
#' @description Score uncertainty for each cluster per level
#' @param cell.uncertainty.per.level uncertainty per cell per level. Can be obtained from `scoreCellUncertaintyPerLevel`
#' @param clusters vector with cluster label for each cell
#' @return Uncertainty per level for each cluster
#' @examples
#'   score_info <- getMarkerScoreInfo(clf_data)
#'   unc_info <- scoreCellUncertaintyPerLevel(ann_by_level, score_info)
#'   unc_per_clust <- scoreClusterUncertaintyPerLevel(unc_info, clusters)
#'
#' @export
scoreClusterUncertaintyPerLevel <- function(cell.uncertainty.per.level, clusters) {
  lapply(cell.uncertainty.per.level, scorePerClusterUncertainty, clusters)
}

filterAnnotationByUncertainty <- function(annotation, scores, score.info, cur.types=unique(annotation), clusters=NULL,
                                          thresholds=c(coverage=0.5, negative=0.5, positive=0.75)) {
  unc.per.cell <- scorePerCellUncertainty(annotation, scores, score.info, cur.types=cur.types)

  if (!is.null(clusters)) {
    unc.per.clust <- scorePerClusterUncertainty(unc.per.cell, clusters)
    for (n in names(thresholds)) {
      annotation[unc.per.clust[[n]][as.character(clusters[names(annotation)])] > thresholds[[n]]] <- NA
    }
  } else {
    for (n in names(thresholds)) {
      annotation[unc.per.cell[[n]][names(annotation)] > thresholds[[n]]] <- NA
    }
  }

  return(annotation)
}
