getAnnotationConfidence <- function(annotation, scores) {
  mapply(function(i,j) scores[i, j], 1:nrow(scores), match(annotation, colnames(scores))) %>%
    setNames(rownames(scores))
}

#' @export
scoreUncertaintyPerCell <- function(annotation, scores.norm, score.info, cur.types=unique(annotation), coverage.max.quantile=0.75) {
  positive <- getAnnotationConfidence(annotation, scores.norm)

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

#' @export
scoreUncertaintyPerCluster <- function(unc.per.cell, clusters) {
  return(unc.per.cell %$% list(
    coverage = coverage %>% split(clusters[names(.)]) %>% sapply(function(x) mean(x > 0.99)),
    positive = positive %>% split(clusters[names(.)]) %>% sapply(mean),
    negative = negative %>% split(clusters[names(.)]) %>% sapply(mean)
  ))
}

filterAnnotationByUncertainty <- function(annotation, scores, score.info, cur.types=unique(annotation), clusters=NULL,
                                          thresholds=c(coverage=0.5, negative=0.5, positive=0.75)) {
  unc.per.cell <- scoreUncertaintyPerCell(annotation, scores, score.info, cur.types=cur.types)

  if (!is.null(clusters)) {
    unc.per.clust <- scoreUncertaintyPerCluster(unc.per.cell, clusters)
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
