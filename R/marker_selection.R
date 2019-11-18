aggregateScoreChangePerGene <- function(d.scores, annotation, marker.type, cell.type, target.type=NULL, balance.cell.types=T) {
  aggr.func <- if (balance.cell.types) Matrix::colMeans else Matrix::colSums

  scores.per.type <- split(names(annotation), annotation) %>%
    sapply(function(cbs) aggr.func(d.scores[cbs,])) %>% `*`(-1)
  scores.per.type[,cell.type] %<>% `*`(-1)
  d.scores <- rowSums(scores.per.type) %>% sort(decreasing=T)
  res <- tibble::tibble(Gene=names(d.scores), Score=d.scores, Type=cell.type, MT=marker.type)
  if (!is.null(target.type)) {
    res$ScoreTarget <- scores.per.type[res$Gene, target.type]
  }

  return(res)
}

estimatePositiveMarkerScoreChange <- function(cell.type, annotation, cm.norm, de.genes, pos.scores, neg.scores, sum.scores, ...) {
  if (length(de.genes) == 0)
    stop("de.genes are empty for ", cell.type)

  sum.scores %<>% pmax(1e-30)

  c.exprs <- as.matrix(cm.norm[, de.genes])

  c.neg.scores <- (1 - neg.scores[, cell.type])
  c.pos.scores <- pos.scores[, cell.type] * c.neg.scores
  ds <- c.exprs * c.neg.scores
  ds <- (ds + c.pos.scores) / (ds + sum.scores) - (c.pos.scores / sum.scores)

  return(aggregateScoreChangePerGene(ds, annotation, marker.type="expressed", cell.type=cell.type, ...))
}

estimateNegativeMarkerScoreChange <- function(cell.type, annotation, cm.norm, de.genes, pos.scores, neg.scores,
                                              sum.scores, max.pos.scores, max.neg.score.per.pos.type=0.1, ...) {
  if (length(de.genes) == 0)
    stop("de.genes are empty for ", cell.type)

  sum.scores %<>% pmax(1e-30)
  max.pos.scores <- pmax(max.pos.scores[,cell.type], 1e-30)
  mask <- (annotation == cell.type)
  pos.ids <- which(mask)
  neg.ids <- which(!mask)

  c.exprs <- as.matrix(cm.norm[, de.genes])

  c.pos.scores <- pos.scores[,cell.type]
  c.pos.scores[pos.ids] <- sum.scores[pos.ids] # In ideal world, all positive score for a cell came from the proper cell type

  new.neg.scores <- estimateNewNegativeScores(c.exprs, max.pos.scores, neg.scores[,cell.type]) %>%
    `dimnames<-`(dimnames(c.exprs))

  old.scores <- c.pos.scores * (1 - neg.scores[,cell.type])
  ds <- (c.pos.scores * (1 - new.neg.scores) - old.scores)

  new.scores <- c.pos.scores * (1 - new.neg.scores)
  ds <- new.scores / pmax(sum.scores + new.scores - old.scores, 1e-30) - (old.scores / pmax(sum.scores, 1e-30))

  res <- aggregateScoreChangePerGene(ds, annotation, marker.type="not_expressed", cell.type=cell.type, ...)
  res$Score[Matrix::colMeans(new.neg.scores[pos.ids, res$Gene]) > max.neg.score.per.pos.type] <- -Inf
  return(res)
}

prepareMarkerScoringInfo <- function(score.infos) {
  return(list(
    sum.scores=sapply(score.infos, `[[`, "scores") %>% rowSums(),
    pos.scores=sapply(score.infos, `[[`, "scores.raw"),
    neg.scores=1 - sapply(score.infos, `[[`, "mult"),
    max.pos.scores=sapply(score.infos, `[[`, "max.positive")
  ))
}

initializeMarkerList <- function(pos.markers.per.type, cm.norm, annotation, debug=F) {
  marker.list <- lapply(pos.markers.per.type, function(x) list(expressed=c(), not_expressed=c()))
  subtypes.left <- names(pos.markers.per.type)
  for (i in 1:length(subtypes.left)) {
    s.info <- lapply(marker.list, getCellTypeScoreInfo, cm.norm) %>% prepareMarkerScoringInfo()

    marker.score <- subtypes.left %>% lapply(function(ct)
      estimatePositiveMarkerScoreChange(ct, annotation, cm.norm, pos.markers.per.type[[ct]], s.info$pos.scores,
                                        s.info$neg.scores, s.info$sum.scores, balance.cell.types=F)
    ) %>% Reduce(rbind, .) %>% .[which.max(.$Score),]

    if (marker.score$Score < 0) {
      marker.score <- subtypes.left %>% lapply(function(ct)
        estimatePositiveMarkerScoreChange(ct, annotation, cm.norm, pos.markers.per.type[[ct]], s.info$pos.scores,
                                          s.info$neg.scores, s.info$sum.scores)
      ) %>% Reduce(rbind, .) %>% .[which.max(.$Score),]

      if (marker.score$Score < 0) {
        warning("Negative score for ", marker.score$Type, ": ", marker.score$Score, ". Consider changing extending marker list.")
      }
    }

    marker.list[[c(marker.score$Type, marker.score$MT)]] %<>% c(marker.score$Gene)
    subtypes.left %<>% setdiff(marker.score$Type)

    if (debug) print(marker.score)
  }

  return(marker.list)
}

preSelectMarkerCandidates <- function(de.info) {
  de.info %<>% lapply(function(x)
    x[order(abs(x$Z), decreasing=T)[1:max(sum(abs(x$Z) > 2), min(nrow(x), 5))], ])

  res <- list(
    positive=lapply(de.info, function(df) df$Gene[df$Z > 0]),
    negative=lapply(de.info, function(df) df$Gene[df$Z < 0])
  )

  if (any(unlist(lapply(res, sapply, length)) == 0)) {
    warning("Some cell types don't have positive or negative markers")
  }

  return(res)
}

updateMarkersPerType <- function(markers.per.type, marker.list=NULL, marker.info=NULL) {
  if (!is.null(marker.list)) {
    for (n in names(marker.list)) {
      markers.per.type$positive[[n]] %<>% setdiff(marker.list[[n]]$expressed)
      for (n2 in (names(markers.per.type$negative) %>% .[. != n])) {
        markers.per.type$negative[[n2]] %<>% union(marker.list[[n]]$expressed)
      }
      markers.per.type$negative[[n]] %<>% setdiff(marker.list[[n]]$not_expressed)
    }
  } else if (!is.null(marker.info)) {
    if (marker.info$MT == "expressed") {
      markers.per.type$positive[[marker.info$Type]] %<>% setdiff(marker.info$Gene)
      for (n in (names(markers.per.type$negative) %>% .[. != marker.info$Type])) {
        markers.per.type$negative[[n]] %<>% union(marker.info$Gene)
      }
    } else {
      markers.per.type$negative[[marker.info$Type]] %<>% setdiff(marker.info$Gene)
    }
  } else {
    stop("Either marker.list or marker.info must be provided")
  }

  return(markers.per.type)
}

selectMarkersPerType <- function(cm.norm, marker.list, annotation, markers.per.type, max.iters=ncol(cm.norm), max.markers.per.type=10, max.uncertainty=0.2, log.step=0, debug=F) {
  score.info <- lapply(marker.list, getCellTypeScoreInfo, cm.norm)
  s.info <- prepareMarkerScoringInfo(score.info)
  mean.confidence <- 0
  mean.neg.score <- Inf

  for (i in 1:max.iters) {
    n.pos.markers <- sapply(marker.list, function(x) length(x$expressed))
    n.neg.markers <- sapply(marker.list, function(x) length(x$not_expressed))

    pos.marker.scores <- which(n.pos.markers <= max.markers.per.type) %>% names() %>% lapply(function(ct)
      estimatePositiveMarkerScoreChange(ct, annotation, cm.norm, markers.per.type$positive[[ct]], s.info$pos.scores, s.info$neg.scores, s.info$sum.scores)
    ) %>% Reduce(rbind, .)

    marker.info <- which(n.neg.markers <= max.markers.per.type) %>% names() %>% lapply(function(ct)
      estimateNegativeMarkerScoreChange(ct, annotation, cm.norm, markers.per.type$negative[[ct]], s.info$pos.scores, s.info$neg.scores, s.info$sum.scores, s.info$max.pos.scores)
    ) %>% Reduce(rbind, .) %>% rbind(pos.marker.scores)

    marker.info %<>% .[which.max(.$Score),]
    if (marker.info$Score < 0)
      break

    marker.list.new <- marker.list
    marker.list.new[[c(marker.info$Type, marker.info$MT)]] %<>% c(marker.info$Gene)

    score.info.new <- lapply(marker.list.new, getCellTypeScoreInfo, cm.norm)
    scores <- lapply(score.info.new, `[[`, "scores") %>% as.data.frame(optional=T) %>% normalizeScores()
    confidence <- getAnnotationConfidence(annotation, scores)
    mean.confidence.new <- mean(confidence)

    if (mean.confidence < mean.confidence.new) {
      mean.confidence <- mean.confidence.new
      score.info <- score.info.new
      marker.list <- marker.list.new
    }

    if ((log.step > 0) && (i %% log.step == 0)) {
      message("Iteration ", i, ", gene ", marker.info$Gene, ". Max uncertainty: ", 1 - min(confidence),
              ", mean uncertainty: ", 1 - mean.confidence, ", min #positive: ", min(n.pos.markers),
              ", min #negative: ", min(n.neg.markers))
    }

    if ((min(min(n.neg.markers), min(n.pos.markers)) >= max.markers.per.type) || (1 - min(confidence)) < max.uncertainty)
      break

    markers.per.type %<>% updateMarkersPerType(marker.info=marker.info)

    if (debug) print(marker.info)
  }

  return(marker.list)
}

## can add small bonus for negative markers, which target a log of negative cells even though there are no positive expression yet
## Looks like we don't have proper negative markers in the DE set now. append_specifisity metrics and sort by negative specificity instead of Z
