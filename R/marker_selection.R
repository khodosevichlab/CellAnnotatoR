aggregateScoreChangePerGene <- function(d.scores, annotation, marker.type, cell.type, target.type=NULL, balance.cell.types=T, self.mult=1) {
  aggr.func <- if (balance.cell.types) Matrix::colMeans else Matrix::colSums

  cbs.per.type <- split(names(annotation), annotation)
  scores.per.type <- cbs.per.type %>%
    lapply(function(cbs) aggr.func(d.scores[cbs,, drop=F])) %>%
    Reduce(cbind, .) %>% `colnames<-`(names(cbs.per.type)) %>% `*`(-1)
  # scores.per.type[,cell.type] %<>% `*`(-self.mult)
  scores.per.type[,cell.type] %<>% `*`(-1)
  if (!is.null(target.type)) {
    scores.per.type[,target.type] %<>% `*`(self.mult)
  }
  d.scores <- rowSums(scores.per.type) %>% sort(decreasing=T)
  res <- tibble::tibble(Gene=names(d.scores), Score=d.scores, Type=cell.type, MT=marker.type)
  if (!is.null(target.type)) {
    res$ScoreTarget <- scores.per.type[res$Gene, target.type]
  }

  return(res)
}

estimatePositiveMarkerScoreChange <- function(cell.type, annotation, cm.norm, de.genes, pos.scores, neg.scores, sum.scores, aggr=T, ...) {
  if (length(de.genes) == 0)
    stop("de.genes are empty for ", cell.type)

  sum.scores %<>% pmax(1e-30)

  c.exprs <- as.matrix(cm.norm[, de.genes, drop=F])

  c.neg.scores <- (1 - neg.scores[, cell.type])
  c.pos.scores <- pos.scores[, cell.type] * c.neg.scores
  ds <- c.exprs * c.neg.scores
  ds <- (ds + c.pos.scores) / (ds + sum.scores) - (c.pos.scores / sum.scores)

  if (!aggr)
    return(ds)

  return(aggregateScoreChangePerGene(ds, annotation, marker.type="expressed", cell.type=cell.type, ...))
}

estimateNegativeMarkerScoreChange <- function(cell.type, annotation, cm.norm, de.genes, pos.scores, neg.scores,
                                              sum.scores, max.pos.scores, max.neg.score.per.pos.type=0.1, aggr=T, ...) {
  if (length(de.genes) == 0)
    stop("de.genes are empty for ", cell.type)

  sum.scores %<>% pmax(1e-30)
  max.pos.scores <- pmax(max.pos.scores[,cell.type], 1e-30)
  mask <- (annotation == cell.type)
  pos.ids <- which(mask)
  neg.ids <- which(!mask)

  c.exprs <- as.matrix(cm.norm[, de.genes, drop=F])

  c.pos.scores <- pos.scores[,cell.type]
  c.pos.scores[pos.ids] <- sum.scores[pos.ids] # In ideal world, all positive score for a cell came from the proper cell type

  new.neg.scores <- estimateNewNegativeScores(c.exprs, max.pos.scores, neg.scores[,cell.type]) %>%
    `dimnames<-`(dimnames(c.exprs))

  old.scores <- c.pos.scores * (1 - neg.scores[,cell.type])
  ds <- (c.pos.scores * (1 - new.neg.scores) - old.scores)

  new.scores <- c.pos.scores * (1 - new.neg.scores)
  ds <- new.scores / pmax(sum.scores + new.scores - old.scores, 1e-30) - (old.scores / pmax(sum.scores, 1e-30))

  if (!aggr)
    return(ds)

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
        warning("Negative score for ", marker.score$Type, ": ", round(marker.score$Score, 4), ". Consider extending marker list.\n")
      }
    }

    marker.list[[c(marker.score$Type, marker.score$MT)]] %<>% c(marker.score$Gene)
    subtypes.left %<>% setdiff(marker.score$Type)

    if (debug) print(marker.score)
  }

  return(marker.list)
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

preSelectMarkersForType <- function(de.df, min.pos.markers=5, max.pos.markers=50, min.pos.specificity=0.2, min.pos.expression.frac=0.1,
                                    pos.expression.frac.weight=0.2, max.neg.expression.frac=0.1) {
  de.pos <- de.df[de.df$Z > 0, ]
  if (sum(de.pos$Specificity > min.pos.specificity) < min.pos.markers) {
    pos.markers <- de.pos$Gene[order(de.pos$Specificity, decreasing=T)[1:min.pos.markers]]
  } else {
    de.pos %<>% .[.$Specificity > min.pos.specificity,]
    if (sum(de.pos$ExpressionFraction > min.pos.expression.frac) < min.pos.markers) {
      pos.markers <- de.pos$Gene[order(de.pos$ExpressionFraction, decreasing=T)[1:min.pos.markers]]
    } else {
      pos.markers <- de.pos %>% .[.$ExpressionFraction > min.pos.expression.frac,] %>%
        .[order(.$Specificity + .$ExpressionFraction * pos.expression.frac.weight, decreasing=T),] %>% .$Gene %>%
        .[1:min(length(.), max.pos.markers)]
    }
  }
  # Negative: ExpressionFraction < 0.1 && Z < 0 && top by specificity (or > 0.95)

  neg.markers <- de.df %>% .[(.$Z < 0) & (.$ExpressionFraction < max.neg.expression.frac), ] %>% .$Gene
  return(list(positive=pos.markers, negative=neg.markers))
}

#' @inheritDotParams preSelectMarkersForType
#' @export
preSelectMarkerCandidates <- function(de.info, ...) {
  markers.per.type <- lapply(de.info, preSelectMarkersForType, ...)
  markers.per.type <- c("positive", "negative") %>% setNames(., .) %>%
    lapply(function(n) lapply(markers.per.type, `[[`, n))

  if (any(unlist(lapply(markers.per.type, sapply, length)) == 0)) {
    warning("Some cell types don't have positive or negative markers")
  }

  return(markers.per.type)
}

#### Version 2

getTopNegativeGenes <- function(pos.gene, cell.type, cm.norm, annotation, markers.per.type, s.info, pos.score.changes, n.neg.genes, score.change.threshold) {
  c.max.scores <- pmax(s.info$max.pos.scores[,cell.type], cm.norm[, pos.gene])
  cm.norm.neg <- cm.norm[, markers.per.type$negative[[cell.type]]]
  neg.scores <- estimateNewNegativeScores(cm.norm.neg, c.max.scores, s.info$neg.scores[,cell.type]) %>%
    `dimnames<-`(dimnames(cm.norm.neg))

  neg.score.info <- ((1 - neg.scores) * pos.score.changes[,pos.gene]) %>%
    aggregateScoreChangePerGene(annotation, "both", cell.type)
  neg.score.base <- neg.score.info %>% .[which.max(.$Score),]
  top.neg.ids <- neg.score.info %$% Gene[order(Score, decreasing=T)[1:n.neg.genes]] %>%
    match(colnames(neg.scores))

  d.score.per.neg <- lapply(top.neg.ids, function(base.id)
    1:ncol(neg.scores) %>%
      lapply(function(i) (1 - pmax(neg.scores[,base.id], neg.scores[,i]))) %>%
      Reduce(cbind, .) %>% `colnames<-`(colnames(neg.scores)) %>%
      `*`(pos.score.changes[,pos.gene]))

  res <- lapply(d.score.per.neg, aggregateScoreChangePerGene, annotation, "both", cell.type) %>%
    lapply(`[`,1,) %>% Reduce(rbind, .) %>% dplyr::rename(NGene1=Gene) %>%
    dplyr::mutate(NGene2=colnames(neg.scores)[top.neg.ids]) %>%  .[which.max(.$Score),]

  if ((res$Score - neg.score.base$Score) / abs(neg.score.base$Score) < score.change.threshold) {
    res <- neg.score.base %>% rename(NGene1=Gene) %>% dplyr::mutate(NGene2=NA)
  }

  return(dplyr::mutate(res, PGene=pos.gene))
}

getNextMarkers <- function(cell.type, cm.norm, annotation, marker.list, markers.per.type, n.pos.genes=3,
                           n.neg.genes=10, score.change.threshold=0.05, verbose=F, n.cores=1) {
  s.info <- lapply(marker.list, getCellTypeScoreInfo, cm.norm) %>% prepareMarkerScoringInfo()
  pos.score.changes <- estimatePositiveMarkerScoreChange(cell.type, annotation, cm.norm, markers.per.type$positive[[cell.type]],
                                                         s.info$pos.scores, s.info$neg.scores, s.info$sum.scores, aggr=F)

  pos.score.changes.aggr <- aggregateScoreChangePerGene(pos.score.changes, annotation, marker.type="expressed", cell.type=cell.type)
  top.pos.genes <- pos.score.changes.aggr %$% Gene[order(Score, decreasing=T)[1:n.pos.genes]]

  pos.score <- pos.score.changes.aggr %>% .[which.max(.$Score),]
  res.score <- plapply(top.pos.genes, getTopNegativeGenes, cell.type, cm.norm, annotation, markers.per.type, s.info,
                       pos.score.changes, n.neg.genes=n.neg.genes, score.change.threshold=score.change.threshold,
                       verbose=verbose, n.cores=max(min(n.cores, n.pos.genes), 1)) %>%
    Reduce(rbind, .) %>% .[which.max(.$Score),]

  if (verbose) {
    message("Neg. score: ", round(res.score$Score, 3), ", Pos. score: ", round(pos.score$Score, 3))
  }

  if ((res.score$Score - pos.score$Score) / abs(pos.score$Score) < score.change.threshold) {
    res.score <- pos.score %>% dplyr::rename(PGene=Gene) %>% dplyr::mutate(NGene1=NA, NGene2=NA)
  }

  if (verbose) {
    message("Score: ", round(res.score$Score, 3), ", PG: ", res.score$PGene, ", NG1:", res.score$NGene1, ", NG2:", res.score$NGene2)
  }

  return(list(expressed=res.score$PGene, not_expressed=unlist(res.score[c("NGene1", "NGene2")]) %>% .[!is.na(.)]))
}

filterMarkerList <- function(gene, marker.list, cell.type, expr.type) {
  marker.list[[cell.type]][[expr.type]] %<>% setdiff(gene)
  return(marker.list)
}

getMeanConfidencePerType <- function(marker.list, cm.norm, annotation) {
  c.scores <- lapply(marker.list, getCellTypeScoreInfo, cm.norm) %>%
    lapply(`[[`, "scores") %>% as.data.frame(optional=T) %>% normalizeScores()
  confidence <- getAnnotationConfidence(annotation, c.scores)
  return(confidence %>% split(annotation[names(.)]) %>% sapply(mean))
}

filterMarkerListByScore <- function(marker.list, cm.norm, annotation, verbose=F, n.cores=1, do.recursive=T, change.threshold=1e-5) {
  mls.filt <- lapply(names(marker.list), function(ct) lapply(names(marker.list[[ct]]), function(et)
    lapply(marker.list[[ct]][[et]], filterMarkerList, marker.list, ct, et))) %>%
    unlist(recursive=F) %>% unlist(recursive=F)

  mean.conf.per.type <- getMeanConfidencePerType(marker.list, cm.norm, annotation)
  conf.per.ml <- plapply(mls.filt, getMeanConfidencePerType, cm.norm, annotation, verbose=verbose, n.cores=n.cores)

  mean.score.per.ml <- sapply(conf.per.ml, mean)
  t.ids <- which((mean.score.per.ml >= mean(mean.conf.per.type)) & (sapply(conf.per.ml, min) >= min(mean.conf.per.type)))
  if (length(t.ids) == 0)
    return(marker.list)

  best.id <- t.ids[which.max(mean.score.per.ml[t.ids])]
  change <- mean.score.per.ml[best.id] - mean(mean.conf.per.type)

  ml.res <- mls.filt[[best.id]]
  if (do.recursive && (change > change.threshold)) {
    ml.res <- filterMarkerListByScore(ml.res, cm.norm, annotation, verbose=verbose, n.cores=n.cores,
                                      do.recursive=do.recursive, change.threshold=change.threshold)
  }

  if (verbose) {
    message("Score improvement: ", mean.score.per.ml[best.id] - mean(mean.conf.per.type))
  }

  return(ml.res)
}

selectMarkersPerType <- function(cm.norm, marker.list, annotation, markers.per.type, max.iters=ncol(cm.norm), max.uncertainty=0.25, verbose=0, max.pos.markers=10, log.step=1, n.cores=1, refinement.period=10, ret.all=T) {
  mean.unc.per.type <- (1 - getMeanConfidencePerType(marker.list, cm.norm, annotation))
  for (i in 1:max.iters) {
    n.markers.per.cell <- sapply(marker.list, function(x) length(x$expressed))
    if (all(n.markers.per.cell >= max.pos.markers))
      break

    cell.type <- mean.unc.per.type %>% .[n.markers.per.cell[names(.)] < max.pos.markers] %>%
      which.max() %>% names()
    m.update <- getNextMarkers(cell.type, cm.norm, annotation, marker.list=marker.list, markers.per.type=markers.per.type,
                               verbose=(verbose > 1), n.cores=n.cores)

    marker.list.new <- marker.list
    marker.list.new[[cell.type]]$expressed %<>% union(m.update$expressed)
    marker.list.new[[cell.type]]$not_expressed %<>% union(m.update$not_expressed)

    markers.per.type %<>% updateMarkersPerType(marker.list=setNames(list(m.update), cell.type))

    mean.unc.per.type.new <- (1 - getMeanConfidencePerType(marker.list.new, cm.norm, annotation))
    # if (mean.unc.per.type.new[cell.type] < mean.unc.per.type[cell.type]) {
    if ((mean(mean.unc.per.type.new) < mean(mean.unc.per.type)) || (max(mean.unc.per.type.new) < max(mean.unc.per.type))) {
      mean.unc.per.type <- mean.unc.per.type.new
      marker.list <- marker.list.new
    }

    if (verbose && (log.step > 0) && (i %% log.step == 0)) {
      message("Iteration ", i, ". Max uncertainty: ", round(max(mean.unc.per.type), 3), ", mean uncertainty: ", round(mean(mean.unc.per.type), 3))
    }

    if (max(mean.unc.per.type) < max.uncertainty)
      break

    if ((refinement.period > 0) && (i %% refinement.period == 0)) {
      if (verbose) message("Refine markers...")
      marker.list %<>% filterMarkerListByScore(cm.norm, annotation, verbose=(verbose > 1), n.cores=n.cores)
    }
  }

  if (refinement.period != 0) {
    if (verbose) message("Refine markers...")
    marker.list %<>% filterMarkerListByScore(cm.norm, annotation, verbose=(verbose > 1), n.cores=n.cores)
  }

  if (ret.all)
    return(list(marker.list=marker.list, markers.per.type=markers.per.type))

  return(marker.list)
}


## can add small bonus for negative markers, which target a log of negative cells even though there are no positive expression yet