estimatePositiveMarkerScoreChange <- function(cell.type, annotation, cm.norm, de.genes, pos.scores, neg.scores, sum.scores, balance.cell.types=T) {
  if (length(de.genes) == 0)
    stop("de.genes are empty for ", cell.type)

  sum.scores %<>% pmax(1e-30)

  c.exprs <- as.matrix(cm.norm[, de.genes])

  c.neg.scores <- (1 - neg.scores[, cell.type])
  c.pos.scores <- pos.scores[, cell.type] * c.neg.scores
  ds <- c.exprs * c.neg.scores
  ds <- (ds + c.pos.scores) / (ds + sum.scores) - (c.pos.scores / sum.scores)

  if (balance.cell.types) {
    scores.per.type <- split(names(annotation), annotation) %>%
      sapply(function(cbs) Matrix::colMeans(ds[cbs,]))
    ds <- scores.per.type[,cell.type]
    scores.per.type[,cell.type] <- 0
    ds <- ds - rowSums(scores.per.type)
  } else {
    mask <- (annotation == cell.type)
    # pos.ids <- which(mask)
    # neg.ids <- which(!mask)
    ds <- Matrix::colSums(ds[mask,]) - Matrix::colSums(ds[!mask,])
  }

  # TODO: also need to account on possible change of positive max values

  ds <- sort(ds, decreasing=T)
  return(tibble::tibble(Gene=names(ds), Score=ds, Type=cell.type, MT="expressed"))
}

estimateNegativeMarkerScoreChange <- function(cell.type, annotation, cm.norm, de.genes, pos.scores,
                                              neg.scores, sum.scores, max.pos.scores, balance.cell.types=T) {
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

  # ds <- ((c.pos.scores + ds) / pmax(sum.scores + ds, 1e-30)) - (c.pos.scores / sum.scores)
  # ds <- Matrix::colSums(ds[pos.ids,]) - Matrix::colSums(ds[neg.ids,])

  new.scores <- c.pos.scores * (1 - new.neg.scores)
  ds <- new.scores / pmax(sum.scores + new.scores - old.scores, 1e-30) - (old.scores / pmax(sum.scores, 1e-30))

  if (balance.cell.types) {
    scores.per.type <- split(names(annotation), annotation) %>%
      sapply(function(cbs) Matrix::colMeans(ds[cbs,, drop=F]))
    ds <- scores.per.type[,cell.type]
    scores.per.type[,cell.type] <- 0
    ds <- ds - rowSums(scores.per.type)
  } else {
    ds <- Matrix::colSums(ds[pos.ids,, drop=F]) - Matrix::colSums(ds[neg.ids,, drop=F])
  }

  # ds <- estimateDNegativeScores(ds, c.pos.scores, sum.scores, mask) %>%
  #   setNames(colnames(ds))

  ds <- sort(ds, decreasing=T)
  return(tibble::tibble(Gene=names(ds), Score=ds, Type=cell.type, MT="not_expressed"))
}

prepareMarkerScoringInfo <- function(marker.list, cm.norm) {
  res <- lapply(marker.list, getCellTypeScoreFull, cm.norm)

  return(list(
    sum.scores=sapply(res, `[[`, "sm") %>% rowSums(),
    pos.scores=sapply(res, `[[`, "scores"),
    neg.scores=1 - sapply(res, `[[`, "mult"),
    max.pos.scores=sapply(res, `[[`, "max.positive")
  ))
}

initializeMarkerList <- function(pos.markers.per.type, cm.norm, debug=F) {
  marker.list <- lapply(pos.markers.per.type, function(x) list(expressed=list(), not_expressed=list()))
  subtypes.left <- names(pos.markers.per.type)
  for (i in 1:length(subtypes.left)) {
    s.info <- prepareMarkerScoringInfo(marker.list, cm.norm)

    marker.score <- subtypes.left %>% lapply(function(ct)
      estimatePositiveMarkerScoreChange(ct, annotation_adj, cm.norm, pos.markers.per.type[[ct]], s.info$pos.scores,
                                        s.info$neg.scores, s.info$sum.scores, balance.cell.types=F)
    ) %>% Reduce(rbind, .) %>% .[which.max(.$Score),]

    if (marker.score$Score < 0) {
      marker.score <- subtypes.left %>% lapply(function(ct)
        estimatePositiveMarkerScoreChange(ct, annotation_adj, cm.norm, pos.markers.per.type[[ct]], s.info$pos.scores,
                                          s.info$neg.scores, s.info$sum.scores)
      ) %>% Reduce(rbind, .) %>% .[which.max(.$Score),]

      if (marker.score$Score < 0) {
        warning("Negative score for ", n, ": ", marker.score$Score, ". Consider changing extending marker list.")
      }
    }

    marker.list[[c(marker.score$Type, marker.score$MT)]] %<>% c(marker.score$Gene)
    subtypes.left %<>% setdiff(marker.score$Type)

    if (debug) print(marker.score)
  }

  return(marker.list)
}

## Now cell types are not balanced by number of cells
## can add small bonus for negative markers, which target a log of negative cells even though there are no positive expression yet
## Looks like we don't have proper negative markers in the DE set now. append_specifisity metrics and sort by negative specificity instead of Z!!!
## Estimate score per marker for all cell types together? Looks complicated.
