#' @useDynLib CellAnnotatoR, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom magrittr %<>% %$% %>%
NULL

## Correct unloading of the library
.onUnload <- function (libpath) {
  library.dynam.unload("CellAnnotatoR", libpath)
}

plapply <- function(..., n.cores=1, verbose=F) {
  if (verbose)
    return(pbapply::pblapply(..., cl=n.cores))

  if (n.cores == 1)
    return(lapply(...))

  return(parallel::mclapply(..., mc.cores=n.cores))
}

expandAnnotationToClusters <- function(scores, clusters) {
  clusters <- droplevels(clusters[rownames(scores)])
  ann.update <- split(1:length(clusters), clusters) %>%
    sapply(function(cids) colSums(scores[cids, ,drop=F]) %>% which.max() %>% names()) %>%
    .[clusters] %>% setNames(names(clusters))

  return(ann.update)
}

annotationFromScores <- function(scores, clusters=NULL) {
  if (!is.null(clusters))
    return(expandAnnotationToClusters(scores, clusters))

  return(colnames(scores)[apply(scores, 1, which.max)] %>% setNames(rownames(scores)))
}

#' @export
getAnnotationConfidence <- function(annotation, scores) {
  mapply(function(i,j) scores[i, j], 1:nrow(scores), match(annotation, colnames(scores))) %>%
    setNames(rownames(scores))
}

#' Run diffusion on graph
#'
#' @param fading fading level for graph diffusion
#' @param fading.const constant in exponent for graph diffusion
#' @param score.fixing.threshold threshold for a label to be considered true
#' @param verbose print progress bar
#' @param tol tolerance for diffusion stopping
diffuseGraph <- function(graph, scores, fading=10, fading.const=0.5, score.fixing.threshold=0.8,
                         verbose=FALSE, max.iters=1000, tol=1e-3) {
  cbs <- igraph::V(graph)$name
  if (length(cbs) == 0)
    return(NULL)

  scores <- as.matrix(scores[cbs, ])
  scores[rowSums(scores) < 1e-5, ] <- 1
  scores %<>% `/`(rowSums(.))

  if (any(is.na(scores)))
    stop("NAs in scores")

  edges <- igraph::as_edgelist(graph)

  is.fixed <- (apply(scores, 1, max) > score.fixing.threshold)

  if (nrow(edges) == 0)
    return(scores)

  edge.weights <- igraph::edge.attributes(graph)$weight
  res <- conos:::smoothMatrixOnGraph(edges, edge.weights, scores, is.label.fixed=is.fixed, max_n_iters=max.iters,
                                     diffusion_fading=fading, diffusion_fading_const=fading.const, verbose=verbose,
                                     tol=tol, normalize=T)
  return(res)
}

#' Assign cell types for each cell based on type scores. Optionally uses `clusters` to expand annotation.
#'
#' @param graph cell graph from Seurat, Pagoda2 or some other tool. Can be either in igraph or adjacency matrix format.
#'    Use `graph=NULL` to skip graph diffusion step and get raw score annotation (useful when debug marker genes).
#' @param clf.data classification data from `getClassificationData`
#' @param scores cell type scores from `getMarkerScoresPerCellType` function. Re-estimated if NULL
#' @param clusters cluster assignment of data. Used to expand annotation on these clusters.
#' @param verbose verbosity level (from 0 to 2)
#' @inheritDotParams diffuseGraph fading fading.const verbose tol score.fixing.threshold
#'
#' @export
assignCellsByScores <- function(graph, clf.data, scores=NULL, clusters=NULL, verbose=0, n.cores=1, ...) {
  if (!is.null(clusters)) {
    clusters <- as.factor(clusters)
  }

  if (is.null(scores)) {
    scores <- getMarkerScoresPerCellType(clf.data)
  }

  if (!is.null(graph)) {
    if ((is(graph, "Matrix") || is(graph, "matrix")) && ncol(graph) == nrow(graph)) {
      graph <- igraph::graph_from_adjacency_matrix(graph, weighted=T)
    } else if (!is(graph, "igraph")) {
      stop("Unknown graph format. Only adjacency matrix or igraph are supported")
    }
  }

  subtypes.per.depth.level <- classificationTreeToDf(clf.data$classification.tree) %$%
    split(., PathLen) %>% lapply(function(df) split(df$Node, df$Parent)) %>% .[order(names(.))]

  c.ann <- rep("root", nrow(scores)) %>% setNames(rownames(scores))
  ann.by.level <- list()
  scores.by.level <- list()
  scores.posterior <- scores

  possible.ann.levels <- c()
  for (pl in 1:length(subtypes.per.depth.level)) {
    if (verbose > 0) message("Level ", pl, "...")

    c.subtypes.per.parent <- subtypes.per.depth.level[[pl]]
    possible.ann.levels %<>% c(unlist(c.subtypes.per.parent)) %>% unique()

    c.parents <- names(c.subtypes.per.parent) %>% setNames(., .)
    cbs.per.typs <- lapply(c.parents, function(p) names(c.ann)[c.ann == p]) %>%
      .[sapply(., length) > 0]
    c.parents %<>% .[names(cbs.per.typs)]

    scores.per.type <- lapply(c.parents, function(p) scores[cbs.per.typs[[p]], c.subtypes.per.parent[[p]]])

    if (!is.null(graph)) {
      scores.per.type <- plapply(c.parents, function(p) {
        diffuseGraph(igraph::induced_subgraph(graph, cbs.per.typs[[p]]),
                     scores=scores.per.type[[p]], verbose=(verbose > 1), ...)
      }, verbose=(verbose > 0), n.cores=n.cores)

      for (cs in scores.per.type) {
        scores.posterior[rownames(cs), colnames(cs)] <- cs
      }
    }

    res.ann <- lapply(scores.per.type, annotationFromScores, clusters) %>% Reduce(c, .)
    c.ann[names(res.ann)] <- res.ann

    level.name <- paste0("l", pl)
    ann.by.level[[level.name]] <- c.ann
    scores.by.level[[level.name]] <- scores.posterior[,possible.ann.levels]
    if (verbose > 0) message("Done")
  }

  return(list(annotation=ann.by.level, scores=scores.by.level))
}


## Score assignment

#' @export
normalizeTfIdfWithFeatures <- function(cm, max.quantile=0.95, max.smooth=1e-10) {
  cm@x <- cm@x / rep(Matrix::colSums(cm), diff(cm@p))
  cm <- Matrix::t(cm)

  max.vals <- split(cm@x, rep(1:(length(cm@p)-1), diff(cm@p)))[paste0(1:ncol(cm))]
  max.vals[is.null(max.vals)] <- c()
  max.vals %<>% sapply(quantile, max.quantile) %>% `+`(max.smooth)

  cm@x <- cm@x / rep(max.vals, diff(cm@p))

  tf.idf <- cm
  idf.weights <- log(1 + nrow(tf.idf) / (Matrix::colSums(tf.idf > 0) + 1))
  tf.idf@x <- tf.idf@x * rep(idf.weights, diff(tf.idf@p))
  return(tf.idf)
}

getCellTypeScoreFull <- function(markers, tf.idf, aggr=T) {
  if (length(markers$expressed) == 0) {
    if (aggr) {
      scores <- setNames(rep(0, nrow(tf.idf)), rownames(tf.idf))
    } else {
      scores <- matrix(0, nrow=nrow(tf.idf), ncol=0, dimnames=list(rownames(tf.idf), character()))
    }

    max.positive.expr <- setNames(rep(0, nrow(tf.idf)), rownames(tf.idf))
  } else {
    m.expressed <- intersect(markers$expressed, colnames(tf.idf))
    c.submat <- tf.idf[, m.expressed, drop=F]
    c.submat.t <- Matrix::t(c.submat)
    scores <- if (aggr) Matrix::colSums(c.submat.t) else c.submat
    max.positive.expr <- apply(c.submat.t, 2, max)
  }

  not.expressed.genes <- intersect(markers$not_expressed, colnames(tf.idf))
  if (length(not.expressed.genes) == 0) {
    score.mult <- setNames(rep(1, nrow(tf.idf)), rownames(tf.idf))
  } else {
    max.negative.expr <- apply(tf.idf[, not.expressed.genes, drop=F], 1, max)
    score.mult <- pmax(max.positive.expr - max.negative.expr, 0) / max.positive.expr
    score.mult[is.na(score.mult)] <- 0
  }

  res <- list(scores=scores, mult=score.mult, max.positive=max.positive.expr, sm=(scores * score.mult))
  return(res)
}

getCellTypeScore <- function(markers, tf.idf, aggr=T, do.multiply=T, check.gene.presence=T) {
  if (length(markers$expressed) == 0) {
    if (aggr) {
      res <- setNames(rep(0, nrow(tf.idf)), rownames(tf.idf))
    } else {
      res <- matrix(0, nrow=nrow(tf.idf), ncol=0, dimnames=list(rownames(tf.idf), character()))
    }

    if (do.multiply)
      return(res)

    score.mult <- setNames(rep(1, nrow(tf.idf)), rownames(tf.idf))
    return(list(scores=res, mult=score.mult))
  }

  m.expressed <- if (check.gene.presence) intersect(markers$expressed, colnames(tf.idf)) else markers$expressed # TODO: Remove this?

  c.submat <- tf.idf[, m.expressed, drop=F]
  c.submat.t <- Matrix::t(c.submat)
  scores <- if (aggr) Matrix::colSums(c.submat.t) else c.submat

  not.expressed.genes <- intersect(markers$not_expressed, colnames(tf.idf))
  if (length(not.expressed.genes) == 0) {
    if (do.multiply)
      return(scores)

    score.mult <- setNames(rep(1, nrow(tf.idf)), rownames(tf.idf))
    return(list(scores=res, mult=score.mult))
  }

  max.positive.expr <- apply(c.submat.t, 2, max)
  max.negative.expr <- apply(tf.idf[, not.expressed.genes, drop=F], 1, max)

  score.mult <- pmax(max.positive.expr - max.negative.expr, 0) / max.positive.expr
  score.mult[is.na(score.mult)] <- 0

  if (!do.multiply)
    return(list(scores=scores, mult=score.mult))

  return(scores * score.mult)
}

#' Return initial scores of each cell type for each cell
#'
#' @param clf classification data from `getClassificationData`
#' @param aggr should individual gene scores be aggregated per cell type? If `FALSE`,
#' returns list of data.frames, showing scores of each gene for each cell.
#' Useful for debugging list of markers.
#'
#' @return data.frame with rows corresponding to cells and columns corresponding to cell types.
#' Values are cell type scores, normalized per level of hierarchy
#'
#' @export
getMarkerScoresPerCellType <- function(clf, aggr=T) {
  res <- lapply(clf$marker.list, getCellTypeScore, clf$cm, aggr=aggr)

  if (!aggr)
    return(lapply(res, as.matrix) %>% lapply(as.data.frame, optional=T))

  clf.nodes <- classificationTreeToDf(clf$classification.tree)
  res %<>% as.data.frame(optional=T)

  for (nodes in split(clf.nodes$Node, clf.nodes$PathLen)) {
    res[rowSums(res[, nodes]) < 1e-10, nodes] <- 1
    res[, nodes]  %<>% `/`(rowSums(.))
  }

  return(res)
}

## Utils

mergeAnnotationToLevel <- function(level, annotation, classification.tree) {
  parent.types <- classificationTreeToDf(classification.tree) %$% Node[PathLen == level] %>% unique()
  if (length(parent.types) == 0) {
    warning("Nothing to merge at level ", level)
    return(annotation)
  }

  type.map <- unique(annotation) %>% setNames(., .)
  for (pt in parent.types) {
    for (st in getAllSubtypes(pt, classification.tree)) {
      type.map[st] = pt
    }
  }

  return(setNames(type.map[annotation], names(annotation)))
}

#' @export
mergeAnnotationByLevels <- function(annotation, classification.tree) {
  max.level <- classificationTreeToDf(classification.tree)$PathLen %>% max()
  anns <- 1:max.level %>% setNames(., paste0("l", .)) %>%
    lapply(mergeAnnotationToLevel, annotation, classification.tree)
  return(anns)
}

classificationTreeToDf <- function(classification.tree) {
  igraph::V(classification.tree)$name %>% .[. != "root"] %>%
    lapply(function(cn) {
      p <- igraph::all_simple_paths(classification.tree, "root", to=cn, mode="out")[[1]]$name
      tibble::tibble(Parent=tail(p, 2)[1], Node=tail(p, 1), PathLen=length(p) - 1)
    }) %>%  Reduce(rbind, .)
}

getAllSubtypes <- function(parent.type, classification.tree, max.depth=NULL) {
  paths <- igraph::dfs(classification.tree, parent.type, neimode="out", unreachable=F, dist=T)
  paths <- if (!is.null(max.depth)) names(which(paths$dist <= max.depth)) else names(paths$order)
  return(paths %>% .[!is.na(.)] %>% .[. != parent.type])
}

## Validation

rateMetricssPerType <- function(ann.res, ann.reference) {
  ann.reference <- ann.reference[!is.na(ann.reference)]
  ann.res <- ann.res[names(ann.reference)]

  res <- unique(ann.reference) %>% lapply(function(v){
    obs.pos <- (ann.res == v)
    real.pos <- (ann.reference == v)
    tpr <- sum(obs.pos & real.pos) / sum(real.pos)
    tnr <- sum(!obs.pos & !real.pos) / sum(!real.pos)
    precision <- sum(obs.pos & real.pos) / pmax(sum(obs.pos), 0.1)
    return(c(TPR=tpr, TNR=tnr, Precision=precision))
  }) %>% as.data.frame() %>% t() %>% magrittr::set_rownames(unique(ann.reference)) %>%
    as.data.frame()
  return(res)
}
