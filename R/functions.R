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

#' @inheritDotParams diffuseGraph fading fading.const verbose tol score.fixing.threshold
diffuseScorePerType <- function(scores.per.type, graph, parents, cbs.per.type, verbose, n.cores=1, ...) {
  plapply(parents, function(p)
    diffuseGraph(igraph::induced_subgraph(graph, cbs.per.type[[p]]),
                 scores=scores.per.type[[p]], verbose=verbose, ...),
    verbose=(verbose > 0), n.cores=n.cores)
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
#' @param score.info cell type scores from `getMarkerScoreInfo` function. Re-estimated if NULL
#' @param clusters cluster assignment of data. Used to expand annotation on these clusters.
#' @param verbose verbosity level (from 0 to 2)
#' @inheritDotParams diffuseScorePerType
#'
#' @export
assignCellsByScores <- function(graph, clf.data, score.info=NULL, clusters=NULL, verbose=0, uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75), ...) {
  if (!is.null(clusters)) {
    clusters <- as.factor(clusters)
  }

  if (is.null(score.info)) {
    score.info <- getMarkerScoreInfo(clf.data)
  }

  scores <- getMarkerScoresPerCellType(clf.data, score.info=score.info)

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
  c.ann.filt <- c.ann
  ann.by.level <- list()
  ann.filt.by.level <- list()
  scores.by.level <- list()
  scores.posterior <- scores

  possible.ann.levels <- c()
  for (pl in 1:length(subtypes.per.depth.level)) {
    if (verbose > 0) message("Level ", pl, "...")

    c.subtypes.per.parent <- subtypes.per.depth.level[[pl]]
    possible.ann.levels %<>% c(unlist(c.subtypes.per.parent)) %>% unique()

    c.parents <- names(c.subtypes.per.parent) %>% setNames(., .)
    cbs.per.type <- lapply(c.parents, function(p) names(c.ann)[c.ann == p]) %>%
      .[sapply(., length) > 0]
    c.parents %<>% .[names(cbs.per.type)]

    scores.per.type <- lapply(c.parents, function(p) scores[cbs.per.type[[p]], c.subtypes.per.parent[[p]]])

    if (!is.null(graph)) {
      scores.per.type %<>% diffuseScorePerType(graph, c.parents, cbs.per.type, verbose=(verbose > 1), ...)

      for (cs in scores.per.type) {
        scores.posterior[rownames(cs), colnames(cs)] <- cs
      }
    }

    res.ann <- lapply(scores.per.type, annotationFromScores, clusters) %>% Reduce(c, .)
    res.ann.filt <- filterAnnotationByUncertainty(res.ann, scores.posterior[,possible.ann.levels], score.info=score.info,
                                                  cur.types=unique(res.ann), clusters=clusters, thresholds=uncertainty.thresholds)

    c.ann[names(res.ann)] <- res.ann
    c.ann.filt[names(res.ann.filt)] %<>% is.na() %>% ifelse(NA, res.ann.filt)

    level.name <- paste0("l", pl)
    ann.by.level[[level.name]] <- c.ann
    ann.filt.by.level[[level.name]] <- c.ann.filt
    scores.by.level[[level.name]] <- scores.posterior[,possible.ann.levels]
    if (verbose > 0) message("Done")
  }

  return(list(annotation=ann.by.level, scores=scores.by.level, annotation.filt=ann.filt.by.level))
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

#' @export
getMarkerScoreInfo <- function(clf.data, ...) {
  lapply(clf.data$marker.list, getCellTypeScoreInfo, clf.data$cm, ...)
}

mergeScores <- function(score.name, scores, aggr.func=c) {
  names(scores[[1]][[score.name]]) %>% setNames(., .) %>% lapply(function(n)
    lapply(scores, `[[`, c(score.name, n)) %>% Reduce(aggr.func, .))
}

#' @export
mergeScoreInfos <- function(score.infos, verbose=F) {
  names(score.infos[[1]]) %>% setNames(., .) %>% lapply(mergeScores, score.infos)
  # names(score.infos[[1]]) %>% setNames(., .) %>%
  #   plapply(function(tn) names(score.infos[[1]][[1]]) %>% setNames(., .) %>% lapply(function(sn)
  #     lapply(score.infos, function(si) si[[tn]][[sn]]) %>% Reduce(c, .)), verbose=verbose)
}

mergeAnnotationInfos <- function(ann.infos) {
  return(list(
    annotation=mergeScores("annotation", ann.infos),
    scores=mergeScores("scores", ann.infos, aggr.func=rbind),
    annotation.filt=mergeScores("annotation.filt", ann.infos)
  ))
}

getCellTypeScoreInfo <- function(markers, tf.idf, aggr=T) {
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

  res <- list(scores.raw=scores, mult=score.mult, max.positive=max.positive.expr, scores=(scores * score.mult))
  return(res)
}

getCellTypeScore <- function(markers, tf.idf, aggr=T, do.multiply=T, check.gene.presence=T) {
  .Deprecated("getCellTypeScoreInfo(...)$scores")
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

normalizeScores <- function(scores, min.val=1e-10) {
  scores[rowSums(scores) < 1e-10,] <- 1
  scores  %<>% `/`(rowSums(.))
  return(scores)
}

#' Return initial scores of each cell type for each cell
#'
#' @param clf classification data from `getClassificationData`
#' @param score.info pre-calculated score info from `getMarkerScoreInfo`
#' @param aggr should individual gene scores be aggregated per cell type? If `FALSE`,
#' returns list of data.frames, showing scores of each gene for each cell.
#' Useful for debugging list of markers.
#'
#' @return data.frame with rows corresponding to cells and columns corresponding to cell types.
#' Values are cell type scores, normalized per level of hierarchy
#'
#' @export
getMarkerScoresPerCellType <- function(clf, score.info=NULL, aggr=T) {
  if (is.null(score.info)) {
    score.info <- getMarkerScoreInfo(clf, aggr=aggr)
  }

  scores <- lapply(score.info, `[[`, "scores")
  if (!aggr)
    return(lapply(scores, as.matrix) %>% lapply(as.data.frame, optional=T))

  clf.nodes <- classificationTreeToDf(clf$classification.tree)
  scores %<>% as.data.frame(optional=T)

  for (nodes in split(clf.nodes$Node, clf.nodes$PathLen)) {
    scores[, nodes]  %<>% normalizeScores()
  }

  return(scores)
}

## Utils

appendHierarchyBranch <- function(branch, parent.name) {
  is.leaf <- (sapply(branch, is.character) & (sapply(branch, length) == 1))
  if (any((sapply(names(branch), nchar) > 0) != !is.leaf))
    stop("Can't parse ", parent.name, " branch: only lists must have names")

  current.nodes <- lapply(c(branch[is.leaf], names(branch[!is.leaf])), function(x)
    list(list(expressed=c(), not_expressed=c(), parent=parent.name)) %>% setNames(x))

  sub.branches <- lapply(names(branch)[!is.leaf], function(n) appendHierarchyBranch(branch[[n]], n))
  return(c(current.nodes, unlist(sub.branches, recursive=F)))
}

getAnnotationPerParent <- function(clf.tree, annotation) {
  classificationTreeToDf(clf.tree) %>% split(.$Parent) %>% lapply(function(df)
    mergeAnnotationToLevel(df$PathLen[1], annotation, clf.tree) %>% .[. %in% df$Node])
}

#' @export
hierarchyToClassificationTree <- function(hierarchy) {appendHierarchyBranch(hierarchy, "root") %>% unlist(recursive=F) %>% createClassificationTree()}

mergeAnnotationToLevel <- function(level, annotation, classification.tree) {
  parent.types <- classificationTreeToDf(classification.tree) %$% Node[PathLen == level] %>% unique()
  if (length(parent.types) == 0) {
    warning("Nothing to merge at level ", level)
    return(annotation)
  }

  if (is.factor(annotation)) {
    annotation <- setNames(as.character(annotation), names(annotation))
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

getAnnotationPerCluster <- function(annotation, clusters) {
  ann.per.clust <- table(annotation, clusters[names(annotation)])
  if (any(colSums(ann.per.clust > 0) != 1))
    stop("Some clusters match to multiple cell types")

  ann.per.clust %<>% apply(2, which.max) %>% rownames(ann.per.clust)[.] %>%
    setNames(colnames(ann.per.clust))

  return(ann.per.clust)
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
