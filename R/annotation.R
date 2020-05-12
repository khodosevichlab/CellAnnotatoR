#' @useDynLib CellAnnotatoR, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom magrittr %<>% %$% %>%
NULL

## Correct unloading of the library
.onUnload <- function (libpath) {
  library.dynam.unload("CellAnnotatoR", libpath)
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

#' Diffuse Score per Type
#' @inheritDotParams diffuseGraph fading fading.const verbose tol score.fixing.threshold
diffuseScorePerType <- function(scores.per.type, graph, parents, cbs.per.type, verbose, n.cores=1, ...) {
  plapply(parents, function(p)
    diffuseGraph(igraph::induced_subgraph(graph, cbs.per.type[[p]]),
                 scores=scores.per.type[[p]], verbose=verbose, ...),
    verbose=(verbose > 0), n.cores=n.cores)
}

#' Diffuse Graph
#' @description Run diffusion on graph
#' @param graph graph to diffuse on
#' @param scores table of scores
#' @param fading fading level for graph diffusion
#' @param fading.const constant in exponent for graph diffusion
#' @param score.fixing.threshold threshold for a label to be considered certain. Such labels can't be changed during diffusion.
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

#' Assign Cells By Scores
#' @description Assign cell types for each cell based on type scores. Optionally uses `clusters` to expand annotation.
#'
#' @param graph cell graph from Seurat, Pagoda2 or some other tool. Can be either in igraph or adjacency matrix format.
#'    Use `graph=NULL` to skip graph diffusion step and get raw score annotation (useful when debug marker genes).
#' @param score.info cell type scores from `getMarkerScoreInfo` function. Re-estimated if NULL
#' @param clusters vector with cluster labels named by cell ids. Used to expand annotation on these clusters.
#' @param verbose verbosity level (from 0 to 2)
#' @param max.depth maximal depth for which annotation is done. Useful during manual marker selection
#' @inheritParams getMarkerScoreInfo
#' @inheritDotParams diffuseScorePerType
#' @return list with parameters:\itemize{
#'   \item{annotation: annotation per level}
#'   \item{scores: assignment scores per level}
#'   \item{annotation.filt: the same as annotation, but cells, which don't pass QC are assigned to NA class}
#' }
#'
#' @examples
#'   clf_data <- getClassificationData(cm, marker_path)
#'   ann_by_level <- assignCellsByScores(graph, clf_data, clusters=clusters)
#'
#' @export
assignCellsByScores <- function(graph, clf.data, score.info=NULL, clusters=NULL, verbose=0, uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75), max.depth=NULL, ...) {
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

    if (length(setdiff(igraph::V(graph)$name, rownames(scores))) > 0)
      stop("Not all cells from the graph are presented in clf.data")

    if (length(setdiff(rownames(scores), igraph::V(graph)$name)) > 0) {
      warning("Not all cells from the clf.data are presented in the graph. Omitting ",
              nrow(scores) - length(igraph::V(graph)$name), " cells")
      scores %<>% .[igraph::V(graph)$name,]
    }
  }

  subtypes.per.depth.level <- classificationTreeToDf(clf.data$classification.tree) %$%
    split(., PathLen) %>% lapply(function(df) split(df$Node, df$Parent)) %>% .[order(names(.))]

  max.depth <- if (is.null(max.depth)) length(subtypes.per.depth.level) else min(max.depth, length(subtypes.per.depth.level))
  c.ann <- rep("root", nrow(scores)) %>% setNames(rownames(scores))
  c.ann.filt <- c.ann
  ann.by.level <- list()
  ann.filt.by.level <- list()
  scores.by.level <- list()
  scores.posterior <- scores

  possible.ann.levels <- c()
  for (pl in 1:max.depth) {
    if (verbose > 0) message("Level ", pl, "...")

    c.subtypes.per.parent <- subtypes.per.depth.level[[pl]]
    possible.ann.levels %<>% c(unlist(c.subtypes.per.parent)) %>% unique()

    c.parents <- names(c.subtypes.per.parent) %>% setNames(., .)
    cbs.per.type <- lapply(c.parents, function(p) names(c.ann)[c.ann == p]) %>%
      .[sapply(., length) > 0]

    if (length(cbs.per.type) == 0) # some subtypes have deeper level of annotation, but non of them is found in the dataset
      break

    c.parents %<>% .[names(cbs.per.type)]

    scores.per.type <- lapply(c.parents, function(p) scores[cbs.per.type[[p]], c.subtypes.per.parent[[p]], drop=F])

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

#' Normalize Total Count with Features
#' @description Normalize `cm` matrix using total counts and then column-wise min-max scaling
#'
#' @param max.quantile quantile to be used for max estimation during scaling
#' @return Normalized matrix of the same shape as `cm`
#' @inheritParams getClassificationData
#'
#' @export
normalizeTCWithFeatures <- function(cm, max.quantile=0.95, max.smooth=1e-10, transpose=T) {
  cm %<>% as("dgCMatrix") %>% Matrix::drop0()

  # Total count normalization (i.e. TF-step)
  cm@x <- cm@x / rep(Matrix::colSums(cm), diff(cm@p))
  cm <- Matrix::t(cm)

  # Factors for min-max gene normalization
  max.vals <- split(cm@x, rep(1:(length(cm@p)-1), diff(cm@p)))[paste0(1:ncol(cm))]
  max.vals[is.null(max.vals)] <- c()
  max.vals %<>% sapply(quantile, max.quantile) %>% `+`(max.smooth) # Robust alternative to maximum

  cm@x <- cm@x / rep(max.vals, diff(cm@p)) # fast way to do columnwise normalization for sparse matrices
  cm@x %<>% pmin(1.0)

  if (!transpose)
    return(cm)

  return(Matrix::t(cm))
}

#' Normalize TF-IDF with Features
#' @description Normalize `cm` matrix using total counts, column-wise min-max scaling and then IDF normalization
#' @inheritDotParams normalizeTCWithFeatures max.quantile max.smooth
#'
#' @export
normalizeTfIdfWithFeatures <- function(cm, ...) {
  tf.idf <- normalizeTCWithFeatures(cm, transpose=F, ...)
  # IDF-factors: log(1 + fraction of expressing cells)
  idf.weights <- log(1 + nrow(tf.idf) / (Matrix::colSums(tf.idf > 0) + 1))
  tf.idf@x <- tf.idf@x * rep(idf.weights, diff(tf.idf@p))
  return(Matrix::t(tf.idf))
}

#' Get Marker Score Info
#' @description estimate info, neccessary for scoring of cells by cell types
#'
#' @param clf.data classification data from `getClassificationData`
#' @inheritDotParams getCellTypeScoreInfo aggr
#' @return List of score info for each cell type. See `CellAnnotatoR:::getCellTypeScoreInfo` for more info
#' @export
getMarkerScoreInfo <- function(clf.data, ...) {
  lapply(clf.data$marker.list, getCellTypeScoreInfo, clf.data$cm, ...)
}

#' Merge Scores
#' @param score.name type of score to merge
#' @param score.infos scoring info per cell per cell type
mergeScores <- function(score.name, score.infos, aggr.func=c) {
  names(score.infos[[1]][[score.name]]) %>% setNames(., .) %>% lapply(function(n)
    lapply(score.infos, `[[`, c(score.name, n)) %>% Reduce(aggr.func, .))
}

#' Merge Score Infos
#' @description merge score information from multiple datasets
#'
#' @inheritParams mergeScores
#' @inheritParams plapply
#' @inherit getMarkerScoreInfo return
#' @examples
#'   clf_datas <- lapply(cms, getClassificationData, marker_path)
#'   score_infos <- lapply(clf_datas, getMarkerScoreInfo)
#'   all_score_info <- mergeScoreInfos(score_infos, verbose=T)
#'
#' @export
mergeScoreInfos <- function(score.infos, verbose=F) {
  names(score.infos[[1]]) %>% setNames(., .) %>% plapply(mergeScores, score.infos, verbose=verbose)
}

mergeAnnotationInfos <- function(ann.infos) {
  return(list(
    annotation=mergeScores("annotation", ann.infos),
    scores=mergeScores("scores", ann.infos, aggr.func=rbind),
    annotation.filt=mergeScores("annotation.filt", ann.infos)
  ))
}

#' Get Cell Type Score Info
#' @description estimate info, neccessary for scoring of cells by cell types for the specified cell type
#'
#' @param markers element of marker list. List with fields `expressed` and `not_expressed`
#' @param tf.idf TF-IDF normalized matrix. Can be obtained with `normalizeTfIdfWithFeatures`
#' @param aggr if scores must be aggregated for the whole cell type or returned for each marker separately
#' @return list with score info:\itemize{
#'   \item{scores.raw: scores from positive markers}
#'   \item{mult: score multiplier, estimated with negative markers}
#'   \item{max.positive: maximal expression of positive markers. Used for estimation of negative scores}
#'   \item{scores: final scores. Equal to `scores * score.mult`}
#' }
getCellTypeScoreInfo <- function(markers, tf.idf, aggr=T) {
  expressed.genes <- intersect(markers$expressed, colnames(tf.idf))
  if (length(expressed.genes) == 0) {
    if (aggr) {
      scores <- setNames(rep(0, nrow(tf.idf)), rownames(tf.idf))
    } else {
      scores <- matrix(0, nrow=nrow(tf.idf), ncol=0, dimnames=list(rownames(tf.idf), character()))
    }

    max.positive.expr <- setNames(rep(0, nrow(tf.idf)), rownames(tf.idf))
  } else {
    c.submat <- tf.idf[, expressed.genes, drop=F]
    c.submat.t <- Matrix::t(c.submat)
    scores <- if (aggr) Matrix::colSums(c.submat.t) else c.submat
    max.positive.expr <- apply(c.submat.t, 2, max)
  }

  not.expressed.genes <- intersect(markers$not_expressed, colnames(tf.idf))
  if (length(not.expressed.genes) == 0) {
    score.mult <- setNames(rep(1, nrow(tf.idf)), rownames(tf.idf))
  } else {
    max.negative.expr <- apply(tf.idf[, not.expressed.genes, drop=F], 1, max)
    # max.negative.expr <- sparseRowMax(tf.idf[, not.expressed.genes, drop=F])
    score.mult <- pmax(max.positive.expr - max.negative.expr, 0) / max.positive.expr
    score.mult[is.na(score.mult)] <- 0
  }

  res <- list(scores.raw=scores, mult=score.mult, max.positive=max.positive.expr, scores=(scores * score.mult))
  return(res)
}

normalizeScores <- function(scores, min.val=1e-10) {
  scores[rowSums(scores) < 1e-10,] <- 1
  scores  %<>% `/`(rowSums(.))
  return(scores)
}

#' Return initial scores of each cell type for each cell
#'
#' @param score.info pre-calculated score info from `getMarkerScoreInfo`
#' @param aggr should individual gene scores be aggregated per cell type? If `FALSE`,
#' returns list of data.frames, showing scores of each gene for each cell.
#' Useful for debugging list of markers.
#' @inheritParams getMarkerScoreInfo
#'
#' @return data.frame with rows corresponding to cells and columns corresponding to cell types.
#'   Values are cell type scores, normalized per level of hierarchy
getMarkerScoresPerCellType <- function(clf.data, score.info=NULL, aggr=T) {
  if (is.null(score.info)) {
    score.info <- getMarkerScoreInfo(clf.data, aggr=aggr)
  }

  scores <- lapply(score.info, `[[`, "scores")
  if (!aggr)
    return(lapply(scores, as.matrix) %>% lapply(as.data.frame, optional=T))

  clf.nodes <- classificationTreeToDf(clf.data$classification.tree)
  scores %<>% as.data.frame(optional=T)

  for (nodes in split(clf.nodes$Node, clf.nodes$PathLen)) {
    scores[, nodes]  %<>% normalizeScores()
  }

  return(scores)
}
