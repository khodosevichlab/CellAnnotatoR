#' @importFrom magrittr %<>% %$% %>%

diffuseSubgraph <- function(graph, annotation, subtype, scores, fading=10, fading.const=0.5,
                            clusters=NULL, verbose=FALSE) {
  cbs <- names(annotation)[annotation == subtype]
  if (length(cbs) == 0)
    return(NULL)

  subgraph <- igraph::induced_subgraph(graph, cbs)

  scores <- as.matrix(scores[cbs, ])
  scores[rowSums(scores) < 1e-5, ] <- 1
  scores %<>% `/`(rowSums(.))

  edges <- igraph::as_edgelist(subgraph)
  edge.weights <- igraph::edge.attributes(subgraph)$weight
  res <- conos:::smooth_count_matrix(edges, edge.weights, scores, max_n_iters=1000,
                                     diffusion_fading=fading, diffusion_fading_const=fading.const,
                                     verbose=verbose, normalize=T)

  if (is.null(clusters)) {
    ann.update <- colnames(res)[apply(res, 1, which.max)] %>% setNames(rownames(res))
  } else {
    clusters <- droplevels(clusters[rownames(res)])
    ann.update <- clusters %>% split(names(.), .) %>%
      sapply(function(cbs) colSums(res[cbs, ,drop=F]) %>% which.max() %>% names()) %>%
      .[clusters] %>% setNames(names(clusters))
  }

  annotation[names(ann.update)] <- ann.update
  return(list(annotation=annotation, scores=res))
}

#' @export
assignCellsByScores <- function(graph, scores, classification.tree, clusters=NULL, ...) {
  if (!is.null(clusters)) {
    clusters <- as.factor(clusters)
  }

  paths <- classificationTreeToDf(classification.tree) %$% split(., PathLen) %>%
    lapply(function(df) split(df$Node, df$Parent)) %>% .[order(names(.))]

  c.ann <- rep("root", length(igraph::V(graph))) %>% setNames(names(igraph::V(graph)))
  ann.by.level <- list()
  scores.by.level <- list()
  scores.posterior <- scores

  possible.ann.levels <- c()
  for (pl in 1:length(paths)) {
    c.path <- paths[[pl]]
    possible.ann.levels %<>% c(unlist(c.path)) %>% unique()

    for (p in names(c.path)) {
      res <- diffuseSubgraph(graph, c.ann, p, scores[, c.path[[p]]], clusters=clusters, ...)
      if (is.null(res))
        next

      cur.type.cbs <- names(c.ann)[c.ann == p]
      scores.posterior[cur.type.cbs, colnames(res$scores)] <- res$scores[cur.type.cbs,]
      c.ann <- res$annotation
    }

    level.name <- paste0("l", pl)
    ann.by.level[[level.name]] <- c.ann
    scores.by.level[[level.name]] <- scores.posterior[,possible.ann.levels]
  }

  return(list(annotation=ann.by.level, scores=scores.by.level))
}


## Score assignment

tfIdfFeatureNorm <- function(cm, max.quantile=0.95, max.smooth=1e-10) {
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

getCellTypeScore <- function(markers, tf.idf, aggr=T, do.multiply=T) {
  c.submat <- tf.idf[, intersect(markers@expressed, colnames(tf.idf)), drop=F]
  c.submat.t <- Matrix::t(c.submat)
  scores <- if (aggr) Matrix::colSums(c.submat.t) else c.submat

  not.expressed.genes <- intersect(markers@not_expressed, colnames(tf.idf))
  if (length(not.expressed.genes) == 0)
    return(scores)

  max.positive.expr <- apply(c.submat.t, 2, max)
  max.negative.expr <- apply(tf.idf[, not.expressed.genes, drop=F], 1, max)

  score.mult <- pmax(max.positive.expr - max.negative.expr, 0) / max.positive.expr
  score.mult[is.na(score.mult)] <- 0

  if (!do.multiply)
    return(list(scores=scores, mult=score.mult))

  return(scores * score.mult)
}

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

#' @export
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

classificationTreeToDf <- function(classification.tree) {
  igraph::V(classification.tree)$name %>% .[. != "root"] %>%
    lapply(function(cn) {
      p <- igraph::all_simple_paths(classification.tree, "root", to=cn, mode="out")[[1]]$name
      tibble::tibble(Parent=tail(p, 2)[1], Node=tail(p, 1), PathLen=length(p) - 1)
    }) %>%  Reduce(rbind, .)
}

getAllSubtypes <- function(parent.type, classification.tree) {
  igraph::dfs(classification.tree, parent.type, neimode="out", unreachable=F)$order %>%
    names() %>% .[!is.na(.)] %>% .[. != parent.type]
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
