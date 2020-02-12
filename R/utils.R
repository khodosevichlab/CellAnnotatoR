#' Parallel Lapply
#' @description parallel, optionally verbose lapply
#' @param n.cores number of cores to use
#' @param verbose show progress bar
plapply <- function(..., n.cores=1, verbose=F) {
  if (verbose)
    return(pbapply::pblapply(..., cl=n.cores))

  if (n.cores == 1)
    return(lapply(...))

  return(parallel::mclapply(..., mc.cores=n.cores))
}

sparseColMax <- function(mtx) {
  max.vals <- rep(0, ncol(mtx))
  facs <- split(mtx@x, rep(1:(length(mtx@p)-1), diff(mtx@p)))
  facs <- sapply(facs, max)
  max.vals[as.integer(names(facs))] <- unlist(facs)
  return(max.vals)
}

sparseRowMax <- function(mtx) {
  max.vals <- rep(0, nrow(mtx))
  facs <- split(mtx@x, mtx@i + 1)
  facs <- sapply(facs, max)
  max.vals[as.integer(names(facs))] <- as.vector(facs)
  return(max.vals)
}

#' Get Annotation Per Parent
#' @description for each cell type get annotation of its subtypes on the next hierarchy level
#'
#' @inheritParams classificationTreeToDf
#' @inheritParams mergeAnnotationToLevel
#' @return list of sub-annotations named by parent types
#' @examples
#'   ann_by_parent <- getAnnotationPerParent(clf_data$classification.tree, annotation)
#' @export
getAnnotationPerParent <- function(classification.tree, annotation) {
  classificationTreeToDf(classification.tree) %>% split(.$Parent) %>% lapply(function(df)
    mergeAnnotationToLevel(df$PathLen[1], annotation, classification.tree) %>% .[. %in% df$Node])
}

#' Classification Tree to DataFrame
#' @param classification.tree cell type hierarchy represented by graph. Part of `clf_data` object from `getClassificationData`
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

rateMetricsPerType <- function(ann.res, ann.reference) {
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
