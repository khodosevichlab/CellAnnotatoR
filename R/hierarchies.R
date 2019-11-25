#' Merge Annotation to Level
#' @param annotation annotation for high-resolution. Cell type names must correspond to nodes in the `classification.tree`
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

#' Merge Annotation By Levels
#' @description merge provided annotation using `classification.tree` to get annotations for different levels of hierarchy
#'
#' @inheritParams classificationTreeToDf
#' @inheritParams mergeAnnotationToLevel
#' @return list of annotations where each element correspond to some hierarchy level
#'
#' @examples
#'   ann_by_level <- mergeAnnotationByLevels(annotation, clf_data$classification.tree)
#'
#' @export
mergeAnnotationByLevels <- function(annotation, classification.tree) {
  max.level <- classificationTreeToDf(classification.tree)$PathLen %>% max()
  anns <- 1:max.level %>% setNames(., paste0("l", .)) %>%
    lapply(mergeAnnotationToLevel, annotation, classification.tree)
  return(anns)
}

simplifyHierarchy <- function(branch) {
  if (!is.list(branch))
    return(branch)

  for (i in 1:length(branch)) {
    while (is.list(branch[[i]]) && length(branch[[i]]) == 1) {
      branch[[i]] <- branch[[i]][[1]]
    }
  }

  if (length(branch) == 1)
    return(branch[[1]])

  return(lapply(branch, simplifyHierarchy))
}

appendHierarchyBranch <- function(branch, parent.name) {
  branch %<>% simplifyHierarchy()
  is.leaf <- (sapply(branch, is.character) & (sapply(branch, length) == 1))

  current.nodes <- c(unlist(branch[is.leaf], use.names=F), names(branch[!is.leaf])) %>%
    setNames(., .) %>% lapply(function(x) list(expressed=c(), not_expressed=c(), parent=parent.name))

  sub.branches <- names(branch)[!is.leaf] %>%
    lapply(function(n) appendHierarchyBranch(branch[[n]], n)) %>% unlist(recursive=F)
  return(c(current.nodes, sub.branches))
}

#' Hierarchy to Classification Tree
#' @description Convert list with cell type hierarchy to classification tree in igraph format
#' @param hierarchy list of cell types, where each element represent a cell type. If cell type has some subtypes it's represented as another list, otherwise it's just a string
#' @inherit createClassificationTree return
#' @examples
#'   hierarchy <- list(
#'     Alveolar=c("AT1 Cell", "AT2 Cell"),
#'     B=c("B Cell", "Ig-producing B cell"),
#'     NK_T=list(`T`=c("Dividing T cells", "T Cell_Cd8b1 high", "Nuocyte"), "NK Cell"),
#'     "Endothelial"
#'   )
#'
#'   clf_tree <- hierarchyToClassificationTree(hierarchy)
#'
#' @export
hierarchyToClassificationTree <- function(hierarchy) {appendHierarchyBranch(hierarchy, "root") %>% createClassificationTree()}

splitClusteringDf <- function(df) {
  if (ncol(df) == 1)
    return(df[[1]])

  return(split(df, df[,1]) %>% lapply(function(x) splitClusteringDf(x[,2:ncol(x)])))
}

#' Derive Hierarchy
#' @description derive hierarchy from the data using hclust
#' @param feature.matrix matrix where rows represent cells and columns represent either genes or some embedded space (e.g. PCA)
#' @param annotation vector with cell type label per cell
#' @param dist.method method for pairwise distance estimation. Either "cor" for correlation distance or any method supported by `dist` function
#' @param max.depth maximal depth of the hierarchy
#' @return list with cell type hierarchy
#' @examples
#'   hierarchy <- deriveHierarchy(pca_mtx, annotation, max.depth=3)
#'   clf_tree <- hierarchyToClassificationTree(hierarchy)
#'   plotTypeHierarchy(clf_tree)
#'
#' @export
deriveHierarchy <- function(feature.matrix, annotation, dist.method="cor", max.depth=2) {
  feature.matrix %<>% as("dgCMatrix") %>% conos:::collapseCellsByType(groups=annotation, min.cell.count=0)

  if (dist.method == "cor") {
    c.dist <- (1 - cor(Matrix::t(feature.matrix))) %>% as.dist()
  } else {
    c.dist <- dist(feature.matrix, method=dist.method)
  }

  clusts <- hclust(c.dist)
  heights <- quantile(clusts$height, seq(0, 1, length.out=max.depth + 1) %>% .[2:(length(.)-1)]) %>% rev()
  clusts %<>% cutree(h=heights) %>% as.matrix()
  for (i in 1:ncol(clusts)) {
    clusts[,i] %<>% paste0("l", i, "_", .)
  }

  clusts %<>% cbind(rownames(.)) %>% `colnames<-`(paste0("l", 1:ncol(.))) %>% tibble::as_tibble()

  return(clusts %>% splitClusteringDf() %>% simplifyHierarchy())
}
