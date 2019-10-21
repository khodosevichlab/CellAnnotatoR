plotAssignmentScores <- function(parent.node, scores, classification.tree, embedding, ...) {
  classificationTreeToDf(classification.tree) %>%
    dplyr::filter(Parent == parent.node) %>% .$Node %>%
    lapply(function(n) conos:::embeddingPlot(embedding, colors=setNames(scores[, n], rownames(scores)),
                                             title=n, color.range=c(0, 1), ...))
}

plotGarnettAssignemnts <- function(embedding, assignment, plot.ambig=T, show.legend=T, mark.groups=F, ...) {
  lapply(colnames(assignment), function(n) {
    groups <- setNames(assignment[[n]], rownames(assignment));
    if (!plot.ambig) {
      groups <- droplevels(groups[!(groups %in% c("None", "Ambiguous"))])
    }

    conos:::embeddingPlot(
      embedding, groups=groups, show.legend=show.legend, mark.groups=mark.groups, title=n,
      plot.theme=ggplot2::theme(legend.background=ggplot2::element_rect(fill=ggplot2::alpha("white", 0.2))), ...)
  })
}

#' @export
plotTypeMarkers <- function(parse.list, cell.type, embedding, count.matrix, show.legend=T, ...) {
  plot.func <- function(gene, title) {
    conos:::embeddingPlot(embedding, colors=count.matrix[,gene], show.legend=show.legend, title=title, ...)
  }

  res <- parse.list[[cell.type]]@expressed %>% intersect(colnames(count.matrix)) %>%
    lapply(function(g) plot.func(g, paste0(g, "+"))) %>%
    c(parse.list[[cell.type]]@not_expressed %>% intersect(colnames(count.matrix)) %>%
        lapply(function(g) plot.func(g, paste0(g, "-"))))

  return(res)
}

#' @export
plotSubtypeMarkers <- function(parent.type, embedding, count.matrix, parse.list, classification.tree, show.legend=F, drop.missing=T, ...) {
  markers <- list()

  for (type in getAllSubtypes(parent.type, classification.tree)) {
    for (g in parse.list[[type]]@expressed) {
      markers[[g]] %<>% c(paste0(type, "+"))
    }

    for (g in parse.list[[type]]@not_expressed) {
      markers[[g]] %<>% c(paste0(type, "-"))
    }
  }

  if (drop.missing) {
    markers <- markers[names(markers) %in% colnames(count.matrix)]
  }

  titles <- mapply(function(n,ts) paste0(n, ": ", ts), names(markers), lapply(markers, paste, collapse=", "))
  mapply(function(gene, title) {
    expr.vec <- if (gene %in% colnames(count.matrix)) count.matrix[,gene] else NA
    conos:::embeddingPlot(embedding, colors=expr.vec, show.legend=show.legend, title=title, ...)
  },names(markers), titles, SIMPLIFY=F
  )
}

#' @export
plotAnnotationByLevels <- function(embedding, annotation.by.level, clusters=NULL, ...) {
  res <- lapply(1:length(annotation.by.level),
                function(i) conos:::embeddingPlot(embedding, groups=annotation.by.level[[i]], title=paste("Level", i), ...))

  if (!is.null(clusters)) {
    res[[length(res) + 1]] <- conos:::embeddingPlot(embedding, groups=clusters, title="Clustering", ...)
  }

  return(res)
}

plotConfidenceByLevels <- function(embedding, annotation.by.level, scores, show.legend=T, ...) {
  if (class(scores) == "data.frame") {
    scores <- lapply(annotation.by.level, function(x) scores)
  }

  conf.per.level <- lapply(names(scores), function(n)
    mapply(function(i,j) scores[[n]][i, j], 1:nrow(scores[[n]]),
           match(annotation.by.level[[n]], colnames(scores[[n]]))) %>%
      setNames(rownames(scores[[n]])))

  res <- lapply(1:length(conf.per.level), function(i)
    conos:::embeddingPlot(embedding, colors=conf.per.level[[i]], title=paste("Level", i),
                          legend.title="Confidence", show.legend=T, color.range=c(0, 1), ...))

  return(res)
}

#' @export
plotTypeHierarchy <- function(classification.tree, layout="slanted", xlims=NULL, font.size=3, ...) {
  if (!requireNamespace("ggtree", quietly=T))
    stop("You need to install package 'ggtree' to be able to plot hierarchical tree. ",
         "`Try devtools::install_github('YuLab-SMU/ggtree')`")

  c.df <- classificationTreeToDf(classification.tree)
  cg <- c.df %$% data.frame(parent=Parent, node=Node) %>% ape::as.phylo()
  cg$edge.length <- rep(1, nrow(cg$edge))

  if (is.null(xlims)) {
    xlims <- c(0, max(c.df$PathLen) + 0.5)
  }

  ggtree::ggtree(cg, layout = layout, ...) +
    ggtree::geom_rootpoint() +
    ggtree::geom_label2(ggplot2::aes(label=c(cg$tip.label, cg$node.label)[node]), size=font.size) +
    ggplot2::xlim(xlims)
}
