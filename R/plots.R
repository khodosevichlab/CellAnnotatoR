#' @export
plotAssignmentScores <- function(embedding, scores, classification.tree, parent.node="root", ...) {
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
plotTypeMarkers <- function(embedding, count.matrix, cell.type, marker.list, show.legend=T, ...) {
  plot.func <- function(gene, title) {
    conos:::embeddingPlot(embedding, colors=count.matrix[,gene], show.legend=show.legend, title=title, ...)
  }

  res <- marker.list[[cell.type]]$expressed %>% intersect(colnames(count.matrix)) %>%
    lapply(function(g) plot.func(g, paste0(g, "+"))) %>%
    c(marker.list[[cell.type]]$not_expressed %>% intersect(colnames(count.matrix)) %>%
        lapply(function(g) plot.func(g, paste0(g, "-"))))

  return(res)
}

#' @export
plotSubtypeMarkers <- function(embedding, count.matrix, parent.type="root", clf.data=NULL, clf.tree=NULL, marker.list=NULL,
                               show.legend=F, max.depth=NULL, drop.missing=T, build.panel=T, n.col=NULL, n.row=NULL, ...) {
  if (is.null(clf.data) && (is.null(clf.tree) || is.null(marker.list)))
    stop("Either clf.data or both clf.tree and marker.list must be provided")

  if (!is.null(clf.data)) {
    clf.tree <- clf.data$classification.tree
    marker.list <- clf.data$marker.list
  }

  markers <- list()

  for (type in getAllSubtypes(parent.type, clf.tree, max.depth=max.depth)) {
    for (g in marker.list[[type]]$expressed) {
      markers[[g]] %<>% c(paste0(type, "(+)"))
    }

    for (g in marker.list[[type]]$not_expressed) {
      markers[[g]] %<>% c(paste0(type, "(-)"))
    }
  }

  if (drop.missing) {
    markers <- markers[names(markers) %in% colnames(count.matrix)]
  }

  titles <- mapply(function(n,ts) paste0(n, ": ", ts), names(markers), lapply(markers, paste, collapse=", "))
  plots <- mapply(function(gene, title) {
    expr.vec <- if (gene %in% colnames(count.matrix)) count.matrix[,gene] else NA
    conos:::embeddingPlot(embedding, colors=expr.vec, show.legend=show.legend, title=title, ...)
  },names(markers), titles, SIMPLIFY=F
  )

  if (!build.panel)
    return(plots)

  return(cowplot::plot_grid(plotlist=plots, ncol=n.col, nrow=n.row))
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
    getAnnotationConfidence(annotation.by.level[[n]], scores[[n]]))

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

plotUncertaintyPerCell <- function(embedding, uncertainty.info, palette=colorRampPalette(c("gray", "#ffeda0", "#fec44f", "#f03b20")), alpha=0.3, ...) {
  names(uncertainty.info) %>% setNames(., .) %>% lapply(function(n)
    conos::embeddingPlot(embedding, colors=uncertainty.info[[n]], alpha=alpha, title=n, palette=palette, ...))
}

plotUncertaintyPerClust <- function(uncertainty.per.clust, annotation=NULL, clusters=NULL, ann.per.clust=NULL, threshold=0.5, text.angle=45, title=NULL) {
  if (!is.null(ann.per.clust) && (is.null(annotation) != is.null(clusters)))
    stop("Either both or none of annotation and clusters must be provided")

  p.df <- tibble::enframe(uncertainty.per.clust, name="Cluster", value="Uncertainty")
  p.aes <- ggplot2::aes()
  if (!is.null(annotation)) {
    ann.per.clust <- getAnnotationPerCluster(annotation, clusters)
  }

  if (!is.null(ann.per.clust)) {
    p.aes <- ggplot2::aes(fill=Type)
    p.df %<>% dplyr::mutate(Type=ann.per.clust[Cluster], Cluster=factor(Cluster, levels=names(sort(ann.per.clust))))
  }

  gg <- ggplot2::ggplot(p.df, ggplot2::aes(x=Cluster, y=Uncertainty)) +
    ggplot2::geom_bar(p.aes, stat="identity") +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=text.angle, hjust=1)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=threshold)) + ggplot2::ylim(0, 1)

  if (!is.null(title)) {
    gg <- gg + ggplot2::ggtitle(title)
  }

  return(gg)
}
