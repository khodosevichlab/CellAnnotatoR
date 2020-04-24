#' Arrange Plots
#' @param plot.list list of plots
#' @param build.panel combine individual plots to a single panel
#' @param n.col number of columns in the panel
#' @param n.row number of rows in the panel
#' @param title.size font size for title of individual plots. Used only if build.panel==T
#' @return panel of plots if `build.panel==T` or `plot.list` otherwise
arrangePlots <- function(plot.list, build.panel, n.col=NULL, n.row=NULL, title.size=10) {
  if (build.panel) {
    p.theme <- ggplot2::theme(plot.title=ggplot2::element_text(size=10, vjust=-1),
                              plot.margin=ggplot2::margin())
    plot.list %<>% lapply(`+`, p.theme)
    return(cowplot::plot_grid(plotlist=plot.list, ncol=n.col, nrow=n.row))
  }

  return(plot.list)
}

#' Plot Assignment Scores
#' @description plot assignment scores on a separate scatterplot for each cell type
#' @param scores assignment scores. Can be obtained with `assignCellsByScores`
#' @param parent.node cell type, which subtypes' scores must be plotted
#' @inheritParams conos::embeddingPlot
#' @inheritParams classificationTreeToDf
#' @inheritParams arrangePlots
#' @inheritDotParams conos::embeddingPlot
#' @inherit arrangePlots return
#' @examples
#'   clf_data <- getClassificationData(cm, marker_path)
#'   ann_by_level <- assignCellsByScores(graph, clf_data)
#'   plotAssignmentScores(t_sne, ann_by_level$scores$l1, clf_data$classification.tree)
#' @export
plotAssignmentScores <- function(embedding, scores, classification.tree, parent.node="root", build.panel=T, n.col=NULL, n.row=NULL, title.size=12, ...) {
  classificationTreeToDf(classification.tree) %>%
    dplyr::filter(Parent == parent.node) %>% .$Node %>%
    lapply(function(n) conos::embeddingPlot(embedding, colors=setNames(scores[, n], rownames(scores)),
                                             title=n, color.range=c(0, 1), ...)) %>%
    arrangePlots(build.panel=build.panel, n.row=n.row, n.col=n.col, title.size=title.size)
}

#' Plot Type Markers
#' @description plot markers for the specified `cell.type`
#' @param count.matrix gene expression matrix with genes by columns and cells by rows
#' @param cell.type cell type for which the markers must be plotted
#' @param marker.list list of markers per cell type. Can be obtained with `parseMarkerFile`
#' @inheritParams conos::embeddingPlot
#' @inheritParams arrangePlots
#' @inheritDotParams conos::embeddingPlot
#' @export
plotTypeMarkers <- function(embedding, count.matrix, cell.type, marker.list, show.legend=T, build.panel=T, n.col=NULL, n.row=NULL, title.size=10, ...) {
  plot.func <- function(gene, title) {
    conos::embeddingPlot(embedding, colors=count.matrix[,gene], show.legend=show.legend, title=title, ...)
  }

  res <- marker.list[[cell.type]]$expressed %>% intersect(colnames(count.matrix)) %>%
    lapply(function(g) plot.func(g, paste0(g, "+"))) %>%
    c(marker.list[[cell.type]]$not_expressed %>% intersect(colnames(count.matrix)) %>%
        lapply(function(g) plot.func(g, paste0(g, "-")))) %>%
    arrangePlots(build.panel=build.panel, n.row=n.row, n.col=n.col, title.size=title.size)

  return(res)
}

#' Extract Markers From Subtypes
#' @param parent.type cell type for which the markers of the subtypes must be plotted
extractMarkersFromSubtypes <- function(parent.type="root", clf.data=NULL, clf.tree=NULL, marker.list=NULL, count.matrix=NULL,
                                       max.depth=NULL, drop.missing=T, marker.type=c("expressed", "not_expressed")) {
  if (length(setdiff(marker.type, c("expressed", "not_expressed"))) > 0)
    stop("Unknown marker.type")

  if (is.null(clf.data) && is.null(marker.list))
    stop("Either clf.data or marker.list must be provided")

  if (!is.null(clf.data)) {
    clf.tree <- clf.data$classification.tree
    marker.list <- clf.data$marker.list
  }

  if (is.null(clf.tree)) {
    clf.tree <- createClassificationTree(marker.list)
  }

  markers <- list()

  for (type in sort(getAllSubtypes(parent.type, clf.tree, max.depth=max.depth))) {
    if ("expressed" %in% marker.type) {
      for (g in marker.list[[type]]$expressed) {
        markers[[g]] %<>% c(paste0(type, "(+)"))
      }
    }

    if ("not_expressed" %in% marker.type) {
      for (g in marker.list[[type]]$not_expressed) {
        markers[[g]] %<>% c(paste0(type, "(-)"))
      }
    }
  }

  if (drop.missing && !is.null(count.matrix)) {
    markers <- markers[names(markers) %in% colnames(count.matrix)]
  }

  return(markers)
}

#' Plot Subtype Markers
#' @description plot markers for all subtypes of the specified `parent.type`
#' @param count.matrix gene expression matrix with genes by columns and cells by rows

#' @param marker.list list of markers per cell type. Can be obtained with `parseMarkerFile`
#' @inheritParams conos::embeddingPlot
#' @inheritParams arrangePlots
#' @inheritParams plotTypeMarkers
#' @inheritParams extractMarkersFromSubtypes
#' @inheritDotParams conos::embeddingPlot
#' @export
plotSubtypeMarkers <- function(embedding, count.matrix, parent.type="root", clf.data=NULL, clf.tree=NULL, marker.list=NULL,
                               show.legend=F, max.depth=NULL, build.panel=T, n.col=NULL, n.row=NULL, marker.type=c("expressed", "not_expressed"),
                               title.size=10, ...) {
  markers <- extractMarkersFromSubtypes(parent.type=parent.type, clf.data=clf.data, clf.tree=clf.tree, marker.list=marker.list,
                                        count.matrix=count.matrix, max.depth=max.depth, drop.missing=T, marker.type=marker.type)

  titles <- mapply(function(n,ts) paste0(n, ": ", ts), names(markers), lapply(markers, paste, collapse=", "))
  plots <- mapply(function(gene, title) {
    expr.vec <- if (gene %in% colnames(count.matrix)) count.matrix[,gene] else NA
    conos::embeddingPlot(embedding, colors=expr.vec, show.legend=show.legend, title=title, ...)
  },names(markers), titles, SIMPLIFY=F
  )

  return(arrangePlots(plots, build.panel=build.panel, n.row=n.row, n.col=n.col, title.size=title.size))
}

#' @export
plotAnnotationByLevels <- function(embedding, annotation.by.level, clusters=NULL, build.panel=T, n.col=NULL, n.row=NULL, title.size=12, ...) {
  res <- lapply(1:length(annotation.by.level),
                function(i) conos::embeddingPlot(embedding, groups=annotation.by.level[[i]], title=paste("Level", i), ...))

  if (!is.null(clusters)) {
    res[[length(res) + 1]] <- conos::embeddingPlot(embedding, groups=clusters, title="Clustering", ...)
  }

  return(arrangePlots(res, build.panel=build.panel, n.row=n.row, n.col=n.col, title.size=title.size))
}

plotConfidenceByLevels <- function(embedding, annotation.by.level, scores, show.legend=T, ...) {
  if (class(scores) == "data.frame") {
    scores <- lapply(annotation.by.level, function(x) scores)
  }

  conf.per.level <- lapply(names(scores), function(n)
    getAnnotationConfidence(annotation.by.level[[n]], scores[[n]]))

  res <- lapply(1:length(conf.per.level), function(i)
    conos::embeddingPlot(embedding, colors=conf.per.level[[i]], title=paste("Level", i),
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

#' Uncertainty scatterplots per cell
#' @param build.panel join plots to single panel
#' @inheritDotParams conos::embeddingPlot
#' @export
plotUncertaintyPerCell <- function(embedding, uncertainty.info, palette=colorRampPalette(c("gray", "#ffeda0", "#fec44f", "#f03b20")), alpha=0.3,
                                   build.panel=T, n.col=length(uncertainty.info), n.row=NULL, title.size=12, ...) {
  names(uncertainty.info) %>% setNames(., .) %>% lapply(function(n)
    conos::embeddingPlot(embedding, colors=uncertainty.info[[n]], alpha=alpha, title=n, palette=palette, ...)) %>%
    arrangePlots(build.panel=build.panel, n.col=n.col, n.row=n.row, title.size=title.size)
}

#' Plot One Uncertainty per Cluster
#' @param text.angle angle of x-axis labels
#' @inheritDotParams sccore:::styleEmbeddingPlot
plotOneUncertaintyPerClust <- function(uncertainty.per.clust, clusters, annotation=NULL, ann.per.clust=NULL, threshold=0.5, text.angle=45, ...) {
  p.df <- tibble::enframe(uncertainty.per.clust, name="Cluster", value="Uncertainty")
  p.aes <- ggplot2::aes()
  if (is.null(ann.per.clust) && !is.null(annotation)) {
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

  gg %<>% sccore:::styleEmbeddingPlot(show.ticks=T, show.labels=T, ...)

  return(gg + ggplot2::labs(x="", y="Uncertainty"))
}

#' Uncertainty barplots per cluster
#' @inheritParams plotOneUncertaintyPerClust
#' @inheritDotParams plotOneUncertaintyPerClust text.angle
#' @export
plotUncertaintyPerClust <- function(uncertainty.per.clust, clusters, annotation=NULL, ann.per.clust=NULL,
                                    thresholds=c(coverage=0.5, negative=0.5, positive=0.75), build.panel=T,
                                    n.col=1, n.row=NULL, adjust.legend=T, rel.legend.width=0.25, title.size=12, ...) {
  ggs <- names(uncertainty.per.clust) %>% setNames(., .) %>% lapply(function(n)
    plotOneUncertaintyPerClust(uncertainty.per.clust[[n]], annotation=annotation, clusters=clusters,
                               ann.per.clust=ann.per.clust, threshold=thresholds[[n]], title=n, ...))

  if (!build.panel || !adjust.legend || is.null(annotation))
    return(arrangePlots(ggs, build.panel=build.panel, n.col=n.col, n.row=n.row, title.size=title.size))

  if (!requireNamespace("ggpubr", quietly=T))
    stop("You need to install package 'ggpubr' to be able to adjust legend")

  ggl <- ggpubr::get_legend(ggs[[1]])
  ggs %<>% lapply(`+`, ggplot2::theme(legend.position="none", plot.margin=ggplot2::margin()))
  ggr <- arrangePlots(ggs, build.panel=build.panel, n.col=n.col, n.row=n.row, title.size=title.size) %>%
    cowplot::plot_grid(ggl, rel_widths=c(1 - rel.legend.width, rel.legend.width))

  return(ggr)
}

#' @export
plotAssignmentConfusion <- function(scores, annotation=NULL, clusters=annotation, ann.per.clust=NULL) {
  if (!is.null(annotation)) {
    ann.per.clust <- getAnnotationPerCluster(annotation, clusters)
  } else if (is.null(ann.per.clust) || is.null(clusters)) {
    stop("Either annotation or both clusters and ann.per.clust must be provided")
  }

  scores %>% dplyr::group_by(Cluster=clusters[rownames(.)]) %>%
    dplyr::summarise_all(dplyr::funs(mean)) %>%
    data.frame(row.names=.$Cluster, check.names=F) %>% .[,2:ncol(.)] %>%
    .[names(sort(ann.per.clust)), sort(colnames(.))] %>%
    pheatmap::pheatmap(cluster_rows=F, cluster_cols=F, annotation_row=data.frame(ann.per.clust))
}

#' Plot expression violin map for given marker list
#' @inheritParams plotExpressionViolinMap
#' @inheritParams plotSubtypeMarkers
#' @export
plotMarkerListViolinMap <- function(count.matrix, annotation, parent.type="root", clf.data=NULL, clf.tree=NULL, marker.list=NULL,
                                    max.depth=NULL, marker.type="expressed", text.angle=45, gene.order=NULL, ...) {
  extractMarkersFromSubtypes(parent.type=parent.type, clf.data=clf.data, clf.tree=clf.tree, marker.list=marker.list,
                             count.matrix=count.matrix, max.depth=max.depth, drop.missing=T, marker.type=marker.type) %>%
    names() %>% plotExpressionViolinMap(count.matrix, annotation, text.angle=text.angle, gene.order=gene.order)
}

#' @export
plotExpressionViolinMap <- function(markers, count.matrix, annotation, text.angle=45, gene.order=NULL) {
  p.df <- lapply(markers, function(g)
    data.frame(Expr=count.matrix[names(annotation),g], Type=annotation, Gene=g)) %>%
    Reduce(rbind, .)

  if (is.logical(gene.order) && gene.order) {
    gene.order <- unique(markers)
  } else {
    gene.order <- NULL
  }

  if (!is.null(gene.order)) {
    p.df %<>% dplyr::mutate(Gene=factor(as.character(Gene), levels=gene.order))
  }

  ggplot2::ggplot(p.df) +
    ggplot2::geom_violin(ggplot2::aes(x=Type, y=Expr, fill=Type), scale="width") +
    ggplot2::facet_grid(rows=ggplot2::vars(Gene)) +
    ggplot2::theme(
      strip.text.y=ggplot2::element_text(angle=0),
      panel.spacing=ggplot2::unit(0, "pt"), panel.grid=ggplot2::element_blank(),
      axis.text.y=ggplot2::element_blank(), axis.ticks=ggplot2::element_blank(),
      axis.text.x=ggplot2::element_text(angle=text.angle, hjust=1), legend.position="none")
}

#' Plot gene expression on cell embedding
#'
#' @param genes vector of gene names
#' @param embedding cell embedding
#' @param cm count matrix with genes as columns
#' @inheritParams arrangePlots
#' @inheritDotParams conos::embeddingPlot
#' @export
plotGeneExpression <- function(genes, embedding, cm, build.panel=T, n.col=NULL, n.row=NULL, title.size=12, ...) {
  res <- genes %>% setNames(., .) %>%
    lapply(function(g) conos::embeddingPlot(embedding, colors=cm[, g], title=g, ...))

  if (length(res) == 1)
    return(res[[1]])

  return(arrangePlots(res, build.panel=build.panel, n.col=n.col, n.row=n.row, title.size=title.size))
}
