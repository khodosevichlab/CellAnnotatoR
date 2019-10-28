unifyGeneIds <- function(cm, data.gene.id.type, marker.gene.id.type, db=NULL, verbose=F) {
  if (marker.gene.id.type == data.gene.id.type)
    return(list(cm=cm, gene.table=data.frame(data=rownames(cm), marker=rownames(cm), stringsAsFactors=F)))

  if (is.null(db))
    stop("db must be provided for marker conversion")

  new.ids <- garnett:::convert_gene_ids(rownames(cm), db, data.gene.id.type, marker.gene.id.type)
  new.ids %<>% .[!is.na(.)] %>% .[!duplicated(.)]

  n.lost <- nrow(cm) - length(new.ids)
  frac.lost <- n.lost / nrow(cm)
  if (frac.lost > 0.7)
    warning("When converting to", marker.gene.id.type, "IDs ", round(frac.lost * 100), "% were lost.",
            "Did you specify the correct gene ID types and the correct db?")

  if (verbose)
    message("After converting CDS to", marker.gene.id.type, "IDs,", nrow(cm) - length(new.ids),
            " (", round(frac.lost * 100), "%) IDs were lost")

  cm <- cm[names(new.ids), ] %>% magrittr::set_rownames(new.ids)
  gene.table <- data.frame(data=names(new.ids), marker=as.vector(new.ids), stringsAsFactors=F)
  return(list(cm=cm, gene.table=gene.table))
}

#' @export
parseMarkerFile <- function(path) {
  marker.list <- path %>% readChar(file.info(.)$size) %>% paste0("\n") %>% garnett:::parse_input() %>% as.list()
  marker.list$name_order <- NULL
  return(marker.list)
}

createClassificationTree <- function(marker.list) {
  parents <- sapply(marker.list, function(pl) if(length(pl@parenttype) == 0) "root" else pl@parenttype)
  tree <- c(parents, names(marker.list)) %>% matrix(ncol=2) %>% t() %>% igraph::make_directed_graph()
  return(tree)
}

#' Get cell type information from the marker file
#'
#' @param cm gene count matrix. May be in raw, TC-normalized or tf-idf-normalized format (in case of tf-idf, `prenormalized` must be set to `T`)
#' @param markers path to the file with marker genes or parsed marker list from `parseMarkerFile` function
#' @param prenormalized is `cm` in tf-idf-normalized format? Default: FALSE.
#' @export
getClassificationData <- function(cm, markers, prenormalized=F, data.gene.id.type="SYMBOL", marker.gene.id.type="SYMBOL", db=NULL, verbose=F) {
  if (!prenormalized) {
    cm <- normalizeTfIdfWithFeatures(cm)
  }

  gi <- unifyGeneIds(cm, data.gene.id.type=data.gene.id.type, marker.gene.id.type=marker.gene.id.type, db=db, verbose=verbose)

  if (is.character(markers)) {
    markers <- parseMarkerFile(markers)
  } else if (!is.list(markers)) {
     stop("Unknown format of markers")
  }

  classification.tree <- createClassificationTree(markers)

  res <- list(cm=gi$cm, classification.tree=classification.tree, gene.table=gi$gene.table, marker.list=markers)
  return(res)
}
