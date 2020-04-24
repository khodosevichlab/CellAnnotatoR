convertGeneIds <- function (genes, db, source.type, target.type) {
  suppressMessages(AnnotationDbi::mapIds(db, keys=genes, column=target.type, source.type))
}

unifyGeneIds <- function(cm, data.gene.id.type, marker.gene.id.type, db=NULL, verbose=F) {
  if (marker.gene.id.type == data.gene.id.type)
    return(list(cm=cm, gene.table=data.frame(data=rownames(cm), marker=rownames(cm), stringsAsFactors=F)))

  if (is.null(db))
    stop("db must be provided for marker conversion")

  new.ids <- convertGeneIds(rownames(cm), db, data.gene.id.type, marker.gene.id.type) %>%
    .[!is.na(.)] %>% .[!duplicated(.)]

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

wrongBlockError <- function(block, line, msg) {
  stop("Can't parse block for type '", block[[1]], "', line: '", line, "'. ", msg)
}

parseBlock <- function(block) {
  cell.type <- gsub(">", "", block[[1]]) %>% trimws()
  marker.info <- list(expressed=c(), not_expressed=c(), expressed_uniq=c())
  for (l in block[2:length(block)]) {
    l.parts <- strsplit(l, ":", fixed=T)[[1]] %>% trimws()
    if (length(l.parts) != 2) {
      if ((nchar(l.parts[[1]]) == nchar(l)) || (length(l.parts) > 2))
        wrongBlockError(block, l, "Exactly one ':' must be present.")

      wrongBlockError(block, l, "Empty entry after the ':' separator.")
    }

    if (l.parts[[1]] == "expressed") {
      genes <- strsplit(l.parts[[2]], ",", fixed=T)[[1]] %>% trimws()
      if (length(genes) == 0)
        wrongBlockError(block, l, "No genes found.")

      marker.info$expressed %<>% union(genes)
    } else if (l.parts[[1]] == "expressed uniq") {
      genes <- strsplit(l.parts[[2]], ",", fixed=T)[[1]] %>% trimws()
      if (length(genes) == 0)
        wrongBlockError(block, l, "No genes found.")

      marker.info$expressed_uniq %<>% union(genes)
    } else if (l.parts[[1]] == "not expressed") {
      genes <- strsplit(l.parts[[2]], ",", fixed=T)[[1]] %>% trimws()
      if (length(genes) == 0)
        wrongBlockError(block, l, "No genes found.")

      marker.info$not_expressed %<>% union(genes)
    } else if (l.parts[[1]] == "subtype of") {
      if (!is.null(marker.info$parent))
        wrongBlockError(block, l, "Only one 'subtype of' line can be presented.")

      marker.info$parent <- l.parts[[2]]
    } else {
      next
    }
  }

  if (is.null(marker.info$parent)) {
    marker.info$parent <- "root"
  }

  if ((length(marker.info$expressed) == 0) && (length(marker.info$expressed_uniq) == 0))
    wrongBlockError(block, l, "No positive markers specified.")

  return(list(marker.info) %>% setNames(cell.type))
}

#' Parse Marker File
#' @description read markers from the markup file
#' @param path path to the markup file
#' @return list of markers for each cell type. List elements have the following info: \itemize{
#'   \item{expressed: vector of positive markers}
#'   \item{not_expressed: vector of negative markers}
#'   \item{parent: parent cell type (`"root"` if there is no parent)}
#' }
#' @export
parseMarkerFile <- function(path, is.text=FALSE) {
  lines <- if (is.text) unlist(strsplit(path, "\n")) else readLines(path)
  markup.lines <- lines %>% .[nchar(.) > 0] %>%
    strsplit("#", fixed=T) %>% sapply(`[[`, 1) %>%
    trimws() %>% .[nchar(.) > 0]

  start.ids <- which(substr(markup.lines, 1, 1) == ">")
  if (length(start.ids) < 2)
    stop("At least two type entries must be presented in the markup file")

  blocks <- mapply(function(s, e) markup.lines[s:(e-1)],
                   start.ids, c(start.ids[2:length(start.ids)], length(markup.lines) + 1),
                   SIMPLIFY=F)

  marker.list <- Reduce(c, lapply(blocks, parseBlock))

  # Process "expressed uniq" tag
  for (n in names(marker.list)) {
    if (length(marker.list[[n]]$expressed_uniq) > 0) {
      for (n2 in names(marker.list)) {
        if ((n2 != n) && (marker.list[[n2]]$parent == marker.list[[n]]$parent)) {
          marker.list[[n2]]$not_expressed %<>% union(marker.list[[n]]$expressed_uniq)
        }
      }
    }

    marker.list[[n]]$expressed %<>% union(marker.list[[n]]$expressed_uniq)
    marker.list[[n]]$expressed_uniq <- NULL
  }

  return(marker.list)
}

#' Create Classification Tree
#' @description creates graph object for cell type hierarchy
#' @param marker.list list of markers per cell type. Can be obtained with `parseMarkerFile`
#' @return igraph graph with the type hierarhcy
#' @examples
#'   markers <- parseMarkerFile(marker_path)
#'   clf_tree <- createClassificationTree(markers)
#'
#' @export
createClassificationTree <- function(marker.list) {
  if (length(marker.list) == 0)
    stop("Empty marker list")
  parents <- sapply(marker.list, function(pl) if(length(pl$parent) == 0) "root" else pl$parent)
  tree <- c(parents, names(marker.list)) %>% matrix(ncol=2) %>% t() %>% igraph::make_directed_graph()

  n.parents <- igraph::as_adj(tree) %>% Matrix::colSums()
  hanging.verts <- which((n.parents == 0) & (names(n.parents) != "root"))
  if (length(hanging.verts) > 0)
    stop("Cell types ", paste0("'", names(hanging.verts), "'") %>% paste(collapse=", "), " are mentioned as parents, but entries are missing")

  return(tree)
}

#' Get Classification data
#' @description prepare information neccessary for cell type classification
#'
#' @param cm gene count matrix with cells by columns and genes by rows. May be in raw, TC-normalized or tf-idf-normalized format (in case of tf-idf, `prenormalized` must be set to `T`)
#' @param markers path to the file with marker genes or parsed marker list from `parseMarkerFile` function
#' @param prenormalized is `cm` in tf-idf-normalized format? Default: FALSE.
#' @export
getClassificationData <- function(cm, markers, prenormalized=F, data.gene.id.type="SYMBOL", marker.gene.id.type="SYMBOL", db=NULL, verbose=F) {
  if (!prenormalized) {
    cm %<>% normalizeTfIdfWithFeatures()
  }

  gi <- unifyGeneIds(cm, data.gene.id.type=data.gene.id.type, marker.gene.id.type=marker.gene.id.type, db=db, verbose=verbose)

  if (is.character(markers)) {
    markers <- parseMarkerFile(markers)
  } else if (!is.list(markers)) {
     stop("Unknown format of markers")
  }

  missed.genes <- c()
  present.genes <- c()

  # Filter markers not presented in the dataset
  for (n in names(markers)) {
    for (n2 in c("expressed", "not_expressed")) {
      missed.genes %<>% union(setdiff(markers[[n]][[n2]], rownames(cm)))
      markers[[n]][[n2]] %<>% intersect(rownames(cm))
      present.genes %<>% union(markers[[n]][[n2]])
    }

    if (length(markers[[n]]$expressed) == 0) {
      subtypes <- names(markers)[sapply(markers, function(x) x$parent == n)]
      warning("Cell type '", n, "' is not presented in the data. Drop it and its ", length(subtypes), " subtypes")
      markers[[n]] <- NULL
      markers[subtypes] <- NULL
    }
  }

  if (length(missed.genes) != 0) {
    warning("Dropped ", length(missed.genes), " genes missed from the dataset. ", length(present.genes), " marker genes left.\n")
  }

  classification.tree <- createClassificationTree(markers)

  res <- list(cm=Matrix::t(gi$cm), classification.tree=classification.tree, gene.table=gi$gene.table, marker.list=markers)
  return(res)
}

### Save markers

markersToMarkup <- function(markers, name) {
  if (!is.null(markers$parent) && (markers$parent == "root")) {
    markers$parent <- NULL
  }

  expr <- paste0("expressed: ", paste0(markers$expressed, collapse=", "), "\n")
  if (!is.null(markers$not_expressed) && length(markers$not_expressed) > 0) {
    not.expr <- paste0("not expressed: ", paste0(markers$not_expressed, collapse=", "), "\n")
  } else {
    not.expr <- ""
  }

  parent <- if (!is.null(markers$parent)) paste0("subtype of: ", markers$parent[1], "\n") else ""

  return(paste0("> ", name, "\n", expr, not.expr, parent, "\n"))
}

#' Marker List to Markup
#' @description Convert marker list to the markup language and optionally save it to a file
#'
#' @param marker.list list of markers per cell type. Can be obtained with `parseMarkerFile`
#' @param file file to save to. If empty string, returns markup text instead of saving
#' @param group.by.parent group cell types in the markup by the parent type
#' @return path to the file if `file == ""` or markupt text otherwise
#' @examples markerListToMarkup(clf_data$marker.list, file="markers.txt")
#' @export
markerListToMarkup <- function(marker.list, file="", group.by.parent=T) {
  if (!group.by.parent) {
    markup.text <- names(marker.list) %>% sort() %>%
      lapply(function(n) markersToMarkup(marker.list[[n]], n)) %>% paste(collapse="")
  } else {
    for (n in names(marker.list)) {
      if (length(marker.list[[n]]$parent) == 0) {
        marker.list[[n]]$parent <- "root"
      }
    }

    ml.per.parent <- marker.list %>% split(sapply(., `[[`, "parent")) %$%
      c(list(root="root"), .[names(.) != "root"])
    markup.text <- names(ml.per.parent) %>% lapply(function(pn)
      paste0("## ", pn, "\n\n", markerListToMarkup(ml.per.parent[[pn]], group.by.parent=F))) %>%
      paste0(collapse="")
  }

  if (!is.null(file) && nchar(file) > 0) {
    cat(markup.text, file=file)
    return(file)
  }

  return(markup.text)
}
