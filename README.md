# CellAnnotatoR
Automated marker-based annotation of cell types

## Installation

This package uses [Garnett](https://cole-trapnell-lab.github.io/garnett/docs/) format of marker gene markup files and 
[Conos](https://github.com/hms-dbmi/conos) label propagaton, so you need to install these packages as dependencies:

```r
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")

if (!requireNamespace("devtools"))
  install.packages("devtools")

BiocManager::install(c("monocle"))
devtools::install_github("cole-trapnell-lab/garnett")
devtools::install_github("hms-dbmi/conos")
```

Next, install CellAnnotatoR:

```r
devtools::install_github("khodosevichlab/CellAnnotatoR")
```

## Usage

**NOTE:** this package is still in the development, and some functionality can be changed.

Assuming that you have path to your marker file in `marker_path`, gene count matrix `cm`, cell graph `graph`, clustering `clusters` and embedding `emb`. 
Examples of graphs are: [Seurat](https://github.com/satijalab/seurat/wiki/Seurat) `so@graphs[[1]]`, 
[Pagoda 2](https://github.com/hms-dbmi/pagoda2) `p2$graphs[[1]]` or [Conos](https://github.com/hms-dbmi/conos) `con$graph`.

```r
clf_data <- getClassificationData(cm, marker_path, data.gene.id.type="SYMBOL", marker.gene.id.type="SYMBOL")

marker_scores <- getMarkerScoresPerCellType(clf_data)
ann_by_level <- assignCellsByScores(graph, marker_scores, clf_data$classification.tree, clusters=clusters)

plotAnnotationByLevels(emb, ann_by_level$annotation, clusters=clusters, size=0.2, font.size=c(2, 4), shuffle.colors=T)
```
