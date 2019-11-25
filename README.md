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
ann_by_level <- assignCellsByScores(graph, clf_data, clusters=clusters)

plotAnnotationByLevels(emb, ann_by_level$annotation, clusters=clusters, size=0.2, font.size=c(2, 4), shuffle.colors=T)
```

See the vignettes for a [Seurat PBMC3k](vignettes/seurat_pbmc3k.md) example and [Conos BM+CB alignment annotation](vignettes/conos_bm_cb.md).

## Creating annotation file

### Extracting markers from a provided annotation

If you have an annotated dataset from the same tissue, but no existing markers [Automated marker selection based on provided annotation](vignettes/mca_marker_selection.md) (MCA Lung data). Please, be aware that marker selection algorithm is under development and will be improved. In case you already have some markers, check the section ["Improving known list of markers"](vignettes/mca_marker_selection.md#improving-known-list-of-markers)

### De-novo annotation

Creating annotation de-novo depends a lot on the type of your data and the more prior knowledge you have about the markers the better.
Here are some sources where you can get some markers to start with:

- [PanglaoDB](https://panglaodb.se/)
- [CellMarker](http://biocc.hrbmu.edu.cn/CellMarker/)
- [ScType](https://sctype.fimm.fi/database.php)
- [Allen Brain Institute interactive datasets](https://portal.brain-map.org/atlases-and-data/rnaseq): very useful in case you work with brain
- [Garnett marker lists](https://cole-trapnell-lab.github.io/garnett/classifiers/): at the moment (Nov 2019) has few datasets and quality is unclear, but hopefully will be improved
- [SCfind](https://scfind.sanger.ac.uk/): service, which helps to identify unknown clusters given a list of their DE markers. At the moment, resolution of most annotations there is relatevely low.

After getting some markers for your data, use [Garnett specification](https://cole-trapnell-lab.github.io/garnett/docs/#constructing-a-marker-file) to create a markup file.
Indeed, the annotation process mostly follows the next workflow:

1. Find marker candidates either with differential expression or using prior knowledge
2. Plot the markers on your data and ensure that they suit your case
3. Add the markers to the file and re-run the classification
4. Check results, find cell types, which are not well-separated or for which you want to increase the annotation resolution
5. Go to step 1 if there is something to improve

One round of this workflow is shown in the [QC vignette](vignettes/mca_qc.md). And here are some more tips for these steps:

#### Step 1

It depends a lot on your problem and packages you use to work with scRNA-seq data. So, few can be mantioned here in general case. Only that Specificity and ROC AUC metrics really help to select good markers (see [Conos walkthrough](https://github.com/hms-dbmi/conos/blob/dev/vignettes/walkthrough.md#cluster-markers) for an example).

#### Step 2

Seurat has its own functions for plotting gene expression, but for general case CellAnnotatoR provides the function `plotGeneExpression(genes, embedding, cm, ...)`. 
It returns list of plots for individual genes. **Note**: matrix `cm` must be transposed, i.e. have genes as columns and cells as rows.

Example with Pagoda 2 object `p2`:

```R
c("Cldn10", "Igfbpl1", "Ccnd2", "Nes") %>% 
  plotGeneExpression(p2$embeddings$PCA$UMAP, p2$counts)
```

If you want to use panel of violinplots instead, you can use `plotExpressionViolinMap(genes, cm.matrix, annotation, ...)`. It suits better for large panels of markers:

```R
c("Cldn10", "Igfbpl1", "Ccnd2", "Nes", "Id4", "Ascl1", "Egfr", "Serpine2", "Dcx", "Tubb3",
  "Slc1a3", "Slc1a2", "Meis2", "Dlx5", "Dlx6") %>% 
  plotExpressionViolinMap(p2$counts, p2$clusters$PCA$leiden)
```

Two functions above allow you to plot markers, which you just want to test on the dataset. But in case you need to plot the markers, which are already in the markup file, two more functions are provided:

- `plotTypeMarkers(embedding, count.matrix, cell.type, marker.list, ...)`: plot markers for a specific `cell.type` from the `marker.list`
- `plotSubtypeMarkers(embedding, count.matrix, parent.type="root", marker.list=NULL, ...)`: plot markers, which separate subtypes within a given cell type `parent.type` and the provided `marker.list`. Setting `parent.type` to `"root"` plots markers for all subtypes. Setting `max.depth` allow to restrict maximal depth of the subbranch, for which markers are plotted.

#### Step 3

Adding markers is trivial and described in the [Garnett specification](https://cole-trapnell-lab.github.io/garnett/docs/#constructing-a-marker-file). However there are some tricks for running classification.

The simplest way to get the annotation is using the following code:

```R
clf_data <- getClassificationData(cm, marker_path)
ann_by_level <- assignCellsByScores(graph, clf_data, clusters=clusters)
```

Though if you need to re-run annotation multiple times (and you probably are) some steps can be optimized.

First, `getClassificationData` performs TF-IDF normalization inside, which doesn't depend on the marker list. So we can do it only once:

```R
cm_norm <- normalizeTfIdfWithFeatures(cm)
clf_data <- getClassificationData(cm_norm, marker_path, prenormalized=T)
```

Second, the most time-consuming step for classification is label propagation on graph. It improves classification quality in many cases, but for getting approximate results we can avoid this. 
To do so, it's enough to pass `NULL` instead of the graph object:

```R
ann_by_level <- assignCellsByScores(NULL, clf_data, clusters=clusters)
```

So, during marker selection it's recommended to re-run the following code every time you update markers:

```R
clf_data <- getClassificationData(cm_norm, marker_path, prenormalized=T)
ann_by_level <- assignCellsByScores(NULL, clf_data, clusters=clusters)
```

#### Step 4

Validation of the results is crucial for high-quality annotation, so this part is described in the [QC vignette](vignettes/mca_qc.md). Here is just a list of functions, which can be useful:

- `plotAssignmentScores(embedding, scores, classification.tree, parent.node)`
- `plotAnnotationByLevels(embedding, annotation.by.level)`
- `plotUncertaintyPerCell(embedding, uncertainty.info)`
- `plotUncertaintyPerClust(uncertainty.per.clust, clusters)`
- `plotAssignmentConfusion(scores, annotation)`
