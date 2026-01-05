# Single-cell RNA-seq analysis: understanding and guidance

This tutorial reflects my personal understanding of scRNA-seq analysis with Seurat (for full details, see [here](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)). Before you begin, please install the necessary packages based on the [instructions](https://github.com/steverozen/single-cell-2025/blob/main/install-R-etc.md). It's also recommended to get familiar with [the basics of single-cell RNA-seq](https://docs.google.com/presentation/d/1FaS36SJNyT10_rJlYw8EqaE7sMSstwHln9xE9ccBJA8/edit?slide=id.g3af029a619d_2_71#slide=id.g3af029a619d_2_71). Special thanks to Dr. Steven Rozen for his guidance and valuable insights during *the Duke Genome Academy workshop*.

## Setup the Seurat Object

**Seurat** is a popular R package for analyzing single-cell RNA-seq data. It provides tools to read, organize, and analyze data, including quality control, dimensionality reduction, and clustering.

We start by loading our data. The `Read10X()` function reads the output from 10X Genomics’ Cell Ranger pipeline and gives us a **UMI (Unique Molecular Identifier) count matrix**. In this matrix, each row is a gene, each column is a cell, and the numbers tell us how many RNA molecules of each gene were detected in each cell.

If your Cell Ranger output is in the newer **.h5 file format**, you can use `Read10X_h5()` instead to read the data into Seurat.

Next, we use this count matrix to create a **Seurat object**. Think of this object as a container that holds both your data (like the count matrix) and your analysis results (like PCA or clustering). Everything for a single-cell dataset is stored in this object, making it easy to work with.

```r
library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "/path/to/filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
```

## Standard pre-processing workflow
Below are the basic pre-processing steps in Seurat, covering quality control, data normalization, scaling, and the identification of highly variable genes.

### 1. QC and selecting cells for further analysis

**Seurat** allows you to easily assess data quality and filter out low-quality cells based on user-defined QC criteria.

- A commonly used QC metric is **the number of genes detected per cell**:
  - Cells with very few detected genes are often low-quality cells or empty droplets.
  - Cells with unusually high numbers of genes may be **doublets** (more than one cell captured together).

- Another important metric is the **total number of RNA molecules** detected in each cell, which is usually strongly correlated with the number of genes.

- The **percentage of mitochondrial gene expression** is also widely used:
  Cells with a high mitochondrial percentage are often low-quality or dying cells.

In Seurat, mitochondrial content is calculated using the PercentageFeatureSet() function.

By default, genes starting with **MT-** are used to define mitochondrial genes.

```r
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
```

### 2. Normalizing the data

After filtering out low-quality cells, the next step is to normalize the data. Seurat uses a default method called LogNormalize, which adjusts gene expression values so that cells can be compared fairly. This method scales each cell by its total RNA counts, applies a fixed scale factor (10,000 by default), and then log-transforms the data.

```r
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#or directly using
pbmc <- NormalizeData(pbmc)
```

### 3. Identification of highly variable features (feature selection)

We then identify genes that vary a lot across cells, meaning they are highly expressed in some cells but low in others. Focusing on these highly variable genes helps capture meaningful biological signals, and in Seurat this is done using the `FindVariableFeatures()` function, which by default selects 2,000 genes for downstream analyses such as PCA.

```r
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
```

### 4. Scaling the data

Next, we scale the data before running dimensionality reduction methods such as PCA. This step centers and standardizes gene expression values so that all genes are on the same scale, preventing highly expressed genes from dominating the analysis and allowing meaningful biological patterns to be detected more easily.

```r
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
```

### 5. Perform linear dimensional reduction

Next, we run **PCA** on the scaled data to reduce the complexity of the dataset and capture the main sources of variation between cells. PCA uses the most variable genes by default and summarizes groups of genes that tend to be expressed together (or oppositely), helping us understand the major biological patterns in the data.

```r
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca") + NoLegend()

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
```

### 6. Determine the 'dimensionality' of the dataset

Seurat groups cells using their **PCA results**, which helps reduce technical noise by summarizing many related genes into a smaller number of meaningful components. To decide how many principal components to keep, we often use an **Elbow plot**, which shows how much variation each PC explains and helps identify a cutoff point (the “elbow”) where most of the biological signal has already been captured.

```r
ElbowPlot(pbmc)
```

### 7. Cluster the cells

This is an important step in scRNA-seq analysis. Seurat uses a **graph-based clustering** method to group similar cells together based on their PCA results. First, cells are connected to their most similar neighbors to build a graph, and then this graph is divided into clusters using algorithms like **Louvain**, with a resolution parameter controlling how many clusters are produced (higher values give more clusters).

```r
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)
```

### 8. Run non-linear dimensional reduction (UMAP/tSNE)

Seurat provides **non-linear dimensionality reduction methods** like **t-SNE** and **UMAP** to visualize single-cell data in 2D. These methods help place similar cells close together, making it easier to explore and interpret clusters identified earlier. The reason for this step is to **see patterns and relationships between cells**, but it’s important to remember that these plots simplify the data and should not be used alone to make biological conclusions.

```r
pbmc <- RunUMAP(pbmc, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")
```

### 9. Finding differentially expressed features (cluster biomarkers)

Seurat can identify **marker genes** that distinguish each cluster by comparing gene expression between one cluster and all others. This step is important because it helps **understand what makes each cell group unique** and can reveal the biological identity or function of the clusters. Functions like `FindAllMarkers()` automate this for all clusters, and Seurat v5 uses the presto package to make this analysis much faster for large datasets. 

Notes: Clustering groups cells (step 7) based on overall expression patterns, while marker genes identify the key genes that distinguish each cluster.

```r
# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)

# tests for differential expression
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

# tools for visualizing marker expression
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))

# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
    "CD8A"))

pbmc.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
```

### 10. Assigning cell type identity to clusters

In practice, if you already know the cell types and their typical marker genes, you can simply name the clusters accordingly. Luckily, for this dataset, we can use well-known marker genes to match each unbiased cluster to its corresponding cell type:

```r
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
    "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
```

Generate plot:

```r
library(ggplot2)
plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
    theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
ggsave(filename = "../output/images/pbmc3k_umap.jpg", height = 7, width = 12, plot = plot, quality = 50)
```

