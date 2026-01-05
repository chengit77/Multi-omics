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












