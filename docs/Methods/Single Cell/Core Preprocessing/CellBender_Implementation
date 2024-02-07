---
title: "CellBender Implementation"
author: "IBDGC Single-Cell Working Group"
date: '2024-02-01'
output: rmarkdown::github_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

<style type="text/css">

body{ /* Normal  */
      font-size: 16px;
      font-face: "bold"
      font-family: "Arial", serif;
}

h1 { /* Header 1 */
  font-size: 28px;
  color: DarkBlue2;
}
h2 { /* Header 2 */
    font-size: 24px;
  color: DarkBlue2;
  font-face: "bold";
}
h3 { /* Header 3 */
  font-size: 24px;
  font-face: "bold";
}
code.r{ /* Code block */
    font-size: 20px;
    font-family: "Arial", serif;
}
pre { /* Code block - determines code spacing between lines */
    font-size: 14px;
}
</style>


# CellBender

This software package was developed to eliminate technical artifacts from single cell data (Fleming et al. 2023). Even with recent progress in standardizing the single cell sequencing, we are still experiencing issues with background noise in count matrices due to the complexity of this technology.

The purpose of CellBender is to take the CellRanger generated unfiltered count matrices to model and remove background noise and biases and improve estimates of gene expression.

## Tutorial links
```{r results="asis", echo = FALSE}
	cat("[GitHub Link](https://github.com/broadinstitute/CellBender)\n")
```

```{r results="asis", echo = FALSE}
	cat("[Documentation Link](https://cellbender.readthedocs.io/en/latest/)")
```


## Module
	remove-background: this removes the counts due to ambient RNA and barcode swapping from raw/unfiltered feature-	barcode matrices in HDF5 format (.h5)

## Installation on terminal/command prompt
  pip install cellbender

## Test dataset: 10K Heart Cells from an E18 mouse (v3 chemistry)

```{r results = "asis", echo = FALSE}
  cat("[Download Test Dataset](https://www.10xgenomics.com/datasets/10-k-heart-cells-from-an-e-18-mouse-v-3-chemistry-3-standard-3-0-0)")
```

## Input for remove-background module
raw_feature_bc_matrix.h5 from the test dataset

## Command: run for each single cell 10x sample
Default setting is encouraged for first run and based on the information from output.html file we can 	add that as setting
	
	cellbender remove-background --input raw_feature_bc_matrix.h5 --output CB_raw_feature_bc_matrix.h5

## Output for remove-background module
	This produces 9 output files as below:
	1. output_report.html: includes plots and warnings/suggestions for improved parameter setting
	2. output.h5: full count matrix as .h5 with background RNA removed
	3. output_filtered.h5: filtered count matrices as .h5 file with background RNA removed and only has the barcodes with >50% posterior probability of containing cells.
	4. output_cell_barcodes.csv
	5. output.log
	6. output_metrics.csv
	7. ckpt.tar.gz
	8. output_posterior.h5

## Use output count matrix in downstream analyses
	The count matrix that ‘remove-background’ module generated can be easily loaded and used for downstream analyses in Seurat

## Make CellBender output compatible with Seurat
	**ptrepack --complevel 5 tiny_output_filtered.h5:/matrix tiny_output_filtered_seurat.h5:/matrix**
	The file tiny_output_filtered_seurat.h5 is now formatted exactly like a CellRanger v3 h5 file, so Seurat can load it

## Load data from the filtered_Seurat.h5 file on to R
	data.file <- 'tiny_output_filtered_seurat.h5'
	data.data <- Read10X_h5(filename = data.file, use.names = TRUE)
	
## Create Seurat object
	obj <- CreateSeuratObject(counts = data.data)
	obj

## Reference
	Stephen J Fleming, Mark D Chaffin, Alessandro Arduini, Amer-Denis Akkad, Eric Banks, John C Marioni, Anthony A Philippakis, Patrick T Ellinor, and Mehrtash Babadi. Unsupervised 	removal of systematic background noise from droplet-based single-cell experiments using CellBender. Nature Methods, 2023. 
```{r results = "asis", echo = FALSE}
  cat("[https://doi.org/10.1038/s41592-023-01943-7](https://doi.org/10.1038/s41592-023-01943-7)")
```



