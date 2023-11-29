
# Cedars Python Pipeline



## outlines

- [x] [preprocess and qc](#preprocess-and-qc)
- [x] [deg analysis](#deg-analysis)
- [x] [integration](#integration)
- [x] [annotation](#annotation)
- [x] [label transfer](#label-transfer)
- [x] [velocity](#velocity)



### preprocess and qc

#### scanpy

- reference:

https://github.com/scverse/scanpy

- [x] read_10x_mtx
- [x] filter
- [x] normalize
- [x] log1p
- [x] regress out
- [x] scale

## fetch the data from ncbi

### Install the SRA Toolkit

``` bash
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz
tar -xzf sratoolkit.current-centos_linux64.tar.gz

```
Sample list:

"GSM6614365_CD-6"

"GSM6614358_UC-5"

"GSM6614352_HC-5"


### Download the data

``` bash
./prefetch PRJNA886695
fastq-dump --split-files SRR15236500
```


### Run the cellranger count

``` bash
wd={work_directory}
for s in $(ls $wd|grep SR|grep fastq|grep "_1")
do
	sample_id=$(echo $s | cut -f 1 -d '_')
	echo $sample_id
	mv "$sample_id"_1.fastq "$sample_id"_S1_L001_R1_001.fastq
	mv "$sample_id"_2.fastq "$sample_id"_S1_L001_R2_001.fastq
	cellranger-7.1.0/bin/cellranger count --id=$sample_id \
                 --transcriptome={tools_directory}/refdata-gex-GRCh38-2020-A \
                 --fastqs={fastq_directory} \
                 --sample=$sample_id
done
```



load the 10x data



``` python
adata_list = []
for matrix_file in [f for f in all_files if 'SRR' in f]:
    sadata = sc.read_10x_mtx(
        'SRR21787570/outs/filtered_feature_bc_matrix',  # the directory with the `.mtx` file
        var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
        cache=True)                              # write a cache file for faster subsequent reading
    
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    
    adata_list.append(sadata)
adata = adata_list[0].concatenate(adata_list[1:], join='outer')

```


filter the data

``` python
adata.var_names_make_unique()
adata = sc.pp.filter_cells(adata, min_genes=200)
adata = sc.pp.filter_genes(adata, min_cells=3)

```


normalize the data

``` python
sc.pp.normalize_total(adata, target_sum=1e4)

```

log1p the data

``` python
sc.pp.log1p(adata)

```

regress out

``` python
sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])

```

scale the data

``` python
sc.pp.scale(adata, max_value=10)

```







### deg analysis

#### scanpy

Find Marker Genes: 

You can use Scanpy to find marker genes that distinguish between different cell types or conditions. This is often a crucial step before differential expression analysis.

``` python

# find marker genes for each group
sc.tl.rank_genes_groups(adata, groupby='group_id', method='wilcoxon')

```

Visualize Marker Genes: 

Visualize the marker genes to get a sense of the differences between the groups.


``` python

sc.pl.rank_genes_groups(adata, n_genes=25)

```



Differential Expression Analysis: 

To perform differential expression analysis, you can use the sc.tl.rank_genes_groups function again, but this time with a different method such as 'wilcoxon' or 'logreg'.

``` python

sc.tl.rank_genes_groups(adata, groupby='condition', method='wilcoxon')

```



### integration (potential batch effect removal)

#### scVI

- reference:

https://docs.scvi-tools.org/en/stable/tutorials/notebooks/scrna/scarches_scvi_tools.html


demo code:

``` python

# We need to perform normalize and log1p before copy the adata to raw
adata.layers["counts"] = adata.X.copy()

# Normalize the total count of each cell to a target value of 1e4
sc.pp.normalize_total(adata, target_sum=1e4)

# Log-transform the data, adding 1 to avoid taking the logarithm of zero counts
sc.pp.log1p(adata)


# preserve the raw count of adata
adata.raw = adata.copy()


# %% [markdown]

logging.info("Identifying highly variable genes")
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    subset=True,
    layer="counts",
    flavor="seurat_v3",
    span=0.8,
    batch_key="sample_id"
)

# Setup the AnnData object for scVI modeling
logging.info("Training scVI model")
scvi.model.SCVI.setup_anndata(
    adata,
    # layer="counts",
    batch_key="sample_id",
    categorical_covariate_keys=["sample_id"],
)
vae = scvi.model.SCVI(adata, n_latent=30, n_layers=2)
# vae = scvi.model.SCVI(adata)
import torch

# vae.to_device(3)
vae.train(batch_size=1024 * 1024 * 2, max_epochs=300, early_stopping=True)

# plot the training curve to identify the optimal epoch
train_elbo = vae.history["elbo_train"][1:]
test_elbo = vae.history["elbo_validation"]
ax = train_elbo.plot()
test_elbo.plot(ax=ax)

logging.info("Finished training scVI model")


if os.path.exists("scVI_model/"):
    shutil.rmtree("scVI_model/")
vae.save("scVI_model/")

logging.info("Finished training scVI model")

adata.obsm["X_scVI"] = vae.get_latent_representation()


sc.pp.neighbors(adata, n_neighbors=20, n_pcs=30, use_rep="X_scVI")

sc.tl.umap(adata, min_dist=0.8)

sc.tl.leiden(adata, resolution=0.8, key_added = "leiden")

sc.pl.umap(
    adata,
    color=["leiden"],
    frameon=False,
    size = 6,
    legend_loc='on data',

)


sc.pl.umap(
    adata,
    color=["sample_id"],
    frameon=False,
    legend_loc='right margin',
)


```


### LISI score calculation

We offer an algorithm to calculate the Local Inverse Simpson’s Index (LISI) for evaluating the mixture quality of cell clusters in a single-cell RNA-seq dataset with respect to certain categorical variables, such as batch, technology, or donor. 

This assessment helps determine if the clusters are well-mixed and integrated across these categories.

This step we could compare the LISI score between the Harmony, scVI, and Seurat v3 integration method.



### integration (potential batch effect removal)
![](images/2023-10-10-17-34-28.png)


### UMAP by sample id
![](images/2023-10-10-17-34-41.png)



### Compute integration metrics
![](images/2023-10-24-22-42-19.png)

![](images/2023-10-24-22-41-44.png)

    KMeans ARI (Adjusted Rand Index):
        The Adjusted Rand Index is a measure of the similarity between two data clusterings. It assesses the quality of clusters by comparing the pairs of data points and measuring whether they are assigned to the same or different clusters in two different clusterings.
        The KMeans ARI specifically applies the Adjusted Rand Index to assess the quality of clusters generated by the K-Means clustering algorithm. It helps determine how well K-Means is able to group data points into meaningful clusters. A higher ARI score indicates better clustering performance.

    Silhouette Score:
        The Silhouette Score is another metric used to evaluate the quality of clusters. It measures how similar an object is to its own cluster (cohesion) compared to other clusters (separation). In other words, it quantifies how well-separated and compact the clusters are.
        The Silhouette Score ranges from -1 to 1, where a higher score indicates that the data points are well-clustered with clear separation, and a lower score suggests overlapping clusters or poorly defined clusters.

    Aggregate Score (Usually, a combination of multiple metrics):
        The Aggregate Score is not a single metric but rather a general term for a combination of multiple evaluation metrics used together to assess the overall quality of clustering. It's often used to provide a more comprehensive view of clustering performance.
        This aggregate score may include metrics such as the Adjusted Rand Index, Silhouette Score, Davies-Bouldin Index, and others. By considering multiple metrics, it becomes possible to gain a more robust understanding of how well a clustering algorithm is performing on a specific dataset.

### annotation

#### TOSICA

- reference:

https://github.com/JackieHanLab/TOSICA


load the library

``` python

import scanpy as sc
import numpy as np
import pandas as pd
import os
import TOSICA

```


get overlapped genes between reference and query adata

``` python



overlapped_genes = list(set(ref_adata.var_names).intersection(set(query_adata.var_names)))
print(len(overlapped_genes))

ref_adata = ref_adata[:,overlapped_genes]

query_adata = query_adata[:,overlapped_genes]

```

train the model based on transformer model


``` python

TOSICA.train(ref_adata, gmt_path='human_gobp', 
             label_name='MajorCelltype',
             epochs=3,
             batch_size = 8 * 16,
             project='hGOBP_demo')


model_weight_path = './hGOBP_demo/model-0.pth'
new_adata = TOSICA.pre(query_adata, 
                       model_weight_path = model_weight_path,
                       project='hGOBP_demo')

```




![](images/2023-10-10-18-11-23.png)

### velocity

#### veloVI

- reference:

https://github.com/YosefLab/velovi


















































