# CellPhoneDB
## Description
CellPhoneDB is a publicly available repository of curated receptors, ligands and their interactions. 
  
## Overview
Subunit architecture is included for both ligands and receptors, representing heteromeric complexes accurately. This is crucial, as cell-cell communication relies on multi-subunit protein complexes that go beyond the binary representation used in most databases and studies.

CellPhoneDB integrates existing datasets that pertain to cellular communication and new manually reviewed information. Databases from which CellPhoneDB gets information are: UniProt, Ensembl, PDB, the IMEx consortium, IUPHAR.

CellPhoneDB can be used to search for a particular ligand/receptor, or interrogate your own single-cell transcriptomics data.
  
## Dependenncies
- conda
- python v3.6 or later

## Code Base
- cellphonedb_inputs.Rmd
- [github](https://github.com/Teichlab/cellphonedb)

## Workflow
### Installing CellphoneDB
NOTE: Works with Python v3.6 or greater. If your default Python interpreter is for `v2.x` (you can check it with `python --version`),, calls to `python/pip` should be substituted by `python3/pip3`.
We highly recommend using an isolated python environment (as described in steps 1 and 2) using [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) or [virtualenv](https://docs.python.org/3/library/venv.html) but you could of course omit these steps and install via `pip` immediately.
1. Create python=>3.6 environment: `conda create -n cpdb python=3.7`
2. Activate environment: `source activate cpdb`
3. Install CellPhoneDB: `pip install cellphonedb` 

### Input data (from seurat object)
1. Run the cellphonedb_inputs.Rmd to prepare inputs from seurat object:
    Input is the Seurat Object .rds
    Outputs are count_matrix.txt and meta_table.txt

### Running CellphoneDB
1. Activate the environment
    `source activate cpdb`
2. Running with statistical methods 
    `cellphonedb method statistical_analysis YOUR_PATH_TO_OUTPUT/meta_table.txt YOUR_PATH_TO_OUTPUT/count_martix.txt`
    - It’s likely we need to add an parameters: `--counts-data`: [ensembl | gene_name | hgnc_symbol] Type of gene identifiers in the counts data
    - Other optional parameters can be found in the github page:[github](https://github.com/Teichlab/cellphonedb)
3. Plotting with statistical results
    - Currently there are two plot types available: dot_plot & heatmap_plot
    - Once you have the needed files (means & pvalues) you can proceed as follows:
    `cellphonedb plot dot_plot`
    `cellphonedb plot heatmap_plot yourmeta.txt`
    **dot_plot**
    This plot type requires ggplot2 R package installed and working; 