# How to use the code
This repository consists of the scripts and data used in order to run the human PPI clustering algorithm. The scripts divided in different repositories labeled in chronological order and should be followed accordingly:

## 01_preprocessing
### `01_collaps_vir_fam.py` 
Collapses all human protein acession numbers related to a virus family in one file per family.

**Necessary input data:**
- `in_file_path` --> Path to the virus dataset. 
- `out_file_path` --> Path to created family files with collapsed human accession numbers. 

**How to run the script:**
 Run `python 01_collaps_vir_fam.py`


## 02_virus_propagation
- `02_run_propagation.sh` --> Runs the virus_propagation.py script (RWR algorithm) for each virus family.

**Necessary input data:**
- `in_data_path` --> Path to virus family directory with all virus family files containing human protein accession numbers. 
- `uniprot_file` --> Path to file with all human acession numbers and their gene name. 
- `network_path` --> Path to netwrok firectory containing all PPI networks. 
- `out_data_path` --> Path to results directory. 

**How to run the script:**
Simply run `bash 02_run propagation.sh`

**Other scripts:**
- `virus_propagation.py` --> Runs the human PPI network algorithm with RWR and is used in the `02_run_propagation.sh`script. 
- `fisher_test_virus.py` --> Finds functionally enriched pathways and is used in the `virus_propagation.py`script.


## 03_cluster_generation:
- `03_1_virus_signature.py` --> Uses the propagation results to create a significance matrix and the clusters.

**Necessary input data:**
- Make sure the paths in the script lead to your specific results from the `02_run_propagation.sh`script.

**How to run the script:**
Run for example `03_1_virus_signature.py 0.3` where 0.3 represents the global frequency threshold, meaning only proteins found in >30% of the results of the virus family propagation results are included. 
 
- `cophenetic_correlation.py` --> Analysis of the cophenetic correlation coefficient in order to decide the number of clusters.  
- `03_2_run_cophenetic_correlation.py` --> Runs the `cophenetic_correlation.py` script. 
- `03_3_cophenetic_correlation_figure.py` --> Creates a figure showing the results from the cophenetic correlation coefficient analysis. 
