## How to use the code
This repository consists of the scripts and data used in order to run the human PPI clustering algorithm. The scripts divided in different repositories labeled in chronological order and should be followed accordingly:

### 01_preprocessing
- 01_collaps_vir_fam.py --> Collapses all human protein acession numbers related to a virus family in one file per family.

### 02_virus_propagation
  - 02_run_propagation.sh --> Runs the virus_propagation.py script (RWR algorithm) for each virus family.
  - virus_propagation.py --> Runs the human PPI network algorithm with RWR. 
  - fisher_test_virus.py --> Script used in virus_propagation.py to find functionally enriched pathways

### 03_cluster_generation:
- 03_1_virus_signature.py
- 03_2_run_cophenetic_correlation.py
- 03_3_cophenetic_correlation_figure.py
- cophenetic_correlation.py
