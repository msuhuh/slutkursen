# Project Overview

This repository contains the code developed for the course **"Applied Bioinformatics."**  
The project was conducted by:  
- Julia Ancker  
- Alexander Blomlöf  
- Minna Sayehban  

**Clients / Co-supervisors:**  
- Ylva Ivarsson  
- Leandro Simonetti  

**Supervisor:**  
- Jan Andersson  

The code in this folder was primarily developed by **Alexander Blomlöf.**

---

## How to Use the Code

### 1. Prepare the Dataset
The code requires the raw dataset provided by the Ivarsson research group at Uppsala University, which is not in the repo.  
Store the dataset in the folder `00_original_data`. In the current version of the code, this pathway was used:  
`src/00_data_folder/00_original_data/20241119_Compiled_viral_results_compiled.xlsx - Sheet1.csv`  

### 2. Run the First Phase
To initiate the analysis, execute the script `first_phase_pipeline.py`.  

### 3. Update for the Second Phase
After running the first phase, update the script `src/03_get_baits_from_clusters/03_03_extract_baits_from_significant_hits.py` based on the output generated in phase one.  
Specifically, the output file:  
`src/00_data_folder/02_02_overlapping_virus_families/analysis_summary_filter_3_pruned.txt`  
shows which virus families belong to which cluster. This information is essential for updating the second-phase script.
Thereafter run 'second_phase_pipeline.py' to complete the pipeline.
The code found in "src/05_data_visualization" can used to visulize the result, but the input folders needs to be updated.

### 4. Use Consistent Randomized Graphs
The results may vary slightly depending on the randomized graphs used to calculate significant clusters.  
To replicate the exact results reported, refer to the **Appendix zip file** attached to the report.

### 5. Adjust Parameters
In all scripts, look for the section labeled `< ### PARAMETERS ### >`.  
Modify the input and output arguments (file paths) in this section as needed.

---

## Contact

**Alexander Blomlöf**  
Email: alexander.blomlof@hotmail.com  
Date: 2025-01-10
