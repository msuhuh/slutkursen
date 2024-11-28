#!/bin/bash

# This script is a pipeline for running the clustering method based on human PPI networks 

# Written by: Minna Sayehban, student at Uppsala University
# E-mail: minna.sayehban.1224@student.uu.se or minna.sayehban@gmail.com
# Date: 28-11-2024


########## INSTALL & LOAD MODULES ###########
#pip install -r requirements.txt


########### 1. GENERATE HUMAN PPI NETWORKS #############

# 1.1 PATHS & DATA SETS 
in_data_path="../data/datasets/virus/family/"
uniprot_file="../data/uniprot_to_gene.tab"
network_path="../networks/TCSS"
out_data_path="../test_results/virus/family/"

# 1.2 RUN RWR ALGORITHM
for virus_family in "$in_data_path"*; do
	python propagation_virus.py "$virus_family" "$in_data_path" "$uniprot_file" "$network_path" "$out_data_path"
	echo "Processed virus family: $(basename "$virus_family") successfully."
done
echo "Script ran successfully! All virus families have been processed and generated human PPI networks."


########### 2. ANALYSE SIGNATURE FEATURES #############

# 2.1 PATHS & DATA SETS 


# 2.2 RUN NMF ALGORITHM









