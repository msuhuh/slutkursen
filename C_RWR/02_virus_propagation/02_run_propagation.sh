#!/bin/bash

# This script is a pipeline for running a mofified RWR algorithm, inspired by
# the one used in the paper 'Large-scale phage-based screening reveals extensive 
# pan-viral mimicry of host short linear motifs' by Mihalič et al. Nature 
# communications 2023.The script runs 02_virus_propagation.py and does it in parallell 
# for efficeny. 

# Written by: Minna Sayehban, student at Uppsala University
# E-mail: minna.sayehban.1224@student.uu.se or minna.sayehban@gmail.com
# Date: 28-11-2024
# Modified: 19-12-2024


########## INSTALL & LOAD MODULES ###########
pip install -r requirements.txt


########### GENERATE HUMAN PPI NETWORKS #############

# PATHS & DATA SETS 
in_data_path="/data/family/"
uniprot_file="/data/uniprot_to_gene.tab"
network_path="/networks/TCSS"
out_data_path="results/family/"

# RUN PROPAGATION SCRIPT IN PARALLELL 
process_family() {
    virus_family="$1"
    python propagation_virus.py "$virus_family" "$in_data_path" "$uniprot_file" "$network_path" "$out_data_path"
    echo "Processed virus family: $(basename "$virus_family") successfully."
}

export -f process_family

find "$in_data_path" -mindepth 1 -maxdepth 1 -type d | parallel -j 4 process_family {}

# echo "Script ran successfully! All virus families have been processed."








