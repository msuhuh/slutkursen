'''
This script collapses the data provided by Ylva Ivarsson group where a file
is created for each virus family, containing the human accession ID for all
human proteins interacting with the virus taxa within that family.

Written by: Minna Sayehban, student at Uppsala University
E-mail: minna.sayehban.1224@student.uu.se or minna.sayehban@gmail.com
Date: 26-11-2024
'''

############### IMPORT MODULES ###############
import pandas as pd
import os


############### LOAD & EXTRACT DATA ###############
in_file_path = "../data/datasets/virus/test_data.csv" 
out_file_path = "../data/datasets/virus/family"  

os.makedirs(out_file_path, exist_ok=True)

data = pd.read_csv(in_file_path)

vir_fam = "virus_taxa_family"    
vir_taxa_id = "virus_taxa_id"          
human_acc = "human_accession"

filtered_data = data[[vir_fam, vir_taxa_id, human_acc]]

grouped_data = filtered_data.groupby(vir_fam)


############# COLLAPS DATA BASED ON VIRUS FAMILY #############
for vir_fam, group in grouped_data:
    file_name = f"{vir_fam.replace(' ', '_').replace('/', '_')}_human_accessions.txt"
    specific_out_file_path = os.path.join(out_file_path, file_name)

    unique_human_accessions = group[human_acc].unique()

    formatted_data = pd.DataFrame({
        'Column1': [1.0] * len(unique_human_accessions),      
        'Human Accession Number': unique_human_accessions,     
        'Column3': [0.05] * len(unique_human_accessions)
    })

    formatted_data.to_csv(specific_out_file_path, index=False, header=False, sep='\t')

print("Files have been created for each unique virus family successfully.")