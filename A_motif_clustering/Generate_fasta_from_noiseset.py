import pandas as pd 

# Script to generate fasta files from the noise dataset created from Generate_noise_set.py

file_path = "test_set_noice_with_clusters.xlsx"  # Replace with your Excel file path
df = pd.read_excel(file_path)

fasta_file = 'test_set_noise_75_clusterid.fa' # Replace with precentage you want

with open(fasta_file, 'w') as fasta:
    # Iterate through each row in the DataFrame
    for _, row in df.iterrows():
        # Extract values from the row
        virus_motif_instance = row['75_percent'] # Replace with precentage you want
        label = row['cluster_id'] # Label with cluster id to know "ground truth" which cluster the original sequence belongs to
        id_value = row['id_value']
        
        # Create the header
        header = f">{id_value} | {1} | {label}\n"
        
        # Write the header and the cleaned sequence to the FASTA file
        fasta.write(header)
        fasta.write(virus_motif_instance + '\n')
    
print(f"FASTA file has been created successfully.")