import pandas as pd


# Load the Excel file
file_path = "20241122_motif_mapped_data.xlsx"  # Replace with your file path
df = pd.read_excel(file_path)

#Name you fasta file, 

fasta_file = 'Data_set_fasta_domain.fa' # Domain choosen as label 

# Open fasta file to write
with open(fasta_file, 'w') as fasta:
    # Iterate through each row in the DataFrame
    for _, row in df.iterrows():
        # Extract values from the row
        virus_motif_instance = row['new_motif'] # Peptide sequence used for clustering
        human_acc = row['human_accession'] # Human accession number can be choosen as label
        id_value = row['Id']  # Id can be choosen as label also included in header
        viral_family = row['virus_taxa_family'] # Virus family name can be choosen as label
        domain = row['human_domain_name'] # Human domain name can be choosen as label        
        
        # Data instances that lack motifs have a '*' instead. These rows are skipped
        if virus_motif_instance != '*':

            # Create the header, first is id value, next is sequence count (always one if clustering algorithm is run with all steps)
            # Last value in header is the label of this sequence, here domain is choosen, this can be changed after preferance
            header = f">{id_value} | {1} | {domain}\n"
            
            # Write the header and the sequence to the FASTA file
            fasta.write(header)
            fasta.write(virus_motif_instance + '\n')

print(f"FASTA file has been created successfully.")
