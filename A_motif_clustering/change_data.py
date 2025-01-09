import pandas as pd

# Function to create a FASTA file from an Excel file
def create_fasta_from_excel(excel_file, output_fasta):
    # Load the Excel file into a pandas DataFrame
    df = pd.read_excel(excel_file)

    # Open the output file for writing the FASTA sequences
    with open(output_fasta, 'w') as fasta_file:
        # Loop through each row of the dataframe
        for index, row in df.iterrows():
            # Extract the required columns (adjust indexing based on 0-indexing in pandas)
            peptide_id = row['Id']  # 'Id' column for peptide ID
            peptide_sequence = row['virus_collapsed_hit']  # 'virus_collapsed_hit' column for peptide sequence
            n_acid = row['nucleic_acid']  # 'nucleic_acid' column for label
            virus_family = row['virus_taxa_family'] 
            human_bait = row['bait']
            # Write the FASTA format entry
            fasta_file.write(f">{peptide_id} | {1} | {n_acid}\n")
            fasta_file.write(f"{peptide_sequence}\n")


# Example usage
path_input = 'juliaancker/Desktop/projekt_bio_inf/input_data.xlsx'


create_fasta_from_excel(path_input, 'output_sequences.fasta')
