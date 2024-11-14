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
            peptide_id = row[0]  # 1st column (ID)
            peptide_sequence = row[10]  # 11th column (peptide sequence)
            label = row[19]  # 20th column (label)

            # Write the FASTA format entry
            fasta_file.write(f">{peptide_id} | {1} | {label}\n")
            fasta_file.write(f"{peptide_sequence}\n")

# Example usage
create_fasta_from_excel('input_data.xlsx', 'output_sequences.fasta')
