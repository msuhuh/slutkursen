import pandas as pd
import random

# Script to create a excell file with motifs with added noise

# Load the Excel file
file_path = "202041125_test_set.xlsx"  # Replace with your Excel file path
df = pd.read_excel(file_path)

motifs = df['new_motif'] # Extract column with motifs

fasta_file = 'test_set_noice_domain.fa'

# Dictionary with amino acid frequencies
aa_freq = {'A':0, 'R':0, 'N' :0, 'D':0, 'Q':0, 'E':0, 'G':0, 'H':0, 'I':0, 'L':0, 'K':0, 'M':0, 'F':0, 'P':0, 'S':0, 'T':0, 'W':0, 'Y':0, 'V':0 }

# Count amino acid frequencies in dataset
for motif in motifs:
    for aa in str(motif):
        aa_freq[aa] += 1

total_frequency = sum(aa_freq.values())  # Total count of amino acids

normalized_freq = {aa: ((freq / total_frequency)) for aa, freq in aa_freq.items()} # Normalise amino acid counts

amino_acids = list(normalized_freq.keys())
frequencies = list(normalized_freq.values())

# Precentages of noise to be added 
precent = [0.25, 0.5, 0.75, 1, 2]

# Create a dictionary to store the modified motifs for each percentage, keep id for each sequence 
modified_motifs = {'id_value': df['Id'], 'new_motif': df['new_motif']}  

# For each percentage, modify the motifs and add to the dictionary
for p in precent:

    modified_column = []  # List to store modified motifs for the current percentage

    # Iterate through each motif
    for _, row in df.iterrows():
        virus_motif = row['new_motif'] # Original motif
        n_added = int(p * len(virus_motif))  # Calculate number of amino acids to add based on precentages 

        # Add random amino acids to random place in the motif
        for _ in range(n_added):
            range_pos = len(virus_motif) - 1 
            pos = random.randint(0, range_pos)
            added_aa = random.choices(amino_acids, weights=frequencies, k=1)[0]

            # add to motif
            virus_motif = virus_motif[:pos] + added_aa + virus_motif[pos:]

        modified_column.append(virus_motif)  # Append the modified motif to the list

    # Add the modified column to the dictionary
    modified_motifs[f"{int(p * 100)}_percent"] = modified_column

# Convert the dictionary into a DataFrame
modified_df = pd.DataFrame(modified_motifs)

# Save the DataFrame to an Excel file
output_file = "test_set_noice.xlsx"
modified_df.to_excel(output_file, index=False)

print(f"Excel file '{output_file}' has been created successfully!")