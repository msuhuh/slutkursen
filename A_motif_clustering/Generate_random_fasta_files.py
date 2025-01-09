import pandas as pd
import random

# Script to generate a 1000 fasta files with randomised motifs

# Load the Excel file
file_path = "20241122_motif_mapped_data.xlsx"  # Replace with your Excel file path
df = pd.read_excel(file_path)

# Extract the motif column in the dataframe
motifs = df['new_motif']

# Create dictionary of amino acid frequency in dataset
aa_freq = {'A':0, 'R':0, 'N' :0, 'D':0, 'Q':0, 'E':0, 'G':0, 'H':0, 'I':0, 'L':0, 'K':0, 'M':0, 'F':0, 'P':0, 'S':0, 'T':0, 'W':0, 'Y':0, 'V':0 }

# Count amount of each amino acid in dataset
for motif in motifs:
    if motif == '*':
        continue
    for aa in str(motif):
        aa_freq[aa] += 1


total_frequency = sum(aa_freq.values())  # Total count of amino acids

# Calculate amino acid precentages 
normalized_freq = {aa: ((freq / total_frequency)) for aa, freq in aa_freq.items()}

# Create list of amino acids and frequencies
amino_acids = list(normalized_freq.keys())
frequencies = list(normalized_freq.values())

# Generate 1000 fasta files with randomised motifs

for i in range(1000):
    fasta_file = f'/Users/juliaancker/Desktop/full_data_random_sets2/{i}R.fa' # Replace with your path to random fasta files
    with open(fasta_file, 'w') as fasta:
                # Iterate through each row in the DataFrame

                for _, row in df.iterrows():
                    # Extract motif from the row
                    motif = row['new_motif']
                    if motif == '*':
                         continue
                    # Create a new motif 
                    virus_motif = ""
                    id_value = row['Id']
                    label = row['human_domain_name']
                    n_added = (len(motif))

                    # Create a new random motif of the same length as the original
                    for n in range(n_added):
                        # Add a random amino acid based on the frequencies of amino acid in teh orignal dataset
                        added_aa = random.choices(amino_acids, weights=frequencies, k=1)[0] 
                        virus_motif += added_aa

                    # Give the random motif the label R for random. 
                    header = f">{id_value} | {1} | 'R' \n"

                    # Keep the original id as it is needed for the bash script run
                    header2 = f">{id_value} | {1} | {id_value} \n"

                    # Add the sequence CCCCC to the original header as there are no C in the original dataset it cannot cluster together with any random cluster and can easialy be identified
                    extra = "CCCCC"
                    
                    # Write the new R header with the random motif and the original header with the CCCCC motif 
                    fasta.write(header)
                    fasta.write(virus_motif + '\n')
                    fasta.write(header2)
                    fasta.write(motif + '\n')



