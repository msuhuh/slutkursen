import pandas as pd

file_mat = 'subs_1_testset.txt' # File to your new substitution matrix

file_path = "202041125_test_set.xlsx"  # Replace with your Excel file path
df = pd.read_excel(file_path)


# Function to read the BLOSUM62 matrix and extract diagonal values
def extract_diagonal_values(filename):
    diagonal_values = []
    
    with open(filename, 'r') as file:

        # Reading the lines from the file
        lines = file.readlines()
        
        # Skip comment lines that start with '#'
        data_lines = [line for line in lines if not line.startswith('#')]
        
        # Loop through the remaining lines to extract diagonal values
        for i, line in enumerate(data_lines[1:], start=0):
            # Split each line into columns
            columns = line.split()
            
            # The first element in each row is the amino acid name (skip it)
            # Extract the diagonal value (same amino acid substitution) at index i
            try:
                value = int(columns[i+1])  # columns[i+1] because columns[0] is the amino acid name
                diagonal_values.append(value)
            except ValueError:
                # If the value is not an integer (for instance, '*'), skip it
                pass
    
    return diagonal_values

# Diagonal_values extractd from function 

diagonal_values = sorted([5, 6, 6, 9, 5, 5, 6, 8, 4, 4, 5, 5, 6, 7, 4, 5, 11, 7, 4])

motifs = df['new_motif']

aa_freq = {'A':0, 'R':0, 'N' :0, 'D':0, 'Q':0, 'E':0, 'G':0, 'H':0, 'I':0, 'L':0, 'K':0, 'M':0, 'F':0, 'P':0, 'S':0, 'T':0, 'W':0, 'Y':0, 'V':0 }

# Calculate the frequency of each amino acid in dataset
for motif in motifs:
    for aa in str(motif):
        aa_freq[aa] += 1

# Calculate the total frequency and normalize
total_frequency = sum(aa_freq.values())  # Total count of amino acids

# Create a dictionary of normalized frequencies
normalized_freq = {aa: ((freq / total_frequency)) for aa, freq in aa_freq.items()}

inverted_values = {aa: 1 / norm_val for aa, norm_val in normalized_freq.items()}

# Re-normalize the inverted values to ensure the sum equals 1
total_inverted = sum(inverted_values.values())
renormalized_inverted = {aa: inv_val / total_inverted for aa, inv_val in inverted_values.items()}

# Check the sum
sum_inverted = sum(renormalized_inverted.values())


# Create a mapping from amino acid names to diagonal values

amino_acids = ['A', 'R', 'N', 'D', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

# Adjust diagonal values based on the frequency distribution
adjusted_diagonal_values = {}

# First, we calculate the total of the diagonal values
total_diag = sum(diagonal_values)

# Sort the amino acids 
sorted_freqs = dict(sorted(renormalized_inverted.items(), key=lambda item: item[1]))

# Adjusted diagonal values based on the frequencies 
i = 0
for keys in sorted_freqs:
    adjusted_diagonal_values[keys] = diagonal_values[i]
    i+= 1

amino_acids = ['A', 'R', 'N', 'D', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

adjusted_diagonal_values1 = {aa: adjusted_diagonal_values[aa] for aa in amino_acids}

print(adjusted_diagonal_values1)