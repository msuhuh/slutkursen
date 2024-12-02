import pandas as pd

# Paths to the two Excel files
file1 = "local_code/data/20241119_Compiled_viral_results_compiled.xlsx"  # New baits
file2 = "local_code/data/data_collapsed_new.xlsx"  # Not with new baits

# Load the files into DataFrames
df1 = pd.read_excel(file1)
df2 = pd.read_excel(file2)

# Get unique accession counts for each file
unique_accessions_file1 = df1['virus_accession'].nunique()
unique_accessions_file2 = df2['virus_accession'].nunique()

# Print the results
print(f"Number of unique accession IDs in {file1}: {unique_accessions_file1}")
print(f"Number of unique accession IDs in {file2}: {unique_accessions_file2}")
