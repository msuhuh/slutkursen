import pandas as pd

# Load the first and second CSV files
matrix = pd.read_csv("local_code\\create_graph_v1\\data_graphs\\domain_virus_family_matrix.csv")
data_with_norm = pd.read_csv("local_code\\data\\proteome_lengths_v2.csv")

# Copy the matrix to create a normalized version
normalized_data = matrix.copy()

# Keep track of unmatched entries
unmatched_entries = []

# Process each row in data_with_norm
for _, norm_row in data_with_norm.iterrows():
    # Extract relevant information
    accession = norm_row["human_accession"]
    virus_family = norm_row["virus_taxa_family"]
    norm_length = norm_row["norm_length"]
    norm_family = norm_row["norm_N_viruses_family"]

    # Check if the accession and virus family exist in the matrix
    if accession in matrix["human_accession"].values and virus_family in matrix.columns:
        # Normalize the corresponding matrix value
        normalized_data.loc[normalized_data["human_accession"] == accession, virus_family] /= (norm_length * norm_family)
    else:
        # Record unmatched entries
        unmatched_entries.append((accession, virus_family))

# Print unmatched entries
print("Unmatched Entries:")
for entry in unmatched_entries:
    print(f"Accession: {entry[0]}, Virus Family: {entry[1]}")

# Save the normalized matrix
normalized_data.to_csv("local_code\\create_graph_v1\\data_graphs\\domain_virus_family_matrix_norm.csv", index=False)
