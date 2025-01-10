import pandas as pd

#### PARAMETERS ####

# Define the input CSV file path
input_csv_path = 'src/00_data_folder/00_original_data/preprocessed_data.csv'

# Save the interaction matrix to a CSV file
intermediate_csv_path = 'src/00_data_folder/01_ouput_graphs/matrices/original_matrix.csv' # non-normalized

# Save mormalized matrix
output_csv_path = 'src/00_data_folder/01_ouput_graphs/matrices/normalized_data.csv'

##########################

# Load the CSV file into a DataFrame
df = pd.read_csv(input_csv_path)

# Create a pivot table to count interactions
interaction_matrix = (
    df.groupby(["human_bait_id", "virus_taxa_family"])
    .size()
    .unstack(fill_value=0)
)

interaction_matrix.to_csv(intermediate_csv_path)

#############################

# Normalize the edges based on amount of hits in each corresponding node.
# Normalized_w_ij = w_ij / (node_degree_i + node_degree_j)

# Find the maximum value in the matrix
max_value = interaction_matrix.values.max()

# Print the maximum value
print("max value before:")
print(max_value)

# Calculate row sums and column sums
row_sums = interaction_matrix.sum(axis=1)
col_sums = interaction_matrix.sum(axis=0)

# Normalize each cell
normalized_matrix = interaction_matrix.copy()
for i in normalized_matrix.index:
    for j in normalized_matrix.columns:
        normalized_matrix.loc[i, j] = (
            interaction_matrix.loc[i, j] / (row_sums[i] + col_sums[j])
        )

# Perform 0-1 normalization (min-max normalization)
min_value = normalized_matrix.values.min()
max_value = normalized_matrix.values.max()
print("max after row and col normalization:")
print(max_value)
normalized_01_matrix = (normalized_matrix - min_value) / (max_value - min_value)

# Print the maximum value after normalization
new_max = normalized_01_matrix.values.max()
print("max after normalization:")
print(new_max)

normalized_01_matrix.to_csv(output_csv_path, index=True)
print(f"Extracted data saved to: {output_csv_path}")