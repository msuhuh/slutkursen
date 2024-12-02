import pandas as pd

# Load the matrix CSV file
matrix = pd.read_csv("local_code\graph_v2_virusFam_and_humanDomain\graph_data\hit_matrix.csv")

# Exclude the 'human_domain' column for normalization calculations
human_domain = matrix["Domain"]  # Preserve the human_domain column
matrix_values = matrix.drop(columns=["Domain"])  # Matrix without human_domain

# Calculate row and column sums
row_sums = matrix_values.sum(axis=1)  # Sum along rows
col_sums = matrix_values.sum(axis=0)  # Sum along columns

# Normalize each cell by its row sum and column sum
normalized_matrix_values = matrix_values.copy()

for i in range(matrix_values.shape[0]):  # Iterate over rows
    for j in matrix_values.columns:  # Iterate over columns
        if row_sums[i] != 0 and col_sums[j] != 0:  # Avoid division by zero
            normalized_matrix_values.loc[i, j] = matrix_values.loc[i, j] / (row_sums[i] * col_sums[j])
        else:
            normalized_matrix_values.loc[i, j] = 0  # Set to 0 if row or column sum is zero

# Min-max normalization to scale between 0 and 1
min_value = normalized_matrix_values.min().min()  # Minimum value in the matrix
max_value = normalized_matrix_values.max().max()  # Maximum value in the matrix

print("old min max")
print(min_value)
print(max_value)

normalized_matrix_values = (normalized_matrix_values - min_value) / (max_value - min_value)

# Min-max normalization to scale between 0 and 1
min_value = normalized_matrix_values.min().min()  # Minimum value in the matrix
max_value = normalized_matrix_values.max().max()  # Maximum value in the matrix

print("new min max")
print(min_value)
print(max_value)

# Combine normalized values with the human_domain column
final_normalized_matrix = pd.concat([human_domain, normalized_matrix_values], axis=1)

# Save the normalized matrix
final_normalized_matrix.to_csv("local_code\graph_v2_virusFam_and_humanDomain\graph_data\HumanDomain_virusFamily_matrix_norm_0_1.csv", index=False)

print("Normalization complete. The normalized matrix (0-1 range) has been saved.")
