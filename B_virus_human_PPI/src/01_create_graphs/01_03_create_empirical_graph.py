import pandas as pd

#### PARAMETERS ####

# Load the normalized matrix (replace with your actual path)
normalized_01_matrix = pd.read_csv('src/00_data_folder/01_ouput_graphs/matrices/normalized_data.csv', index_col=0)

# Define the output file path
ncol_file_path = 'src/00_data_folder/01_ouput_graphs/normalized_graph_ncol.txt'

####################

# Create the NCOL format data
edges = []
for virus_family in normalized_01_matrix.columns:
    for bait in normalized_01_matrix.index:
        edge_weight = normalized_01_matrix.loc[bait, virus_family]
        if edge_weight > 0:  # Include only non-zero edge weights
            edges.append(f"{virus_family} {bait} {edge_weight}")



# Save to a text file in NCOL format
with open(ncol_file_path, 'w') as f:
    for edge in edges:
        f.write(edge + '\n')

print(f"NCOL file saved to: {ncol_file_path}")