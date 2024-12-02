import pandas as pd

# Load the normalized matrix
#data = pd.read_csv("local_code\create_graph_v1\data_graphs\domain_virus_family_matrix_norm_0_1.csv")

data = pd.read_csv("local_code\create_graph_v1\data_graphs\\not_correct\domain_virus_family_matrix.csv")

# Prepare the edge list
edges = []

# Iterate through rows (human_accession) and columns (virus families)
for _, row in data.iterrows():
    human_accession = row["human_accession"]  # Get the human accession
    for virus_family in data.columns[1:]:  # Exclude the human_accession column
        weight = row[virus_family]
        if weight > 0:  # Only include edges with a positive weight
            edges.append((human_accession, virus_family, weight))

# Convert the edges list to a DataFrame
edge_list = pd.DataFrame(edges, columns=["source", "target", "weight"])

# Save the edge list to a file in NCOL format (tab-separated)
edge_list.to_csv("local_code\create_graph_v1\data_graphs\hits_as_weights.txt", sep="\t", index=False, header=False)

print("Edge list prepared and saved.")
