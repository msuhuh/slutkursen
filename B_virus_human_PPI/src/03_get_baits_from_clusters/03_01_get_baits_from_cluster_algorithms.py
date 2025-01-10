import os
import igraph as ig
import random
import win32file
win32file._setmaxstdio(8192)
import pandas as pd

def get_results(file_path):

    # Read the file content into a structured format
    with open(file_path, 'r') as file:
        raw_data = file.readlines()

    # Parse raw data into a structured dictionary
    parsed_data = {}
    current_tuple_size = None

    for line in raw_data:
        line = line.strip()
        if not line:
            continue

        if line.startswith("Shared"):
            current_tuple_size = line.split()[1].replace("-tuples:", "").strip()
            parsed_data[current_tuple_size] = []
        elif current_tuple_size:
            virus_families, algorithms = line.split(";")
            virus_families = virus_families.strip("[] ")
            algorithms = algorithms.strip("[] ").split(", ")
            parsed_data[current_tuple_size].append((virus_families, algorithms))

    # Convert parsed data into a tabular DataFrame for visualization
    parsed_rows = []
    for tuple_size, entries in parsed_data.items():
        for virus_families, algorithms in entries:
            parsed_rows.append({"Tuple Size": tuple_size, "Virus Families": virus_families, "Algorithms": ", ".join(algorithms)})

    df = pd.DataFrame(parsed_rows)
    #print(df)
    return df

def create_virus_bait_matrix(graph, algorithm, virus_names, bait_names, output_file):
    """
    Creates a matrix of virus families (columns) and unique baits (rows), recording adjacency counts from the graph.

    Parameters:
        graph: igraph.Graph
            The input graph containing nodes and edges.
        algorithm: Callable
            The community detection algorithm to be applied on the graph.
        virus_names: list
            List of known virus names.
        bait_names: list
            List of known bait names.
        output_file: str
            Path to the output CSV file.

    Returns:
        pd.DataFrame: The virus-bait adjacency matrix.
    """
    # Run the community detection algorithm
    communities = algorithm(graph)

    # Collect all virus families and unique baits
    virus_families = {graph.vs[node]["name"] for node in range(graph.vcount()) if graph.vs[node]["name"] in virus_names}
    unique_baits = {graph.vs[node]["name"] for node in range(graph.vcount()) if graph.vs[node]["name"] in bait_names}
    
    # Initialize adjacency tracking
    adjacency_counts = {bait: {virus: 0 for virus in virus_families} for bait in unique_baits}

    # Populate adjacency counts
    for community in communities:
        community_nodes = [graph.vs[node]["name"] for node in community]
        
        viruses = set(node for node in community_nodes if node in virus_names)
        baits = set(node for node in community_nodes if node in bait_names)

        for bait in baits:
            bait_index = graph.vs.find(name=bait).index
            for virus in viruses:
                virus_index = graph.vs.find(name=virus).index
                edge = graph.get_eid(bait_index, virus_index, directed=False, error=False)
                
                if edge != -1:  # Check if an edge exists
                    weight = graph.es[edge]["weight"] if "weight" in graph.es[edge].attributes() else 1
                    adjacency_counts[bait][virus] += weight

    # Convert the matrix to a DataFrame
    matrix_df = pd.DataFrame.from_dict(adjacency_counts, orient="index").fillna(0)

    # Save the DataFrame to a CSV file
    matrix_df.to_csv(output_file)

    return matrix_df

def get_viruses_baits_matrices(
        original_path,
        output_file,
        txt_virus_families,
        txt_baits,
        algorithm,
        seed
    ):


    # Set seed for reproducibility
    if seed is not None:
        random.seed(seed)

    # Check paths
    if not os.path.exists(original_path):
        raise FileNotFoundError(f"Original graph file not found: {original_path}")
    if not os.path.exists(txt_virus_families):
        raise FileNotFoundError(f"File not found: {txt_virus_families}")
    if not os.path.exists(txt_baits):
        raise FileNotFoundError(f"File not found: {txt_baits}")

    # Load unique virus and bait names
    with open(txt_virus_families, "r") as f:
        virus_names = set(f.read().strip().split(","))

    with open(txt_baits, "r") as f:
        bait_names = set(f.read().strip().split(","))

    # Load the original graph and its cluster pairs
    original_g = ig.Graph.Read_Ncol(original_path, weights=True, directed=False)
    # Get the results
    df = create_virus_bait_matrix(original_g, algorithm, virus_names, bait_names, output_file)
    return df


############## Main script #################
# Set a seed for reproducibility
seed = 2024

######### PARAMETERS ###############
original_network = "src/00_data_folder/01_ouput_graphs/normalized_graph_ncol.txt"
txt_virus_families = "src/00_data_folder/00_original_data/unique_virus_fam_names.txt"
txt_baits = "src/00_data_folder/00_original_data/unique_bait_names.txt"
output_root_folder = "src/00_data_folder/03_01_get_baits_from_clusters"

algorithms = [
        (lambda g: g.community_multilevel(weights=g.es["weight"]), "louvain"),

        (lambda g: g.community_fastgreedy(weights=g.es["weight"]).as_clustering(), "fastgreedy"),

        (lambda g: g.community_leading_eigenvector(weights=g.es["weight"]), "eigenvector"),

        (lambda g: g.community_walktrap(steps=2, weights=g.es["weight"]).as_clustering(), "walktrap2"),

        (lambda g: g.community_leiden(objective_function='modularity', weights=g.es["weight"],
                                         resolution=1, beta=0.01, initial_membership=None, n_iterations=-1, node_weights=None),
                                           "leiden_mod_res1"),
        
        (lambda g: g.community_leiden(objective_function='CPM', weights=g.es["weight"],
                                         resolution=0.01, beta=0.01, initial_membership=None, n_iterations=-1, node_weights=None),
                                           "leiden_CPM_res0.01"),

        (lambda g: g.community_infomap(edge_weights=g.es["weight"], trials=10), "infomap10"),
        
              ]

###############################

df_list = []
for algorithm, algorithm_name in algorithms:
    print("Starting alogrithm: ", algorithm_name)
    output_file = os.path.join(output_root_folder, algorithm_name + '.csv')

    # Compare original with random graphs
    df_out = get_viruses_baits_matrices(
        original_path=original_network,
        output_file=output_file,
        txt_virus_families=txt_virus_families,
        txt_baits=txt_baits,
        algorithm=algorithm,
        seed=seed
    )

    df_list.append(df_out)


# Align and combine DataFrames
combined_df = None

for df in df_list:
    if combined_df is None:
        combined_df = df.copy()  # Initialize the combined DataFrame
    else:
        # Align the current DataFrame with the combined DataFrame
        df = df.reindex_like(combined_df)
        # Add the current DataFrame to the combined one
        combined_df = combined_df.add(df, fill_value=0)

# Convert any non-zero values to 1 (binary presence)
#combined_df = combined_df.applymap(lambda x: 1 if x > 0 else 0)

# Save the combined DataFrame
combined_output_file = os.path.join(output_root_folder, "combined_virus_bait_matrix.csv")
combined_df.to_csv(combined_output_file)

print("Combined DataFrame saved to:", combined_output_file)