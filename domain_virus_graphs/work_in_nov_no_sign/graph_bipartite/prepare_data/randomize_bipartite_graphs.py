import os
import numpy as np
from igraph import Graph
from collections import Counter

import win32file
win32file._setmaxstdio(4096)

def generate_and_save_random_bipartite_graphs(input_file, output_folder, num_random_graphs=2000, seed=None):
    """
    Generate random bipartite graphs while preserving degree sequence and original node names.
    
    Args:
        input_file (str): Path to the original graph file in .ncol format.
        output_folder (str): Folder to save the generated random graphs.
        num_random_graphs (int): Number of random graphs to generate.
        seed (int, optional): Seed for reproducibility.
    """
    # Set random seed for reproducibility
    if seed is not None:
        np.random.seed(seed)
    
    # Create output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)
    
    # Load the original graph
    original_graph = Graph.Read_Ncol(input_file, weights=True, directed=False)
    
    # Load virus names to determine partitions
    with open(r"local_code\create_graph_v1\data_graphs\unique_virus_fam_names.txt", "r") as file:
        virus_names = set(file.read().strip().split(","))

    # Classify nodes into two partitions
    partition_virus_families = set(v["name"] for v in original_graph.vs if v["name"] in virus_names)
    partition_human_domains = set(v["name"] for v in original_graph.vs if v["name"] not in virus_names)

    # Extract degree sequences for both partitions
    degrees_virus_families = original_graph.degree(list(partition_virus_families))
    degrees_human_domains = original_graph.degree(list(partition_human_domains))

    # Extract edges and weights
    original_edges = set(
        (original_graph.vs[edge.source]["name"], original_graph.vs[edge.target]["name"], edge["weight"])
        for edge in original_graph.es
    )
    total_original_edges = len(original_edges)
    original_weights = [edge[2] for edge in original_edges]

    print(f"Loaded bipartite graph with {len(partition_virus_families)} virus families and {len(partition_human_domains)} human domains.")
    print(f"The original graph has {total_original_edges} edges.")

    # Counter to track edge matches
    edge_occurrence_counter = Counter()

    # Function to generate a random bipartite graph
    def generate_random_bipartite_graph(degrees_virus, degrees_human):
        # Create stubs (nodes with repeated entries based on their degree)
        stubs_virus = []
        for idx, degree in enumerate(degrees_virus):
            stubs_virus.extend([idx] * degree)

        stubs_human = []
        for idx, degree in enumerate(degrees_human):
            stubs_human.extend([idx + len(degrees_virus)] * degree)

        # Shuffle stubs to create random edges
        np.random.shuffle(stubs_virus)
        np.random.shuffle(stubs_human)

        # Create edges while ensuring no duplicates
        edges = set()  # Use a set to ensure uniqueness
        for v, h in zip(stubs_virus, stubs_human):
            edges.add((v, h))  # Add edge as a tuple (source, target)

        # Convert edges back to a list for graph construction
        random_graph = Graph(edges=list(edges), directed=False)

        # Assign vertex types for bipartition
        random_graph.vs["type"] = [1] * len(degrees_virus) + [0] * len(degrees_human)

        return random_graph

    # Generate and save random bipartite graphs
    for i in range(num_random_graphs):
        # Generate random graph
        random_graph = generate_random_bipartite_graph(degrees_virus_families, degrees_human_domains)

        # Assign original node names to the random graph
        random_graph.vs["name"] = list(partition_virus_families) + list(partition_human_domains)

        # Shuffle weights and assign to the random graph
        shuffled_weights = np.random.permutation(original_weights)
        for idx, edge in enumerate(random_graph.es):
            edge["weight"] = shuffled_weights[idx]

        # Check that the random graph is simple
        assert random_graph.is_simple(), f"Random graph {i+1} contains multi-edges!"

        # Get edges in the random graph
        random_edges = set(
            (random_graph.vs[edge.source]["name"], random_graph.vs[edge.target]["name"], edge["weight"])
            for edge in random_graph.es
        )

        # Count identical edges
        identical_edges = original_edges.intersection(random_edges)
        for edge in identical_edges:
            edge_occurrence_counter[edge] += 1

        # Check if the random graph is bipartite
        is_bipartite = random_graph.is_bipartite()
        random_partition_virus_families = {v["name"] for v in random_graph.vs if v["type"] == 1}
        random_partition_human_domains = {v["name"] for v in random_graph.vs if v["type"] == 0}

        if not is_bipartite:
            print(f"Random graph {i+1} is NOT bipartite!")
        elif (random_partition_virus_families != partition_virus_families or
              random_partition_human_domains != partition_human_domains):
            print(f"Random graph {i+1} has incorrect bipartition sets!")
        else:
            print(f"Random graph {i+1} is valid and bipartite.")

        # Save the random graph to a .txt file
        output_file = os.path.join(output_folder, f"{i+1}.txt")
        random_graph.write_ncol(output_file)

    print(f"All {num_random_graphs} random bipartite networks have been saved in the '{output_folder}' folder.")
    print(f"Edge occurrence counts across all random graphs:\n")
    for edge, count in edge_occurrence_counter.items():
        print(f"Edge {edge}: {count} occurrences")

# Define input and output paths
input_file_1 = "local_code\\graph_v2_virusFam_and_humanDomain\\graph_data\\HumanDomain_virusFamily_laplacian.txt"
output_folder_1 = "local_code\\graph_bipartite\\data\\random_graphs_laplacian"

# Call the function
generate_and_save_random_bipartite_graphs(input_file=input_file_1, 
                                          output_folder=output_folder_1, 
                                          num_random_graphs=1000, 
                                          seed=2024)


input_file_2 = "local_code\\graph_v2_virusFam_and_humanDomain\\graph_data\\HumanDomain_virusFamily_matrix_norm_0_1.txt"
output_folder_2 = "local_code\\graph_bipartite\\data\\random_graphs_norm"
# Call the function
generate_and_save_random_bipartite_graphs(input_file=input_file_2, 
                                          output_folder=output_folder_2, 
                                          num_random_graphs=1000, 
                                          seed=2024)
