import os
import numpy as np
from igraph import Graph

import win32file
win32file._setmaxstdio(2048)

def generate_and_save_random_graphs(input_file, output_folder, num_random_graphs=500, seed=None):
    """
    Generate random graphs while preserving degree sequence and original node names.
    
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
    
    # Extract degree sequence, weights, and original edges with node names
    degree_sequence = original_graph.degree()
    original_weights = [edge["weight"] for edge in original_graph.es]
    original_edges = set((original_graph.vs[edge.source]["name"], 
                          original_graph.vs[edge.target]["name"], 
                          edge["weight"]) for edge in original_graph.es)
    total_original_edges = len(original_edges)

    print(f"Loaded graph with {len(original_graph.vs)} nodes and {total_original_edges} edges.")

    # Function to generate a single random graph
    def generate_random_graph(degree_sequence, weights):
        random_graph = Graph.Degree_Sequence(degree_sequence, method="vl")
        np.random.shuffle(weights)
        random_graph.es["weight"] = weights[:len(random_graph.es)]
        return random_graph

    # Generate and save random graphs
    for i in range(num_random_graphs):
        # Generate random graph
        random_graph = generate_random_graph(degree_sequence, original_weights)
        
        # Assign original node names to the random graph
        random_graph.vs["name"] = original_graph.vs["name"]
        
        # Compare edges including weights
        random_edges = set((random_graph.vs[edge.source]["name"], 
                            random_graph.vs[edge.target]["name"], 
                            edge["weight"]) for edge in random_graph.es)
        identical_edges = len(original_edges.intersection(random_edges))
        
        # Calculate percentage of identical edges
        percentage_identical = (identical_edges / total_original_edges) * 100
        
        # Print the results
        print(f"Random Network {i+1}: {identical_edges} identical edges ({percentage_identical:.2f}%) to the original graph.")
        
        # Save the random graph to a .txt file
        output_file = os.path.join(output_folder, f"{i+1}.txt")
        random_graph.write_ncol(output_file)
    
    print(f"All {num_random_graphs} random networks have been saved in the '{output_folder}' folder.")


# Define input and output paths
input_file_1 = "local_code\\graph_v2_virusFam_and_humanDomain\\graph_data\\HumanDomain_virusFamily_laplacian.txt"
output_folder_1 = "local_code\\graph_v2_virusFam_and_humanDomain\\random_graphs\\laplacian_norm"

input_file_2 = "local_code\\graph_v2_virusFam_and_humanDomain\\graph_data\\HumanDomain_virusFamily_matrix_norm_0_1.txt"
output_folder_2 = "local_code\\graph_v2_virusFam_and_humanDomain\\random_graphs\\norm"

# Call the function
generate_and_save_random_graphs(input_file=input_file_1, 
                                output_folder=output_folder_1, 
                                num_random_graphs=2000, 
                                seed=2024)

#generate_and_save_random_graphs(input_file=input_file_2, 
#                                output_folder=output_folder_2, 
#                                num_random_graphs=2000, 
#                                seed=2024)