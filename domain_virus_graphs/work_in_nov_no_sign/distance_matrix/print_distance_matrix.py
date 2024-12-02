import igraph as ig
import numpy as np

# Load the original graph
def load_graph(file_path):
    # Read graph from NCOL format
    g = ig.Graph.Read_Ncol(file_path, weights=True, directed=False)
    normalize_weights(g)  # Normalize edge weights and compute distances
    return g

# Normalize edge weights and compute "distance_weight"
def normalize_weights(g):
    weights = g.es["weight"]
    min_weight = min(weights)
    max_weight = max(weights)
    if max_weight - min_weight == 0:
        return
    normalized = [(w - min_weight) / (max_weight - min_weight) for w in weights]
    g.es["normalized_weight"] = normalized
    epsilon = 1e-6
    g.es["distance_weight"] = [(1 - w) + epsilon for w in normalized]

# Compute and return the shortest path distance matrix
def compute_distance_matrix(g):
    # Compute shortest path distances
    distance_matrix = np.array(g.distances(weights="distance_weight"))

    # Handle disconnected nodes (assign a large distance or another placeholder)
    disconnected_value = np.max(distance_matrix[np.isfinite(distance_matrix)]) + 1
    distance_matrix[distance_matrix == np.inf] = disconnected_value

    return distance_matrix

# Print the distance matrix
def print_distance_matrix(matrix):
    print("Shortest Path Distance Matrix:")
    for i in range(matrix.shape[0]):

        if i in [100]:
            print(matrix[i])
    print(matrix.shape)

# Main function
def main():
    # Specify the path to the graph file
    graph_file1 = "local_code\graph_v2_virusFam_and_humanDomain\graph_data\HumanDomain_virusFamily_laplacian.txt"  # Replace with your file path
    graph_file2 = "local_code\graph_v2_virusFam_and_humanDomain\graph_data\HumanDomain_virusFamily_matrix_norm_0_1.txt"

    graphs = [graph_file1, graph_file2]
    
    for graph_file in graphs:
        # Load the graph
        graph = load_graph(graph_file)

        # Compute the distance matrix
        distance_matrix = compute_distance_matrix(graph)

        # Print the distance matrix
        print_distance_matrix(distance_matrix)

# Run the script
if __name__ == "__main__":
    main()
