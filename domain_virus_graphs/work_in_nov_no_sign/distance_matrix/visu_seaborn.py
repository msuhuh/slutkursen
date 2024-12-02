import igraph as ig
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

# Load the original graph
def load_graph(file_path):
    g = ig.Graph.Read_Ncol(file_path, weights=True, directed=False)
    normalize_weights(g)
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
    distance_matrix = np.array(g.distances(weights="distance_weight"))
    disconnected_value = np.max(distance_matrix[np.isfinite(distance_matrix)]) + 1
    distance_matrix[distance_matrix == np.inf] = disconnected_value
    return distance_matrix

# Create a correlation plot
def plot_correlation(distance_matrix, title, output_file):
    df = pd.DataFrame(distance_matrix)
    correlation_matrix = df.corr()

    # Fix for potential heatmap issues
    plt.figure(figsize=(12, 10))  # Adjusted size for better readability
    ax = sns.heatmap(
        correlation_matrix,
        cmap="coolwarm",
        annot=False,
        fmt=".2f",
        xticklabels=False,  # Prevent overlapping for 145x145
        yticklabels=False,
        cbar=True,
    )
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.show()

# Main function
def main():
    graph_files = [
        "local_code/graph_v2_virusFam_and_humanDomain/graph_data/HumanDomain_virusFamily_laplacian.txt",
        "local_code/graph_v2_virusFam_and_humanDomain/graph_data/HumanDomain_virusFamily_matrix_norm_0_1.txt",
    ]
    for graph_file in graph_files:
        graph = load_graph(graph_file)
        distance_matrix = compute_distance_matrix(graph)
        title = f"Correlation Heatmap for {graph_file.split('/')[-1]}"
        output_file = f"{graph_file.split('/')[-1].split('.')[0]}_correlation_heatmap.png"
        plot_correlation(distance_matrix, title, output_file)

if __name__ == "__main__":
    main()
