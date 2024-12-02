import igraph as ig
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.manifold import MDS

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

# Create a heatmap of the distance matrix
def plot_distance_heatmap(distance_matrix, title, output_file):
    plt.figure(figsize=(12, 10))
    sns.heatmap(
        distance_matrix,
        cmap="viridis",
        square=True,
        cbar=True,
        xticklabels=False,  # Turn off tick labels for readability
        yticklabels=False,
    )
    plt.title(title)
    plt.xlabel("Nodes")
    plt.ylabel("Nodes")
    plt.tight_layout()
    plt.savefig(output_file)
    plt.show()

# Create an MDS plot
def plot_mds(distance_matrix, title, output_file):
    mds = MDS(n_components=2, dissimilarity="precomputed", random_state=42)
    mds_coords = mds.fit_transform(distance_matrix)
    plt.figure(figsize=(10, 8))
    plt.scatter(mds_coords[:, 0], mds_coords[:, 1], s=50, cmap="viridis")
    plt.title(title)
    plt.xlabel("MDS Dimension 1")
    plt.ylabel("MDS Dimension 2")
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

        # Plot distance heatmap
        heatmap_title = f"Distance Heatmap for {graph_file.split('/')[-1]}"
        heatmap_output = f"{graph_file.split('/')[-1].split('.')[0]}_distance_heatmap.png"
        plot_distance_heatmap(distance_matrix, heatmap_title, heatmap_output)

        # Plot MDS visualization
        mds_title = f"MDS Plot for {graph_file.split('/')[-1]}"
        mds_output = f"{graph_file.split('/')[-1].split('.')[0]}_mds_plot.png"
        plot_mds(distance_matrix, mds_title, mds_output)

if __name__ == "__main__":
    main()
