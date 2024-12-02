import os
import igraph as ig
from collections import defaultdict
import win32file

# Increase the file limit for Windows
win32file._setmaxstdio(4096)

# Function to normalize weights
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

# Function to get node pairs in the same cluster using leading eigenvector clustering
def get_cluster_pairs(graph, clusters=None):
    communities = graph.community_leading_eigenvector(clusters=clusters, weights=graph.es["normalized_weight"])
    cluster_pairs = set()
    for community in communities:
        for i in range(len(community)):
            for j in range(i + 1, len(community)):
                node1 = graph.vs[community[i]]["name"]
                node2 = graph.vs[community[j]]["name"]
                cluster_pairs.add(tuple(sorted((node1, node2))))
    return cluster_pairs

# Function to compare the original graph with random graphs
def compare_original_with_random_clusters(original_path, random_folder, output_file_base):
    # Load the original graph
    original_g = ig.Graph.Read_Ncol(original_path, weights=True, directed=False)
    normalize_weights(original_g)

    # Test for various cluster arguments
    for cluster_arg in [None]: #+ list(range(2, 11)):
        print(f"Processing cluster argument: {cluster_arg}")
        output_file = f"{output_file_base}_clusters_{cluster_arg}.txt"

        # Get node names and cluster pairs from the original graph
        original_pairs = get_cluster_pairs(original_g, clusters=cluster_arg)
        print(f"Original pairs (clusters={cluster_arg}): {len(original_pairs)} pairs")

        # Initialize pair counts for only original graph pairs
        pair_counts = defaultdict(int)

        # Load random graphs and track pair occurrences
        random_files = [
            os.path.join(random_folder, f)
            for f in os.listdir(random_folder)
            if f.endswith(".txt")
        ]

        for i, random_path in enumerate(random_files):
            if i % 100 == 0:
                print(f"Processing random graph {i + 1}/{len(random_files)}...")
            if i == 2000:
                break
            random_g = ig.Graph.Read_Ncol(random_path, weights=True, directed=False)
            normalize_weights(random_g)

            random_pairs = get_cluster_pairs(random_g, clusters=cluster_arg)

            # Update counts for pairs in random clusters that exist in the original pairs
            for pair in random_pairs:
                if pair in original_pairs:
                    pair_counts[pair] += 1

        # Calculate p-values and filter unique pairs
        unique_pairs = {
            pair: count
            for pair, count in pair_counts.items()
            if count < 100
        }

        # Write results to file
        with open(output_file, "w") as file:
            file.write(f"Unique pairs with p-value < 0.05 (less than 100 occurrences in random graphs):\n")
            for pair, count in unique_pairs.items():
                file.write(f"Pair {pair}: Random occurrences = {count}\n")
        print(f"Results for clusters={cluster_arg} written to {output_file}")

# Main script
original_network1 = "local_code/graph_v2_virusFam_and_humanDomain/graph_data/HumanDomain_virusFamily_laplacian.txt"
#random_networks_folder1 = "local_code/graph_bipartite/data/random_graphs_laplacian"
random_networks_folder1 = "local_code/graph_v2_virusFam_and_humanDomain/random_graphs/laplacian_norm"
output_file1_base = "local_code/local_clustering_pval/results/leading_eigenvector/laplacian_unique_pairs"

original_network2 = "local_code/graph_v2_virusFam_and_humanDomain/graph_data/HumanDomain_virusFamily_matrix_norm_0_1.txt"
#random_networks_folder2 = "local_code/graph_bipartite/data/random_graphs_norm"
random_networks_folder2 = "local_code/graph_v2_virusFam_and_humanDomain/random_graphs/norm"
output_file2_base = "local_code/local_clustering_pval/results/leading_eigenvector/normalized_unique_pairs"

os.makedirs(output_file1_base, exist_ok=True)
os.makedirs(output_file2_base, exist_ok=True)

# Compare original with random graphs for different cluster arguments
compare_original_with_random_clusters(
    original_path=original_network1,
    random_folder=random_networks_folder1,
    output_file_base=output_file1_base
)

compare_original_with_random_clusters(
    original_path=original_network2,
    random_folder=random_networks_folder2,
    output_file_base=output_file2_base
)
