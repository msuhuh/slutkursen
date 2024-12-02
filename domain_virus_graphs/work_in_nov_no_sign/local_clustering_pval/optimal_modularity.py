import os
import igraph as ig
from collections import defaultdict

import win32file
# Open N files
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

# Function to get node pairs in the same cluster
def get_cluster_pairs(graph, algorithm):
    dendrogram = algorithm(graph)
    communities = dendrogram.as_clustering()
    cluster_pairs = set()
    for community in communities:
        for i in range(len(community)):
            for j in range(i + 1, len(community)):
                node1 = graph.vs[community[i]]["name"]
                node2 = graph.vs[community[j]]["name"]
                cluster_pairs.add(tuple(sorted((node1, node2))))
    return cluster_pairs

# Function to compare the original graph with random graphs
def compare_original_with_random_clusters(original_path, random_folder, output_file, algorithm):
    # Load the original graph
    original_g = ig.Graph.Read_Ncol(original_path, weights=True, directed=False)
    normalize_weights(original_g)
    
    # Get node names and cluster pairs from the original graph
    original_pairs = get_cluster_pairs(original_g, algorithm)

    print(original_pairs)
    
    # Initialize pair counts for only original graph pairs
    pair_counts = defaultdict(int)

    # Load random graphs and track pair occurrences
    random_files = [
        os.path.join(random_folder, f)
        for f in os.listdir(random_folder)
        if f.endswith(".txt")
    ]

    for i, random_path in enumerate(random_files):
        print(i)
        if i % 100 == 0:
            print(f"Processing random graph {i + 1}/{len(random_files)}...")
        if i == 1000:
            break
        random_g = ig.Graph.Read_Ncol(random_path, weights=True, directed=False)
        normalize_weights(random_g)

        random_pairs = get_cluster_pairs(random_g, algorithm)

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
        file.write(f"Unique pairs with p-value < 0.01 (less than 10 occurrences in random graphs):\n")
        for pair, count in unique_pairs.items():
            file.write(f"Pair {pair}: Random occurrences = {count}\n")
    print(f"Results written to {output_file}")

# Main script
original_network1 = "local_code/graph_v2_virusFam_and_humanDomain/graph_data/HumanDomain_virusFamily_laplacian.txt"
random_networks_folder1 = "local_code/graph_bipartite/data/random_graphs_laplacian"
output_file1 = "local_code\local_clustering_pval/laplacian_unique_pairs.txt"

original_network2 = "local_code/graph_v2_virusFam_and_humanDomain/graph_data/HumanDomain_virusFamily_matrix_norm_0_1.txt"
random_networks_folder2 = "local_code/graph_bipartite/data/random_graphs_norm"
output_file2 = "local_code\local_clustering_pval/optimal_modularity_normalized.txt"

# Define the clustering algorithm
algorithm = lambda g: g.community_optimal_modularity(weights=g.es["normalized_weight"])

# Compare original with random graphs
compare_original_with_random_clusters(
    original_path=original_network1,
    random_folder=random_networks_folder1,
    output_file=output_file1,
    algorithm=algorithm
)

compare_original_with_random_clusters(
    original_path=original_network2,
    random_folder=random_networks_folder2,
    output_file=output_file2,
    algorithm=algorithm
)
