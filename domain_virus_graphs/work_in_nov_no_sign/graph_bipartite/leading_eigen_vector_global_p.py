import igraph as ig
import os
import numpy as np
from sklearn.metrics import silhouette_score, davies_bouldin_score
import networkx as nx
from win32 import win32file
import random

# Set seeds for reproducibility
np.random.seed(2024)
random.seed(2024)
ig.set_random_number_generator(random.Random(2024))

# Open N files
win32file._setmaxstdio(8192)

# Define metrics functions
def calculate_conductance(g, communities):
    nx_graph = nx.Graph()
    nx_graph.add_nodes_from(range(len(g.vs)))
    nx_graph.add_weighted_edges_from((edge.source, edge.target, edge["normalized_weight"]) for edge in g.es)

    conductances = []
    for community in communities:
        community_nodes = set(community)
        internal_edges = sum(
            nx_graph[u][v]["weight"]
            for u in community_nodes
            for v in community_nodes
            if nx_graph.has_edge(u, v)
        )
        external_edges = sum(
            nx_graph[u][v]["weight"]
            for u in community_nodes
            for v in nx_graph.neighbors(u)
            if v not in community_nodes
        )
        total_edges = internal_edges + external_edges
        if total_edges > 0:
            conductance = external_edges / total_edges
            conductances.append(conductance)

    avg_conductance = np.mean(conductances) if conductances else None
    return avg_conductance


def calculate_modularity(g, communities):
    return g.modularity(communities, weights=g.es["normalized_weight"])


def calculate_silhouette_score_from_graph(g, communities):
    # Compute shortest path distances for all pairs of nodes
    shortest_paths = np.array(g.distances(weights="distance_weight"))
    
    # Handle disconnected nodes (assign large distances)
    disconnected_value = np.max(shortest_paths[np.isfinite(shortest_paths)]) + 1
    shortest_paths[shortest_paths == np.inf] = disconnected_value

    # Assign cluster labels
    labels = np.zeros(len(g.vs), dtype=int)
    for i, community in enumerate(communities):
        labels[community] = i

    # Calculate silhouette score using precomputed distances
    return silhouette_score(shortest_paths, labels, metric="precomputed")


def calculate_davies_bouldin_index_from_graph(g, communities):
    # Compute shortest path distances for all pairs of nodes
    shortest_paths = np.array(g.distances(weights="distance_weight"))
    
    # Assign cluster labels
    labels = np.zeros(len(g.vs), dtype=int)
    for i, community in enumerate(communities):
        labels[community] = i

    # Calculate Davies-Bouldin Index using precomputed distances
    return davies_bouldin_score(shortest_paths, labels)


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


def compute_metrics(g, clusters=None):

    # Perform community detection using the leading eigenvector method
    communities = g.community_leading_eigenvector(clusters=clusters, weights=g.es["normalized_weight"])

    modularity = calculate_modularity(g, communities)
    silhouette = calculate_silhouette_score_from_graph(g, communities)
    davies_bouldin_index = calculate_davies_bouldin_index_from_graph(g, communities)
    conductance = calculate_conductance(g, communities)
    return {
        "modularity": modularity,
        "silhouette": silhouette,
        "davies_bouldin_index": davies_bouldin_index,
        "conductance": conductance,
    }


# Load networks and compute comparisons
def compare_original_with_random(original_path, random_folder, output_file_base):
    # Load original network
    original_g = ig.Graph.Read_Ncol(original_path, weights=True, directed=False)
    normalize_weights(original_g)

    # Test for various cluster arguments
    for cluster_arg in [None] + list(range(2, 11)):
        print("Current cluster arg: ", cluster_arg)
        output_file = f"{output_file_base}_clusters_{cluster_arg}.txt"
        original_scores = compute_metrics(original_g, clusters=cluster_arg)

        # Load random networks
        random_files = [
            os.path.join(random_folder, f)
            for f in os.listdir(random_folder)
            if f.endswith(".txt")
        ]
        better_counts = {metric: 0 for metric in original_scores}

        i = 0
        for random_path in random_files:
            i += 1
            if i % 500 == 0:
                print(i)
            random_g = ig.Graph.Read_Ncol(random_path, weights=True, directed=False)
            normalize_weights(random_g)
            random_scores = compute_metrics(random_g, clusters=cluster_arg)

            for metric in original_scores:
                if metric == "davies_bouldin_index":  # Lower is better
                    if random_scores[metric] < original_scores[metric]:
                        better_counts[metric] += 1
                else:  # Higher is better
                    if random_scores[metric] > original_scores[metric]:
                        better_counts[metric] += 1

        # Write results to file
        with open(output_file, "w") as file:
            file.write(f"Original Graph Scores (clusters={cluster_arg}):\n")
            for metric, score in original_scores.items():
                file.write(f"{metric.capitalize()}: {score:.4f}\n")

            file.write("\nComparison Results:\n")
            for metric, count in better_counts.items():
                file.write(f"Metric '{metric.capitalize()}' - Random networks scored better {count} times.\n")


# Main script
original_network1 = "local_code/graph_v2_virusFam_and_humanDomain/graph_data/HumanDomain_virusFamily_laplacian.txt"
random_networks_folder1 = "local_code\\graph_bipartite\\data\\random_graphs_laplacian"
output_file1_base = "local_code\\graph_bipartite\\results\\leading_eigenvector\\laplacian_global"

original_network2 = "local_code\graph_v2_virusFam_and_humanDomain\graph_data\HumanDomain_virusFamily_matrix_norm_0_1.txt"
random_networks_folder2 = "local_code\\graph_bipartite\\data\\random_graphs_norm"
output_file2_base = "local_code\\graph_bipartite\\results\\leading_eigenvector\\normalized_global"

print("hi")
os.makedirs(output_file1_base, exist_ok=True)
os.makedirs(output_file2_base, exist_ok=True)

#compare_original_with_random(
#    original_path=original_network1,
#    random_folder=random_networks_folder1,
#    output_file_base=output_file1_base,
    #algorithm=lambda g: g.community_leading_eigenvector()
#)

compare_original_with_random(
    original_path=original_network2,
    random_folder=random_networks_folder2,
    output_file_base=output_file2_base,
    #algorithm=lambda g: g.community_leading_eigenvector()
)
