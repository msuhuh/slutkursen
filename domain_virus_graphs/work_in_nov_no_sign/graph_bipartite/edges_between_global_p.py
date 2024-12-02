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
win32file._setmaxstdio(4096)

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
    distances = np.zeros((len(g.vs), len(g.vs)))
    for edge in g.es:
        distances[edge.source, edge.target] = edge["distance_weight"]
        distances[edge.target, edge.source] = edge["distance_weight"]

    labels = np.zeros(len(g.vs), dtype=int)
    for i, community in enumerate(communities):
        labels[community] = i

    return silhouette_score(distances, labels, metric="precomputed")


def calculate_davies_bouldin_index_from_graph(g, communities):
    distances = np.zeros((len(g.vs), len(g.vs)))
    for edge in g.es:
        distances[edge.source, edge.target] = edge["distance_weight"]
        distances[edge.target, edge.source] = edge["distance_weight"]

    labels = np.zeros(len(g.vs), dtype=int)
    for i, community in enumerate(communities):
        labels[community] = i

    return davies_bouldin_score(distances, labels)


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


def compute_metrics(g, algorithm):
    communities = algorithm(g)
    #print(len(communities))
    # Print the size of each community
    for i, community in enumerate(communities):
        print(f"Community {i + 1}: {len(community)} nodes")
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
def compare_original_with_random(original_path, random_folder, output_file, algorithm):
    # Load original network
    original_g = ig.Graph.Read_Ncol(original_path, weights=True, directed=False)
    normalize_weights(original_g)
    original_scores = compute_metrics(original_g, algorithm)

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
        if i % 200 == 0:
            print(i)
        random_g = ig.Graph.Read_Ncol(random_path, weights=True, directed=False)
        normalize_weights(random_g)
        random_scores = compute_metrics(random_g, algorithm)
        if i == 1000:
            break

        for metric in original_scores:
            if metric == "davies_bouldin_index":  # Lower is better
                if random_scores[metric] < original_scores[metric]:
                    better_counts[metric] += 1
            else:  # Higher is better
                if random_scores[metric] > original_scores[metric]:
                    better_counts[metric] += 1

    # Write results to file
    with open(output_file, "w") as file:
        file.write(f"Original Graph Scores:\n")
        for metric, score in original_scores.items():
            file.write(f"{metric.capitalize()}: {score:.4f}\n")

        file.write("\nComparison Results:\n")
        for metric, count in better_counts.items():
            file.write(f"Metric '{metric.capitalize()}' - Random networks scored better {count} times.\n")


# Main script
original_network1 = "local_code/graph_v2_virusFam_and_humanDomain/graph_data/HumanDomain_virusFamily_laplacian.txt"
random_networks_folder1 = "local_code\\graph_bipartite\\data\\random_graphs_laplacian"
output_file1 = "local_code\\graph_bipartite\\results\\edge_betwenness\\laplacian_global.txt"

original_network2 = "local_code\graph_v2_virusFam_and_humanDomain\graph_data\HumanDomain_virusFamily_matrix_norm_0_1.txt"
random_networks_folder2 = "local_code\\graph_bipartite\\data\\random_graphs_norm"
output_file2 = "local_code\\graph_bipartite\\results\\edge_betwenness\\normalized_global.txt"

#compare_original_with_random(
 #   original_path=original_network1,
#    random_folder=random_networks_folder1,
#    output_file=output_file1,
#    algorithm=lambda g: g.community_edge_betweenness(directed=False, weights=g.es["distance_weight"]).as_clustering(n=5),
#)

compare_original_with_random(
    original_path=original_network2,
    random_folder=random_networks_folder2,
    output_file=output_file2,
    algorithm=lambda g: g.community_edge_betweenness(directed=False, weights=g.es["distance_weight"]).as_clustering(n=5),
)

