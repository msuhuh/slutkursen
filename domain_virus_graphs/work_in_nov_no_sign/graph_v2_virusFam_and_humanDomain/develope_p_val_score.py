import igraph as ig
import csv
import json
import os

from sklearn.metrics import pairwise_distances
import numpy as np

np.random.seed(2024)

# Load virus family names
with open(r"local_code/create_graph_v1/data_graphs/unique_virus_fam_names.txt", "r") as file:
    virus_names = set(file.read().strip().split(","))  # Virus families are comma-separated

# Load RNA/DNA classification
with open(r"local_code/data/virus_fam_2_nucleic_acid.json", "r") as file:
    virus_type_map = json.load(file)  # Example format: {"VirusFamily1": "RNA", "VirusFamily2": "DNA"}
    



def calculate_silhouette_score(g, communities):
    # Create a distance matrix using distance weights
    distances = np.zeros((len(g.vs), len(g.vs)))
    for edge in g.es:
        source = edge.source
        target = edge.target
        distances[source, target] = edge["distance_weight"]
        distances[target, source] = edge["distance_weight"]  # Symmetric matrix

    # Compute silhouette score for each node
    node_silhouette_scores = []
    for i, community in enumerate(communities):
        for node in community:
            # Nodes in the same cluster
            same_cluster = [n for n in community if n != node]
            a_i = np.mean([distances[node, n] for n in same_cluster]) if same_cluster else 0
            
            # Nodes in the nearest cluster
            nearest_distances = []
            for j, other_community in enumerate(communities):
                if i != j:  # Skip the current cluster
                    other_distances = [distances[node, n] for n in other_community]
                    nearest_distances.append(np.mean(other_distances) if other_distances else float('inf'))
            b_i = min(nearest_distances) if nearest_distances else 0

            # Calculate silhouette score for the node
            s_i = (b_i - a_i) / max(a_i, b_i) if max(a_i, b_i) > 0 else 0
            node_silhouette_scores.append(s_i)

    # Compute overall and per-cluster silhouette scores
    overall_silhouette_score = np.mean(node_silhouette_scores)
    cluster_silhouette_scores = []
    for i, community in enumerate(communities):
        cluster_scores = [
            node_silhouette_scores[node]
            for node in community
        ]
        cluster_silhouette_scores.append(np.mean(cluster_scores) if cluster_scores else 0)

    return overall_silhouette_score, cluster_silhouette_scores


# Normalize weights to [0, 1]
def normalize_weights(g):
    weights = g.es["weight"]
    min_weight = min(weights)
    max_weight = max(weights)
    if max_weight - min_weight == 0:
        print("Warning: All weights are identical; skipping normalization.")
        return
    normalized = [(w - min_weight) / (max_weight - min_weight) for w in weights]
    g.es["normalized_weight"] = normalized
    epsilon = 1e-6
    g.es["distance_weight"] = [(1 - w) + epsilon for w in normalized]  # Inverted for distance

# Define a function to analyze clusters and save to CSV
def analyze_and_save_clusters(g, communities, algorithm_name, params, output_dir):
    cluster_data = []  # To store cluster information for CSV

    for cluster_id, cluster in enumerate(communities):
        virus_count = 0
        human_count = 0
        rna_count = 0
        dna_count = 0
        virus_in_cluster = []
        human_in_cluster = []

        for node in cluster:
            node_name = g.vs[node]["name"]
            if node_name in virus_names:
                virus_count += 1
                virus_in_cluster.append(node_name)
                if virus_type_map.get(node_name, None) == "RNA":
                    rna_count += 1
                elif virus_type_map.get(node_name, None) == "DNA":
                    dna_count += 1
            else:
                human_count += 1
                human_in_cluster.append(node_name)

        # Calculate DNA ratio
        dna_ratio = dna_count / (rna_count + dna_count) if (rna_count + dna_count) > 0 else 0

        # Save data for this cluster
        cluster_data.append({
            "Cluster ID": cluster_id,
            "N_RNA": rna_count,
            "N_DNA": dna_count,
            "DNA_Ratio": dna_ratio,
            "Virus_Families": ", ".join(virus_in_cluster),
            "Domains": ", ".join(human_in_cluster),
        })

    # Save all clusters to a CSV file
    csv_file = os.path.join(output_dir, f"{algorithm_name}_clusters_{params}.csv")
    os.makedirs(output_dir, exist_ok=True)
    with open(csv_file, "w", newline="") as csvfile:
        fieldnames = ["Cluster ID", "N_RNA", "N_DNA", "DNA_Ratio", "Virus_Families", "Domains"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(cluster_data)

    print(f"Cluster data saved to {csv_file}")

# Apply multiple clustering algorithms
def test_clustering_algorithms(g_path):
    g = ig.Graph.Read_Ncol(g_path, weights=True, directed=False)
    normalize_weights(g)  # Normalize weights

    algorithms = {
        "Fast Greedy": lambda: g.community_fastgreedy(weights=g.es["distance_weight"]).as_clustering(),
    }

    # Extract network name
    network_name = os.path.splitext(os.path.basename(g_path))[0]

    for algorithm_name, algorithm_fn in algorithms.items():
        try:
            communities = algorithm_fn()
            if communities is not None:
                output_dir = os.path.join(
                    "local_code/graph_v2_virusFam_and_humanDomain/results_brute_force",
                    algorithm_name.lower().replace(" ", "_")
                )
                params = f"{network_name}_distance_weight"
                analyze_and_save_clusters(g, communities, algorithm_name, params, output_dir)

                # Calculate silhouette score
                overall_score, cluster_scores = calculate_silhouette_score(g, communities)
                print(f"Silhouette Score for {algorithm_name} (Overall): {overall_score:.4f}")
                print(f"Silhouette Scores per Cluster: {cluster_scores}")
        except Exception as e:
            print(f"Error with {algorithm_name} and weights=distance_weight: {e}")

# List of networks to test
networks = [
    "local_code/graph_v2_virusFam_and_humanDomain/graph_data/HumanDomain_virusFamily_laplacian.txt",
    "local_code\\graph_v2_virusFam_and_humanDomain\\random_graphs\\laplacian_norm\\1.txt",
    "local_code\\graph_v2_virusFam_and_humanDomain\\random_graphs\\laplacian_norm\\2.txt",
    "local_code\\graph_v2_virusFam_and_humanDomain\\random_graphs\\laplacian_norm\\3.txt",

]

# Apply clustering algorithms to each network
for graph_path in networks:
    print(f"\nTesting graph: {graph_path}")
    test_clustering_algorithms(graph_path)
