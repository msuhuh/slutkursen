import igraph as ig
import csv
import json
import os

# Load virus family names
with open(r"local_code/create_graph_v1/data_graphs/unique_virus_fam_names.txt", "r") as file:
    virus_names = set(file.read().strip().split(","))  # Virus families are comma-separated

# Load RNA/DNA classification
with open(r"local_code/data/virus_fam_2_nucleic_acid.json", "r") as file:
    virus_type_map = json.load(file)  # Example format: {"VirusFamily1": "RNA", "VirusFamily2": "DNA"}

# Normalize weights to [0, 1]
def normalize_weights(g):
    weights = g.es["weight"]
    min_weight = min(weights)
    max_weight = max(weights)
    normalized = [(w - min_weight) / (max_weight - min_weight) for w in weights]
    g.es["normalized_weight"] = normalized
    g.es["distance_weight"] = [1 - w for w in normalized]  # Inverted for distance

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
        "Edge Betweenness": lambda: g.community_edge_betweenness(weights="distance_weight").as_clustering(),
        "Fast Greedy": lambda: g.community_fastgreedy(weights="normalized_weight").as_clustering(),
        "Louvain": lambda: g.community_multilevel(weights="normalized_weight"),
        "Spinglass": lambda: g.community_spinglass(weights="distance_weight") if g.is_connected() else None,
        "Walktrap": lambda: g.community_walktrap(weights="normalized_weight").as_clustering(),
        "Leading Eigenvector": lambda: g.community_leading_eigenvector(weights="normalized_weight"),
    }

    # Extract network name
    network_name = os.path.basename(g_path).split(".")[0]

    for algorithm_name, algorithm_fn in algorithms.items():
        print(f"\nTesting algorithm: {algorithm_name}")

        try:
            communities = algorithm_fn()
            if communities is not None:
                output_dir = "local_code/graph_v2_virusFam_and_humanDomain/results_brute_force"
                params = f"{network_name}"
                analyze_and_save_clusters(g, communities, algorithm_name, params, output_dir)
        except Exception as e:
            print(f"Error with {algorithm_name}: {e}")

# List of networks to test
networks = [
    "local_code/graph_v2_virusFam_and_humanDomain/graph_data/hit_matrix.txt",
    "local_code/graph_v2_virusFam_and_humanDomain/graph_data/HumanDomain_virusFamily_laplacian.txt",
    "local_code/graph_v2_virusFam_and_humanDomain/graph_data/HumanDomain_virusFamily_matrix_norm_0_1.txt",
]

# Apply clustering algorithms to each network
for graph_path in networks:
    print(f"\nTesting graph: {graph_path}")
    test_clustering_algorithms(graph_path)
