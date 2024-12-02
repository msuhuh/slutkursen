import igraph as ig
import csv
import json

# Load virus family names
with open(r"local_code/create_graph_v1/data_graphs/unique_virus_fam_names.txt", "r") as file:
    virus_names = set(file.read().strip().split(","))  # Virus families are comma-separated

# Load RNA/DNA classification
with open(r"local_code/data/virus_fam_2_nucleic_acid.json", "r") as file:
    virus_type_map = json.load(file)  # Example format: {"VirusFamily1": "RNA", "VirusFamily2": "DNA"}

# Define a function to analyze clusters and save to CSV
def analyze_and_save_clusters(g, communities, title="Communities"):
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
    csv_file = f"{title.replace(' ', '_')}_clusters.csv"
    with open(csv_file, "w", newline="") as csvfile:
        fieldnames = ["Cluster ID", "N_RNA", "N_DNA", "DNA_Ratio", "Virus_Families", "Domains"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(cluster_data)

    print(f"Cluster data saved to {csv_file}")

# Main function to apply clustering algorithms
def test_clustering_algorithms(g_path):
    g = ig.Graph.Read_Ncol(g_path, weights=True, directed=False)

    # Perform clustering (e.g., Fast Greedy)
    dendrogram = g.community_fastgreedy(weights="weight")
    communities = dendrogram.as_clustering()

    network_name = ""
    if "hit" in g_path:
        network_name += "_non_normalized"
    elif "laplacian" in g_path:
        network_name += "_laplacian_normalized"
    elif "norm" in g_path:
        network_name += "_normalized"


    # Analyze clusters and save results
    analyze_and_save_clusters(g, communities, title="Fast Greedy Clustering" + network_name)

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
