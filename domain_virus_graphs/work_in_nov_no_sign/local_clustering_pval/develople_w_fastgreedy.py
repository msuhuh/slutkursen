import igraph as ig

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

# Function to load virus families from the classification file
def load_virus_families(classification_file):
    with open(classification_file, "r") as file:
        # Split the line into a list of virus family names
        virus_families = set(file.readline().strip().split(","))
    return virus_families

# Function to detect communities and print cluster statistics
def detect_communities_and_print_statistics(graph_path, classification_file):
    # Load the graph
    g = ig.Graph.Read_Ncol(graph_path, weights=True, directed=False)
    
    # Normalize weights
    normalize_weights(g)
    
    # Load virus family classifications
    virus_families = load_virus_families(classification_file)
    
    # Apply the fast_greedy community detection algorithm
    dendrogram = g.community_fastgreedy(weights=g.es["normalized_weight"])
    communities = dendrogram.as_clustering()
    
    # Print statistics for each community
    for i, community in enumerate(communities):
        community_names = [g.vs[node]["name"] for node in community]
        
        # Classify nodes
        virus_count = sum(1 for node in community_names if node in virus_families)
        domain_count = len(community_names) - virus_count  # Remaining nodes are domains
        
        print(f"Community {i + 1}:")
        print(f"  Virus Families: {virus_count}")
        print(f"  Human Protein Domains: {domain_count}")

# Paths to the graph and classification files
graph_path1 = "local_code/graph_v2_virusFam_and_humanDomain/graph_data/HumanDomain_virusFamily_laplacian.txt"
graph_path2 = "local_code/graph_v2_virusFam_and_humanDomain/graph_data/HumanDomain_virusFamily_matrix_norm_0_1.txt"
classification_file = "local_code\\data\\unique_virus_fam_names.txt"

graphs = [graph_path1, graph_path2]

# Detect communities and print statistics
for graph_path in graphs:
    print(f"Processing graph: {graph_path}")
    detect_communities_and_print_statistics(graph_path, classification_file)
    print("\n")
