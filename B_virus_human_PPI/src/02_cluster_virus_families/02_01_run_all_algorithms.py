import os
import igraph as ig
import random
import win32file
win32file._setmaxstdio(8192)

def compare_original_clusters_w_random(cluster_counts, random_lists_w_viruss_fam):
    """
    Compare original clusters (using sets) with random lists of virus families.
    
    Parameters:
    - cluster_counts: List of tuples, where each tuple contains a cluster (set of virus families) and a count (int).
    - random_lists_w_viruss_fam: List of lists, where each list contains virus families from random graphs.
    
    Returns:
    - Updated cluster_counts with incremented counts for clusters matched in random lists.
    """
    for random_list in random_lists_w_viruss_fam:
        random_set = set(random_list)
        
        # Check each cluster against the current random set
        for idx, (cluster_set, count) in enumerate(cluster_counts):
            if cluster_set.issubset(random_set):  # Subset check
                cluster_counts[idx] = (cluster_set, count + 1)
    
    return cluster_counts


def get_whole_clusters_only_families(graph, algorithm, virus_names):
    """
    Extract clusters of viruses from a graph.
    
    Parameters:
    - graph: Graph object with nodes.
    - algorithm: Community detection algorithm function.
    - virus_names: Set of virus names to filter.
    
    Returns:
    - List of lists: Each sublist contains virus names from one cluster.
    """
    communities = algorithm(graph)
    clusters = []  # To store clusters of viruses
    
    for community in communities:
        # Extract virus names in the current community
        virus_cluster = [
            graph.vs[node]["name"] for node in community if graph.vs[node]["name"] in virus_names
        ]
        #print(virus_cluster)
        # Only include clusters that have at least 2 viruses
        if len(virus_cluster) > 1 :
            clusters.append(virus_cluster)
    
    return clusters

def compare_original_with_random_clusters_n_times(original_path, pre_loaded_rand_graphs, output_file_no_suffix, txt_virus_families, iterations, algorithm, seed=None):
    """Compare node pairs clustered together in the original graph with random graphs."""
    
    # Set seed for reproducibility
    if seed is not None:
        random.seed(seed)

    # Define threshold and number of random files for comparison
    threshold = 0.01
    n_random_files = len(pre_loaded_rand_graphs)

    # Check paths
    if not os.path.exists(original_path):
        raise FileNotFoundError(f"Original graph file not found: {original_path}")
    if not os.path.exists(txt_virus_families):
        raise FileNotFoundError(f"Virus families txt file not found: {txt_virus_families}")

    # Load unique virus and bait names
    with open(txt_virus_families, "r") as f:
        virus_names = set(f.read().strip().split(","))

    # Load the original graph and its cluster pairs
    original_g = ig.Graph.Read_Ncol(original_path, weights=True, directed=False)
    original_lists_w_viruss_fam = get_whole_clusters_only_families(original_g, algorithm, virus_names)

    assert len(original_lists_w_viruss_fam) != 0, "No cluster fams were found"
    print(f"{len(original_lists_w_viruss_fam)} clusters with at least 2 virues were found")

    # Repeat the N random graphs X times to ensure that there is not noise in the result.
    for iter_number in range(iterations):
        
        print("Starting iter: ", iter_number)

        # Initialize cluster counts with 0 using sets for faster comparisons
        cluster_counts = [(set(cluster), 0) for cluster in original_lists_w_viruss_fam]
        # Compare original graph agaisnt random graphs
        for random_graph in (pre_loaded_rand_graphs):
            random_lists_w_viruss_fam = get_whole_clusters_only_families(random_graph, algorithm, virus_names)
            cluster_counts = compare_original_clusters_w_random(cluster_counts, random_lists_w_viruss_fam)

        
        # Save significant result for current iteration of N random graphs
        output_file = output_file_no_suffix + str(iter_number)+'.txt'

        # Open the output file for writing
        with open(output_file, "w") as file:
            # Write header information
            file.write(f"p value was set to: {threshold}.\n")
            file.write(f"{n_random_files} random graphs were analyzed, threshold became: {int(threshold * n_random_files)}\n")
            file.write("Cluster results (TRUE/FALSE based on threshold, followed by the cluster and count):\n\n")
            
            # Process each cluster
            for cluster, count in cluster_counts:
                # Determine if the count meets the threshold
                meets_threshold = count < int(threshold * n_random_files)
                status = "TRUE" if meets_threshold else "FALSE"
                
                # Write the result to the file
                cluster_names = ", ".join(cluster)  # Convert cluster set to a comma-separated string
                file.write(f"{status}: Cluster = [{cluster_names}], Count = {count}\n")


############## Main script #################
# Set a seed for reproducibility
seed = 2024

######### Params ###############
original_network = "src/00_data_folder/01_ouput_graphs/normalized_graph_ncol.txt"
random_networks_folder = "src/00_data_folder/01_ouput_graphs/random_graphs"
output_root_folder = "src/00_data_folder/02_01_virus_outputs"
txt_virus_families = "src/00_data_folder/00_original_data/unique_virus_fam_names.txt"

# READ ME
# See also p-val threshold in function, easy to change if needed.
# Out commented algorithms were tested but did not yield result.
# Other algorithms outside the outcommented ones were also tested. Mostly, arguments were vaired
# Some algorithms were discareded due to the time they took to run.

################################

iterations_per_algo = 5 #Number of times to repeats experiment

algorithms = [

        # Leiden
        (lambda g: g.community_leiden(objective_function='modularity', weights=g.es["weight"],
                                         resolution=1, beta=0.01, initial_membership=None, n_iterations=-1, node_weights=None),
                                           "leiden_mod_res1_"),
        #(lambda g: g.community_leiden(objective_function='modularity', weights=g.es["weight"],
        #                                 resolution=0.5, beta=0.01, initial_membership=None, n_iterations=-1, node_weights=None),
        #                                   "leiden_mod_res0.5_"),
        (lambda g: g.community_leiden(objective_function='CPM', weights=g.es["weight"],
                                         resolution=0.01, beta=0.01, initial_membership=None, n_iterations=-1, node_weights=None),
                                           "leiden_CPM_res0.01_"),
        #(lambda g: g.community_leiden(objective_function='CPM', weights=g.es["weight"],
        #                                 resolution=0.001, beta=0.01, initial_membership=None, n_iterations=-1, node_weights=None),
        #                                   "leiden_CPM_res0.001_"),
                                
        ######################

        (lambda g: g.community_multilevel(weights=g.es["weight"]), "louvain_"),
        (lambda g: g.community_fastgreedy(weights=g.es["weight"]).as_clustering(), "fastgreedy_"),
        (lambda g: g.community_leading_eigenvector(weights=g.es["weight"]), "eigenvector_"),
        (lambda g: g.community_infomap(edge_weights=g.es["weight"], trials=10), "infomap10_"),
        #(lambda g: g.community_label_propagation(weights=g.es["weight"]), "labelprop_"),
        (lambda g: g.community_walktrap(steps=2, weights=g.es["weight"]).as_clustering(), "walktrap2_"),
        #(lambda g: g.community_walktrap(steps=3, weights=g.es["weight"]).as_clustering(), "walktrap3_"), 
        #(lambda g: g.community_edge_betweenness(clusters=2, directed=False, weights=g.es["weight"]).as_clustering(), "edges2_"),

        ]

###############################

# Load the pregenerated graphs
if not os.path.exists(random_networks_folder):
        raise FileNotFoundError(f"Random graphs folder not found: {random_networks_folder}")

# Load all random graphs, only once. Windows issue with accessing many files:
# This is also the steps that takes the most time.
random_graphs_paths = [
    os.path.join(random_networks_folder, f) for f in os.listdir(random_networks_folder) if f.endswith(".txt")
]

# Shorten the amount of graphs for testing (if needed)
#random_graphs_paths = random_graphs_paths[:100]


list_of_random_graphs = []
for i, random_path in enumerate(random_graphs_paths):
    if i % 500 == 0:
        print(i, "in loading_graphs")
    list_of_random_graphs.append(ig.Graph.Read_Ncol(random_path, weights=True, directed=False))

for algorithm, algorithm_name in algorithms:
    print("Starting alogrithm: ", algorithm_name)
    subfolder = os.path.join(output_root_folder, algorithm_name)
    os.makedirs(subfolder, exist_ok=True)
    output_file_no_suffix = os.path.join(subfolder, algorithm_name)

    # Compare original with random graphs
    compare_original_with_random_clusters_n_times(
        original_path=original_network,
        pre_loaded_rand_graphs=list_of_random_graphs,
        output_file_no_suffix=output_file_no_suffix,
        txt_virus_families=txt_virus_families,
        iterations = iterations_per_algo,
        algorithm=algorithm,
        seed=seed
    )
