from igraph import Graph
import numpy as np

np.random.seed(2024)

def check_bipartition(path_to_txt):
    # Load the graph from the edge list
    graph = Graph.Read_Ncol(path_to_txt, weights=True, directed=False)

    # Load the virus names from a text file
    with open(r"local_code\create_graph_v1\data_graphs\unique_virus_fam_names.txt", "r") as file:
        virus_names = set(file.read().strip().split(","))  # Read and split by commas

    # Create a "type" attribute based on whether the vertex is a virus family or not
    graph.vs["type"] = [1 if v["name"] in virus_names else 0 for v in graph.vs]

    # Check if the graph is bipartite
    is_bipartite = graph.is_bipartite()
    if not is_bipartite:
        print(f"Graph loaded from {path_to_txt} is NOT bipartite.")
        return

    # Separate nodes into partitions
    partition_0 = [v["name"] for v in graph.vs if v["type"] == 0]
    partition_1 = [v["name"] for v in graph.vs if v["type"] == 1]

    # Validate bipartition
    viruses_in_partition_0 = set(partition_0).intersection(virus_names)
    viruses_in_partition_1 = set(partition_1).intersection(virus_names)

    print(f"Graph loaded from {path_to_txt} is bipartite.")
    if viruses_in_partition_0 and viruses_in_partition_1:
        print("ERROR: Virus families are split across both partitions.")
    else:
        correct_partition = "partition_0" if viruses_in_partition_0 else "partition_1"
        print(f"Virus families are correctly in {correct_partition}.")

    human_domains_in_partition_0 = set(partition_0) - viruses_in_partition_0
    human_domains_in_partition_1 = set(partition_1) - viruses_in_partition_1

    if human_domains_in_partition_0 and human_domains_in_partition_1:
        print("ERROR: Human domains are split across both partitions.")
    else:
        correct_partition = "partition_0" if human_domains_in_partition_0 else "partition_1"
        print(f"Human domains are correctly in {correct_partition}.")

    # Print nodes in each partition
    print(f"Nodes in partition 0 ({len(partition_0)} nodes): {partition_0}")
    print("")
    print(f"Nodes in partition 1 ({len(partition_1)} nodes): {partition_1}")
    print("")

networks = [
    "local_code\graph_v2_virusFam_and_humanDomain\graph_data\HumanDomain_virusFamily_laplacian.txt",
    #"local_code\graph_v2_virusFam_and_humanDomain\graph_data\HumanDomain_virusFamily_matrix_norm_0_1.txt"
]

for graph_path in networks:
    check_bipartition(graph_path)
