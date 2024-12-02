import igraph as ig
from igraph import Graph
import matplotlib.pyplot as plt

# Define a function to visualize communities
def visualize_communities(g, communities, title="Communities"):
    num_communities = len(communities)
    print(f"Number of communities ({title}):", num_communities)

    # Assign colors to vertices based on community membership
    palette = ig.RainbowPalette(n=num_communities)
    g.vs["color"] = [palette.get(communities.membership[i]) for i in range(len(g.vs))]

    fig, ax = plt.subplots()
    ig.plot(
        g,
        palette=palette,
        edge_width=1,
        target=ax,
        vertex_size=20,
        vertex_label=None,
    )

    # Create a custom color legend
    legend_handles = []
    for i in range(num_communities):
        handle = ax.scatter(
            [], [],
            s=100,
            facecolor=palette.get(i),
            edgecolor="k",
            label=f"Community {i}",
        )
        legend_handles.append(handle)
    ax.legend(
        handles=legend_handles,
        title=title,
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
    )
    plt.show()

# Define a function to apply clustering algorithms with varying parameters
def test_clustering_algorithms(g_path):
    g = Graph.Read_Ncol(g_path, weights=True, directed=False)

    # Original weights
    g.es["weight_original"] = g.es["weight"]

    # Invert weights: max - current (min becomes 0)
    max_weight = max(g.es["weight"])
    g.es["weight_inverted"] = [max_weight - w for w in g.es["weight"]]

    # Loop over original and inverted weights
    for weight_type in ["weight_original", "weight_inverted"]:
        print(f"\nTesting with {weight_type}...")

        # Loop over the desired cluster counts (None and 2 to 10)
        for cluster_count in [None] + list(range(2, 11)):
            print(f"\nParameters: weight = {weight_type}, clusters = {cluster_count}")

            try:
                communities = g.community_edge_betweenness(
                    weights=weight_type,
                    clusters=cluster_count
                ).as_clustering()

                num_communities = len(communities)
                print(f"Result: Number of communities = {num_communities}")
                #visualize_communities(g, communities, title=f"{weight_type}, {cluster_count}")
            except Exception as e:
                print(f"Error with parameters weight = {weight_type}, clusters = {cluster_count}: {e}")

# List of networks to test
networks = [
    "local_code/create_graph_v1/data_graphs/graph_0_1_v1.txt",
    "local_code/create_graph_v1/data_graphs/laplacian_normalized_edge_weights.txt",
    "local_code/create_graph_v1/data_graphs/hits_as_weights.txt",
]

# Apply clustering algorithms to each network
for graph in networks:
    print(f"\nTesting graph: {graph}")
    test_clustering_algorithms(graph)
