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
    #plt.show()

# Define a function to apply clustering algorithms
def test_clustering_algorithms(g_path):
    g = Graph.Read_Ncol(g_path, weights=True, directed=False)

    # Define clustering methods to test
    clustering_methods = {
        "Edge Betweenness": g.community_edge_betweenness().as_clustering(), # Weights are used to calculate weighted shortest paths, so they are interpreted as distances.
        "Fast Greedy": g.community_fastgreedy().as_clustering(), # A larger edge weight means a stronger connection for this function.
        #"Label Propagation": g.community_label_propagation(), #A larger edge weight means a stronger connection for this function.
        "Louvain": g.community_multilevel(),
        "Spinglass": g.community_spinglass() if g.is_connected() else None,
        "Walktrap": g.community_walktrap().as_clustering(),
        "Leading Eigenvector": g.community_leading_eigenvector(),
    }

    # Apply and visualize each method
    for method_name, communities in clustering_methods.items():
        if communities is not None:
            print(f"\nApplying {method_name} method...")
            visualize_communities(g, communities, title=method_name)
        else:
            print(f"{method_name} method is not applicable to this graph.")

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
