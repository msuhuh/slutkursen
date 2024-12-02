from igraph import Graph, plot
import matplotlib.pyplot as plt
import numpy as np

np.random.seed(2024)


def visu_raw_graph(path_to_txt):

    # Load the graph from the edge list
    graph = Graph.Read_Ncol(path_to_txt, weights=True, directed=False)


    # Load the virus names from a text file
    with open(r"local_code\create_graph_v1\data_graphs\unique_virus_fam_names.txt", "r") as file:
        virus_names = set(file.read().strip().split(","))  # Read and split by commas

    # Classify nodes as either virus family or human accession
    node_colors = []
    node_labels = []
    for node in graph.vs["name"]:
        if node in virus_names:  # Check if the node is in the virus names list
            node_colors.append("red")  # Virus families are red
            node_labels.append(node)  # Add the label for virus families
        else:
            node_colors.append("blue")  # Human accessions are blue
            node_labels.append("")  # No label for human accessions

    # Set colors for visualization
    graph.vs["color"] = node_colors

    # Set edge colors based on weights (higher weights = more red)
    max_weight = max(graph.es["weight"])
    min_weight = min(graph.es["weight"])
    edge_colors = [
        (1, 0, 0, (weight - min_weight) / (max_weight - min_weight))  # RGBA for red with alpha
        for weight in graph.es["weight"]
    ]

    # Plot the graph
    layout = graph.layout("fr")  # Fruchterman-Reingold layout

    fig, ax = plt.subplots(figsize=(10, 10))

    # Draw edges with colors
    for edge, color in zip(graph.get_edgelist(), edge_colors):
        start, end = layout[edge[0]], layout[edge[1]]
        ax.plot(
            [start[0], end[0]],
            [start[1], end[1]],
            color=color,
            linewidth=2,
            alpha=color[3],  # Use alpha channel for transparency
        )

    # Draw nodes with colors
    x, y = zip(*layout)
    ax.scatter(
        x,
        y,
        c=node_colors,
        s=100,
        edgecolors="k",
        zorder=2,
    )

    # Add labels only for virus families
    for i, (x_coord, y_coord) in enumerate(layout):
        if node_labels[i]:  # Only if a label exists
            ax.text(
                x_coord,
                y_coord,
                node_labels[i],
                fontsize=8,
                ha="right",
                color="black",
            )

    # Remove axes for better visualization
    ax.axis("off")
    plt.show()


networks = [
    "local_code\create_graph_v1\data_graphs\graph_0_1_v1.txt",
    "local_code\create_graph_v1\data_graphs\laplacian_normalized_edge_weights.txt",
    "local_code\create_graph_v1\data_graphs\hits_as_weights.txt"
]

for graph in networks:
    visu_raw_graph(graph)