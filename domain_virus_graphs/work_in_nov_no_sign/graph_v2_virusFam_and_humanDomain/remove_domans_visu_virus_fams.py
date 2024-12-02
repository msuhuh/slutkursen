from igraph import Graph, plot
import matplotlib.pyplot as plt
import numpy as np

np.random.seed(2024)


def visualize_modified_graph(path_to_txt, virus_names_file):
    # Load the graph from the edge list
    graph = Graph.Read_Ncol(path_to_txt, weights=True, directed=False)

    # Load the virus names from a text file
    with open(virus_names_file, "r") as file:
        virus_names = set(file.read().strip().split(","))  # Read and split by commas

    # Classify nodes as either virus family (red) or human accession (blue)
    blue_nodes = []
    red_nodes = set()
    for node in graph.vs["name"]:
        if node in virus_names:  # Virus family
            red_nodes.add(node)
        else:  # Human accession
            blue_nodes.append(node)

    # Create a new graph by connecting red nodes that were linked via blue nodes
    new_edges = set()
    for blue_node in blue_nodes:
        # Find all red nodes connected to this blue node
        neighbors = graph.neighbors(blue_node, mode="all")
        red_neighbors = [
            graph.vs[n]["name"] for n in neighbors if graph.vs[n]["name"] in red_nodes
        ]

        # Add edges between all red neighbors of this blue node
        for i in range(len(red_neighbors)):
            for j in range(i + 1, len(red_neighbors)):
                new_edges.add((red_neighbors[i], red_neighbors[j]))

    # Create a new graph with only red nodes and new edges
    red_node_names = [v["name"] for v in graph.vs if v["name"] in red_nodes]
    subgraph = Graph()
    subgraph.add_vertices(red_node_names)
    subgraph.add_edges(list(new_edges))  # Add edges using names directly

    # Print node degrees
    for vertex in subgraph.vs:
        print(f"Node {vertex['name']}: Degree {subgraph.degree(vertex)}")
    # Visualize the resulting graph
    layout = subgraph.layout("fr")  # Fruchterman-Reingold layout

    fig, ax = plt.subplots(figsize=(10, 10))

    # Draw edges
    for edge in subgraph.es:
        start, end = layout[edge.source], layout[edge.target]
        ax.plot(
            [start[0], end[0]],
            [start[1], end[1]],
            color="black",
            linewidth=1,
        )

    # Draw nodes
    x, y = zip(*layout)
    ax.scatter(
        x,
        y,
        c="red",  # All nodes are now red
        s=100,
        edgecolors="k",
        zorder=2,
    )

    # Add labels for virus families
    for i, (x_coord, y_coord) in enumerate(layout):
        ax.text(
            x_coord,
            y_coord,
            subgraph.vs["name"][i],
            fontsize=8,
            ha="right",
            color="black",
        )

    # Remove axes for better visualization
    ax.axis("off")
    plt.show()


# Example usage
virus_names_file = r"local_code\create_graph_v1\data_graphs\unique_virus_fam_names.txt"
path_to_graph = r"local_code\graph_v2_virusFam_and_humanDomain\graph_data\hit_matrix.txt"

visualize_modified_graph(path_to_graph, virus_names_file)
