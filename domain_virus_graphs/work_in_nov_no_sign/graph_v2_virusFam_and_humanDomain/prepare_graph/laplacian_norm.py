from igraph import Graph, plot
import matplotlib.pyplot as plt
import numpy as np

np.random.seed(2024)

# Load the graph from the edge list
graph = Graph.Read_Ncol(r"local_code\graph_v2_virusFam_and_humanDomain\graph_data\HumanDomain_virusFamily_matrix_norm_0_1.txt", weights=True, directed=False)

# Calculate degrees for all nodes
degrees = graph.degree()

# Perform Laplacian normalization for edge weights
laplacian_weights = []
for edge in graph.es:
    source, target = edge.tuple  # Get the source and target node indices
    weight = edge["weight"]  # Original weight
    d_i = degrees[source]  # Degree of source node
    d_j = degrees[target]  # Degree of target node
    normalized_weight = weight / np.sqrt(d_i * d_j) if d_i > 0 and d_j > 0 else 0
    laplacian_weights.append(normalized_weight)

# Update graph edge weights with normalized values
graph.es["weight"] = laplacian_weights

# Save Laplacian-normalized weights to a text file
with open(r"local_code\graph_v2_virusFam_and_humanDomain\graph_data\HumanDomain_virusFamily_laplacian.txt", "w") as file:
    for edge, norm_weight in zip(graph.get_edgelist(), laplacian_weights):
        source_name = graph.vs[edge[0]]["name"]
        target_name = graph.vs[edge[1]]["name"]
        file.write(f"{source_name}\t{target_name}\t{norm_weight}\n")

# Visualization of the graph (optional)
layout = graph.layout("fr")  # Fruchterman-Reingold layout

fig, ax = plt.subplots(figsize=(10, 10))

# Draw edges with Laplacian-normalized weights as transparency (alpha)
max_weight = max(laplacian_weights)
min_weight = min(laplacian_weights)
edge_colors = [
    (1, 0, 0, (weight - min_weight) / (max_weight - min_weight))  # Red with alpha based on weight
    for weight in laplacian_weights
]

for edge, color in zip(graph.get_edgelist(), edge_colors):
    start, end = layout[edge[0]], layout[edge[1]]
    ax.plot(
        [start[0], end[0]],
        [start[1], end[1]],
        color=color,
        linewidth=2,
        alpha=color[3],  # Use alpha channel for transparency
    )

# Draw nodes with fixed sizes
x, y = zip(*layout)
ax.scatter(
    x,
    y,
    c="blue",
    s=100,
    edgecolors="k",
    zorder=2,
)

# Remove axes for better visualization
ax.axis("off")
plt.title("Graph with Laplacian-Normalized Edge Weights")
plt.show()
