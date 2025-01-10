import os
import random

##### PARAMETERS #####

# Configuration
input_file = "src/00_data_folder/01_ouput_graphs/normalized_graph_ncol.txt"
output_folder = "src/00_data_folder/01_ouput_graphs/random_graphs"
os.makedirs(output_folder, exist_ok=True)

num_graphs = 5000

#####################


# Ensure the output folder exists
os.makedirs(output_folder, exist_ok=True)

# Read the input file
with open(input_file, "r") as file:
    lines = file.readlines()

# Parse the file content
edges = [line.strip().split() for line in lines[1:]]
edge_weights = [edge[2] for edge in edges]

# Generate randomized graphs
for i in range(num_graphs):
    if i % 500 == 0:
        print(i) # flag

    random.shuffle(edge_weights)
    randomized_edges = [
        f"{edge[0]} {edge[1]} {weight}" for edge, weight in zip(edges, edge_weights)
    ]

    # Save the randomized graph
    output_file = os.path.join(output_folder, f"{i + 1}.txt")
    with open(output_file, "w") as out_file:
        out_file.write("\n".join(randomized_edges) + "\n")

print(f"Generated {num_graphs} randomized graphs in '{output_folder}'")