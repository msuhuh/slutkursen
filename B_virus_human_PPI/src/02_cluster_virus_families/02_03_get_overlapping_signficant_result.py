import re
from collections import defaultdict
from itertools import combinations


############### PARAMETETRS ##############

# Input and output files
input_file = "src/00_data_folder/02_02_overalpping_virus_families/summary.txt"
intermediate_file = "src/00_data_folder/02_02_overalpping_virus_families/analysis_summary_filter_3.txt"
output_file_pruned = "src/00_data_folder/02_02_overalpping_virus_families/analysis_summary_filter_3_pruned.txt"

##########################################

# Regex to parse clusters and subfolders
subfolder_pattern = re.compile(r"Cluster algorithm: (.+)")
cluster_pattern = re.compile(r"Cluster: \[(.+?)\], TRUE Percentage: (\d+\.\d+)%")

# Dictionary to store clusters by subfolder
subfolder_clusters = {}

# Read the input file and extract clusters by subfolder
current_subfolder = None
with open(input_file, 'r') as file:
    for line in file:
        subfolder_match = subfolder_pattern.match(line.strip())
        if subfolder_match:
            current_subfolder = subfolder_match.group(1)
            subfolder_clusters[current_subfolder] = set()
        elif current_subfolder:
            cluster_match = cluster_pattern.match(line.strip())
            if cluster_match:
                cluster = tuple(sorted(cluster_match.group(1).split(", ")))  # Sort for consistency
                percentage = float(cluster_match.group(2))
                if percentage == 100.0:  # Include only clusters with 100% TRUE Percentage
                    subfolder_clusters[current_subfolder].add(cluster)

# Analyze shared pairs, triplets, etc., across subfolders
shared_combinations = defaultdict(lambda: defaultdict(list))

# Generate combinations and check for sharing
for subfolder, clusters in subfolder_clusters.items():
    for cluster in clusters:
        for size in range(2, len(cluster) + 1):
            for comb in combinations(cluster, size):
                shared_combinations[size][comb].append(subfolder)

# Sort the combinations by the number of subfolders they are shared among
sorted_combinations = {
    size: sorted(combinations_dict.items(), key=lambda x: len(x[1]), reverse=True)
    for size, combinations_dict in shared_combinations.items()
}

# Write results to the output file
with open(intermediate_file, 'w') as file:
    for size, combinations_list in sorted_combinations.items():
        file.write(f"\nShared {size}-tuples:\n")
        for comb, subfolders in combinations_list:
            if len(subfolders) >= 3:  # Only include combinations shared among at least 3 subfolders
                comb_str = ", ".join(comb)
                subfolders_str = ", ".join(subfolders)
                file.write(f"Combination: [{comb_str}] is shared among Cluster algorithms: {subfolders_str}\n")

print(f"Shared cluster analysis results have been saved to {intermediate_file}")


# Prune data
# Dictionary to store parsed combinations
combination_map = defaultdict(list)

# Read and parse the intermediate file
with open(intermediate_file, 'r') as file:
    current_size = None
    for line in file:
        if line.startswith("Shared"):
            current_size = int(line.split()[1].split("-")[0])
        elif line.startswith("Combination: ["):
            comb_part = line.split("is shared among Cluster algorithms:")
            comb = tuple(comb_part[0].replace("Combination: [", "").replace("]", "").strip().split(", "))
            algorithms = comb_part[1].strip().split(", ")
            combination_map[current_size].append((comb, algorithms))

# Filter subsets if identical algorithms exist in larger clusters
for size, combinations_list in reversed(sorted(combination_map.items())):
    for row_index, (current_virus_found, current_algos_used) in enumerate(combinations_list):
        for smaller_size in range(size - 1, 0, -1):
            if smaller_size in combination_map:
                combination_map[smaller_size] = [
                    (virus, algos)
                    for virus, algos in combination_map[smaller_size]
                    if not (set(virus).issubset(current_virus_found) and algos == current_algos_used)
                ]

# Write the filtered results to the new output file
with open(output_file_pruned, 'w') as file:
    for size, combinations_list in sorted(combination_map.items()):
        file.write(f"\nShared {size}-tuples:\n")
        for comb, algorithms in combinations_list:
            file.write(f"[{', '.join(comb)}]; [{', '.join(algorithms)}]\n")

print(f"Grouped and filtered cluster analysis results have been saved to {output_file_pruned}")