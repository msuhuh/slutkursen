import os
import re
from collections import defaultdict


######### PARAMETERS ###############
# Parent folder containing the subfolders (e.g., fastgreedy, etc.)
parent_folder_path = "src/00_data_folder/02_01_virus_outputs"
output_file = "src/00_data_folder/02_02_overalpping_virus_families/summary.txt"

###################################

# Regex to parse each line of the file
pattern = re.compile(r"(TRUE|FALSE): Cluster = \[(.+?)\], Count = (\d+)")

# Dictionary to store cluster counts for each subfolder
results = {}

# Iterate over all subfolders in the parent folder
for subfolder_name in os.listdir(parent_folder_path):
    subfolder_path = os.path.join(parent_folder_path, subfolder_name)
    if os.path.isdir(subfolder_path):  # Only process directories
        cluster_true_counts = defaultdict(int)
        total_files = 0

        # Process all .txt files in the subfolder
        for file_name in os.listdir(subfolder_path):
            if file_name.endswith(".txt"):
                total_files += 1
                file_path = os.path.join(subfolder_path, file_name)
                with open(file_path, 'r') as file:
                    for line in file:
                        match = pattern.match(line.strip())
                        if match:
                            status, cluster, count = match.groups()
                            cluster = tuple(cluster.split(", "))  # Convert to a tuple for consistent hashing
                            if status == "TRUE":
                                cluster_true_counts[cluster] += 1

        # Convert counts to percentages
        cluster_true_percentages = {
            cluster: (true_count / total_files) * 100
            for cluster, true_count in cluster_true_counts.items()
        }

        # Store results for this subfolder
        results[subfolder_name] = cluster_true_percentages

# Write the aggregated results to a new .txt file
with open(output_file, 'w') as file:
    file.write("Cluster TRUE Percentages by Subfolder:\n")
    for subfolder_name, cluster_true_percentages in results.items():
        file.write(f"\nCluster algorithm: {subfolder_name}\n")
        for cluster, percentage in cluster_true_percentages.items():
            cluster_str = ", ".join(cluster)  # Convert cluster tuple back to a string
            file.write(f"Cluster: [{cluster_str}], TRUE Percentage: {percentage:.2f}%\n")

print(f"Results have been saved to {output_file}")
