import json
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os, glob


def get_json_paths(folder_path):
    """
    Retrieves all JSON file paths from the specified folder.

    Parameters:
        folder_path (str): Path to the folder containing JSON files.

    Returns:
        list: A list of file paths to all JSON files in the folder.
    """
    # Use glob to find all JSON files in the folder
    json_files = glob.glob(os.path.join(folder_path, "*.json"))
    return json_files

# Define the folder path
folder_path = "z_new_pipeline/stringDB/filtered"

# Get all JSON file paths from the folder
paths = get_json_paths(folder_path)

# Initialize a DataFrame to store the data
heatmap_data = pd.DataFrame()

# Iterate through the JSON files
for i, path in enumerate(paths):
    with open(path, 'r') as json_file:
        data = json.load(json_file)
        for row in data:
            description = row["description"]
            if description == "Diseases of signal transduction by growth factor receptors and second messengers":
                description = "Growth factor receptors"
                
            p_value = float(row["p_value"])
            #if p_value < 0.0000001:
            # Populate the DataFrame with descriptions and p-values
            heatmap_data.loc[description, f"Cluster {i+1}"] = p_value

# Transform p-values to -log10(p_value) for better visualization
log_heatmap_data = -np.log10(heatmap_data)

# Keep only the 5 most significant hits for each cluster
filtered_data = pd.DataFrame()
for cluster in log_heatmap_data.columns:
    cluster_data = log_heatmap_data[cluster].dropna().sort_values(ascending=False).head(5)
    filtered_data = pd.concat([filtered_data, cluster_data], axis=1)

# Plot the heatmap
plt.figure(figsize=(10, 8))
plt.imshow(filtered_data, aspect='auto', cmap='viridis', origin='lower')

# Add labels and title
cbar = plt.colorbar(label='-log10(p-value)      (yellow is better)')
plt.xticks(range(filtered_data.shape[1]), filtered_data.columns, rotation=45, ha='right')
plt.yticks(range(filtered_data.shape[0]), filtered_data.index)
plt.title("Clusters and Pathways")


# Show the plot
plt.tight_layout()
plt.show()