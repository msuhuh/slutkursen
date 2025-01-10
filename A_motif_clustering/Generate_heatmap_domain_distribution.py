import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Path to the Excel file
sequence_path = "/Users/juliaancker/Desktop/projekt_bioinf/results/test_set_domain_stats.xlsx" # Replace with your path

# Load the data
df = pd.read_excel(sequence_path)

# Define non-domain columns, these are the non domain columns from the Hammock algorithm output
non_domain_columns = ['cluster_id', 'New_name', 'Stat scores', 'main_sequence', 'sum']

# Dynamically detect domain columns
domain_columns = [col for col in df.columns if col not in non_domain_columns]

# Filter out domains that are not in any cluster 
active_domains = [col for col in domain_columns if df[col].sum() > 0]

# Subset the data with relevant columns
heatmap_data = df.set_index('cluster_id')[active_domains]

# Normalize data within each cluster (row-wise normalization) to show precentages of each domain in a cluster
row_normalized_data = heatmap_data.div(heatmap_data.sum(axis=1), axis=0).fillna(0) * 100  # Convert to percentages

# Replace zero values with NaN to display them as white in the heatmap
row_normalized_data.replace(0, np.nan, inplace=True)


# Create the heatmap
plt.figure(figsize=(len(row_normalized_data.index) * 0.5, 10))  # Adjust width dynamically
sns.heatmap(
    row_normalized_data.T,  # Transpose so domains are on the y-axis and clusters on the x-axis
    cmap="inferno_r",  # Choose a color map
    cbar=True,
    linewidths=0.5,
    linecolor="grey",
    annot=False,  # Remove the numbers from the cells
    mask=row_normalized_data.T.isna()  # Use a mask for NaN values (white cells)
)

# Add titles and labels
plt.title("Domain distribution", fontsize=16)
plt.xlabel("Clusters", fontsize=12)
plt.ylabel("Domains", fontsize=12)

# Rotate x-axis labels for better readability
plt.xticks(rotation=45, ha="right", fontsize=10)
plt.yticks(fontsize=10)

# Save the heatmap
plt.tight_layout()
plt.savefig("heatmapdomain.png") # Replace with your desired path to output file
plt.show()
