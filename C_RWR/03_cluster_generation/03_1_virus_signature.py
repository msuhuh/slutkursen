'''
This script is a modified version of the virus_signature.py file from 
the paper 'Large-scale phage-based screening reveals extensive pan-viral 
mimicry of host short linear motifs' by Mihalič et al. Nature communications 2023, 
in order to run the algorithm for each virus family instead.

Modified by: Minna Sayehban, student at Uppsala University
E-mail: minna.sayehban.1224@student.uu.se or minna.sayehban@gmail.com
Date: 28-11-2024
'''

############### IMPORT & LOAD NECESSARY MODULES ################
import nimfa
from sklearn.preprocessing import QuantileTransformer
import sys, os
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from collections import Counter


# Check for command-line arguments
if len(sys.argv) < 2:
    print("Usage: python virus_signature.py <global_frequency_threshold>")
    sys.exit(1)

# Parse the threshold from command-line arguments
threshold = float(sys.argv[1])  # Example: python virus_signature.py 0.3
print(f"Using threshold: {threshold}")

# Load PPI network
net = nx.read_weighted_edgelist("/networks/TCSS.txt")

sns.set(font_scale=0.6)
plt.rcParams['font.family'] = 'Arial'

# Load UniProt to gene name mapping
def load_uniprot_to_gene_mapping():
    uniprot_to_gene = {}
    f1 = open("/data/uniprot_to_gene.tab", "r")
    seq = f1.readline()
    while seq != "":
        seq = seq.strip().split("\t")
        uniprot_to_gene[seq[0]] = seq[1].split(";")[0].strip()
        seq = f1.readline()
    return uniprot_to_gene

# Load gene mapping
uniprot_to_gene = load_uniprot_to_gene_mapping()

# Base folder containing the results of the RWR propagation for each virus family
basefolder = "/results/family/"
family_list = os.listdir(basefolder)

# Initialize data structures to store RWR results, seed nodes, etc.
virus_net = {}
virus_net_rwr = {}
virus_rwr = {}
all_nodes = set()
virus_seed = {}
rwr_values = []
seed_nodes = set()

# Load the starting data for each virus family
for family in family_list:
    virus_seed[family] = []

    # Read start seed nodes
    f1 = open(os.path.join(basefolder, family, "start_seeds.txt"))
    seq = f1.readline()
    while seq != "":
        seq = seq.strip().split("\t")
        virus_seed[family].append(seq[0])
        seed_nodes.add(seq[0])
        all_nodes.add(seq[0])
        seq = f1.readline()

    # Read p-values (overrepresented nodes)
    f1 = open(os.path.join(basefolder, family, "significance_hits.txt"))
    seq = f1.readline()
    virus_rwr[family] = {}
    temp = []
    while seq != "":
        seq = seq.strip().split("\t")
        if float(seq[1]) > 990:  # This threshold determines "overrepresentation"
            temp.append(seq[0])
            all_nodes.add(seq[0])
        seq = f1.readline()
    virus_net[family] = temp + virus_seed[family]

    # Read RWR results
    f1 = open(os.path.join(basefolder, family, "rwr_scores.txt"))
    seq = f1.readline()
    while seq != "":
        seq = seq.strip().split("\t")
        if seq[0] in virus_net[family]:
            virus_rwr[family][seq[0]] = float(seq[1])
        else:
            rwr_values.append(float(seq[1]))
        seq = f1.readline()

    # Summing RWR values for each family
    for protein_id in virus_net[family]:
        if protein_id in virus_rwr[family]:
            if protein_id in virus_net_rwr:
                virus_net_rwr[protein_id][family] += virus_rwr[family][protein_id]
            else:
                virus_net_rwr[protein_id] = {family: 0.0 for family in family_list}
                virus_net_rwr[protein_id][family] = virus_rwr[family][protein_id]
        else:
            print(f"Warning: Protein {protein_id} in virus_net[{family}] has no RWR score. Skipping.")


# Calculate Z-scores and select highly hijacked proteins
mean = np.mean(rwr_values)
std = np.std(rwr_values)
f1 = open("/results/most_hijacked_proteins.txt", "w")
f1.write("uniprot\tzscore\tis_seed\n")
for protein_id in virus_net_rwr:
    # Normalize by the number of virus families the protein is involved with
    virus_net_rwr[protein_id] = {k: v / len(virus_net) for k, v in virus_net_rwr[protein_id].items()}
    temp = [value for value in virus_net_rwr[protein_id].values() if value > 0]
    
    # Proteins must be present in more than 7 families and have z-score > 2.326
    if len(temp) > 7:
        num = np.mean(temp) - mean
        den = std
        if (num / den) > 2.326:
            f1.write(protein_id + "\t" + str(num / den) + "\t" + str(protein_id in seed_nodes) + "\n")
f1.close()

# Generate signature matrix
column = {}
col = []
for i, protein_id in enumerate(sorted(all_nodes)):
    column[protein_id] = i
    col.append(protein_id)

mat_family_global = {}
for family in virus_net:
    temp = np.zeros(len(column))
    for protein_id in virus_net[family]:
        if protein_id in virus_rwr[family]:  # Check if RWR score exists
            temp[column[protein_id]] = virus_rwr[family][protein_id]
        else:
            print(f"Warning: Protein {protein_id} in family {family} has no RWR score. Skipping.")
            temp[column[protein_id]] = 0  # Assign a default value (e.g., 0)

    if family in mat_family_global:
        mat_family_global[family] = np.add(mat_family_global[family], temp)
    else:
        mat_family_global[family] = temp

# Write signature matrix to file
f1 = open("/results/signature_matrix.txt", "w")
fam = []
mat = []
for family in mat_family_global:
    f1.write(family + "\t" + "\t".join(map(str, mat_family_global[family])) + "\n")
    mat.append(np.array(mat_family_global[family]))
    fam.append(family)
f1.close()

# The signature matrix
mat = np.array(mat)
mat = mat.T

# Quantile normalize the matrix
qt = QuantileTransformer(n_quantiles=min(263, len(mat)), random_state=0)
mat = qt.fit_transform(mat)

# Perform the NMF
rank = min(6, mat.shape[1])
nmf = nimfa.Nmf(mat, seed='nndsvd', max_iter=500, rank=rank)
nmf_fit = nmf()
w = nmf_fit.basis()
h = nmf_fit.coef()
h = np.asarray(h)

# Extract clusters by assigning the maximum latent factor
clusters = {}
for i, h_val in enumerate(h.T):
    max_factor = np.argmax(h_val)
    if max_factor in clusters:
        clusters[max_factor].append(fam[i])
    else:
        clusters[max_factor] = [fam[i]]

# Write cluster assignment for each family
f1 = open("/results/cluster_family.txt", "w")
proteins_in_clusters = {}
for cluster_id in clusters:
    proteins_in_clusters[cluster_id] = set()
    f1.write(str(cluster_id + 1) + "\t" + ", ".join(clusters[cluster_id]) + "\n")
    for family in clusters[cluster_id]:
        proteins_in_clusters[cluster_id].update(virus_net[family])
f1.close()

# Parsing Fisher's Results for Annotations
annotation_family_virus = {}
all_annot = set()
for cluster_id, cluster_families in clusters.items():
    annotation_family_virus[cluster_id] = {}

    for family in cluster_families:
        fisher_file = os.path.join(basefolder, family, "fisher", "RTfisher.txt")
        if not os.path.exists(fisher_file):
            print(f"Warning: Fisher's results file not found for {family}.")
            continue

        with open(fisher_file, "r") as f:
            seq = f.readline()
            while seq:
                fields = seq.strip().split("\t")
                reactome_id = fields[2]
                proteins = fields[3:]  # Assuming the proteins are listed after the third column
                all_annot.add(reactome_id)

                if family not in annotation_family_virus[cluster_id]:
                    annotation_family_virus[cluster_id][family] = Counter()
                annotation_family_virus[cluster_id][family].update([reactome_id])
                seq = f.readline()

# Generate Data for CSV
data = []
row = []
for cluster_id, family_annotations in annotation_family_virus.items():
    temp = {}
    row.append(cluster_id + 1)

    for reactome_id in all_annot:
        global_count = 0  # Count how many families in the cluster have this pathway

        for family, annotations in family_annotations.items():
            if reactome_id in annotations:
                global_count += 1

        # Compute global frequency
        global_freq = global_count / len(family_annotations)  # Fraction of families in the cluster with this pathway

        # Apply the user-specified threshold
        if global_freq >= threshold:
            temp[reactome_id] = global_freq

    data.append({"cluster": f"C{cluster_id + 1}", "pathways": temp})

# Prepare the data for CSV
# Collect all pathways across clusters (to ensure consistency in rows)
all_pathways = set()
for cluster_data in data:
    all_pathways.update(cluster_data["pathways"].keys())
all_pathways = sorted(all_pathways)  # Sort pathways alphabetically or as needed

# Prepare a DataFrame with pathways as rows and clusters as columns
cluster_ids = [f"C{cluster_id + 1}" for cluster_id in range(len(annotation_family_virus))]
df_data = {cluster: [] for cluster in cluster_ids}

for pathway in all_pathways:
    for cluster_data in data:
        cluster_id = cluster_data["cluster"]
        df_data[cluster_id].append(cluster_data["pathways"].get(pathway, None))  # None for missing values

# Create a DataFrame
df = pd.DataFrame(df_data, index=all_pathways)
df.index.name = "Pathway"

# Convert the threshold to a string with a fixed number of decimal places (e.g., 1 decimal place)
threshold_str = str(threshold).replace('.', '_')  # Replacing the dot with an underscore for filename compatibility

# Create the output CSV filename with the threshold included
output_csv = f"/results/c_heatmap_data_threshold_{threshold_str}.csv"

# Save DataFrame to CSV
df.to_csv(output_csv, na_rep="")  # Empty cell for pathways not present in a cluster

print(f"Data saved to {output_csv}")



