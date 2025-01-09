# Script to make a heatmap of pathways and clusteres
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


filepath = 'Cluster_heatmap_data.csv' #Replace with you filepath to the pathway cluster data

df = pd.read_csv(filepath) # Read csv file

# Pathways found by all algorithms sorted manually into categories, replace with your pathways

pathway_sorted = {
    "Immune system": [
        "Immune System",
        "Immune system Adaptive Immune system",
        "C-type lectin receptors (CLRs)",
        "Fcgamma receptor (FCGR) dependent phagocytosis",
        "Antigen processing: Ubiquitination & Proteasome degradation",
        "DAP12 signaling",
        "FCGR activation"
    ],
    "Cell cycle": [
        "Cell cycle",
        "Mitotic Spindle Checkpoint",
        "Amplification of signal from the kinetochores",
        "Amplification of signal from unattached kinetochores via a MAD2 inhibitory signal",
        "EML4 and NUDC in mitotic spindle formation",
        "Resolution of Sister Chromatid Cohesion"
    ],
    "Gene expression": [
        "Gene expression",
        "YAP1- and WWTR1 (TAZ)-stimulated gene expression",
        "Transcriptional regulation by RUNX2",
        "Transcriptional regulation by RUNX3",
        "RUNX2 regulates bone development",
        "RUNX3 regulates YAP1-mediated transcription",
        "RUNX2 regulates osteoblast differentiation"
    ],
    "Disease": [
        "Disease",
        "Infectious disease",
        "Diseases of signal transduction by growth factor receptors and second messengers",
        "Oncogenic MAPK signaling",
        "Listeria monocytogenes entry into host cells",
        "InlB-mediated entry of Listeria monocytogenes into host cell",
        "Host Interactions of HIV factors",
        "The role of Nef in HIV-1 replication and disease pathogenesis",
        "Nef-mediates down modulation of cell surface receptors by recruiting them to clathrin adapters",
        "Nef Mediated CD4 Down-regulation",
        "Nef Mediated CD8 Down-regulation",
        "Parasite infection",
        "FCGR3A-mediated phagocytosis",
        "Leishmania phagocytosis",
        "FCGR3A-mediated phagocytosis"
    ],
    "Developmental Biology": [
        "Developmental Biology",
        "Nervous system development",
        "Axon guidance",
        "Signaling by ROBO receptors",
        "EPH-Ephrin signaling",
        "EPH-ephrin mediated repulsion of cells",
        "EPHA-mediated growth cone collapse",
        "L1CAM interactions",
        "Recycling pathway of L1"
    ],
    "Vesicle-mediated transport": [
        "Vesicle-mediated transport",
        "Membane trafficking",
        "Clathrin-mediated endocytosis",
        "Cargo recognition for clathrin-mediated endocytosis",
        "TBC/RABGAPs (Homo sapiens)"
    ],
    "Metabolism of RNA": [
        "Metabolism of RNA",
        "Processing of Capped Intron-Containing Pre-mRNA",
        "mRNA splicing",
        "mRNA Splicing - Major Pathway"
    ],
    "Cell - Cell Communication": [
        "Cell - Cell Communication",
        "Nephrin family interactions"
    ],
    "Neuronal System": [
        "Neuronal System",
        "Protein-Protein interactions at synapses",
        "Neurexins and neuroligins"
    ],
    "Autophagy": [
        "Autophagy",
        "Macroautophagy",
        "PINK1-PRKN Mediated Mitophagy",
        "Receptor Mediated Mitophagy"
    ],
    "Hemostasis": [
        "Hemostasis",
        "Platelet activation, signaling and aggregation"
    ],
    "Sensory Perception": [
        "Sensory processing of sound by inner hair cells of the cochlea",
        "Sensory processing of sound by outer hair cells of the cochlea"
    ],
    "Transport of small molecules": [
        "Transport of small molecules",
        "VLDLR internalisation and degradation",
        "LDL clearance"
    ],
    "Signal Transduction": [
        "Signal Transduction",
        "WNT5A-dependent internalization of FZD2, FZD5 and ROR2",
        "WNT5A-dependent internalization of FZD4",
        "Signaling by Rho GTPases",
        "Miro GTPases and RHOBTB3",
        "Signaling by Rho GTPases",
        "RHO GTPase cycle",
        "CDC42 GTPase cycle",
        "RAC1 GTPase cycle",
        "RHOV GTPase cycle",
        "RHOBTB2 GTPase cycle",
        "RHO GTPase Effectors",
        "RHO GTPases Activate Formins",
        "RHO GTPases Activate WASPs and WAVEs",
        "Signaling by Receptor Tyrosine Kinases",
        "Signaling by ERBB4",
        "Downregulation of ERBB4 signaling",
        "EGFR downregulation",
        "Signaling by VEGF",
        "VEGFA-VEGFR2 Pathway",
        "Signaling by NTRKs",
        "Signaling by NTRK1 (TRKA)",
        "Retrograde neurotrophin signalling",
        "Signaling by MET",
        "Signaling by ERBB2",
        "Signaling by EGFR",
        "Regulation of KIT signaling",
        "Signaling by Hippo",
        "Intracellular signaling by second messengers",
        "Regulation of PTEN localization",
        "RAF/MAP kinase cascade"
    ]
}


df = df.drop(columns=['A15', 'A4', 'A5', 'A7', 'A13', 'A10', 'A12'], errors='ignore') # Drop columns if neede

threshold = 0.1 # Set threshold for values to show

filtered_df = df.loc[(df.iloc[:, 1:] > threshold).any(axis=1)]  # Filter values below a threshold


heatmap_data = filtered_df 


# Build an ordered DataFrame based on the sorted pathway dictionary
ordered_rows = []

for category, pathways in pathway_sorted.items():

    # Check which pathways from dictionary that exists in dataframe
    existing_pathways = [heatmap_data[heatmap_data['Pathways'] == pathway].iloc[0] for pathway in pathways if not heatmap_data[heatmap_data['Pathways'] == pathway].empty]
    # Add the header row
    if existing_pathways:
        # Add the header row first
        header_row = pd.Series([category.upper()] + [np.nan] * (heatmap_data.shape[1] - 1), index=heatmap_data.columns)
        ordered_rows.append(header_row)

        # Add the pathways under the header
        ordered_rows.extend(existing_pathways)

# Reverse order of rows to get them in correct order
ordered_rows.reverse()

# Combine the ordered rows into a DataFrame
ordered_df = pd.DataFrame(ordered_rows).reset_index(drop=True)

# Prepare heatmap data
heatmap_data = ordered_df.drop(columns=['Pathways'], errors='ignore')

# Mask values that are 0 so they appear as white
masked_data = np.ma.masked_where(heatmap_data.isna() | (heatmap_data == 0), heatmap_data)

# Create colormap
cmap = plt.cm.viridis.reversed()
cmap.set_bad(color='white')  # Set NaN (headers) to white

# Plot the heatmap
plt.figure(figsize=(12, 10))
plt.imshow(masked_data, aspect='auto', cmap=cmap, origin='lower')

# Add colorbar
cbar = plt.colorbar(label='Enrichment score')

# Set axis labels
plt.xticks(range(heatmap_data.shape[1]), heatmap_data.columns, rotation=45, ha='right')
y_labels = ordered_df['Pathways']
y_tick_labels = []

for label in y_labels:
    if label.isupper():  # Bold for uppercase headers
        formatted_label = label.replace(' ', r'\ ')  # Replace spaces with LaTeX spaces
        y_tick_labels.append(f"$\\bf{{{formatted_label}}}$")
    else:
        y_tick_labels.append(label)  # Non-uppercase labels remain as is

plt.yticks(range(ordered_df.shape[0]), y_tick_labels, va='center')
plt.title("Clusters and Pathways")

# Adjust layout and show the plot
plt.tight_layout()
plt.savefig('A_heatmap_sorted_pathways.png', dpi=300, bbox_inches='tight')
plt.show()