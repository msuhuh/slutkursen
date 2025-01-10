import pandas as pd
import os

def extract_rows_based_on_threshold(input_csv, viruses, output_txt):
    """
    Extract row names from a CSV based on overlapping values for multiple columns with a threshold.

    Parameters:
        input_csv (str): Path to the input CSV file.
        viruses (list): List of virus names (column names).
        output_txt (str): Path to save the output text file.

    Returns:
        None
    """
    # Load the CSV file
    df = pd.read_csv(input_csv)

    # Create a mask for the threshold condition across all specified viruses
    condition = (df[viruses] > 0).all(axis=1)

    # Extract the row names (first column) where the condition is True
    row_names = df.loc[condition, df.columns[0]]

    # Save the row names to a text file
    with open(output_txt, 'w') as file:
        for row_name in row_names:
            file.write(f"{row_name}\n")

    print(f"Extracted row names saved to {output_txt}")

########### PARAMETERS #################

# Sort of bad implementaton, sorry:
# The virus family clusters must be in a nested list, see example below.
# The virus families names cen easiest be seen in:
# <src\00_data_folder\02_02_overalpping_virus_families\analysis_summary_filter_3_pruned.txt>

large_list = [

['Caliciviridae', 'Picornaviridae'],
['Coronaviridae', 'Paramyxoviridae'],
['Orthomyxoviridae', 'Paramyxoviridae'],
['Anelloviridae', 'Nairoviridae'],
['Hepeviridae', 'Retroviridae'],
['Kolmioviridae', 'Papillomaviridae'],
['Caliciviridae', 'Paramyxoviridae'],
['Arenaviridae', 'Pneumoviridae'],
['Astroviridae', 'Togaviridae'],
['Astroviridae', 'Papillomaviridae'],
['Coronaviridae', 'Orthomyxoviridae', 'Paramyxoviridae'],
['Anelloviridae', 'Hepeviridae', 'Nairoviridae', 'Retroviridae'],
['Caliciviridae', 'Coronaviridae', 'Paramyxoviridae', 'Picornaviridae'],
['Anelloviridae', 'Hepeviridae', 'Nairoviridae', 'Retroviridae', 'Tobaniviridae'],
['Adenoviridae', 'Filoviridae', 'Herpesviridae', 'Parvoviridae', 'Reoviridae', 'Rhabdoviridae'],
['Adenoviridae', 'Filoviridae', 'Herpesviridae', 'Parvoviridae', 'Poxviridae', 'Reoviridae', 'Rhabdoviridae']
]

# Data from clustering
input_csv_file = "src/00_data_folder/03_02_human_accession_ids/virus_accession_matrix.csv"
output_folder = "src/00_data_folder/03_03_significant_baits"

#####################################################################


for viruses_to_check in large_list:
    name = ""
    for virus in viruses_to_check:
        name = name + virus + "_"

    subdir = os.path.join(output_folder, name)
    os.makedirs(subdir, exist_ok=True)

    filename = ".txt"
    output_file = os.path.join(subdir, filename)
    

    # Extract
    extract_rows_based_on_threshold(
        input_csv=input_csv_file,
        viruses=viruses_to_check,
        output_txt=output_file
    )
