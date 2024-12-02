import pandas as pd
import json

# Load the CSV file
def collapse_matrix(csv_file_path, outpath):
    df = pd.read_csv(csv_file_path)

    # Load the JSON file
    json_file_path = "local_code\data\human_accession_to_human_domain.json"
    with open(json_file_path, 'r') as f:
        protein_to_domain = json.load(f)

    # Map human_accession to domains
    df['Domain'] = df['human_accession'].map(protein_to_domain)

    # Drop rows where Domain is NaN (i.e., proteins not in the JSON mapping)
    df = df.dropna(subset=['Domain'])

    # Group by Domain and sum the values
    collapsed_df = df.groupby('Domain').sum(numeric_only=True)

    # Display or save the result
    collapsed_df.reset_index(inplace=True)
    print(collapsed_df)
    collapsed_df.to_csv(outpath, index=False)

collapse_matrix("local_code\create_graph_v1\data_graphs\\not_correct\domain_virus_family_matrix.csv",
                "local_code\graph_v2_virusFam_and_humanDomain\graph_data\hit_matrix.csv")