import pandas as pd
import json

def rename_csv_rows(input_csv, json_file, output_csv):
    """
    Rename rows in a CSV file based on a JSON mapping and save a new CSV.

    Parameters:
        input_csv (str): Path to the input CSV file.
        json_file (str): Path to the JSON file with the mapping.
        output_csv (str): Path to save the output CSV file.
        column_to_rename (str): Column name in the CSV to rename based on the JSON.

    Returns:
        None
    """
    # Load the CSV file
    df = pd.read_csv(input_csv)

    duplicates = df[df.duplicated(subset=df.columns[0], keep=False)]
    if not duplicates.empty:
        print("Duplicate rows after remapping:")
        print(duplicates)
    else:
        print("No duplicates in original")

    # Load the JSON mapping
    with open(json_file, 'r') as file:
        id_to_accession_mapping = json.load(file)

    # Rename rows in the specified column using the JSON mapping
    df.iloc[:, 0] = df.iloc[:, 0].map(id_to_accession_mapping).fillna(df.iloc[:, 0])

    duplicates = df[df.duplicated(subset=df.columns[0], keep=False)]
    if not duplicates.empty:
        print("Ampunt of duplicate rows after remapping:")
        print(len(duplicates))

        # Collapse rows that share the same name after renaming
    if not duplicates.empty:
        print("Collapsing duplicate rows...")
        df = df.groupby(df.columns[0], as_index=False).sum()
        print("Duplicates collapsed successfully.")

    # Ensure no duplicates remain
    assert df.iloc[:, 0].is_unique, "Duplicate rows still exist after collapsing!"


    # Save the updated DataFrame to a new CSV
    df.to_csv(output_csv, index=False)

    print(f"Renamed CSV saved to {output_csv}")




##### PARAMETERS #########
input_csv_file = "src/00_data_folder/03_01_get_baits_from_clusters/combined_virus_bait_matrix.csv"
json_mapping_file = "src/00_data_folder/00_original_data/map_baits_to_accession.json"
output_csv_file = "src/00_data_folder/03_02_human_accession_ids/virus_accession_matrix.csv"
##########################


rename_csv_rows(
    input_csv=input_csv_file,
    json_file=json_mapping_file,
    output_csv=output_csv_file,
)