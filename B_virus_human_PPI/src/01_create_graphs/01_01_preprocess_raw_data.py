import pandas as pd


#### PARAMETER #####

# Define the input CSV file path
input_csv_path = 'src/00_data_folder/00_original_data/20241119_Compiled_viral_results_compiled.xlsx - Sheet1.csv'
output_csv_path = 'src/00_data_folder/00_original_data/preprocessed_data.csv'

#################

# Load the CSV file into a DataFrame
df = pd.read_csv(input_csv_path)

# Extract the desired columns
extracted_columns = df[["human_accession", "human_bait_id", "virus_accession", "virus_taxa_sp", "virus_taxa_family", "nucleic_acid"]]

# Save to a new CSV file
extracted_columns.to_csv(output_csv_path, index=False)
print(f"Extracted data saved to: {output_csv_path}")