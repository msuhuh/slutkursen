import pandas as pd

# Paths for input and output files
path_input = '/Users/juliaancker/Desktop/projekt_bioinf/data_formated2.xlsx'  # Adjust extension if needed
path_output = '/Users/juliaancker/Desktop/projekt_bioinf/filtered_output.xlsx'

# Load the input Excel file into a DataFrame
df = pd.read_excel(path_input)

# Remove duplicate rows based on the specified columns
df_filtered = df.drop_duplicates(subset=['viral_hit_collapsed', 'human_bait_id', 'virus_taxa_id'])

# Save the filtered DataFrame to a new Excel file
df_filtered.to_excel(path_output, index=False)

print(f"Filtered data has been saved to {path_output}")
