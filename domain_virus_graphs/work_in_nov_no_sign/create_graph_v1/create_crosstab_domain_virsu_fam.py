import pandas as pd

# Load your data from Excel
file_path = "local_code\data\Compiled_viral_results_compiled.xlsx"  # Replace with your actual file path
df = pd.read_excel(file_path)

# Create a crosstab to count occurrences of each (domain, virus family) pair
matrix = pd.crosstab(df['human_accession'], df['virus_taxa_family'])

# Save the resulting matrix to a new Excel file
output_path = "local_code\create_graph_v1\data_graphs\domain_virus_family_matrix.csv"
matrix.to_csv(output_path)

print(f"Matrix saved to {output_path}")