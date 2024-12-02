import pandas as pd

# Load your Excel file
file_path = "local_code\data\Compiled_viral_results_compiled.xlsx"  # Replace with your actual file name

# Read the Excel file
df = pd.read_excel(file_path)

# Create the mapping dictionary
mapping_dict = dict(zip(df['human_accession'], df['human_domain_name']))

# Save the mapping dictionary to a file (JSON format)
output_path = "local_code\data\human_accession_to_human_domain.json"  # Replace with your desired output file name
with open(output_path, "w") as f:
    import json
    json.dump(mapping_dict, f, indent=4)

print(f"Mapping dictionary saved to {output_path}")
