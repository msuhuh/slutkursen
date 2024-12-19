import pandas as pd

# Script to create statistics of the Data

# Load the Excel file
file_path = '202041125_test_set.xlsx'  # Replace with your file path
df = pd.read_excel(file_path)
 
domain = df['human_domain_name'] # Name of human domain
viral_family = df['virus_taxa_family'] # Name of viral family
n_acid = df['nucleic_acid'] # DNA / RNA

# Function to make data counts the 

def count_data(column):
    count_dict = {}
    for name in column:
        if name in count_dict:
            count_dict[name] += 1
        else:
            count_dict[name] = 1
    return count_dict

# Statistict for each column

domain_count = count_data(domain)
viral_family_count = count_data(viral_family)
n_acid_count = count_data(n_acid)

domain_df = pd.DataFrame(domain_count.items(), columns=['Domain Name', 'Domain Count'])
viral_family_df = pd.DataFrame(viral_family_count.items(), columns=['Viral Family', 'Family Count'])
n_acid_df = pd.DataFrame(n_acid_count.items(), columns=['Nucleic Acid', 'Nucleic Acid Count'])

# Combine into a single DataFrame with potential mismatched rows
combined_df = pd.concat([domain_df, viral_family_df, n_acid_df], axis=1)

# Save to Excel

output_file = 'viral_data_counts.xlsx' # Replace with your output file path
combined_df.to_excel(output_file, index=False)


print(f"Data saved to {output_file}")

