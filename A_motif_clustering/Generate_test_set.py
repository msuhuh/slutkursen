import pandas as pd
 
input_path = '20241122_motif_mapped_data.xlsx' # Replace with your path
output_path = '202041125_test_set.xlsx' # Replace with your path 

# Function to filter the data and extract the values with highest confidence
def filter_excel_file(input_file, output_file):
    try:
        # Load the Excel file
        df = pd.read_excel(input_file)

        # Filter rows based on the conditions
        condition = (df['max_confidence'] >= 3) & (df['new_motif'] != '*') & (df['virus_motif_instances'].isna() | (df['virus_motif_instances'] == ''))
        
        # Combine the two conditions
        filtered_df = df[condition]

        # Save the filtered rows to a new Excel file
        filtered_df.to_excel(output_file, index=False)

        print(f"Filtered data saved to {output_file}")

    except Exception as e:
        print(f"An error occurred: {e}")

filter_excel_file(input_path, output_path)
