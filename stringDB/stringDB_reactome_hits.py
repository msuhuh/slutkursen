import requests  # python -m pip install requests
import json
import os
import time

# Define folders
"""
intermediate_folder --> store raw hits from StringDB
res_folder --> Filtered hits on reactome and p_val < 0.01

input_folder --> root folder with subfolders in it and txt in subfolders.
Ex: root;
    subfolderA;
        file1.txt --> where human accession is seperated by new row 
        file2.txt

    SubfolderB;
        fileX.txt
        fileY.txt

OUTPUT: 
file output will be named "subfolder_filename_filtered.json
file output can be changed in line 93
"""
intermediate_folder = "stringDB_api/intermediate_jsons"
res_folder = "stringDB_api/RCTM_filtered_jsons"
input_folder = "using_bait_fam_matrix/txt_with_res"

# Ensure output folders exist
os.makedirs(intermediate_folder, exist_ok=True)
os.makedirs(res_folder, exist_ok=True)

#################################################

"""
API
"""
string_api_url = "https://version-12-0.string-db.org/api"
output_format = "json"
method = "enrichment"
request_url = "/".join([string_api_url, output_format, method])

#################################################

"""
Initiate paths
"""
# Initialize list for storing text file paths
txt_paths = []

# Iterate through subfolders and files
for subfolder, _, files in os.walk(input_folder):
    for file in files:
        if file.endswith(".txt"):  # Check if the file has a .txt extension
            file_path = os.path.join(subfolder, file)
            txt_paths.append((subfolder, file_path))  # Append subfolder and file path

####################################################

"""
Start StringDB retrival
"""

# Start loop per input
for subfolder, txt in txt_paths:
    time.sleep(1) # Avoid crashing the website
    
    human_accessions = []
    with open(txt, 'r') as file:
        # Read the entire file content, split by newlines, and store it as a list
        file_content = file.read().splitlines()
        human_accessions.extend(file_content)  # Add lines to the list (flatten structure)

    if len(human_accessions) > 0: # Ensure some data actually exists
        # Print the resulting list
        print(f"Processing file: {txt}")
        #print(human_accessions)

        # Prepare parameters for STRING API
        params = {
            "identifiers": "%0d".join(human_accessions),  # your protein
            "species": 9606,  # species NCBI identifier
            "caller_identity": "www.awesome_app.org",  # your app name --> StringDB example
        }

        # Call STRING API
        response = requests.post(request_url, data=params)
        data = json.loads(response.text)

        # Generate unique output filenames
        subfolder_name = os.path.basename(subfolder)  # Get the name of the subfolder
        base_filename = os.path.basename(txt).replace(".txt", "")  # Strip the extension
        unique_filename = f"{subfolder_name}_{base_filename}"

        # Save intermediate JSON
        intermediate_json = os.path.join(intermediate_folder, f"{unique_filename}_intermediate.json")
        with open(intermediate_json, 'w') as json_file:
            json.dump(data, json_file, indent=4)  # Save intermediate JSON
        #print(f"Intermediate JSON data has been saved to {intermediate_json}")

        # Filter results
        filtered_data = [row for row in data if row["category"] == "RCTM" and float(row["p_value"]) < 0.01]

        # Data without significant hits will result in empty json files.
        # Save filtered JSON
        res_json = os.path.join(res_folder, f"{unique_filename}_filtered.json")
        with open(res_json, 'w') as json_file:
            json.dump(filtered_data, json_file, indent=4)  # Save filtered JSON
        #print(f"Filtered JSON data has been saved to {res_json}")