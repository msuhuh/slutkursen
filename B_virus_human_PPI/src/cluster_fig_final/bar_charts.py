import os
import matplotlib.pyplot as plt
import re

# Custom sorting function to sort files numerically
def numeric_sort_key(path):
    # Extract the file name without the folder path
    file_name = os.path.basename(path)
    # Extract the numeric part using a regular expression
    numbers = re.findall(r'\d+', file_name)
    # Convert to integers for sorting, default to 0 if no numbers found
    return int(numbers[0]) if numbers else 0

def read_txt_file(file_path, use_counts=True):
    """
    Reads a .txt file with virus family names and optional counts.

    Parameters:
        file_path (str): Path to the .txt file.
        use_counts (bool): Whether to use counts from the file if provided.

    Returns:
        list: A list of tuples (family_name, count).
    """
    virus_list = []

    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()
            if "," in line and use_counts:
                family_name, count = line.split(",")
                virus_list.append((family_name.strip(), int(count.strip())))
            else:
                virus_list.append((line.strip(), 1))

    return virus_list

def map_viruses_to_classifications(virus_list, classification_mapping):
    """
    Maps virus families from the input list to their classifications.

    Parameters:
        virus_list (list): List of tuples (family_name, count).
        classification_mapping (dict): Mapping of virus families to classifications.

    Returns:
        dict: A dictionary with family names and their full classification details.
    """
    mapped_viruses = {}

    for family_name, count in virus_list:
        classification = classification_mapping.get(family_name, None)
        if classification:
            mapped_viruses[family_name] = {
                "count": count,
                **classification
            }
        else:
            mapped_viruses[family_name] = {
                "count": count,
                "type": "No match",
                "acid": "No match",
                "ds_ss": "No match",
                "plus_minus": "No match",
                "env": "No match"
            }

    return mapped_viruses

def create_classification_mapping():
    """
    Creates a mapping of virus families to their classifications and attributes.

    Returns:
        dict: A dictionary with family names as keys and attributes as values.
    """
    classification_mapping = {

        "Togaviridae": {"type": "(+) ssRNA env", "acid": "RNA", "ds_ss": "ss", "plus_minus": "+", "env": "Enveloped"},
        "Tobaniviridae": {"type": "(+) ssRNA env", "acid": "RNA", "ds_ss": "ss", "plus_minus": "+", "env": "Enveloped"},
        "Rhabdoviridae": {"type": "(-) ssRNA env", "acid": "RNA", "ds_ss": "ss", "plus_minus": "-", "env": "Enveloped"},
        "Retroviridae": {"type": "(+) ssRNA env", "acid": "RNA", "ds_ss": "ss", "plus_minus": "+", "env": "Enveloped"},
        "Reoviridae": {"type": "dsRNA non-env", "acid": "RNA", "ds_ss": "ds", "plus_minus": None, "env": "Non-enveloped"},
        "Poxviridae": {"type": "dsDNA env", "acid": "DNA", "ds_ss": "ds", "plus_minus": None, "env": "Enveloped"},
        "Polyomaviridae": {"type": "dsDNA non-env", "acid": "DNA", "ds_ss": "ds", "plus_minus": None, "env": "Non-enveloped"},
        "Pneumoviridae": {"type": "(-) ssRNA env", "acid": "RNA", "ds_ss": "ss", "plus_minus": "-", "env": "Enveloped"},
        "Picornaviridae": {"type": "(+) ssRNA non-env", "acid": "RNA", "ds_ss": "ss", "plus_minus": "+", "env": "Non-enveloped"},
        "Picobirnaviridae": {"type": "dsRNA non-env", "acid": "RNA", "ds_ss": "ds", "plus_minus": None, "env": "Non-enveloped"},
        "Phenuiviridae": {"type": "(-) ssRNA env", "acid": "RNA", "ds_ss": "ss", "plus_minus": "-", "env": "Enveloped"},
        "Peribunyaviridae": {"type": "(-) ssRNA env", "acid": "RNA", "ds_ss": "ss", "plus_minus": "-", "env": "Enveloped"},
        "Parvoviridae": {"type": "ssDNA non-env", "acid": "DNA", "ds_ss": "ss", "plus_minus": None, "env": "Non-enveloped"},
        "Paramyxoviridae": {"type": "(-) ssRNA env", "acid": "RNA", "ds_ss": "ss", "plus_minus": "-", "env": "Enveloped"},
        "Papillomaviridae": {"type": "dsDNA non-env", "acid": "DNA", "ds_ss": "ds", "plus_minus": None, "env": "Non-enveloped"},
        "Orthomyxoviridae": {"type": "(-) ssRNA env", "acid": "RNA", "ds_ss": "ss", "plus_minus": "-", "env": "Enveloped"},
        "Nodaviridae": {"type": "(+) ssRNA non-env", "acid": "RNA", "ds_ss": "ss", "plus_minus": "+", "env": "Non-enveloped"},
        "Nairoviridae": {"type": "(-) ssRNA env", "acid": "RNA", "ds_ss": "ss", "plus_minus": "-", "env": "Enveloped"},
        "Matonaviridae": {"type": "(+) ssRNA env", "acid": "RNA", "ds_ss": "ss", "plus_minus": "+", "env": "Enveloped"},
        "Kolmioviridae": {"type": "(-) ssRNA circular", "acid": "RNA", "ds_ss": "ss", "plus_minus": "-", "env": "Non-enveloped"},
        "Herpesviridae": {"type": "dsDNA env", "acid": "DNA", "ds_ss": "ds", "plus_minus": None, "env": "Enveloped"},
        "Hepeviridae": {"type": "(+) ssRNA non-env", "acid": "RNA", "ds_ss": "ss", "plus_minus": "+", "env": "Non-enveloped"},
        "Hepadnaviridae": {"type": "dsDNA env", "acid": "DNA", "ds_ss": "ds", "plus_minus": None, "env": "Enveloped"},
        "Hantaviridae": {"type": "(-) ssRNA env", "acid": "RNA", "ds_ss": "ss", "plus_minus": "-", "env": "Enveloped"},
        "Flaviviridae": {"type": "(+) ssRNA env", "acid": "RNA", "ds_ss": "ss", "plus_minus": "+", "env": "Enveloped"},
        "Filoviridae": {"type": "(-) ssRNA env", "acid": "RNA", "ds_ss": "ss", "plus_minus": "-", "env": "Enveloped"},
        "Coronaviridae": {"type": "(+) ssRNA env", "acid": "RNA", "ds_ss": "ss", "plus_minus": "+", "env": "Enveloped"},
        "Caliciviridae": {"type": "(+) ssRNA non-env", "acid": "RNA", "ds_ss": "ss", "plus_minus": "+", "env": "Non-enveloped"},
        "Bornaviridae": {"type": "(-) ssRNA env", "acid": "RNA", "ds_ss": "ss", "plus_minus": "-", "env": "Enveloped"},
        "Birnaviridae": {"type": "dsRNA non-env", "acid": "RNA", "ds_ss": "ds", "plus_minus": None, "env": "Non-enveloped"},
        "Astroviridae": {"type": "(+) ssRNA non-env", "acid": "RNA", "ds_ss": "ss", "plus_minus": "+", "env": "Non-enveloped"},
        "Arteriviridae": {"type": "(+) ssRNA env", "acid": "RNA", "ds_ss": "ss", "plus_minus": "+", "env": "Enveloped"},
        "Arenaviridae": {"type": "(-) ssRNA env", "acid": "RNA", "ds_ss": "ss", "plus_minus": "-", "env": "Enveloped"},
        "Anelloviridae": {"type": "ssDNA non-env", "acid": "DNA", "ds_ss": "ss", "plus_minus": None, "env": "Non-enveloped"},
        "Adenoviridae": {"type": "dsDNA non-env", "acid": "DNA", "ds_ss": "ds", "plus_minus": None, "env": "Non-enveloped"}
    }
    return classification_mapping

def get_color_acid(acid_type):
    print(acid_type)
    "Return hexcode based on the acid type"
    assert (acid_type=="DNA" or acid_type=="RNA"), f"{acid_type} was sent to get_color_acid. Only DNA or RNA are accepted"

    if acid_type == "DNA":
        return "FFCCB6"
    else:
        return "A2E1DB"


def get_color_env(env_type):
    print(env_type)
    "Return hexcode based on the acid type"
    assert (env_type=="Enveloped" or env_type=="Non-enveloped"), f"{env_type} was sent to get_color_acid. Only DNA or RNA are accepted"

    if env_type == "Enveloped":
        return "AACCB6"
    else:
        return "FFE1DF"



def get_virus_type(virus_type):
    """
    No asserting here, too lazy.
    """
    color = None

    # Virus types    

    if virus_type == "dsDNA non-env":
        color = "B22222"  # Crimson Red

    elif virus_type == "(-) ssRNA env":
        color = "228B22"  # Forest Green

    elif virus_type == "dsDNA env":
        color = "000000"  # Black

    elif virus_type == "dsRNA non-env":
        color = "9932CC"  # Dark Orchid

    elif virus_type == "(+) ssRNA env":
        color = "8B4513"  # Saddle Brown

    elif virus_type == "(+) ssRNA non-env":
        color = "FF4500"  # Orange Red (distinct and vibrant)

    elif virus_type == "(-) ssRNA circular":
        color = "DAA520"  # Goldenrod

    elif virus_type == "ssDNA non-env":
        color = "4682B4"  # Steel Blue

    else:
        print(f"Virus type not found. Got: {virus_type}, returning 'None'. ")

    return color


def create_multi_row_pie_charts(txt_paths, classification_mapping, use_counts):
    """
    Generate rows of pie charts for multiple input files (clusters), compactly displaying the data.

    Parameters:
        txt_paths (list): List of file paths containing virus data.
        classification_mapping (dict): Dictionary mapping virus families to classifications.
        use_counts (bool): Whether to use counts from the input files.
    """
    # Prepare a list to store data for each file
    rows_data = []

    for file_path in txt_paths:
        # Read the input file
        virus_list = read_txt_file(file_path, use_counts)

        # Map the viruses to classifications
        mapped_viruses = map_viruses_to_classifications(virus_list, classification_mapping)

        # Prepare data for pie charts
        acid_types = {}
        virus_types = {}
        env_types = {}
        family_names = []
        total_count = 0

        for family, details in mapped_viruses.items():
            acid_label = details["acid"]
            virus_type_label = details["type"]
            env_label = details["env"]

            family_names.append(family)
            total_count += details["count"]

            # Acid type distribution
            acid_types[acid_label] = acid_types.get(acid_label, 0) + details["count"]

            # Virus type distribution
            virus_types[virus_type_label] = virus_types.get(virus_type_label, 0) + details["count"]

            # Envelope status distribution
            env_types[env_label] = env_types.get(env_label, 0) + details["count"]

        # Save the data for this file
        cluster_name = os.path.splitext(os.path.basename(file_path))[0]  # Extract cluster name
        rows_data.append((cluster_name, family_names, virus_types, acid_types, env_types, total_count))

    # Create the main figure
    num_rows = len(rows_data)
    fig_height = max(2.5 * num_rows, 12)  # Adjust figure height dynamically
    #print(fig_height)
    #a=b
    fig, axes = plt.subplots(
        num_rows, 5, figsize=(12, fig_height),  # Add an extra column
        gridspec_kw={'width_ratios': [0, 1, 1, 1, 0]}  # Narrow column for N=N
    )

    # Handle single-row case
    if num_rows == 1:
        axes = [axes]

    
    # Populate each row
    for row_idx, (cluster_name, family_names, virus_types, acid_types, env_types, total_count) in enumerate(rows_data):
       

        # Bar chart data
        for ax, (title, data, color_func) in zip(
            axes[row_idx][1:4],  # Columns 1 to 3 for bar charts
            [
                ("Virus type distribution", virus_types, get_virus_type),
                ("Nucleic acid", acid_types, get_color_acid),
                ("Enveloped,\nNon-enveloped", env_types, get_color_env),
            ],
        ):
            if not data:
                ax.text(0.5, 0.5, "No Data", ha='center', va='center', fontsize=12)
                ax.axis('off')
                continue

            # Extract keys and values for bar chart
            keys = list(data.keys())
            values = list(data.values())

            # Generate colors for bars
            colors = [f"#{color_func(key)}" for key in keys]

            # Create the bar chart
            ax.bar(keys, values, color=colors)

            # Set chart title and labels
            ax.set_title(title, fontsize=12)
            #ax.set_xlabel('Categories', fontsize=10)
            ax.set_ylabel('Count', fontsize=12)

            # Adjust label rotation for better readability
            ax.tick_params(axis='x', labelsize=0, rotation=90)
            ax.tick_params(axis='y', labelsize=10)

    # Adjust layout for sufficient spacing
    fig.subplots_adjust(
    top=0.995,
    bottom=0.005,
    left=0.05,  # Move left boundary closer to the charts
    right=0.995,
    hspace=0.064,
    wspace=0.0  # Minimize column spacing
    )

    plt.tight_layout()
    #plt.show()
    plt.savefig('output.png', bbox_inches='tight')
    #plt.savefig()

# Example usage

input_folder = "cluster_fig_final/input_folder_all_virus_types"
use_counts = True
txt_paths = []


# Iterate through subfolders and files
for _, _, files in os.walk(input_folder):
    for file in files:
        if file.endswith(".txt"):  # Check if the file has a .txt extension
            file_path = os.path.join(input_folder, file)
            txt_paths.append(file_path)  # Append subfolder and file path


# Sort the paths numerically
txt_paths.sort(key=numeric_sort_key)
# Create the classification mapping
classification_mapping = create_classification_mapping()
create_multi_row_pie_charts(txt_paths, classification_mapping, use_counts)