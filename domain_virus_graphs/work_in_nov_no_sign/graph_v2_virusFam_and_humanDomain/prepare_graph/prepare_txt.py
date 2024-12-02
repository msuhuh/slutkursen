import pandas as pd

def prepare_txt(csv_path, out_path):
    data = pd.read_csv(csv_path)

    # Prepare the edge list
    edges = []

    # Iterate through rows (Domain) and columns (virus families)
    for _, row in data.iterrows():
        Domain = row["Domain"]  # Get the human domain
        Domain = Domain.replace(" ", "_")
        for virus_family in data.columns[1:]:  # Exclude the Domain column
            weight = row[virus_family]
            if weight != 0:  # Only include edges with a non 0 weight
                edges.append((Domain, virus_family, weight))

    # Convert the edges list to a DataFrame
    edge_list = pd.DataFrame(edges, columns=["source", "target", "weight"])

    # Save the edge list to a file in NCOL format (tab-separated)
    edge_list.to_csv(out_path, sep="\t", index=False, header=False)

    print("Edge list prepared and saved.")

input_mats = ["local_code\graph_v2_virusFam_and_humanDomain\graph_data\HumanDomain_virusFamily_matrix_norm_0_1.csv",
              "local_code\graph_v2_virusFam_and_humanDomain\graph_data\hit_matrix.csv"]
output_txts = ["local_code\graph_v2_virusFam_and_humanDomain\graph_data\HumanDomain_virusFamily_matrix_norm_0_1.txt",
               "local_code\graph_v2_virusFam_and_humanDomain\graph_data\hit_matrix.txt"]


for i in range(len(input_mats)):
    prepare_txt(input_mats[i], output_txts[i])
