import pandas as pd
from Bio import Entrez, SeqIO
import csv

# Set your email for NCBI API (required by NCBI)
Entrez.email = "alexander.blomlof@hotmail.com"

def fetch_proteome(accession_number):
    """
    Fetch proteome sequence from NCBI using accession number.
    """
    try:
        handle = Entrez.efetch(db="protein", id=accession_number, rettype="fasta", retmode="text")
        sequences = list(SeqIO.parse(handle, "fasta"))
        if not sequences:
            raise ValueError("No sequences found.")
        return sequences
    except Exception as e:
        print(f"Failed to retrieve {accession_number}: {str(e)}")
        return None

def calculate_proteome_length(sequences):
    """
    Calculate total proteome length from a list of sequences.
    """
    return sum(len(seq.seq) for seq in sequences)

def process_accessions(file_path):
    """
    Process accession numbers from an Excel file to fetch proteomes and calculate their lengths.
    """
    # Load Excel file and extract unique accession numbers
    df = pd.read_excel(file_path)
    unique_accessions = df['virus_accession'].drop_duplicates().tolist()

    proteome_lengths = {}
    for accession in unique_accessions:
        print(f"Processing accession: {accession}")
        sequences = fetch_proteome(accession)
        if sequences:
            total_length = calculate_proteome_length(sequences)
            proteome_lengths[accession] = total_length
        else:
            proteome_lengths[accession] = "Failed to retrieve"

    return proteome_lengths

def save_to_csv(results, output_file):
    """
    Save results to a CSV file.
    """
    with open(output_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Accession Number", "Proteome Length"])
        for accession, length in results.items():
            writer.writerow([accession, length])

# File path to the Excel file
excel_file = "local_code\data\data_collapsed_new.xlsx"  # Replace with your file path
output_csv = "local_code\data\proteome_lengths.csv"

# Process the file and fetch proteome lengths
results = process_accessions(excel_file)

# Save results to a CSV file
save_to_csv(results, output_csv)

print(f"Proteome lengths saved to {output_csv}.")
