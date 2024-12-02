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

def process_accessions(new_accessions):
    """
    Process new accession numbers to fetch proteomes and calculate their lengths.
    """
    proteome_lengths = {}
    for accession in new_accessions:
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
    with open(output_file, "a", newline="") as csvfile:  # Append mode
        writer = csv.writer(csvfile)
        for accession, length in results.items():
            writer.writerow([accession, length])

# File paths
excel_file = "local_code\data\Compiled_viral_results_compiled.xlsx"  
existing_csv = "local_code\data\proteome_lengths_v2.csv"  # Existing CSV with processed results
output_csv = existing_csv  # Output will overwrite the existing CSV

# Load new accession numbers from Excel
df_excel = pd.read_excel(excel_file)
new_accessions = set(df_excel['virus_accession'].drop_duplicates())

# Load existing accession numbers from CSV
try:
    df_existing = pd.read_csv(existing_csv)
    existing_accessions = set(df_existing['Accession Number'])
except FileNotFoundError:
    # If the CSV does not exist, there are no existing accessions
    existing_accessions = set()

# Filter out already processed accession numbers
accessions_to_process = new_accessions - existing_accessions

print(f"Total new accessions to process: {len(accessions_to_process)}")

# Process only new accession numbers
if accessions_to_process:
    results = process_accessions(accessions_to_process)
    save_to_csv(results, output_csv)
else:
    print("No new accession numbers to process.")

print(f"Proteome lengths updated in {output_csv}.")
