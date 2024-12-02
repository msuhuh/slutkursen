import os
import pandas as pd
from Bio import Entrez, SeqIO

# Set your email for NCBI API (required by NCBI)
Entrez.email = "alexander.blomlof@hotmail.com"

def fetch_and_save_fasta(accession_number, output_folder):
    """
    Fetch proteome sequence from NCBI using accession number and save it as a FASTA file.
    """
    try:
        handle = Entrez.efetch(db="protein", id=accession_number, rettype="fasta", retmode="text")
        sequences = list(SeqIO.parse(handle, "fasta"))
        if not sequences:
            raise ValueError("No sequences found.")
        
        # Save sequences to a FASTA file in the specified folder
        fasta_path = os.path.join(output_folder, f"{accession_number}.fasta")
        with open(fasta_path, "w") as fasta_file:
            SeqIO.write(sequences, fasta_file, "fasta")
        
        print(f"Saved {accession_number} to {fasta_path}")
    except Exception as e:
        print(f"Failed to retrieve {accession_number}: {str(e)}")

def process_accessions(file_path, output_folder):
    """
    Process accession numbers from an Excel file and fetch FASTA files.
    """
    # Create output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)

    # Load Excel file and extract unique accession numbers
    df = pd.read_excel(file_path)
    unique_accessions = df['virus_accession'].drop_duplicates().tolist()

    for accession in unique_accessions:
        print(f"Processing accession: {accession}")
        fetch_and_save_fasta(accession, output_folder)

# File paths
excel_file = "local_code/data/data_collapsed_new.xlsx"  # Replace with your file path
output_folder = "local_code/data/fasta_files"  # Folder to store FASTA files

# Process the file and fetch FASTA files
process_accessions(excel_file, output_folder)

print(f"FASTA files have been saved to {output_folder}.")
