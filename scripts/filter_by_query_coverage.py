import os
import argparse
from Bio import Entrez, SeqIO

def fetch_rna_sequences(term, taxid, email, output_dir, retmax=500):
    """
    Fetch RNA sequences from NCBI based on a given term and TaxID.

    Parameters:
    - term (str): Search term for the gene or feature of interest.
    - taxid (str): NCBI Taxonomy ID for the desired organism.
    - email (str): Email address for identification to NCBI.
    - output_dir (str): Directory to save downloaded sequences.
    - retmax (int): Maximum number of records to retrieve per request.

    Returns:
    None
    """
    Entrez.email = email

    # Step 1: Search for RNA sequences using the specified term and taxid
    search_query=f'txid{taxid}[Organism] AND {term}[All Fields]'
    search_handle = Entrez.esearch(db="gene", term=search_query, retmax=retmax)
    search_results = Entrez.read(search_handle)
    id_list = search_results["IdList"]

    if not id_list:
        print(f"No RNA sequences found for the term '{term}' in TaxID {taxid}")
        return

    # Step 2: Fetch and save the RNA sequences
    for sequence_id in id_list:
        try:
            # Fetch the sequence
            seq_handle = Entrez.efetch(db="nucleotide", id=sequence_id, rettype="fasta", retmode="text")
            sequence = SeqIO.read(seq_handle, "fasta")

            # Save the sequence to a file
            output_file = os.path.join(output_dir, f"{sequence_id}_{sequence.description}.fasta")
            SeqIO.write(sequence, output_file, "fasta")
            print(f"Downloaded {sequence_id} - {sequence.description}")

        except Exception as e:
            print(f"Error fetching RNA sequence {sequence_id}: {str(e)}")

def parse_arguments():
    """
    Parse command-line arguments.

    Returns:
    argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(description='Fetch RNA sequences from NCBI.')
    parser.add_argument('--term', type=str, required=True, help='Search term for the gene or feature of interest.')
    parser.add_argument('--taxid', type=str, required=True, help='NCBI Taxonomy ID for the desired organism.')
    parser.add_argument('--email', type=str, required=True, help='Your email address for identification to NCBI.')
    parser.add_argument('--output_dir', type=str, default='output_sequences', help='Directory to save downloaded sequences.')
    parser.add_argument('--retmax', type=int, default=500, help='Maximum number of records to retrieve per request.')
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    fetch_rna_sequences(args.term, args.taxid, args.email, args.output_dir, args.retmax)

