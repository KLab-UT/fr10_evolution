import argparse
from Bio import Entrez

def download_sequences(email, database, start_id, end_id, output_file):
    """
    Download sequences from NCBI within the specified ID range and save them to a file.

    Parameters:
    - email (str): Your email address for NCBI.
    - database (str): NCBI database name (e.g., "nuccore" for nucleotide sequences).
    - start_id (int): Starting ID of the sequence range.
    - end_id (int): Ending ID of the sequence range.
    - output_file (str): Output file to save the downloaded sequences.

    Returns:
    None
    """
    # Provide your email address to NCBI
    Entrez.email = email

    # Generate a list of IDs within the specified range
    id_range = list(range(start_id, end_id + 1))

    # Convert the list of IDs to a comma-separated string
    id_list = ",".join(map(str, id_range))

    try:
        # Fetch the sequences from NCBI
        handle = Entrez.efetch(db=database, id=id_list, rettype="fasta", retmode="text")
        sequences = handle.read()
        handle.close()

        # Save the sequences to the specified output file
        with open(output_file, "w") as file:
            file.write(sequences)

        print(f"Sequences downloaded successfully and saved to {output_file}")

    except Exception as e:
        print(f"An error occurred: {e}")

def parse_arguments():
    parser = argparse.ArgumentParser(description="Download sequences from NCBI within a specified ID range.")
    parser.add_argument("email", help="Your email address for NCBI.")
    parser.add_argument("database", help="NCBI database name (e.g., 'nuccore' for nucleotide sequences).")
    parser.add_argument("start_id", type=int, help="Starting ID of the sequence range.")
    parser.add_argument("end_id", type=int, help="Ending ID of the sequence range.")
    parser.add_argument("output_file", help="Output file to save the downloaded sequences.")
    return parser.parse_args()

if __name__ == "__main__":
    # Example usage:
    # python download_sequences.py your_email@example.com nuccore 1 10 downloaded_sequences.fasta
    args = parse_arguments()
    download_sequences(args.email, args.database, args.start_id, args.end_id, args.output_file)

