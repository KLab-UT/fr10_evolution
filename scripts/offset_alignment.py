import argparse
from Bio import SeqIO

def add_gaps_to_sequences(input_file, output_file, num_gaps):
    sequences = []

    # Read the input FASTA file
    with open(input_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            # Add gaps to the start of each sequence
            new_sequence = "-" * num_gaps + str(record.seq)
            sequences.append((record.id, new_sequence))

    # Write the modified sequences to the output FASTA file
    with open(output_file, "w") as handle:
        for seq_id, seq in sequences:
            handle.write(f">{seq_id}\n{seq}\n")

def main():
    parser = argparse.ArgumentParser(description="Add a specified number of gaps to the start of every sequence in a FASTA file.")
    parser.add_argument("input_file", help="Path to the input FASTA file")
    parser.add_argument("output_file", help="Path to the output FASTA file")
    parser.add_argument("num_gaps", type=int, help="Number of gaps to add to the start of each sequence")

    args = parser.parse_args()

    add_gaps_to_sequences(args.input_file, args.output_file, args.num_gaps)
    print(f"Gaps added to sequences. Output written to {args.output_file}")

if __name__ == "__main__":
    main()

