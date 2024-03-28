from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import sys

def calculate_percent_identity(alignment):
    total_positions = alignment.get_alignment_length()
    num_sequences = len(alignment)

    identical_positions = 0
    for i in range(total_positions):
        column = alignment[:, i]
        if len(set(column)) == 1 and '-' not in column:
            identical_positions += 1

    percent_identity = (identical_positions / total_positions) * 100
    return percent_identity

def main():
    # Replace 'your_alignment_file.fasta' with the path to your multiple sequence alignment file
    alignment_file = sys.argv[1]

    try:
        alignment = AlignIO.read(alignment_file, 'fasta')
    except FileNotFoundError:
        print(f"Error: The file '{alignment_file}' was not found.")
        return
    except Exception as e:
        print(f"An error occurred: {e}")
        return

    percent_identity = calculate_percent_identity(alignment)
    print(f"Percent Identity: {percent_identity:.2f}%")

if __name__ == "__main__":
    main()

