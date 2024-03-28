'''script for extracting a specified sequence from fasta file.  Take fasta path, output path, seq id, start and end of desired sequence, saves to output file'''
import sys
from Bio import SeqIO

def extract_sequence(input_file, output_file, sequence_id, start, end):
    with open(input_file, "r") as in_file, open(output_file, "w") as out_file:
        for record in SeqIO.parse(in_file, "fasta"):
            if record.id == sequence_id:
                sequence = record.seq[start - 1:end]  # Adjust for 0-based indexing
                out_file.write(f">{record.id}_{start}-{end}\n")
                out_file.write(str(sequence) + "\n")

if __name__ == '__main__':
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    sequence_id = sys.argv[3]
    start = int(sys.argv[4])  # Start position (1-based)
    end = int(sys.argv[5])    # End position (1-based)
    extract_sequence(input_file, output_file, sequence_id, start, end)
#python3 pull_seq.py /scratch/general/nfs1/utu_4310/fr10_evolution_wd/fnas_dbsGCA_028564925.1_aRanSyl1.merge_genomic.fna test.txt CM053038.1 30256486 30260748
