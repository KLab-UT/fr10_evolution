import argparse
from Bio import SeqIO
from collections import Counter

def get_hit_consensus(hits_fasta, query_fasta, out_file):
    '''
    Create a consensus sequence with positions that agree in hits_fasta alignment, positions that disagree are used from the query provided that the nucleotide is found in at least one of the hit seqs at that position
    '''
    query_record = SeqIO.read(query_fasta, "fasta")
    hit_records = list(SeqIO.parse(hits_fasta, "fasta"))

    query_seq = str(query_record.seq)
    hit_seqs = [str(record.seq) for record in hit_records]

    hit_consensus = ""
    for i in range(len(hit_seqs[0])):
        # make a list of the nucleotide at position i for each hit sequence
        nts = [hit[i] for hit in hit_seqs if hit[i] != "-"]

        if not nts:  # all "-" at this position
            hit_consensus += "-"  # or any other default nucleotide you want
        elif all(nt == nts[0] for nt in nts if nt != "-"):  # check if nucleotides agree, ignoring "-"
            hit_consensus += nts[0]
        elif "-" in nts:  # check if there is at least one non-"-" nucleotide
            non_dash_nts = [nt for nt in nts if nt != "-"]
            hit_consensus += Counter(non_dash_nts).most_common(1)[0][0]
        elif query_seq[i] in nts:  # check if query nucleotide is in nts
            hit_consensus += query_seq[i]
        else:  # nucleotides disagree
            most_common_hit_nt = Counter(nts).most_common(1)[0][0]
            hit_consensus += most_common_hit_nt

    # Write the hit_consensus sequence to the output file
    with open(out_file, "w") as output_handle:
        output_handle.write(f">{hit_records[0].id}_consensus\n")
        output_handle.write(hit_consensus)

def main():
    parser = argparse.ArgumentParser(description="Generate a hit consensus sequence.")
    parser.add_argument("hits_fasta", help="Path to the hits fasta file")
    parser.add_argument("query_fasta", help="Path to the query fasta file")
    parser.add_argument("out_file", help="Path to the output consensus fasta file")

    args = parser.parse_args()
    get_hit_consensus(args.hits_fasta, args.query_fasta, args.out_file)
    print(f"Consensus sequence written to {args.out_file}.")

if __name__ == "__main__":
    main()
