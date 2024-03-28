from Bio import SeqIO, Entrez
import json
import argparse
import subprocess
import os

def make_amino_acid_alignment(gene, ids, output_dir):
    """
    make amino acid alignment named after gene by fetching all of the sequences identified in the ids list
    """
    seqs = []
    problem_ids =[]
    for seq_id in ids:
        try:
            handle = Entrez.efetch(db="protein", id=seq_id, rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            seqs.append(record)
        except:
            print(f"problem fetching: {seq_id}")
            problem_ids.append(seq_id)

    with open("problem_seqs.txt", "a") as problems:
        for id in problem_ids:
            problems.write(f"{id}\n")

    fout = os.path.join(output_dir, f"{gene}.fa")
    SeqIO.write(seqs, fout, "fasta")

def align_with_query(gene_fasta, query_fasta, output_dir):
    """
    Align sequences from query FASTA to gene FASTA using MAFFT and store the result in output directory.

    Args:
        gene_fasta (str): Path to the original gene FASTA file.
        query_fasta (str): Path to the query FASTA file containing sequences to align.
        output_dir (str): Directory to store the aligned sequences.

    Returns:
        None
    """
    # Read gene FASTA file
    gene_sequences = SeqIO.to_dict(SeqIO.parse(gene_fasta, "fasta"))

    # Read query sequences
    query_sequences = SeqIO.parse(query_fasta, "fasta")

    # Combine gene and query sequences
    for query_record in query_sequences:
        gene_sequences[query_record.id] = query_record

    # Write combined sequences to a temporary FASTA file
    temp_combined_fasta = os.path.join(output_dir, "temp.combined.fasta")
    with open(temp_combined_fasta, "w") as temp_handle:
        SeqIO.write(gene_sequences.values(), temp_handle, "fasta")

    # Call MAFFT to align sequences and write to output directory
    aligned_fasta = os.path.join(output_dir, f"{os.path.splitext(os.path.basename(gene_fasta))[0]}_aligned.fa")
    command = f"mafft --auto {temp_combined_fasta} > {aligned_fasta}"
    subprocess.run(command, shell=True, check=True)

    estimate_tree(aligned_fasta)

def estimate_tree(alignment_fasta):
    """ 
    Use iqtree to estimate a phylogeny and update the tree descriptions
    """

    # Estimate tree
    iqtree_command = f"iqtree2 -s {alignment_fasta} -B 1000 -redo"
    subprocess.run(iqtree_command, shell=True, check=True)
    contree_path = f"{alignment_fasta}.contree"

    # Update descriptions
    desc_tree_path = f"{contree_path.split('.')[0]}.contree"
    decriptions_command = f"python3 get_tree_descriptions.py {contree_path} {alignment_fasta} {desc_tree_path}"
    subprocess.run(decriptions_command, shell=True, check=True)

def parse_genes(gene_ids_json):
    genes_dict = json.load(gene_ids_json)
    return genes_dict

def main():
    parser = argparse.ArgumentParser(description="Fetch protein sequences and create amino acid alignment.")
    parser.add_argument("gene_json", help="JSON file containing gene information and sequence IDs.")
    parser.add_argument("query_fasta", help="fasta containg query sequences to add to alignments")
    parser.add_argument("email", help="email to use NCBI Entrez")
    parser.add_argument("output_dir", help="directory to save output files")
    args = parser.parse_args()

    Entrez.email = args.email

    with open(args.gene_json, "r") as gene_ids_json:
        genes_dict = parse_genes(gene_ids_json)

    # Create the output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    for gene, ids in genes_dict.items():
        print(f"Gene: {gene} ids: {ids}")
        make_amino_acid_alignment(gene, ids, args.output_dir)

        gene_fasta = os.path.join(args.output_dir, f"{gene}.fa")
        align_with_query(gene_fasta, args.query_fasta, args.output_dir)
        print(f"Alignment created!: {os.path.join(args.output_dir, f'{gene}_aligned.fa')}")

if __name__ == "__main__":
    main()
