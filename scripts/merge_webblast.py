from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.Align.AlignInfo import SummaryInfo
from Bio.SeqRecord import SeqRecord
import sys
import argparse

def consensus_by_species(input_file, output_file):
    """
    Generate a consensus sequence for each species in the input FASTA file.

    Parameters:
        input_file (str): Path to the input FASTA file containing sequences.
        output_file (str): Path to the output FASTA file to store consensus sequences.

    Returns:
        None
    """
    records = list(SeqIO.parse(input_file, 'fasta'))
    species_records = {}

    for record in records:
        seq_species = " ".join(record.description.split()[-2:])
        seq_range = record.description.split()[1]
        if seq_species not in species_records:
            species_records[seq_species] = [record]
            print(seq_species)
        else:
            species_records[seq_species].append(record)

    merged_records = []
    for seq_species, record_list in species_records.items():
        temp_input_file = f"temp_input_{seq_species}.fasta"
        SeqIO.write(record_list, temp_input_file, "fasta")
        alignment = AlignIO.read(temp_input_file, "fasta")

        summary_info = SummaryInfo(alignment)
        consensus_sequence = summary_info.dumb_consensus(threshold=0.51)
        species_ids = []
        for record in record_list:
            species_ids.append(record.id)

        consensus_id = "+".join(species_ids)
        print(f"consensus id: {consensus_id}")
        consensus_description = "| " + get_ranges(record_list) + " | " + seq_species
        consensus_record = SeqRecord(consensus_sequence, id=consensus_id, description=consensus_description)
        merged_records.append(consensus_record)

    SeqIO.write(merged_records, output_file, 'fasta')
    replace_x_with_gap(output_file, output_file)

def get_ranges(record_list):
    """
    return str the ranges of each record in the list where the range comes after :
    in the seqID
    """
    ranges=[]
    for record in record_list:
        seq_range = record.description.split()[1]
        ranges.append(seq_range)
    return ":".join(ranges)


def replace_x_with_gap(input_file, output_file):
    """
    Replace 'X' characters with '-' in sequences and write to a new FASTA file.

    Parameters:
        input_file (str): Path to the input FASTA file.
        output_file (str): Path to the output FASTA file with 'X' replaced by '-'.

    Returns:
        None
    """
    records = list(SeqIO.parse(input_file, "fasta"))
    for record in records:
        modified_seq = str(record.seq).replace("X", "-")
        record.seq = Seq(modified_seq)

    with open(output_file, 'w') as output_handle:
        SeqIO.write(records, output_handle, "fasta")


def consensus_by_id(input_file, output_file):
    """
    Generate a consensus sequence for each unique sequence ID in the input FASTA file.

    Parameters:
        input_file (str): Path to the input FASTA file containing sequences.
        output_file (str): Path to the output FASTA file to store consensus sequences.

    Returns:
        None
    """
    records = list(SeqIO.parse(input_file, 'fasta'))
    id_records = {}

    for record in records:
        seq_id = record.id
        seq_range = record.description.split()[1]
        if seq_id not in id_records:
            id_records[seq_id] = [record]
            print(seq_id)
        else:
            id_records[seq_id].append(record)

    merged_records = []
    for seq_id, record_list in id_records.items():
        temp_input_file = f"temp_input_{seq_id}.fasta"
        SeqIO.write(record_list, temp_input_file, "fasta")
        alignment = AlignIO.read(temp_input_file, "fasta")

        summary_info = SummaryInfo(alignment)
        consensus_sequence = summary_info.dumb_consensus(threshold=0.51)

        species = " ".join(record_list[0].description.split()[3:5])
        print(species)
        consensus_id = seq_id
        consensus_description = "| " + get_ranges(record_list) + " | " + species
        consensus_record = SeqRecord(consensus_sequence, id=consensus_id, description=consensus_description)
        merged_records.append(consensus_record)

    SeqIO.write(merged_records, output_file, 'fasta')
    replace_x_with_gap(output_file, output_file)

def get_ranges(record_list):
    """
    return str the ranges of each record in the list where the range comes after :
    in the seqID
    """
    ranges=[]
    for record in record_list:
        seq_range = record.description.split()[1]
        ranges.append(seq_range)
    return ":".join(ranges)


def replace_x_with_gap(input_file, output_file):
    """
    Replace 'X' characters with '-' in sequences and write to a new FASTA file.

    Parameters:
        input_file (str): Path to the input FASTA file.
        output_file (str): Path to the output FASTA file with 'X' replaced by '-'.

    Returns:
        None
    """
    records = list(SeqIO.parse(input_file, "fasta"))
    for record in records:
        modified_seq = str(record.seq).replace("X", "-")
        record.seq = Seq(modified_seq)

    with open(output_file, 'w') as output_handle:
        SeqIO.write(records, output_handle, "fasta")


def pad_seqs(input_file, output_file):
    """
    Pad sequences with '-' to have equal length and write to a new FASTA file.

    Parameters:
        input_file (str): Path to the input FASTA file.
        output_file (str): Path to the output FASTA file with padded sequences.

    Returns:
        None
    """
    records = list(SeqIO.parse(input_file, "fasta"))
    max_length = max(len(record.seq) for record in records)
    print(f"max seq length: {max_length}")
    padded_records = []
    for record in records:
        padding_length = max_length - len(record.seq)
        padded_seq = record.seq + Seq('-' * padding_length)
        padded_record = SeqRecord(padded_seq, id=record.id, description=record.description)
        padded_records.append(padded_record)

    SeqIO.write(padded_records, output_file, "fasta")


def remove_query(input_file, query_id):
    """
    Remove sequences containing the specified query ID and write to the same FASTA file.

    Parameters:
        input_file (str): Path to the input FASTA file.
        query_id (str): Identifier for the query sequence.

    Returns:
        None
    """
    records_to_write = []
    records = list(SeqIO.parse(input_file, 'fasta'))
    for record in records:
        if query_id not in record.id:
            records_to_write.append(record)
    SeqIO.write(records_to_write, input_file, 'fasta')


def main():
    parser = argparse.ArgumentParser(description="Generate consensus sequences from a FASTA file.")
    parser.add_argument("input_file", help="Path to the input FASTA file containing sequences.")
    parser.add_argument("output_file", help="Path to the output FASTA file to store consensus sequences.")
    parser.add_argument("--query_id", help="Identifier for the query sequence.")
    parser.add_argument("merge_by", choices=["id", "species"], help="Method for merging sequences.")

    args = parser.parse_args()

    padded_output_file = args.input_file.rsplit('.')[0] + '_padded.fasta'
    pad_seqs(args.input_file, padded_output_file)

    if args.query_id:
        remove_query(padded_output_file, args.query_id)

    if args.merge_by == "id":
        consensus_by_id(padded_output_file, args.output_file)
    else:
        consensus_by_species(padded_output_file, args.output_file)

if __name__ == '__main__':
    main()
