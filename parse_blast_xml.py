"""parse_blast_xml.py
processes complete sequences from webblast xml results intedted to be used with blastn and tblastn for creating nucleuotide alignments
"""
import argparse
import xml.etree.ElementTree as ET
import subprocess
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import Counter

class BlastHit:
    """
    Manage blast hits from ET element
    """
    def __init__(self, hit_element, blast_type):
        self.accession = hit_element.find("Hit_accession").text
        self.definition = hit_element.find("Hit_def").text
        self.length = int(hit_element.find("Hit_len").text)

        if "TSA" in self.definition:
            self.species = " ".join(self.definition.split()[:3])
        else:
            self.species = " ".join(self.definition.split()[:2])

        print(f"species for {self.accession}: {self.species}")
        self.hsps = []

        # Iterate over Hsp elements and create Hsp objects
        for hsp_element in hit_element.findall("Hit_hsps/Hsp"):
            if blast_type == "blastn":
                hsp = self.create_blastn_hsp_object(hsp_element)
            elif blast_type == "tblastn":
                hsp = self.create_tblastn_hsp_object(hsp_element)
            elif blast_type == "tblastx":
                hsp = self.create_tblastx_hsp_object(hsp_element)
            else:
                raise ValueError(f"Unsupported blast type: {blast_type}")
            self.log_hsp(hsp)
            self.hsps.append(hsp)

    def create_blastn_hsp_object(self, hsp_element):
        """
        Create Hsp object for blastn hits from Hsp XML element
        """
        hit_start=int(hsp_element.find("Hsp_hit-from").text)
        hit_end=int(hsp_element.find("Hsp_hit-to").text)
        query_start=int(hsp_element.find("Hsp_query-from").text)
        query_end=int(hsp_element.find("Hsp_query-to").text)
        hseq = fetch_sequence(self.accession, hit_start, hit_end)

        return Hsp(
            bit_score=float(hsp_element.find("Hsp_bit-score").text),
            score=int(hsp_element.find("Hsp_score").text),
            evalue=float(hsp_element.find("Hsp_evalue").text),
            query_start=query_start,
            query_end=query_end,
            hit_start=hit_start,
            hit_end=hit_end,
            query_frame=int(hsp_element.find("Hsp_query-frame").text),
            hit_frame=int(hsp_element.find("Hsp_hit-frame").text),
            identity=int(hsp_element.find("Hsp_identity").text),
            positive=int(hsp_element.find("Hsp_positive").text),
            gaps=int(hsp_element.find("Hsp_gaps").text),
            align_len=int(hsp_element.find("Hsp_align-len").text),
            qseq=hsp_element.find("Hsp_qseq").text,
            hseq=hseq,
            midline=hsp_element.find("Hsp_midline").text,
        )

    def create_tblastn_hsp_object(self, hsp_element):
        hit_start=int(hsp_element.find("Hsp_hit-from").text)
        hit_end=int(hsp_element.find("Hsp_hit-to").text)
        query_start=3*int(hsp_element.find("Hsp_query-from").text)
        query_end=3*int(hsp_element.find("Hsp_query-to").text)
        hseq = fetch_sequence(self.accession, hit_start, hit_end)

        return Hsp(
            bit_score=float(hsp_element.find("Hsp_bit-score").text),
            score=int(hsp_element.find("Hsp_score").text),
            evalue=float(hsp_element.find("Hsp_evalue").text),
            query_start=query_start,
            query_end=query_end,
            hit_start=hit_start,
            hit_end=hit_end,
            query_frame=int(hsp_element.find("Hsp_query-frame").text),
            hit_frame=int(hsp_element.find("Hsp_hit-frame").text),
            identity=int(hsp_element.find("Hsp_identity").text),
            positive=int(hsp_element.find("Hsp_positive").text),
            gaps=int(hsp_element.find("Hsp_gaps").text),
            align_len=int(hsp_element.find("Hsp_align-len").text),
            qseq=hsp_element.find("Hsp_qseq").text,
            hseq=hseq,
            midline=hsp_element.find("Hsp_midline").text,
        )

    def create_tblastx_hsp_object(self, hsp_element):
        hit_start=int(hsp_element.find("Hsp_hit-from").text)
        hit_end=int(hsp_element.find("Hsp_hit-to").text)
        query_start=int(hsp_element.find("Hsp_query-from").text)
        query_end=int(hsp_element.find("Hsp_query-to").text)
        hseq = fetch_sequence(self.accession, hit_start, hit_end)
        return Hsp(
            bit_score=float(hsp_element.find("Hsp_bit-score").text),
            score=int(hsp_element.find("Hsp_score").text),
            evalue=float(hsp_element.find("Hsp_evalue").text),
            query_start=query_start,
            query_end=query_end,
            hit_start=hit_start,
            hit_end=hit_end,
            query_frame=int(hsp_element.find("Hsp_query-frame").text),
            hit_frame=int(hsp_element.find("Hsp_hit-frame").text),
            identity=int(hsp_element.find("Hsp_identity").text),
            positive=int(hsp_element.find("Hsp_positive").text),
            gaps=int(hsp_element.find("Hsp_gaps").text),
            align_len=int(hsp_element.find("Hsp_align-len").text),
            qseq=hsp_element.find("Hsp_qseq").text,
            hseq=hseq,
            midline=hsp_element.find("Hsp_midline").text,
        )
    
    
    def log_hsp(self, hsp_object):
        '''
        log start and stop of where query is hit for hsp object in csv format
        '''
        with open("hit_coverage.csv", "a") as fout:
            formatted_species = self.species.replace(" ","_")
            fout.write(f"{self.accession},{formatted_species},{hsp_object.query_start},{hsp_object.query_end}\n")





    def get_fasta_sequences(self, query_length):
        """
        Get FASTA-formatted sequences for each HSP's hit sequence
        """
        fasta_sequences = []
        for hsp in self.hsps:
            fasta_sequences.append(hsp.get_fasta_sequence(self.accession, self.species, query_length))
        return fasta_sequences


class Hsp:
    """
    Manage HSP information
    An HSP is a high scoring pair representing a specific aligngment within a BLAST hit of a sequence
    """
    def __init__(self, bit_score, score, evalue, query_start, query_end,
                 hit_start, hit_end, query_frame, hit_frame, identity,
                 positive, gaps, align_len, qseq, hseq, midline):
        self.bit_score = bit_score
        self.score = score
        self.evalue = evalue
        self.query_start = query_start
        self.query_end = query_end
        self.hit_start = hit_start
        self.hit_end = hit_end
        self.query_frame = query_frame
        self.hit_frame = hit_frame
        self.identity = identity
        self.positive = positive
        self.gaps = gaps
        self.align_len = align_len
        self.qseq = qseq
        self.hseq = hseq
        self.midline = midline
        if self.hit_start>self.hit_end or self.hit_frame<0: #get reverse complement
            self.hseq=reverse_complement(self.hseq)

    def get_fasta_sequence(self, hit_accession, species, query_length):
        """
        Get FASTA-formatted sequence for HSP hits positioned at where it hit in the query
        """
        delta_query = self.query_end - self.query_start + 1 
        print(f"hit sequence: {self.hseq} query start: {self.query_start}  query end: {self.query_end} Hit seq length: {len(self.hseq)} delta query: {delta_query}")

        if delta_query != len(self.hseq):
            print("error! query delta different from hseq length")


        seq = "-" * (self.query_start - 1) + self.hseq
        seq += ("-" * (query_length-len(seq)))


        print(f"formatted hit sequence: {seq} /n length: {len(seq)}")

        header = f">{hit_accession}"
        description = f"{self.hit_start}-{self.hit_end} | {species}"
        seq_record = SeqRecord(Seq(seq), id=hit_accession, description=description)
        return seq_record

def reverse_complement(dna_sequence):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', "-":"-"}
    reverse_sequence = reversed(dna_sequence)
    complement_sequence = [complement_dict[base] for base in reverse_sequence]
    return ''.join(complement_sequence)

def get_sequence_length(query_fasta):
    """
    Get the length of the sequence in a FASTA file.

    Parameters:
    - fasta_file (str): Path to the FASTA file.

    Returns:
    - int: Length of the sequence.
    """
    with open(query_fasta, "r") as file:
        lines = file.readlines()

    # Remove newline characters and check if the file is not empty
    lines = [line.strip() for line in lines if line.strip()]
    if not lines:
        raise ValueError("FASTA file is empty.")

    # Assume the first line is the header, and concatenate the rest as the sequence
    sequence = "".join(lines[1:])

    return len(sequence)

def fetch_sequence(accession_number, start, end):
    """
    Fetch a specified range of sequence from NCBI for a given accession number.

    Parameters:
    - accession_number (str): Accession number to fetch.
    - start (int): Start position of the sequence.
    - end (int): End position of the sequence.

    Returns:
    - sequence (str): Fetched sequence in FASTA format.
    """
    try:
        # Fetch the sequence from NCBI for the specified accession number and range
        print(f"fetching sequence for: {accession_number}")
        handle = Entrez.efetch(db="nuccore", id=accession_number, rettype="fasta", retmode="text", seq_start=start, seq_stop=end)
        record = SeqIO.read(handle, "fasta")
        handle.close()

        nucleotides = {"A", "T", "G", "C"}
        sequence = str(record.seq)

        # Replace characters not in nucleotides or '-' with '-'
        sanitized_sequence = ''.join(c if c in nucleotides or c == '-' else '-' for c in sequence)

        return sanitized_sequence

    except Exception as e:
        print(f"An error occurred: {e}")
        return None

def fetch_query(query_fasta, start, end):
    """
    Fetch a specified range from a query sequence in a FASTA file.

    Parameters:
    - query_fasta (str): Path to the FASTA file containing the query sequence.
    - start (int): Start position of the desired range.
    - end (int): End position of the desired range.

    Returns:
    - sequence (str): Fetched sequence in the specified range.
    """
    try:
        # Read the query sequence from the FASTA file
        with open(query_fasta, "r") as file:
            records = list(SeqIO.parse(file, "fasta"))

        if not records:
            raise ValueError("No sequences found in the FASTA file.")

        # Assume there is only one sequence in the file
        query_sequence = records[0].seq

        # Fetch the specified range from the query sequence
        sequence_range = query_sequence[start - 1:end]

        return str(sequence_range)

    except Exception as e:
        print(f"An error occurred: {e}")
        return None

def get_ranges(record_list):
    """
    return str the ranges of each record in the list where the range is the first string
    after the seqID
    """
    ranges=[]
    for record in record_list:
        seq_range = record.description.split()[1]
        ranges.append(seq_range)
    return ":".join(ranges)

def pad_seqs(hit_seqs):
    max_length = max(len(seq) for seq in hit_seqs)  # Find the maximum length among all strings
    padded_seqs = [seq + "--" * (max_length - len(seq)) for seq in hit_seqs]
    return padded_seqs



def get_hit_consensus(hits_fasta, query_fasta, blast_type, db):
    '''
    Create a consensus sequence with positions that agree in hits_fasta alignment, positions that disagree are used from the query provided that the nucleotide is found in at least one of the hit seqs at that position
    '''
    query_record = SeqIO.read(query_fasta, "fasta")
    hit_records = list(SeqIO.parse(hits_fasta, "fasta"))

    query_seq = str(query_record.seq)
    hit_seqs = [str(record.seq) for record in hit_records]
    hit_seqs = pad_seqs(hit_seqs)
    
    hit_ranges=get_ranges(hit_records)
    hit_species=" ".join(hit_records[0].description.split()[-2:])

    if len(hit_seqs) <= 1:
        hit_description=f" | {hit_ranges} | {hit_species} {blast_type} {db}"
    if len(hit_seqs) > 1:
        hit_description = f" | {hit_ranges} | {hit_species} {blast_type} {db} consensus"
    print(f"num of seqs: {len(hit_seqs)} hit description: {hit_description}")
    hit_consensus = ""
    for i in range(len(hit_seqs[0])):
        # make a list of the nucleotide at position i for each hit sequence
        nts = [hit[i] for hit in hit_seqs if hit[i]!= "-"]

        if not nts:  # all "-" at this position
            hit_consensus += "-"
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

    hit_consensus_record = SeqRecord(Seq(hit_consensus), id=hit_records[0].id, description=hit_description)
    return hit_consensus_record


def get_species_by_assembly_ids(assembly_ids):
    with open("assemblies_searched.txt", "a") as output_file:
        # Iterate through each assembly ID in the list
        for assembly_id in assembly_ids:
            # Use the efetch utility to retrieve information about the assembly
            print(assembly_id)


            try:
                handle = Entrez.efetch(db="nucleotide", id=assembly_id, retmode="xml")
            except:
                assembly_id += "0"
                handle = Entrez.efetch(db="nucleotide", id=assembly_id, retmode="xml")
            

            record = handle.read()
            print(record)

            # Parse the XML response
            tree = ET.fromstring(record)

            # Extract species information
            species_element = tree.find(".//GBSeq_organism")
            species = species_element.text if species_element is not None else 'Species information not available'


            # Write the assembly ID and species to the file
            output_file.write(f"{assembly_id}, {species}\n")
            print(f"Assembly ID {assembly_id}: {species}")
        output_file.write(f"\n")

def extract_assembly_ids_from_string(BlastOutput_db):
    # Split by space
    components = BlastOutput_db.split()

    # Extract assembly IDs from each component that contains "WGS_VDB://"
    assembly_ids = [component.split("WGS_VDB://")[1] for component in components if "WGS_VDB://" in component]

    #Format ids for retrival
    ids = []
    for assembly_id in assembly_ids:
        id = assembly_id[:-2]
        id+="00000000"
        ids.append(id)

    return ids

def parse_blast(tree_root, query_fasta):
    """
    Creates a nucleotide alignment of blast hits by positioning sequences where they hit the query in the blast search
    """
    blast_type = tree_root.find(".//BlastOutput_program").text

    # Determine database type searched
    if "TSA" in tree_root.find(".//Hit_def").text:
        db = "TSA"
    else:
        db = "WGS"

    #Get assembly accesion and species searched
    BlastOutput_db = tree_root.find(".//BlastOutput_db").text
    assembly_ids=extract_assembly_ids_from_string(BlastOutput_db)
    get_species_by_assembly_ids(assembly_ids)

    output_alignment = f"{query_fasta.rsplit('.')[0]}_{blast_type}_{db}.fasta"
    print(f"output alignment: {output_alignment}")

    #Process blast hits
    blast_hits = []
    for hit_element in tree_root.findall(".//Hit"):
        blast_hit = BlastHit(hit_element, blast_type)
        blast_hits.append(blast_hit)

    # Write FASTA sequences to the output alignment file
    sequences = []
    query_length = get_sequence_length(query_fasta)
    for blast_hit in blast_hits:
        hit_fout=blast_hit.accession+"_hits.fasta"
        fasta_sequences = blast_hit.get_fasta_sequences(query_length)
        SeqIO.write(fasta_sequences, hit_fout, "fasta")
        consensus=get_hit_consensus(hit_fout, query_fasta, db, blast_type)
        sequences.append(consensus)
    print(f"Writing consensus sequences to {output_alignment}")
    SeqIO.write(sequences, output_alignment,"fasta")



def main(xml_file, query_fasta):
    # Parse XML file
    tree = ET.parse(xml_file)
    root = tree.getroot()


    # Process blast hits
    parse_blast(root, query_fasta)




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process BLAST XML results and create nucleotide alignments.")
    parser.add_argument("xml_file", help="Path to the BLAST XML file")
    parser.add_argument("query_fasta", help="Path to the query FASTA file")
    parser.add_argument("email", help="email address to use NCBI Entrez")
    args = parser.parse_args()
    Entrez.email = args.email
    main(args.xml_file, args.query_fasta)


